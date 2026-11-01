"""The ratio harness: interleaved HEAD-vs-reference CPU pairs per workload, written as a ``perf-ratio/1`` session.

Two arm worktrees (the pinned reference and HEAD) are proved, rebuilt and
settled, then every workload runs in alternating pairs (AB, BA, AB, ...; A =
reference), each arm child in a fresh process of its own interpreter, with the
1-minute load read around every pair. Both arms are bracketed around the whole
session, so a tree, library or native build that moves mid-session is refused.
Load is recorded, never gated, and no ratio is gated here.

Usage: ``python -m tools.perf_ratio [--set NAME | --workload NAME ...] [--pairs N] [--inject-cpu-seconds ARM=S]``.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import Path
from typing import Callable
from typing import Optional

from tools.perf_ratio_arms import DEFAULT_HEAD_WORKTREE
from tools.perf_ratio_arms import DEFAULT_REF_WORKTREE
from tools.perf_ratio_arms import PINNED_REF_COMMIT
from tools.perf_ratio_arms import ROOT
from tools.perf_ratio_arms import Arm
from tools.perf_ratio_arms import ArmBracket
from tools.perf_ratio_arms import arm_python
from tools.perf_ratio_arms import checkout_head
from tools.perf_ratio_arms import live_bracket
from tools.perf_ratio_arms import open_arm
from tools.perf_ratio_arms import rebuild_arm
from tools.perf_ratio_arms import require_distinct_arms
from tools.perf_ratio_arms import require_settled
from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import LoadOneMinute
from tools.perf_ratio_record import PairRun
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import WorkloadRatio
from tools.perf_ratio_workloads import DEFAULT_SET_NAME
from tools.perf_ratio_workloads import WORKLOAD_SETS
from tools.perf_ratio_workloads import ChildRun
from tools.perf_ratio_workloads import Family
from tools.perf_ratio_workloads import InjectionUnsupportedError
from tools.perf_ratio_workloads import Workload
from tools.perf_ratio_workloads import child_command
from tools.perf_ratio_workloads import parse_child
from tools.perf_ratio_workloads import workload_named
from tools.pipeline_timing_1_1_base import editable_build_dir_of
from tools.pipeline_timing_1_1_base import require_measured_tree
from tools.pipeline_timing_1_1_base import run_in_tree
from tools.run_native_identity import track_native_in
from tools.run_provenance import file_set_sha256
from tools.run_provenance import machine_identity
from tools.run_provenance import utc_now

# Pairs per workload (spec: N = 5 pairs, so min and median each rest on five runs per arm).
DEFAULT_PAIRS = 5
PERF_RATIO_SESSION_PATH = Path("benchmarks") / "results" / "perf_ratio_session.json"
# Every source whose code shapes a ratio number: the harness itself, the reused tree/rebuild proofs and
# provenance readers, and the native probe's whole import closure from this tree (``PROBE_SNIPPET`` puts
# this tree's ``tools`` and ``benchmarks`` first on each arm's sys.path; a test measures that closure).
HARNESS_SOURCES = (
    Path("tools") / "perf_ratio.py",
    Path("tools") / "perf_ratio_record.py",
    Path("tools") / "perf_ratio_workloads.py",
    Path("tools") / "perf_ratio_arms.py",
    Path("tools") / "perf_ratio_history.py",
    Path("tools") / "pipeline_timing_1_1_base.py",
    Path("tools") / "__init__.py",
    Path("tools") / "run_provenance.py",
    Path("tools") / "run_native_identity.py",
    Path("tools") / "held_native_identity.py",
    Path("tools") / "held_source_closure.py",
    Path("benchmarks") / "__init__.py",
    Path("benchmarks") / "errors.py",
    Path("benchmarks") / "performance_record.py",
)
# Separates the arm from the seconds in ``--inject-cpu-seconds ARM=S``.
INJECTION_SEPARATOR = "="
# Injected CPU of an arm the injection does not name: no spin.
NO_INJECTION = CpuSeconds(0.0)

# Runs one workload once in one arm's fresh child and returns its measurement.
Launch = Callable[[ArmName, Workload], ChildRun]
# Reads the 1-minute load average.
LoadReader1 = Callable[[], LoadOneMinute]


class NondeterministicArmError(Exception):
    """Two runs of one workload within one arm report different output digests."""


class InvalidInjectionError(Exception):
    """``--inject-cpu-seconds`` is not ``ARM=S`` with ARM in {ref, head} and S finite and positive."""


class HeadWorktreeIsSharedError(Exception):
    """``--head-worktree`` resolves to the orchestrator's tree or the ref worktree, which a HEAD checkout would move."""


class HarnessSourceMissingError(Exception):
    """A file of `HARNESS_SOURCES` is missing from this tree, so the harness digest cannot be read."""


class HarnessChangedDuringRunError(Exception):
    """The harness sources' digest at the end of the session differs from the one read before anything ran."""


class NoPairsError(Exception):
    """The session asks for fewer than one pair per workload, so no min, median or ratio is defined."""


class DuplicateWorkloadError(Exception):
    """A workload is named twice in one session, which a ``perf-ratio/1`` record refuses."""


def harness_sha256() -> str:
    """SHA-256 over `HARNESS_SOURCES` in this tree.

    Raises:
        HarnessSourceMissingError: A source file does not exist (each missing one is named).
    """
    missing = [path.as_posix() for path in HARNESS_SOURCES if not (ROOT / path).is_file()]
    if missing:
        raise HarnessSourceMissingError(f"the harness digest covers {', '.join(missing)}, missing from {ROOT}.")
    return file_set_sha256(ROOT, HARNESS_SOURCES)


def require_unchanged_harness(before: str) -> None:
    """Refuse a harness whose sources changed since *before* was read.

    The orchestrator imports the harness at start, and every arm's native probe
    imports this tree's probe code at settle, bracket start and bracket finish,
    so an edit anywhere in that span would give a digest of code that did not run.

    Raises:
        HarnessChangedDuringRunError: The digest now differs from *before*.
        HarnessSourceMissingError: A source file was removed.
    """
    after = harness_sha256()
    if after != before:
        raise HarnessChangedDuringRunError(f"the harness sources changed during the session: digest {before} at start, {after} at end.")


def load_one_minute() -> LoadOneMinute:
    """The 1-minute load average of this machine."""
    return LoadOneMinute(os.getloadavg()[0])


def pair_order(index: int) -> tuple[ArmName, ArmName]:
    """The arm order of pair *index* (0-based): AB, BA, AB, ... with A the reference."""
    return (ArmName.HEAD, ArmName.REF) if index % 2 else (ArmName.REF, ArmName.HEAD)


def require_plan(workloads: Sequence[Workload], pairs: int) -> None:
    """Refuse a session plan whose record could not be built.

    Raises:
        NoPairsError: *pairs* is below one.
        DuplicateWorkloadError: A workload appears more than once in *workloads*.
    """
    if pairs < 1:
        raise NoPairsError(f"a session runs at least one pair per workload, not {pairs}.")
    names = [workload.name for workload in workloads]
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise DuplicateWorkloadError(f"each workload is measured once per session; named more than once: {', '.join(repeated)}.")


def measure_workload(workload: Workload, *, pairs: int, launch: Launch, load_reader: LoadReader1) -> WorkloadRatio:
    """Run *workload* in *pairs* interleaved pairs, reading the load before and after every pair.

    Args:
        workload: The workload.
        pairs: How many pairs; pair *i* runs the arms in `pair_order` (i).
        launch: Runs one arm's child once.
        load_reader: Reads the 1-minute load.

    Returns:
        The workload's pairs and each arm's (single) output digest.

    Raises:
        NoPairsError: *pairs* is below one.
        NondeterministicArmError: A run's digest differs from its arm's first run (refused at that run).
        NonPositiveCpuError: A run reports non-positive CPU.
    """
    require_plan((workload,), pairs)
    digests: dict[ArmName, str] = {}
    measured: list[PairRun] = []
    for index in range(pairs):
        order = pair_order(index)
        load_before = load_reader()
        runs: dict[ArmName, ChildRun] = {}
        for arm in order:
            run = launch(arm, workload)
            first = digests.setdefault(arm, run.digest)
            if run.digest != first:
                raise NondeterministicArmError(f"{workload.name}: the {arm.value} arm's pair {index} digest {run.digest!r} differs from its first run's {first!r}.")
            runs[arm] = run
        load_after = load_reader()
        measured.append(PairRun.build(order=order, ref_cpu=runs[ArmName.REF].cpu_seconds, head_cpu=runs[ArmName.HEAD].cpu_seconds, load_before=load_before, load_after=load_after))
    return WorkloadRatio(
        workload=workload.name,
        family=workload.family.value,
        pairs=tuple(measured),
        ref_digest=digests[ArmName.REF],
        head_digest=digests[ArmName.HEAD],
    )


def measure_session(
    workloads: Sequence[Workload],
    *,
    ref_commit: str,
    brackets: Mapping[ArmName, ArmBracket],
    pairs: int,
    launch: Launch,
    load_reader: LoadReader1,
    injected: Optional[Injection],
    harness: str,
) -> PerfRatioSession:
    """Measure every workload inside both arms' brackets.

    Order: the start time, start both brackets (ref, head), measure the
    workloads, finish both brackets (ref, head), then read the end time and
    `machine_identity`.

    Args:
        workloads: The workloads, measured in this order.
        ref_commit: Full commit id of the pinned reference.
        brackets: Each arm's bracket (`live_bracket` for a real arm).
        pairs: Pairs per workload.
        launch: Runs one arm's child once.
        load_reader: Reads the 1-minute load.
        injected: The audit's injected CPU, None for an honest session.
        harness: Digest of the harness sources (`harness_sha256` for a real session).

    Returns:
        The unchecked session; `PerfRatioSession.build` checks it.

    Raises:
        NoPairsError: *pairs* is below one (before any bracket starts).
        DuplicateWorkloadError: A workload is listed twice (before any bracket starts).
        NondeterministicArmError: Runs within one arm disagree on their digest.
    """
    require_plan(workloads, pairs)
    started_utc = utc_now()
    finishers = {arm: brackets[arm]() for arm in ArmName}
    measured = tuple(measure_workload(workload, pairs=pairs, launch=launch, load_reader=load_reader) for workload in workloads)
    arms: dict[ArmName, ArmIdentity] = {arm: finishers[arm]() for arm in ArmName}
    ended_utc = utc_now()
    host, cpu_model = machine_identity()
    return PerfRatioSession(
        ref_commit=ref_commit,
        arms=arms,
        workloads=measured,
        host=host,
        cpu_model=cpu_model,
        started_utc=started_utc,
        ended_utc=ended_utc,
        harness_sha256=harness,
        injected=injected,
    )


def arm_launcher(arms: Mapping[ArmName, Arm], injected: Optional[Injection]) -> Launch:
    """The launch that runs a workload's child in the named arm's tree, interpreter and pinned environment.

    Only the arm *injected* names receives its spin; the other arm's child spins nothing.
    """

    def inject_for(arm: Arm) -> CpuSeconds:
        return injected.seconds if injected is not None and injected.arm is arm.name else NO_INJECTION

    def launch(name: ArmName, workload: Workload) -> ChildRun:
        arm = arms[name]
        return parse_child(workload, run_in_tree(arm.root, arm.env)(child_command(workload, arm.python, inject_for(arm))))

    return launch


def parse_injection(text: str) -> Injection:
    """The injection ``ARM=S`` names: ARM one of ``ref`` / ``head``, S finite and positive CPU seconds.

    Raises:
        InvalidInjectionError: Any other form.
    """
    arm_text, separator, seconds_text = text.partition(INJECTION_SEPARATOR)
    arms = ", ".join(arm.value for arm in ArmName)
    if not separator:
        raise InvalidInjectionError(f"--inject-cpu-seconds is ARM=S, not {text!r}.")
    try:
        arm = ArmName(arm_text)
    except ValueError as error:
        raise InvalidInjectionError(f"--inject-cpu-seconds names an arm ({arms}), not {arm_text!r}.") from error
    try:
        seconds = float(seconds_text)
    except ValueError as error:
        raise InvalidInjectionError(f"--inject-cpu-seconds seconds are a number, not {seconds_text!r}.") from error
    if not (math.isfinite(seconds) and seconds > 0.0):
        raise InvalidInjectionError(f"--inject-cpu-seconds seconds are finite and positive, not {seconds_text!r}.")
    return Injection(arm=arm, seconds=CpuSeconds(seconds))


def require_injectable(workloads: Sequence[Workload], injected: Optional[Injection]) -> None:
    """Refuse an injection over any workload whose clock the harness does not own (every family but GT).

    Raises:
        InjectionUnsupportedError: *injected* is set and a workload is not GT.
    """
    if injected is None:
        return
    unsupported = [workload.name for workload in workloads if workload.family is not Family.GT]
    if unsupported:
        raise InjectionUnsupportedError(f"CPU injection needs a harness-owned clock, which only the GT family has; refused for {', '.join(unsupported)}.")


def require_head_worktree_apart(head_root: Path, ref_root: Path) -> None:
    """Refuse a HEAD worktree that is the orchestrator's tree or the ref worktree (compared resolved).

    Raises:
        HeadWorktreeIsSharedError: *head_root* resolves to `ROOT` or to *ref_root*.
    """
    head = head_root.resolve()
    for label, other in (("the orchestrator's tree", ROOT), ("the ref worktree", ref_root)):
        if head == other.resolve():
            raise HeadWorktreeIsSharedError(f"--head-worktree {head_root} is {label} ({other.resolve()}); checking out HEAD there would move it.")


def summary_lines(session: PerfRatioSession) -> list[str]:
    """One line per workload (ratio, pair spread, same_output, load range), then the session's wall time."""
    lines = []
    for workload in session.workloads:
        ratios = [pair.ratio for pair in workload.pairs]
        low, high = workload.load_range
        lines.append(
            f"{workload.workload}: ratio {workload.ratio:.3f}x (pair median {workload.pair_ratio_median:.3f}x, "
            f"pairs {min(ratios):.3f}-{max(ratios):.3f}), same_output {workload.same_output}, load {low:.1f}-{high:.1f}"
        )
    lines.append(f"wall {session.wall_seconds:.0f} s")
    return lines


def write_session(session: PerfRatioSession, path: Path) -> PerfRatioSession:
    """Check *session* by `PerfRatioSession.build` on its JSON, write it to *path*, and return the checked session.

    Raises:
        MalformedRatioRecordError: The session is malformed.
        InconsistentRatioRecordError: A derived value disagrees, or ``ref_commit`` is not the ref arm's rev.
    """
    checked = PerfRatioSession.build(session.to_json())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(checked.dumps(), encoding="utf-8")
    return checked


def _orchestrator_head() -> str:
    return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m tools.perf_ratio", description="Measure HEAD-over-reference CPU ratios in interleaved pairs.")
    parser.add_argument("--ref-worktree", type=Path, default=DEFAULT_REF_WORKTREE, help="the clean reference worktree")
    parser.add_argument("--ref-rev", default=PINNED_REF_COMMIT, help="the full commit id the reference worktree must be at")
    parser.add_argument("--head-worktree", type=Path, default=DEFAULT_HEAD_WORKTREE, help="the clean worktree checked out to --head-rev")
    parser.add_argument("--head-rev", default=None, help="the HEAD arm's revision (default: this tree's HEAD)")
    parser.add_argument("--set", choices=sorted(WORKLOAD_SETS), default=DEFAULT_SET_NAME, help="the workload set")
    parser.add_argument("--workload", action="append", default=None, metavar="NAME", help="a workload to measure (repeatable; replaces --set)")
    parser.add_argument("--pairs", type=int, default=DEFAULT_PAIRS, help="interleaved pairs per workload")
    parser.add_argument("--inject-cpu-seconds", default=None, metavar="ARM=S", help="audit: spin S CPU seconds inside ARM's measured span (GT only)")
    parser.add_argument("--write", type=Path, default=PERF_RATIO_SESSION_PATH, help="where to write the session record")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry: measure both arms and write the ``perf-ratio/1`` session.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    args = _parser().parse_args(argv)
    workloads = tuple(workload_named(name) for name in args.workload) if args.workload else WORKLOAD_SETS[args.set]
    ref_root: Path = args.ref_worktree
    head_root: Path = args.head_worktree

    # 1. The injection is checked before anything else; then the refusals that need no git.
    injected = None if args.inject_cpu_seconds is None else parse_injection(args.inject_cpu_seconds)
    require_injectable(workloads, injected)
    require_plan(workloads, args.pairs)
    require_head_worktree_apart(head_root, ref_root)
    # The harness digest is read before anything runs; step 5 records it and the end re-reads it.
    harness = harness_sha256()

    # 2. Move the HEAD worktree to the HEAD revision.
    head_rev = checkout_head(head_root, args.head_rev if args.head_rev is not None else _orchestrator_head())

    # 3. Prove both arms, then that they are distinct.
    ref = open_arm(ArmName.REF, ref_root, arm_python(ref_root), args.ref_rev)
    head = open_arm(ArmName.HEAD, head_root, arm_python(head_root), head_rev)
    require_distinct_arms(ref, head)
    arms = {ArmName.REF: ref, ArmName.HEAD: head}

    # 4. Rebuild both arms (pinning this tree's build so it is never rebuilt), then settle both native builds.
    orchestrator_build_dir = editable_build_dir_of((sys.executable,), ROOT, os.environ)
    for arm in arms.values():
        print(f"rebuilt the {arm.name.value} arm: {' '.join(rebuild_arm(arm, orchestrator_build_dir))}")
    for arm in arms.values():
        require_settled(track_native_in(arm.python, arm.root, arm.env))

    # 5. Measure inside both arms' brackets.
    session = measure_session(
        workloads,
        ref_commit=ref.rev,
        brackets={name: live_bracket(arm) for name, arm in arms.items()},
        pairs=args.pairs,
        launch=arm_launcher(arms, injected),
        load_reader=load_one_minute,
        injected=injected,
        harness=harness,
    )
    require_unchanged_harness(harness)

    # 6. Both trees are still their proved clean revisions.
    for arm in arms.values():
        require_measured_tree(arm.root, arm.rev)

    # 7. Check, write and summarise.
    write = args.write if args.write.is_absolute() else ROOT / args.write
    checked = write_session(session, write)
    print(f"wrote {write}")
    for line in summary_lines(checked):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
