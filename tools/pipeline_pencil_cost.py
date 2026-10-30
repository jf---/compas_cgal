"""Measure the tangent pencil's own decision cost `C_pencil` per ``pipeline`` gate case (Task 1.2, ruling D6).

The Task 1.2 time gate holds each ``pipeline`` case's generation median to the
slowest 1.1 generator run plus `C_pencil`, the pencil's decision cost measured
in the same session as `tools.pipeline_timing`. This tool measures it the way
`tools.pipeline_timing` measures generation time: every ``pipeline`` gate case
not excluded by the admissions, one case at a time, each of `RUNS_PER_CASE`
runs in a fresh interpreter (``python -m tools.pipeline_pencil_cost --case
KEY``). In the child, `pencil_instrumentation` wraps each of the pencil's
decisions (`TangentPencil.build`, `RampOut.build`, `RampOut.next_member`) and
the first cut (`generate.enter_first_cut`) in its own pass-through wall clock
and call counter.

The wrappers decide nothing: arguments, return values and exceptions pass
through untouched. That is checked, not assumed: every run's ops digest must
equal the uninstrumented digest of the same case in
`benchmarks/results/pipeline_timing.json`, and the runs must agree with each
other, or the tool fails loud. What `C_pencil` contains is checked too: each
callable must be reached, and the counts must match the first cut the
generator reports (one build of each, one first cut, one ``next_member`` per
emitted member plus the root's), so a wrapper that is dropped, renamed or
re-routed fails loud instead of under-reading `C_pencil`.

The 1/5/15-minute load average is read before each case and after it. A
measurement refuses to start or finish a case above the quiet limit
(`tools.held_generation_clock.QUIET_LOAD_AVERAGE_MAX`, raising its
`MachineNotQuietError` with this tool's own hint); ``--indicative`` measures
at any load and marks the record indicative, which yields no time-gate verdict.

Run as ``pixi run pipeline-pencil-cost`` right after ``pixi run
pipeline-timing``; it writes `benchmarks/results/pipeline_pencil_cost.json`
and its rendering `benchmarks/results/pipeline_pencil_cost.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Iterator
from typing import Mapping
from typing import Sequence

from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket
from compas_cgal.pipeline import generate
from compas_cgal.pipeline.pencil import RampOut
from compas_cgal.pipeline.pencil import TangentPencil
from compas_cgal.pipeline.result import PipelineResult
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import LoadAverage
from tools.held_generation_clock import MachineNotQuietError
from tools.pipeline_pencil_cost_results import DECIDE_CALLABLES
from tools.pipeline_pencil_cost_results import ENTER_FIRST_CUT
from tools.pipeline_pencil_cost_results import PENCIL_GENERATOR
from tools.pipeline_pencil_cost_results import PIPELINE_PENCIL_COST_PATH
from tools.pipeline_pencil_cost_results import PIPELINE_PENCIL_COST_TABLE_PATH
from tools.pipeline_pencil_cost_results import RAMP_OUT_BUILD
from tools.pipeline_pencil_cost_results import RAMP_OUT_NEXT_MEMBER
from tools.pipeline_pencil_cost_results import TANGENT_PENCIL_BUILD
from tools.pipeline_pencil_cost_results import WRAPPED_CALLABLES
from tools.pipeline_pencil_cost_results import CasePencilCost
from tools.pipeline_pencil_cost_results import PencilRun
from tools.pipeline_pencil_cost_results import PipelinePencilCosts
from tools.pipeline_pencil_cost_results import is_pencil_key
from tools.pipeline_pencil_cost_results import pencil_cases
from tools.pipeline_pencil_cost_results import render_pipeline_pencil_cost_table
from tools.pipeline_pencil_cost_results import require_call_composition
from tools.pipeline_pencil_cost_results import rooted
from tools.pipeline_timing import CaseProcessFailedError
from tools.pipeline_timing import MissingCaseResultError
from tools.pipeline_timing import NondeterministicCaseError
from tools.pipeline_timing import gate_timing_keys
from tools.pipeline_timing import ops_digest
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import STORED_SECONDS_DECIMALS
from tools.pipeline_timing_results import PipelineTimings
from tools.pipeline_timing_results import format_load
from tools.pipeline_timing_results import load_pipeline_timings
from tools.pipeline_timing_results import stored_load
from tools.run_native_identity import NativeIdentity
from tools.run_native_identity import NativeTracker
from tools.run_native_identity import track_native_here
from tools.run_provenance import RunProvenance
from tools.run_provenance import TreeIdentity
from tools.run_provenance import TreeReader
from tools.run_provenance import file_set_sha256
from tools.run_provenance import interpreter_identity_here
from tools.run_provenance import machine_identity
from tools.run_provenance import require_unchanged_tree
from tools.run_provenance import sha256_file
from tools.run_provenance import tree_identity
from tools.run_provenance import utc_now

ROOT = Path(__file__).resolve().parents[1]
# The child prints its result on one stdout line after this prefix, so library
# output on stdout cannot be mistaken for it.
RESULT_LINE_PREFIX = "PIPELINE_PENCIL_COST_RESULT "
# What to do instead when the machine is loaded, named in the refusal.
INDICATIVE_HINT = (
    "Re-run `pixi run pipeline-timing` and `pixi run pipeline-pencil-cost` back to back on a quiet machine, or measure indicatively "
    "(no time-gate verdict): python -m tools.pipeline_pencil_cost --indicative."
)

LoadReader = Callable[[], LoadAverage]
# Builds a run's provenance from its UTC start and end.
Describer = Callable[[TreeIdentity, NativeIdentity, str, str], RunProvenance]
# The sources of this tool, digested into the provenance.
PENCIL_TOOL_SOURCES = (
    Path("tools") / "pipeline_pencil_cost.py",
    Path("tools") / "pipeline_pencil_cost_results.py",
    Path("tools") / "run_provenance.py",
    Path("tools") / "run_native_identity.py",
)
PIXI_LOCK = Path("pixi.lock")


class UnknownPencilCaseError(Exception):
    """A case key that is not a ``pipeline`` gate case."""


class InstrumentationNotReachedError(Exception):
    """A wrapped callable never ran: the pencil's decisions or the first cut moved out from under the wrappers."""


class NotAPipelineResultError(Exception):
    """The ``pipeline`` gate generator returned something other than a `PipelineResult`, so its ramp-out cannot be read."""


class PencilDigestMismatchError(Exception):
    """An instrumented run's ops digest differs from the uninstrumented digest of the same case."""


class MissingTimingReferenceError(Exception):
    """The timing document holds no uninstrumented digest for a case to measure."""


@dataclass(frozen=True)
class PencilCaseRun:
    """One instrumented generation of one case, as its child process reports it.

    Attributes:
        run: The three wall clocks.
        calls: Wrapped callable (`WRAPPED_CALLABLES`) to its call count.
        ramp_out_members: Pencil members the first cut emitted.
        operations: The number of operations generated.
        digest: The ops digest (`tools.pipeline_timing.ops_digest`).
    """

    run: PencilRun
    calls: Mapping[str, int]
    ramp_out_members: int
    operations: int
    digest: str


Launcher = Callable[[str], PencilCaseRun]


@dataclass
class CallClock:
    """Accumulated wall-clock seconds and call count of one wrapped callable."""

    seconds: float = 0.0
    calls: int = 0


@dataclass
class PencilClocks:
    """One clock per wrapped callable, filled by one instrumented generation."""

    clocks: dict[str, CallClock] = field(default_factory=lambda: {name: CallClock() for name in WRAPPED_CALLABLES})

    @property
    def decide_s(self) -> float:
        """`C_pencil`: the summed seconds of the pencil's decisions."""
        return sum(self.clocks[name].seconds for name in DECIDE_CALLABLES)

    @property
    def first_cut_s(self) -> float:
        """The summed seconds of the first cut."""
        return self.clocks[ENTER_FIRST_CUT].seconds

    @property
    def calls(self) -> dict[str, int]:
        """Wrapped callable to its call count."""
        return {name: clock.calls for name, clock in self.clocks.items()}


def _timed(function: Callable[..., Any], clock: CallClock) -> Callable[..., Any]:
    """*function* unchanged in arguments, result and exceptions, its wall clock added to *clock*."""

    def wrapper(*args: Any, **kwargs: Any) -> Any:
        started = time.perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            clock.seconds += time.perf_counter() - started
            clock.calls += 1

    return wrapper


def _wrapped(original: Any, clock: CallClock) -> Any:
    """The timed replacement of a class or module attribute, keeping a classmethod a classmethod."""
    if isinstance(original, classmethod):
        return classmethod(_timed(original.__func__, clock))
    return _timed(original, clock)


@contextmanager
def pencil_instrumentation() -> Iterator[PencilClocks]:
    """Wrap each of the pencil's decisions and the first cut in its own pass-through clock for the duration of the block.

    The class attributes and the module attribute are restored on exit,
    whether the block raises or not. No decision nests in another
    (`TangentPencil.build` calls only the unwrapped `TangentPencil.on_grid`),
    so `PencilClocks.decide_s` does not double-count.

    Yields:
        The clocks the block's generation fills.
    """
    clocks = PencilClocks()
    targets: list[tuple[object, str, str]] = [
        (TangentPencil, "build", TANGENT_PENCIL_BUILD),
        (RampOut, "build", RAMP_OUT_BUILD),
        (RampOut, "next_member", RAMP_OUT_NEXT_MEMBER),
        (generate, "enter_first_cut", ENTER_FIRST_CUT),
    ]
    originals = [(owner, attribute, vars(owner)[attribute]) for owner, attribute, _name in targets]
    try:
        for (owner, attribute, original), (_owner, _attribute, name) in zip(originals, targets):
            setattr(owner, attribute, _wrapped(original, clocks.clocks[name]))
        yield clocks
    finally:
        for owner, attribute, original in originals:
            setattr(owner, attribute, original)


def run_case_here(key: str) -> PencilCaseRun:
    """Generate *key*'s case once in this process under `pencil_instrumentation`.

    Args:
        key: A ``pipeline`` gate case key.

    Returns:
        The run.

    Raises:
        UnknownPencilCaseError: *key* is not a ``pipeline`` gate case.
        NotAPipelineResultError: The generator's result carries no ramp-out record.
        InstrumentationNotReachedError: The generation never reached a wrapped callable.
        CallCompositionError: The call counts do not match the first cut the generator reports.
    """
    cases = pencil_cases()
    if key not in cases:
        raise UnknownPencilCaseError(f"{key!r} is not a {PENCIL_GENERATOR} gate case.")
    cap, pocket = cases[key]
    spec = gate_pocket(pocket, tea_cap_deg=GateCapDegrees(cap))
    with pencil_instrumentation() as clocks:
        started = time.perf_counter()
        result = GATE_GENERATORS[PENCIL_GENERATOR](spec)
        elapsed = time.perf_counter() - started
    if not isinstance(result, PipelineResult):
        raise NotAPipelineResultError(f"{key}: the {PENCIL_GENERATOR} generator returned {type(result).__name__}, not PipelineResult.")
    calls = clocks.calls
    unreached = [name for name, count in calls.items() if count == 0]
    if unreached:
        raise InstrumentationNotReachedError(f"{key}: {unreached} never ran; every wrapped callable must be reached ({calls}).")
    members = len(result.ramp_out.members)
    require_call_composition(key, calls, members)
    return PencilCaseRun(
        run=PencilRun(generation_s=elapsed, decide_s=clocks.decide_s, first_cut_s=clocks.first_cut_s),
        calls=calls,
        ramp_out_members=members,
        operations=len(result.operations),
        digest=ops_digest(result),
    )


def result_payload(run: PencilCaseRun) -> dict[str, object]:
    """The JSON object a child prints after `RESULT_LINE_PREFIX` for *run*."""
    return {
        "generation_s": run.run.generation_s,
        "decide_s": run.run.decide_s,
        "first_cut_s": run.run.first_cut_s,
        "calls": dict(run.calls),
        "ramp_out_members": run.ramp_out_members,
        "operations": run.operations,
        "digest": run.digest,
    }


def launch_fresh_process(key: str) -> PencilCaseRun:
    """Run *key*'s case once in a fresh interpreter and read its result line.

    The child inherits this process's environment (the editable build the pixi
    task selects) and runs from the repository root.

    Args:
        key: A ``pipeline`` gate case key.

    Returns:
        The child's run.

    Raises:
        CaseProcessFailedError: The child exited non-zero; its stderr is quoted.
        MissingCaseResultError: The child printed no result line.
    """
    command = [sys.executable, "-m", "tools.pipeline_pencil_cost", "--case", key]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise CaseProcessFailedError(f"{key}: {' '.join(command)} exited {completed.returncode}:\n{completed.stderr}")
    lines = [line for line in completed.stdout.splitlines() if line.startswith(RESULT_LINE_PREFIX)]
    if len(lines) != 1:
        raise MissingCaseResultError(f"{key}: the child printed {len(lines)} result lines, not one.")
    payload = json.loads(lines[0][len(RESULT_LINE_PREFIX) :])
    return PencilCaseRun(
        run=PencilRun(generation_s=float(payload["generation_s"]), decide_s=float(payload["decide_s"]), first_cut_s=float(payload["first_cut_s"])),
        calls={str(name): int(count) for name, count in payload["calls"].items()},
        ramp_out_members=int(payload["ramp_out_members"]),
        operations=int(payload["operations"]),
        digest=str(payload["digest"]),
    )


def _stored(run: PencilRun) -> PencilRun:
    return PencilRun(
        generation_s=round(run.generation_s, STORED_SECONDS_DECIMALS),
        decide_s=round(run.decide_s, STORED_SECONDS_DECIMALS),
        first_cut_s=round(run.first_cut_s, STORED_SECONDS_DECIMALS),
    )


def cost_case(key: str, reference_digest: str, launch: Launcher, read_load: LoadReader) -> CasePencilCost:
    """Run *key* `RUNS_PER_CASE` times through *launch*, consecutively; require agreement with each other and with *reference_digest*.

    Args:
        key: A ``pipeline`` gate case key.
        reference_digest: The case's uninstrumented ops digest.
        launch: Runs the case once, in a fresh process for a real measurement.
        read_load: Reads the load average; called before the first run and after the last.

    Returns:
        The case's cost.

    Raises:
        MachineNotQuietError: *read_load* refuses a loaded machine.
        NondeterministicCaseError: Two runs disagree on the digest, the
            operation count, the member count or a call count.
        PencilDigestMismatchError: The runs' digest is not *reference_digest*:
            the instrumentation changed a decision, or the timing document is
            stale against the code.
    """
    load_before = read_load()
    runs = [launch(key) for _run in range(RUNS_PER_CASE)]
    load_after = read_load()
    produced = [(run.digest, run.operations, run.ramp_out_members, tuple(sorted(run.calls.items()))) for run in runs]
    if len(set(produced)) != 1:
        raise NondeterministicCaseError(f"{key}: the runs produced different (digest, operations, members, calls): {produced}.")
    first = runs[0]
    if first.digest != reference_digest:
        raise PencilDigestMismatchError(f"{key}: instrumented digest {first.digest} is not the uninstrumented {reference_digest} of {PIPELINE_TIMING_PATH}.")
    return CasePencilCost(
        runs=tuple(_stored(run.run) for run in runs),
        calls={name: first.calls[name] for name in WRAPPED_CALLABLES},
        ramp_out_members=first.ramp_out_members,
        operations=first.operations,
        digest=first.digest,
        load_before=stored_load(load_before),
        load_after=stored_load(load_after),
    )


def read_this_tree() -> TreeIdentity:
    """This tree's identity (`ROOT`)."""
    return tree_identity(ROOT)


def describe_this_run(tree: TreeIdentity, native: NativeIdentity, started_utc: str, ended_utc: str) -> RunProvenance:
    """The provenance of a pencil-cost run of this tree (`ROOT`) in this interpreter, between *started_utc* and *ended_utc*."""
    rev, dirty, patch = tree
    environment, python_version, stock_2 = interpreter_identity_here()
    host, cpu_model = machine_identity()
    return RunProvenance(
        rev=rev,
        dirty=dirty,
        patch_sha256=patch,
        timing_tool_sha256=file_set_sha256(ROOT, PENCIL_TOOL_SOURCES),
        pixi_lock_sha256=sha256_file(ROOT / PIXI_LOCK),
        environment=environment,
        python_version=python_version,
        stock_2_sha256=stock_2,
        native_sha256=native.native_sha256,
        native_sources_sha256=native.native_sources_sha256,
        host=host,
        cpu_model=cpu_model,
        started_utc=started_utc,
        ended_utc=ended_utc,
    )


def measure_pipeline_pencil_costs(
    keys: Sequence[str],
    *,
    excluded: Sequence[str],
    reference: Mapping[str, str],
    timing_sha256: str,
    launch: Launcher = launch_fresh_process,
    load_reader: LoadReader = os.getloadavg,
    indicative: bool = False,
    describe: Describer = describe_this_run,
    read_tree: TreeReader = read_this_tree,
    track_native: NativeTracker = track_native_here(ROOT),
) -> PipelinePencilCosts:
    """Measure *keys* serially, one case at a time, in order.

    Args:
        keys: The ``pipeline`` case keys to measure.
        excluded: The ``pipeline`` case keys the admissions exclude, recorded unmeasured.
        reference: Case key to its uninstrumented ops digest.
        timing_sha256: SHA-256 of the timing document *reference* came from.
        launch: Runs one case once.
        load_reader: Reads the 1/5/15-minute load average.
        indicative: Measure at any load and mark the record indicative.
        describe: Builds the run's provenance from its tree identity and UTC start and end.
        read_tree: Reads the measured tree's identity; read at the start and the end, which must agree.
        track_native: Starts the native bracket after the start read; its end, after the end read, yields the run's native identity.

    Returns:
        The costs (unvalidated; `PipelinePencilCosts.build` checks coverage).

    Raises:
        MissingTimingReferenceError: A key has no reference digest.
        TreeChangedDuringRunError: The tree's identity changed during the run.
        NativeIdentityError: The native build changed during the run, is stale, or was not compiled from this tree
            (``tools.held_native_identity``; ``NativeSourcesMismatchError`` from ``tools.run_native_identity``).
        MachineNotQuietError: Not *indicative*, and the load was above the
            quiet limit before or after a case.
    """
    missing = [key for key in keys if key not in reference]
    if missing:
        raise MissingTimingReferenceError(f"{PIPELINE_TIMING_PATH} holds no digest for {missing}; run `pixi run pipeline-timing` first.")

    def read_load() -> LoadAverage:
        load = load_reader()
        if not indicative and load[0] > QUIET_LOAD_AVERAGE_MAX:
            raise MachineNotQuietError(
                f"1-minute load average {format_load(load)} (1 / 5 / 15) is above the quiet limit {QUIET_LOAD_AVERAGE_MAX}: "
                f"the pencil cost counts toward the time gate only from a quiet machine. {INDICATIVE_HINT}"
            )
        return load

    tree_at_start, started_utc = read_tree(), utc_now()
    finish_native = track_native()
    cases = {key: cost_case(key, reference[key], launch, read_load) for key in keys}
    ended_utc, tree_at_end = utc_now(), read_tree()
    require_unchanged_tree(tree_at_start, tree_at_end, "the pencil tree")
    native = finish_native()
    return PipelinePencilCosts(
        timing_sha256=timing_sha256, indicative=indicative, cases=cases, excluded=tuple(excluded), provenance=describe(tree_at_end, native, started_utc, ended_utc)
    )


def pencil_cost_keys() -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Every ``pipeline`` gate case key to measure and the ones the committed admissions exclude, in gate-case order."""
    to_time, excluded = gate_timing_keys()
    return tuple(key for key in to_time if is_pencil_key(key)), tuple(key for key in excluded if is_pencil_key(key))


def reference_digests(timings: PipelineTimings) -> dict[str, str]:
    """Case key to uninstrumented ops digest, for every case *timings* timed."""
    return {key: case.digest for key, case in timings.cases.items()}


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry: measure every ``pipeline`` case and write the document and its table, or run one case as a child.

    Relative ``--timing``, ``--write`` and ``--table`` paths resolve against the repository root.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(prog="python -m tools.pipeline_pencil_cost", description="Measure the pencil's decision cost per pipeline case.")
    parser.add_argument("--case", default=None, metavar="KEY", help="child mode: run one case once in this process and print its result line")
    parser.add_argument("--timing", type=Path, default=PIPELINE_TIMING_PATH, help="the uninstrumented timing document whose digests every run must equal")
    parser.add_argument("--write", type=Path, default=PIPELINE_PENCIL_COST_PATH)
    parser.add_argument("--table", type=Path, default=PIPELINE_PENCIL_COST_TABLE_PATH)
    parser.add_argument("--indicative", action="store_true", help="measure at any load; the record is marked indicative and yields no verdict")
    args = parser.parse_args(argv)
    if args.case is not None:
        print(RESULT_LINE_PREFIX + json.dumps(result_payload(run_case_here(args.case))))
        return 0
    timing_path, write, table = rooted(args.timing), rooted(args.write), rooted(args.table)
    timing_sha256 = hashlib.sha256(timing_path.read_bytes()).hexdigest()
    keys, excluded = pencil_cost_keys()
    measured = measure_pipeline_pencil_costs(
        keys, excluded=excluded, reference=reference_digests(load_pipeline_timings(timing_path)), timing_sha256=timing_sha256, indicative=args.indicative
    )
    costs = PipelinePencilCosts.build(measured.to_json())
    write.parent.mkdir(parents=True, exist_ok=True)
    table.parent.mkdir(parents=True, exist_ok=True)
    write.write_text(costs.dumps())
    table.write_text(render_pipeline_pencil_cost_table(costs))
    print(f"wrote {write} and {table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
