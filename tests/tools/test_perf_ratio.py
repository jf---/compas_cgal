"""Orchestration of the ratio harness: pair order, the session bracket, injection and the CLI, on scripted launches and fake arms."""

from __future__ import annotations

import itertools
import json
import os
import subprocess
import sys
from collections.abc import Iterator
from collections.abc import Sequence
from pathlib import Path
from typing import Callable
from typing import NoReturn
from typing import Optional

import pytest

from tests.tools.perf_ratio_fake_arms import FAKE_NATIVE
from tests.tools.perf_ratio_fake_arms import fake_arm
from tests.tools.perf_ratio_fake_arms import fake_arm_tree
from tests.tools.perf_ratio_fake_arms import fake_bracket
from tests.tools.perf_ratio_fake_arms import fake_interpreter
from tests.tools.perf_ratio_fake_arms import fake_tracker
from tools import perf_ratio
from tools.perf_ratio import HARNESS_SOURCES
from tools.perf_ratio import DuplicateWorkloadError
from tools.perf_ratio import HarnessChangedDuringRunError
from tools.perf_ratio import HarnessSourceMissingError
from tools.perf_ratio import HeadWorktreeIsSharedError
from tools.perf_ratio import InvalidInjectionError
from tools.perf_ratio import NondeterministicArmError
from tools.perf_ratio import NoPairsError
from tools.perf_ratio import arm_launcher
from tools.perf_ratio import measure_session
from tools.perf_ratio import measure_workload
from tools.perf_ratio import pair_order
from tools.perf_ratio import parse_injection
from tools.perf_ratio import summary_lines
from tools.perf_ratio import write_session
from tools.perf_ratio_arms import ROOT
from tools.perf_ratio_arms import Arm
from tools.perf_ratio_arms import ArmBracket
from tools.perf_ratio_arms import arm_python
from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import LoadOneMinute
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session
from tools.perf_ratio_workloads import WORKLOAD_SETS
from tools.perf_ratio_workloads import ChildRun
from tools.perf_ratio_workloads import InjectionUnsupportedError
from tools.perf_ratio_workloads import Workload
from tools.perf_ratio_workloads import workload_named
from tools.run_native_identity import PROBE_SNIPPET
from tools.run_provenance import TreeIdentity
from tools.run_provenance import tree_identity

REF, HEAD = ArmName.REF, ArmName.HEAD
# Pairs per workload in the scripted tests (the spec's N).
PAIRS = 5
# The pair (0-based) whose HEAD digest the nondeterminism test changes.
CHANGED_PAIR = 2
# A literal harness digest: the real one covers Task 5's file, which does not exist yet.
FAKE_HARNESS = "e" * 64
# The digest a harness edited mid-run reads at the end.
EDITED_HARNESS = "f" * 64
FAKE_REF_REV = "1" * 40
FAKE_HEAD_REV = "2" * 40
# CPU seconds the end-to-end fake arms report; HEAD is twice the reference.
FAKE_REF_CPU_S = 1.0
FAKE_HEAD_CPU_S = 2.0
FAKE_RATIO = FAKE_HEAD_CPU_S / FAKE_REF_CPU_S
FAKE_DIGEST = "fake-digest"
# The 1-minute load the scripted reader reports, counting up per read so each read is distinguishable.
FIRST_LOAD = 1.0
GT_WORKLOAD = WORKLOAD_SETS["gt-cap40"][0]
POLYGON_WORKLOAD = workload_named("polygon/figure5_refined")


def _identity(rev: str) -> ArmIdentity:
    return ArmIdentity.of(
        tree=TreeIdentity(rev, False, None),
        native=FAKE_NATIVE,
        interpreter=fake_interpreter(),
        pixi_lock_sha256="d" * 64,
        compas_cgal_file=f"/{rev}/src/compas_cgal/__init__.py",
        stock_2_origin=f"/{rev}/build/_stock_2.so",
        editable_build_dir=f"/{rev}/build",
        interpreter_prefix=f"/{rev}/.env",
    )


class Script:
    """A recording launch, load reader and pair of brackets that share one event log."""

    def __init__(self, digest_of: Callable[[ArmName, Workload, int], str], cpu_of: Callable[[ArmName], float] = lambda arm: 1.0 if arm is REF else 1.5) -> None:
        self.events: list[str] = []
        self.calls: dict[ArmName, int] = {REF: 0, HEAD: 0}
        self.digest_of = digest_of
        self.cpu_of = cpu_of
        self.loads = itertools.count(FIRST_LOAD)

    def launch(self, arm: ArmName, workload: Workload) -> ChildRun:
        self.events.append(f"{arm.value} {workload.name}")
        index = self.calls[arm]
        self.calls[arm] += 1
        return ChildRun(CpuSeconds(self.cpu_of(arm)), self.digest_of(arm, workload, index))

    def load(self) -> LoadOneMinute:
        self.events.append("load")
        return LoadOneMinute(float(next(self.loads)))

    def bracket(self, arm: ArmName, rev: str) -> ArmBracket:
        def start() -> Callable[[], ArmIdentity]:
            self.events.append(f"start {arm.value}")

            def finish() -> ArmIdentity:
                self.events.append(f"finish {arm.value}")
                return _identity(rev)

            return finish

        return start

    def session(self, workloads: Sequence[Workload], injected: Optional[Injection] = None, pairs: int = PAIRS) -> PerfRatioSession:
        return measure_session(
            workloads,
            ref_commit=FAKE_REF_REV,
            brackets={REF: self.bracket(REF, FAKE_REF_REV), HEAD: self.bracket(HEAD, FAKE_HEAD_REV)},
            pairs=pairs,
            launch=self.launch,
            load_reader=self.load,
            injected=injected,
            harness=FAKE_HARNESS,
        )


def _same(arm: ArmName, workload: Workload, index: int) -> str:
    return f"digest {workload.name}"


def _fail(*args: object, **kwargs: object) -> NoReturn:
    raise AssertionError(f"ran before the refusal: {args!r}")


def _no_git(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make checkout_head and every subprocess fail if called."""
    monkeypatch.setattr(perf_ratio, "checkout_head", _fail)
    monkeypatch.setattr(subprocess, "run", _fail)


def test_pair_order_alternates_starting_with_the_reference() -> None:
    assert [pair_order(i) for i in range(5)] == [(REF, HEAD), (HEAD, REF), (REF, HEAD), (HEAD, REF), (REF, HEAD)]


def test_each_pair_runs_in_its_order_with_load_read_around_it() -> None:
    script = Script(_same)
    result = measure_workload(GT_WORKLOAD, pairs=PAIRS, launch=script.launch, load_reader=script.load)
    expected: list[str] = []
    for index in range(PAIRS):
        first, second = pair_order(index)
        expected += ["load", f"{first.value} {GT_WORKLOAD.name}", f"{second.value} {GT_WORKLOAD.name}", "load"]
    assert script.events == expected
    assert len(result.pairs) == PAIRS
    assert [pair.order for pair in result.pairs] == [pair_order(i) for i in range(PAIRS)]
    # Each pair's loads are the reads directly around it.
    assert [(pair.load_before, pair.load_after) for pair in result.pairs] == [(FIRST_LOAD + 2 * i, FIRST_LOAD + 2 * i + 1) for i in range(PAIRS)]


def test_runs_that_disagree_within_an_arm_are_refused() -> None:
    def changes(arm: ArmName, workload: Workload, index: int) -> str:
        return "changed" if arm is HEAD and index == CHANGED_PAIR else "same"

    script = Script(changes)
    with pytest.raises(NondeterministicArmError):
        measure_workload(GT_WORKLOAD, pairs=PAIRS, launch=script.launch, load_reader=script.load)


def test_arms_with_different_output_are_recorded_and_flagged() -> None:
    script = Script(lambda arm, workload, index: f"{arm.value} output")
    session = PerfRatioSession.build(script.session([GT_WORKLOAD]).to_json())
    (workload,) = session.workloads
    assert workload.same_output is False
    assert (workload.ref_digest, workload.head_digest) == ("ref output", "head output")


def test_the_session_brackets_both_arms_around_every_workload() -> None:
    workloads = WORKLOAD_SETS["gt-cap40"]
    script = Script(_same)
    script.session(workloads)
    launches = [event for event in script.events if event.startswith(("ref ", "head "))]
    assert len(launches) == 2 * PAIRS * len(workloads)
    measured = [event for event in script.events if event != "load"]
    assert measured == ["start ref", "start head", *launches, "finish ref", "finish head"]


def _fake_tree(path: Path, cpu_seconds: float) -> tuple[Path, str]:
    path.mkdir()
    return path, fake_arm_tree(path, cpu_seconds=cpu_seconds, digest=FAKE_DIGEST)


def _fake_live_bracket(arm: Arm) -> ArmBracket:
    return fake_bracket(arm, trees=(tree_identity(arm.root) for _ in itertools.count()))


def test_end_to_end_on_two_fake_arms_writes_a_valid_record(tmp_path: Path) -> None:
    ref_root, ref_rev = _fake_tree(tmp_path / "ref", FAKE_REF_CPU_S)
    head_root, head_rev = _fake_tree(tmp_path / "head", FAKE_HEAD_CPU_S)
    arms = {REF: fake_arm(REF, ref_root, ref_rev), HEAD: fake_arm(HEAD, head_root, head_rev)}
    workloads = WORKLOAD_SETS["default"]
    session = measure_session(
        workloads,
        ref_commit=ref_rev,
        brackets={name: _fake_live_bracket(arm) for name, arm in arms.items()},
        pairs=2,
        launch=arm_launcher(arms, None),
        load_reader=perf_ratio.load_one_minute,
        injected=None,
        harness=FAKE_HARNESS,
    )
    path = tmp_path / "out" / "session.json"
    write_session(session, path)
    loaded = load_session(path)
    assert len(loaded.workloads) == len(workloads) == 10
    assert [workload.ratio for workload in loaded.workloads] == [FAKE_RATIO] * len(workloads)
    assert all(workload.same_output for workload in loaded.workloads)
    assert (loaded.arms[REF].rev, loaded.arms[HEAD].rev) == (ref_rev, head_rev)
    assert loaded.injected is None


def test_the_launcher_injects_only_the_named_arm(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    arms = {name: fake_arm(name, tmp_path / name.value, rev) for name, rev in ((REF, FAKE_REF_REV), (HEAD, FAKE_HEAD_REV))}
    seen: list[tuple[str, tuple[str, ...], float]] = []

    def command(workload: Workload, python: Sequence[str], inject: CpuSeconds) -> tuple[str, ...]:
        seen.append((workload.name, tuple(python), inject))
        return ("child",)

    def runner(root: Path, env: object) -> Callable[[Sequence[str]], str]:
        return lambda argv: f"{root.name}"

    def parse(workload: Workload, stdout: str) -> ChildRun:
        return ChildRun(CpuSeconds(1.0), stdout)

    monkeypatch.setattr(perf_ratio, "child_command", command)
    monkeypatch.setattr(perf_ratio, "run_in_tree", runner)
    monkeypatch.setattr(perf_ratio, "parse_child", parse)
    launch = arm_launcher(arms, Injection(HEAD, CpuSeconds(0.5)))
    assert launch(REF, GT_WORKLOAD).digest == "ref"
    assert launch(HEAD, GT_WORKLOAD).digest == "head"
    assert seen == [(GT_WORKLOAD.name, arms[REF].python, 0.0), (GT_WORKLOAD.name, arms[HEAD].python, 0.5)]


@pytest.mark.parametrize("text", ["head=abc", "both=1", "head=-1", "head=0", "head=nan", "head", "head=inf", "head=1=2", "=1", "HEAD=1"])
def test_malformed_injection_is_refused(text: str) -> None:
    with pytest.raises(InvalidInjectionError):
        parse_injection(text)


@pytest.mark.parametrize("text", ["head", "ref1"])
def test_an_injection_without_its_separator_names_the_form(text: str) -> None:
    with pytest.raises(InvalidInjectionError, match="is ARM=S"):
        parse_injection(text)


def test_a_well_formed_injection_parses() -> None:
    assert parse_injection("ref=0.25") == Injection(REF, CpuSeconds(0.25))
    assert parse_injection("head=2") == Injection(HEAD, CpuSeconds(2.0))


def test_injection_with_a_non_gt_workload_is_refused_before_anything_runs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _no_git(monkeypatch)
    nonexistent = tmp_path / "missing"
    with pytest.raises(InjectionUnsupportedError):
        perf_ratio.main(["--workload", POLYGON_WORKLOAD.name, "--inject-cpu-seconds", "head=1", "--head-worktree", str(nonexistent)])


def test_a_malformed_injection_on_the_command_line_is_refused_before_anything_runs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _no_git(monkeypatch)
    with pytest.raises(InvalidInjectionError):
        perf_ratio.main(["--set", "gt-cap40", "--inject-cpu-seconds", "head=0", "--head-worktree", str(tmp_path / "missing")])


@pytest.mark.parametrize("which", ["orchestrator", "ref", "ref-alias"])
def test_a_head_worktree_that_is_the_orchestrator_or_the_ref_tree_is_refused(which: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ref = tmp_path / "ref"
    ref.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(ref)
    head = {"orchestrator": ROOT, "ref": tmp_path / "ref", "ref-alias": alias}[which]
    _no_git(monkeypatch)
    with pytest.raises(HeadWorktreeIsSharedError):
        perf_ratio.main(["--ref-worktree", str(ref), "--head-worktree", str(head), "--head-rev", FAKE_HEAD_REV])


@pytest.mark.parametrize("pairs", ["0", "-1"])
def test_a_session_without_pairs_is_refused_before_anything_runs(pairs: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _no_git(monkeypatch)
    with pytest.raises(NoPairsError):
        perf_ratio.main(["--pairs", pairs, "--head-worktree", str(tmp_path / "missing"), "--head-rev", FAKE_HEAD_REV])
    script = Script(_same)
    with pytest.raises(NoPairsError):
        script.session([GT_WORKLOAD], pairs=int(pairs))
    assert script.events == []


def test_a_workload_named_twice_is_refused_before_anything_runs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _no_git(monkeypatch)
    with pytest.raises(DuplicateWorkloadError):
        perf_ratio.main(["--workload", GT_WORKLOAD.name, "--workload", GT_WORKLOAD.name, "--head-worktree", str(tmp_path / "missing"), "--head-rev", FAKE_HEAD_REV])
    script = Script(_same)
    with pytest.raises(DuplicateWorkloadError):
        script.session([GT_WORKLOAD, GT_WORKLOAD])
    assert script.events == []


def test_an_injected_session_carries_the_injection() -> None:
    injected = Injection(HEAD, CpuSeconds(1.0))
    session = Script(_same).session([GT_WORKLOAD], injected=injected)
    assert session.injected == injected
    assert json.loads(session.dumps())["injected"] == {"arm": "head", "seconds": 1.0}
    assert PerfRatioSession.build(session.to_json()).injected == injected


def test_the_harness_digest_covers_the_probe_code() -> None:
    assert Path("tools/held_native_identity.py") in HARNESS_SOURCES
    assert Path("tools/run_native_identity.py") in HARNESS_SOURCES
    assert {path for path in HARNESS_SOURCES if path.name.startswith("perf_ratio")} == {
        Path("tools") / f"perf_ratio{suffix}.py" for suffix in ("", "_record", "_workloads", "_arms", "_history")
    }


# The call the native probe makes once its tree is on sys.path; the closure test swaps it for a module dump.
PROBE_CALL = "from tools.run_native_identity import probe_main; probe_main(sys.argv[1:])"
# Prints the source file of every module the probe's imports loaded, one per line.
MODULE_DUMP = "import tools.run_native_identity; print(chr(10).join(sorted({getattr(m, '__file__', None) or '' for m in list(sys.modules.values())})))"
# The packages the probe puts on sys.path from the orchestrator's tree (``compas_cgal`` resolves through the arm's own finder).
PROBE_PACKAGES = frozenset({"tools", "benchmarks"})


def test_the_harness_digest_covers_every_orchestrator_module_the_probe_imports(tmp_path: Path) -> None:
    assert PROBE_CALL in PROBE_SNIPPET
    env = {name: value for name, value in os.environ.items() if name != "PYTHONPATH"}
    dump = subprocess.run([sys.executable, "-c", PROBE_SNIPPET.replace(PROBE_CALL, MODULE_DUMP)], cwd=tmp_path, env=env, capture_output=True, text=True, check=True)
    loaded = {Path(line).relative_to(ROOT) for line in dump.stdout.splitlines() if line and Path(line).is_relative_to(ROOT)}
    imported = {path for path in loaded if path.parts[0] in PROBE_PACKAGES}
    assert Path("tools") / "run_native_identity.py" in imported
    assert imported <= set(HARNESS_SOURCES), sorted(path.as_posix() for path in imported - set(HARNESS_SOURCES))


def test_summary_lines_name_each_workload_then_the_wall_time() -> None:
    script = Script(_same)
    session = script.session(WORKLOAD_SETS["gt-cap40"][:1])
    (workload,) = session.workloads
    low, high = workload.load_range
    lines = summary_lines(session)
    assert lines[0] == (
        f"{workload.workload}: ratio {workload.ratio:.3f}x (pair median {workload.pair_ratio_median:.3f}x, pairs "
        f"{min(p.ratio for p in workload.pairs):.3f}-{max(p.ratio for p in workload.pairs):.3f}), same_output True, load {low:.1f}-{high:.1f}"
    )
    assert lines[1:] == [f"wall {session.wall_seconds:.0f} s"]


def _patch_main(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, harness_digests: Iterator[str]) -> tuple[list[str], list[str], Path, str, str]:
    """Patch main's steps over two fake arms to record into one event log; return (events, argv, out, ref rev, head rev)."""
    ref_root, ref_rev = _fake_tree(tmp_path / "ref", FAKE_REF_CPU_S)
    head_root, head_rev = _fake_tree(tmp_path / "head", FAKE_HEAD_CPU_S)
    events: list[str] = []

    def record(name: str, function: Callable[..., object], label: Callable[..., str] = lambda *a, **k: "") -> None:
        def wrapper(*args: object, **kwargs: object) -> object:
            events.append(f"{name}{label(*args, **kwargs)}")
            return function(*args, **kwargs)

        monkeypatch.setattr(perf_ratio, name, wrapper)

    def open_fake(name: ArmName, root: Path, python: Sequence[str], rev: str) -> Arm:
        assert tuple(python) == arm_python(root)
        return fake_arm(name, root, rev)

    def tracker(python: Sequence[str], root: Path, env: object) -> object:
        events.append(f"track_native_in {root.name}")
        return fake_tracker()

    steps = (
        "parse_injection",
        "require_injectable",
        "require_plan",
        "require_head_worktree_apart",
        "checkout_head",
        "require_distinct_arms",
        "measure_session",
        "require_unchanged_harness",
        "write_session",
    )
    for name in steps:
        record(name, getattr(perf_ratio, name))
    record("harness_sha256", lambda: next(harness_digests))
    record("open_arm", open_fake, lambda name, *rest: f" {name.value}")
    record("rebuild_arm", lambda arm, build: ("rebuild", str(build)), lambda arm, build: f" {arm.name.value}")
    record("require_settled", lambda track: track()(), lambda track: "")
    record("require_measured_tree", perf_ratio.require_measured_tree, lambda root, rev: f" {root.name}")
    record("live_bracket", _fake_live_bracket, lambda arm: f" {arm.name.value}")
    monkeypatch.setattr(perf_ratio, "track_native_in", tracker)
    monkeypatch.setattr(perf_ratio, "editable_build_dir_of", lambda python, root, env: tmp_path / "orchestrator-build")
    out = tmp_path / "session.json"
    argv = ["--ref-worktree", str(ref_root), "--ref-rev", ref_rev, "--head-worktree", str(head_root), "--head-rev", head_rev]
    argv += ["--set", "gt-cap40", "--pairs", "1", "--inject-cpu-seconds", "head=0.01", "--write", str(out)]
    return events, argv, out, ref_rev, head_rev


def test_main_runs_its_steps_in_order(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    events, argv, out, ref_rev, head_rev = _patch_main(tmp_path, monkeypatch, itertools.repeat(FAKE_HARNESS))
    assert perf_ratio.main(argv) == 0
    assert events == [
        "parse_injection",
        "require_injectable",
        "require_plan",
        "require_head_worktree_apart",
        "harness_sha256",
        "checkout_head",
        "open_arm ref",
        "open_arm head",
        "require_distinct_arms",
        "rebuild_arm ref",
        "rebuild_arm head",
        "track_native_in ref",
        "require_settled",
        "track_native_in head",
        "require_settled",
        "live_bracket ref",
        "live_bracket head",
        "measure_session",
        # measure_session re-checks its plan, then measure_workload once per workload.
        *["require_plan"] * (1 + len(WORKLOAD_SETS["gt-cap40"])),
        "require_unchanged_harness",
        "harness_sha256",
        "require_measured_tree ref",
        "require_measured_tree head",
        "write_session",
    ]
    loaded = load_session(out)
    assert loaded.harness_sha256 == FAKE_HARNESS
    assert loaded.injected == Injection(HEAD, CpuSeconds(0.01))
    assert (loaded.ref_commit, loaded.arms[HEAD].rev) == (ref_rev, head_rev)
    assert [workload.ratio for workload in loaded.workloads] == [FAKE_RATIO] * len(WORKLOAD_SETS["gt-cap40"])
    printed = capsys.readouterr().out.splitlines()
    assert printed[-len(loaded.workloads) - 1 :] == summary_lines(loaded)


def test_a_harness_edited_during_the_run_is_refused_and_nothing_is_written(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    events, argv, out, _, _ = _patch_main(tmp_path, monkeypatch, iter([FAKE_HARNESS, EDITED_HARNESS]))
    with pytest.raises(HarnessChangedDuringRunError):
        perf_ratio.main(argv)
    assert events[-2:] == ["require_unchanged_harness", "harness_sha256"]
    assert not out.exists()


def _harness_tree(root: Path) -> None:
    for path in HARNESS_SOURCES:
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        (root / path).write_text(f"# {path}\n")


def test_an_edited_harness_source_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _harness_tree(tmp_path)
    monkeypatch.setattr(perf_ratio, "ROOT", tmp_path)
    before = perf_ratio.harness_sha256()
    perf_ratio.require_unchanged_harness(before)
    (tmp_path / HARNESS_SOURCES[-1]).write_text("# edited\n")
    with pytest.raises(HarnessChangedDuringRunError):
        perf_ratio.require_unchanged_harness(before)


def test_a_missing_harness_source_is_refused_by_name_before_anything_runs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _harness_tree(tmp_path)
    missing = Path("tools") / "perf_ratio_history.py"
    (tmp_path / missing).unlink()
    monkeypatch.setattr(perf_ratio, "ROOT", tmp_path)
    with pytest.raises(HarnessSourceMissingError, match="perf_ratio_history.py"):
        perf_ratio.harness_sha256()
    _no_git(monkeypatch)
    with pytest.raises(HarnessSourceMissingError):
        perf_ratio.main(["--set", "gt-cap40", "--head-worktree", str(tmp_path / "missing"), "--head-rev", FAKE_HEAD_REV])
