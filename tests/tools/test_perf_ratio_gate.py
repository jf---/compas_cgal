"""The ratio gate CLI: check, calibrate, promote.

Every session is synthetic (`synthetic_session`), every record lives under ``tmp_path``
(``root=tmp_path``), the `Measure` is a fake and promotion's checkout is a fake: no
timing, no arms, no harness subprocess.
"""

from __future__ import annotations

import dataclasses
import hashlib
import subprocess
import sys
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path
from typing import Optional

import pytest
from tools.perf_history_store import REPOSITORY
from tools.perf_ratio import DEFAULT_PAIRS
from tools.perf_ratio_arms import DEFAULT_REF_WORKTREE
from tools.perf_ratio_gate import GATE_SCRATCH
from tools.perf_ratio_gate import CalibrationExistsError
from tools.perf_ratio_gate import CalibrationHarnessMismatchError
from tools.perf_ratio_gate import CalibrationMachineMismatchError
from tools.perf_ratio_gate import MissingRegistryError
from tools.perf_ratio_gate import PromotedReferenceMismatchError
from tools.perf_ratio_gate import PromotionRecordsExistError
from tools.perf_ratio_gate import PromotionRefusedError
from tools.perf_ratio_gate import RebaseWithoutHarnessChangeError
from tools.perf_ratio_gate import harness_command
from tools.perf_ratio_gate import main
from tools.perf_ratio_gate import promote as promote_records
from tools.perf_ratio_gate import subprocess_measure
from tools.perf_ratio_gate_chain import AnchorValueMismatchError
from tools.perf_ratio_gate_decide import UNIT_ANCHOR
from tools.perf_ratio_gate_decide import Verdict
from tools.perf_ratio_gate_decide import decide
from tools.perf_ratio_gate_registry import CALIBRATION_SESSIONS
from tools.perf_ratio_gate_registry import InvalidCalibrationSessionError
from tools.perf_ratio_gate_registry import TooFewCalibrationSessionsError
from tools.perf_ratio_gate_registry import UncalibratedMachineError
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import derive_margins
from tools.perf_ratio_gate_registry import load_registry
from tools.perf_ratio_gate_registry import promotion_dir
from tools.perf_ratio_gate_registry import write_registry
from tools.perf_ratio_gate_verdict import ExtensionRequest
from tools.perf_ratio_gate_verdict import InjectedSessionError
from tools.perf_ratio_gate_verdict import Measure
from tools.perf_ratio_gate_verdict import RatioRegressionError
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session

from tests.tools.perf_ratio_sessions import CPU_MODEL
from tests.tools.perf_ratio_sessions import HARNESS
from tests.tools.perf_ratio_sessions import OTHER_REV
from tests.tools.perf_ratio_sessions import REF
from tests.tools.perf_ratio_sessions import synthetic_session
from tests.tools.perf_ratio_sessions import write_record

GT = "gt/square"
OMEGA = "omega/monstera"
WORKLOADS = (GT, OMEGA)
# The HEAD rev of a gated session (a change against REF).
HEAD = "fedcba9876543210fedcba9876543210fedcba98"
# A second calibrated machine.
OTHER_CPU = "Intel Xeon W-3245"
# HEAD CPU seconds per A/A session (reference 1.0 s): worst deviation 0.02, so the margin is 2 * 0.02 ~ 0.04.
AA_HEAD_CPU = (1.0, 1.01, 0.99, 1.005, 1.02)
# A second A/A series (worst deviation 0.01, margin ~ 0.02), for a re-calibration.
AA_HEAD_CPU_QUIET = (1.0, 1.005, 0.995, 1.01, 0.999)
# HEAD CPU seconds (reference 1.0 s) well inside, below and above a ~0.04 margin.
PASSING, IMPROVED, REGRESSED = 1.0, 0.9, 1.2
# The harness digest of an edited harness (any 64-hex string other than HARNESS).
NEW_HARNESS = "e" * 64
# A workload only the edited harness measures.
NEW_WORKLOAD = "x/new"
# One audit injection on the HEAD arm.
INJECTION = Injection(arm=ArmName.HEAD, seconds=CpuSeconds(0.5))


def pairs(head_cpu: float, count: int = DEFAULT_PAIRS) -> list[tuple[float, float]]:
    """*count* identical pairs at reference 1.0 s and HEAD *head_cpu* s."""
    return [(1.0, head_cpu)] * count


def aa_records(
    root: Path,
    head_cpus: Sequence[float] = AA_HEAD_CPU,
    *,
    cpu_model: str = CPU_MODEL,
    series: str = "aa",
    workloads: Sequence[str] = WORKLOADS,
    **overrides: object,
) -> list[Path]:
    """One A/A record per *head_cpus* entry under *root*/<series>-<k>/, every one of *workloads* at that ratio."""
    paths = []
    for index, head_cpu in enumerate(head_cpus, start=1):
        options: dict[str, object] = {"cpu_model": cpu_model}
        options.update(overrides)
        session = synthetic_session({name: pairs(head_cpu) for name in workloads}, **options)  # type: ignore[arg-type]
        paths.append(write_record(session, root / cpu_model.replace(" ", "-") / f"{series}-{index}" / "perf_ratio_session.json"))
    return paths


def gated(head_cpus: Mapping[str, float], **overrides: object) -> PerfRatioSession:
    """A REF -> HEAD session, workload *w* at HEAD CPU ``head_cpus[w]``."""
    options: dict[str, object] = {"head_rev": HEAD}
    options.update(overrides)
    return synthetic_session({name: pairs(head_cpu) for name, head_cpu in head_cpus.items()}, **options)  # type: ignore[arg-type]


class FakeMeasure:
    """A `Measure` that records each request and answers its n-th with every workload at HEAD CPU ``head_cpus[n]``."""

    def __init__(self, *head_cpus: float, cpu_model: str = CPU_MODEL) -> None:
        self.requests: list[ExtensionRequest] = []
        self.head_cpus = head_cpus
        self.cpu_model = cpu_model

    def __call__(self, request: ExtensionRequest) -> PerfRatioSession:
        self.requests.append(request)
        head_cpu = self.head_cpus[(len(self.requests) - 1) % len(self.head_cpus)]
        return synthetic_session(
            {name: pairs(head_cpu, request.pairs) for name in request.workloads},
            ref_commit=request.ref_commit,
            head_rev=request.head_rev,
            injected=request.injected,
            cpu_model=self.cpu_model,
        )


class FakeCheckout:
    """A checkout that records each call and answers with *answer* (default: the rev it was given)."""

    def __init__(self, answer: Optional[str] = None) -> None:
        self.calls: list[tuple[Path, str]] = []
        self.answer = answer

    def __call__(self, root: Path, rev: str) -> str:
        self.calls.append((root, rev))
        return rev if self.answer is None else self.answer


def no_measure(request: ExtensionRequest) -> PerfRatioSession:
    raise AssertionError(f"no measurement was expected, got {request}")


def no_checkout(root: Path, rev: str) -> str:
    raise AssertionError(f"no checkout was expected, got {root} {rev}")


def no_cpu_model() -> str:
    raise AssertionError("the machine identity was not expected to be read")


def calibrated(tmp_path: Path) -> Path:
    """A registry at REF under HARNESS with CPU_MODEL calibrated from `AA_HEAD_CPU`; returns its path."""
    path = tmp_path / "perf_ratio_gate.json"
    paths = aa_records(tmp_path / "records")
    assert main(["calibrate", "--records", *map(str, paths), "--registry", str(path)], root=tmp_path, measure=no_measure, checkout=no_checkout) == 0
    return path


def check(tmp_path: Path, registry: Path, session: PerfRatioSession, *extra: str, measure: Measure = no_measure) -> int:
    record = write_record(session, tmp_path / "gate" / "perf_ratio_session.json")
    return main(["check", str(record), "--registry", str(registry), *extra], root=tmp_path, measure=measure, checkout=no_checkout)


def promote(tmp_path: Path, registry: Path, session: PerfRatioSession, checkout: Callable[[Path, str], str], measure: Measure = no_measure) -> int:
    record = write_record(session, tmp_path / "gate" / "perf_ratio_session.json")
    return main(["promote", str(record), "--registry", str(registry)], root=tmp_path, measure=measure, checkout=checkout)


# --- harness_command / subprocess_measure ------------------------------------------------------------


def test_harness_command_shape() -> None:
    write = Path("/scratch/ext-1.json")
    honest = ExtensionRequest(ref_commit=REF, head_rev=HEAD, workloads=(GT, OMEGA), pairs=5, injected=None)
    assert harness_command(honest, write) == (
        sys.executable, "-m", "tools.perf_ratio", "--ref-rev", REF, "--head-rev", HEAD,
        "--workload", GT, "--workload", OMEGA, "--pairs", "5", "--write", "/scratch/ext-1.json",
    )  # fmt: skip
    audit = dataclasses.replace(honest, workloads=(GT,), pairs=3, injected=INJECTION)
    assert harness_command(audit, write) == (
        sys.executable, "-m", "tools.perf_ratio", "--ref-rev", REF, "--head-rev", HEAD,
        "--workload", GT, "--pairs", "3", "--write", "/scratch/ext-1.json", "--inject-cpu-seconds", "head=0.5",
    )  # fmt: skip


def test_gate_scratch_is_under_the_gitignored_build_dir() -> None:
    assert GATE_SCRATCH == REPOSITORY / "build" / "ratio-gate"


def test_subprocess_measure_runs_the_harness_and_reads_each_new_record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    runs: list[tuple[tuple[str, ...], bool]] = []
    request = ExtensionRequest(ref_commit=REF, head_rev=HEAD, workloads=(GT,), pairs=5, injected=None)

    def fake_run(command: Sequence[str], *, check: bool, cwd: Path) -> None:
        runs.append((tuple(command), check))
        assert cwd == REPOSITORY
        write_record(gated({GT: 1.0 + len(runs) / 100}), Path(command[command.index("--write") + 1]))

    monkeypatch.setattr(subprocess, "run", fake_run)
    write_record(gated({GT: 3.0}), tmp_path / "ext-1.json")  # a stale record from an earlier gate run
    measure = subprocess_measure(tmp_path)
    first, second = measure(request), measure(request)
    assert runs == [(harness_command(request, tmp_path / "ext-2.json"), True), (harness_command(request, tmp_path / "ext-3.json"), True)]
    assert first == load_session(tmp_path / "ext-2.json") and second == load_session(tmp_path / "ext-3.json")
    assert first.workloads[0].ratio == pytest.approx(1.01) and second.workloads[0].ratio == pytest.approx(1.02)


# --- check -------------------------------------------------------------------------------------------


def test_check_pass_exits_zero_and_prints_each_workload(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    capsys.readouterr()
    assert check(tmp_path, registry, gated({GT: PASSING, OMEGA: 1.01}, head_digest="e")) == 0
    margin = load_registry(registry, root=tmp_path).margins_for(CPU_MODEL)[GT]
    lines = capsys.readouterr().out.splitlines()
    assert lines == [
        f"{GT}: PASS rho 1.000 A_w 1.000 A_w*rho_w 1.000 margin {margin:.3f} pairs {DEFAULT_PAIRS} changed output",
        f"{OMEGA}: PASS rho 1.010 A_w 1.000 A_w*rho_w 1.010 margin {margin:.3f} pairs {DEFAULT_PAIRS} changed output",
        "ratio gate: 2 workloads: 2 PASS, 0 IMPROVED, 0 REGRESSION",
    ]


def test_check_regression_raises_with_workload_rho_margin_pairs(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    capsys.readouterr()
    measure = FakeMeasure(REGRESSED)
    # One message form whatever the anchor: A_w*rho_w against 1 + m_w, with A_w, rho_w and the pairs.
    with pytest.raises(RatioRegressionError, match=rf"{OMEGA}: A_w\*rho_w 1\.2 > 1 \+ m_w 0\.04\d* \(A_w 1\.0, rho_w 1\.2\) over 20 pairs") as raised:
        check(tmp_path, registry, gated({GT: IMPROVED, OMEGA: REGRESSED}), measure=measure)
    assert [verdict.workload for verdict in raised.value.report.regressions] == [OMEGA]
    assert len(measure.requests) == 3
    out = capsys.readouterr().out
    assert f"{OMEGA}: REGRESSION rho 1.200 A_w 1.000 A_w*rho_w 1.200 margin 0.040 pairs 20\n" in out
    assert "ratio gate: 2 workloads: 0 PASS, 1 IMPROVED, 1 REGRESSION" in out
    assert "promotion available" not in out


def test_check_improved_prints_the_promotion_hint(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    capsys.readouterr()
    record = write_record(gated({GT: IMPROVED, OMEGA: PASSING}), tmp_path / "gate" / "perf_ratio_session.json")
    assert main(["check", str(record), "--registry", str(registry)], root=tmp_path, measure=no_measure, checkout=no_checkout) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == f"{GT}: IMPROVED rho 0.900 A_w 1.000 A_w*rho_w 0.900 margin 0.040 pairs {DEFAULT_PAIRS}"
    assert lines[-1] == f"promotion available: pixi run perf-ratio-promote {record}"


def test_check_audit_decides_an_injected_session_and_never_offers_promotion(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    session = gated({GT: IMPROVED}, injected=INJECTION)
    with pytest.raises(InjectedSessionError):
        check(tmp_path, registry, session)
    capsys.readouterr()
    assert check(tmp_path, registry, session, "--audit") == 0
    out = capsys.readouterr().out
    assert "ratio gate (audit): 1 workloads: 0 PASS, 1 IMPROVED, 0 REGRESSION" in out
    assert "promotion available" not in out


def test_check_never_offers_promotion_of_a_dirty_head(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    capsys.readouterr()
    assert check(tmp_path, registry, dirty_head(gated({GT: IMPROVED, OMEGA: PASSING}))) == 0
    assert "promotion available" not in capsys.readouterr().out


def test_check_never_offers_promotion_of_the_reference_itself(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = calibrated(tmp_path)
    capsys.readouterr()
    assert check(tmp_path, path, gated({GT: IMPROVED, OMEGA: PASSING}, head_rev=REF)) == 0
    assert "promotion available" not in capsys.readouterr().out


def test_check_on_an_uncalibrated_machine_is_refused(tmp_path: Path) -> None:
    registry = calibrated(tmp_path)
    with pytest.raises(UncalibratedMachineError):
        check(tmp_path, registry, gated({GT: PASSING, OMEGA: PASSING}, cpu_model=OTHER_CPU))


def test_check_without_a_registry_is_refused(tmp_path: Path) -> None:
    with pytest.raises(MissingRegistryError):
        check(tmp_path, tmp_path / "absent.json", gated({GT: PASSING, OMEGA: PASSING}))


# --- calibrate --records -----------------------------------------------------------------------------


def test_calibrate_from_records_writes_a_loadable_registry(tmp_path: Path) -> None:
    registry = load_registry(calibrated(tmp_path), root=tmp_path)
    sessions = [load_session(path) for path in aa_records(tmp_path / "records")]
    assert registry.ref_commit == REF and registry.harness_sha256 == HARNESS
    assert registry.workloads == tuple(sorted(WORKLOADS))
    assert dict(registry.margins_for(CPU_MODEL)) == derive_margins(sessions)


def test_calibrate_replaces_only_its_cpu(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    other = aa_records(tmp_path / "records", cpu_model=OTHER_CPU)
    assert main(["calibrate", "--records", *map(str, other), "--registry", str(path)], root=tmp_path, measure=no_measure, checkout=no_checkout) == 0
    before = load_registry(path, root=tmp_path)
    quiet = aa_records(tmp_path / "records", AA_HEAD_CPU_QUIET, series="quiet")
    assert main(["calibrate", "--records", *map(str, quiet), "--registry", str(path)], root=tmp_path, measure=no_measure, checkout=no_checkout) == 0
    after = load_registry(path, root=tmp_path)
    (other_before,) = [c for c in before.calibrations if c.cpu_model == OTHER_CPU]
    (other_after,) = [c for c in after.calibrations if c.cpu_model == OTHER_CPU]
    assert other_after == other_before
    assert [c.cpu_model for c in after.calibrations] == sorted([CPU_MODEL, OTHER_CPU])
    assert dict(after.margins_for(CPU_MODEL)) == derive_margins([load_session(p) for p in quiet])
    assert dict(after.margins_for(CPU_MODEL)) != dict(before.margins_for(CPU_MODEL))


@pytest.mark.parametrize("override", [{"ref_commit": OTHER_REV, "head_rev": OTHER_REV}, {"harness": "g" * 64}], ids=["ref", "harness"])
def test_calibrate_without_a_registry_refuses_records_that_disagree(tmp_path: Path, override: dict[str, object]) -> None:
    paths = aa_records(tmp_path / "records")
    (odd,) = aa_records(tmp_path / "records", (1.0,), series="odd", **override)  # type: ignore[arg-type]
    with pytest.raises(InvalidCalibrationSessionError):
        main(
            ["calibrate", "--records", *map(str, [*paths[:-1], odd]), "--registry", str(tmp_path / "perf_ratio_gate.json")],
            root=tmp_path,
            measure=no_measure,
            checkout=no_checkout,
        )
    assert not (tmp_path / "perf_ratio_gate.json").exists()


def test_calibrate_into_a_registry_refuses_records_at_another_reference(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    elsewhere = aa_records(tmp_path / "records", cpu_model=OTHER_CPU, ref_commit=OTHER_REV, head_rev=OTHER_REV)
    with pytest.raises(InvalidCalibrationSessionError):
        main(["calibrate", "--records", *map(str, elsewhere), "--registry", str(path)], root=tmp_path, measure=no_measure, checkout=no_checkout)
    assert path.read_bytes() == before


# --- calibrate --rebase-harness ----------------------------------------------------------------------


def calibrate(tmp_path: Path, registry: Path, records: Sequence[Path], *extra: str) -> int:
    return main(["calibrate", "--records", *map(str, records), "--registry", str(registry), *extra], root=tmp_path, measure=no_measure, checkout=no_checkout)


def two_machine_registry(tmp_path: Path) -> Path:
    """`calibrated`, plus an OTHER_CPU calibration under HARNESS."""
    path = calibrated(tmp_path)
    assert calibrate(tmp_path, path, aa_records(tmp_path / "records", cpu_model=OTHER_CPU)) == 0
    return path


def test_calibrate_refuses_a_new_harness_without_rebase_and_names_the_flag(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    edited = aa_records(tmp_path / "records", cpu_model=OTHER_CPU, harness=NEW_HARNESS)
    with pytest.raises(CalibrationHarnessMismatchError, match="--rebase-harness"):
        calibrate(tmp_path, path, edited)
    assert path.read_bytes() == before


def test_rebase_harness_takes_the_records_harness_and_workloads_and_drops_every_other_calibration(tmp_path: Path) -> None:
    path = two_machine_registry(tmp_path)
    before = load_registry(path, root=tmp_path)
    rebased = aa_records(tmp_path / "records", series="rebased", harness=NEW_HARNESS, workloads=(GT, NEW_WORKLOAD))
    assert calibrate(tmp_path, path, rebased, "--rebase-harness") == 0
    after = load_registry(path, root=tmp_path)
    assert after.ref_commit == before.ref_commit == REF
    assert after.harness_sha256 == NEW_HARNESS
    assert after.workloads == tuple(sorted((GT, NEW_WORKLOAD)))
    assert [c.cpu_model for c in after.calibrations] == [CPU_MODEL]
    assert dict(after.margins_for(CPU_MODEL)) == derive_margins([load_session(p) for p in rebased])
    with pytest.raises(UncalibratedMachineError):
        after.margins_for(OTHER_CPU)


def test_rebase_harness_refuses_records_at_another_reference(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    elsewhere = aa_records(tmp_path / "records", series="elsewhere", harness=NEW_HARNESS, ref_commit=OTHER_REV, head_rev=OTHER_REV)
    with pytest.raises(InvalidCalibrationSessionError, match=REF):
        calibrate(tmp_path, path, elsewhere, "--rebase-harness")
    assert path.read_bytes() == before


def test_rebase_harness_refuses_records_under_the_registry_harness(tmp_path: Path) -> None:
    path = two_machine_registry(tmp_path)
    before = path.read_bytes()
    with pytest.raises(RebaseWithoutHarnessChangeError):
        calibrate(tmp_path, path, aa_records(tmp_path / "records", AA_HEAD_CPU_QUIET, series="quiet"), "--rebase-harness")
    assert path.read_bytes() == before


def test_rebase_harness_needs_a_registry(tmp_path: Path) -> None:
    path = tmp_path / "perf_ratio_gate.json"
    with pytest.raises(MissingRegistryError):
        calibrate(tmp_path, path, aa_records(tmp_path / "records", harness=NEW_HARNESS), "--rebase-harness")
    assert not path.exists()


def test_calibrate_sessions_rebase_harness_measures_at_the_reference_and_rebases(tmp_path: Path) -> None:
    path = tmp_path / "perf_ratio_gate.json"
    write_registry(build_registry(ref_commit=REF, harness_sha256=NEW_HARNESS, calibrations=(), workloads=WORKLOADS), path, root=tmp_path)
    measure = FakeMeasure(*AA_HEAD_CPU)
    argv = ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--rebase-harness", "--registry", str(path)]
    assert main(argv, root=tmp_path, measure=measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL, tree_harness=lambda: HARNESS) == 0
    assert {(request.ref_commit, request.head_rev) for request in measure.requests} == {(REF, REF)}
    registry = load_registry(path, root=tmp_path)
    assert registry.harness_sha256 == HARNESS and registry.ref_commit == REF
    assert [c.cpu_model for c in registry.calibrations] == [CPU_MODEL]


# --- calibrate --sessions ----------------------------------------------------------------------------


def uncalibrated_registry(tmp_path: Path, ref_commit: str = REF) -> Path:
    """A calibration-less registry at *ref_commit* under HARNESS over `WORKLOADS` (a freshly promoted one)."""
    path = tmp_path / "perf_ratio_gate.json"
    write_registry(build_registry(ref_commit=ref_commit, harness_sha256=HARNESS, calibrations=(), workloads=WORKLOADS), path, root=tmp_path)
    return path


def test_calibrate_after_promotion_refuses_records_at_the_old_reference(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path, ref_commit=HEAD)
    before = path.read_bytes()
    with pytest.raises(InvalidCalibrationSessionError):
        main(["calibrate", "--records", *map(str, aa_records(tmp_path / "records")), "--registry", str(path)], root=tmp_path, measure=no_measure, checkout=no_checkout)
    assert path.read_bytes() == before


def sessions_target(tmp_path: Path, cpu_model: str = CPU_MODEL) -> Path:
    return tmp_path / calibration_dir(REF, cpu_model).relative_to(REPOSITORY)


def test_calibrate_sessions_runs_k_aa_sessions_and_calibrates(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    measure = FakeMeasure(*AA_HEAD_CPU, 1.015)
    sessions = CALIBRATION_SESSIONS + 1
    assert (
        main(
            ["calibrate", "--sessions", str(sessions), "--registry", str(path)],
            root=tmp_path,
            measure=measure,
            checkout=no_checkout,
            cpu_model=lambda: CPU_MODEL,
            tree_harness=lambda: HARNESS,
        )
        == 0
    )
    assert measure.requests == [ExtensionRequest(ref_commit=REF, head_rev=REF, workloads=tuple(sorted(WORKLOADS)), pairs=DEFAULT_PAIRS, injected=None)] * sessions
    target = sessions_target(tmp_path)
    records = [target / f"aa-{k}" / "perf_ratio_session.json" for k in range(1, sessions + 1)]
    assert sorted(target.iterdir()) == sorted(record.parent for record in records)
    registry = load_registry(path, root=tmp_path)
    (calibration,) = registry.calibrations
    assert [record.path for record in calibration.sessions] == [record.relative_to(tmp_path).as_posix() for record in records]
    assert dict(calibration.margins) == derive_margins([load_session(record) for record in records])


def test_calibrate_sessions_refuses_an_existing_target_before_any_run(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    sessions_target(tmp_path).mkdir(parents=True)
    before = path.read_bytes()
    with pytest.raises(CalibrationExistsError):
        main(
            ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(path)],
            root=tmp_path,
            measure=no_measure,
            checkout=no_checkout,
            cpu_model=lambda: CPU_MODEL,
            tree_harness=lambda: HARNESS,
        )
    assert path.read_bytes() == before


def test_calibrate_sessions_refuses_fewer_than_k_sessions(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    with pytest.raises(TooFewCalibrationSessionsError):
        main(
            ["calibrate", "--sessions", str(CALIBRATION_SESSIONS - 1), "--registry", str(path)],
            root=tmp_path,
            measure=no_measure,
            checkout=no_checkout,
            cpu_model=lambda: CPU_MODEL,
            tree_harness=lambda: HARNESS,
        )
    assert not sessions_target(tmp_path).exists()


def test_calibrate_sessions_needs_a_registry(tmp_path: Path) -> None:
    with pytest.raises(MissingRegistryError):
        main(
            ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(tmp_path / "absent.json")],
            root=tmp_path,
            measure=no_measure,
            checkout=no_checkout,
            cpu_model=no_cpu_model,
            tree_harness=lambda: HARNESS,
        )


def test_calibrate_sessions_refuses_a_session_from_another_machine(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    before = path.read_bytes()
    measure = FakeMeasure(1.0, cpu_model=OTHER_CPU)
    with pytest.raises(CalibrationMachineMismatchError):
        main(
            ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(path)],
            root=tmp_path,
            measure=measure,
            checkout=no_checkout,
            cpu_model=lambda: CPU_MODEL,
            tree_harness=lambda: HARNESS,
        )
    assert len(measure.requests) == 1
    assert not sessions_target(tmp_path).exists()
    assert path.read_bytes() == before


def test_calibrate_sessions_refuses_an_edited_tree_harness_before_any_run(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    before = path.read_bytes()
    measure = FakeMeasure(*AA_HEAD_CPU)
    argv = ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(path)]
    with pytest.raises(CalibrationHarnessMismatchError, match="--rebase-harness"):
        main(argv, root=tmp_path, measure=measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL, tree_harness=lambda: NEW_HARNESS)
    assert measure.requests == []
    assert not sessions_target(tmp_path).exists()
    assert path.read_bytes() == before


def test_calibrate_sessions_refuses_rebase_under_an_unchanged_tree_harness_before_any_run(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    before = path.read_bytes()
    measure = FakeMeasure(*AA_HEAD_CPU)
    argv = ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--rebase-harness", "--registry", str(path)]
    with pytest.raises(RebaseWithoutHarnessChangeError, match="without it"):
        main(argv, root=tmp_path, measure=measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL, tree_harness=lambda: HARNESS)
    assert measure.requests == []
    assert not sessions_target(tmp_path).exists()
    assert path.read_bytes() == before


class FailingMeasure(FakeMeasure):
    """A `FakeMeasure` whose *fail_at*-th request fails as a harness that exits non-zero."""

    def __init__(self, fail_at: int, *head_cpus: float) -> None:
        super().__init__(*head_cpus)
        self.fail_at = fail_at

    def __call__(self, request: ExtensionRequest) -> PerfRatioSession:
        if len(self.requests) + 1 == self.fail_at:
            self.requests.append(request)
            raise subprocess.CalledProcessError(1, ("perf_ratio",))
        return super().__call__(request)


def test_a_failed_sessions_run_leaves_no_calibration_directory_and_reruns(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    before = path.read_bytes()
    argv = ["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(path)]
    failing = FailingMeasure(2, *AA_HEAD_CPU)
    with pytest.raises(subprocess.CalledProcessError):
        main(argv, root=tmp_path, measure=failing, checkout=no_checkout, cpu_model=lambda: CPU_MODEL, tree_harness=lambda: HARNESS)
    assert len(failing.requests) == 2
    assert not sessions_target(tmp_path).exists()
    assert not (tmp_path / calibration_dir(REF, CPU_MODEL).relative_to(REPOSITORY).parents[1]).exists()
    assert path.read_bytes() == before
    rerun = FakeMeasure(*AA_HEAD_CPU)
    assert main(argv, root=tmp_path, measure=rerun, checkout=no_checkout, cpu_model=lambda: CPU_MODEL, tree_harness=lambda: HARNESS) == 0
    assert len(rerun.requests) == CALIBRATION_SESSIONS
    assert sorted(child.name for child in sessions_target(tmp_path).iterdir()) == [f"aa-{k}" for k in range(1, CALIBRATION_SESSIONS + 1)]
    (calibration,) = load_registry(path, root=tmp_path).calibrations
    assert calibration.cpu_model == CPU_MODEL


# --- promote -----------------------------------------------------------------------------------------


def dirty_head(session: PerfRatioSession) -> PerfRatioSession:
    head = dataclasses.replace(session.arms[ArmName.HEAD], dirty=True, patch_sha256="p" * 64)
    return dataclasses.replace(session, arms={ArmName.REF: session.arms[ArmName.REF], ArmName.HEAD: head})


def test_promote_moves_the_reference_and_clears_calibrations(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = calibrated(tmp_path)
    before = load_registry(path, root=tmp_path)
    checkout = FakeCheckout()
    capsys.readouterr()
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: PASSING}), checkout) == 0
    assert checkout.calls == [(DEFAULT_REF_WORKTREE, HEAD)]
    promoted = load_registry(path, root=tmp_path)
    assert promoted.ref_commit == HEAD
    assert promoted.calibrations == ()
    assert promoted.harness_sha256 == before.harness_sha256 and promoted.workloads == before.workloads
    with pytest.raises(UncalibratedMachineError):
        promoted.margins_for(CPU_MODEL)
    assert f"next: pixi run perf-ratio-calibrate --sessions {CALIBRATION_SESSIONS}" in capsys.readouterr().out


PROMOTION_REFUSALS: dict[str, tuple[Callable[[], PerfRatioSession], type[Exception], str]] = {
    "regression": (lambda: gated({GT: IMPROVED, OMEGA: REGRESSED}), PromotionRefusedError, "REGRESSION"),
    "no-improvement": (lambda: gated({GT: PASSING, OMEGA: PASSING}), PromotionRefusedError, "no IMPROVED"),
    "dirty-head": (lambda: dirty_head(gated({GT: IMPROVED, OMEGA: PASSING})), PromotionRefusedError, "dirty"),
    "injected": (lambda: gated({GT: IMPROVED, OMEGA: PASSING}, injected=INJECTION), InjectedSessionError, "audit session"),
    "head-is-reference": (lambda: gated({GT: IMPROVED, OMEGA: PASSING}, head_rev=REF), PromotionRefusedError, "is the reference"),
}


@pytest.mark.parametrize("reason", sorted(PROMOTION_REFUSALS))
def test_promotion_refusals_leave_the_registry_alone(tmp_path: Path, reason: str) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    build, error, message = PROMOTION_REFUSALS[reason]
    with pytest.raises(error, match=message):
        promote(tmp_path, path, build(), no_checkout, measure=FakeMeasure(REGRESSED))
    assert path.read_bytes() == before


def test_a_promotion_refusal_names_every_reason(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    with pytest.raises(PromotionRefusedError, match=r"no IMPROVED.*; .*dirty"):
        promote(tmp_path, path, dirty_head(gated({GT: PASSING, OMEGA: PASSING})), no_checkout)


def test_a_checkout_that_lands_elsewhere_leaves_the_registry_alone(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    with pytest.raises(PromotedReferenceMismatchError):
        promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: PASSING}), FakeCheckout(answer=OTHER_REV))
    assert path.read_bytes() == before


# --- promote: the anchor chain -----------------------------------------------------------------------

# The second and third promoted commits of a chain (REF -> HEAD -> HEAD2 -> HEAD3).
HEAD2 = "2222222222222222222222222222222222222222"
HEAD3 = "3333333333333333333333333333333333333333"
# HEAD CPU seconds (reference 1.0 s) of a non-improved workload that drifts inside a ~0.04 margin.
DRIFT = 1.03
# Red against a ~0.04 margin at first; its extension's pairs merge down to a PASS ratio.
RED_FIRST, EXTENDED_TO = 1.05, 1.01


def calibrate_at(tmp_path: Path, registry: Path, ref: str) -> None:
    """Calibrate CPU_MODEL at *ref* from `AA_HEAD_CPU` (margin ~ 0.04), its records in their own series."""
    assert calibrate(tmp_path, registry, aa_records(tmp_path / "records", series=f"aa-{ref[:8]}", ref_commit=ref, head_rev=ref)) == 0


def promotion_target(tmp_path: Path, from_ref: str, to_ref: str) -> Path:
    return tmp_path / promotion_dir(from_ref, to_ref).relative_to(REPOSITORY)


def gate_record(tmp_path: Path) -> Path:
    """Where `promote` and `check` write the session they gate."""
    return tmp_path / "gate" / "perf_ratio_session.json"


def test_promote_carries_a_drifted_anchor_and_resets_an_improved_one(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = calibrated(tmp_path)
    before = load_registry(path, root=tmp_path)
    capsys.readouterr()
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}), FakeCheckout()) == 0
    promoted = load_registry(path, root=tmp_path)
    # Carried: A' = A * rho exactly (A = 1 at a fresh registry); reset: the improved workload's anchor is 1.
    assert dict(promoted.anchors_exact()) == {GT: Fraction(1), OMEGA: Fraction(DRIFT)}
    assert promoted.ref_commit == HEAD and promoted.anchor_ref == REF and promoted.calibrations == ()
    (link,) = promoted.anchor_links
    assert (link.from_ref, link.to_ref) == (REF, HEAD)
    assert link.calibration == before.calibrations[0]
    target = promotion_target(tmp_path, REF, HEAD)
    assert [record.path for record in link.sessions] == [(target / "gate" / "perf_ratio_session.json").relative_to(tmp_path).as_posix()]
    assert (target / "gate" / "perf_ratio_session.json").read_bytes() == gate_record(tmp_path).read_bytes()
    out = capsys.readouterr().out.splitlines()
    assert f"  {GT}: A_w 1.000 -> 1.000 (reset)" in out
    assert f"  {OMEGA}: A_w 1.000 -> 1.030 (carried)" in out


def test_promote_persists_the_extensions_and_anchors_on_the_extended_ratio(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    measure = FakeMeasure(EXTENDED_TO)
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: RED_FIRST}), FakeCheckout(), measure=measure) == 0
    assert [request.workloads for request in measure.requests] == [(OMEGA,)]
    promoted = load_registry(path, root=tmp_path)
    (link,) = promoted.anchor_links
    target = promotion_target(tmp_path, REF, HEAD)
    assert [record.path for record in link.sessions] == [(target / name / "perf_ratio_session.json").relative_to(tmp_path).as_posix() for name in ("gate", "ext-1")]
    assert sorted(child.name for child in target.iterdir()) == ["ext-1", "gate"]
    (extension,) = load_session(target / "ext-1" / "perf_ratio_session.json").workloads
    assert extension.workload == OMEGA and len(extension.pairs) == 5
    # The extended run's ratio is min head / min ref over all ten pairs: EXTENDED_TO, not RED_FIRST.
    assert promoted.anchors_exact()[OMEGA] == Fraction(EXTENDED_TO)


def test_the_promoted_registry_reloads_to_the_registry_promote_returned(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    record = write_record(gated({GT: IMPROVED, OMEGA: DRIFT}), gate_record(tmp_path))
    returned = promote_records(record, registry_path=path, root=tmp_path, measure=no_measure, checkout=FakeCheckout())
    assert load_registry(path, root=tmp_path) == returned
    assert path.read_text(encoding="utf-8") == returned.dumps()


def test_a_second_promotion_composes_on_the_first(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}), FakeCheckout()) == 0
    calibrate_at(tmp_path, path, HEAD)
    # Against the anchors {GT: 1, OMEGA: 1.03}: GT drifts 1.02 (carried), OMEGA at 0.9 is 0.927 < 1 - m (reset).
    assert promote(tmp_path, path, gated({GT: 1.02, OMEGA: IMPROVED}, ref_commit=HEAD, head_rev=HEAD2), FakeCheckout()) == 0
    promoted = load_registry(path, root=tmp_path)
    assert promoted.ref_commit == HEAD2 and promoted.anchor_ref == REF
    assert [(link.from_ref, link.to_ref) for link in promoted.anchor_links] == [(REF, HEAD), (HEAD, HEAD2)]
    assert dict(promoted.anchors_exact()) == {GT: Fraction(1.02), OMEGA: Fraction(1)}


def test_a_drift_creeping_over_two_promotions_is_a_regression_against_its_anchor(tmp_path: Path) -> None:
    """The creep the landed gate (9ae1cdd2) let through, caught.

    OMEGA never improves: it drifts 1.03 at the first promotion and 1.009 at the second, and its
    third session measures 1.03 again. Every one of those ratios is inside ``1 + m_w`` on its
    own, so on the landed gate's semantics (``rho_w`` against ``1 + m_w``, no anchor) the third
    session PASSES and the drift compounds without bound. Against its anchor
    ``A_w = 1.03 * 1.009`` the third session is ``A_w * rho_w ~ 1.07 > 1 + m_w``: REGRESSION. A
    second +1.03 drift is already refused at the second promotion.
    """
    path = calibrated(tmp_path)
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}), FakeCheckout()) == 0
    calibrate_at(tmp_path, path, HEAD)
    after_first = path.read_bytes()
    with pytest.raises(PromotionRefusedError, match=f"REGRESSION on {OMEGA}"):
        promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}, ref_commit=HEAD, head_rev=HEAD2), no_checkout, measure=FakeMeasure(DRIFT))
    assert path.read_bytes() == after_first
    second = 1.009
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: second}, ref_commit=HEAD, head_rev=HEAD2), FakeCheckout()) == 0
    calibrate_at(tmp_path, path, HEAD2)
    registry = load_registry(path, root=tmp_path)
    assert registry.anchors_exact()[OMEGA] == Fraction(DRIFT) * Fraction(second)
    third = gated({GT: PASSING, OMEGA: DRIFT}, ref_commit=HEAD2, head_rev=HEAD3)
    margin = registry.margins_for(CPU_MODEL)[OMEGA]
    (omega_run,) = [run for run in third.workloads if run.workload == OMEGA]
    # The landed gate's decision: rho_w alone, no anchor.
    assert decide(omega_run, margin, UNIT_ANCHOR) is Verdict.PASS
    with pytest.raises(RatioRegressionError, match=rf"{OMEGA}: A_w\*rho_w 1\.07\d* > 1 \+ m_w 0\.04\d* \(A_w 1\.039\d*, rho_w 1\.03\) over 20 pairs"):
        check(tmp_path, path, third, measure=FakeMeasure(DRIFT))


def test_promotion_refusals_name_the_anchor(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    with pytest.raises(PromotionRefusedError, match=r"REGRESSION on omega/monstera \(A_w\*rho_w > 1 \+ m_w\)") as raised:
        promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: REGRESSED}), no_checkout, measure=FakeMeasure(REGRESSED))
    assert "never moves to a slower commit" not in str(raised.value)
    with pytest.raises(PromotionRefusedError, match=r"no IMPROVED workload; promotion needs a certified improvement \(A_w\*rho_w < 1 - m_w\)"):
        promote(tmp_path, path, gated({GT: PASSING, OMEGA: PASSING}), no_checkout)


class RaisingCheckout(FakeCheckout):
    """A checkout that fails as ``git checkout`` exiting non-zero."""

    def __call__(self, root: Path, rev: str) -> str:
        self.calls.append((root, rev))
        raise subprocess.CalledProcessError(128, ("git", "checkout", rev))


FAILED_CHECKOUTS: dict[str, tuple[Callable[[], FakeCheckout], type[Exception]]] = {
    "lands-elsewhere": (lambda: FakeCheckout(answer=OTHER_REV), PromotedReferenceMismatchError),
    "raises": (RaisingCheckout, subprocess.CalledProcessError),
}


@pytest.mark.parametrize("failure", sorted(FAILED_CHECKOUTS))
def test_a_failed_checkout_leaves_the_registry_and_a_rerun_finishes_without_measuring(tmp_path: Path, failure: str) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    session = gated({GT: IMPROVED, OMEGA: RED_FIRST})
    build_checkout, error = FAILED_CHECKOUTS[failure]
    with pytest.raises(error):
        promote(tmp_path, path, session, build_checkout(), measure=FakeMeasure(EXTENDED_TO))
    assert path.read_bytes() == before
    target = promotion_target(tmp_path, REF, HEAD)
    assert sorted(child.name for child in target.iterdir()) == ["ext-1", "gate"]
    # The rerun replays the persisted extension instead of measuring a new one, so it cites the same records.
    checkout = FakeCheckout()
    assert promote(tmp_path, path, session, checkout, measure=no_measure) == 0
    assert checkout.calls == [(DEFAULT_REF_WORKTREE, HEAD)]
    promoted = load_registry(path, root=tmp_path)
    (link,) = promoted.anchor_links
    assert [record.path for record in link.sessions] == [(target / name / "perf_ratio_session.json").relative_to(tmp_path).as_posix() for name in ("gate", "ext-1")]
    assert promoted.anchors_exact()[OMEGA] == Fraction(EXTENDED_TO)


def test_promotion_records_of_another_session_are_refused_before_measuring(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    target = promotion_target(tmp_path, REF, HEAD)
    other = write_record(gated({GT: IMPROVED, OMEGA: PASSING}), target / "gate" / "perf_ratio_session.json").read_bytes()
    # An extension the promoted session's gate would accept: without the identity check the gate would run on it.
    extension = write_record(gated({OMEGA: EXTENDED_TO}), target / "ext-1" / "perf_ratio_session.json").read_bytes()
    capsys.readouterr()
    with pytest.raises(PromotionRecordsExistError) as raised:
        promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: RED_FIRST}), no_checkout, measure=no_measure)
    # Refused before the gate runs: no verdict is printed (checked before the message), and the foreign records are untouched.
    assert capsys.readouterr().out == ""
    assert "gate/perf_ratio_session.json is not the record" in str(raised.value)
    assert path.read_bytes() == before
    assert (target / "gate" / "perf_ratio_session.json").read_bytes() == other
    assert (target / "ext-1" / "perf_ratio_session.json").read_bytes() == extension


def test_promotion_records_beyond_the_ones_promote_would_write_are_refused(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    session = gated({GT: IMPROVED, OMEGA: PASSING})
    with pytest.raises(PromotedReferenceMismatchError):
        promote(tmp_path, path, session, FakeCheckout(answer=OTHER_REV))
    target = promotion_target(tmp_path, REF, HEAD)
    write_record(gated({OMEGA: EXTENDED_TO}), target / "ext-1" / "perf_ratio_session.json")
    with pytest.raises(PromotionRecordsExistError, match=r"other records than this promotion writes \(ext-1/perf_ratio_session.json differ\)"):
        promote(tmp_path, path, session, no_checkout)
    assert path.read_bytes() == before


def test_a_rerun_whose_gate_asks_for_more_extensions_than_were_persisted_is_refused(tmp_path: Path) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    session = gated({GT: IMPROVED, OMEGA: RED_FIRST})
    target = promotion_target(tmp_path, REF, HEAD)
    write_record(session, target / "gate" / "perf_ratio_session.json")
    with pytest.raises(PromotionRecordsExistError, match="no further extension record"):
        promote(tmp_path, path, session, no_checkout, measure=no_measure)
    assert path.read_bytes() == before


def promoted_once(tmp_path: Path) -> Path:
    """`calibrated`, promoted REF -> HEAD (GT reset, OMEGA carried at `DRIFT`): one link, no calibrations."""
    path = calibrated(tmp_path)
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}), FakeCheckout()) == 0
    return path


def test_calibrate_records_keeps_the_anchor_chain(tmp_path: Path) -> None:
    path = promoted_once(tmp_path)
    before = load_registry(path, root=tmp_path)
    calibrate_at(tmp_path, path, HEAD)
    after = load_registry(path, root=tmp_path)
    assert [c.cpu_model for c in after.calibrations] == [CPU_MODEL]
    assert (after.anchor_ref, after.anchor_links) == (before.anchor_ref, before.anchor_links) == (REF, before.anchor_links)
    assert len(after.anchor_links) == 1
    assert dict(after.anchors_exact()) == dict(before.anchors_exact()) == {GT: Fraction(1), OMEGA: Fraction(DRIFT)}


def test_rebase_harness_keeps_the_anchor_chain(tmp_path: Path) -> None:
    path = promoted_once(tmp_path)
    calibrate_at(tmp_path, path, HEAD)
    before = load_registry(path, root=tmp_path)
    rebased = aa_records(tmp_path / "records", series="rebased", harness=NEW_HARNESS, ref_commit=HEAD, head_rev=HEAD)
    assert calibrate(tmp_path, path, rebased, "--rebase-harness") == 0
    after = load_registry(path, root=tmp_path)
    assert after.harness_sha256 == NEW_HARNESS and after.ref_commit == HEAD
    assert (after.anchor_ref, after.anchor_links) == (before.anchor_ref, before.anchor_links)
    assert len(after.anchor_links) == 1
    assert dict(after.anchors_exact()) == {GT: Fraction(1), OMEGA: Fraction(DRIFT)}


# --- promote: the carry and reset rules, persisted (anchor-t3 fix wave) ---------------------------


def test_a_carried_anchor_is_persisted_as_its_anchor_times_rho(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """The carry rule ``A'_w = A_w * rho_w``, read back from the written registry.

    OMEGA enters at ``A_w = 1.03`` and measures ``rho_w = 0.95``: improved on ``rho_w`` alone
    (``< 1 - m_w``), but ``A_w * rho_w ~ 0.979`` is not, so the drift is carried, not reset.
    """
    path = promoted_once(tmp_path)
    calibrate_at(tmp_path, path, HEAD)
    rho = 0.95
    capsys.readouterr()
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: rho}, ref_commit=HEAD, head_rev=HEAD2), FakeCheckout()) == 0
    persisted = load_registry(path, root=tmp_path).anchors_exact()
    assert persisted[OMEGA] == Fraction(DRIFT) * Fraction(rho)
    assert f"  {OMEGA}: A_w 1.030 -> 0.978 (carried)" in capsys.readouterr().out.splitlines()


def test_a_reset_anchor_is_persisted_as_one(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """The reset rule ``A'_w = 1`` on a certified improvement (``A_w * rho_w < 1 - m_w``), read back.

    OMEGA enters at ``A_w = 1.03`` and measures 0.9: ``A_w * rho_w = 0.927``, so its anchor is 1,
    neither ``rho_w`` nor ``A_w * rho_w``. GT (``A_w = 1``, flat) is carried at 1.
    """
    path = promoted_once(tmp_path)
    calibrate_at(tmp_path, path, HEAD)
    capsys.readouterr()
    assert promote(tmp_path, path, gated({GT: PASSING, OMEGA: IMPROVED}, ref_commit=HEAD, head_rev=HEAD2), FakeCheckout()) == 0
    persisted = load_registry(path, root=tmp_path).anchors_exact()
    assert persisted[OMEGA] == Fraction(1)
    assert persisted[GT] == Fraction(1)
    out = capsys.readouterr().out.splitlines()
    assert f"  {OMEGA}: A_w 1.030 -> 1.000 (reset)" in out
    assert f"  {GT}: A_w 1.000 -> 1.000 (carried)" in out


# --- promote: the record is read once (review M1) -------------------------------------------------


def test_a_record_rewritten_during_the_gate_persists_the_bytes_that_were_decided(tmp_path: Path) -> None:
    """Promote reads *record* once; the gate, the persisted copy and the rerun check all use those bytes.

    A `Measure` rewrites the record mid-gate (GT no longer improved). The promotion still decides,
    persists and cites the bytes read at the start, and leaves no orphan promotion directory or staging.
    """
    path = calibrated(tmp_path)
    record = write_record(gated({GT: IMPROVED, OMEGA: RED_FIRST}), gate_record(tmp_path))
    decided = record.read_bytes()
    inner = FakeMeasure(EXTENDED_TO)

    def rewriting_measure(request: ExtensionRequest) -> PerfRatioSession:
        write_record(gated({GT: PASSING, OMEGA: RED_FIRST}), record)
        return inner(request)

    promoted = promote_records(record, registry_path=path, root=tmp_path, measure=rewriting_measure, checkout=FakeCheckout())
    assert record.read_bytes() != decided
    target = promotion_target(tmp_path, REF, HEAD)
    assert (target / "gate" / "perf_ratio_session.json").read_bytes() == decided
    (link,) = promoted.anchor_links
    assert [cited.path for cited in link.sessions] == [(target / name / "perf_ratio_session.json").relative_to(tmp_path).as_posix() for name in ("gate", "ext-1")]
    assert link.sessions[0].record_sha256 == hashlib.sha256(decided).hexdigest()
    # Decided on the bytes read first: GT improved (reset), OMEGA extended to EXTENDED_TO (carried).
    assert dict(promoted.anchors_exact()) == {GT: Fraction(1), OMEGA: Fraction(EXTENDED_TO)}
    assert load_registry(path, root=tmp_path) == promoted
    # No orphan: the only promotion directory is the cited one, and no staging directory is left.
    assert [child.name for child in target.parent.iterdir()] == [target.name]
    scratch = tmp_path / GATE_SCRATCH.relative_to(REPOSITORY)
    assert list(scratch.glob("promotion-staging-*")) == []


# --- promote: the anchor cross-check (review M2) --------------------------------------------------


def test_the_anchor_cross_check_refuses_a_report_the_chain_does_not_replay(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A CLI carry that disagrees with the chain is `AnchorValueMismatchError` before checkout.

    The records are already staged when the guard fires; they stay (a rerun reproduces the error).
    """
    path = calibrated(tmp_path)
    before = path.read_bytes()
    monkeypatch.setattr("tools.perf_ratio_gate.carry_anchor", lambda anchor, ratio, margin: UNIT_ANCHOR)
    with pytest.raises(AnchorValueMismatchError, match="re-derives the anchors"):
        promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}), no_checkout)
    assert path.read_bytes() == before
    assert sorted(child.name for child in promotion_target(tmp_path, REF, HEAD).iterdir()) == ["gate"]


# --- promote: coverage gaps (review M6) -----------------------------------------------------------


def test_a_failed_registry_write_leaves_the_registry_and_a_rerun_finishes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = calibrated(tmp_path)
    before = path.read_bytes()
    session = gated({GT: IMPROVED, OMEGA: RED_FIRST})

    def failing_write(*args: object, **kwargs: object) -> None:
        raise OSError("disk full")

    with monkeypatch.context() as patch:
        patch.setattr("tools.perf_ratio_gate.write_registry", failing_write)
        with pytest.raises(OSError, match="disk full"):
            promote(tmp_path, path, session, FakeCheckout(), measure=FakeMeasure(EXTENDED_TO))
    assert path.read_bytes() == before
    target = promotion_target(tmp_path, REF, HEAD)
    assert sorted(child.name for child in target.iterdir()) == ["ext-1", "gate"]
    checkout = FakeCheckout()
    assert promote(tmp_path, path, session, checkout, measure=no_measure) == 0
    assert checkout.calls == [(DEFAULT_REF_WORKTREE, HEAD)]
    promoted = load_registry(path, root=tmp_path)
    (link,) = promoted.anchor_links
    assert [cited.path for cited in link.sessions] == [(target / name / "perf_ratio_session.json").relative_to(tmp_path).as_posix() for name in ("gate", "ext-1")]
    assert dict(promoted.anchors_exact()) == {GT: Fraction(1), OMEGA: Fraction(EXTENDED_TO)}


@pytest.mark.parametrize("cpu", [CPU_MODEL, OTHER_CPU])
def test_the_link_cites_the_promoting_machines_calibration(tmp_path: Path, cpu: str) -> None:
    path = two_machine_registry(tmp_path)
    calibrations = {calibration.cpu_model: calibration for calibration in load_registry(path, root=tmp_path).calibrations}
    assert sorted(calibrations) == sorted([CPU_MODEL, OTHER_CPU])
    assert promote(tmp_path, path, gated({GT: IMPROVED, OMEGA: DRIFT}, cpu_model=cpu), FakeCheckout()) == 0
    (link,) = load_registry(path, root=tmp_path).anchor_links
    assert link.calibration == calibrations[cpu]


def test_rebase_harness_with_a_changed_workload_set_keeps_the_chain_and_reanchors(tmp_path: Path) -> None:
    """Spec: a rebase keeps the chain; an added workload starts at 1, a removed one drops its anchor."""
    path = promoted_once(tmp_path)
    calibrate_at(tmp_path, path, HEAD)
    before = load_registry(path, root=tmp_path)
    rebased = aa_records(tmp_path / "records", series="rebased", workloads=(OMEGA, NEW_WORKLOAD), harness=NEW_HARNESS, ref_commit=HEAD, head_rev=HEAD)
    assert calibrate(tmp_path, path, rebased, "--rebase-harness") == 0
    after = load_registry(path, root=tmp_path)
    assert after.harness_sha256 == NEW_HARNESS and sorted(after.workloads) == sorted([OMEGA, NEW_WORKLOAD])
    assert (after.anchor_ref, after.anchor_links) == (before.anchor_ref, before.anchor_links)
    assert dict(after.anchors_exact()) == {OMEGA: Fraction(DRIFT), NEW_WORKLOAD: Fraction(1)}
