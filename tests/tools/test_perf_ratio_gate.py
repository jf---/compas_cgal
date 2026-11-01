"""The ratio gate CLI: check, calibrate, promote.

Every session is synthetic (`synthetic_session`), every record lives under ``tmp_path``
(``root=tmp_path``), the `Measure` is a fake and promotion's checkout is a fake: no
timing, no arms, no harness subprocess.
"""

from __future__ import annotations

import dataclasses
import subprocess
import sys
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
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
from tools.perf_ratio_gate import PromotionRefusedError
from tools.perf_ratio_gate import RebaseWithoutHarnessChangeError
from tools.perf_ratio_gate import harness_command
from tools.perf_ratio_gate import main
from tools.perf_ratio_gate import subprocess_measure
from tools.perf_ratio_gate_registry import CALIBRATION_SESSIONS
from tools.perf_ratio_gate_registry import InvalidCalibrationSessionError
from tools.perf_ratio_gate_registry import TooFewCalibrationSessionsError
from tools.perf_ratio_gate_registry import UncalibratedMachineError
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import derive_margins
from tools.perf_ratio_gate_registry import load_registry
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
        f"{GT}: PASS rho 1.000 margin {margin:.3f} pairs {DEFAULT_PAIRS} changed output",
        f"{OMEGA}: PASS rho 1.010 margin {margin:.3f} pairs {DEFAULT_PAIRS} changed output",
        "ratio gate: 2 workloads: 2 PASS, 0 IMPROVED, 0 REGRESSION",
    ]


def test_check_regression_raises_with_workload_rho_margin_pairs(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    capsys.readouterr()
    measure = FakeMeasure(REGRESSED)
    with pytest.raises(RatioRegressionError, match=rf"{OMEGA}: rho_w 1\.2 > 1 \+ m_w 0\.04\d* over 20 pairs") as raised:
        check(tmp_path, registry, gated({GT: IMPROVED, OMEGA: REGRESSED}), measure=measure)
    assert [verdict.workload for verdict in raised.value.report.regressions] == [OMEGA]
    assert len(measure.requests) == 3
    out = capsys.readouterr().out
    assert f"{OMEGA}: REGRESSION rho 1.200 margin 0.040 pairs 20\n" in out
    assert "ratio gate: 2 workloads: 0 PASS, 1 IMPROVED, 1 REGRESSION" in out
    assert "promotion available" not in out


def test_check_improved_prints_the_promotion_hint(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    capsys.readouterr()
    record = write_record(gated({GT: IMPROVED, OMEGA: PASSING}), tmp_path / "gate" / "perf_ratio_session.json")
    assert main(["check", str(record), "--registry", str(registry)], root=tmp_path, measure=no_measure, checkout=no_checkout) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == f"{GT}: IMPROVED rho 0.900 margin 0.040 pairs {DEFAULT_PAIRS}"
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
    assert main(argv, root=tmp_path, measure=measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL) == 0
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
    assert main(["calibrate", "--sessions", str(sessions), "--registry", str(path)], root=tmp_path, measure=measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL) == 0
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
        main(["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(path)], root=tmp_path, measure=no_measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL)
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
        )


def test_calibrate_sessions_refuses_a_session_from_another_machine(tmp_path: Path) -> None:
    path = uncalibrated_registry(tmp_path)
    before = path.read_bytes()
    measure = FakeMeasure(1.0, cpu_model=OTHER_CPU)
    with pytest.raises(CalibrationMachineMismatchError):
        main(["calibrate", "--sessions", str(CALIBRATION_SESSIONS), "--registry", str(path)], root=tmp_path, measure=measure, checkout=no_checkout, cpu_model=lambda: CPU_MODEL)
    assert len(measure.requests) == 1
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
        main(argv, root=tmp_path, measure=failing, checkout=no_checkout, cpu_model=lambda: CPU_MODEL)
    assert len(failing.requests) == 2
    assert not sessions_target(tmp_path).exists()
    assert not (tmp_path / calibration_dir(REF, CPU_MODEL).relative_to(REPOSITORY).parents[1]).exists()
    assert path.read_bytes() == before
    rerun = FakeMeasure(*AA_HEAD_CPU)
    assert main(argv, root=tmp_path, measure=rerun, checkout=no_checkout, cpu_model=lambda: CPU_MODEL) == 0
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
