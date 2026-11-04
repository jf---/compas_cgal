"""The ratio gate's A/A calibration citations: margins re-derived from cited records, and their named refusals."""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pytest
from tools.perf_ratio_gate_citations import CALIBRATION_SESSIONS
from tools.perf_ratio_gate_citations import SAFETY
from tools.perf_ratio_gate_citations import Calibration
from tools.perf_ratio_gate_citations import CalibrationRecordDigestError
from tools.perf_ratio_gate_citations import CalibrationRecordMissingError
from tools.perf_ratio_gate_citations import CalibrationRecordOutsideRootError
from tools.perf_ratio_gate_citations import CitedRecord
from tools.perf_ratio_gate_citations import DuplicateCalibrationSessionError
from tools.perf_ratio_gate_citations import InvalidCalibrationSessionError
from tools.perf_ratio_gate_citations import MarginMismatchError
from tools.perf_ratio_gate_citations import TooFewCalibrationSessionsError
from tools.perf_ratio_gate_citations import calibration_of
from tools.perf_ratio_gate_citations import derive_margins
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import Margin
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import load_registry
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session
from tools.run_provenance import sha256_file

from tests.tools.perf_ratio_sessions import CPU_MODEL
from tests.tools.perf_ratio_sessions import HARNESS
from tests.tools.perf_ratio_sessions import OTHER_REV
from tests.tools.perf_ratio_sessions import REF
from tests.tools.perf_ratio_sessions import V_PAIRS
from tests.tools.perf_ratio_sessions import W_PAIRS
from tests.tools.perf_ratio_sessions import _aa_sessions
from tests.tools.perf_ratio_sessions import _dirty
from tests.tools.perf_ratio_sessions import _dump_unchecked
from tests.tools.perf_ratio_sessions import _ratios
from tests.tools.perf_ratio_sessions import _rewrite
from tests.tools.perf_ratio_sessions import _write_all
from tests.tools.perf_ratio_sessions import _written
from tests.tools.perf_ratio_sessions import cite
from tests.tools.perf_ratio_sessions import synthetic_session
from tests.tools.perf_ratio_sessions import write_record


# A harness digest other than HARNESS, for a calibration session that ran another harness.
OLD_HARNESS = "g" * 64


def _cited(root: Path, paths: Sequence[Path]) -> tuple[CitedRecord, ...]:
    return tuple(CitedRecord(path=path.relative_to(root).as_posix(), record_sha256=sha256_file(path)) for path in paths)


def test_margins_are_safety_times_worst_aa_deviation() -> None:
    expected = SAFETY * max(abs(head / ref - 1.0) for ref, head in W_PAIRS)
    margins = derive_margins(_aa_sessions())
    assert margins["w"] == expected
    assert margins["w"] == pytest.approx(2 * 0.02)
    assert SAFETY == 2.0
    assert margins["v"] == SAFETY * abs(min(head for _, head in V_PAIRS) / min(ref for ref, _ in V_PAIRS) - 1.0)


def test_a_re_digested_registry_with_a_wrong_margin_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)

    def widen(document: dict[str, Any]) -> None:
        document["calibrations"][0]["margins"]["w"] *= 2.0

    _rewrite(path, widen, redigest=True)
    with pytest.raises(MarginMismatchError):
        load_registry(path, root=tmp_path)


def test_a_regenerated_calibration_record_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)
    regenerated = synthetic_session({"w": ((100.0, 100.5),), "v": V_PAIRS})
    write_record(regenerated, tmp_path / "apple-m1-max" / "aa-1" / "perf_ratio_session.json")
    with pytest.raises(CalibrationRecordDigestError):
        load_registry(path, root=tmp_path)


def test_a_deleted_calibration_record_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)
    (tmp_path / "apple-m1-max" / "aa-2" / "perf_ratio_session.json").unlink()
    with pytest.raises(CalibrationRecordMissingError):
        load_registry(path, root=tmp_path)


INVALID_SESSIONS: Mapping[str, Callable[[], PerfRatioSession]] = {
    "injected": lambda: synthetic_session(_ratios(4), injected=Injection(arm=ArmName.HEAD, seconds=CpuSeconds(1.0))),
    "head-at-another-rev": lambda: synthetic_session(_ratios(4), head_rev=OTHER_REV),
    "ref-at-another-rev": lambda: synthetic_session(_ratios(4), ref_commit=OTHER_REV, head_rev=OTHER_REV),
    "dirty": lambda: _dirty(synthetic_session(_ratios(4))),
    "harness-differs": lambda: synthetic_session(_ratios(4), harness="g" * 64),
    "cpu-differs": lambda: synthetic_session(_ratios(4), cpu_model="Intel Xeon"),
    "same-output-false": lambda: synthetic_session(_ratios(4), head_digest="e"),
    "missing-workload": lambda: synthetic_session({"w": (W_PAIRS[4],)}),
}


@pytest.mark.parametrize("reason", sorted(INVALID_SESSIONS))
def test_an_invalid_calibration_session_is_refused(tmp_path: Path, reason: str) -> None:
    sessions = [*_aa_sessions(CALIBRATION_SESSIONS - 1), INVALID_SESSIONS[reason]()]
    with pytest.raises(InvalidCalibrationSessionError):
        calibration_of(_write_all(tmp_path, sessions), ref_commit=REF, harness_sha256=HARNESS, root=tmp_path)


@pytest.mark.parametrize("reason", sorted(INVALID_SESSIONS))
def test_a_registry_citing_an_invalid_session_is_refused_on_load(tmp_path: Path, reason: str) -> None:
    sessions = [*_aa_sessions(CALIBRATION_SESSIONS - 1), INVALID_SESSIONS[reason]()]
    paths = _write_all(tmp_path, sessions)
    workloads = ("v", "w")
    margins = {name: Margin(1.0) for name in workloads}
    calibration = Calibration(cpu_model=CPU_MODEL, sessions=_cited(tmp_path, paths), margins=margins)
    registry = GateRegistry(ref_commit=REF, harness_sha256=HARNESS, safety=SAFETY, workloads=workloads, calibrations=(calibration,))
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(InvalidCalibrationSessionError):
        load_registry(path, root=tmp_path)


def test_a_registry_whose_sessions_all_miss_a_workload_is_refused_on_load(tmp_path: Path) -> None:
    sessions = [synthetic_session({"w": (pair,)}) for pair in W_PAIRS]
    paths = _write_all(tmp_path, sessions)
    calibration = Calibration(cpu_model=CPU_MODEL, sessions=_cited(tmp_path, paths), margins=derive_margins(sessions))
    registry = GateRegistry(ref_commit=REF, harness_sha256=HARNESS, safety=SAFETY, workloads=("v", "w"), calibrations=(calibration,))
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(InvalidCalibrationSessionError):
        load_registry(path, root=tmp_path)


def test_derive_margins_refuses_sessions_over_different_workloads() -> None:
    sessions = [*_aa_sessions(CALIBRATION_SESSIONS - 1), synthetic_session({"w": (W_PAIRS[4],)})]
    with pytest.raises(InvalidCalibrationSessionError):
        derive_margins(sessions)


def test_fewer_than_k_sessions_are_refused(tmp_path: Path) -> None:
    paths = _write_all(tmp_path, _aa_sessions(CALIBRATION_SESSIONS - 1))
    with pytest.raises(TooFewCalibrationSessionsError):
        calibration_of(paths, ref_commit=REF, harness_sha256=HARNESS, root=tmp_path)


def test_a_registry_citing_fewer_than_k_sessions_is_refused_on_load(tmp_path: Path) -> None:
    sessions = _aa_sessions(CALIBRATION_SESSIONS - 1)
    paths = _write_all(tmp_path, sessions)
    calibration = Calibration(cpu_model=CPU_MODEL, sessions=_cited(tmp_path, paths), margins=derive_margins(sessions))
    registry = build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration])
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(TooFewCalibrationSessionsError):
        load_registry(path, root=tmp_path)


def test_derive_margins_refuses_no_sessions() -> None:
    with pytest.raises(TooFewCalibrationSessionsError):
        derive_margins([])


def test_one_record_cited_twice_is_refused(tmp_path: Path) -> None:
    paths = _write_all(tmp_path, _aa_sessions())
    with pytest.raises(DuplicateCalibrationSessionError):
        calibration_of([*paths[:-1], paths[0]], ref_commit=REF, harness_sha256=HARNESS, root=tmp_path)


def test_a_byte_copy_of_one_record_is_refused(tmp_path: Path) -> None:
    paths = _write_all(tmp_path, _aa_sessions())
    paths[-1].write_bytes(paths[0].read_bytes())
    with pytest.raises(DuplicateCalibrationSessionError):
        calibration_of(paths, ref_commit=REF, harness_sha256=HARNESS, root=tmp_path)


def test_a_registry_citing_one_record_twice_is_refused_on_load(tmp_path: Path) -> None:
    sessions = _aa_sessions()
    paths = _write_all(tmp_path, sessions)
    cited = _cited(tmp_path, [*paths[:-1], paths[0]])
    calibration = Calibration(cpu_model=CPU_MODEL, sessions=cited, margins=derive_margins([*sessions[:-1], sessions[0]]))
    registry = build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration])
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(DuplicateCalibrationSessionError):
        load_registry(path, root=tmp_path)


def test_a_record_outside_the_root_is_refused(tmp_path: Path) -> None:
    paths = _write_all(tmp_path / "outside", _aa_sessions())
    with pytest.raises(CalibrationRecordOutsideRootError):
        calibration_of(paths, ref_commit=REF, harness_sha256=HARNESS, root=tmp_path / "root")


def test_a_registry_citing_a_path_outside_the_root_is_refused_on_load(tmp_path: Path) -> None:
    root = tmp_path / "root"
    path = _written(root)

    def escape(document: dict[str, Any]) -> None:
        document["calibrations"][0]["sessions"][0]["path"] = "../escaped/perf_ratio_session.json"

    _rewrite(path, escape, redigest=True)
    with pytest.raises(CalibrationRecordOutsideRootError):
        load_registry(path, root=root)


def test_a_registry_citing_an_absolute_path_is_refused_on_load(tmp_path: Path) -> None:
    path = _written(tmp_path)

    def absolute(document: dict[str, Any]) -> None:
        cited = document["calibrations"][0]["sessions"][0]
        cited["path"] = (tmp_path / cited["path"]).resolve().as_posix()

    _rewrite(path, absolute, redigest=True)
    with pytest.raises(CalibrationRecordOutsideRootError):
        load_registry(path, root=tmp_path)


def test_a_registry_calibration_harness_refusal_names_the_registrys_harness(tmp_path: Path) -> None:
    sessions = [*_aa_sessions(CALIBRATION_SESSIONS - 1), synthetic_session(_ratios(4), harness=OLD_HARNESS)]
    paths = _write_all(tmp_path, sessions)
    with pytest.raises(InvalidCalibrationSessionError, match=rf"not the registry's {HARNESS}"):
        calibration_of(paths, ref_commit=REF, harness_sha256=HARNESS, root=tmp_path)
    calibration = Calibration(cpu_model=CPU_MODEL, sessions=tuple(cite(tmp_path, path) for path in paths), margins=derive_margins(sessions))
    registry = GateRegistry(ref_commit=REF, harness_sha256=HARNESS, safety=SAFETY, workloads=("v", "w"), calibrations=(calibration,))
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(InvalidCalibrationSessionError, match=rf"not the registry's {HARNESS}"):
        load_registry(path, root=tmp_path)


def test_the_committed_calibration_records_calibrate() -> None:
    paths = sorted(calibration_dir(REF, CPU_MODEL).glob("aa-*/perf_ratio_session.json"))
    assert len(paths) == CALIBRATION_SESSIONS
    harness = load_session(paths[0]).harness_sha256
    calibration = calibration_of(paths, ref_commit=REF, harness_sha256=harness)
    assert calibration.cpu_model == CPU_MODEL
    assert calibration.margins == derive_margins([load_session(path) for path in paths])
    assert all(not record.path.startswith("/") for record in calibration.sessions)
