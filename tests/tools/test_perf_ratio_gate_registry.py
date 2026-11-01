"""The ratio gate registry: margins re-derived from cited A/A records, and its named refusals."""

from __future__ import annotations

import dataclasses
import hashlib
import json
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pytest
from tools.perf_ratio_gate_registry import CALIBRATION_ROOT
from tools.perf_ratio_gate_registry import CALIBRATION_SESSIONS
from tools.perf_ratio_gate_registry import SAFETY
from tools.perf_ratio_gate_registry import Calibration
from tools.perf_ratio_gate_registry import CalibrationRecordDigestError
from tools.perf_ratio_gate_registry import CalibrationRecordMissingError
from tools.perf_ratio_gate_registry import CalibrationRecordOutsideRootError
from tools.perf_ratio_gate_registry import CalibrationWorkloadMismatchError
from tools.perf_ratio_gate_registry import CitedRecord
from tools.perf_ratio_gate_registry import DuplicateCalibrationError
from tools.perf_ratio_gate_registry import DuplicateCalibrationSessionError
from tools.perf_ratio_gate_registry import DuplicateRegistryWorkloadError
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import InvalidCalibrationSessionError
from tools.perf_ratio_gate_registry import Margin
from tools.perf_ratio_gate_registry import MalformedRegistryError
from tools.perf_ratio_gate_registry import MarginMismatchError
from tools.perf_ratio_gate_registry import MissingWorkloadsError
from tools.perf_ratio_gate_registry import RegistryDigestMismatchError
from tools.perf_ratio_gate_registry import SafetyMismatchError
from tools.perf_ratio_gate_registry import TooFewCalibrationSessionsError
from tools.perf_ratio_gate_registry import UncalibratedMachineError
from tools.perf_ratio_gate_registry import UnsluggableCpuModelError
from tools.perf_ratio_gate_registry import UnsortedRegistryWorkloadsError
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibrated_slug
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import calibration_of
from tools.perf_ratio_gate_registry import derive_margins
from tools.perf_ratio_gate_registry import load_registry
from tools.perf_ratio_gate_registry import write_registry
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
from tests.tools.perf_ratio_sessions import synthetic_session
from tests.tools.perf_ratio_sessions import write_record

# (ref_cpu, head_cpu) of each A/A session's single "w" pair: ratios 1.00, 0.98, 1.01, 1.00, 1.005.
W_PAIRS = ((100.0, 100.0), (100.0, 98.0), (100.0, 101.0), (150.0, 150.0), (200.0, 201.0))
# Pairs of the second workload "v", the same in every session (two pairs, ratio 0.99).
V_PAIRS = ((50.0, 49.5), (51.0, 52.0))


def _ratios(index: int) -> dict[str, Sequence[tuple[float, float]]]:
    return {"w": (W_PAIRS[index],), "v": V_PAIRS}


def _aa_sessions(count: int = CALIBRATION_SESSIONS) -> list[PerfRatioSession]:
    return [synthetic_session(_ratios(index)) for index in range(count)]


def _write_all(root: Path, sessions: Sequence[PerfRatioSession], cpu: str = "apple-m1-max") -> list[Path]:
    return [write_record(session, root / cpu / f"aa-{index}" / "perf_ratio_session.json") for index, session in enumerate(sessions)]


def _registry(root: Path) -> GateRegistry:
    calibration = calibration_of(_write_all(root, _aa_sessions()), ref_commit=REF, harness_sha256=HARNESS, root=root)
    return build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration])


def _written(root: Path) -> Path:
    path = root / "perf_ratio_gate.json"
    write_registry(_registry(root), path, root=root)
    return path


def _content_digest(document: Mapping[str, Any]) -> str:
    content = {key: value for key, value in document.items() if key != "registry_sha256"}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _rewrite(path: Path, edit: Callable[[dict[str, Any]], None], *, redigest: bool) -> None:
    document = json.loads(path.read_text(encoding="utf-8"))
    edit(document)
    if redigest:
        document["registry_sha256"] = _content_digest(document)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _dump_unchecked(registry: GateRegistry, path: Path) -> Path:
    """Write *registry* bypassing `write_registry`'s verification, as a hand edit would."""
    path.write_text(registry.dumps(), encoding="utf-8")
    return path


def _cited(root: Path, paths: Sequence[Path]) -> tuple[CitedRecord, ...]:
    return tuple(CitedRecord(path=path.relative_to(root).as_posix(), record_sha256=sha256_file(path)) for path in paths)


def test_margins_are_safety_times_worst_aa_deviation() -> None:
    expected = SAFETY * max(abs(head / ref - 1.0) for ref, head in W_PAIRS)
    margins = derive_margins(_aa_sessions())
    assert margins["w"] == expected
    assert margins["w"] == pytest.approx(2 * 0.02)
    assert SAFETY == 2.0
    assert margins["v"] == SAFETY * abs(min(head for _, head in V_PAIRS) / min(ref for ref, _ in V_PAIRS) - 1.0)


def test_a_registry_round_trips_and_re_derives_on_load(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    path = tmp_path / "perf_ratio_gate.json"
    write_registry(registry, path, root=tmp_path)
    loaded = load_registry(path, root=tmp_path)
    assert loaded == registry
    assert loaded.registry_sha256 == registry.registry_sha256
    assert json.loads(path.read_text(encoding="utf-8"))["registry_sha256"] == registry.registry_sha256
    assert path.read_text(encoding="utf-8") == registry.dumps() == loaded.dumps()
    assert loaded.workloads == ("v", "w")
    assert loaded.safety == SAFETY
    assert loaded.margins_for(CPU_MODEL) == derive_margins(_aa_sessions())
    assert [record.path for record in loaded.calibrations[0].sessions] == [f"apple-m1-max/aa-{index}/perf_ratio_session.json" for index in range(5)]


def test_an_edited_registry_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)
    text = path.read_text(encoding="utf-8")
    margin = repr(_registry(tmp_path).margins_for(CPU_MODEL)["w"])
    assert text.count(margin) == 1
    edited = margin[:-1] + ("1" if margin[-1] != "1" else "2")
    path.write_text(text.replace(margin, edited), encoding="utf-8")
    with pytest.raises(RegistryDigestMismatchError):
        load_registry(path, root=tmp_path)


def test_a_re_digested_registry_with_a_wrong_margin_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)

    def widen(document: dict[str, Any]) -> None:
        document["calibrations"][0]["margins"]["w"] *= 2.0

    _rewrite(path, widen, redigest=True)
    with pytest.raises(MarginMismatchError):
        load_registry(path, root=tmp_path)


def test_a_re_digested_registry_with_another_safety_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)

    def loosen(document: dict[str, Any]) -> None:
        document["safety"] = 3.0

    _rewrite(path, loosen, redigest=True)
    with pytest.raises(SafetyMismatchError):
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


def _dirty(session: PerfRatioSession) -> PerfRatioSession:
    head = dataclasses.replace(session.arms[ArmName.HEAD], dirty=True, patch_sha256="p" * 64)
    return dataclasses.replace(session, arms={ArmName.REF: session.arms[ArmName.REF], ArmName.HEAD: head})


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


def test_an_uncalibrated_cpu_model_is_refused(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    with pytest.raises(UncalibratedMachineError):
        registry.margins_for("Intel Xeon")


def test_a_calibration_less_registry_round_trips_and_refuses_every_cpu(tmp_path: Path) -> None:
    registry = build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[], workloads=["w", "v"])
    path = tmp_path / "perf_ratio_gate.json"
    write_registry(registry, path, root=tmp_path)
    loaded = load_registry(path, root=tmp_path)
    assert loaded == registry
    assert loaded.calibrations == ()
    assert loaded.workloads == ("v", "w")
    with pytest.raises(UncalibratedMachineError):
        loaded.margins_for(CPU_MODEL)


def test_a_calibration_less_registry_needs_its_workloads() -> None:
    with pytest.raises(MissingWorkloadsError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[])
    with pytest.raises(MissingWorkloadsError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[], workloads=[])


def test_repeated_registry_workloads_are_refused() -> None:
    with pytest.raises(DuplicateRegistryWorkloadError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[], workloads=["w", "v", "w"])


def test_given_workloads_must_agree_with_the_calibrations(tmp_path: Path) -> None:
    calibration = _registry(tmp_path).calibrations[0]
    assert build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration], workloads=["w", "v"]).workloads == ("v", "w")
    with pytest.raises(CalibrationWorkloadMismatchError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration], workloads=["w"])


def test_calibrations_over_different_workloads_are_refused(tmp_path: Path) -> None:
    calibration = _registry(tmp_path).calibrations[0]
    other = Calibration(cpu_model="Intel Xeon", sessions=calibration.sessions, margins={"w": calibration.margins["w"]})
    with pytest.raises(CalibrationWorkloadMismatchError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration, other])


def test_two_calibrations_for_one_cpu_are_refused(tmp_path: Path) -> None:
    calibration = _registry(tmp_path).calibrations[0]
    with pytest.raises(DuplicateCalibrationError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=[calibration, calibration])


def test_a_registry_with_two_calibrations_for_one_cpu_is_refused_on_load(tmp_path: Path) -> None:
    calibration = _registry(tmp_path).calibrations[0]
    registry = GateRegistry(ref_commit=REF, harness_sha256=HARNESS, safety=SAFETY, workloads=("v", "w"), calibrations=(calibration, calibration))
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(DuplicateCalibrationError):
        load_registry(path, root=tmp_path)


def test_write_registry_refuses_a_registry_that_would_not_load(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    calibration = registry.calibrations[0]
    wrong = dataclasses.replace(calibration, margins={name: Margin(margin * 2.0) for name, margin in calibration.margins.items()})
    path = tmp_path / "perf_ratio_gate.json"
    with pytest.raises(MarginMismatchError):
        write_registry(dataclasses.replace(registry, calibrations=(wrong,)), path, root=tmp_path)
    assert not path.exists()


# A calibration-less registry's workload set, each shape the shared check refuses, and the refusal.
BAD_WORKLOAD_SHAPES: Mapping[str, tuple[tuple[str, ...], type[Exception]]] = {
    "unsorted": (("w", "v"), UnsortedRegistryWorkloadsError),
    "repeated": (("v", "v", "w"), DuplicateRegistryWorkloadError),
    "empty": ((), MissingWorkloadsError),
}


@pytest.mark.parametrize("shape", sorted(BAD_WORKLOAD_SHAPES))
def test_write_registry_refuses_a_badly_shaped_workload_set(tmp_path: Path, shape: str) -> None:
    workloads, error = BAD_WORKLOAD_SHAPES[shape]
    registry = GateRegistry(ref_commit=REF, harness_sha256=HARNESS, safety=SAFETY, workloads=workloads, calibrations=())
    path = tmp_path / "perf_ratio_gate.json"
    with pytest.raises(error):
        write_registry(registry, path, root=tmp_path)
    assert not path.exists()
    assert not path.with_name(f"{path.name}.writing").exists()


@pytest.mark.parametrize("shape", sorted(BAD_WORKLOAD_SHAPES))
def test_a_registry_with_a_badly_shaped_workload_set_is_refused_on_load(tmp_path: Path, shape: str) -> None:
    workloads, error = BAD_WORKLOAD_SHAPES[shape]
    registry = GateRegistry(ref_commit=REF, harness_sha256=HARNESS, safety=SAFETY, workloads=workloads, calibrations=())
    path = _dump_unchecked(registry, tmp_path / "perf_ratio_gate.json")
    with pytest.raises(error):
        load_registry(path, root=tmp_path)


MALFORMED_EDITS: Mapping[str, Callable[[dict[str, Any]], None]] = {
    "schema": lambda document: document.update(schema="perf-ratio-gate/0"),
    "extra-key": lambda document: document.update(extra=1),
    "missing-key": lambda document: document.pop("ref_commit"),
    "safety-string": lambda document: document.update(safety="2.0"),
    "margin-not-a-number": lambda document: document["calibrations"][0]["margins"].update(w=True),
    "cited-not-an-object": lambda document: document["calibrations"][0]["sessions"].append("x"),
}


@pytest.mark.parametrize("edit", sorted(MALFORMED_EDITS))
def test_a_malformed_registry_is_refused(tmp_path: Path, edit: str) -> None:
    path = _written(tmp_path)
    _rewrite(path, MALFORMED_EDITS[edit], redigest=True)
    with pytest.raises(MalformedRegistryError):
        load_registry(path, root=tmp_path)


def test_a_registry_that_is_not_json_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "perf_ratio_gate.json"
    path.write_text("{", encoding="utf-8")
    with pytest.raises(MalformedRegistryError):
        load_registry(path, root=tmp_path)


def test_calibration_dir_names_ref_and_cpu() -> None:
    assert calibration_dir(REF, "Apple M1 Max") == CALIBRATION_ROOT / "15fbed70" / "apple-m1-max"


def test_calibrated_slug_collapses_runs_of_other_characters() -> None:
    assert calibrated_slug("Apple M1 Max") == "apple-m1-max"
    assert calibrated_slug("  Intel(R) Xeon(R) CPU @ 2.20GHz ") == "intel-r-xeon-r-cpu-2-20ghz"
    with pytest.raises(UnsluggableCpuModelError):
        calibrated_slug(" () ")


def test_the_committed_calibration_records_calibrate() -> None:
    paths = sorted(calibration_dir(REF, CPU_MODEL).glob("aa-*/perf_ratio_session.json"))
    assert len(paths) == CALIBRATION_SESSIONS
    harness = load_session(paths[0]).harness_sha256
    calibration = calibration_of(paths, ref_commit=REF, harness_sha256=harness)
    assert calibration.cpu_model == CPU_MODEL
    assert calibration.margins == derive_margins([load_session(path) for path in paths])
    assert all(not record.path.startswith("/") for record in calibration.sessions)
