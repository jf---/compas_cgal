"""The ratio gate registry file: round-trip, self-digest, schema, shape and write refusals."""

from __future__ import annotations

import dataclasses
import json
from collections.abc import Callable
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest
from tools.perf_ratio_gate_registry import CALIBRATION_ROOT
from tools.perf_ratio_gate_registry import CALIBRATE_TASK
from tools.perf_ratio_gate_registry import REGISTRY_SCHEMA
from tools.perf_ratio_gate_registry import SAFETY
from tools.perf_ratio_gate_registry import Calibration
from tools.perf_ratio_gate_registry import CalibrationWorkloadMismatchError
from tools.perf_ratio_gate_registry import DuplicateCalibrationError
from tools.perf_ratio_gate_registry import DuplicateRegistryWorkloadError
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import Margin
from tools.perf_ratio_gate_registry import MalformedRegistryError
from tools.perf_ratio_gate_registry import MarginMismatchError
from tools.perf_ratio_gate_registry import MIGRATION_COMMAND
from tools.perf_ratio_gate_registry import MissingWorkloadsError
from tools.perf_ratio_gate_registry import RegistryDigestMismatchError
from tools.perf_ratio_gate_registry import SafetyMismatchError
from tools.perf_ratio_gate_registry import SupersededRegistrySchemaError
from tools.perf_ratio_gate_registry import UncalibratedMachineError
from tools.perf_ratio_gate_registry import UnsluggableCpuModelError
from tools.perf_ratio_gate_registry import UnsortedRegistryWorkloadsError
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibrated_slug
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import derive_margins
from tools.perf_ratio_gate_registry import load_registry
from tools.perf_ratio_gate_registry import write_registry

from tests.tools.perf_ratio_sessions import CPU_MODEL
from tests.tools.perf_ratio_sessions import HARNESS
from tests.tools.perf_ratio_sessions import REF
from tests.tools.perf_ratio_sessions import _aa_sessions
from tests.tools.perf_ratio_sessions import _registry
from tests.tools.perf_ratio_sessions import _written
from tests.tools.perf_ratio_sessions import _rewrite
from tests.tools.perf_ratio_sessions import _dump_unchecked


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


def test_a_re_digested_registry_with_another_safety_is_refused(tmp_path: Path) -> None:
    path = _written(tmp_path)

    def loosen(document: dict[str, Any]) -> None:
        document["safety"] = 3.0

    _rewrite(path, loosen, redigest=True)
    with pytest.raises(SafetyMismatchError):
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


def test_the_schema_is_two() -> None:
    assert REGISTRY_SCHEMA == "perf-ratio-gate/2"


def test_a_superseded_schema_one_registry_names_the_migration(tmp_path: Path) -> None:
    path = _written(tmp_path)

    def downgrade(document: dict[str, Any]) -> None:
        document["schema"] = "perf-ratio-gate/1"
        for key in ("anchor_ref", "anchor_links", "anchors"):
            del document[key]

    _rewrite(path, downgrade, redigest=True)
    with pytest.raises(SupersededRegistrySchemaError, match=r"perf-ratio-calibrate --records") as caught:
        load_registry(path, root=tmp_path)
    # One calibrate call per machine: calibration_of refuses records from mixed cpu_models.
    assert "per machine" in str(caught.value)
    assert MIGRATION_COMMAND in str(caught.value)


def test_the_migration_command_is_the_calibrate_task() -> None:
    assert MIGRATION_COMMAND == f"{CALIBRATE_TASK} --records"
