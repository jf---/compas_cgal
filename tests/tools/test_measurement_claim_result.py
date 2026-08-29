from __future__ import annotations

# ruff: noqa: E501  # Strict-mypy fixture literals stay whole and independently auditable.

import copy
import datetime
import hashlib
import importlib
import importlib.util
import json
import pathlib
import subprocess
import sys
from typing import Any
from typing import Callable
from typing import Optional

import pytest


UTC = datetime.timezone.utc
EXTRACTION_COMMIT = "eec665c1df1cd8d1e98dd9dd1001b5984e17a703"
HISTORY_COMMIT = "29050b01e656ea7bf577b18f7bb50a04ff9a23c9"
CASE_ORDER = [
    "radial-station",
    "radial-subdivisions",
    "radial-floor",
    "radial-margin",
    "advance-placement",
    "advance-probe-count",
]


def _artifact() -> Any:
    return importlib.import_module("tools.measurement_artifact")


def _claims() -> Any:
    return importlib.import_module("tools.measurement_claim_result")


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    return subprocess.run(["git", "-C", str(repository), *arguments], check=True, capture_output=True).stdout


def _repository(tmp_path: pathlib.Path) -> pathlib.Path:
    repository = tmp_path / "repository"
    repository.mkdir()
    _git(repository, "init", "-q")
    (repository / "pixi.lock").write_bytes(b"claim lock\n")
    _git(repository, "add", "pixi.lock")
    _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "lock")
    return repository.resolve()


def _public_call(width: float, height: float, cap: float, *, climb: bool = True) -> dict[str, object]:
    return {
        "polygon": [[0.0, 0.0, 0.0], [width, 0.0, 0.0], [width, height, 0.0], [0.0, height, 0.0]],
        "holes": [],
        "tool_diameter": 2.0,
        "tea_cap": cap,
        "guide_step": 0.025,
        "max_advance": 1.0,
        "radial_clearance": 0.002,
        "climb": climb,
        "cut_z": 0.0,
        "clearance_z": 2.0,
        "max_passes": 1000,
        "samples_per_radian": 10.0,
    }


def _radial_provenance() -> dict[str, object]:
    return {
        "policy": "radial-known-reporting-driven-selection/v1",
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": [
            {
                "symbol": "compas_cgal.engagement_radial_toolpath._least_bad_rung",
                "value_source": "compas_cgal._stock_2.engagement_at[1]",
                "effect": "forced-radius-selection",
            },
            {
                "symbol": "compas_cgal.engagement_radial_toolpath._largest_admissible_radius",
                "value_source": "compas_cgal.engagement_radial_toolpath._GentlestRung.peak",
                "effect": "refined-scan-control",
            },
        ],
    }


def _advance_provenance() -> dict[str, object]:
    return {
        "policy": "advance-native-cap-reporting-observation/v1",
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": [],
    }


def _counts(observations: int, exceeded: int) -> dict[str, int]:
    return {"observations": observations, "accepted": observations - exceeded, "exceeded": exceeded}


def _radial_audit() -> dict[str, object]:
    return {
        "phase": "entry-angle",
        "probe_offsets": [360.0 * index / 16 for index in range(16)],
        "includes_entry_phase": True,
        "adds_separate_entry_probe": False,
        "excludes_chain_entry_circles": True,
    }


def _payload(source_commit: str) -> dict[str, object]:
    generator_angles = [360.0 * index / 32 for index in range(32)]
    radial = _radial_provenance()
    advance = _advance_provenance()
    station_reporting = {
        "maximal_radius": 0.5156,
        "rung_6_radius": 0.2156,
        "rung_6_peak": 61.3,
        "rung_7_radius": 0.1656,
        "rung_7_peak": 5.9,
        "refined_band_min_radius": 0.1719,
        "refined_band_max_radius": 0.2123,
        "forced_peak": 107.0,
        "rescued_peak": 59.0,
        "angle_unit": "degree",
        "length_unit": "mm",
    }
    subdivision_reporting = [
        {"subdivisions": value, "worst_peak": 88.6, "circles_over_cap": count, "angle_unit": "degree"} for value, count in zip([1, 2, 4, 8, 16], [12, 12, 8, 8, 8])
    ]
    floor_reporting = [{"floor_steps": value, "worst_peak": 88.6, "circles_over_cap": count, "angle_unit": "degree"} for value, count in zip([0.25, 0.5, 1.0], [8, 8, 12])]
    margin_values = [1.25, 1.4, 1.5, 1.75, 2.0]
    margin_reporting = [
        {
            "refinement_margin": value,
            "worst_peak": peak,
            "circles_over_cap": count,
            "cutting_length": length,
            "angle_unit": "degree",
            "length_unit": "mm",
        }
        for value, peak, count, length in zip(margin_values, [54.5, 61.1, 76.3, 76.3, 110.0], [40, 23, 18, 17, 11], [413.0, 593.0, 661.0, 704.0, 835.0])
    ]
    placement_keys = [(80.0, True), (80.0, False), (100.0, True), (100.0, False)]
    placement_reporting = [
        {
            "cap": cap,
            "climb": climb,
            "selected_circles": 1,
            "positions_over_cap": 1,
            "worst_peak": 101.7,
            "worst_peak_offset": -135.0 if climb else 135.0,
            "old_probe_peak": 39.1,
            "offset_bin_counts": [
                {
                    "offset": angle if angle <= 180.0 else angle - 360.0,
                    "positions_over_cap": 1 if index == 0 else 0,
                    "angle_unit": "degree",
                }
                for index, angle in enumerate(generator_angles)
            ],
            "angle_unit": "degree",
        }
        for cap, climb in placement_keys
    ]
    labels = [("3-old", 3), ("8", 8), ("12", 12), ("16", 16), ("24", 24), ("32", 32), ("40", 40), ("48", 48)]
    count_reporting = [
        {"label": label, "probe_count": count, "cap": cap, "worst_peak": 43.4 if cap == 40.0 else 81.7, "angle_unit": "degree"} for label, count in labels for cap in [40.0, 80.0]
    ]
    cases: list[dict[str, object]] = [
        {
            "case": "radial-station",
            "source_claim_ids": ["MC-001", "MC-002", "MC-003"],
            "config": {
                "public_call": _public_call(20.0, 12.0, 60.0),
                "subdivisions": 8,
                "floor_steps": 0.5,
                "refinement_margin": 1.4,
                "generator_probe_angles": generator_angles,
                "ladder_step": 0.025,
                "ladder_span": 1.0,
                "ladder_rungs": 40,
                "max_radial_sweeps": 40,
            },
            "reconstruction": {
                "stock_model": "generator-faithful-radial-pre-bridge/v1",
                "target_centre": [18.482, 10.482],
                "centre_decimal_places": 3,
                "target_maximal_radius": 0.5156,
                "radius_decimal_places": 4,
                "occurrence_count": 1,
                "coarse_step": 0.05,
                "coarse_rungs": [6, 7],
                "refined_radius_sequence": [0.2123, 0.1719],
            },
            "native_sampled_decisions": {
                "rung_6_cap_exceeded": True,
                "rung_6_cuts_material": True,
                "rung_7_cap_exceeded": False,
                "rung_7_cuts_material": False,
                "refined_candidates": [
                    {"radius": 0.2123, "cap_exceeded": False, "cuts_material": True},
                    {"radius": 0.1719, "cap_exceeded": False, "cuts_material": True},
                ],
            },
            "reporting_values": station_reporting,
            "selection_decision_provenance": radial,
            "continuous_certificate": None,
        },
        {
            "case": "radial-subdivisions",
            "source_claim_ids": ["MC-004", "MC-005", "MC-006"],
            "config": {
                "public_call": _public_call(20.0, 12.0, 60.0),
                "subdivision_values": [1, 2, 4, 8, 16],
                "floor_steps": 0.5,
                "refinement_margin": 1.4,
                "generator_probe_angles": generator_angles,
                "audit": _radial_audit(),
            },
            "reconstruction": {
                "stock_model": "generator-faithful-radial-pre-bridge/v1",
                "audit_position_count": 16,
                "audit_includes_entry_phase": True,
                "audit_adds_separate_entry_probe": False,
                "excludes_chain_entry_circles": True,
                "non_entry_circle_counts": [244, 244, 244, 244, 244],
            },
            "native_sampled_decisions": [
                {"subdivisions": value, "selected_circles": 244, "observations": _counts(3904, count)} for value, count in zip([1, 2, 4, 8, 16], [12, 12, 8, 8, 8])
            ],
            "reporting_values": subdivision_reporting,
            "selection_decision_provenance": radial,
            "continuous_certificate": None,
        },
        {
            "case": "radial-floor",
            "source_claim_ids": ["MC-007"],
            "config": {
                "public_call": _public_call(20.0, 12.0, 60.0),
                "floor_values": [0.25, 0.5, 1.0],
                "subdivisions": 8,
                "refinement_margin": 1.4,
                "generator_probe_angles": generator_angles,
                "audit": _radial_audit(),
            },
            "reconstruction": {
                "stock_model": "generator-faithful-radial-pre-bridge/v1",
                "audit_position_count": 16,
                "audit_includes_entry_phase": True,
                "audit_adds_separate_entry_probe": False,
                "excludes_chain_entry_circles": True,
                "non_entry_circle_counts": [240, 241, 242],
            },
            "native_sampled_decisions": [
                {"floor_steps": value, "selected_circles": selected, "observations": _counts(selected * 16, count)}
                for value, selected, count in zip([0.25, 0.5, 1.0], [240, 241, 242], [8, 8, 12])
            ],
            "reporting_values": floor_reporting,
            "selection_decision_provenance": radial,
            "continuous_certificate": None,
        },
        {
            "case": "radial-margin",
            "source_claim_ids": ["MC-008"],
            "config": {
                "public_call": _public_call(6.0, 4.0, 40.0),
                "margin_values": margin_values,
                "subdivisions": 8,
                "floor_steps": 0.5,
                "generator_probe_angles": generator_angles,
                "audit": _radial_audit(),
            },
            "reconstruction": {
                "stock_model": "generator-faithful-radial-pre-bridge/v1",
                "audit_position_count": 16,
                "audit_includes_entry_phase": True,
                "audit_adds_separate_entry_probe": False,
                "excludes_chain_entry_circles": True,
                "non_entry_circle_counts": [100, 101, 102, 103, 104],
            },
            "native_sampled_decisions": [
                {"refinement_margin": value, "selected_circles": selected, "observations": _counts(selected * 16, count)}
                for value, selected, count in zip(margin_values, [100, 101, 102, 103, 104], [40, 23, 18, 17, 11])
            ],
            "reporting_values": margin_reporting,
            "selection_decision_provenance": radial,
            "continuous_certificate": None,
        },
        {
            "case": "advance-placement",
            "source_claim_ids": ["MC-009"],
            "config": {
                "public_calls": [_public_call(20.0, 12.0, cap, climb=climb) for cap, climb in placement_keys],
                "generator_probe_angles": [-60.0, 0.0, 60.0],
                "generator_prepends_entry_probe": True,
                "audit": {
                    "phase": "advance-direction",
                    "probe_offsets": generator_angles,
                    "excludes_entry_probe": True,
                },
            },
            "reconstruction": {
                "stock_model": "generator-pre-bridge/v1",
                "selected_circle_policy": "accepted-non-forced/v1",
                "original_calls_per_wrapper": 1,
                "generator_prepends_entry_probe": True,
                "audit_excludes_entry_probe": True,
                "selected_circle_count": 4,
            },
            "native_sampled_decisions": [{"cap": cap, "climb": climb, "selected_circles": 1, "observations": _counts(32, 1)} for cap, climb in placement_keys],
            "reporting_values": placement_reporting,
            "selection_decision_provenance": advance,
            "continuous_certificate": None,
        },
        {
            "case": "advance-probe-count",
            "source_claim_ids": ["MC-010"],
            "config": {
                "public_calls": [_public_call(20.0, 12.0, cap) for cap in [40.0, 80.0]],
                "generator_configurations": [
                    {
                        "label": label,
                        "probe_count": count,
                        "probe_angles": [-60.0, 0.0, 60.0] if label == "3-old" else [360.0 * index / count for index in range(count)],
                    }
                    for label, count in labels
                ],
                "generator_prepends_entry_probe": True,
                "audit": {
                    "phase": "advance-half-step",
                    "probe_offsets": [360.0 * (index + 0.5) / 60 for index in range(60)],
                    "excludes_entry_probe": True,
                },
            },
            "reconstruction": {
                "stock_model": "generator-pre-bridge/v1",
                "selected_circle_policy": "accepted-non-forced/v1",
                "original_calls_per_wrapper": 1,
                "generator_prepends_entry_probe": True,
                "audit_excludes_entry_probe": True,
                "selected_circle_count": 16,
            },
            "native_sampled_decisions": [
                {"label": label, "probe_count": count, "cap": cap, "selected_circles": 1, "observations": _counts(60, 1)} for label, count in labels for cap in [40.0, 80.0]
            ],
            "reporting_values": count_reporting,
            "selection_decision_provenance": advance,
            "continuous_certificate": None,
        },
    ]
    claims = [
        {
            "claim_id": "MC-001",
            "case": "radial-station",
            "disposition": "historical",
            "reason": "station measured",
            "selection_decision_provenance": radial,
            "evidence": {"occurrence_count": 1, "station_centre": [18.482, 10.482], "maximal_radius": 0.5156, "length_unit": "mm"},
        },
        {
            "claim_id": "MC-002",
            "case": "radial-station",
            "disposition": "historical",
            "reason": "coarse rungs measured",
            "selection_decision_provenance": radial,
            "evidence": {
                "coarse_step": 0.05,
                "rung_6_radius": 0.2156,
                "rung_6_peak": 61.3,
                "rung_6_cuts_material": True,
                "rung_7_radius": 0.1656,
                "rung_7_peak": 5.9,
                "rung_7_cuts_material": False,
                "angle_unit": "degree",
                "length_unit": "mm",
            },
        },
        {
            "claim_id": "MC-003",
            "case": "radial-station",
            "disposition": "historical",
            "reason": "refined band measured",
            "selection_decision_provenance": radial,
            "evidence": {
                "refined_band_min_radius": 0.1719,
                "refined_band_max_radius": 0.2123,
                "forced_peak": 107.0,
                "rescued_peak": 59.0,
                "angle_unit": "degree",
                "length_unit": "mm",
            },
        },
        {
            "claim_id": "MC-004",
            "case": "radial-subdivisions",
            "disposition": "historical",
            "reason": "audit phase recorded",
            "selection_decision_provenance": radial,
            "evidence": {"audit_position_count": 16, "audit_phase": "entry-angle", "includes_entry_phase": True, "adds_separate_entry_probe": False},
        },
        {
            "claim_id": "MC-005",
            "case": "radial-subdivisions",
            "disposition": "historical",
            "reason": "subdivision rows recorded",
            "selection_decision_provenance": radial,
            "evidence": {"non_entry_circle_count": 244, "rows": subdivision_reporting},
        },
        {
            "claim_id": "MC-006",
            "case": "radial-subdivisions",
            "disposition": "historical",
            "reason": "configuration absent",
            "selection_decision_provenance": radial,
            "evidence": {"history_commit": HISTORY_COMMIT, "missing_configuration": "refinement-without-reporting-ranking"},
        },
        {
            "claim_id": "MC-007",
            "case": "radial-floor",
            "disposition": "historical",
            "reason": "floor rows recorded",
            "selection_decision_provenance": radial,
            "evidence": {"audit_position_count": 16, "rows": floor_reporting},
        },
        {
            "claim_id": "MC-008",
            "case": "radial-margin",
            "disposition": "historical",
            "reason": "margin claims corrected",
            "selection_decision_provenance": radial,
            "evidence": {"rows": margin_reporting, "baseline_available": False, "open_ended_gate_claim_removed": True, "reporting_selection_disclosed": True},
        },
        {
            "claim_id": "MC-009",
            "case": "advance-placement",
            "disposition": "historical",
            "reason": "placement rows recorded",
            "selection_decision_provenance": advance,
            "evidence": {"selected_circle_count": 4, "audit_position_count": 32, "rows": placement_reporting},
        },
        {
            "claim_id": "MC-010",
            "case": "advance-probe-count",
            "disposition": "historical",
            "reason": "unsupported costs removed",
            "selection_decision_provenance": advance,
            "evidence": {"audit_position_count": 60, "rows": count_reporting, "timing_claim_removed": True, "relative_cost_claim_removed": True},
        },
    ]
    return {
        "schema_version": "measurement-claim-payload/v1",
        "batch": "generator",
        "extraction_commit": EXTRACTION_COMMIT,
        "source_commit": source_commit,
        "case_order": CASE_ORDER,
        "cases": cases,
        "claims": claims,
    }


def _input_payload(payload: dict[str, object]) -> dict[str, object]:
    cases = payload["cases"]
    assert type(cases) is list
    return {
        "extraction_commit": payload["extraction_commit"],
        "case_order": payload["case_order"],
        "case_inputs": [
            {
                "case": case["case"],
                "source_claim_ids": case["source_claim_ids"],
                "config": case["config"],
                "selection_decision_provenance": case["selection_decision_provenance"],
            }
            for case in cases
        ],
    }


def _typed_cases(source_commit: str) -> list[Any]:
    cases = copy.deepcopy(_payload(source_commit)["cases"])
    for case in cases:
        config = case["config"]
        calls = config.get("public_calls", [config.get("public_call")])
        for call in calls:
            if call is not None:
                call["polygon"] = tuple(tuple(point) for point in call["polygon"])
    cases[0]["reconstruction"]["target_centre"] = tuple(cases[0]["reconstruction"]["target_centre"])
    return cases


def _source(repository: pathlib.Path, commit: str) -> Any:
    artifact = _artifact()
    lock = _git(repository, "show", f"{commit}:pixi.lock")
    return artifact.SourceSnapshot.build(repository=repository, commit=commit, pixi_lock_sha256=hashlib.sha256(lock).hexdigest())


def _first_nested_key_paths(value: object) -> list[tuple[object, ...]]:
    paths: list[tuple[object, ...]] = []
    seen: set[tuple[tuple[str, ...], str]] = set()

    def visit(item: object, path: tuple[object, ...]) -> None:
        if type(item) is dict:
            mapping = item
            shape = tuple(sorted(mapping))
            for key, child in mapping.items():
                identity = (shape, key)
                if path and identity not in seen:
                    seen.add(identity)
                    paths.append((*path, key))
                visit(child, (*path, key))
        elif type(item) is list:
            for index, child in enumerate(item):
                visit(child, (*path, index))

    visit(value, ())
    return paths


NESTED_KEY_PATHS = _first_nested_key_paths(_payload("a" * 40))


def _delete_path(value: object, path: tuple[object, ...]) -> None:
    parent = value
    for component in path[:-1]:
        parent = parent[component]
    del parent[path[-1]]


def _set_path(value: object, path: tuple[object, ...], replacement: object) -> None:
    parent = value
    for component in path[:-1]:
        parent = parent[component]
    parent[path[-1]] = replacement


def _write_artifact(
    repository: pathlib.Path,
    *,
    hidden: bool = False,
    payload_transform: Optional[Callable[[bytes], bytes]] = None,
    input_transform: Optional[Callable[[dict[str, object]], dict[str, object]]] = None,
    argv: Optional[tuple[str, ...]] = None,
    input_version: str = "generator-measurement-claim-input/v1",
    result_version: str = "generator-measurement-claim-result/v1",
    stamp_transform: Optional[Callable[[bytes], bytes]] = None,
) -> tuple[pathlib.Path, pathlib.PurePosixPath, dict[str, object]]:
    artifact = _artifact()
    commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    payload = _payload(commit)
    payload_bytes = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    if payload_transform is not None:
        payload_bytes = payload_transform(payload_bytes)
    semantic_input = _input_payload(payload)
    if input_transform is not None:
        semantic_input = input_transform(semantic_input)
    started = datetime.datetime(2026, 8, 28, 23, 59, 59, 123456, tzinfo=UTC)
    envelope = artifact.build_envelope(
        artifact_kind=artifact.ArtifactKind("generator-measurement-claims/v1"),
        source=_source(repository, commit),
        started=started,
        finished=started + datetime.timedelta(seconds=2),
        argv=argv or ("/python", "-m", "tools.measurement_claim_probes", "run-generator", "--all"),
        input_version=artifact.IdentityVersion(input_version),
        input_payload=semantic_input,
        result_version=artifact.IdentityVersion(result_version),
        payloads={"generator-claims.json": payload_bytes},
    )
    input_digest = envelope["input_identity"]["sha256"]
    logical_name = f"2026-08-28-{commit[:12]}-generator-{input_digest[:12]}"
    actual_name = f".{logical_name}.stage-fixture" if hidden else logical_name
    result = repository / "benchmarks" / "measurement_claim_results" / actual_name
    result.mkdir(parents=True)
    (result / "generator-claims.json").write_bytes(payload_bytes)
    artifact.write_envelope(result, envelope)
    stamp_path = result / artifact.STAMP_NAME
    if stamp_transform is not None:
        stamp_path.write_bytes(stamp_transform(stamp_path.read_bytes()))
    canonical = pathlib.PurePosixPath("benchmarks", "measurement_claim_results", logical_name)
    return result, canonical, payload


def test_measurement_claim_result_module_exists() -> None:
    assert importlib.util.find_spec("tools.measurement_claim_result") is not None


def test_claim_artifact_identity_constants_are_single_source() -> None:
    module = _claims()
    assert module.ARTIFACT_KIND == "generator-measurement-claims/v1"
    assert module.INPUT_VERSION == "generator-measurement-claim-input/v1"
    assert module.RESULT_VERSION == "generator-measurement-claim-result/v1"
    assert module.PAYLOAD_NAME == "generator-claims.json"
    assert module.EXTRACTION_COMMIT == EXTRACTION_COMMIT


def test_compose_generator_payload_owns_claims_and_exact_semantic_projection() -> None:
    source_commit = "a" * 40
    cases = _typed_cases(source_commit)
    payload = _claims().compose_generator_payload(_artifact().GitObjectId(source_commit), cases)
    assert payload["case_order"] == CASE_ORDER
    assert payload["cases"] == cases
    assert [claim["claim_id"] for claim in payload["claims"]] == [f"MC-{index:03d}" for index in range(1, 11)]
    assert [claim["disposition"] for claim in payload["claims"]] == [
        "re-earned",
        "re-earned",
        "re-earned",
        "corrected",
        "re-earned",
        "historical",
        "re-earned",
        "corrected",
        "corrected",
        "corrected",
    ]
    assert "omitted forward-peak histogram" in payload["claims"][8]["reason"]
    assert _claims().generator_semantic_input(payload) == _input_payload(payload)
    decoded = _claims()._decode(json.dumps(payload, allow_nan=False).encode("utf-8"), "payload")
    assert _claims().validate_generator_payload(decoded) == decoded


def test_compose_generator_payload_fails_conditional_claims_closed() -> None:
    source_commit = "a" * 40
    cases = _typed_cases(source_commit)
    cases[0]["reporting_values"]["rung_6_peak"] = 61.4
    cases[1]["reporting_values"][0]["worst_peak"] = 88.7
    cases[2]["reporting_values"][0]["circles_over_cap"] = 9
    cases[2]["native_sampled_decisions"][0]["observations"] = _counts(3904, 9)
    payload = _claims().compose_generator_payload(_artifact().GitObjectId(source_commit), cases)
    assert [payload["claims"][index]["disposition"] for index in (0, 1, 2)] == ["historical"] * 3
    assert payload["claims"][4]["disposition"] == "corrected"
    assert payload["claims"][6]["disposition"] == "corrected"


def test_validate_generator_payload_accepts_complete_concrete_payload(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    payload = _payload(commit)
    assert _claims().validate_generator_payload(payload) == payload


@pytest.mark.parametrize(
    ("mutate", "fragment"),
    [
        (lambda payload: {**payload, "batch": "benchmark"}, "batch"),
        (lambda payload: {**payload, "extra": True}, "keys"),
        (lambda payload: {**payload, "case_order": list(reversed(CASE_ORDER))}, "case_order"),
        (lambda payload: {**payload, "claims": payload["claims"][:-1]}, "claims"),
    ],
)
def test_validate_generator_payload_rejects_root_contract_damage(
    tmp_path: pathlib.Path,
    mutate: Callable[[dict[str, object]], dict[str, object]],
    fragment: str,
) -> None:
    repository = _repository(tmp_path)
    commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match=fragment):
        _claims().validate_generator_payload(mutate(_payload(commit)))


def test_validate_generator_payload_rejects_nested_discriminant_and_nonfinite(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    payload = _payload(commit)
    cases = payload["cases"]
    assert type(cases) is list
    cases[0]["config"]["public_call"]["tool_diameter"] = float("inf")
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="finite"):
        _claims().validate_generator_payload(payload)


@pytest.mark.parametrize("path", NESTED_KEY_PATHS, ids=lambda path: ".".join(map(str, path)))
def test_validate_generator_payload_rejects_every_nested_missing_key(path: tuple[object, ...]) -> None:
    payload = copy.deepcopy(_payload("a" * 40))
    _delete_path(payload, path)
    with pytest.raises(_claims().MeasurementClaimError):
        _claims().validate_generator_payload(payload)


@pytest.mark.parametrize(
    ("path", "replacement"),
    [
        (("cases", 0, "config", "subdivisions"), 9),
        (("cases", 1, "config", "audit", "phase"), "advance-direction"),
        (("cases", 2, "config", "floor_values"), [1.0, 0.5, 0.25]),
        (("cases", 3, "config", "public_call", "tea_cap"), 60.0),
        (("cases", 4, "config", "audit", "excludes_entry_probe"), False),
        (("cases", 5, "config", "generator_configurations", 0, "probe_count"), 8),
        (("cases", 0, "reconstruction", "stock_model"), "generator-pre-bridge/v1"),
        (("cases", 1, "reconstruction", "audit_position_count"), 17),
        (("cases", 1, "reconstruction", "non_entry_circle_counts"), [244, 244, 244, 244, 245]),
        (("cases", 2, "reconstruction", "non_entry_circle_counts"), [240, 241]),
        (("cases", 4, "reconstruction", "selected_circle_policy"), "forced/v1"),
        (("cases", 0, "native_sampled_decisions", "refined_candidates", 0, "radius"), 0.2),
        (("cases", 1, "native_sampled_decisions", 0, "subdivisions"), 2),
        (("cases", 2, "native_sampled_decisions", 0, "floor_steps"), 0.5),
        (("cases", 2, "native_sampled_decisions", 1, "selected_circles"), 240),
        (("cases", 2, "native_sampled_decisions", 1, "observations"), _counts(3840, 8)),
        (("cases", 3, "native_sampled_decisions", 0, "refinement_margin"), 1.4),
        (("cases", 4, "native_sampled_decisions", 0, "cap"), 100.0),
        (("cases", 5, "native_sampled_decisions", 0, "label"), "8"),
        (("cases", 0, "reporting_values", "angle_unit"), "radian"),
        (("cases", 1, "reporting_values", 0, "circles_over_cap"), 0),
        (("cases", 1, "reporting_values", 0, "circles_over_cap"), 13),
        (("cases", 1, "reporting_values", 0, "circles_over_cap"), 245),
        (("cases", 3, "reporting_values", 0, "length_unit"), "metre"),
        (("cases", 4, "reporting_values", 0, "offset_bin_counts", 0, "angle_unit"), "radian"),
        (("cases", 1, "selection_decision_provenance", "native_cap_decision_site"), "stock[1]"),
        (("cases", 5, "selection_decision_provenance", "reporting_driven_decision_sites"), [{}]),
        (("claims", 0, "reason"), "unsafe|reason"),
        (("claims", 3, "disposition"), "re-earned"),
        (("claims", 5, "disposition"), "corrected"),
        (("claims", 7, "disposition"), "re-earned"),
        (("claims", 5, "evidence", "missing_configuration"), "present"),
        (("claims", 9, "disposition"), "re-earned"),
    ],
)
def test_validate_generator_payload_rejects_nested_value_contract_damage(path: tuple[object, ...], replacement: object) -> None:
    payload = copy.deepcopy(_payload("a" * 40))
    _set_path(payload, path, replacement)
    with pytest.raises(_claims().MeasurementClaimError):
        _claims().validate_generator_payload(payload)


def test_validate_generator_payload_rejects_noncanonical_placement_bin_centres() -> None:
    payload = copy.deepcopy(_payload("a" * 40))
    for row in payload["cases"][4]["reporting_values"]:
        row["offset_bin_counts"][0]["offset"] = -140.0
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="offset"):
        _claims().validate_generator_payload(payload)


def test_validate_generator_payload_rejects_reordered_decision_sites() -> None:
    payload = copy.deepcopy(_payload("a" * 40))
    sites = payload["cases"][0]["selection_decision_provenance"]["reporting_driven_decision_sites"]
    sites.reverse()
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="must equal"):
        _claims().validate_generator_payload(payload)


def test_validate_generator_payload_distinguishes_radial_position_and_circle_counts() -> None:
    payload = copy.deepcopy(_payload("a" * 40))
    case = payload["cases"][1]
    case["reconstruction"]["non_entry_circle_counts"] = [2, 2, 2, 2, 2]
    payload["claims"][4]["evidence"]["non_entry_circle_count"] = 2
    for native, report in zip(case["native_sampled_decisions"], case["reporting_values"]):
        native["selected_circles"], native["observations"] = 2, _counts(32, 16)
        report["circles_over_cap"] = 1
    assert _claims().validate_generator_payload(payload) == payload


@pytest.mark.parametrize("hidden", [False, True])
def test_validate_claim_artifact_accepts_final_and_owned_stage(tmp_path: pathlib.Path, hidden: bool) -> None:
    repository = _repository(tmp_path)
    result, canonical, payload = _write_artifact(repository, hidden=hidden)
    validated, envelope, started, artifact_directory = _claims().validate_claim_artifact(result)
    assert validated == payload
    assert str(envelope.commit) == payload["source_commit"]
    assert started == datetime.datetime(2026, 8, 28, 23, 59, 59, 123456, tzinfo=UTC)
    assert artifact_directory == canonical


@pytest.mark.parametrize(
    "transform",
    [
        lambda raw: raw.replace(b"{", b'{"schema_version":"measurement-claim-payload/v1",', 1),
        lambda raw: b"\xff" + raw,
        lambda raw: raw.replace(b'"tool_diameter":2.0', b'"tool_diameter":NaN', 1),
        lambda raw: raw.replace(b'"tool_diameter":2.0', b'"tool_diameter":1e999', 1),
    ],
)
def test_validate_claim_artifact_rejects_strict_json_damage(tmp_path: pathlib.Path, transform: Callable[[bytes], bytes]) -> None:
    repository = _repository(tmp_path)
    result, _, _ = _write_artifact(repository, payload_transform=transform)
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError):
        _claims().validate_claim_artifact(result)


@pytest.mark.parametrize(
    "transform",
    [
        lambda raw: raw.replace(b"{", b'{"envelope_version":"measurement-artifact-envelope/v1",', 1),
        lambda raw: b"\xff" + raw,
        lambda raw: raw.replace(b'"dirty": false', b'"dirty": NaN', 1),
        lambda raw: raw.replace(b'"dirty": false', b'"dirty": 1e999', 1),
    ],
)
def test_validate_claim_artifact_rejects_strict_stamp_json_damage(tmp_path: pathlib.Path, transform: Callable[[bytes], bytes]) -> None:
    repository = _repository(tmp_path)
    result, _, _ = _write_artifact(repository, stamp_transform=transform)
    with pytest.raises(_artifact().MeasurementArtifactError):
        _claims().validate_claim_artifact(result)


def test_validate_claim_artifact_calls_common_before_family_parsing(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    repository = _repository(tmp_path)
    result, _, _ = _write_artifact(repository)
    events: list[str] = []
    common = _artifact().validate_envelope
    decode = _claims()._decode

    def recording_common(*args: object, **kwargs: object) -> Any:
        events.append("common")
        return common(*args, **kwargs)

    def recording_decode(data: bytes, field: str) -> object:
        events.append(f"decode:{field}")
        return decode(data, field)

    monkeypatch.setattr(_artifact(), "validate_envelope", recording_common)
    monkeypatch.setattr(_claims(), "_decode", recording_decode)
    _claims().validate_claim_artifact(result)
    assert events == ["common", "decode:stamp", "decode:generator-claims.json"]


def test_validate_claim_artifact_rejects_changed_stamp_bytes(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    repository = _repository(tmp_path)
    result, _, _ = _write_artifact(repository)
    original = _artifact().validate_envelope

    def mutate_after_validation(*args: object, **kwargs: object) -> Any:
        validated = original(*args, **kwargs)
        stamp = result / _artifact().STAMP_NAME
        stamp.write_bytes(stamp.read_bytes() + b" ")
        return validated

    monkeypatch.setattr(_artifact(), "validate_envelope", mutate_after_validation)
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="changed"):
        _claims().validate_claim_artifact(result)


def test_validate_claim_artifact_rejects_wrong_root_and_symlink(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    result, _, _ = _write_artifact(repository)
    wrong_root = repository / "benchmarks" / "other_results"
    wrong_root.mkdir()
    moved = wrong_root / result.name
    result.rename(moved)
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="directly below"):
        _claims().validate_claim_artifact(moved)
    moved.rename(result)
    link = result.with_name(f"{result.name}-link")
    link.symlink_to(result, target_is_directory=True)
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="symlink"):
        _claims().validate_claim_artifact(link)


def test_validate_claim_artifact_calls_common_validation_once_and_rejects_changed_bytes(
    tmp_path: pathlib.Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repository = _repository(tmp_path)
    result, _, _ = _write_artifact(repository)
    original = _artifact().validate_envelope
    calls = 0

    def mutate_after_validation(*args: object, **kwargs: object) -> Any:
        nonlocal calls
        calls += 1
        validated = original(*args, **kwargs)
        (result / "generator-claims.json").write_bytes(b"{}")
        return validated

    monkeypatch.setattr(_artifact(), "validate_envelope", mutate_after_validation)
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError, match="changed"):
        _claims().validate_claim_artifact(result)
    assert calls == 1


@pytest.mark.parametrize(
    "damage",
    ["wrong-name", "empty-stage", "wrong-argv", "wrong-input-version", "wrong-result-version", "wrong-input"],
)
def test_validate_claim_artifact_rejects_authenticated_family_contract_damage(tmp_path: pathlib.Path, damage: str) -> None:
    repository = _repository(tmp_path)
    kwargs: dict[str, object] = {}
    if damage == "wrong-argv":
        kwargs["argv"] = ("/python", "-m", "tools.measurement_claim_probes", "run-generator", "--case", "radial-station")
    elif damage == "wrong-input-version":
        kwargs["input_version"] = "wrong/v1"
    elif damage == "wrong-result-version":
        kwargs["result_version"] = "wrong/v1"
    elif damage == "wrong-input":
        kwargs["input_transform"] = lambda value: {**value, "case_order": list(reversed(CASE_ORDER))}
    result, _, _ = _write_artifact(repository, hidden=damage == "empty-stage", **kwargs)
    if damage == "wrong-name":
        renamed = result.with_name("wrong-name")
        result.rename(renamed)
        result = renamed
    elif damage == "empty-stage":
        logical_name = result.name.removeprefix(".").split(".stage-", 1)[0]
        renamed = result.with_name(f".{logical_name}.stage-")
        result.rename(renamed)
        result = renamed
    with pytest.raises(_claims().InvalidMeasurementClaimPayloadError):
        _claims().validate_claim_artifact(result)


def test_concrete_schema_passes_positive_and_negative_strict_mypy_fixtures(tmp_path: pathlib.Path) -> None:
    positive = tmp_path / "positive_claim_schema.py"
    positive.write_text(
        """from tools.measurement_artifact import GitObjectId
import tools.measurement_claim_result as m
source = GitObjectId("a" * 40)
point: m.WorldPointMillimetres = (m.WorldMillimetres(0.0), m.WorldMillimetres(0.0), m.WorldMillimetres(0.0))
rectangle: m.WorldRectangleMillimetres = (point, point, point, point)
call: m.GeneratorPublicCallPayload = {"polygon": rectangle, "holes": [], "tool_diameter": m.Millimetres(2.0), "tea_cap": m.Degrees(60.0), "guide_step": m.ToolDiameters(0.025), "max_advance": m.ToolDiameters(1.0), "radial_clearance": m.Millimetres(0.002), "climb": True, "cut_z": m.WorldMillimetres(0.0), "clearance_z": m.WorldMillimetres(2.0), "max_passes": 1, "samples_per_radian": m.SamplesPerRadian(10.0)}
radial_p: m.RadialSelectionDecisionProvenancePayload = {"policy": "radial-known-reporting-driven-selection/v1", "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]", "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]", "reporting_driven_decision_sites": []}
advance_p: m.AdvanceSelectionDecisionProvenancePayload = {"policy": "advance-native-cap-reporting-observation/v1", "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]", "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]", "reporting_driven_decision_sites": []}
radial_audit: m.RadialAuditConfigPayload = {"phase": "entry-angle", "probe_offsets": [m.Degrees(0.0)], "includes_entry_phase": True, "adds_separate_entry_probe": False, "excludes_chain_entry_circles": True}
placement_audit: m.AdvancePlacementAuditConfigPayload = {"phase": "advance-direction", "probe_offsets": [m.Degrees(0.0)], "excludes_entry_probe": True}
count_audit: m.AdvanceProbeCountAuditConfigPayload = {"phase": "advance-half-step", "probe_offsets": [m.Degrees(0.0)], "excludes_entry_probe": True}
station_cfg: m.RadialStationConfigPayload = {"public_call": call, "subdivisions": 8, "floor_steps": m.CoarseSteps(0.5), "refinement_margin": m.DimensionlessRatio(1.4), "generator_probe_angles": [m.Degrees(0.0)], "ladder_step": m.ToolDiameters(0.025), "ladder_span": m.ToolDiameters(1.0), "ladder_rungs": 40, "max_radial_sweeps": 40}
sub_cfg: m.RadialSubdivisionsConfigPayload = {"public_call": call, "subdivision_values": [8], "floor_steps": m.CoarseSteps(0.5), "refinement_margin": m.DimensionlessRatio(1.4), "generator_probe_angles": [m.Degrees(0.0)], "audit": radial_audit}
floor_cfg: m.RadialFloorConfigPayload = {"public_call": call, "floor_values": [m.CoarseSteps(0.5)], "subdivisions": 8, "refinement_margin": m.DimensionlessRatio(1.4), "generator_probe_angles": [m.Degrees(0.0)], "audit": radial_audit}
margin_cfg: m.RadialMarginConfigPayload = {"public_call": call, "margin_values": [m.DimensionlessRatio(1.4)], "subdivisions": 8, "floor_steps": m.CoarseSteps(0.5), "generator_probe_angles": [m.Degrees(0.0)], "audit": radial_audit}
placement_cfg: m.AdvancePlacementConfigPayload = {"public_calls": [call], "generator_probe_angles": [m.Degrees(0.0)], "generator_prepends_entry_probe": True, "audit": placement_audit}
probe_cfg: m.AdvanceProbeCountConfigPayload = {"public_calls": [call], "generator_configurations": [{"label": "8", "probe_count": 8, "probe_angles": [m.Degrees(0.0)]}], "generator_prepends_entry_probe": True, "audit": count_audit}
station_recon: m.RadialStationReconstructionPayload = {"stock_model": "generator-faithful-radial-pre-bridge/v1", "target_centre": (m.WorldMillimetres(0.0), m.WorldMillimetres(0.0)), "centre_decimal_places": 3, "target_maximal_radius": m.Millimetres(0.5), "radius_decimal_places": 4, "occurrence_count": 1, "coarse_step": m.Millimetres(0.05), "coarse_rungs": [6], "refined_radius_sequence": [m.Millimetres(0.2)]}
radial_recon: m.RadialSweepReconstructionPayload = {"stock_model": "generator-faithful-radial-pre-bridge/v1", "audit_position_count": 1, "audit_includes_entry_phase": True, "audit_adds_separate_entry_probe": False, "excludes_chain_entry_circles": True, "non_entry_circle_counts": [1]}
advance_recon: m.AdvanceReconstructionPayload = {"stock_model": "generator-pre-bridge/v1", "selected_circle_policy": "accepted-non-forced/v1", "original_calls_per_wrapper": 1, "generator_prepends_entry_probe": True, "audit_excludes_entry_probe": True, "selected_circle_count": 1}
counts: m.NativeObservationCountsPayload = {"observations": 1, "accepted": 1, "exceeded": 0}
station_native: m.RadialStationNativeDecisionsPayload = {"rung_6_cap_exceeded": False, "rung_6_cuts_material": True, "rung_7_cap_exceeded": None, "rung_7_cuts_material": None, "refined_candidates": [{"radius": m.Millimetres(0.2), "cap_exceeded": False, "cuts_material": True}]}
station_report: m.RadialStationReportingPayload = {"maximal_radius": m.Millimetres(0.5), "rung_6_radius": m.Millimetres(0.2), "rung_6_peak": m.Degrees(50.0), "rung_7_radius": None, "rung_7_peak": None, "refined_band_min_radius": m.Millimetres(0.1), "refined_band_max_radius": m.Millimetres(0.2), "forced_peak": m.Degrees(70.0), "rescued_peak": m.Degrees(50.0), "angle_unit": "degree", "length_unit": "mm"}
sub_native: list[m.RadialSubdivisionNativeRowPayload] = [{"subdivisions": 8, "selected_circles": 1, "observations": counts}]
sub_report: list[m.RadialSubdivisionReportingRowPayload] = [{"subdivisions": 8, "worst_peak": m.Degrees(50.0), "circles_over_cap": 0, "angle_unit": "degree"}]
floor_native: list[m.RadialFloorNativeRowPayload] = [{"floor_steps": m.CoarseSteps(0.5), "selected_circles": 1, "observations": counts}]
floor_report: list[m.RadialFloorReportingRowPayload] = [{"floor_steps": m.CoarseSteps(0.5), "worst_peak": m.Degrees(50.0), "circles_over_cap": 0, "angle_unit": "degree"}]
margin_native: list[m.RadialMarginNativeRowPayload] = [{"refinement_margin": m.DimensionlessRatio(1.4), "selected_circles": 1, "observations": counts}]
margin_report: list[m.RadialMarginReportingRowPayload] = [{"refinement_margin": m.DimensionlessRatio(1.4), "worst_peak": m.Degrees(50.0), "circles_over_cap": 0, "cutting_length": m.Millimetres(1.0), "angle_unit": "degree", "length_unit": "mm"}]
placement_native: list[m.AdvancePlacementNativeRowPayload] = [{"cap": m.Degrees(80.0), "climb": True, "selected_circles": 1, "observations": counts}]
placement_report: list[m.AdvancePlacementReportingRowPayload] = [{"cap": m.Degrees(80.0), "climb": True, "selected_circles": 1, "positions_over_cap": 0, "worst_peak": m.Degrees(50.0), "worst_peak_offset": m.Degrees(0.0), "old_probe_peak": m.Degrees(40.0), "offset_bin_counts": [{"offset": m.Degrees(0.0), "positions_over_cap": 0, "angle_unit": "degree"}], "angle_unit": "degree"}]
probe_native: list[m.AdvanceProbeCountNativeRowPayload] = [{"label": "8", "probe_count": 8, "cap": m.Degrees(80.0), "selected_circles": 1, "observations": counts}]
probe_report: list[m.AdvanceProbeCountReportingRowPayload] = [{"label": "8", "probe_count": 8, "cap": m.Degrees(80.0), "worst_peak": m.Degrees(50.0), "angle_unit": "degree"}]
station: m.RadialStationCasePayload = {"case": "radial-station", "source_claim_ids": ["MC-001", "MC-002", "MC-003"], "config": station_cfg, "reconstruction": station_recon, "native_sampled_decisions": station_native, "reporting_values": station_report, "selection_decision_provenance": radial_p, "continuous_certificate": None}
sub: m.RadialSubdivisionsCasePayload = {"case": "radial-subdivisions", "source_claim_ids": ["MC-004", "MC-005", "MC-006"], "config": sub_cfg, "reconstruction": radial_recon, "native_sampled_decisions": sub_native, "reporting_values": sub_report, "selection_decision_provenance": radial_p, "continuous_certificate": None}
floor: m.RadialFloorCasePayload = {"case": "radial-floor", "source_claim_ids": ["MC-007"], "config": floor_cfg, "reconstruction": radial_recon, "native_sampled_decisions": floor_native, "reporting_values": floor_report, "selection_decision_provenance": radial_p, "continuous_certificate": None}
margin: m.RadialMarginCasePayload = {"case": "radial-margin", "source_claim_ids": ["MC-008"], "config": margin_cfg, "reconstruction": radial_recon, "native_sampled_decisions": margin_native, "reporting_values": margin_report, "selection_decision_provenance": radial_p, "continuous_certificate": None}
placement: m.AdvancePlacementCasePayload = {"case": "advance-placement", "source_claim_ids": ["MC-009"], "config": placement_cfg, "reconstruction": advance_recon, "native_sampled_decisions": placement_native, "reporting_values": placement_report, "selection_decision_provenance": advance_p, "continuous_certificate": None}
probe: m.AdvanceProbeCountCasePayload = {"case": "advance-probe-count", "source_claim_ids": ["MC-010"], "config": probe_cfg, "reconstruction": advance_recon, "native_sampled_decisions": probe_native, "reporting_values": probe_report, "selection_decision_provenance": advance_p, "continuous_certificate": None}
station_input: m.RadialStationCaseInputPayload = {"case": "radial-station", "source_claim_ids": ["MC-001", "MC-002", "MC-003"], "config": station_cfg, "selection_decision_provenance": radial_p}
sub_input: m.RadialSubdivisionsCaseInputPayload = {"case": "radial-subdivisions", "source_claim_ids": ["MC-004", "MC-005", "MC-006"], "config": sub_cfg, "selection_decision_provenance": radial_p}
floor_input: m.RadialFloorCaseInputPayload = {"case": "radial-floor", "source_claim_ids": ["MC-007"], "config": floor_cfg, "selection_decision_provenance": radial_p}
margin_input: m.RadialMarginCaseInputPayload = {"case": "radial-margin", "source_claim_ids": ["MC-008"], "config": margin_cfg, "selection_decision_provenance": radial_p}
placement_input: m.AdvancePlacementCaseInputPayload = {"case": "advance-placement", "source_claim_ids": ["MC-009"], "config": placement_cfg, "selection_decision_provenance": advance_p}
probe_input: m.AdvanceProbeCountCaseInputPayload = {"case": "advance-probe-count", "source_claim_ids": ["MC-010"], "config": probe_cfg, "selection_decision_provenance": advance_p}
e1: m.MC001EvidencePayload = {"occurrence_count": 1, "station_centre": station_recon["target_centre"], "maximal_radius": station_report["maximal_radius"], "length_unit": "mm"}
e2: m.MC002EvidencePayload = {"coarse_step": m.Millimetres(0.05), "rung_6_radius": station_report["rung_6_radius"], "rung_6_peak": station_report["rung_6_peak"], "rung_6_cuts_material": True, "rung_7_radius": None, "rung_7_peak": None, "rung_7_cuts_material": None, "angle_unit": "degree", "length_unit": "mm"}
e3: m.MC003EvidencePayload = {"refined_band_min_radius": station_report["refined_band_min_radius"], "refined_band_max_radius": station_report["refined_band_max_radius"], "forced_peak": station_report["forced_peak"], "rescued_peak": station_report["rescued_peak"], "angle_unit": "degree", "length_unit": "mm"}
e4: m.MC004EvidencePayload = {"audit_position_count": 1, "audit_phase": "entry-angle", "includes_entry_phase": True, "adds_separate_entry_probe": False}
e5: m.MC005EvidencePayload = {"non_entry_circle_count": 1, "rows": sub_report}
e6: m.MC006EvidencePayload = {"history_commit": source, "missing_configuration": "refinement-without-reporting-ranking"}
e7: m.MC007EvidencePayload = {"audit_position_count": 1, "rows": floor_report}
e8: m.MC008EvidencePayload = {"rows": margin_report, "baseline_available": False, "open_ended_gate_claim_removed": True, "reporting_selection_disclosed": True}
e9: m.MC009EvidencePayload = {"selected_circle_count": 1, "audit_position_count": 1, "rows": placement_report}
e10: m.MC010EvidencePayload = {"audit_position_count": 1, "rows": probe_report, "timing_claim_removed": True, "relative_cost_claim_removed": True}
c1: m.MC001ClaimPayload = {"claim_id": "MC-001", "case": "radial-station", "disposition": "re-earned", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e1}
c2: m.MC002ClaimPayload = {"claim_id": "MC-002", "case": "radial-station", "disposition": "re-earned", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e2}
c3: m.MC003ClaimPayload = {"claim_id": "MC-003", "case": "radial-station", "disposition": "re-earned", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e3}
c4: m.MC004ClaimPayload = {"claim_id": "MC-004", "case": "radial-subdivisions", "disposition": "corrected", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e4}
c5: m.MC005ClaimPayload = {"claim_id": "MC-005", "case": "radial-subdivisions", "disposition": "corrected", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e5}
c6: m.MC006ClaimPayload = {"claim_id": "MC-006", "case": "radial-subdivisions", "disposition": "historical", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e6}
c7: m.MC007ClaimPayload = {"claim_id": "MC-007", "case": "radial-floor", "disposition": "re-earned", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e7}
c8: m.MC008ClaimPayload = {"claim_id": "MC-008", "case": "radial-margin", "disposition": "corrected", "reason": "r", "selection_decision_provenance": radial_p, "evidence": e8}
c9: m.MC009ClaimPayload = {"claim_id": "MC-009", "case": "advance-placement", "disposition": "corrected", "reason": "r", "selection_decision_provenance": advance_p, "evidence": e9}
c10: m.MC010ClaimPayload = {"claim_id": "MC-010", "case": "advance-probe-count", "disposition": "corrected", "reason": "r", "selection_decision_provenance": advance_p, "evidence": e10}
case_order: list[m.GeneratorCase] = ["radial-station", "radial-subdivisions", "radial-floor", "radial-margin", "advance-placement", "advance-probe-count"]
cases: list[m.GeneratorCasePayload] = [station, sub, floor, margin, placement, probe]
inputs: list[m.GeneratorCaseInputPayload] = [station_input, sub_input, floor_input, margin_input, placement_input, probe_input]
claims: list[m.GeneratorClaimRecord] = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10]
generator_input: m.GeneratorClaimInputPayload = {"extraction_commit": source, "case_order": case_order, "case_inputs": inputs}
generator_result: m.GeneratorClaimPayload = {"schema_version": "measurement-claim-payload/v1", "batch": "generator", "extraction_commit": source, "source_commit": source, "case_order": case_order, "cases": cases, "claims": claims}
""",
        encoding="utf-8",
    )
    negative = tmp_path / "negative_claim_schema.py"
    negative.write_text(
        """from typing import Dict
from tools.measurement_claim_result import AdvancePlacementCasePayload, AdvancePlacementConfigPayload
from tools.measurement_claim_result import GeneratorClaimPayload, Millimetres, RadialStationCasePayload, RadialStationConfigPayload
from tools.measurement_claim_result import WorldMillimetres, WorldPointMillimetres
wrong_frame: Millimetres = WorldMillimetres(0.0)
wrong_arity: WorldPointMillimetres = (WorldMillimetres(0.0), WorldMillimetres(0.0))
station_config: RadialStationConfigPayload; wrong_config: AdvancePlacementConfigPayload = station_config
station_case: RadialStationCasePayload; wrong_case: AdvancePlacementCasePayload = station_case
erased: Dict[str, object] = {}; escaped: GeneratorClaimPayload = erased
""",
        encoding="utf-8",
    )
    repository, command = pathlib.Path(__file__).resolve().parents[2], [sys.executable, "-m", "mypy", "--strict", "--warn-unused-ignores"]
    positive_result = subprocess.run([*command, str(positive)], cwd=repository, check=False, capture_output=True, text=True)
    negative_result = subprocess.run([*command, str(negative)], cwd=repository, check=False, capture_output=True, text=True)
    assert positive_result.returncode == 0 and negative_result.returncode != 0, positive_result.stdout + positive_result.stderr
    assert all(fragment in negative_result.stdout for fragment in "Incompatible types;tuple[World;RadialStationConfig;RadialStationCase;dict[str".split(";"))
