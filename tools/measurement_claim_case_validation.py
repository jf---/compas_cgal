"""Canonical configuration and structural validation for generator cases."""

from __future__ import annotations

import copy
import math
from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import cast

from tools import measurement_claim_json
from tools.measurement_claim_errors import InvalidMeasurementClaimConfigError
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_errors import UnknownMeasurementClaimCaseError
from tools.measurement_claim_schema import GeneratorCase
from tools.measurement_claim_schema import GeneratorCasePayload

CASE_ORDER: List[GeneratorCase] = [
    "radial-station",
    "radial-subdivisions",
    "radial-floor",
    "radial-margin",
    "advance-placement",
    "advance-probe-count",
]
CASE_CLAIMS = {
    "radial-station": ["MC-001", "MC-002", "MC-003"],
    "radial-subdivisions": ["MC-004", "MC-005", "MC-006"],
    "radial-floor": ["MC-007"],
    "radial-margin": ["MC-008"],
    "advance-placement": ["MC-009"],
    "advance-probe-count": ["MC-010"],
}
_PLACEMENT_BIN_CENTRES = [angle if angle <= 180.0 else angle - 360.0 for angle in (360.0 * index / 32 for index in range(32))]
_COMMON_KEYS = (
    "case",
    "source_claim_ids",
    "config",
    "reconstruction",
    "native_sampled_decisions",
    "reporting_values",
    "selection_decision_provenance",
    "continuous_certificate",
)


def fail_payload(field: str, detail: str) -> None:
    """Raise the named payload error for one invalid field.

    Args:
        field: Field path that owns the violation.
        detail: Human-readable contract violation.

    Raises:
        InvalidMeasurementClaimPayloadError: Always.
    """
    raise InvalidMeasurementClaimPayloadError(f"{field}: {detail}")


def validate_object(value: object, keys: Sequence[str], field: str) -> Dict[str, object]:
    """Validate an exact JSON object with one exact key set.

    Args:
        value: Candidate decoded JSON value.
        keys: Required object keys.
        field: Field path used in diagnostics.

    Returns:
        The validated mutable object view.

    Raises:
        InvalidMeasurementClaimPayloadError: The value or keys differ.
    """
    if type(value) is not dict:
        fail_payload(field, "must be an exact JSON object")
    result = cast(Dict[str, object], value)
    if set(result) != set(keys):
        fail_payload(field, f"keys are {tuple(result)}, expected {tuple(keys)}")
    return result


def validate_array(value: object, field: str) -> List[object]:
    """Validate one exact decoded JSON array.

    Args:
        value: Candidate decoded JSON value.
        field: Field path used in diagnostics.

    Returns:
        The validated mutable list view.

    Raises:
        InvalidMeasurementClaimPayloadError: The value is not an exact list.
    """
    if type(value) is not list:
        fail_payload(field, "must be an exact JSON array")
    return cast(List[object], value)


def validate_integer(value: object, field: str) -> int:
    """Validate one non-negative exact integer.

    Args:
        value: Candidate decoded JSON value.
        field: Field path used in diagnostics.

    Returns:
        The validated integer.

    Raises:
        InvalidMeasurementClaimPayloadError: The value is negative or not an int.
    """
    if type(value) is not int or value < 0:
        fail_payload(field, "must be a non-negative exact integer")
    return cast(int, value)


def validate_float(value: object, field: str, *, optional: bool = False) -> Optional[float]:
    """Validate one finite exact float, optionally allowing null.

    Args:
        value: Candidate decoded JSON value.
        field: Field path used in diagnostics.
        optional: Whether ``None`` is accepted.

    Returns:
        The validated float, or ``None`` when allowed.

    Raises:
        InvalidMeasurementClaimPayloadError: The value is not finite float data.
    """
    if optional and value is None:
        return None
    if type(value) is not float or not math.isfinite(value):
        fail_payload(field, "must be a finite exact float")
    return cast(float, value)


def _boolean(value: object, field: str, *, optional: bool = False) -> Optional[bool]:
    if optional and value is None:
        return None
    if type(value) is not bool:
        fail_payload(field, "must be an exact boolean")
    return cast(bool, value)


def validate_literal(value: object, expected: object, field: str) -> None:
    """Require exact runtime type and value equality.

    Args:
        value: Candidate decoded JSON value.
        expected: Canonical literal.
        field: Field path used in diagnostics.

    Raises:
        InvalidMeasurementClaimPayloadError: Type or value differs.
    """
    if type(value) is not type(expected) or value != expected:
        fail_payload(field, f"must equal {expected!r}")


def validate_same(value: object, expected: object, field: str) -> None:
    """Recursively require exact JSON-tree equality.

    Args:
        value: Candidate decoded JSON tree.
        expected: Canonical JSON tree.
        field: Field path used in diagnostics.

    Raises:
        InvalidMeasurementClaimPayloadError: Shape, type, or value differs.
    """
    if type(value) is not type(expected):
        fail_payload(field, f"type differs from canonical {type(expected).__name__}")
    if type(value) is dict:
        actual_object = cast(Dict[str, object], value)
        expected_object = cast(Dict[str, object], expected)
        if set(actual_object) != set(expected_object):
            fail_payload(field, "keys differ from canonical configuration")
        for key in expected_object:
            validate_same(actual_object[key], expected_object[key], f"{field}.{key}")
    elif type(value) is list:
        actual_list = cast(List[object], value)
        expected_list = cast(List[object], expected)
        if len(actual_list) != len(expected_list):
            fail_payload(field, "length differs from canonical configuration")
        for index, item in enumerate(expected_list):
            validate_same(actual_list[index], item, f"{field}[{index}]")
    elif value != expected:
        fail_payload(field, f"must equal {expected!r}")


def validate_finite_tree(value: object, field: str) -> None:
    """Reject non-finite floats anywhere in a decoded tree.

    Args:
        value: Candidate decoded JSON tree.
        field: Field path used in diagnostics.

    Raises:
        InvalidMeasurementClaimPayloadError: A non-finite float is present.
    """
    measurement_claim_json.validate_finite(value, field, InvalidMeasurementClaimPayloadError)


def _public_call(width: float, height: float, cap: float, climb: bool = True) -> Dict[str, object]:
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


def _radial_audit() -> Dict[str, object]:
    return {
        "phase": "entry-angle",
        "probe_offsets": [360.0 * index / 16 for index in range(16)],
        "includes_entry_phase": True,
        "adds_separate_entry_probe": False,
        "excludes_chain_entry_circles": True,
    }


def _expected_config(case: str) -> Dict[str, object]:
    generator_angles = [360.0 * index / 32 for index in range(32)]
    if case == "radial-station":
        return {
            "public_call": _public_call(20.0, 12.0, 60.0),
            "subdivisions": 8,
            "floor_steps": 0.5,
            "refinement_margin": 1.4,
            "generator_probe_angles": generator_angles,
            "ladder_step": 0.025,
            "ladder_span": 1.0,
            "ladder_rungs": 40,
            "max_radial_sweeps": 40,
        }
    if case == "radial-subdivisions":
        return {
            "public_call": _public_call(20.0, 12.0, 60.0),
            "subdivision_values": [1, 2, 4, 8, 16],
            "floor_steps": 0.5,
            "refinement_margin": 1.4,
            "generator_probe_angles": generator_angles,
            "audit": _radial_audit(),
        }
    if case == "radial-floor":
        return {
            "public_call": _public_call(20.0, 12.0, 60.0),
            "floor_values": [0.25, 0.5, 1.0],
            "subdivisions": 8,
            "refinement_margin": 1.4,
            "generator_probe_angles": generator_angles,
            "audit": _radial_audit(),
        }
    if case == "radial-margin":
        return {
            "public_call": _public_call(6.0, 4.0, 40.0),
            "margin_values": [1.25, 1.4, 1.5, 1.75, 2.0],
            "subdivisions": 8,
            "floor_steps": 0.5,
            "generator_probe_angles": generator_angles,
            "audit": _radial_audit(),
        }
    if case == "advance-placement":
        return {
            "public_calls": [_public_call(20.0, 12.0, cap, climb) for cap in (80.0, 100.0) for climb in (True, False)],
            "generator_probe_angles": [-60.0, 0.0, 60.0],
            "generator_prepends_entry_probe": True,
            "audit": {"phase": "advance-direction", "probe_offsets": generator_angles, "excludes_entry_probe": True},
        }
    if case == "advance-probe-count":
        labels = [("3-old", 3), ("8", 8), ("12", 12), ("16", 16), ("24", 24), ("32", 32), ("40", 40), ("48", 48)]
        return {
            "public_calls": [_public_call(20.0, 12.0, cap) for cap in (40.0, 80.0)],
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
        }
    raise UnknownMeasurementClaimCaseError(case)


def _expected_provenance(radial: bool) -> Dict[str, object]:
    sites: List[object] = []
    policy = "advance-native-cap-reporting-observation/v1"
    if radial:
        policy = "radial-known-reporting-driven-selection/v1"
        sites = [
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
        ]
    return {
        "policy": policy,
        "native_cap_decision_site": "compas_cgal._stock_2.engagement_at[2]",
        "engagement_reporting_value_site": "compas_cgal._stock_2.engagement_at[1]",
        "reporting_driven_decision_sites": sites,
    }


def _counts(value: object, field: str) -> Dict[str, object]:
    counts = validate_object(value, ("observations", "accepted", "exceeded"), field)
    observations = validate_integer(counts["observations"], f"{field}.observations")
    accepted = validate_integer(counts["accepted"], f"{field}.accepted")
    exceeded = validate_integer(counts["exceeded"], f"{field}.exceeded")
    if observations != accepted + exceeded:
        fail_payload(field, "observations must equal accepted + exceeded")
    return counts


def _station_case(case: Dict[str, object], field: str) -> None:
    reconstruction = validate_object(
        case["reconstruction"],
        (
            "stock_model",
            "target_centre",
            "centre_decimal_places",
            "target_maximal_radius",
            "radius_decimal_places",
            "occurrence_count",
            "coarse_step",
            "coarse_rungs",
            "refined_radius_sequence",
        ),
        f"{field}.reconstruction",
    )
    validate_literal(reconstruction["stock_model"], "generator-faithful-radial-pre-bridge/v1", f"{field}.reconstruction.stock_model")
    centre = validate_array(reconstruction["target_centre"], f"{field}.reconstruction.target_centre")
    if len(centre) != 2:
        fail_payload(f"{field}.reconstruction.target_centre", "must contain two world coordinates")
    for index, value in enumerate(centre):
        validate_float(value, f"{field}.reconstruction.target_centre[{index}]")
    for name in ("centre_decimal_places", "radius_decimal_places", "occurrence_count"):
        validate_integer(reconstruction[name], f"{field}.reconstruction.{name}")
    for name in ("target_maximal_radius", "coarse_step"):
        validate_float(reconstruction[name], f"{field}.reconstruction.{name}")
    coarse_rungs = validate_array(reconstruction["coarse_rungs"], f"{field}.reconstruction.coarse_rungs")
    for index, value in enumerate(coarse_rungs):
        validate_integer(value, f"{field}.reconstruction.coarse_rungs[{index}]")
    sequence = validate_array(reconstruction["refined_radius_sequence"], f"{field}.reconstruction.refined_radius_sequence")
    for index, value in enumerate(sequence):
        validate_float(value, f"{field}.reconstruction.refined_radius_sequence[{index}]")

    native = validate_object(
        case["native_sampled_decisions"],
        ("rung_6_cap_exceeded", "rung_6_cuts_material", "rung_7_cap_exceeded", "rung_7_cuts_material", "refined_candidates"),
        f"{field}.native_sampled_decisions",
    )
    for name in ("rung_6_cap_exceeded", "rung_6_cuts_material", "rung_7_cap_exceeded", "rung_7_cuts_material"):
        _boolean(native[name], f"{field}.native_sampled_decisions.{name}", optional=True)
    refined = validate_array(native["refined_candidates"], f"{field}.native_sampled_decisions.refined_candidates")
    if len(refined) != len(sequence):
        fail_payload(field, "refined candidates must align with refined radius sequence")
    for index, value in enumerate(refined):
        row = validate_object(value, ("radius", "cap_exceeded", "cuts_material"), f"{field}.native_sampled_decisions.refined_candidates[{index}]")
        validate_float(row["radius"], f"{field}.native_sampled_decisions.refined_candidates[{index}].radius")
        validate_literal(row["radius"], sequence[index], f"{field}.native_sampled_decisions.refined_candidates[{index}].radius")
        _boolean(row["cap_exceeded"], f"{field}.native_sampled_decisions.refined_candidates[{index}].cap_exceeded")
        _boolean(row["cuts_material"], f"{field}.native_sampled_decisions.refined_candidates[{index}].cuts_material")

    reporting = validate_object(
        case["reporting_values"],
        (
            "maximal_radius",
            "rung_6_radius",
            "rung_6_peak",
            "rung_7_radius",
            "rung_7_peak",
            "refined_band_min_radius",
            "refined_band_max_radius",
            "forced_peak",
            "rescued_peak",
            "angle_unit",
            "length_unit",
        ),
        f"{field}.reporting_values",
    )
    for name in (
        "maximal_radius",
        "rung_6_radius",
        "rung_6_peak",
        "rung_7_radius",
        "rung_7_peak",
        "refined_band_min_radius",
        "refined_band_max_radius",
        "forced_peak",
        "rescued_peak",
    ):
        validate_float(reporting[name], f"{field}.reporting_values.{name}", optional=True)
    validate_literal(reporting["angle_unit"], "degree", f"{field}.reporting_values.angle_unit")
    validate_literal(reporting["length_unit"], "mm", f"{field}.reporting_values.length_unit")


def _sweep_case(case: Dict[str, object], field: str, discriminator: str, report_keys: Sequence[str]) -> None:
    reconstruction = validate_object(
        case["reconstruction"],
        ("stock_model", "audit_position_count", "audit_includes_entry_phase", "audit_adds_separate_entry_probe", "excludes_chain_entry_circles", "non_entry_circle_counts"),
        f"{field}.reconstruction",
    )
    validate_literal(reconstruction["stock_model"], "generator-faithful-radial-pre-bridge/v1", f"{field}.reconstruction.stock_model")
    audit_count = validate_integer(reconstruction["audit_position_count"], f"{field}.reconstruction.audit_position_count")
    selected_counts = [
        validate_integer(value, f"{field}.reconstruction.non_entry_circle_counts[{index}]")
        for index, value in enumerate(validate_array(reconstruction["non_entry_circle_counts"], f"{field}.reconstruction.non_entry_circle_counts"))
    ]
    validate_literal(reconstruction["audit_includes_entry_phase"], True, f"{field}.reconstruction.audit_includes_entry_phase")
    validate_literal(reconstruction["audit_adds_separate_entry_probe"], False, f"{field}.reconstruction.audit_adds_separate_entry_probe")
    validate_literal(reconstruction["excludes_chain_entry_circles"], True, f"{field}.reconstruction.excludes_chain_entry_circles")
    config = cast(Dict[str, object], case["config"])
    validate_literal(audit_count, len(cast(List[object], cast(Dict[str, object], config["audit"])["probe_offsets"])), f"{field}.reconstruction.audit_position_count")
    native_rows = validate_array(case["native_sampled_decisions"], f"{field}.native_sampled_decisions")
    report_rows = validate_array(case["reporting_values"], f"{field}.reporting_values")
    expected_values = cast(List[object], config[{"subdivisions": "subdivision_values", "floor_steps": "floor_values", "refinement_margin": "margin_values"}[discriminator]])
    if len(selected_counts) != len(expected_values) or len(native_rows) != len(expected_values) or len(report_rows) != len(expected_values):
        fail_payload(field, "reconstruction/native/reporting rows must align with configuration")
    for index, expected in enumerate(expected_values):
        selected_count = selected_counts[index]
        native = validate_object(native_rows[index], (discriminator, "selected_circles", "observations"), f"{field}.native_sampled_decisions[{index}]")
        validate_literal(native[discriminator], expected, f"{field}.native_sampled_decisions[{index}].{discriminator}")
        validate_literal(
            validate_integer(native["selected_circles"], f"{field}.native_sampled_decisions[{index}].selected_circles"),
            selected_count,
            f"{field}.native_sampled_decisions[{index}].selected_circles",
        )
        observations = _counts(native["observations"], f"{field}.native_sampled_decisions[{index}].observations")
        validate_literal(observations["observations"], selected_count * audit_count, f"{field}.native_sampled_decisions[{index}].observations.observations")
        report = validate_object(report_rows[index], report_keys, f"{field}.reporting_values[{index}]")
        validate_literal(report[discriminator], expected, f"{field}.reporting_values[{index}].{discriminator}")
        validate_float(report["worst_peak"], f"{field}.reporting_values[{index}].worst_peak")
        circles_over_cap = validate_integer(report["circles_over_cap"], f"{field}.reporting_values[{index}].circles_over_cap")
        exceeded_positions = cast(int, observations["exceeded"])
        if circles_over_cap > selected_count or circles_over_cap > exceeded_positions or (circles_over_cap == 0) != (exceeded_positions == 0):
            fail_payload(f"{field}.reporting_values[{index}].circles_over_cap", "must be sound for selected circles and exceeded audit positions")
        validate_literal(report["angle_unit"], "degree", f"{field}.reporting_values[{index}].angle_unit")
        if discriminator == "refinement_margin":
            validate_float(report["cutting_length"], f"{field}.reporting_values[{index}].cutting_length")
            validate_literal(report["length_unit"], "mm", f"{field}.reporting_values[{index}].length_unit")


def _advance_reconstruction(case: Dict[str, object], field: str, audit_count: int) -> Tuple[Dict[str, object], int]:
    reconstruction = validate_object(
        case["reconstruction"],
        ("stock_model", "selected_circle_policy", "original_calls_per_wrapper", "generator_prepends_entry_probe", "audit_excludes_entry_probe", "selected_circle_count"),
        f"{field}.reconstruction",
    )
    validate_literal(reconstruction["stock_model"], "generator-pre-bridge/v1", f"{field}.reconstruction.stock_model")
    validate_literal(reconstruction["selected_circle_policy"], "accepted-non-forced/v1", f"{field}.reconstruction.selected_circle_policy")
    validate_literal(reconstruction["original_calls_per_wrapper"], 1, f"{field}.reconstruction.original_calls_per_wrapper")
    validate_literal(reconstruction["generator_prepends_entry_probe"], True, f"{field}.reconstruction.generator_prepends_entry_probe")
    validate_literal(reconstruction["audit_excludes_entry_probe"], True, f"{field}.reconstruction.audit_excludes_entry_probe")
    selected = validate_integer(reconstruction["selected_circle_count"], f"{field}.reconstruction.selected_circle_count")
    config = cast(Dict[str, object], case["config"])
    validate_literal(audit_count, len(cast(List[object], cast(Dict[str, object], config["audit"])["probe_offsets"])), f"{field}.audit count")
    return reconstruction, selected


def _advance_placement(case: Dict[str, object], field: str) -> None:
    _, total_selected = _advance_reconstruction(case, field, 32)
    config = cast(Dict[str, object], case["config"])
    calls = cast(List[object], config["public_calls"])
    native_rows = validate_array(case["native_sampled_decisions"], f"{field}.native_sampled_decisions")
    reporting_rows = validate_array(case["reporting_values"], f"{field}.reporting_values")
    if len(native_rows) != len(calls) or len(reporting_rows) != len(calls):
        fail_payload(field, "placement rows must align with public calls")
    selected_sum = 0
    for index, call_value in enumerate(calls):
        call = cast(Dict[str, object], call_value)
        native = validate_object(native_rows[index], ("cap", "climb", "selected_circles", "observations"), f"{field}.native_sampled_decisions[{index}]")
        report = validate_object(
            reporting_rows[index],
            ("cap", "climb", "selected_circles", "positions_over_cap", "worst_peak", "worst_peak_offset", "old_probe_peak", "offset_bin_counts", "angle_unit"),
            f"{field}.reporting_values[{index}]",
        )
        for name in ("cap", "climb"):
            validate_literal(native[name], call["tea_cap" if name == "cap" else name], f"{field}.native_sampled_decisions[{index}].{name}")
            validate_literal(report[name], native[name], f"{field}.reporting_values[{index}].{name}")
        selected = validate_integer(native["selected_circles"], f"{field}.native_sampled_decisions[{index}].selected_circles")
        selected_sum += selected
        validate_literal(
            validate_integer(report["selected_circles"], f"{field}.reporting_values[{index}].selected_circles"), selected, f"{field}.reporting_values[{index}].selected_circles"
        )
        observations = _counts(native["observations"], f"{field}.native_sampled_decisions[{index}].observations")
        validate_literal(observations["observations"], selected * 32, f"{field}.native_sampled_decisions[{index}].observations.observations")
        over_cap = validate_integer(report["positions_over_cap"], f"{field}.reporting_values[{index}].positions_over_cap")
        validate_literal(over_cap, observations["exceeded"], f"{field}.reporting_values[{index}].positions_over_cap")
        for name in ("worst_peak", "worst_peak_offset", "old_probe_peak"):
            validate_float(report[name], f"{field}.reporting_values[{index}].{name}")
        validate_literal(report["angle_unit"], "degree", f"{field}.reporting_values[{index}].angle_unit")
        bins = validate_array(report["offset_bin_counts"], f"{field}.reporting_values[{index}].offset_bin_counts")
        offsets: List[object] = []
        bin_sum = 0
        for bin_index, bin_value in enumerate(bins):
            bin_row = validate_object(bin_value, ("offset", "positions_over_cap", "angle_unit"), f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}]")
            validate_float(bin_row["offset"], f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}].offset")
            offsets.append(bin_row["offset"])
            bin_sum += validate_integer(bin_row["positions_over_cap"], f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}].positions_over_cap")
            validate_literal(bin_row["angle_unit"], "degree", f"{field}.reporting_values[{index}].offset_bin_counts[{bin_index}].angle_unit")
        validate_same(offsets, _PLACEMENT_BIN_CENTRES, f"{field}.reporting_values[{index}].offset_bin_counts offsets")
        validate_literal(bin_sum, over_cap, f"{field}.reporting_values[{index}].offset_bin_counts")
    validate_literal(selected_sum, total_selected, f"{field}.reconstruction.selected_circle_count")


def _advance_probe_count(case: Dict[str, object], field: str) -> None:
    _, total_selected = _advance_reconstruction(case, field, 60)
    config = cast(Dict[str, object], case["config"])
    configurations = cast(List[object], config["generator_configurations"])
    calls = cast(List[object], config["public_calls"])
    expected = [
        (cast(Dict[str, object], configuration)["label"], cast(Dict[str, object], configuration)["probe_count"], cast(Dict[str, object], call)["tea_cap"])
        for configuration in configurations
        for call in calls
    ]
    native_rows = validate_array(case["native_sampled_decisions"], f"{field}.native_sampled_decisions")
    reporting_rows = validate_array(case["reporting_values"], f"{field}.reporting_values")
    if len(native_rows) != len(expected) or len(reporting_rows) != len(expected):
        fail_payload(field, "probe-count rows must align with configurations and caps")
    selected_sum = 0
    for index, discriminants in enumerate(expected):
        native = validate_object(native_rows[index], ("label", "probe_count", "cap", "selected_circles", "observations"), f"{field}.native_sampled_decisions[{index}]")
        report = validate_object(reporting_rows[index], ("label", "probe_count", "cap", "worst_peak", "angle_unit"), f"{field}.reporting_values[{index}]")
        for name, expected_value in zip(("label", "probe_count", "cap"), discriminants):
            validate_literal(native[name], expected_value, f"{field}.native_sampled_decisions[{index}].{name}")
            validate_literal(report[name], expected_value, f"{field}.reporting_values[{index}].{name}")
        selected = validate_integer(native["selected_circles"], f"{field}.native_sampled_decisions[{index}].selected_circles")
        selected_sum += selected
        observations = _counts(native["observations"], f"{field}.native_sampled_decisions[{index}].observations")
        validate_literal(observations["observations"], selected * 60, f"{field}.native_sampled_decisions[{index}].observations.observations")
        validate_float(report["worst_peak"], f"{field}.reporting_values[{index}].worst_peak")
        validate_literal(report["angle_unit"], "degree", f"{field}.reporting_values[{index}].angle_unit")
    validate_literal(selected_sum, total_selected, f"{field}.reconstruction.selected_circle_count")


def validate_case(value: object, expected_case: str, index: int) -> Dict[str, object]:
    """Validate one generator case against its canonical case contract.

    Args:
        value: Candidate decoded case.
        expected_case: Canonical case discriminator.
        index: Position in the ordered generator batch.

    Returns:
        The validated case object.

    Raises:
        InvalidMeasurementClaimConfigError: Canonical configuration differs.
        InvalidMeasurementClaimPayloadError: Result structure or values differ.
    """
    field = f"cases[{index}]"
    case = validate_object(value, _COMMON_KEYS, field)
    validate_literal(case["case"], expected_case, f"{field}.case")
    validate_same(case["source_claim_ids"], CASE_CLAIMS[expected_case], f"{field}.source_claim_ids")
    try:
        validate_same(case["config"], _expected_config(expected_case), f"{field}.config")
    except InvalidMeasurementClaimPayloadError as exc:
        raise InvalidMeasurementClaimConfigError(str(exc)) from exc
    radial = expected_case.startswith("radial-")
    validate_same(case["selection_decision_provenance"], _expected_provenance(radial), f"{field}.selection_decision_provenance")
    validate_literal(case["continuous_certificate"], None, f"{field}.continuous_certificate")
    if expected_case == "radial-station":
        _station_case(case, field)
    elif expected_case == "radial-subdivisions":
        _sweep_case(case, field, "subdivisions", ("subdivisions", "worst_peak", "circles_over_cap", "angle_unit"))
    elif expected_case == "radial-floor":
        _sweep_case(case, field, "floor_steps", ("floor_steps", "worst_peak", "circles_over_cap", "angle_unit"))
    elif expected_case == "radial-margin":
        _sweep_case(case, field, "refinement_margin", ("refinement_margin", "worst_peak", "circles_over_cap", "cutting_length", "angle_unit", "length_unit"))
    elif expected_case == "advance-placement":
        _advance_placement(case, field)
    else:
        _advance_probe_count(case, field)
    return case


def validate_generator_cases(value: object) -> List[GeneratorCasePayload]:
    """Validate the one ordered six-case generator collection."""
    case_values = validate_array(value, "cases")
    if len(case_values) != len(CASE_ORDER):
        fail_payload("cases", "must contain exactly six cases")
    for index, case_name in enumerate(CASE_ORDER):
        validate_case(case_values[index], case_name, index)
    return cast(List[GeneratorCasePayload], case_values)


def _wire_tuple(value: object, arity: int, field: str) -> object:
    if type(value) is not tuple:
        return value
    items = cast(Tuple[object, ...], value)
    if len(items) != arity:
        fail_payload(field, f"tuple must contain exactly {arity} coordinates")
    return list(items)


def _normalize_public_call_wire(value: object, field: str) -> None:
    if type(value) is not dict:
        return
    call = cast(Dict[str, object], value)
    polygon = _wire_tuple(call.get("polygon"), 4, f"{field}.polygon")
    if type(polygon) is list:
        points = cast(List[object], polygon)
        polygon = [_wire_tuple(point, 3, f"{field}.polygon[{index}]") for index, point in enumerate(points)]
        call["polygon"] = polygon
    holes = call.get("holes")
    if type(holes) is list:
        for hole_index, hole in enumerate(cast(List[object], holes)):
            if type(hole) is list:
                points = cast(List[object], hole)
                points[:] = [_wire_tuple(point, 3, f"{field}.holes[{hole_index}][{index}]") for index, point in enumerate(points)]


def generator_cases_wire_view(cases: Sequence[GeneratorCasePayload]) -> List[GeneratorCasePayload]:
    """Project typed coordinate tuples to their exact JSON-array wire shape."""
    wire = copy.deepcopy(list(cases))
    for case_index, case in enumerate(wire):
        config = cast(Dict[str, object], case.get("config")) if type(case.get("config")) is dict else None
        if config is not None:
            public_call = config.get("public_call")
            if public_call is not None:
                _normalize_public_call_wire(public_call, f"cases[{case_index}].config.public_call")
            public_calls = config.get("public_calls")
            if type(public_calls) is list:
                for call_index, call in enumerate(cast(List[object], public_calls)):
                    _normalize_public_call_wire(call, f"cases[{case_index}].config.public_calls[{call_index}]")
        if case.get("case") == "radial-station" and type(case.get("reconstruction")) is dict:
            reconstruction = cast(Dict[str, object], case["reconstruction"])
            reconstruction["target_centre"] = _wire_tuple(
                reconstruction.get("target_centre"),
                2,
                f"cases[{case_index}].reconstruction.target_centre",
            )
    return wire
