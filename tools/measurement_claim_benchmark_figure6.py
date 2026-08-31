"""Strict JSON grammar for the fixed Figure-6 benchmark result."""

from __future__ import annotations

import math
from typing import Dict
from typing import List
from typing import Tuple
from typing import cast

from tools.measurement_claim_benchmark_semantic_input import FIGURE6_CAPS
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_COMPARISON_SPACING
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_FINE_SPACING
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SPACINGS
from tools.measurement_claim_case_validation import fail_payload
from tools.measurement_claim_case_validation import validate_array
from tools.measurement_claim_case_validation import validate_float
from tools.measurement_claim_case_validation import validate_integer
from tools.measurement_claim_case_validation import validate_literal
from tools.measurement_claim_case_validation import validate_object
from tools.measurement_claim_case_validation import validate_same
from tools.measurement_claim_schema import Degrees
from tools.measurement_claim_schema import ToolDiameters

_ROOT_KEYS = ("pocket", "engagement_measured_at_cap_deg", "points", "spacing_trials")
_POCKET_KEYS = ("name", "family", "tool_diameter", "params", "holes")
_POINT_KEYS = (
    "cap_deg",
    "controlled_length",
    "controlled_cut_motions",
    "controlled_entry_cuts",
    "controlled_max_tea_deg",
    "controlled_max_tea_after_entry_deg",
    "controlled_exceedances_after_entry",
    "controlled_meets_cap",
    "mathsm_spacing_tool_diameters",
    "mathsm_length",
    "mathsm_max_tea_after_entry_deg",
    "length_ratio",
)
_TRIAL_KEYS = ("spacing_tool_diameters", "length", "cut_motions", "entry_cuts", "max_tea_deg", "max_tea_after_entry_deg")


def _validate_holes(value: object) -> List[object]:
    holes = validate_array(value, "figure6.json.pocket.holes")
    for ring_index, ring_value in enumerate(holes):
        ring = validate_array(ring_value, f"figure6.json.pocket.holes[{ring_index}]")
        for point_index, point_value in enumerate(ring):
            field = f"figure6.json.pocket.holes[{ring_index}][{point_index}]"
            point = validate_array(point_value, field)
            if len(point) != 3:
                fail_payload(field, "must contain exactly three coordinates")
            for coordinate_index, coordinate in enumerate(point):
                if type(coordinate) not in (int, float) or not math.isfinite(cast(float, coordinate)):
                    fail_payload(f"{field}[{coordinate_index}]", "must be finite numeric")
    return holes


def validate_figure6_json(figure6_payload: object) -> Tuple[str, Degrees, Degrees, int, int]:
    """Validate the fixed Figure-6 JSON and select the MC-013 observations.

    Args:
        figure6_payload: Decoded raw Figure-6 JSON.

    Returns:
        Pocket name, selected angles, point count, and trial count.

    Raises:
        InvalidMeasurementClaimPayloadError: Structure, config, or selected
            measurement ordering violates the fixed contract.
    """
    root = validate_object(figure6_payload, _ROOT_KEYS, "figure6.json")
    pocket = validate_object(root["pocket"], _POCKET_KEYS, "figure6.json.pocket")
    validate_literal(pocket["name"], "rect_20x12", "figure6.json.pocket.name")
    validate_literal(pocket["family"], "analytic", "figure6.json.pocket.family")
    validate_literal(pocket["tool_diameter"], 2.0, "figure6.json.pocket.tool_diameter")
    validate_same(pocket["params"], {"width": 20.0, "height": 12.0}, "figure6.json.pocket.params")
    holes = _validate_holes(pocket["holes"])
    if holes:
        fail_payload("figure6.json.pocket.holes", "fixed Figure-6 config requires exactly []")
    validate_literal(root["engagement_measured_at_cap_deg"], 180.0, "figure6.json.engagement_measured_at_cap_deg")
    points = validate_array(root["points"], "figure6.json.points")
    if len(points) != len(FIGURE6_CAPS):
        fail_payload("figure6.json.points", "requires the exact cap count")
    for index, (item, cap) in enumerate(zip(points, FIGURE6_CAPS)):
        point = validate_object(item, _POINT_KEYS, f"figure6.json.points[{index}]")
        validate_literal(point["cap_deg"], cap, f"figure6.json.points[{index}].cap_deg")
        for name in ("controlled_length", "controlled_max_tea_deg", "controlled_max_tea_after_entry_deg"):
            validate_float(point[name], f"figure6.json.points[{index}].{name}")
        for name in ("controlled_cut_motions", "controlled_entry_cuts", "controlled_exceedances_after_entry"):
            validate_integer(point[name], f"figure6.json.points[{index}].{name}")
        if type(point["controlled_meets_cap"]) is not bool:
            fail_payload(f"figure6.json.points[{index}].controlled_meets_cap", "must be an exact boolean")
        for name in ("mathsm_spacing_tool_diameters", "mathsm_length", "mathsm_max_tea_after_entry_deg", "length_ratio"):
            validate_float(point[name], f"figure6.json.points[{index}].{name}", optional=True)
    trials = validate_array(root["spacing_trials"], "figure6.json.spacing_trials")
    if len(trials) != len(FIGURE6_SPACINGS):
        fail_payload("figure6.json.spacing_trials", "requires the exact spacing count")
    selected: Dict[ToolDiameters, Degrees] = {}
    for index, (item, spacing) in enumerate(zip(trials, FIGURE6_SPACINGS)):
        trial = validate_object(item, _TRIAL_KEYS, f"figure6.json.spacing_trials[{index}]")
        validate_literal(trial["spacing_tool_diameters"], spacing, f"figure6.json.spacing_trials[{index}].spacing_tool_diameters")
        for name in ("length", "max_tea_deg", "max_tea_after_entry_deg"):
            validate_float(trial[name], f"figure6.json.spacing_trials[{index}].{name}")
        validate_integer(trial["cut_motions"], f"figure6.json.spacing_trials[{index}].cut_motions")
        validate_integer(trial["entry_cuts"], f"figure6.json.spacing_trials[{index}].entry_cuts")
        if spacing in (FIGURE6_FINE_SPACING, FIGURE6_COMPARISON_SPACING):
            selected[spacing] = Degrees(cast(float, trial["max_tea_after_entry_deg"]))
    fine = selected[FIGURE6_FINE_SPACING]
    comparison = selected[FIGURE6_COMPARISON_SPACING]
    if fine <= comparison:
        fail_payload("figure6.json.spacing_trials", "requires unique 0.025 > 0.1 selected measurements")
    return cast(str, pocket["name"]), fine, comparison, len(points), len(trials)
