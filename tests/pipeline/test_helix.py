"""The exact helix model of the helical entry (Task 1.2, brief Step 1 and the C10 footprint).

Every decision here is exact: the turns rule, the ramp verdict, the footprint
radius floor and the footprint's incidence and tangency are `Fraction` or
native exact predicates. The helix length is REPORTING only.
"""

from __future__ import annotations

import math
from itertools import takewhile
from fractions import Fraction

import pytest
from compas.geometry import Circle
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.gate import gate_pocket
from benchmarks.gate import pipeline
from compas_cgal.incident_loop import INCIDENCE_CHART_DENOMINATOR
from compas_cgal.incident_loop import INCIDENCE_GRID_EXPONENT
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.helix import Helix
from compas_cgal.helix import HelixDescentError
from compas_cgal.pipeline.pencil import HelixFootprintExceedsRootError
from compas_cgal.helix import HelixTurnsError
from compas_cgal.pipeline.pencil import helix_footprint
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

GRID = 2**INCIDENCE_GRID_EXPONENT
GATE_POCKETS = ("rect_12x8", "rect_20x12", "L_shape")
GATE_CAPS_DEG = (120.0, 40.0)

# Probe P3 (W/task-1.2-probe-P34-report.md): depth 4 mm, footprint radius 1 mm.
P3_DEPTH_MM = 4.0
P3_TURNS_BY_ANGLE_DEG = {1.0: 37, 2.0: 19, 3.0: 13, 5.0: 8}
# P3 reported helix length at 2 deg, mm, to the probe's printed 0.01 mm.
P3_LENGTH_AT_2_DEG_MM = 119.45
P3_LENGTH_PRINT_HALF_ULP_MM = 0.005

# Probe P4 witness: chart (0, 16) with r = 1 is the first chart-0 numerator on
# which round() puts the footprint radius above the tool radius.
WITNESS_CHART = (0, 16)


def _root(radius: float, phase: float = 0.0, centre: tuple[float, float] = (0.0, 0.0)) -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY](Millimetre(centre[0]), Millimetre(centre[1])), GuideRadius(Millimetre(radius)), Radian(phase))


def _plane(bottom: float, top: float) -> CutPlane:
    return CutPlane.build(CutZ.build(bottom), ClearanceZ.build(top))


def _unit_footprint() -> IncidentLoop:
    footprint = helix_footprint(_root(3.0), ToolRadius.build(1.0))
    assert Fraction(footprint.radius.value) == 1
    return footprint


def _with_turns(helix: Helix, turns: int) -> Helix:
    return Helix(footprint=helix.footprint, top_z=helix.top_z, bottom_z=helix.bottom_z, turns=turns, clockwise=helix.clockwise)


@pytest.mark.parametrize(("angle_deg", "turns"), sorted(P3_TURNS_BY_ANGLE_DEG.items()))
def test_the_p3_table_turns_hold_at_turns_and_fail_one_below(angle_deg: float, turns: int) -> None:
    limit = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(angle_deg))
    helix = Helix.build(_unit_footprint(), _plane(0.0, P3_DEPTH_MM), limit, clockwise=False)
    assert helix.turns == turns
    assert helix.ramp_within(limit)
    assert not _with_turns(helix, turns - 1).ramp_within(limit)
    assert helix.pitch == Fraction(P3_DEPTH_MM) / turns


# Integer-ratio boundary: K * rho_h divides the descent exactly, so the pitch
# at the least turn count EQUALS the bound. K values are exact binary doubles.
INTEGER_RATIO_CASES = ((1.0, 4), (0.5, 8))


@pytest.mark.parametrize(("pitch_per_radius", "turns"), INTEGER_RATIO_CASES)
def test_an_integer_turn_ratio_holds_at_equality_and_fails_one_below(pitch_per_radius: float, turns: int) -> None:
    limit = HelixPitchPerRadius.build(pitch_per_radius)
    footprint = _unit_footprint()
    assert P3_DEPTH_MM / (pitch_per_radius * footprint.radius.value) == turns
    helix = Helix.build(footprint, _plane(0.0, P3_DEPTH_MM), limit, clockwise=False)
    assert helix.turns == turns
    assert helix.pitch == limit.exact * Fraction(footprint.radius.value)
    assert helix.ramp_within(limit)
    assert not _with_turns(helix, turns - 1).ramp_within(limit)


@settings(max_examples=200, deadline=None)
@given(
    radius=st.floats(min_value=0.05, max_value=200.0),
    phase=st.floats(min_value=-math.pi, max_value=math.pi),
    cut_z=st.floats(min_value=-50.0, max_value=50.0),
    depth=st.floats(min_value=1e-3, max_value=100.0),
    angle_deg=st.floats(min_value=1e-3, max_value=89.9),
)
def test_turns_are_the_least_integer_meeting_the_exact_ramp_surrogate(radius: float, phase: float, cut_z: float, depth: float, angle_deg: float) -> None:
    footprint = _root(radius, phase)
    limit = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(angle_deg))
    helix = Helix.build(footprint, _plane(cut_z, cut_z + depth), limit, clockwise=True)
    assert helix.ramp_within(limit)
    assert helix.turns == 1 or not _with_turns(helix, helix.turns - 1).ramp_within(limit)


def test_exit_xy_is_the_entry_xy_exactly() -> None:
    footprint = _root(2.5, 1.0, (3.25, -1.5))
    limit = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(2.0))
    helix = Helix.build(footprint, _plane(0.0, 4.0), limit, clockwise=False)
    assert (Fraction(helix.exit_xy.x), Fraction(helix.exit_xy.y)) == (Fraction(footprint.entry.x), Fraction(footprint.entry.y))
    assert helix.entry_xy == helix.exit_xy == footprint.entry


def test_the_helix_spans_the_cut_plane_and_reports_its_length() -> None:
    limit = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(2.0))
    helix = Helix.build(_unit_footprint(), _plane(0.0, P3_DEPTH_MM), limit, clockwise=False)
    assert (helix.bottom_z.value, helix.top_z.value) == (0.0, P3_DEPTH_MM)
    assert math.isclose(helix.ramp_length(), P3_LENGTH_AT_2_DEG_MM, rel_tol=0.0, abs_tol=P3_LENGTH_PRINT_HALF_ULP_MM)
    assert helix.length() == helix.ramp_length() + 2.0 * math.pi * helix.footprint.radius.value
    circle = helix.footprint_circle()
    assert (circle.radius, tuple(circle.frame.point)) == (helix.footprint.radius.value, (helix.footprint.centre.x, helix.footprint.centre.y, 0.0))


@pytest.mark.parametrize("turns", [0, -1, True])
def test_a_helix_needs_a_positive_integer_turn_count(turns: int) -> None:
    footprint = _unit_footprint()
    with pytest.raises(HelixTurnsError):
        Helix(footprint=footprint, top_z=ClearanceZ.build(4.0), bottom_z=CutZ.build(0.0), turns=turns, clockwise=False)


def test_a_helix_that_does_not_descend_is_refused() -> None:
    with pytest.raises(HelixDescentError):
        Helix(footprint=_unit_footprint(), top_z=ClearanceZ.build(0.0), bottom_z=CutZ.build(0.0), turns=1, clockwise=False)


def test_the_helix_footprint_radius_is_floored_onto_the_root_grid() -> None:
    chart, numerator = WITNESS_CHART
    root = _root(3.0, 2.0 * math.atan2(numerator, INCIDENCE_CHART_DENOMINATOR))
    assert (root.chart, root.chart_numerator) == WITNESS_CHART
    tool = ToolRadius.build(1.0)
    footprint = helix_footprint(root, tool)
    r = Fraction(tool.value)
    rho_h = Fraction(footprint.radius.value)
    assert rho_h <= r
    # The witness is live: round() on this chart would exceed r.
    h = _chart_hypotenuse(root)
    assert Fraction(round(r * GRID / h) * h, GRID) > r
    assert rho_h == Fraction((r * GRID) // h * h, GRID)
    # The footprint is a pencil member of the root: same entry, internally tangent.
    assert footprint.certifies_incidence()
    assert footprint.entry == root.entry
    dx = Fraction(root.centre.x) - Fraction(footprint.centre.x)
    dy = Fraction(root.centre.y) - Fraction(footprint.centre.y)
    assert dx * dx + dy * dy == (Fraction(root.radius.value) - rho_h) ** 2


def _chart_hypotenuse(root: IncidentLoop) -> int:
    """The chart's reduced hypotenuse ``m``, read from the exact entry direction."""
    dx = Fraction(root.entry.x) - Fraction(root.centre.x)
    dy = Fraction(root.entry.y) - Fraction(root.centre.y)
    rho = Fraction(root.radius.value)
    cos, sin = dx / rho, dy / rho
    return math.lcm(cos.denominator, sin.denominator)


def test_a_root_narrower_than_the_footprint_is_refused() -> None:
    with pytest.raises(HelixFootprintExceedsRootError):
        helix_footprint(_root(0.5), ToolRadius.build(1.0))


def _gate_root(cap_deg: float, pocket: str) -> tuple[IncidentLoop, ToolRadius, Stock, Stock, Helix]:
    """The FIRST_CUT root (the last loop entered at the helix's entry, the pencil's end) and the emitted helix."""
    spec = gate_pocket(pocket, tea_cap_deg=cap_deg)
    operations = pipeline(spec).operations
    index = next(k for k, op in enumerate(operations) if op.operation is OperationType.HELIX)
    helix = operations[index].geometry
    assert isinstance(helix, Helix)
    pencil = list(takewhile(lambda op: isinstance(op.geometry, Circle) and op.entry_xy == operations[index].entry_xy, operations[index + 1 :]))
    assert pencil, "no root loop follows the helix"
    root_op = pencil[-1]
    assert root_op.entry_xy is not None
    circle = root_op.geometry
    centre = (float(circle.frame.point[0]), float(circle.frame.point[1]))
    ex, ey = root_op.entry_xy
    root = _root(float(circle.radius), math.atan2(ey - centre[1], ex - centre[0]), centre)
    assert (root.entry.x, root.entry.y, root.radius.value) == (ex, ey, circle.radius), "gate root chart not reproduced"
    return root, ToolRadius.build(spec.tool_radius), Stock(spec.polygon, list(spec.holes)), Stock(spec.polygon, list(spec.holes)), helix


@pytest.mark.parametrize("cap_deg", GATE_CAPS_DEG)
@pytest.mark.parametrize("pocket", GATE_POCKETS)
def test_a_full_turn_at_the_footprint_is_the_exact_bore_disk(cap_deg: float, pocket: str) -> None:
    root, tool, sweep, disk, helix = _gate_root(cap_deg, pocket)
    footprint = helix_footprint(root, tool)
    assert footprint == helix.footprint, "the emitted helix is not the root's floored footprint"
    bore_radius = Fraction(footprint.radius.value) + Fraction(tool.value)
    assert Fraction(float(bore_radius)) == bore_radius, "bore radius is not an exact double"
    sweep.subtract_circle_sweep_local(footprint.centre.x, footprint.centre.y, footprint.radius.value, tool.value)
    disk.subtract_disk(footprint.centre.x, footprint.centre.y, float(bore_radius))
    assert sweep.exactly_equals(disk)
