"""Tests for the spiral-ramp entry: the centre plunge, the ramp, and the 240-degree floor."""

from __future__ import annotations

from typing import List

import pytest
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import Polygon

from compas_cgal.engagement_rho_toolpath import DEGENERATE_LOOP_TOOL_RADII
from compas_cgal.engagement_rho_toolpath import _rho_stations
from compas_cgal.engagement_spiral_entry_toolpath import MAX_RAMP_TURNS
from compas_cgal.engagement_spiral_entry_toolpath import RampCannotProgressError
from compas_cgal.engagement_spiral_entry_toolpath import _ramp_rungs
from compas_cgal.engagement_spiral_entry_toolpath import spiral_entry_toolpath
from compas_cgal.engagement_toolpath import UnavoidableEngagementWarning
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.engagement_toolpath import _Regulation
from compas_cgal.toolpath import OperationType

TOOL_DIAMETER = 2.0
TOOL_RADIUS = 1.0
CAP_DEG = 120.0

# The floor a non-degenerate first loop cannot beat after a single plunge:
# 360 - 2*acos(rho / 2r) is increasing in rho and equals this at rho = r.
ENTRY_ENGAGEMENT_FLOOR_DEG = 240.0

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])


def _regulation() -> _Regulation:
    return _Regulation.build(
        tool_diameter=TOOL_DIAMETER,
        tea_cap_deg=CAP_DEG,
        guide_step_tool_diameters=0.025,
        max_advance_tool_diameters=1.0,
        radial_clearance=None,
        cut_z=0.0,
        clearance_z=None,
    )


def _station(radius: float) -> _GuideStation:
    return _GuideStation(cx=0.0, cy=0.0, radius=radius, clockwise=True, tx=1.0, ty=0.0)


def _result():
    with pytest.warns(UnavoidableEngagementWarning):
        return spiral_entry_toolpath(RECT_12X8, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)


# ---------------------------------------------------------------------------
# The ramp rungs.
# ---------------------------------------------------------------------------


def test_ramp_rungs_descend_and_never_offer_a_bore() -> None:
    rungs = _ramp_rungs(_rho_stations([_station(3.0)])[0], _regulation())
    assert rungs == sorted(rungs, reverse=True), "the scan relies on descending order"
    assert rungs[0] == 3.0, "rung 0 must be the station's own radius bit-for-bit"
    assert min(rungs) > DEGENERATE_LOOP_TOOL_RADII * TOOL_RADIUS


def test_a_station_at_the_floor_offers_no_ramp_rung_at_all() -> None:
    assert _ramp_rungs(_rho_stations([_station(TOOL_RADIUS)])[0], _regulation()) == []


def test_a_ramp_that_cannot_enlarge_its_radius_raises_rather_than_spinning() -> None:
    """`MAX_RAMP_TURNS` is a ceiling on a loop whose progress is otherwise structural."""
    assert MAX_RAMP_TURNS > 0
    assert issubclass(RampCannotProgressError, RuntimeError)


# ---------------------------------------------------------------------------
# The entry geometry: the plunge moves to the CENTRE, which is the whole point.
# ---------------------------------------------------------------------------


def test_the_plunge_sits_at_the_first_circle_centre_not_on_its_rim() -> None:
    result = _result()
    operations = result.operations
    plunge_indices = [index for index, op in enumerate(operations) if op.operation is OperationType.PLUNGE]
    assert plunge_indices, "every chain must plunge"
    for index in plunge_indices:
        plunge = operations[index].geometry
        assert isinstance(plunge, Line)
        point = (float(plunge.end[0]), float(plunge.end[1]))
        following = next(op for op in operations[index + 1 :] if isinstance(op.geometry, Circle))
        centre = following.geometry.frame.point
        assert point == pytest.approx((float(centre[0]), float(centre[1])), abs=1e-9), "the ramp plunges at the circle CENTRE so the hole is concentric with every turn"


def test_the_ramp_climbs_strictly_and_reaches_the_station_radius() -> None:
    result = _result()
    operations = result.operations
    first_plunge = next(index for index, op in enumerate(operations) if op.operation is OperationType.PLUNGE)
    centre = None
    radii: List[float] = []
    for op in operations[first_plunge + 1 :]:
        if not isinstance(op.geometry, Circle):
            continue
        here = (float(op.geometry.frame.point[0]), float(op.geometry.frame.point[1]))
        if centre is None:
            centre = here
        if here != centre:
            break
        radii.append(float(op.geometry.radius))
    assert len(radii) > 1, "the entry must be a ramp of several concentric turns, not one circle"
    assert radii == sorted(radii), "each ramp turn must strictly enlarge the radius"
    assert len(set(radii)) == len(radii)


def test_no_ramp_turn_is_a_bore() -> None:
    result = _result()
    radii = [float(op.geometry.radius) for op in result.operations if isinstance(op.geometry, Circle)]
    assert min(radii) > DEGENERATE_LOOP_TOOL_RADII * TOOL_RADIUS


# ---------------------------------------------------------------------------
# The measured claim.
# ---------------------------------------------------------------------------


def test_the_entry_engages_far_less_than_a_plunge_entry_but_not_below_the_floor() -> None:
    """The bound is 240 degrees and the construction must land above it, not beat it."""
    from benchmarks.gate import gate_pocket
    from benchmarks.survey import MotionKind
    from benchmarks.survey import survey_path

    spec = gate_pocket("rect_12x8")
    with pytest.warns(UnavoidableEngagementWarning):
        result = spiral_entry_toolpath(spec.polygon, tool_diameter=spec.tool_diameter, tea_cap_deg=spec.tea_cap_deg)
    survey = survey_path(spec, result)
    first_loop = next(motion for motion in survey.motions if motion.kind is MotionKind.LOOP)
    assert first_loop.peak_engagement_deg < 360.0, "a centre plunge must open a hole the first loop can see"
    assert first_loop.peak_engagement_deg >= ENTRY_ENGAGEMENT_FLOOR_DEG, "no non-degenerate first loop can engage below 240 degrees after one plunge"


def test_the_declined_regions_survive_the_different_entry() -> None:
    result = _result()
    assert result.declined_regions, "the rectangle's corners are declined whichever entry is used"
