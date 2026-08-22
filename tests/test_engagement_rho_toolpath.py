"""Tests for the rho-regulated generator: the scan direction, the entry geometry, the floor."""

from __future__ import annotations

import math
from typing import List

import pytest
from compas.geometry import Circle
from compas.geometry import Polygon

from compas_cgal import engagement_rho_toolpath as rho_module
from compas_cgal.engagement_rho_toolpath import DEGENERATE_LOOP_TOOL_RADII
from compas_cgal.engagement_rho_toolpath import MAX_CLEARANCE_SLOPE
from compas_cgal.engagement_rho_toolpath import MAX_LOOP_RADIUS_STEP_TOOL_RADII
from compas_cgal.engagement_rho_toolpath import ImplausibleClearanceSlopeError
from compas_cgal.engagement_rho_toolpath import _entry_normal
from compas_cgal.engagement_rho_toolpath import _largest_admissible_advance
from compas_cgal.engagement_rho_toolpath import _rho_stations
from compas_cgal.engagement_rho_toolpath import rho_regulated_toolpath
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.engagement_toolpath import _Regulation

TOOL_DIAMETER = 2.0
TOOL_RADIUS = 1.0
CAP_DEG = 120.0

# A chain running along +X whose radius climbs at this rate per unit of travel.
# 0.7071 is the slope of a 45-degree corner bisector, the case the rectangle
# corner chains actually present.
BISECTOR_SLOPE = math.sqrt(0.5)

# Guide station spacing used by the synthetic chains, in model units.
SYNTHETIC_SPACING = 0.05

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


def _straight_chain(count: int, slope: float, first_radius: float, spacing: float = SYNTHETIC_SPACING) -> List[_GuideStation]:
    """A chain along +X whose gouge-free radius climbs at a constant *slope*."""
    return [
        _GuideStation(
            cx=index * spacing,
            cy=0.0,
            radius=first_radius + slope * index * spacing,
            clockwise=True,
            tx=1.0,
            ty=0.0,
        )
        for index in range(count)
    ]


# ---------------------------------------------------------------------------
# The scan direction. The team-lead brief calls this out explicitly: the bug this
# family of searches has already had once was an INVERTED rung choice.
# ---------------------------------------------------------------------------


def test_the_advance_scan_returns_the_furthest_admissible_station(monkeypatch) -> None:
    """A non-monotone pass/fail pattern: only a downward scan finds the furthest pass.

    Admissible set is ``{2, 5}`` out of a window reaching 6. The furthest is 5.
    A bisection on ``[1, 6]`` probes 3 (fail), then 1..2, and returns 2. An
    inverted scan returns 2 as well. Only a scan from the top returns 5.
    """
    stations = _rho_stations(_straight_chain(8, slope=0.0, first_radius=2.0))
    admissible = {2, 5}
    seen: List[int] = []

    def fake(_stock, guide_station, _advance, _tool_radius, _cap_ratio) -> bool:
        index = round(guide_station.cx / SYNTHETIC_SPACING)
        seen.append(index)
        return index in admissible

    monkeypatch.setattr(rho_module, "_station_is_admissible", fake)
    index, forced = _largest_admissible_advance(None, stations, 0, 6, _regulation())

    assert index == 5, f"expected the furthest admissible station 5, got {index}"
    assert forced is False
    # The scan must have started at the far end, not the near one.
    assert seen[0] == 6, f"the scan must begin at the furthest candidate; it began at {seen[0]}"


def test_the_advance_scan_reports_forced_when_nothing_is_admissible(monkeypatch) -> None:
    stations = _rho_stations(_straight_chain(8, slope=0.0, first_radius=2.0))
    monkeypatch.setattr(rho_module, "_station_is_admissible", lambda *_args, **_kwargs: False)
    index, forced = _largest_admissible_advance(None, stations, 0, 6, _regulation())
    assert forced is True
    assert index == 1, "a forced advance is the minimum one, never a jump"


def test_the_advance_scan_never_steps_over_the_radius_bound(monkeypatch) -> None:
    """The furthest station is admissible by the cap but out of radius reach.

    The spacing here is deliberately far wider than the shipped guide's. At the
    shipped defaults the bound CANNOT bind: a clearance function is 1-Lipschitz,
    so over the widest advance the search may consider --
    ``MAX_ADVANCE_TOOL_DIAMETERS * tool_diameter``, which is ``2 * r`` -- the
    radius cannot change by more than ``2 * r``, which is the bound itself. The
    rule is therefore a guard that binds only on a coarser guide or a wider
    advance window, and this test is what keeps it working for that case.
    """
    regulation = _regulation()
    bound = MAX_LOOP_RADIUS_STEP_TOOL_RADII * regulation.tool_radius
    spacing = 0.5
    slope = 0.9
    stations = _rho_stations(_straight_chain(8, slope=slope, first_radius=2.0, spacing=spacing))
    assert abs(stations[6].radius - stations[0].radius) > bound, "the fixture must place station 6 out of reach"
    monkeypatch.setattr(rho_module, "_station_is_admissible", lambda *_args, **_kwargs: True)
    index, forced = _largest_admissible_advance(None, stations, 0, 6, regulation)
    assert forced is False
    assert index < 6, "the scan must refuse the furthest station on the radius bound alone"
    assert abs(stations[index].radius - stations[0].radius) <= bound


# ---------------------------------------------------------------------------
# The entry geometry.
# ---------------------------------------------------------------------------


def test_a_constant_radius_chain_keeps_the_guide_normal_entry() -> None:
    """At zero slope the external-tangent entry IS the advance-only generator's entry."""
    station = _GuideStation(cx=0.0, cy=0.0, radius=2.0, clockwise=True, tx=1.0, ty=0.0)
    entry_x, entry_y = _entry_normal(station, 0.0)
    guide_entry_x, guide_entry_y = station.entry
    assert entry_x == pytest.approx((guide_entry_x - station.cx) / station.radius, abs=1e-12)
    assert entry_y == pytest.approx((guide_entry_y - station.cy) / station.radius, abs=1e-12)


@pytest.mark.parametrize("clockwise", [True, False])
@pytest.mark.parametrize("slope", [0.0, 0.25, BISECTOR_SLOPE, -0.25, -BISECTOR_SLOPE])
def test_the_bridge_is_tangent_to_both_circles_at_any_radius_difference(slope: float, clockwise: bool) -> None:
    """The whole point of the external-tangent entry: G1 at both ends of every bridge."""
    chain = _straight_chain(6, slope=slope, first_radius=2.0)
    if not clockwise:
        chain = [_GuideStation(cx=s.cx, cy=s.cy, radius=s.radius, clockwise=False, tx=s.tx, ty=s.ty) for s in chain]
    stations = _rho_stations(chain)
    for earlier, later in zip(stations[1:-2], stations[2:-1]):
        start_x, start_y = earlier.entry
        end_x, end_y = later.entry
        length = math.hypot(end_x - start_x, end_y - start_y)
        assert length > 0.0
        chord = ((end_x - start_x) / length, (end_y - start_y) / length)
        for station in (earlier, later):
            tangent = station.entry_tangent
            dot = chord[0] * tangent[0] + chord[1] * tangent[1]
            assert dot == pytest.approx(1.0, abs=1e-9), f"junction is not G1: dot={dot} at slope={slope} clockwise={clockwise}"


def test_a_slope_a_clearance_function_cannot_have_is_refused() -> None:
    station = _GuideStation(cx=0.0, cy=0.0, radius=2.0, clockwise=True, tx=1.0, ty=0.0)
    with pytest.raises(ImplausibleClearanceSlopeError):
        _entry_normal(station, MAX_CLEARANCE_SLOPE * 1.5)


def test_the_entry_direction_stays_a_unit_vector_across_the_slope_range() -> None:
    station = _GuideStation(cx=0.0, cy=0.0, radius=2.0, clockwise=True, tx=1.0, ty=0.0)
    for slope in (0.0, 0.1, 0.5, BISECTOR_SLOPE, 0.99, MAX_CLEARANCE_SLOPE):
        entry_x, entry_y = _entry_normal(station, slope)
        assert math.hypot(entry_x, entry_y) == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# End to end.
# ---------------------------------------------------------------------------


def _loop_radii(result) -> List[float]:
    return [float(op.geometry.radius) for op in result.operations if isinstance(op.geometry, Circle)]


def test_no_emitted_machining_circle_is_a_bore() -> None:
    """The degeneracy floor, asserted on the pocket whose corners break the older generators."""
    with pytest.warns(rho_module.UnavoidableEngagementWarning):
        result = rho_regulated_toolpath(RECT_12X8, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)
    radii = _loop_radii(result)
    assert radii, "the pocket must produce machining circles"
    floor = DEGENERATE_LOOP_TOOL_RADII * TOOL_RADIUS
    assert min(radii) > floor, f"emitted a bore: smallest radius {min(radii)} against a floor of {floor}"


def test_consecutive_machining_circles_keep_their_annuli_overlapping() -> None:
    with pytest.warns(rho_module.UnavoidableEngagementWarning):
        result = rho_regulated_toolpath(RECT_12X8, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)
    radii = _loop_radii(result)
    bound = MAX_LOOP_RADIUS_STEP_TOOL_RADII * TOOL_RADIUS
    steps = [abs(later - earlier) for earlier, later in zip(radii, radii[1:])]
    assert max(steps) <= bound, f"radius step {max(steps)} exceeds the annulus-overlap bound {bound}"


def test_the_generator_emits_the_documented_operation_vocabulary() -> None:
    from compas_cgal.toolpath import OperationType

    with pytest.warns(rho_module.UnavoidableEngagementWarning):
        result = rho_regulated_toolpath(RECT_12X8, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)
    kinds = {op.operation for op in result.operations}
    assert OperationType.PLUNGE in kinds
    assert OperationType.CUT in kinds
    assert OperationType.RETRACT in kinds
    assert result.polyline.shape[1] == 3


# ---------------------------------------------------------------------------
# The declined-region contract. Silent under-cut is the one failure a roughing
# generator must not have, so the material it refuses is DATA, not a warning.
# ---------------------------------------------------------------------------


def test_declined_regions_name_the_material_the_generator_refused() -> None:
    with pytest.warns(rho_module.UnavoidableEngagementWarning):
        result = rho_regulated_toolpath(RECT_12X8, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)
    regions = result.declined_regions
    assert regions, "the rectangle's four corners are below the floor and must be declared"
    for region in regions:
        assert region.station_count > 0
        assert region.largest_gouge_free_radius <= DEGENERATE_LOOP_TOOL_RADII * TOOL_RADIUS, "a declined run must be below the floor by definition"
        assert len(region.first_center) == 2 and len(region.last_center) == 2


def test_a_declined_region_is_a_maximal_run_not_a_span_hiding_machined_stations() -> None:
    """Two separated sub-threshold runs on one chain must come back as two regions.

    The stations are built directly rather than through `_rho_stations`: the radius
    pattern this needs is a step function, whose slope no clearance function could
    have, and `_entry_normal` rightly refuses it. Region detection reads only the
    radius and the centre, so it is exercised on exactly what it consumes.
    """
    tiny = 0.5 * TOOL_RADIUS
    stations = [
        rho_module._RhoStation(
            cx=index * SYNTHETIC_SPACING,
            cy=0.0,
            radius=(tiny if index in (0, 1, 7, 8) else 3.0),
            clockwise=True,
            tx=1.0,
            ty=0.0,
            wx=0.0,
            wy=1.0,
        )
        for index in range(9)
    ]
    regions = rho_module._declined_regions(stations, path_index=3, tool_radius=TOOL_RADIUS)
    assert len(regions) == 2
    assert [region.station_count for region in regions] == [2, 2]
    assert all(region.path_index == 3 for region in regions)


def test_a_fully_machined_chain_declines_nothing() -> None:
    stations = _rho_stations(_straight_chain(6, slope=0.0, first_radius=3.0))
    assert rho_module._declined_regions(stations, path_index=0, tool_radius=TOOL_RADIUS) == []
