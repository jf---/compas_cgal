"""Narrow-station rules: ports of the rho and chain-ordered tests, with injected conditions.

Ported from `tests/test_engagement_rho_toolpath.py:72-183` and
`9fa6eea8^:tests/test_engagement_ordered_toolpath.py:121-173`. Monkeypatches
became injected condition tuples (Ruling 1.1-R13).
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from compas.geometry import Polygon
from hypothesis import given
from hypothesis import strategies as st

from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import EmptyLadderWindowError
from compas_cgal.pipeline.errors import ImplausibleClearanceSlopeError
from compas_cgal.pipeline.errors import SubFloorStationError
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.narrow import MAX_CLEARANCE_SLOPE
from compas_cgal.pipeline.narrow import MAX_LOOP_RADIUS_STEP_TOOL_RADII
from compas_cgal.pipeline.narrow import Advance
from compas_cgal.pipeline.narrow import ChainSplit
from compas_cgal.pipeline.narrow import NarrowCandidate
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.narrow import RadiusStepCondition
from compas_cgal.pipeline.narrow import advance_ladder
from compas_cgal.pipeline.narrow import advance_window
from compas_cgal.pipeline.narrow import clearance_slope
from compas_cgal.pipeline.narrow import radius_step_admissible
from compas_cgal.pipeline.narrow import entry_normal
from compas_cgal.pipeline.narrow import narrow_stations
from compas_cgal.pipeline.narrow import orientations
from compas_cgal.pipeline.narrow import station_loop
from compas_cgal.pipeline.result import DeclineCause
from compas_cgal.pipeline.result import DeclinedRun
from compas_cgal.pipeline.result import ForcedAdvance
from compas_cgal.stock import Stock
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)

# The slope of a 45-degree corner bisector, which the rectangle corner chains present.
BISECTOR_SLOPE = math.sqrt(0.5)

# Synthetic guide spacing, mm.
SYNTHETIC_SPACING = 0.05

# Test-side G1 tolerance on unit-vector dot products (a measurement, not a decision).
G1_DOT_TOLERANCE = 1e-9
UNIT_TOLERANCE = 1e-12

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])


def _station(cx: float, radius: float, clockwise: bool = True, cy: float = 0.0, tx: float = 1.0, ty: float = 0.0) -> GuideStation:
    return GuideStation(
        centre=Point2[WorldXY](Millimetre(cx), Millimetre(cy)),
        radius=GuideRadius.build(radius),
        tangent=Direction2[WorldXY].build(tx, ty),
        clockwise=clockwise,
    )


def _straight_chain(count: int, slope: float, first_radius: float, spacing: float = SYNTHETIC_SPACING, clockwise: bool = True, origin: int = 0) -> GuideChain:
    """A chain along +X whose gouge-free radius climbs at a constant *slope*."""
    return GuideChain(
        stations=tuple(_station(index * spacing, first_radius + slope * index * spacing, clockwise) for index in range(count)),
        origin=origin,
    )


def _oriented(chain: GuideChain) -> OrientedChain:
    oriented = OrientedChain.build(chain, TOOL)
    assert oriented is not None
    return oriented


def _stock() -> Stock:
    return Stock(RECT_12X8)


# ---------------------------------------------------------------------------
# The advance scan.
# ---------------------------------------------------------------------------


def test_the_advance_ladder_returns_the_furthest_admissible_station() -> None:
    chain = _oriented(_straight_chain(8, slope=0.0, first_radius=2.0))
    seen: list[int] = []

    def admissible(_stock: Stock, _origin: NarrowCandidate, candidate: NarrowCandidate) -> bool:
        seen.append(candidate.index)
        return candidate.index in {2, 5}

    advance = advance_ladder(_stock(), chain, 0, 6, (admissible,))
    assert advance.index == 5
    assert advance.forced is None
    assert seen[0] == 6, "the scan must begin at the furthest candidate"


def test_the_advance_ladder_records_a_forced_minimum_advance_when_nothing_is_admissible() -> None:
    chain = _oriented(_straight_chain(8, slope=0.0, first_radius=2.0, origin=3))

    def refuse(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
        return False

    advance = advance_ladder(_stock(), chain, 0, 6, (refuse,))
    assert advance.index == 1, "a forced advance is the minimum one, never a jump"
    assert advance.forced == ForcedAdvance(chain=3, reversed_walk=False, origin_station=0, station=1)


def test_the_advance_ladder_never_steps_over_the_radius_bound() -> None:
    """Spacing far wider than the shipped guide's, where the bound binds."""
    bound = MAX_LOOP_RADIUS_STEP_TOOL_RADII * TOOL.value
    chain = _oriented(_straight_chain(8, slope=0.9, first_radius=2.0, spacing=0.5))
    assert abs(chain.stations[6].radius.value - chain.stations[0].radius.value) > bound, "the fixture must place station 6 out of reach"

    def accept(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
        return True

    advance = advance_ladder(_stock(), chain, 0, 6, (RadiusStepCondition(TOOL), accept))
    assert advance.forced is None
    assert advance.index < 6
    assert abs(chain.stations[advance.index].radius.value - chain.stations[0].radius.value) <= bound


def test_the_forced_advance_keeps_the_radius_bound() -> None:
    # Stations 1 and 2 are out of radius reach of station 0; 3 is back in reach.
    # Spacing 3 mm keeps every slope within the 1-Lipschitz bound.
    stations = [_station(0.0, 2.0), _station(3.0, 4.5), _station(6.0, 4.4), _station(9.0, 3.9), _station(12.0, 3.9)]
    chain = _oriented(GuideChain(stations=tuple(stations), origin=0))

    def refuse(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
        return False

    advance = advance_ladder(_stock(), chain, 0, 4, (refuse,))
    assert advance.index == 3
    assert advance.forced is not None


def test_an_inadmissible_bridge_refuses_a_candidate_whose_circle_complies() -> None:
    chain = _oriented(_straight_chain(8, slope=0.0, first_radius=2.0))

    def circle(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
        return True

    def bridge(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
        return False

    assert advance_ladder(_stock(), chain, 0, 6, (circle, bridge)).forced is not None


def test_the_ladder_never_offers_a_station_below_the_floor() -> None:
    # Station 2 dips to the floor inside the chain: it is never a candidate.
    stations = [_station(0.0, 2.0), _station(1.0, 2.0), _station(2.0, 1.0), _station(3.0, 2.0)]
    chain = _oriented(GuideChain(stations=tuple(stations), origin=0))
    seen: list[int] = []

    def record(_stock: Stock, _origin: NarrowCandidate, candidate: NarrowCandidate) -> bool:
        seen.append(candidate.index)
        return False

    advance = advance_ladder(_stock(), chain, 0, 3, (record,))
    assert seen == [3, 1]
    assert advance.index == 1


# ---------------------------------------------------------------------------
# The entry geometry.
# ---------------------------------------------------------------------------


def test_a_constant_radius_chain_keeps_the_guide_normal_entry() -> None:
    station = _station(0.0, 2.0)
    wx, wy = entry_normal(station, 0.0)
    # Clockwise on +X: the material-side guide normal is (-ty, tx) = (0, 1).
    assert (wx, wy) == pytest.approx((0.0, 1.0), abs=UNIT_TOLERANCE)


@pytest.mark.parametrize("clockwise", [True, False])
@pytest.mark.parametrize("slope", [0.0, 0.25, BISECTOR_SLOPE, -0.25, -BISECTOR_SLOPE])
def test_the_bridge_is_tangent_to_both_circles_at_any_radius_difference(slope: float, clockwise: bool) -> None:
    stations = narrow_stations(_straight_chain(6, slope=slope, first_radius=2.0, clockwise=clockwise))
    for earlier, later in zip(stations[1:-2], stations[2:-1]):
        (sx, sy), (ex, ey) = earlier.entry, later.entry
        length = math.hypot(ex - sx, ey - sy)
        assert length > 0.0
        chord = ((ex - sx) / length, (ey - sy) / length)
        for station in (earlier, later):
            tx, ty = station.entry_tangent
            assert chord[0] * tx + chord[1] * ty == pytest.approx(1.0, abs=G1_DOT_TOLERANCE)


def test_a_slope_a_clearance_function_cannot_have_is_refused() -> None:
    with pytest.raises(ImplausibleClearanceSlopeError):
        narrow_stations(_straight_chain(4, slope=MAX_CLEARANCE_SLOPE * 1.5, first_radius=2.0))


def test_the_entry_direction_stays_a_unit_vector_across_the_slope_range() -> None:
    for slope in (0.0, 0.1, 0.5, BISECTOR_SLOPE, 0.99, MAX_CLEARANCE_SLOPE):
        wx, wy = entry_normal(_station(0.0, 2.0), slope)
        assert math.hypot(wx, wy) == pytest.approx(1.0, abs=UNIT_TOLERANCE)


# ---------------------------------------------------------------------------
# Orientation and the floor.
# ---------------------------------------------------------------------------


def test_a_chain_is_offered_in_both_walk_directions() -> None:
    offered = orientations(_straight_chain(8, slope=0.0, first_radius=3.0, origin=2), TOOL)
    assert len(offered) == 2
    assert {chain.reversed_walk for chain in offered} == {False, True}
    assert all(chain.origin == 2 for chain in offered)
    forward, backward = offered
    assert forward.entry_station.centre != backward.entry_station.centre, "the two directions begin at opposite ends"


def test_a_chain_with_no_trochoidal_station_is_offered_in_no_direction() -> None:
    assert orientations(_straight_chain(6, slope=0.0, first_radius=0.5), TOOL) == ()


def test_a_sub_floor_station_is_never_authored() -> None:
    station = narrow_stations(_straight_chain(1, slope=0.0, first_radius=1.0))[0]
    with pytest.raises(SubFloorStationError):
        station_loop(station, TOOL)


def test_the_floor_is_re_decided_on_the_authored_radius() -> None:
    # Above r as a station, but the authored radius (quantum 2^-40 mm at this
    # axis phase) rounds onto r exactly: the loop would be a bore, so it is declined.
    radius = math.nextafter(1.0, 2.0)
    chain = _straight_chain(1, slope=0.0, first_radius=radius)
    assert chain.stations[0].radius.value > TOOL.value
    assert OrientedChain.build(chain, TOOL) is None


@given(radius=st.floats(1.0, 1.0 + 2e-5))
def test_the_floor_is_decided_on_the_station_then_on_the_authored_radius(radius: float) -> None:
    chain = _straight_chain(1, slope=0.0, first_radius=radius)
    station = narrow_stations(chain)[0]
    expected = radius > TOOL.value and station_loop(station, TOOL).radius.value > TOOL.value
    assert (OrientedChain.build(chain, TOOL) is not None) == expected


def test_declined_runs_are_maximal_and_in_guide_order() -> None:
    radii = [0.5, 0.8, 2.0, 2.0, 0.9, 2.0, 0.7]
    chain = GuideChain(stations=tuple(_station(2.0 * index, value) for index, value in enumerate(radii)), origin=4)
    forward, backward = orientations(chain, TOOL)
    expected = ((4, 0, 1), (4, 4, 4), (4, 6, 6))
    for oriented in (forward, backward):
        assert tuple((run.chain, run.first_station, run.last_station) for run in oriented.declined) == expected
    assert (forward.first, forward.last) == (2, 5)


# ---------------------------------------------------------------------------
# Rulings 1.1-B-boundinvariant, B-unbridgeable, B-slope-exact, B-minors.
# ---------------------------------------------------------------------------


def _accept(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
    return True


def test_injecting_accept_alone_still_refuses_a_bound_breaking_candidate() -> None:
    """B-boundinvariant: the radius-step bound is structural; a caller cannot drop it."""
    bound = MAX_LOOP_RADIUS_STEP_TOOL_RADII * TOOL.value
    chain = _oriented(_straight_chain(8, slope=0.9, first_radius=2.0, spacing=0.5))
    advance = advance_ladder(_stock(), chain, 0, 6, (_accept,))
    assert isinstance(advance, Advance)
    assert abs(chain.stations[advance.index].radius.value - chain.stations[0].radius.value) <= bound


def test_a_floor_dip_longer_than_the_window_splits_the_chain_instead_of_raising() -> None:
    """B-unbridgeable: no floor-clearing station in the window -> a split, a gap record, a re-entry station."""
    window = 3
    radii = [2.0, 2.0, 0.5, 0.5, 0.5, 0.5, 2.0, 2.0]
    chain = _oriented(GuideChain(stations=tuple(_station(1.6 * index, value) for index, value in enumerate(radii)), origin=5))
    split = advance_ladder(_stock(), chain, 1, 1 + window, (_accept,))
    assert split == ChainSplit(origin=1, resume=6, gap=DeclinedRun(chain=5, first_station=1, last_station=6, cause=DeclineCause.RADIUS_STEP))


def test_no_station_keeping_the_bound_splits_instead_of_breaking_it() -> None:
    """B-unbridgeable: rho's forced origin+1 would open an uncut patch; the port splits there."""
    # Stations 1-3 all sit more than 2r above station 0; spacing 3 mm keeps the slope plausible.
    stations = [_station(0.0, 1.5), _station(3.0, 4.0), _station(6.0, 4.2), _station(9.0, 4.4)]
    chain = _oriented(GuideChain(stations=tuple(stations), origin=0))
    split = advance_ladder(_stock(), chain, 0, 3, (_accept,))
    assert isinstance(split, ChainSplit)
    assert (split.origin, split.resume) == (0, 1)
    assert not radius_step_admissible(chain.stations[0].radius, chain.stations[1].radius, TOOL)


def test_the_radius_step_bound_admits_its_closed_end() -> None:
    """B-minors: |rho' - rho| == 2r is admissible (``<=``), one ulp more is not."""
    assert radius_step_admissible(GuideRadius.build(1.0), GuideRadius.build(3.0), TOOL)
    assert not radius_step_admissible(GuideRadius.build(1.0), GuideRadius.build(math.nextafter(3.0, 4.0)), TOOL)


def test_an_empty_ladder_window_is_refused_by_name() -> None:
    chain = _oriented(_straight_chain(4, slope=0.0, first_radius=2.0))
    with pytest.raises(EmptyLadderWindowError):
        advance_ladder(_stock(), chain, 2, 2, (_accept,))
    with pytest.raises(EmptyLadderWindowError):
        advance_ladder(_stock(), chain, chain.last, chain.last + 1, (_accept,))


def test_forced_advances_and_declined_runs_share_guide_order() -> None:
    """B-minors: one station numbering for the ledger, whatever the walk direction."""
    forward = _straight_chain(8, slope=0.0, first_radius=2.0, origin=3)
    backward = OrientedChain.build(forward.reversed(), TOOL)
    assert backward is not None

    def refuse(_stock: Stock, _origin: NarrowCandidate, _candidate: NarrowCandidate) -> bool:
        return False

    advance = advance_ladder(_stock(), backward, 0, 6, (refuse,))
    assert isinstance(advance, Advance)
    assert advance.forced == ForcedAdvance(chain=3, reversed_walk=True, origin_station=7, station=6)


@given(
    ra=st.floats(0.5, 5.0),
    rb=st.floats(0.5, 5.0),
    dx=st.floats(-0.2, 0.2),
    dy=st.floats(-0.2, 0.2),
)
def test_the_slope_guard_is_decided_exactly(ra: float, rb: float, dx: float, dy: float) -> None:
    """B-slope-exact: refused iff (rho_b - rho_a)^2 > |c_b - c_a|^2 over the rationals of the doubles."""
    before, after = _station(0.0, ra), _station(dx, rb, cy=dy)
    rise = Fraction(rb) - Fraction(ra)
    implausible = rise * rise > Fraction(dx) ** 2 + Fraction(dy) ** 2
    if implausible:
        with pytest.raises(ImplausibleClearanceSlopeError):
            clearance_slope(before, after)
    else:
        clearance_slope(before, after)


def test_the_window_is_rhos_floor_division() -> None:
    """B-window: 1.0 // 0.025 == 39 stations, kept verbatim from the frozen generators."""
    assert advance_window(0.025) == 39


def test_the_slope_guard_admits_its_exact_bound_and_refuses_one_ulp_past_it() -> None:
    """B-slope-exact: centres (0, 0) -> (0.75, 1) are 1.25 apart exactly; a rise of 1.25 is the bound itself."""
    before = _station(0.0, 1.0)
    assert clearance_slope(before, _station(0.75, 2.25, cy=1.0)) == MAX_CLEARANCE_SLOPE
    with pytest.raises(ImplausibleClearanceSlopeError):
        clearance_slope(before, _station(0.75, math.nextafter(2.25, 3.0), cy=1.0))
