"""The arrangement-free rim engine ``AnnulusView2`` against the arrangement engine it must equal (unit U1 of the actual-stock refinement design).

Every test compares with ``Stock2`` replayed from the same part polygon and the
same closed annuli, or pins a verdict derived by hand. ``station_matches_stock``
and ``circle_matches_stock`` compare exact run endpoints, stronger than the
reported doubles. Integer fixtures reach the measure-zero contacts (tangency,
rim overlap, rim through a vertex); generic doubles are the witness that the
engine does not stall on them.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _stock_2

CAP_RATIO = _stock_2.cap_chord_ratio(math.radians(80.0))
LOOP_PARAMETERS = [k / 45.0 for k in range(45)]
# Star directions of the Hypothesis polygons: any positive radii give a simple polygon.
STAR_DIRECTIONS = ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))
# A generic-double decision that stalls is the failure this bound exists for (seconds, wall, under any load).
GENERIC_DECISION_WALL_S = 30.0
# Two runs equal at their exact endpoints report spans that differ only by to_double of different one-root representations (radians).
REPORTED_SPAN_AGREEMENT_RAD = 1e-9


def rows(points: list[tuple[float, float]]) -> np.ndarray:
    return np.array([(x, y, 0.0) for x, y in points], dtype=np.float64)


def annulus_rows(annuli: list[tuple[float, float, float]]) -> np.ndarray:
    return np.array(annuli, dtype=np.float64).reshape(-1, 3)


def square(half: float) -> np.ndarray:
    return rows([(-half, -half), (half, -half), (half, half), (-half, half)])


def replay(polygon: np.ndarray, annuli: list[tuple[float, float, float]], tool_radius: float) -> _stock_2.Stock2:
    stock = _stock_2.Stock2(polygon, [])
    for cx, cy, rho in annuli:
        stock.subtract_circle_sweep(cx, cy, rho, tool_radius)
    return stock


def view_of(polygon: np.ndarray, annuli: list[tuple[float, float, float]], tool_radius: float) -> _stock_2.AnnulusView2:
    return _stock_2.AnnulusView2.build(polygon, annulus_rows(annuli), tool_radius)


def assert_station_equal(
    polygon: np.ndarray, annuli: list[tuple[float, float, float]], tool_radius: float, x: float, y: float, cap: float = CAP_RATIO, gap: float = 0.0
) -> tuple[float, float, bool]:
    stock = replay(polygon, annuli, tool_radius)
    view = view_of(polygon, annuli, tool_radius)
    engine = view.station(0, len(annuli), x, y, cap, gap)
    reference = _stock_2.engagement_at(stock, x, y, tool_radius, cap, gap)
    assert engine[2] == reference[2]
    assert engine[0] == pytest.approx(reference[0], abs=REPORTED_SPAN_AGREEMENT_RAD)
    assert engine[1] == pytest.approx(reference[1], abs=REPORTED_SPAN_AGREEMENT_RAD)
    assert view.station_matches_stock(stock, 0, len(annuli), x, y, cap, gap)
    return engine


# ---- tangency ----------------------------------------------------------------------


def test_rim_tangent_to_outer_circle_from_outside_splits_nothing() -> None:
    total, max_run, exceeded = assert_station_equal(square(20.0), [(3.0, 0.0, 1.0)], 1.0, 0.0, 0.0)
    assert total == pytest.approx(2.0 * math.pi)
    assert max_run == pytest.approx(2.0 * math.pi)
    assert exceeded


def test_rim_tangent_to_inner_circle_inside_the_hole_splits_nothing() -> None:
    total, _max_run, exceeded = assert_station_equal(square(20.0), [(0.0, 0.0, 4.0)], 1.0, 2.0, 0.0)
    assert total == pytest.approx(2.0 * math.pi)
    assert exceeded


def test_rim_tangent_to_outer_circle_from_inside_the_annulus_has_no_material() -> None:
    total, _max_run, exceeded = assert_station_equal(square(20.0), [(0.0, 0.0, 3.0)], 1.0, 3.0, 0.0)
    assert total == 0.0
    assert not exceeded


def test_rim_tangent_to_a_wall_from_inside_merges_across_the_contact() -> None:
    total, max_run, exceeded = assert_station_equal(square(10.0), [], 1.0, 9.0, 0.0)
    assert total == pytest.approx(2.0 * math.pi)
    assert max_run == pytest.approx(2.0 * math.pi)
    assert exceeded


def test_rim_tangent_to_a_wall_from_outside_has_no_material() -> None:
    total, _max_run, exceeded = assert_station_equal(square(10.0), [], 1.0, 11.0, 0.0)
    assert total == 0.0
    assert not exceeded


def test_two_tangent_annuli_merge_at_a_shared_rim_point() -> None:
    assert_station_equal(square(20.0), [(3.0, 0.0, 1.0), (-3.0, 0.0, 1.0), (0.0, 3.0, 1.0)], 1.0, 0.0, 0.0)


# ---- rim overlap with an annulus boundary ------------------------------------------


def test_rim_overlapping_the_inner_circle_is_closed_annulus_not_material() -> None:
    total, _max_run, exceeded = assert_station_equal(square(20.0), [(0.0, 0.0, 2.0)], 1.0, 0.0, 0.0)
    assert total == 0.0
    assert not exceeded


def test_rim_overlapping_the_outer_circle_is_closed_annulus_not_material() -> None:
    total, _max_run, exceeded = assert_station_equal(square(20.0), [(0.0, 0.0, 0.0)], 1.0, 0.0, 0.0)
    assert total == 0.0
    assert not exceeded


def test_rim_overlap_cut_by_a_wall_has_no_material() -> None:
    total, _max_run, exceeded = assert_station_equal(square(10.0), [(9.5, 0.0, 2.0)], 1.0, 9.5, 0.0)
    assert total == 0.0
    assert not exceeded


def test_guide_radius_equal_to_tool_radius_is_a_full_disk() -> None:
    total, _max_run, exceeded = assert_station_equal(square(20.0), [(0.0, 0.0, 1.0)], 1.0, 0.25, 0.0)
    assert total == 0.0
    assert not exceeded


def test_duplicate_annuli_are_one_region() -> None:
    assert_station_equal(square(20.0), [(0.5, 0.0, 1.5), (0.5, 0.0, 1.5)], 1.0, 0.0, 1.0)


# ---- rim through vertices ----------------------------------------------------------


def test_rim_through_a_convex_vertex() -> None:
    assert_station_equal(rows([(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]), [], 1.0, 4.0, 1.0)


def test_rim_through_a_reflex_vertex_keeps_material_on_both_sides() -> None:
    polygon = rows([(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (0.0, 0.0), (-5.0, 5.0)])
    assert_station_equal(polygon, [], 1.0, 0.0, -1.0)


# ---- station 0 at a wall (the population of review M3) -----------------------------


@pytest.mark.parametrize(
    ("entry_x", "expected_exceeded"),
    [
        pytest.param(math.nextafter(9.0, 10.0), False, id="crosses_the_wall_by_one_ulp_and_splits_the_run"),
        pytest.param(9.0, True, id="tangent_keeps_one_run"),
        pytest.param(math.nextafter(9.0, 0.0), True, id="clear_of_the_wall_by_one_ulp_keeps_one_run"),
    ],
)
def test_sub_ulp_wall_crossing_is_real_stock(entry_x: float, expected_exceeded: bool) -> None:
    # A disk of squared radius 1 + 1 + 2cos(50 deg) about (entry - 1, 0) leaves about +-50 deg of rim: 100 deg merged
    # (over the 80 deg cap), 50 deg per side when the wall event splits it (ruling M3: no tolerance, no near-tangent class).
    rho = 0.8125
    annuli = [(entry_x - 1.0, 0.0, rho)]
    _total, _max_run, exceeded = assert_station_equal(square(10.0), annuli, 1.0, entry_x, 0.0)
    assert exceeded is expected_exceeded


def test_generic_double_wall_station_decides_quickly_and_equals_the_arrangement() -> None:
    # A rotated wall with generic doubles: the near-tangent wall event of real stations, with no integer structure.
    angle = 0.7390851332151607
    ux, uy = math.cos(angle), math.sin(angle)
    half = math.e * 3.0
    polygon = rows([(ux * a - uy * b + math.pi, uy * a + ux * b + 1.0 / 3.0) for a, b in ((-half, -half), (half, -half), (half, half), (-half, half))])
    radius = math.sqrt(2.0) / 2.0
    entry_x = ux * (half - radius) + math.pi
    entry_y = uy * (half - radius) + 1.0 / 3.0
    annuli = [(entry_x - ux * radius, entry_y - uy * radius, radius * 0.8125)]
    start = time.perf_counter()
    assert_station_equal(polygon, annuli, radius, entry_x, entry_y)
    assert time.perf_counter() - start < GENERIC_DECISION_WALL_S


# ---- the material domain -----------------------------------------------------------


def test_view_on_part_refuses_where_the_offset_polygon_accepts() -> None:
    # Only the r-wide wall band carries material at this station: a disk removes the rim's left part, and the rest of the rim is outside the offset.
    annuli = [(8.0, 0.0, 0.5)]
    part = view_of(square(10.0), annuli, 1.0).station(0, 1, 9.0, 0.0, CAP_RATIO, 0.0)
    offset = view_of(square(9.0), annuli, 1.0).station(0, 1, 9.0, 0.0, CAP_RATIO, 0.0)
    assert part[2]
    assert not offset[2]
    assert offset[0] == 0.0


# ---- the view range ----------------------------------------------------------------


def test_range_selects_the_emission_window() -> None:
    # The station's rim lies inside the first annulus (radii 2..4) and clear of the second.
    annuli = [(0.0, 0.0, 3.0), (15.0, 15.0, 3.0)]
    view = view_of(square(20.0), annuli, 1.0)
    assert view.size == 2
    assert view.station(0, 0, 3.0, 0.0, CAP_RATIO, 0.0)[2]
    assert not view.station(0, 1, 3.0, 0.0, CAP_RATIO, 0.0)[2]
    assert view.station(1, 2, 3.0, 0.0, CAP_RATIO, 0.0)[2]


@pytest.mark.parametrize(("first", "end"), [(1, 0), (0, 2), (3, 3)])
def test_invalid_range_raises(first: int, end: int) -> None:
    view = view_of(square(20.0), [(0.0, 0.0, 0.0)], 1.0)
    with pytest.raises(_stock_2.InvalidAnnulusViewRangeError):
        view.station(first, end, 0.5, 0.0, CAP_RATIO, 0.0)


# ---- input contract: where the arrangement engine has no verdict, both raise --------


def test_non_simple_polygon_raises_on_both_engines() -> None:
    bowtie = rows([(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)])
    with pytest.raises(ValueError):
        _stock_2.Stock2(bowtie, [])
    with pytest.raises(_stock_2.InvalidAnnulusViewPolygonError):
        view_of(bowtie, [], 1.0)


def test_too_few_vertices_raise_on_both_engines() -> None:
    degenerate = rows([(0.0, 0.0), (1.0, 1.0)])
    with pytest.raises(ValueError):
        _stock_2.Stock2(degenerate, [])
    with pytest.raises(_stock_2.InvalidAnnulusViewPolygonError):
        view_of(degenerate, [], 1.0)


def test_non_finite_annulus_raises_on_both_engines() -> None:
    stock = _stock_2.Stock2(square(5.0), [])
    with pytest.raises(ValueError):
        stock.subtract_circle_sweep(math.nan, 0.0, 1.0, 1.0)
    with pytest.raises(_stock_2.NonFiniteAnnulusViewInputError):
        view_of(square(5.0), [(math.nan, 0.0, 1.0)], 1.0)


@pytest.mark.parametrize(("rho", "tool"), [(-1.0, 1.0), (1.0, 0.0), (1.0, -1.0)])
def test_invalid_radii_raise_on_both_engines(rho: float, tool: float) -> None:
    stock = _stock_2.Stock2(square(5.0), [])
    with pytest.raises(ValueError):
        stock.subtract_circle_sweep(0.0, 0.0, rho, tool)
    with pytest.raises(_stock_2.InvalidAnnulusViewRadiiError):
        view_of(square(5.0), [(0.0, 0.0, rho)], tool)


def test_non_finite_station_raises_on_both_engines() -> None:
    stock = _stock_2.Stock2(square(5.0), [])
    view = view_of(square(5.0), [], 1.0)
    with pytest.raises(ValueError):
        _stock_2.engagement_at(stock, math.inf, 0.0, 1.0, CAP_RATIO)
    with pytest.raises(ValueError):
        view.station(0, 0, math.inf, 0.0, CAP_RATIO, 0.0)


# ---- whole circles: stations, own prefix, entry -------------------------------------


def assert_circle_equal(
    polygon: np.ndarray, annuli: list[tuple[float, float, float]], tool_radius: float, circle: tuple[float, float, float, float, float], clockwise: bool = False
) -> list[tuple[float, float, bool]]:
    stock = replay(polygon, annuli, tool_radius)
    view = view_of(polygon, annuli, tool_radius)
    cx, cy, rho, ex, ey = circle
    engine = view.circle_prefixes(0, len(annuli), cx, cy, rho, ex, ey, clockwise, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    reference = _stock_2.engagement_circle_prefixes(stock, cx, cy, rho, ex, ey, clockwise, LOOP_PARAMETERS, tool_radius, CAP_RATIO, 0.0)
    assert len(engine) == len(reference) == len(LOOP_PARAMETERS)
    for ours, theirs in zip(engine, reference):
        assert ours[2] == theirs[2]
        assert ours[0] == pytest.approx(theirs[0], abs=REPORTED_SPAN_AGREEMENT_RAD)
        assert ours[1] == pytest.approx(theirs[1], abs=REPORTED_SPAN_AGREEMENT_RAD)
    decision = view.circle_exceeded(0, len(annuli), cx, cy, rho, ex, ey, clockwise, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    assert decision == any(sample[2] for sample in reference)
    assert decision == _stock_2.engagement_circle_exceeded(stock, cx, cy, rho, ex, ey, clockwise, LOOP_PARAMETERS, tool_radius, CAP_RATIO, 0.0)
    assert all(verdict and runs for verdict, runs in view.circle_matches_stock(stock, 0, len(annuli), cx, cy, rho, ex, ey, clockwise, LOOP_PARAMETERS, CAP_RATIO, 0.0))
    return engine


def test_circle_after_a_predecessor_equals_the_arrangement() -> None:
    assert_circle_equal(square(30.0), [(0.0, 0.0, 3.0)], 1.0, (1.0, 0.0, 3.0, 4.0, 0.0))


def test_clockwise_circle_equals_the_arrangement() -> None:
    assert_circle_equal(square(30.0), [(0.0, 0.0, 3.0)], 1.0, (1.0, 0.0, 3.0, 4.0, 0.0), clockwise=True)


def test_circle_on_a_wall_equals_the_arrangement() -> None:
    assert_circle_equal(square(10.0), [(5.0, 0.0, 3.0)], 1.0, (5.5, 0.0, 3.5, 9.0, 0.0))


def test_station_zero_is_the_entry_not_a_point_on_the_circle() -> None:
    # The entry is tangent to the wall, 0.5 past the guide circle. Only the entry reads the wall band; the same material is clear of a rim on the circle.
    annuli = [(8.0, 0.0, 0.5)]
    circle = (5.0, 0.0, 3.5, 9.0, 0.0)
    samples = assert_circle_equal(square(10.0), annuli, 1.0, circle)
    at_entry = view_of(square(10.0), annuli, 1.0).station(0, 1, 9.0, 0.0, CAP_RATIO, 0.0)
    assert samples[0] == pytest.approx(at_entry)
    on_circle = view_of(square(10.0), annuli, 1.0).station(0, 1, 8.5, 0.0, CAP_RATIO, 0.0)
    assert samples[0][1] != pytest.approx(on_circle[1], abs=1e-3)


# ---- H1: random polygons and annulus sequences ---------------------------------------

integer_annuli = st.lists(st.tuples(st.integers(-9, 9), st.integers(-9, 9), st.integers(0, 6)), max_size=6)
half_integer_annuli = st.lists(st.tuples(st.integers(-18, 18), st.integers(-18, 18), st.integers(0, 12)), max_size=6)
generic = st.floats(min_value=-9.0, max_value=9.0, allow_nan=False, allow_infinity=False, allow_subnormal=False)
generic_annuli = st.lists(st.tuples(generic, generic, st.floats(min_value=0.0, max_value=6.0, allow_subnormal=False)), max_size=6)


def star(radii: list[float]) -> np.ndarray:
    return rows([(dx * radius, dy * radius) for (dx, dy), radius in zip(STAR_DIRECTIONS, radii)])


integer_radii = st.lists(st.integers(1, 12), min_size=8, max_size=8)
generic_radii = st.lists(st.floats(min_value=0.5, max_value=12.0, allow_subnormal=False), min_size=8, max_size=8)


@settings(max_examples=250, deadline=None)
@given(integer_radii, integer_annuli, st.integers(1, 3), st.integers(-12, 12), st.integers(-12, 12))
def test_h1_integer_fixtures_station_equals_the_arrangement(radii: list[int], annuli: list[tuple[int, int, int]], tool: int, x: int, y: int) -> None:
    assert_station_equal(star([float(radius) for radius in radii]), [(float(a), float(b), float(c)) for a, b, c in annuli], float(tool), float(x), float(y))


@settings(max_examples=200, deadline=None)
@given(integer_radii, half_integer_annuli, st.integers(-24, 24), st.integers(-24, 24))
def test_h1_half_integer_fixtures_station_equals_the_arrangement(radii: list[int], annuli: list[tuple[int, int, int]], x: int, y: int) -> None:
    assert_station_equal(star([float(radius) for radius in radii]), [(a / 2.0, b / 2.0, c / 2.0) for a, b, c in annuli], 1.0, x / 2.0, y / 2.0)


@settings(max_examples=200, deadline=None)
@given(generic_radii, generic_annuli, st.floats(min_value=0.25, max_value=3.0, allow_subnormal=False), generic, generic)
def test_h1_generic_doubles_station_equals_the_arrangement(radii: list[float], annuli: list[tuple[float, float, float]], tool: float, x: float, y: float) -> None:
    assert_station_equal(star(radii), annuli, tool, x, y)


@settings(max_examples=60, deadline=None)
@given(integer_radii, integer_annuli, st.integers(-9, 9), st.integers(-9, 9), st.integers(1, 6), st.booleans())
def test_h1_integer_fixtures_circle_equals_the_arrangement(radii: list[int], annuli: list[tuple[int, int, int]], x: int, y: int, rho: int, clockwise: bool) -> None:
    circle = (float(x), float(y), float(rho), float(x + rho), float(y))
    assert_circle_equal(star([float(radius) for radius in radii]), [(float(a), float(b), float(c)) for a, b, c in annuli], 1.0, circle, clockwise)


@settings(max_examples=40, deadline=None)
@given(
    generic_radii,
    generic_annuli,
    generic,
    generic,
    st.floats(min_value=0.5, max_value=6.0, allow_subnormal=False),
    st.floats(min_value=0.0, max_value=6.28, allow_subnormal=False),
)
def test_h1_generic_doubles_circle_equals_the_arrangement(radii: list[float], annuli: list[tuple[float, float, float]], x: float, y: float, rho: float, phase: float) -> None:
    entry_x, entry_y = x + rho * math.cos(phase), y + rho * math.sin(phase)
    assume((entry_x, entry_y) != (x, y))
    assert_circle_equal(star(radii), annuli, 1.0, (x, y, rho, entry_x, entry_y))


@settings(max_examples=120, deadline=None)
@given(integer_radii, integer_annuli, st.integers(-12, 12), st.integers(-12, 12), st.floats(min_value=0.0, max_value=4.0, allow_subnormal=False))
def test_h1_gap_closure_equals_the_arrangement(radii: list[int], annuli: list[tuple[int, int, int]], x: int, y: int, gap: float) -> None:
    assert_station_equal(star([float(radius) for radius in radii]), [(float(a), float(b), float(c)) for a, b, c in annuli], 1.0, float(x), float(y), gap=gap)


# ---- large coordinates: the witness construction's exact refinement ------------------


def test_cells_finer_than_the_double_steering_resolve_exactly() -> None:
    # Far from the origin, to_double of an event relative to the centre resolves about 1e-9 rad while two annulus boundaries a few ulps apart cut a thinner cell.
    far = 2.0**40
    outer = 2.0
    near = math.nextafter(outer, 3.0)
    polygon = square(2.0**41)
    annuli = [(far + 3.0, far, outer - 1.0), (far + 3.0, far + 2.0**-12, near - 1.0)]
    assert_station_equal(polygon, annuli, 1.0, far, far)
    assert_station_equal(polygon, annuli, 1.0, math.nextafter(far, far * 2.0), far)
