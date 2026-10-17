"""Exact arc predicates: the radius tolerance and the swept-direction test."""

from fractions import Fraction

from hypothesis import given
from hypothesis import strategies as st

from compas_cgal.post.arc_geometry import Vector
from compas_cgal.post.arc_geometry import radius_within_tolerance
from compas_cgal.post.arc_geometry import sweep_contains
from compas_cgal.post.grid import quantise
from compas_cgal.post.quantities import DecimalPlaces

TOLERANCE = Fraction(2, 1000)
DECIMALS = DecimalPlaces(4)


def test_radius_gap_exactly_at_tolerance_passes() -> None:
    assert radius_within_tolerance(Fraction(1), (1 + TOLERANCE) ** 2, TOLERANCE)


def test_radius_gap_beyond_tolerance_fails() -> None:
    assert not radius_within_tolerance(Fraction(1), (1 + TOLERANCE + Fraction(1, 10**9)) ** 2, TOLERANCE)
    assert not radius_within_tolerance((1 + TOLERANCE + Fraction(1, 10**9)) ** 2, Fraction(1), TOLERANCE)


def test_squared_difference_form_would_reject_a_rounding_level_arc() -> None:
    # r = 9 mm, ends 0.0005 mm apart in radius: |rs^2 - re^2| = 0.009 >> tol^2 = 4e-6, yet (rs - re)^2 <= tol^2.
    rs2, re2 = Fraction(9) ** 2, (Fraction(9) + Fraction(5, 10000)) ** 2
    assert abs(rs2 - re2) > TOLERANCE**2
    assert radius_within_tolerance(rs2, re2, TOLERANCE)


def _on_circle(t: Fraction, radius: Fraction) -> Vector:
    """An exact point of the circle of *radius* about the origin (rational parametrisation)."""
    denominator = 1 + t * t
    return Vector(radius * (1 - t * t) / denominator, radius * 2 * t / denominator)


coordinates = st.fractions(min_value=-200, max_value=200, max_denominator=10**6)
parameters = st.fractions(min_value=-50, max_value=50, max_denominator=10**4)
radii = st.fractions(min_value=Fraction(1, 10), max_value=100, max_denominator=10**4)


@given(coordinates, coordinates, radii, parameters, parameters)
def test_quantised_exact_arcs_pass(cx: Fraction, cy: Fraction, radius: Fraction, ts: Fraction, te: Fraction) -> None:
    centre = Vector(quantise(cx, DECIMALS), quantise(cy, DECIMALS))
    start, end = _on_circle(ts, radius), _on_circle(te, radius)
    qs = Vector(quantise(cx + start.x, DECIMALS), quantise(cy + start.y, DECIMALS))
    qe = Vector(quantise(cx + end.x, DECIMALS), quantise(cy + end.y, DECIMALS))
    rs2 = (qs.x - centre.x) ** 2 + (qs.y - centre.y) ** 2
    re2 = (qe.x - centre.x) ** 2 + (qe.y - centre.y) ** 2
    assert radius_within_tolerance(rs2, re2, TOLERANCE)


# Smallest radius step the boundary test straddles, mm: far below the 1e-4 mm quantum.
BOUNDARY_STEP = Fraction(1, 10**12)


@given(coordinates, coordinates, radii, parameters, parameters)
def test_tolerance_boundary_is_exact(cx: Fraction, cy: Fraction, radius: Fraction, ts: Fraction, te: Fraction) -> None:
    centre = Vector(cx, cy)
    start = _on_circle(ts, radius)
    rs2 = start.x**2 + start.y**2
    for grown, inside in ((radius + TOLERANCE, True), (radius + TOLERANCE + BOUNDARY_STEP, False), (radius - TOLERANCE, True), (radius - TOLERANCE - BOUNDARY_STEP, False)):
        end = _on_circle(te, grown)
        moved = Vector(centre.x + end.x, centre.y + end.y)
        re2 = (moved.x - centre.x) ** 2 + (moved.y - centre.y) ** 2
        assert radius_within_tolerance(rs2, re2, TOLERANCE) is inside, (grown, inside)


EAST, NORTH, WEST, SOUTH = Vector(Fraction(1), Fraction(0)), Vector(Fraction(0), Fraction(1)), Vector(Fraction(-1), Fraction(0)), Vector(Fraction(0), Fraction(-1))


def test_sweep_quarter_turns() -> None:
    assert sweep_contains(EAST, NORTH, Vector(Fraction(1), Fraction(1)), clockwise=False)
    assert not sweep_contains(EAST, NORTH, WEST, clockwise=False)
    assert sweep_contains(EAST, NORTH, WEST, clockwise=True)
    assert sweep_contains(EAST, NORTH, SOUTH, clockwise=True)


def test_sweep_half_turns() -> None:
    assert sweep_contains(EAST, WEST, NORTH, clockwise=False)
    assert not sweep_contains(EAST, WEST, SOUTH, clockwise=False)
    assert sweep_contains(EAST, WEST, SOUTH, clockwise=True)
    assert not sweep_contains(EAST, WEST, NORTH, clockwise=True)
