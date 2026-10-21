"""The batched offset contact projection decides what the per-contact native loop decided.

Oracle: ``_circle_geometry_2.project_boundary_contact`` for the side, the tie
and the bound; an independent ``Fraction`` projection for the parameter, which
must be the correctly rounded double of the exact rational (the native
``to_double`` of a lazy Epeck value is not).
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from hypothesis import event
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_offset_contact_projection import project_bounded_offset_contacts
from benchmarks.held_offset_contact_projection import project_offset_contacts
from compas_cgal import _circle_geometry_2
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

_LATTICE = st.integers(min_value=-12, max_value=12).map(lambda value: value / 4)
_GENERIC = st.floats(min_value=-6.0, max_value=6.0, allow_nan=False)
_BOUND = st.one_of(st.sampled_from((0.0, 0.25, 0.5, 1.0, 2.0, 100.0)), st.floats(min_value=0.0, max_value=8.0, allow_nan=False))


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


@st.composite
def _star_ring(draw: st.DrawFn) -> tuple[tuple[float, float], ...]:
    """A simple ring: strictly increasing angles around the origin, positive radii."""
    count = draw(st.integers(min_value=3, max_value=12))
    angles = sorted(set(draw(st.lists(st.integers(min_value=0, max_value=359), min_size=count, max_size=count))))
    if len(angles) < 3 or any(later - earlier >= 180 for earlier, later in zip(angles, (*angles[1:], angles[0] + 360))):
        angles = [0, 120, 240]
    radii = draw(st.lists(st.one_of(st.sampled_from((1.0, 2.0, 4.0)), st.floats(min_value=0.5, max_value=5.0)), min_size=len(angles), max_size=len(angles)))
    return tuple((radius * math.cos(math.radians(angle)), radius * math.sin(math.radians(angle))) for angle, radius in zip(angles, radii))


@st.composite
def _lattice_rectangle(draw: st.DrawFn) -> tuple[tuple[float, float], ...]:
    """Axis-aligned lattice rectangles make exact distance ties common."""
    left = draw(_LATTICE)
    bottom = draw(_LATTICE)
    width = draw(st.integers(min_value=1, max_value=16)) / 4
    height = draw(st.integers(min_value=1, max_value=16)) / 4
    return ((left, bottom), (left + width, bottom), (left + width, bottom + height), (left, bottom + height))


def _exact_parameter(ring: tuple[tuple[float, float], ...], side: int, query: tuple[float, float]) -> Fraction:
    (ax, ay), (bx, by) = ring[side], ring[(side + 1) % len(ring)]
    edge_x, edge_y = Fraction(bx) - Fraction(ax), Fraction(by) - Fraction(ay)
    along = (Fraction(query[0]) - Fraction(ax)) * edge_x + (Fraction(query[1]) - Fraction(ay)) * edge_y
    return min(Fraction(1), max(Fraction(0), along / (edge_x * edge_x + edge_y * edge_y)))


def _native(ring: tuple[tuple[float, float], ...], query: tuple[float, float], bound: float) -> ProjectionBoundarySite | str:
    try:
        side, _ = _circle_geometry_2.project_boundary_contact(list(ring), query, bound)
    except _circle_geometry_2.BoundaryProjectionDistanceError:
        return "beyond bound"
    except _circle_geometry_2.AmbiguousBoundaryProjectionError:
        return "ambiguous"
    rounded = Fraction(float(_exact_parameter(ring, side, query)))
    return ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), rounded, len(ring))


def _batched(ring: tuple[tuple[float, float], ...], query: tuple[float, float], bound: float) -> ProjectionBoundarySite | str:
    try:
        (site,) = project_offset_contacts((_point(*query),), tuple(_point(*vertex) for vertex in ring), Millimetre(bound))
    except _circle_geometry_2.AmbiguousBoundaryProjectionError:
        return "ambiguous"
    return "beyond bound" if site is None else site


@settings(max_examples=1500, deadline=None)
@given(
    ring=st.one_of(_star_ring(), _lattice_rectangle()),
    query=st.tuples(st.one_of(_LATTICE, _GENERIC), st.one_of(_LATTICE, _GENERIC)),
    bound=_BOUND,
)
def test_batched_projection_matches_native_side_tie_bound_with_correctly_rounded_parameter(
    ring: tuple[tuple[float, float], ...],
    query: tuple[float, float],
    bound: float,
) -> None:
    expected = _native(ring, query, bound)
    event(expected if isinstance(expected, str) else "site")
    assert _batched(ring, query, bound) == expected


@settings(max_examples=300, deadline=None)
@given(ring=_star_ring(), queries=st.lists(st.tuples(_GENERIC, _GENERIC), min_size=1, max_size=20), bound=_BOUND)
def test_one_batch_equals_per_contact_queries(
    ring: tuple[tuple[float, float], ...],
    queries: list[tuple[float, float]],
    bound: float,
) -> None:
    boundary = tuple(_point(*vertex) for vertex in ring)
    try:
        batch = project_offset_contacts(tuple(_point(*query) for query in queries), boundary, Millimetre(bound))
    except _circle_geometry_2.AmbiguousBoundaryProjectionError:
        assert any(_batched(ring, query, bound) == "ambiguous" for query in queries)
        return
    assert batch == tuple(project_offset_contacts((_point(*query),), boundary, Millimetre(bound))[0] for query in queries)


def test_parameter_is_correctly_rounded_where_the_native_view_is_not() -> None:
    """Native returns ``to_double`` of a lazy value, one ulp off here; the batch rounds correctly."""
    ring = ((0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0))
    query = (1.4872182438574346, 0.18530151526760788)
    exact = Fraction(query[0]) / 3
    side, native_parameter = _circle_geometry_2.project_boundary_contact(list(ring), query, 1.0)
    assert side == 0
    assert native_parameter != float(exact)
    (site,) = project_offset_contacts((_point(*query),), tuple(_point(*vertex) for vertex in ring), Millimetre(1.0))
    assert site is not None
    assert site.parameter == Fraction(float(exact))


def test_contract_violations_raise_the_native_error_types() -> None:
    square = tuple(_point(*vertex) for vertex in ((0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)))
    with pytest.raises(_circle_geometry_2.AmbiguousBoundaryProjectionError):
        project_offset_contacts((_point(2.0, 2.0),), square, Millimetre(3.0))
    assert project_offset_contacts((_point(2.0, 2.0),), square, Millimetre(1.0)) == (None,)
    with pytest.raises(_circle_geometry_2.InvalidBoundaryProjectionError):
        project_offset_contacts((_point(2.0, 2.0),), square[:2], Millimetre(1.0))
    with pytest.raises(_circle_geometry_2.InvalidBoundaryProjectionError):
        project_offset_contacts((_point(2.0, 2.0),), square, Millimetre(-1.0))


def _native_loop(ring: tuple[tuple[float, float], ...], queries: list[tuple[float, float]], bound: float) -> object:
    sites = []
    try:
        for query in queries:
            side, _ = _circle_geometry_2.project_boundary_contact(list(ring), query, bound)
            rounded = Fraction(float(_exact_parameter(ring, side, query)))
            sites.append(ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), rounded, len(ring)))
    except (_circle_geometry_2.BoundaryProjectionDistanceError, _circle_geometry_2.AmbiguousBoundaryProjectionError) as error:
        return type(error).__name__
    return tuple(sites)


@settings(max_examples=400, deadline=None)
@given(
    ring=st.one_of(_star_ring(), _lattice_rectangle()),
    queries=st.lists(st.tuples(st.one_of(_LATTICE, _GENERIC), st.one_of(_LATTICE, _GENERIC)), min_size=1, max_size=12),
    bound=_BOUND,
)
def test_bounded_batch_fails_like_the_first_failing_native_call(
    ring: tuple[tuple[float, float], ...],
    queries: list[tuple[float, float]],
    bound: float,
) -> None:
    boundary = tuple(_point(*vertex) for vertex in ring)
    try:
        batched: object = project_bounded_offset_contacts(tuple(_point(*query) for query in queries), boundary, Millimetre(bound))
    except (_circle_geometry_2.BoundaryProjectionDistanceError, _circle_geometry_2.AmbiguousBoundaryProjectionError) as error:
        batched = type(error).__name__
    assert batched == _native_loop(ring, queries, bound)
