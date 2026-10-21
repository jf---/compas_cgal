"""The integer offset-point lerp reports exactly what the Fraction lerp reported.

The oracle is the pre-change ``_canonical_offset_point``, kept verbatim: it
injects both side endpoints with ``Fraction.from_float``, lerps exactly and
rounds once with ``float``. The integer path must give the same correctly
rounded double bit for bit (``==`` plus sign of zero), not approximately.
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_boundary_path import AmbiguousFigure5BoundaryPathError
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_boundary_path import _dyadic_lerp
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

SMALLEST_SUBNORMAL = math.ulp(0.0)
# Coordinates stay far from overflow: the lerp of two finite doubles is finite.
COORDINATE_LIMIT = 1e300

_SUBNORMAL = st.integers(min_value=-(2**52) + 1, max_value=2**52 - 1).map(lambda count: count * SMALLEST_SUBNORMAL)
_COORDINATE = st.one_of(
    st.floats(min_value=-COORDINATE_LIMIT, max_value=COORDINATE_LIMIT, allow_nan=False, allow_infinity=False),
    st.integers(min_value=-64, max_value=64).map(lambda value: value / 8),
    _SUBNORMAL,
    st.sampled_from((0.0, -0.0)),
)
_DENOMINATOR = st.one_of(st.integers(min_value=1, max_value=2**64), st.integers(min_value=0, max_value=1100).map(lambda power: 2**power))
_PARAMETER = st.one_of(
    st.sampled_from((Fraction(0), Fraction(1), Fraction(1, 2))),
    _DENOMINATOR.flatmap(lambda denominator: st.integers(min_value=0, max_value=denominator).map(lambda numerator: Fraction(numerator, denominator))),
)


def _fraction_lerp(start: float, end: float, parameter: Fraction) -> float:
    exact_start = Fraction.from_float(start)
    exact_end = Fraction.from_float(end)
    return float(exact_start + parameter * (exact_end - exact_start))


def _old_canonical_offset_point(
    component: tuple[Point2[WorldXY], ...],
    site: ProjectionBoundarySite,
) -> Point2[WorldXY]:
    """Verbatim pre-change ``_canonical_offset_point`` (helpers inlined)."""
    side_index = int(site.segment_id)
    if side_index >= len(component):
        raise AmbiguousFigure5BoundaryPathError("A Figure 5 hypothesis names no radius-one offset side.")
    first = component[side_index]
    second = component[(side_index + 1) % len(component)]
    start = Fraction.from_float(float(first.x)), Fraction.from_float(float(first.y))
    end = Fraction.from_float(float(second.x)), Fraction.from_float(float(second.y))
    point = (
        start[0] + site.parameter * (end[0] - start[0]),
        start[1] + site.parameter * (end[1] - start[1]),
    )
    return Point2[WorldXY].build(float(point[0]), float(point[1]))


def _same_double(first: float, second: float) -> bool:
    return first == second and math.copysign(1.0, first) == math.copysign(1.0, second)


@settings(max_examples=2000, deadline=None)
@given(_COORDINATE, _COORDINATE, _PARAMETER)
@example(0.0, 1.0, Fraction(0))
@example(0.0, 1.0, Fraction(1))
@example(-0.0, -0.0, Fraction(1, 3))
@example(SMALLEST_SUBNORMAL, -SMALLEST_SUBNORMAL, Fraction(1, 2))
@example(SMALLEST_SUBNORMAL, 2 * SMALLEST_SUBNORMAL, Fraction(1, 2))
@example(1.0, 1.0 + math.ulp(1.0), Fraction(1, 2))
@example(1.0, 1.0 + math.ulp(1.0), Fraction(3, 4))
@example(-COORDINATE_LIMIT, COORDINATE_LIMIT, Fraction(1, 3))
def test_dyadic_lerp_equals_the_fraction_lerp(start: float, end: float, parameter: Fraction) -> None:
    assert _same_double(_dyadic_lerp(start, end, parameter), _fraction_lerp(start, end, parameter))


@st.composite
def _sites(draw: st.DrawFn) -> tuple[tuple[Point2[WorldXY], ...], ProjectionBoundarySite]:
    vertices = tuple(Point2[WorldXY].build(draw(_COORDINATE), draw(_COORDINATE)) for _ in range(draw(st.integers(min_value=3, max_value=6))))
    segment = draw(st.integers(min_value=0, max_value=len(vertices) - 1))
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(segment), draw(_PARAMETER), len(vertices))
    return vertices, site


@settings(max_examples=500, deadline=None)
@given(_sites())
def test_canonical_offset_point_equals_the_pre_change_fraction_path(case: tuple[tuple[Point2[WorldXY], ...], ProjectionBoundarySite]) -> None:
    component, site = case
    new = _canonical_offset_point(component, site)
    old = _old_canonical_offset_point(component, site)
    assert repr(new) == repr(old)
    assert _same_double(float(new.x), float(old.x)) and _same_double(float(new.y), float(old.y))


def test_canonical_offset_point_rejects_a_missing_side() -> None:
    component = (Point2[WorldXY].build(0.0, 0.0), Point2[WorldXY].build(1.0, 0.0), Point2[WorldXY].build(0.0, 1.0))
    site = ProjectionBoundarySite(ProjectionBoundarySegmentId(3), Fraction(1, 2), None)
    with pytest.raises(AmbiguousFigure5BoundaryPathError):
        _canonical_offset_point(component, site)
