"""The tier-1b candidate filter decides exactly, without libm: ``|c - m| <= rho + 2 r + h``.

``m`` is the sub-piece midpoint and ``h`` its half-length (the filter's
documented criterion). The filter must return exactly the prior circles
meeting it, on the gate's grid; the oracle here is rational arithmetic
(squaring, no square root). ``hypot`` is
not correctly rounded on every platform (RS: ``np.hypot(3.6163379554499997,
48.30545762211739)`` is one ulp above ``math.hypot``, whose square already
exceeds the exact ``x² + y²``), so a filter built on it cannot be exact.
"""

from __future__ import annotations

import math
from fractions import Fraction

from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_connector_engagement import _CandidateIndex
from benchmarks.held_connector_engagement import _DyadicGrid
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
RS_OFFSET = (3.6163379554499997, 48.30545762211739)


def _p(x: float, y: float) -> Point2[WorldXY]:
    point: Point2[WorldXY] = Point2[WorldXY].build(x, y)
    return point


def _circle(index: int, centre: tuple[float, float], rho: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(index, 10_000), 4)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, _p(*centre), _p(*centre), GuideRadius.build(rho))


def _oracle(circles: list[tuple[tuple[float, float], float]], a: tuple[float, float], b: tuple[float, float]) -> list[int]:
    ax, ay, bx, by = (Fraction(v) for v in (*a, *b))
    mx, my = (ax + bx) / 2, (ay + by) / 2
    half_sq = ((bx - ax) ** 2 + (by - ay) ** 2) / 4
    kept = []
    for k in range(len(circles) - 1, -1, -1):
        (cx, cy), rho = circles[k]
        reach = Fraction(rho) + 2 * Fraction(float(TOOL.value))
        excess = (Fraction(cx) - mx) ** 2 + (Fraction(cy) - my) ** 2 - reach**2 - half_sq
        if excess <= 0 or excess**2 <= 4 * reach**2 * half_sq:
            kept.append(k)
    return kept


def _near(circles: list[tuple[tuple[float, float], float]], a: tuple[float, float], b: tuple[float, float]) -> list[int]:
    emitted = tuple(_circle(k, centre, rho) for k, (centre, rho) in enumerate(circles))
    grid = _DyadicGrid([float(TOOL.value), *(v for centre, rho in circles for v in (*centre, rho)), *a, *b])
    index = _CandidateIndex(emitted, TOOL, grid)
    return index.near((grid.length(a[0]), grid.length(a[1])), (grid.length(b[0]), grid.length(b[1])), range(len(circles) - 1, -1, -1))


def _boundary_circles(end: tuple[float, float], offset: tuple[float, float]) -> list[tuple[tuple[float, float], float]]:
    """Circles at ``end + offset`` whose reach is the float nearest ``|offset|``, nudged by -3..3 ulps."""
    centre = (end[0] + offset[0], end[1] + offset[1])
    distance = math.hypot(*offset)
    circles = []
    for nudge in range(-3, 4):
        reach = distance
        for _ in range(abs(nudge)):
            reach = math.nextafter(reach, math.inf if nudge > 0 else -math.inf)
        circles.append((centre, reach - 2 * float(TOOL.value)))
    return circles


def test_rs_hypot_witness_is_decided_exactly() -> None:
    a, b = (0.0, 0.0), (0.0, 2.0**-30)
    circles = _boundary_circles(a, RS_OFFSET)
    assert _near(circles, a, b) == _oracle(circles, a, b)


@settings(max_examples=300, deadline=None)
@example(end=(0.0, 0.0), offset=RS_OFFSET, step=(0.0, 2.0**-30))
@given(
    end=st.tuples(st.floats(-60.0, 60.0), st.floats(-60.0, 60.0)),
    offset=st.tuples(st.floats(-50.0, 50.0), st.floats(-50.0, 50.0)).filter(lambda v: math.hypot(*v) > 2.5),
    step=st.tuples(st.floats(-1.0, 1.0), st.floats(-1.0, 1.0)).filter(lambda v: math.hypot(*v) > 1e-6),
)
def test_candidate_filter_equals_the_exact_midpoint_criterion(end: tuple[float, float], offset: tuple[float, float], step: tuple[float, float]) -> None:
    a = end
    b = (end[0] + step[0], end[1] + step[1])
    if a == b:
        return
    length = math.hypot(*offset) / math.hypot(*step)
    along = (step[0] * length, step[1] * length)  # on the piece's extension beyond b: |c - m| = h + |along|
    circles = _boundary_circles(a, offset) + _boundary_circles(b, offset) + _boundary_circles(b, along)
    assert _near(circles, a, b) == _oracle(circles, a, b)
