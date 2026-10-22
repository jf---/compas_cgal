"""Tier-1b neighbour prefilter: a provable superset of every circle the candidate filter can return.

``_ConnectorGate`` puts only the emitted circles ``_tier_1b_neighbours`` keeps
on the exact grid; ``_CandidateIndex.near`` then asks, per dyadic sub-piece,
binary64 questions about them. The prefilter must keep every circle that any
sub-piece's binary64 test returns, with no slack constant: the fixtures put
circles exactly on the reach boundary, on the piece's extension line, where
the un-slacked whole-piece binary64 test drops a circle a sub-piece needs.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import benchmarks.held_figure5_engagement_refinement as refinement
from benchmarks.held_connector_engagement import MAX_BISECTION_DEPTH
from benchmarks.held_connector_engagement import _at
from benchmarks.held_connector_engagement import _CandidateIndex
from benchmarks.held_connector_engagement import _DyadicGrid
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import _ConnectorGate
from benchmarks.held_figure5_engagement_refinement import _tier_1b_neighbours
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
FULL = 1 << MAX_BISECTION_DEPTH
# Found by build/neighbour-prefilter/search.py: the circle sits on the
# extension of the piece beyond b at exactly reach (one ulp nudged); the
# sub-piece [1 - 2**-9, 1] (depth 9, sharing b) returns it, the un-slacked
# whole-piece binary64 test does not.
WITNESS_A = (15.982389175107969, -7.370062128614549)
WITNESS_B = (15.932203634019697, -8.005242413715974)
WITNESS_CENTRE = (15.660104537837267, -11.449102492558177)
WITNESS_RHO = 1.4545926186435227
WITNESS_DEPTH = 9


def _p(x: float, y: float) -> Point2[WorldXY]:
    point: Point2[WorldXY] = Point2[WorldXY].build(x, y)
    return point


def _circle(index: int, centre: tuple[float, float], rho: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(index, 10_000), 4)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, _p(*centre), _p(*centre), GuideRadius.build(rho))


def _reach(rho: float) -> float:
    return rho + 2 * float(TOOL.value)


def _sound(circles: list[tuple[tuple[float, float], float]], samples: list[tuple[float, float]]) -> set[int]:
    centres = np.array([centre for centre, _ in circles], dtype=np.float64).reshape(-1, 2)
    reach = np.array([_reach(rho) for _, rho in circles], dtype=np.float64)
    return {int(k) for k in np.nonzero(_tier_1b_neighbours(centres, reach, samples))[0]}


def _naive(circles: list[tuple[tuple[float, float], float]], samples: list[tuple[float, float]]) -> set[int]:
    """The pre-fix whole-piece binary64 test with the slack constant at zero."""
    centres = np.array([centre for centre, _ in circles], dtype=np.float64).reshape(-1, 2)
    reach = np.array([_reach(rho) for _, rho in circles], dtype=np.float64)
    near = np.zeros(len(circles), dtype=bool)
    for (ax, ay), (bx, by) in zip(samples, samples[1:]):
        half = math.hypot(bx - ax, by - ay) / 2
        near |= np.hypot(centres[:, 0] - (ax + bx) / 2, centres[:, 1] - (ay + by) / 2) <= reach + half
    return {int(k) for k in np.nonzero(near)[0]}


def _sub_pieces(extra: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Every dyadic sub-piece sharing an end with its piece, at every depth, plus ``extra``."""
    ends = [(FULL - (FULL >> depth), FULL) for depth in range(MAX_BISECTION_DEPTH + 1)]
    starts = [(0, FULL >> depth) for depth in range(MAX_BISECTION_DEPTH + 1)]
    return ends + starts + extra


def _candidates(circles: list[tuple[tuple[float, float], float]], samples: list[tuple[float, float]], sub_pieces: list[tuple[int, int]]) -> set[int]:
    """Every circle ``_CandidateIndex.near`` returns for some sub-piece, on the gate's own grid."""
    emitted = tuple(_circle(k, centre, rho) for k, (centre, rho) in enumerate(circles))
    grid = _DyadicGrid([float(TOOL.value), *(v for centre, rho in circles for v in (*centre, rho)), *(v for point in samples for v in point)])
    index = _CandidateIndex(emitted, TOOL, grid)
    prior = range(len(circles) - 1, -1, -1)
    found: set[int] = set()
    vertices = [(grid.length(x), grid.length(y)) for x, y in samples]
    for a, b in zip(vertices, vertices[1:]):
        if a == b:
            continue
        for start, end in sub_pieces:
            found.update(index.near(_at(a, b, start), _at(a, b, end), prior))
    return found


def test_the_slack_constant_is_gone() -> None:
    assert not hasattr(refinement, "TIER_1B_PREFILTER_SLACK_TOOL_RADII")


def test_boundary_witness_naive_drops_sound_keeps() -> None:
    circles = [(WITNESS_CENTRE, WITNESS_RHO)]
    samples = [WITNESS_A, WITNESS_B]
    needed = _candidates(circles, samples, [(FULL - (FULL >> WITNESS_DEPTH), FULL)])
    # The exact filter rule |c - m| <= rho + 2 r + h drops the witness on that
    # sub-piece: its squared excess (|c - m|² - R² - h²)² - 4 R² h² is 6.5e-17 > 0
    # (build/oc/witness_check.py). The old binary64 near kept it only by rounding.
    assert needed <= {0}
    assert _naive(circles, samples) == set()
    assert _sound(circles, samples) == {0}


def test_connector_gate_nearby_keeps_the_boundary_witness() -> None:
    boundary = (_p(-100.0, -100.0), _p(100.0, -100.0), _p(100.0, 100.0), _p(-100.0, 100.0))
    gate = _ConnectorGate(boundary, TOOL, EngagementCap.build(math.radians(80.0)))
    witness = _circle(0, WITNESS_CENTRE, WITNESS_RHO)
    last = _circle(1, (90.0, 90.0), 1.0)
    transition = Figure5BoundaryTransition((_p(*WITNESS_A), _p(*WITNESS_B)), Millimetre(math.dist(WITNESS_A, WITNESS_B)), False)
    assert gate._nearby((witness, last), transition) == (witness, last)


@st.composite
def _boundary_fixture(draw: st.DrawFn) -> tuple[list[tuple[tuple[float, float], float]], list[tuple[float, float]], list[tuple[int, int]]]:
    """2-4 samples; one circle per piece on the reach boundary beyond either end, nudged by up to 3 ulps."""
    coordinate = st.floats(-60.0, 60.0, allow_nan=False)
    first = (draw(coordinate), draw(coordinate))
    samples = [first]
    for _ in range(draw(st.integers(1, 3))):
        angle = draw(st.floats(0.0, 2 * math.pi))
        length = draw(st.floats(1e-3, 8.0))
        x, y = samples[-1]
        samples.append((x + length * math.cos(angle), y + length * math.sin(angle)))
    circles: list[tuple[tuple[float, float], float]] = []
    for a, b in zip(samples, samples[1:]):
        rho = draw(st.floats(0.25, 6.0))
        start, end = (a, b) if draw(st.booleans()) else (b, a)
        length = math.dist(start, end)
        if length == 0.0:
            continue
        ux, uy = (end[0] - start[0]) / length, (end[1] - start[1]) / length
        reach = _reach(rho)
        centre = [end[0] + reach * ux, end[1] + reach * uy]
        for axis in (0, 1):
            for _ in range(draw(st.integers(0, 3))):
                centre[axis] = math.nextafter(centre[axis], draw(st.sampled_from([math.inf, -math.inf])))
        circles.append(((centre[0], centre[1]), rho))
    if not circles:
        circles.append((first, 1.0))
    low = draw(st.integers(0, FULL - 1))
    high = draw(st.integers(low + 1, FULL))
    return circles, samples, [(low, high)]


@settings(max_examples=300, deadline=None)
@given(_boundary_fixture())
def test_prefilter_is_a_superset_of_every_sub_piece_candidate(fixture: tuple[list[tuple[tuple[float, float], float]], list[tuple[float, float]], list[tuple[int, int]]]) -> None:
    circles, samples, extra = fixture
    assert _candidates(circles, samples, _sub_pieces(extra)) <= _sound(circles, samples)


# RS's input: np.hypot is one ulp above math.hypot here, whose square already
# exceeds the exact x² + y²; the filter must not depend on it.
RS_OFFSET = (3.6163379554499997, 48.30545762211739)


def test_near_keeps_exactly_the_circles_within_reach_plus_half_length_on_rs_input() -> None:
    """``near`` keeps exactly the circles with ``|c - m| <= rho + 2 r + h``, decided in rationals by squaring."""
    # a piece far shorter than an ulp of the reach, so h does not mask the boundary
    a, b = (0.0, 0.0), (0.0, 2.0**-60)
    centre = RS_OFFSET
    circles: list[tuple[tuple[float, float], float]] = []
    for nudge in range(-3, 4):
        reach = math.hypot(*RS_OFFSET)
        for _ in range(abs(nudge)):
            reach = math.nextafter(reach, math.inf if nudge > 0 else -math.inf)
        circles.append((centre, reach - 2 * float(TOOL.value)))
    ax, ay, bx, by = (Fraction(v) for v in (*a, *b))
    mx, my = (ax + bx) / 2, (ay + by) / 2
    half_sq = ((bx - ax) ** 2 + (by - ay) ** 2) / 4
    expected = set()
    for k, ((cx, cy), rho) in enumerate(circles):
        reach_exact = Fraction(rho) + 2 * Fraction(float(TOOL.value))
        excess = (Fraction(cx) - mx) ** 2 + (Fraction(cy) - my) ** 2 - reach_exact**2 - half_sq
        if excess <= 0 or excess**2 <= 4 * reach_exact**2 * half_sq:
            expected.add(k)
    assert 0 < len(expected) < len(circles), "the nudges straddle the exact boundary"
    assert _candidates(circles, [a, b], [(0, FULL)]) == expected
