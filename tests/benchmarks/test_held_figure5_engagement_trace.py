"""Contract tests of the refinement trace that the splice's exact early stop reads.

``EngagementTrace.last_influenced`` must flag exactly the decisions a circle
could change: a tier-1b decision iff the gate's own filter keeps the circle
for one of its pieces, a convex-corner floor hypothesis iff the circle is
within rho + 2 r of the vertex (exact, boundary included). Design:
docs/held_coverage_preserving_selection.md, step 3.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_engagement_refinement import circle_on_side
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_engagement_trace import EngagementFootprint
from benchmarks.held_figure5_engagement_trace import EngagementRefinementKey
from benchmarks.held_figure5_engagement_trace import EngagementTrace
from benchmarks.held_figure5_engagement_trace import PieceBounds
from benchmarks.held_figure5_engagement_trace import _tier_1b_neighbours
from benchmarks.held_figure5_engagement_trace import tier_1b_reach
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80))
# The tooth of test_held_connector_repair.py: vertex 2 = (0.5, 0.0625) is convex with an
# exterior angle of about 92.8 degrees, over the 80-degree cap, so it is a floor vertex.
TOOTH = tuple(Point2[WorldXY].build(x, y) for x, y in ((-3.0, 0.0), (0.25, 0.0), (0.5, 0.0625), (-0.5, 3.375), (-3.0, 3.375)))
KEY = EngagementRefinementKey(TOOTH, TOOL, CAP, 16, False, False)
# Below the binary64 spacing of 3.5 (2**-51 * 2): the next double above the exact boundary centre.
BEYOND = math.nextafter(3.5, math.inf)


def _circle(x: float, y: float, radius: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(1, 2), len(TOOTH))
    contact = Point2[WorldXY].build(x, y - radius)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(0), site, site, Point2[WorldXY].build(x, y), contact, GuideRadius.build(radius))


def test_a_floor_footprint_is_influenced_exactly_within_rho_plus_two_r_of_its_vertex() -> None:
    """|c - V| = 3 = rho + 2 r exactly: influenced (the floor hypothesis reads it); one ulp further: not."""
    trace = EngagementTrace(KEY, (EngagementFootprint(5, vertex=2),))
    assert trace.last_influenced((_circle(3.5, 0.0625, 1.0),), 0) == 5
    assert trace.last_influenced((_circle(BEYOND, 0.0625, 1.0),), 0) == -1
    assert trace.last_influenced((_circle(3.5, 0.0625, 1.0),), 6) == -1, "a decision before `since` is not a later decision"


def test_the_last_influenced_position_is_the_latest_of_both_kinds() -> None:
    near_piece = PieceBounds.of([(0.0, 0.0), (0.1, 0.0)])
    trace = EngagementTrace(
        KEY, (EngagementFootprint(3, pieces=near_piece), EngagementFootprint(7, vertex=2), EngagementFootprint(9, pieces=PieceBounds.of([(-90.0, 50.0), (-89.0, 50.0)])))
    )
    assert trace.last_influenced((_circle(0.5, 1.0, 1.0),), 0) == 7
    assert trace.last_influenced((_circle(-2.0, 1.0, 0.5),), 0) == 3, "beyond rho + 2 r of the vertex, kept by the near piece only"


@settings(max_examples=300, deadline=None)
@given(
    st.floats(-20, 20, allow_nan=False),
    st.floats(-20, 20, allow_nan=False),
    st.floats(0.125, 4, allow_nan=False),
    st.lists(st.tuples(st.floats(-20, 20, allow_nan=False), st.floats(-20, 20, allow_nan=False)), min_size=2, max_size=5),
)
def test_tier_1b_influence_is_the_gates_own_filter(x: float, y: float, radius: float, samples: list[tuple[float, float]]) -> None:
    """Replay identity: a circle influences a tier-1b footprint iff ``_tier_1b_neighbours`` keeps it for that connector."""
    circle = _circle(x, y, radius)
    trace = EngagementTrace(KEY, (EngagementFootprint(4, pieces=PieceBounds.of(samples)),))
    kept = _tier_1b_neighbours(np.array([[x, y]]), np.array([tier_1b_reach(circle, TOOL)]), samples)
    assert (trace.last_influenced((circle,), 0) == 4) == bool(kept[0])


def _on(side: int, parameter: Fraction, radius: float, ordinal: int) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, len(TOOTH))
    contact = _canonical_offset_point(TOOTH, site)
    raw = Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(ordinal), site, site, contact, contact, GuideRadius.build(radius))
    return circle_on_side(raw, site, contact, GuideRadius.build(radius), TOOTH)


def test_a_real_refinement_records_floor_footprints_and_its_trace_compares_by_value() -> None:
    """The tooth: the R = 1.8 circle tangent at (0, 0) leaves vertex 2's wedge, so the arrival side 1 -> side 2 is gated, not refused."""
    sources = (_on(0, Fraction(12, 13), 1.8, 0), _on(1, Fraction(1, 2), 0.25, 1), _on(2, Fraction(7, 8), 0.5, 2))
    first = refine_figure5_engagement(sources, TOOTH, TOOL, CAP)
    second = refine_figure5_engagement(sources, TOOTH, TOOL, CAP)
    floors = [footprint for footprint in first.trace.footprints if footprint.vertex is not None]
    assert floors and all(footprint.vertex == 2 for footprint in floors)
    assert first.trace.last_influenced((first.circles[0],), 0) >= floors[-1].position
    assert first.trace == second.trace and hash(first.trace) == hash(second.trace)
    pieces = PieceBounds.of([(0.0, 0.0), (1.0, 0.5)])
    assert pieces == PieceBounds.of([(0.0, 0.0), (1.0, 0.5)]) and pieces != PieceBounds.of([(0.0, 0.0), (1.0, 0.25)])
    with pytest.raises(ValueError):
        pieces.half[0] = 0.0
