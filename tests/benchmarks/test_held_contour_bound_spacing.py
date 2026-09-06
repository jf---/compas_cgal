"""Ordered selection must append every accepted circle to the native contour."""

import math
from fractions import Fraction

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal import _stock_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
TOOL = ToolRadius.build(1)
CAP = EngagementCap.build(math.radians(80))


def _circle(index: int) -> Figure5CounterclockwiseCircle:
    x = 2 + index / 100
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x / 10), 4)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, Point2[WorldXY].build(x, 1), Point2[WorldXY].build(x, 0), GuideRadius.build(1))


def _verify_replay(circles: tuple[Figure5CounterclockwiseCircle, ...]) -> None:
    first = circles[0]
    contour = _stock_2.HeldDiskContour2((float(first.center.x), float(first.center.y)), float(first.radius.value), float(TOOL.value))
    for circle in circles[1:]:
        center = float(circle.center.x), float(circle.center.y)
        assert not contour.engagement_bound(center, float(circle.radius.value), CAP.chord_ratio)[1]
        contour.append(center, float(circle.radius.value))


def test_thinning_and_refinement_replay_under_the_same_evolving_contour() -> None:
    from benchmarks.held_contour_bound_spacing import select_contour_bound_path

    circles = tuple(_circle(index) for index in range(101))
    selected, refined, bounds = select_contour_bound_path(circles, ((GuideRunId(0),),) * len(circles), BOUNDARY, TOOL, CAP)
    assert selected[0] == 0 and selected[-1] == len(circles) - 1
    assert len(refined.circles) < len(circles) / 2
    assert len(bounds) == len(refined.circles) - 1
    assert refined.original_indices[0] == 0
    assert refined.original_indices[-1] == len(refined.circles) - 1
    _verify_replay(refined.circles)


def test_wide_gap_is_subdivided_before_any_over_bound_motion_is_emitted() -> None:
    from benchmarks.held_contour_bound_spacing import select_contour_bound_path

    circles = (_circle(0), _circle(600))
    selected, refined, bounds = select_contour_bound_path(circles, ((GuideRunId(0),), (GuideRunId(1),)), BOUNDARY, TOOL, CAP)
    assert selected == (0, 1)
    assert len(refined.circles) > 2
    assert len(bounds) == len(refined.circles) - 1
    _verify_replay(refined.circles)


def test_last_representative_of_each_source_run_is_retained() -> None:
    from benchmarks.held_contour_bound_spacing import select_contour_bound_path

    circles = tuple(_circle(index) for index in range(5))
    sources = ((GuideRunId(0),), (GuideRunId(1),), (GuideRunId(1), GuideRunId(2)), (GuideRunId(0),), (GuideRunId(0),))
    selected, refined, _ = select_contour_bound_path(circles, sources, BOUNDARY, TOOL, CAP)
    assert 2 in selected
    assert {run for index in selected for run in sources[index]} == {GuideRunId(0), GuideRunId(1), GuideRunId(2)}
    assert len(refined.original_indices) == len(selected)


def test_required_dense_sources_are_refined_without_being_thinned() -> None:
    from benchmarks.held_contour_bound_spacing import select_contour_bound_path

    circles = tuple(_circle(index) for index in range(5))
    selected, refined, _ = select_contour_bound_path(
        circles,
        ((GuideRunId(0),),) * len(circles),
        BOUNDARY,
        TOOL,
        CAP,
        required_sources=frozenset(range(len(circles))),
    )
    assert selected == tuple(range(len(circles)))
    assert len(refined.original_indices) == len(circles)
    _verify_replay(refined.circles)
