"""Boundary-order selection carries one predecessor and preserves run evidence."""

import math
from fractions import Fraction

import pytest

from benchmarks.held_boundary_order_spacing import InvalidBoundarySpacingInputError
from benchmarks.held_boundary_order_spacing import select_boundary_ordered_sources
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import maximum_predecessor_engagement
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


def test_boundary_spacing_removes_dense_stations_without_losing_endpoints() -> None:
    circles = tuple(_circle(index) for index in range(101))
    selected = select_boundary_ordered_sources(circles, ((GuideRunId(0),),) * len(circles), TOOL, CAP)
    assert selected[0] == 0 and selected[-1] == len(circles) - 1
    assert len(selected) < len(circles) / 2
    assert tuple(sorted(set(selected))) == selected
    for a, b in zip(selected, selected[1:]):
        assert float(maximum_predecessor_engagement(_paper_candidate(circles[a]), _paper_candidate(circles[b]), TOOL)) <= float(CAP.theta)


def test_last_opportunity_retains_a_disconnected_source_family() -> None:
    circles = tuple(_circle(index) for index in range(5))
    sources = ((GuideRunId(0),), (GuideRunId(1),), (GuideRunId(1), GuideRunId(2)), (GuideRunId(0),), (GuideRunId(0),))
    selected = select_boundary_ordered_sources(circles, sources, TOOL, CAP)
    assert 2 in selected
    assert {run for index in selected for run in sources[index]} == {GuideRunId(0), GuideRunId(1), GuideRunId(2)}


def test_wide_gap_reaches_refinement_without_silently_forcing_final_motion() -> None:
    circles = (_circle(0), _circle(600))
    selected = select_boundary_ordered_sources(circles, ((GuideRunId(0),), (GuideRunId(1),)), TOOL, CAP)
    result = refine_figure5_engagement(tuple(circles[index] for index in selected), BOUNDARY, TOOL, CAP)
    assert len(result.circles) > 2
    for a, b in zip(result.circles, result.circles[1:]):
        assert float(maximum_predecessor_engagement(_paper_candidate(a), _paper_candidate(b), TOOL)) <= float(CAP.theta)


def test_missing_source_ledger_fails_loudly() -> None:
    with pytest.raises(InvalidBoundarySpacingInputError):
        select_boundary_ordered_sources((_circle(0), _circle(1)), ((GuideRunId(0),),), TOOL, CAP)
