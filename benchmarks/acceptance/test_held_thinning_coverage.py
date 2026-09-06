"""Regression gate: engagement-safe thinning must retain required swept stock."""

import math
from fractions import Fraction

import pytest

from benchmarks.held_contour_bound_spacing import select_contour_bound_path
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_stock_replay import replay_circle_stock
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


@pytest.mark.parametrize("scale", (0.125, 1.0, 8.0))
def test_thinning_preserves_material_reached_within_one_source_run(scale: float) -> None:
    boundary = tuple(Point2[WorldXY].build(x * scale, y * scale) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
    circles = []
    for index, x in enumerate((2.0, 2.2, 2.4)):
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x / 10), 4)
        circles.append(
            Figure5CounterclockwiseCircle(
                GuideRunId(0),
                GuideRunStationOrdinal(index),
                site,
                site,
                Point2[WorldXY].build(x * scale, scale),
                Point2[WorldXY].build(x * scale, 0),
                GuideRadius.build(scale),
            )
        )
    tool = ToolRadius.build(scale)
    _, result, _ = select_contour_bound_path(tuple(circles), ((GuideRunId(0),),) * 3, boundary, tool, EngagementCap.build(math.radians(80)))
    # The middle circle reaches this interior point. Both retained neighbors
    # miss it; bottom-boundary connectors also cannot reach its y coordinate.
    witness = (2.2 * scale, 2.995 * scale)
    source = replay_circle_stock(boundary, tuple(circles), tool)
    thinned = replay_circle_stock(boundary, result.circles, tool)
    assert not source.contains(*witness)
    assert not thinned.contains(*witness), "Engagement-safe thinning discarded required material coverage."
