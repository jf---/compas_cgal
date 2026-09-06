"""Coverage comparisons include the material cut by boundary connectors."""

from fractions import Fraction

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_stock_replay import replay_circle_stock
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def test_longer_connector_clears_material_between_retained_circles() -> None:
    from benchmarks.held_motion_coverage import replay_circle_connector_stock

    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ((-10, -10), (10, -10), (10, 10), (-10, 10)))
    circles = []
    for index, x in enumerate((0, 6)):
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(index), 4)
        circles.append(
            Figure5CounterclockwiseCircle(
                GuideRunId(0),
                GuideRunStationOrdinal(index),
                site,
                site,
                Point2[WorldXY].build(x, 0),
                Point2[WorldXY].build(x, -1),
                GuideRadius.build(1),
            )
        )
    transition = Figure5BoundaryTransition((circles[0].contact_point, circles[1].contact_point), Millimetre(6), False)
    tool = ToolRadius.build(1)
    circle_only = replay_circle_stock(boundary, tuple(circles), tool)
    complete = replay_circle_connector_stock(boundary, tuple(circles), (transition,), tool)
    assert circle_only.contains(3, -1)
    assert not complete.contains(3, -1)
    assert complete.contains(3, 1)
