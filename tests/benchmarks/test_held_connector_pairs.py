"""Connector pairing of open and closed (seam-connected) Held circle paths."""

import math
from fractions import Fraction

import pytest

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_motion_coverage import InvalidCoverageMotionError
from benchmarks.held_motion_coverage import connector_pairs
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


def _circles(xs: tuple[float, ...]) -> tuple[Figure5CounterclockwiseCircle, ...]:
    circles = []
    for index, x in enumerate(xs):
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x / 10), 4)
        circles.append(
            Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, Point2[WorldXY].build(x, 1), Point2[WorldXY].build(x, 0), GuideRadius.build(1))
        )
    return tuple(circles)


def _connect(a: Figure5CounterclockwiseCircle, b: Figure5CounterclockwiseCircle) -> Figure5BoundaryTransition:
    return Figure5BoundaryTransition((a.contact_point, b.contact_point), Millimetre(abs(float(b.contact_point.x) - float(a.contact_point.x))), False)


def test_closed_path_pairs_the_seam_last_to_first() -> None:
    circles = _circles((2, 4, 6))
    transitions = (_connect(circles[0], circles[1]), _connect(circles[1], circles[2]), _connect(circles[2], circles[0]))
    pairs = connector_pairs(circles, transitions)
    assert [(a, b) for a, b, _ in pairs] == [(circles[0], circles[1]), (circles[1], circles[2]), (circles[2], circles[0])]


def test_open_path_pairs_successors_only() -> None:
    circles = _circles((2, 4, 6))
    transitions = (_connect(circles[0], circles[1]), _connect(circles[1], circles[2]))
    assert len(connector_pairs(circles, transitions)) == 2


@pytest.mark.parametrize("count", (0, 1, 4))
def test_connector_count_must_be_open_or_closed(count: int) -> None:
    circles = _circles((2, 4, 6))
    transitions = tuple(_connect(circles[0], circles[0]) for _ in range(count))
    with pytest.raises(InvalidCoverageMotionError, match="connectors"):
        connector_pairs(circles, transitions)


def test_seam_connector_must_end_at_the_first_contact() -> None:
    circles = _circles((2, 4, 6))
    transitions = (_connect(circles[0], circles[1]), _connect(circles[1], circles[2]), _connect(circles[2], circles[1]))
    with pytest.raises(InvalidCoverageMotionError, match="endpoints"):
        connector_pairs(circles, transitions)


def test_thinning_keeps_the_seam_connector() -> None:
    from benchmarks.held_coverage_thinning import thin_covered_contour_path

    circles = _circles((2, 2.2, 2.4))
    seam = _connect(circles[2], circles[0])
    transitions = (_connect(circles[0], circles[1]), _connect(circles[1], circles[2]), seam)
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
    selected, path, bounds, _ = thin_covered_contour_path(circles, ((GuideRunId(0),),) * 3, transitions, boundary, ToolRadius.build(1), EngagementCap.build(math.radians(80)))
    assert len(path.transitions) == len(path.circles) == len(selected)
    assert path.transitions[-1] is seam
    assert len(bounds) == len(path.circles) - 1
