"""Full-motion coverage rejects residuals independently of inspection grids."""

from fractions import Fraction

import pytest

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _circle(x: float, radius: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(0), 4)
    return Figure5CounterclockwiseCircle(
        GuideRunId(0),
        GuideRunStationOrdinal(0),
        site,
        site,
        Point2[WorldXY].build(x, 0),
        Point2[WorldXY].build(x, -radius),
        GuideRadius.build(radius),
    )


def _square(half_width: float) -> tuple[Point2[WorldXY], ...]:
    return tuple(Point2[WorldXY].build(x, y) for x, y in ((-half_width, -half_width), (half_width, -half_width), (half_width, half_width), (-half_width, half_width)))


@pytest.mark.parametrize("scale", (0.125, 1.0, 8.0))
def test_small_uncut_island_is_rejected_continuously(scale: float) -> None:
    from benchmarks.held_exact_motion_coverage import replay_exact_motion_coverage

    # Inner radius 0.03125*scale: much smaller than a quarter-tool grid.
    residual = replay_exact_motion_coverage(_square(scale), (_circle(0, 1.03125 * scale),), (), ToolRadius.build(scale))
    assert not residual.is_empty()
    assert residual.contains(0, 0)


def test_circle_sweep_can_clear_entire_target_without_invented_preclear() -> None:
    from benchmarks.held_exact_motion_coverage import replay_exact_motion_coverage

    residual = replay_exact_motion_coverage(_square(1), (_circle(0, 0.5),), (), ToolRadius.build(1))
    assert residual.is_empty()


def test_connector_contributes_to_exact_coverage() -> None:
    from benchmarks.held_exact_motion_coverage import replay_exact_motion_coverage

    circles = (_circle(-3, 1), _circle(3, 1))
    transition = Figure5BoundaryTransition((circles[0].contact_point, circles[1].contact_point), Millimetre(6), False)
    residual = replay_exact_motion_coverage(_square(5), circles, (transition,), ToolRadius.build(1))
    assert not residual.contains(0, -1)
    assert residual.contains(0, 1)


def test_identical_incomplete_paths_are_not_absolute_coverage() -> None:
    from benchmarks.held_exact_motion_coverage import replay_exact_motion_coverage

    residual = replay_exact_motion_coverage(_square(1), (_circle(0, 1.5),), (), ToolRadius.build(1))
    assert residual.is_subset_of(residual)
    assert not residual.is_empty()
