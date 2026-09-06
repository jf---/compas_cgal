"""Circle deletion preserves material and the actual connector chain."""

import math
from dataclasses import replace
from fractions import Fraction

import pytest

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_motion_coverage import replay_circle_connector_stock
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _source(xs: tuple[float, ...], scale: float = 1.0) -> tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[Figure5BoundaryTransition, ...]]:
    circles = []
    for index, x in enumerate(xs):
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
    transitions = tuple(
        Figure5BoundaryTransition((a.contact_point, b.contact_point), Millimetre(abs(float(b.contact_point.x) - float(a.contact_point.x))), a.contact_point == b.contact_point)
        for a, b in zip(circles, circles[1:])
    )
    return tuple(circles), transitions


def _boundary(scale: float = 1.0) -> tuple[Point2[WorldXY], ...]:
    return tuple(Point2[WorldXY].build(x * scale, y * scale) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))


@pytest.mark.parametrize("scale", (0.125, 1.0, 8.0))
def test_required_circle_is_retained_at_every_scale(scale: float) -> None:
    from benchmarks.held_coverage_thinning import thin_covered_contour_path

    circles, transitions = _source((2, 2.2, 2.4), scale)
    selected, path, _, _ = thin_covered_contour_path(
        circles, ((GuideRunId(0),),) * 3, transitions, _boundary(scale), ToolRadius.build(scale), EngagementCap.build(math.radians(80))
    )
    assert 1 in selected
    after = replay_circle_connector_stock(_boundary(scale), path.circles, path.transitions, ToolRadius.build(scale))
    assert not after.contains(2.2 * scale, 2.995 * scale)


def test_redundant_circles_are_deleted_without_rebuilding_connectors() -> None:
    from benchmarks.held_coverage_thinning import thin_covered_contour_path

    circles, _ = _source((2, 2, 2, 2))
    q = circles[0].contact_point
    detour = Point2[WorldXY].build(5, 0)
    original = (
        Figure5BoundaryTransition((q, detour, q), Millimetre(6), False),
        Figure5BoundaryTransition((q, q), Millimetre(0), True),
        Figure5BoundaryTransition((q, q), Millimetre(0), True),
    )
    selected, path, _, _ = thin_covered_contour_path(circles, ((GuideRunId(0),),) * 4, original, _boundary(), ToolRadius.build(1), EngagementCap.build(math.radians(80)))
    assert selected == (0, 3)
    assert path.transitions[0].samples == (q, detour, q, q, q)
    assert path.transitions[0].boundary_progress == Millimetre(6)
    assert not path.transitions[0].within_fan
    before = replay_circle_connector_stock(_boundary(), circles, original, ToolRadius.build(1))
    after = replay_circle_connector_stock(_boundary(), path.circles, path.transitions, ToolRadius.build(1))
    assert after.exactly_equals(before)


def test_incomplete_dense_source_is_not_claimed_to_be_complete() -> None:
    from benchmarks.held_coverage_thinning import thin_covered_contour_path

    circles, transitions = _source((2, 2, 2))
    _, path, _, _ = thin_covered_contour_path(circles, ((GuideRunId(0),),) * 3, transitions, _boundary(), ToolRadius.build(1), EngagementCap.build(math.radians(80)))
    after = replay_circle_connector_stock(_boundary(), path.circles, path.transitions, ToolRadius.build(1))
    assert after.contains(8, 8)


def test_pipeline_refines_before_coverage_preserving_thinning() -> None:
    from benchmarks.held_coverage_path import build_coverage_preserving_path

    circles, _ = _source((2, 2.2, 2.4))
    dense, path, bounds, stock = build_coverage_preserving_path(
        circles,
        ((GuideRunId(0),),) * 3,
        _boundary(),
        _boundary(),
        ToolRadius.build(1),
        EngagementCap.build(math.radians(80)),
    )
    assert len(dense.circles) >= 3
    assert len(bounds) == len(path.circles) - 1
    assert not stock.contains(2.2, 2.995)


def test_distinct_source_sites_with_equal_contacts_keep_both_endpoints() -> None:
    from benchmarks.held_figure5_boundary_path import _side_lengths
    from benchmarks.held_figure5_boundary_path import _transition

    circles, _ = _source((5, 5))
    # Distinct source parameters can project to the same stored double contact.
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(576460752303423489, 1152921504606846976), 4)
    following = replace(circles[1], boundary_site=site)
    transition = _transition(_boundary(), _side_lengths(_boundary()), circles[0], following)
    assert transition.samples == (circles[0].contact_point, following.contact_point)
