"""Held predecessor decisions consume exact native circle proposals."""

import math

import numpy as np
import pytest

from compas_cgal import _circle_geometry_2 as native
from compas_cgal import _coverage_2 as coverage


def _circle(cx: float, cy: float, radius: float, tool: float = 1.0) -> native.BoundaryNormalCircleProposal2:
    """A rectangle's bottom normal constructs the requested circle natively."""
    clearance = 2 * radius + tool
    bottom = cy - radius - tool
    half_width = 2 * clearance
    owner = native.BoundaryNormalCircle2(
        [
            (cx - half_width, bottom),
            (cx + half_width, bottom),
            (cx + half_width, bottom + 2 * clearance),
            (cx - half_width, bottom + 2 * clearance),
        ]
    )
    return owner.query(0, 0.5, tool)


def _stationary() -> native.BoundaryNormalCircleProposal2:
    polygon = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
    cycle = coverage.build_center_boundary_cycle(np.array([(x, y, 0.0) for x, y in polygon]), [], 1.0)
    # Select the known analytic corner while passing its native contact owner.
    contact = next(primitive.start for primitive in cycle.primitives if primitive.start_mm == (1.0, 1.0))
    return coverage.boundary_circle_at_contact(native.BoundaryNormalCircle2(polygon), contact, 1.0)


def test_figure4b_angle_and_exact_cap_boundary() -> None:
    previous, current = _circle(0, 0, 2), _circle(1, 0, 2)
    angle, exceeded = native.boundary_circle_engagement(previous, current, 2.5)
    assert angle == pytest.approx(math.acos(-0.25))
    assert not exceeded
    assert native.boundary_circle_engagement(previous, current, math.nextafter(2.5, 0.0))[1]


def test_figure4d_uses_corrected_outer_disk_intersection() -> None:
    previous, current = _circle(0, 0, 2), _circle(2, 0, 0.5)
    angle, exceeded = native.boundary_circle_engagement(previous, current, 2.0)
    # Q=31/6, H=-17/6; cos(theta)=2*(289/744)-1=-83/372.
    assert angle == pytest.approx(math.acos(-83 / 372))
    assert exceeded


def test_contained_and_identical_successors_remove_no_new_material() -> None:
    previous = _circle(0, 0, 2)
    assert native.boundary_circle_engagement(previous, _circle(0.5, 0, 1), 0.0) == (0.0, False)
    assert native.boundary_circle_engagement(previous, previous, 0.0) == (0.0, False)


def test_concentric_growth_has_nonzero_engagement() -> None:
    angle, exceeded = native.boundary_circle_engagement(_circle(0, 0, 1), _circle(0, 0, 2), 2.0)
    assert angle == pytest.approx(math.acos(-0.25))
    assert exceeded


def test_eq4_gap_is_rejected_but_tangency_is_full_slot() -> None:
    previous = _circle(0, 0, 1)
    with pytest.raises(native.BoundaryCircleSpacingError):
        native.boundary_circle_engagement(previous, _circle(3, 0, 1), 4.0)
    angle, exceeded = native.boundary_circle_engagement(previous, _circle(2, 0, 1), 4.0)
    assert angle == pytest.approx(math.pi)
    assert not exceeded
    assert native.boundary_circle_engagement(previous, _circle(2, 0, 1), math.nextafter(4.0, 0.0))[1]


def test_stationary_successor_requires_already_cleared_disk() -> None:
    stationary = _stationary()
    assert stationary.is_stationary
    assert native.boundary_circle_engagement(_circle(2, 2, 2), stationary, 0.0) == (0.0, False)
    with pytest.raises(native.UncoveredStationaryCircleError):
        native.boundary_circle_engagement(_circle(20, 20, 1), stationary, 4.0)


def test_stationary_predecessor_supplies_its_declared_cleared_disk() -> None:
    stationary = _stationary()
    current = _circle(1, 1, 1)
    angle, exceeded = native.boundary_circle_engagement(stationary, current, 3.0)
    assert angle == pytest.approx(math.acos(-0.5))
    assert not exceeded


@pytest.mark.parametrize("ratio", [-1.0, 5.0, math.nan, math.inf])
def test_invalid_surrogate_has_named_error(ratio: float) -> None:
    circle = _circle(0, 0, 1)
    with pytest.raises(native.InvalidBoundaryEngagementCapError):
        native.boundary_circle_engagement(circle, circle, ratio)


def test_tool_radius_mismatch_is_rejected_before_containment() -> None:
    with pytest.raises(native.BoundaryCircleToolMismatchError):
        native.boundary_circle_engagement(_circle(0, 0, 2), _circle(0, 0, 1, 0.5), 4.0)


@pytest.mark.parametrize("previous_radius,distance", [(0.5, 1.75), (1.0, 1.5)])
def test_unresolved_closest_contact_reports_conservative_pi(previous_radius: float, distance: float) -> None:
    previous, current = _circle(0, 0, previous_radius), _circle(distance, 0, 0.25)
    angle, exceeded = native.boundary_circle_engagement(previous, current, 2.0)
    assert angle == pytest.approx(math.pi)
    assert exceeded


@pytest.mark.parametrize("scale", [1 / 1024, 1024.0])
def test_binary_scale_preserves_exact_threshold(scale: float) -> None:
    previous = _circle(0, 0, 2 * scale, scale)
    current = _circle(scale, 0, 2 * scale, scale)
    assert not native.boundary_circle_engagement(previous, current, 2.5)[1]
    assert native.boundary_circle_engagement(previous, current, math.nextafter(2.5, 0.0))[1]
