"""Whole-orbit leading-rim upper bound from native cleared-contour clearance."""

import math

import pytest

from compas_cgal import _stock_2


def test_single_disk_bound_matches_radial_clearance_equation() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 2, 1)
    angle, exceeded = contour.engagement_bound((1, 0), 2, 2)
    assert angle == pytest.approx(math.acos(-0.25))
    assert exceeded


def test_joint_disk_union_increases_clearance_beyond_either_disk() -> None:
    contour = _stock_2.HeldDiskContour2((-1, 0), 2, 1)
    contour.append((1, 0), 2)
    # Closest union boundary is (0, +/-sqrt(8)), a trimmed-arc endpoint.
    angle, exceeded = contour.engagement_bound((0, 0), 2, 1)
    assert angle == pytest.approx(math.acos(0.75))
    assert not exceeded


def test_clockwise_contact_shortcut_cannot_authorize_unsafe_spacing() -> None:
    contour = _stock_2.HeldDiskContour2((3, -1), 0.1, 1)
    contour.append((0, 0), 2)
    angle, exceeded = contour.engagement_bound((1, 0), 2, 2)
    # CW contact + Eq7 gives ~85.37 degrees. At q=(3,0), however, an
    # untouched connected forward-rim arc spans acos(-1/6) > 99 degrees.
    assert angle >= math.acos(-1 / 6)
    assert exceeded


def test_small_orbit_bound_does_not_assume_the_figure4d_shortcut() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 2, 1)
    angle, exceeded = contour.engagement_bound((2, 0), 0.5, 2)
    assert angle == pytest.approx(math.acos(-0.25))
    assert exceeded


def test_hole_boundary_limits_clearance_inside_another_disk() -> None:
    contour = _stock_2.HeldDiskContour2((0, 2), 0.8, 1)
    for center in ((2, 0), (0, -2), (-2, 0)):
        contour.append(center, 0.8)
    # C=(0,1) is inside the upper disk, 0.8 from the central hole boundary.
    angle, exceeded = contour.engagement_bound((0, 1), 0.5, 2)
    assert angle == pytest.approx(math.acos(-0.61))
    assert exceeded


def test_contained_and_outside_orbits_have_explicit_bounds() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 2, 1)
    assert contour.engagement_bound((0, 0), 2, 2) == (0, False)
    assert contour.engagement_bound((4, 0), 2, 2) == (math.pi, True)
    assert contour.engagement_bound((3, 0), 2, 2) == (math.pi, True)


def test_cap_comparison_is_exact_at_the_surrogate_boundary() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 1, 1)
    assert not contour.engagement_bound((1, 0), 1, 3)[1]
    assert contour.engagement_bound((1, 0), 1, math.nextafter(3, 0))[1]


def test_bound_is_available_when_latest_disk_is_fully_covered() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 4, 1)
    contour.append((0, 0), 1)
    assert contour.engagement_bound((1, 0), 2, 2) == (0, False)
