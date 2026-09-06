"""Held Section 3.1: clockwise contact correction on the union of outer disks."""

import math

import pytest

from compas_cgal import _stock_2


def test_uncovered_contact_is_unchanged() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 2, 1)
    point, moved = contour.contact_toward((1, 1))
    assert point == pytest.approx((3 / math.sqrt(2), 3 / math.sqrt(2)))
    assert not moved


@pytest.mark.parametrize("direction", [(1, 0), (1, 1), (1, -1)])
def test_covered_contact_moves_clockwise_to_exposed_intersection(direction: tuple[float, float]) -> None:
    contour = _stock_2.HeldDiskContour2((2, 0), 2, 1)
    contour.append((0, 0), 2)
    point, moved = contour.contact_toward(direction)
    assert moved
    assert point == pytest.approx((1, -math.sqrt(8)))


def test_clockwise_scan_passes_intersection_hidden_by_an_older_disk() -> None:
    contour = _stock_2.HeldDiskContour2((2, 0), 2, 1)
    contour.append((0, -2), 2)
    contour.append((0, 0), 2)
    point, moved = contour.contact_toward((1, 0))
    assert moved
    assert point == pytest.approx((-math.sqrt(8), -1))


def test_tangent_contact_remains_on_exposed_boundary() -> None:
    contour = _stock_2.HeldDiskContour2((6, 0), 2, 1)
    contour.append((0, 0), 2)
    assert contour.contact_toward((1, 0)) == ((3, 0), False)


def test_fully_covered_predecessor_has_no_exposed_contact() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 4, 1)
    contour.append((0, 0), 2)
    with pytest.raises(_stock_2.NoExposedPredecessorArcError):
        contour.contact_toward((1, 0))


def test_repeated_equal_disk_preserves_exposed_contact() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 2, 1)
    contour.append((0, 0), 2)
    assert contour.contact_toward((1, 0)) == ((3, 0), False)


def test_outer_radius_is_added_inside_cgal() -> None:
    contour = _stock_2.HeldDiskContour2((-0.3, 0), 0.1, 0.2)
    point, moved = contour.contact_toward((0, 0))
    expected = 2.0**-55
    # Approximate reporting allows four ulps for two conversions and division.
    # Both zero and the doubled endpoint from a rounded radius sum fail this.
    assert point[0] == pytest.approx(expected, rel=0, abs=4 * math.ulp(expected))
    assert point[1] == 0
    assert not moved


def test_invalid_append_leaves_predecessor_and_union_unchanged() -> None:
    contour = _stock_2.HeldDiskContour2((0, 0), 2, 1)
    with pytest.raises(_stock_2.InvalidHeldContourInputError):
        contour.append((0, 0), -1)
    with pytest.raises(_stock_2.InvalidHeldContourInputError):
        contour.append((math.nan, 0), 2)
    assert contour.contact_toward((1, 0)) == ((3, 0), False)
    with pytest.raises(_stock_2.UndefinedPredecessorDirectionError):
        contour.contact_toward((0, 0))


@pytest.mark.parametrize("scale", [1e155, 1e-160])
def test_coordinate_reporting_does_not_convert_squared_scale_to_double(scale: float) -> None:
    contour = _stock_2.HeldDiskContour2((-scale, 0), scale / 2, scale / 2)
    point, moved = contour.contact_toward((0, scale))
    assert all(math.isfinite(coordinate) for coordinate in point)
    # Dimensionless comparison keeps the small-scale witness above pytest's
    # default absolute tolerance. This is approximate coordinate reporting.
    assert (point[0] / scale, point[1] / scale) == pytest.approx((-1 + 1 / math.sqrt(2), 1 / math.sqrt(2)))
    assert not moved
