"""Consumer witnesses for the native exact decisions used by Held placement."""

import math

import pytest

from compas_cgal import _circle_geometry_2 as geometry


def test_disk_containment_includes_tangency_but_not_a_positive_gap() -> None:
    assert geometry.disk_contains_disk((0, 0), 2, (1, 0), 1)
    assert not geometry.disk_contains_disk((0, 0), 2, (math.nextafter(1, math.inf), 0), 1)
    assert not geometry.disk_contains_disk((0, 0), 1, (0, 0), 2)


def test_near_tangent_intersection_has_positive_height() -> None:
    _, height_squared = geometry.swept_disk_intersection(
        (16.847395307512397, 28.211893255992088),
        5.758981134680233,
        (16.85127854510552, 28.22374119808311),
        5.746513044770512,
        1,
    )
    assert height_squared > 0


def test_intersection_reports_known_geometry_and_rejects_disjoint_disks() -> None:
    assert geometry.swept_disk_intersection((0, 0), 1, (2, 0), 1, 1) == pytest.approx((-1, 3))
    with pytest.raises(geometry.NoCircleIntersectionError):
        geometry.swept_disk_intersection((0, 0), 1, (5, 0), 1, 1)
    with pytest.raises(geometry.NoCircleIntersectionError):
        geometry.swept_disk_intersection((0, 0), 1, (0, 0), 1, 1)


def test_orientation_preserves_small_nonzero_turn() -> None:
    assert geometry.orientation((0, 0), (1, 1), (2, math.nextafter(2, math.inf))) == 1
    assert geometry.orientation((0, 0), (1, 1), (2, 2)) == 0
    assert geometry.orientation((0, 0), (2, 2), (1, math.nextafter(1, -math.inf))) == -1


def test_invalid_native_input_fails_loudly() -> None:
    with pytest.raises(geometry.InvalidCircleGeometryError):
        geometry.disk_contains_disk((math.nan, 0), 2, (0, 0), 1)
    with pytest.raises(geometry.InvalidCircleGeometryError):
        geometry.disk_contains_disk((0, 0), -1, (0, 0), 1)


def test_offset_projection_uses_geometry_and_canonicalizes_shared_vertex() -> None:
    boundary = [(1.0, 1.0), (9.0, 1.0), (9.0, 9.0), (1.0, 9.0)]
    side, parameter = geometry.project_boundary_contact(boundary, (9, 5), 0.01)
    assert (side, parameter) == (1, 0.5)
    assert geometry.project_boundary_contact(boundary, (9, 9), 0.01) == (2, 0)
    with pytest.raises(geometry.AmbiguousBoundaryProjectionError):
        geometry.project_boundary_contact(boundary, (5, 5), 5)
    with pytest.raises(geometry.BoundaryProjectionDistanceError):
        geometry.project_boundary_contact(boundary, (10, 5), 0.01)


def test_corrected_engagement_uses_exact_squared_chord_geometry() -> None:
    assert geometry.corrected_engagement_cosine_squared((0, 0), 2, (1, 0), 2, 1) == pytest.approx(121 / 156)
    with pytest.raises(geometry.NoCircleIntersectionError):
        geometry.corrected_engagement_cosine_squared((0, 0), 1, (5, 0), 1, 1)
