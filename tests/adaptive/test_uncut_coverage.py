"""Continuous coverage of explicitly uncut design and reachable regions."""

import numpy as np
import pytest

from compas_cgal import _coverage_2


def _square(half_width: float) -> np.ndarray:
    return np.array(
        [
            [-half_width, -half_width, 0.0],
            [half_width, -half_width, 0.0],
            [half_width, half_width, 0.0],
            [-half_width, half_width, 0.0],
        ]
    )


def test_uncut_circle_retains_center_until_actually_swept() -> None:
    target = _coverage_2.ExactRegion2.from_polygon(_square(3.0), [])
    coverage = _coverage_2.Coverage2.from_uncut(target)
    assert coverage.residual().exactly_equals(target)
    assert coverage.accumulated_sweeps().is_empty()
    assert coverage.sweep_records == []
    coverage.add_full_circle_sweep(0.0, 0.0, 2.0, 0.0, 1.0)
    assert coverage.residual().contains(0.0, 0.0)
    assert not coverage.residual().contains(2.0, 0.0)
    assert not coverage.residual_is_empty()
    coverage.add_full_circle_sweep(0.0, 0.0, 1.0, 0.0, 4.0)
    assert coverage.residual_is_empty()
    assert target.contains(0.0, 0.0)
    assert coverage.exact_residual_relation()


def test_diagonal_connector_removes_true_capsule_edge_sliver() -> None:
    target = _coverage_2.ExactRegion2.from_polygon(_square(3.0), [])
    coverage = _coverage_2.Coverage2.from_uncut(target)
    coverage.add_segment_sweep(-2.0, -2.0, 2.0, 2.0, 1.0)
    # Squared distance = 2 * 0.7071**2 < 1, outside the conservative quad.
    assert not coverage.residual().contains(0.7071, -0.7071)
    assert coverage.residual().contains(0.708, -0.708)


def test_polygon_holes_remain_excluded_from_uncut_target() -> None:
    target = _coverage_2.ExactRegion2.from_polygon(_square(3.0), [_square(1.0)])
    assert not target.contains(0.0, 0.0)
    assert target.contains(2.0, 0.0)
    assert _coverage_2.Coverage2.from_uncut(target).residual().exactly_equals(target)


def test_uncut_reachable_target_and_invalid_role() -> None:
    domain = _coverage_2.ReachableDomain2(_square(3.0), [], 1.0)
    target = domain.reachable_material()
    assert _coverage_2.Coverage2.from_uncut(target).residual().exactly_equals(target)
    with pytest.raises(_coverage_2.CoverageTransitionError, match="design or reachable"):
        _coverage_2.Coverage2.from_uncut(domain.center_domain())


def test_invalid_polygon_fails_without_offset_construction() -> None:
    with pytest.raises(_coverage_2.InvalidReachableDomainInputError):
        _coverage_2.ExactRegion2.from_polygon(np.zeros((3, 3)), [])


def test_actual_plunge_removes_only_its_cutter_disk() -> None:
    target = _coverage_2.ExactRegion2.from_polygon(_square(3.0), [])
    coverage = _coverage_2.Coverage2.from_uncut(target)
    coverage.add_disk_sweep(0.0, 0.0, 1.0)
    assert not coverage.residual().contains(0.0, 0.0)
    assert not coverage.residual().contains(0.5, 0.5)
    assert coverage.residual().contains(1.0, 1.0)
    assert len(coverage.sweep_records) == 1
    assert coverage.exact_residual_relation()
    coverage.add_disk_sweep(0.0, 0.0, 5.0)
    assert coverage.residual_is_empty()


@pytest.mark.parametrize("center_x,radius", [(float("nan"), 1.0), (0.0, 0.0), (0.0, -1.0), (0.0, float("inf"))])
def test_invalid_plunge_preserves_uncut_target(center_x: float, radius: float) -> None:
    target = _coverage_2.ExactRegion2.from_polygon(_square(3.0), [])
    coverage = _coverage_2.Coverage2.from_uncut(target)
    with pytest.raises(_coverage_2.InvalidCoverageGeometryError):
        coverage.add_disk_sweep(center_x, 0.0, radius)
    assert coverage.residual().exactly_equals(target)
    assert coverage.sweep_records == []
