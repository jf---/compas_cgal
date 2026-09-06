"""Residual-only exact replay agrees with the chronological coverage engine."""

import numpy as np
import pytest

from compas_cgal import _coverage_2


def _target() -> _coverage_2.ExactRegion2:
    return _coverage_2.ExactRegion2.from_polygon(np.array([[-3, -3, 0], [3, -3, 0], [3, 3, 0], [-3, 3, 0]], dtype=float), [])


def test_mixed_replay_matches_existing_exact_coverage_and_preserves_target() -> None:
    target = _target()
    circles = np.array([[0.0, 0.0, 2.0]])
    segments = np.array([[-2.0, -2.0, 2.0, 2.0]])
    disks = np.array([[1.0, -1.0]])
    expected = _coverage_2.Coverage2.from_uncut(target)
    expected.add_full_circle_sweep(0.0, 0.0, 2.0, 0.0, 1.0)
    expected.add_segment_sweep(-2.0, -2.0, 2.0, 2.0, 1.0)
    expected.add_disk_sweep(1.0, -1.0, 1.0)
    actual = _coverage_2.remaining_material(target, circles, segments, disks, 1.0)
    assert actual.exactly_equals(expected.residual())
    assert target.exactly_equals(_target())


def test_subgrid_central_island_is_not_rounded_away() -> None:
    target = _target()
    # The annulus leaves a radius 1/1024 mm island, far below a 0.25 mm grid.
    actual = _coverage_2.remaining_material(target, np.array([[0.0, 0.0, 1.0009765625]]), np.empty((0, 4)), np.empty((0, 2)), 1.0)
    assert actual.contains(0.0, 0.0)
    assert not actual.contains(0.01, 0.0)
    assert not actual.is_empty()


def test_empty_motion_returns_unchanged_region() -> None:
    target = _target()
    assert _coverage_2.remaining_material(target, np.empty((0, 3)), np.empty((0, 4)), np.empty((0, 2)), 1.0).exactly_equals(target)


@pytest.mark.parametrize("invalid", ["circle", "segment", "disk", "radius", "shape"])
def test_invalid_tail_rejected_even_when_first_circle_clears_target(invalid: str) -> None:
    target = _target()
    circles = np.array([[0.0, 0.0, 1.0]])
    segments = np.empty((0, 4))
    disks = np.empty((0, 2))
    radius = 5.0
    if invalid == "circle":
        circles = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
    elif invalid == "segment":
        segments = np.zeros((1, 4))
    elif invalid == "disk":
        disks = np.array([[float("nan"), 0.0]])
    elif invalid == "radius":
        radius = float("inf")
    else:
        circles = np.empty((0, 2))
    with pytest.raises(_coverage_2.InvalidCoverageGeometryError):
        _coverage_2.remaining_material(target, circles, segments, disks, radius)
    assert target.exactly_equals(_target())


def test_actual_disk_can_complete_removal() -> None:
    assert _coverage_2.remaining_material(_target(), np.empty((0, 3)), np.empty((0, 4)), np.array([[0.0, 0.0]]), 5.0).is_empty()


def test_residual_role_cannot_be_reinterpreted_as_fresh_target() -> None:
    residual = _coverage_2.Coverage2.from_uncut(_target()).residual()
    with pytest.raises(_coverage_2.CoverageTransitionError, match="design or reachable"):
        _coverage_2.remaining_material(residual, np.empty((0, 3)), np.empty((0, 4)), np.empty((0, 2)), 1.0)
