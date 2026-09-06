"""Local circle deletion preserves material without assuming a cleared seed."""

import numpy as np
import pytest

from compas_cgal import _stock_2


BOUNDARY = np.array([[-10, -10, 0], [10, -10, 0], [10, 10, 0], [-10, 10, 0]], dtype=float)
NO_CONNECTOR = np.zeros((2, 2))


def test_unique_bulge_retains_current_circle() -> None:
    target = _stock_2.Stock2(BOUNDARY, [])
    assert not target.can_remove_circle((-1, 0, 2), (0, 0, 2), (1, 0, 2), np.array([[-1.0, 2.0], [1.0, 2.0]]), 1.0)
    assert target.exactly_equals(_stock_2.Stock2(BOUNDARY, []))


def test_preserved_connector_can_cover_removed_circle_material() -> None:
    boundary = np.array([[-0.25, 1.75, 0], [0.25, 1.75, 0], [0.25, 2.25, 0], [-0.25, 2.25, 0]])
    target = _stock_2.Stock2(boundary, [])
    previous, current, successor = (-5, 2, 3), (0, 2, 0.5), (5, 2, 3)
    assert not target.can_remove_circle(previous, current, successor, NO_CONNECTOR, 1.0)
    assert target.can_remove_circle(previous, current, successor, np.array([[-2.0, 2.0], [2.0, 2.0]]), 1.0)
    assert target.exactly_equals(_stock_2.Stock2(boundary, []))


def test_outer_disks_cannot_hide_uncut_annulus_hole() -> None:
    target = _stock_2.Stock2(BOUNDARY, [])
    assert not target.can_remove_circle((0, 0, 3), (0, 0, 0.5), (0, 0, 3), NO_CONNECTOR, 1.0)


def test_identical_retained_sweep_authorizes_deletion() -> None:
    target = _stock_2.Stock2(BOUNDARY, [])
    assert target.can_remove_circle((0, 0, 2), (0, 0, 2), (4, 0, 1), NO_CONNECTOR, 1.0)


@pytest.mark.parametrize("invalid", ["shape", "short", "connector_nan", "circle_nan", "negative_guide", "zero_tool"])
def test_invalid_input_rejected_before_empty_local_region(invalid: str) -> None:
    target = _stock_2.Stock2(BOUNDARY, [])
    previous = (30.0, 30.0, 1.0)
    current = (30.0, 30.0, 1.0)  # Entirely outside target: no clipped material.
    connector = NO_CONNECTOR
    tool_radius = 1.0
    if invalid == "shape":
        connector = np.zeros((2, 3))
    elif invalid == "short":
        connector = np.zeros((1, 2))
    elif invalid == "connector_nan":
        connector = np.array([[0.0, 0.0], [float("nan"), 0.0]])
    elif invalid == "circle_nan":
        previous = (float("nan"), 0.0, 1.0)
    elif invalid == "negative_guide":
        previous = (0.0, 0.0, -1.0)
    else:
        tool_radius = 0.0
    with pytest.raises(_stock_2.InvalidCircleRemovalInputError):
        target.can_remove_circle(previous, current, (0.0, 0.0, 1.0), connector, tool_radius)
    assert target.exactly_equals(_stock_2.Stock2(BOUNDARY, []))


def test_permanent_prefix_material_supports_later_deletion() -> None:
    target = _stock_2.Stock2(BOUNDARY, [])
    previous, current, successor = (-5, 0, 1), (0, 0, 2), (5, 0, 1)
    assert not target.can_remove_circle(previous, current, successor, NO_CONNECTOR, 1.0)
    target.subtract_circle_sweep(0.0, 0.0, 2.0, 1.0)
    before = target.clone()
    assert target.can_remove_circle(previous, current, successor, NO_CONNECTOR, 1.0)
    assert target.exactly_equals(before)
