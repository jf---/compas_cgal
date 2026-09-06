"""Circle-sweep inputs reach exact stock without Python radius arithmetic."""

import numpy as np
import pytest

from compas_cgal import _stock_2


BOUNDARY = np.array([[-10, -10, 0], [10, -10, 0], [10, 10, 0], [-10, 10, 0]], dtype=np.float64)


def test_circle_sweep_keeps_the_uncut_central_island() -> None:
    stock = _stock_2.Stock2(BOUNDARY, [])
    stock.subtract_circle_sweep(0, 0, 3, 1)
    assert stock.contains(0, 0)
    assert not stock.contains(3, 0)
    assert stock.contains(5, 0)


def test_circle_sweep_with_larger_tool_clears_the_center() -> None:
    stock = _stock_2.Stock2(BOUNDARY, [])
    stock.subtract_circle_sweep(0, 0, 1, 2)
    assert not stock.contains(0, 0)
    assert not stock.contains(2, 0)
    assert stock.contains(4, 0)


def test_circle_sweep_adds_radii_exactly_after_injection() -> None:
    stock = _stock_2.Stock2(BOUNDARY, [])
    rounded = stock.clone()
    stock.subtract_circle_sweep(-0.3, 0, 0.1, 0.2)
    rounded.subtract_annulus(-0.3, 0, 0, 0.1 + 0.2)
    # The exact outer endpoint is 2**-55 mm; the rounded one is 2**-54 mm.
    # Their binary64 midpoint witnesses excess removal, away from boundaries.
    probe_x = 3 * 2.0**-56
    assert stock.contains(probe_x, 0)
    assert not rounded.contains(probe_x, 0)


def test_invalid_circle_sweep_does_not_mutate_stock() -> None:
    stock = _stock_2.Stock2(BOUNDARY, [])
    before = stock.clone()
    with pytest.raises(_stock_2.InvalidAnnulusRadiiError):
        stock.subtract_circle_sweep(0, 0, -1, 1)
    with pytest.raises(_stock_2.NonFiniteAnnulusInputError):
        stock.subtract_circle_sweep(float("nan"), 0, 1, 1)
    assert stock.exactly_equals(before)
