import math

import numpy as np
import pytest

from compas_cgal import _stock_2

SQUARE = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=np.float64)
ISLAND = np.array([[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]], dtype=np.float64)


def test_stock_init_square():
    stock = _stock_2.Stock2(SQUARE, [])
    assert stock.contains(5.0, 5.0)
    assert not stock.contains(-1.0, 5.0)
    assert not stock.contains(11.0, 5.0)
    assert not stock.is_empty()


def test_stock_init_with_island():
    stock = _stock_2.Stock2(SQUARE, [ISLAND])
    assert stock.contains(2.0, 2.0)
    assert not stock.contains(5.0, 5.0)  # inside island = void


def test_stock_rejects_non_simple():
    bowtie = np.array([[0, 0, 0], [6, 0, 0], [0, 4, 0], [6, 4, 0]], dtype=np.float64)
    with pytest.raises(ValueError, match="simple"):
        _stock_2.Stock2(bowtie, [])


def test_subtract_disk_clears_point():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not stock.contains(5.0, 5.0)
    assert not stock.contains(5.9, 5.0)   # inside disk
    assert stock.contains(6.1, 5.0)       # outside disk
    assert stock.contains(2.0, 2.0)


def test_subtract_capsule_geometry():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(2.0, 5.0, 8.0, 5.0, 0.5)
    assert not stock.contains(5.0, 5.0)       # on the axis
    assert not stock.contains(5.0, 5.4)       # inside half-width
    assert stock.contains(5.0, 5.6)           # outside half-width
    assert not stock.contains(8.4, 5.0)       # inside end cap
    assert stock.contains(8.6, 5.0)           # outside end cap


def test_overlapping_subtractions_are_regularized():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(2.0, 5.0, 8.0, 5.0, 0.5)
    stock.subtract_capsule(2.0, 5.2, 8.0, 5.2, 0.5)  # heavy overlap
    stock.subtract_disk(5.0, 5.0, 0.5)               # fully inside cleared
    assert not stock.contains(5.0, 5.3)
    assert stock.contains(5.0, 6.0)


def test_subtract_degenerate_capsule_is_disk():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(5.0, 5.0, 5.0, 5.0, 1.0)
    assert not stock.contains(5.5, 5.0)
    assert stock.contains(6.5, 5.0)


def test_subtract_arc_sweep_full_circle():
    stock = _stock_2.Stock2(SQUARE, [])
    # tool r=0.4 swept around full circle of radius 2 centered at (5,5)
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, 0.4)
    assert not stock.contains(7.0, 5.0)        # on guide circle
    assert not stock.contains(5.0, 7.0)        # opposite side of circle
    assert not stock.contains(7.3, 5.0)        # within tool band (outer)
    assert not stock.contains(6.7, 5.0)        # within tool band (inner)
    assert stock.contains(5.0, 5.0)            # annulus center survives
    assert stock.contains(7.6, 5.0)            # outside band + slack


def test_subtract_arc_sweep_quarter_arc():
    stock = _stock_2.Stock2(SQUARE, [])
    # CCW quarter arc from (7,5) to (5,7) about (5,5), tool r=0.4
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 5.0, 7.0, False, 0.4)
    mid = (5.0 + 2.0 * math.cos(math.pi / 4), 5.0 + 2.0 * math.sin(math.pi / 4))
    assert not stock.contains(*mid)
    assert stock.contains(5.0, 3.0)            # untouched opposite quadrant
    assert stock.contains(3.0, 5.0)


def test_arc_sweep_under_covers_never_over():
    """Conservatism: nothing farther than tool_radius from the guide arc is removed."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, 0.4)
    rng = np.random.default_rng(3)
    for _ in range(400):
        x, y = rng.uniform(0.2, 9.8, 2)
        d_guide = abs(math.hypot(x - 5.0, y - 5.0) - 2.0)
        if d_guide > 0.4:               # strictly outside the true sweep
            assert stock.contains(x, y), (x, y, d_guide)
