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
