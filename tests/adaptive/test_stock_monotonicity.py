import math

import numpy as np

from compas_cgal import _stock_2

SQUARE = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]])
LOWER_HALF = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 5.0, 0.0], [0.0, 5.0, 0.0]])
LOWER_THREE_QUARTERS = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 5.25, 0.0], [0.0, 5.25, 0.0]])


def _cap_exceeded(stock: _stock_2.Stock2, cap_radians: float) -> bool:
    cap_ratio = _stock_2.cap_chord_ratio(cap_radians)
    return _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio)[2]


def test_stock_inclusion_and_cap_fixture_is_non_vacuous() -> None:
    smaller = _stock_2.Stock2(LOWER_HALF, [])
    larger = _stock_2.Stock2(SQUARE, [])

    assert smaller.is_subset_of(larger)
    assert not _cap_exceeded(smaller, math.pi)
    assert _cap_exceeded(larger, math.pi)


def test_exact_cap_predicate_is_monotone_under_stock_inclusion() -> None:
    smaller = _stock_2.Stock2(LOWER_THREE_QUARTERS, [])
    larger = _stock_2.Stock2(SQUARE, [])

    assert smaller.is_subset_of(larger)
    assert _cap_exceeded(smaller, math.pi)
    assert _cap_exceeded(larger, math.pi)
