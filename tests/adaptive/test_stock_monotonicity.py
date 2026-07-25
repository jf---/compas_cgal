import math

import numpy as np

from compas_cgal import _stock_2

SQUARE = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]])
LOWER_HALF = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 5.0, 0.0], [0.0, 5.0, 0.0]])
LOWER_THREE_QUARTERS = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 5.25, 0.0], [0.0, 5.25, 0.0]])


def _cap_exceeded(
    stock: _stock_2.Stock2,
    cap_ratio: float,
    gap_close_ratio: float = 0.0,
) -> bool:
    return _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio, gap_close_ratio)[2]


def _assert_cap_implication(
    smaller: _stock_2.Stock2,
    larger: _stock_2.Stock2,
    cap_ratio: float,
    gap_close_ratio: float = 0.0,
) -> None:
    assert smaller.is_subset_of(larger)
    assert _cap_exceeded(smaller, cap_ratio, gap_close_ratio)
    assert _cap_exceeded(larger, cap_ratio, gap_close_ratio)


def _slotted_lower_half() -> _stock_2.Stock2:
    stock = _stock_2.Stock2(LOWER_HALF, [])
    stock.subtract_capsule(5.0, 4.30, 5.0, 4.53, 0.05)
    return stock


def test_stock_inclusion_and_cap_fixture_is_non_vacuous() -> None:
    smaller = _stock_2.Stock2(LOWER_HALF, [])
    larger = _stock_2.Stock2(SQUARE, [])

    assert smaller.is_subset_of(larger)
    assert not _cap_exceeded(smaller, _stock_2.cap_chord_ratio(math.pi))
    assert _cap_exceeded(larger, _stock_2.cap_chord_ratio(math.pi))


def test_exact_cap_predicate_is_monotone_under_stock_inclusion() -> None:
    smaller = _stock_2.Stock2(LOWER_THREE_QUARTERS, [])
    larger = _stock_2.Stock2(SQUARE, [])

    _assert_cap_implication(smaller, larger, _stock_2.cap_chord_ratio(math.pi))


def test_cap_monotonicity_for_multiple_runs() -> None:
    _assert_cap_implication(
        _slotted_lower_half(),
        _stock_2.Stock2(SQUARE, []),
        _stock_2.cap_chord_ratio(math.radians(60.0)),
    )


def test_cap_monotonicity_for_tangent_and_coincident_boundaries() -> None:
    tangent = _stock_2.Stock2(SQUARE, [])
    tangent.subtract_disk(5.0, 5.25, 0.25)
    _assert_cap_implication(
        tangent,
        _stock_2.Stock2(SQUARE, []),
        _stock_2.cap_chord_ratio(math.pi),
    )
    _assert_cap_implication(
        _stock_2.Stock2(LOWER_HALF, []),
        _stock_2.Stock2(SQUARE, []),
        3.0,
    )


def test_cap_monotonicity_for_seam_merged_run() -> None:
    seam = _stock_2.Stock2(SQUARE, [])
    seam.subtract_capsule(3.0, 0.0, 3.0, 10.0, 2.0)

    _assert_cap_implication(
        seam,
        _stock_2.Stock2(SQUARE, []),
        _stock_2.cap_chord_ratio(math.pi),
    )


def test_cap_monotonicity_with_nonzero_gap_closure_pessimism() -> None:
    _assert_cap_implication(
        _slotted_lower_half(),
        _stock_2.Stock2(SQUARE, []),
        _stock_2.cap_chord_ratio(math.radians(135.0)),
        _stock_2.cap_chord_ratio(math.radians(30.0)),
    )
