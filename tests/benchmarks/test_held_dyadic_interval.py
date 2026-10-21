"""Dyadic interval enclosures always contain the exact value (the filter never decides wrongly)."""

from __future__ import annotations

from fractions import Fraction

from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_dyadic_interval import FRACTION_BITS
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_one_root import OneRootNumber

_SCALE = Fraction(1, 2**FRACTION_BITS)
_values = st.fractions(min_value=-4096, max_value=4096, max_denominator=2**80)
_positive = st.fractions(min_value=Fraction(1, 2**40), max_value=4096, max_denominator=2**80)


def _contains(box: Interval, value: Fraction) -> bool:
    return box.lo * _SCALE <= value <= box.hi * _SCALE


@settings(max_examples=300, deadline=None)
@given(_values, _values)
def test_arithmetic_encloses_the_exact_result(first: Fraction, second: Fraction) -> None:
    a, b = Interval.of(first), Interval.of(second)
    assert _contains(a + b, first + second)
    assert _contains(a - b, first - second)
    assert _contains(a * b, first * second)
    assert _contains(a.square(), first * first)
    assert _contains(a.absolute(), abs(first))
    if second != 0 and b.sign() is not None:
        assert _contains(a / b, first / second)


@settings(max_examples=300, deadline=None)
@given(_positive, _values, _values)
def test_roots_enclose_the_exact_root(radicand: Fraction, rational: Fraction, coefficient: Fraction) -> None:
    root = Interval.sqrt_of(radicand)
    assert (root.lo * _SCALE) ** 2 <= radicand <= (root.hi * _SCALE) ** 2
    assert (Interval.of(radicand).sqrt().lo * _SCALE) ** 2 <= radicand <= (Interval.of(radicand).sqrt().hi * _SCALE) ** 2
    number = OneRootNumber.build(rational, coefficient, radicand)
    box = Interval.of_one_root(number)
    assert number.compare(box.lo * _SCALE) >= 0 and number.compare(box.hi * _SCALE) <= 0


@settings(max_examples=300, deadline=None)
@given(_values)
def test_a_decided_sign_is_the_exact_sign(value: Fraction) -> None:
    decided = Interval.of(value).sign()
    if decided is not None:
        assert decided == (value > 0) - (value < 0)
    else:
        assert abs(value) <= _SCALE
