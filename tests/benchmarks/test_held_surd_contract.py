"""Contract of the exact surd ``rational + radical * sqrt(radicand)`` behind the placement predicates.

Characterises ``_Surd`` independently of its representation (F5 lever L2 proposed a
``NamedTuple``; this contract is what any faster representation must keep): exact arithmetic
and sign, one shared radicand, immutability, value equality, the canonical
``repr`` (it reaches path digests through frozen dataclasses), and a fail-loud
operator surface: no ordering, no equality with bare tuples, no sequence
repetition or concatenation.
"""

from __future__ import annotations

import math
from decimal import Decimal
from decimal import localcontext
from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st

from benchmarks.held_standard_placement import MixedRadicandError
from benchmarks.held_standard_placement import _Surd

# Coefficient bound for the property tests: 62-bit integers, the size of the
# dyadic-scaled coordinates the placement predicates square (2**62 > the 2**50
# scale of Figure 5 coordinates with headroom).
COEFFICIENT_BITS = 62
_coefficients = st.integers(min_value=-(2**COEFFICIENT_BITS), max_value=2**COEFFICIENT_BITS)
_radicands = st.integers(min_value=0, max_value=2**COEFFICIENT_BITS)
# Decimal digits for the sign oracle: a non-zero ``a + b*sqrt(c)`` with 62-bit
# integers is at least ``1 / |a - b*sqrt(c)|`` >= 2**-126 in magnitude, i.e.
# about 38 decimal digits below the 2**126 operand scale; 120 digits dwarf that.
ORACLE_DIGITS = 120


def _oracle_sign(rational: int, radical: int, radicand: int) -> int:
    root = math.isqrt(radicand)
    if root * root == radicand:
        exact = rational + radical * root
        return (exact > 0) - (exact < 0)
    with localcontext() as context:
        context.prec = ORACLE_DIGITS
        value = Decimal(rational) + Decimal(radical) * Decimal(radicand).sqrt()
    return (value > 0) - (value < 0)


@given(_coefficients, _coefficients, _coefficients, _coefficients, _radicands)
def test_arithmetic_is_exact_over_one_radicand(a: int, b: int, c: int, d: int, radicand: int) -> None:
    first, second = _Surd(a, b, radicand), _Surd(c, d, radicand)
    assert first + second == _Surd(a + c, b + d, radicand)
    assert first - second == _Surd(a - c, b - d, radicand)
    assert first * second == _Surd(a * c + b * d * radicand, a * d + b * c, radicand)
    assert first.scaled(c) == _Surd(a * c, b * c, radicand)
    assert first.plus(c) == _Surd(a + c, b, radicand)


@given(_coefficients, _coefficients, _radicands)
def test_sign_is_the_exact_sign(rational: int, radical: int, radicand: int) -> None:
    assert _Surd(rational, radical, radicand).sign() == _oracle_sign(rational, radical, radicand)


@given(_coefficients, st.integers(min_value=-(2**31), max_value=2**31))
def test_sign_of_an_exact_zero_on_a_perfect_square_radicand(root: int, radical: int) -> None:
    assert _Surd(-radical * abs(root), radical, root * root).sign() == 0


def test_mixed_radicands_raise_by_name() -> None:
    with pytest.raises(MixedRadicandError):
        _Surd(1, 1, 2) + _Surd(1, 1, 3)
    with pytest.raises(MixedRadicandError):
        _Surd(1, 1, 2) - _Surd(1, 1, 3)
    with pytest.raises(MixedRadicandError):
        _Surd(1, 1, 2) * _Surd(1, 1, 3)


def test_fields_are_immutable() -> None:
    surd = _Surd(1, 2, 3)
    with pytest.raises(AttributeError):
        surd.rational = 5  # type: ignore[misc]


def test_value_equality_hash_and_repr() -> None:
    assert _Surd(1, 2, 3) == _Surd(1, 2, 3)
    assert _Surd(1, 2, 3) != _Surd(1, 2, 4)
    assert hash(_Surd(1, 2, 3)) == hash(_Surd(1, 2, 3))
    assert repr(_Surd(1, -2, 3)) == "_Surd(rational=1, radical=-2, radicand=3)"


def test_a_bare_tuple_is_never_equal_to_a_surd() -> None:
    # ``Any``: the checks below are deliberately ill-typed; the contract is the runtime refusal.
    surd: Any = _Surd(1, 2, 3)
    assert surd != (1, 2, 3)
    assert (1, 2, 3) != surd


def test_surds_have_no_order() -> None:
    first: Any = _Surd(1, 2, 3)
    second: Any = _Surd(2, 0, 3)
    with pytest.raises(TypeError):
        first < second  # noqa: B015
    with pytest.raises(TypeError):
        first <= second  # noqa: B015
    with pytest.raises(TypeError):
        first > second  # noqa: B015
    with pytest.raises(TypeError):
        first >= second  # noqa: B015
    with pytest.raises(TypeError):
        (1, 2, 3) < second  # noqa: B015


def test_no_sequence_repetition_or_concatenation() -> None:
    surd: Any = _Surd(1, 2, 3)
    with pytest.raises(TypeError):
        2 * surd  # noqa: B018
    with pytest.raises((TypeError, AttributeError)):
        surd * 2  # noqa: B018
    with pytest.raises(TypeError):
        (1, 2) + surd  # noqa: B018
    with pytest.raises((TypeError, AttributeError)):
        surd + (1, 2, 3)  # noqa: B018
