"""Property tests: BoundaryArclength decisions against an independent exact oracle.

The oracle shares nothing with the implementation: it clears denominators to
integer sums ``sum(n * sqrt(m))``, factors each ``m`` into a square-free prime
set, and decides the sign by eliminating one prime at a time in the field tower
(``A + B*sqrt(p)`` with ``sign(A**2 - p*B**2)``). No enclosures, no folding by
pairwise square-class tests.
"""

from __future__ import annotations

import math
from fractions import Fraction

from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_boundary_arclength import BoundaryArclength

Terms = tuple[tuple[Fraction, Fraction], ...]  # (radicand > 0, coefficient)
Tower = dict[frozenset[int], int]  # square-free prime set -> integer coefficient

MAX_RADICAND_NUMERATOR = 72  # Keeps trial-division factoring and tower size small; still mixes many square classes.
MAX_RADICAND_DENOMINATOR = 12
MAX_COEFFICIENT = 24
MAX_SQUARE_FACTOR = 6  # k in sqrt(k**2 * r) == k * sqrt(r): the folded-class rewrites.

radicands = st.fractions(min_value=Fraction(1, MAX_RADICAND_DENOMINATOR), max_value=MAX_RADICAND_NUMERATOR, max_denominator=MAX_RADICAND_DENOMINATOR)
coefficients = st.fractions(min_value=-MAX_COEFFICIENT, max_value=MAX_COEFFICIENT, max_denominator=MAX_RADICAND_DENOMINATOR)
term_lists = st.lists(st.tuples(radicands, coefficients), min_size=0, max_size=5).map(tuple)


def _prime_factors(value: int) -> list[int]:
    factors, divisor = [], 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            factors.append(divisor)
            value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    return factors


def _tower(terms: Terms) -> Tower:
    """Positive integer multiple of ``sum(c * sqrt(r))`` as ``sum(n * sqrt(prod(primes)))``, square-free keys."""
    # c*sqrt(a/b) == c/b * sqrt(a*b); scale every term by one positive common denominator.
    integer_terms = [(coefficient / radicand.denominator, radicand.numerator * radicand.denominator) for radicand, coefficient in terms]
    scale = math.lcm(1, *(coefficient.denominator for coefficient, _ in integer_terms))
    tower: Tower = {}
    for coefficient, radicand in integer_terms:
        outside, primes = 1, set[int]()
        for prime in _prime_factors(radicand):
            if prime in primes:
                primes.remove(prime)
                outside *= prime
            else:
                primes.add(prime)
        key = frozenset(primes)
        tower[key] = tower.get(key, 0) + int(coefficient * scale) * outside
    return {key: value for key, value in tower.items() if value != 0}


def _product(left: Tower, right: Tower) -> Tower:
    out: Tower = {}
    for left_key, left_value in left.items():
        for right_key, right_value in right.items():
            key = left_key ^ right_key
            out[key] = out.get(key, 0) + left_value * right_value * math.prod(left_key & right_key)
    return {key: value for key, value in out.items() if value != 0}


def _oracle_sign(tower: Tower) -> int:
    """Exact sign of ``sum(n * sqrt(prod(key)))`` by eliminating the largest prime."""
    tower = {key: value for key, value in tower.items() if value != 0}
    if not tower:
        return 0
    primes = set().union(*tower)
    if not primes:
        value = tower[frozenset()]
        return (value > 0) - (value < 0)
    prime = max(primes)
    rest = {key: value for key, value in tower.items() if prime not in key}
    factor = {key - {prime}: value for key, value in tower.items() if prime in key}
    rest_sign, factor_sign = _oracle_sign(rest), _oracle_sign(factor)
    if rest_sign == 0 or factor_sign == 0 or rest_sign == factor_sign:
        return rest_sign or factor_sign
    # Opposite signs: |rest| vs |factor|*sqrt(prime). Equality is impossible (sqrt(prime) is not in the subfield).
    difference = _product(rest, rest)
    for key, value in _product(factor, factor).items():
        difference[key] = difference.get(key, 0) - prime * value
    return rest_sign * _oracle_sign(difference)


def _total(terms: Terms) -> BoundaryArclength:
    total = BoundaryArclength.zero()
    for radicand, coefficient in terms:
        total = total + BoundaryArclength.root(radicand, coefficient)
    return total


def _square_class_rewrite(terms: Terms, factors: list[int]) -> Terms:
    """The same value with each ``c*sqrt(r)`` written as ``(c/k)*sqrt(k**2*r)``: a different representative of its class."""
    return tuple((radicand * k * k, coefficient / k) for (radicand, coefficient), k in zip(terms, factors))


@settings(max_examples=400, deadline=None)
@given(term_lists)
def test_sign_matches_exact_oracle(terms: Terms) -> None:
    assert _total(terms).sign() == _oracle_sign(_tower(terms))


@settings(max_examples=300, deadline=None)
@given(term_lists, st.lists(st.integers(1, MAX_SQUARE_FACTOR), min_size=5, max_size=5))
def test_folded_class_cancellation_is_exact_zero(terms: Terms, factors: list[int]) -> None:
    difference = _total(terms) - _total(_square_class_rewrite(terms, factors))
    assert difference.terms == ()
    assert difference.sign() == 0
    assert _oracle_sign(_tower(terms + tuple((r, -c) for r, c in _square_class_rewrite(terms, factors)))) == 0


@settings(max_examples=300, deadline=None)
@given(term_lists, term_lists, st.lists(st.integers(1, MAX_SQUARE_FACTOR), min_size=5, max_size=5))
def test_partial_folded_tie_leaves_the_remainder_sign(kept: Terms, cancelled: Terms, factors: list[int]) -> None:
    """``kept + cancelled - rewrite(cancelled)`` must decide exactly like ``kept`` alone (ties folded across classes)."""
    rewritten = tuple((radicand, -coefficient) for radicand, coefficient in _square_class_rewrite(cancelled, factors))
    mixed = _total(kept + cancelled + rewritten)
    assert mixed.sign() == _oracle_sign(_tower(kept))
    assert mixed.sign() == _total(kept).sign()


@settings(max_examples=300, deadline=None)
@given(term_lists, coefficients, radicands)
def test_scaled_and_per_root_signs_match_oracle(terms: Terms, factor: Fraction, radicand: Fraction) -> None:
    total = _total(terms)
    assert total.scaled(factor).sign() == _oracle_sign(_tower(tuple((r, c * factor) for r, c in terms)))
    # Dividing by sqrt(radicand) > 0 preserves the sign.
    assert total.per_root(radicand).sign() == _oracle_sign(_tower(terms))


@settings(max_examples=300, deadline=None)
@given(term_lists)
def test_nearest_binary64_is_correctly_rounded(terms: Terms) -> None:
    """The value lies within half an ulp of the returned double on both sides (oracle-checked)."""
    nearest = _total(terms).nearest_binary64()
    exact = Fraction(nearest)
    below, above = Fraction(math.nextafter(nearest, -math.inf)), Fraction(math.nextafter(nearest, math.inf))
    lower_midpoint, upper_midpoint = (below + exact) / 2, (exact + above) / 2
    assert _oracle_sign(_tower(terms + ((Fraction(1), -lower_midpoint),))) >= 0
    assert _oracle_sign(_tower(terms + ((Fraction(1), -upper_midpoint),))) <= 0


MIN_CANCELLATION_POWER = 40  # (sqrt(r) - isqrt(r))**k < 0.5**40: magnitudes far below the 2**-64 first enclosure width.
MAX_CANCELLATION_POWER = 160


@settings(max_examples=200, deadline=None)
@given(
    st.integers(2, MAX_RADICAND_NUMERATOR).filter(lambda r: math.isqrt(r) ** 2 != r),
    st.integers(MIN_CANCELLATION_POWER, MAX_CANCELLATION_POWER),
    st.sampled_from((-1, 1)),
)
def test_near_cancellation_sign_matches_oracle(radicand: int, power: int, orientation: int) -> None:
    """``orientation * (sqrt(r) - isqrt(r))**k`` expanded to ``A + B*sqrt(r)``: huge coefficients, tiny value.

    Random sums never come within the first enclosure width of zero; these do, so an
    unsound enclosure endpoint (a bound off by one root's width) decides them wrongly.
    """
    floor = math.isqrt(radicand)
    rational, irrational = 1, 0  # (A + B*sqrt(r)) * (sqrt(r) - floor) == (B*r - A*floor) + (A - B*floor)*sqrt(r)
    for _ in range(power):
        rational, irrational = irrational * radicand - rational * floor, rational - irrational * floor
    terms: Terms = ((Fraction(1), Fraction(orientation * rational)), (Fraction(radicand), Fraction(orientation * irrational)))
    assert _oracle_sign(_tower(terms)) == orientation
    assert _total(terms).sign() == orientation


def test_oracle_decides_known_values() -> None:
    """Guard the oracle itself on hand-checked identities."""
    assert _oracle_sign(_tower(((Fraction(8), Fraction(1)), (Fraction(2), Fraction(-2))))) == 0  # sqrt8 == 2*sqrt2
    assert _oracle_sign(_tower(((Fraction(2), Fraction(1)), (Fraction(3), Fraction(1)), (Fraction(10), Fraction(-1))))) == -1  # 3.146 < 3.162
    assert _oracle_sign(_tower(((Fraction(1, 2), Fraction(2)), (Fraction(2), Fraction(-1))))) == 0  # 2*sqrt(1/2) == sqrt2
    assert _oracle_sign(_tower(((Fraction(5), Fraction(1)), (Fraction(1), Fraction(-2))))) == 1
