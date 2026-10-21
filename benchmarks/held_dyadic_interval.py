"""Rigorous dyadic interval enclosures: the filter stage of exact predicates in Python.

An :class:`Interval` encloses a real value in ``[lo, hi] · 2^-FRACTION_BITS``
with integer ends. Every operation rounds its ends outward with integer floor
and ceiling (``>>`` floors in Python), so an interval always contains the exact
value; no floating-point arithmetic and no error derivation is involved. A
sign is decided by the filter only when the interval excludes zero; otherwise
the caller falls back to the exact rational / one-root predicate. This is the
interval-then-exact filter the native kernels use, with the same guarantee:
the filter never decides wrongly, it only declines.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.held_one_root import OneRootNumber

# Fractional bits of the integer interval ends: 2^-96 mm, so enclosures of
# corpus quantities (|x| < 2^12 mm, squares < 2^24 mm²) stay tight far below
# any exact tie the filter must leave to the exact stage.
FRACTION_BITS = 96
_ONE = 1 << FRACTION_BITS


def _floor_div(numerator: int, denominator: int) -> int:
    return numerator // denominator


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


@dataclass(frozen=True, slots=True)
class Interval:
    """The closed interval ``[lo, hi] · 2^-FRACTION_BITS`` (``lo <= hi``)."""

    lo: int
    hi: int

    @classmethod
    def of(cls, value: Fraction | int | float) -> Interval:
        """Tightest enclosure of a rational (a double is a rational)."""
        exact = Fraction(value)
        scaled = exact.numerator * _ONE
        return cls(_floor_div(scaled, exact.denominator), _ceil_div(scaled, exact.denominator))

    @classmethod
    def sqrt_of(cls, value: Fraction) -> Interval:
        """Enclosure of ``√value`` for a rational ``value >= 0``."""
        numerator, denominator = value.numerator, value.denominator
        # √(p/q) = √(p·q)/q; at scale 2^K: isqrt(p·q·4^K)/q, rounded outward.
        root = math.isqrt(numerator * denominator * _ONE * _ONE)
        return cls(_floor_div(root, denominator), _ceil_div(root + 1, denominator))

    @classmethod
    def of_one_root(cls, value: OneRootNumber) -> Interval:
        """Enclosure of ``a + b·√c``."""
        if value.is_rational:
            return cls.of(value.rational)
        return cls.of(value.rational) + cls.of(value.coefficient) * cls.sqrt_of(value.radicand)

    def __add__(self, other: Interval) -> Interval:
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other: Interval) -> Interval:
        return Interval(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self) -> Interval:
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other: Interval) -> Interval:
        products = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return Interval(min(products) >> FRACTION_BITS, -((-max(products)) >> FRACTION_BITS))

    def __truediv__(self, other: Interval) -> Interval:
        """Enclosure of ``x / y`` for a divisor enclosure that excludes zero.

        Raises:
            ZeroDivisionError: If the divisor enclosure contains zero.
        """
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("Interval divisor encloses zero.")
        quotients = [(numerator << FRACTION_BITS, denominator) for numerator in (self.lo, self.hi) for denominator in (other.lo, other.hi)]
        return Interval(min(_floor_div(n, d) if d > 0 else _floor_div(-n, -d) for n, d in quotients), max(_ceil_div(n, d) if d > 0 else _ceil_div(-n, -d) for n, d in quotients))

    def midpoint(self) -> float:
        """A binary64 inside the enclosure's hull (a proposal, never a decision)."""
        return float(Fraction(self.lo + self.hi, 2 * _ONE))

    def square(self) -> Interval:
        """Enclosure of ``x²`` (non-negative)."""
        low, high = abs(self.lo), abs(self.hi)
        if self.lo <= 0 <= self.hi:
            return Interval(0, -((-(max(low, high) ** 2)) >> FRACTION_BITS))
        small, large = min(low, high), max(low, high)
        return Interval((small * small) >> FRACTION_BITS, -((-(large * large)) >> FRACTION_BITS))

    def sqrt(self) -> Interval:
        """Enclosure of ``√x`` over the non-negative part."""
        low, high = max(self.lo, 0), max(self.hi, 0)
        return Interval(math.isqrt(low * _ONE), math.isqrt(high * _ONE) + 1)

    def absolute(self) -> Interval:
        """Enclosure of ``|x|``."""
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return Interval(0, max(-self.lo, self.hi))

    def sign(self) -> int | None:
        """``1`` or ``-1`` when the interval excludes zero, ``None`` when the filter declines."""
        if self.lo > 0:
            return 1
        if self.hi < 0:
            return -1
        return None


def minimum(intervals: list[Interval]) -> Interval:
    """Enclosure of the minimum of enclosed values."""
    return Interval(min(item.lo for item in intervals), min(item.hi for item in intervals))


def maximum(intervals: list[Interval]) -> Interval:
    """Enclosure of the maximum of enclosed values."""
    return Interval(max(item.lo for item in intervals), max(item.hi for item in intervals))
