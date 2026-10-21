"""Exact boundary arclength: sums ``sum(c * sqrt(r))`` over positive rational radicands.

Square-class folding keeps the representation canonical, so zero is decided
structurally and nonzero signs by certified ``isqrt`` enclosures; no float
enters a decision.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Self

FIRST_ENCLOSURE_BITS = 64  # Initial binary digits of each certified square-root enclosure; doubled until a decision separates.


def _rational_square_root(value: Fraction) -> Fraction | None:
    """The exact rational square root of ``value >= 0``, or None when it is irrational."""
    numerator, denominator = math.isqrt(value.numerator), math.isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Fraction(numerator, denominator)
    return None


def _root_enclosure(radicand: Fraction, bits: int) -> tuple[Fraction, Fraction]:
    """Rational bounds ``lo <= sqrt(radicand) <= hi`` with ``hi - lo <= 2**-bits / denominator``."""
    scaled = radicand.numerator * radicand.denominator << (2 * bits)
    root = math.isqrt(scaled)
    scale = radicand.denominator << bits
    return Fraction(root, scale), Fraction(root + (root * root != scaled), scale)


@dataclass(frozen=True)
class BoundaryArclength:
    """Exact arclength ``sum(c * sqrt(r))`` along the polygonal Ω, ``r > 0`` rational.

    Folded so no two radicands differ by a rational-square factor; rational
    terms carry radicand 1. Square roots of positive rationals pairwise in
    distinct square classes are linearly independent over Q (Besicovitch),
    so the sum is zero iff it has no terms, and a nonzero sum's sign is
    settled by certified ``isqrt`` enclosures refined until they exclude 0.
    """

    terms: tuple[tuple[Fraction, Fraction], ...]

    @classmethod
    def zero(cls) -> Self:
        return cls(())

    @classmethod
    def root(cls, radicand: Fraction, coefficient: Fraction) -> Self:
        """``coefficient * sqrt(radicand)`` for ``radicand >= 0``."""
        return cls.zero()._folded(((radicand, coefficient),))

    def _folded(self, extra: tuple[tuple[Fraction, Fraction], ...]) -> Self:
        classes = dict(self.terms)
        for radicand, coefficient in extra:
            if radicand == 0 or coefficient == 0:
                continue
            square_root = _rational_square_root(radicand)
            if square_root is not None:
                radicand, coefficient = Fraction(1), coefficient * square_root
            for representative in classes:
                ratio = _rational_square_root(radicand * representative)
                if ratio is not None:
                    # sqrt(radicand) == sqrt(radicand * representative) / representative * sqrt(representative)
                    classes[representative] += coefficient * ratio / representative
                    break
            else:
                classes[radicand] = coefficient
        return type(self)(tuple((radicand, coefficient) for radicand, coefficient in classes.items() if coefficient != 0))

    def __add__(self, other: BoundaryArclength) -> BoundaryArclength:
        return self._folded(other.terms)

    def __sub__(self, other: BoundaryArclength) -> BoundaryArclength:
        return self._folded(tuple((radicand, -coefficient) for radicand, coefficient in other.terms))

    def scaled(self, factor: Fraction) -> BoundaryArclength:
        return self.zero()._folded(tuple((radicand, coefficient * factor) for radicand, coefficient in self.terms))

    def per_root(self, radicand: Fraction) -> BoundaryArclength:
        """``self / sqrt(radicand)`` for ``radicand > 0``: ``c*sqrt(r)/sqrt(s) == (c/s)*sqrt(r*s)``."""
        return self.zero()._folded(tuple((term * radicand, coefficient / radicand) for term, coefficient in self.terms))

    def rational(self) -> Fraction | None:
        """The exact value when it is rational (only the radicand-1 class), else None."""
        if not self.terms:
            return Fraction(0)
        if len(self.terms) == 1 and self.terms[0][0] == 1:
            return self.terms[0][1]
        return None

    def enclosure(self, bits: int) -> tuple[Fraction, Fraction]:
        low = high = Fraction(0)
        for radicand, coefficient in self.terms:
            root_low, root_high = _root_enclosure(radicand, bits)
            low += coefficient * (root_low if coefficient > 0 else root_high)
            high += coefficient * (root_high if coefficient > 0 else root_low)
        return low, high

    def sign(self) -> int:
        exact = self.rational()
        if exact is not None:
            return (exact > 0) - (exact < 0)
        bits = FIRST_ENCLOSURE_BITS
        while True:
            low, high = self.enclosure(bits)
            if low > 0:
                return 1
            if high < 0:
                return -1
            bits *= 2

    def nearest_binary64(self) -> float:
        """The correctly rounded binary64 value: enclosures refine until both ends round alike."""
        exact = self.rational()
        if exact is not None:
            return float(exact)
        bits = FIRST_ENCLOSURE_BITS
        while True:
            low, high = self.enclosure(bits)
            # Rounding boundaries are rational, so an irrational value is never on one.
            if float(low) == float(high):
                return float(low)
            bits *= 2
