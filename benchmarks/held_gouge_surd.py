"""Exact numbers of one quadratic extension over the integers, ``(a + b·√K) / D``.

This is the Python mirror of CGAL's one-root ``Sqrt_extension``: ``a``, ``b``,
``K >= 0`` and ``D > 0`` are integers, and ``sign`` is decided exactly by sign
analysis and one squaring. Every gouge decision in the Held gouge audit is the
sign of one such number, so no decision ever depends on a tolerance.
Integers (not ``Fraction``) carry the arithmetic: the audit puts every input
double on one common dyadic grid, where all coordinates are integers.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from decimal import localcontext

from benchmarks.errors import BenchmarkError

# Decimal digits used to REPORT an exact surd as a double (never to decide).
# A depth is r - dist, a difference of two O(100) mm values; the deepest
# rounding-scale gouge seen is 2e-17 mm, so ~20 digits of cancellation plus 17
# for the double need ~40; 120 keeps proof-level slack for finer inputs.
REPORT_DIGITS = 120


class SurdError(BenchmarkError):
    """Base error of the exact quadratic-surd arithmetic."""


class InvalidSurdError(SurdError):
    """A surd was built over a negative radicand or a non-positive denominator."""


def sign_surd(a: int, b: int, k: int) -> int:
    """Exact ``sign(a + b·√k)`` for integers, ``k >= 0``."""
    if b == 0 or k == 0:
        return (a > 0) - (a < 0)
    sb = 1 if b > 0 else -1
    if a == 0:
        return sb
    sa = 1 if a > 0 else -1
    if sa == sb:
        return sa
    # a and b·√k have opposite signs: |a| against |b|·√k, compared squared.
    a2, b2k = a * a, b * b * k
    if a2 == b2k:
        return 0
    return sa if a2 > b2k else sb


@dataclass(frozen=True, slots=True)
class Surd:
    """The exact real number ``(rational + coefficient·√radicand) / denominator``."""

    rational: int
    coefficient: int
    radicand: int
    denominator: int

    @classmethod
    def build(cls, rational: int, coefficient: int, radicand: int, denominator: int) -> Surd:
        """Validated factory.

        Raises:
            InvalidSurdError: ``radicand < 0`` or ``denominator <= 0``.
        """
        if radicand < 0 or denominator <= 0:
            raise InvalidSurdError(f"radicand {radicand} must be >= 0 and denominator {denominator} > 0")
        return cls(rational, coefficient, radicand, denominator)

    @classmethod
    def of(cls, numerator: int, denominator: int = 1) -> Surd:
        """The rational ``numerator / denominator`` (``denominator > 0``)."""
        return cls(numerator, 0, 0, denominator)

    def sign(self) -> int:
        """Exact sign in ``{-1, 0, 1}``."""
        return sign_surd(self.rational, self.coefficient, self.radicand)

    def compare_square(self, t: int) -> int:
        """Exact ``sign(self - t²)``."""
        return sign_surd(self.rational - t * t * self.denominator, self.coefficient, self.radicand)

    def to_decimal(self) -> Decimal:
        """Reporting-only decimal evaluation at ``REPORT_DIGITS`` digits."""
        with localcontext() as context:
            context.prec = REPORT_DIGITS
            value = Decimal(self.rational)
            if self.coefficient != 0 and self.radicand != 0:
                value += Decimal(self.coefficient) * Decimal(self.radicand).sqrt()
            return value / Decimal(self.denominator)
