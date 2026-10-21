"""Exact one-root numbers ``a + b√c`` and points over them, for the Held centre domain Ω*.

One-root numbers are ``a + b√c`` over rationals; comparisons across different
roots are decided exactly by squaring (``sign(u + w√β)`` with ``u, w ∈ ℚ(√α)``
reduces to ``sign(u)``, ``sign(w)`` and ``sign(u² − β·w²)``), never by
``float``. Arithmetic across different roots is refused
(:class:`CrossRootArithmeticError`). Design: docs/held_exact_centre_domain.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Self
from typing import TypeAlias

from benchmarks.held_centre_domain_errors import CentreDomainError
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import Sense
from compas_cgal.units import WorldXY

_Scalar: TypeAlias = "OneRootNumber | Fraction | int"


class NegativeRadicandError(CentreDomainError):
    """A one-root number was given a negative radicand."""


class CrossRootArithmeticError(CentreDomainError):
    """Arithmetic was attempted between numbers of different quadratic extensions."""


class MixedRootPointError(CentreDomainError):
    """A one-root point's coordinates do not share one quadratic extension."""


# ---------------------------------------------------------------------------
# One-root numbers
# ---------------------------------------------------------------------------


def _is_rational_square(value: Fraction) -> bool:
    numerator, denominator = value.numerator, value.denominator
    return math.isqrt(numerator) ** 2 == numerator and math.isqrt(denominator) ** 2 == denominator


def _rational_sqrt(value: Fraction) -> Fraction:
    return Fraction(math.isqrt(value.numerator), math.isqrt(value.denominator))


def _sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


@dataclass(frozen=True, eq=False)
class OneRootNumber:
    """Exact ``a + b·√c`` with rational ``a, b, c``, ``c >= 0``.

    The representation is folded: ``b == 0`` iff ``c == 0``, and ``c`` is never a
    rational square. Equality and order are exact value comparisons (also across
    different radicands), so the type is deliberately unhashable.
    """

    rational: Fraction
    coefficient: Fraction
    radicand: Fraction

    __hash__ = None  # type: ignore[assignment]

    @classmethod
    def build(cls, rational: Fraction | int, coefficient: Fraction | int = 0, radicand: Fraction | int = 0) -> Self:
        """Validate and fold ``a + b√c``.

        Args:
            rational: The rational part ``a``.
            coefficient: The coefficient ``b`` of the root.
            radicand: The radicand ``c >= 0``.

        Returns:
            The folded number.

        Raises:
            NegativeRadicandError: If ``c < 0``.
        """
        a, b, c = Fraction(rational), Fraction(coefficient), Fraction(radicand)
        if c < 0:
            raise NegativeRadicandError(f"Radicand {c} is negative.")
        if b == 0 or c == 0:
            return cls(a, Fraction(0), Fraction(0))
        if _is_rational_square(c):
            return cls(a + b * _rational_sqrt(c), Fraction(0), Fraction(0))
        return cls(a, b, c)

    @classmethod
    def of(cls, value: Fraction | int) -> Self:
        """The rational ``value`` as a one-root number."""
        return cls(Fraction(value), Fraction(0), Fraction(0))

    @property
    def is_rational(self) -> bool:
        """Whether the number is rational (folded root is zero)."""
        return self.radicand == 0

    def sign(self) -> int:
        """Exact sign: ``-1``, ``0`` or ``1``."""
        sign_a, sign_b = _sign(self.rational), _sign(self.coefficient)
        if sign_b == 0:
            return sign_a
        if sign_a * sign_b >= 0:
            return sign_b if sign_a == 0 else sign_a
        return sign_a * _sign(self.rational * self.rational - self.coefficient * self.coefficient * self.radicand)

    def common_radicand(self, other: OneRootNumber) -> Fraction:
        """The shared radicand of two numbers of one extension.

        Raises:
            CrossRootArithmeticError: If the two numbers lie in different extensions.
        """
        if self.radicand == 0:
            return other.radicand
        if other.radicand == 0 or other.radicand == self.radicand:
            return self.radicand
        raise CrossRootArithmeticError(f"√{self.radicand} and √{other.radicand} are different extensions.")

    def __add__(self, other: _Scalar) -> OneRootNumber:
        right = _lift(other)
        radicand = self.common_radicand(right)
        return OneRootNumber.build(self.rational + right.rational, self.coefficient + right.coefficient, radicand)

    def __radd__(self, other: _Scalar) -> OneRootNumber:
        return self + other

    def __neg__(self) -> OneRootNumber:
        return OneRootNumber(-self.rational, -self.coefficient, self.radicand)

    def __sub__(self, other: _Scalar) -> OneRootNumber:
        return self + (-_lift(other))

    def __rsub__(self, other: _Scalar) -> OneRootNumber:
        return _lift(other) - self

    def __mul__(self, other: _Scalar) -> OneRootNumber:
        right = _lift(other)
        radicand = self.common_radicand(right)
        rational = self.rational * right.rational + self.coefficient * right.coefficient * radicand
        coefficient = self.rational * right.coefficient + self.coefficient * right.rational
        return OneRootNumber.build(rational, coefficient, radicand)

    def __rmul__(self, other: _Scalar) -> OneRootNumber:
        return self * other

    def compare(self, other: _Scalar) -> int:
        """Exact three-way comparison, also across different extensions.

        Returns:
            ``-1``, ``0`` or ``1`` as ``self`` is smaller, equal or larger.
        """
        right = _lift(other)
        if self.radicand == 0 or right.radicand == 0 or self.radicand == right.radicand:
            return (self - right).sign()
        difference = OneRootNumber.build(self.rational - right.rational, self.coefficient, self.radicand)
        return _sign_of_sum(difference, OneRootNumber.of(-right.coefficient), right.radicand)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (OneRootNumber, Fraction, int)):
            return NotImplemented
        return self.compare(other) == 0

    def __lt__(self, other: _Scalar) -> bool:
        return self.compare(other) < 0

    def __le__(self, other: _Scalar) -> bool:
        return self.compare(other) <= 0

    def __gt__(self, other: _Scalar) -> bool:
        return self.compare(other) > 0

    def __ge__(self, other: _Scalar) -> bool:
        return self.compare(other) >= 0

    def __float__(self) -> float:
        """Nearest-ish binary64 (reporting only, never a decision)."""
        return float(self.rational) + float(self.coefficient) * math.sqrt(float(self.radicand))

    def __repr__(self) -> str:
        if self.is_rational:
            return f"OneRootNumber({self.rational})"
        return f"OneRootNumber({self.rational} + {self.coefficient}·√{self.radicand})"


def _lift(value: _Scalar) -> OneRootNumber:
    if isinstance(value, OneRootNumber):
        return value
    return OneRootNumber.of(Fraction(value))


def _sign_of_sum(u: OneRootNumber, w: OneRootNumber, beta: Fraction) -> int:
    """Exact sign of ``u + w·√β`` with ``u, w`` in one extension ``ℚ(√α)`` and rational ``β >= 0``."""
    sign_u, sign_w = u.sign(), w.sign()
    if beta == 0 or sign_w == 0:
        return sign_u
    if sign_u == 0 or sign_u == sign_w:
        return sign_w if sign_u == 0 else sign_u
    return sign_u * (u * u - w * w * beta).sign()


@dataclass(frozen=True)
class _TwoRootLinear:
    """``A + B√α + C√β + D√α√β``: the products of numbers from two extensions, for signs only."""

    a: Fraction
    b: Fraction
    c: Fraction
    d: Fraction
    alpha: Fraction
    beta: Fraction

    @classmethod
    def product(cls, left: OneRootNumber, right: OneRootNumber) -> _TwoRootLinear:
        return cls(
            left.rational * right.rational,
            left.coefficient * right.rational,
            left.rational * right.coefficient,
            left.coefficient * right.coefficient,
            left.radicand,
            right.radicand,
        )

    def plus(self, other: _TwoRootLinear, scale: int = 1) -> _TwoRootLinear:
        """``self + scale·other``; both must use the same two extensions (or a rational one)."""
        alpha = self.alpha or other.alpha
        beta = self.beta or other.beta
        if other.alpha not in (0, alpha) or other.beta not in (0, beta):
            raise CrossRootArithmeticError(f"(√{self.alpha}, √{self.beta}) and (√{other.alpha}, √{other.beta}) are different extension pairs.")
        return _TwoRootLinear(self.a + scale * other.a, self.b + scale * other.b, self.c + scale * other.c, self.d + scale * other.d, alpha, beta)

    def sign(self) -> int:
        """Exact sign by grouping ``(A + B√α) + √β·(C + D√α)``."""
        u = OneRootNumber.build(self.a, self.b, self.alpha)
        w = OneRootNumber.build(self.c, self.d, self.alpha)
        return _sign_of_sum(u, w, self.beta)


# ---------------------------------------------------------------------------
# One-root points
# ---------------------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class OneRootPoint2:
    """A point whose coordinates lie in one quadratic extension ``ℚ(√δ)`` (world XY, mm)."""

    x: OneRootNumber
    y: OneRootNumber

    __hash__ = None  # type: ignore[assignment]

    @classmethod
    def build(cls, x: OneRootNumber, y: OneRootNumber) -> Self:
        """Validate that both coordinates share one extension.

        Raises:
            MixedRootPointError: If the coordinates lie in different extensions.
        """
        if x.radicand != 0 and y.radicand != 0 and x.radicand != y.radicand:
            raise MixedRootPointError(f"Coordinates in √{x.radicand} and √{y.radicand}.")
        return cls(x, y)

    @classmethod
    def from_exact(cls, point: ExactPoint2[WorldXY]) -> Self:
        """Embed a rational point."""
        return cls(OneRootNumber.of(point.x), OneRootNumber.of(point.y))

    @classmethod
    def rational(cls, x: Fraction | int, y: Fraction | int) -> Self:
        """A rational point from two rationals."""
        return cls(OneRootNumber.of(x), OneRootNumber.of(y))

    @property
    def radicand(self) -> Fraction:
        """The shared radicand (``0`` for a rational point)."""
        return self.x.radicand if self.x.radicand != 0 else self.y.radicand

    def minus(self, origin: ExactPoint2[WorldXY]) -> tuple[OneRootNumber, OneRootNumber]:
        """The exact vector ``self - origin`` (same extension)."""
        return self.x - origin.x, self.y - origin.y

    def compare_xy(self, other: OneRootPoint2) -> int:
        """Exact lexicographic (x, then y) comparison, across extensions."""
        by_x = self.x.compare(other.x)
        return by_x if by_x != 0 else self.y.compare(other.y)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OneRootPoint2):
            return NotImplemented
        return self.compare_xy(other) == 0

    def as_floats(self) -> tuple[float, float]:
        """Binary64 view (reporting only)."""
        return float(self.x), float(self.y)

    def __repr__(self) -> str:
        return f"OneRootPoint2({self.x!r}, {self.y!r})"


def orientation(centre: ExactPoint2[WorldXY], first: OneRootPoint2, second: OneRootPoint2) -> int:
    """Exact sign of ``cross(first - centre, second - centre)``, across extensions."""
    first_x, first_y = first.minus(centre)
    second_x, second_y = second.minus(centre)
    return _TwoRootLinear.product(first_x, second_y).plus(_TwoRootLinear.product(first_y, second_x), -1).sign()


def _dot_sign(centre: ExactPoint2[WorldXY], first: OneRootPoint2, second: OneRootPoint2) -> int:
    first_x, first_y = first.minus(centre)
    second_x, second_y = second.minus(centre)
    return _TwoRootLinear.product(first_x, second_x).plus(_TwoRootLinear.product(first_y, second_y)).sign()


# ---------------------------------------------------------------------------
# Order on a circle
# ---------------------------------------------------------------------------


def compare_on_arc(centre: ExactPoint2[WorldXY], sense: Sense, start: OneRootPoint2, first: OneRootPoint2, second: OneRootPoint2) -> int:
    """Exact order of two points of one circle by swept angle from ``start`` in ``sense``."""
    first_half, second_half = _half_turn(centre, sense, start, first), _half_turn(centre, sense, start, second)
    if first_half != second_half:
        return -1 if first_half < second_half else 1
    return -sense * orientation(centre, first, second)


def _half_turn(centre: ExactPoint2[WorldXY], sense: Sense, start: OneRootPoint2, point: OneRootPoint2) -> int:
    """``0`` for swept angle in ``[0, π)``, ``1`` for ``[π, 2π)``."""
    side = sense * orientation(centre, start, point)
    if side > 0:
        return 0
    if side < 0:
        return 1
    return 0 if _dot_sign(centre, start, point) > 0 else 1
