"""Independent exact reference for the Section 2.3 standard-model maximum engagement.

Test oracle only. It shares no code with ``exact_predecessor_engagement``: the
model is re-derived from its geometric construction in an explicit local frame
and evaluated over ``Fraction`` coordinates in the field Q(sqrt(D)), with
divisions taken in the field instead of multiplied through. The Fig. 4(d)
correction constructs the swept-disk intersection point I and the corrected
tool centre q explicitly and takes the half-angle from a dot product, where the
code under test expands a law-of-cosines identity over integers.

Frame: origin at the candidate centre c, +x from the predecessor centre p to
c, so p = (-d, 0) with d = |c - p| = sqrt(D). R_p = rho_prev + r is the
predecessor swept radius and R_c = rho + r the candidate's.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import PredecessorEngagementCosine
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius


class ReferenceFragmentationError(ArithmeticError):
    """The corrected construction has no real intersection point or chord direction."""


class ReferenceCosineRangeError(ArithmeticError):
    """A law-of-cosines branch produced a cosine outside [-1, 1]: no real tool-circle intersection."""


class ReferenceRadicandError(ArithmeticError):
    """Field arithmetic across two different radicands."""


def _rational_sqrt(value: Fraction) -> Fraction | None:
    """Exact square root of a non-negative rational, or None when it is irrational."""
    numerator, denominator = math.isqrt(value.numerator), math.isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Fraction(numerator, denominator)
    return None


@dataclass(frozen=True)
class QuadraticNumber:
    """``a + b sqrt(D)`` with rational a, b; D > 0 is never a rational square (collapsed to Q)."""

    a: Fraction
    b: Fraction
    radicand: Fraction

    @staticmethod
    def field(radicand: Fraction) -> tuple[QuadraticNumber, QuadraticNumber]:
        """Return ``(1, sqrt(radicand))`` of the field, collapsing a rational root into Q."""
        root = _rational_sqrt(radicand)
        if root is not None:
            return QuadraticNumber(Fraction(1), Fraction(0), Fraction(0)), QuadraticNumber(root, Fraction(0), Fraction(0))
        return QuadraticNumber(Fraction(1), Fraction(0), radicand), QuadraticNumber(Fraction(0), Fraction(1), radicand)

    def _lift(self, other: QuadraticNumber | Fraction | int) -> QuadraticNumber:
        if isinstance(other, QuadraticNumber):
            if other.radicand != self.radicand:
                raise ReferenceRadicandError("Operands live in different quadratic fields.")
            return other
        return QuadraticNumber(Fraction(other), Fraction(0), self.radicand)

    def __add__(self, other: QuadraticNumber | Fraction | int) -> QuadraticNumber:
        o = self._lift(other)
        return QuadraticNumber(self.a + o.a, self.b + o.b, self.radicand)

    __radd__ = __add__

    def __neg__(self) -> QuadraticNumber:
        return QuadraticNumber(-self.a, -self.b, self.radicand)

    def __sub__(self, other: QuadraticNumber | Fraction | int) -> QuadraticNumber:
        return self + -self._lift(other)

    def __rsub__(self, other: Fraction | int) -> QuadraticNumber:
        return self._lift(other) - self

    def __mul__(self, other: QuadraticNumber | Fraction | int) -> QuadraticNumber:
        o = self._lift(other)
        return QuadraticNumber(self.a * o.a + self.b * o.b * self.radicand, self.a * o.b + self.b * o.a, self.radicand)

    __rmul__ = __mul__

    def __truediv__(self, other: QuadraticNumber | Fraction | int) -> QuadraticNumber:
        o = self._lift(other)
        # D is not a rational square, so the norm vanishes only for zero.
        norm = o.a * o.a - o.b * o.b * self.radicand
        if norm == 0:
            raise ZeroDivisionError("Division by zero in Q(sqrt(D)).")
        return self * QuadraticNumber(o.a / norm, -o.b / norm, self.radicand)

    def sign(self) -> int:
        """Exact sign: like-signed parts decide; unlike-signed parts compare a^2 with b^2 D."""
        sa, sb = (self.a > 0) - (self.a < 0), (self.b > 0) - (self.b < 0)
        if sb == 0:
            return sa
        if sa == 0 or sa == sb:
            return sb
        difference = self.a * self.a - self.b * self.b * self.radicand
        return sa if difference > 0 else -sa if difference < 0 else 0


@dataclass(frozen=True)
class ReferenceEngagement:
    """Exact cos(theta) of the standard-model maximum engagement, in Q(sqrt(D))."""

    cosine: QuadraticNumber
    squared_spacing: Fraction
    branch: str  # model branch that produced the cosine, for coverage events

    def within(self, cap: EngagementCap) -> bool:
        """theta <= cap on the EngagementCap surrogate: cos theta >= 1 - chord_ratio / 2."""
        return (self.cosine - (1 - Fraction(float(cap.chord_ratio)) / 2)).sign() >= 0

    def reaches(self, floor: EngagementCap) -> bool:
        """theta >= floor on the surrogate: cos theta <= 1 - chord_ratio / 2."""
        return (self.cosine - (1 - Fraction(float(floor.chord_ratio)) / 2)).sign() <= 0

    def equals(self, cosine: PredecessorEngagementCosine) -> bool:
        """Exact value equality with the code-under-test's integer-scaled representation.

        Its radicand is D times the square of the common power-of-two scale,
        so sqrt(code radicand) = k sqrt(D) with k = sqrt(code radicand / D) rational.
        """
        numerator, denominator = cosine.numerator, Fraction(cosine.denominator)
        value = self.cosine - Fraction(numerator.rational) / denominator
        if numerator.radical != 0 and numerator.radicand != 0:
            if self.squared_spacing == 0:
                return False
            scale = _rational_sqrt(Fraction(numerator.radicand) / self.squared_spacing)
            if scale is None:
                return False
            _, root = QuadraticNumber.field(self.squared_spacing)
            value = value - root * (Fraction(numerator.radical) * scale / denominator)
        return value.sign() == 0


def _injected(candidate: PaperCircleCandidate) -> tuple[Fraction, Fraction, Fraction]:
    return Fraction(float(candidate.center.x)), Fraction(float(candidate.center.y)), Fraction(float(candidate.guide_radius.value))


def _law_of_cosines(cosine: QuadraticNumber, squared_spacing: Fraction, branch: str) -> ReferenceEngagement:
    if (cosine - 1).sign() > 0 or (cosine + 1).sign() < 0:
        raise ReferenceCosineRangeError("Tool circle does not reach the bounding circle.")
    return ReferenceEngagement(cosine, squared_spacing, branch)


def reference_engagement(predecessor: PaperCircleCandidate, candidate: PaperCircleCandidate, tool_radius: ToolRadius) -> ReferenceEngagement:
    """Standard-model maximum engagement, derived geometrically in the local frame.

    Raises:
        ReferenceFragmentationError: Fig. 4(d) has no real I or no chord direction.
        ReferenceCosineRangeError: A law-of-cosines branch has no real intersection.
    """
    px, py, rho_prev = _injected(predecessor)
    cx, cy, rho = _injected(candidate)
    r = Fraction(float(tool_radius.value))
    squared_spacing = (cx - px) ** 2 + (cy - py) ** 2
    one, d = QuadraticNumber.field(squared_spacing)
    outer_prev, outer = rho_prev + r, rho + r
    # Candidate guide disk inside the predecessor's: nothing new is swept.
    if rho_prev >= rho and squared_spacing <= (rho_prev - rho) ** 2:
        return ReferenceEngagement(one, squared_spacing, "contained")

    def full_slot(branch: str) -> ReferenceEngagement:
        return ReferenceEngagement(-one, squared_spacing, branch)

    if squared_spacing == 0:
        if rho - r >= outer_prev:
            return full_slot("concentric full slot")
        # Tool at T = (rho, 0); its circle meets |X| = R_p at X_x; angle from +x at T.
        crossing_x = (outer_prev * outer_prev - r * r + rho * rho) / (2 * rho)
        return _law_of_cosines(one * ((crossing_x - rho) / r), squared_spacing, "concentric")
    # Eq. 4: the swept disks do not overlap.
    if (2 * r - rho + rho_prev - d).sign() <= 0:
        return full_slot("no overlap")
    # B = (b, 0): where the predecessor swept circle crosses the +x axis beyond c.
    b = outer_prev - d
    if b.sign() <= 0:
        return full_slot("b <= 0")
    # Tool centre q on the candidate circle whose tool circle passes through B.
    q_x = (b * b + rho * rho - r * r) / (2 * b)
    q_y_squared = rho * rho - q_x * q_x
    if q_y_squared.sign() < 0:
        return full_slot("no contact station")
    # Outermost tool point W = (R_c / rho) q; is it outside the predecessor swept disk?
    w_x, w_y_squared = q_x * (outer / rho), q_y_squared * (outer / rho) ** 2
    if ((w_x + d) * (w_x + d) + w_y_squared - outer_prev * outer_prev).sign() > 0:
        # Eq. 7: tool at T = (rho, 0); the engaged arc ends where its circle reaches |X| = b.
        boundary_crossing_x = (b * b + rho * rho - r * r) / (2 * rho)
        return _law_of_cosines((boundary_crossing_x - rho) / r, squared_spacing, "equation 7")
    # Fig. 4(d): I = (x_I, sqrt(beta)) on both swept circles; q = (rho / R_c) I.
    x_i = one * (outer_prev * outer_prev - outer * outer - squared_spacing) / (2 * d)
    beta = outer * outer - x_i * x_i
    if beta.sign() < 0:
        raise ReferenceFragmentationError("Swept circles do not intersect.")
    scale = rho / outer
    # Vectors as (x, y / sqrt(beta)); dot(u, v) = u_x v_x + beta u_y v_y.
    to_p = (-d - x_i * scale, -scale * one)
    to_i = (x_i * (1 - scale), (1 - scale) * one)
    dot = to_p[0] * to_i[0] + beta * to_p[1] * to_i[1]
    p_squared = to_p[0] * to_p[0] + beta * to_p[1] * to_p[1]
    i_squared = to_i[0] * to_i[0] + beta * to_i[1] * to_i[1]
    if p_squared.sign() == 0:
        raise ReferenceFragmentationError("Corrected tool centre coincides with the predecessor centre.")
    half_cosine_squared = dot * dot / (p_squared * i_squared)
    return ReferenceEngagement(2 * half_cosine_squared - 1, squared_spacing, "figure 4(d)")
