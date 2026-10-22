"""The single-sweep engagement bound decided over a sweep's realised squared-distance range, any sweep shape.

Doctrine and proof: ``docs/held_connector_engagement.md`` (Lemma, and Lemma
for arc motion). The bound of one prior annulus ``k`` at a tool centre ``p``
depends on ``p`` only through ``D = |p - c_k|²``: the outer arc is
non-decreasing in ``D`` and the hole arc is largest at the smallest
``kappa_h``, which is convex in ``d`` with its minimum at ``d*² = r² - h²``.
The tool centre moves continuously over a connected sweep, so ``D`` takes
every value of ``[D_min, D_max]`` and no other. A sweep type therefore only
has to name a finite set of REALISED values of ``D`` that contains
``D_min`` and ``D_max`` (:meth:`Sweep.squared_distances`); the decision over
that set is exact and is the same for lines and arcs:

- the outer test at every candidate (it is monotone, so this is the test at
  ``D_max``);
- the hole test at every candidate (realised values, so no candidate can
  refuse what the true range would certify), and at ``d*`` when some
  candidate lies below ``d*²`` and some above.

Candidates are one-root numbers ``X / Y`` (``Y > 0``) whose numerator and
denominator share one quadratic extension, so every test is a sign of a
same-extension one-root number (:class:`benchmarks.held_one_root.OneRootNumber`).
Candidates of one sweep may live in different extensions; they are never
combined arithmetically, only each tested on its own. Every sign is first
asked of the repo's rigorous dyadic interval enclosure
(:class:`benchmarks.held_dyadic_interval.Interval`, outward-rounded integer
ends): it decides only when the enclosure excludes zero, otherwise the exact
one-root sign decides (:func:`filtered_sign`). The filter never decides
wrongly; it only declines. Degrees are reporting only
(:meth:`RangeEngagementLemma.bound_deg`).
"""

from __future__ import annotations

import enum
import math
from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Protocol
from typing import Self

from benchmarks.held_connector_engagement import CHORD_RATIO_MAX
from benchmarks.held_connector_engagement import ConnectorEngagementError
from benchmarks.held_connector_engagement import InvalidConnectorEngagementInputError
from benchmarks.held_connector_engagement import hole_arc_deg
from benchmarks.held_connector_engagement import outer_arc_deg
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


class ExactSweepError(ConnectorEngagementError):
    """Base failure of an exact connector sweep's contract."""


class DegenerateSweepError(ExactSweepError):
    """A line sweep has equal endpoints, or an arc sweep a non-positive radius."""


class MixedRootLineSweepError(ExactSweepError):
    """A line sweep's endpoints lie in different extensions and no rational supporting line was given."""


class PointOffSweepSupportError(ExactSweepError):
    """A sweep endpoint is not exactly on the declared supporting line or circle."""


class DiscontinuousConnectorError(ExactSweepError):
    """Consecutive sweeps of one connector do not share their endpoint exactly."""


class SweepKind(enum.Enum):
    """What a sweep is in the path, for refusal reporting; the decision does not read it."""

    LINE = "line"
    ARC = "arc"
    OFFSET_LINE = "offset_line"
    OFFSET_ARC = "offset_arc"
    JUNCTION_ARC = "junction_arc"
    CIRCLE_JUNCTION = "circle_junction"


def filtered_sign(enclosure: Interval, exact: Callable[[], int]) -> int:
    """The sign of a value: its interval enclosure's when that excludes zero, else the exact sign."""
    sign = enclosure.sign()
    return exact() if sign is None else sign


@dataclass(frozen=True)
class PointBox:
    """Interval enclosures of an exact point's coordinates (filter only)."""

    x: Interval
    y: Interval

    @classmethod
    def of(cls, point: OneRootPoint2) -> Self:
        return cls(Interval.of_one_root(point.x), Interval.of_one_root(point.y))

    @classmethod
    def of_exact(cls, point: ExactPoint2[WorldXY]) -> Self:
        return cls(Interval.of(point.x), Interval.of(point.y))


def points_apart(first: OneRootPoint2, second: OneRootPoint2) -> bool:
    """Exact ``first != second``: disjoint interval enclosures decide, else exact comparison."""
    one, other = PointBox.of(first), PointBox.of(second)
    if (one.x - other.x).sign() is not None or (one.y - other.y).sign() is not None:
        return True
    return first != second


class SquaredDistance:
    """A realised squared distance ``X / Y`` with ``Y > 0``, ``X`` and ``Y`` in one extension.

    Carries interval enclosures of ``X`` and ``Y`` (the filter) and computes
    the exact one-root pair only when a filter declines.
    """

    __slots__ = ("_exact", "_pair", "denominator_box", "numerator_box")

    def __init__(self, numerator_box: Interval, denominator_box: Interval, exact: Callable[[], tuple[OneRootNumber, OneRootNumber]]) -> None:
        self.numerator_box = numerator_box
        self.denominator_box = denominator_box
        self._exact = exact
        self._pair: tuple[OneRootNumber, OneRootNumber] | None = None

    @classmethod
    def of(cls, value: OneRootNumber) -> SquaredDistance:
        """``value / 1``."""
        return cls(Interval.of_one_root(value), Interval.of(1), lambda: (value, OneRootNumber.of(1)))

    @property
    def numerator(self) -> OneRootNumber:
        return self.pair[0]

    @property
    def denominator(self) -> OneRootNumber:
        return self.pair[1]

    @property
    def pair(self) -> tuple[OneRootNumber, OneRootNumber]:
        """The exact ``(X, Y)``, computed once."""
        if self._pair is None:
            self._pair = self._exact()
        return self._pair

    def compare(self, value: Fraction) -> int:
        """Exact sign of ``X / Y - value``."""
        return filtered_sign(self.numerator_box - Interval.of(value) * self.denominator_box, lambda: (self.numerator - value * self.denominator).sign())

    def reported(self) -> float:
        """Binary64 view (reporting only)."""
        return float(self.numerator) / float(self.denominator)


@dataclass(frozen=True)
class PriorCircle:
    """One cut circle as the lemma reads it: exact centre and exact guide radius ``rho > 0``."""

    centre: ExactPoint2[WorldXY]
    guide_radius: Fraction
    box: PointBox = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "box", PointBox.of_exact(self.centre))

    @classmethod
    def from_binary64(cls, centre: Point2[WorldXY], guide_radius: float) -> Self:
        """Inject a binary64 circle exactly."""
        return cls(ExactPoint2[WorldXY].from_binary64(float(centre.x), float(centre.y)), Fraction(guide_radius))


class Sweep(Protocol):
    """A connected motion of the tool centre with exact endpoints."""

    @property
    def start(self) -> OneRootPoint2: ...

    @property
    def end(self) -> OneRootPoint2: ...

    @property
    def kind(self) -> SweepKind: ...

    def squared_distances(self, circle: PriorCircle) -> tuple[SquaredDistance, ...]:
        """Realised values of ``|p - c|²`` over the sweep (``c`` the circle's centre), including its minimum and maximum."""
        ...

    def split(self) -> tuple[Self, Self] | None:
        """Two sweeps of the same type whose union is this sweep, or ``None`` when no exact split point is found."""
        ...


def squared_distance(point: OneRootPoint2, box: PointBox, circle: PriorCircle) -> SquaredDistance:
    """``|point - c|²`` (in ``point``'s extension), with ``box`` enclosing ``point``."""

    def exact() -> tuple[OneRootNumber, OneRootNumber]:
        dx, dy = point.minus(circle.centre)
        return dx * dx + dy * dy, OneRootNumber.of(1)

    return SquaredDistance((box.x - circle.box.x).square() + (box.y - circle.box.y).square(), Interval.of(1), exact)


@dataclass(frozen=True)
class RangeEngagementLemma:
    """The single-sweep bound with its exact decision over realised squared distances.

    ``tool_radius`` is the exact tool radius ``r``; ``cos_sq_half_cap`` is
    ``C = 1 - s / 4`` for the cap's binary64 chord surrogate ``s``.
    """

    tool_radius: Fraction
    cos_sq_half_cap: Fraction

    @classmethod
    def build(cls, tool_radius: Fraction, chord_ratio: float) -> Self:
        """Inject an exact tool radius and a binary64 chord surrogate ``s``.

        Raises:
            InvalidConnectorEngagementInputError: ``s`` is not a finite binary64
                in ``(0, 4]`` or ``r`` is not positive.
        """
        if not math.isfinite(chord_ratio) or not 0.0 < chord_ratio <= CHORD_RATIO_MAX:
            raise InvalidConnectorEngagementInputError(f"chord surrogate must lie in (0, 4], got {chord_ratio!r}")
        if tool_radius <= 0:
            raise InvalidConnectorEngagementInputError(f"tool radius must be positive, got {tool_radius}")
        return cls(Fraction(tool_radius), 1 - Fraction(chord_ratio) / 4)

    def outer_within(self, distance: SquaredDistance, outer_radius: Fraction) -> bool:
        """Outer arc within the cap at ``D = X / Y``: ``N_O <= 0`` and ``N_O² >= 4 C D r²`` (times ``Y²``)."""
        return self._arc_within(distance, self.tool_radius * self.tool_radius - outer_radius * outer_radius, outer=True)

    def _arc_within(self, distance: SquaredDistance, offset: Fraction, *, outer: bool) -> bool:
        """``N = X + offset·Y`` on the arc's side of zero and ``N² >= 4 C X Y r²``, each sign filtered then exact."""
        c = self.cos_sq_half_cap
        r_sq = self.tool_radius * self.tool_radius
        x_box, y_box = distance.numerator_box, distance.denominator_box
        n_box = x_box + Interval.of(offset) * y_box

        def exact_n() -> int:
            return (distance.numerator + offset * distance.denominator).sign()

        n_sign = filtered_sign(n_box, exact_n)
        if (outer and n_sign > 0) or (not outer and n_sign < 0):
            return False
        weight = 4 * c.numerator * r_sq

        def exact_square() -> int:
            x, y = distance.pair
            n = x + offset * y
            return (n * n * c.denominator - weight * x * y).sign()

        square_box = n_box.square() * Interval.of(c.denominator) - Interval.of(weight) * x_box * y_box
        return filtered_sign(square_box, exact_square) >= 0

    def hole_within(self, distance: SquaredDistance, inner_radius: Fraction) -> bool:
        """Hole arc within the cap at ``D = X / Y > 0``; ``D = 0`` is refused (no certificate, never a verdict)."""
        if filtered_sign(distance.numerator_box, lambda: distance.numerator.sign()) == 0:
            return False
        return self._arc_within(distance, self.tool_radius * self.tool_radius - inner_radius * inner_radius, outer=False)

    def range_within(self, distances: tuple[SquaredDistance, ...], circle: PriorCircle) -> bool:
        """Exact: every engaged run over a sweep whose realised ``D`` set contains ``distances``' min and max is within the cap.

        ``distances`` must be realised values that include the minimum and the
        maximum of ``D`` over the sweep (the sweep type's contract).
        """
        r, rho = self.tool_radius, circle.guide_radius
        if not all(self.outer_within(distance, rho + r) for distance in distances):
            return False
        if rho <= r:
            return True
        if all(distance.compare(rho * rho) >= 0 for distance in distances):
            return True
        inner = rho - r
        if not all(self.hole_within(distance, inner) for distance in distances):
            return False
        star = r * r - inner * inner
        if star > 0 and any(distance.compare(star) < 0 for distance in distances) and any(distance.compare(star) > 0 for distance in distances):
            c = self.cos_sq_half_cap
            return star * c.denominator >= c.numerator * r * r
        return True

    def certifies(self, sweep: Sweep, circle: PriorCircle) -> bool:
        """Exact: ``circle``'s annulus alone bounds every engaged run along ``sweep`` within the cap."""
        return self.range_within(sweep.squared_distances(circle), circle)

    def bound_deg(self, sweep: Sweep, circle: PriorCircle) -> float:
        """REPORTING ONLY: the lemma's bound over the sweep, degrees."""
        r, rho = float(self.tool_radius), float(circle.guide_radius)
        values = [distance.reported() for distance in sweep.squared_distances(circle)]
        d_min, d_max = math.sqrt(max(min(values), 0.0)), math.sqrt(max(values))
        bound = outer_arc_deg(d_max, r, rho + r)
        if rho > r and d_min < rho:
            h = rho - r
            distances = [d_min, d_max]
            if r > h and d_min * d_min < r * r - h * h < d_max * d_max:
                distances.append(math.sqrt(r * r - h * h))
            bound = max(bound, *(hole_arc_deg(d, r, h) for d in distances))
        return bound
