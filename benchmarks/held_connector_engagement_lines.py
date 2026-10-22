"""Exact straight sweeps for the connector-engagement certificate: rational or one-root endpoints, no rounding.

Along a straight sweep ``[a, b]`` the squared distance to a centre ``c`` is a
convex quadratic in the parameter, so its maximum is at an end and its
minimum at an end or at the foot of ``c`` when the foot is strictly inside.
Those are the realised candidates (:mod:`benchmarks.held_connector_engagement_range`).

Two exact representations, chosen by the endpoints:

- **one extension** (both ends rational, or both in one ``Q(√α)``): the foot
  test and ``cross² / |b - a|²`` stay in that extension;
- **two extensions** (ends in ``Q(√α)`` and ``Q(√β)``, ``α != β``): the
  sweep must declare its rational supporting line (anchor, direction), as a
  ∂Ω* offset line does. The foot parameter and the line distance are
  rational; each end's parameter is compared with the rational foot
  parameter in its own extension. Without a support the sweep is refused at
  build (:class:`MixedRootLineSweepError`): cross-root arithmetic is never
  written.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Self

from benchmarks.held_connector_engagement_range import DegenerateSweepError
from benchmarks.held_connector_engagement_range import MixedRootLineSweepError
from benchmarks.held_connector_engagement_range import PointBox
from benchmarks.held_connector_engagement_range import PointOffSweepSupportError
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import SquaredDistance
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_connector_engagement_range import filtered_sign
from benchmarks.held_connector_engagement_range import points_apart
from benchmarks.held_connector_engagement_range import squared_distance
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import WorldXY

_HALF = Fraction(1, 2)


@dataclass(frozen=True)
class RationalLine:
    """The rational line ``anchor + t·direction`` (``direction != 0``)."""

    anchor: ExactPoint2[WorldXY]
    direction: tuple[Fraction, Fraction]

    @classmethod
    def build(cls, anchor: ExactPoint2[WorldXY], direction: tuple[Fraction, Fraction]) -> Self:
        """Validate a non-zero direction.

        Raises:
            DegenerateSweepError: ``direction`` is the zero vector.
        """
        if direction[0] == 0 and direction[1] == 0:
            raise DegenerateSweepError("a supporting line needs a non-zero direction")
        return cls(anchor, (Fraction(direction[0]), Fraction(direction[1])))

    def parameter(self, point: OneRootPoint2) -> OneRootNumber:
        """``(point - anchor)·direction`` (``point``'s extension; not divided by ``|direction|²``)."""
        dx, dy = point.minus(self.anchor)
        return self.direction[0] * dx + self.direction[1] * dy

    def contains(self, point: OneRootPoint2) -> bool:
        """Whether ``point`` is exactly on the line."""
        dx, dy = point.minus(self.anchor)
        return (dx * self.direction[1] - dy * self.direction[0]).sign() == 0


def _common(first: OneRootPoint2, second: OneRootPoint2) -> bool:
    return first.radicand == 0 or second.radicand == 0 or first.radicand == second.radicand


@dataclass(frozen=True)
class LineSweep:
    """The tool centre moving straight from ``start`` to ``end`` (exact; ``support`` required across extensions)."""

    start: OneRootPoint2
    end: OneRootPoint2
    support: RationalLine | None = None
    kind: SweepKind = SweepKind.LINE
    start_box: PointBox = field(init=False, repr=False, compare=False)
    end_box: PointBox = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_box", PointBox.of(self.start))
        object.__setattr__(self, "end_box", PointBox.of(self.end))

    @classmethod
    def build(cls, start: OneRootPoint2, end: OneRootPoint2, support: RationalLine | None = None, kind: SweepKind = SweepKind.LINE) -> Self:
        """Validate a non-degenerate straight sweep.

        Raises:
            DegenerateSweepError: ``start == end``.
            MixedRootLineSweepError: the ends lie in different extensions and
                no ``support`` is given.
            PointOffSweepSupportError: an end is not exactly on ``support``.
        """
        if not points_apart(start, end):
            raise DegenerateSweepError(f"a line sweep needs distinct ends, got {start!r} twice")
        if support is not None:
            for point in (start, end):
                if not support.contains(point):
                    raise PointOffSweepSupportError(f"{point!r} is not on its supporting line")
        elif not _common(start, end):
            raise MixedRootLineSweepError(f"{start!r} and {end!r} lie in different extensions and no rational supporting line was given")
        return cls(start, end, support, kind)

    def squared_distances(self, circle: PriorCircle) -> tuple[SquaredDistance, ...]:
        """Both ends, and the foot of the circle's centre when it is strictly inside."""
        ends = (squared_distance(self.start, self.start_box, circle), squared_distance(self.end, self.end_box, circle))
        foot = self._foot_in_extension(circle) if _common(self.start, self.end) else self._foot_on_support(self._required_support(), circle)
        return ends if foot is None else (*ends, foot)

    def _required_support(self) -> RationalLine:
        if self.support is None:
            raise MixedRootLineSweepError(f"{self.start!r} and {self.end!r} lie in different extensions and no rational supporting line was given")
        return self.support

    def _foot_in_extension(self, circle: PriorCircle) -> SquaredDistance | None:
        centre, box = circle.centre, circle.box
        ex_box, ey_box = self.end_box.x - self.start_box.x, self.end_box.y - self.start_box.y
        wx_box, wy_box = box.x - self.start_box.x, box.y - self.start_box.y
        dot_box = wx_box * ex_box + wy_box * ey_box

        def vectors() -> tuple[OneRootNumber, OneRootNumber, OneRootNumber, OneRootNumber]:
            return self.end.x - self.start.x, self.end.y - self.start.y, centre.x - self.start.x, centre.y - self.start.y

        def exact_dot() -> int:
            ex, ey, wx, wy = vectors()
            return (wx * ex + wy * ey).sign()

        if filtered_sign(dot_box, exact_dot) <= 0:
            return None
        length_box = ex_box.square() + ey_box.square()

        def exact_beyond() -> int:
            ex, ey, wx, wy = vectors()
            return (ex * ex + ey * ey - (wx * ex + wy * ey)).sign()

        if filtered_sign(length_box - dot_box, exact_beyond) <= 0:
            return None

        def exact_foot() -> tuple[OneRootNumber, OneRootNumber]:
            ex, ey, wx, wy = vectors()
            cross = wx * ey - wy * ex
            return cross * cross, ex * ex + ey * ey

        return SquaredDistance((wx_box * ey_box - wy_box * ex_box).square(), length_box, exact_foot)

    def _foot_on_support(self, support: RationalLine, circle: PriorCircle) -> SquaredDistance | None:
        ux, uy = support.direction
        wx, wy = circle.centre.minus(support.anchor)
        foot = ux * wx + uy * wy
        anchor_x, anchor_y, ux_box, uy_box, foot_box = Interval.of(support.anchor.x), Interval.of(support.anchor.y), Interval.of(ux), Interval.of(uy), Interval.of(foot)

        def side(point: OneRootPoint2, box: PointBox) -> int:
            parameter_box = ux_box * (box.x - anchor_x) + uy_box * (box.y - anchor_y)
            return filtered_sign(parameter_box - foot_box, lambda: (support.parameter(point) - foot).sign())

        if side(self.start, self.start_box) * side(self.end, self.end_box) >= 0:
            return None
        length_sq = ux * ux + uy * uy
        distance = (wx * wx + wy * wy) * length_sq - foot * foot
        return SquaredDistance(Interval.of(distance), Interval.of(length_sq), lambda: (OneRootNumber.of(distance), OneRootNumber.of(length_sq)))

    def split(self) -> tuple[LineSweep, LineSweep] | None:
        """Halve at the exact midpoint (one extension), or at a rational support point strictly inside (two extensions).

        Across extensions the rational parameter is proposed in binary64 and
        accepted only when it lies strictly between the two ends exactly;
        otherwise ``None`` (the caller refuses the sweep by name).
        """
        if _common(self.start, self.end):
            middle = OneRootPoint2.build(_HALF * (self.start.x + self.end.x), _HALF * (self.start.y + self.end.y))
        else:
            support = self._required_support()
            low, high = support.parameter(self.start), support.parameter(self.end)
            proposal = Fraction((float(low) + float(high)) / 2)
            if (low - proposal).sign() * (high - proposal).sign() >= 0:
                return None
            ux, uy = support.direction
            scale = proposal / (ux * ux + uy * uy)
            middle = OneRootPoint2.rational(support.anchor.x + scale * ux, support.anchor.y + scale * uy)
        return LineSweep(self.start, middle, self.support, self.kind), LineSweep(middle, self.end, self.support, self.kind)
