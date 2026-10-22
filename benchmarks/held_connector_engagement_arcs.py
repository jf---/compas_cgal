"""Exact arc sweeps for the connector-engagement certificate: the tool centre on a rational circle of ∂Ω*.

Lemma for arc motion (proof in ``docs/held_connector_engagement.md``): on the
circle ``(q, R)`` with ``L = |c - q| > 0``, the squared distance to a centre
``c`` is ``D(φ) = R² + L² - 2RL·cos(φ - φ_c)``, where ``φ_c`` is the direction
of ``c - q``. Its only critical points are the NEAR point (direction
``c - q``, ``D = (L - R)²``) and the FAR point (direction ``q - c``,
``D = (L + R)²``). Over an arc from ``start`` to ``end`` the extremes of ``D``
are therefore among the two ends and whichever of NEAR / FAR lies on the arc;
for ``L = 0``, ``D = R²`` everywhere.

Exactness:

- ``D`` at an end lies in that end's extension; the two ends may be in
  different extensions (they are only tested one at a time).
- ``(L ± R)² = β + R² ± 2R√β`` with ``β = |c - q|²`` rational: one new
  extension ``Q(√β)``, its own candidate, never mixed with the ends.
- Whether NEAR / FAR lies on the arc is an order question on the circle about
  the rational directions ``±(c - q)``. The order (:meth:`ArcSweep.order`,
  the same rule as :func:`benchmarks.held_one_root.compare_on_arc`) reads
  only signs of cross and dot products about ``q``, so the rational points
  ``c`` and ``2q - c`` stand for NEAR and FAR (same directions). Each sign is
  asked of the interval filter first and otherwise decided by
  :func:`benchmarks.held_one_root.orientation`, the exact mixed-radical
  decomposition (a dot product is the orientation against the point turned a
  quarter about ``q``).
- Split points for tier 1b are RATIONAL points of the circle,
  ``q + R·u`` with ``u`` an exactly unit Pythagorean direction (half-angle
  tangent on the ``2^-56`` grid, :func:`benchmarks.held_declared_boundary_emission.choose_direction`).
  Binary64 only proposes the direction; that it lies strictly inside the arc
  is decided exactly, and a failed proposal refuses the split by name.
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Self

from benchmarks.held_connector_engagement_range import DegenerateSweepError
from benchmarks.held_connector_engagement_range import PointBox
from benchmarks.held_connector_engagement_range import PointOffSweepSupportError
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import SquaredDistance
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_connector_engagement_range import filtered_sign
from benchmarks.held_connector_engagement_range import points_apart
from benchmarks.held_connector_engagement_range import squared_distance
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import Sense
from benchmarks.held_declared_boundary_emission import choose_direction
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from benchmarks.held_one_root import orientation
from compas_cgal.units import WorldXY


def on_circle(centre: ExactPoint2[WorldXY], radius: Fraction, point: OneRootPoint2) -> bool:
    """Whether ``point`` is exactly on the circle ``(centre, radius)``."""
    dx, dy = point.minus(centre)
    return (dx * dx + dy * dy).compare(radius * radius) == 0


@dataclass(frozen=True)
class ArcSweep:
    """The tool centre moving on the circle ``(centre, radius)`` from ``start`` to ``end`` in ``sense``.

    ``start == end`` is a full turn (a zero-length arc is never a sweep).
    ``full`` is that exact decision, made once at construction
    (:func:`benchmarks.held_connector_engagement_range.points_apart`: interval
    filter, then the exact cross-root comparison).
    """

    centre: ExactPoint2[WorldXY]
    radius: Fraction
    sense: Sense
    start: OneRootPoint2
    end: OneRootPoint2
    kind: SweepKind = SweepKind.ARC
    centre_box: PointBox = field(init=False, repr=False, compare=False)
    start_box: PointBox = field(init=False, repr=False, compare=False)
    end_box: PointBox = field(init=False, repr=False, compare=False)
    full: bool = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "centre_box", PointBox.of_exact(self.centre))
        object.__setattr__(self, "start_box", PointBox.of(self.start))
        object.__setattr__(self, "end_box", PointBox.of(self.end))
        object.__setattr__(self, "full", not points_apart(self.start, self.end))

    @classmethod
    def build(cls, centre: ExactPoint2[WorldXY], radius: Fraction, sense: Sense, start: OneRootPoint2, end: OneRootPoint2, kind: SweepKind = SweepKind.ARC) -> Self:
        """Validate the radius and that both ends lie exactly on the circle.

        Raises:
            DegenerateSweepError: ``radius <= 0`` or ``sense`` is not ``±1``.
            PointOffSweepSupportError: an end is not exactly on the circle.
        """
        if radius <= 0 or sense not in (1, -1):
            raise DegenerateSweepError(f"an arc sweep needs a positive radius and sense ±1, got {radius}, {sense!r}")
        for point in (start, end):
            if not on_circle(centre, Fraction(radius), point):
                raise PointOffSweepSupportError(f"{point!r} is not on the circle about {centre!r} of radius {radius}")
        return cls(centre, Fraction(radius), sense, start, end, kind)

    def _orientation(self, first: OneRootPoint2, first_box: PointBox, second: OneRootPoint2, second_box: PointBox) -> int:
        q = self.centre_box
        box = (first_box.x - q.x) * (second_box.y - q.y) - (first_box.y - q.y) * (second_box.x - q.x)
        return filtered_sign(box, lambda: orientation(self.centre, first, second))

    def _dot(self, first: OneRootPoint2, first_box: PointBox, second: OneRootPoint2, second_box: PointBox) -> int:
        q = self.centre_box
        box = (first_box.x - q.x) * (second_box.x - q.x) + (first_box.y - q.y) * (second_box.y - q.y)

        def exact() -> int:
            dx, dy = second.minus(self.centre)
            turned = OneRootPoint2.build(self.centre.x - dy, self.centre.y + dx)
            return orientation(self.centre, first, turned)

        return filtered_sign(box, exact)

    def _half_turn(self, point: OneRootPoint2, box: PointBox) -> int:
        """``0`` for a swept angle from ``start`` in ``[0, π)``, ``1`` for ``[π, 2π)``."""
        side = self.sense * self._orientation(self.start, self.start_box, point, box)
        if side != 0:
            return 0 if side > 0 else 1
        return 0 if self._dot(self.start, self.start_box, point, box) > 0 else 1

    @functools.cached_property
    def _end_half_turn(self) -> int:
        """:meth:`_half_turn` of ``end``: exact, of immutable data, computed once."""
        return self._half_turn(self.end, self.end_box)

    def order(self, point: OneRootPoint2) -> int:
        """Exact order of the direction of ``point`` against ``end``, by swept angle from ``start`` in ``sense``.

        ``-1`` / ``0`` / ``1`` as ``point``'s direction comes before, with, or
        after ``end``'s; ``point`` need only be off the centre (the rule of
        :func:`benchmarks.held_one_root.compare_on_arc`, sign-filtered).
        """
        box = PointBox.of(point)
        first, second = self._half_turn(point, box), self._end_half_turn
        if first != second:
            return -1 if first < second else 1
        return -self.sense * self._orientation(point, box, self.end, self.end_box)

    def _reaches(self, direction: OneRootPoint2) -> bool:
        """Whether the ray from the centre through ``direction`` meets the arc (ends included)."""
        return self.full or self.order(direction) <= 0

    def squared_distances(self, circle: PriorCircle) -> tuple[SquaredDistance, ...]:
        """Both ends, and NEAR / FAR ``(L ∓ R)²`` when their directions lie on the arc."""
        centre = circle.centre
        candidates = [squared_distance(self.start, self.start_box, circle), squared_distance(self.end, self.end_box, circle)]
        vx, vy = centre.minus(self.centre)
        if vx == 0 and vy == 0:
            return tuple(candidates)
        beta, radius = vx * vx + vy * vy, self.radius
        if self._reaches(OneRootPoint2.rational(centre.x, centre.y)):
            candidates.append(SquaredDistance.of(OneRootNumber.build(beta + radius * radius, -2 * radius, beta)))
        if self._reaches(OneRootPoint2.rational(self.centre.x - vx, self.centre.y - vy)):
            candidates.append(SquaredDistance.of(OneRootNumber.build(beta + radius * radius, 2 * radius, beta)))
        return tuple(candidates)

    def split(self) -> tuple[ArcSweep, ArcSweep] | None:
        """Halve at a rational point of the circle strictly inside the arc, or ``None`` when the proposal is not inside."""
        cx, cy = float(self.centre.x), float(self.centre.y)
        sx, sy = self.start.as_floats()
        ex, ey = self.end.as_floats()
        begin = math.atan2(sy - cy, sx - cx)
        sweep = math.tau if self.full else (self.sense * (math.atan2(ey - cy, ex - cx) - begin)) % math.tau
        middle_angle = begin + self.sense * sweep / 2
        ux, uy = choose_direction(math.cos(middle_angle), math.sin(middle_angle)).vector
        middle = OneRootPoint2.rational(self.centre.x + self.radius * ux, self.centre.y + self.radius * uy)
        if middle == self.start or middle == self.end:
            return None
        if not self.full and compare_on_arc(self.centre, self.sense, self.start, middle, self.end) >= 0:
            return None
        return ArcSweep(self.centre, self.radius, self.sense, self.start, middle, self.kind), ArcSweep(self.centre, self.radius, self.sense, middle, self.end, self.kind)
