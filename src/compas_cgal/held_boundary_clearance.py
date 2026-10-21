"""Exact distance-admissible nearest-boundary sites on a closed polygonal ring."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property

from compas_cgal import _stock_2
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import SquaredMillimetre
from compas_cgal.units import WorldXY


def _scaled(numerator: int, denominator: int, exponent: int) -> Fraction:
    """Exact ``numerator / denominator * 2**exponent``; Fraction reduces it."""
    if exponent >= 0:
        return Fraction(numerator << exponent, denominator)
    return Fraction(numerator, denominator << -exponent)


def _scaled_double(numerator: int, denominator: int, exponent: int) -> float:
    """Correctly rounded ``numerator / denominator * 2**exponent``.

    Python's ``int / int`` rounds the exact quotient once, as
    ``float(Fraction)`` does (it is the same division on the reduced pair), so
    both give the same double for the same rational.

    Raises:
        OverflowError: the quotient exceeds the binary64 range.
    """
    if exponent >= 0:
        return (numerator << exponent) / denominator
    return numerator / (denominator << -exponent)


@dataclass(frozen=True, eq=False)
class BoundaryClearanceSite:
    """One distance-admissible boundary site, every value exact.

    The native batch reports the values on one integer grid, unreduced: with
    ``L = length_squared > 0`` and ``e = exponent``, the footpoint is
    ``(foot_x, foot_y) * 2**e / L`` and the squared distance is
    ``distance * 2**(2e) / L``. Exact ``Fraction`` values are built only on
    request; the doubles are correctly rounded integer quotients, and
    ``squared_distance_exceeds`` decides in integers. Equality is exact value
    equality.

    Attributes:
        segment: Canonical ring segment index; parameter one is reported as
            parameter zero on the following segment.
        parameter: Exact clamped projection parameter in ``[0, 1)``.
        foot_x: Footpoint x numerator on the grid ``2**e / L``, world XY mm.
        foot_y: Footpoint y numerator on the grid ``2**e / L``, world XY mm.
        distance: Squared-distance numerator on the grid ``2**(2e) / L``, mm².
        length_squared: The positive common denominator ``L``.
        exponent: The binary grid exponent ``e``.
    """

    segment: int
    parameter: Fraction
    foot_x: int
    foot_y: int
    distance: int
    length_squared: int
    exponent: int

    @cached_property
    def footpoint(self) -> tuple[ExactMillimetre, ExactMillimetre]:
        """Exact nearest point on the segment, world XY."""
        return (
            ExactMillimetre(_scaled(self.foot_x, self.length_squared, self.exponent)),
            ExactMillimetre(_scaled(self.foot_y, self.length_squared, self.exponent)),
        )

    @cached_property
    def squared_distance(self) -> SquaredMillimetre:
        """Exact squared station-to-segment distance."""
        return SquaredMillimetre(_scaled(self.distance, self.length_squared, 2 * self.exponent))

    @cached_property
    def footpoint_doubles(self) -> tuple[Millimetre, Millimetre]:
        """The correctly rounded footpoint, world XY; reporting only.

        Raises:
            OverflowError: a coordinate exceeds the binary64 range.
        """
        return (
            Millimetre(_scaled_double(self.foot_x, self.length_squared, self.exponent)),
            Millimetre(_scaled_double(self.foot_y, self.length_squared, self.exponent)),
        )

    @cached_property
    def squared_distance_double(self) -> float:
        """The correctly rounded squared distance in mm²; reporting only.

        Raises:
            OverflowError: the value exceeds the binary64 range.
        """
        return _scaled_double(self.distance, self.length_squared, 2 * self.exponent)

    def squared_distance_exceeds(self, bound: SquaredMillimetre) -> bool:
        """Decide ``squared_distance > bound`` exactly, in integers.

        ``distance * 2**(2e) / L > p / q`` with ``L, q > 0`` is
        ``distance * q * 2**(2e) > p * L``.
        """
        left = self.distance * bound.denominator
        right = bound.numerator * self.length_squared
        shift = 2 * self.exponent
        if shift >= 0:
            return left << shift > right
        return left > right << -shift

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BoundaryClearanceSite):
            return NotImplemented
        return self.segment == other.segment and self.parameter == other.parameter and self.footpoint == other.footpoint and self.squared_distance == other.squared_distance

    def __hash__(self) -> int:
        return hash((self.segment, self.parameter, self.footpoint, self.squared_distance))


def _site(record: tuple[int, str, str, str, str, str, int]) -> BoundaryClearanceSite:
    """Decode one native record: parse the hex integers, reduce only the parameter."""
    segment, along, length_squared, foot_x, foot_y, distance, exponent = record
    denominator = int(length_squared, 16)
    return BoundaryClearanceSite(
        segment,
        Fraction(int(along, 16), denominator),
        int(foot_x, 16),
        int(foot_y, 16),
        int(distance, 16),
        denominator,
        exponent,
    )


class HeldBoundaryClearanceIndex:
    """Spatially indexed exact clearance query over one closed ring in world XY mm.

    Decides, for each station ``p``, the set of segments ``i`` with
    ``sqrt(d_i) - sqrt(min_j d_j) <= gap`` exactly, where ``d_i`` is the exact
    squared distance from ``p`` to segment ``i``. Coordinates and the gap
    enter the kernel once by exact double injection.
    """

    def __init__(self, native: _stock_2.HeldBoundaryClearance2) -> None:
        self._native = native

    @classmethod
    def build(cls, ring: Sequence[Point2[WorldXY]]) -> HeldBoundaryClearanceIndex:
        """Index a closed ring; the last vertex connects back to the first.

        Raises:
            InvalidHeldBoundaryClearanceInputError: fewer than two vertices.
            DegenerateHeldBoundarySegmentError: two consecutive vertices coincide.
        """
        return cls(_stock_2.HeldBoundaryClearance2([(float(point.x), float(point.y)) for point in ring]))

    @property
    def segment_count(self) -> int:
        return int(self._native.segment_count)

    def distance_admissible_sites(
        self,
        stations: Sequence[Point2[WorldXY]],
        gap: Millimetre,
    ) -> tuple[tuple[BoundaryClearanceSite, ...], ...]:
        """Return every station's admissible sites, sorted by (segment, parameter).

        Raises:
            InvalidHeldBoundaryClearanceInputError: negative or non-finite gap.
        """
        rows = self._native.distance_admissible_sites([(float(point.x), float(point.y)) for point in stations], float(gap))
        return tuple(tuple(_site(record) for record in row) for row in rows)
