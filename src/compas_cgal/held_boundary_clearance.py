"""Exact distance-admissible nearest-boundary sites on a closed polygonal ring."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

from compas_cgal import _stock_2
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import SquaredMillimetre
from compas_cgal.units import WorldXY


@dataclass(frozen=True)
class BoundaryClearanceSite:
    """One distance-admissible boundary site, every value exact.

    Attributes:
        segment: Canonical ring segment index; parameter one is reported as
            parameter zero on the following segment.
        parameter: Exact clamped projection parameter in ``[0, 1)``.
        footpoint: Exact nearest point on the segment, world XY.
        squared_distance: Exact squared station-to-segment distance.
    """

    segment: int
    parameter: Fraction
    footpoint: tuple[ExactMillimetre, ExactMillimetre]
    squared_distance: SquaredMillimetre


def _scaled(numerator: int, denominator: int, exponent: int) -> Fraction:
    """Exact ``numerator / denominator * 2**exponent``; Fraction reduces it."""
    if exponent >= 0:
        return Fraction(numerator << exponent, denominator)
    return Fraction(numerator, denominator << -exponent)


def _site(record: tuple[int, str, str, str, str, str, int]) -> BoundaryClearanceSite:
    segment, along, length_squared, foot_x, foot_y, distance, exponent = record
    denominator = int(length_squared, 16)
    return BoundaryClearanceSite(
        segment,
        Fraction(int(along, 16), denominator),
        (
            ExactMillimetre(_scaled(int(foot_x, 16), denominator, exponent)),
            ExactMillimetre(_scaled(int(foot_y, 16), denominator, exponent)),
        ),
        SquaredMillimetre(_scaled(int(distance, 16), denominator, 2 * exponent)),
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
