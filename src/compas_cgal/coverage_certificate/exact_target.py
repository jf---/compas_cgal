"""The exact machinable target T = Omega* (+) D_r, as the certificate takes it.

Omega* is the exact centre domain of a declared (Pythagorean) boundary
(docs/held_exact_centre_domain.md). The certificate needs three things from it,
all exact:

- the declared boundary itself (rational segments with unit rational directions,
  rational arcs): its supports are target cut supports, and P-membership plus
  dist(x, dP) >= r decide membership in Omega*;
- every dOmega* piece (offset line, offset circle, junction r-circle) with its
  one-root end vertices: its band joins T;
- the vertex disks D_r(w) at the convex vertices w (one-root centres).

This module validates the types and builds the wire rows; the geometry lives in
``coverage_certificate_2/exact_target_2.h``. Building an ``ExactCentreDomainTarget``
from the benchmarks' Omega* construction is the adapter
``benchmarks.held_coverage_target.exact_coverage_target``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TypeAlias

from compas_cgal.coverage_certificate.target import OneRootMillimetre
from compas_cgal.coverage_certificate.target import VertexDisk
from compas_cgal.units import ExactMillimetre

ExactXY = tuple[ExactMillimetre, ExactMillimetre]
RationalText = tuple[str, str]
PrimitiveRow = tuple[bool, int, list[RationalText]]
PieceRow = tuple[bool, int, list[RationalText], list[RationalText], list[RationalText]]


class InvalidExactTargetError(ValueError):
    """An exact-target part is not exact (a non-Fraction), not unit where a unit is required, or empty."""


def _text(value: Fraction) -> RationalText:
    return str(value.numerator), str(value.denominator)


def _require_fractions(*values: object) -> None:
    for value in values:
        if not isinstance(value, Fraction):
            raise InvalidExactTargetError(f"exact-target parts must be Fractions, got {type(value).__name__}")


def _require_unit(vector: tuple[Fraction, Fraction], name: str) -> None:
    _require_fractions(*vector)
    if vector[0] * vector[0] + vector[1] * vector[1] != 1:
        raise InvalidExactTargetError(f"{name} {vector} is not an exactly unit vector")


def _require_sense(sense: int) -> None:
    if sense not in (1, -1) or isinstance(sense, bool):
        raise InvalidExactTargetError(f"sense must be +1 or -1, got {sense!r}")


@dataclass(frozen=True)
class OneRootPoint:
    """A point whose coordinates are rational + coefficient * sqrt(radicand), one radicand for both."""

    x: OneRootMillimetre
    y: OneRootMillimetre
    radicand: Fraction

    @classmethod
    def build(cls, x: OneRootMillimetre, y: OneRootMillimetre, radicand: Fraction) -> OneRootPoint:
        """Raises:
        InvalidExactTargetError: a part is not a Fraction or the radicand is negative.
        """
        _require_fractions(x.rational, x.coefficient, y.rational, y.coefficient, radicand)
        if radicand < 0:
            raise InvalidExactTargetError(f"one-root radicand must be non-negative, got {radicand}")
        return cls(x, y, radicand)

    def wire(self) -> list[RationalText]:
        """(x0, x1, y0, y1, radicand) as decimal text."""
        return [_text(v) for v in (self.x.rational, self.x.coefficient, self.y.rational, self.y.coefficient, self.radicand)]


@dataclass(frozen=True)
class DeclaredSegment:
    """Declared boundary line ``start + lambda direction``, lambda in [0, length]; the pocket is on its left."""

    start: ExactXY
    end: ExactXY
    direction: tuple[Fraction, Fraction]
    length: ExactMillimetre

    @classmethod
    def build(cls, start: ExactXY, direction: tuple[Fraction, Fraction], length: ExactMillimetre) -> DeclaredSegment:
        """Raises:
        InvalidExactTargetError: non-Fraction parts, a non-unit direction or a non-positive length.
        """
        _require_fractions(*start, length)
        _require_unit(direction, "segment direction")
        if length <= 0:
            raise InvalidExactTargetError(f"segment length must be positive, got {length}")
        end = (ExactMillimetre(start[0] + length * direction[0]), ExactMillimetre(start[1] + length * direction[1]))
        return cls(start, end, direction, length)

    def wire(self) -> PrimitiveRow:
        values = (*self.start, *self.end, *self.direction, self.length)
        return True, 1, [_text(v) for v in values]


@dataclass(frozen=True)
class DeclaredCircleArc:
    """Declared boundary arc from ``start`` to ``end`` about ``centre``, radius ``radius``, sense +1 (CCW) or -1."""

    start: ExactXY
    end: ExactXY
    centre: ExactXY
    radius: ExactMillimetre
    sense: int

    @classmethod
    def build(cls, start: ExactXY, end: ExactXY, centre: ExactXY, radius: ExactMillimetre, sense: int) -> DeclaredCircleArc:
        """Raises:
        InvalidExactTargetError: non-Fraction parts, a non-positive radius, a bad sense, or ends off the circle.
        """
        _require_fractions(*start, *end, *centre, radius)
        _require_sense(sense)
        if radius <= 0:
            raise InvalidExactTargetError(f"arc radius must be positive, got {radius}")
        for point in (start, end):
            if (point[0] - centre[0]) ** 2 + (point[1] - centre[1]) ** 2 != radius * radius:
                raise InvalidExactTargetError(f"arc end {point} is not on its circle")
        return cls(start, end, centre, radius, sense)

    def wire(self) -> PrimitiveRow:
        values = (*self.start, *self.end, *self.centre, self.radius)
        return False, self.sense, [_text(v) for v in values]


@dataclass(frozen=True)
class OffsetLinePiece:
    """A dOmega* piece on the offset line n.x = level (n the unit interior normal), travelled along (n_y, -n_x)."""

    normal: tuple[Fraction, Fraction]
    level: ExactMillimetre
    source: OneRootPoint
    target: OneRootPoint

    @classmethod
    def build(cls, normal: tuple[Fraction, Fraction], level: ExactMillimetre, source: OneRootPoint, target: OneRootPoint) -> OffsetLinePiece:
        """Raises:
        InvalidExactTargetError: a non-unit normal or a non-Fraction level.
        """
        _require_unit(normal, "offset-line normal")
        _require_fractions(level)
        return cls(normal, level, source, target)

    def wire(self) -> PieceRow:
        return True, 1, [_text(v) for v in (*self.normal, self.level)], self.source.wire(), self.target.wire()


@dataclass(frozen=True)
class OffsetCirclePiece:
    """A dOmega* piece on the circle (centre, radius), offsetting a primitive of radius ``boundary_radius`` (0: a junction)."""

    centre: ExactXY
    radius: ExactMillimetre
    boundary_radius: ExactMillimetre
    sense: int
    source: OneRootPoint
    target: OneRootPoint

    @classmethod
    def build(cls, centre: ExactXY, radius: ExactMillimetre, boundary_radius: ExactMillimetre, sense: int, source: OneRootPoint, target: OneRootPoint) -> OffsetCirclePiece:
        """Raises:
        InvalidExactTargetError: non-Fraction parts, R <= 0, rho < 0, or a bad sense.
        """
        _require_fractions(*centre, radius, boundary_radius)
        _require_sense(sense)
        if radius <= 0 or boundary_radius < 0:
            raise InvalidExactTargetError(f"circle piece needs R > 0 and rho >= 0, got R = {radius}, rho = {boundary_radius}")
        return cls(centre, radius, boundary_radius, sense, source, target)

    def wire(self) -> PieceRow:
        values = (*self.centre, self.radius, self.boundary_radius)
        return False, self.sense, [_text(v) for v in values], self.source.wire(), self.target.wire()


DeclaredPrimitive: TypeAlias = DeclaredSegment | DeclaredCircleArc
CentreDomainPiece: TypeAlias = OffsetLinePiece | OffsetCirclePiece


@dataclass(frozen=True)
class ExactCentreDomainTarget:
    """T = Omega* (+) D_r: the declared boundary, every dOmega* piece and the convex-vertex disks.

    Omega* must be connected: the certificate's interior-anchor rule needs int T
    connected (``exact_certificate_2.h``), which holds for one component of
    Omega*. The benchmarks' ``CentreDomain.build`` guarantees it (a neck that
    splits Omega* is a hard error there), and ``build`` checks the count again.
    """

    boundary: tuple[DeclaredPrimitive, ...]
    pieces: tuple[CentreDomainPiece, ...]
    vertex_disks: tuple[VertexDisk, ...]

    @classmethod
    def build(
        cls, boundary: Sequence[DeclaredPrimitive], pieces: Sequence[CentreDomainPiece], vertex_disks: Sequence[VertexDisk], *, components: int
    ) -> ExactCentreDomainTarget:
        """Raises:
        InvalidExactTargetError: an empty boundary or piece list, or Omega* not one component.
        """
        if not boundary or not pieces:
            raise InvalidExactTargetError("the exact target needs a declared boundary and at least one dOmega* piece")
        if components != 1:
            raise InvalidExactTargetError(f"Omega* must be one component for the interior-anchor rule, got {components}")
        return cls(tuple(boundary), tuple(pieces), tuple(vertex_disks))

    def wire(self) -> tuple[list[PrimitiveRow], list[PieceRow], list[list[RationalText]]]:
        """The native rows: boundary, pieces, vertex disks."""
        return [p.wire() for p in self.boundary], [p.wire() for p in self.pieces], [d.wire() for d in self.vertex_disks]
