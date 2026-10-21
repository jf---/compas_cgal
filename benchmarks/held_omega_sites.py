"""Exact sites on ∂Ω* for Held path generation: order keys, points and inward normals.

A site is a decided piece of ∂Ω* plus a rational probe whose orthogonal
projection onto the piece's support is the site (``None`` probe: the piece's
source vertex). Projection preserves the order data exactly without a root:

- on an offset line the position is the rational parameter ``u·(probe - anchor)``;
- on a circle the position is the direction of the rational vector
  ``probe - centre`` (the site lies on that ray), ordered by its half-turn from
  the piece's source and then by orientation, exactly as
  :func:`benchmarks.held_one_root.compare_on_arc`.

So sites sort by a key of rationals and exact comparisons; the one-root point
(``c + ρ'(probe - c)/√D`` on a circle) is built only when a consumer needs it.
Every sign is decided by a dyadic interval filter first
(:mod:`benchmarks.held_dyadic_interval`) and exactly when the filter declines.
The order is :meth:`benchmarks.held_exact_centre_domain.CentreDomainSite.compare`
restricted to one CCB (contract-tested).
"""

from __future__ import annotations

import enum
import functools
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Self
from typing import TypeAlias

from benchmarks.errors import BenchmarkError
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import CentreDomainPiece
from benchmarks.held_exact_centre_domain import CentreDomainSite
from benchmarks.held_exact_centre_domain import DecidedContact
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import orientation
from compas_cgal.units import WorldXY

_Vector: TypeAlias = tuple[Fraction, Fraction]
_RootVector: TypeAlias = tuple[OneRootNumber, OneRootNumber]
_IntervalVector: TypeAlias = tuple[Interval, Interval]


class OmegaSiteError(BenchmarkError):
    """Base error of exact ∂Ω* sites."""


class MultiLoopCentreDomainError(OmegaSiteError):
    """Ω* has holes: one cyclic ∂Ω* order is not defined (the Held workloads have none)."""


class ProbeAtCircleCentreError(OmegaSiteError):
    """A circle site's probe is the circle's centre, so it names no direction."""


class SiteOutsidePieceError(OmegaSiteError):
    """A materialised site is not before its piece's target (a decision-contract breach)."""


def _enclose(vector: _RootVector) -> _IntervalVector:
    return Interval.of_one_root(vector[0]), Interval.of_one_root(vector[1])


def _cross_sign(left: _RootVector, left_box: _IntervalVector, right: _Vector) -> int:
    """Filtered exact ``sign(left × right)`` for a one-root ``left`` and a rational ``right``."""
    value = left_box[0] * Interval.of(right[1]) - left_box[1] * Interval.of(right[0])
    decided = value.sign()
    return decided if decided is not None else (left[0] * right[1] - left[1] * right[0]).sign()


def _dot_sign(left: _RootVector, left_box: _IntervalVector, right: _Vector) -> int:
    """Filtered exact ``sign(left · right)``."""
    value = left_box[0] * Interval.of(right[0]) + left_box[1] * Interval.of(right[1])
    decided = value.sign()
    return decided if decided is not None else (left[0] * right[0] + left[1] * right[1]).sign()


@dataclass(frozen=True)
class PieceFrame:
    """Per-piece data the site predicates reuse, exact and enclosed.

    On a line: the source and target parameters along ``u`` (first entries). On
    a circle: the source and target vectors from the centre, and the target's
    half-turn from the source (``1`` for a closed circle).
    """

    index: int
    piece: CentreDomainPiece
    closed: bool
    source: _RootVector
    source_box: _IntervalVector
    target: _RootVector
    target_box: _IntervalVector
    target_half: int
    source_point_box: _IntervalVector
    origin_box: _IntervalVector
    axis_box: _IntervalVector
    radius_box: Interval

    @classmethod
    def build(cls, index: int, piece: CentreDomainPiece) -> Self:
        support = piece.support
        closed = piece.source == piece.target
        zero = OneRootNumber.of(0)
        point_box = _enclose((piece.source.x, piece.source.y))
        if isinstance(support, OffsetLineSupport):
            source, target = (support.parameter(piece.source), zero), (support.parameter(piece.target), zero)
            anchor = (Interval.of(support.anchor.x), Interval.of(support.anchor.y))
            direction = (Interval.of(support.direction[0]), Interval.of(support.direction[1]))
            return cls(index, piece, closed, source, _enclose(source), target, _enclose(target), 0, point_box, anchor, direction, Interval.of(0))
        source, target = piece.source.minus(support.centre), piece.target.minus(support.centre)
        # Distinct points of one circle are never in one direction from its
        # centre, so a collinear target is the antipode (half-turn 1).
        side = support.sense * orientation(support.centre, piece.source, piece.target)
        half = 1 if closed or side <= 0 else 0
        centre = (Interval.of(support.centre.x), Interval.of(support.centre.y))
        return cls(
            index, piece, closed, source, _enclose(source), target, _enclose(target), half, point_box, centre, (Interval.of(0), Interval.of(0)), Interval.of(support.radius)
        )

    def half_turn(self, vector: _Vector) -> int:
        """``0`` for a swept angle from the source in ``[0, π)``, ``1`` for ``[π, 2π)`` (filtered exact)."""
        support = self.piece.support
        assert isinstance(support, OffsetCircleSupport)
        side = support.sense * _cross_sign(self.source, self.source_box, vector)
        if side != 0:
            return 0 if side > 0 else 1
        return 0 if _dot_sign(self.source, self.source_box, vector) > 0 else 1

    def before_target(self, vector: _Vector, half: int) -> bool:
        """Whether the ray along ``vector`` (half-turn ``half``) meets the arc strictly before its target."""
        support = self.piece.support
        assert isinstance(support, OffsetCircleSupport)
        if self.closed or half != self.target_half:
            return self.closed or half < self.target_half
        # Same half-turn: vector < target iff σ·cross(vector, target) > 0.
        return support.sense * -_cross_sign(self.target, self.target_box, vector) > 0

    def compare_along(self, along: Fraction, end: _RootVector, box: _IntervalVector) -> int:
        """Filtered exact ``sign(along - end parameter)`` on a line piece."""
        decided = (Interval.of(along) - box[0]).sign()
        return decided if decided is not None else -end[0].compare(along)


@dataclass(frozen=True)
class OmegaLoop:
    """Ω*'s one CCB with per-piece frames."""

    domain: CentreDomain
    pieces: tuple[CentreDomainPiece, ...]
    frames: tuple[PieceFrame, ...]
    convex_sources: tuple[bool, ...]

    @classmethod
    def build(cls, domain: CentreDomain) -> Self:
        """Frames of every piece of the outer CCB, and which piece sources are convex vertices.

        Raises:
            MultiLoopCentreDomainError: If Ω* has a hole.
        """
        (component,) = domain.components
        if component.holes:
            raise MultiLoopCentreDomainError(f"Ω* has {len(component.holes)} hole(s); the Held path follows one loop.")
        pieces = component.outer
        convex = []
        for index, piece in enumerate(pieces):
            incoming = pieces[index - 1]
            first, second = incoming.support.travel_tangent(piece.source), piece.support.travel_tangent(piece.source)
            # The loop turns left (convex, the face on the left) at the source.
            convex.append(len(pieces) > 1 and (first[0] * second[1] - first[1] * second[0]).sign() > 0)
        return cls(domain, pieces, tuple(PieceFrame.build(index, piece) for index, piece in enumerate(pieces)), tuple(convex))


class _Placement(enum.IntEnum):
    """Where on its piece a site lies: at the source vertex, or interior."""

    SOURCE = 0
    INTERIOR = 1


@dataclass(frozen=True)
class _Direction:
    """A rational direction about a circle centre, ordered within one half-turn in the travel sense."""

    vector: _Vector
    sense: int
    box: _IntervalVector

    def compare(self, other: _Direction) -> int:
        """Exact three-way order within one half-turn: first < second iff σ·cross(first, second) > 0."""
        decided = (self.box[0] * other.box[1] - self.box[1] * other.box[0]).sign()
        if decided is None:
            value = self.vector[0] * other.vector[1] - self.vector[1] * other.vector[0]
            decided = (value > 0) - (value < 0)
        return -self.sense * decided


def _three_way(first: Fraction | int, second: Fraction | int) -> int:
    return (first > second) - (first < second)


@functools.total_ordering
@dataclass(frozen=True, eq=False)
class SiteKey:
    """Exact sort key of a site on Ω*'s one loop: piece, source-or-interior, then position on the piece.

    ``along`` is the rational line parameter; ``half`` and ``direction`` order a
    circle site (half-turn from the piece's source, then orientation).
    """

    piece: int
    placement: _Placement
    along: Fraction | None = None
    half: int | None = None
    direction: _Direction | None = None

    def compare(self, other: SiteKey) -> int:
        """Exact three-way order along the loop (one exact decision per call)."""
        if self.piece != other.piece:
            return _three_way(self.piece, other.piece)
        if self.placement is not other.placement:
            return _three_way(self.placement, other.placement)
        if self.placement is _Placement.SOURCE:
            return 0
        if self.along is not None and other.along is not None:
            return _three_way(self.along, other.along)
        assert self.half is not None and other.half is not None and self.direction is not None and other.direction is not None
        if self.half != other.half:
            return _three_way(self.half, other.half)
        return self.direction.compare(other.direction)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SiteKey):
            return NotImplemented
        return self.compare(other) == 0

    def __lt__(self, other: SiteKey) -> bool:
        return self.compare(other) < 0

    __hash__ = None  # type: ignore[assignment]


@dataclass(frozen=True, eq=False)
class OmegaSite:
    """An exact point of ∂Ω*: piece ``frame.index``, at its source or projected from ``probe``."""

    frame: PieceFrame
    probe: ExactPoint2[WorldXY] | None
    key: SiteKey
    _point: list[OneRootPoint2] = field(default_factory=list, repr=False)

    @property
    def index(self) -> int:
        return self.frame.index

    @property
    def piece(self) -> CentreDomainPiece:
        return self.frame.piece

    @classmethod
    def at_source(cls, frame: PieceFrame) -> Self:
        """The source vertex of the frame's piece."""
        return cls(frame, None, SiteKey(frame.index, _Placement.SOURCE))

    @classmethod
    def project(cls, frame: PieceFrame, probe: ExactPoint2[WorldXY]) -> Self:
        """The orthogonal projection of ``probe`` onto the piece's support (normalised to the source when equal).

        Raises:
            ProbeAtCircleCentreError: If the probe is a circle piece's centre.
            SiteOutsidePieceError: If the projection is not strictly before the piece's target.
        """
        support = frame.piece.support
        if isinstance(support, OffsetLineSupport):
            along = support.direction[0] * (probe.x - support.anchor.x) + support.direction[1] * (probe.y - support.anchor.y)
            versus_source = frame.compare_along(along, frame.source, frame.source_box)
            if versus_source == 0:
                return cls.at_source(frame)
            if versus_source < 0 or frame.compare_along(along, frame.target, frame.target_box) >= 0:
                raise SiteOutsidePieceError(f"Parameter {float(along)} is outside line piece {frame.index}.")
            return cls(frame, probe, SiteKey(frame.index, _Placement.INTERIOR, along=along))
        vector = probe.minus(support.centre)
        if vector[0] == 0 and vector[1] == 0:
            raise ProbeAtCircleCentreError(f"Probe {probe} is the centre of piece {frame.index}.")
        half = frame.half_turn(vector)
        if half == 0 and _cross_sign(frame.source, frame.source_box, vector) == 0:
            return cls.at_source(frame)
        if not frame.before_target(vector, half):
            raise SiteOutsidePieceError(f"Probe {probe} projects at or beyond the target of arc piece {frame.index}.")
        return cls(
            frame, probe, SiteKey(frame.index, _Placement.INTERIOR, half=half, direction=_Direction(vector, support.sense, (Interval.of(vector[0]), Interval.of(vector[1]))))
        )

    @classmethod
    def decided(cls, loop: OmegaLoop, contact: DecidedContact) -> Self:
        """The site of a native contact decision."""
        frame = loop.frames[contact.piece.index]
        if contact.is_vertex:
            return cls.at_source(frame)
        return cls.project(frame, ExactPoint2[WorldXY].from_binary64(*contact.query))

    @property
    def at_vertex(self) -> bool:
        """Whether the site is its piece's source vertex."""
        return self.probe is None

    def same_position(self, other: OmegaSite) -> bool:
        """Exact equality of the two sites' positions on ∂Ω*."""
        return self.key == other.key

    def point(self) -> OneRootPoint2:
        """The exact one-root point, materialised once (checked on the support)."""
        if not self._point:
            self._point.append(self._materialise())
        return self._point[0]

    def _materialise(self) -> OneRootPoint2:
        if self.probe is None:
            return self.piece.source
        support = self.piece.support
        if isinstance(support, OffsetLineSupport):
            along = self.key.along
            assert along is not None
            foot = support.anchor.offset(support.direction, along)
            point = OneRootPoint2.rational(foot.x, foot.y)
        else:
            dx, dy = self.probe.minus(support.centre)
            radial_squared = dx * dx + dy * dy
            scale = support.radius / radial_squared
            point = OneRootPoint2.build(OneRootNumber.build(support.centre.x, scale * dx, radial_squared), OneRootNumber.build(support.centre.y, scale * dy, radial_squared))
        return CentreDomainSite.build(self.piece, point).point

    def inward_normal(self) -> _RootVector:
        """Unit normal of ∂Ω* into Ω* at the site (the outgoing piece's at a vertex).

        On a line it is the rational ``n̂``; on a circle ``-σ (q - c) / ρ'``,
        in the site's extension.
        """
        support = self.piece.support
        if isinstance(support, OffsetLineSupport):
            return OneRootNumber.of(support.normal[0]), OneRootNumber.of(support.normal[1])
        dx, dy = self.point().minus(support.centre)
        scale = Fraction(-support.sense) / support.radius
        return scale * dx, scale * dy

    def point_box(self) -> _IntervalVector:
        """Enclosure of the site's point, computed from the probe without roots."""
        frame = self.frame
        if self.probe is None:
            return frame.source_point_box
        probe = (Interval.of(self.probe.x), Interval.of(self.probe.y))
        if isinstance(self.piece.support, OffsetLineSupport):
            along = self.key.along
            assert along is not None
            step = Interval.of(along)
            return frame.origin_box[0] + step * frame.axis_box[0], frame.origin_box[1] + step * frame.axis_box[1]
        wx, wy = probe[0] - frame.origin_box[0], probe[1] - frame.origin_box[1]
        length = (wx.square() + wy.square()).sqrt()
        if length.lo <= 0:
            # A probe within the enclosure resolution of the centre: enclose the exact point instead.
            point = self.point()
            return Interval.of_one_root(point.x), Interval.of_one_root(point.y)
        return frame.origin_box[0] + frame.radius_box * wx / length, frame.origin_box[1] + frame.radius_box * wy / length

    def squared_distance_to(self, origin: ExactPoint2[WorldXY]) -> OneRootNumber:
        """Exact ``|q - origin|²`` (one-root on a circle: ``|c - o|² + ρ'² + 2ρ'(c - o)·w/√D``)."""
        support = self.piece.support
        if self.probe is None or isinstance(support, OffsetLineSupport):
            dx, dy = self.point().minus(origin)
            return dx * dx + dy * dy
        wx, wy = self.probe.minus(support.centre)
        radial = wx * wx + wy * wy
        cx, cy = support.centre.minus(origin)
        rho = support.radius
        return OneRootNumber.build(cx * cx + cy * cy + rho * rho, 2 * rho * (cx * wx + cy * wy) / radial, radial)


def single_loop_pieces(domain: CentreDomain) -> tuple[CentreDomainPiece, ...]:
    """The pieces of Ω*'s one CCB, in CCW order.

    Raises:
        MultiLoopCentreDomainError: If Ω* has a hole.
    """
    return OmegaLoop.build(domain).pieces


def sort_sites(sites: list[OmegaSite]) -> list[int]:
    """Indices of ``sites`` in exact ∂Ω* order (stable for equal positions)."""
    return sorted(range(len(sites)), key=lambda index: sites[index].key)
