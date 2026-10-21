"""Exact Held centre domain Ω* = P ⊖ D_r on the declared (Pythagorean) boundary.

Design: docs/held_exact_centre_domain.md. User rulings (2026-10-17):

- the declared boundary (:mod:`benchmarks.held_declared_boundary`) IS the
  boundary; "exact" means exact with respect to it;
- the certified object is the exact path (exact ∂Ω* connectors, binary64
  circles);
- a neck narrower than or equal to 2r, which splits Ω* into more than one
  component, is a hard named error (:class:`NeckSplitsCentreDomainError`);
- the coverage certificate gains the one-root-centre ``VERTEX_DISK`` member
  (:class:`VertexDisk`).

Pure rational parts live beside this module: one-root numbers and points in
:mod:`benchmarks.held_one_root`, every support of ∂Ω* and the local junction
vertices in :mod:`benchmarks.held_centre_domain_supports`. This module holds the
domain itself: pieces, sites and their exact order, the identity, and the
native seam.

What is NOT here is the arrangement: the aggregated join of the members of
``∂P ⊕ D_r`` and the difference ``P ∖ U``, which trims the local pieces against
each other globally (necks, far interactions, arcs tighter than the tool). It is
native (``src/held_centre_domain_2.cpp``, the CORE-free extension
``_held_centre_domain_2``; wire layer :mod:`benchmarks.held_centre_domain_native`)
and sits behind :func:`construct_centre_domain_arrangement`. Its pieces are
matched back to the rational supports here by exact support keys.

Precondition: the declared cycle is simple and counter-clockwise (the
``ReferenceBoundary`` contract), so the interior is on the left of every
primitive. It is checked exactly in Python before the native call
(:func:`benchmarks.held_centre_domain_supports.validate_simple_counter_clockwise`),
because the native construction does not check it yet.

Unreachable attribution is by unreached boundary: a boundary point b is
reached iff its inward centre b + r·n̂ lies in Ω*, so the reached part of a
primitive is the contact-map image of the ∂Ω* pieces it owns plus the points
the vertex disks touch, and its complement is exactly the boundary that
unreachable material of P ∖ T touches (see :class:`UnreachableAttribution`).
"""

from __future__ import annotations

import enum
import functools
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Any
from typing import Self
from typing import TypeAlias

from benchmarks.held_centre_domain_errors import CentreDomainError
from benchmarks.held_centre_domain_native import AMBIGUOUS
from benchmarks.held_centre_domain_native import BEYOND_GAP
from benchmarks.held_centre_domain_native import NativeCircle
from benchmarks.held_centre_domain_native import NativeLine
from benchmarks.held_centre_domain_native import NativePiece
from benchmarks.held_centre_domain_native import build_native_centre_domain
from benchmarks.held_centre_domain_native import native_decide_binary64
from benchmarks.held_centre_domain_native import native_locate
from benchmarks.held_centre_domain_native import native_nearest
from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_centre_domain_supports import CentreDomainSupport
from benchmarks.held_centre_domain_supports import CentreDomainSupports
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_centre_domain_supports import PrimitiveBox
from benchmarks.held_centre_domain_supports import declared_primitive_box
from benchmarks.held_centre_domain_supports import on_declared_arc
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredBoundaryIdentity
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import DeclaredPrimitive
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import WorldXY

# Version of the centre-domain contract embedded in every identity hash.
CENTRE_DOMAIN_CONTRACT_VERSION = 1
# Name of the native construction component the arrangement seam calls.
NATIVE_CENTRE_DOMAIN_COMPONENT = "_held_centre_domain_2.CentreDomain2"
# Width of the x-strips (in tool radii) that bucket primitives for the
# vertex-disk touch search; only prunes, every candidate is decided exactly.
TOUCH_STRIP_TOOL_RADII = 4
# Bits of the rational bracket of a one-root vertex coordinate used to bucket
# it; only prunes (the bracket is exact: integer square roots).
VERTEX_BRACKET_BITS = 40


class PointNotOnPieceError(CentreDomainError):
    """A site position was requested for a point that is not exactly on the piece's support."""


class UnmatchedCentreDomainPieceError(CentreDomainError):
    """A native ∂Ω* piece lies on no rational support, or travels against its support (a construction contract breach)."""


class ContactBeyondGapError(CentreDomainError):
    """The nearest ∂Ω* site of a hypothesis contact is farther than the admissible gap."""


class AmbiguousContactProjectionError(CentreDomainError):
    """Two distinct ∂Ω* points are exactly nearest to a hypothesis contact."""


class EmptyCentreDomainError(CentreDomainError):
    """Ω* has no component: the tool fits nowhere in the pocket."""


class NeckSplitsCentreDomainError(CentreDomainError):
    """A neck no wider than 2r splits Ω* into more than one component (user ruling: hard error)."""

    def __init__(self, message: str, component_count: int) -> None:
        super().__init__(message)
        self.component_count = component_count


class InconsistentUnreachableAttributionError(CentreDomainError):
    """Unreached intervals of adjacent primitives disagree at their shared junction (an attribution contract breach)."""


# ---------------------------------------------------------------------------
# Pieces, sites and order
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CentreDomainPiece:
    """One curve of ∂Ω*, from ``source`` to ``target`` in CCW order along its CCB.

    ``owners`` is the full owner set of the support (coincident supports, e.g.
    collinear lines, carry every owner).
    """

    component: int
    ccb: int
    index: int
    owners: frozenset[BoundaryFeature]
    support: CentreDomainSupport
    source: OneRootPoint2
    target: OneRootPoint2

    @property
    def kind(self) -> PieceKind:
        """The piece kind, from its support."""
        return self.support.kind


@dataclass(frozen=True)
class CentreDomainSite:
    """A point on piece ``piece`` of CCB ``ccb`` of component ``component``.

    The order is lexicographic in ``(component, ccb, piece)`` then exact position
    along the piece (:func:`compare_on_piece`).
    """

    component: int
    ccb: int
    piece: CentreDomainPiece
    point: OneRootPoint2

    @classmethod
    def build(cls, piece: CentreDomainPiece, point: OneRootPoint2) -> Self:
        """Validate that ``point`` lies exactly on the piece's support.

        Raises:
            PointNotOnPieceError: If it does not.
        """
        if not piece.support.contains(point):
            raise PointNotOnPieceError(f"{point!r} is not on piece {piece.index} ({piece.kind.value}).")
        return cls(piece.component, piece.ccb, piece, point)

    def compare(self, other: CentreDomainSite) -> int:
        """Exact three-way order of two sites."""
        mine = (self.component, self.ccb, self.piece.index)
        theirs = (other.component, other.ccb, other.piece.index)
        if mine != theirs:
            return -1 if mine < theirs else 1
        return compare_on_piece(self.piece, self.point, other.point)

    def contact(self) -> OneRootPoint2:
        """The boundary contact of this site (a rational affine image, same extension)."""
        return self.piece.support.contact(self.point)


@dataclass(frozen=True, eq=False)
class DecidedContact:
    """A binary64 contact whose nearest ∂Ω* piece is decided exactly; its one-root site is built on demand.

    ``is_vertex`` means the nearest point is the source of ``piece``; otherwise it
    is interior to it, i.e. the orthogonal projection of the contact onto the
    piece's support (rational on a line, ``c + R(q - c)/√D`` on a circle).
    """

    query: tuple[float, float]
    piece: CentreDomainPiece
    is_vertex: bool

    def site(self) -> CentreDomainSite:
        """Materialise the exact site (pure ``Fraction`` work, no native call).

        Raises:
            PointNotOnPieceError: If the projection is not on the piece's support (a contract breach).
        """
        if self.is_vertex:
            return CentreDomainSite.build(self.piece, self.piece.source)
        return CentreDomainSite.build(self.piece, _projection(self.piece.support, ExactPoint2[WorldXY].from_binary64(*self.query)))


@dataclass(frozen=True, eq=False)
class ContactDecisions:
    """One batch projection: per contact its decision (``None`` beyond the gap), and the exact-stage count."""

    contacts: tuple[DecidedContact | None, ...]
    exact_fallbacks: int


def _projection(support: CentreDomainSupport, query: ExactPoint2[WorldXY]) -> OneRootPoint2:
    """Orthogonal projection of ``query`` onto the support (the query is not a circle's centre)."""
    if isinstance(support, OffsetLineSupport):
        along = support.direction[0] * (query.x - support.anchor.x) + support.direction[1] * (query.y - support.anchor.y)
        foot = support.anchor.offset(support.direction, along)
        return OneRootPoint2.rational(foot.x, foot.y)
    dx, dy = query.minus(support.centre)
    radial_squared = dx * dx + dy * dy
    scale = support.radius / radial_squared
    return OneRootPoint2.build(OneRootNumber.build(support.centre.x, scale * dx, radial_squared), OneRootNumber.build(support.centre.y, scale * dy, radial_squared))


def compare_on_piece(piece: CentreDomainPiece, first: OneRootPoint2, second: OneRootPoint2) -> int:
    """Exact order of two points of ``piece`` along its travel direction.

    A line orders by the exact parameter ``u·(q - anchor)``; an arc orders by
    orientation about its centre from ``piece.source`` in the travel sense (the
    arc-prefix idiom: half-turn classification, then orientation, no angles).
    """
    support = piece.support
    if isinstance(support, OffsetLineSupport):
        return support.parameter(first).compare(support.parameter(second))
    return compare_on_arc(support.centre, support.sense, piece.source, first, second)


# ---------------------------------------------------------------------------
# The domain
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CentreDomainComponent:
    """One face of Ω*: an outer CCB (CCW) and its holes (CW), each a closed piece chain."""

    outer: tuple[CentreDomainPiece, ...]
    holes: tuple[tuple[CentreDomainPiece, ...], ...]


@dataclass(frozen=True)
class VertexDisk:
    """Coverage-target member D_r(w) at a convex non-G1 vertex w of Ω* (one-root centre).

    ``owners`` are the owners of the two pieces meeting at w: a convex boundary
    kink, or a crossing of non-adjacent offsets (an arc tighter than the tool
    in between, or a neck), so no single junction names every such vertex.
    """

    centre: OneRootPoint2
    radius: ExactMillimetre
    owners: frozenset[BoundaryFeature]


class UnreachableCause(enum.Enum):
    """Why a stretch of ∂P is unreached by T = Ω* ⊕ D_r."""

    CONVEX_KINK_CUSP = "convex_kink_cusp"
    TIGHT_CONVEX_ARC = "tight_convex_arc"
    NECK = "neck"


@dataclass(frozen=True)
class UnreachableInterval:
    """A maximal unreached sub-interval of declared primitive ``primitive``, ``start`` → ``end`` in travel order.

    The ends are exact points of the primitive; ``start_closed`` / ``end_closed``
    say whether the end point itself is unreached. A junction point that is
    unreached closes the interval on both primitives that meet there.
    """

    primitive: BoundaryFeature
    start: OneRootPoint2
    end: OneRootPoint2
    start_closed: bool
    end_closed: bool


@dataclass(frozen=True)
class UnreachableAttribution:
    """One connected stretch of unreached ∂P (its intervals in boundary order) and its cause.

    Completeness (lemma): if x ∈ P ∖ T and b is its nearest point of ∂P, the
    inward centre b + r·n̂ is not in Ω* (else x would lie in its disk), so b is
    unreached. Every part of P ∖ T therefore touches one of these stretches,
    and the stretches are exact: a boundary point is reached iff it is the
    contact of a point of Ω*, i.e. lies in the closed contact-map image of an
    owned ∂Ω* piece, is touched by a vertex disk, or is an arc of ρ = r whose
    centre is in Ω* (decided by the native ``locate``).

    Cause: a stretch is *local* when the Ω* points that reach its two ends
    share one vertex w (the stretch lies under the single vertex disk D_r(w)).
    A local stretch containing a convex arc with ρ < r is
    :attr:`UnreachableCause.TIGHT_CONVEX_ARC`; else one containing a convex
    kink is :attr:`UnreachableCause.CONVEX_KINK_CUSP`; ``owners`` are those
    causal features. Every other stretch, e.g. a dead-end slot no wider than
    2r, is :attr:`UnreachableCause.NECK`, owned by its primitives and causal
    features.
    """

    cause: UnreachableCause
    owners: frozenset[BoundaryFeature]
    intervals: tuple[UnreachableInterval, ...]


class PointLocation(enum.Enum):
    """Exact location of a point relative to Ω*."""

    INTERIOR = "interior"
    BOUNDARY = "boundary"
    EXTERIOR = "exterior"


@dataclass(frozen=True)
class CentreDomainArrangement:
    """What the native construction returns, matched to supports: faces, attributed unreachable parts, vertex disks."""

    components: tuple[CentreDomainComponent, ...]
    unreachable: tuple[UnreachableAttribution, ...]
    vertex_disks: tuple[VertexDisk, ...]
    native: Any = field(compare=False, repr=False)


@dataclass(frozen=True)
class CentreDomainIdentity:
    """Content address of Ω*: contract version, native component, boundary identity, tool radius."""

    contract_version: int
    native_component: str
    boundary: DeclaredBoundaryIdentity
    tool_radius: ExactMillimetre

    @property
    def sha256(self) -> str:
        """Hash of every identity field."""
        document = {
            "contract_version": self.contract_version,
            "native_component": self.native_component,
            "declared_sha256": self.boundary.declared_sha256,
            "source_sha256": self.boundary.source_sha256,
            "tool_radius": f"{self.tool_radius.numerator}/{self.tool_radius.denominator}",
        }
        return hashlib.sha256(json.dumps(document, sort_keys=True).encode()).hexdigest()


# ---------------------------------------------------------------------------
# Matching native pieces to supports (exact support keys)
# ---------------------------------------------------------------------------

_LineKey = tuple[Fraction, Fraction, Fraction]
_CircleKey = tuple[Fraction, Fraction, Fraction]


def _line_key(a: Fraction, b: Fraction, c: Fraction) -> _LineKey:
    """``a·x + b·y + c = 0`` scaled so its first non-zero coefficient is 1 (exact, sign-free)."""
    leading = a if a != 0 else b
    return a / leading, b / leading, c / leading


def _support_index(supports: CentreDomainSupports) -> tuple[dict[_LineKey, list[OffsetLineSupport]], dict[_CircleKey, list[OffsetCircleSupport]]]:
    lines: dict[_LineKey, list[OffsetLineSupport]] = {}
    for line in supports.offset_lines:
        lines.setdefault(_line_key(line.normal[0], line.normal[1], -line.level), []).append(line)
    circles: dict[_CircleKey, list[OffsetCircleSupport]] = {}
    for circle in (*supports.offset_arcs, *supports.junction_arcs):
        circles.setdefault((circle.centre.x, circle.centre.y, circle.radius_squared), []).append(circle)
    return lines, circles


def _matched_piece(
    native: NativePiece,
    position: tuple[int, int, int],
    index: tuple[dict[_LineKey, list[OffsetLineSupport]], dict[_CircleKey, list[OffsetCircleSupport]]],
) -> CentreDomainPiece:
    """Match a native piece to every support on its exact key that it travels along.

    A coincident support of the opposite sense (e.g. an offset arc of ρ = 2r
    centred on a reflex junction, whose r-circle is the same circle travelled
    clockwise) is not an owner of this piece; a piece no candidate travels with
    is a contract breach.
    """
    lines, circles = index
    support = native.support
    candidates: list[OffsetLineSupport] | list[OffsetCircleSupport]
    if isinstance(support, NativeLine):
        candidates = lines.get(_line_key(support.a, support.b, support.c), [])
    else:
        assert isinstance(support, NativeCircle)
        candidates = circles.get((support.centre_x, support.centre_y, support.radius_squared), [])
    if not candidates:
        raise UnmatchedCentreDomainPieceError(f"∂Ω* piece {position} lies on {support}, which offsets no declared feature.")
    travelling: list[CentreDomainSupport] = [candidate for candidate in candidates if _travels_along(native, candidate)]
    if not travelling:
        raise UnmatchedCentreDomainPieceError(f"∂Ω* piece {position} travels against every support on its key ({[candidate.owner for candidate in candidates]}).")
    owners = frozenset(candidate.owner for candidate in travelling)
    return CentreDomainPiece(position[0], position[1], position[2], owners, travelling[0], native.source, native.target)


def _travels_along(native: NativePiece, candidate: CentreDomainSupport) -> bool:
    if isinstance(candidate, OffsetLineSupport):
        return candidate.parameter(native.target).compare(candidate.parameter(native.source)) > 0
    return native.sense == candidate.sense


# ---------------------------------------------------------------------------
# Vertex disks
# ---------------------------------------------------------------------------


def _cross(left: tuple[OneRootNumber, OneRootNumber], right: tuple[OneRootNumber, OneRootNumber]) -> int:
    return (left[0] * right[1] - left[1] * right[0]).sign()


def _dot(left: tuple[OneRootNumber, OneRootNumber], right: tuple[OneRootNumber, OneRootNumber]) -> int:
    return (left[0] * right[0] + left[1] * right[1]).sign()


def _vertex_disks(components: tuple[CentreDomainComponent, ...], tool_radius: ExactMillimetre) -> tuple[VertexDisk, ...]:
    """D_r(w) at every convex vertex w of every cycle: a left turn, or a π-turn cusp.

    The face is on the left. A π-turn (``cross = 0``, ``dot < 0``) is always a
    convex spike of Ω*: the complement of Ω* is a union of open r-disks, so it
    satisfies an exterior r-ball condition and has no inward cusp.
    """
    disks: list[VertexDisk] = []
    for component in components:
        for cycle in (component.outer, *component.holes):
            if len(cycle) == 1:
                continue
            for incoming, outgoing in zip(cycle, (*cycle[1:], cycle[0])):
                vertex = incoming.target
                arriving, leaving = incoming.support.travel_tangent(vertex), outgoing.support.travel_tangent(vertex)
                turn = _cross(arriving, leaving)
                if turn > 0 or (turn == 0 and _dot(arriving, leaving) < 0):
                    disks.append(VertexDisk(vertex, tool_radius, incoming.owners | outgoing.owners))
    return tuple(disks)


# ---------------------------------------------------------------------------
# Unreachable attribution: reached boundary per primitive, exact
# ---------------------------------------------------------------------------

# The Ω* points that reach an end of a reached interval (empty when the end is
# a primitive end reached from a piece's interior).
_Witness: TypeAlias = tuple[OneRootPoint2, ...]


@dataclass
class _Reach:
    """A closed reached interval ``[low, high]`` of one primitive, its order keys, and the Ω* points that reach its ends."""

    low: OneRootPoint2
    high: OneRootPoint2
    low_key: OneRootNumber | None
    high_key: OneRootNumber | None
    low_witness: _Witness
    high_witness: _Witness


@dataclass(frozen=True)
class _Gap:
    """An unreached interval of primitive ``position`` and the witnesses of the reached ends that bound it."""

    position: int
    interval: UnreachableInterval
    start_witness: _Witness
    end_witness: _Witness


class _PrimitiveTrack:
    """Exact order along one declared primitive, and the reached intervals collected on it.

    Order is by an exact scalar key where one exists: ``u·(p - start)`` on a
    line, and on an arc of sweep ≤ π the chord projection ``(E - S)·(p - S)``,
    which is strictly increasing along such an arc (its derivative is
    ``ρ·|E - S|·cos(φ - θ/2) > 0`` for ``φ ∈ (0, θ)``). Keys of points of one
    extension are same-root numbers, compared across roots by squaring. A
    major arc (sweep > π) orders by :func:`compare_on_arc`.
    """

    def __init__(self, primitive: DeclaredPrimitive) -> None:
        self.primitive = primitive
        self.start = OneRootPoint2.from_exact(primitive.start)
        self.end = OneRootPoint2.from_exact(primitive.end)
        self.feature = BoundaryFeature(BoundaryFeatureKind.ARC if isinstance(primitive, DeclaredArc) else BoundaryFeatureKind.LINE, primitive.source_index)
        self.axis: tuple[Fraction, Fraction] | None = primitive.direction.vector if isinstance(primitive, DeclaredLine) else _minor_arc_chord(primitive)
        self.start_key = self.key(self.start)
        self.end_key = self.key(self.end)
        self.reached: list[_Reach] = []

    def key(self, point: OneRootPoint2) -> OneRootNumber | None:
        """The exact order key of a point of the primitive's support (``None`` on a major arc)."""
        if self.axis is None:
            return None
        dx, dy = point.minus(self.primitive.start)
        return self.axis[0] * dx + self.axis[1] * dy

    def order(self, first: OneRootPoint2, first_key: OneRootNumber | None, second: OneRootPoint2, second_key: OneRootNumber | None) -> int:
        """Exact order of two points of the primitive along its travel direction."""
        if first_key is not None and second_key is not None:
            return first_key.compare(second_key)
        primitive = self.primitive
        assert isinstance(primitive, DeclaredArc)
        return compare_on_arc(primitive.centre, primitive.sense, self.start, first, second)

    def within(self, point: OneRootPoint2) -> bool:
        """Whether a point of the primitive's support lies within its extent."""
        primitive = self.primitive
        if isinstance(primitive, DeclaredLine):
            along = self.key(point)
            assert along is not None
            return along.sign() >= 0 and along.compare(primitive.length) <= 0
        return on_declared_arc(primitive, point)

    def add(self, low: OneRootPoint2, high: OneRootPoint2, low_witness: _Witness, high_witness: _Witness) -> None:
        self.reached.append(_Reach(low, high, self.key(low), self.key(high), low_witness, high_witness))

    def add_line_image(self, support: OffsetLineSupport, piece: CentreDomainPiece) -> None:
        """The contact image of a piece on this line's offset, clipped to the primitive."""
        assert isinstance(self.primitive, DeclaredLine)
        length = self.primitive.length
        first, last = support.parameter(piece.source), support.parameter(piece.target)  # = u·(contact - start)
        if last.sign() < 0 or first.compare(length) > 0:
            return
        low, low_key, low_witness = (support.contact(piece.source), first, (piece.source,)) if first.sign() >= 0 else (self.start, self.start_key, ())
        high, high_key, high_witness = (support.contact(piece.target), last, (piece.target,)) if last.compare(length) <= 0 else (self.end, self.end_key, ())
        self.reached.append(_Reach(low, high, low_key, high_key, low_witness, high_witness))

    def add_arc_image(self, support: OffsetCircleSupport, piece: CentreDomainPiece) -> None:
        """The contact image of a piece on this arc's offset circle, intersected with the arc (≤ 2 intervals)."""
        primitive = self.primitive
        assert isinstance(primitive, DeclaredArc)
        if piece.source == piece.target:
            self.reached.append(_Reach(self.start, self.end, self.start_key, self.end_key, (), ()))
            return
        first, last = support.contact(piece.source), support.contact(piece.target)
        first_in, last_in = self.within(first), self.within(last)
        if first_in:
            if last_in and self.order(first, self.key(first), last, self.key(last)) <= 0:
                self.add(first, last, (piece.source,), (piece.target,))
                return
            self.add(first, self.end, (piece.source,), ())
            if last_in:
                self.add(self.start, last, (), (piece.target,))
            return
        # The image enters the arc only through its start, if it passes it.
        if compare_on_arc(primitive.centre, primitive.sense, first, self.start, last) <= 0:
            self.add(self.start, last if last_in else self.end, (), (piece.target,) if last_in else ())

    def reaches_start(self) -> list[_Reach]:
        """The reached intervals that contain the primitive's start."""
        return [reach for reach in self.reached if self.order(reach.low, reach.low_key, self.start, self.start_key) == 0]

    def reaches_end(self) -> list[_Reach]:
        """The reached intervals that contain the primitive's end."""
        return [reach for reach in self.reached if self.order(reach.high, reach.high_key, self.end, self.end_key) == 0]

    def gaps(self, position: int) -> list[_Gap]:
        """The complement of the merged reached intervals in ``[start, end]``, with the bounding witnesses."""
        merged = self._merged()
        if not merged:
            return [_Gap(position, UnreachableInterval(self.feature, self.start, self.end, True, True), (), ())]
        gaps: list[_Gap] = []
        first, last = merged[0], merged[-1]
        if self.order(self.start, self.start_key, first.low, first.low_key) < 0:
            gaps.append(_Gap(position, UnreachableInterval(self.feature, self.start, first.low, True, False), (), first.low_witness))
        for before, after in itertools.pairwise(merged):
            gaps.append(_Gap(position, UnreachableInterval(self.feature, before.high, after.low, False, False), before.high_witness, after.low_witness))
        if self.order(last.high, last.high_key, self.end, self.end_key) < 0:
            gaps.append(_Gap(position, UnreachableInterval(self.feature, last.high, self.end, False, True), last.high_witness, ()))
        return gaps

    def _by_low(self, first: _Reach, second: _Reach) -> int:
        return self.order(first.low, first.low_key, second.low, second.low_key)

    def _merged(self) -> list[_Reach]:
        ordered = sorted(self.reached, key=functools.cmp_to_key(self._by_low))
        merged: list[_Reach] = []
        for reach in ordered:
            if not merged or self.order(reach.low, reach.low_key, merged[-1].high, merged[-1].high_key) > 0:
                merged.append(_Reach(reach.low, reach.high, reach.low_key, reach.high_key, reach.low_witness, reach.high_witness))
                continue
            last = merged[-1]
            if self.order(reach.low, reach.low_key, last.low, last.low_key) == 0:
                last.low_witness += reach.low_witness
            order = self.order(reach.high, reach.high_key, last.high, last.high_key)
            if order > 0:
                last.high, last.high_key, last.high_witness = reach.high, reach.high_key, reach.high_witness
            elif order == 0:
                last.high_witness += reach.high_witness
        return merged


def _minor_arc_chord(arc: DeclaredArc) -> tuple[Fraction, Fraction] | None:
    """``E - S`` when the arc sweeps at most π (its centre is not strictly on its own side of the chord), else ``None``."""
    chord = arc.end.minus(arc.start)
    to_centre = arc.centre.minus(arc.start)
    return chord if arc.sense * _cross_rational(chord, to_centre) >= 0 else None


def _attribute_unreachable(
    boundary: DeclaredBoundary,
    supports: CentreDomainSupports,
    components: tuple[CentreDomainComponent, ...],
    handle: Any,
) -> tuple[UnreachableAttribution, ...]:
    """Every unreached stretch of ∂P, exactly (see :class:`UnreachableAttribution`)."""
    primitives = boundary.primitives
    count = len(primitives)
    radius = supports.tool_radius
    tracks = [_PrimitiveTrack(primitive) for primitive in primitives]
    positions = {track.feature: position for position, track in enumerate(tracks)}
    owned_supports: tuple[CentreDomainSupport, ...] = (*supports.offset_lines, *supports.offset_arcs)
    owned = {support.owner: support for support in owned_supports}
    junction_hits: list[list[_Witness]] = [[] for _ in range(count)]
    pieces = [piece for cycle in _cycles(components) for piece in cycle]
    # 1. Contact-map images of owned pieces; junction-arc pieces reach their junction.
    for piece in pieces:
        for owner in piece.owners:
            if owner.kind is BoundaryFeatureKind.JUNCTION:
                junction_hits[owner.index].append((piece.source, piece.target))
                continue
            track, support = tracks[positions[owner]], owned[owner]
            if isinstance(support, OffsetLineSupport):
                track.add_line_image(support, piece)
            else:
                track.add_arc_image(support, piece)
    # 2. Points touched by vertex disks away from the owned images (three-way ties at a vertex).
    incident = [(piece.source, piece.owners | previous.owners) for cycle in _cycles(components) for previous, piece in zip((cycle[-1], *cycle[:-1]), cycle)]
    _add_vertex_touches(tracks, owned, incident, radius, junction_hits)
    # 3. An arc of ρ = r is wholly reached iff its centre (a point of ∂Ω* then) is in Ω*.
    for track in tracks:
        primitive = track.primitive
        if isinstance(primitive, DeclaredArc) and primitive.sense == 1 and primitive.radius == radius and native_locate(handle, primitive.centre) >= 0:
            centre = (OneRootPoint2.from_exact(primitive.centre),)
            track.reached.append(_Reach(track.start, track.end, track.start_key, track.end_key, centre, centre))
    # 4. A junction reached from either side, or by its own arc or a touch, is reached on both primitives.
    for junction in range(count):
        before, after = tracks[junction - 1], tracks[junction]
        witnesses = [point for hit in junction_hits[junction] for point in hit]
        reached = bool(junction_hits[junction])
        for reach in before.reaches_end():
            reached, witnesses = True, [*witnesses, *reach.high_witness]
        for reach in after.reaches_start():
            reached, witnesses = True, [*witnesses, *reach.low_witness]
        if reached:
            before.reached.append(_Reach(before.end, before.end, before.end_key, before.end_key, tuple(witnesses), tuple(witnesses)))
            after.reached.append(_Reach(after.start, after.start, after.start_key, after.start_key, tuple(witnesses), tuple(witnesses)))
    gaps = [gap for position, track in enumerate(tracks) for gap in track.gaps(position)]
    return tuple(_attributed(boundary, supports, stretch) for stretch in _stretches(boundary.case_name, gaps, count))


def _cycles(components: tuple[CentreDomainComponent, ...]) -> list[tuple[CentreDomainPiece, ...]]:
    return [cycle for component in components for cycle in (component.outer, *component.holes)]


def _add_vertex_touches(
    tracks: list[_PrimitiveTrack],
    owned: dict[BoundaryFeature, CentreDomainSupport],
    vertices: list[tuple[OneRootPoint2, frozenset[BoundaryFeature]]],
    radius: ExactMillimetre,
    junction_hits: list[list[_Witness]],
) -> None:
    """Boundary points at distance exactly r from a vertex w of ∂Ω* (w is at distance r from ∂P).

    Such a point is a foot b with w = b + r·n̂ (w on the primitive's offset, its
    contact within the primitive) or a junction at distance r. The feet on the
    owners of the two pieces meeting at w are their image ends already. Candidates are
    pruned by exact boxes in x-strips, then by exact dyadic-interval
    enclosures of the two equalities (a sound filter: an enclosure that
    excludes the value proves inequality); survivors are decided exactly.
    """
    count = len(tracks)
    strip = TOUCH_STRIP_TOOL_RADII * radius
    boxes = [_grown(declared_primitive_box(track.primitive), radius) for track in tracks]
    strips: dict[int, list[int]] = {}
    for position, box in enumerate(boxes):
        for key in range(math.floor(box[0] / strip), math.floor(box[2] / strip) + 1):
            strips.setdefault(key, []).append(position)
    squared_radius = radius * radius
    ends = [(_dyadic_point(track.primitive.start), _dyadic_point(track.primitive.end)) for track in tracks]
    for vertex, owners in vertices:
        x_low, x_high = _bracket(vertex.x)
        y_low, y_high = _bracket(vertex.y)
        candidates = {position for key in range(math.floor(x_low / strip), math.floor(x_high / strip) + 1) for position in strips.get(key, [])}
        for position in sorted(candidates):
            box = boxes[position]
            if box[0] > x_high or box[2] < x_low or box[1] > y_high or box[3] < y_low:
                continue
            track = tracks[position]
            support = owned.get(track.feature)
            if support is not None and track.feature not in owners and _may_lie_on(support, (x_low, x_high), (y_low, y_high)) and support.contains(vertex):
                contact = support.contact(vertex)
                if track.within(contact):
                    track.add(contact, contact, (vertex,), (vertex,))
            for junction, end, enclosed in ((position, track.primitive.start, ends[position][0]), ((position + 1) % count, track.primitive.end, ends[position][1])):
                if not _contains(_squared_distance_enclosure(enclosed, (x_low, x_high), (y_low, y_high)), squared_radius):
                    continue
                dx, dy = vertex.minus(end)
                if (dx * dx + dy * dy).compare(squared_radius) == 0:
                    junction_hits[junction].append((vertex,))


_Enclosure: TypeAlias = tuple[Fraction, Fraction]


def _contains(enclosure: _Enclosure, value: Fraction) -> bool:
    return enclosure[0] <= value <= enclosure[1]


def _product(left: _Enclosure, right: _Enclosure) -> _Enclosure:
    products = (left[0] * right[0], left[0] * right[1], left[1] * right[0], left[1] * right[1])
    return min(products), max(products)


def _dyadic(value: Fraction) -> _Enclosure:
    scale = 1 << VERTEX_BRACKET_BITS
    return Fraction(math.floor(value * scale), scale), Fraction(math.ceil(value * scale), scale)


def _dyadic_point(point: ExactPoint2[WorldXY]) -> tuple[_Enclosure, _Enclosure]:
    """Outward dyadic enclosure of a rational point, so the prefilter works on short rationals."""
    return _dyadic(point.x), _dyadic(point.y)


def _squared_distance_enclosure(centre: tuple[_Enclosure, _Enclosure], x: _Enclosure, y: _Enclosure) -> _Enclosure:
    """Exact enclosure of ``|w - c|²`` for ``w`` in the box ``x × y`` and ``c`` in the box ``centre``."""
    total = (Fraction(0), Fraction(0))
    for low, high in ((x[0] - centre[0][1], x[1] - centre[0][0]), (y[0] - centre[1][1], y[1] - centre[1][0])):
        square = (Fraction(0) if low <= 0 <= high else min(low * low, high * high), max(low * low, high * high))
        total = (total[0] + square[0], total[1] + square[1])
    return total


def _may_lie_on(support: CentreDomainSupport, x: _Enclosure, y: _Enclosure) -> bool:
    """False only if no point of the box ``x × y`` lies on the support (exact enclosure)."""
    if isinstance(support, OffsetLineSupport):
        along_x, along_y = _product((support.normal[0], support.normal[0]), x), _product((support.normal[1], support.normal[1]), y)
        return along_x[0] + along_y[0] <= support.level <= along_x[1] + along_y[1]
    return _contains(_squared_distance_enclosure(_dyadic_point(support.centre), x, y), support.radius_squared)


def _grown(box: PrimitiveBox, margin: Fraction) -> PrimitiveBox:
    return box[0] - margin, box[1] - margin, box[2] + margin, box[3] + margin


def _bracket(value: OneRootNumber) -> tuple[Fraction, Fraction]:
    """Exact dyadic ``[low, high]`` (denominator 2^VERTEX_BRACKET_BITS) around a one-root number; prunes only.

    The square root is bracketed by integer square roots and the ends rounded
    outward to the dyadic grid, so later enclosure arithmetic stays on short
    rationals.
    """
    scale = 1 << VERTEX_BRACKET_BITS
    if value.is_rational:
        low = high = value.rational
    else:
        radicand = value.radicand
        root_low = Fraction(math.isqrt(radicand.numerator * radicand.denominator * scale * scale), radicand.denominator * scale)
        root_high = root_low + Fraction(1, radicand.denominator * scale)
        ends = (value.rational + value.coefficient * root_low, value.rational + value.coefficient * root_high)
        low, high = min(ends), max(ends)
    return Fraction(math.floor(low * scale), scale), Fraction(math.ceil(high * scale), scale)


def _stretches(case_name: str, gaps: list[_Gap], count: int) -> list[list[_Gap]]:
    """Group the gaps into maximal runs connected through unreached junctions, in boundary order."""
    if not gaps:
        return []
    size = len(gaps)
    links = []
    for index, gap in enumerate(gaps):
        following = gaps[(index + 1) % size]
        if not gap.interval.end_closed:
            links.append(False)
            continue
        if following.position != (gap.position + 1) % count or not following.interval.start_closed:
            raise InconsistentUnreachableAttributionError(f"{case_name}: junction {(gap.position + 1) % count} is unreached on primitive {gap.position} only.")
        links.append(True)
    if all(links):
        return [gaps]
    first = next(index for index in range(size) if not links[index - 1])
    runs: list[list[int]] = []
    for offset in range(size):
        index = (first + offset) % size
        if not links[index - 1]:
            runs.append([])
        runs[-1].append(index)
    return [[gaps[index] for index in run] for run in sorted(runs)]


def _attributed(boundary: DeclaredBoundary, supports: CentreDomainSupports, stretch: list[_Gap]) -> UnreachableAttribution:
    primitives = boundary.primitives
    count = len(primitives)
    junctions = [(gap.position + 1) % count for gap in stretch if gap.interval.end_closed]
    kinks = {_feature(junction, False) for junction in junctions if _cross_rational(primitives[junction - 1].end_tangent, primitives[junction].start_tangent) > 0}
    tight = {
        gap.interval.primitive
        for gap in stretch
        if isinstance(primitive := primitives[gap.position], DeclaredArc) and primitive.sense == 1 and primitive.radius < supports.tool_radius
    }
    start_witness, end_witness = stretch[0].start_witness, stretch[-1].end_witness
    local = any(first == second for first in start_witness for second in end_witness)
    intervals = tuple(gap.interval for gap in stretch)
    if local and tight:
        return UnreachableAttribution(UnreachableCause.TIGHT_CONVEX_ARC, frozenset(tight | kinks), intervals)
    if local and kinks:
        return UnreachableAttribution(UnreachableCause.CONVEX_KINK_CUSP, frozenset(kinks), intervals)
    return UnreachableAttribution(UnreachableCause.NECK, frozenset(tight | kinks | {interval.primitive for interval in intervals}), intervals)


def _feature(index: int, is_arc: bool) -> BoundaryFeature:
    return BoundaryFeature(BoundaryFeatureKind.ARC if is_arc else BoundaryFeatureKind.JUNCTION, index)


def _cross_rational(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def construct_centre_domain_arrangement(boundary: DeclaredBoundary, supports: CentreDomainSupports) -> CentreDomainArrangement:
    """Native seam: Ω* = P ∖ join(line bands, arc bands, junction disks) in Gps circle-segment traits.

    The component count is checked before anything queries the native handle
    (an empty native domain must never be queried). Every native piece is
    matched to the rational supports by an exact support key; the owner set is
    every support on that key it travels along.

    Raises:
        EmptyCentreDomainError: If Ω* is empty.
        NeckSplitsCentreDomainError: If Ω* has more than one component.
        UnmatchedCentreDomainPieceError: If a piece lies on no support or travels against every one.
    """
    native = build_native_centre_domain(boundary, supports.tool_radius)
    count = len(native.components)
    if count == 0:
        raise EmptyCentreDomainError(f"{boundary.case_name}: the tool fits nowhere.")
    if count > 1:
        raise NeckSplitsCentreDomainError(f"{boundary.case_name}: a neck of width <= 2r splits Ω* into {count} components.", count)
    index = _support_index(supports)
    components: list[CentreDomainComponent] = []
    for component_index, cycles in enumerate(native.components):
        matched = tuple(
            tuple(_matched_piece(piece, (component_index, ccb_index, piece_index), index) for piece_index, piece in enumerate(cycle)) for ccb_index, cycle in enumerate(cycles)
        )
        components.append(CentreDomainComponent(matched[0], matched[1:]))
    frozen = tuple(components)
    unreachable = _attribute_unreachable(boundary, supports, frozen, native.handle)
    return CentreDomainArrangement(frozen, unreachable, _vertex_disks(frozen, supports.tool_radius), native.handle)


# ---------------------------------------------------------------------------
# The domain
# ---------------------------------------------------------------------------

_LOCATIONS = {1: PointLocation.INTERIOR, 0: PointLocation.BOUNDARY, -1: PointLocation.EXTERIOR}


@dataclass(frozen=True)
class CentreDomain:
    """Ω* = {x ∈ P : dist(x, ∂P) >= r}, exact with respect to the declared boundary."""

    supports: CentreDomainSupports
    components: tuple[CentreDomainComponent, ...]
    unreachable: tuple[UnreachableAttribution, ...]
    vertex_disks: tuple[VertexDisk, ...]
    identity: CentreDomainIdentity
    native: Any = field(compare=False, repr=False)

    @classmethod
    def build(cls, boundary: DeclaredBoundary, tool_radius: Fraction) -> Self:
        """Construct Ω* exactly.

        Args:
            boundary: The declared counter-clockwise boundary.
            tool_radius: Tool radius r, a positive binary64 injected exactly.

        Returns:
            The single-component centre domain.

        Raises:
            NonSimpleDeclaredBoundaryError: If the declared cycle is not simple.
            ClockwiseDeclaredBoundaryError: If the declared cycle is clockwise.
            EmptyCentreDomainError: If Ω* is empty.
            NeckSplitsCentreDomainError: If Ω* has more than one component.
            UnmatchedCentreDomainPieceError: If a native piece matches no support.
        """
        supports = CentreDomainSupports.build(boundary, tool_radius)
        arrangement = construct_centre_domain_arrangement(boundary, supports)
        identity = CentreDomainIdentity(CENTRE_DOMAIN_CONTRACT_VERSION, NATIVE_CENTRE_DOMAIN_COMPONENT, boundary.identity, supports.tool_radius)
        return cls(supports, arrangement.components, arrangement.unreachable, arrangement.vertex_disks, identity, arrangement.native)

    @functools.cached_property
    def pieces(self) -> tuple[CentreDomainPiece, ...]:
        """Every piece of every CCB, in component, CCB and CCW order."""
        return tuple(piece for component in self.components for ccb in (component.outer, *component.holes) for piece in ccb)

    def locate(self, point: ExactPoint2[WorldXY]) -> PointLocation:
        """Exact location of a rational point (interior, on ∂Ω*, exterior), decided by the native arrangement."""
        return _LOCATIONS[native_locate(self.native, point)]

    def squared_distance_to(self, point: ExactPoint2[WorldXY]) -> OneRootNumber:
        """Exact squared distance from a rational point to Ω* (zero inside); one-root in ``|p - c|²``."""
        if self.locate(point) is not PointLocation.EXTERIOR:
            return OneRootNumber.of(0)
        return native_nearest(self.native, point).squared_distance

    def project_contact(self, point: ExactPoint2[WorldXY], gap: ExactMillimetre) -> CentreDomainSite:
        """Project a hypothesis contact to its nearest ∂Ω* site, rejecting beyond ``gap``.

        Native and exact: squared distances are one-root numbers (to an offset
        circle ``(|p - c| - R)²``) compared across roots by squaring; a
        certified box per piece only prunes. A vertex is reported at the source
        of the piece that starts there, so a vertex has one site.

        Raises:
            ContactBeyondGapError: If the nearest site is farther than ``gap``.
            AmbiguousContactProjectionError: If two distinct points are exactly nearest.
        """
        nearest = native_nearest(self.native, point)
        if nearest.squared_distance.compare(gap * gap) > 0:
            raise ContactBeyondGapError(f"Nearest ∂Ω* site of {point} is farther than the gap {gap}.")
        if nearest.ambiguous:
            raise AmbiguousContactProjectionError(f"Distinct ∂Ω* points are exactly nearest to {point}.")
        return CentreDomainSite.build(self.pieces[nearest.piece], nearest.point)

    def decide_contacts(self, points: list[tuple[float, float]], gap: ExactMillimetre) -> ContactDecisions:
        """Decide, in one native call, the piece each binary64 contact projects onto; ``None`` marks beyond ``gap``.

        The decision is the same exact decision as :meth:`project_contact`
        (vertices reported on the piece they start). Sites are lazy: each
        :class:`DecidedContact` builds its exact one-root site only when asked.

        Raises:
            AmbiguousContactProjectionError: If any contact has two distinct nearest points.
        """
        pieces = self.pieces
        batch = native_decide_binary64(self.native, points, gap)
        decided: list[DecidedContact | None] = []
        for index, decision in enumerate(batch.decisions):
            if decision.status == BEYOND_GAP:
                decided.append(None)
                continue
            if decision.status == AMBIGUOUS:
                raise AmbiguousContactProjectionError(f"Distinct ∂Ω* points are exactly nearest to contact {index} {points[index]}.")
            decided.append(DecidedContact(points[index], pieces[decision.piece], decision.is_vertex))
        return ContactDecisions(tuple(decided), batch.exact_fallbacks)

    def project_contacts(self, points: list[tuple[float, float]], gap: ExactMillimetre) -> tuple[DecidedContact | None, ...]:
        """Project binary64 hypothesis contacts: one native call for all, exact sites materialised lazily.

        Raises:
            AmbiguousContactProjectionError: If any contact has two distinct nearest points.
        """
        return self.decide_contacts(points, gap).contacts
