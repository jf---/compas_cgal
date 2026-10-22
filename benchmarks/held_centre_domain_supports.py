"""Rational supports of ∂Ω* for the Held centre domain, from the declared boundary alone.

Every support of ∂Ω* (offset lines, offset circles of radius ``ρ ± r``,
r-circles at reflex junctions), the raw offset endpoints, the local one-root
vertex at every convex junction, and the contact map. Precondition, checked
exactly by :func:`validate_simple_counter_clockwise`: the declared cycle is
simple and counter-clockwise, so the interior is on the left of every
primitive. Design: docs/held_exact_centre_domain.md.
"""

from __future__ import annotations

import enum
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction
from typing import Self
from typing import TypeAlias

from benchmarks.held_centre_domain_errors import CentreDomainError
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredBoundaryIdentity
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import DeclaredPrimitive
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import Sense
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import SquaredMillimetre
from compas_cgal.units import WorldXY

_Vector: TypeAlias = tuple[Fraction, Fraction]
# Exact axis-aligned box (xmin, ymin, xmax, ymax).
PrimitiveBox: TypeAlias = tuple[Fraction, Fraction, Fraction, Fraction]
# The four axis directions, in which an arc may reach beyond its endpoint box.
_AXES: tuple[_Vector, ...] = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)), (Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1)))


class InvalidCentreDomainToolRadiusError(CentreDomainError):
    """The tool radius is not a positive binary64."""


class ReversedJunctionError(CentreDomainError):
    """A junction reverses direction exactly (a cusp): the offset is undefined there."""


class UnresolvedJunctionVertexError(CentreDomainError):
    """The two supports at a convex junction do not cross exactly once in the convex orientation."""


class NoSupportError(CentreDomainError):
    """A boundary feature has no support on ∂Ω* (a tight convex arc, or a non-reflex junction)."""


class NoJunctionVertexError(CentreDomainError):
    """A junction has no local vertex (it is reflex, or its vertex is deferred to the native arrangement)."""


class ClockwiseDeclaredBoundaryError(CentreDomainError):
    """The declared cycle runs clockwise: the interior would be on the right of every primitive."""


class NonSimpleDeclaredBoundaryError(CentreDomainError):
    """Two primitives of the declared cycle cross, touch or overlap away from their shared junction."""


# ---------------------------------------------------------------------------
# Supports
# ---------------------------------------------------------------------------


class PieceKind(enum.Enum):
    """What a ∂Ω* piece offsets."""

    OFFSET_LINE = "offset_line"
    OFFSET_ARC = "offset_arc"
    JUNCTION_ARC = "junction_arc"


class BoundaryFeatureKind(enum.Enum):
    """A feature of the declared boundary: a primitive or a junction."""

    LINE = "line"
    ARC = "arc"
    JUNCTION = "junction"


@dataclass(frozen=True)
class BoundaryFeature:
    """Primitive ``index`` (LINE / ARC) or junction ``index`` (the start of primitive ``index``)."""

    kind: BoundaryFeatureKind
    index: int


_KIND_ORDER: dict[BoundaryFeatureKind, int] = {kind: position for position, kind in enumerate(BoundaryFeatureKind)}


def owner_order_key(feature: BoundaryFeature) -> tuple[int, int]:
    """The exact canonical order of boundary features: kind in declaration order (LINE, ARC, JUNCTION), then index."""
    return _KIND_ORDER[feature.kind], feature.index


class OwnerSet(frozenset[BoundaryFeature]):
    """A set of boundary features that iterates and prints in canonical order (:func:`owner_order_key`).

    A plain ``frozenset`` iterates and prints in hash-table order, and the
    enum member hashes its name with the per-process string hash seed, so
    digests of identical Ω* paths differed across processes and any output
    built by iterating owners could too. Set semantics (equality with any
    set, ``in``, hashing) are the frozenset's; only order is fixed. Build it
    with :func:`canonical_owners`; set operators return plain frozensets, so
    a union is canonicalised again by passing it through that factory.
    """

    __slots__ = ()

    def __iter__(self) -> Iterator[BoundaryFeature]:
        return iter(sorted(super().__iter__(), key=owner_order_key))

    def __repr__(self) -> str:
        return f"OwnerSet({tuple(self)!r})"


def canonical_owners(features: Iterable[BoundaryFeature]) -> OwnerSet:
    """The distinct ``features`` as an :class:`OwnerSet` (canonical iteration and repr)."""
    return OwnerSet(features)


@dataclass(frozen=True)
class OffsetLineSupport:
    """The rational line ``n̂·x = level`` at distance r inside declared line ``owner``.

    Travel is along the exactly unit direction ``u``; ``n̂ = u^⊥`` is the interior
    normal; ``raw_source = anchor`` and ``raw_target = anchor + length·u`` are the
    offsets of the primitive's ends before any trimming.
    """

    owner: BoundaryFeature
    anchor: ExactPoint2[WorldXY]
    direction: _Vector
    normal: _Vector
    level: ExactMillimetre
    length: ExactMillimetre
    tool_radius: ExactMillimetre

    @property
    def kind(self) -> PieceKind:
        """Always :attr:`PieceKind.OFFSET_LINE`."""
        return PieceKind.OFFSET_LINE

    @property
    def raw_source(self) -> ExactPoint2[WorldXY]:
        """Offset of the primitive's start."""
        return self.anchor

    @property
    def raw_target(self) -> ExactPoint2[WorldXY]:
        """Offset of the primitive's end."""
        return self.anchor.offset(self.direction, self.length)

    def contains(self, point: OneRootPoint2) -> bool:
        """Whether ``point`` is exactly on the support line."""
        return (self.normal[0] * point.x + self.normal[1] * point.y).compare(self.level) == 0

    def parameter(self, point: OneRootPoint2) -> OneRootNumber:
        """Exact position ``u·(point - anchor)`` along the line (same extension as ``point``)."""
        dx, dy = point.minus(self.anchor)
        return self.direction[0] * dx + self.direction[1] * dy

    def travel_tangent(self, point: OneRootPoint2) -> tuple[OneRootNumber, OneRootNumber]:
        """Travel direction at ``point`` (constant ``u``)."""
        return OneRootNumber.of(self.direction[0]), OneRootNumber.of(self.direction[1])

    def contact(self, centre: OneRootPoint2) -> OneRootPoint2:
        """Boundary contact ``q - r·n̂`` of a centre ``q`` on this support (adds no root)."""
        return OneRootPoint2.build(centre.x - self.tool_radius * self.normal[0], centre.y - self.tool_radius * self.normal[1])


@dataclass(frozen=True)
class OffsetCircleSupport:
    """A rational circle of ∂Ω*: an offset arc (radius ``ρ ∓ r``) or a junction r-circle.

    ``boundary_radius`` is the owner's radius ``ρ`` for an offset arc and ``0`` for
    a junction (the junction is a degenerate circle), so the contact map
    ``c + (ρ / R)(q - c)`` covers both. ``sense`` is the travel sense along ∂Ω*
    (``1`` counter-clockwise); a junction arc is always travelled clockwise.
    """

    owner: BoundaryFeature
    centre: ExactPoint2[WorldXY]
    radius: ExactMillimetre
    boundary_radius: ExactMillimetre
    sense: Sense
    raw_source: ExactPoint2[WorldXY]
    raw_target: ExactPoint2[WorldXY]

    @property
    def kind(self) -> PieceKind:
        """:attr:`PieceKind.JUNCTION_ARC` for a junction owner, else :attr:`PieceKind.OFFSET_ARC`."""
        return PieceKind.JUNCTION_ARC if self.owner.kind is BoundaryFeatureKind.JUNCTION else PieceKind.OFFSET_ARC

    @property
    def radius_squared(self) -> SquaredMillimetre:
        """Exact squared radius."""
        return SquaredMillimetre(self.radius * self.radius)

    def contains(self, point: OneRootPoint2) -> bool:
        """Whether ``point`` is exactly on the support circle."""
        dx, dy = point.minus(self.centre)
        return (dx * dx + dy * dy).compare(self.radius_squared) == 0

    def travel_tangent(self, point: OneRootPoint2) -> tuple[OneRootNumber, OneRootNumber]:
        """Travel direction at ``point``, scaled by the radius: ``σ·(point - c)^⊥``."""
        dx, dy = point.minus(self.centre)
        return -self.sense * dy, self.sense * dx

    def contact(self, centre: OneRootPoint2) -> OneRootPoint2:
        """Boundary contact ``c + (ρ / R)(q - c)`` of a centre ``q`` on this support (adds no root)."""
        ratio = self.boundary_radius / self.radius
        dx, dy = centre.minus(self.centre)
        return OneRootPoint2.build(ratio * dx + self.centre.x, ratio * dy + self.centre.y)


CentreDomainSupport: TypeAlias = "OffsetLineSupport | OffsetCircleSupport"


class JunctionVertexKind(enum.Enum):
    """How the two supports meet at a junction."""

    CONVEX = "convex"
    G1 = "g1"


@dataclass(frozen=True)
class JunctionVertex:
    """The local vertex of the offset at a convex or exactly-G1 junction.

    It is the crossing of the incoming and outgoing supports in the convex
    orientation. It is a vertex of Ω* unless the global arrangement trims it.
    """

    junction: int
    kind: JunctionVertexKind
    point: OneRootPoint2
    incoming: BoundaryFeature
    outgoing: BoundaryFeature


@dataclass(frozen=True)
class CentreDomainSupports:
    """Every rational support of ∂Ω* and the local vertices, from the declared boundary alone.

    ``deferred_junctions`` are convex junctions next to a primitive without a
    support (an arc tighter than the tool); their vertex is global and native.
    """

    boundary: DeclaredBoundaryIdentity
    tool_radius: ExactMillimetre
    offset_lines: tuple[OffsetLineSupport, ...]
    offset_arcs: tuple[OffsetCircleSupport, ...]
    junction_arcs: tuple[OffsetCircleSupport, ...]
    tight_convex_arcs: tuple[int, ...]
    junction_vertices: tuple[JunctionVertex, ...]
    deferred_junctions: tuple[int, ...]

    @classmethod
    def build(cls, boundary: DeclaredBoundary, tool_radius: Fraction) -> Self:
        """Derive every support and local vertex exactly.

        Args:
            boundary: The declared counter-clockwise boundary.
            tool_radius: Tool radius r, a positive binary64 injected exactly.

        Returns:
            The supports.

        Raises:
            InvalidCentreDomainToolRadiusError: If r is not a positive binary64.
            NonSimpleDeclaredBoundaryError: If the cycle crosses, touches or overlaps itself.
            ClockwiseDeclaredBoundaryError: If the cycle is clockwise.
            ReversedJunctionError: If a junction reverses direction exactly.
            UnresolvedJunctionVertexError: If a convex junction's supports do not cross as required.
        """
        radius = _validated_tool_radius(tool_radius)
        validate_simple_counter_clockwise(boundary)
        primitives = boundary.primitives
        supports: dict[int, CentreDomainSupport] = {}
        lines: list[OffsetLineSupport] = []
        arcs: list[OffsetCircleSupport] = []
        tight: list[int] = []
        for primitive in primitives:
            if isinstance(primitive, DeclaredLine):
                line = _offset_line(primitive, radius)
                supports[primitive.source_index] = line
                lines.append(line)
                continue
            arc = _offset_arc(primitive, radius)
            if arc is None:
                tight.append(primitive.source_index)
                continue
            supports[primitive.source_index] = arc
            arcs.append(arc)
        junction_arcs: list[OffsetCircleSupport] = []
        vertices: list[JunctionVertex] = []
        deferred: list[int] = []
        for index, outgoing in enumerate(primitives):
            incoming = primitives[index - 1]
            junction = outgoing.start
            incoming_tangent, outgoing_tangent = incoming.end_tangent, outgoing.start_tangent
            turn = _cross(incoming_tangent, outgoing_tangent)
            if turn == 0 and _dot(incoming_tangent, outgoing_tangent) < 0:
                raise ReversedJunctionError(f"{boundary.case_name}: junction {index} reverses direction.")
            incoming_normal, outgoing_normal = _left(incoming_tangent), _left(outgoing_tangent)
            if turn < 0:
                junction_arcs.append(
                    OffsetCircleSupport(
                        BoundaryFeature(BoundaryFeatureKind.JUNCTION, index),
                        junction,
                        radius,
                        ExactMillimetre(Fraction(0)),
                        -1,
                        junction.offset(incoming_normal, radius),
                        junction.offset(outgoing_normal, radius),
                    )
                )
                continue
            incoming_support = supports.get(incoming.source_index)
            outgoing_support = supports.get(outgoing.source_index)
            if incoming_support is None or outgoing_support is None:
                deferred.append(index)
                continue
            if turn == 0:
                point = OneRootPoint2.from_exact(junction.offset(incoming_normal, radius))
                vertices.append(JunctionVertex(index, JunctionVertexKind.G1, point, incoming_support.owner, outgoing_support.owner))
                continue
            point = _convex_crossing(boundary.case_name, index, incoming_support, outgoing_support)
            vertices.append(JunctionVertex(index, JunctionVertexKind.CONVEX, point, incoming_support.owner, outgoing_support.owner))
        return cls(
            boundary.identity,
            radius,
            tuple(lines),
            tuple(arcs),
            tuple(junction_arcs),
            tuple(tight),
            tuple(vertices),
            tuple(deferred),
        )

    def support_of(self, owner: BoundaryFeature) -> CentreDomainSupport:
        """The support offsetting ``owner``.

        Raises:
            NoSupportError: If ``owner`` has no support.
        """
        candidates: tuple[CentreDomainSupport, ...] = (*self.offset_lines, *self.offset_arcs, *self.junction_arcs)
        for support in candidates:
            if support.owner == owner:
                return support
        raise NoSupportError(f"{owner} has no support on ∂Ω*.")

    def junction_vertex(self, junction: int) -> JunctionVertex:
        """The local vertex at ``junction``.

        Raises:
            NoJunctionVertexError: If the junction is reflex or deferred.
        """
        for vertex in self.junction_vertices:
            if vertex.junction == junction:
                return vertex
        raise NoJunctionVertexError(f"Junction {junction} has no local vertex.")


def _validated_tool_radius(tool_radius: Fraction) -> ExactMillimetre:
    radius = Fraction(tool_radius)
    if radius <= 0 or Fraction(float(radius)) != radius:
        raise InvalidCentreDomainToolRadiusError(f"Tool radius {radius} is not a positive binary64.")
    return ExactMillimetre(radius)


def _offset_line(line: DeclaredLine, radius: ExactMillimetre) -> OffsetLineSupport:
    direction = line.direction.vector
    normal = _left(direction)
    anchor = line.start.offset(normal, radius)
    level = ExactMillimetre(normal[0] * anchor.x + normal[1] * anchor.y)
    return OffsetLineSupport(BoundaryFeature(BoundaryFeatureKind.LINE, line.source_index), anchor, direction, normal, level, line.length, radius)


def _offset_arc(arc: DeclaredArc, radius: ExactMillimetre) -> OffsetCircleSupport | None:
    offset_radius = arc.radius - arc.sense * radius
    if offset_radius <= 0:
        return None
    return OffsetCircleSupport(
        BoundaryFeature(BoundaryFeatureKind.ARC, arc.source_index),
        arc.centre,
        ExactMillimetre(offset_radius),
        arc.radius,
        arc.sense,
        arc.start.offset(_left(arc.start_tangent), radius),
        arc.end.offset(_left(arc.end_tangent), radius),
    )


def _convex_crossing(case_name: str, junction: int, incoming: CentreDomainSupport, outgoing: CentreDomainSupport) -> OneRootPoint2:
    candidates = _support_intersections(incoming, outgoing)
    convex = [point for point in candidates if _cross_one_root(incoming.travel_tangent(point), outgoing.travel_tangent(point)) > 0]
    if len(convex) != 1:
        raise UnresolvedJunctionVertexError(f"{case_name}: junction {junction} has {len(convex)} convex crossings among {len(candidates)} intersections.")
    return convex[0]


def _support_intersections(first: CentreDomainSupport, second: CentreDomainSupport) -> tuple[OneRootPoint2, ...]:
    if isinstance(first, OffsetLineSupport) and isinstance(second, OffsetLineSupport):
        return _line_line(first, second)
    if isinstance(first, OffsetLineSupport) and isinstance(second, OffsetCircleSupport):
        return _line_circle(first.anchor, first.direction, second.centre, second.radius_squared)
    if isinstance(first, OffsetCircleSupport) and isinstance(second, OffsetLineSupport):
        return _line_circle(second.anchor, second.direction, first.centre, first.radius_squared)
    assert isinstance(first, OffsetCircleSupport) and isinstance(second, OffsetCircleSupport)
    return _circle_circle(first, second)


def _line_line(first: OffsetLineSupport, second: OffsetLineSupport) -> tuple[OneRootPoint2, ...]:
    determinant = _cross(first.normal, second.normal)
    if determinant == 0:
        return ()
    x = (first.level * second.normal[1] - second.level * first.normal[1]) / determinant
    y = (first.normal[0] * second.level - second.normal[0] * first.level) / determinant
    return (OneRootPoint2.rational(x, y),)


def _line_circle(anchor: ExactPoint2[WorldXY], direction: _Vector, centre: ExactPoint2[WorldXY], radius_squared: Fraction) -> tuple[OneRootPoint2, ...]:
    """Points ``anchor + λ·d`` on the circle: ``|d|²λ² + 2Bλ + C = 0``."""
    relative = anchor.minus(centre)
    a = _dot(direction, direction)
    b = _dot(direction, relative)
    c = _dot(relative, relative) - radius_squared
    discriminant = b * b - a * c
    if discriminant < 0:
        return ()
    base = (anchor.x - b / a * direction[0], anchor.y - b / a * direction[1])
    scale = 1 / a
    points = []
    for branch in (1, -1) if discriminant != 0 else (1,):
        x = OneRootNumber.build(base[0], branch * scale * direction[0], discriminant)
        y = OneRootNumber.build(base[1], branch * scale * direction[1], discriminant)
        points.append(OneRootPoint2.build(x, y))
    return tuple(points)


def _circle_circle(first: OffsetCircleSupport, second: OffsetCircleSupport) -> tuple[OneRootPoint2, ...]:
    return _circle_circle_points(first.centre, first.radius_squared, second.centre, second.radius_squared)


def _circle_circle_points(first: ExactPoint2[WorldXY], first_squared: Fraction, second: ExactPoint2[WorldXY], second_squared: Fraction) -> tuple[OneRootPoint2, ...]:
    """Radical line ``2(c2 - c1)·x = |c2|² - |c1|² + R1² - R2²`` met with the first circle (none if concentric)."""
    d = second.minus(first)
    if d[0] == 0 and d[1] == 0:
        return ()
    level = (_dot((second.x, second.y), (second.x, second.y)) - _dot((first.x, first.y), (first.x, first.y))) + (first_squared - second_squared)
    t = (level - 2 * _dot(d, (first.x, first.y))) / (2 * _dot(d, d))
    foot = first.offset(d, t)
    return _line_circle(foot, _left(d), first, first_squared)


def _cross_one_root(left: tuple[OneRootNumber, OneRootNumber], right: tuple[OneRootNumber, OneRootNumber]) -> int:
    return (left[0] * right[1] - left[1] * right[0]).sign()


def _left(vector: _Vector) -> _Vector:
    return -vector[1], vector[0]


def _dot(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _cross(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


# ---------------------------------------------------------------------------
# Precondition: a simple counter-clockwise cycle (exact)
# ---------------------------------------------------------------------------


def declared_primitive_box(primitive: DeclaredPrimitive) -> PrimitiveBox:
    """Exact bounding box of a declared primitive: its ends, and an arc's axis extremes inside its sweep."""
    xs, ys = [primitive.start.x, primitive.end.x], [primitive.start.y, primitive.end.y]
    if isinstance(primitive, DeclaredArc):
        for axis in _AXES:
            extreme = primitive.centre.offset(axis, primitive.radius)
            if _on_declared_arc_rational(primitive, extreme):
                xs.append(extreme.x)
                ys.append(extreme.y)
    return min(xs), min(ys), max(xs), max(ys)


def validate_simple_counter_clockwise(boundary: DeclaredBoundary) -> None:
    """Refuse a declared cycle that is not simple, or not counter-clockwise, exactly.

    Simplicity: every pair of primitives whose exact boxes meet shares no point
    except the junction of adjacent primitives (crossings, touches and
    collinear or co-circular overlaps are refused). Orientation: at the
    lexicographically lowest point of the cycle (min y, then min x) a simple
    cycle turns left iff it is counter-clockwise; there an arc's own sense, or
    the exact turn ``cross(T_in, T_out)`` of the junction, decides it (a G1
    junction at the bottom is CCW iff it travels towards +x).

    Raises:
        ReversedJunctionError: If a junction reverses direction exactly.
        NonSimpleDeclaredBoundaryError: If two primitives meet away from a shared junction.
        ClockwiseDeclaredBoundaryError: If the cycle is clockwise.
    """
    primitives = boundary.primitives
    for index, outgoing in enumerate(primitives):
        incoming = primitives[index - 1]
        if _cross(incoming.end_tangent, outgoing.start_tangent) == 0 and _dot(incoming.end_tangent, outgoing.start_tangent) < 0:
            raise ReversedJunctionError(f"{boundary.case_name}: junction {index} reverses direction.")
    _validate_simple(boundary)
    if not _counter_clockwise(primitives):
        raise ClockwiseDeclaredBoundaryError(f"{boundary.case_name}: the declared cycle is clockwise.")


def _validate_simple(boundary: DeclaredBoundary) -> None:
    primitives = boundary.primitives
    count = len(primitives)
    boxes = [declared_primitive_box(primitive) for primitive in primitives]
    order = sorted(range(count), key=lambda index: boxes[index][0])
    for position, first in enumerate(order):
        for second in order[position + 1 :]:
            if boxes[second][0] > boxes[first][2]:
                break
            if boxes[second][1] > boxes[first][3] or boxes[first][1] > boxes[second][3]:
                continue
            low, high = min(first, second), max(first, second)
            shared = [primitives[high].start] if high == low + 1 else []
            if low == 0 and high == count - 1:
                shared.append(primitives[0].start)
            allowed = [OneRootPoint2.from_exact(junction) for junction in shared]
            if len(shared) == 1:
                common = _common_points_through(primitives[low], primitives[high], shared[0])
            else:
                common = _common_points(primitives[low], primitives[high])
            if common is None:
                raise NonSimpleDeclaredBoundaryError(f"{boundary.case_name}: primitives {low} and {high} overlap.")
            for point in common:
                if not any(point == junction for junction in allowed):
                    raise NonSimpleDeclaredBoundaryError(f"{boundary.case_name}: primitives {low} and {high} meet at {point.as_floats()}.")


def _common_points_through(first: DeclaredPrimitive, second: DeclaredPrimitive, shared: ExactPoint2[WorldXY]) -> tuple[OneRootPoint2, ...] | None:
    """:func:`_common_points` of two primitives through the rational point ``shared``, in rationals only.

    The second intersection of a line through a rational point of a rational
    circle is rational (``J + t·u``, ``t = -2u·(J - c)``), and so is the second
    intersection of two circles through J (J reflected in their centre line),
    so adjacent primitives need no square root.
    """
    if isinstance(first, DeclaredLine) and isinstance(second, DeclaredLine):
        if _cross(first.direction.vector, second.direction.vector) != 0:
            return (OneRootPoint2.from_exact(shared),)
        return _line_line_common(first, second)
    if isinstance(first, DeclaredArc) and isinstance(second, DeclaredLine):
        first, second = second, first
    candidates: list[ExactPoint2[WorldXY]] = [shared]
    if isinstance(first, DeclaredLine):
        assert isinstance(second, DeclaredArc)
        u = first.direction.vector
        along = -2 * _dot(u, shared.minus(second.centre))
        other = shared.offset(u, along)
        if other != shared and _on_line(first, OneRootPoint2.from_exact(other)) and _on_declared_arc_rational(second, other):
            candidates.append(other)
        return tuple(OneRootPoint2.from_exact(point) for point in candidates)
    assert isinstance(second, DeclaredArc)
    if first.centre == second.centre:
        return _arc_arc_common(first, second)
    axis = second.centre.minus(first.centre)
    relative = shared.minus(first.centre)
    foot = _dot(relative, axis) / _dot(axis, axis)
    other = ExactPoint2[WorldXY].build(2 * (first.centre.x + foot * axis[0]) - shared.x, 2 * (first.centre.y + foot * axis[1]) - shared.y)
    if other != shared and _on_declared_arc_rational(first, other) and _on_declared_arc_rational(second, other):
        candidates.append(other)
    return tuple(OneRootPoint2.from_exact(point) for point in candidates)


def _common_points(first: DeclaredPrimitive, second: DeclaredPrimitive) -> tuple[OneRootPoint2, ...] | None:
    """Every common point of two primitives, or ``None`` for an overlap of positive length."""
    if isinstance(first, DeclaredLine) and isinstance(second, DeclaredLine):
        return _line_line_common(first, second)
    if isinstance(first, DeclaredLine) and isinstance(second, DeclaredArc):
        return _line_arc_common(first, second)
    if isinstance(first, DeclaredArc) and isinstance(second, DeclaredLine):
        return _line_arc_common(second, first)
    assert isinstance(first, DeclaredArc) and isinstance(second, DeclaredArc)
    return _arc_arc_common(first, second)


def _line_line_common(first: DeclaredLine, second: DeclaredLine) -> tuple[OneRootPoint2, ...] | None:
    u, v = first.direction.vector, second.direction.vector
    offset = second.start.minus(first.start)
    determinant = _cross(u, v)
    if determinant != 0:
        along_first, along_second = _cross(offset, v) / determinant, _cross(offset, u) / determinant
        if 0 <= along_first <= first.length and 0 <= along_second <= second.length:
            point = first.start.offset(u, along_first)
            return (OneRootPoint2.rational(point.x, point.y),)
        return ()
    if _cross(offset, u) != 0:
        return ()
    ends = (_dot(u, offset), _dot(u, second.end.minus(first.start)))
    low, high = max(Fraction(0), min(ends)), min(Fraction(first.length), max(ends))
    if low < high:
        return None
    if low == high:
        point = first.start.offset(u, low)
        return (OneRootPoint2.rational(point.x, point.y),)
    return ()


def _on_line(line: DeclaredLine, point: OneRootPoint2) -> bool:
    dx, dy = point.minus(line.start)
    along = line.direction.vector[0] * dx + line.direction.vector[1] * dy
    return along.sign() >= 0 and along.compare(line.length) <= 0


def on_declared_arc(arc: DeclaredArc, point: OneRootPoint2) -> bool:
    """Whether a point of the arc's circle lies on the arc (ends included), by one exact chord-side sign.

    An arc of sweep ≤ π is the part of its circle on its own side of the chord
    ``S → E``: ``σ·cross(E - S, p - S) <= 0``. A major arc is the circle minus
    the open minor arc ``E → S``: ``σ·cross(S - E, p - E) >= 0``.
    """
    chord = arc.end.minus(arc.start)
    if arc.sense * _cross(chord, arc.centre.minus(arc.start)) >= 0:
        dx, dy = point.minus(arc.start)
        return arc.sense * (chord[0] * dy - chord[1] * dx).sign() <= 0
    dx, dy = point.minus(arc.end)
    return arc.sense * (chord[1] * dx - chord[0] * dy).sign() >= 0


def _on_declared_arc_rational(arc: DeclaredArc, point: ExactPoint2[WorldXY]) -> bool:
    """:func:`on_declared_arc` for a rational point, in ``Fraction`` arithmetic."""
    chord = arc.end.minus(arc.start)
    if arc.sense * _cross(chord, arc.centre.minus(arc.start)) >= 0:
        return arc.sense * _cross(chord, point.minus(arc.start)) <= 0
    return arc.sense * _cross(point.minus(arc.end), chord) >= 0


def _line_arc_common(line: DeclaredLine, arc: DeclaredArc) -> tuple[OneRootPoint2, ...]:
    crossings = _line_circle(line.start, line.direction.vector, arc.centre, arc.radius_squared)
    return tuple(point for point in crossings if _on_line(line, point) and on_declared_arc(arc, point))


def _arc_arc_common(first: DeclaredArc, second: DeclaredArc) -> tuple[OneRootPoint2, ...] | None:
    if first.centre != second.centre:
        crossings = _circle_circle_points(first.centre, first.radius_squared, second.centre, second.radius_squared)
        return tuple(point for point in crossings if on_declared_arc(first, point) and on_declared_arc(second, point))
    if first.radius != second.radius:
        return ()
    centre = first.centre
    first_low, first_high = _counter_clockwise_ends(first)
    second_low, second_high = _counter_clockwise_ends(second)
    if compare_on_arc(centre, 1, first_low, second_low, first_high) < 0 or compare_on_arc(centre, 1, second_low, first_low, second_high) < 0:
        return None
    return (*(point for point in (first_high,) if point == second_low), *(point for point in (second_high,) if point == first_low))


def _counter_clockwise_ends(arc: DeclaredArc) -> tuple[OneRootPoint2, OneRootPoint2]:
    start, end = OneRootPoint2.from_exact(arc.start), OneRootPoint2.from_exact(arc.end)
    return (start, end) if arc.sense == 1 else (end, start)


def _counter_clockwise(primitives: tuple[DeclaredPrimitive, ...]) -> bool:
    """Turn direction at the lexicographically lowest point of a simple cycle (see the validator)."""
    lowest: tuple[Fraction, Fraction] | None = None
    counter_clockwise = True
    for index, primitive in enumerate(primitives):
        start = primitive.start
        if lowest is None or (start.y, start.x) < lowest:
            lowest = (start.y, start.x)
            incoming, outgoing = primitives[index - 1].end_tangent, primitive.start_tangent
            turn = _cross(incoming, outgoing)
            counter_clockwise = turn > 0 if turn != 0 else outgoing[0] > 0
        if not isinstance(primitive, DeclaredArc):
            continue
        bottom = primitive.centre.offset((Fraction(0), Fraction(-1)), primitive.radius)
        strictly_inside = bottom != primitive.start and bottom != primitive.end and _on_declared_arc_rational(primitive, bottom)
        assert lowest is not None
        if strictly_inside and (bottom.y, bottom.x) < lowest:
            lowest = (bottom.y, bottom.x)
            counter_clockwise = primitive.sense == 1
    return counter_clockwise
