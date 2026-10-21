"""Wire layer between the declared boundary and the native Ω* arrangement (``_held_centre_domain_2``).

Rationals cross as ``(numerator, denominator)`` decimal text and one-root
numbers as ``(rational, coefficient, radicand)`` triples of those; no float
crosses in either direction. This module only encodes the declared boundary,
calls the native construction and decodes its boundary cycles into
:class:`NativePiece` records. Matching pieces to supports is the domain's job
(:mod:`benchmarks.held_exact_centre_domain`). Design:
docs/held_exact_centre_domain.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any
from typing import TypeAlias

from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal import _held_centre_domain_2
from compas_cgal.units import WorldXY

_RationalWire: TypeAlias = tuple[str, str]
_PrimitiveWire: TypeAlias = tuple[int, int, list[_RationalWire]]

# Primitive kinds of the native wire (held_centre_domain_2.h, PrimitiveKind).
LINE_KIND = 0
ARC_KIND = 1
# Nearest-point kinds (held_centre_domain_query_2.h, NearestKind).
VERTEX_KIND = 0
INTERIOR_KIND = 1
# Batch projection statuses (held_centre_domain_2_bindings.cpp).
PROJECTED = 0
AMBIGUOUS = 1
BEYOND_GAP = 2


@dataclass(frozen=True)
class NativeLine:
    """A line support ``a·x + b·y + c = 0`` as the arrangement stores it (not normalised)."""

    a: Fraction
    b: Fraction
    c: Fraction


@dataclass(frozen=True)
class NativeCircle:
    """A circle support: rational centre and squared radius."""

    centre_x: Fraction
    centre_y: Fraction
    radius_squared: Fraction


@dataclass(frozen=True)
class NativePiece:
    """One maximal curve of a boundary cycle on one support, travelled ``source → target``.

    ``sense`` is the travel sense on a circle (``1`` CCW, ``-1`` CW), ``0`` on a line.
    The face of Ω* lies to the left.
    """

    support: NativeLine | NativeCircle
    sense: int
    source: OneRootPoint2
    target: OneRootPoint2


@dataclass(frozen=True)
class NativeCentreDomain:
    """The decoded arrangement: per component, the outer cycle then its holes; plus the live native handle."""

    components: tuple[tuple[tuple[NativePiece, ...], ...], ...]
    member_count: int
    construction_seconds: float
    handle: Any


def _wire(value: Fraction) -> _RationalWire:
    value = Fraction(value)
    return str(value.numerator), str(value.denominator)


def _rational(wire: _RationalWire) -> Fraction:
    return Fraction(int(wire[0]), int(wire[1]))


def _one_root(wire: tuple[_RationalWire, _RationalWire, _RationalWire]) -> OneRootNumber:
    return OneRootNumber.build(_rational(wire[0]), _rational(wire[1]), _rational(wire[2]))


def _point(wire: Any) -> OneRootPoint2:
    return OneRootPoint2.build(_one_root(wire[0]), _one_root(wire[1]))


def _piece(wire: Any) -> NativePiece:
    circular, support, sense, source, target = wire
    values = tuple(_rational(value) for value in support)
    decoded: NativeLine | NativeCircle = NativeCircle(*values) if circular else NativeLine(*values)
    return NativePiece(decoded, int(sense), _point(source), _point(target))


def encode_boundary(boundary: DeclaredBoundary) -> list[_PrimitiveWire]:
    """The declared cycle on the native wire (lines carry their exact interior normal)."""
    primitives: list[_PrimitiveWire] = []
    for primitive in boundary.primitives:
        ends = [*_point_wire(primitive.start), *_point_wire(primitive.end)]
        if isinstance(primitive, DeclaredLine):
            normal = primitive.direction.left_normal
            primitives.append((LINE_KIND, 0, [*ends, _wire(normal[0]), _wire(normal[1])]))
            continue
        primitives.append((ARC_KIND, primitive.sense, [*ends, *_point_wire(primitive.centre), _wire(primitive.radius)]))
    return primitives


def _point_wire(point: ExactPoint2[WorldXY]) -> tuple[_RationalWire, _RationalWire]:
    return _wire(point.x), _wire(point.y)


def build_native_centre_domain(boundary: DeclaredBoundary, tool_radius: Fraction) -> NativeCentreDomain:
    """Construct Ω* natively and decode its boundary cycles.

    Raises:
        compas_cgal._held_centre_domain_2.InvalidCentreDomainInputError: Malformed declared input.
        compas_cgal._held_centre_domain_2.BrokenCentreDomainChainError: An extracted cycle does not chain.
    """
    handle = _held_centre_domain_2.CentreDomain2.build(encode_boundary(boundary), _wire(tool_radius))
    components = tuple(tuple(tuple(_piece(piece) for piece in cycle) for cycle in component) for component in handle.components())
    return NativeCentreDomain(components, int(handle.member_count), float(handle.construction_seconds), handle)


@dataclass(frozen=True, eq=False)
class NativeNearest:
    """Exact nearest ∂Ω* point of a query: on piece ``piece`` (global index), at its source for a vertex."""

    ambiguous: bool
    is_vertex: bool
    piece: int
    squared_distance: OneRootNumber
    point: OneRootPoint2


@dataclass(frozen=True)
class NativeContactDecision:
    """One batch projection decision: ``status`` (:data:`PROJECTED`, :data:`AMBIGUOUS`, :data:`BEYOND_GAP`).

    The nearest point is the source of piece ``piece`` when ``is_vertex``, else
    interior to it (its coordinates come from :func:`native_nearest`).
    """

    status: int
    is_vertex: bool
    piece: int


def native_nearest(handle: Any, point: ExactPoint2[WorldXY]) -> NativeNearest:
    """Exact nearest boundary point of a rational point (native, exact, box-pruned)."""
    status, kind, piece, squared, x, y = handle.nearest(_wire(point.x), _wire(point.y))
    location = OneRootPoint2.build(_one_root(x), _one_root(y))
    return NativeNearest(status == AMBIGUOUS, kind == VERTEX_KIND, int(piece), _one_root(squared), location)


@dataclass(frozen=True)
class NativeContactDecisions:
    """One batch call's decisions and the number of predicates its interval stage left to the exact stage."""

    decisions: tuple[NativeContactDecision, ...]
    exact_fallbacks: int


def native_decide_binary64(handle: Any, points: list[tuple[float, float]], gap: Fraction) -> NativeContactDecisions:
    """Decide the nearest ∂Ω* point of binary64 points (injected exactly) in one native call.

    The exact-fallback count is read from the extension's module-level
    ``last_exact_fallbacks()`` immediately after the call and returned with
    the decisions, so no caller reads the global. Stopgap: the binding should
    return the count with the rows (queued for the native slot).
    """
    rows = handle.project_binary64([point[0] for point in points], [point[1] for point in points], _wire(gap))
    exact_fallbacks = int(_held_centre_domain_2.last_exact_fallbacks())
    decisions = tuple(NativeContactDecision(int(status), kind == VERTEX_KIND, int(piece)) for status, kind, piece in rows)
    return NativeContactDecisions(decisions, exact_fallbacks)


def native_locate(handle: Any, point: ExactPoint2[WorldXY]) -> int:
    """Exact location of a rational point: ``1`` interior, ``0`` boundary, ``-1`` exterior."""
    return int(handle.locate(_wire(point.x), _wire(point.y)))
