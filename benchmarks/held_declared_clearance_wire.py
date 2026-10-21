"""Wire encoding for the native declared-boundary clearance oracle (``_held_declared_clearance_2``).

Rationals cross as ``(numerator, denominator)`` decimal text, one-root numbers
as ``(rational, coefficient, radicand)`` triples of those, points as pairs of
one-root numbers (the ``_held_centre_domain_2`` convention). Binary64 values
cross as floats and are injected exactly on the C++ side; this module refuses
any value that is not exactly a binary64 rather than rounding it.

Pure Python: it never imports the extension, so it is testable on its own.
"""

from __future__ import annotations

from collections.abc import Iterable
from collections.abc import Sequence
from fractions import Fraction
from typing import TypeAlias

from benchmarks.errors import BenchmarkError
from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import WorldXY

RationalWire: TypeAlias = tuple[str, str]
RootWire: TypeAlias = tuple[RationalWire, RationalWire, RationalWire]
PointWire: TypeAlias = tuple[RootWire, RootWire]
PrimitiveWire: TypeAlias = tuple[int, int, list[RationalWire]]
ElementWire: TypeAlias = tuple[int, int, list[int], PointWire, PointWire, list[RationalWire]]
Binary64Circle: TypeAlias = tuple[float, float, float]

# Kinds of the native wire (src/held_declared_clearance_2/wire_2.h).
LINE_KIND = 0
ARC_KIND = 1


class ClearanceWireError(BenchmarkError):
    """Base error of the native clearance wire."""


class NonBinary64ClearanceInputError(ClearanceWireError):
    """A value sent on a binary64 lane is not exactly a binary64."""


def rational_wire(value: Fraction | int) -> RationalWire:
    """``(numerator, denominator)`` decimal text of an exact rational."""
    exact = Fraction(value)
    return str(exact.numerator), str(exact.denominator)


def root_wire(value: OneRootNumber) -> RootWire:
    """``(a, b, c)`` of ``a + b·√c``, as the Python oracle folds it."""
    return rational_wire(value.rational), rational_wire(value.coefficient), rational_wire(value.radicand)


def point_wire(point: OneRootPoint2) -> PointWire:
    """Both coordinates of a one-root point."""
    return root_wire(point.x), root_wire(point.y)


def primitive_wires(boundary: DeclaredBoundary) -> list[PrimitiveWire]:
    """Every declared primitive, in index order.

    A line is ``(0, 0, [sx, sy, ex, ey, ux, uy, length])``; an arc is
    ``(1, sense, [sx, sy, ex, ey, cx, cy, radius])``.
    """
    wires: list[PrimitiveWire] = []
    for primitive in boundary.primitives:
        ends = [rational_wire(primitive.start.x), rational_wire(primitive.start.y), rational_wire(primitive.end.x), rational_wire(primitive.end.y)]
        if isinstance(primitive, DeclaredLine):
            ux, uy = primitive.direction.vector
            wires.append((LINE_KIND, 0, [*ends, rational_wire(ux), rational_wire(uy), rational_wire(primitive.length)]))
        else:
            arc: DeclaredArc = primitive
            wires.append((ARC_KIND, arc.sense, [*ends, rational_wire(arc.centre.x), rational_wire(arc.centre.y), rational_wire(arc.radius)]))
    return wires


def element_wire(element: SegmentElement | ArcElement, owners: Sequence[int]) -> ElementWire:
    """One exact element and the declared primitives it may offset (each checked natively)."""
    start, end = point_wire(element.start), point_wire(element.end)
    if isinstance(element, SegmentElement):
        line = [
            rational_wire(element.direction[0]),
            rational_wire(element.direction[1]),
            rational_wire(element.normal[0]),
            rational_wire(element.normal[1]),
            rational_wire(element.level),
        ]
        return (LINE_KIND, 0, list(owners), start, end, line)
    circle = [rational_wire(element.centre.x), rational_wire(element.centre.y), rational_wire(element.radius)]
    return (ARC_KIND, element.sense, list(owners), start, end, circle)


def offset_owners(owners: Iterable[BoundaryFeature]) -> tuple[int, ...]:
    """The declared LINE and ARC owners of a ∂Ω* piece, ascending (junction owners have no lemma)."""
    return tuple(sorted(feature.index for feature in owners if feature.kind in (BoundaryFeatureKind.LINE, BoundaryFeatureKind.ARC)))


def binary64(value: Fraction) -> float:
    """The binary64 equal to ``value``.

    Raises:
        NonBinary64ClearanceInputError: If ``value`` is not exactly a binary64.
    """
    converted = float(value)
    if Fraction(converted) != value:
        raise NonBinary64ClearanceInputError(f"{value} is not exactly a binary64.")
    return converted


def binary64_circle(element: CircleElement) -> Binary64Circle:
    """``(x, y, radius)`` of a machining circle whose values are binary64.

    Raises:
        NonBinary64ClearanceInputError: If a value is not exactly a binary64.
    """
    return binary64(element.centre.x), binary64(element.centre.y), binary64(element.radius)


def binary64_point(point: ExactPoint2[WorldXY]) -> tuple[float, float]:
    """``(x, y)`` of a point whose coordinates are binary64.

    Raises:
        NonBinary64ClearanceInputError: If a coordinate is not exactly a binary64.
    """
    return binary64(point.x), binary64(point.y)
