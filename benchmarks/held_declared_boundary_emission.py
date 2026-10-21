"""Re-emit stored Held reference primitives as a declared (Pythagorean) boundary.

Float ``atan2`` / ``tan`` / ``hypot`` only CHOOSE the binary64 / dyadic
parameters (the half-angle tangent ``t`` and the line length ``L``); no exact
value is ever derived through a float square root.

Junction rule (docs/held_exact_centre_domain.md): an arc ends exactly on its
published binary64 end; every primitive starts at its predecessor's exact end;
a line ends at ``start + L·u``. The walk begins right after an arc, so the
last primitive is an arc that ends on the walk's first published junction and
the cycle closes exactly. Stored starts other than the walk's first are not
read: a stored gap between a primitive's end and its successor's start is
repaired, and the certificate's move bound accounts for it.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from fractions import Fraction

from benchmarks.held_declared_boundary import QUARTER_TURNS
from benchmarks.held_declared_boundary import T_GRID_DENOMINATOR
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredBoundaryError
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import DeclaredPrimitive
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_declared_boundary import Sense
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

QUARTER_TURN_RADIANS = math.pi / 2


class UnanchoredDeclaredCycleError(DeclaredBoundaryError):
    """The stored cycle has no arc, so no published junction anchors exact closure."""


class DegenerateStoredDirectionError(DeclaredBoundaryError):
    """A stored line or arc start tangent has zero binary64 length."""


def choose_direction(x: float, y: float) -> PythagoreanDirection:
    """Choose the grid direction nearest the float direction ``(x, y)``.

    Args:
        x: Binary64 x component (any non-zero length).
        y: Binary64 y component.

    Returns:
        The Pythagorean direction with ``t`` rounded to the ``2^-56`` grid.

    Raises:
        DegenerateStoredDirectionError: If ``(x, y)`` is the zero vector.
    """
    if x == 0.0 and y == 0.0:
        raise DegenerateStoredDirectionError("Cannot choose a direction for the zero vector.")
    angle = math.atan2(y, x)
    turns = round(angle / QUARTER_TURN_RADIANS)
    residual = angle - turns * QUARTER_TURN_RADIANS
    t = Fraction(round(math.tan(residual / 2.0) * T_GRID_DENOMINATOR), T_GRID_DENOMINATOR)
    return PythagoreanDirection.build(t, turns % QUARTER_TURNS)


def stored_primitives_sha256(primitives: Sequence[ReferencePrimitive]) -> str:
    """Hash the stored primitives by the exact hex of every binary64 they carry.

    Args:
        primitives: Stored reference primitives.

    Returns:
        SHA-256 hex digest of their canonical serialisation.
    """
    records: list[list[str]] = []
    for primitive in primitives:
        points = (primitive.start, primitive.end) if isinstance(primitive, ReferenceLine) else (primitive.start, primitive.end, primitive.centre)
        record = ["line" if isinstance(primitive, ReferenceLine) else "arc"]
        record += [float(value).hex() for point in points for value in (point.x, point.y)]
        if isinstance(primitive, ReferenceArc):
            record.append(float(primitive.sweep).hex())
        records.append(record)
    return hashlib.sha256(json.dumps(records, separators=(",", ":")).encode()).hexdigest()


def declared_boundary_from_reference(case_name: str, boundary: ReferenceBoundary) -> DeclaredBoundary:
    """Re-emit a validated stored boundary.

    Args:
        case_name: Name of the reference case.
        boundary: Stored boundary (closed under its float convention).

    Returns:
        The exactly closed declared boundary.
    """
    return reemit_primitives(case_name, boundary.primitives)


def reemit_primitives(case_name: str, primitives: Sequence[ReferencePrimitive]) -> DeclaredBoundary:
    """Re-emit stored primitives as an exactly closed declared boundary.

    The stored sequence need not satisfy the stored closure check; a stored
    gap is repaired here and bounded by the certificate's move bound.

    Args:
        case_name: Name of the reference case.
        primitives: Stored primitives in cycle order.

    Returns:
        The declared boundary, primitive ``i`` re-emitting stored primitive ``i``.

    Raises:
        UnanchoredDeclaredCycleError: If no primitive is an arc.
    """
    count = len(primitives)
    anchor = next((index for index in range(count) if isinstance(primitives[index - 1], ReferenceArc)), None)
    if anchor is None:
        raise UnanchoredDeclaredCycleError(f"{case_name}: a cycle without arcs has no published junction to close on.")
    emitted: dict[int, DeclaredPrimitive] = {}
    current = _exact(primitives[anchor].start)
    for offset in range(count):
        index = (anchor + offset) % count
        primitive = primitives[index]
        declared: DeclaredPrimitive
        if isinstance(primitive, ReferenceLine):
            declared = _reemit_line(primitive, current, index)
        else:
            declared = _reemit_arc(primitive, current, index)
        emitted[index] = declared
        current = declared.end
    return DeclaredBoundary.build(case_name, tuple(emitted[index] for index in range(count)), stored_primitives_sha256(primitives))


def _reemit_line(line: ReferenceLine, start: ExactPoint2[WorldXY], index: int) -> DeclaredLine:
    dx = float(line.end.x) - float(line.start.x)
    dy = float(line.end.y) - float(line.start.y)
    length = math.hypot(dx, dy)
    if length == 0.0:
        raise DegenerateStoredDirectionError(f"Stored line {index} has zero binary64 length.")
    return DeclaredLine.build(start, choose_direction(dx / length, dy / length), Fraction(length), index)


def _reemit_arc(arc: ReferenceArc, start: ExactPoint2[WorldXY], index: int) -> DeclaredArc:
    sense: Sense = 1 if float(arc.sweep) > 0.0 else -1
    radial_x = float(arc.start.x) - float(arc.centre.x)
    radial_y = float(arc.start.y) - float(arc.centre.y)
    tangent_x, tangent_y = -sense * radial_y, sense * radial_x
    length = math.hypot(tangent_x, tangent_y)
    if length == 0.0:
        raise DegenerateStoredDirectionError(f"Stored arc {index} has a zero start radius.")
    return DeclaredArc.build(start, _exact(arc.end), choose_direction(tangent_x / length, tangent_y / length), sense, index)


def _exact(point: Point2[WorldXY]) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].from_binary64(float(point.x), float(point.y))
