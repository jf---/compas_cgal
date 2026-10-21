"""Exact number census of the stored Held reference boundaries.

Promoted from ``build/held-exact-omega/boundary_number_census.py`` (the design
input of docs/held_exact_centre_domain.md, "What the published data gives").
Every binary64 is injected exactly as a Fraction; every count below is an
exact decision. Fields suffixed ``_float`` are reporting only.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferencePrimitive
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

DEFAULT_OUTPUT = Path("benchmarks/results/held_boundary_census.json")
# Reporting bins of junction turn magnitude, radians (upper edges; the last bin is open).
TURN_BIN_EDGES_RAD = (1e-12, 1e-9, 1e-6, 1e-3, math.radians(1.0))

_Vector = tuple[Fraction, Fraction]


@dataclass(frozen=True)
class StoredBoundaryCensus:
    """Exact census of one stored boundary."""

    case: str
    primitives: int
    lines: int
    arcs: int
    convex_arcs: int
    concave_arcs: int
    lines_with_rational_length: int
    arcs_with_rational_start_radius: int
    arcs_with_end_exactly_on_start_circle: int
    arc_end_radius_gap_max_float_mm: float
    arcs_with_radius_below_tool_radius: int
    junctions_exact_g1: int
    junctions_left_turn: int
    junctions_right_turn: int
    junctions_antiparallel: int
    junction_turn_bins_float: tuple[int, ...]
    certified_reconstruction_deviation_mm: float


def census(case: HeldReferenceCase) -> StoredBoundaryCensus:
    """Census of one case's stored primitives.

    Args:
        case: Loaded reference case.

    Returns:
        The census.
    """
    primitives = case.boundary.primitives
    arcs = [primitive for primitive in primitives if isinstance(primitive, ReferenceArc)]
    lines = [primitive for primitive in primitives if not isinstance(primitive, ReferenceArc)]
    tool_radius_squared = Fraction(float(case.tool_radius.value)) ** 2
    rational_lengths = sum(_is_rational_square(_squared(_minus(_exact(line.end), _exact(line.start)))) for line in lines)
    rational_radii = on_circle = below_tool = 0
    end_gap = 0.0
    for arc in arcs:
        centre = _exact(arc.centre)
        start_squared = _squared(_minus(_exact(arc.start), centre))
        end_squared = _squared(_minus(_exact(arc.end), centre))
        rational_radii += _is_rational_square(start_squared)
        on_circle += start_squared == end_squared
        below_tool += start_squared < tool_radius_squared
        end_gap = max(end_gap, abs(math.sqrt(float(end_squared)) - math.sqrt(float(start_squared))))
    g1 = left = right = antiparallel = 0
    bins = [0] * (len(TURN_BIN_EDGES_RAD) + 1)
    for incoming, outgoing in zip((primitives[-1], *primitives[:-1]), primitives):
        tangent_in = _tangents(incoming)[1]
        tangent_out = _tangents(outgoing)[0]
        turn = _cross(tangent_in, tangent_out)
        if turn > 0:
            left += 1
        elif turn < 0:
            right += 1
        elif _dot(tangent_in, tangent_out) > 0:
            g1 += 1
        else:
            antiparallel += 1
        magnitude = abs(math.atan2(float(turn), float(_dot(tangent_in, tangent_out))))
        bins[next((index for index, edge in enumerate(TURN_BIN_EDGES_RAD) if magnitude <= edge), len(TURN_BIN_EDGES_RAD))] += 1
    return StoredBoundaryCensus(
        case=case.name,
        primitives=len(primitives),
        lines=len(lines),
        arcs=len(arcs),
        convex_arcs=sum(float(arc.sweep) > 0.0 for arc in arcs),
        concave_arcs=sum(float(arc.sweep) < 0.0 for arc in arcs),
        lines_with_rational_length=rational_lengths,
        arcs_with_rational_start_radius=rational_radii,
        arcs_with_end_exactly_on_start_circle=on_circle,
        arc_end_radius_gap_max_float_mm=end_gap,
        arcs_with_radius_below_tool_radius=below_tool,
        junctions_exact_g1=g1,
        junctions_left_turn=left,
        junctions_right_turn=right,
        junctions_antiparallel=antiparallel,
        junction_turn_bins_float=tuple(bins),
        certified_reconstruction_deviation_mm=float(case.reconstruction.deviation_upper_bound),
    )


def _tangents(primitive: ReferencePrimitive) -> tuple[_Vector, _Vector]:
    start, end = _exact(primitive.start), _exact(primitive.end)
    if not isinstance(primitive, ReferenceArc):
        chord = _minus(end, start)
        return chord, chord
    centre = _exact(primitive.centre)
    sense = 1 if float(primitive.sweep) > 0.0 else -1
    start_radial, end_radial = _minus(start, centre), _minus(end, centre)
    return (-sense * start_radial[1], sense * start_radial[0]), (-sense * end_radial[1], sense * end_radial[0])


def _is_rational_square(value: Fraction) -> bool:
    return value >= 0 and math.isqrt(value.numerator) ** 2 == value.numerator and math.isqrt(value.denominator) ** 2 == value.denominator


def _exact(point: Point2[WorldXY]) -> _Vector:
    return Fraction(float(point.x)), Fraction(float(point.y))


def _minus(left: _Vector, right: _Vector) -> _Vector:
    return left[0] - right[0], left[1] - right[1]


def _squared(vector: _Vector) -> Fraction:
    return vector[0] * vector[0] + vector[1] * vector[1]


def _dot(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _cross(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def main() -> None:
    """Write the census of all four cases as JSON and print one line per case."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    records = [asdict(census(load_held_reference_case(name))) for name in CANONICAL_CASE_NAMES]
    arguments.output.write_text(json.dumps(records, indent=1) + "\n")
    for record in records:
        print(json.dumps(record))


if __name__ == "__main__":
    main()
