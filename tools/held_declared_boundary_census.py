"""Re-emit, re-certify, measure and time the declared Held boundaries.

Promoted from ``build/held-exact-omega/pythagorean_reemission_tgrid56.py``,
now driven by :mod:`benchmarks.held_declared_boundary_emission` and
:mod:`benchmarks.held_declared_boundary_certificate`. Exact quantities are
written as ``numerator/denominator`` text beside a float for reading.
Timing excludes the case load (which re-runs the publisher reconstruction);
it is the minimum over ``--repeats`` runs of re-emission and certification.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path

from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary_certificate import DeclaredBoundaryCertificate
from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc

DEFAULT_OUTPUT = Path("benchmarks/results/held_declared_boundary.json")
DEFAULT_REPEATS = 5


def _bits(value: Fraction) -> int:
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _number_sizes(declared: DeclaredBoundary, tool_radius: Fraction) -> dict[str, int]:
    sizes = {"junction": 0, "direction": 0, "arc_centre": 0, "arc_radius": 0, "offset_radius_squared": 0}
    for primitive in declared.primitives:
        sizes["junction"] = max(sizes["junction"], _bits(primitive.start.x), _bits(primitive.start.y))
        sizes["direction"] = max(sizes["direction"], *(_bits(component) for component in primitive.start_tangent))
        if isinstance(primitive, DeclaredArc):
            offset_radius = primitive.radius - tool_radius if primitive.sense > 0 else primitive.radius + tool_radius
            sizes["arc_centre"] = max(sizes["arc_centre"], _bits(primitive.centre.x), _bits(primitive.centre.y))
            sizes["arc_radius"] = max(sizes["arc_radius"], _bits(primitive.radius))
            sizes["offset_radius_squared"] = max(sizes["offset_radius_squared"], _bits(offset_radius * offset_radius))
    return sizes


def _float_displacements(case: HeldReferenceCase, declared: DeclaredBoundary) -> dict[str, float]:
    """The prototype's reporting proxies (centre and radius moves are not boundary moves)."""
    line_end = arc_centre = arc_radius = 0.0
    for primitive in declared.primitives:
        stored = case.boundary.primitives[primitive.source_index]
        if isinstance(primitive, DeclaredLine) or not isinstance(stored, ReferenceArc):
            line_end = max(line_end, math.dist(primitive.end.as_floats(), (float(stored.end.x), float(stored.end.y))))
            continue
        stored_centre = (float(stored.centre.x), float(stored.centre.y))
        arc_centre = max(arc_centre, math.dist(primitive.centre.as_floats(), stored_centre))
        arc_radius = max(arc_radius, abs(float(primitive.radius) - math.dist((float(stored.start.x), float(stored.start.y)), stored_centre)))
    return {"line_end_move_mm": line_end, "arc_centre_move_mm": arc_centre, "arc_radius_change_mm": arc_radius}


def measure(case: HeldReferenceCase, repeats: int) -> dict[str, object]:
    """Re-emit and certify one case, timing both.

    Args:
        case: Loaded reference case.
        repeats: Number of timed runs (the minimum is reported).

    Returns:
        A JSON-ready record.
    """
    emission_s = certification_s = math.inf
    for _ in range(repeats):
        started = time.perf_counter()
        declared = declared_boundary_from_reference(case.name, case.boundary)
        emitted = time.perf_counter()
        certificate = DeclaredBoundaryCertificate.build(case, declared)
        certified = time.perf_counter()
        emission_s = min(emission_s, emitted - started)
        certification_s = min(certification_s, certified - emitted)
    worst = max(certificate.moves, key=lambda move: move.bound)
    census = declared.junction_census()
    return {
        "case": case.name,
        "primitives": len(declared.primitives),
        "closes_exactly": True,
        "junctions": {"exact_g1": census.exact_g1, "left": census.left_turns, "right": census.right_turns, "reversals": census.reversals},
        "max_move_upper_mm": float(certificate.maximum_move),
        "max_move_upper_exact": _text(certificate.maximum_move),
        "max_move_primitive": worst.source_index,
        "max_move_kind": worst.kind,
        "stored_bound_mm": float(certificate.stored_bound),
        "declared_bound_mm": float(certificate.declared_bound),
        "declared_bound_exact": _text(certificate.declared_bound),
        "reconstruction_limit_mm": float(certificate.reconstruction_limit),
        "max_bits": _number_sizes(declared, Fraction(float(case.tool_radius.value))),
        "float_displacement_proxies": _float_displacements(case, declared),
        "declared_sha256": declared.identity.declared_sha256,
        "source_sha256": declared.identity.source_sha256,
        "certificate_sha256": certificate.identity_sha256,
        "emission_s": emission_s,
        "certification_s": certification_s,
    }


def main() -> None:
    """Measure all four cases, write JSON and print one line per case."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    arguments = parser.parse_args()
    records = [measure(load_held_reference_case(name), arguments.repeats) for name in CANONICAL_CASE_NAMES]
    arguments.output.write_text(json.dumps(records, indent=1) + "\n")
    for record in records:
        print(json.dumps(record))


if __name__ == "__main__":
    main()
