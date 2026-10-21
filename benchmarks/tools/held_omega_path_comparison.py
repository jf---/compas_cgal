"""Old (polygon Ω) versus new (Ω*) Held paths: exact gouge audits, counts and generation time.

For each workload both generators run in this process, one after the other.
Both paths are audited by one exact oracle against the declared boundary
(:mod:`benchmarks.held_omega_gouge_audit`); the polygon path is also audited
by the stored-boundary audit it was published under
(:mod:`benchmarks.held_gouge_audit`). Writes one JSON record per workload.

Usage: ``python -m benchmarks.tools.held_omega_path_comparison --output build/omega-paths/comparison.json [--workloads ...]``.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import pathlib
import time
from fractions import Fraction
from typing import Any

from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_gouge_workloads import audit_generated
from benchmarks.held_omega_gouge_audit import audit_polygon_path
from benchmarks.held_omega_toolpath import OMEGA_WORKLOADS
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.held_omega_toolpath import omega_metrics
from benchmarks.held_toolpath_benchmark import _metrics
from benchmarks.held_toolpath_benchmark import generate


def _audit_summary(report: HeldGougeAuditReport) -> dict[str, Any]:
    return {
        "gouging_circles": len(report.gouging_circles),
        "gouging_connectors": len(report.gouging_connectors),
        "touching_circles": len(report.touching_circles),
        "touching_connectors": len(report.touching_connectors),
        "max_circle_depth_mm": report.max_depth_mm("circle"),
        "max_connector_depth_mm": report.max_depth_mm("connector"),
        "gouge_primitive_kinds": sorted({item.primitive_kind for item in report.gouges}),
        "audit_seconds": report.seconds,
    }


def compare(workload: str) -> dict[str, Any]:
    """Generate both paths of one workload and audit them."""
    old = generate(workload)
    declared = declared_boundary_from_reference(old.case.name, old.case.boundary)
    clearance = DeclaredClearance.build(declared, Fraction(float(old.case.tool_radius.value)))
    old_declared = audit_polygon_path(clearance, [(circle.center, float(circle.radius.value)) for circle in old.circles], [transition.samples for transition in old.transitions])
    old_stored = audit_generated(old)
    started = time.perf_counter()
    new = generate_omega(workload)
    wall = time.perf_counter() - started
    new_audit = new.audit()
    return {
        "workload": workload,
        "old": {
            "metrics": dataclasses.asdict(_metrics(old)),
            "stage_seconds": old.stage_seconds,
            "generation_seconds": old.generation_seconds,
            "audit_declared": _audit_summary(old_declared),
            "audit_stored": _audit_summary(old_stored),
        },
        "new": {
            "metrics": dataclasses.asdict(omega_metrics(new)),
            "stage_seconds": new.stage_seconds,
            "generation_seconds": new.generation_seconds,
            "generation_wall_seconds": wall,
            "audit_declared": _audit_summary(new_audit),
            "omega_pieces": len(new.loop.pieces),
            "canonical_candidates": new.initial.canonical_candidate_count,
            "placement_lanes": new.initial.placement_lane_count,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--workloads", nargs="+", default=list(OMEGA_WORKLOADS), choices=list(OMEGA_WORKLOADS))
    args = parser.parse_args()
    records = []
    for workload in args.workloads:
        record = compare(workload)
        records.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"RECORD {args.output}")


if __name__ == "__main__":
    main()
