"""Native versus Python declared-boundary gouge audits on the five Held workloads: decisions, depths and time.

For each workload the polygon-Ω path and the Ω* path are generated once, then
each is audited by the Python oracle and by the native oracle (native build
time included in its audit time). Decisions are compared element for element:
(element, index, clearance, primitive, primitive kind, segment) must be
identical; reported depths are compared bitwise and their mismatches counted
(reporting only). The 1-minute load average is recorded around every timing.

The ``generation`` section generates each Ω* path with both containment
oracles in one process, alternating which goes first, ``--generation-repeats``
times each: the paths must be byte-identical (SHA-256 of the ``repr`` of the
initial path, the refined circles and the connectors; a ``repr`` carrying a
memory address is refused) and the record holds the same-process time ratio.

Usage: ``python -m benchmarks.tools.held_declared_clearance_differential --output build/native-clearance/differential.json [--workloads ...]``;
add ``--sections audit generation [--generation-repeats N]`` for the generation A/B.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import re
import statistics
import time
from fractions import Fraction
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.held_clearance_backend import ClearanceBackend
from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_declared_clearance_native import NativeDeclaredClearance
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_gouge_audit import audit_polygon_path
from benchmarks.held_omega_gouge_audit_native import audit_omega_path_native
from benchmarks.held_omega_gouge_audit_native import audit_polygon_path_native
from benchmarks.held_omega_toolpath import OMEGA_WORKLOADS
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.held_toolpath_benchmark import generate

# ``object.__repr__`` and friends: "<module.Class object at 0x10a2b3c4d>" (not canonical).
_MEMORY_ADDRESS = re.compile(r" at 0x[0-9a-fA-F]+>")


class NonCanonicalReprError(BenchmarkError):
    """A digested ``repr`` carries a memory address, so equal paths could digest differently."""


def _decisions(report: HeldGougeAuditReport) -> list[tuple[object, ...]]:
    return [(item.element, item.index, int(item.clearance), item.primitive, item.primitive_kind, item.segment) for item in report.elements]


def _summary(report: HeldGougeAuditReport, seconds: float, load: tuple[float, float]) -> dict[str, Any]:
    return {
        "gouging_circles": len(report.gouging_circles),
        "gouging_connectors": len(report.gouging_connectors),
        "touching_circles": len(report.touching_circles),
        "touching_connectors": len(report.touching_connectors),
        "seconds": seconds,
        "load_1min_before_after": load,
    }


def _timed(run: Any) -> tuple[HeldGougeAuditReport, float, tuple[float, float]]:
    before = os.getloadavg()[0]
    started = time.perf_counter()
    report = run()
    seconds = time.perf_counter() - started
    return report, seconds, (before, os.getloadavg()[0])


def _compare(python: HeldGougeAuditReport, native: HeldGougeAuditReport) -> dict[str, Any]:
    expected, found = _decisions(python), _decisions(native)
    depth_mismatches = sum(1 for left, right in zip(python.elements, native.elements, strict=False) if left.depth_mm != right.depth_mm)
    mismatched = [(left, right) for left, right in zip(expected, found, strict=False) if left != right]
    return {
        "identical_decisions": expected == found,
        "listed_python": len(expected),
        "listed_native": len(found),
        "first_mismatches": [list(map(str, pair)) for pair in mismatched[:5]],
        "depth_bit_mismatches": depth_mismatches,
    }


def differential(workload: str) -> dict[str, Any]:
    """Generate both paths of one workload and audit each with both oracles."""
    old = generate(workload)
    declared = declared_boundary_from_reference(old.case.name, old.case.boundary)
    radius = Fraction(float(old.case.tool_radius.value))
    python_oracle = DeclaredClearance.build(declared, radius)
    circles = [(circle.center, float(circle.radius.value)) for circle in old.circles]
    polylines = [transition.samples for transition in old.transitions]
    old_python, old_python_s, old_python_load = _timed(lambda: audit_polygon_path(python_oracle, circles, polylines, ClearanceBackend.PYTHON))
    old_native, old_native_s, old_native_load = _timed(lambda: audit_polygon_path_native(NativeDeclaredClearance.build(declared, radius), circles, polylines))
    new = generate_omega(workload)
    connectors = new.connectors.connectors
    new_python, new_python_s, new_python_load = _timed(lambda: audit_omega_path(new.clearance, new.circles, connectors, ClearanceBackend.PYTHON))
    new_native, new_native_s, new_native_load = _timed(
        lambda: audit_omega_path_native(NativeDeclaredClearance.build(new.declared, new.clearance.tool_radius), new.circles, connectors)
    )
    return {
        "workload": workload,
        "polygon": {
            "python": _summary(old_python, old_python_s, old_python_load),
            "native": _summary(old_native, old_native_s, old_native_load),
            "comparison": _compare(old_python, old_native),
        },
        "omega": {
            "python": _summary(new_python, new_python_s, new_python_load),
            "native": _summary(new_native, new_native_s, new_native_load),
            "comparison": _compare(new_python, new_native),
        },
    }


def _digest(value: object) -> str:
    text = repr(value)
    if _MEMORY_ADDRESS.search(text):
        raise NonCanonicalReprError(f"repr of {type(value).__name__} carries a memory address")
    return hashlib.sha256(text.encode()).hexdigest()


def path_digests(path: OmegaGeneratedToolpath) -> dict[str, str]:
    """SHA-256 of the ``repr`` of the initial path, the refined circles and the connectors."""
    return {"initial": _digest(path.initial), "circles": _digest(path.circles), "connectors": _digest(path.connectors)}


def generation_ab(workload: str, repeats: int) -> dict[str, Any]:
    """Generate one Ω* path ``repeats`` times per containment oracle in this process, alternating the order."""
    runs: list[dict[str, Any]] = []
    for repeat in range(repeats):
        order = (ClearanceBackend.PYTHON, ClearanceBackend.NATIVE) if repeat % 2 == 0 else (ClearanceBackend.NATIVE, ClearanceBackend.PYTHON)
        for backend in order:
            before = os.getloadavg()[0]
            started = time.perf_counter()
            path = generate_omega(workload, backend)
            seconds = time.perf_counter() - started
            load = (before, os.getloadavg()[0])
            runs.append(
                {
                    "backend": backend.value,
                    "seconds": seconds,
                    "stage_seconds": path.stage_seconds,
                    "load_1min_before_after": load,
                    "circles": len(path.circles),
                    "digests": path_digests(path),
                }
            )
    medians = {backend.value: statistics.median(run["seconds"] for run in runs if run["backend"] == backend.value) for backend in ClearanceBackend}
    digests = {json.dumps(run["digests"], sort_keys=True) for run in runs}
    return {
        "repeats": repeats,
        "byte_identical": len(digests) == 1,
        "median_seconds": medians,
        "python_over_native": medians["python"] / medians["native"],
        "runs": runs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--workloads", nargs="+", default=list(OMEGA_WORKLOADS), choices=list(OMEGA_WORKLOADS))
    parser.add_argument("--sections", nargs="+", default=["audit"], choices=["audit", "generation"])
    parser.add_argument("--generation-repeats", type=int, default=3)
    args = parser.parse_args()
    records = []
    for workload in args.workloads:
        record: dict[str, Any] = differential(workload) if "audit" in args.sections else {"workload": workload}
        if "generation" in args.sections:
            record["generation"] = generation_ab(workload, args.generation_repeats)
        records.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"RECORD {args.output}")


if __name__ == "__main__":
    main()
