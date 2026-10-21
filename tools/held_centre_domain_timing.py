"""Time the exact Held centre domain Ω*: construction and contact projection, one fresh process per sample.

For every canonical Held case, each of :data:`SAMPLES` samples runs in a fresh
interpreter (``python -m tools.held_centre_domain_timing --case NAME``) and
measures, wall clock and process CPU:

- **construction**: ``CentreDomain.build`` on the declared boundary (native
  arrangement plus support matching), and the native part alone;
- **projection**: the exact native decision of every distance-admissible guide
  hypothesis contact (``CentreDomain.decide_contacts``), one batch call.

The hypotheses are built before the clock starts. The result records the
median of each figure, the decision counts (projected, beyond the gap) and the
number of predicates the interval stage left to the exact stage. Every sample
must report identical counts: the decisions are deterministic.

Budgets (docs/held_exact_centre_domain.md, "Time budget"): Monstera
construction ≤ 2 s, Monstera projection of all hypotheses ≤ 1 s.

Run as ``pixi run held-centre-domain-timing``; writes
``benchmarks/results/held_centre_domain_timing.json``.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import _held_centre_domain_2
from compas_cgal.units import ExactMillimetre

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = ROOT / "benchmarks" / "results" / "held_centre_domain_timing.json"
# Fresh-process samples per case; the median is reported.
SAMPLES = 3
# The child prints its result on one stdout line after this prefix.
RESULT_LINE_PREFIX = "HELD_CENTRE_DOMAIN_TIMING "
# Stored seconds are rounded to milliseconds (reporting only).
SECONDS_DECIMALS = 3
# Counts every sample of one case must agree on.
DETERMINISTIC_KEYS = ("pieces", "members", "hypotheses", "projected", "beyond_gap", "exact_predicates")


class NondeterministicTimingSampleError(Exception):
    """Two samples of one case disagree on a decision count."""


class TimingProcessFailedError(Exception):
    """A child timing process failed."""


def measure_case(name: str) -> dict[str, Any]:
    """Measure one case in this process."""
    case = load_held_reference_case(name)
    declared = declared_boundary_from_reference(name, case.boundary)
    radius = Fraction(float(case.tool_radius.value))
    guide = build_held_reference_raw_guide(case)
    points = [
        (float(hypothesis.contact_point.x), float(hypothesis.contact_point.y))
        for run in guide.runs
        for station in run.stations
        for hypothesis in build_distance_admissible_hypotheses(guide, station)
    ]
    gap = ExactMillimetre(Fraction(float(guide.site_budget.admissible_distance_gap)))
    wall, cpu = time.perf_counter(), time.process_time()
    domain = CentreDomain.build(declared, radius)
    build_wall, build_cpu = time.perf_counter() - wall, time.process_time() - cpu
    wall, cpu = time.perf_counter(), time.process_time()
    decided = domain.decide_contacts(points, gap)
    project_wall, project_cpu = time.perf_counter() - wall, time.process_time() - cpu
    return {
        "case": name,
        "pieces": len(domain.pieces),
        "members": int(domain.native.member_count),
        "construction_wall_s": build_wall,
        "construction_cpu_s": build_cpu,
        "construction_native_s": float(domain.native.construction_seconds),
        "hypotheses": len(points),
        "projected": sum(piece is not None for piece in decided.contacts),
        "beyond_gap": sum(piece is None for piece in decided.contacts),
        "exact_predicates": int(_held_centre_domain_2.last_exact_fallbacks()),
        "projection_wall_s": project_wall,
        "projection_cpu_s": project_cpu,
    }


def _sample(name: str) -> dict[str, Any]:
    environment = dict(os.environ)
    completed = subprocess.run([sys.executable, "-m", "tools.held_centre_domain_timing", "--case", name], cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise TimingProcessFailedError(f"{name}: exit {completed.returncode}\n{completed.stderr[-2000:]}")
    lines = [line for line in completed.stdout.splitlines() if line.startswith(RESULT_LINE_PREFIX)]
    if len(lines) != 1:
        raise TimingProcessFailedError(f"{name}: expected one result line, got {len(lines)}")
    result: dict[str, Any] = json.loads(lines[0][len(RESULT_LINE_PREFIX) :])
    return result


def measure_all() -> dict[str, Any]:
    """Every case, :data:`SAMPLES` fresh-process samples each, serially."""
    cases = []
    for name in CANONICAL_CASE_NAMES:
        samples = [_sample(name) for _ in range(SAMPLES)]
        for key in DETERMINISTIC_KEYS:
            if len({sample[key] for sample in samples}) != 1:
                raise NondeterministicTimingSampleError(f"{name}: samples disagree on {key}: {[sample[key] for sample in samples]}")
        summary: dict[str, Any] = {key: samples[0][key] for key in ("case", *DETERMINISTIC_KEYS)}
        for key in ("construction_wall_s", "construction_cpu_s", "construction_native_s", "projection_wall_s", "projection_cpu_s"):
            summary[key] = round(statistics.median(sample[key] for sample in samples), SECONDS_DECIMALS)
        cases.append(summary)
    return {"samples": SAMPLES, "load_average_1m_at_end": round(os.getloadavg()[0], 1), "cases": cases}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CANONICAL_CASE_NAMES)
    arguments = parser.parse_args()
    if arguments.case is not None:
        print(RESULT_LINE_PREFIX + json.dumps(measure_case(arguments.case)))
        return
    result = measure_all()
    OUTPUT_PATH.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
