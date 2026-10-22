"""Results of the Ω* connector repair stage for docs/held_connector_engagement.md.

Two modes:

- ``--measure`` (minutes): generates every Ω* workload
  (:func:`benchmarks.held_omega_toolpath.generate_omega`), whose last stage,
  ``connector_repair``, gates every connector exactly and repairs each refused
  one with ∂Ω* midpoint circles
  (:func:`benchmarks.held_omega_connector_repair.repair_omega_connectors`). It
  records, per workload: the refusals before (the stage's own gate on the
  NOTHING-stage circles) and after (a fresh exact assessment of the emitted
  path, which must refuse nothing and give the stage's tiers), the circles
  added per repaired connector, the stage seconds split into gate and
  insertion pass, and the exact gouge audit of the path before and after
  repair. The raw record goes to
  ``benchmarks/results/held_omega_connector_repair.json``, then the table is
  rendered.
- default: renders ``benchmarks/results/held_omega_connector_repair.md`` from
  that JSON without measuring.

Both files are written by this tool; do not edit them. Times are seconds on a
shared machine with the 1-minute load average beside each; only the
same-process ratio of the stage to the rest of generation is meaningful.
"""

from __future__ import annotations

import argparse
import collections
import gc
import json
import os
import pathlib
import time
from collections.abc import Sequence
from dataclasses import asdict
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_connector_engagement import assess_omega_connectors
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_toolpath import OMEGA_WORKLOADS
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.performance_record import REPOSITORY
from benchmarks.performance_record import build_identity
from benchmarks.performance_record import machine_identity

RESULTS = REPOSITORY / "benchmarks" / "results"
JSON_PATH = RESULTS / "held_omega_connector_repair.json"
MARKDOWN_PATH = RESULTS / "held_omega_connector_repair.md"


class RepairRecordContractError(BenchmarkError):
    """The emitted path still has a refused connector, or a fresh assessment disagrees with the repair stage's verdicts."""


def check_repair_contract(workload: str, *, refused_after: int, fresh_equal: bool) -> None:
    """Fail loudly unless the emitted path refuses nothing and the fresh tiers equal the stage's.

    Raises:
        RepairRecordContractError: Either condition is breached.
    """
    if refused_after:
        raise RepairRecordContractError(f"{workload}: {refused_after} connectors refused after repair")
    if not fresh_equal:
        raise RepairRecordContractError(f"{workload}: a fresh assessment of the emitted path disagrees with the repair verdicts")


def _gouge(report: HeldGougeAuditReport) -> list[int]:
    return [len(report.gouging_circles), len(report.gouging_connectors)]


def measure_workload(workload: str) -> dict[str, Any]:
    """Generate one workload; record its repair stage, a fresh post-gate and the gouge audits."""
    gc.collect()
    load_before = os.getloadavg()[0]
    wall, cpu = time.perf_counter(), time.process_time()
    path = generate_omega(workload)
    generation_wall, generation_cpu = time.perf_counter() - wall, time.process_time() - cpu
    repair, before = path.repair, path.nothing.circles
    fresh = assess_omega_connectors(path)
    fresh_equal = [v.tier for v in fresh.verdicts] == [v.tier for v in repair.verdicts]
    check_repair_contract(workload, refused_after=len(fresh.refused), fresh_equal=fresh_equal)
    unrepaired = OmegaConnectors.build(path.loop.pieces, before)
    stage = path.stage_seconds["connector_repair"]
    return {
        "workload": workload,
        "cap_deg": float(path.case.tea_cap),
        "circles_before": len(before),
        "circles_after": len(path.circles),
        "refused_before": len(repair.actions),
        "refused_after": len(fresh.refused),
        "repaired_connectors": [action.connector for action in repair.actions],
        "seam_repaired": any(action.seam for action in repair.actions),
        "inserted": repair.inserted,
        "inserted_per_action": [action.inserted for action in repair.actions],
        "first_refusal_kinds": dict(collections.Counter(action.refusals[0].kind.value for action in repair.actions)),
        "tiers_after": {tier.value: fresh.count(tier) for tier in ConnectorTier},
        "fresh_tiers_equal_repair_verdicts": fresh_equal,
        "gouge_before": _gouge(audit_omega_path(path.clearance, before, unrepaired.connectors)),
        "gouge_after": _gouge(path.audit()),
        "connector_repair_seconds": stage,
        "repair_gate_seconds": repair.gate_seconds,
        "repair_pass_seconds": repair.repair_seconds,
        "generation_seconds": path.generation_seconds,
        "stage_over_rest": stage / (path.generation_seconds - stage),
        "generation_wall_seconds": generation_wall,
        "generation_cpu_seconds": generation_cpu,
        "load_1min_before_after": [load_before, os.getloadavg()[0]],
    }


def measure(workloads: Sequence[str]) -> dict[str, Any]:
    """Every workload in one process, sequentially; return the record."""
    build, machine = build_identity(), machine_identity()
    rows = []
    for workload in workloads:
        rows.append(measure_workload(workload))
        print(f"repair {workload} done", flush=True)
    return {"build": asdict(build), "machine": asdict(machine), "workloads": rows}


def _kinds(counts: dict[str, int]) -> str:
    return ", ".join(f"{kind} {count}" for kind, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))) or "—"


def _per_action(inserted: list[int]) -> str:
    counts = collections.Counter(inserted)
    return ", ".join(f"{count}×{n}" for n, count in sorted(counts.items())) or "—"


def render(record: dict[str, Any]) -> str:
    """The Markdown table of one record."""
    build = record["build"]
    lines = [
        f"Recorded {build['recorded_at']} at commit `{build['commit'][:10]}`{' (dirty tree)' if build['dirty'] else ''} on {record['machine']['cpu_brand']}.",
        "",
        "**Ω\\* connector repair.** Refused before: the repair stage's exact gate on the NOTHING-stage circles. Refused after: a fresh",
        "exact assessment of the emitted path (its tiers equal the stage's verdicts). Inserted per repaired connector: count × circles.",
        "Gouge: gouging circles / connectors by the exact audit. Stage s = gate s + pass s (+ bookkeeping); stage / rest is the",
        "same-process ratio to the other generation stages.",
        "",
        "| workload | circles before → after | refused before → after | inserted | per repaired connector | first refused sweep | seam repaired "
        "| gouge before / after | stage s | gate s | pass s | gen s | stage / rest | load |",
        "|---|---|---|---:|---|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in record["workloads"]:
        load = row["load_1min_before_after"]
        gb, ga = row["gouge_before"], row["gouge_after"]
        lines.append(
            f"| `{row['workload']}` | {row['circles_before']} → {row['circles_after']} | {row['refused_before']} → {row['refused_after']} | {row['inserted']} "
            f"| {_per_action(row['inserted_per_action'])} | {_kinds(row['first_refusal_kinds'])} | {'yes' if row['seam_repaired'] else 'no'} "
            f"| {gb[0]}/{gb[1]} / {ga[0]}/{ga[1]} | {row['connector_repair_seconds']:.2f} | {row['repair_gate_seconds']:.2f} | {row['repair_pass_seconds']:.3f} "
            f"| {row['generation_seconds']:.2f} | {row['stage_over_rest']:.2f} | {load[0]:.0f}–{load[1]:.0f} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measure", action="store_true", help="re-measure (minutes) before rendering")
    parser.add_argument("--json", type=pathlib.Path, default=JSON_PATH)
    parser.add_argument("--markdown", type=pathlib.Path, default=MARKDOWN_PATH)
    parser.add_argument("--workloads", nargs="+", default=list(OMEGA_WORKLOADS), choices=list(OMEGA_WORKLOADS))
    args = parser.parse_args()
    if args.measure:
        record = measure(args.workloads)
        args.json.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    record = json.loads(args.json.read_text(encoding="utf-8"))
    args.markdown.write_text(render(record), encoding="utf-8")
    print(f"WROTE {args.json} {args.markdown}")


if __name__ == "__main__":
    main()
