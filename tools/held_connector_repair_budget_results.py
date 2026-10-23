"""Cost of connector repair on the Held polygon workloads, for docs/held_connector_engagement.md.

Two measurements, one record:

- **Circle budget.** Every polygon workload is generated twice in one
  process: as published (``benchmarks.held_toolpath_benchmark.generate``:
  refinement with the connector gate, NOTHING-circle removal, connectors,
  post-gate) and as an ablation with the connector gate off (the refinement
  accepts a transition on its circle pair alone, as before the gate was
  wired, and the post-gate is replaced by a plain assessment that counts the
  connectors left without a certificate). Everything else is the same code.
  The difference in emitted circles is what connector repair costs; its
  budget is ``CIRCLE_BUDGET`` of the ablation's circles, decided exactly on
  integers.
- **Tier-2 replay cost.** The first ``TIER2_PROBE_CONNECTORS`` connectors of
  Figure 5 are replayed in emission order on both bracketing stocks, with
  the calls of ``tools.held_connector_engagement_tier2.replay_full``, timing
  each native operation kind. The mean costs give a prefix-rate estimate of
  a full tier-2 replay of Monstera's published path. Per-operation cost
  depends on the arrangement's size and on the load, so the estimate is
  the prefix rate scaled, not a forecast.

Modes: ``--measure`` (a few minutes) writes
``benchmarks/results/held_connector_repair_budget.json`` and renders; the
default re-renders the Markdown from that JSON. Both files are generated;
do not edit them. Times are seconds on a shared machine with the 1-minute
load average beside them.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
import pathlib
import time
from collections.abc import Iterator
from collections.abc import Sequence
from contextlib import contextmanager
from dataclasses import asdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Any
from unittest import mock

import numpy as np

import benchmarks.held_figure5_engagement_refinement as refinement_module
import benchmarks.held_toolpath_benchmark as toolpath_module
from benchmarks.errors import BenchmarkError
from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import moving_pieces
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_connector_gate import _ConnectorGate
from benchmarks.held_toolpath_benchmark import GENERATORS
from benchmarks.held_toolpath_benchmark import generate
from benchmarks.performance_record import REPOSITORY
from benchmarks.performance_record import build_identity
from benchmarks.performance_record import machine_identity
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.pipeline.admissibility import motion_parameters
from compas_cgal.units import ToolRadius
from tools.held_native_identity import capture_current_native
from tools.held_native_identity import native_identity_line
from tools.held_native_identity import native_record_fields
from tools.held_source_closure import capture
from tools.held_source_closure import digest_contents
from tools.held_source_closure import source_closure

RESULTS = REPOSITORY / "benchmarks" / "results"
JSON_PATH = RESULTS / "held_connector_repair_budget.json"
MARKDOWN_PATH = RESULTS / "held_connector_repair_budget.md"
# Budget for circles added by connector refinement, as a fraction of the
# ablation's circle count: set before the gate was wired (docs page, Cost
# model and budget). Exact, decided on integers.
CIRCLE_BUDGET = Fraction(2, 100)
# Connectors of the Figure 5 prefix the tier-2 cost probe replays: the
# design-stage windowing probe's prefix length.
TIER2_PROBE_WORKLOAD = "figure5_refined"
TIER2_PROBE_CONNECTORS = 200
# The workload the full-replay estimate is stated for (the largest).
TIER2_ESTIMATE_WORKLOAD = "figure8_monstera_boundary_qualified_corners"
# No gap closure for station verdicts (as the tier-2 tool).
NO_GAP_CLOSURE = 0.0
MS_PER_S = 1000.0


class BudgetRecordContractError(BenchmarkError):
    """The published path or the ablation breaks a premise of the comparison."""


@dataclass
class GateCalls:
    """How often the ablation's swapped connector gate was asked to decide: evidence that the swap reached the refinement."""

    certifies: int = 0
    sharp_vertex: int = 0


def _pair_cap_only_gate(calls: GateCalls) -> type[_ConnectorGate]:
    """Ablation gate class: accepts every connector, never refuses a sharp first arrival, counts its calls into ``calls``."""

    class _PairCapOnlyGate(_ConnectorGate):
        def certifies(self, emitted: Sequence[Figure5CounterclockwiseCircle], transition: Figure5BoundaryTransition) -> bool:
            calls.certifies += 1
            return True

        def sharp_vertex(self, emitted: Sequence[Figure5CounterclockwiseCircle], successor: Figure5CounterclockwiseCircle) -> int | None:
            calls.sharp_vertex += 1
            return None

    return _PairCapOnlyGate


@contextmanager
def connector_gate_off(left: list[int], calls: GateCalls) -> Iterator[None]:
    """Generate with the connector gate off; ``left`` receives each path's count of uncertified connectors, ``calls`` the swapped gate's calls."""

    def assess(circles: tuple[Figure5CounterclockwiseCircle, ...], connectors: Figure5BoundaryConnectors, tool: ToolRadius, cap: EngagementCap) -> None:
        left.append(len(ConnectorEngagementAssessment.build(circles, connectors, tool, cap).refused))

    with mock.patch.object(refinement_module, "_ConnectorGate", _pair_cap_only_gate(calls)), mock.patch.object(toolpath_module, "certify_connectors", assess):
        yield


def within_budget(added: int, ablation: int) -> bool:
    """Exact: ``added <= CIRCLE_BUDGET * ablation``."""
    return added * CIRCLE_BUDGET.denominator <= CIRCLE_BUDGET.numerator * ablation


def measure_budget(workload: str) -> dict[str, Any]:
    """One workload, published and with the connector gate off, in this process."""
    gc.collect()
    load = os.getloadavg()[0]
    published = generate(workload)
    left: list[int] = []
    calls = GateCalls()
    with connector_gate_off(left, calls):
        ablation = generate(workload)
    if len(left) != 1:
        raise BudgetRecordContractError(f"{workload}: the ablation ran its post-gate {len(left)} times")
    if calls.certifies == 0 or calls.sharp_vertex == 0:
        raise BudgetRecordContractError(f"{workload}: the refinement never called the swapped connector gate ({calls}); the ablation measured the published path")
    added = len(published.circles) - len(ablation.circles)
    if added < 0:
        raise BudgetRecordContractError(f"{workload}: the connector gate removed {-added} circles")
    return {
        "workload": workload,
        "circles_published": len(published.circles),
        "circles_gate_off": len(ablation.circles),
        "added": added,
        "within_budget": within_budget(added, len(ablation.circles)),
        "gate_off_left_uncertified": left[0],
        "gate_off_calls": asdict(calls),
        "generation_seconds_published": published.generation_seconds,
        "generation_seconds_gate_off": ablation.generation_seconds,
        "load_1min_before_after": [load, os.getloadavg()[0]],
    }


def _segment(a: Any, b: Any) -> tuple[float, float, float, float]:
    return float(a[0]), float(a[1]), float(b[0]), float(b[1])


def measure_tier2_rates(workload: str = TIER2_PROBE_WORKLOAD, connectors: int = TIER2_PROBE_CONNECTORS) -> dict[str, Any]:
    """Mean seconds per native operation over a replayed prefix, both stocks, as ``replay_full`` calls them."""
    path = generate(workload)
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    r = float(path.case.tool_radius.value)
    parameters = motion_parameters()
    region = np.array([(float(p.x), float(p.y), 0.0) for p in path.case.projection.points], dtype=np.float64)
    under, over = _stock_2.Stock2(region, []), _stock_2.Stock2(region, [])
    seconds = {"circle_under": 0.0, "circle_over": 0.0, "stations": 0.0, "piece_under": 0.0, "piece_over": 0.0}
    circles = pieces = stations = 0
    load = os.getloadavg()[0]
    for index, transition in enumerate(path.transitions[:connectors]):
        circle = path.circles[index]
        cx, cy, rho = float(circle.center.x), float(circle.center.y), float(circle.radius.value)
        started = time.perf_counter()
        under.subtract_circle_sweep_local(cx, cy, rho, r)
        seconds["circle_under"] += time.perf_counter() - started
        started = time.perf_counter()
        over.subtract_circle_sweep(cx, cy, rho, r)
        seconds["circle_over"] += time.perf_counter() - started
        circles += 1
        for a, b in moving_pieces(transition):
            x0, y0, x1, y1 = _segment(a, b)
            started = time.perf_counter()
            rows = _stock_2.engagement_line_prefixes(under, x0, y0, x1, y1, parameters, r, cap.chord_ratio, NO_GAP_CLOSURE)
            rows_over = _stock_2.engagement_line_prefixes(over, x0, y0, x1, y1, parameters, r, cap.chord_ratio, NO_GAP_CLOSURE)
            seconds["stations"] += time.perf_counter() - started
            stations += len(rows) + len(rows_over)
            started = time.perf_counter()
            under.subtract_line_sweep_local(x0, y0, x1, y1, r)
            seconds["piece_under"] += time.perf_counter() - started
            started = time.perf_counter()
            over.subtract_line_sweep_over(x0, y0, x1, y1, r)
            seconds["piece_over"] += time.perf_counter() - started
            pieces += 1
    return {
        "workload": workload,
        "connectors": min(connectors, len(path.transitions)),
        "circles": circles,
        "pieces": pieces,
        "stations_both_stocks": stations,
        "stations_per_piece": len(parameters),
        "seconds": seconds,
        "mean_ms": {
            "circle_under": MS_PER_S * seconds["circle_under"] / circles,
            "circle_over": MS_PER_S * seconds["circle_over"] / circles,
            "piece_under": MS_PER_S * seconds["piece_under"] / pieces,
            "piece_over": MS_PER_S * seconds["piece_over"] / pieces,
            "station": MS_PER_S * seconds["stations"] / stations,
        },
        "load_1min_before_after": [load, os.getloadavg()[0]],
    }


def tier2_estimate(rates: dict[str, Any], circles: int, pieces: int) -> dict[str, Any]:
    """REPORTING: prefix-rate seconds of a full replay of ``circles`` circles and ``pieces`` pieces on both stocks."""
    ms = rates["mean_ms"]
    subtraction = (circles * (ms["circle_under"] + ms["circle_over"]) + pieces * (ms["piece_under"] + ms["piece_over"])) / MS_PER_S
    stations = 2 * rates["stations_per_piece"] * pieces * ms["station"] / MS_PER_S
    return {"circles": circles, "pieces": pieces, "subtraction_seconds": subtraction, "station_seconds": stations, "total_seconds": subtraction + stations}


def measure(workloads: Sequence[str]) -> dict[str, Any]:
    """The whole record in one process; the build identity is read before and after and must agree."""
    before = build_identity()
    loaded = capture(source_closure())
    native_run = capture_current_native()
    rows = []
    for workload in workloads:
        rows.append(measure_budget(workload))
        print(f"budget {workload} done", flush=True)
    rates = measure_tier2_rates()
    print("tier-2 rates done", flush=True)
    target = generate(TIER2_ESTIMATE_WORKLOAD)
    pieces = sum(len(moving_pieces(transition)) for transition in target.transitions)
    estimate = {"workload": TIER2_ESTIMATE_WORKLOAD, **tier2_estimate(rates, len(target.circles), pieces)}
    if capture(loaded) != loaded:
        raise BudgetRecordContractError("the measured sources changed during the run; re-measure once edits settle")
    closure = source_closure()
    return {
        "build": asdict(before),
        "source_sha256": digest_contents(capture(closure)),
        "sources": list(closure),
        **native_record_fields(native_run),
        "machine": asdict(machine_identity()),
        "circle_budget": [CIRCLE_BUDGET.numerator, CIRCLE_BUDGET.denominator],
        "workloads": rows,
        "tier2_rates": rates,
        "tier2_estimate": estimate,
    }


def _percent(added: int, base: int) -> str:
    return f"{100.0 * added / base:+.2f} %"


def render(record: dict[str, Any]) -> str:
    """The Markdown fragment the docs page includes."""
    build = record["build"]
    budget = Fraction(*record["circle_budget"])
    lines = [
        "<!-- Generated by `pixi run held-connector-repair-budget` (tools/held_connector_repair_budget_results.py --measure); do not edit. -->",
        "",
        f"Recorded {build['recorded_at']} at commit `{build['commit'][:10]}`{' (dirty tree)' if build['dirty'] else ''}, "
        f"measured sources sha256 `{record['source_sha256'][:12]}`, {native_identity_line(record)}, on {record['machine']['cpu_brand']}.",
        "",
        "**Circles added by connector refinement.** Each workload generated as published and with the connector gate off (pair cap only, "
        "as before the gate was wired; NOTHING removal and connectors unchanged). Left uncertified: connectors of the gate-off path that "
        f"tiers 1 + 1b refuse, which the gate repairs. Budget: added ≤ {float(budget) * 100:g} % of the gate-off circles (exact).",
        "",
        "| workload | circles, gate off | circles, published | added | added / gate off | within budget | left uncertified, gate off | gen s (published / gate off) | load |",
        "|---|---:|---:|---:|---:|---|---:|---|---|",
    ]
    for row in record["workloads"]:
        load = row["load_1min_before_after"]
        lines.append(
            f"| `{row['workload']}` | {row['circles_gate_off']:,} | {row['circles_published']:,} | {row['added']:+d} | {_percent(row['added'], row['circles_gate_off'])} "
            f"| {'yes' if row['within_budget'] else '**NO**'} | {row['gate_off_left_uncertified']} "
            f"| {row['generation_seconds_published']:.1f} / {row['generation_seconds_gate_off']:.1f} | {load[0]:.0f}–{load[1]:.0f} |"
        )
    rates, estimate = record["tier2_rates"], record["tier2_estimate"]
    ms, load = rates["mean_ms"], rates["load_1min_before_after"]
    lines += [
        "",
        f"**Tier-2 replay cost.** First {rates['connectors']} connectors of `{rates['workload']}` replayed on both stocks "
        f"({rates['circles']} circles, {rates['pieces']} pieces, {rates['stations_both_stocks']:,} stations), load {load[0]:.0f}–{load[1]:.0f}. Mean ms per operation:",
        "",
        "| circle, UNDER (local) | circle, OVER (global) | piece, UNDER (local) | piece, OVER | station (one stock) |",
        "|---:|---:|---:|---:|---:|",
        f"| {ms['circle_under']:.2f} | {ms['circle_over']:.2f} | {ms['piece_under']:.2f} | {ms['piece_over']:.2f} | {ms['station']:.3f} |",
        "",
        f"Prefix-rate estimate of a full tier-2 replay of `{estimate['workload']}` ({estimate['circles']:,} circles, {estimate['pieces']:,} pieces, "
        f"{rates['stations_per_piece']} stations per piece on each stock): subtractions {estimate['subtraction_seconds']:.0f} s + stations "
        f"{estimate['station_seconds']:.0f} s ≈ {estimate['total_seconds']:.0f} s. Per-operation cost depends on the arrangement's size and on the "
        "load, so this is the prefix rate scaled, not a forecast.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--measure", action="store_true", help="re-measure (a few minutes) before rendering")
    parser.add_argument("--json", type=pathlib.Path, default=JSON_PATH)
    parser.add_argument("--markdown", type=pathlib.Path, default=MARKDOWN_PATH)
    parser.add_argument("--workloads", nargs="+", default=list(GENERATORS), choices=list(GENERATORS))
    args = parser.parse_args(argv)
    if args.measure:
        record = measure(args.workloads)
        args.json.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    record = json.loads(args.json.read_text(encoding="utf-8"))
    args.markdown.write_text(render(record), encoding="utf-8")
    print(args.markdown.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
