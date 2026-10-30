"""The Omega* coverage record: the budget verdict is the production gate's, the diagnostic carries none, stale records are refused.

The verdict applies to the production workload (the Omega*-built Monstera
path, ``tools.coverage_certificate_budget_gate``), never to the diagnostic
polygonal-path run against T(Omega*), which exists to prove non-coverage.
Records of an older schema are refused, never rendered (strict render).
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest

from benchmarks.performance_record import BuildIdentity
from benchmarks.performance_record import MachineIdentity
from tools.coverage_certificate_budget_gate import BUDGET_RUNS
from tools.coverage_certificate_budget_gate import BUDGET_S
from tools.coverage_certificate_budget_gate import PRODUCTION_WORKLOAD
from tools.coverage_certificate_budget_gate import BudgetRun
from tools.coverage_certificate_budget_gate import CoverageVerdict
from tools.coverage_certificate_budget_gate import ProductionBudget
from tools.coverage_certificate_omega_star import SCHEMA
from tools.coverage_certificate_omega_star import ClusterRecord
from tools.coverage_certificate_omega_star import OmegaStarRecord
from tools.coverage_certificate_omega_star import StaleOmegaStarRecordError
from tools.coverage_certificate_omega_star import load_record
from tools.coverage_certificate_omega_star import render_markdown
from tools.held_generation_clock import MachineNotQuietError

QUIET = 1.0
LOADED = 450.0
# A wall time far above the budget: a wall-time gate would say OVER.
SLOW_WALL_S = 9.0


def _run(cpu_s: float, load: float = QUIET) -> BudgetRun:
    return BudgetRun(cpu_s=cpu_s, wall_s=SLOW_WALL_S, load_1min=load, certificate_load_1min=(load, load), verdict=CoverageVerdict.COVERED, witnesses=0, edges=1_210_313)


def _quiet_budget() -> ProductionBudget:
    return ProductionBudget.build(tuple(_run(c) for c in (4.9, 5.2, 5.3, 5.1, 5.0)), indicative=False)


def _record(budget: ProductionBudget, diagnostic_cpu_s: float = 9.0) -> OmegaStarRecord:
    cluster = ClusterRecord(
        centroid_mm=(45.45, 61.08),
        count=3,
        kinds=("ARC",),
        label="A",
        cluster_class="interior A",
        witnesses_by_feature=(0, 3, 0, 0),
        deepest_feature="convex arc",
        deepest_primitive=7,
        feature_radius_mm=2.5,
        depth_mm=1e-3,
        target_depth_mm=1e-4,
    )
    return OmegaStarRecord(
        generated_by="test",
        schema=SCHEMA,
        identity=MachineIdentity("arm64", "Darwin", "Apple M1 Max", "3.12"),
        build=BuildIdentity("0" * 40, False, "0.0", "2026-10-22T00:00:00+00:00"),
        budget=budget,
        cuts=4428,
        connectors=4148,
        target_pieces=487,
        vertex_disks=138,
        one_root_vertex_disks=138,
        target_build_s=0.9,
        attribution_s=1.0,
        region_count=1,
        piece_count=1,
        edge_count=1,
        witnesses=3,
        undecided=0,
        runs_wall_s=(diagnostic_cpu_s, diagnostic_cpu_s, diagnostic_cpu_s),
        runs_cpu_s=(diagnostic_cpu_s, diagnostic_cpu_s, diagnostic_cpu_s),
        runs_load=(QUIET, QUIET, QUIET),
        median_wall_s=diagnostic_cpu_s,
        median_cpu_s=diagnostic_cpu_s,
        native_phase_s=(1.0, 1.0, 1.0, 1.0),
        polygon_target_median_wall_s=2.4,
        polygon_target_witnesses=8,
        counters={"edges": 1},
        clusters=(cluster,),
    )


def test_an_indicative_record_renders_no_verdict() -> None:
    budget = ProductionBudget.build(tuple(_run(4.0, load=LOADED) for _ in range(BUDGET_RUNS)), indicative=True)
    text = render_markdown(_record(budget))
    assert "INDICATIVE" in text
    assert "within the" not in text
    assert "OVER the" not in text


def test_the_verdict_names_the_production_workload_and_the_diagnostic_carries_none() -> None:
    text = render_markdown(_record(_quiet_budget(), diagnostic_cpu_s=9.0))
    assert text.count(f"within the {BUDGET_S:g} s budget") == 1
    assert "OVER the" not in text
    verdict_line = next(line for line in text.splitlines() if f"within the {BUDGET_S:g} s budget" in line)
    assert PRODUCTION_WORKLOAD in verdict_line
    assert "process CPU" in verdict_line
    diagnostic_lines = [line for line in text.splitlines() if "polygonal centre domain" in line]
    assert diagnostic_lines
    assert all("Diagnostic" in line and "not the budget gate" in line for line in diagnostic_lines)


def test_a_record_round_trips_through_json(tmp_path: Path) -> None:
    record = _record(_quiet_budget())
    path = tmp_path / "record.json"
    path.write_text(record.to_json(), encoding="utf-8")
    assert load_record(path) == record


def test_a_schema_2_record_is_refused_not_rendered(tmp_path: Path) -> None:
    raw = json.loads(_record(_quiet_budget()).to_json())
    raw["schema"] = "coverage-certificate-omega-star/2"
    del raw["budget"]
    raw["budget_s"] = BUDGET_S
    path = tmp_path / "old.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(StaleOmegaStarRecordError):
        load_record(path)


def test_loading_re_validates_the_budget(tmp_path: Path) -> None:
    record = _record(_quiet_budget())
    loud = dataclasses.replace(record.budget, runs=tuple(dataclasses.replace(r, load_1min=LOADED) for r in record.budget.runs))
    path = tmp_path / "loud.json"
    path.write_text(dataclasses.replace(record, budget=loud).to_json(), encoding="utf-8")
    with pytest.raises(MachineNotQuietError):
        load_record(path)
