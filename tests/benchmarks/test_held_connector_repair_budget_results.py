"""Contract of tools/held_connector_repair_budget_results.py: the ablation, the exact budget, the tier-2 probe, the rendering."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import pytest

import benchmarks.held_figure5_engagement_refinement as refinement_module
import benchmarks.held_toolpath_benchmark as toolpath_module
import tools.held_connector_repair_budget_results as budget_tool
from tools.held_connector_repair_budget_results import BudgetRecordContractError
from tools.held_connector_repair_budget_results import GateCalls
from tools.held_connector_repair_budget_results import connector_gate_off
from tools.held_connector_repair_budget_results import measure_budget
from tools.held_connector_repair_budget_results import measure_tier2_rates
from tools.held_connector_repair_budget_results import render
from tools.held_connector_repair_budget_results import tier2_estimate
from tools.held_connector_repair_budget_results import within_budget
from tools.held_native_identity import EditableBuild
from tools.held_native_identity import capture_native
from tools.held_native_identity import native_closure
from tools.held_native_identity import native_fields
from tools.held_source_closure import capture
from tools.held_source_closure import digest_contents
from tools.held_source_closure import source_closure


@pytest.fixture(scope="module")
def figure5_row() -> dict[str, Any]:
    return measure_budget("figure5_refined")


@pytest.fixture(scope="module")
def rates() -> dict[str, Any]:
    return measure_tier2_rates(connectors=3)


def test_budget_is_exact_at_the_boundary() -> None:
    assert within_budget(2, 100)
    assert not within_budget(3, 100)
    assert within_budget(0, 1)
    assert not within_budget(1, 49)


def test_the_ablation_restores_the_published_generator() -> None:
    gate, post_gate = refinement_module._ConnectorGate, toolpath_module.certify_connectors
    left: list[int] = []
    with connector_gate_off(left, GateCalls()):
        assert refinement_module._ConnectorGate is not gate
        assert toolpath_module.certify_connectors is not post_gate
    assert refinement_module._ConnectorGate is gate and toolpath_module.certify_connectors is post_gate


def test_figure5_row_compares_published_with_the_gate_off(figure5_row: dict[str, Any]) -> None:
    row = figure5_row
    assert row["circles_published"] == row["circles_gate_off"] + row["added"]
    assert row["added"] >= 0
    assert row["within_budget"] == within_budget(row["added"], row["circles_gate_off"])
    # The published path passes the post-gate, so a gate-off path that tiers 1 + 1b refuse must have needed circles.
    assert row["gate_off_left_uncertified"] == 0 or row["added"] > 0


def test_figure5_ablation_ran_through_the_swapped_gate(figure5_row: dict[str, Any]) -> None:
    # figure5_refined adds no circles, so only the call counts show that the swapped gate decided its transitions.
    assert figure5_row["gate_off_calls"]["certifies"] > 0
    assert figure5_row["gate_off_calls"]["sharp_vertex"] > 0


def test_a_swap_the_refinement_never_reaches_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    # Patch a stand-in module: the refinement keeps the published gate, the swapped gate is never called.
    monkeypatch.setattr(budget_tool, "refinement_module", SimpleNamespace(_ConnectorGate=refinement_module._ConnectorGate))
    with pytest.raises(BudgetRecordContractError, match="swapped connector gate"):
        measure_budget("figure5_refined")


def test_tier2_probe_counts_every_operation(rates: dict[str, Any]) -> None:
    assert rates["connectors"] == rates["circles"] == 3
    assert rates["pieces"] >= 1
    assert rates["stations_both_stocks"] == 2 * rates["stations_per_piece"] * rates["pieces"]
    assert all(value > 0 for value in rates["mean_ms"].values())


def test_estimate_scales_the_prefix_rates(rates: dict[str, Any]) -> None:
    one = tier2_estimate(rates, 1, 1)
    ten = tier2_estimate(rates, 10, 10)
    assert ten["total_seconds"] == pytest.approx(10 * one["total_seconds"])
    assert one["total_seconds"] == pytest.approx(one["subtraction_seconds"] + one["station_seconds"])


def test_render_writes_one_row_per_workload_and_the_identity(figure5_row: dict[str, Any], rates: dict[str, Any]) -> None:
    build = EditableBuild.locate()
    record = {
        "build": {"recorded_at": "t", "commit": "0123456789ab", "dirty": True},
        "source_sha256": digest_contents(capture(source_closure())),
        **native_fields(capture_native(build, native_closure(build))),
        "machine": {"cpu_brand": "cpu"},
        "circle_budget": [2, 100],
        "workloads": [figure5_row],
        "tier2_rates": rates,
        "tier2_estimate": {"workload": "w", **tier2_estimate(rates, 100, 100)},
    }
    text = render(record)
    assert text.count("| `figure5_refined` |") == 1
    assert "(dirty tree)" in text and "`0123456789`" in text
    assert record["source_sha256"][:12] in text
    assert record["native_sha256"][:12] in text
    assert "do not edit" in text.splitlines()[0]
