"""Contract of tools/held_omega_connector_repair_results.py: one generated Figure 5 row, and the rendering."""

from __future__ import annotations

from typing import Any

import pytest

from tools.held_omega_connector_repair_results import RepairRecordContractError
from tools.held_omega_connector_repair_results import check_repair_contract
from tools.held_omega_connector_repair_results import measure_workload
from tools.held_omega_connector_repair_results import render


@pytest.fixture(scope="module")
def figure5_row() -> dict[str, Any]:
    return measure_workload("figure5_refined")


def test_figure5_row_counts_the_repair(figure5_row: dict[str, Any]) -> None:
    row = figure5_row
    assert row["refused_before"] == len(row["repaired_connectors"]) >= 1
    assert row["refused_after"] == 0
    assert row["circles_after"] == row["circles_before"] + row["inserted"]
    assert row["inserted"] == sum(row["inserted_per_action"]) >= row["refused_before"]
    assert row["gouge_before"] == [0, 0] and row["gouge_after"] == [0, 0]
    assert row["fresh_tiers_equal_repair_verdicts"] is True


def test_figure5_row_splits_the_stage_seconds(figure5_row: dict[str, Any]) -> None:
    row = figure5_row
    assert 0 < row["repair_gate_seconds"] + row["repair_pass_seconds"] <= row["connector_repair_seconds"] < row["generation_seconds"]
    assert row["stage_over_rest"] == pytest.approx(row["connector_repair_seconds"] / (row["generation_seconds"] - row["connector_repair_seconds"]))


def test_contract_check_fails_loud_on_a_refusal_left_after_repair() -> None:
    check_repair_contract("w", refused_after=0, fresh_equal=True)
    with pytest.raises(RepairRecordContractError, match="w"):
        check_repair_contract("w", refused_after=1, fresh_equal=True)
    with pytest.raises(RepairRecordContractError, match="w"):
        check_repair_contract("w", refused_after=0, fresh_equal=False)


def test_render_writes_one_row_per_workload(figure5_row: dict[str, Any]) -> None:
    record = {"build": {"recorded_at": "t", "commit": "0123456789ab", "dirty": True}, "machine": {"cpu_brand": "cpu"}, "workloads": [figure5_row]}
    text = render(record)
    assert text.count("| `figure5_refined` |") == 1
    assert f"| {figure5_row['refused_before']} → 0 |" in text
    assert "(dirty tree)" in text
