"""Contract of tools/held_omega_connector_gate_results.py after the connector repair stage: emitted and pre-repair gates."""

from __future__ import annotations

from typing import Any

import pytest

from tools.held_omega_connector_gate_results import GateRecordContractError
from tools.held_omega_connector_gate_results import check_gate_contract
from tools.held_omega_connector_gate_results import measure_workload
from tools.held_omega_connector_gate_results import render


@pytest.fixture(scope="module")
def figure5_row() -> dict[str, Any]:
    return measure_workload("figure5_refined")


def test_emitted_path_refuses_nothing_and_the_pre_repair_gate_names_the_repaired_connectors(figure5_row: dict[str, Any]) -> None:
    row, pre = figure5_row, figure5_row["pre_repair"]
    assert row["tiers"]["refused"] == 0 and row["refused_connectors"] == []
    assert pre["tiers"]["refused"] == len(pre["refused_connectors"]) >= 1
    assert pre["refused_connectors"] == row["repaired_connectors"]
    assert row["circles"] == pre["circles"] + row["repair_inserted"]
    assert pre["certified"] == pre["connectors"] - len(pre["refused_connectors"])
    assert pre["refused_sweeps_deg"][0] > row["cap_deg"]


def test_contract_check_fails_loud() -> None:
    check_gate_contract("w", emitted_refused=0, pre_refused=[3], repaired=[3], fresh_equal=True)
    for bad in (
        {"emitted_refused": 1, "pre_refused": [3], "repaired": [3], "fresh_equal": True},
        {"emitted_refused": 0, "pre_refused": [3], "repaired": [4], "fresh_equal": True},
        {"emitted_refused": 0, "pre_refused": [3], "repaired": [3], "fresh_equal": False},
    ):
        with pytest.raises(GateRecordContractError, match="w"):
            check_gate_contract("w", **bad)


def test_render_adds_the_pre_repair_table(figure5_row: dict[str, Any]) -> None:
    record = {"build": {"recorded_at": "t", "commit": "0123456789ab", "dirty": False}, "machine": {"cpu_brand": "cpu"}, "workloads": [figure5_row]}
    text = render(record)
    assert text.count("| `figure5_refined` |") == 3
    pre = figure5_row["pre_repair"]
    assert f"| `figure5_refined` before repair | {pre['certified']} / {pre['connectors']} |" in text
