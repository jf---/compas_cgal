"""Contract of tools/held_omega_coverage_results.py: Figure 5 against T(Ω*), the region definition, the rendering."""

from __future__ import annotations

from typing import Any

import pytest

from tools.held_omega_coverage_results import RegionAccountingError
from tools.held_omega_coverage_results import check_region_accounting
from tools.held_omega_coverage_results import measure_workload
from tools.held_omega_coverage_results import render


@pytest.fixture(scope="module")
def figure5_row() -> dict[str, Any]:
    return measure_workload("figure5_refined")


def test_figure5_is_covered_with_one_region_per_circle_and_per_omega_piece(figure5_row: dict[str, Any]) -> None:
    row = figure5_row
    assert (row["verdict"], row["witnesses"]) == ("COVERED", 0)
    assert row["regions"] == row["circles"] + row["omega_pieces"]
    assert row["connectors"] == row["circles"]
    assert row["sub_pieces"] > row["omega_pieces"]
    assert row["certificate_cpu_seconds"] > 0 and row["certificate_wall_seconds"] > 0


def test_region_accounting_fails_loud_on_the_sub_piece_count() -> None:
    check_region_accounting(1360, 1256, 104)
    with pytest.raises(RegionAccountingError):
        check_region_accounting(2616, 1256, 104)


def test_render_writes_one_row_per_workload(figure5_row: dict[str, Any]) -> None:
    uncovered = {**figure5_row, "workload": "u", "verdict": "UNCOVERED", "witnesses": 2, "witnesses_by_kind": {"ANNULUS_OUTER": 2}, "regions": None, "edges": None}
    record = {"build": {"recorded_at": "t", "commit": "0123456789ab", "dirty": False}, "machine": {"cpu_brand": "cpu"}, "workloads": [figure5_row, uncovered]}
    text = render(record)
    assert text.count("| `figure5_refined` | COVERED | 0 |") == 1
    assert "| `u` | UNCOVERED | 2 (ANNULUS_OUTER 2) |" in text
