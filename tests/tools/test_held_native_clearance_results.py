"""Contract of tools/held_native_clearance_results.py and its committed record.

The record ``benchmarks/results/held_native_clearance.{json,md}`` backs the
claims of docs/held_native_clearance.md: native and Python oracles make
identical decisions on every workload and both paths, the Ω* paths audit to
0 gouges, and Ω* generation is byte-identical under either oracle. These tests
re-render the Markdown from the JSON without measuring and check the JSON
against those claims.
"""

from __future__ import annotations

import json
import statistics
import sys
from typing import Any

import pytest

from benchmarks.held_omega_toolpath import OMEGA_WORKLOADS
from benchmarks.performance_record import REPOSITORY
from tools.held_native_clearance_results import GENERATION_WORKLOADS
from tools.held_native_clearance_results import JSON_PATH
from tools.held_native_clearance_results import MARKDOWN_PATH
from tools.held_native_clearance_results import PATHS
from tools.held_native_clearance_results import main
from tools.held_native_clearance_results import render

DOC_PATH = REPOSITORY / "docs" / "held_native_clearance.md"
BACKENDS = ("python", "native")
COUNT_KEYS = ("gouging_circles", "gouging_connectors", "touching_circles", "touching_connectors")


@pytest.fixture(scope="module")
def record() -> dict[str, Any]:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def test_render_reproduces_the_committed_markdown_byte_for_byte(record: dict[str, Any]) -> None:
    assert render(record).encode("utf-8") == MARKDOWN_PATH.read_bytes()


def test_the_render_mode_cli_writes_the_committed_markdown(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "held_native_clearance.md"
    monkeypatch.setattr(sys, "argv", ["held_native_clearance_results", "--json", str(JSON_PATH), "--markdown", str(out)])
    main()
    assert out.read_bytes() == MARKDOWN_PATH.read_bytes()


def test_the_page_includes_the_rendered_record() -> None:
    assert '--8<-- "benchmarks/results/held_native_clearance.md"' in DOC_PATH.read_text(encoding="utf-8")


def test_record_schema(record: dict[str, Any]) -> None:
    assert set(record) == {"audits", "build", "generation", "machine"}
    assert {"commit", "dirty", "recorded_at"} <= set(record["build"])
    assert "cpu_brand" in record["machine"]
    for audit in record["audits"]:
        assert set(audit) == {"workload", *(key for key, _ in PATHS)}
        for key, _ in PATHS:
            part = audit[key]
            assert set(part) == {"python", "native", "comparison"}
            for backend in BACKENDS:
                assert set(COUNT_KEYS) <= set(part[backend])
                assert all(isinstance(part[backend][name], int) and part[backend][name] >= 0 for name in COUNT_KEYS)
                assert part[backend]["seconds"] > 0
                assert len(part[backend]["load_1min_before_after"]) == 2
    for item in record["generation"]:
        assert set(item) == {"workload", "repeats", "byte_identical", "median_seconds", "python_over_native", "runs"}


def test_every_workload_is_audited_on_both_paths(record: dict[str, Any]) -> None:
    assert [audit["workload"] for audit in record["audits"]] == list(OMEGA_WORKLOADS)


def test_native_and_python_decisions_are_identical(record: dict[str, Any]) -> None:
    for audit in record["audits"]:
        for key, _ in PATHS:
            part = audit[key]
            comparison = part["comparison"]
            assert comparison["identical_decisions"] is True, (audit["workload"], key)
            assert comparison["depth_bit_mismatches"] == 0
            assert comparison["first_mismatches"] == []
            assert comparison["listed_python"] == comparison["listed_native"]
            assert {name: part["python"][name] for name in COUNT_KEYS} == {name: part["native"][name] for name in COUNT_KEYS}


def test_omega_paths_audit_to_zero_gouges(record: dict[str, Any]) -> None:
    for audit in record["audits"]:
        native = audit["omega"]["native"]
        assert (native["gouging_circles"], native["gouging_connectors"]) == (0, 0), audit["workload"]


def test_generation_ab_is_byte_identical_and_its_summary_matches_its_runs(record: dict[str, Any]) -> None:
    assert [item["workload"] for item in record["generation"]] == list(GENERATION_WORKLOADS)
    for item in record["generation"]:
        runs = item["runs"]
        assert [run["backend"] for run in runs].count("python") == item["repeats"]
        assert [run["backend"] for run in runs].count("native") == item["repeats"]
        digests = {json.dumps(run["digests"], sort_keys=True) for run in runs}
        assert item["byte_identical"] is (len(digests) == 1)
        assert item["byte_identical"] is True, item["workload"]
        assert len({run["circles"] for run in runs}) == 1
        medians = {backend: statistics.median(run["seconds"] for run in runs if run["backend"] == backend) for backend in BACKENDS}
        assert item["median_seconds"] == medians
        assert item["python_over_native"] == medians["python"] / medians["native"]
