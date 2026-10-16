"""Stage acceptance rows select cases by fnmatch pattern (Rulings 1.1-C15, 1.1-C15-amend), and the stage-1 rows 1.1 adds hold."""

from __future__ import annotations

import json
from pathlib import Path

from benchmarks.stage_acceptance import ALL_CASES
from benchmarks.stage_acceptance import STAGE_ACCEPTANCE
from benchmarks.stage_acceptance import TestRow
from benchmarks.stage_acceptance import VectorRow
from benchmarks.stage_acceptance import matching_cases
from tests.benchmarks.test_stage_exit import _test_exists
from tests.benchmarks.test_stage_exit import _vector_row_failures
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import GateVectors
from tools.gate_vectors import load_gate_vectors

ROOT = Path(__file__).resolve().parents[2]

TOY_CASES = (
    "cap-120-default/engagement_controlled/rect_12x8",
    "cap-120-default/pipeline/rect_12x8",
    "cap-120-default/pipeline/L_shape",
    "cap-40-attribution/pipeline/rect_12x8",
    "cap-40-attribution/engagement_controlled/L_shape",
)


def test_a_pattern_row_selects_exactly_the_matching_keys() -> None:
    assert matching_cases("cap-120-default/pipeline/*", TOY_CASES) == ("cap-120-default/pipeline/rect_12x8", "cap-120-default/pipeline/L_shape")
    assert matching_cases("*/pipeline/rect_12x8", TOY_CASES) == ("cap-120-default/pipeline/rect_12x8", "cap-40-attribution/pipeline/rect_12x8")
    assert matching_cases(ALL_CASES, TOY_CASES) == TOY_CASES
    assert matching_cases("cap-120-default/pipeline/L_shape", TOY_CASES) == ("cap-120-default/pipeline/L_shape",)


def test_a_pattern_row_over_a_two_column_document_binds_only_its_cases() -> None:
    raw = json.loads((ROOT / GATE_VECTORS_PATH).read_text())
    base = next(iter(raw["cases"].values()))
    raw["cases"] = {
        "cap-120-default/engagement_controlled/rect_12x8": {**base, "redundant_operations": "7"},
        "cap-120-default/pipeline/rect_12x8": {**base, "redundant_operations": "0"},
    }
    raw["admissions"] = []
    vectors = GateVectors.build(raw)
    assert _vector_row_failures(vectors, VectorRow("cap-120-default/pipeline/*", "redundant_operations", "==", "0", "pattern")) == []
    assert _vector_row_failures(vectors, VectorRow("cap-120-default/*/rect_12x8", "redundant_operations", "==", "0", "both")) == [
        "cap-120-default/engagement_controlled/rect_12x8/redundant_operations = 7, required == 0 (both)"
    ]


def test_a_pattern_matching_no_case_is_a_failure_and_not_a_vacuous_pass() -> None:
    vectors = load_gate_vectors(ROOT / GATE_VECTORS_PATH)
    assert _vector_row_failures(vectors, VectorRow("cap-120-default/no_such_generator/*", "redundant_operations", "==", "0", "typo")) == [
        "cap-120-default/no_such_generator/* matches no case (typo)"
    ]


def test_the_stage_1_rows_task_1_1_adds_hold() -> None:
    """Ruling 1.1-C15: the stage-1 rows are checked here until the stage's exit test lands (1.7)."""
    vectors = load_gate_vectors(ROOT / GATE_VECTORS_PATH)
    rows = STAGE_ACCEPTANCE["stage-1"]
    failures = [failure for row in rows if isinstance(row, VectorRow) for failure in _vector_row_failures(vectors, row, stage="stage-1")]
    assert failures == []
    assert all(_test_exists(row.nodeid) for row in rows if isinstance(row, TestRow))
    assert {row.criterion for row in rows if isinstance(row, VectorRow)} >= {"material_plunges", "redundant_operations", "degenerate_loops", "tangent_breaks"}
