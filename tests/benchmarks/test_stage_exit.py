"""One exit-gate test per product-pipeline stage, reading `benchmarks.stage_acceptance`."""

from __future__ import annotations

import ast
import copy
import json
import re
from decimal import Decimal
from pathlib import Path

import pytest

from benchmarks.stage_acceptance import ALL_CASES
from benchmarks.stage_acceptance import STAGE_ACCEPTANCE
from benchmarks.stage_acceptance import TestRow
from benchmarks.stage_acceptance import UnknownAcceptanceCriterionError
from benchmarks.stage_acceptance import VectorRow
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import GateVectors
from tools.gate_vectors import load_gate_vectors

ROOT = Path(__file__).resolve().parents[2]
RED_MANIFEST = ROOT / "docs" / "red_manifest.json"

RELATIONS = {
    "<=": lambda value, bound: value <= bound,
    "==": lambda value, bound: value == bound,
    "<": lambda value, bound: value < bound,
}


def _test_exists(nodeid: str) -> bool:
    path, function = nodeid.split("::")
    module = ast.parse((ROOT / path).read_text())
    return any(isinstance(node, ast.FunctionDef) and node.name == function for node in module.body)


def _junit_identity(nodeid: str) -> str:
    path, function = nodeid.split("::")
    return f"{Path(path).with_suffix('').as_posix().replace('/', '.')}::{function}"


def _vector_row_failures(vectors: GateVectors, row: VectorRow) -> list[str]:
    """Check one row against the current vectors by criterion name and bound (Ruling 0.3-exit-unpinned).

    Raises:
        UnknownAcceptanceCriterionError: The vectors carry no criterion of the
            row's name; a bump that changes a row's criterion edits the row in
            the same commit.
    """
    if row.criterion not in {criterion.name for criterion in vectors.criteria}:
        raise UnknownAcceptanceCriterionError(f"acceptance row criterion {row.criterion!r} is not in the gate vectors (rules {vectors.rules_version}); edit the row ({row.why}).")
    failures: list[str] = []
    cases = list(vectors.cases) if row.case == ALL_CASES else [row.case]
    for case in cases:
        value = vectors.value(case, row.criterion)
        if not RELATIONS[row.relation](value, Decimal(row.bound)):
            failures.append(f"{case}/{row.criterion} = {value}, required {row.relation} {row.bound} ({row.why})")
    return failures


def _assert_stage_exits(stage: str, vectors: GateVectors) -> None:
    manifest = [re.compile(entry["match"]) for entry in json.loads(RED_MANIFEST.read_text())["expected_red"]]
    failures: list[str] = []
    for row in STAGE_ACCEPTANCE[stage]:
        if isinstance(row, TestRow):
            if not _test_exists(row.nodeid):
                failures.append(f"missing acceptance test {row.nodeid} ({row.why})")
            identity = _junit_identity(row.nodeid)
            if any(pattern.match(identity) or pattern.match(f"{identity}[") for pattern in manifest):
                failures.append(f"acceptance test {row.nodeid} is a declared red ({row.why})")
            continue
        assert isinstance(row, VectorRow)
        failures.extend(_vector_row_failures(vectors, row))
    assert failures == []


def test_stage_0_exits() -> None:
    """Stage 0: the gate measures physics (0.1), sweeps are exact (0.2), rules and ratchet (0.3)."""
    _assert_stage_exits("stage-0", load_gate_vectors(ROOT / GATE_VECTORS_PATH))


def test_a_vector_row_naming_a_criterion_the_vectors_lack_fails_by_name() -> None:
    """Ruling 0.3-exit-unpinned: rows bind by criterion name; a vanished name is loud, never a skip."""
    vectors = load_gate_vectors(ROOT / GATE_VECTORS_PATH)
    row = VectorRow(ALL_CASES, "no_such_criterion", "==", "0", "a row whose criterion a bump renamed")
    with pytest.raises(UnknownAcceptanceCriterionError, match="no_such_criterion"):
        _vector_row_failures(vectors, row)


def test_the_exit_rows_are_not_pinned_to_a_rules_version() -> None:
    """Ruling 0.3-exit-unpinned: the rows read the current vectors whatever their rules version."""
    raw = json.loads((ROOT / GATE_VECTORS_PATH).read_text())
    bumped = copy.deepcopy(raw)
    bumped["rules_version"] = "99.0"
    bumped["rebaseline"] = {
        "from": raw["rules_version"],
        "approval": "Ruling 99.0-X: synthetic",
        "criteria": [next(name for name, spec in raw["criteria"].items() if spec["direction"] == "lower")],
        "retired": [],
    }
    _assert_stage_exits("stage-0", GateVectors.build(bumped))
