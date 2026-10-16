"""One exit-gate test per product-pipeline stage, reading `benchmarks.stage_acceptance`."""

from __future__ import annotations

import ast
import copy
import json
import re
from decimal import Decimal
from pathlib import Path

import pytest

from benchmarks.gate import FROZEN_PAIR
from benchmarks.stage_acceptance import ALL_CASES
from benchmarks.stage_acceptance import STAGE_ACCEPTANCE
from benchmarks.stage_acceptance import TestRow
from benchmarks.stage_acceptance import UnknownAcceptanceCriterionError
from benchmarks.stage_acceptance import VectorRow
from benchmarks.stage_acceptance import matching_cases
from tools.gate_admissions import EXEMPTABLE_STAGE
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


def _vector_row_failures(vectors: GateVectors, row: VectorRow, *, stage: str | None = None) -> list[str]:
    """Check one row against the current vectors by criterion name and bound (Ruling 0.3-exit-unpinned).

    Raises:
        UnknownAcceptanceCriterionError: The vectors carry no criterion of the
            row's name; a bump that changes a row's criterion edits the row in
            the same commit.
    """
    if row.criterion not in {criterion.name for criterion in vectors.criteria}:
        raise UnknownAcceptanceCriterionError(f"acceptance row criterion {row.criterion!r} is not in the gate vectors (rules {vectors.rules_version}); edit the row ({row.why}).")
    failures: list[str] = []
    # G4 (Ruling 0.3-C2-amend): a stage-0 all-case row skips exactly the
    # (case, criterion) pairs an admission's exempt_rows lists; nothing else,
    # and no other stage's rows.
    exempt = vectors.exempt_rows if stage == EXEMPTABLE_STAGE and row.case == ALL_CASES else {}
    matched = matching_cases(row.case, tuple(vectors.cases))
    if not matched:
        return [f"{row.case} matches no case ({row.why})"]
    cases = [case for case in matched if row.criterion not in exempt.get(case, frozenset())]
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
        failures.extend(_vector_row_failures(vectors, row, stage=stage))
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


def _frozen_pair_document() -> dict:
    """The committed vectors restricted to the frozen pair's cases, with no admission: a base the exemption tests add to."""
    raw = json.loads((ROOT / GATE_VECTORS_PATH).read_text())
    raw["cases"] = {key: cells for key, cells in raw["cases"].items() if key.split("/")[1] in FROZEN_PAIR}
    raw["admissions"] = []
    return raw


def test_an_exempt_row_skips_exactly_its_admitted_case_and_criterion() -> None:
    """G4 (Ruling 0.3-C2-amend): an admission's ``exempt_rows`` lifts named (case, criterion) pairs; every other row binds."""
    raw = _frozen_pair_document()
    sibling = next(iter(raw["cases"]))
    exempt = "cap-120-default/rho_regulated/L_shape"
    raw["cases"][exempt] = {**raw["cases"][sibling], "misplaced_plunges": "1"}
    raw["admissions"] = [{"ruling": "0.3-C2-amend", "cases": [exempt], "reference": {"rho_regulated": None}, "exempt_rows": {exempt: ["misplaced_plunges"]}, "excluded": {}}]
    misplaced = VectorRow(ALL_CASES, "misplaced_plunges", "==", "0", "rule (e)")
    unsafe = VectorRow(ALL_CASES, "unsafe_rapids", "==", "0", "rapids")
    assert _vector_row_failures(GateVectors.build(raw), misplaced, stage="stage-0") == []
    raw["cases"][exempt]["unsafe_rapids"] = "1"
    assert _vector_row_failures(GateVectors.build(raw), unsafe, stage="stage-0") == [f"{exempt}/unsafe_rapids = 1, required == 0 (rapids)"]
    raw["cases"][sibling]["misplaced_plunges"] = "1"
    assert _vector_row_failures(GateVectors.build(raw), misplaced, stage="stage-0") == [f"{sibling}/misplaced_plunges = 1, required == 0 (rule (e))"]


def test_exemptions_apply_only_to_the_stage_0_rows() -> None:
    """An admission's exempt_rows lifts stage-0 '*' rows; a later stage's '*' row on the same criterion still binds."""
    raw = _frozen_pair_document()
    exempt = "cap-120-default/rho_regulated/L_shape"
    raw["cases"][exempt] = {**raw["cases"][next(iter(raw["cases"]))], "misplaced_plunges": "1"}
    raw["admissions"] = [{"ruling": "0.3-C2-amend", "cases": [exempt], "reference": {"rho_regulated": None}, "exempt_rows": {exempt: ["misplaced_plunges"]}, "excluded": {}}]
    row = VectorRow(ALL_CASES, "misplaced_plunges", "==", "0", "a later stage")
    vectors = GateVectors.build(raw)
    assert _vector_row_failures(vectors, row, stage="stage-0") == []
    assert _vector_row_failures(vectors, row, stage="stage-1") == [f"{exempt}/misplaced_plunges = 1, required == 0 (a later stage)"]
    assert _vector_row_failures(vectors, VectorRow(exempt, "misplaced_plunges", "==", "0", "named"), stage="stage-0") == [f"{exempt}/misplaced_plunges = 1, required == 0 (named)"]
