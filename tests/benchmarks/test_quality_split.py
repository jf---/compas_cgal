"""The H1 split of `benchmarks.quality_observations` keeps every moved name importable from its old path."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import benchmarks.quality_binding as quality_binding
import benchmarks.quality_criteria as quality_criteria
import benchmarks.quality_observations as quality_observations
from tools.gate_ratchet import GATE_PATHS

REPO_ROOT = Path(__file__).resolve().parents[2]

# Every moved name an existing importer reads from `benchmarks.quality_observations`.
CRITERIA_REEXPORTS = (
    "CRITERION_NAMES",
    "EVIDENCE_BY_CRITERION",
    "CountCriterion",
    "DegreesCriterion",
    "FractionCriterion",
    "MeasuredStep",
    "OperationPair",
    "PathQualityAttribution",
    "ToolRadiusMultipleCriterion",
    "_build_record",
)
BINDING_REEXPORTS = (
    "HelixHasNoMotionKindError",
    "_snapshot_cut_height",
    "_snapshot_replay_category",
    "snapshot_motion_kind",
)


@pytest.mark.parametrize("name", CRITERIA_REEXPORTS)
def test_criteria_names_resolve_to_one_object(name: str) -> None:
    assert getattr(quality_observations, name) is getattr(quality_criteria, name)


@pytest.mark.parametrize("name", BINDING_REEXPORTS)
def test_binding_names_resolve_to_one_object(name: str) -> None:
    assert getattr(quality_observations, name) is getattr(quality_binding, name)


def test_assessment_stays_with_its_reductions() -> None:
    assert quality_observations.PathQualityAssessment.__module__ == "benchmarks.quality_observations"


@pytest.mark.parametrize("path", ["benchmarks/quality_criteria.py", "benchmarks/quality_binding.py"])
def test_split_modules_are_gate_inputs(path: str) -> None:
    assert path in GATE_PATHS


def _benchmark_imports(module_path: Path) -> set[str]:
    tree = ast.parse(module_path.read_text())
    return {node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom) and node.module is not None and node.module.startswith("benchmarks.quality")}


def test_split_layers_are_acyclic() -> None:
    # binding <- criteria <- observations: no module imports a layer above it.
    assert _benchmark_imports(REPO_ROOT / "benchmarks/quality_binding.py") == set()
    assert _benchmark_imports(REPO_ROOT / "benchmarks/quality_criteria.py") == {"benchmarks.quality_binding"}
