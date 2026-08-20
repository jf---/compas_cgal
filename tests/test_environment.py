from __future__ import annotations

import tomllib
from pathlib import Path


def test_pixi_manifest_declares_the_tasks_this_repo_documents():
    """CLAUDE.md mandates `pixi run baseline|pytest|affected|lint`; the manifest must provide them."""
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())
    tasks = manifest["tool"]["pixi"]["tasks"]
    assert {"baseline", "pytest", "affected", "lint"} <= set(tasks)


def test_editable_rebuild_only_imports_modules_this_branch_builds():
    """The rebuild probe must not reference branch-B extensions (_containment_2 etc.)."""
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())
    probe = manifest["tool"]["pixi"]["tasks"]["_editable-rebuild"]
    for absent in ("_containment_2", "_coverage_2", "_medial_axis_2"):
        assert absent not in probe
    assert "_stock_2" in probe and "_toolpath" in probe
