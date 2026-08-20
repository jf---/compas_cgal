from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from scikit_build_core.settings.skbuild_read_settings import SettingsReader

from scripts.pytest_editable import PROBE_MODULE
from scripts.pytest_editable import EditableBuildDirNotFound
from scripts.pytest_editable import find_editable_build_dir


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


def test_strict_config_is_relaxed_for_pixi_only():
    """pip and cibuildwheel must keep scikit-build-core's typo detection armed.

    `strict-config = false` disables validation of this repo's own `[tool.scikit-build]` table
    and of `[project]` metadata, so a typo in e.g. `wheel.py-api` would yield a silently wrong
    wheel. It is gated behind pixi's `PIXI_PROJECT_NAME` so only pixi's own build relaxes it.
    """
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())

    def strict_config_under(env):
        return SettingsReader(manifest, {}, state="editable", env=env).settings.strict_config

    assert strict_config_under({}) is True, "pip/cibuildwheel builds must keep strict-config enabled"
    assert strict_config_under({"PIXI_PROJECT_NAME": "compas_cgal"}) is False, "pixi builds must relax it"


def test_missing_editable_finder_raises_a_named_error():
    """A failed build-dir lookup must name what it looked for, not raise bare StopIteration."""
    with pytest.raises(EditableBuildDirNotFound) as excinfo:
        find_editable_build_dir("compas_cgal._definitely_not_built")

    message = str(excinfo.value)
    assert "compas_cgal._definitely_not_built" in message
    assert "SKBUILD_EDITABLE_SKIP" in message


def test_editable_build_dir_resolves_to_a_real_directory():
    """The path exported as SKBUILD_EDITABLE_SKIP must be the finder's actual build directory."""
    assert Path(find_editable_build_dir(PROBE_MODULE)).is_dir()
