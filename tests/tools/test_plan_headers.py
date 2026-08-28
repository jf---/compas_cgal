from __future__ import annotations

import pathlib
import subprocess
import sys

from tools.plan_headers import missing_headers

PROJECT_ROOT = pathlib.Path(__file__).parents[2]


def _plan(plans_dir: pathlib.Path, name: str, text: str) -> None:
    (plans_dir / name).write_text(text)


def _run_cli(
    plans_dir: pathlib.Path,
    *arguments: str,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.plan_headers",
            "--plans-dir",
            str(plans_dir),
            *arguments,
        ],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_a_stamped_plan_passes(tmp_path: pathlib.Path) -> None:
    _plan(tmp_path, "a.md", "# A\n\n> **status: landed** — evidence.\n\nbody\n")

    assert missing_headers(tmp_path) == []


def test_an_unstamped_plan_is_named(tmp_path: pathlib.Path) -> None:
    _plan(tmp_path, "b.md", "# B\n\nbody with no header\n")

    assert [plan.name for plan in missing_headers(tmp_path)] == ["b.md"]


def test_a_header_after_the_first_twelve_lines_does_not_count(
    tmp_path: pathlib.Path,
) -> None:
    _plan(tmp_path, "c.md", "# C\n" + "\n" * 11 + "> **status: landed**\n")

    assert [plan.name for plan in missing_headers(tmp_path)] == ["c.md"]


def test_an_explicit_function_prefix_exempts_a_plan(tmp_path: pathlib.Path) -> None:
    _plan(tmp_path, "owned-a.md", "# Owned\n")
    _plan(tmp_path, "unowned.md", "# Unowned\n")

    assert [plan.name for plan in missing_headers(tmp_path, allow_prefixes=("owned-",))] == ["unowned.md"]


def test_cli_builtin_prefix_exempts_the_auditor_fleet(
    tmp_path: pathlib.Path,
) -> None:
    _plan(
        tmp_path,
        "2026-08-23-auditor-convergence-p0-reconciliation.md",
        "# Auditor-owned plan\n",
    )

    result = _run_cli(tmp_path)

    assert result.returncode == 0
    assert result.stdout == ""
    assert result.stderr == ""


def test_cli_prefixes_are_additive_and_keep_the_builtin_exemption(
    tmp_path: pathlib.Path,
) -> None:
    _plan(
        tmp_path,
        "2026-08-23-auditor-convergence-p4-enforcement-evidence.md",
        "# Auditor-owned plan\n",
    )
    _plan(tmp_path, "alpha-plan.md", "# Alpha\n")
    _plan(tmp_path, "beta-plan.md", "# Beta\n")
    _plan(tmp_path, "gamma-plan.md", "# Gamma\n")

    result = _run_cli(
        tmp_path,
        "--allow-prefix",
        "alpha-",
        "--allow-prefix",
        "beta-",
    )

    assert result.returncode == 1
    assert result.stdout == f"missing status header: {tmp_path / 'gamma-plan.md'}\n"
    assert result.stderr == ""
