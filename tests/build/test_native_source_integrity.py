from __future__ import annotations

import re
import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).parents[2]
NATIVE_SOURCE_VERIFIER = REPOSITORY_ROOT / "cmake" / "VerifyNativeSource.cmake"
DIGEST_PATTERN = re.compile(r"source-tree SHA256: ([0-9a-f]{64})")


def _write_fixture(source_dir: Path) -> None:
    include_dir = source_dir / "include"
    include_dir.mkdir(parents=True)
    (source_dir / ".locked-metadata").write_bytes(b"metadata\n")
    (include_dir / "exact.h").write_bytes(b"exact predicates\n")
    (include_dir / "kernel.h").write_bytes(b"exact constructions\n")


def _run_verifier(
    source_dir: Path,
    *,
    expected_digest: str | None = None,
    observed_digest_file: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    command = [
        "cmake",
        "-DNATIVE_PACKAGE=fixture",
        f"-DNATIVE_SOURCE_DIR={source_dir}",
        "-DNATIVE_PRINT_DIGEST=ON",
    ]
    if expected_digest is not None:
        command.append(f"-DNATIVE_EXPECTED_DIGEST={expected_digest}")
    if observed_digest_file is not None:
        command.append(f"-DNATIVE_OBSERVED_DIGEST_FILE={observed_digest_file}")
    command.extend(["-P", str(NATIVE_SOURCE_VERIFIER)])
    return subprocess.run(command, check=False, capture_output=True, text=True)


def _digest(source_dir: Path) -> str:
    result = _run_verifier(source_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    match = DIGEST_PATTERN.search(result.stdout + result.stderr)
    assert match is not None
    return match.group(1)


def test_independent_trees_have_identical_digest(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    _write_fixture(first)
    _write_fixture(second)

    assert _digest(first) == _digest(second)


def test_metadata_changes_do_not_change_digest(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    _write_fixture(source_dir)
    expected_digest = _digest(source_dir)
    header = source_dir / "include" / "exact.h"

    header.touch()
    header.chmod(0o600)

    assert _digest(source_dir) == expected_digest


def _modify_file(source_dir: Path) -> None:
    (source_dir / "include" / "exact.h").write_bytes(b"rounded predicates\n")


def _add_file(source_dir: Path) -> None:
    (source_dir / "unexpected.h").write_bytes(b"unexpected\n")


def _remove_file(source_dir: Path) -> None:
    (source_dir / "include" / "exact.h").unlink()


def _rename_file(source_dir: Path) -> None:
    (source_dir / "include" / "exact.h").rename(
        source_dir / "include" / "inexact.h"
    )


@pytest.mark.parametrize(
    "mutate",
    [_modify_file, _add_file, _remove_file, _rename_file],
    ids=["modify", "add", "remove", "rename"],
)
def test_content_mutations_fail_loudly(
    tmp_path: Path,
    mutate: Callable[[Path], None],
) -> None:
    source_dir = tmp_path / "source"
    _write_fixture(source_dir)
    expected_digest = _digest(source_dir)
    mutate(source_dir)

    result = _run_verifier(source_dir, expected_digest=expected_digest)
    diagnostic = result.stdout + result.stderr

    assert result.returncode != 0
    assert "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=fixture" in diagnostic
    assert f"expected={expected_digest}" in diagnostic
    assert "observed=" in diagnostic


def test_file_symlink_fails_loudly(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    _write_fixture(source_dir)
    (source_dir / "linked.h").symlink_to(source_dir / "include" / "exact.h")

    result = _run_verifier(source_dir)

    assert result.returncode != 0
    assert "expected=regular-file-tree observed=symlink:linked.h" in (
        result.stdout + result.stderr
    )


def test_root_symlink_fails_loudly(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    _write_fixture(source_dir)
    source_link = tmp_path / "source-link"
    source_link.symlink_to(source_dir, target_is_directory=True)

    result = _run_verifier(source_link)
    diagnostic = result.stdout + result.stderr

    assert result.returncode != 0
    assert "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=fixture" in diagnostic
    assert "expected=source-directory" in diagnostic
    assert "observed=root-symlink" in diagnostic


def test_observed_digest_is_persisted(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    _write_fixture(source_dir)
    expected_digest = _digest(source_dir)
    digest_file = tmp_path / "build" / "fixture.sha256"

    result = _run_verifier(
        source_dir,
        expected_digest=expected_digest,
        observed_digest_file=digest_file,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert digest_file.read_text(encoding="ascii") == f"{expected_digest}\n"
