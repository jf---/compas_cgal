from __future__ import annotations

import hashlib
import io
import re
import shutil
import subprocess
import tarfile
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).parents[2]
NATIVE_SOURCE_VERIFIER = REPOSITORY_ROOT / "cmake" / "VerifyNativeSource.cmake"
NATIVE_DEPENDENCIES = REPOSITORY_ROOT / "cmake" / "NativeDependencies.cmake"
DIGEST_PATTERN = re.compile(r"source-tree SHA256: ([0-9a-f]{64})")
FIXTURE_FILES = {
    ".locked-metadata": b"metadata\n",
    "include/exact.h": b"exact predicates\n",
    "include/kernel.h": b"exact constructions\n",
}
# Same version range as the real CMakeLists.txt: its upper bound is what sets CMP0135 NEW
# (extract with --touch), which the configure-depends on the source directory relies on.
CMAKE_MINIMUM_REQUIRED = "3.15...3.26"


@pytest.fixture
def scratch() -> Iterator[Path]:
    # A fresh configure, under the repository's gitignored build/ (never /tmp).
    root = REPOSITORY_ROOT / "build" / "cmake-configure-test" / uuid.uuid4().hex
    root.mkdir(parents=True)
    yield root
    shutil.rmtree(root)


def _write_archive(archive: Path) -> str:
    with tarfile.open(archive, "w:gz") as bundle:
        for relative_path, contents in FIXTURE_FILES.items():
            entry = tarfile.TarInfo(f"fixture-1.0/{relative_path}")
            entry.size = len(contents)
            bundle.addfile(entry, io.BytesIO(contents))
    return hashlib.sha256(archive.read_bytes()).hexdigest()


def _source_tree_digest(reference_dir: Path) -> str:
    for relative_path, contents in FIXTURE_FILES.items():
        path = reference_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)
    result = subprocess.run(
        [
            "cmake",
            "-DNATIVE_PACKAGE=fixture",
            f"-DNATIVE_SOURCE_DIR={reference_dir}",
            "-DNATIVE_PRINT_DIGEST=ON",
            "-P",
            str(NATIVE_SOURCE_VERIFIER),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    match = DIGEST_PATTERN.search(result.stdout + result.stderr)
    assert match is not None
    return match.group(1)


def _run(command: list[str]) -> str:
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    return output


def _has_target(build_dir: Path, target: str) -> bool:
    result = subprocess.run(
        ["ninja", "-C", str(build_dir), "-t", "query", target],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def test_sources_landing_after_configure_reconfigure_and_settle(scratch: Path) -> None:
    """A fresh worktree: external/ is empty at configure and populated by the first build.

    The second build must re-run CMake (the source directory is a configure input), take the
    populated branch, verify at configure time and drop the always-dirty verify target, so a
    third invocation has no work to do.
    """
    archive = scratch / "fixture-1.0.tar.gz"
    archive_sha256 = _write_archive(archive)
    source_tree_sha256 = _source_tree_digest(scratch / "reference")
    source_dir = scratch / "external" / "fixture"
    project_dir = scratch / "project"
    project_dir.mkdir()
    (project_dir / "CMakeLists.txt").write_text(
        f"""\
cmake_minimum_required(VERSION {CMAKE_MINIMUM_REQUIRED})
project(native_source_reconfigure NONE)
include(ExternalProject)
add_custom_target(external_downloads ALL)
include("{NATIVE_DEPENDENCIES.as_posix()}")
native_add_locked_dependency(
    NAME fixture
    SOURCE_DIR "{source_dir.as_posix()}"
    VERSION "1.0"
    URL "{archive.as_uri()}"
    ARCHIVE_SHA256 "{archive_sha256}"
    SOURCE_TREE_SHA256 "{source_tree_sha256}"
)
""",
        encoding="utf-8",
    )
    build_dir = scratch / "build"

    _run(["cmake", "-S", str(project_dir), "-B", str(build_dir), "-G", "Ninja"])
    assert _has_target(build_dir, "verify_fixture_source")

    _run(["ninja", "-C", str(build_dir)])
    assert (source_dir / "include" / "exact.h").read_bytes() == FIXTURE_FILES["include/exact.h"]

    second_build = _run(["ninja", "-C", str(build_dir)])
    assert "Re-running CMake" in second_build
    assert "VerifyNativeSource.cmake" not in second_build
    assert not _has_target(build_dir, "verify_fixture_source")

    settle = _run(["ninja", "-C", str(build_dir), "-n"])
    assert "no work to do" in settle
