"""The native libraries a timing run measured: tracked from its start to its end, identified in its provenance."""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

from tools.held_native_identity import NativeBuildChangedError
from tools.run_native_identity import NativeIdentity
from tools.run_native_identity import NativeSourcesMismatchError
from tools.run_native_identity import native_sources_sha256
from tools.run_native_identity import require_built_from
from tools.run_native_identity import require_unchanged_natives
from tools.run_native_identity import track_native_here
from tools.run_native_identity import track_native_in

ROOT = Path(__file__).resolve().parents[2]
SHA256_HEX_CHARS = 64


def _tree(root: Path) -> Path:
    subprocess.run(["git", "init", "-q"], cwd=root, check=True, capture_output=True)
    (root / ".gitignore").write_text("build/\n")
    (root / "src").mkdir()
    (root / "src" / "a.cpp").write_text("int a;")
    (root / "src" / "b.h").write_text("int b();")
    (root / "src" / "wrapper.py").write_text("python, not native")
    (root / "CMakeLists.txt").write_text("project(x)")
    return root


def test_the_native_sources_digest_covers_cpp_headers_and_cmake_only(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    digest = native_sources_sha256(root)
    (root / "src" / "wrapper.py").write_text("edited")
    (root / "src" / "build").mkdir()
    (root / "src" / "build" / "generated.cpp").write_text("ignored by .gitignore")
    assert native_sources_sha256(root) == digest
    (root / "src" / "b.h").write_text("int b(int);")
    assert native_sources_sha256(root) != digest
    (root / "src" / "b.h").write_text("int b();")
    (root / "CMakeLists.txt").write_text("project(y)")
    assert native_sources_sha256(root) != digest


def test_a_build_from_other_native_sources_than_the_measured_tree_is_refused(tmp_path: Path) -> None:
    (tmp_path / "build_tree").mkdir()
    (tmp_path / "measured").mkdir()
    build, measured = _tree(tmp_path / "build_tree"), _tree(tmp_path / "measured")
    require_built_from(build, measured)
    (build / "src" / "parked").mkdir()
    (build / "src" / "parked" / "never_compiled.cpp").write_text("only in the build's tree")
    require_built_from(build, measured)
    (build / "src" / "a.cpp").write_text("int a = 1;")
    with pytest.raises(NativeSourcesMismatchError, match="a.cpp"):
        require_built_from(build, measured)


def test_an_extension_rewritten_between_the_two_reads_is_refused() -> None:
    at_start = {"build/x/_stock_2.so": hashlib.sha256(b"loaded").hexdigest()}
    require_unchanged_natives(at_start, {"build/x/_stock_2.so": b"loaded", "build/x/_lazy.so": b"loaded later"})
    with pytest.raises(NativeBuildChangedError):
        require_unchanged_natives(at_start, {"build/x/_stock_2.so": b"rebuilt"})
    with pytest.raises(NativeBuildChangedError):
        require_unchanged_natives(at_start, {})


def test_this_interpreter_s_native_run_is_identified_end_to_end() -> None:
    finish = track_native_here(ROOT)()
    identity = finish()
    assert isinstance(identity, NativeIdentity)
    assert len(identity.native_sha256) == len(identity.native_sources_sha256) == SHA256_HEX_CHARS
    assert identity.native_sources_sha256 == native_sources_sha256(ROOT)


def test_the_probe_in_another_interpreter_agrees_with_the_in_process_bracket() -> None:
    """The 1.1 base's bracket runs as a probe in the measured interpreter; on this interpreter it must give the in-process identity."""
    probed = track_native_in((sys.executable,), ROOT, dict(os.environ))()()
    assert probed == track_native_here(ROOT)()()
