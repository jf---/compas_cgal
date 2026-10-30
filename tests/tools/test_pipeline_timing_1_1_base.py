"""The 1.1 base wrapper: times the unmodified 1.1 tree through its own child mode, reads the load, lifts to pipeline-timing/2."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import pytest

from tests.tools.test_run_provenance import provenance
from tools.held_native_identity import NativeBuildChangedError
from tools.run_native_identity import NativeFinisher
from tools.run_native_identity import NativeIdentity
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import LoadAverage
from tools.held_generation_clock import MachineNotQuietError
from tools.pipeline_timing import CaseProcessFailedError
from tools.pipeline_timing import CaseRun
from tools.pipeline_timing import NondeterministicCaseError
from tools.pipeline_timing_1_1_base import BASE_RESULT_LINE_PREFIX
from tools.pipeline_timing_1_1_base import BASE_TOOL_SOURCES
from tools.pipeline_timing_1_1_base import REBUILD_TASK
from tools.pipeline_timing_1_1_base import BaseRebuildFailedError
from tools.pipeline_timing_1_1_base import rebuild_base_build
from tools.pipeline_timing_1_1_base import EXPECTED_1_1_COMMIT
from tools.pipeline_timing_1_1_base import ForeignLibraryError
from tools.pipeline_timing_1_1_base import base_tool_sha256
from tools.pipeline_timing_1_1_base import editable_build_dir_of
from tools.pipeline_timing_1_1_base import expected_base_tool_sha256
from tools.pipeline_timing_1_1_base import library_origin_of
from tools.pipeline_timing_1_1_base import require_library_in_tree
from tools.pipeline_timing_1_1_base import MeasuredTreeError
from tools.pipeline_timing_1_1_base import base_child_env
from tools.pipeline_timing_1_1_base import base_keys
from tools.pipeline_timing_1_1_base import base_launcher
from tools.pipeline_timing_1_1_base import measure_1_1_base
from tools.pipeline_timing_1_1_base import require_measured_tree
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import SCHEMA
from tools.pipeline_timing_results import PipelineTimings
from tools.run_provenance import RunProvenance
from tools.run_provenance import TreeChangedDuringRunError
from tools.run_provenance import TreeIdentity

ROOT = Path(__file__).resolve().parents[2]
PYTHON = ("pixi", "run", "python")
KEYS = ("cap-120-default/pipeline/rect_12x8", "cap-120-default/pipeline/L_shape")
EXCLUDED = ("cap-40-attribution/pipeline/L_shape",)
DIGEST = "0123456789abcdef"
QUIET: LoadAverage = (0.5, 0.75, 1.25)
LOADED: LoadAverage = (QUIET_LOAD_AVERAGE_MAX + 1.0, 1.0, 1.0)


TREE = (EXPECTED_1_1_COMMIT, False, None)


def _tree() -> TreeIdentity:
    return TREE


NATIVE = NativeIdentity("1" * 64, "2" * 64)


def _track_native() -> NativeFinisher:
    return lambda: NATIVE


def _native_build_changed() -> NativeFinisher:
    def finish() -> NativeIdentity:
        raise NativeBuildChangedError("a native extension changed on disk during the run")

    return finish


def _describe(tree: TreeIdentity, native: NativeIdentity, started_utc: str, ended_utc: str) -> RunProvenance:
    return provenance(rev=tree[0], native_sha256=native.native_sha256, native_sources_sha256=native.native_sources_sha256, started_utc=started_utc, ended_utc=ended_utc)


def _no_rebuild() -> None:
    return None


def _result_line(seconds: float = 0.5, digest: str = DIGEST) -> str:
    return f'noise\n{BASE_RESULT_LINE_PREFIX}{{"generation_s": {seconds}, "operations": 7, "digest": "{digest}"}}\n'


def test_the_1_1_tool_prints_the_result_prefix_the_wrapper_parses() -> None:
    source = subprocess.run(["git", "show", f"{EXPECTED_1_1_COMMIT}:tools/pipeline_timing.py"], capture_output=True, text=True, check=True).stdout
    assert f'RESULT_LINE_PREFIX = "{BASE_RESULT_LINE_PREFIX}"' in source
    assert "--case" in source


def test_the_case_keys_come_from_the_1_1_tree_itself() -> None:
    commands: list[Sequence[str]] = []

    def run(command: Sequence[str]) -> str:
        commands.append(command)
        return f"[{list(KEYS)!r}, {list(EXCLUDED)!r}]\n".replace("'", '"')

    assert base_keys(run, PYTHON) == (KEYS, EXCLUDED)
    (command,) = commands
    assert tuple(command[: len(PYTHON)]) == PYTHON and command[len(PYTHON)] == "-c"
    assert "from tools.pipeline_timing import gate_timing_keys" in command[-1]


def test_each_run_is_the_1_1_tools_own_child_mode() -> None:
    commands: list[Sequence[str]] = []

    def run(command: Sequence[str]) -> str:
        commands.append(command)
        return _result_line()

    run_ = base_launcher(run, PYTHON)(KEYS[0])
    assert list(commands[0]) == [*PYTHON, "-m", "tools.pipeline_timing", "--case", KEYS[0]]
    assert (run_.generation_s, run_.digest) == (0.5, DIGEST)


def test_the_base_is_timed_serially_with_its_load_and_provenance_and_lifted_to_the_load_schema() -> None:
    launched: list[str] = []

    def launch(key: str) -> CaseRun:
        launched.append(key)
        return base_launcher(lambda _command: _result_line(), PYTHON)(key)

    timings = measure_1_1_base(
        KEYS, excluded=EXCLUDED, launch=launch, load_reader=lambda: QUIET, describe=_describe, rebuild=_no_rebuild, track_native=_track_native, read_tree=_tree
    )
    assert launched == [key for key in KEYS for _run in range(RUNS_PER_CASE)]
    assert (timings.load_start, timings.load_end) == (QUIET, QUIET)
    assert timings.provenance is not None and timings.provenance.rev == EXPECTED_1_1_COMMIT
    assert timings.to_json()["schema"] == SCHEMA
    assert isinstance(timings, PipelineTimings)


def test_a_case_on_a_loaded_machine_is_refused_before_it_runs() -> None:
    loads = iter((QUIET, QUIET, LOADED))
    launched: list[str] = []

    def launch(key: str) -> CaseRun:
        launched.append(key)
        return base_launcher(lambda _command: _result_line(), PYTHON)(key)

    with pytest.raises(MachineNotQuietError) as refused:
        measure_1_1_base(
            KEYS, excluded=EXCLUDED, launch=launch, load_reader=lambda: next(loads), describe=_describe, rebuild=_no_rebuild, track_native=_track_native, read_tree=_tree
        )
    assert launched == [KEYS[0]] * RUNS_PER_CASE
    assert "--indicative" in str(refused.value)


def test_an_indicative_base_run_measures_under_load() -> None:
    timings = measure_1_1_base(
        KEYS,
        excluded=EXCLUDED,
        launch=base_launcher(lambda _command: _result_line(), PYTHON),
        load_reader=lambda: LOADED,
        describe=_describe,
        rebuild=_no_rebuild,
        track_native=_track_native,
        read_tree=_tree,
        indicative=True,
    )
    assert timings.load_start == LOADED


def test_runs_that_disagree_are_refused() -> None:
    digests = iter((DIGEST, DIGEST, "fedcba9876543210"))
    with pytest.raises(NondeterministicCaseError):
        measure_1_1_base(
            KEYS[:1],
            excluded=(),
            launch=base_launcher(lambda _command: _result_line(digest=next(digests)), PYTHON),
            load_reader=lambda: QUIET,
            describe=_describe,
            rebuild=_no_rebuild,
            track_native=_track_native,
            read_tree=_tree,
        )


def _repo(tmp_path: Path) -> str:
    def git(*args: str) -> None:
        subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    (tmp_path / "f.txt").write_text("one")
    git("add", "f.txt")
    git("commit", "-q", "-m", "one")
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout.strip()


def test_a_tree_at_another_commit_is_refused(tmp_path: Path) -> None:
    _repo(tmp_path)
    with pytest.raises(MeasuredTreeError):
        require_measured_tree(tmp_path, EXPECTED_1_1_COMMIT)


def test_a_modified_tree_is_refused(tmp_path: Path) -> None:
    rev = _repo(tmp_path)
    require_measured_tree(tmp_path, rev)
    (tmp_path / "f.txt").write_text("patched")
    with pytest.raises(MeasuredTreeError):
        require_measured_tree(tmp_path, rev)


def test_children_run_the_1_1_tree_with_its_own_build_pinned(tmp_path: Path) -> None:
    build = tmp_path / "build" / "editable"
    env = base_child_env(tmp_path, {"PYTHONPATH": ".", "SKBUILD_EDITABLE_SKIP": "/elsewhere", "HOME": "/home"}, build)
    assert env["PYTHONPATH"] == str(tmp_path)
    assert env["SKBUILD_EDITABLE_SKIP"] == str(build)
    assert env["HOME"] == "/home"
    unpinned = base_child_env(tmp_path, {"SKBUILD_EDITABLE_SKIP": "/elsewhere"}, None)
    assert "SKBUILD_EDITABLE_SKIP" not in unpinned


def _this_interpreter_in(root: Path) -> dict[str, str]:
    """This interpreter's child env in *root*, its own build pinned (so no import triggers a rebuild)."""
    return base_child_env(root, os.environ, Path(os.environ["SKBUILD_EDITABLE_SKIP"]))


def test_this_interpreters_library_is_refused_as_the_1_1_tree_of_another_directory(tmp_path: Path) -> None:
    origin = library_origin_of((sys.executable,), tmp_path, _this_interpreter_in(tmp_path))
    with pytest.raises(ForeignLibraryError):
        require_library_in_tree(origin, tmp_path)


def test_this_interpreters_library_lies_in_this_tree() -> None:
    origin = library_origin_of((sys.executable,), ROOT, _this_interpreter_in(ROOT))
    require_library_in_tree(origin, ROOT)
    assert ROOT in origin.compas_cgal_file.parents and ROOT in origin.stock_2_origin.parents


def test_the_editable_build_dir_is_read_without_importing_the_library() -> None:
    build = editable_build_dir_of((sys.executable,), ROOT, base_child_env(ROOT, os.environ, None))
    assert build == Path(os.environ["SKBUILD_EDITABLE_SKIP"])


def test_the_base_tool_digest_is_recomputable_from_the_1_1_commit(tmp_path: Path) -> None:
    for source in BASE_TOOL_SOURCES:
        target = tmp_path / source
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(subprocess.run(["git", "show", f"{EXPECTED_1_1_COMMIT}:{source.as_posix()}"], capture_output=True, check=True).stdout)
    assert base_tool_sha256(tmp_path) == expected_base_tool_sha256(EXPECTED_1_1_COMMIT)


def test_a_base_tree_that_changes_during_the_run_is_refused() -> None:
    trees = iter((TREE, (EXPECTED_1_1_COMMIT, True, "c" * 64)))
    with pytest.raises(TreeChangedDuringRunError):
        measure_1_1_base(
            KEYS,
            excluded=EXCLUDED,
            launch=base_launcher(lambda _command: _result_line(), PYTHON),
            load_reader=lambda: QUIET,
            describe=_describe,
            rebuild=_no_rebuild,
            track_native=_track_native,
            read_tree=lambda: next(trees),
        )


def test_the_1_1_build_is_rebuilt_by_its_own_task_before_the_first_case(tmp_path: Path) -> None:
    events: list[str] = []

    def run(command: Sequence[str]) -> str:
        events.append(" ".join(command))
        return ""

    def launch(key: str) -> CaseRun:
        events.append(key)
        return base_launcher(lambda _command: _result_line(), PYTHON)(key)

    def read_tree() -> TreeIdentity:
        events.append("tree")
        return TREE

    measure_1_1_base(
        KEYS,
        excluded=EXCLUDED,
        launch=launch,
        load_reader=lambda: QUIET,
        describe=_describe,
        rebuild=lambda: rebuild_base_build(tmp_path, run),
        track_native=_track_native,
        read_tree=read_tree,
    )
    assert events[0] == " ".join(("pixi", "run", "--manifest-path", str(tmp_path / "pyproject.toml"), REBUILD_TASK))
    assert events[1] == "tree" and events[2] == KEYS[0]
    assert sum(event.endswith(REBUILD_TASK) for event in events) == 1


def test_the_rebuild_reports_the_command_it_ran(tmp_path: Path) -> None:
    assert rebuild_base_build(tmp_path, lambda _command: "") == ("pixi", "run", "--manifest-path", str(tmp_path / "pyproject.toml"), REBUILD_TASK)


def test_a_failing_rebuild_is_refused_before_any_case_runs(tmp_path: Path) -> None:
    launched: list[str] = []

    def run(command: Sequence[str]) -> str:
        raise CaseProcessFailedError(f"{' '.join(command)} exited 1")

    def launch(key: str) -> CaseRun:
        launched.append(key)
        return base_launcher(lambda _command: _result_line(), PYTHON)(key)

    with pytest.raises(BaseRebuildFailedError):
        measure_1_1_base(
            KEYS,
            excluded=EXCLUDED,
            launch=launch,
            load_reader=lambda: QUIET,
            describe=_describe,
            rebuild=lambda: rebuild_base_build(tmp_path, run),
            track_native=_track_native,
            read_tree=_tree,
        )
    assert launched == []


def test_a_1_1_native_build_that_changes_during_the_base_run_is_refused() -> None:
    """The base runs its own build: its identity is held to itself, start to end, never to the 1.2 build's."""
    with pytest.raises(NativeBuildChangedError):
        measure_1_1_base(
            KEYS,
            excluded=EXCLUDED,
            launch=base_launcher(lambda _command: _result_line(), PYTHON),
            load_reader=lambda: QUIET,
            describe=_describe,
            rebuild=_no_rebuild,
            track_native=_native_build_changed,
            read_tree=_tree,
        )


def test_the_1_1_native_run_starts_after_the_rebuild(tmp_path: Path) -> None:
    events: list[str] = []

    def track() -> NativeFinisher:
        events.append("native start")
        return lambda: NATIVE

    def rebuild() -> None:
        events.append("rebuild")

    measured = measure_1_1_base(
        KEYS,
        excluded=EXCLUDED,
        launch=base_launcher(lambda _command: _result_line(), PYTHON),
        load_reader=lambda: QUIET,
        describe=_describe,
        rebuild=rebuild,
        track_native=track,
        read_tree=_tree,
    )
    assert events == ["rebuild", "native start"]
    assert measured.provenance is not None and measured.provenance.native_sha256 == NATIVE.native_sha256
