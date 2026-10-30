"""Time the Task 1.1 tree as the generation-time gate's base, without modifying it (orchestrator ruling on RV27 I1, option b).

The gate (`tools.pipeline_time_gate`) needs the 1.1 tree's serial timing
measured in the same session as the 1.2 timing, with its load and provenance.
The 1.1 tree's own timing tool writes ``pipeline-timing/1``, which records
neither. This wrapper therefore runs at HEAD and drives the UNMODIFIED 1.1
worktree through the tool's own interfaces:

- the case keys and admissions-excluded keys come from the 1.1 tree's
  ``tools.pipeline_timing.gate_timing_keys``, run in its interpreter;
- every run is the 1.1 tool's own child mode (``python -m tools.pipeline_timing
  --case KEY``), in a fresh 1.1 interpreter, so the generator call and its
  timer are the 1.1 code;
- the wrapper reads the load before the first case, before every case (the
  child mode is a per-case hook) and after the last, and refuses a case on a
  loaded machine unless ``--indicative``;
- it lifts the result to ``pipeline-timing/2``, with load and provenance. The
  provenance records the 1.1 tree's commit (it must be `EXPECTED_1_1_COMMIT`, or
  the one named explicitly, and clean), the digest of the 1.1 timing tool
  together with this wrapper, the 1.1 tree's ``pixi.lock``, and the 1.1
  interpreter's environment, Python version and ``_stock_2`` digest.

The measured code is never patched: a modified or wrong tree is refused
before anything runs (`MeasuredTreeError`), and so is an interpreter whose
``compas_cgal`` or ``_stock_2`` does not resolve inside the worktree
(`ForeignLibraryError`): a HEAD interpreter would time the 1.1 tool against the
1.2 library. The children run with the 1.1 tree's own editable build pinned
(``SKBUILD_EDITABLE_SKIP``), exactly as the 1.1 ``pipeline-timing`` task runs
its tool, so no child rebuilds or checks the build at import. Because the pin
drops that task's ``_editable-rebuild`` dependency, the wrapper runs the 1.1
tree's own ``_editable-rebuild`` task once, in the 1.1 tree, before the first
case (`rebuild_base_build`), and prints the command it ran: a stale 1.1 build
would otherwise be timed silently. A failing rebuild refuses the run
(`BaseRebuildFailedError`). The rebuild's environment pins this tree's build,
so even a misresolved interpreter cannot rebuild it.

Run it right before ``pixi run pipeline-timing`` and ``pixi run
pipeline-pencil-cost``, on the same quiet machine::

    python -m tools.pipeline_timing_1_1_base --worktree ../compas_cgal_prs-1.1 \\
        --python "pixi run --no-install python"
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from typing import Mapping
from typing import Optional
from typing import Sequence

from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import MachineNotQuietError
from tools.pipeline_timing import CaseProcessFailedError
from tools.pipeline_timing import CaseRun
from tools.pipeline_timing import Describer
from tools.pipeline_timing import Launcher
from tools.pipeline_timing import LoadReader
from tools.pipeline_timing import MissingCaseResultError
from tools.pipeline_timing import time_case
from tools.pipeline_timing_results import PipelineTimings
from tools.pipeline_timing_results import format_load
from tools.pipeline_timing_results import stored_load
from tools.run_native_identity import NativeIdentity
from tools.run_native_identity import NativeTracker
from tools.run_native_identity import track_native_in
from tools.run_provenance import RunProvenance
from tools.run_provenance import TreeIdentity
from tools.run_provenance import TreeReader
from tools.run_provenance import bytes_set_sha256
from tools.run_provenance import file_set_sha256
from tools.run_provenance import interpreter_identity_of
from tools.run_provenance import machine_identity
from tools.run_provenance import require_unchanged_tree
from tools.run_provenance import sha256_file
from tools.run_provenance import tree_identity
from tools.run_provenance import utc_now

ROOT = Path(__file__).resolve().parents[1]
# The Task 1.1 tree: the last commit touching pipeline_timing.json before any Task 1.2 pipeline source commit.
EXPECTED_1_1_COMMIT = "7b001d473d6e94b79fc8ff3cb852ecde1ed67733"
PIPELINE_TIMING_1_1_SESSION_PATH = Path("benchmarks") / "results" / "pipeline_timing_1_1_session.json"
# The 1.1 tool's result-line prefix (checked against its source at EXPECTED_1_1_COMMIT by the tests).
BASE_RESULT_LINE_PREFIX = "PIPELINE_TIMING_RESULT "
# The 1.1 tool's sources (in the measured tree) and this wrapper's (at HEAD), digested into the provenance.
BASE_TOOL_SOURCES = (Path("tools") / "pipeline_timing.py", Path("tools") / "pipeline_timing_results.py")
WRAPPER_SOURCES = (Path("tools") / "pipeline_timing_1_1_base.py", Path("tools") / "run_provenance.py", Path("tools") / "run_native_identity.py")
PIXI_LOCK = Path("pixi.lock")
# Prints the editable build directory of the interpreter that runs it, as the 1.1 `pipeline-timing` task
# reads it (scikit-build-core's redirecting finder), without importing compas_cgal.
EDITABLE_BUILD_DIR_SNIPPET = "import sys; print(next(f.path for f in sys.meta_path if hasattr(f, 'known_wheel_files') and 'compas_cgal._stock_2' in f.known_wheel_files))"
# Prints {compas_cgal_file, stock_2_origin} for the interpreter that runs it.
LIBRARY_ORIGIN_SNIPPET = (
    "import importlib.util, json; import compas_cgal; spec = importlib.util.find_spec('compas_cgal._stock_2'); "
    "print(json.dumps({'compas_cgal_file': compas_cgal.__file__, 'stock_2_origin': spec.origin if spec else None}))"
)
# Environment variables of this tree's interpreter that must not leak into the 1.1 tree's children.
LEAKING_ENV = ("PYTHONPATH", "SKBUILD_EDITABLE_SKIP")
# Prints the 1.1 tree's (keys to time, excluded keys) as JSON, by its own function.
KEYS_SNIPPET = "import json; from tools.pipeline_timing import gate_timing_keys; print(json.dumps(gate_timing_keys()))"
# The 1.1 tree's own pixi task that brings its editable build up to date (an unpinned import of every
# extension); its `pipeline-timing` task depends on it.
REBUILD_TASK = "_editable-rebuild"
PYPROJECT = Path("pyproject.toml")
INDICATIVE_HINT = "Re-run on a quiet machine, or measure without a verdict: python -m tools.pipeline_timing_1_1_base --indicative ..."

Runner = Callable[[Sequence[str]], str]


class MeasuredTreeError(Exception):
    """The worktree to time is not the expected 1.1 commit, or differs from it: its timing would not be the 1.1 tree's."""


class ForeignLibraryError(Exception):
    """The 1.1 interpreter's ``compas_cgal``, ``_stock_2`` or editable build lies outside the 1.1 worktree: it would time another library."""


class BaseRebuildFailedError(Exception):
    """The 1.1 tree's own rebuild task failed: its build is not known to match its sources, so its timing could be a stale library's."""


@dataclass(frozen=True)
class LibraryOrigin:
    """Where an interpreter's library resolves.

    Attributes:
        compas_cgal_file: ``compas_cgal.__file__``.
        stock_2_origin: ``find_spec("compas_cgal._stock_2").origin``.
    """

    compas_cgal_file: Path
    stock_2_origin: Path


def require_measured_tree(root: Path, expected_rev: str) -> None:
    """Refuse a worktree at *root* that is not at *expected_rev* (full commit id) or is modified.

    Raises:
        MeasuredTreeError: Another commit, or tracked or untracked changes.
    """
    rev, dirty, _patch = tree_identity(root)
    if rev != expected_rev:
        raise MeasuredTreeError(f"{root} is at {rev}, not the 1.1 tree {expected_rev}.")
    if dirty:
        raise MeasuredTreeError(f"{root} differs from {rev}: the 1.1 tree is timed unmodified.")


def base_child_env(root: Path, environ: Mapping[str, str], editable_build_dir: Optional[Path]) -> dict[str, str]:
    """*environ* for a child in the 1.1 tree at *root*: this tree's import path and build pin removed, the 1.1 tree on the path, its build pinned.

    Args:
        root: The 1.1 worktree.
        environ: The environment to start from.
        editable_build_dir: The 1.1 tree's editable build directory to pin (`editable_build_dir_of`), or None to pin none.
    """
    env = {name: value for name, value in environ.items() if name not in LEAKING_ENV}
    env["PYTHONPATH"] = str(root)
    if editable_build_dir is not None:
        env["SKBUILD_EDITABLE_SKIP"] = str(editable_build_dir)
    return env


def editable_build_dir_of(python: Sequence[str], root: Path, env: Mapping[str, str]) -> Path:
    """The editable build directory of the interpreter *python* in *root*, read without importing the library.

    Raises:
        subprocess.CalledProcessError: The interpreter has no editable compas_cgal install.
    """
    completed = subprocess.run([*python, "-c", EDITABLE_BUILD_DIR_SNIPPET], cwd=root, env=dict(env), capture_output=True, text=True, check=True)
    return Path(completed.stdout.strip().splitlines()[-1])


def library_origin_of(python: Sequence[str], root: Path, env: Mapping[str, str]) -> LibraryOrigin:
    """Where ``compas_cgal`` and ``_stock_2`` resolve for the interpreter *python* run in *root* with *env*.

    Raises:
        ForeignLibraryError: ``_stock_2`` does not resolve at all.
    """
    completed = subprocess.run([*python, "-c", LIBRARY_ORIGIN_SNIPPET], cwd=root, env=dict(env), capture_output=True, text=True, check=True)
    payload = json.loads(completed.stdout.strip().splitlines()[-1])
    if payload["stock_2_origin"] is None:
        raise ForeignLibraryError(f"compas_cgal._stock_2 does not resolve for {' '.join(python)} in {root}.")
    return LibraryOrigin(compas_cgal_file=Path(payload["compas_cgal_file"]).resolve(), stock_2_origin=Path(payload["stock_2_origin"]).resolve())


def require_library_in_tree(origin: LibraryOrigin, root: Path, editable_build_dir: Optional[Path] = None) -> None:
    """Refuse a library (and editable build directory) that does not lie inside the worktree *root*.

    Raises:
        ForeignLibraryError: A path resolves outside *root*.
    """
    tree = root.resolve()
    paths = {"compas_cgal": origin.compas_cgal_file, "_stock_2": origin.stock_2_origin}
    if editable_build_dir is not None:
        paths["editable build"] = editable_build_dir.resolve()
    outside = {name: str(path) for name, path in paths.items() if tree not in path.parents}
    if outside:
        raise ForeignLibraryError(f"{outside} lie outside the 1.1 worktree {tree}: the run would time another library.")


def run_in_tree(root: Path, env: Mapping[str, str]) -> Runner:
    """A runner that executes a command in the 1.1 tree at *root* with *env* and returns its stdout.

    Raises (from the runner):
        CaseProcessFailedError: The command exited non-zero; its stderr is quoted.
    """

    def run(command: Sequence[str]) -> str:
        completed = subprocess.run(list(command), cwd=root, env=dict(env), capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise CaseProcessFailedError(f"{' '.join(command)} exited {completed.returncode} in {root}:\n{completed.stderr}")
        return completed.stdout

    return run


def rebuild_base_build(root: Path, run: Runner) -> tuple[str, ...]:
    """Bring the 1.1 tree's editable build at *root* up to date with its own `REBUILD_TASK`, through *run*, and return the command run.

    Raises:
        BaseRebuildFailedError: The task failed (*run* raised `CaseProcessFailedError`).
    """
    command = ("pixi", "run", "--manifest-path", str(root / PYPROJECT), REBUILD_TASK)
    try:
        run(command)
    except CaseProcessFailedError as error:
        raise BaseRebuildFailedError(f"the 1.1 tree's {REBUILD_TASK} failed; its build cannot be trusted to match its sources:\n{error}") from error
    return command


def base_keys(run: Runner, python: Sequence[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """The 1.1 tree's keys to time and admissions-excluded keys, from its own `gate_timing_keys`."""
    to_time, excluded = json.loads(run([*python, "-c", KEYS_SNIPPET]).strip().splitlines()[-1])
    return tuple(str(key) for key in to_time), tuple(str(key) for key in excluded)


def base_launcher(run: Runner, python: Sequence[str]) -> Launcher:
    """A launcher running one case once through the 1.1 tool's child mode, in a fresh 1.1 interpreter."""

    def launch(key: str) -> CaseRun:
        stdout = run([*python, "-m", "tools.pipeline_timing", "--case", key])
        lines = [line for line in stdout.splitlines() if line.startswith(BASE_RESULT_LINE_PREFIX)]
        if len(lines) != 1:
            raise MissingCaseResultError(f"{key}: the 1.1 child printed {len(lines)} result lines, not one.")
        payload = json.loads(lines[0][len(BASE_RESULT_LINE_PREFIX) :])
        return CaseRun(generation_s=float(payload["generation_s"]), operations=int(payload["operations"]), digest=str(payload["digest"]))

    return launch


def measure_1_1_base(
    keys: Sequence[str],
    *,
    excluded: Sequence[str],
    launch: Launcher,
    load_reader: LoadReader,
    describe: Describer,
    rebuild: Callable[[], object],
    read_tree: TreeReader,
    track_native: NativeTracker,
    indicative: bool = False,
) -> PipelineTimings:
    """Time *keys* serially through *launch*, reading the load before each case, and lift the result to ``pipeline-timing/2``.

    *rebuild* runs once, before the tree, time and load are first read and
    before the first case (`rebuild_base_build` in the 1.1 tree). *track_native*
    then starts the bracket of the 1.1 build (after the rebuild, which writes
    it); its end yields the base's own native identity, held to itself from
    start to end and never compared with the 1.2 build's.

    Raises:
        BaseRebuildFailedError: *rebuild* failed; no case ran.
        MachineNotQuietError: Not *indicative*, and the 1-minute load before a case was above `QUIET_LOAD_AVERAGE_MAX`.
        NondeterministicCaseError: A case's runs disagree.
        TreeChangedDuringRunError: The tree's identity changed during the run.
        NativeIdentityError: The 1.1 build changed during the run or is stale (``tools.held_native_identity``).
    """
    rebuild()
    tree_at_start, started_utc, load_start = read_tree(), utc_now(), load_reader()
    finish_native = track_native()
    cases = {}
    for key in keys:
        load = load_reader()
        if not indicative and load[0] > QUIET_LOAD_AVERAGE_MAX:
            raise MachineNotQuietError(f"{key}: load {format_load(load)} (1 / 5 / 15) is above the quiet limit {QUIET_LOAD_AVERAGE_MAX}. {INDICATIVE_HINT}")
        cases[key] = time_case(key, launch)
    load_end, ended_utc, tree_at_end = load_reader(), utc_now(), read_tree()
    require_unchanged_tree(tree_at_start, tree_at_end, "the 1.1 tree")
    native = finish_native()
    return PipelineTimings(
        cases=cases, excluded=tuple(excluded), load_start=stored_load(load_start), load_end=stored_load(load_end), provenance=describe(tree_at_end, native, started_utc, ended_utc)
    )


def base_tool_sha256(root: Path) -> str:
    """The timing-tool digest of a base run: the 1.1 tool's sources in *root* together with this wrapper's sources at HEAD."""
    return _combined_tool_sha256(file_set_sha256(root, BASE_TOOL_SOURCES))


def expected_base_tool_sha256(rev: str) -> str:
    """`base_tool_sha256` as it must read for a run of the unmodified tree at *rev*, from its committed sources (``git show``)."""
    contents = {path: subprocess.run(["git", "show", f"{rev}:{path.as_posix()}"], cwd=ROOT, capture_output=True, check=True).stdout for path in BASE_TOOL_SOURCES}
    return _combined_tool_sha256(bytes_set_sha256(contents))


def _combined_tool_sha256(base_sources_sha256: str) -> str:
    return hashlib.sha256(bytes.fromhex(base_sources_sha256) + bytes.fromhex(file_set_sha256(ROOT, WRAPPER_SOURCES))).hexdigest()


def base_describer(root: Path, python: Sequence[str], env: Mapping[str, str], editable_build_dir: Path) -> Describer:
    """The provenance of a wrapped 1.1 run of the tree at *root* in the interpreter *python* with *env*; re-checks the library's origin."""

    def describe(tree: TreeIdentity, native: NativeIdentity, started_utc: str, ended_utc: str) -> RunProvenance:
        require_library_in_tree(library_origin_of(python, root, env), root, editable_build_dir)
        rev, dirty, patch = tree
        environment, python_version, stock_2 = interpreter_identity_of(python, root, dict(env))
        host, cpu_model = machine_identity()
        return RunProvenance(
            rev=rev,
            dirty=dirty,
            patch_sha256=patch,
            timing_tool_sha256=base_tool_sha256(root),
            pixi_lock_sha256=sha256_file(root / PIXI_LOCK),
            environment=environment,
            python_version=python_version,
            stock_2_sha256=stock_2,
            native_sha256=native.native_sha256,
            native_sources_sha256=native.native_sources_sha256,
            host=host,
            cpu_model=cpu_model,
            started_utc=started_utc,
            ended_utc=ended_utc,
        )

    return describe


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry: time the 1.1 worktree and write its ``pipeline-timing/2`` record.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(prog="python -m tools.pipeline_timing_1_1_base", description="Time the unmodified 1.1 tree as the time gate's base.")
    parser.add_argument("--worktree", type=Path, required=True, help="a clean worktree of the 1.1 tree with its own built environment")
    parser.add_argument("--python", required=True, help="the 1.1 worktree's interpreter command, e.g. 'pixi run --no-install python'")
    parser.add_argument("--expected-rev", default=EXPECTED_1_1_COMMIT, help="the full commit id the worktree must be at")
    parser.add_argument("--indicative", action="store_true", help="measure at any load; the gate will render it as no verdict")
    parser.add_argument("--write", type=Path, default=PIPELINE_TIMING_1_1_SESSION_PATH)
    args = parser.parse_args(argv)
    root = args.worktree.resolve()
    python = tuple(shlex.split(args.python))
    require_measured_tree(root, args.expected_rev)
    editable_build_dir = editable_build_dir_of(python, root, base_child_env(root, os.environ, None))
    env = base_child_env(root, os.environ, editable_build_dir)
    require_library_in_tree(library_origin_of(python, root, env), root, editable_build_dir)
    run = run_in_tree(root, env)
    keys, excluded = base_keys(run, python)
    # The rebuild runs unpinned for the 1.1 build but pins this tree's own build.
    this_build_dir = editable_build_dir_of((sys.executable,), ROOT, os.environ)
    rebuild_run = run_in_tree(root, base_child_env(root, os.environ, this_build_dir))

    def rebuild() -> None:
        print(f"rebuilt the 1.1 build: {' '.join(rebuild_base_build(root, rebuild_run))}")

    measured = measure_1_1_base(
        keys,
        excluded=excluded,
        launch=base_launcher(run, python),
        load_reader=os.getloadavg,
        describe=base_describer(root, python, env, editable_build_dir),
        rebuild=rebuild,
        read_tree=lambda: tree_identity(root),
        track_native=track_native_in(python, root, env),
        indicative=args.indicative,
    )
    require_measured_tree(root, args.expected_rev)
    timings = PipelineTimings.build(measured.to_json())
    write = args.write if args.write.is_absolute() else ROOT / args.write
    write.parent.mkdir(parents=True, exist_ok=True)
    write.write_text(timings.dumps())
    print(f"wrote {write}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
