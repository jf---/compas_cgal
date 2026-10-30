"""The native libraries a timing run measured: held unchanged from its start to its end, identified in its provenance.

The timing tools time fresh child processes, and every child loads compas_cgal's
extensions from the editable build. Neither the tree identity
(``tools.run_provenance.tree_identity``: ``.so`` files are ignored) nor a hash
read once at the end can tell that a rebuild replaced them during the run. So
the run is bracketed with ``tools.held_native_identity``, the native-identity
mechanism the other timing records use:

- **Start** (`start_native_run`): import every extension the editable finder
  maps, then `capture_current_native` reads their bytes, refuses any written
  since this process started and refuses a stale build.
- **End** (`finish_native_run`): `native_record_fields` reads them again,
  refuses any change since the start (`NativeBuildChangedError`) or any write
  since the process start (`NativeExtensionReplacedAfterStartError`) and
  yields ``native_sha256``.
- **Sources** (`require_built_from`): the build must have been compiled from
  the measured tree's C++ sources. Under the export protocol the ``.so`` comes
  from the worktree's build while the ``.py`` comes from the export, so a
  worktree with C++ edits that the export lacks is refused
  (`NativeSourcesMismatchError`). The record keeps the measured tree's
  `native_sources_sha256`, the native identity the time gate compares, never
  a commit id.

The 1.1 base runs in another tree's interpreter, so its bracket runs there, as
a probe (`probe_main`, launched by `track_native_in`) that imports this tree's
module ahead of the 1.1 tree's. The probe's start is the process start that the
end check holds the files to.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import pathlib
import subprocess
import sys
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from typing import NamedTuple
from typing import Optional

from tools.held_native_identity import EXTENSION_SUFFIXES
from tools.held_native_identity import NATIVE_PACKAGE_PREFIX
from tools.held_native_identity import EditableBuild
from tools.held_native_identity import NativeBuildChangedError
from tools.held_native_identity import NativeRun
from tools.held_native_identity import capture_current_native
from tools.held_native_identity import capture_native
from tools.held_native_identity import native_closure
from tools.held_native_identity import native_fields
from tools.held_native_identity import native_record_fields
from tools.held_native_identity import refuse_replaced_after_start
from tools.held_native_identity import require_current
from tools.run_provenance import file_set_sha256
from tools.run_provenance import not_ignored

ROOT = pathlib.Path(__file__).resolve().parents[1]
# The native sources: every C++ source and header under src/, and the build definition.
NATIVE_SOURCE_DIRECTORY = pathlib.Path("src")
NATIVE_SOURCE_SUFFIXES = (".cpp", ".h", ".hpp")
NATIVE_BUILD_DEFINITION = pathlib.Path("CMakeLists.txt")
# The probe's two modes, run in another tree's interpreter.
PROBE_START = "start"
PROBE_FINISH = "finish"
# Runs `probe_main` from THIS tree in another interpreter: this tree's ``tools`` and
# ``benchmarks`` go first on sys.path; compas_cgal still resolves through that
# interpreter's own editable finder (this tree keeps it under src/, not at the top).
PROBE_SNIPPET = f"import sys; sys.path.insert(0, {str(ROOT)!r}); from tools.run_native_identity import probe_main; probe_main(sys.argv[1:])"


class NativeSourcesMismatchError(Exception):
    """The native build was compiled from other C++ sources than the measured tree's."""


class NativeProbeError(Exception):
    """The native-identity probe in another interpreter failed or printed no identity."""


class NativeIdentity(NamedTuple):
    """What a run's native code was.

    Attributes:
        native_sha256: ``held_native_identity`` digest over every loaded extension (path and bytes).
        native_sources_sha256: `native_sources_sha256` of the measured tree.
    """

    native_sha256: str
    native_sources_sha256: str


# Called at the end of a run: the run's native identity, or a refusal.
NativeFinisher = Callable[[], NativeIdentity]
# Called at the start of a run (after any rebuild): starts the bracket, returns its end.
NativeTracker = Callable[[], NativeFinisher]


def native_source_paths(root: pathlib.Path) -> list[pathlib.Path]:
    """The native sources under *root* (relative), sorted: C++ under src/ and CMakeLists.txt, minus what *root*'s ``.gitignore`` files ignore."""
    sources = [path.relative_to(root) for path in (root / NATIVE_SOURCE_DIRECTORY).rglob("*") if path.is_file() and path.suffix in NATIVE_SOURCE_SUFFIXES]
    return sorted(not_ignored(root, [*sources, NATIVE_BUILD_DEFINITION]), key=lambda path: path.as_posix())


def native_sources_sha256(root: pathlib.Path) -> str:
    """`file_set_sha256` of `native_source_paths`."""
    return file_set_sha256(root, native_source_paths(root))


def require_built_from(build_root: pathlib.Path, measured_root: pathlib.Path) -> None:
    """Refuse a build whose source tree differs from the measured tree in any of the measured tree's native sources.

    Files only the build's tree holds are not compared: one the build compiles
    needs a ``CMakeLists.txt`` change, which is compared.

    Raises:
        NativeSourcesMismatchError: A native source is missing from, or differs in, the build's tree.
    """
    differing = [
        path.as_posix()
        for path in native_source_paths(measured_root)
        if not (build_root / path).is_file() or (build_root / path).read_bytes() != (measured_root / path).read_bytes()
    ]
    if differing:
        raise NativeSourcesMismatchError(
            f"the native build in {build_root} was not compiled from the measured tree {measured_root}: {differing} differ; rebuild from the measured sources"
        )


def require_unchanged_natives(at_start: Mapping[str, str], now: Mapping[str, bytes]) -> None:
    """Refuse when an extension read at the start (path to SHA-256) is missing now or has other bytes.

    Raises:
        NativeBuildChangedError: An extension changed on disk during the run.
    """
    changed = sorted(path for path, sha in at_start.items() if path not in now or hashlib.sha256(now[path]).hexdigest() != sha)
    if changed:
        raise NativeBuildChangedError(f"native extensions changed on disk during the run: {changed}; re-measure once the build settles")


def _import_every_extension(build: EditableBuild) -> None:
    for name, path in sorted(build.wheel_files.items()):
        if name.startswith(NATIVE_PACKAGE_PREFIX) and str(path).endswith(EXTENSION_SUFFIXES):
            importlib.import_module(name)


def start_native_run() -> NativeRun:
    """Import every extension of this interpreter's editable build, then start the ``held_native_identity`` bracket."""
    _import_every_extension(EditableBuild.locate())
    return capture_current_native()


def finish_native_run(run: NativeRun, measured_root: pathlib.Path) -> NativeIdentity:
    """End the bracket started by `start_native_run`; the build must be compiled from *measured_root*'s native sources."""
    native = native_record_fields(run)["native_sha256"]
    require_built_from(run.build.source_directory, measured_root)
    return NativeIdentity(str(native), native_sources_sha256(measured_root))


def track_native_here(measured_root: pathlib.Path) -> NativeTracker:
    """The bracket of this interpreter's native build, for a run measuring *measured_root*."""

    def start() -> NativeFinisher:
        run = start_native_run()
        return lambda: finish_native_run(run, measured_root)

    return start


def _probe(python: Sequence[str], cwd: pathlib.Path, env: Optional[Mapping[str, str]], *arguments: str) -> dict[str, object]:
    completed = subprocess.run([*python, "-c", PROBE_SNIPPET, *arguments], cwd=cwd, env=None if env is None else dict(env), capture_output=True, text=True, check=False)
    lines = completed.stdout.strip().splitlines()
    if completed.returncode != 0 or not lines:
        raise NativeProbeError(f"the native-identity probe ({arguments[0]}) in {cwd} failed ({completed.returncode}): {completed.stderr.strip()[-2000:]}")
    payload = json.loads(lines[-1])
    if not isinstance(payload, dict):
        raise NativeProbeError(f"the native-identity probe printed {lines[-1]!r}, not an object")
    return payload


def track_native_in(python: Sequence[str], root: pathlib.Path, env: Optional[Mapping[str, str]]) -> NativeTracker:
    """The bracket of the native build that the interpreter *python* (run in *root* with *env*) loads, measuring *root*."""

    def start() -> NativeFinisher:
        started = _probe(python, root, env, PROBE_START)

        def finish() -> NativeIdentity:
            identity = _probe(python, root, env, PROBE_FINISH, json.dumps(started), str(root))
            return NativeIdentity(str(identity["native_sha256"]), str(identity["native_sources_sha256"]))

        return finish

    return start


def probe_main(argv: Sequence[str]) -> None:
    """The probe, run in the measured interpreter: ``start`` prints the bracket's start; ``finish START ROOT`` checks it and prints the identity."""
    build = EditableBuild.locate()
    _import_every_extension(build)
    if argv[0] == PROBE_START:
        run = capture_current_native()
        print(json.dumps({"started_ns": run.started_ns, "natives": {path: hashlib.sha256(data).hexdigest() for path, data in run.loaded.items()}}))
        return
    if argv[0] != PROBE_FINISH or len(argv) != 3:
        raise NativeProbeError(f"probe arguments {list(argv)}: expected {PROBE_START!r} or {PROBE_FINISH!r} START ROOT")
    started = json.loads(argv[1])
    at_start: dict[str, str] = dict(started["natives"])
    closure = native_closure(build)
    contents = capture_native(build, tuple(sorted(set(closure) | set(at_start))))
    require_unchanged_natives(at_start, contents)
    refuse_replaced_after_start(build, tuple(contents), int(started["started_ns"]))
    require_current(build, tuple(contents))
    measured_root = pathlib.Path(argv[2])
    require_built_from(build.source_directory, measured_root)
    print(json.dumps({"native_sha256": native_fields(contents)["native_sha256"], "native_sources_sha256": native_sources_sha256(measured_root)}))


if __name__ == "__main__":
    probe_main(sys.argv[1:])
