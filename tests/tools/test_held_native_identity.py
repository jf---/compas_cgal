"""Contract of tools/held_native_identity.py: a record names the native extensions it loaded, and refuses a stale build."""

from __future__ import annotations

from collections.abc import Iterator
import dataclasses
import hashlib
import os
import pathlib
import shutil
import subprocess
import sys
import time
import types
import uuid

import pytest

import compas_cgal._stock_2  # noqa: F401  (a native module the closure must name)
from benchmarks.performance_record import REPOSITORY
import tools.held_connector_repair_budget_results as budget_tool
import tools.held_native_identity as native_identity
import tools.held_splice_ab_results as splice_tool
from tools.held_native_identity import AmbiguousEditableFinderError
from tools.held_native_identity import CMakeCacheMissingError
from tools.held_native_identity import CMakeHomeDirectoryUnsetError
from tools.held_native_identity import EditableBuild
from tools.held_native_identity import EditableBuildDirectoryUnknownError
from tools.held_native_identity import EditableBuildOutsideSourceTreeError
from tools.held_native_identity import MakeProgramUnsetError
from tools.held_native_identity import NativeBuildChangedError
from tools.held_native_identity import NativeExtensionReplacedAfterStartError
from tools.held_native_identity import NinjaMissingError
from tools.held_native_identity import NoEditableFinderError
from tools.held_native_identity import NotNinjaBuildError
from tools.held_native_identity import NativeBuildStaleError
from tools.held_native_identity import NativeIdentityMissingError
from tools.held_native_identity import NativeModuleOutsideEditableBuildError
from tools.held_native_identity import NativeRun
from tools.held_native_identity import ProcessStartUnavailableError
from tools.held_native_identity import capture_current_native
from tools.held_native_identity import capture_native
from tools.held_native_identity import native_closure
from tools.held_native_identity import native_fields
from tools.held_native_identity import native_identity_line
from tools.held_native_identity import native_record_fields
from tools.held_native_identity import process_start_ns
from tools.held_native_identity import refuse_replaced_after_start
from tools.held_native_identity import require_current
from tools.held_native_identity import stale_reasons
from tools.held_source_closure import digest_contents

EXTENSION = "_probe.abi3.so"
MODULE = "compas_cgal._probe"
# The synthetic build's directory inside its source tree, and its install copy (scikit-build-core's layout).
BUILD_DIRECTORY = pathlib.Path("build") / "cp-synthetic"
INSTALLED = pathlib.Path("install") / "platlib" / "compas_cgal" / EXTENSION
BUILD_NINJA = """rule configure
  command = touch build.ninja
  generator = 1
build build.ninja: configure CMakeLists.txt
rule copy
  command = cp $in $out
build _probe.abi3.so: copy probe.cpp
"""
# How far (ns) a stale input's mtime is set past its output's: one second, far above any file-system timestamp granularity.
STALE_STEP_NS = 1_000_000_000


def _synthetic_build(source: pathlib.Path, *, newer: str | None) -> tuple[EditableBuild, str]:
    """A one-target Ninja build in ``source``/BUILD_DIRECTORY with a CMake-like regeneration rule, really built once and installed.

    ``source`` is the build's source tree (its ``CMAKE_HOME_DIRECTORY``); the
    returned path is relative to it. With ``newer`` naming ``probe.cpp`` (the
    target's source) or ``CMakeLists.txt`` (the manifest's input), that file is
    then made newer than what was built from it.
    """
    root = source / BUILD_DIRECTORY
    build = EditableBuild(directory=root, source_directory=source, wheel_files={MODULE: root / INSTALLED}, make_program=EditableBuild.locate().make_program)
    root.mkdir(parents=True)
    (root / "CMakeLists.txt").write_bytes(b"configure")
    (root / "build.ninja").write_text(BUILD_NINJA, encoding="utf-8")
    manifest = (root / "build.ninja").stat().st_mtime_ns
    os.utime(root / "CMakeLists.txt", ns=(manifest - STALE_STEP_NS, manifest - STALE_STEP_NS))
    (root / "probe.cpp").write_bytes(b"binary")
    subprocess.run([str(build.make_program), "-C", str(root)], capture_output=True, check=True)
    installed = root / INSTALLED
    installed.parent.mkdir(parents=True)
    shutil.copyfile(root / EXTENSION, installed)
    if newer is not None:
        output = root / ("build.ninja" if newer == "CMakeLists.txt" else EXTENSION)
        later = output.stat().st_mtime_ns + STALE_STEP_NS
        os.utime(root / newer, ns=(later, later))
    return build, installed.relative_to(source).as_posix()


def _load_only_the_probe(monkeypatch: pytest.MonkeyPatch, build: EditableBuild) -> None:
    """The module sees one loaded compas_cgal extension: the synthetic probe, from the file the build maps it to."""
    probe = types.ModuleType(MODULE)
    probe.__file__ = str(build.wheel_files[MODULE])
    monkeypatch.setattr(native_identity, "sys", types.SimpleNamespace(modules={MODULE: probe}, meta_path=sys.meta_path))


@pytest.fixture
def scratch() -> Iterator[pathlib.Path]:
    # A synthetic build's source tree, under the repository's gitignored build/ (never /tmp).
    root = REPOSITORY / "build" / "native-identity-tests" / uuid.uuid4().hex
    yield root
    shutil.rmtree(root)


def test_the_native_closure_names_every_loaded_extension_relative_to_the_builds_source_tree_and_sorted() -> None:
    build = EditableBuild.locate()
    assert build.directory.is_relative_to(build.source_directory)
    closure = native_closure(build)
    assert list(closure) == sorted(set(closure))
    stock = build.wheel_files["compas_cgal._stock_2"].resolve().relative_to(build.source_directory).as_posix()
    assert stock in closure
    contents = capture_native(build, closure)
    assert hashlib.sha256(contents[stock]).hexdigest() == hashlib.sha256(pathlib.Path(sys.modules["compas_cgal._stock_2"].__file__).read_bytes()).hexdigest()


def test_the_native_digest_changes_with_any_extension_byte() -> None:
    build = EditableBuild.locate()
    contents = capture_native(build, native_closure(build))
    path = next(iter(contents))
    changed = {**contents, path: contents[path] + b"\0"}
    assert digest_contents(changed) != digest_contents(contents)


@pytest.mark.parametrize(
    "filename",
    [pathlib.Path("/opt/elsewhere/compas_cgal/_stray.abi3.so"), REPOSITORY / "build" / "cp39-old" / "compas_cgal" / "_stray.abi3.so"],
    ids=["outside the repository", "an older build directory in the repository"],
)
def test_a_loaded_extension_that_is_not_the_editable_builds_file_is_refused(monkeypatch: pytest.MonkeyPatch, filename: pathlib.Path) -> None:
    stray = types.ModuleType("compas_cgal._stray")
    stray.__file__ = str(filename)
    monkeypatch.setitem(sys.modules, "compas_cgal._stray", stray)
    with pytest.raises(NativeModuleOutsideEditableBuildError, match="_stray.*but the editable finder maps it to None"):
        native_closure(EditableBuild.locate())


def test_an_extension_the_finder_maps_outside_the_build_directory_is_refused(monkeypatch: pytest.MonkeyPatch, scratch: pathlib.Path) -> None:
    # The finder's own map is not trusted to stay inside the build: a .so outside the build directory is not that build's output.
    build, _ = _synthetic_build(scratch, newer=None)
    outside = scratch / "elsewhere" / EXTENSION
    outside.parent.mkdir()
    outside.write_bytes(b"binary")
    build = dataclasses.replace(build, wheel_files={MODULE: outside})
    _load_only_the_probe(monkeypatch, build)
    with pytest.raises(NativeModuleOutsideEditableBuildError, match=f"{MODULE} is loaded from {outside}, outside the editable build directory"):
        native_closure(build)


def test_the_native_identity_is_rooted_at_the_builds_source_tree_not_the_python_sources(monkeypatch: pytest.MonkeyPatch, scratch: pathlib.Path) -> None:
    # The export layout (build/commit-check/cp13/run.sh): the .py come from an export tree (REPOSITORY), the .so from an
    # editable build configured from another tree. The whole run guard must pass, with paths relative to the build's tree.
    build, path = _synthetic_build(scratch, newer=None)
    assert build.source_directory != REPOSITORY
    _load_only_the_probe(monkeypatch, build)
    monkeypatch.setattr(EditableBuild, "locate", classmethod(lambda cls: build))
    monkeypatch.setattr(native_identity, "process_start_ns", time.time_ns)
    run = capture_current_native()
    fields = native_record_fields(run)
    assert path == (BUILD_DIRECTORY / INSTALLED).as_posix()
    assert list(fields["natives"]) == [path]
    assert fields["natives"][path] == hashlib.sha256(b"binary").hexdigest()


def test_a_target_newer_than_its_sources_is_current(scratch: pathlib.Path) -> None:
    build, path = _synthetic_build(scratch, newer=None)
    log = (build.directory / ".ninja_log").read_bytes()
    assert stale_reasons(build, (path,)) == ()
    require_current(build, (path,))
    # A dry run: Ninja's log is untouched.
    assert (build.directory / ".ninja_log").read_bytes() == log


def test_a_source_newer_than_the_target_is_stale(scratch: pathlib.Path) -> None:
    build, path = _synthetic_build(scratch, newer="probe.cpp")
    assert stale_reasons(build, (path,))
    with pytest.raises(NativeBuildStaleError, match="probe"):
        require_current(build, (path,))


def test_an_install_copy_that_differs_from_the_build_tree_is_stale(scratch: pathlib.Path) -> None:
    build, path = _synthetic_build(scratch, newer=None)
    (scratch / path).write_bytes(b"older binary")
    with pytest.raises(NativeBuildStaleError, match="differs from the build tree"):
        require_current(build, (path,))


def test_an_installed_extension_no_target_builds_is_stale(scratch: pathlib.Path) -> None:
    build, path = _synthetic_build(scratch, newer=None)
    leftover = scratch / path
    orphan = leftover.with_name("_leftover.abi3.so")
    orphan.write_bytes(b"binary")
    with pytest.raises(NativeBuildStaleError, match="no build target"):
        require_current(build, (orphan.relative_to(scratch).as_posix(),))


def test_a_rebuild_during_the_run_is_refused(monkeypatch: pytest.MonkeyPatch, scratch: pathlib.Path) -> None:
    build, path = _synthetic_build(scratch, newer=None)
    _load_only_the_probe(monkeypatch, build)
    run = NativeRun(build=build, started_ns=time.time_ns(), loaded=capture_native(build, (path,)))
    (scratch / path).write_bytes(b"rebuilt")
    with pytest.raises(NativeBuildChangedError):
        native_record_fields(run)


def test_the_record_fields_name_each_extension_and_the_digest() -> None:
    build = EditableBuild.locate()
    before = capture_native(build, native_closure(build))
    fields = native_fields(before)
    assert list(fields["natives"]) == sorted(before)
    assert fields["native_sha256"] == digest_contents(before)
    assert all(fields["natives"][path] == hashlib.sha256(before[path]).hexdigest() for path in before)
    assert fields["native_sha256"][:12] in native_identity_line(fields)


def test_a_record_without_the_native_identity_does_not_render() -> None:
    with pytest.raises(NativeIdentityMissingError):
        native_identity_line({"source_sha256": "f" * 64})


def test_a_pending_regeneration_of_the_build_manifest_is_stale(scratch: pathlib.Path) -> None:
    # Ninja's dry run returns early at the regeneration ("Re-running CMake..."), before it looks at any target.
    build, path = _synthetic_build(scratch, newer="CMakeLists.txt")
    manifest = (build.directory / "build.ninja").stat().st_mtime_ns
    with pytest.raises(NativeBuildStaleError, match="build.ninja"):
        require_current(build, (path,))
    assert (build.directory / "build.ninja").stat().st_mtime_ns == manifest


def test_the_process_start_precedes_now() -> None:
    assert 0 < process_start_ns() <= time.time_ns()


# A child that reads its own process start, then (after a stdin line) checks one file against it with the module's refusal.
_START_PROBE = "import tools.held_native_identity as m; print(m.process_start_ns())"
_REFUSAL_PROBE = """
import pathlib, sys, types
import tools.held_native_identity as m
started = m.process_start_ns()
print(started, flush=True)
sys.stdin.readline()
m.refuse_replaced_after_start(types.SimpleNamespace(source_directory=pathlib.Path(sys.argv[1])), (sys.argv[2],), started)
"""
# Spawns per range check: a whole-second start passes only when a spawn lands on a second boundary, so all five must hold.
START_SPAWNS = 5
# The second fraction (ns) before which a chmod+spawn starts, so both land in one wall-clock second (a spawn takes ms).
SAME_SECOND_BEFORE_NS = 500_000_000


def _child_start_ns() -> tuple[int, int, int]:
    before = time.time_ns()
    printed = subprocess.run([sys.executable, "-c", _START_PROBE], cwd=REPOSITORY, capture_output=True, text=True, check=True).stdout
    after = time.time_ns()
    return before, int(printed), after


def test_a_child_process_start_lies_between_the_clock_readings_around_its_spawn() -> None:
    for _ in range(START_SPAWNS):
        before, started, after = _child_start_ns()
        if sys.platform == "darwin":
            # The kernel keeps the start in microseconds: truncated to one, it cannot precede the microsecond of `before`.
            assert before // native_identity.NS_PER_US * native_identity.NS_PER_US <= started <= after
        else:
            # The ps path: the start floored to the whole second.
            assert started % native_identity.NS_PER_S == 0
            assert started <= after


def _refusal_child(extension: pathlib.Path) -> subprocess.Popen[str]:
    return subprocess.Popen(
        [sys.executable, "-c", _REFUSAL_PROBE, str(extension.parent), extension.name],
        cwd=REPOSITORY,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _extension_file(scratch: pathlib.Path) -> pathlib.Path:
    scratch.mkdir(parents=True)
    extension = scratch / EXTENSION
    extension.write_bytes(b"loaded image")
    return extension


def test_an_extension_written_before_a_process_started_in_the_same_second_is_accepted_where_the_start_is_sub_second(scratch: pathlib.Path) -> None:
    # Session aa-2: the no-op editable rebuild bumped the .so's ctime in the second the probe started, strictly before it.
    extension = _extension_file(scratch)
    while time.time_ns() % native_identity.NS_PER_S >= SAME_SECOND_BEFORE_NS:
        pass
    os.chmod(extension, 0o600)
    child = _refusal_child(extension)
    out, err = child.communicate("\n")
    assert extension.stat().st_ctime_ns // native_identity.NS_PER_S == int(out.split()[0]) // native_identity.NS_PER_S, "not one second"
    if sys.platform == "darwin":
        assert child.returncode == 0, err
    else:
        # The ps path floors the start: a same-second write is refused by contract (the conservative direction).
        assert "NativeExtensionReplacedAfterStartError" in err


def test_an_extension_written_after_a_process_started_is_refused_by_that_process(scratch: pathlib.Path) -> None:
    # Defect injection for the sub-second start: the same child must refuse a write made after it read its start.
    extension = _extension_file(scratch)
    child = _refusal_child(extension)
    assert child.stdout is not None
    child.stdout.readline()
    os.chmod(extension, 0o600)
    _, err = child.communicate("\n")
    assert child.returncode != 0
    assert "NativeExtensionReplacedAfterStartError" in err


def test_an_extension_written_after_the_process_started_is_refused(scratch: pathlib.Path) -> None:
    # The loaded image may then differ from the bytes on disk: they must not be recorded as what ran.
    build, path = _synthetic_build(scratch, newer=None)
    written = (scratch / path).stat()
    latest = max(written.st_ctime_ns, written.st_mtime_ns)
    with pytest.raises(NativeExtensionReplacedAfterStartError, match="probe"):
        refuse_replaced_after_start(build, (path,), min(written.st_ctime_ns, written.st_mtime_ns) - 1)
    # A tie is refused: on a whole-second file system a write inside the (floored) start second carries the start second.
    with pytest.raises(NativeExtensionReplacedAfterStartError, match="probe"):
        refuse_replaced_after_start(build, (path,), latest)
    refuse_replaced_after_start(build, (path,), latest + 1)


def test_the_run_start_capture_refuses_an_extension_newer_than_the_process(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(native_identity, "process_start_ns", lambda: 0)
    with pytest.raises(NativeExtensionReplacedAfterStartError):
        capture_current_native()


class _Sentinel(Exception):
    """Raised by the patched native capture: the tool reached it."""


class _WorkloadStarted(Exception):
    """Raised by the patched first workload: the tool measured before its native check."""


def _refuse(*_: object, **__: object) -> None:
    raise _Sentinel


def _started(*_: object, **__: object) -> None:
    raise _WorkloadStarted


def test_the_splice_tool_checks_the_native_build_before_any_workload(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(splice_tool, "capture_current_native", _refuse)
    monkeypatch.setattr(splice_tool, "load_held_reference_case", _started)
    with pytest.raises(_Sentinel):
        splice_tool.measure()


def test_the_budget_tool_checks_the_native_build_before_any_workload(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(budget_tool, "capture_current_native", _refuse)
    monkeypatch.setattr(budget_tool, "measure_budget", _started)
    with pytest.raises(_Sentinel):
        budget_tool.measure(["figure5_refined"])


class _Finder:
    """A stand-in for scikit-build-core's editable finder."""

    def __init__(self, path: object) -> None:
        self.path = path
        self.known_wheel_files = {native_identity.EDITABLE_PROBE_MODULE: "x.abi3.so"}


def _cache(directory: pathlib.Path, text: str) -> pathlib.Path:
    directory.mkdir(parents=True)
    (directory / "CMakeCache.txt").write_text(text, encoding="utf-8")
    return directory


NINJA_CACHE = "CMAKE_GENERATOR:INTERNAL=Ninja\n"


def _configured_cache(directory: pathlib.Path, home: pathlib.Path | None) -> pathlib.Path:
    """A Ninja cache with a make program on disk and, when given, ``home`` as the source tree it was configured from."""
    make_program = f"CMAKE_MAKE_PROGRAM:FILEPATH={EditableBuild.locate().make_program}\n"
    return _cache(directory, NINJA_CACHE + make_program + ("" if home is None else f"CMAKE_HOME_DIRECTORY:INTERNAL={home}\n"))


def test_no_editable_finder_is_refused() -> None:
    with pytest.raises(NoEditableFinderError):
        EditableBuild.from_meta_path([object()])


def test_two_editable_finders_are_refused() -> None:
    with pytest.raises(AmbiguousEditableFinderError):
        EditableBuild.from_meta_path([_Finder("a"), _Finder("b")])


def test_a_finder_without_a_build_directory_is_refused() -> None:
    with pytest.raises(EditableBuildDirectoryUnknownError):
        EditableBuild.from_meta_path([_Finder(None)])


def test_a_build_without_a_cmake_cache_is_refused(scratch: pathlib.Path) -> None:
    scratch.mkdir(parents=True)
    with pytest.raises(CMakeCacheMissingError):
        EditableBuild.from_meta_path([_Finder(str(scratch))])


def test_a_non_ninja_build_is_refused(scratch: pathlib.Path) -> None:
    _cache(scratch, "CMAKE_GENERATOR:INTERNAL=Unix Makefiles\nCMAKE_MAKE_PROGRAM:FILEPATH=/usr/bin/make\n")
    with pytest.raises(NotNinjaBuildError):
        EditableBuild.from_meta_path([_Finder(str(scratch))])


def test_a_cache_without_a_make_program_is_refused(scratch: pathlib.Path) -> None:
    _cache(scratch, NINJA_CACHE)
    with pytest.raises(MakeProgramUnsetError):
        EditableBuild.from_meta_path([_Finder(str(scratch))])


def test_a_missing_ninja_is_refused(scratch: pathlib.Path) -> None:
    _cache(scratch, NINJA_CACHE + f"CMAKE_MAKE_PROGRAM:FILEPATH={scratch / 'no-ninja'}\n")
    with pytest.raises(NinjaMissingError):
        EditableBuild.from_meta_path([_Finder(str(scratch))])


def test_a_cache_without_a_source_tree_is_refused(scratch: pathlib.Path) -> None:
    _configured_cache(scratch / BUILD_DIRECTORY, None)
    with pytest.raises(CMakeHomeDirectoryUnsetError):
        EditableBuild.from_meta_path([_Finder(str(scratch / BUILD_DIRECTORY))])


def test_a_build_outside_its_source_tree_is_refused(scratch: pathlib.Path) -> None:
    # Native paths are recorded relative to the source tree, so the build (and its install copies) must lie inside it.
    _configured_cache(scratch / BUILD_DIRECTORY, scratch / "elsewhere")
    with pytest.raises(EditableBuildOutsideSourceTreeError, match="elsewhere"):
        EditableBuild.from_meta_path([_Finder(str(scratch / BUILD_DIRECTORY))])


def test_the_source_tree_is_the_one_the_build_was_configured_from(scratch: pathlib.Path) -> None:
    _configured_cache(scratch / BUILD_DIRECTORY, scratch)
    build = EditableBuild.from_meta_path([_Finder(str(scratch / BUILD_DIRECTORY))])
    assert (build.directory, build.source_directory) == ((scratch / BUILD_DIRECTORY).resolve(), scratch.resolve())


def test_a_graph_target_that_was_never_built_is_stale(scratch: pathlib.Path) -> None:
    build, path = _synthetic_build(scratch, newer=None)
    (build.directory / EXTENSION).unlink()
    assert any("never built" in reason for reason in stale_reasons(build, (path,)))


def test_an_unparseable_ps_process_start_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(native_identity.subprocess, "run", lambda *_, **__: subprocess.CompletedProcess([], 0, stdout="not a date\n", stderr=""))
    with pytest.raises(ProcessStartUnavailableError, match="not a date"):
        native_identity.ps_process_start_ns()


def _order_probe(monkeypatch: pytest.MonkeyPatch, calls: list[str]) -> None:
    """Record the order of byte reads, start-time checks and Ninja checks in the module."""
    read, times = native_identity.capture_native, native_identity.refuse_replaced_after_start

    def capture(build: EditableBuild, paths: object) -> dict[str, bytes]:
        calls.append("bytes")
        return read(build, paths)

    def check_times(build: EditableBuild, paths: object, started_ns: int) -> None:
        calls.append("times")
        times(build, paths, started_ns)

    def check_current(build: object, paths: object) -> None:
        # Recorded, not run: the shared build may be stale while other lanes edit; the verdict is tested on synthetic builds.
        calls.append("ninja")

    monkeypatch.setattr(native_identity, "capture_native", capture)
    monkeypatch.setattr(native_identity, "refuse_replaced_after_start", check_times)
    monkeypatch.setattr(native_identity, "require_current", check_current)
    monkeypatch.setattr(native_identity, "process_start_ns", lambda: time.time_ns())


def test_the_run_start_reads_the_bytes_before_it_checks_their_times(monkeypatch: pytest.MonkeyPatch) -> None:
    # A replacement before the read shows a later ctime; one after it leaves the held bytes those of the loaded image.
    calls: list[str] = []
    _order_probe(monkeypatch, calls)
    run = capture_current_native()
    assert calls == ["bytes", "times", "ninja"]
    assert set(run.loaded) == set(native_closure(run.build))


def test_the_run_end_reads_the_bytes_once_and_records_those(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    _order_probe(monkeypatch, calls)
    run = capture_current_native()
    calls.clear()
    fields = native_record_fields(run)
    assert calls == ["bytes", "times", "ninja"]
    assert fields == native_fields(run.loaded)


def test_an_extension_loaded_during_the_run_is_checked_against_the_process_start() -> None:
    # Nothing loaded at the start: every extension in the end-of-run closure was loaded lazily, after a start at the epoch.
    run = NativeRun(build=EditableBuild.locate(), started_ns=0, loaded={})
    with pytest.raises(NativeExtensionReplacedAfterStartError):
        native_record_fields(run)


def test_the_ps_process_start_is_the_whole_second_at_or_before_the_kernel_start() -> None:
    # Both paths read one process: on macOS the kernel start floored to the second is what ps prints.
    started = native_identity.ps_process_start_ns()
    assert started % native_identity.NS_PER_S == 0
    assert started <= time.time_ns()
    if sys.platform == "darwin":
        assert started == native_identity.kernel_process_start_ns() // native_identity.NS_PER_S * native_identity.NS_PER_S
