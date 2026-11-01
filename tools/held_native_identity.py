"""Identity of the native extensions a measuring process loaded, and the refusal of a stale editable build.

A record's commit and dirty flag name the C++ sources, but the counts are
decided by the compiled extensions the process actually loaded. Through
``pixi run <task>`` the editable install rebuilds them first; through
``build/pyrun.sh`` (``SKBUILD_EDITABLE_SKIP``) or a direct ``measure()`` call
it does not, so a stale ``.so`` would be invisible in the record.

Two parts:

- **Identity.** Every loaded ``compas_cgal.*`` extension module must be the
  file the scikit-build-core editable finder maps it to, and must lie inside
  that build's directory (else a stray install is measured). Its path
  relative to the build's own source tree (the ``CMAKE_HOME_DIRECTORY`` it
  was configured from) and its bytes are hashed with
  ``held_source_closure.digest_contents`` into ``native_sha256``; the
  per-file SHA-256s go in ``natives``. The root is the build's tree, not the
  one the Python sources come from: a ``.so``'s provenance is its build, a
  ``.py``'s is the source closure digest, and under the export protocol
  (``build/commit-check/cp13/run.sh``: ``.py`` from an export, ``.so`` from
  the worktree's editable build) the two trees differ by design.
- **Currency.** scikit-build-core records no source hash, and the build
  directory holds no stamp of the repository's C++ sources (its
  ``native-source-digests`` cover the vendored CGAL/Boost/Eigen only), so a
  content comparison against ``git ls-files src`` has nothing to compare
  with. The sound check available is the build system's own: a Ninja dry
  run (``-n``, which writes neither ``.ninja_log`` nor ``.ninja_deps``) of
  each loaded extension's build target, over the dependency graph CMake
  generated from the build's source tree and the compiler's header depfiles. It must report no work,
  and the installed copy must equal the build tree's target byte for byte.
  This is exactly the criterion the editable rebuild of the pixi tasks
  applies, and it shares that criterion's blind spots, each a false
  *fresh* (a stale ``.so`` accepted as current): a source edit carrying an
  mtime older than its target (an mtime-preserving copy such as ``cp -p``,
  ``rsync -a``, ``tar x`` or a backup restore; a clock step back), and a
  change of ``[tool.scikit-build] cmake.define`` or the build type in
  ``pyproject.toml``, which does not invalidate ``build.ninja``. In those
  cases the record still names the bytes that ran (the identity stays
  true; only the currency verdict is wrong). A regeneration of
  ``build.ninja`` pending (``CMakeLists.txt`` newer than the last
  configure) stops the dry run at "Re-running CMake" and is reported as
  stale, the conservative direction.
- **Loaded image.** Hashing reads the file on disk, but the process runs the
  image it mapped at import. So the bytes are read FIRST, then every
  extension whose ctime or mtime is at or after the process start is
  refused (it may have been replaced after it was loaded): a replacement
  before the read shows a later time, and one after the read leaves the
  held bytes those of the pre-start file, the loaded image. At the end of
  the run the closure (now with the lazily loaded extensions) is read once,
  compared with the start-of-run bytes, checked against the same start, and
  the record is built from exactly those bytes. A tie is refused. On macOS
  the start is the kernel's own (``sysctl`` ``KERN_PROC_PID``, the
  ``p_starttime`` timeval), kept in microseconds and truncated to one, so a
  file written in the microsecond before the start is refused (the
  conservative direction) and one written earlier in the start second is
  accepted. Elsewhere the start is read from ``ps`` in whole seconds,
  rounded down, so a file written in the start second is refused, on a
  whole-second file system too. Conservative refusals (re-run): a write in
  the second (``ps``) or microsecond (kernel) before the start, and a
  ctime-only change (chmod, an xattr, a hard link). False accept: on Linux
  ``ps`` derives the start as wall-clock now minus the monotonic elapsed
  time, so a forward clock step between the start and the check reads the
  start late. ``strptime``'s ``%a``/``%b`` follow this process's
  ``LC_TIME``, which matches ``ps``'s ``LC_ALL=C`` output unless something
  has called ``setlocale``; a mismatch raises, it never misreads.
"""

from __future__ import annotations

import ctypes
import hashlib
import importlib.machinery
import os
import pathlib
import subprocess
import sys
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from datetime import timezone
from typing import Any
from typing import cast

from benchmarks.errors import BenchmarkError
from tools.held_source_closure import digest_contents

# Loaded modules whose native code can decide a measured number.
NATIVE_PACKAGE_PREFIX = "compas_cgal."
# The editable install's generated finder: the module whose meta-path finder knows this extension is the one in use.
EDITABLE_PROBE_MODULE = "compas_cgal._stock_2"
EXTENSION_SUFFIXES = tuple(importlib.machinery.EXTENSION_SUFFIXES)
CMAKE_CACHE = "CMakeCache.txt"
CMAKE_GENERATOR_KEY = "CMAKE_GENERATOR:INTERNAL"
CMAKE_MAKE_PROGRAM_KEY = "CMAKE_MAKE_PROGRAM"
# The top-level source tree the build was configured from (CMAKE_SOURCE_DIR of that configure; build.ninja's paths point into it).
CMAKE_HOME_DIRECTORY_KEY = "CMAKE_HOME_DIRECTORY"
NINJA_GENERATOR = "Ninja"
# Ninja's dry-run verdict when every requested target is newer than all its inputs.
NINJA_UP_TO_DATE = "ninja: no work to do."
NINJA_ENTERING = "ninja: Entering directory"
# Process start time from ps, printed in UTC under the C locale so that it parses without DST or locale ambiguity.
PS_START_COMMAND = ("ps", "-o", "lstart=", "-p")
PS_ENVIRONMENT = {"TZ": "UTC0", "LC_ALL": "C"}
PS_START_FORMAT = "%a %b %d %H:%M:%S %Y"
NS_PER_S = 1_000_000_000
NS_PER_US = 1_000
# macOS <sys/sysctl.h>: the MIB {CTL_KERN, KERN_PROC, KERN_PROC_PID, pid} returns that process's struct kinfo_proc.
CTL_KERN = 1
KERN_PROC = 14
KERN_PROC_PID = 1


class NativeIdentityError(BenchmarkError):
    """Base of the native-identity failure modes."""


class NativeModuleOutsideEditableBuildError(NativeIdentityError):
    """A loaded compas_cgal extension is not the file the editable build maps it to (a stray or stale install)."""


class NativeBuildStaleError(NativeIdentityError):
    """The editable rebuild would change a loaded extension: a source is newer, or the install copy is not the build's."""


class NoEditableFinderError(NativeIdentityError):
    """No scikit-build-core editable finder maps compas_cgal's extensions: compas_cgal is not an editable install here."""


class AmbiguousEditableFinderError(NativeIdentityError):
    """Several editable finders map compas_cgal's extensions: which build is loaded is undecided."""


class EditableBuildDirectoryUnknownError(NativeIdentityError):
    """The editable finder names no build directory."""


class CMakeCacheMissingError(NativeIdentityError):
    """The editable build directory holds no CMakeCache.txt: it was never configured, or it is not the build."""


class CMakeHomeDirectoryUnsetError(NativeIdentityError):
    """The CMake cache names no CMAKE_HOME_DIRECTORY: the source tree the build was configured from is unknown."""


class EditableBuildOutsideSourceTreeError(NativeIdentityError):
    """The editable build directory is not inside the source tree it was configured from; native paths could not be rooted there."""


class NotNinjaBuildError(NativeIdentityError):
    """The editable build was not generated for Ninja; only Ninja's dependency graph is checked."""


class MakeProgramUnsetError(NativeIdentityError):
    """The CMake cache names no CMAKE_MAKE_PROGRAM."""


class NinjaMissingError(NativeIdentityError):
    """The Ninja the build was configured with is not on disk."""


class NativeBuildToolError(NativeIdentityError):
    """Ninja failed on the build's own graph (a missing input with no rule, a malformed manifest)."""


class ProcessStartUnavailableError(NativeIdentityError):
    """The kernel (sysctl, on macOS) or ps (elsewhere) did not report this process's start time."""


class NativeExtensionReplacedAfterStartError(NativeIdentityError):
    """A loaded extension's file was written after this process started: the image that runs may not be these bytes."""


class NativeBuildChangedError(NativeIdentityError):
    """An extension file loaded at the start changed on disk during the run (a rebuild mid-measurement)."""


class NativeIdentityMissingError(NativeIdentityError):
    """A record carries no native extension identity; it predates the field and must be re-measured."""


@dataclass(frozen=True)
class EditableBuild:
    """The scikit-build-core editable build: its directory, the source tree it was configured from, its module-to-installed-file map, its Ninja.

    The native identity is rooted at ``source_directory`` (the build's own
    ``CMAKE_HOME_DIRECTORY``), never at the tree the Python sources come from:
    under the export protocol (``build/commit-check/cp13/run.sh``) the ``.py``
    are an export's and the ``.so`` the worktree build's, two identities.
    """

    directory: pathlib.Path
    source_directory: pathlib.Path
    wheel_files: Mapping[str, pathlib.Path]
    make_program: pathlib.Path

    @classmethod
    def locate(cls) -> EditableBuild:
        """The editable build this interpreter imports compas_cgal from (the finder ``build/pyrun.sh`` pins)."""
        return cls.from_meta_path(sys.meta_path)

    @classmethod
    def from_meta_path(cls, meta_path: Sequence[object]) -> EditableBuild:
        """The editable build of the one finder in ``meta_path`` that maps compas_cgal's extensions."""
        finders = [finder for finder in meta_path if EDITABLE_PROBE_MODULE in getattr(finder, "known_wheel_files", {})]
        if not finders:
            raise NoEditableFinderError(f"no editable finder maps {EDITABLE_PROBE_MODULE}; run through pixi (build/pyrun.sh or a task)")
        if len(finders) > 1:
            raise AmbiguousEditableFinderError(f"{len(finders)} editable finders map {EDITABLE_PROBE_MODULE}")
        finder = finders[0]
        path = getattr(finder, "path", None)
        if path is None:
            raise EditableBuildDirectoryUnknownError(f"the editable finder {finder!r} names no build directory")
        directory = pathlib.Path(path).resolve()
        cache = _cmake_cache(directory)
        if cache.get(CMAKE_GENERATOR_KEY) != NINJA_GENERATOR:
            raise NotNinjaBuildError(f"{directory} is not a Ninja build ({cache.get(CMAKE_GENERATOR_KEY)!r}); only Ninja's graph is checked")
        if CMAKE_MAKE_PROGRAM_KEY not in cache:
            raise MakeProgramUnsetError(f"{directory / CMAKE_CACHE} names no {CMAKE_MAKE_PROGRAM_KEY}")
        make_program = pathlib.Path(cache[CMAKE_MAKE_PROGRAM_KEY])
        if not make_program.is_file():
            raise NinjaMissingError(f"{make_program} (the build's {CMAKE_MAKE_PROGRAM_KEY}) is not on disk")
        if CMAKE_HOME_DIRECTORY_KEY not in cache:
            raise CMakeHomeDirectoryUnsetError(f"{directory / CMAKE_CACHE} names no {CMAKE_HOME_DIRECTORY_KEY}")
        source_directory = pathlib.Path(cache[CMAKE_HOME_DIRECTORY_KEY]).resolve()
        if not directory.is_relative_to(source_directory):
            raise EditableBuildOutsideSourceTreeError(f"the editable build {directory} is not inside its source tree {source_directory}")
        known: Mapping[str, str] = cast("Mapping[str, str]", getattr(finder, "known_wheel_files", {}))
        wheel_files = {name: pathlib.Path(wheel_file).resolve() for name, wheel_file in known.items()}
        return cls(directory=directory, source_directory=source_directory, wheel_files=wheel_files, make_program=make_program)


def _cmake_cache(directory: pathlib.Path) -> dict[str, str]:
    """``KEY:TYPE=value`` entries of the build's CMake cache, keyed by ``KEY`` (and ``KEY:TYPE`` for the generator)."""
    cache = directory / CMAKE_CACHE
    if not cache.is_file():
        raise CMakeCacheMissingError(f"{cache} does not exist; the editable build was never configured there")
    entries = {}
    for line in cache.read_text(encoding="utf-8").splitlines():
        if line.startswith(("#", "//")) or "=" not in line:
            continue
        key, value = line.split("=", 1)
        entries[key] = value
        entries[key.split(":", 1)[0]] = value
    return entries


def native_closure(build: EditableBuild) -> tuple[str, ...]:
    """Paths, relative to the build's source tree, of every loaded ``compas_cgal`` extension module, sorted; each must be the editable build's file.

    Call it after the measured code has run: an extension imported lazily
    during the run is only in the closure once it is loaded.
    """
    paths = set()
    for name, module in list(sys.modules.items()):
        filename = getattr(module, "__file__", None)
        if not name.startswith(NATIVE_PACKAGE_PREFIX) or filename is None or not filename.endswith(EXTENSION_SUFFIXES):
            continue
        path = pathlib.Path(filename).resolve()
        mapped = build.wheel_files.get(name)
        if mapped != path:
            raise NativeModuleOutsideEditableBuildError(f"{name} is loaded from {path}, but the editable finder maps it to {mapped}")
        if not path.is_relative_to(build.directory):
            raise NativeModuleOutsideEditableBuildError(f"{name} is loaded from {path}, outside the editable build directory {build.directory}")
        paths.add(path.relative_to(build.source_directory).as_posix())
    return tuple(sorted(paths))


def capture_native(build: EditableBuild, paths: Sequence[str]) -> dict[str, bytes]:
    """The current bytes of each extension path (relative to the build's source tree)."""
    return {path: (build.source_directory / path).read_bytes() for path in paths}


def _ninja(build: EditableBuild, *arguments: str) -> subprocess.CompletedProcess[str]:
    command = [str(build.make_program), "-C", str(build.directory), *arguments]
    if not build.make_program.is_file():
        raise NinjaMissingError(f"{build.make_program} is not on disk")
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise NativeBuildToolError(f"{' '.join(command)} failed ({result.returncode}): {result.stderr.strip() or result.stdout.strip()}")
    return result


def stale_reasons(build: EditableBuild, paths: Sequence[str]) -> tuple[str, ...]:
    """Why the editable rebuild would change the extensions at ``paths``; empty when they are current."""
    targets: dict[str, list[str]] = {}
    for line in _ninja(build, "-t", "targets", "all").stdout.splitlines():
        target = line.rsplit(": ", 1)[0]
        targets.setdefault(pathlib.PurePosixPath(target).name, []).append(target)
    reasons, requested = [], []
    for path in paths:
        installed = build.source_directory / path
        matches = targets.get(installed.name, [])
        if len(matches) != 1:
            reasons.append(f"{path}: no build target produces it" if not matches else f"{path}: ambiguous build targets {matches}")
            continue
        target = matches[0]
        requested.append(target)
        if not (build.directory / target).is_file():
            reasons.append(f"{path}: its build target {target} was never built")
        elif (build.directory / target).read_bytes() != installed.read_bytes():
            reasons.append(f"{path}: differs from the build tree's {target} (the install step is not current)")
    if requested:
        dry_run = _ninja(build, "-n", "-d", "explain", *requested)
        lines = [line for line in dry_run.stdout.splitlines() if not line.startswith(NINJA_ENTERING)]
        if lines != [NINJA_UP_TO_DATE]:
            explained = [line for line in dry_run.stderr.splitlines() if line.strip()]
            reasons.append(f"ninja would rebuild {requested}: " + "; ".join(explained + lines))
    return tuple(reasons)


def require_current(build: EditableBuild, paths: Sequence[str]) -> None:
    """Raise unless the extensions at ``paths`` are what the editable rebuild would produce now."""
    reasons = stale_reasons(build, paths)
    if reasons:
        raise NativeBuildStaleError("the native build is stale; run the pixi task (it rebuilds) instead of build/pyrun.sh: " + " | ".join(reasons))


class _Timeval(ctypes.Structure):
    """macOS 64-bit ``struct timeval``: ``time_t tv_sec`` (8 bytes), ``__darwin_suseconds_t tv_usec`` (int32, then 4 padding)."""

    _fields_ = (("tv_sec", ctypes.c_int64), ("tv_usec", ctypes.c_int32))


# struct kinfo_proc opens with struct extern_proc kp_proc, which opens with the union p_un whose member
# __p_starttime is a struct timeval: the start sits at byte 0 (arm64 and x86_64 alike; <sys/sysctl.h>, <sys/proc.h>).
KINFO_PROC_STARTTIME_OFFSET = 0


def kernel_process_start_ns() -> int:
    """This process's start time (ns since the epoch) from the macOS kernel, at its microsecond precision.

    Raises:
        ProcessStartUnavailableError: sysctl failed (with its errno) or returned no record for this pid.
    """
    libc = ctypes.CDLL(None, use_errno=True)
    sysctl = libc.sysctl
    sysctl.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_uint, ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p, ctypes.c_size_t)
    sysctl.restype = ctypes.c_int
    mib = (ctypes.c_int * 4)(CTL_KERN, KERN_PROC, KERN_PROC_PID, os.getpid())
    size = ctypes.c_size_t(0)
    # First call sizes the record (the kernel's sizeof(struct kinfo_proc)), so no struct size is hard-coded here.
    if sysctl(mib, len(mib), None, ctypes.byref(size), None, 0) != 0:
        number = ctypes.get_errno()
        raise ProcessStartUnavailableError(f"sysctl KERN_PROC_PID sizing failed: errno {number} ({os.strerror(number)})")
    record = ctypes.create_string_buffer(size.value)
    if sysctl(mib, len(mib), record, ctypes.byref(size), None, 0) != 0:
        number = ctypes.get_errno()
        raise ProcessStartUnavailableError(f"sysctl KERN_PROC_PID failed: errno {number} ({os.strerror(number)})")
    if size.value < KINFO_PROC_STARTTIME_OFFSET + ctypes.sizeof(_Timeval):
        raise ProcessStartUnavailableError(f"sysctl KERN_PROC_PID returned {size.value} bytes, no kinfo_proc for pid {os.getpid()}")
    started = _Timeval.from_buffer_copy(record, KINFO_PROC_STARTTIME_OFFSET)
    return int(started.tv_sec) * NS_PER_S + int(started.tv_usec) * NS_PER_US


def process_start_ns() -> int:
    """This process's start time (ns since the epoch): the kernel's microseconds on macOS, ``ps``'s whole second elsewhere."""
    if sys.platform == "darwin":
        return kernel_process_start_ns()
    return ps_process_start_ns()


def ps_process_start_ns() -> int:
    """This process's start time (ns since the epoch), rounded down to the whole second ``ps`` reports.

    Raises:
        ProcessStartUnavailableError: ps failed, printed nothing, or printed a start that does not parse.
    """
    command = [*PS_START_COMMAND, str(os.getpid())]
    result = subprocess.run(command, capture_output=True, text=True, check=False, env={**os.environ, **PS_ENVIRONMENT})
    if result.returncode != 0 or not result.stdout.strip():
        raise ProcessStartUnavailableError(f"{' '.join(command)} failed ({result.returncode}): {result.stderr.strip()}")
    printed = " ".join(result.stdout.split())
    try:
        started = datetime.strptime(printed, PS_START_FORMAT).replace(tzinfo=timezone.utc)
    except ValueError as error:
        raise ProcessStartUnavailableError(f"{' '.join(command)} printed {printed!r}, not {PS_START_FORMAT!r}") from error
    return int(started.timestamp()) * NS_PER_S


def refuse_replaced_after_start(build: EditableBuild, paths: Sequence[str], started_ns: int) -> None:
    """Raise if any extension file at ``paths`` was written (ctime or mtime) at or after ``started_ns`` (a tie is refused)."""
    replaced = []
    for path in paths:
        status = (build.source_directory / path).stat()
        if max(status.st_ctime_ns, status.st_mtime_ns) >= started_ns:
            replaced.append(path)
    if replaced:
        raise NativeExtensionReplacedAfterStartError(f"written after this process started, so possibly not the loaded image: {replaced}; re-run")


@dataclass(frozen=True)
class NativeRun:
    """The start of a measured run: the editable build, the process start, the bytes of every extension loaded then."""

    build: EditableBuild
    started_ns: int
    loaded: Mapping[str, bytes]


def capture_current_native() -> NativeRun:
    """Before the run: read the loaded extensions' bytes, then refuse any written since the process start, then any stale."""
    build = EditableBuild.locate()
    closure = native_closure(build)
    started_ns = process_start_ns()
    loaded = capture_native(build, closure)
    refuse_replaced_after_start(build, closure, started_ns)
    require_current(build, closure)
    return NativeRun(build=build, started_ns=started_ns, loaded=loaded)


def native_fields(contents: Mapping[str, bytes]) -> dict[str, Any]:
    """The record fields: ``native_sha256`` over the (path, bytes) pairs and each extension's own SHA-256."""
    return {
        "native_sha256": digest_contents(contents),
        "natives": {path: hashlib.sha256(contents[path]).hexdigest() for path in sorted(contents)},
    }


def native_record_fields(run: NativeRun) -> dict[str, Any]:
    """After the run: one read of every loaded extension, unchanged since the start, none written since the process start, current; its fields."""
    closure = native_closure(run.build)
    contents = capture_native(run.build, tuple(sorted(set(closure) | set(run.loaded))))
    if {path: contents[path] for path in run.loaded} != dict(run.loaded):
        raise NativeBuildChangedError("a native extension changed on disk during the run; re-measure once the build settles")
    refuse_replaced_after_start(run.build, tuple(contents), run.started_ns)
    require_current(run.build, tuple(contents))
    return native_fields(contents)


def native_identity_line(record: Mapping[str, Any]) -> str:
    """The rendered identity: the digest prefix and the module count."""
    if "native_sha256" not in record or "natives" not in record:
        raise NativeIdentityMissingError("the record has no native extension identity; re-measure it through its pixi task")
    return f"native extensions sha256 `{record['native_sha256'][:12]}` ({len(record['natives'])} modules)"
