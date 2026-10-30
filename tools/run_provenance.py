"""Where a timing record came from: tree, tool, environment, machine and wall-clock interval.

A timing record (``pipeline-timing/2``, the pencil-cost record) is evidence
for the Task 1.2 generation-time gate only when the gate can check that its
records describe the trees they claim, ran in one environment on one machine,
and ran back to back. `RunProvenance` is that record. The fields are:

- the tree: its full commit, whether the working tree differed from it, and,
  if it did, a digest of the difference (tracked diff plus untracked files,
  without the gate session's own output records, `GATE_OUTPUT_RECORDS`). An
  export (a directory inside a repository that is not that repository's top
  level, e.g. a landing export under the ignored ``build/``) is identified by
  its content instead: the containing repository's HEAD as the parent commit,
  dirty, and a digest of every file in it (`ExportTreeIdentity`);
- the timing tool's source digest;
- the environment: ``pixi.lock``, the pixi environment name, the Python version
  and the compiled ``compas_cgal._stock_2`` extension's digest;
- the native code: ``native_sha256`` over every loaded extension, held
  unchanged from the run's start to its end, and ``native_sources_sha256`` of
  the measured tree's C++ sources (``tools.run_native_identity``);
- the machine: its host name and CPU model;
- the UTC start and end of the run.

The interpreter facts are read by `INTERPRETER_IDENTITY_SNIPPET`, a standalone
program that is run in-process here (`interpreter_identity_here`) and in a
subprocess for a tree whose tools predate this module
(`tools.pipeline_timing_1_1_base`). It only reads the environment and hashes the
extension file; it never touches the code being measured.
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from datetime import timezone
from pathlib import Path
from typing import Callable
from typing import Mapping
from typing import NamedTuple
from typing import Optional
from typing import Sequence

# Hex characters of a full git commit id (SHA-1).
GIT_REV_HEX_CHARS = 40
# Hex characters of a full SHA-256.
SHA256_HEX_CHARS = 64
# The timespec of a recorded UTC time: whole seconds are below every interval the gate checks.
UTC_TIMESPEC = "seconds"
_REV_PATTERN = re.compile(rf"[0-9a-f]{{{GIT_REV_HEX_CHARS}}}")
_SHA256_PATTERN = re.compile(rf"[0-9a-f]{{{SHA256_HEX_CHARS}}}")
_TEXT_FIELDS = ("environment", "python_version", "host", "cpu_model")
_SHA_FIELDS = ("timing_tool_sha256", "pixi_lock_sha256", "stock_2_sha256", "native_sha256", "native_sources_sha256")
_KEYS = frozenset({"rev", "dirty", "patch_sha256", *_SHA_FIELDS, *_TEXT_FIELDS, "started_utc", "ended_utc"})
# The gate session's own output records (repository-relative). `tree_identity` leaves them out: they
# are what the session writes, not the code it measures, and the documented sequence (timing, then
# pencil cost, then gate) rewrites the tracked timing records after the timing run's end read.
# Equal to the tools' default output paths (checked by the tests).
_RESULTS = Path("benchmarks") / "results"
GATE_OUTPUT_RECORDS = (
    _RESULTS / "pipeline_timing.json",
    _RESULTS / "pipeline_timing.md",
    _RESULTS / "pipeline_timing_1_1_session.json",
    _RESULTS / "pipeline_pencil_cost.json",
    _RESULTS / "pipeline_pencil_cost.md",
    _RESULTS / "pipeline_time_gate.md",
)
# The pathspec selecting the whole tree except `GATE_OUTPUT_RECORDS`, for ``git diff`` and ``git ls-files``.
_IDENTITY_PATHSPEC = (".", *(f":(top,exclude){path.as_posix()}" for path in GATE_OUTPUT_RECORDS))
# Prints one JSON line {environment, python_version, stock_2_sha256} for the interpreter that runs it.
INTERPRETER_IDENTITY_SNIPPET = """
import hashlib, importlib.util, json, os, sys
spec = importlib.util.find_spec("compas_cgal._stock_2")
if spec is None or spec.origin is None:
    raise SystemExit("compas_cgal._stock_2 is not importable in this interpreter")
environment = os.environ.get("PIXI_ENVIRONMENT_NAME")
if not environment:
    raise SystemExit("PIXI_ENVIRONMENT_NAME is unset: run through pixi")
with open(spec.origin, "rb") as extension:
    stock_2 = hashlib.sha256(extension.read()).hexdigest()
print(json.dumps({"environment": environment, "python_version": sys.version.split()[0], "stock_2_sha256": stock_2}))
"""


# Prefixed to an export's file-set digest, so it can never equal a worktree's difference digest:
# the two kinds of digest hash different things and must never pair.
EXPORT_DIGEST_DOMAIN = b"run-provenance/export-content/1\0"


class TreeIdentity(NamedTuple):
    """A git worktree's identity, as `tree_identity` reads it.

    Attributes:
        rev: Full commit id of the worktree's HEAD.
        dirty: Whether the working tree differs from *rev*.
        patch_sha256: Digest of that difference, None when clean.
    """

    rev: str
    dirty: bool
    patch_sha256: Optional[str]


class ExportTreeIdentity(TreeIdentity):
    """An export's identity: its content, not a commit.

    *rev* is the containing repository's HEAD when read (the parent the export was
    cut from, as far as git can say), *dirty* is always True, and *patch_sha256* is
    the domain-separated digest of every file in the export except
    `GATE_OUTPUT_RECORDS` and what its ``.gitignore`` files ignore. A commit in the
    containing repository moves *rev* and nothing else; it is not a change of the
    measured tree.
    """

    __slots__ = ()


TreeReader = Callable[[], TreeIdentity]


class TreeChangedDuringRunError(Exception):
    """The measured tree's identity at the end of a run differs from its start: the children imported mixed code."""


class MalformedProvenanceError(Exception):
    """A provenance record breaks its schema."""


class MachineIdentityError(Exception):
    """This machine's CPU model cannot be read."""


def _parse_utc(name: str, value: object) -> datetime:
    if not isinstance(value, str):
        raise MalformedProvenanceError(f"{name} is an ISO 8601 UTC time string.")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise MalformedProvenanceError(f"{name} {value!r} is not ISO 8601.") from error
    if parsed.utcoffset() is None:
        raise MalformedProvenanceError(f"{name} {value!r} carries no UTC offset.")
    return parsed


@dataclass(frozen=True)
class RunProvenance:
    """Where one timing run came from. Raw construction is unchecked; `build` validates.

    Attributes:
        rev: The full commit id of the measured tree.
        dirty: Whether the working tree differed from *rev*.
        patch_sha256: SHA-256 of that difference (`tree_identity`), None when clean.
        timing_tool_sha256: `file_set_sha256` of the timing tool's sources.
        pixi_lock_sha256: SHA-256 of the measured tree's ``pixi.lock``.
        environment: The pixi environment name.
        python_version: The interpreter version.
        stock_2_sha256: SHA-256 of the compiled ``compas_cgal._stock_2`` extension.
        native_sha256: Digest of every compas_cgal extension the run loaded, unchanged start to end (``tools.run_native_identity``).
        native_sources_sha256: Digest of the measured tree's native sources (``tools.run_native_identity.native_sources_sha256``).
        host: The host name.
        cpu_model: The CPU model.
        started_utc: The run's start, ISO 8601 with offset.
        ended_utc: The run's end, likewise.
    """

    rev: str
    dirty: bool
    patch_sha256: Optional[str]
    timing_tool_sha256: str
    pixi_lock_sha256: str
    environment: str
    python_version: str
    stock_2_sha256: str
    native_sha256: str
    native_sources_sha256: str
    host: str
    cpu_model: str
    started_utc: str
    ended_utc: str

    @classmethod
    def build(cls, raw: object) -> RunProvenance:
        """Validate a parsed JSON provenance.

        Raises:
            MalformedProvenanceError: Wrong keys, a malformed rev or digest, a
                dirty flag that disagrees with the patch digest, an empty text
                field, a time without offset, or an end before the start.
        """
        if not isinstance(raw, dict) or set(raw) != _KEYS:
            raise MalformedProvenanceError(f"a provenance has exactly the keys {sorted(_KEYS)}.")
        rev = raw["rev"]
        if not isinstance(rev, str) or _REV_PATTERN.fullmatch(rev) is None:
            raise MalformedProvenanceError(f"rev is the full {GIT_REV_HEX_CHARS}-character commit id, not {rev!r}.")
        dirty, patch = raw["dirty"], raw["patch_sha256"]
        if not isinstance(dirty, bool):
            raise MalformedProvenanceError("dirty is a boolean.")
        if dirty != (patch is not None) or (patch is not None and (not isinstance(patch, str) or _SHA256_PATTERN.fullmatch(patch) is None)):
            raise MalformedProvenanceError("patch_sha256 is a SHA-256 exactly when the tree is dirty.")
        for name in _SHA_FIELDS:
            value = raw[name]
            if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
                raise MalformedProvenanceError(f"{name} is {SHA256_HEX_CHARS} lowercase hex characters.")
        for name in _TEXT_FIELDS:
            if not isinstance(raw[name], str) or not raw[name]:
                raise MalformedProvenanceError(f"{name} is a non-empty string.")
        started, ended = _parse_utc("started_utc", raw["started_utc"]), _parse_utc("ended_utc", raw["ended_utc"])
        if ended < started:
            raise MalformedProvenanceError(f"ended_utc {raw['ended_utc']} precedes started_utc {raw['started_utc']}.")
        return cls(**{name: raw[name] for name in _KEYS})

    @property
    def started(self) -> datetime:
        """The start as an aware datetime."""
        return _parse_utc("started_utc", self.started_utc)

    @property
    def ended(self) -> datetime:
        """The end as an aware datetime."""
        return _parse_utc("ended_utc", self.ended_utc)

    def to_json(self) -> dict[str, object]:
        """The JSON-ready record."""
        return {
            "rev": self.rev,
            "dirty": self.dirty,
            "patch_sha256": self.patch_sha256,
            "timing_tool_sha256": self.timing_tool_sha256,
            "pixi_lock_sha256": self.pixi_lock_sha256,
            "environment": self.environment,
            "python_version": self.python_version,
            "stock_2_sha256": self.stock_2_sha256,
            "native_sha256": self.native_sha256,
            "native_sources_sha256": self.native_sources_sha256,
            "host": self.host,
            "cpu_model": self.cpu_model,
            "started_utc": self.started_utc,
            "ended_utc": self.ended_utc,
        }


def utc_now() -> str:
    """The current UTC time, ISO 8601 with offset, to whole seconds."""
    return datetime.now(timezone.utc).isoformat(timespec=UTC_TIMESPEC)


def bytes_set_sha256(contents: Mapping[Path, bytes]) -> str:
    """SHA-256 over each (path, bytes) of *contents*, in sorted path order, as its path then its bytes."""
    digest = hashlib.sha256()
    for path in sorted(contents, key=lambda path: path.as_posix()):
        digest.update(path.as_posix().encode() + b"\0")
        digest.update(contents[path] + b"\0")
    return digest.hexdigest()


def file_set_sha256(root: Path, paths: Sequence[Path]) -> str:
    """`bytes_set_sha256` of each of *paths* (relative to *root*) and its file bytes."""
    return bytes_set_sha256({path: (root / path).read_bytes() for path in paths})


def require_unchanged_tree(before: TreeIdentity, after: TreeIdentity, label: str) -> None:
    """Refuse a run whose tree identity changed between its start (*before*) and its end (*after*).

    Two export identities agree when their content digests do (the containing
    repository's HEAD may move during the run); a worktree identity agrees only
    with an equal worktree identity (there a HEAD move is a change of the tree).
    A worktree and an export identity never agree.

    Raises:
        TreeChangedDuringRunError: The identities differ.
    """
    exports = (isinstance(before, ExportTreeIdentity), isinstance(after, ExportTreeIdentity))
    if exports == (True, True):
        unchanged = before.patch_sha256 == after.patch_sha256
    else:
        unchanged = exports == (False, False) and tuple(before) == tuple(after)
    if not unchanged:
        raise TreeChangedDuringRunError(f"{label} changed during the run: {before} at the start, {after} at the end.")


def sha256_file(path: Path) -> str:
    """SHA-256 of *path*'s bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(root: Path, *args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=root, capture_output=True, check=True).stdout


def not_ignored(root: Path, paths: Sequence[Path]) -> list[Path]:
    """*paths* (relative to *root*) minus those *root*'s ``.gitignore`` files ignore, *root* taken as its own work tree.

    Raises:
        subprocess.CalledProcessError: ``git check-ignore`` failed (*root* is not inside a git repository).
    """
    listed = b"".join(path.as_posix().encode() + b"\0" for path in paths)
    # --no-index: by pattern only, whether or not a path is tracked.
    checked = subprocess.run(["git", "--work-tree=.", "check-ignore", "--no-index", "-z", "--stdin"], cwd=root, input=listed, capture_output=True)
    if checked.returncode not in (0, 1):  # 1: nothing ignored
        raise subprocess.CalledProcessError(checked.returncode, checked.args, checked.stdout, checked.stderr)
    ignored = {Path(path) for path in checked.stdout.decode().split("\0") if path}
    return [path for path in paths if path not in ignored]


def _export_files(root: Path) -> list[Path]:
    """Every file under the export *root* except `GATE_OUTPUT_RECORDS` and what its ``.gitignore`` files ignore."""
    candidates = sorted(path.relative_to(root) for path in root.rglob("*") if path.is_file())
    return [path for path in not_ignored(root, candidates) if path not in GATE_OUTPUT_RECORDS]


def tree_identity(root: Path) -> TreeIdentity:
    """The tree at *root*: (full commit id, dirty, difference digest or None).

    In a worktree (*root* is its repository's top level) the difference digest
    is SHA-256 over ``git diff HEAD --binary`` and every untracked, not-ignored
    file (path and bytes), both without `GATE_OUTPUT_RECORDS`: any edit to the
    measured code, committed or not, changes it, and writing the session's own
    output records does not. An export (*root* inside a repository but not its
    top level) is identified by its content (`ExportTreeIdentity`): git sees
    none of its files, so a worktree reading would name the containing HEAD as
    clean while the export holds other code.

    Raises:
        subprocess.CalledProcessError: *root* is not inside a git work tree, or (an export) ``git check-ignore`` failed.
    """
    rev = _git(root, "rev-parse", "HEAD").decode().strip()
    toplevel = Path(_git(root, "rev-parse", "--show-toplevel").decode().strip())
    if toplevel.resolve() != root.resolve():
        content = hashlib.sha256(EXPORT_DIGEST_DOMAIN)
        content.update(bytes.fromhex(file_set_sha256(root, _export_files(root))))
        return ExportTreeIdentity(rev, True, content.hexdigest())
    diff = _git(root, "diff", "HEAD", "--binary", "--", *_IDENTITY_PATHSPEC)
    untracked = sorted(path for path in _git(root, "ls-files", "--others", "--exclude-standard", "-z", "--", *_IDENTITY_PATHSPEC).decode().split("\0") if path)
    if not diff and not untracked:
        return TreeIdentity(rev, False, None)
    digest = hashlib.sha256(diff)
    digest.update(bytes.fromhex(file_set_sha256(root, [Path(path) for path in untracked])))
    return TreeIdentity(rev, True, digest.hexdigest())


def machine_identity() -> tuple[str, str]:
    """(host name, CPU model) of this machine.

    Raises:
        MachineIdentityError: The CPU model cannot be read (neither ``sysctl`` nor ``/proc/cpuinfo``).
    """
    host = platform.node()
    if sys.platform == "darwin":
        model = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=True).stdout.strip()
    else:
        cpuinfo = Path("/proc/cpuinfo")
        names = [line.split(":", 1)[1].strip() for line in cpuinfo.read_text().splitlines() if line.startswith("model name")] if cpuinfo.exists() else []
        model = names[0] if names else ""
    if not host or not model:
        raise MachineIdentityError(f"host {host!r}, CPU model {model!r}: both must be readable.")
    return host, model


def _parse_interpreter_identity(text: str) -> tuple[str, str, str]:
    payload = json.loads(text)
    return str(payload["environment"]), str(payload["python_version"]), str(payload["stock_2_sha256"])


def interpreter_identity_here() -> tuple[str, str, str]:
    """(pixi environment, Python version, ``_stock_2`` SHA-256) of this interpreter, by `INTERPRETER_IDENTITY_SNIPPET`."""
    completed = subprocess.run([sys.executable, "-c", INTERPRETER_IDENTITY_SNIPPET], capture_output=True, text=True, check=True)
    return _parse_interpreter_identity(completed.stdout)


def interpreter_identity_of(python: Sequence[str], cwd: Path, env: Optional[dict[str, str]] = None) -> tuple[str, str, str]:
    """The same facts for the interpreter *python* (a command prefix) run in *cwd*, e.g. another tree's pixi environment."""
    completed = subprocess.run([*python, "-c", INTERPRETER_IDENTITY_SNIPPET], cwd=cwd, env=env, capture_output=True, text=True, check=True)
    return _parse_interpreter_identity(completed.stdout.strip().splitlines()[-1])
