"""Turn record bytes into history entries: from a file on disk (``ingest``) or from every git version of every record (``backfill``).

``backfill`` walks ``git log --reverse`` of ``benchmarks/results`` on the
checked-out branch and reads every added or modified ``*.json`` with
``git show <rev>:<path>``, oldest first. A listed non-timing record is counted
and named in the result (never silently skipped); a record with no adapter
raises ``UnknownRecordFamilyError`` and stops the backfill.
"""

from __future__ import annotations

import hashlib
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from pathlib import PurePosixPath
from typing import Optional

from tools.perf_history_adapters import ADAPTER_VERSIONS
from tools.perf_history_adapters import MissingSiblingRecordError
from tools.perf_history_adapters import NotATimingRecordError
from tools.perf_history_adapters import SiblingReader
from tools.perf_history_adapters import adapt
from tools.perf_history_adapters import family_of
from tools.perf_history_entry import PerfEntry
from tools.perf_history_entry import RecordSource
from tools.perf_history_store import REPOSITORY

RESULTS_DIRECTORY = PurePosixPath("benchmarks/results")
RECORD_SUFFIX = ".json"
# Separates commits in the ``git log`` output backfill parses.
_COMMIT_MARK = "\x00"


class GitReadError(Exception):
    """A git command the backfill needs failed."""


class RecordOutsideRepositoryError(Exception):
    """``ingest`` was given a file outside the repository: its record path would not name it."""


def record_entries(record_path: str, data: bytes, sibling: SiblingReader, committed_in: Optional[str] = None, committed_at: Optional[str] = None) -> tuple[PerfEntry, ...]:
    """The history entries of one record version.

    Raises:
        Everything ``tools.perf_history_adapters.adapt`` raises.
    """
    adapted = adapt(record_path, data, sibling)
    family = family_of(record_path)
    source = RecordSource(
        family=family,
        record_path=record_path,
        record_sha256=hashlib.sha256(data).hexdigest(),
        committed_in=committed_in,
        committed_at=committed_at,
        adapter_version=ADAPTER_VERSIONS[family],
    )
    return tuple(PerfEntry.build(source, adapted.provenance, measurement) for measurement in adapted.measurements)


def file_entries(path: Path, root: Optional[Path] = None) -> tuple[PerfEntry, ...]:
    """The entries of the record file at ``path`` as it is on disk now; its record path is relative to ``root`` (default: the repository).

    Raises:
        RecordOutsideRepositoryError: ``path`` is not under ``root``.
        Everything ``record_entries`` raises.
    """
    root = REPOSITORY if root is None else root
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(root.resolve()).as_posix()
    except ValueError as error:
        raise RecordOutsideRepositoryError(f"{path} is not under {root}.") from error

    def sibling(name: str) -> bytes:
        candidate = resolved.parent / name
        if not candidate.is_file():
            raise MissingSiblingRecordError(f"{relative} needs {name} beside it.")
        return candidate.read_bytes()

    return record_entries(relative, resolved.read_bytes(), sibling)


def _git(root: Path, *args: str) -> bytes:
    completed = subprocess.run(["git", *args], cwd=root, capture_output=True, check=False)
    if completed.returncode != 0:
        raise GitReadError(f"git {' '.join(args)}: {completed.stderr.decode('utf-8', 'replace').strip()}")
    return completed.stdout


@dataclass(frozen=True)
class RecordVersion:
    """One committed version of one record file."""

    rev: str
    committed_at: str
    path: str


def record_versions(root: Path = REPOSITORY, results: PurePosixPath = RESULTS_DIRECTORY) -> tuple[RecordVersion, ...]:
    """Every added or modified ``*.json`` under ``results`` on the checked-out branch, oldest commit first.

    ``results`` is taken from the repository top (``:(top)`` pathspec), as ``git
    show <rev>:<path>`` reads it, so the result does not depend on where inside
    the repository ``root`` is (an export under ``build/`` reads the history of
    the repository it sits in).

    Raises:
        GitReadError: ``git log`` failed.
    """
    log = _git(root, "log", "--reverse", "--no-renames", "--diff-filter=AM", "--format=%x00%H%x09%cI", "--name-only", "--", f":(top){results.as_posix()}").decode("utf-8")
    versions: list[RecordVersion] = []
    for block in log.split(_COMMIT_MARK)[1:]:
        lines = [line for line in block.splitlines() if line]
        rev, committed_at = lines[0].split("\t")
        versions.extend(RecordVersion(rev, committed_at, path) for path in sorted(lines[1:]) if path.endswith(RECORD_SUFFIX))
    return tuple(versions)


@dataclass(frozen=True)
class BackfillResult:
    """What the backfill read.

    Attributes:
        entries: Every entry of every timing-record version, oldest commit first.
        timing_versions: The record versions the entries came from.
        non_timing_versions: Versions of listed non-timing records, not ingested.
    """

    entries: tuple[PerfEntry, ...]
    timing_versions: tuple[RecordVersion, ...]
    non_timing_versions: tuple[RecordVersion, ...]


def backfill_entries(versions: Sequence[RecordVersion], root: Path = REPOSITORY) -> BackfillResult:
    """The entries of every timing-record version in ``versions``.

    Raises:
        UnknownRecordFamilyError: a version's file has no adapter and is not listed as non-timing.
        GitReadError: ``git show`` failed.
        Everything ``record_entries`` raises.
    """
    entries: list[PerfEntry] = []
    timing: list[RecordVersion] = []
    non_timing: list[RecordVersion] = []
    for version in versions:
        try:
            family_of(version.path)
        except NotATimingRecordError:
            non_timing.append(version)
            continue
        directory = PurePosixPath(version.path).parent

        def sibling(name: str, rev: str = version.rev, directory: PurePosixPath = directory) -> bytes:
            try:
                return _git(root, "show", f"{rev}:{(directory / name).as_posix()}")
            except GitReadError as error:
                raise MissingSiblingRecordError(f"{rev}: {(directory / name).as_posix()}: {error}") from error

        data = _git(root, "show", f"{version.rev}:{version.path}")
        entries.extend(record_entries(version.path, data, sibling, committed_in=version.rev, committed_at=version.committed_at))
        timing.append(version)
    return BackfillResult(tuple(entries), tuple(timing), tuple(non_timing))
