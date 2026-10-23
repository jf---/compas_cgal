"""The committed, append-only performance history: ``benchmarks/results/perf_history.jsonl``.

One ``PerfEntry`` or ``Rebaseline`` per line. Lines are only ever appended:
``append_entries`` and ``append_rebaseline`` write in append mode and then
re-read the file to verify the bytes that were there before are unchanged
(``AppendOnlyViolationError`` otherwise). An entry whose ``identity`` (the
same record read by the same adapter version) is already present is a
duplicate and is not written again, so re-ingesting a record is a no-op.
Loading re-derives every line's digest and quiet verdict (an edited line is
refused) and refuses a repeated digest or identity.

A wrong historic entry is never rewritten. Its adapter is corrected and its
family's ``ADAPTER_VERSIONS`` entry bumped; ``backfill`` then appends the
corrected entries (a new identity), and readers keep only the newest adapter
version's entries per record version
(``tools.perf_history_baseline.current_entries``).
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from tools.perf_history_entry import MalformedEntryError
from tools.perf_history_entry import PerfEntry
from tools.perf_history_rebaseline import REBASELINE_SCHEMA
from tools.perf_history_rebaseline import Rebaseline

REPOSITORY = Path(__file__).resolve().parents[1]
HISTORY_PATH = REPOSITORY / "benchmarks" / "results" / "perf_history.jsonl"


class AppendOnlyViolationError(Exception):
    """The history's existing bytes changed while appending."""


class DuplicateHistoryEntryError(Exception):
    """Two history lines carry the same ``entry_sha256``."""


class UnterminatedHistoryError(Exception):
    """The history does not end in a newline: its last line was cut off."""


@dataclass(frozen=True)
class AppendResult:
    """What one append did.

    Attributes:
        added: Entries written, in order.
        duplicates: Entries already present (by ``entry_sha256``), not written.
    """

    added: tuple[PerfEntry, ...]
    duplicates: tuple[PerfEntry, ...]


@dataclass(frozen=True)
class History:
    """The history file's lines: measurement entries and rebaselines, each in file order."""

    entries: tuple[PerfEntry, ...]
    rebaselines: tuple[Rebaseline, ...]


def _line(text: str) -> Union[PerfEntry, Rebaseline]:
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as error:
        raise MalformedEntryError(f"a history line is JSON: {error}") from error
    if isinstance(raw, dict) and raw.get("schema") == REBASELINE_SCHEMA:
        return Rebaseline.from_json(raw)
    return PerfEntry.from_json_line(text)


def load_history(path: Path = HISTORY_PATH) -> History:
    """Every line of the history at ``path``, in file order; empty when the file does not exist.

    Raises:
        UnterminatedHistoryError: the file does not end in a newline.
        DuplicateHistoryEntryError: an entry digest or identity repeats.
        MalformedEntryError, MalformedRebaselineError: a line is neither an entry nor a rebaseline.
        EntryDigestMismatchError: an entry line was edited.
        RebaselineTargetError: a rebaseline names an entry the history does not hold.
    """
    if not path.exists():
        return History((), ())
    text = path.read_text(encoding="utf-8")
    if text and not text.endswith("\n"):
        raise UnterminatedHistoryError(f"{path}: the last line is not newline-terminated.")
    lines = [_line(line) for line in text.splitlines()]
    entries = tuple(line for line in lines if isinstance(line, PerfEntry))
    rebaselines = tuple(line for line in lines if isinstance(line, Rebaseline))
    digests: set[str] = set()
    identities: set[tuple[str, str, str, str, int]] = set()
    for entry in entries:
        if entry.entry_sha256 in digests or entry.identity in identities:
            raise DuplicateHistoryEntryError(f"{path}: entry {entry.entry_sha256} ({entry.identity}) appears twice.")
        digests.add(entry.entry_sha256)
        identities.add(entry.identity)
    for rebaseline in rebaselines:
        rebaseline.target(entries)
    return History(entries, rebaselines)


def _append_text(path: Path, text: str) -> None:
    before = path.read_bytes() if path.exists() else b""
    if text:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
    after = path.read_bytes() if path.exists() else b""
    if not after.startswith(before):
        raise AppendOnlyViolationError(f"{path}: the existing {len(before)} bytes changed during the append.")


def append_entries(entries: Iterable[PerfEntry], path: Path = HISTORY_PATH) -> AppendResult:
    """Append the entries not yet in the history, verifying the existing bytes are untouched.

    An entry is already present when an entry with the same ``identity`` (record
    path, record SHA-256, workload, stage, adapter version) is: the same record
    read by the same adapter version. A corrected adapter (a higher version)
    therefore appends new entries beside the old ones, never over them.

    Raises:
        AppendOnlyViolationError: the bytes before the append are not a prefix of the bytes after.
        (and everything ``load_history`` raises for the existing file)
    """
    present = {entry.identity for entry in load_history(path).entries}
    added: list[PerfEntry] = []
    duplicates: list[PerfEntry] = []
    for entry in entries:
        if entry.identity in present:
            duplicates.append(entry)
            continue
        present.add(entry.identity)
        added.append(entry)
    _append_text(path, "".join(entry.to_json_line() for entry in added))
    return AppendResult(tuple(added), tuple(duplicates))


def append_rebaseline(rebaseline: Rebaseline, path: Path = HISTORY_PATH) -> None:
    """Append ``rebaseline``, verifying it names an entry of the history and the existing bytes are untouched.

    Raises:
        RebaselineTargetError: the history does not hold the accepted entry.
        AppendOnlyViolationError: the existing bytes changed.
    """
    rebaseline.target(load_history(path).entries)
    _append_text(path, rebaseline.to_json_line())
