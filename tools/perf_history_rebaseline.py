"""An explicit acceptance of a slower baseline: one appended ``rebaseline`` line of the history.

The baseline of a key is the best (fastest) quiet entry, so a slower result
never becomes the baseline by being recorded. Accepting one takes a
``Rebaseline``: it names the accepted quiet entry by its ``entry_sha256``, says
why, and when. From then on ``tools.perf_history_baseline.baselines`` only
considers the key's quiet entries measured at or after the accepted entry
(compared as instants), and takes the best of those.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import TYPE_CHECKING
from typing import Any

from tools.perf_history_entry import PerfEntry
from tools.perf_history_entry import parse_instant

if TYPE_CHECKING:
    from tools.perf_history_store import History

REBASELINE_SCHEMA = "perf-history-rebaseline/1"
_KEYS = ("schema", "workload", "stage", "machine", "entry_sha256", "reason", "recorded_at")


class RebaselineError(Exception):
    """Base of the rebaseline's failure modes."""


class RebaselineTargetError(RebaselineError):
    """The named entry is not in the history, or is not quiet: only a quiet entry can be a baseline."""


class RebaselineKeyMismatchError(RebaselineError):
    """A rebaseline's (workload, stage, machine) is not the key of the entry it names."""


class RebaselineReasonError(RebaselineError):
    """The rebaseline gives no reason."""


class MalformedRebaselineError(RebaselineError):
    """A history line is not a well-formed rebaseline, or its digest is not its content's."""


@dataclass(frozen=True)
class Rebaseline:
    """The accepted baseline of one (workload, stage, machine). Raw construction is unchecked; ``build`` validates.

    Attributes:
        key: (workload, stage, machine), taken from the accepted entry.
        entry_sha256: The accepted quiet entry.
        reason: Why the slower baseline is accepted.
        recorded_at: When it was accepted (ISO 8601 with offset).
        rebaseline_sha256: SHA-256 of the line's content.
    """

    key: tuple[str, str, str]
    entry_sha256: str
    reason: str
    recorded_at: str
    rebaseline_sha256: str

    @classmethod
    def build(cls, history: History, entry_sha256: str, reason: str, recorded_at: str) -> Rebaseline:
        """Raises:
        RebaselineTargetError: no entry of ``history`` has ``entry_sha256``, or it is indicative.
        RebaselineReasonError: ``reason`` is blank.
        NaiveTimestampError: ``recorded_at`` has no UTC offset.
        """
        matches = [entry for entry in history.entries if entry.entry_sha256 == entry_sha256]
        if not matches:
            raise RebaselineTargetError(f"no history entry has entry_sha256 {entry_sha256}.")
        target = matches[0]
        if not target.quiet:
            raise RebaselineTargetError(f"{target.key}: entry {entry_sha256} is indicative ({target.indicative_reason}); only a quiet entry can be a baseline.")
        if not reason.strip():
            raise RebaselineReasonError("a rebaseline says why the slower baseline is accepted.")
        parse_instant(recorded_at)
        return cls._sealed(target.key, entry_sha256, reason, recorded_at)

    @classmethod
    def _sealed(cls, key: tuple[str, str, str], entry_sha256: str, reason: str, recorded_at: str) -> Rebaseline:
        return cls(key, entry_sha256, reason, recorded_at, _digest(_content(key, entry_sha256, reason, recorded_at)))

    def target(self, entries: tuple[PerfEntry, ...]) -> PerfEntry:
        """The accepted entry among ``entries``.

        Raises:
            RebaselineTargetError: it is not there (the history lost it).
            RebaselineKeyMismatchError: it is there under another (workload, stage, machine).
        """
        for entry in entries:
            if entry.entry_sha256 == self.entry_sha256:
                if entry.key != self.key:
                    raise RebaselineKeyMismatchError(f"rebaseline for {self.key} names entry {self.entry_sha256} of {entry.key}.")
                return entry
        raise RebaselineTargetError(f"{self.key}: the accepted entry {self.entry_sha256} is not in the history.")

    def to_json_line(self) -> str:
        document = dict(_content(self.key, self.entry_sha256, self.reason, self.recorded_at))
        document["rebaseline_sha256"] = self.rebaseline_sha256
        return json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"

    @classmethod
    def from_json(cls, raw: object) -> Rebaseline:
        """Raises:
        MalformedRebaselineError: ``raw`` is not a rebaseline line, or its digest is not its content's.
        NaiveTimestampError: its time has no UTC offset.
        """
        if not isinstance(raw, dict) or raw.get("schema") != REBASELINE_SCHEMA or set(raw) != {*_KEYS, "rebaseline_sha256"}:
            raise MalformedRebaselineError(f"a rebaseline line is an object of schema {REBASELINE_SCHEMA!r} with the keys {sorted({*_KEYS, 'rebaseline_sha256'})}.")
        for name in (*_KEYS, "rebaseline_sha256"):
            if not isinstance(raw[name], str):
                raise MalformedRebaselineError(f"{name} is a string, not {raw[name]!r}.")
        parse_instant(raw["recorded_at"])
        rebaseline = cls._sealed((raw["workload"], raw["stage"], raw["machine"]), raw["entry_sha256"], raw["reason"], raw["recorded_at"])
        if rebaseline.rebaseline_sha256 != raw["rebaseline_sha256"]:
            raise MalformedRebaselineError(f"{rebaseline.key}: the stored rebaseline digest is not its content's (the line was edited).")
        return rebaseline


def _content(key: tuple[str, str, str], entry_sha256: str, reason: str, recorded_at: str) -> dict[str, Any]:
    workload, stage, machine = key
    return {"schema": REBASELINE_SCHEMA, "workload": workload, "stage": stage, "machine": machine, "entry_sha256": entry_sha256, "reason": reason, "recorded_at": recorded_at}


def _digest(content: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
