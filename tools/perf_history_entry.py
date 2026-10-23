"""One performance-history entry: a measured (workload, stage), where it came from, and its content address.

An entry is what ``tools.perf_history_adapters`` extracts from one timing
record (a ``Measurement`` plus the record's ``RecordProvenance``), stamped
with the record it was read from (``RecordSource``). Its ``entry_sha256`` is
the SHA-256 of everything the record determines (``content_fields``); where the
record was later found in git (``committed_in``, ``committed_at``) is context
only and is left out, so the same record bytes give the same entry whichever
way they are ingested, and the history dedupes on that hash.

Quiet is derived, never passed: an entry is quiet exactly when the 1-minute
load average was recorded before and after the measurement, both are at most
``tools.held_generation_clock.QUIET_LOAD_AVERAGE_MAX``, and the record itself
does not mark the run indicative. Every other entry is indicative and says why
(``indicative_reason``).
"""

from __future__ import annotations

import enum
import hashlib
import json
import math
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from typing import Any
from typing import NewType
from typing import Optional

from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX

# Seconds of one clock (process CPU or wall), as a record reports them.
Seconds = NewType("Seconds", float)
# The 1-minute load average (runnable processes) read before or after a measurement.
LoadAverage1m = NewType("LoadAverage1m", float)

# The machine key of an entry whose record names no CPU.
UNRECORDED_MACHINE = "unrecorded"
# The schema tag every new history line carries: /1 plus the adapter version that read the record.
ENTRY_SCHEMA = "perf-history-entry/2"
# Lines written before adapter versions existed; read as adapter version 1, never written again.
LEGACY_ENTRY_SCHEMA = "perf-history-entry/1"
# The adapter version of a /1 line.
LEGACY_ADAPTER_VERSION = 1


class PerfHistoryEntryError(Exception):
    """Base of the entry's failure modes."""


class InvalidClockError(PerfHistoryEntryError):
    """A clock summary has no runs, a negative, non-finite or non-numeric time, or minimum above maximum."""


class NoClockError(PerfHistoryEntryError):
    """A measurement carries neither a CPU nor a wall clock."""


class MalformedEntryError(PerfHistoryEntryError):
    """A history line is not a well-formed entry of ``ENTRY_SCHEMA``."""


class NaiveTimestampError(PerfHistoryEntryError):
    """A recorded time carries no UTC offset, so it cannot be ordered against the others."""


class EntryDigestMismatchError(PerfHistoryEntryError):
    """A history line's ``entry_sha256`` is not the digest of its content: the line was edited."""


class ClockStatistic(enum.Enum):
    """What a clock summary's minimum and maximum are over.

    RUNS: the record lists every run; minimum and maximum are over them.
    SINGLE: the record holds one run.
    REPORTED_MEDIAN: the record kept only the median of ``runs`` runs; minimum = maximum = that median.
    SUM_OF_SINGLES: the sum of single runs of ``runs`` sub-items (one total).
    SUM_OF_MEDIANS: the sum of per-item medians over ``runs`` repeats each (one total).
    """

    RUNS = "runs"
    SINGLE = "single"
    REPORTED_MEDIAN = "reported_median"
    SUM_OF_SINGLES = "sum_of_singles"
    SUM_OF_MEDIANS = "sum_of_medians"


@dataclass(frozen=True)
class ClockSummary:
    """One clock of one measured stage. Raw construction is unchecked; use the factories.

    Attributes:
        minimum: The compared value (interference only adds time).
        maximum: Reported beside it; equal to ``minimum`` unless ``statistic`` is RUNS.
        runs: How many runs (or summed items) the value covers.
        statistic: What ``minimum`` and ``maximum`` are over.
    """

    minimum: Seconds
    maximum: Seconds
    runs: int
    statistic: ClockStatistic

    @classmethod
    def build(cls, minimum: float, maximum: float, runs: int, statistic: ClockStatistic) -> ClockSummary:
        """Raises:
        InvalidClockError: a time is negative, non-finite or not a number, ``runs`` < 1, or minimum > maximum.
        """
        for value in (minimum, maximum):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise InvalidClockError(f"a clock reads a finite non-negative number of seconds, not {value!r}.")
        if isinstance(runs, bool) or not isinstance(runs, int) or runs < 1:
            raise InvalidClockError(f"a clock covers at least one run, not {runs!r}.")
        if minimum > maximum:
            raise InvalidClockError(f"minimum {minimum!r} exceeds maximum {maximum!r}.")
        if statistic is not ClockStatistic.RUNS and minimum != maximum:
            raise InvalidClockError(f"a {statistic.value} clock holds one value, not {minimum!r}..{maximum!r}.")
        return cls(Seconds(float(minimum)), Seconds(float(maximum)), runs, statistic)

    @classmethod
    def over_runs(cls, seconds: Sequence[float]) -> ClockSummary:
        """Minimum and maximum over every listed run (one run is SINGLE).

        Raises:
            InvalidClockError: ``seconds`` is empty or holds an invalid time.
        """
        if not seconds:
            raise InvalidClockError("a clock over runs needs at least one run.")
        if len(seconds) == 1:
            return cls.single(seconds[0])
        return cls.build(min(seconds), max(seconds), len(seconds), ClockStatistic.RUNS)

    @classmethod
    def single(cls, seconds: float) -> ClockSummary:
        """One run."""
        return cls.build(seconds, seconds, 1, ClockStatistic.SINGLE)

    @classmethod
    def reported_median(cls, seconds: float, runs: int) -> ClockSummary:
        """The record kept only the median of ``runs`` runs."""
        return cls.build(seconds, seconds, runs, ClockStatistic.REPORTED_MEDIAN)

    @classmethod
    def total(cls, seconds: float, items: int, statistic: ClockStatistic) -> ClockSummary:
        """A sum over ``items`` sub-items (``SUM_OF_SINGLES`` or ``SUM_OF_MEDIANS``)."""
        return cls.build(seconds, seconds, items, statistic)

    @property
    def spread_ratio(self) -> float:
        """Maximum over minimum (1 for one value); the run-to-run spread the regression threshold is derived from."""
        if self.minimum == self.maximum:
            return 1.0
        return self.maximum / self.minimum if self.minimum > 0 else math.inf

    def to_json(self) -> dict[str, Any]:
        return {"minimum": self.minimum, "maximum": self.maximum, "runs": self.runs, "statistic": self.statistic.value}

    @classmethod
    def from_json(cls, raw: object) -> ClockSummary:
        """Raises:
        MalformedEntryError: ``raw`` is not a clock object.
        InvalidClockError: its values are invalid.
        """
        if not isinstance(raw, dict) or set(raw) != {"minimum", "maximum", "runs", "statistic"}:
            raise MalformedEntryError(f"a clock is an object with minimum, maximum, runs and statistic, not {raw!r}.")
        try:
            statistic = ClockStatistic(raw["statistic"])
        except ValueError as error:
            raise MalformedEntryError(f"unknown clock statistic {raw['statistic']!r}.") from error
        return cls.build(raw["minimum"], raw["maximum"], raw["runs"], statistic)


@dataclass(frozen=True)
class RecordProvenance:
    """What a timing record says about the tree and machine it was measured on; every field may be unrecorded.

    Attributes:
        commit: The measured commit.
        dirty: Whether the measured tree differed from ``commit``.
        source_sha256: The record's source-closure digest.
        native_sha256: The record's native-extension digest.
        machine: The CPU brand, or ``UNRECORDED_MACHINE``.
        recorded_at: When the record was measured (ISO 8601), as the record says.
    """

    commit: Optional[str] = None
    dirty: Optional[bool] = None
    source_sha256: Optional[str] = None
    native_sha256: Optional[str] = None
    machine: str = UNRECORDED_MACHINE
    recorded_at: Optional[str] = None


@dataclass(frozen=True)
class Measurement:
    """One measured (workload, stage) as one record reports it.

    Attributes:
        workload: The measured workload (case) name.
        stage: The measured stage, qualified by the record family that measured it,
            so harnesses that count different work never share a baseline.
        cpu: Process CPU, when the record has it.
        wall: Wall clock, when the record has it.
        load_before: 1-minute load average before the measurement, when recorded.
        load_after: 1-minute load average after it, when recorded.
        record_indicative: Why the record itself marks the run indicative, else None.
    """

    workload: str
    stage: str
    cpu: Optional[ClockSummary]
    wall: Optional[ClockSummary]
    load_before: Optional[LoadAverage1m]
    load_after: Optional[LoadAverage1m]
    record_indicative: Optional[str] = None


@dataclass(frozen=True)
class RecordSource:
    """The record an entry was read from.

    Attributes:
        family: The record family (``tools.perf_history_adapters.FAMILIES`` key).
        record_path: The record's repository-relative path.
        record_sha256: SHA-256 of the record's bytes.
        committed_in: The git commit the bytes were read from (backfill), else None. Context only.
        committed_at: That commit's committer date, else None. Context only.
        adapter_version: The version of the family's adapter that read the record
            (``tools.perf_history_adapters.ADAPTER_VERSIONS``). A corrected adapter
            bumps it; the corrected entries are appended beside the old ones and
            supersede them (``tools.perf_history_baseline.current_entries``).
    """

    family: str
    record_path: str
    record_sha256: str
    committed_in: Optional[str] = None
    committed_at: Optional[str] = None
    adapter_version: int = LEGACY_ADAPTER_VERSION


def parse_instant(text: str) -> datetime:
    """An ISO 8601 time with a UTC offset, as an aware datetime (compared as an instant, never as text).

    Raises:
        NaiveTimestampError: ``text`` has no offset or is not ISO 8601.
    """
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError as error:
        raise NaiveTimestampError(f"{text!r} is not an ISO 8601 time.") from error
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise NaiveTimestampError(f"{text!r} has no UTC offset.")
    return parsed


def indicative_reason(load_before: Optional[float], load_after: Optional[float], record_indicative: Optional[str]) -> Optional[str]:
    """Why a measurement is indicative, or None when it is quiet (``QUIET_LOAD_AVERAGE_MAX``, inclusive)."""
    if record_indicative is not None:
        return record_indicative
    if load_before is None and load_after is None:
        return "no load recorded"
    if load_before is None:
        return "no load recorded before the measurement"
    if load_after is None:
        return "no load recorded after the measurement"
    peak = max(load_before, load_after)
    if peak > QUIET_LOAD_AVERAGE_MAX:
        return f"1-minute load {peak:.2f} > {QUIET_LOAD_AVERAGE_MAX}"
    return None


def _optional_load(value: object) -> Optional[LoadAverage1m]:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise MalformedEntryError(f"a load average is a finite non-negative number or null, not {value!r}.")
    return LoadAverage1m(float(value))


def _optional_text(name: str, value: object) -> Optional[str]:
    if value is not None and not isinstance(value, str):
        raise MalformedEntryError(f"{name} is a string or null, not {value!r}.")
    return value


@dataclass(frozen=True)
class PerfEntry:
    """One line of the performance history. Raw construction is unchecked; ``build`` and ``from_json`` validate.

    Attributes:
        source: The record the entry was read from.
        provenance: What that record says about its tree and machine.
        measurement: The measured stage.
        quiet: Derived by ``indicative_reason``; only quiet entries set a baseline.
        indicative_reason: Why the entry is indicative, else None.
        entry_sha256: SHA-256 of ``content_fields`` (integrity; a repeated digest is refused).
        schema: ``ENTRY_SCHEMA``, or ``LEGACY_ENTRY_SCHEMA`` for a line read from before adapter versions.
    """

    source: RecordSource
    provenance: RecordProvenance
    measurement: Measurement
    quiet: bool
    indicative_reason: Optional[str]
    entry_sha256: str
    schema: str = ENTRY_SCHEMA

    @classmethod
    def build(cls, source: RecordSource, provenance: RecordProvenance, measurement: Measurement, schema: str = ENTRY_SCHEMA) -> PerfEntry:
        """Raises:
        NoClockError: the measurement has neither clock.
        MalformedEntryError: a legacy-schema entry with an adapter version other than 1, or an unknown schema.
        NaiveTimestampError: ``recorded_at`` or ``committed_at`` has no UTC offset or is not ISO 8601.
        """
        if measurement.cpu is None and measurement.wall is None:
            raise NoClockError(f"{source.record_path}: {measurement.workload} / {measurement.stage} has neither a CPU nor a wall clock.")
        if schema not in (ENTRY_SCHEMA, LEGACY_ENTRY_SCHEMA):
            raise MalformedEntryError(f"unknown entry schema {schema!r}.")
        if schema == LEGACY_ENTRY_SCHEMA and source.adapter_version != LEGACY_ADAPTER_VERSION:
            raise MalformedEntryError(f"a {LEGACY_ENTRY_SCHEMA} entry is adapter version {LEGACY_ADAPTER_VERSION}.")
        for stamp in (provenance.recorded_at, source.committed_at):
            if stamp is not None:
                parse_instant(stamp)  # an unorderable time never enters the history (RV40 N1)
        reason = indicative_reason(measurement.load_before, measurement.load_after, measurement.record_indicative)
        content = _content_fields(source, provenance, measurement, reason, schema)
        return cls(source, provenance, measurement, reason is None, reason, _digest(content), schema)

    @property
    def identity(self) -> tuple[str, str, str, str, int]:
        """What was read and how: (record path, record SHA-256, workload, stage, adapter version); the dedupe key."""
        return (self.source.record_path, self.source.record_sha256, self.measurement.workload, self.measurement.stage, self.source.adapter_version)

    @property
    def key(self) -> tuple[str, str, str]:
        """The baseline key: (workload, stage, machine)."""
        return (self.measurement.workload, self.measurement.stage, self.provenance.machine)

    @property
    def gated_clock(self) -> ClockSummary:
        """The compared clock: CPU when recorded, else wall."""
        clock = self.measurement.cpu if self.measurement.cpu is not None else self.measurement.wall
        assert clock is not None  # build refuses an entry with neither
        return clock

    @property
    def gated_metric(self) -> str:
        """``cpu`` or ``wall``: which clock ``gated_clock`` is."""
        return "cpu" if self.measurement.cpu is not None else "wall"

    @property
    def when(self) -> Optional[datetime]:
        """The instant of the measurement: the record's ``recorded_at``, else the commit date it was found in, else None.

        Raises:
            NaiveTimestampError: the time has no UTC offset.
        """
        text = self.provenance.recorded_at or self.source.committed_at
        return None if text is None else parse_instant(text)

    @property
    def when_text(self) -> str:
        """The time as recorded (for display), or empty."""
        return self.provenance.recorded_at or self.source.committed_at or ""

    def to_json_line(self) -> str:
        """The entry as one canonical JSON line, newline-terminated."""
        document = dict(_content_fields(self.source, self.provenance, self.measurement, self.indicative_reason, self.schema))
        document["quiet"] = self.quiet
        document["committed_in"] = self.source.committed_in
        document["committed_at"] = self.source.committed_at
        document["entry_sha256"] = self.entry_sha256
        return json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"

    @classmethod
    def from_json_line(cls, line: str) -> PerfEntry:
        """Raises:
        MalformedEntryError: the line is not an ``ENTRY_SCHEMA`` or ``LEGACY_ENTRY_SCHEMA`` entry.
        EntryDigestMismatchError: its ``entry_sha256`` or ``quiet`` is not what its content gives.
        InvalidClockError: a clock is invalid.
        """
        try:
            raw = json.loads(line)
        except json.JSONDecodeError as error:
            raise MalformedEntryError(f"a history line is JSON: {error}") from error
        if not isinstance(raw, dict) or raw.get("schema") not in (ENTRY_SCHEMA, LEGACY_ENTRY_SCHEMA):
            raise MalformedEntryError(f"a history entry line is an object of schema {ENTRY_SCHEMA!r} or {LEGACY_ENTRY_SCHEMA!r}.")
        schema = str(raw["schema"])
        expected = set(_content_keys(schema)) | {"quiet", "committed_in", "committed_at", "entry_sha256"}
        if set(raw) != expected:
            raise MalformedEntryError(f"a history line has exactly the keys {sorted(expected)}, not {sorted(raw)}.")
        for name in ("family", "record_path", "record_sha256", "workload", "stage", "machine", "entry_sha256"):
            if not isinstance(raw[name], str):
                raise MalformedEntryError(f"{name} is a string, not {raw[name]!r}.")
        if raw["dirty"] is not None and not isinstance(raw["dirty"], bool):
            raise MalformedEntryError(f"dirty is a boolean or null, not {raw['dirty']!r}.")
        adapter_version = raw.get("adapter_version", LEGACY_ADAPTER_VERSION)
        if isinstance(adapter_version, bool) or not isinstance(adapter_version, int) or adapter_version < LEGACY_ADAPTER_VERSION:
            raise MalformedEntryError(f"adapter_version is a positive integer, not {adapter_version!r}.")
        source = RecordSource(
            family=raw["family"],
            record_path=raw["record_path"],
            record_sha256=raw["record_sha256"],
            committed_in=_optional_text("committed_in", raw["committed_in"]),
            committed_at=_optional_text("committed_at", raw["committed_at"]),
            adapter_version=adapter_version,
        )
        provenance = RecordProvenance(
            commit=_optional_text("commit", raw["commit"]),
            dirty=raw["dirty"],
            source_sha256=_optional_text("source_sha256", raw["source_sha256"]),
            native_sha256=_optional_text("native_sha256", raw["native_sha256"]),
            machine=raw["machine"],
            recorded_at=_optional_text("recorded_at", raw["recorded_at"]),
        )
        measurement = Measurement(
            workload=raw["workload"],
            stage=raw["stage"],
            cpu=None if raw["cpu"] is None else ClockSummary.from_json(raw["cpu"]),
            wall=None if raw["wall"] is None else ClockSummary.from_json(raw["wall"]),
            load_before=_optional_load(raw["load_before"]),
            load_after=_optional_load(raw["load_after"]),
            record_indicative=_optional_text("record_indicative", raw["record_indicative"]),
        )
        entry = cls.build(source, provenance, measurement, schema)
        if entry.entry_sha256 != raw["entry_sha256"] or entry.quiet != raw["quiet"] or entry.indicative_reason != raw["indicative_reason"]:
            raise EntryDigestMismatchError(f"{source.record_path}: {measurement.workload} / {measurement.stage}: the stored digest or quiet verdict is not its content's.")
        return entry


_LEGACY_CONTENT_KEYS = (
    "schema",
    "family",
    "record_path",
    "record_sha256",
    "workload",
    "stage",
    "cpu",
    "wall",
    "load_before",
    "load_after",
    "record_indicative",
    "indicative_reason",
    "commit",
    "dirty",
    "source_sha256",
    "native_sha256",
    "machine",
    "recorded_at",
)
_CONTENT_KEYS = (*_LEGACY_CONTENT_KEYS, "adapter_version")


def _content_keys(schema: str) -> tuple[str, ...]:
    return _LEGACY_CONTENT_KEYS if schema == LEGACY_ENTRY_SCHEMA else _CONTENT_KEYS


def _content_fields(source: RecordSource, provenance: RecordProvenance, measurement: Measurement, reason: Optional[str], schema: str) -> Mapping[str, Any]:
    values: dict[str, Any] = {
        "schema": schema,
        "family": source.family,
        "record_path": source.record_path,
        "record_sha256": source.record_sha256,
        "workload": measurement.workload,
        "stage": measurement.stage,
        "cpu": None if measurement.cpu is None else measurement.cpu.to_json(),
        "wall": None if measurement.wall is None else measurement.wall.to_json(),
        "load_before": measurement.load_before,
        "load_after": measurement.load_after,
        "record_indicative": measurement.record_indicative,
        "indicative_reason": reason,
        "commit": provenance.commit,
        "dirty": provenance.dirty,
        "source_sha256": provenance.source_sha256,
        "native_sha256": provenance.native_sha256,
        "machine": provenance.machine,
        "recorded_at": provenance.recorded_at,
    }
    if schema == ENTRY_SCHEMA:
        values["adapter_version"] = source.adapter_version
    assert tuple(values) == _content_keys(schema)
    return values


def _digest(content: Mapping[str, Any]) -> str:
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
