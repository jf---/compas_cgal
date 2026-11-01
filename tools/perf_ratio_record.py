"""The ``perf-ratio/1`` record: one session of interleaved HEAD-vs-reference CPU runs.

A session pairs two arms (the pinned reference tree and the HEAD tree), each
identified the same way, and records per workload every pair of in-child
process-CPU runs with the 1-minute load around it. The headline per workload
is ``min(head) / min(ref)``; the median of per-pair ratios sits beside it as
the spread witness.

`PerfRatioSession.to_json` stores the derived values (min, median, ratio,
pair-ratio median, same_output) beside the runs for readers.
`PerfRatioSession.build` recomputes each one from the runs and refuses any
mismatch by exact equality: the same float operations on the same inputs give
the same bits, and JSON round-trips a float exactly, so no tolerance applies.
"""

from __future__ import annotations

import json
import math
import statistics
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any
from typing import NewType
from typing import Optional

from tools.perf_history_json import MalformedTimingRecordError
from tools.perf_history_json import field
from tools.perf_history_json import items
from tools.perf_history_json import number
from tools.perf_history_json import obj
from tools.perf_history_json import optional_text
from tools.perf_history_json import text
from tools.run_native_identity import NativeIdentity
from tools.run_provenance import TreeIdentity

PERF_RATIO_SCHEMA = "perf-ratio/1"

# Characters of the reference commit id in a stage name: this repository's ``git log --oneline``
# abbreviation, which git sizes automatically from the object count (8 here); a stable choice, not a git default.
REF_SHORT_CHARS = 8
# The stage name's prefix before the reference's short rev.
STAGE_PREFIX = "ratio-vs-"

# Process CPU of one arm's child, in seconds (in-child, never wall clock).
CpuSeconds = NewType("CpuSeconds", float)
# A dimensionless quotient of HEAD CPU over reference CPU.
Ratio = NewType("Ratio", float)
# The 1-minute load average, as ``os.getloadavg()[0]`` reports it.
LoadOneMinute = NewType("LoadOneMinute", float)


class MalformedRatioRecordError(Exception):
    """A ratio record has the wrong shape, schema or field type."""


class InconsistentRatioRecordError(Exception):
    """A stored derived value (min, median, ratio, pair-ratio median, same_output) differs from its recomputation,
    or ``ref_commit`` is not the ref arm's measured rev."""


class NonPositiveCpuError(Exception):
    """A run's CPU seconds are not positive, so no ratio over it is defined."""


class ArmName(str, Enum):
    """The two arms of a ratio session."""

    REF = "ref"
    HEAD = "head"


# The two orders a pair may run in: each arm exactly once.
_PAIR_ORDERS = frozenset({(ArmName.REF, ArmName.HEAD), (ArmName.HEAD, ArmName.REF)})


@dataclass(frozen=True)
class ArmIdentity:
    """Which tree, interpreter and native build one arm measured.

    Attributes:
        rev: Full commit id of the arm's tree.
        dirty: Whether the arm's working tree differed from *rev*.
        patch_sha256: Digest of that difference, None when clean.
        native_sha256: Digest over every native extension the arm's child loaded.
        native_sources_sha256: Digest of the arm tree's C++ sources.
        environment: The arm's pixi environment name.
        python_version: The arm interpreter's version.
        stock_2_sha256: SHA-256 of the arm's ``_stock_2`` extension.
        pixi_lock_sha256: SHA-256 of the arm tree's ``pixi.lock``.
        compas_cgal_file: Where the arm's child imported ``compas_cgal`` from.
        stock_2_origin: Where the arm's child loaded ``_stock_2`` from.
        editable_build_dir: The arm's editable native build directory.
        interpreter_prefix: ``sys.prefix`` of the arm's interpreter.
    """

    rev: str
    dirty: bool
    patch_sha256: Optional[str]
    native_sha256: str
    native_sources_sha256: str
    environment: str
    python_version: str
    stock_2_sha256: str
    pixi_lock_sha256: str
    compas_cgal_file: str
    stock_2_origin: str
    editable_build_dir: str
    interpreter_prefix: str

    @classmethod
    def of(
        cls,
        *,
        tree: TreeIdentity,
        native: NativeIdentity,
        interpreter: tuple[str, str, str],
        pixi_lock_sha256: str,
        compas_cgal_file: str,
        stock_2_origin: str,
        editable_build_dir: str,
        interpreter_prefix: str,
    ) -> ArmIdentity:
        """An arm identity from the shared identity readers.

        Args:
            tree: `tree_identity` of the arm's tree.
            native: The arm child's `NativeIdentity`.
            interpreter: ``(environment, python_version, stock_2_sha256)`` from `interpreter_identity_of`.
            pixi_lock_sha256: SHA-256 of the arm tree's ``pixi.lock``.
            compas_cgal_file: Where the arm's child imported ``compas_cgal`` from.
            stock_2_origin: Where the arm's child loaded ``_stock_2`` from.
            editable_build_dir: The arm's editable native build directory.
            interpreter_prefix: ``sys.prefix`` of the arm's interpreter.

        Returns:
            The arm identity.
        """
        environment, python_version, stock_2_sha256 = interpreter
        return cls(
            rev=tree.rev,
            dirty=tree.dirty,
            patch_sha256=tree.patch_sha256,
            native_sha256=native.native_sha256,
            native_sources_sha256=native.native_sources_sha256,
            environment=environment,
            python_version=python_version,
            stock_2_sha256=stock_2_sha256,
            pixi_lock_sha256=pixi_lock_sha256,
            compas_cgal_file=compas_cgal_file,
            stock_2_origin=stock_2_origin,
            editable_build_dir=editable_build_dir,
            interpreter_prefix=interpreter_prefix,
        )

    def to_json(self) -> dict[str, object]:
        """The arm identity as a JSON object."""
        record: dict[str, object] = {name: getattr(self, name) for name in _ARM_TEXT_FIELDS}
        record["dirty"] = self.dirty
        record["patch_sha256"] = self.patch_sha256
        return record

    @classmethod
    def build(cls, raw: object, where: str = "arm") -> ArmIdentity:
        """A checked arm identity from its parsed JSON object (`to_json`'s shape).

        Raises:
            MalformedRatioRecordError: Wrong keys or field type, or a patch digest present exactly when the tree is clean.
        """
        try:
            return _build_arm(raw, where)
        except MalformedTimingRecordError as error:
            raise MalformedRatioRecordError(str(error)) from error


_ARM_TEXT_FIELDS = (
    "rev",
    "native_sha256",
    "native_sources_sha256",
    "environment",
    "python_version",
    "stock_2_sha256",
    "pixi_lock_sha256",
    "compas_cgal_file",
    "stock_2_origin",
    "editable_build_dir",
    "interpreter_prefix",
)


@dataclass(frozen=True)
class PairRun:
    """One interleaved pair: each arm's child run once, in *order*. Raw construction is unchecked; `build` validates.

    Attributes:
        order: Which arm ran first.
        ref_cpu: The reference child's process CPU.
        head_cpu: The HEAD child's process CPU.
        load_before: 1-minute load before the pair.
        load_after: 1-minute load after the pair.
    """

    order: tuple[ArmName, ArmName]
    ref_cpu: CpuSeconds
    head_cpu: CpuSeconds
    load_before: LoadOneMinute
    load_after: LoadOneMinute

    @classmethod
    def build(
        cls,
        *,
        order: tuple[ArmName, ArmName],
        ref_cpu: CpuSeconds,
        head_cpu: CpuSeconds,
        load_before: LoadOneMinute,
        load_after: LoadOneMinute,
    ) -> PairRun:
        """A checked pair.

        Raises:
            NonPositiveCpuError: Either arm's CPU is not positive (NaN included).
            MalformedRatioRecordError: *order* is not each arm once, or a CPU or load is not finite.
        """
        if tuple(order) not in _PAIR_ORDERS:
            raise MalformedRatioRecordError(f"pair order {order!r} does not run each arm exactly once.")
        for name, cpu in ((ArmName.REF, ref_cpu), (ArmName.HEAD, head_cpu)):
            if not cpu > 0.0:
                raise NonPositiveCpuError(f"{name.value} run reports {cpu!r} CPU seconds; a ratio needs positive CPU.")
        for label, value in (("ref_cpu", ref_cpu), ("head_cpu", head_cpu), ("load_before", load_before), ("load_after", load_after)):
            if not math.isfinite(value):
                raise MalformedRatioRecordError(f"pair {label} is finite, not {value!r}.")
        return cls(
            order=(ArmName(order[0]), ArmName(order[1])),
            ref_cpu=CpuSeconds(float(ref_cpu)),
            head_cpu=CpuSeconds(float(head_cpu)),
            load_before=LoadOneMinute(float(load_before)),
            load_after=LoadOneMinute(float(load_after)),
        )

    @property
    def ratio(self) -> Ratio:
        """HEAD CPU over reference CPU for this pair."""
        return Ratio(self.head_cpu / self.ref_cpu)

    def to_json(self) -> dict[str, object]:
        """The pair as a JSON object."""
        return {
            "order": [arm.value for arm in self.order],
            "ref_cpu": self.ref_cpu,
            "head_cpu": self.head_cpu,
            "load_before": self.load_before,
            "load_after": self.load_after,
        }


@dataclass(frozen=True)
class WorkloadRatio:
    """All pairs of one workload, and the two arms' output digests.

    Attributes:
        workload: The workload's name.
        family: The workload's family (which clock owns its number).
        pairs: The interleaved pairs, in run order.
        ref_digest: Digest of the reference arm's output.
        head_digest: Digest of the HEAD arm's output.
    """

    workload: str
    family: str
    pairs: tuple[PairRun, ...]
    ref_digest: str
    head_digest: str

    def __post_init__(self) -> None:
        """Refuse a workload without pairs, whose min, median and ratio are undefined.

        Raises:
            MalformedRatioRecordError: *pairs* is empty.
        """
        if not self.pairs:
            raise MalformedRatioRecordError(f"workload {self.workload!r} holds at least one pair.")

    @property
    def ref_min(self) -> CpuSeconds:
        """The fastest reference run."""
        return CpuSeconds(min(pair.ref_cpu for pair in self.pairs))

    @property
    def ref_median(self) -> CpuSeconds:
        """The median reference run."""
        return CpuSeconds(statistics.median(pair.ref_cpu for pair in self.pairs))

    @property
    def head_min(self) -> CpuSeconds:
        """The fastest HEAD run."""
        return CpuSeconds(min(pair.head_cpu for pair in self.pairs))

    @property
    def head_median(self) -> CpuSeconds:
        """The median HEAD run."""
        return CpuSeconds(statistics.median(pair.head_cpu for pair in self.pairs))

    @property
    def ratio(self) -> Ratio:
        """The headline: ``head_min / ref_min``."""
        return Ratio(self.head_min / self.ref_min)

    @property
    def pair_ratio_median(self) -> Ratio:
        """The median of per-pair ratios, the spread witness beside `ratio`."""
        return Ratio(statistics.median(pair.ratio for pair in self.pairs))

    @property
    def same_output(self) -> bool:
        """Whether both arms produced the same output digest (flagged, never refused)."""
        return self.ref_digest == self.head_digest

    @property
    def load_range(self) -> tuple[LoadOneMinute, LoadOneMinute]:
        """(lowest, highest) 1-minute load seen around any pair."""
        loads = [load for pair in self.pairs for load in (pair.load_before, pair.load_after)]
        return LoadOneMinute(min(loads)), LoadOneMinute(max(loads))

    @classmethod
    def build(cls, raw: object, where: str = "workload") -> WorkloadRatio:
        """A checked workload from its parsed JSON object (`to_json`'s shape).

        Raises:
            MalformedRatioRecordError: Wrong keys or field type, or no pairs.
            InconsistentRatioRecordError: A stored derived value differs from its recomputation from the pairs.
            NonPositiveCpuError: A run's CPU is not positive.
        """
        try:
            return _build_workload(raw, where)
        except MalformedTimingRecordError as error:
            raise MalformedRatioRecordError(str(error)) from error

    def to_json(self) -> dict[str, object]:
        """The workload as a JSON object, derived values stored beside the runs."""
        return {
            "workload": self.workload,
            "family": self.family,
            "pairs": [pair.to_json() for pair in self.pairs],
            "ref_digest": self.ref_digest,
            "head_digest": self.head_digest,
            "ref_min": self.ref_min,
            "ref_median": self.ref_median,
            "head_min": self.head_min,
            "head_median": self.head_median,
            "ratio": self.ratio,
            "pair_ratio_median": self.pair_ratio_median,
            "same_output": self.same_output,
        }


@dataclass(frozen=True)
class Injection:
    """A deliberate CPU spin added inside one arm's measured span (the audit).

    Attributes:
        arm: The arm whose child spun.
        seconds: Process CPU seconds spun.
    """

    arm: ArmName
    seconds: CpuSeconds

    def to_json(self) -> dict[str, object]:
        """The injection as a JSON object."""
        return {"arm": self.arm.value, "seconds": self.seconds}


@dataclass(frozen=True)
class PerfRatioSession:
    """One ``perf-ratio/1`` session. Raw construction is unchecked; `build` validates a parsed record.

    Attributes:
        ref_commit: Full commit id of the pinned reference.
        arms: Identity of each arm, keyed by `ArmName`; both arms present.
        workloads: Per-workload pairs, in run order; names unique.
        host: Host name of the machine.
        cpu_model: CPU model of the machine.
        started_utc: Session start, ISO 8601 with a UTC offset.
        ended_utc: Session end, ISO 8601 with a UTC offset.
        harness_sha256: Digest of the harness sources.
        injected: The audit's injected CPU, None for an honest session.
    """

    ref_commit: str
    arms: Mapping[ArmName, ArmIdentity]
    workloads: tuple[WorkloadRatio, ...]
    host: str
    cpu_model: str
    started_utc: str
    ended_utc: str
    harness_sha256: str
    injected: Optional[Injection]

    @property
    def stage(self) -> str:
        """The history stage name: `stage_of` the reference commit."""
        return stage_of(self.ref_commit)

    @property
    def wall_seconds(self) -> float:
        """Session duration in wall seconds (reporting only, never a measurement)."""
        return (_parse_utc(self.ended_utc, "ended_utc") - _parse_utc(self.started_utc, "started_utc")).total_seconds()

    def to_json(self) -> dict[str, object]:
        """The session as a JSON object."""
        return {
            "schema": PERF_RATIO_SCHEMA,
            "ref_commit": self.ref_commit,
            "arms": arms_to_json(self.arms),
            "workloads": [workload.to_json() for workload in self.workloads],
            "host": self.host,
            "cpu_model": self.cpu_model,
            "started_utc": self.started_utc,
            "ended_utc": self.ended_utc,
            "harness_sha256": self.harness_sha256,
            "injected": None if self.injected is None else self.injected.to_json(),
        }

    def dumps(self) -> str:
        """The session as deterministic JSON text (sorted keys, two-space indent, trailing newline)."""
        return json.dumps(self.to_json(), indent=2, sort_keys=True) + "\n"

    @classmethod
    def build(cls, raw: object) -> PerfRatioSession:
        """A checked session from a parsed ``perf-ratio/1`` record.

        Args:
            raw: The parsed JSON record.

        Returns:
            The session; ``build(json.loads(s.dumps())) == s``.

        Raises:
            MalformedRatioRecordError: Wrong shape, schema or field type.
            InconsistentRatioRecordError: A stored derived value differs from its recomputation, or
                ``ref_commit`` is not the ref arm's rev.
            NonPositiveCpuError: A run's CPU is not positive.
        """
        try:
            return _build_session(raw)
        except MalformedTimingRecordError as error:
            raise MalformedRatioRecordError(str(error)) from error


def stage_of(ref_commit: str) -> str:
    """The history stage name of a reference commit: ``ratio-vs-<first REF_SHORT_CHARS characters>``."""
    return f"{STAGE_PREFIX}{ref_commit[:REF_SHORT_CHARS]}"


def arms_to_json(arms: Mapping[ArmName, ArmIdentity]) -> dict[str, object]:
    """Both arms' identities as one JSON object keyed by arm name."""
    return {arm.value: arms[arm].to_json() for arm in ArmName}


def build_arms(raw: object, where: str = "arms") -> dict[ArmName, ArmIdentity]:
    """Both arms' checked identities from their parsed JSON object (`arms_to_json`'s shape).

    Raises:
        MalformedRatioRecordError: Not an object of exactly both arms, or an arm is malformed.
    """
    try:
        return _build_arms(raw, where)
    except MalformedTimingRecordError as error:
        raise MalformedRatioRecordError(str(error)) from error


def load_session(path: Path) -> PerfRatioSession:
    """The checked session stored at *path*.

    Raises:
        MalformedRatioRecordError: The file is not JSON, or not a ``perf-ratio/1`` record.
        InconsistentRatioRecordError: A stored derived value differs from its recomputation.
        NonPositiveCpuError: A run's CPU is not positive.
    """
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise MalformedRatioRecordError(f"{path} is not JSON: {error}") from error
    return PerfRatioSession.build(raw)


def _parse_utc(value: str, where: str) -> datetime:
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise MalformedRatioRecordError(f"{where} {value!r} is not ISO 8601.") from error
    if parsed.utcoffset() is None:
        raise MalformedRatioRecordError(f"{where} {value!r} carries no UTC offset.")
    return parsed


def _json_number(value: object, where: str) -> float:
    # The shared accessor also accepts decimal strings; this record writes only JSON numbers.
    if isinstance(value, str):
        raise MalformedRatioRecordError(f"{where} is a JSON number, not the string {value!r}.")
    return number(value, where)


def _boolean(value: object, where: str) -> bool:
    if not isinstance(value, bool):
        raise MalformedRatioRecordError(f"{where} is a boolean, not {value!r}.")
    return value


def _arm_name(value: object, where: str) -> ArmName:
    name = text(value, where)
    try:
        return ArmName(name)
    except ValueError as error:
        raise MalformedRatioRecordError(f"{where} names an arm ({', '.join(a.value for a in ArmName)}), not {name!r}.") from error


def _require_keys(record: Mapping[str, Any], expected: frozenset[str], where: str) -> None:
    keys = frozenset(record)
    if keys != expected:
        missing = sorted(expected - keys)
        extra = sorted(keys - expected)
        raise MalformedRatioRecordError(f"{where} keys differ: missing {missing}, unexpected {extra}.")


_ARM_KEYS = frozenset(_ARM_TEXT_FIELDS) | {"dirty", "patch_sha256"}
_PAIR_KEYS = frozenset({"order", "ref_cpu", "head_cpu", "load_before", "load_after"})
_WORKLOAD_RUN_KEYS = frozenset({"workload", "family", "pairs", "ref_digest", "head_digest"})
_WORKLOAD_DERIVED_KEYS = ("ref_min", "ref_median", "head_min", "head_median", "ratio", "pair_ratio_median", "same_output")
_WORKLOAD_KEYS = _WORKLOAD_RUN_KEYS | frozenset(_WORKLOAD_DERIVED_KEYS)
_SESSION_KEYS = frozenset({"schema", "ref_commit", "arms", "workloads", "host", "cpu_model", "started_utc", "ended_utc", "harness_sha256", "injected"})
_INJECTION_KEYS = frozenset({"arm", "seconds"})


def _build_arm(raw: object, where: str) -> ArmIdentity:
    record = obj(raw, where)
    _require_keys(record, _ARM_KEYS, where)
    values = {name: text(field(record, name, where), f"{where}.{name}") for name in _ARM_TEXT_FIELDS}
    dirty = _boolean(field(record, "dirty", where), f"{where}.dirty")
    patch_sha256 = optional_text(record, "patch_sha256", where)
    if dirty != (patch_sha256 is not None):
        raise MalformedRatioRecordError(f"{where}: a patch digest is present exactly when the tree is dirty (dirty={dirty}, patch_sha256={patch_sha256!r}).")
    return ArmIdentity(dirty=dirty, patch_sha256=patch_sha256, **values)


def _build_arms(raw: object, where: str) -> dict[ArmName, ArmIdentity]:
    record = obj(raw, where)
    _require_keys(record, frozenset(arm.value for arm in ArmName), where)
    return {arm: _build_arm(record[arm.value], f"{where}.{arm.value}") for arm in ArmName}


def _build_pair(raw: object, where: str) -> PairRun:
    record = obj(raw, where)
    _require_keys(record, _PAIR_KEYS, where)
    order_items = items(field(record, "order", where), f"{where}.order")
    if len(order_items) != len(ArmName):
        raise MalformedRatioRecordError(f"{where}.order names {len(ArmName)} arms, not {len(order_items)}.")
    first, second = (_arm_name(item, f"{where}.order[{index}]") for index, item in enumerate(order_items))
    return PairRun.build(
        order=(first, second),
        ref_cpu=CpuSeconds(_json_number(field(record, "ref_cpu", where), f"{where}.ref_cpu")),
        head_cpu=CpuSeconds(_json_number(field(record, "head_cpu", where), f"{where}.head_cpu")),
        load_before=LoadOneMinute(_json_number(field(record, "load_before", where), f"{where}.load_before")),
        load_after=LoadOneMinute(_json_number(field(record, "load_after", where), f"{where}.load_after")),
    )


def _build_workload(raw: object, where: str) -> WorkloadRatio:
    record = obj(raw, where)
    _require_keys(record, _WORKLOAD_KEYS, where)
    pair_items = items(field(record, "pairs", where), f"{where}.pairs")
    workload = WorkloadRatio(
        workload=text(field(record, "workload", where), f"{where}.workload"),
        family=text(field(record, "family", where), f"{where}.family"),
        pairs=tuple(_build_pair(item, f"{where}.pairs[{index}]") for index, item in enumerate(pair_items)),
        ref_digest=text(field(record, "ref_digest", where), f"{where}.ref_digest"),
        head_digest=text(field(record, "head_digest", where), f"{where}.head_digest"),
    )
    for name in _WORKLOAD_DERIVED_KEYS:
        stored_raw = field(record, name, where)
        recomputed = getattr(workload, name)
        stored: object = _boolean(stored_raw, f"{where}.{name}") if isinstance(recomputed, bool) else _json_number(stored_raw, f"{where}.{name}")
        if stored != recomputed:
            raise InconsistentRatioRecordError(f"{where}.{name} is stored as {stored!r} but the runs give {recomputed!r}.")
    return workload


def _build_injection(record: Mapping[str, Any], where: str) -> Optional[Injection]:
    raw = record.get("injected")
    if raw is None:
        return None
    injected = obj(raw, f"{where}.injected")
    _require_keys(injected, _INJECTION_KEYS, f"{where}.injected")
    seconds = _json_number(field(injected, "seconds", f"{where}.injected"), f"{where}.injected.seconds")
    if not seconds > 0.0:
        raise MalformedRatioRecordError(f"{where}.injected.seconds is positive, not {seconds!r}.")
    return Injection(arm=_arm_name(field(injected, "arm", f"{where}.injected"), f"{where}.injected.arm"), seconds=CpuSeconds(seconds))


def _build_session(raw: object) -> PerfRatioSession:
    where = "perf-ratio session"
    record = obj(raw, where)
    schema = text(field(record, "schema", where), f"{where}.schema")
    if schema != PERF_RATIO_SCHEMA:
        raise MalformedRatioRecordError(f"{where}.schema is {PERF_RATIO_SCHEMA!r}, not {schema!r}.")
    _require_keys(record, _SESSION_KEYS, where)

    arms = _build_arms(field(record, "arms", where), f"{where}.arms")

    workload_items = items(field(record, "workloads", where), f"{where}.workloads")
    if not workload_items:
        raise MalformedRatioRecordError(f"{where}.workloads holds at least one workload.")
    workloads = tuple(_build_workload(item, f"{where}.workloads[{index}]") for index, item in enumerate(workload_items))
    names = [workload.workload for workload in workloads]
    if len(set(names)) != len(names):
        raise MalformedRatioRecordError(f"{where}.workloads names each workload once, not {names}.")

    started_utc = text(field(record, "started_utc", where), f"{where}.started_utc")
    ended_utc = text(field(record, "ended_utc", where), f"{where}.ended_utc")
    _parse_utc(started_utc, f"{where}.started_utc")
    _parse_utc(ended_utc, f"{where}.ended_utc")

    ref_commit = text(field(record, "ref_commit", where), f"{where}.ref_commit")
    if ref_commit != arms[ArmName.REF].rev:
        raise InconsistentRatioRecordError(f"{where}.ref_commit {ref_commit!r} is not the ref arm's measured rev {arms[ArmName.REF].rev!r}.")

    return PerfRatioSession(
        ref_commit=ref_commit,
        arms=arms,
        workloads=workloads,
        host=text(field(record, "host", where), f"{where}.host"),
        cpu_model=text(field(record, "cpu_model", where), f"{where}.cpu_model"),
        started_utc=started_utc,
        ended_utc=ended_utc,
        harness_sha256=text(field(record, "harness_sha256", where), f"{where}.harness_sha256"),
        injected=_build_injection(record, where),
    )
