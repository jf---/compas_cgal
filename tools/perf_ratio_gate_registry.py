"""The ratio gate registry: ``benchmarks/results/perf_ratio_gate.json``.

The registry holds the reference commit the gate decides against, the harness
digest every gated session must carry, the gated workload set, and per CPU model
one A/A calibration: the margins ``m_w`` and the cited ``perf-ratio/1`` records
they come from (spec §Margins: ``m_w = SAFETY * max_k |rho_{w,k} - 1|``).

The registry is content-addressed. Its ``registry_sha256`` is the SHA-256 of its
canonical JSON content, and each cited record carries its file's SHA-256. On load
nothing stored is trusted: the self-digest is re-derived, every cited record is
re-hashed, re-read and re-validated as an honest A/A session of this registry,
and every margin is re-derived from those records and compared with ``==``
(the same float operations on the same inputs give the same bits, and JSON
round-trips a float exactly, so no tolerance applies). `write_registry` runs
the same verification before it writes, so it never writes a registry that
would not load.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import NewType
from typing import Optional

from tools.perf_history_json import MalformedTimingRecordError
from tools.perf_history_json import field
from tools.perf_history_json import items
from tools.perf_history_json import number
from tools.perf_history_json import obj
from tools.perf_history_json import text
from tools.perf_history_store import REPOSITORY
from tools.perf_ratio_record import REF_SHORT_CHARS
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session
from tools.run_provenance import sha256_file

REGISTRY_SCHEMA = "perf-ratio-gate/1"
REGISTRY_PATH = REPOSITORY / "benchmarks" / "results" / "perf_ratio_gate.json"
CALIBRATION_ROOT = REPOSITORY / "benchmarks" / "results" / "perf_ratio_calibration"
# Dimensionless factor on the worst observed A/A deviation |rho - 1| (spec §Margins): it doubles it.
SAFETY = 2.0
# K, the number of A/A sessions a calibration cites at least (spec §Margins: K = 5).
CALIBRATION_SESSIONS = 5

# A dimensionless allowance on a workload's ratio: the gate passes ``rho_w <= 1 + m_w``.
Margin = NewType("Margin", float)

# Runs of characters a CPU model slug replaces with one hyphen.
_SLUG_SEPARATORS = re.compile(r"[^a-z0-9]+")


class MalformedRegistryError(Exception):
    """The registry file is not JSON, or has the wrong schema, keys or field types."""


class RegistryDigestMismatchError(Exception):
    """The stored ``registry_sha256`` is not the digest of the registry's content."""


class SafetyMismatchError(Exception):
    """The registry's ``safety`` is not this module's `SAFETY`."""


class CalibrationRecordDigestError(Exception):
    """A cited calibration record's bytes no longer have the cited ``record_sha256``."""


class CalibrationRecordMissingError(Exception):
    """A cited calibration record does not exist."""


class CalibrationRecordOutsideRootError(Exception):
    """A calibration record lies outside the root its cited path is relative to."""


class MarginMismatchError(Exception):
    """A stored margin (or margin set) differs from the one its cited records re-derive."""


class InvalidCalibrationSessionError(Exception):
    """A calibration session is not an honest A/A session of this registry and machine."""


class TooFewCalibrationSessionsError(Exception):
    """A calibration cites fewer than `CALIBRATION_SESSIONS` sessions."""


class DuplicateCalibrationSessionError(Exception):
    """A calibration cites one record (by path or by content) more than once."""


class DuplicateCalibrationError(Exception):
    """A registry holds two calibrations for one ``cpu_model``."""


class CalibrationWorkloadMismatchError(Exception):
    """A calibration's workload set is not the registry's."""


class MissingWorkloadsError(Exception):
    """A registry's workload set is empty (a calibration-less registry is built without one)."""


class DuplicateRegistryWorkloadError(Exception):
    """A registry's workload set names one workload more than once."""


class UnsortedRegistryWorkloadsError(Exception):
    """A registry's workload set is not in sorted order."""


class UncalibratedMachineError(Exception):
    """The registry holds no calibration for a ``cpu_model``."""


class UnsluggableCpuModelError(Exception):
    """A CPU model has no letter or digit, so it has no calibration directory name."""


@dataclass(frozen=True)
class CitedRecord:
    """One calibration session record, cited by path and content.

    Attributes:
        path: The record's path relative to the registry root (`REPOSITORY`), POSIX separators.
        record_sha256: SHA-256 of the record file's bytes.
    """

    path: str
    record_sha256: str

    def to_json(self) -> dict[str, object]:
        """The citation as a JSON object."""
        return {"path": self.path, "record_sha256": self.record_sha256}


@dataclass(frozen=True)
class Calibration:
    """One machine's A/A calibration. Raw construction is unchecked; `calibration_of` builds a verified one.

    Attributes:
        cpu_model: The CPU model every cited session ran on.
        sessions: The cited A/A session records, in citation order.
        margins: Workload -> its margin, as `derive_margins` gives it from the cited sessions.
    """

    cpu_model: str
    sessions: tuple[CitedRecord, ...]
    margins: Mapping[str, Margin]

    def to_json(self) -> dict[str, object]:
        """The calibration as a JSON object."""
        return {
            "cpu_model": self.cpu_model,
            "sessions": [record.to_json() for record in self.sessions],
            "margins": dict(self.margins),
        }


@dataclass(frozen=True)
class GateRegistry:
    """The gate's reference, harness, workload set and per-machine margins.

    Raw construction is unchecked; `build_registry` and `load_registry` build checked ones.

    Attributes:
        ref_commit: Full commit id of the reference the gate decides against.
        harness_sha256: The harness digest every calibration and gated session carries.
        safety: The margin factor the margins were derived with (`SAFETY`).
        workloads: The gated workload names, sorted.
        calibrations: One calibration per CPU model, sorted by ``cpu_model``.
    """

    ref_commit: str
    harness_sha256: str
    safety: float
    workloads: tuple[str, ...]
    calibrations: tuple[Calibration, ...]

    def margins_for(self, cpu_model: str) -> Mapping[str, Margin]:
        """The margins calibrated on *cpu_model*.

        Raises:
            UncalibratedMachineError: No calibration names *cpu_model*.
        """
        for calibration in self.calibrations:
            if calibration.cpu_model == cpu_model:
                return calibration.margins
        calibrated = ", ".join(repr(calibration.cpu_model) for calibration in self.calibrations) or "none"
        raise UncalibratedMachineError(f"no A/A calibration for cpu_model {cpu_model!r} (calibrated: {calibrated}); run an A/A calibration on that machine first.")

    def content(self) -> dict[str, object]:
        """The registry's JSON content, without its digest."""
        return {
            "schema": REGISTRY_SCHEMA,
            "ref_commit": self.ref_commit,
            "harness_sha256": self.harness_sha256,
            "safety": self.safety,
            "workloads": list(self.workloads),
            "calibrations": [calibration.to_json() for calibration in self.calibrations],
        }

    @property
    def registry_sha256(self) -> str:
        """SHA-256 of the canonical JSON content (sorted keys, ``(",", ":")`` separators)."""
        return _digest(self.content())

    def dumps(self) -> str:
        """The registry as deterministic JSON text (sorted keys, two-space indent, trailing newline), digest included."""
        document = self.content()
        document["registry_sha256"] = self.registry_sha256
        return json.dumps(document, indent=2, sort_keys=True) + "\n"


def calibrated_slug(cpu_model: str) -> str:
    """The directory name of *cpu_model*: lowercase, runs of non-``[a-z0-9]`` as one ``-``, stripped of ``-``.

    Raises:
        UnsluggableCpuModelError: *cpu_model* has no ASCII letter or digit.
    """
    slug = _SLUG_SEPARATORS.sub("-", cpu_model.lower()).strip("-")
    if not slug:
        raise UnsluggableCpuModelError(f"cpu_model {cpu_model!r} has no letter or digit to name a calibration directory.")
    return slug


def calibration_dir(ref_commit: str, cpu_model: str) -> Path:
    """Where the A/A records of *ref_commit* on *cpu_model* live: ``CALIBRATION_ROOT/<ref short>/<cpu slug>``."""
    return CALIBRATION_ROOT / ref_commit[:REF_SHORT_CHARS] / calibrated_slug(cpu_model)


def derive_margins(sessions: Sequence[PerfRatioSession], safety: float = SAFETY) -> dict[str, Margin]:
    """Per workload ``m_w = safety * max_k |rho_{w,k} - 1|`` over *sessions*.

    Args:
        sessions: The A/A sessions; all over one workload set.
        safety: The factor on the worst deviation.

    Returns:
        Workload -> its margin.

    Raises:
        TooFewCalibrationSessionsError: *sessions* is empty.
        InvalidCalibrationSessionError: The sessions are not all over one workload set.
    """
    if not sessions:
        raise TooFewCalibrationSessionsError("margins are derived from at least one A/A session, not none.")
    workloads = _workload_set(sessions[0])
    for index, session in enumerate(sessions):
        if _workload_set(session) != workloads:
            raise InvalidCalibrationSessionError(f"A/A session {index} measures {sorted(_workload_set(session))}, not the first session's {sorted(workloads)}.")
    ratios = [{workload.workload: workload.ratio for workload in session.workloads} for session in sessions]
    return {name: Margin(safety * max(abs(by_name[name] - 1.0) for by_name in ratios)) for name in sorted(workloads)}


def calibration_of(record_paths: Sequence[Path], *, ref_commit: str, harness_sha256: str, root: Path = REPOSITORY) -> Calibration:
    """The calibration the A/A records at *record_paths* give, every record verified.

    The calibration's ``cpu_model`` and workload set are the first record's; every
    other record must match them.

    Args:
        record_paths: The ``perf-ratio/1`` A/A session records, at least `CALIBRATION_SESSIONS`.
        ref_commit: The reference commit both arms of every session ran at.
        harness_sha256: The harness digest every session carries.
        root: The directory cited paths are relative to (`REPOSITORY`).

    Returns:
        The calibration, citing each record by its path under *root* and its SHA-256.

    Raises:
        TooFewCalibrationSessionsError: Fewer than `CALIBRATION_SESSIONS` records.
        CalibrationRecordOutsideRootError: A record lies outside *root*.
        DuplicateCalibrationSessionError: One record (by path or content) is given twice.
        InvalidCalibrationSessionError: A session is not an honest A/A session at *ref_commit* under
            *harness_sha256* on the first session's machine and workload set.
    """
    _require_session_count(len(record_paths), "the calibration")
    resolved_root = root.resolve()
    cited = []
    for path in record_paths:
        resolved = path.resolve()
        if not resolved.is_relative_to(resolved_root):
            raise CalibrationRecordOutsideRootError(f"calibration record {path} lies outside {root}; a registry cites records under its root only.")
        cited.append(CitedRecord(path=resolved.relative_to(resolved_root).as_posix(), record_sha256=sha256_file(resolved)))
    _require_distinct(cited, "the calibration")
    sessions = [load_session(path) for path in record_paths]
    cpu_model = sessions[0].cpu_model
    workloads = _workload_set(sessions[0])
    for record, session in zip(cited, sessions):
        _require_aa_session(session, ref_commit=ref_commit, harness_sha256=harness_sha256, cpu_model=cpu_model, workloads=workloads, where=record.path)
    return Calibration(cpu_model=cpu_model, sessions=tuple(cited), margins=derive_margins(sessions))


def build_registry(
    *,
    ref_commit: str,
    harness_sha256: str,
    calibrations: Sequence[Calibration],
    workloads: Optional[Sequence[str]] = None,
) -> GateRegistry:
    """A registry over *calibrations*, its workload set theirs.

    Args:
        ref_commit: Full commit id of the reference.
        harness_sha256: The harness digest.
        calibrations: At most one calibration per CPU model; may be empty (a freshly promoted reference).
        workloads: The gated workload names; required when *calibrations* is empty, and when given
            beside calibrations it must be their workload set.

    Returns:
        The registry, workloads sorted and calibrations sorted by ``cpu_model``, `SAFETY` as its safety.

    Raises:
        DuplicateCalibrationError: Two calibrations name one ``cpu_model``.
        MissingWorkloadsError: No calibrations and no (or an empty) workload set.
        DuplicateRegistryWorkloadError: *workloads* names one workload twice.
        CalibrationWorkloadMismatchError: A calibration's workload set is not the registry's.
    """
    _require_one_calibration_per_cpu(calibrations)
    if workloads is not None:
        names = tuple(sorted(workloads))
    elif calibrations:
        names = tuple(sorted(calibrations[0].margins))
    else:
        names = ()
    _require_workload_shape(names)
    for calibration in calibrations:
        _require_calibration_workloads(calibration, names)
    return GateRegistry(
        ref_commit=ref_commit,
        harness_sha256=harness_sha256,
        safety=SAFETY,
        workloads=names,
        calibrations=tuple(sorted(calibrations, key=lambda calibration: calibration.cpu_model)),
    )


def load_registry(path: Path = REGISTRY_PATH, *, root: Path = REPOSITORY) -> GateRegistry:
    """The registry at *path*, every stored value re-derived.

    Checks, in order: schema and shape; the self-digest; safety; the workload set (non-empty, sorted,
    no repeats); one calibration per CPU; then per calibration: K, each cited record's digest, each session's validity, the margin re-derivation.

    Args:
        path: The registry file.
        root: The directory cited record paths are relative to (`REPOSITORY`).

    Returns:
        The verified registry.

    Raises:
        MalformedRegistryError: Not JSON, or the wrong schema, keys or field types.
        RegistryDigestMismatchError: The stored ``registry_sha256`` is not the content's.
        SafetyMismatchError: The stored safety is not `SAFETY`.
        MissingWorkloadsError: The workload set is empty.
        UnsortedRegistryWorkloadsError: The workload set is not sorted.
        DuplicateRegistryWorkloadError: The workload set names one workload twice.
        DuplicateCalibrationError: Two calibrations name one ``cpu_model``.
        TooFewCalibrationSessionsError: A calibration cites fewer than `CALIBRATION_SESSIONS` records.
        DuplicateCalibrationSessionError: A calibration cites one record twice.
        CalibrationRecordOutsideRootError: A cited path leaves *root*.
        CalibrationRecordMissingError: A cited record does not exist.
        CalibrationRecordDigestError: A cited record's bytes changed.
        InvalidCalibrationSessionError: A cited session is not an honest A/A session of this registry.
        MarginMismatchError: A stored margin set differs from the re-derived one.
    """
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise MalformedRegistryError(f"{path} is not JSON: {error}") from error
    try:
        registry, stored_digest = _parse(raw, str(path))
    except MalformedTimingRecordError as error:
        raise MalformedRegistryError(str(error)) from error
    if stored_digest != registry.registry_sha256:
        raise RegistryDigestMismatchError(f"{path}: stored registry_sha256 {stored_digest} is not its content's {registry.registry_sha256}.")
    _verify(registry, root)
    return registry


def write_registry(registry: GateRegistry, path: Path = REGISTRY_PATH, *, root: Path = REPOSITORY) -> None:
    """Write *registry* to *path* after verifying it as `load_registry` would.

    The text is written to a sibling file and moved over *path*, so *path* is
    never left half-written.

    Args:
        registry: The registry to write.
        path: The registry file.
        root: The directory cited record paths are relative to (`REPOSITORY`).

    Raises:
        SafetyMismatchError: As `load_registry`; *path* is then untouched.
        MissingWorkloadsError: As `load_registry`.
        UnsortedRegistryWorkloadsError: As `load_registry`.
        DuplicateRegistryWorkloadError: As `load_registry`.
        DuplicateCalibrationError: As `load_registry`.
        TooFewCalibrationSessionsError: As `load_registry`.
        DuplicateCalibrationSessionError: As `load_registry`.
        CalibrationRecordOutsideRootError: As `load_registry`.
        CalibrationRecordMissingError: As `load_registry`.
        CalibrationRecordDigestError: As `load_registry`.
        InvalidCalibrationSessionError: As `load_registry`.
        MarginMismatchError: As `load_registry`.
    """
    _verify(registry, root)
    staging = path.with_name(f"{path.name}.writing")
    staging.write_text(registry.dumps(), encoding="utf-8")
    staging.replace(path)


def _verify(registry: GateRegistry, root: Path) -> None:
    if registry.safety != SAFETY:
        raise SafetyMismatchError(f"registry safety is {registry.safety!r}, not this gate's SAFETY {SAFETY!r}.")
    _require_workload_shape(registry.workloads)
    _require_one_calibration_per_cpu(registry.calibrations)
    workloads = frozenset(registry.workloads)
    for calibration in registry.calibrations:
        where = f"calibration {calibration.cpu_model!r}"
        _require_session_count(len(calibration.sessions), where)
        _require_distinct(calibration.sessions, where)
        record_paths = [_cited_path(record, root) for record in calibration.sessions]
        for record, record_path in zip(calibration.sessions, record_paths):
            actual = sha256_file(record_path)
            if actual != record.record_sha256:
                raise CalibrationRecordDigestError(f"{where}: {record.path} now has SHA-256 {actual}, not the cited {record.record_sha256}; the record was rewritten.")
        sessions = [load_session(record_path) for record_path in record_paths]
        for record, session in zip(calibration.sessions, sessions):
            _require_aa_session(
                session,
                ref_commit=registry.ref_commit,
                harness_sha256=registry.harness_sha256,
                cpu_model=calibration.cpu_model,
                workloads=workloads,
                where=record.path,
            )
        derived = derive_margins(sessions, registry.safety)
        if dict(calibration.margins) != derived:
            raise MarginMismatchError(f"{where}: stored margins {dict(calibration.margins)} are not the ones its cited records give, {derived}.")


def _require_workload_shape(workloads: Sequence[str]) -> None:
    # The one shape check build, write and load share: non-empty, sorted, no repeats.
    if not workloads:
        raise MissingWorkloadsError("a registry names a non-empty workload set; a calibration-less registry is built with workloads.")
    if len(set(workloads)) != len(workloads):
        raise DuplicateRegistryWorkloadError(f"a registry names each workload once, not {list(workloads)}.")
    if list(workloads) != sorted(workloads):
        raise UnsortedRegistryWorkloadsError(f"a registry's workloads are sorted, not {list(workloads)}.")


def _require_session_count(count: int, where: str) -> None:
    if count < CALIBRATION_SESSIONS:
        raise TooFewCalibrationSessionsError(f"{where} cites {count} A/A sessions; a calibration needs at least {CALIBRATION_SESSIONS}.")


def _require_distinct(cited: Sequence[CitedRecord], where: str) -> None:
    for attribute in ("path", "record_sha256"):
        values = [getattr(record, attribute) for record in cited]
        repeated = sorted({value for value in values if values.count(value) > 1})
        if repeated:
            raise DuplicateCalibrationSessionError(f"{where} cites one record more than once (same {attribute}: {', '.join(repeated)}).")


def _require_one_calibration_per_cpu(calibrations: Sequence[Calibration]) -> None:
    models = [calibration.cpu_model for calibration in calibrations]
    repeated = sorted({model for model in models if models.count(model) > 1})
    if repeated:
        raise DuplicateCalibrationError(f"a registry holds one calibration per cpu_model; more than one for {', '.join(repr(model) for model in repeated)}.")


def _require_calibration_workloads(calibration: Calibration, workloads: Sequence[str]) -> None:
    if frozenset(calibration.margins) != frozenset(workloads):
        raise CalibrationWorkloadMismatchError(f"calibration {calibration.cpu_model!r} covers {sorted(calibration.margins)}, not the registry's workloads {list(workloads)}.")


def _require_aa_session(
    session: PerfRatioSession,
    *,
    ref_commit: str,
    harness_sha256: str,
    cpu_model: str,
    workloads: frozenset[str],
    where: str,
) -> None:
    if session.injected is not None:
        raise InvalidCalibrationSessionError(f"{where} is an audit session (injected {session.injected.seconds} s on {session.injected.arm.value}), not an A/A calibration.")
    for arm in ArmName:
        identity = session.arms[arm]
        if identity.rev != ref_commit:
            raise InvalidCalibrationSessionError(f"{where}: the {arm.value} arm ran at {identity.rev}, not the reference {ref_commit}; an A/A session runs both arms there.")
        if identity.dirty:
            raise InvalidCalibrationSessionError(f"{where}: the {arm.value} arm's tree was dirty (patch {identity.patch_sha256}).")
    if session.harness_sha256 != harness_sha256:
        raise InvalidCalibrationSessionError(f"{where} ran harness {session.harness_sha256}, not the registry's {harness_sha256}.")
    if session.cpu_model != cpu_model:
        raise InvalidCalibrationSessionError(f"{where} ran on {session.cpu_model!r}, not the calibration's {cpu_model!r}.")
    changed = [workload.workload for workload in session.workloads if not workload.same_output]
    if changed:
        raise InvalidCalibrationSessionError(f"{where}: the arms' outputs differ on {', '.join(changed)}; A/A arms must agree.")
    if _workload_set(session) != workloads:
        raise InvalidCalibrationSessionError(f"{where} measures {sorted(_workload_set(session))}, not the workloads {sorted(workloads)}.")


def _workload_set(session: PerfRatioSession) -> frozenset[str]:
    return frozenset(workload.workload for workload in session.workloads)


def _cited_path(record: CitedRecord, root: Path) -> Path:
    resolved_root = root.resolve()
    if Path(record.path).is_absolute():
        raise CalibrationRecordOutsideRootError(f"cited path {record.path} is absolute; a registry cites paths relative to its root.")
    resolved = (resolved_root / record.path).resolve()
    if not resolved.is_relative_to(resolved_root):
        raise CalibrationRecordOutsideRootError(f"cited path {record.path} leaves the root {root}.")
    if not resolved.is_file():
        raise CalibrationRecordMissingError(f"cited calibration record {record.path} does not exist under {root}.")
    return resolved


_REGISTRY_KEYS = frozenset({"schema", "ref_commit", "harness_sha256", "safety", "workloads", "calibrations", "registry_sha256"})
_CALIBRATION_KEYS = frozenset({"cpu_model", "sessions", "margins"})
_CITED_KEYS = frozenset({"path", "record_sha256"})


def _parse(raw: object, where: str) -> tuple[GateRegistry, str]:
    record = obj(raw, where)
    schema = text(field(record, "schema", where), f"{where}.schema")
    if schema != REGISTRY_SCHEMA:
        raise MalformedRegistryError(f"{where}.schema is {REGISTRY_SCHEMA!r}, not {schema!r}.")
    _require_keys(record, _REGISTRY_KEYS, where)
    workloads = tuple(text(item, f"{where}.workloads[{index}]") for index, item in enumerate(items(field(record, "workloads", where), f"{where}.workloads")))
    calibration_items = items(field(record, "calibrations", where), f"{where}.calibrations")
    registry = GateRegistry(
        ref_commit=text(field(record, "ref_commit", where), f"{where}.ref_commit"),
        harness_sha256=text(field(record, "harness_sha256", where), f"{where}.harness_sha256"),
        safety=_json_number(field(record, "safety", where), f"{where}.safety"),
        workloads=workloads,
        calibrations=tuple(_parse_calibration(item, f"{where}.calibrations[{index}]") for index, item in enumerate(calibration_items)),
    )
    return registry, text(field(record, "registry_sha256", where), f"{where}.registry_sha256")


def _parse_calibration(raw: object, where: str) -> Calibration:
    record = obj(raw, where)
    _require_keys(record, _CALIBRATION_KEYS, where)
    margins = obj(field(record, "margins", where), f"{where}.margins")
    session_items = items(field(record, "sessions", where), f"{where}.sessions")
    return Calibration(
        cpu_model=text(field(record, "cpu_model", where), f"{where}.cpu_model"),
        sessions=tuple(_parse_cited(item, f"{where}.sessions[{index}]") for index, item in enumerate(session_items)),
        margins={name: Margin(_json_number(value, f"{where}.margins.{name}")) for name, value in margins.items()},
    )


def _parse_cited(raw: object, where: str) -> CitedRecord:
    record = obj(raw, where)
    _require_keys(record, _CITED_KEYS, where)
    return CitedRecord(
        path=text(field(record, "path", where), f"{where}.path"),
        record_sha256=text(field(record, "record_sha256", where), f"{where}.record_sha256"),
    )


def _json_number(value: object, where: str) -> float:
    # The shared accessor also accepts decimal strings; this registry writes only JSON numbers.
    if isinstance(value, str):
        raise MalformedRegistryError(f"{where} is a JSON number, not the string {value!r}.")
    return number(value, where)


def _require_keys(record: Mapping[str, Any], expected: frozenset[str], where: str) -> None:
    keys = frozenset(record)
    if keys != expected:
        raise MalformedRegistryError(f"{where} keys differ: missing {sorted(expected - keys)}, unexpected {sorted(keys - expected)}.")


def _canonical(document: Mapping[str, object]) -> str:
    return json.dumps(document, sort_keys=True, separators=(",", ":"))


def _digest(content: Mapping[str, object]) -> str:
    return hashlib.sha256(_canonical(content).encode("utf-8")).hexdigest()
