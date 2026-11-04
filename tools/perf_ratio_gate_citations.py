"""The ratio gate's A/A calibration citations: cited records, margins, and their verification.

A calibration cites ``perf-ratio/1`` A/A session records by path (relative to the
registry root) and SHA-256, and carries the margins they give (spec §Margins:
``m_w = SAFETY * max_k |rho_{w,k} - 1|``, over at least `CALIBRATION_SESSIONS`
sessions). `calibration_of` builds a verified calibration from records on disk;
`_verify_calibration` re-verifies a stored one, live or an anchor link's: K,
distinct citations, each record's digest, each session an honest A/A session,
and the margins re-derived and compared with ``==``. `_load_cited` is the one load of a
cited record, for calibrations and anchor links alike: a record that has its cited digest
but does not load is a refusal of the cited evidence, not an infra failure. The registry
(`tools.perf_ratio_gate_registry`) and the anchor chain
(`tools.perf_ratio_gate_chain`) build on this module.
"""

from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from tools.perf_history_store import REPOSITORY
from tools.perf_ratio_gate_decide import Margin
from tools.perf_ratio_gate_errors import GateRefusal
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import InconsistentRatioRecordError
from tools.perf_ratio_record import MalformedRatioRecordError
from tools.perf_ratio_record import NonPositiveCpuError
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session
from tools.run_provenance import sha256_file

# Dimensionless factor on the worst observed A/A deviation |rho - 1| (spec §Margins): it doubles it.
SAFETY = 2.0
# K, the number of A/A sessions a calibration cites at least (spec §Margins: K = 5).
CALIBRATION_SESSIONS = 5

# Whose harness a calibration is checked against, as an A/A refusal names it: the registry's
# (live calibration) or an anchor link's gate session's (the calibration the link cites).
HarnessOwner = Literal["registry's", "gate session's"]

# What the harness loader `load_session` raises on a record that does not load. Its module is a HARNESS_SOURCE,
# so these cannot derive from `GateRefusal`; at a citation they are wrapped in one (`_load_cited`).
# UnicodeDecodeError: bytes that are not UTF-8, which `load_session` decodes before it parses JSON.
RECORD_LOAD_ERRORS: tuple[type[Exception], ...] = (MalformedRatioRecordError, InconsistentRatioRecordError, NonPositiveCpuError, UnicodeDecodeError)

# The ruling every unloadable-citation refusal states: committed evidence that stopped loading is the harness change's
# defect, never the gate's to recover (CI-N3 follow-up).
UNLOADABLE_RULING = (
    "The defect is in the harness change that broke loading: that change must keep committed records loadable or ship a migration; the gate has no recovery path of its own."
)


class CalibrationRecordDigestError(GateRefusal):
    """A cited calibration record's bytes no longer have the cited ``record_sha256``."""


class CalibrationRecordMissingError(GateRefusal):
    """A cited calibration record does not exist."""


class CalibrationRecordOutsideRootError(GateRefusal):
    """A calibration record lies outside the root its cited path is relative to."""


class CalibrationRecordUnloadableError(GateRefusal):
    """A cited calibration record has its cited digest but does not load as a ``perf-ratio/1`` session.

    The message names the cited path, the loader's error and `UNLOADABLE_RULING`: the harness change that broke
    loading owns the defect, and the gate has no recovery path of its own.
    """


class MarginMismatchError(GateRefusal):
    """A stored margin (or margin set) differs from the one its cited records re-derive."""


class InvalidCalibrationSessionError(GateRefusal):
    """A calibration session is not an honest A/A session of this registry and machine."""


class TooFewCalibrationSessionsError(GateRefusal):
    """A calibration cites fewer than `CALIBRATION_SESSIONS` sessions."""


class DuplicateCalibrationSessionError(GateRefusal):
    """A calibration cites one record (by path or by content) more than once."""


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
        _require_aa_session(session, ref_commit=ref_commit, harness_sha256=harness_sha256, harness_owner="registry's", cpu_model=cpu_model, workloads=workloads, where=record.path)
    return Calibration(cpu_model=cpu_model, sessions=tuple(cited), margins=derive_margins(sessions))


def _verify_calibration(calibration: Calibration, *, ref_commit: str, harness_sha256: str, harness_owner: HarnessOwner, workloads: frozenset[str], root: Path, where: str) -> None:
    # One calibration, live or a link's: K, distinct records, digests, honest A/A sessions, margins re-derived.
    # *harness_owner* names whose harness *harness_sha256* is (the registry's, or a link's gate session's).
    _require_session_count(len(calibration.sessions), where)
    _require_distinct(calibration.sessions, where)
    record_paths = [_cited_path(record, root, where) for record in calibration.sessions]
    for record, record_path in zip(calibration.sessions, record_paths):
        actual = sha256_file(record_path)
        if actual != record.record_sha256:
            raise CalibrationRecordDigestError(f"{where}: {record.path} now has SHA-256 {actual}, not the cited {record.record_sha256}; the record was rewritten.")
    sessions = [_load_cited(record, record_path, where) for record, record_path in zip(calibration.sessions, record_paths)]
    for record, session in zip(calibration.sessions, sessions):
        _require_aa_session(
            session,
            ref_commit=ref_commit,
            harness_sha256=harness_sha256,
            harness_owner=harness_owner,
            cpu_model=calibration.cpu_model,
            workloads=workloads,
            where=record.path,
        )
    derived = derive_margins(sessions, SAFETY)
    if dict(calibration.margins) != derived:
        raise MarginMismatchError(f"{where}: stored margins {dict(calibration.margins)} are not the ones its cited records give, {derived}.")


def _load_cited(record: CitedRecord, path: Path, where: str, *, unloadable: type[GateRefusal] = CalibrationRecordUnloadableError) -> PerfRatioSession:
    # The one load of a record cited by *where*, its digest already verified. Bytes the citation pins that no longer
    # load (a harness parser change, a hand edit) refuse the cited evidence: a rerun reads the same bytes, so this is
    # never infra. *unloadable* names the citing kind's refusal; the loader's error is chained.
    try:
        return load_session(path)
    except RECORD_LOAD_ERRORS as error:
        cause = f"{type(error).__name__}: {str(error).rstrip('.')}"
        raise unloadable(f"{where}: cited record {record.path} has its cited digest but does not load: {cause}. {UNLOADABLE_RULING}") from error


def _require_session_count(count: int, where: str) -> None:
    if count < CALIBRATION_SESSIONS:
        raise TooFewCalibrationSessionsError(f"{where} cites {count} A/A sessions; a calibration needs at least {CALIBRATION_SESSIONS}.")


def _require_distinct(cited: Sequence[CitedRecord], where: str, *, error: type[Exception] = DuplicateCalibrationSessionError) -> None:
    # No record cited twice, by path or by content; *error* names the citing kind's refusal.
    for attribute in ("path", "record_sha256"):
        values = [getattr(record, attribute) for record in cited]
        repeated = sorted({value for value in values if values.count(value) > 1})
        if repeated:
            raise error(f"{where} cites one record more than once (same {attribute}: {', '.join(repeated)}).")


def _require_aa_session(
    session: PerfRatioSession,
    *,
    ref_commit: str,
    harness_sha256: str,
    harness_owner: HarnessOwner,
    cpu_model: str,
    workloads: frozenset[str],
    where: str,
) -> None:
    # *harness_owner* names whose harness *harness_sha256* is, for the refusal: "registry's" or "gate session's".
    if session.injected is not None:
        raise InvalidCalibrationSessionError(f"{where} is an audit session (injected {session.injected.seconds} s on {session.injected.arm.value}), not an A/A calibration.")
    for arm in ArmName:
        identity = session.arms[arm]
        if identity.rev != ref_commit:
            raise InvalidCalibrationSessionError(f"{where}: the {arm.value} arm ran at {identity.rev}, not the reference {ref_commit}; an A/A session runs both arms there.")
        if identity.dirty:
            raise InvalidCalibrationSessionError(f"{where}: the {arm.value} arm's tree was dirty (patch {identity.patch_sha256}).")
    if session.harness_sha256 != harness_sha256:
        raise InvalidCalibrationSessionError(f"{where} ran harness {session.harness_sha256}, not the {harness_owner} {harness_sha256}.")
    if session.cpu_model != cpu_model:
        raise InvalidCalibrationSessionError(f"{where} ran on {session.cpu_model!r}, not the calibration's {cpu_model!r}.")
    changed = [workload.workload for workload in session.workloads if not workload.same_output]
    if changed:
        raise InvalidCalibrationSessionError(f"{where}: the arms' outputs differ on {', '.join(changed)}; A/A arms must agree.")
    if _workload_set(session) != workloads:
        raise InvalidCalibrationSessionError(f"{where} measures {sorted(_workload_set(session))}, not the workloads {sorted(workloads)}.")


def _workload_set(session: PerfRatioSession) -> frozenset[str]:
    return frozenset(workload.workload for workload in session.workloads)


def _cited_path(
    record: CitedRecord,
    root: Path,
    where: str,
    *,
    outside: type[Exception] = CalibrationRecordOutsideRootError,
    missing: type[Exception] = CalibrationRecordMissingError,
) -> Path:
    # Root containment and existence of one record cited by *where*; *outside* and *missing* name the citing kind's refusals.
    resolved_root = root.resolve()
    if Path(record.path).is_absolute():
        raise outside(f"{where}: cited path {record.path} is absolute; a registry cites paths relative to its root.")
    resolved = (resolved_root / record.path).resolve()
    if not resolved.is_relative_to(resolved_root):
        raise outside(f"{where}: cited path {record.path} leaves the root {root}.")
    if not resolved.is_file():
        raise missing(f"{where}: cited record {record.path} does not exist under {root}.")
    return resolved
