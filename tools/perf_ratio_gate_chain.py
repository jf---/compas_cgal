"""The ratio gate's anchor chain: promotion links, the reset rule, and the link replay.

Per workload the registry carries an anchor ``A_w`` (spec amendment 2026-11-01):
the reference's CPU time over that of the workload's anchor commit. Each promotion
appends an `AnchorLink` citing the promoted gate session, its extension sessions
and the calibration that supplied ``m_w``. `derive_anchors` replays the links in
order with the reset rule (`carry_anchor`) in exact `fractions.Fraction`
arithmetic over the recorded doubles; `link_ratios` replays one link, every
cited record verified, and re-decides it with the gate's own `decide`
(`tools.perf_ratio_gate_decide`): a link replays only if promote would have
accepted it. The calibration a link cites is verified by
`tools.perf_ratio_gate_citations`; chain continuity against a registry's
reference is checked by `tools.perf_ratio_gate_registry`.

Not every refusal is an `AnchorChainMismatchError`. `link_ratios` and
`derive_anchors` refuse with:

- an `AnchorChainMismatchError` subclass for every check of the link's own
  records and its replay: citations, digests, the gate session, extension
  targets, the pair cap and promotability;
- `ExtensionMismatchError` (`tools.perf_ratio_gate_merge`) when the merge
  refuses an extension;
- the calibration errors of `tools.perf_ratio_gate_citations` for the link's
  calibration: `TooFewCalibrationSessionsError`,
  `DuplicateCalibrationSessionError`, `CalibrationRecordOutsideRootError`,
  `CalibrationRecordMissingError`, `CalibrationRecordDigestError`,
  `InvalidCalibrationSessionError` and `MarginMismatchError`;
- the record-load errors of `tools.perf_ratio_record.load_session` for a cited
  record that is not a valid ``perf-ratio/1`` record.

A caller that has to catch every chain refusal catches all four groups.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from tools.perf_history_store import REPOSITORY
from tools.perf_ratio_gate_citations import Calibration
from tools.perf_ratio_gate_citations import CitedRecord
from tools.perf_ratio_gate_citations import _cited_path
from tools.perf_ratio_gate_citations import _require_distinct
from tools.perf_ratio_gate_citations import _verify_calibration
from tools.perf_ratio_gate_citations import _workload_set
from tools.perf_ratio_gate_decide import MAX_PAIRS
from tools.perf_ratio_gate_decide import UNIT_ANCHOR
from tools.perf_ratio_gate_decide import Anchor
from tools.perf_ratio_gate_decide import ExactRatio
from tools.perf_ratio_gate_decide import Margin
from tools.perf_ratio_gate_decide import Verdict
from tools.perf_ratio_gate_decide import decide
from tools.perf_ratio_gate_decide import exact_ratio
from tools.perf_ratio_gate_merge import Merge
from tools.perf_ratio_gate_merge import merge_extension
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import WorkloadRatio
from tools.perf_ratio_record import load_session
from tools.run_provenance import sha256_file


class AnchorChainMismatchError(Exception):
    """The anchor chain does not re-derive the registry's anchors; each failure mode is a subclass."""


class AnchorRecordMissingError(AnchorChainMismatchError):
    """A record cited by an anchor link does not exist."""


class AnchorRecordDigestError(AnchorChainMismatchError):
    """A record cited by an anchor link no longer has the cited ``record_sha256``."""


class AnchorRecordOutsideRootError(AnchorChainMismatchError):
    """A record cited by an anchor link lies outside the registry root."""


class AnchorChainBrokenError(AnchorChainMismatchError):
    """A link does not start where the chain stands: the first at ``anchor_ref``, each later one at the previous ``to_ref``."""


class AnchorChainEndError(AnchorChainMismatchError):
    """The chain does not end at the registry's reference (an empty chain: ``anchor_ref`` is not ``ref_commit``)."""


class AnchorValueMismatchError(AnchorChainMismatchError):
    """The stored anchors are not the ones the chain re-derives (value or workload set)."""


class InvalidAnchorSessionError(AnchorChainMismatchError):
    """A link's cited sessions are not the honest record of promoting ``from_ref`` to ``to_ref``."""


class AnchorLinkSelfLoopError(AnchorChainMismatchError):
    """A link promotes a commit to itself (``from_ref == to_ref``): the anchor would move while the reference does not."""


class DuplicateAnchorRecordError(AnchorChainMismatchError):
    """A link cites one record (by path or by content) more than once."""


class RegressedAnchorLinkError(AnchorChainMismatchError):
    """A link's (extended) session leaves a workload red, ``A_w * rho_w > 1 + m_w``: promote would have refused it."""


class UnimprovedAnchorLinkError(AnchorChainMismatchError):
    """A link's (extended) session has no certified improvement, ``A_w * rho_w < 1 - m_w``: promote would have refused it."""


class UnrequestedExtensionError(AnchorChainMismatchError):
    """A link's extension targets a workload that was not red at that point: the gate never requested it."""


class ExtensionPairCapError(AnchorChainMismatchError):
    """A link's workload holds more than `MAX_PAIRS` pairs after its extensions: the gate stops extending there."""


@dataclass(frozen=True)
class AnchorLink:
    """One promotion of the reference, cited by its records. Raw construction is unchecked.

    Attributes:
        from_ref: The reference the promoted gate session measured against.
        to_ref: The promoted commit, the session's HEAD arm.
        sessions: The promoted gate session, then its extension sessions in merge order.
        calibration: The A/A calibration at *from_ref* that supplied the session machine's ``m_w``.
    """

    from_ref: str
    to_ref: str
    sessions: tuple[CitedRecord, ...]
    calibration: Calibration

    def to_json(self) -> dict[str, object]:
        """The link as a JSON object."""
        return {
            "from_ref": self.from_ref,
            "to_ref": self.to_ref,
            "sessions": [record.to_json() for record in self.sessions],
            "calibration": self.calibration.to_json(),
        }


def carry_anchor(anchor: Anchor, ratio: ExactRatio, margin: Margin) -> Anchor:
    """The reset rule: a workload's anchor after a promotion whose (extended) session measured *ratio*.

    ``A' = A * rho`` when ``A * rho >= 1 - m``: the new reference did not certifiably beat the
    anchor, so the drift is carried. Otherwise ``A' = 1``: the new reference is a certified
    improvement and becomes the anchor. Exact: *margin* enters as the rational its double is.

    Args:
        anchor: ``A_w`` before the promotion.
        ratio: The promoted session's ``rho_w`` (`exact_ratio`).
        margin: The ``m_w`` the promoted session was gated with.

    Returns:
        ``A'_w``.
    """
    carried = anchor * ratio
    if carried >= 1 - Fraction(margin):
        return Anchor(carried)
    return UNIT_ANCHOR


def link_ratios(link: AnchorLink, anchors: Mapping[str, Anchor], *, root: Path = REPOSITORY, merge: Merge = merge_extension) -> dict[str, ExactRatio]:
    """Each workload's ``rho_w`` in *link*'s promoted session, its extensions replayed, every record verified.

    The replay re-decides the link as the gate did (`decide`, under *anchors* and the link
    calibration's margins): each extension, in merge order, must target a workload that was red
    at that point; no workload may end with more than `MAX_PAIRS` pairs; and the extended session
    must be promotable: no workload red, at least one certified improvement.

    Args:
        link: The link.
        anchors: Workload -> ``A_w`` before the link; a workload absent here enters at `UNIT_ANCHOR`.
        root: The directory cited paths are relative to (`REPOSITORY`).
        merge: Merges an extension into a run (`merge_extension`, the gate's own).

    Returns:
        Workload of the gate session -> its exact ratio (`exact_ratio`) after every extension of it.

    Raises:
        InvalidAnchorSessionError: No session is cited; the gate session is injected, dirty, not
            ``from_ref -> to_ref``, or not on the calibration's machine; an extension measures
            other than one workload of the gate session.
        DuplicateAnchorRecordError: The link cites one record (by path or content) twice.
        UnrequestedExtensionError: An extension targets a workload that was not red at that point.
        ExtensionPairCapError: A workload holds more than `MAX_PAIRS` pairs after its extensions.
        RegressedAnchorLinkError: A workload ends with ``A_w * rho_w > 1 + m_w``.
        UnimprovedAnchorLinkError: No workload ends with ``A_w * rho_w < 1 - m_w``.
        AnchorRecordOutsideRootError: A cited path leaves *root*.
        AnchorRecordMissingError: A cited record does not exist.
        AnchorRecordDigestError: A cited record's bytes changed.
        ExtensionMismatchError: *merge* refuses an extension (`merge_extension`).
        TooFewCalibrationSessionsError: As `load_registry`, for the link's calibration.
        DuplicateCalibrationSessionError: As `load_registry`, for the link's calibration.
        CalibrationRecordOutsideRootError: As `load_registry`, for the link's calibration.
        CalibrationRecordMissingError: As `load_registry`, for the link's calibration.
        CalibrationRecordDigestError: As `load_registry`, for the link's calibration.
        InvalidCalibrationSessionError: The link's calibration is not honest A/A at ``from_ref`` under
            the session's harness, on its machine, over its workloads.
        MarginMismatchError: The link's calibration margins do not re-derive.
    """
    where = f"anchor link {link.from_ref} -> {link.to_ref}"
    if not link.sessions:
        raise InvalidAnchorSessionError(f"{where} cites no session; a link cites its promoted gate session first.")
    paths = [_cited_path(record, root, where, outside=AnchorRecordOutsideRootError, missing=AnchorRecordMissingError) for record in link.sessions]
    _require_distinct(link.sessions, where, error=DuplicateAnchorRecordError)
    for record, path in zip(link.sessions, paths):
        actual = sha256_file(path)
        if actual != record.record_sha256:
            raise AnchorRecordDigestError(f"{where}: {record.path} now has SHA-256 {actual}, not the cited {record.record_sha256}; the record was rewritten.")
    gate_session, *extensions = [load_session(path) for path in paths]
    _require_promoted_session(gate_session, link, where)
    _verify_calibration(
        link.calibration,
        ref_commit=link.from_ref,
        harness_sha256=gate_session.harness_sha256,
        harness_owner="gate session's",
        workloads=_workload_set(gate_session),
        root=root,
        where=f"{where}: calibration {link.calibration.cpu_model!r}",
    )
    runs: dict[str, WorkloadRatio] = {run.workload: run for run in gate_session.workloads}
    margins = link.calibration.margins

    def verdict(name: str) -> Verdict:
        return decide(runs[name], margins[name], anchors.get(name, UNIT_ANCHOR))

    for record, extension in zip(link.sessions[1:], extensions):
        names = [run.workload for run in extension.workloads]
        if len(names) != 1 or names[0] not in runs:
            raise InvalidAnchorSessionError(f"{where}: extension {record.path} measures {names}, not one workload of the gate session {sorted(runs)}.")
        (name,) = names
        before = verdict(name)
        if before is not Verdict.REGRESSION:
            raise UnrequestedExtensionError(f"{where}: extension {record.path} extends {name}, {before.value} at that point; the gate extends only red workloads.")
        runs[name] = merge(runs[name], gate_session, extension)
    _require_promotable(runs, verdict, where)
    return {name: exact_ratio(run) for name, run in runs.items()}


def _require_promotable(runs: Mapping[str, WorkloadRatio], verdict: Callable[[str], Verdict], where: str) -> None:
    # The extended session is one promote accepts: every run within MAX_PAIRS, none red, one improved.
    overfull = sorted(name for name, run in runs.items() if len(run.pairs) > MAX_PAIRS)
    if overfull:
        raise ExtensionPairCapError(f"{where}: {', '.join(overfull)} hold more than MAX_PAIRS {MAX_PAIRS} pairs; the gate stops extending there.")
    verdicts = {name: verdict(name) for name in runs}
    regressed = sorted(name for name, decided in verdicts.items() if decided is Verdict.REGRESSION)
    if regressed:
        raise RegressedAnchorLinkError(f"{where}: {', '.join(regressed)} end red (A_w*rho_w > 1 + m_w); promote refuses a regressed session.")
    if Verdict.IMPROVED not in verdicts.values():
        raise UnimprovedAnchorLinkError(f"{where}: no workload is a certified improvement (A_w*rho_w < 1 - m_w); promote needs one.")


def derive_anchors(links: Sequence[AnchorLink], workloads: Sequence[str], *, root: Path = REPOSITORY, merge: Merge = merge_extension) -> dict[str, Anchor]:
    """Each workload's exact anchor ``A_w`` after *links*, by the reset rule (`carry_anchor`).

    Every workload starts at 1. Each link replays its records (`link_ratios`) and carries or
    resets the anchor of every workload its session measured, with that link's calibrated
    margin. A workload a link did not measure has left the registry: its anchor is dropped,
    and a workload (re-)entering later starts at 1 at the link that first measures it.

    Args:
        links: The chain, in order.
        workloads: The registry's workloads: the anchors returned.
        root: The directory cited paths are relative to (`REPOSITORY`).
        merge: Merges an extension into a run (`merge_extension`).

    Returns:
        Workload of *workloads* -> its anchor; 1 for a workload no link carried.

    Raises:
        InvalidAnchorSessionError: As `link_ratios`.
        DuplicateAnchorRecordError: As `link_ratios`.
        UnrequestedExtensionError: As `link_ratios`.
        ExtensionPairCapError: As `link_ratios`.
        RegressedAnchorLinkError: As `link_ratios`.
        UnimprovedAnchorLinkError: As `link_ratios`.
        AnchorRecordOutsideRootError: As `link_ratios`.
        AnchorRecordMissingError: As `link_ratios`.
        AnchorRecordDigestError: As `link_ratios`.
        ExtensionMismatchError: As `link_ratios`.
        TooFewCalibrationSessionsError: As `link_ratios`.
        DuplicateCalibrationSessionError: As `link_ratios`.
        CalibrationRecordOutsideRootError: As `link_ratios`.
        CalibrationRecordMissingError: As `link_ratios`.
        CalibrationRecordDigestError: As `link_ratios`.
        InvalidCalibrationSessionError: As `link_ratios`.
        MarginMismatchError: As `link_ratios`.
        MalformedRatioRecordError: A cited record is not JSON or not a ``perf-ratio/1`` record (`load_session`).
        InconsistentRatioRecordError: A cited record's stored derived value differs from its recomputation (`load_session`).
        NonPositiveCpuError: A cited record holds a run whose CPU is not positive (`load_session`).
    """
    anchors: dict[str, Anchor] = {}
    for link in links:
        ratios = link_ratios(link, anchors, root=root, merge=merge)
        anchors = {name: carry_anchor(anchors.get(name, UNIT_ANCHOR), ratio, link.calibration.margins[name]) for name, ratio in ratios.items()}
    return {name: anchors.get(name, UNIT_ANCHOR) for name in workloads}


def _require_promoted_session(session: PerfRatioSession, link: AnchorLink, where: str) -> None:
    # The cited gate session is the honest record of promoting link.from_ref to link.to_ref.
    if session.injected is not None:
        raise InvalidAnchorSessionError(f"{where}: the gate session is an audit session (injected {session.injected.seconds} s on {session.injected.arm.value}).")
    if session.ref_commit != link.from_ref:
        raise InvalidAnchorSessionError(f"{where}: the gate session measured against {session.ref_commit}, not the link's from_ref.")
    head = session.arms[ArmName.HEAD].rev
    if head != link.to_ref:
        raise InvalidAnchorSessionError(f"{where}: the gate session's HEAD arm ran at {head}, not the link's to_ref.")
    for arm in ArmName:
        if session.arms[arm].dirty:
            raise InvalidAnchorSessionError(f"{where}: the gate session's {arm.value} arm tree was dirty (patch {session.arms[arm].patch_sha256}).")
    if session.cpu_model != link.calibration.cpu_model:
        raise InvalidAnchorSessionError(f"{where}: the gate session ran on {session.cpu_model!r}, not the link calibration's {link.calibration.cpu_model!r}.")
