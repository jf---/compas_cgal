"""The ratio gate CLI: ``python -m tools.perf_ratio_gate check | calibrate | promote``.

- ``check RECORD`` gates a ``perf-ratio/1`` session against the registry (`gate`) and
  prints one verdict per workload. A regression raises `RatioRegressionError`.
- ``calibrate --records PATH...`` builds one machine's A/A calibration from those
  records and merges it into the registry. ``calibrate --sessions K`` first measures K
  A/A sessions at the registry's reference, through the same `Measure` seam. Records
  under another harness are refused unless ``--rebase-harness`` is given: the registry then
  takes their harness and workload set and drops every other calibration (all void under
  the new harness), at the same reference.
- ``promote RECORD`` moves the reference to the session's HEAD (spec amendment 2026-11-01,
  the anchor chain): only from a clean HEAD, other than the reference, whose session passes
  against every workload's anchor (``A_w * rho_w <= 1 + m_w``) with at least one IMPROVED
  workload. The gate session and its extensions are committed under `promotion_dir`, and the
  promoted registry appends the `AnchorLink` citing them: each workload's anchor is carried
  (``A'_w = A_w * rho_w``) or, for an IMPROVED workload, reset to 1 (`carry_anchor`). It holds
  no calibrations, so every check is `UncalibratedMachineError` until an A/A calibration runs
  at the new reference. Calibrating keeps the chain, a harness rebase included.

Red workloads are extended through `Measure`; in production that runs the harness
(`harness_command`) as a subprocess, its records under `GATE_SCRATCH`.
"""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import tempfile
from collections.abc import Callable
from collections.abc import Iterator
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Optional

from tools.perf_history_store import REPOSITORY
from tools.perf_ratio import DEFAULT_PAIRS
from tools.perf_ratio import harness_sha256
from tools.perf_ratio_arms import DEFAULT_REF_WORKTREE
from tools.perf_ratio_arms import checkout_head
from tools.perf_ratio_gate_chain import AnchorLink
from tools.perf_ratio_gate_chain import AnchorValueMismatchError
from tools.perf_ratio_gate_chain import carry_anchor
from tools.perf_ratio_gate_citations import CitedRecord
from tools.perf_ratio_gate_decide import Anchor
from tools.perf_ratio_gate_decide import ExactRatio
from tools.perf_ratio_gate_registry import CALIBRATE_TASK
from tools.perf_ratio_gate_registry import CALIBRATION_SESSIONS
from tools.perf_ratio_gate_registry import REGISTRY_PATH
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import TooFewCalibrationSessionsError
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import calibration_dir
from tools.perf_ratio_gate_registry import calibration_of
from tools.perf_ratio_gate_registry import load_registry
from tools.perf_ratio_gate_registry import promotion_dir
from tools.perf_ratio_gate_registry import write_registry
from tools.perf_ratio_gate_verdict import ExtensionRequest
from tools.perf_ratio_gate_verdict import GateReport
from tools.perf_ratio_gate_verdict import Measure
from tools.perf_ratio_gate_verdict import RatioRegressionError
from tools.perf_ratio_gate_verdict import Verdict
from tools.perf_ratio_gate_verdict import anchors_for
from tools.perf_ratio_gate_verdict import gate
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import MalformedRatioRecordError
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session
from tools.run_provenance import machine_identity
from tools.run_provenance import sha256_file

# Extension and calibration-run records before they are committed (``build/`` is gitignored).
GATE_SCRATCH = REPOSITORY / "build" / "ratio-gate"
# The file name of one A/A session record inside its ``aa-<k>`` directory (the landed layout).
SESSION_RECORD_NAME = "perf_ratio_session.json"
# Prefix of the k-th A/A session's directory under a calibration directory, k = 1..K.
AA_DIR_PREFIX = "aa-"
# Prefix of the scratch directory a ``calibrate --sessions`` run stages its records in before the rename.
CALIBRATION_STAGING_PREFIX = "calibration-staging-"
# The directory of a promotion's gate session record inside its `promotion_dir`.
GATE_RECORD_DIR = "gate"
# Prefix of the n-th extension session's directory inside a `promotion_dir`, n = 1.. in merge order.
EXTENSION_DIR_PREFIX = "ext-"
# Prefix of the scratch directory a promotion stages its records in before the rename.
PROMOTION_STAGING_PREFIX = "promotion-staging-"
# The promote pixi task the printed next step names (the calibrate task is the registry's `CALIBRATE_TASK`).
PROMOTE_TASK = "pixi run perf-ratio-promote"

# Reads this machine's CPU model, as the harness records it in ``cpu_model``.
CpuModelReader = Callable[[], str]
# Reads the SHA-256 of the harness this tree runs, as the harness records it in ``harness_sha256``.
HarnessReader = Callable[[], str]
# Moves the clean worktree at the path to the commit the rev names; returns its full id.
Checkout = Callable[[Path, str], str]


class MissingRegistryError(Exception):
    """The registry file does not exist where a subcommand needs one."""


class CalibrationExistsError(Exception):
    """``calibrate --sessions`` would write into a calibration directory that already exists."""


class CalibrationMachineMismatchError(Exception):
    """A measured A/A session ran on another CPU model than the one its directory is named after."""


class CalibrationHarnessMismatchError(Exception):
    """Calibration records ran another harness than the registry's, and ``--rebase-harness`` was not given."""


class RebaseWithoutHarnessChangeError(Exception):
    """``--rebase-harness`` was given, but the records ran the registry's own harness."""


class PromotionRefusedError(Exception):
    """A session does not qualify for promotion: a regression, no improvement, a dirty HEAD, or HEAD is the reference."""


class PromotionRecordsExistError(Exception):
    """The promotion's record directory exists with other records than the ones this promotion writes."""


class PromotedReferenceMismatchError(Exception):
    """The reference worktree's checkout landed on another commit than the promoted HEAD."""


def harness_command(request: ExtensionRequest, write: Path) -> tuple[str, ...]:
    """The harness command line that measures *request* and writes its session to *write*.

    Args:
        request: What to measure: reference, HEAD rev, workloads, pairs, and the injection to replay.
        write: Where the harness writes the session record.

    Returns:
        ``python -m tools.perf_ratio`` with ``--ref-rev``, ``--head-rev``, one ``--workload`` per
        workload, ``--pairs``, ``--write`` and, for an audit request, ``--inject-cpu-seconds``.
    """
    command = (
        sys.executable,
        "-m",
        "tools.perf_ratio",
        "--ref-rev",
        request.ref_commit,
        "--head-rev",
        request.head_rev,
        *itertools.chain.from_iterable(("--workload", workload) for workload in request.workloads),
        "--pairs",
        str(request.pairs),
        "--write",
        str(write),
    )
    if request.injected is None:
        return command
    return (*command, "--inject-cpu-seconds", f"{request.injected.arm.value}={request.injected.seconds!r}")


def subprocess_measure(scratch: Path) -> Measure:
    """A `Measure` that runs the harness for each request and reads the session it wrote.

    Each run writes ``scratch / f"ext-{n}.json"`` for the next n (counting up from 1) whose
    file does not exist yet, so a record from an earlier run is never read back as this one.

    Args:
        scratch: The directory the harness writes its records into; created on first use.

    Returns:
        The measure. It raises ``subprocess.CalledProcessError`` when the harness fails.
    """
    numbers = itertools.count(1)

    def measure(request: ExtensionRequest) -> PerfRatioSession:
        scratch.mkdir(parents=True, exist_ok=True)
        write = next(path for path in (scratch / f"ext-{n}.json" for n in numbers) if not path.exists())
        subprocess.run(harness_command(request, write), check=True, cwd=REPOSITORY)
        return load_session(write)

    return measure


def report_lines(report: GateReport) -> list[str]:
    """One line per workload verdict, then the summary line.

    Each workload line names ``rho_w``, its anchor ``A_w`` and the decided product ``A_w*rho_w``
    (reporting floats; the verdict was decided exactly), ``m_w`` and the pair count.
    """
    lines = []
    for verdict in report.verdicts:
        changed = " changed output" if verdict.changed_output else ""
        lines.append(
            f"{verdict.workload}: {verdict.verdict.value} rho {verdict.ratio:.3f} A_w {verdict.anchor:.3f} A_w*rho_w {verdict.cumulative_ratio:.3f} "
            f"margin {verdict.margin:.3f} pairs {verdict.pairs}{changed}"
        )
    counts = ", ".join(f"{sum(1 for verdict in report.verdicts if verdict.verdict is kind)} {kind.value}" for kind in (Verdict.PASS, Verdict.IMPROVED, Verdict.REGRESSION))
    mode = " (audit)" if report.audit else ""
    lines.append(f"ratio gate{mode}: {len(report.verdicts)} workloads: {counts}")
    return lines


def promotion_refusals(session: PerfRatioSession, report: GateReport, ref_commit: str) -> tuple[str, ...]:
    """Why *session*, gated as *report*, may not be promoted; empty when it may.

    Args:
        session: The gated session.
        report: Its gate report.
        ref_commit: The registry's reference commit.

    Returns:
        One message per reason: a regression, no improvement, a dirty HEAD tree, a HEAD that is the
        reference itself (an A/A session's IMPROVED is noise), an audit-mode report.
    """
    reasons = []
    if report.regressions:
        reasons.append(
            f"REGRESSION on {', '.join(verdict.workload for verdict in report.regressions)} (A_w*rho_w > 1 + m_w); "
            "a promotion carries every workload's drift into its anchor, so none may end above its margin"
        )
    if not report.improvements:
        reasons.append("no IMPROVED workload; promotion needs a certified improvement (A_w*rho_w < 1 - m_w)")
    head = session.arms[ArmName.HEAD]
    if head.dirty:
        reasons.append(f"the HEAD arm's tree was dirty (patch {head.patch_sha256}); promotion needs a clean, committed HEAD")
    if head.rev == ref_commit:
        reasons.append(f"the HEAD arm is the reference {ref_commit} itself; promotion moves the reference to another commit")
    if report.audit:
        reasons.append("an audit-mode report (injected or partial session) is never a verdict on HEAD")
    return tuple(reasons)


def _load(path: Path, root: Path) -> GateRegistry:
    if not path.is_file():
        raise MissingRegistryError(f"no ratio gate registry at {path}; build one with {CALIBRATE_TASK} --records PATH...")
    return load_registry(path, root=root)


def check(record: Path, *, registry_path: Path, root: Path, measure: Measure, audit: bool) -> int:
    """Gate the session at *record* and print its verdicts.

    Args:
        record: The ``perf-ratio/1`` session record.
        registry_path: The registry file.
        root: The directory the registry's cited records are relative to.
        measure: Measures extensions of red workloads.
        audit: Audit mode: an injected or partial session is decided, never offered for promotion.

    Returns:
        0: no workload regressed.

    Raises:
        MissingRegistryError: No registry at *registry_path*.
        RatioRegressionError: A workload stayed above ``1 + m_w`` after extension.
    """
    registry = _load(registry_path, root)
    session = load_session(record)
    report = gate(session, registry, measure, audit=audit)
    for line in report_lines(report):
        print(line)
    if not promotion_refusals(session, report, registry.ref_commit):
        print(f"promotion available: {PROMOTE_TASK} {record}")
    if report.regressions:
        raise RatioRegressionError(report)
    return 0


def calibrate_records(record_paths: Sequence[Path], *, registry_path: Path, root: Path, rebase_harness: bool = False) -> GateRegistry:
    """Calibrate the records' machine and merge it into the registry at *registry_path*.

    The new calibration replaces that ``cpu_model``'s and keeps every other. With no registry
    yet, the reference and harness are the first record's, and every record must agree
    (`calibration_of` refuses one that does not). Into an existing registry, records under
    another harness are refused unless *rebase_harness*: then the registry takes the first
    record's harness and the records' workload set, keeps its reference (every record must
    still be at it), and drops every other calibration, since all are void under a new harness.
    An existing registry's anchor chain (``anchor_ref``, ``anchor_links``) is kept either way, so a
    recalibration never resets an anchor.

    Args:
        record_paths: The A/A session records, at least `CALIBRATION_SESSIONS`, all on one machine.
        registry_path: The registry file, written in place.
        root: The directory cited records are relative to.
        rebase_harness: Rebase an existing registry onto the records' (edited) harness.

    Returns:
        The registry written.

    Raises:
        MissingRegistryError: *rebase_harness* and no registry at *registry_path*.
        CalibrationHarnessMismatchError: A record ran another harness than the registry's, without
            *rebase_harness*.
        RebaseWithoutHarnessChangeError: *rebase_harness*, but the records ran the registry's harness.
        TooFewCalibrationSessionsError: As `calibration_of`.
        InvalidCalibrationSessionError: As `calibration_of`: a record is not an honest A/A session
            at the reference, under the harness, on the first record's machine and workload set.
        CalibrationWorkloadMismatchError: The records' workload set is not the registry's.
    """
    existing = _load(registry_path, root) if rebase_harness or registry_path.exists() else None
    if existing is None:
        first = load_session(record_paths[0])
        ref_commit, harness_sha256, workloads = first.ref_commit, first.harness_sha256, None
    elif rebase_harness:
        harness_sha256 = load_session(record_paths[0]).harness_sha256
        if harness_sha256 == existing.harness_sha256:
            raise RebaseWithoutHarnessChangeError(
                f"--rebase-harness: the records ran the registry's own harness {harness_sha256}; calibrate without it, which keeps the other machines' calibrations."
            )
        ref_commit, workloads = existing.ref_commit, None
    else:
        stale = [str(path) for path in record_paths if load_session(path).harness_sha256 != existing.harness_sha256]
        if stale:
            raise CalibrationHarnessMismatchError(
                f"{', '.join(stale)} ran another harness than the registry's {existing.harness_sha256}. After a harness edit, re-run the A/A "
                f"sessions at the reference {existing.ref_commit} and calibrate with --rebase-harness, which voids every other calibration."
            )
        ref_commit, harness_sha256, workloads = existing.ref_commit, existing.harness_sha256, existing.workloads
    calibration = calibration_of(record_paths, ref_commit=ref_commit, harness_sha256=harness_sha256, root=root)
    kept = [] if existing is None or rebase_harness else [other for other in existing.calibrations if other.cpu_model != calibration.cpu_model]
    # A recalibration (a harness rebase included) keeps the anchor chain: rho_w is a same-harness
    # ratio, so links measured under another harness still compose (spec amendment 2026-11-01).
    anchor_ref = None if existing is None else existing.chain_start
    anchor_links = () if existing is None else existing.anchor_links
    registry = build_registry(
        ref_commit=ref_commit,
        harness_sha256=harness_sha256,
        calibrations=[*kept, calibration],
        workloads=workloads,
        anchor_ref=anchor_ref,
        anchor_links=anchor_links,
        root=root,
    )
    write_registry(registry, registry_path, root=root)
    if existing is not None and rebase_harness:
        dropped = ", ".join(repr(other.cpu_model) for other in existing.calibrations) or "none"
        print(f"harness rebased: {existing.harness_sha256} -> {harness_sha256}; dropped calibrations: {dropped}.")
    print(f"calibrated {calibration.cpu_model!r} at {ref_commit} from {len(record_paths)} A/A sessions:")
    for workload, margin in sorted(calibration.margins.items()):
        print(f"  {workload}: margin {margin:.3f}")
    return registry


def calibrate_sessions(
    sessions: int,
    *,
    registry_path: Path,
    root: Path,
    measure: Measure,
    cpu_model: CpuModelReader,
    tree_harness: HarnessReader,
    rebase_harness: bool = False,
) -> GateRegistry:
    """Measure *sessions* A/A sessions at the registry's reference, then calibrate from them.

    Session k is written to ``calibration_dir(ref, cpu)/aa-<k>/perf_ratio_session.json`` (under
    *root*), k = 1..*sessions*. Every refusal but a session from another machine fires before
    anything is measured, the two harness refusals included: *measure* runs this tree's harness,
    so its digest (*tree_harness*) is what the measured records will carry. All K sessions are measured before any record is written; the records
    are staged under `GATE_SCRATCH` and moved into place with one rename, so a failed run leaves
    no calibration directory and a rerun starts afresh. The scratch and the results tree are
    both inside the repository, so the rename stays on one filesystem; if they ever were not,
    the rename raises ``OSError`` (EXDEV) and the target is still not created.

    Args:
        sessions: K, at least `CALIBRATION_SESSIONS`.
        registry_path: The registry file; it must exist (it names the reference and workloads).
        root: The directory the calibration directory and cited records live under.
        measure: Measures each A/A session.
        cpu_model: Reads this machine's CPU model.
        tree_harness: Reads the digest of the harness *measure* runs (this tree's).
        rebase_harness: Calibrate the measured records as `calibrate_records` with *rebase_harness*.

    Returns:
        The registry written.

    Raises:
        TooFewCalibrationSessionsError: *sessions* is below `CALIBRATION_SESSIONS`.
        MissingRegistryError: No registry at *registry_path*.
        CalibrationHarnessMismatchError: Before measuring: the tree's harness is not the registry's,
            without *rebase_harness*. After measuring: as `calibrate_records`.
        RebaseWithoutHarnessChangeError: Before measuring: *rebase_harness*, but the tree's harness
            is the registry's. After measuring: as `calibrate_records`.
        CalibrationExistsError: The calibration directory already exists.
        CalibrationMachineMismatchError: A measured session ran on another CPU model.
        subprocess.CalledProcessError: A harness run failed (from the production `Measure`).
    """
    if sessions < CALIBRATION_SESSIONS:
        raise TooFewCalibrationSessionsError(f"--sessions {sessions}: a calibration needs at least {CALIBRATION_SESSIONS} A/A sessions.")
    registry = _load(registry_path, root)
    tree = tree_harness()
    if rebase_harness and tree == registry.harness_sha256:
        raise RebaseWithoutHarnessChangeError(
            f"--rebase-harness: this tree runs the registry's own harness {tree}; calibrate without it, which keeps the other machines' calibrations."
        )
    if not rebase_harness and tree != registry.harness_sha256:
        raise CalibrationHarnessMismatchError(
            f"this tree runs harness {tree}, not the registry's {registry.harness_sha256}. After a harness edit, calibrate with "
            f"--rebase-harness: the A/A sessions run at the reference {registry.ref_commit} and every other calibration is voided."
        )
    cpu = cpu_model()
    target = root / calibration_dir(registry.ref_commit, cpu).relative_to(REPOSITORY)
    if target.exists():
        raise CalibrationExistsError(f"{target} exists; an A/A calibration of {cpu!r} at {registry.ref_commit} was already run there. Calibrate from its records with --records.")
    request = ExtensionRequest(ref_commit=registry.ref_commit, head_rev=registry.ref_commit, workloads=registry.workloads, pairs=DEFAULT_PAIRS, injected=None)
    measured = []
    for index in range(1, sessions + 1):
        session = measure(request)
        if session.cpu_model != cpu:
            raise CalibrationMachineMismatchError(f"A/A session {index} ran on {session.cpu_model!r}, not this machine's {cpu!r} that {target} is named after.")
        measured.append(session)
        print(f"A/A session {index}/{sessions} measured")
    # Stage the K records in the gitignored scratch, then move them into place in one rename: a run
    # that fails at any point before the rename leaves nothing under the committed results tree.
    scratch = root / GATE_SCRATCH.relative_to(REPOSITORY)
    scratch.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=CALIBRATION_STAGING_PREFIX, dir=scratch))
    for index, session in enumerate(measured, start=1):
        record = staging / f"{AA_DIR_PREFIX}{index}" / SESSION_RECORD_NAME
        record.parent.mkdir()
        record.write_text(session.dumps(), encoding="utf-8")
    target.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(target)
    written = [target / f"{AA_DIR_PREFIX}{index}" / SESSION_RECORD_NAME for index in range(1, sessions + 1)]
    print(f"wrote {sessions} A/A session records to {target}")
    return calibrate_records(written, registry_path=registry_path, root=root, rebase_harness=rebase_harness)


@dataclass(frozen=True)
class AnchorMove:
    """One workload's anchor across a promotion (spec amendment 2026-11-01).

    Attributes:
        workload: The workload's name.
        before: ``A_w`` before the promotion.
        after: ``A'_w`` (`carry_anchor`): ``A_w * rho_w``, or 1 on a reset.
        reset: Whether the workload was IMPROVED (``A_w * rho_w < 1 - m_w``), so its anchor reset to 1
            rather than carried; the reset condition of `carry_anchor` is exactly `decide`'s IMPROVED.
    """

    workload: str
    before: Anchor
    after: Anchor
    reset: bool


def anchor_moves(report: GateReport, anchors: Mapping[str, Anchor]) -> tuple[AnchorMove, ...]:
    """Each workload's anchor move if *report*'s session is promoted, in report order.

    Args:
        report: The promoted session's gate report; each ``rho_w`` is the extended run's.
        anchors: Workload -> its exact ``A_w`` the report was decided under (`anchors_for`).

    Returns:
        One move per verdict, by the reset rule (`carry_anchor`) under the verdict's ``m_w``.
    """
    return tuple(
        AnchorMove(
            workload=verdict.workload,
            before=anchors[verdict.workload],
            # The verdict's ratio is the run's recorded double, so this is the rational `exact_ratio` reads.
            after=carry_anchor(anchors[verdict.workload], ExactRatio(Fraction(verdict.ratio)), verdict.margin),
            reset=verdict.verdict is Verdict.IMPROVED,
        )
        for verdict in report.verdicts
    )


def promote(record: Path, *, registry_path: Path, root: Path, measure: Measure, checkout: Checkout) -> GateRegistry:
    """Move the reference to the HEAD of the session at *record*, appending its anchor link.

    After every refusal and before the reference worktree moves, the gate session (its record's
    bytes, verbatim) and every extension the gate merged (`GateReport.extensions`, merge order) are
    written under ``promotion_dir(ref, head)`` as ``gate/`` and ``ext-<n>/`` session records: staged
    under `GATE_SCRATCH` and moved into place with one rename, as `calibrate_sessions` does. The
    promoted registry appends the `AnchorLink` that cites them with this machine's calibration,
    keeps ``anchor_ref``, and re-derives every anchor from the chain; those must be the anchors the
    gate report gives (`anchor_moves`).

    *record* is read once: the session is parsed from those bytes, and the gate, the persisted
    ``gate/`` copy (and its citation) and the rerun's identity check all use them. A record
    rewritten while the gate runs (extensions can take minutes) does not change what is decided
    or persisted.

    A failed checkout or registry write leaves the registry unchanged and the records in place. A
    rerun on the same record finishes the promotion: when the promotion directory exists, its gate
    record must be byte-identical to *record*, the gate replays the persisted extensions in order
    instead of measuring new ones (so the rerun decides on the evidence already written), and the
    directory must then hold exactly the files this promotion would write, byte for byte. Anything
    else is another promotion's records and is refused. The directory is not bound to the
    calibration it was written under: a rerun after a same-reference recalibration of this machine
    replays the written extensions under the new calibration's margins, and the link cites the new
    calibration (`link_ratios` re-verifies that replay on every load). This is no more than a fresh
    ``promote`` under the new calibration could do with the same session.

    `AnchorValueMismatchError` fires after the records are staged, so they stay in place; a rerun
    replays the same evidence and reproduces the error.

    Args:
        record: The gate session to promote.
        registry_path: The registry file, rewritten in place.
        root: The directory the registry's cited records (and the promotion records) are relative to.
        measure: Measures extensions of red workloads; unused on a rerun (see above).
        checkout: Moves the reference worktree (`DEFAULT_REF_WORKTREE`) to the promoted commit.

    Returns:
        The promoted registry: ``ref_commit`` the session's HEAD, the same harness, workloads and
        ``anchor_ref``, the chain plus the new link, no calibrations.

    Raises:
        MissingRegistryError: No registry at *registry_path*.
        MalformedRatioRecordError: *record* is not JSON, or not a ``perf-ratio/1`` record.
        InconsistentRatioRecordError: As `PerfRatioSession.build`.
        NonPositiveCpuError: As `PerfRatioSession.build`.
        InjectedSessionError: The session is an audit session (`require_gateable`).
        ReferenceMismatchError: As `gate`.
        HarnessMismatchError: As `gate`.
        UncalibratedMachineError: As `gate`.
        UnregisteredWorkloadError: As `gate`.
        IncompleteSessionError: As `gate`.
        UnanchoredWorkloadError: As `gate`.
        ExtensionMismatchError: As `gate`.
        PromotionRecordsExistError: The promotion directory exists with other records than this
            promotion writes (checked before measuring for the gate record).
        PromotionRefusedError: A regression, no IMPROVED workload, a dirty HEAD tree, or HEAD is the reference.
        RegressedAnchorLinkError: As `build_registry` (the new link fails its replay).
        UnimprovedAnchorLinkError: As `build_registry`.
        InvalidAnchorSessionError: As `build_registry`.
        AnchorValueMismatchError: The chain re-derives other anchors than the gate report gives.
        PromotedReferenceMismatchError: The checkout landed on another commit.
        OSError: Staging or renaming the promotion records, or writing the registry, failed.
    """
    registry = _load(registry_path, root)
    data, session = _read_session(record)
    head_rev = session.arms[ArmName.HEAD].rev
    target = root / promotion_dir(registry.ref_commit, head_rev).relative_to(REPOSITORY)
    resuming = target.exists()
    if resuming:
        persisted_gate = target / GATE_RECORD_DIR / SESSION_RECORD_NAME
        if not persisted_gate.is_file() or persisted_gate.read_bytes() != data:
            raise PromotionRecordsExistError(
                f"{target} exists, but its {GATE_RECORD_DIR}/{SESSION_RECORD_NAME} is not the record {record}: it holds another promotion's records. "
                "Move it aside to promote afresh."
            )
        measure = _replay_measure(target)
    report = gate(session, registry, measure, audit=False)
    for line in report_lines(report):
        print(line)
    refusals = promotion_refusals(session, report, registry.ref_commit)
    if refusals:
        raise PromotionRefusedError(f"{record} is not promoted: {'; '.join(refusals)}.")
    moves = anchor_moves(report, anchors_for(registry, tuple(verdict.workload for verdict in report.verdicts)))
    files = _promotion_files(data, report.extensions)
    if resuming:
        persisted = _files_under(target)
        differing = sorted(name for name in persisted.keys() | files.keys() if persisted.get(name) != files.get(name))
        if differing:
            raise PromotionRecordsExistError(f"{target} holds other records than this promotion writes ({', '.join(differing)} differ). Move it aside to promote afresh.")
    else:
        _stage_promotion(files, target, root)
    (calibration,) = [candidate for candidate in registry.calibrations if candidate.cpu_model == session.cpu_model]
    link = AnchorLink(from_ref=registry.ref_commit, to_ref=head_rev, sessions=tuple(_cite(target / name, root) for name in files), calibration=calibration)
    promoted = build_registry(
        ref_commit=head_rev,
        harness_sha256=registry.harness_sha256,
        calibrations=(),
        workloads=registry.workloads,
        anchor_ref=registry.chain_start,
        anchor_links=(*registry.anchor_links, link),
        root=root,
    )
    expected = {move.workload: move.after for move in moves}
    if dict(promoted.anchors_exact()) != expected:
        raise AnchorValueMismatchError(f"the promoted chain re-derives the anchors {dict(promoted.anchors_exact())}, not the gate report's {expected}; nothing was checked out.")
    moved = checkout(DEFAULT_REF_WORKTREE, head_rev)
    if moved != head_rev:
        raise PromotedReferenceMismatchError(
            f"{DEFAULT_REF_WORKTREE} was checked out at {moved}, not the promoted {head_rev}; the registry is unchanged and the promotion records stay in {target}, "
            f"so re-running promote on {record} finishes it."
        )
    write_registry(promoted, registry_path, root=root)
    print(f"reference promoted: {registry.ref_commit} -> {head_rev}; the registry holds no calibrations.")
    print(f"anchors (A_w -> A'_w), records in {target}:")
    for move in moves:
        print(f"  {move.workload}: A_w {float(move.before):.3f} -> {float(move.after):.3f} ({'reset' if move.reset else 'carried'})")
    print(f"next: {CALIBRATE_TASK} --sessions {CALIBRATION_SESSIONS}")
    return promoted


def _read_session(record: Path) -> tuple[bytes, PerfRatioSession]:
    # The one read of a promoted record: its bytes, and the session parsed from those same bytes (as `load_session`).
    data = record.read_bytes()
    try:
        raw = json.loads(data.decode("utf-8"))
    except json.JSONDecodeError as error:
        raise MalformedRatioRecordError(f"{record} is not JSON: {error}") from error
    return data, PerfRatioSession.build(raw)


def _promotion_files(gate_bytes: bytes, extensions: Sequence[PerfRatioSession]) -> dict[str, bytes]:
    # Path inside the promotion directory -> its bytes, in citation order: the gate record verbatim, then each extension.
    files = {f"{GATE_RECORD_DIR}/{SESSION_RECORD_NAME}": gate_bytes}
    for index, extension in enumerate(extensions, start=1):
        files[f"{EXTENSION_DIR_PREFIX}{index}/{SESSION_RECORD_NAME}"] = extension.dumps().encode("utf-8")
    return files


def _files_under(directory: Path) -> dict[str, bytes]:
    return {path.relative_to(directory).as_posix(): path.read_bytes() for path in sorted(directory.rglob("*")) if path.is_file()}


def _stage_promotion(files: Mapping[str, bytes], target: Path, root: Path) -> None:
    # Stage in the gitignored scratch, then one rename: a failure before it leaves no promotion directory.
    scratch = root / GATE_SCRATCH.relative_to(REPOSITORY)
    scratch.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=PROMOTION_STAGING_PREFIX, dir=scratch))
    for name, data in files.items():
        path = staging / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(target)


def _persisted_extensions(target: Path) -> Iterator[PerfRatioSession]:
    index = 1
    path = target / f"{EXTENSION_DIR_PREFIX}{index}" / SESSION_RECORD_NAME
    while path.is_file():
        yield load_session(path)
        index += 1
        path = target / f"{EXTENSION_DIR_PREFIX}{index}" / SESSION_RECORD_NAME


def _replay_measure(target: Path) -> Measure:
    # A rerun's `Measure`: the extensions persisted in *target*, in merge order, each the one the gate requests.
    persisted = _persisted_extensions(target)

    def measure(request: ExtensionRequest) -> PerfRatioSession:
        session = next(persisted, None)
        if session is None:
            raise PromotionRecordsExistError(
                f"{target} holds no further extension record, but the gate requests an extension of {', '.join(request.workloads)}: "
                "these are not this promotion's records. Move it aside to promote afresh."
            )
        measured = tuple(run.workload for run in session.workloads)
        if measured != request.workloads:
            raise PromotionRecordsExistError(
                f"{target}: the next persisted extension measures {list(measured)}, but the gate requests {list(request.workloads)}: "
                "these are not this promotion's records. Move it aside to promote afresh."
            )
        return session

    return measure


def _cite(path: Path, root: Path) -> CitedRecord:
    return CitedRecord(path=path.relative_to(root).as_posix(), record_sha256=sha256_file(path))


def _this_cpu_model() -> str:
    _host, cpu = machine_identity()
    return cpu


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m tools.perf_ratio_gate", description="Gate HEAD-over-reference CPU ratios against A/A-calibrated margins.")
    commands = parser.add_subparsers(dest="command", required=True)
    check_parser = commands.add_parser("check", help="gate a perf-ratio/1 session")
    check_parser.add_argument("record", type=Path, help="the session record")
    check_parser.add_argument("--audit", action="store_true", help="decide an injected or partial session (never a verdict on HEAD)")
    calibrate_parser = commands.add_parser("calibrate", help="calibrate this machine's margins from A/A sessions")
    source = calibrate_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--records", type=Path, nargs="+", metavar="PATH", help="A/A session records of one machine")
    source.add_argument("--sessions", type=int, metavar="K", help=f"measure K >= {CALIBRATION_SESSIONS} A/A sessions at the registry's reference first")
    calibrate_parser.add_argument(
        "--rebase-harness",
        action="store_true",
        help="after a harness edit: take the records' harness and workload set, drop every other calibration (same reference)",
    )
    promote_parser = commands.add_parser("promote", help="move the reference to an improving session's HEAD")
    promote_parser.add_argument("record", type=Path, help="the session record")
    for sub in (check_parser, calibrate_parser, promote_parser):
        sub.add_argument("--registry", type=Path, default=REGISTRY_PATH, help="the registry file")
    return parser


def main(
    argv: Optional[Sequence[str]] = None,
    *,
    measure: Optional[Measure] = None,
    checkout: Optional[Checkout] = None,
    root: Path = REPOSITORY,
    cpu_model: Optional[CpuModelReader] = None,
    tree_harness: Optional[HarnessReader] = None,
) -> int:
    """Command-line entry.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.
        measure: Measures extension and A/A sessions; None runs the harness (`subprocess_measure`
            over `GATE_SCRATCH`).
        checkout: Moves the reference worktree on promotion; None is `checkout_head`.
        root: The directory cited records and calibration directories live under (`REPOSITORY`).
        cpu_model: Reads this machine's CPU model for ``calibrate --sessions``; None reads it as the
            harness does (`machine_identity`).
        tree_harness: Reads the digest of the harness ``calibrate --sessions`` runs; None is
            `harness_sha256` over this tree, the harness the subprocess `Measure` runs.

    Returns:
        The process exit status: 0. Every refusal, and a regression, is raised.
    """
    args = _parser().parse_args(argv)
    chosen_measure = subprocess_measure(GATE_SCRATCH) if measure is None else measure
    registry_path: Path = args.registry
    if args.command == "check":
        return check(args.record, registry_path=registry_path, root=root, measure=chosen_measure, audit=args.audit)
    if args.command == "calibrate":
        if args.records is not None:
            calibrate_records(args.records, registry_path=registry_path, root=root, rebase_harness=args.rebase_harness)
        else:
            calibrate_sessions(
                args.sessions,
                registry_path=registry_path,
                root=root,
                measure=chosen_measure,
                cpu_model=_this_cpu_model if cpu_model is None else cpu_model,
                tree_harness=harness_sha256 if tree_harness is None else tree_harness,
                rebase_harness=args.rebase_harness,
            )
        return 0
    promote(args.record, registry_path=registry_path, root=root, measure=chosen_measure, checkout=checkout_head if checkout is None else checkout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
