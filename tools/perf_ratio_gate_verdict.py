"""The ratio gate verdict: each workload's anchored ratio ``A_w * rho_w`` against its calibrated margin ``m_w``.

`gate` first refuses a session it cannot decide (`require_gateable`), then decides
every workload in session order (`decide`) against its anchor ``A_w`` (spec
amendment 2026-11-01, the anchor chain; `registry.anchors_exact`). A workload in the
red zone ``A_w * rho_w > 1 + m_w`` is not yet a regression: the gate asks `Measure`
for more pairs of that workload alone, merges them into the run (`merge_extension`)
and decides again, until the verdict is not red or the workload holds `MAX_PAIRS`
pairs (spec §The hardware answer: a busy runner costs time, not false reds). The
merged ratio is the record's own ``min head / min ref`` over every pair measured.
Every accepted extension session is kept on the report, in merge order, so a
promotion can persist the evidence its anchors are derived from.

`decide`, `Verdict` and `MAX_PAIRS` live in `tools.perf_ratio_gate_decide`, below this
module and the registry, so the registry's chain replay re-decides a promoted link
with the gate's own function; they are re-exported here. Decisions are exact.

Changed output is reported beside the verdict, never refused (spec amendment
2026-11-01): the gate decides ``rho_w`` whether or not the arms' outputs agree.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Optional

from tools.perf_ratio_gate_decide import MAX_PAIRS as MAX_PAIRS
from tools.perf_ratio_gate_decide import UNIT_ANCHOR
from tools.perf_ratio_gate_decide import Anchor
from tools.perf_ratio_gate_decide import Margin
from tools.perf_ratio_gate_decide import Verdict as Verdict
from tools.perf_ratio_gate_decide import anchored
from tools.perf_ratio_gate_decide import decide as decide
from tools.perf_ratio_gate_errors import GateRefusal
from tools.perf_ratio_gate_merge import ExtensionMismatchError
from tools.perf_ratio_gate_merge import merge_extension
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import Ratio

# Pairs per extension batch: one default batch (perf_ratio.DEFAULT_PAIRS).
EXTENSION_PAIRS = 5


class InjectedSessionError(GateRefusal):
    """An audit session (one with an injected spin) is gated outside audit mode."""


class ReferenceMismatchError(GateRefusal):
    """The session's reference commit is not the registry's."""


class HarnessMismatchError(GateRefusal):
    """The session ran another harness than the one the registry was calibrated with."""


class UnregisteredWorkloadError(GateRefusal):
    """A session workload has no margin in the registry."""


class IncompleteSessionError(GateRefusal):
    """Outside audit mode, the session's workload set is not the registry's."""


class UnanchoredWorkloadError(GateRefusal):
    """A gated workload has no anchor ``A_w`` in a registry whose anchor chain is not empty."""


@dataclass(frozen=True)
class ExtensionRequest:
    """More pairs of the named workloads, measured exactly as the base session was.

    Attributes:
        ref_commit: The base session's reference commit.
        head_rev: The base session's HEAD arm rev.
        workloads: The workloads to measure (one: the red workload).
        pairs: Interleaved pairs per workload.
        injected: The base session's audit injection, replayed; None for an honest session.
    """

    ref_commit: str
    head_rev: str
    workloads: tuple[str, ...]
    pairs: int
    injected: Optional[Injection]


# Measures an extension: runs the harness for the request and returns its session.
Measure = Callable[[ExtensionRequest], PerfRatioSession]


@dataclass(frozen=True)
class WorkloadVerdict:
    """One workload's final verdict.

    Attributes:
        workload: The workload's name.
        verdict: The decision on the final (possibly extended) run.
        ratio: The final run's ``rho_w``.
        margin: The workload's calibrated ``m_w``.
        pairs: The pairs the final run holds.
        changed_output: Whether the arms' output digests differ (reported, not refused).
        anchor: The workload's anchor ``A_w`` as the nearest double (reporting; decided exactly).
        cumulative_ratio: ``A_w * rho_w`` as the nearest double to the exact product (reporting).
    """

    workload: str
    verdict: Verdict
    ratio: Ratio
    margin: Margin
    pairs: int
    changed_output: bool
    anchor: float
    cumulative_ratio: float


@dataclass(frozen=True)
class GateReport:
    """Every workload's verdict, in session order.

    Attributes:
        verdicts: One verdict per session workload.
        audit: Whether the session was gated in audit mode.
        extensions: Every extension session the gate merged, in merge order (workload by workload,
            in session order); a promotion persists them beside the gate session.
    """

    verdicts: tuple[WorkloadVerdict, ...]
    audit: bool
    extensions: tuple[PerfRatioSession, ...]

    @property
    def regressions(self) -> tuple[WorkloadVerdict, ...]:
        """The verdicts that are `Verdict.REGRESSION`."""
        return tuple(verdict for verdict in self.verdicts if verdict.verdict is Verdict.REGRESSION)

    @property
    def improvements(self) -> tuple[WorkloadVerdict, ...]:
        """The verdicts that are `Verdict.IMPROVED`."""
        return tuple(verdict for verdict in self.verdicts if verdict.verdict is Verdict.IMPROVED)


class RatioRegressionError(Exception):
    """A gated session regressed: some workload stayed above ``1 + m_w`` after extension.

    Attributes:
        report: The gate report whose regressions this error names.
    """

    def __init__(self, report: GateReport) -> None:
        """Name each regressed workload of *report* with its ``rho_w`` (anchored, if ``A_w != 1``), ``m_w`` and pair count."""
        self.report = report
        named = "; ".join(_regression_line(verdict) for verdict in report.regressions)
        super().__init__(f"ratio regression: {named}")


def _regression_line(verdict: WorkloadVerdict) -> str:
    """One regressed workload in one form, ``A_w = 1`` included: ``A_w*rho_w``, ``1 + m_w``, ``A_w``, ``rho_w``, pairs (reporting floats)."""
    return f"{verdict.workload}: A_w*rho_w {verdict.cumulative_ratio!r} > 1 + m_w {verdict.margin!r} (A_w {verdict.anchor!r}, rho_w {verdict.ratio!r}) over {verdict.pairs} pairs"


def anchors_for(registry: GateRegistry, workloads: tuple[str, ...]) -> dict[str, Anchor]:
    """The exact anchor ``A_w`` of each of *workloads*.

    A workload the registry holds no anchor for is at `UNIT_ANCHOR` only when the chain is empty
    (a raw registry with no links: nothing was ever carried). On a non-empty chain every gated
    workload must carry its derived anchor (`derive_anchors` covers the registry's workloads).

    Args:
        registry: The gate registry.
        workloads: The workloads to anchor.

    Returns:
        Workload -> its exact anchor.

    Raises:
        UnanchoredWorkloadError: The chain is not empty and a workload has no anchor.
    """
    exact = registry.anchors_exact()
    if registry.anchor_links:
        missing = [name for name in workloads if name not in exact]
        if missing:
            raise UnanchoredWorkloadError(f"no anchor for {', '.join(missing)} in a registry with {len(registry.anchor_links)} anchor link(s); anchored: {sorted(exact)}.")
    return {name: exact.get(name, UNIT_ANCHOR) for name in workloads}


def require_gateable(session: PerfRatioSession, registry: GateRegistry, *, audit: bool) -> Mapping[str, Margin]:
    """The margin of each of *session*'s workloads, or the reason *session* cannot be gated.

    Args:
        session: The session to gate.
        registry: The gate registry.
        audit: Audit mode: an injected session and a partial workload set are allowed.

    Returns:
        Session workload -> its margin on the session's machine.

    Raises:
        InjectedSessionError: *session* is injected and *audit* is false.
        ReferenceMismatchError: *session*'s reference is not the registry's.
        HarnessMismatchError: *session*'s harness digest is not the registry's.
        UncalibratedMachineError: The registry holds no calibration for *session*'s CPU model.
        UnregisteredWorkloadError: A session workload has no margin.
        IncompleteSessionError: *audit* is false and the workload set is not the registry's.
    """
    if session.injected is not None and not audit:
        raise InjectedSessionError(
            f"session is an audit session (injected {session.injected.seconds} s on {session.injected.arm.value}); gate it with audit mode, never as a verdict on HEAD."
        )
    if session.ref_commit != registry.ref_commit:
        raise ReferenceMismatchError(f"session measured against reference {session.ref_commit}, not the registry's {registry.ref_commit}.")
    if session.harness_sha256 != registry.harness_sha256:
        raise HarnessMismatchError(
            f"session ran harness {session.harness_sha256}, not the registry's calibrated {registry.harness_sha256}; "
            f"every margin is void under an edited harness. Re-run the A/A sessions at the reference {registry.ref_commit} under this harness, "
            "then pixi run perf-ratio-calibrate --rebase-harness."
        )
    calibrated = registry.margins_for(session.cpu_model)
    names = [workload.workload for workload in session.workloads]
    unregistered = [name for name in names if name not in calibrated]
    if unregistered:
        raise UnregisteredWorkloadError(f"no margin for {', '.join(unregistered)} on {session.cpu_model!r}; registered: {sorted(calibrated)}.")
    if not audit and frozenset(names) != frozenset(registry.workloads):
        missing = sorted(frozenset(registry.workloads) - frozenset(names))
        raise IncompleteSessionError(
            f"session measures {sorted(names)}, not the registry's {list(registry.workloads)} (missing {missing}); a gate session measures every workload."
        )
    return {name: calibrated[name] for name in names}


def gate(session: PerfRatioSession, registry: GateRegistry, measure: Measure, *, audit: bool = False) -> GateReport:
    """Decide every workload of *session*, extending each red one up to `MAX_PAIRS` pairs.

    Args:
        session: The gate session.
        registry: The gate registry.
        measure: Measures an extension request; called only for red workloads.
        audit: Audit mode (see `require_gateable`).

    Returns:
        The report, one verdict per session workload in session order.

    Raises:
        InjectedSessionError: As `require_gateable`.
        ReferenceMismatchError: As `require_gateable`.
        HarnessMismatchError: As `require_gateable`.
        UncalibratedMachineError: As `require_gateable`.
        UnregisteredWorkloadError: As `require_gateable`.
        IncompleteSessionError: As `require_gateable`.
        UnanchoredWorkloadError: As `anchors_for`.
        ExtensionMismatchError: An extension is refused by `merge_extension`, or holds another
            pair count than requested.
    """
    margins = require_gateable(session, registry, audit=audit)
    anchors = anchors_for(registry, tuple(margins))
    verdicts = []
    extensions: list[PerfRatioSession] = []
    for base in session.workloads:
        margin = margins[base.workload]
        anchor = anchors[base.workload]
        run = base
        verdict = decide(run, margin, anchor)
        while verdict is Verdict.REGRESSION and len(run.pairs) < MAX_PAIRS:
            request = ExtensionRequest(
                ref_commit=session.ref_commit,
                head_rev=session.arms[ArmName.HEAD].rev,
                workloads=(base.workload,),
                pairs=min(EXTENSION_PAIRS, MAX_PAIRS - len(run.pairs)),
                injected=session.injected,
            )
            extension = measure(request)
            merged = merge_extension(run, session, extension)
            added = len(merged.pairs) - len(run.pairs)
            if added != request.pairs:
                raise ExtensionMismatchError(f"extension of {base.workload} holds {added} pairs, not the requested {request.pairs}.")
            extensions.append(extension)
            run = merged
            verdict = decide(run, margin, anchor)
        verdicts.append(
            WorkloadVerdict(
                workload=base.workload,
                verdict=verdict,
                ratio=run.ratio,
                margin=margin,
                pairs=len(run.pairs),
                changed_output=not run.same_output,
                anchor=float(anchor),
                cumulative_ratio=float(anchored(run, anchor)),
            )
        )
    return GateReport(verdicts=tuple(verdicts), audit=audit, extensions=tuple(extensions))
