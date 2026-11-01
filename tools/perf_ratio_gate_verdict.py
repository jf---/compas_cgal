"""The ratio gate verdict: each workload's ``rho_w`` against its calibrated margin ``m_w``.

`gate` first refuses a session it cannot decide (`require_gateable`), then decides
every workload in session order (`decide`). A workload in the red zone
``rho_w > 1 + m_w`` is not yet a regression: the gate asks `Measure` for more pairs
of that workload alone, merges them into the run (`merge_extension`) and decides
again, until the verdict is not red or the workload holds `MAX_PAIRS` pairs (spec
§The hardware answer: a busy runner costs time, not false reds). The merged ratio
is the record's own ``min head / min ref`` over every pair measured.

Changed output is reported beside the verdict, never refused (spec amendment
2026-11-01): the gate decides ``rho_w`` whether or not the arms' outputs agree.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import Margin
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import Ratio
from tools.perf_ratio_record import WorkloadRatio

# Pairs per extension batch: one default batch (perf_ratio.DEFAULT_PAIRS).
EXTENSION_PAIRS = 5
# Cap on a workload's pairs: 4x the calibrated N, three extensions; bounds a red workload's extra wall
# at 15 pairs (omega monstera, the costliest, ~27 s/pair on M1 Max -> ~7 min).
MAX_PAIRS = 20


class InjectedSessionError(Exception):
    """An audit session (one with an injected spin) is gated outside audit mode."""


class ReferenceMismatchError(Exception):
    """The session's reference commit is not the registry's."""


class HarnessMismatchError(Exception):
    """The session ran another harness than the one the registry was calibrated with."""


class UnregisteredWorkloadError(Exception):
    """A session workload has no margin in the registry."""


class IncompleteSessionError(Exception):
    """Outside audit mode, the session's workload set is not the registry's."""


class ExtensionMismatchError(Exception):
    """An extension session did not measure what its base session and request measured."""


class Verdict(str, Enum):
    """One workload's gate verdict."""

    PASS = "PASS"
    IMPROVED = "IMPROVED"
    REGRESSION = "REGRESSION"


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
    """

    workload: str
    verdict: Verdict
    ratio: Ratio
    margin: Margin
    pairs: int
    changed_output: bool


@dataclass(frozen=True)
class GateReport:
    """Every workload's verdict, in session order.

    Attributes:
        verdicts: One verdict per session workload.
        audit: Whether the session was gated in audit mode.
    """

    verdicts: tuple[WorkloadVerdict, ...]
    audit: bool

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
        """Name each regressed workload of *report* with its ``rho_w``, ``m_w`` and pair count."""
        self.report = report
        named = "; ".join(f"{verdict.workload}: rho_w {verdict.ratio!r} > 1 + m_w {verdict.margin!r} over {verdict.pairs} pairs" for verdict in report.regressions)
        super().__init__(f"ratio regression: {named}")


def decide(run: WorkloadRatio, margin: Margin) -> Verdict:
    """The verdict on *run*'s ``rho_w`` against *margin*; both bounds pass.

    Args:
        run: The workload's pairs.
        margin: The workload's calibrated ``m_w``.

    Returns:
        REGRESSION when ``rho_w > 1 + m_w``, IMPROVED when ``rho_w < 1 - m_w``, else PASS.
    """
    if run.ratio > 1.0 + margin:
        return Verdict.REGRESSION
    if run.ratio < 1.0 - margin:
        return Verdict.IMPROVED
    return Verdict.PASS


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


def merge_extension(base: WorkloadRatio, base_session: PerfRatioSession, extension: PerfRatioSession) -> WorkloadRatio:
    """*base* with *extension*'s pairs of that workload appended, each pair as measured.

    Args:
        base: The workload's run so far.
        base_session: The session *base* comes from.
        extension: The extension session: *base*'s workload alone, measured as *base_session* was.

    Returns:
        The merged run: *base*'s pairs then the extension's, *base*'s name, family and digests.

    Raises:
        ExtensionMismatchError: The extension's arms, harness, CPU model, host, reference or
            injection differ from *base_session*'s, it measures other workloads than *base*'s,
            or its output digests (or family) differ from *base*'s.
    """
    for arm in ArmName:
        if extension.arms[arm] != base_session.arms[arm]:
            raise ExtensionMismatchError(f"extension's {arm.value} arm is {extension.arms[arm]}, not the base session's {base_session.arms[arm]}; it measured another build.")
    for label, measured, expected in (
        ("harness_sha256", extension.harness_sha256, base_session.harness_sha256),
        ("cpu_model", extension.cpu_model, base_session.cpu_model),
        ("host", extension.host, base_session.host),
        ("ref_commit", extension.ref_commit, base_session.ref_commit),
        ("injected", extension.injected, base_session.injected),
    ):
        if measured != expected:
            raise ExtensionMismatchError(f"extension's {label} is {measured!r}, not the base session's {expected!r}.")
    names = tuple(workload.workload for workload in extension.workloads)
    if names != (base.workload,):
        raise ExtensionMismatchError(f"extension measures {list(names)}, not exactly the requested {base.workload!r}.")
    (run,) = extension.workloads
    for label, measured, expected in (
        ("family", run.family, base.family),
        ("ref_digest", run.ref_digest, base.ref_digest),
        ("head_digest", run.head_digest, base.head_digest),
    ):
        if measured != expected:
            raise ExtensionMismatchError(f"extension's {base.workload} {label} is {measured!r}, not the base run's {expected!r}.")
    return WorkloadRatio(workload=base.workload, family=base.family, pairs=base.pairs + run.pairs, ref_digest=base.ref_digest, head_digest=base.head_digest)


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
        ExtensionMismatchError: An extension is refused by `merge_extension`, or holds another
            pair count than requested.
    """
    margins = require_gateable(session, registry, audit=audit)
    verdicts = []
    for base in session.workloads:
        margin = margins[base.workload]
        run = base
        verdict = decide(run, margin)
        while verdict is Verdict.REGRESSION and len(run.pairs) < MAX_PAIRS:
            request = ExtensionRequest(
                ref_commit=session.ref_commit,
                head_rev=session.arms[ArmName.HEAD].rev,
                workloads=(base.workload,),
                pairs=min(EXTENSION_PAIRS, MAX_PAIRS - len(run.pairs)),
                injected=session.injected,
            )
            merged = merge_extension(run, session, measure(request))
            added = len(merged.pairs) - len(run.pairs)
            if added != request.pairs:
                raise ExtensionMismatchError(f"extension of {base.workload} holds {added} pairs, not the requested {request.pairs}.")
            run = merged
            verdict = decide(run, margin)
        verdicts.append(WorkloadVerdict(workload=base.workload, verdict=verdict, ratio=run.ratio, margin=margin, pairs=len(run.pairs), changed_output=not run.same_output))
    return GateReport(verdicts=tuple(verdicts), audit=audit)
