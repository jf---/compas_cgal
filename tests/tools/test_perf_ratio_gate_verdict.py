"""The ratio gate verdict: decide rho_w against m_w, extend a red workload, refuse what is not gateable.

The registry in these tests is a raw `GateRegistry` with chosen margins and anchors (the
bypass-safe constructor): verdict logic does not depend on how either was derived.
"""

from __future__ import annotations

import ast
import dataclasses
import math
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path
from typing import Optional

import pytest
from tools.perf_ratio import DEFAULT_PAIRS
from tools import perf_ratio_gate_decide
from tools import perf_ratio_gate_verdict
from tools.perf_ratio_gate_decide import exact_ratio
from tools.perf_ratio_gate_merge import ExtensionMismatchError
from tools.perf_ratio_gate_merge import merge_extension
from tools.perf_ratio_gate_registry import SAFETY
from tools.perf_ratio_gate_registry import UNIT_ANCHOR
from tools.perf_ratio_gate_registry import Anchor
from tools.perf_ratio_gate_registry import AnchorLink
from tools.perf_ratio_gate_registry import Calibration
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import Margin
from tools.perf_ratio_gate_registry import UncalibratedMachineError
from tools.perf_ratio_gate_verdict import EXTENSION_PAIRS
from tools.perf_ratio_gate_verdict import MAX_PAIRS
from tools.perf_ratio_gate_verdict import ExtensionRequest
from tools.perf_ratio_gate_verdict import GateReport
from tools.perf_ratio_gate_verdict import HarnessMismatchError
from tools.perf_ratio_gate_verdict import IncompleteSessionError
from tools.perf_ratio_gate_verdict import InjectedSessionError
from tools.perf_ratio_gate_verdict import RatioRegressionError
from tools.perf_ratio_gate_verdict import ReferenceMismatchError
from tools.perf_ratio_gate_verdict import UnanchoredWorkloadError
from tools.perf_ratio_gate_verdict import UnregisteredWorkloadError
from tools.perf_ratio_gate_verdict import Verdict
from tools.perf_ratio_gate_verdict import WorkloadVerdict
from tools.perf_ratio_gate_verdict import decide
from tools.perf_ratio_gate_verdict import gate
from tools.perf_ratio_gate_verdict import require_gateable
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import PerfRatioSession

from tests.tools.perf_ratio_sessions import CPU_MODEL
from tests.tools.perf_ratio_sessions import HARNESS
from tests.tools.perf_ratio_sessions import OTHER_REV
from tests.tools.perf_ratio_sessions import REF
from tests.tools.perf_ratio_sessions import synthetic_session

GT = "gt/square"
OMEGA = "omega/monstera"
WORKLOADS = (GT, OMEGA)
# A dimensionless margin; 1 + M and 1 - M are the decided bounds. 0.04 is not dyadic, so 1 + M is a rounded
# double, but the boundary test is exact: its ratio is (1.0 + M) / 1.0, the same float decide builds.
M = 0.04
# The HEAD rev of a session that measures a change (not an A/A session).
HEAD = "fedcba9876543210fedcba9876543210fedcba98"
# One audit injection on the HEAD arm.
INJECTION = Injection(arm=ArmName.HEAD, seconds=CpuSeconds(0.5))


def registry(
    margins: Optional[Mapping[str, float]] = None,
    anchors: Optional[Mapping[str, Fraction]] = None,
    links: tuple[AnchorLink, ...] = (),
) -> GateRegistry:
    """A raw registry at REF under HARNESS, CPU_MODEL calibrated with *margins* (M on every workload by default).

    *anchors* are the exact ``A_w`` (none: an empty chain); *links* the chain they came from (raw, never replayed here).
    """
    chosen = {name: M for name in WORKLOADS} if margins is None else margins
    return GateRegistry(
        ref_commit=REF,
        harness_sha256=HARNESS,
        safety=SAFETY,
        workloads=tuple(sorted(chosen)),
        calibrations=(Calibration(cpu_model=CPU_MODEL, sessions=(), margins={name: Margin(value) for name, value in chosen.items()}),),
        anchor_links=links,
        exact_anchors={} if anchors is None else {name: Anchor(value) for name, value in anchors.items()},
    )


def pairs(ref_cpu: float, head_cpu: float, count: int = EXTENSION_PAIRS) -> list[tuple[float, float]]:
    """*count* identical (ref, head) pairs."""
    return [(ref_cpu, head_cpu)] * count


def session(ratios: Mapping[str, Sequence[tuple[float, float]]], **overrides: object) -> PerfRatioSession:
    """A session at REF -> HEAD; *overrides* go to `synthetic_session`."""
    options: dict[str, object] = {"head_rev": HEAD}
    options.update(overrides)
    return synthetic_session(ratios, **options)  # type: ignore[arg-type]


def passing(**overrides: object) -> PerfRatioSession:
    """Every registered workload at ratio 1.0."""
    return session({name: pairs(1.0, 1.0) for name in WORKLOADS}, **overrides)


class FakeMeasure:
    """A `Measure` that records each request and answers it with the next batch's (ref, head) pair."""

    def __init__(self, batches: Sequence[tuple[float, float]], **overrides: object) -> None:
        self.requests: list[ExtensionRequest] = []
        self.returned: list[PerfRatioSession] = []
        self.batches = list(batches)
        self.overrides = overrides

    def __call__(self, request: ExtensionRequest) -> PerfRatioSession:
        self.requests.append(request)
        ref_cpu, head_cpu = self.batches[len(self.requests) - 1]
        options: dict[str, object] = {"ref_commit": request.ref_commit, "head_rev": request.head_rev, "injected": request.injected}
        options.update(self.overrides)
        extension = synthetic_session({name: pairs(ref_cpu, head_cpu, request.pairs) for name in request.workloads}, **options)  # type: ignore[arg-type]
        self.returned.append(extension)
        return extension


def no_measure(request: ExtensionRequest) -> PerfRatioSession:
    raise AssertionError(f"no extension was expected, got {request}")


def verdict_of(report: GateReport, workload: str) -> WorkloadVerdict:
    (found,) = [verdict for verdict in report.verdicts if verdict.workload == workload]
    return found


def test_an_extension_batch_is_one_default_batch() -> None:
    assert EXTENSION_PAIRS == DEFAULT_PAIRS
    assert MAX_PAIRS == 4 * DEFAULT_PAIRS


# --- decide ---------------------------------------------------------------------------------------


def test_decide_boundaries_are_exact() -> None:
    # The bounds are the exact rationals 1 +- m, m the double M is. 0.04 is not dyadic, so the doubles
    # 1.0 + M and 1.0 - M round away from them: 1.0 + M lies 5 / 2**57 above the exact upper bound and
    # 1.0 - M as far below the exact lower one. An exact gate decides them REGRESSION and IMPROVED; the
    # ulp next to each, on the inside, is the last ratio that passes.
    margin = Margin(M)
    upper = 1.0 + M
    lower = 1.0 - M
    assert Fraction(upper) - (1 + Fraction(M)) == Fraction(5, 2**57)
    assert Fraction(lower) - (1 - Fraction(M)) == -Fraction(5, 2**57)

    def decided(head_cpu: float) -> Verdict:
        run = session({GT: [(1.0, head_cpu)]}).workloads[0]
        assert run.ratio == head_cpu
        return decide(run, margin, UNIT_ANCHOR)

    assert decided(math.nextafter(upper, 0.0)) is Verdict.PASS
    assert decided(upper) is Verdict.REGRESSION
    assert decided(math.nextafter(lower, math.inf)) is Verdict.PASS
    assert decided(lower) is Verdict.IMPROVED
    assert decided(1.0) is Verdict.PASS


# A double anchor and two adjacent double ratios straddling the exact bound 1 + 0.02 (the double 0.02),
# where IEEE products and the exact Fraction products disagree (both checked in the test):
#   A * RHO_ABOVE:  exact 3.0e-17 above the bound (REGRESSION); the float product rounds onto 1.0 + 0.02 (PASS).
#   A * RHO_BELOW:  exact 8.4e-17 below the bound (PASS); the float product, read exactly, lies above it.
BOUNDARY_ANCHOR = 1.03
BOUNDARY_MARGIN = 0.02
RHO_ABOVE = 0.9902912621359223
RHO_BELOW = 0.9902912621359222


def test_decide_on_the_exact_anchored_product_not_a_float_one() -> None:
    anchor = Anchor(Fraction(BOUNDARY_ANCHOR))
    margin = Margin(BOUNDARY_MARGIN)
    bound = 1 + Fraction(BOUNDARY_MARGIN)
    assert RHO_BELOW == math.nextafter(RHO_ABOVE, 0.0)
    # Above: exact says over the bound, the float product says not.
    assert anchor * Fraction(RHO_ABOVE) > bound
    assert not BOUNDARY_ANCHOR * RHO_ABOVE > 1.0 + BOUNDARY_MARGIN
    # Below: exact says within, the float product (even compared exactly) says over.
    assert anchor * Fraction(RHO_BELOW) < bound
    assert Fraction(BOUNDARY_ANCHOR * RHO_BELOW) > bound

    def decided(ratio: float) -> Verdict:
        run = session({GT: [(1.0, ratio)]}).workloads[0]
        assert run.ratio == ratio
        return decide(run, margin, anchor)

    assert decided(RHO_ABOVE) is Verdict.REGRESSION
    assert decided(RHO_BELOW) is Verdict.PASS


def test_an_anchor_turns_a_flat_ratio_into_a_regression() -> None:
    run = session({GT: pairs(1.0, 1.0)}).workloads[0]
    assert decide(run, Margin(0.02), UNIT_ANCHOR) is Verdict.PASS
    assert decide(run, Margin(0.02), Anchor(Fraction(103, 100))) is Verdict.REGRESSION
    assert decide(run, Margin(0.03), Anchor(Fraction(96, 100))) is Verdict.IMPROVED


# --- gate: the sequential extension -----------------------------------------------------------------


def test_a_passing_session_needs_no_measurement() -> None:
    measure = FakeMeasure([])
    report = gate(passing(), registry(), measure)
    assert measure.requests == []
    assert [verdict.verdict for verdict in report.verdicts] == [Verdict.PASS, Verdict.PASS]
    assert [verdict.workload for verdict in report.verdicts] == list(WORKLOADS)
    assert report.regressions == ()
    assert report.improvements == ()
    assert report.audit is False


def test_a_noisy_red_workload_extends_and_passes() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.0)})
    assert base.workloads[0].ratio == 1.10
    # The extension's HEAD ran faster: merged rho = min head / min ref = 1.01 / 1.0.
    measure = FakeMeasure([(1.0, 1.01)])
    report = gate(base, registry(), measure)
    verdict = verdict_of(report, GT)
    assert verdict.verdict is Verdict.PASS
    assert verdict.ratio == 1.01
    assert verdict.pairs == 10
    assert len(measure.requests) == 1
    request = measure.requests[0]
    assert request.workloads == (GT,)
    assert request.pairs == 5
    assert request.ref_commit == REF
    assert request.head_rev == HEAD
    assert request.injected is None
    assert verdict_of(report, OMEGA).pairs == 5


def test_a_persistent_regression_stops_at_max_pairs() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.0)})
    measure = FakeMeasure([(1.0, 1.10)] * 10)
    report = gate(base, registry(), measure)
    verdict = verdict_of(report, GT)
    assert verdict.verdict is Verdict.REGRESSION
    assert verdict.pairs == MAX_PAIRS
    assert len(measure.requests) == 3
    assert [request.pairs for request in measure.requests] == [5, 5, 5]
    assert report.regressions == (verdict,)


def test_the_last_extension_is_cut_to_the_cap() -> None:
    base = session({GT: pairs(1.0, 1.10, 7), OMEGA: pairs(1.0, 1.0)})
    measure = FakeMeasure([(1.0, 1.10)] * 10)
    verdict = verdict_of(gate(base, registry(), measure), GT)
    assert [request.pairs for request in measure.requests] == [5, 5, 3]
    assert verdict.pairs == MAX_PAIRS


def test_an_improvement_is_not_extended() -> None:
    base = session({GT: pairs(1.0, 0.90), OMEGA: pairs(1.0, 1.0)})
    report = gate(base, registry(), no_measure)
    assert verdict_of(report, GT).verdict is Verdict.IMPROVED
    assert report.improvements == (verdict_of(report, GT),)


def test_the_merged_run_keeps_every_pair_as_measured() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.0)})
    extension = session({GT: pairs(1.0, 1.01)})
    merged = merge_extension(base.workloads[0], base, extension)
    assert merged.pairs == base.workloads[0].pairs + extension.workloads[0].pairs
    assert merged.ref_digest == base.workloads[0].ref_digest
    assert merged.head_digest == base.workloads[0].head_digest
    assert merged.workload == GT
    assert merged.family == base.workloads[0].family


def test_the_extension_replays_the_injection() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.0)}, injected=INJECTION)
    measure = FakeMeasure([(1.0, 1.10)] * 3)
    report = gate(base, registry(), measure, audit=True)
    assert report.audit is True
    assert measure.requests
    assert all(request.injected == INJECTION for request in measure.requests)
    assert verdict_of(report, GT).verdict is Verdict.REGRESSION


# --- gate: the anchor chain -------------------------------------------------------------------------


def test_an_anchored_regression_extends_and_an_extension_to_098_passes() -> None:
    anchor = Fraction(103, 100)
    chosen = registry({GT: 0.02, OMEGA: 0.02}, anchors={GT: anchor, OMEGA: Fraction(1)})
    base = session({GT: pairs(1.0, 1.0), OMEGA: pairs(1.0, 1.0)})
    # The extension's HEAD ran faster: merged rho = 0.98 / 1.0; 1.03 * 0.98 = 1.0094 <= 1.02.
    measure = FakeMeasure([(1.0, 0.98)])
    report = gate(base, chosen, measure)
    verdict = verdict_of(report, GT)
    assert [request.workloads for request in measure.requests] == [(GT,)]
    assert verdict.verdict is Verdict.PASS
    assert verdict.ratio == 0.98
    assert verdict.pairs == 10
    assert verdict.anchor == 1.03
    assert verdict.cumulative_ratio == float(anchor * Fraction(0.98))
    assert report.extensions == tuple(measure.returned)
    assert verdict_of(report, OMEGA).verdict is Verdict.PASS
    assert verdict_of(report, OMEGA).pairs == 5


def test_an_anchored_regression_that_persists_stays_red() -> None:
    chosen = registry({GT: 0.02, OMEGA: 0.02}, anchors={GT: Fraction(103, 100), OMEGA: Fraction(1)})
    base = session({GT: pairs(1.0, 1.0), OMEGA: pairs(1.0, 1.0)})
    measure = FakeMeasure([(1.0, 1.0)] * 3)
    report = gate(base, chosen, measure)
    verdict = verdict_of(report, GT)
    assert verdict.verdict is Verdict.REGRESSION
    assert verdict.pairs == MAX_PAIRS
    assert verdict.ratio == 1.0
    assert verdict.cumulative_ratio == 1.03
    assert report.regressions == (verdict,)
    message = str(RatioRegressionError(report))
    assert "A_w 1.03" in message
    assert "A_w*rho_w 1.03" in message


def test_an_anchor_below_one_turns_a_flat_ratio_into_an_improvement() -> None:
    chosen = registry({GT: 0.03, OMEGA: 0.03}, anchors={GT: Fraction(96, 100), OMEGA: Fraction(1)})
    report = gate(passing(), chosen, no_measure)
    assert verdict_of(report, GT).verdict is Verdict.IMPROVED
    assert verdict_of(report, GT).cumulative_ratio == 0.96
    assert verdict_of(report, OMEGA).verdict is Verdict.PASS
    assert report.improvements == (verdict_of(report, GT),)
    assert report.extensions == ()


def test_an_empty_chain_keeps_the_unanchored_verdicts() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 0.90)})
    measure = FakeMeasure([(1.0, 1.01)])
    report = gate(base, registry(), measure)
    assert [(verdict.verdict, verdict.anchor) for verdict in report.verdicts] == [(Verdict.PASS, 1.0), (Verdict.IMPROVED, 1.0)]
    assert [verdict.cumulative_ratio for verdict in report.verdicts] == [verdict.ratio for verdict in report.verdicts]
    assert report.extensions == tuple(measure.returned)


def test_every_extension_is_kept_on_the_report_in_merge_order() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.10)})
    # GT: two red batches then a passing one; OMEGA: one passing batch.
    measure = FakeMeasure([(1.0, 1.10), (1.0, 1.10), (1.0, 1.01), (1.0, 1.01)])
    report = gate(base, registry(), measure)
    assert [request.workloads for request in measure.requests] == [(GT,), (GT,), (GT,), (OMEGA,)]
    assert len(report.extensions) == 4
    assert all(kept is returned for kept, returned in zip(report.extensions, measure.returned))
    assert [verdict.verdict for verdict in report.verdicts] == [Verdict.PASS, Verdict.PASS]


def test_a_workload_without_an_anchor_on_a_nonempty_chain_is_refused() -> None:
    link = AnchorLink(from_ref=OTHER_REV, to_ref=REF, sessions=(), calibration=Calibration(cpu_model=CPU_MODEL, sessions=(), margins={}))
    chained = registry(anchors={GT: Fraction(1)}, links=(link,))
    with pytest.raises(UnanchoredWorkloadError, match=OMEGA):
        gate(passing(), chained, no_measure)


# --- merge_extension refusals -----------------------------------------------------------------------


def red_base(**overrides: object) -> PerfRatioSession:
    return session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.0)}, **overrides)


def family_changed(extension: PerfRatioSession) -> PerfRatioSession:
    """*extension* with its one workload's family renamed (name and digests untouched)."""
    (run,) = extension.workloads
    return dataclasses.replace(extension, workloads=(dataclasses.replace(run, family="other"),))


def rebuilt(arm: ArmName, measure: FakeMeasure) -> Callable[[ExtensionRequest], PerfRatioSession]:
    """*measure* with *arm* alone rebuilt: its ``native_sha256`` changed, its rev and the other arm untouched."""

    def measured(request: ExtensionRequest) -> PerfRatioSession:
        extension = measure(request)
        arms = dict(extension.arms)
        arms[arm] = dataclasses.replace(arms[arm], native_sha256="m" * 64)
        return dataclasses.replace(extension, arms=arms)

    return measured


def test_an_extension_from_another_build_is_refused() -> None:
    base = red_base()
    # The realistic case: HEAD rebuilt at the same rev, REF unchanged.
    with pytest.raises(ExtensionMismatchError, match="head arm"):
        gate(base, registry(), rebuilt(ArmName.HEAD, FakeMeasure([(1.0, 1.01)])))
    with pytest.raises(ExtensionMismatchError, match="ref arm"):
        gate(base, registry(), rebuilt(ArmName.REF, FakeMeasure([(1.0, 1.01)])))
    with pytest.raises(ExtensionMismatchError, match="injected"):
        gate(base, registry(), FakeMeasure([(1.0, 1.01)], injected=INJECTION))


@pytest.mark.parametrize(
    ("field", "extension"),
    [
        ("harness", lambda: session({GT: pairs(1.0, 1.01)}, harness="g" * 64)),
        ("cpu_model", lambda: session({GT: pairs(1.0, 1.01)}, cpu_model="Apple M2")),
        ("host", lambda: dataclasses.replace(session({GT: pairs(1.0, 1.01)}), host="other-host")),
        # ref_commit alone (arms untouched), so the arm check cannot mask it.
        ("ref_commit", lambda: dataclasses.replace(session({GT: pairs(1.0, 1.01)}), ref_commit=OTHER_REV)),
        ("head arm", lambda: session({GT: pairs(1.0, 1.01)}, head_rev=OTHER_REV)),
        ("head_digest", lambda: session({GT: pairs(1.0, 1.01)}, head_digest="e")),
        ("ref_digest", lambda: session({GT: pairs(1.0, 1.01)}, ref_digest="e")),
        # Same family and digests as GT, so only the workload-name check can refuse it.
        ("workload", lambda: session({"gt/other": pairs(1.0, 1.01)})),
        ("workload", lambda: session({GT: pairs(1.0, 1.01), OMEGA: pairs(1.0, 1.0)})),
        ("family", lambda: family_changed(session({GT: pairs(1.0, 1.01)}))),
    ],
    ids=["harness", "cpu-model", "host", "ref-commit", "head-arm", "head-digest", "ref-digest", "other-workload", "extra-workload", "family"],
)
def test_a_mismatched_extension_is_refused(field: str, extension: object) -> None:
    base = red_base()
    with pytest.raises(ExtensionMismatchError):
        merge_extension(base.workloads[0], base, extension())  # type: ignore[operator]


def test_an_extension_with_another_pair_count_is_refused() -> None:
    base = red_base()

    def short(request: ExtensionRequest) -> PerfRatioSession:
        return session({GT: pairs(1.0, 1.01, request.pairs - 1)})

    with pytest.raises(ExtensionMismatchError, match="pairs"):
        gate(base, registry(), short)


# --- require_gateable refusals ----------------------------------------------------------------------


def test_a_session_from_another_harness_is_refused() -> None:
    with pytest.raises(HarnessMismatchError):
        gate(passing(harness="g" * 64), registry(), no_measure)


def test_a_harness_mismatch_names_the_recovery() -> None:
    with pytest.raises(HarnessMismatchError, match=r"A/A.*at the reference.*perf-ratio-calibrate --rebase-harness"):
        gate(passing(harness="g" * 64), registry(), no_measure)


def test_a_session_against_another_reference_is_refused() -> None:
    with pytest.raises(ReferenceMismatchError):
        gate(passing(ref_commit=OTHER_REV), registry(), no_measure)


def test_an_injected_session_is_refused_outside_audit() -> None:
    with pytest.raises(InjectedSessionError):
        gate(passing(injected=INJECTION), registry(), no_measure)
    report = gate(passing(injected=INJECTION), registry(), no_measure, audit=True)
    assert [verdict.verdict for verdict in report.verdicts] == [Verdict.PASS, Verdict.PASS]


def test_an_uncalibrated_machine_is_refused() -> None:
    with pytest.raises(UncalibratedMachineError):
        gate(passing(cpu_model="Apple M2"), registry(), no_measure)


def test_a_partial_session_is_refused_outside_audit() -> None:
    partial = session({GT: pairs(1.0, 1.0)})
    with pytest.raises(IncompleteSessionError):
        gate(partial, registry(), no_measure)
    report = gate(partial, registry(), no_measure, audit=True)
    assert [(verdict.workload, verdict.verdict) for verdict in report.verdicts] == [(GT, Verdict.PASS)]


def test_an_unregistered_workload_is_refused() -> None:
    extra = session({GT: pairs(1.0, 1.0), OMEGA: pairs(1.0, 1.0), "gt/new": pairs(1.0, 1.0)})
    with pytest.raises(UnregisteredWorkloadError, match="gt/new"):
        gate(extra, registry(), no_measure)
    with pytest.raises(UnregisteredWorkloadError):
        gate(extra, registry(), no_measure, audit=True)


def test_refusals_are_checked_in_order() -> None:
    # Injected, at another reference, another harness, an uncalibrated machine: the first check wins.
    worst = passing(injected=INJECTION, ref_commit=OTHER_REV, harness="g" * 64, cpu_model="Apple M2")
    with pytest.raises(InjectedSessionError):
        require_gateable(worst, registry(), audit=False)
    with pytest.raises(ReferenceMismatchError):
        require_gateable(worst, registry(), audit=True)
    with pytest.raises(HarnessMismatchError):
        require_gateable(dataclasses.replace(worst, ref_commit=REF), registry(), audit=True)
    with pytest.raises(UncalibratedMachineError):
        require_gateable(passing(cpu_model="Apple M2", harness=HARNESS), registry(), audit=True)


def test_require_gateable_returns_the_session_workloads_margins() -> None:
    margins = require_gateable(session({GT: pairs(1.0, 1.0)}), registry({GT: 0.04, OMEGA: 0.25}), audit=True)
    assert dict(margins) == {GT: 0.04}


# --- changed output and the regression error ---------------------------------------------------------


def test_changed_output_is_reported_not_refused() -> None:
    changed = session({GT: pairs(1.0, 1.0), OMEGA: pairs(1.0, 0.5)}, head_digest="e")
    report = gate(changed, registry(), no_measure)
    assert all(verdict.changed_output for verdict in report.verdicts)
    assert verdict_of(report, GT).verdict is Verdict.PASS
    assert verdict_of(report, OMEGA).verdict is Verdict.IMPROVED
    assert not any(verdict.changed_output for verdict in gate(passing(), registry(), no_measure).verdicts)


def test_the_regression_error_names_workload_ratio_margin_and_pairs() -> None:
    base = session({GT: pairs(1.0, 1.10), OMEGA: pairs(1.0, 1.0)})
    report = gate(base, registry(), FakeMeasure([(1.0, 1.10)] * 3))
    error = RatioRegressionError(report)
    assert error.report is report
    message = str(error)
    assert GT in message
    assert OMEGA not in message
    assert "1.1" in message
    assert "0.04" in message
    assert str(MAX_PAIRS) in message


def test_decide_lives_below_the_registry_and_the_verdict() -> None:
    # The registry's chain replay and the verdict share one decide; the shared module imports neither.
    tree = ast.parse(Path(perf_ratio_gate_decide.__file__).read_text(encoding="utf-8"))
    imported = {node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)}
    imported |= {alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names}
    assert not {name for name in imported if name and name.startswith(("tools.perf_ratio_gate_registry", "tools.perf_ratio_gate_verdict"))}
    assert perf_ratio_gate_verdict.decide is perf_ratio_gate_decide.decide
    assert perf_ratio_gate_verdict.Verdict is perf_ratio_gate_decide.Verdict
    # The verdict re-exports, never redefines, what moved down (an int's identity proves nothing, so read the source).
    verdict_tree = ast.parse(Path(perf_ratio_gate_verdict.__file__).read_text(encoding="utf-8"))
    defined = {
        target.id
        for node in verdict_tree.body
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        for target in (node.targets if isinstance(node, ast.Assign) else [node.target])
        if isinstance(target, ast.Name)
    }
    defined |= {node.name for node in verdict_tree.body if isinstance(node, (ast.FunctionDef, ast.ClassDef))}
    assert not defined & {"MAX_PAIRS", "Margin", "Anchor", "UNIT_ANCHOR", "Verdict", "decide", "anchored", "exact_ratio"}


def test_exact_ratio_is_the_record_ratio_double_read_exactly() -> None:
    # 1 / 3 is not a double: the record's ratio is its rounding, and that rounded double is the decided rho_w.
    run = session({GT: [(3.0, 1.0)]}).workloads[0]
    assert exact_ratio(run) == Fraction(run.ratio) == Fraction(1.0 / 3.0)
    assert exact_ratio(run) != Fraction(1, 3)
