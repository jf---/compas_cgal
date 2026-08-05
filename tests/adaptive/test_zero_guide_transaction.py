"""Atomic transaction contracts for exact zero-guide MAT advances."""

import math
from dataclasses import dataclass
from dataclasses import fields
from dataclasses import replace
from typing import NoReturn
from typing import Self

import pytest

from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.candidates import enumerate_zero_guide_link_candidates
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.coverage import SweepWitness
from compas_cgal.adaptive.errors import InvalidZeroGuideTransactionError
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.errors import StaleZeroGuideTransactionError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import MotionWitness
from compas_cgal.adaptive.motion_certificate import SweptPrefixMotionWitness
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckOwnerId
from compas_cgal.adaptive.operation import NeckTraversalOrientation
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import WidthClassId
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import ZeroGuideLinkTransaction
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from tests.adaptive.task13f_fixture import Task13FFixture


@dataclass(frozen=True)
class _ZeroGuideTransactionFixture:
    """Real post-route-0 state and both owned zero-guide candidate variants."""

    evaluator: CandidateEvaluator
    state: GenerationState
    traversal_before: TraversalCursorState
    candidate: ZeroGuideLinkCandidate
    alternate_candidate: ZeroGuideLinkCandidate
    foreign_run_candidate: ZeroGuideLinkCandidate

    @classmethod
    def build(cls) -> Self:
        task13f = Task13FFixture.build()
        physical, route_zero, _ = advance_active_candidate_family(
            evaluator=task13f.evaluator,
            physical=task13f.physical,
            traversal=task13f.traversal,
        )
        route_one = route_zero.activate_next()
        if route_one.pending_transit is not None:
            raise AssertionError("the proved width-2 arm must be outside a causal neck transit")
        active = route_one.active_cursor
        span = MiddleCurveSpan.build(
            axis=route_one.authority.axis,
            cursor_before=active.cursor,
            cursor_limit=active.terminal_cursor,
        )
        full_cap = FullCapDecision.build(
            user_cap=task13f.evaluator.user_cap,
            effective_cap=task13f.evaluator.user_cap,
        )
        candidates = enumerate_zero_guide_link_candidates(
            span=span,
            policy=task13f.evaluator.candidate_policy,
            neck_scope=route_one.neck_scope,
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=True,
        )
        foreign_run = next(run for run in route_one.authority.axis.zero_guide_inventory.runs if run.edge_id != candidates[0].zero_guide_run.edge_id)
        foreign_samples = tuple(sample for sample in route_one.authority.axis.samples if sample.edge_id == foreign_run.edge_id)
        foreign_span = MiddleCurveSpan.build(
            axis=route_one.authority.axis,
            cursor_before=foreign_samples[0],
            cursor_limit=foreign_samples[-1],
        )
        foreign = enumerate_zero_guide_link_candidates(
            span=foreign_span,
            policy=task13f.evaluator.candidate_policy,
            neck_scope=route_one.neck_scope,
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=True,
        )[0]
        return cls(
            evaluator=task13f.evaluator,
            state=physical,
            traversal_before=TraversalCursorState.before(
                candidates[0].traversal_decision,
            ),
            candidate=candidates[0],
            alternate_candidate=candidates[1],
            foreign_run_candidate=foreign,
        )


@pytest.fixture(scope="module")
def zero_guide() -> _ZeroGuideTransactionFixture:
    """Build the costly native radius-1 continuation authority once per worker."""
    return _ZeroGuideTransactionFixture.build()


def _state_snapshot(state: GenerationState) -> tuple[object, ...]:
    stock = state.fork_stock()
    coverage = state.clone_coverage()
    return (
        state.digest,
        state.stock_boundary_digest,
        state.stock_lineage_digest,
        tuple(witness.digest for witness in stock.lineage),
        coverage.certificate.canonical_bytes,
        tuple(witness.digest for witness in coverage.lineage),
        state.phase_point,
        state.traversal,
        state.operations,
        state.passages,
    )


def test_zero_guide_evaluation_and_commit_are_ordered_and_byte_identical(
    zero_guide: _ZeroGuideTransactionFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bind the real arm segment through containment, TEA, and both mutations.

    Task 5 closes the gap between a proved spatial-only candidate and an
    accepted physical cut. The trial must certify against pre-depletion stock,
    update stock before coverage, and reproduce the identical transaction when
    committed from the same authenticated cursor.
    """
    observed: list[str] = []
    real_containment = GougeContainment.certify_segment
    real_certify = MotionCertifier.certify_swept_prefix_segment
    real_deplete = Stock2Area.deplete
    real_add_sweep = CoverageLedger.add_sweep

    def record_containment(
        owner: GougeContainment,
        motion: ExactSegmentMotion,
        tool_radius: ToolRadius,
    ) -> SegmentContainmentCertificate:
        observed.append("containment")
        return real_containment(owner, motion, tool_radius)

    def record_certify(
        owner: MotionCertifier,
        *,
        operation_index: int,
        motion: ExactSegmentMotion,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> SweptPrefixMotionWitness:
        observed.append("tea")
        return real_certify(
            owner,
            operation_index=operation_index,
            motion=motion,
            user_cap=user_cap,
            effective_cap=effective_cap,
        )

    def record_deplete(
        owner: Stock2Area,
        depletion: ExactSegmentMotion,
        tool_radius: ToolRadius,
        policy: DepletionPolicy,
    ) -> DepletionWitness:
        observed.append("depletion")
        return real_deplete(owner, depletion, tool_radius, policy)

    def record_add_sweep(
        owner: CoverageLedger,
        motion: ExactSegmentMotion,
        tool_radius: ToolRadius,
    ) -> SweepWitness:
        observed.append("coverage")
        return real_add_sweep(owner, motion, tool_radius)

    monkeypatch.setattr(GougeContainment, "certify_segment", record_containment)
    monkeypatch.setattr(
        MotionCertifier,
        "certify_swept_prefix_segment",
        record_certify,
    )
    monkeypatch.setattr(Stock2Area, "deplete", record_deplete)
    monkeypatch.setattr(CoverageLedger, "add_sweep", record_add_sweep)
    parent = _state_snapshot(zero_guide.state)

    transaction = zero_guide.evaluator.evaluate_zero_guide_from_cursor(
        zero_guide.state,
        zero_guide.traversal_before,
        zero_guide.candidate,
    )

    assert observed == ["containment", "tea", "depletion", "coverage"]
    assert type(transaction) is ZeroGuideLinkTransaction
    assert type(transaction.segment_witness.operation) is AdvanceSegmentOperation
    assert transaction.segment_witness.operation.motion.start == zero_guide.state.phase_point
    assert transaction.segment_witness.operation.motion.end == zero_guide.candidate.target
    assert transaction.segment_witness.motion_witness.stock_lineage_digest == zero_guide.state.stock_lineage_digest
    assert b"zero-guide-link-transaction-v1" in transaction.canonical_bytes
    assert _state_snapshot(zero_guide.state) == parent

    observed.clear()
    child = zero_guide.evaluator.commit_zero_guide_from_cursor(
        zero_guide.state,
        zero_guide.traversal_before,
        transaction,
    )

    assert observed == ["containment", "tea", "depletion", "coverage"]
    assert child.digest == transaction.result_state_digest
    assert child.operations[-1] == transaction.segment_witness.operation
    assert child.phase_point == zero_guide.candidate.target
    assert child.traversal == transaction.traversal_after
    assert _state_snapshot(zero_guide.state) == parent


class _InjectedStageFailure(RuntimeError):
    """Sentinel proving that one failed trial cannot leak forked mutations."""


def test_every_zero_guide_proof_stage_fails_atomically(
    zero_guide: _ZeroGuideTransactionFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Preserve every parent owner when any ordered proof stage aborts.

    Containment and TEA fail before mutation; depletion and coverage fail after
    progressively more work on trial-local forks. All four paths must leave
    stock, coverage, phase, cursor, operation history, and passage state intact.
    """
    parent = _state_snapshot(zero_guide.state)
    stages = (
        (GougeContainment, "certify_segment", "containment"),
        (MotionCertifier, "certify_swept_prefix_segment", "tea"),
        (Stock2Area, "deplete", "depletion"),
        (CoverageLedger, "add_sweep", "coverage"),
    )

    for owner, method_name, stage in stages:

        def fail_stage(*_arguments: object, **_keywords: object) -> NoReturn:
            raise _InjectedStageFailure(stage)

        with monkeypatch.context() as stage_patch:
            stage_patch.setattr(owner, method_name, fail_stage)
            with pytest.raises(_InjectedStageFailure, match=stage):
                zero_guide.evaluator.evaluate_zero_guide_from_cursor(
                    zero_guide.state,
                    zero_guide.traversal_before,
                    zero_guide.candidate,
                )
        assert _state_snapshot(zero_guide.state) == parent


def test_zero_guide_transaction_rejects_cross_wires_and_stale_commit_roots(
    zero_guide: _ZeroGuideTransactionFixture,
) -> None:
    """Reject every independently mutable proof, candidate, and commit root.

    Local semantic fields fail during transaction construction. Opaque parent
    and child state digests are authenticated at commit, where their preimages
    are available; a stale cursor is likewise distinguished from malformed
    transaction evidence.
    """
    transaction = zero_guide.evaluator.evaluate_zero_guide_from_cursor(
        zero_guide.state,
        zero_guide.traversal_before,
        zero_guide.candidate,
    )
    witness = transaction.segment_witness
    operation = witness.operation
    if type(operation) is not AdvanceSegmentOperation:
        raise AssertionError("zero-guide witness must own the advancing operation")

    target_candidate = ZeroGuideLinkCandidate.build(
        zero_guide_run=zero_guide.candidate.zero_guide_run,
        policy=zero_guide.candidate.policy,
        spatial_progress=zero_guide.candidate.spatial_progress,
        spatial_levels=zero_guide.candidate.spatial_levels,
        target=Point2[WorldXY].build(
            zero_guide.candidate.target.x,
            zero_guide.candidate.target.y + 0.125,
        ),
        cursor_limit_identity=zero_guide.candidate.cursor_limit_identity,
        neck_scope=zero_guide.candidate.neck_scope,
        effective_cap_decision=zero_guide.candidate.effective_cap_decision,
        traversal_decision=zero_guide.candidate.traversal_decision,
    )
    foreign_cap = EngagementCap.build(math.pi / 2.0)
    cap_candidate = ZeroGuideLinkCandidate.build(
        zero_guide_run=zero_guide.candidate.zero_guide_run,
        policy=zero_guide.candidate.policy,
        spatial_progress=zero_guide.candidate.spatial_progress,
        spatial_levels=zero_guide.candidate.spatial_levels,
        target=zero_guide.candidate.target,
        cursor_limit_identity=zero_guide.candidate.cursor_limit_identity,
        neck_scope=zero_guide.candidate.neck_scope,
        effective_cap_decision=FullCapDecision.build(
            user_cap=foreign_cap,
            effective_cap=foreign_cap,
        ),
        traversal_decision=zero_guide.candidate.traversal_decision,
    )
    oriented_candidate = ZeroGuideLinkCandidate.build(
        zero_guide_run=zero_guide.candidate.zero_guide_run,
        policy=zero_guide.candidate.policy,
        spatial_progress=zero_guide.candidate.spatial_progress,
        spatial_levels=zero_guide.candidate.spatial_levels,
        target=zero_guide.candidate.target,
        cursor_limit_identity=zero_guide.candidate.cursor_limit_identity,
        neck_scope=OrientedNeckScope.build(
            neck_owner_id=NeckOwnerId(b"foreign-neck-owner"),
            orientation=NeckTraversalOrientation.FORWARD,
        ),
        effective_cap_decision=NeckCapDecision.build(
            neck_evidence_digest=b"\x91" * 32,
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.FIRST_PASS_COMPLETE,
            user_cap=zero_guide.evaluator.user_cap,
            effective_cap=zero_guide.evaluator.user_cap,
        ),
        traversal_decision=zero_guide.candidate.traversal_decision,
    )
    cross_wired_operation = AdvanceSegmentOperation.build(
        motion=operation.motion,
        cut_z=operation.cut_z,
        neck_scope=operation.neck_scope,
        effective_cap_decision=operation.effective_cap_decision,
        traversal_decision=zero_guide.alternate_candidate.traversal_decision,
    )
    stock_cross_wire = replace(
        witness,
        motion_witness=replace(
            witness.motion_witness,
            stock_lineage_digest=b"\xa1" * 32,
        ),
    )
    coverage_cross_wire = replace(
        witness,
        sweep_witness=replace(
            witness.sweep_witness,
            parent_lineage=witness.sweep_witness.parent_lineage + (b"\xa2" * 32,),
        ),
    )
    operation_cross_wire = replace(
        witness,
        operation=cross_wired_operation,
    )
    relabelled_generic_motion = MotionWitness(
        witness.motion_witness.operation_index,
        witness.motion_witness.operation_kind,
        witness.motion_witness.motion,
        witness.motion_witness.user_cap_bytes,
        witness.motion_witness.effective_cap_bytes,
        b"foreign-segment-strategy-v1",
        witness.motion_witness.stock_lineage_digest,
        witness.motion_witness.event_trace_digest,
        "certified",
        witness.motion_witness.event_cell_count,
        0,
    )
    with pytest.raises(InvalidReplayTraceError, match="swept-prefix"):
        replace(
            witness,
            motion_witness=relabelled_generic_motion,
        )
    bypassed_replay = object.__new__(ReplayLateralWitness)
    for field in fields(witness):
        object.__setattr__(
            bypassed_replay,
            field.name,
            getattr(witness, field.name),
        )
    object.__setattr__(
        bypassed_replay,
        "motion_witness",
        relabelled_generic_motion,
    )
    with pytest.raises(InvalidZeroGuideTransactionError, match="swept-prefix"):
        replace(
            transaction,
            segment_witness=bypassed_replay,
        )
    semantic_mutations = (
        {"candidate": zero_guide.foreign_run_candidate},
        {"candidate": target_candidate},
        {"candidate": cap_candidate},
        {"candidate": oriented_candidate},
        {"segment_witness": stock_cross_wire},
        {"segment_witness": coverage_cross_wire},
        {"segment_witness": operation_cross_wire},
        {"traversal_after": zero_guide.traversal_before},
        {"passage_after": object()},
    )
    for mutation in semantic_mutations:
        with pytest.raises(InvalidZeroGuideTransactionError):
            replace(transaction, **mutation)  # type: ignore[arg-type]

    coverage_parents = witness.sweep_witness.parent_lineage
    if not coverage_parents:
        raise AssertionError("the real post-route-0 fixture must own coverage history")
    opaque_coverage_cross_wire = replace(
        witness,
        sweep_witness=replace(
            witness.sweep_witness,
            parent_lineage=(b"\xa3" * 32, *coverage_parents[1:]),
        ),
    )
    structurally_valid_cross_wire = replace(
        transaction,
        segment_witness=opaque_coverage_cross_wire,
    )
    with pytest.raises(InvalidZeroGuideTransactionError, match="replay"):
        zero_guide.evaluator.commit_zero_guide_from_cursor(
            zero_guide.state,
            zero_guide.traversal_before,
            structurally_valid_cross_wire,
        )

    stale_parent = replace(
        transaction,
        parent_state_digest=b"\xb1" * 32,
    )
    with pytest.raises(StaleZeroGuideTransactionError, match="parent"):
        zero_guide.evaluator.commit_zero_guide_from_cursor(
            zero_guide.state,
            zero_guide.traversal_before,
            stale_parent,
        )

    stale_cursor = replace(
        zero_guide.traversal_before,
        cursor=CursorIdentity(b"stale-zero-guide-cursor"),
    )
    with pytest.raises(StaleZeroGuideTransactionError, match="cursor"):
        zero_guide.evaluator.commit_zero_guide_from_cursor(
            zero_guide.state,
            stale_cursor,
            transaction,
        )

    foreign_child = replace(
        transaction,
        result_state_digest=b"\xb2" * 32,
    )
    with pytest.raises(InvalidZeroGuideTransactionError, match="replay"):
        zero_guide.evaluator.commit_zero_guide_from_cursor(
            zero_guide.state,
            zero_guide.traversal_before,
            foreign_child,
        )
