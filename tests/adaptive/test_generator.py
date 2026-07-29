"""Contracts for exact global adaptive-generation orchestration."""

from unittest.mock import Mock

import pytest

from compas_cgal import _continuous_tea_2
import compas_cgal.adaptive.generator as generator_module
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.generator import GenerationContinuation
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.generator import generate_exact_adaptive_continuation
from compas_cgal.adaptive.generator import materialize_active_candidate_family
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.traversal import MatTraversalState
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from tests.adaptive.test_acceptance import _BranchFixture
from tests.adaptive.test_acceptance import _build_branch_fixture
from tests.adaptive.task13f_fixture import Task13FFixture


@pytest.fixture(scope="module")
def branch() -> _BranchFixture:
    """Build the exact terminal-launch branch-switch fixture once."""
    return _build_branch_fixture()


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build the exact radius-1 one-root continuation fixture once."""
    return Task13FFixture.build()


def _forward_limit_identities(
    branch: _BranchFixture,
) -> tuple[bytes, ...]:
    active = branch.traversal.active_cursor
    samples = branch.traversal.authority.sample_index.samples_by_edge[active.route_step.edge_id]
    return tuple(sample.cursor_identity for sample in samples[1 : 1 + branch.identity.traversal_policy.forward_window])


def test_active_family_materializes_each_forward_span_once(
    branch: _BranchFixture,
) -> None:
    """Build one complete invariant-ordered family over the finite window."""
    family = materialize_active_candidate_family(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
    )

    limit_identities = _forward_limit_identities(branch)
    assert family
    assert family == branch.identity.candidate_policy.order_candidates(
        family,
        key=lambda candidate: candidate.order_key,
    )
    assert {candidate.cursor_limit_identity for candidate in family} == set(limit_identities)
    assert branch.candidate in family


def test_active_family_commits_one_winner_and_builds_one_span_family_each(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Commit a real exact branch proof while counting structural builds."""
    accepted = branch.evaluator.evaluate_from_cursor(
        branch.physical,
        TraversalCursorState.before(
            branch.candidate.traversal_decision,
        ),
        branch.candidate,
    )
    original_enumeration = generator_module.enumerate_middle_curve_candidates
    counted_enumeration = Mock(wraps=original_enumeration)

    def select_known_exact_winner(
        *,
        evaluator: CandidateEvaluator,
        physical: GenerationState,
        traversal: MatTraversalState,
        candidates: tuple[MiddleCurveCandidate, ...],
    ) -> CandidateTransaction:
        assert branch.candidate in candidates
        return accepted

    monkeypatch.setattr(
        generator_module,
        "enumerate_middle_curve_candidates",
        counted_enumeration,
    )
    monkeypatch.setattr(
        generator_module,
        "evaluate_first_feasible_candidate",
        select_known_exact_winner,
    )

    physical_after, traversal_after, commit = advance_active_candidate_family(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
    )
    continuation = GenerationContinuation.build(
        launch_transaction=branch.launch_transaction,
        physical=physical_after,
        traversal=traversal_after,
        commits=(commit,),
    )

    assert counted_enumeration.call_count == len(
        _forward_limit_identities(branch),
    )
    assert type(commit) is TraversalCommit
    assert commit.transaction == accepted
    assert physical_after.digest == commit.physical_child_digest
    assert not traversal_after.active_cursor.terminal
    assert continuation.physical == physical_after
    assert continuation.traversal.active_route_index == 1
    assert continuation.commits == (commit,)


def test_real_active_family_stops_at_unresolved_exact_event(
    branch: _BranchFixture,
) -> None:
    """Refuse a later feasible candidate after a higher-ranked proof gap."""
    with pytest.raises(
        UnresolvedMotionEventError,
        match="unresolved exact event",
    ):
        advance_active_candidate_family(
            evaluator=branch.evaluator,
            physical=branch.physical,
            traversal=branch.traversal,
        )


def test_task13f_fourth_trial_accepts_mixed_seam_with_inactive_incidence(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reach the real post-link one-root circle through exact trial order.

    The fixed family has sixteen cells. Its first three trials are proved
    gouges. The fourth link certifies and leaves four simultaneous boundary
    vertices at the following circle's phase seam: two enter, one leaves, and
    one is tangent. The neighboring active sets are incomparable but do not
    exhaust the exact incidences. The circle is feasible under the pi cap, so
    its accepted transaction must bind the independently reproduced mixed
    native trace instead of rejecting its inactive tangent evidence.
    """
    family = materialize_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=task13f.physical,
        traversal=task13f.traversal,
    )

    assert len(family) == 16
    fourth = family[3]
    original_evaluate = CandidateEvaluator.evaluate_from_cursor
    calls: list[bytes] = []
    gouges: list[bytes] = []

    def track_real_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        calls.append(candidate.identity)
        try:
            return original_evaluate(
                self,
                state,
                traversal_before,
                candidate,
            )
        except GougeContainmentError:
            gouges.append(candidate.identity)
            raise

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_from_cursor",
        track_real_trial,
    )
    physical_after, traversal_after, commit = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=task13f.physical,
        traversal=task13f.traversal,
    )
    transaction = commit.transaction

    assert calls == [candidate.identity for candidate in family[:4]]
    assert gouges == [candidate.identity for candidate in family[:3]]
    assert transaction.candidate == fourth
    assert physical_after.digest == (commit.physical_child_digest)
    assert traversal_after.active_cursor.terminal

    motion = transaction.candidate.motion
    phase_point = Point2[WorldXY].build(
        motion.center.x + motion.phase_vector.x,
        motion.center.y + motion.phase_vector.y,
    )
    post_link_stock = task13f.physical.fork_stock()
    post_link_stock.deplete(
        ExactSegmentMotion.build(
            task13f.physical.phase_point,
            phase_point,
        ),
        task13f.identity.tool_radius,
        task13f.identity.depletion_policy,
    )
    verdict, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        post_link_stock.raw,
        motion.center.x,
        motion.center.y,
        motion.phase_vector.x,
        motion.phase_vector.y,
        motion.clockwise,
        task13f.identity.tool_radius.value,
        task13f.identity.user_cap.chord_ratio,
    )
    seam_fibres = tuple(fibre for fibre in trace.partition.fibres if fibre.seam_id)
    assert len(seam_fibres) == 1
    seam = seam_fibres[0]

    def physical_key(item: object) -> tuple[object, ...]:
        return (
            item.kind,
            item.feature_id,
            item.support_id,
            item.trim_id,
            item.vertex_id,
            item.endpoint_role,
        )

    left_incidence = {physical_key(branch.physical_incidence) for branch in seam.left_active_branches}
    right_incidence = {physical_key(branch.physical_incidence) for branch in seam.right_active_branches}
    active_incidence = left_incidence | right_incidence
    witnessed_incidence = {physical_key(witness) for witness in seam.local_event_witnesses if witness.kind == "endpoint-order"}

    assert verdict == "certified"
    assert trace.canonical_digest == (transaction.circle_witness.motion_witness.event_trace_digest)
    assert left_incidence - right_incidence
    assert right_incidence - left_incidence
    assert active_incidence < witnessed_incidence
    assert (seam.ccw_direction, seam.cw_direction) == (
        "mixed",
        "mixed",
    )
    assert (
        _continuous_tea_2.verify_event_partition(
            trace.partition,
        ).verdict.name
        == "CERTIFIED"
    )


def test_generation_commits_launch_before_propagating_route_uncertainty(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Authenticate and activate the launch root before any route trial."""
    observed: list[tuple[bytes, bytes]] = []

    def unresolved_route(
        *,
        evaluator: CandidateEvaluator,
        physical: GenerationState,
        traversal: MatTraversalState,
    ) -> tuple[GenerationState, MatTraversalState, TraversalCommit]:
        observed.append(
            (
                physical.digest,
                traversal.digest,
            )
        )
        raise UnresolvedMotionEventError(
            "route event partition remains unresolved",
        )

    monkeypatch.setattr(
        generator_module,
        "advance_active_candidate_family",
        unresolved_route,
    )

    with pytest.raises(
        UnresolvedMotionEventError,
        match="route event partition",
    ):
        generate_exact_adaptive_continuation(
            initial_evaluator=branch.initial_evaluator,
            evaluator=branch.evaluator,
            seeded_traversal=branch.seeded_traversal,
            launch_transaction=branch.launch_transaction,
        )

    assert observed == [
        (
            branch.physical.digest,
            branch.traversal.digest,
        )
    ]
