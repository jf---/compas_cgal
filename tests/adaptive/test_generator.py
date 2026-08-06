"""Contracts for exact global adaptive-generation orchestration."""

import hashlib
from copy import copy
from unittest.mock import Mock

import pytest

from compas_cgal import _continuous_tea_2
import compas_cgal.adaptive.generator as generator_module
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import TraversalCandidate
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidCandidateFamilyError
from compas_cgal.adaptive.errors import InvalidTraversalCommitError
from compas_cgal.adaptive.errors import NoFeasibleCandidateError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.generator import GenerationContinuation
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.generator import generate_exact_adaptive_continuation
from compas_cgal.adaptive.generator import materialize_active_candidate_family
from compas_cgal.adaptive.medial_axis import MatZeroGuideRun
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.transaction import ZeroGuideLinkTransaction
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


def test_task13f_route_one_accepts_exact_zero_guide_link(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Machine the width-2 arm by its exact MAT centerline segment.

    Route 0 must retain its established fourth-circle winner. Activating route
    1 then selects only the native-proved zero-guide family, stops dispatch at
    its first accepted member, and commits one advancing segment whose endpoint
    becomes the authoritative physical phase point.
    """
    route_zero_family = materialize_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=task13f.physical,
        traversal=task13f.traversal,
    )
    physical_zero, traversal_zero, route_zero_commit = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=task13f.physical,
        traversal=task13f.traversal,
    )
    traversal_one = traversal_zero.activate_next()
    route_one_family = materialize_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=physical_zero,
        traversal=traversal_one,
    )
    active_run = traversal_one.authority.axis.zero_guide_run_by_edge_id.get(
        traversal_one.active_cursor.route_step.edge_id,
    )
    assert active_run is not None
    with pytest.raises(
        InvalidCandidateFamilyError,
        match="active MAT proof variant",
    ):
        generator_module.evaluate_first_feasible_candidate(
            evaluator=task13f.evaluator,
            physical=physical_zero,
            traversal=traversal_one,
            candidates=route_zero_family,
        )
    other_run = next(run for run in traversal_one.authority.axis.zero_guide_inventory.runs if run != active_run)
    foreign_run = MatZeroGuideRun.build(
        edge_id=active_run.edge_id,
        mat_certificate=traversal_one.authority.axis.mat_certificate,
        native_certificate=other_run.native_certificate,
    )
    owned_candidate = route_one_family[0]
    assert type(owned_candidate) is ZeroGuideLinkCandidate
    foreign_candidate = ZeroGuideLinkCandidate.build(
        zero_guide_run=foreign_run,
        policy=owned_candidate.policy,
        spatial_progress=owned_candidate.spatial_progress,
        spatial_levels=owned_candidate.spatial_levels,
        target=owned_candidate.target,
        cursor_limit_identity=owned_candidate.cursor_limit_identity,
        neck_scope=owned_candidate.neck_scope,
        effective_cap_decision=owned_candidate.effective_cap_decision,
        traversal_decision=owned_candidate.traversal_decision,
    )
    with pytest.raises(
        InvalidCandidateFamilyError,
        match="owned native proof bytes",
    ):
        generator_module.evaluate_first_feasible_candidate(
            evaluator=task13f.evaluator,
            physical=physical_zero,
            traversal=traversal_one,
            candidates=(foreign_candidate, *route_one_family[1:]),
        )
    with pytest.raises(
        InvalidTraversalCommitError,
        match="owned native proof bytes",
    ):
        generator_module.evaluate_traversal_candidate(
            evaluator=task13f.evaluator,
            physical=physical_zero,
            traversal=traversal_one,
            candidate=foreign_candidate,
        )
    attempts: list[TraversalCandidate] = []
    real_evaluate = CandidateEvaluator.evaluate_zero_guide_from_cursor

    def track_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: ZeroGuideLinkCandidate,
    ) -> ZeroGuideLinkTransaction:
        attempts.append(candidate)
        return real_evaluate(
            self,
            state,
            traversal_before,
            candidate,
        )

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_zero_guide_from_cursor",
        track_trial,
    )

    physical_after, traversal_after, route_one_commit = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=physical_zero,
        traversal=traversal_one,
    )

    assert type(route_zero_commit.transaction) is CandidateTransaction
    assert len(route_one_family) == 36
    assert all(
        type(candidate) is ZeroGuideLinkCandidate and candidate.zero_guide_run == active_run and candidate.zero_guide_run.native_certificate == active_run.native_certificate
        for candidate in route_one_family
    )
    assert attempts
    assert all(type(candidate) is ZeroGuideLinkCandidate for candidate in attempts)
    assert len(attempts) == 1
    assert type(route_one_commit.transaction) is ZeroGuideLinkTransaction
    assert route_one_commit.transaction.candidate == route_one_family[0] == attempts[0]
    assert len(physical_after.operations) == len(physical_zero.operations) + 1
    assert physical_after.operations[-1].motion.end == physical_after.phase_point
    assert traversal_after.active_cursor != traversal_one.active_cursor

    transaction = route_one_commit.transaction
    foreign_transaction = ZeroGuideLinkTransaction.build(
        parent_state_digest=transaction.parent_state_digest,
        candidate=foreign_candidate,
        segment_witness=transaction.segment_witness,
        traversal_after=transaction.traversal_after,
        passage_after=transaction.passage_after,
        result_state_digest=transaction.result_state_digest,
    )
    with pytest.raises(
        InvalidTraversalCommitError,
        match="owned native proof bytes",
    ):
        generator_module.commit_traversal_candidate(
            evaluator=task13f.evaluator,
            physical=physical_zero,
            traversal=traversal_one,
            transaction=foreign_transaction,
        )
    hidden_prefix_child = copy(physical_after)
    object.__setattr__(
        hidden_prefix_child,
        "operations",
        (
            *physical_zero.operations,
            physical_zero.operations[-1],
            transaction.segment_witness.operation,
        ),
    )
    with pytest.raises(
        InvalidTraversalCommitError,
        match="exact transaction suffix",
    ):
        TraversalCommit.build(
            physical_before=physical_zero,
            traversal_before=traversal_one,
            transaction=transaction,
            physical_after=hidden_prefix_child,
            traversal_after=traversal_after,
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


def test_task13f_full_continuation(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Expose the first exact boundary after both established Task 13F routes.

    The bounded generator must rediscover the mixed-seam circle and zero-guide
    segment from the authenticated launch root. Route 2 then exposes a physical
    discontinuity in DFS edge-discovery order: every direct link from the
    horizontal leaf gouges, although reversing the accepted route-1 segment
    restores six containment-valid link/circle pairs. This contract preserves
    the route, cursor, family, disposition, and geometric evidence needed by
    the next inter-route-transit stage.
    """
    observed: list[
        tuple[
            GenerationState,
            MatTraversalState,
            tuple[TraversalCandidate, ...],
        ]
    ] = []
    real_materialize = generator_module.materialize_active_candidate_family

    def track_family(
        *,
        evaluator: CandidateEvaluator,
        physical: GenerationState,
        traversal: MatTraversalState,
    ) -> tuple[TraversalCandidate, ...]:
        family = real_materialize(
            evaluator=evaluator,
            physical=physical,
            traversal=traversal,
        )
        observed.append((physical, traversal, family))
        return family

    monkeypatch.setattr(
        generator_module,
        "materialize_active_candidate_family",
        track_family,
    )

    with pytest.raises(
        NoFeasibleCandidateError,
        match=("finite candidate family exhausted.*attempts=56; cap=0; gouge=56; degenerate-link=0"),
    ) as raised:
        generate_exact_adaptive_continuation(
            initial_evaluator=task13f.initial_evaluator,
            evaluator=task13f.evaluator,
            seeded_traversal=task13f.seeded_traversal,
            launch_transaction=task13f.launch_transaction,
        )

    assert [traversal.active_route_index for _, traversal, _ in observed] == [
        0,
        1,
        2,
    ]
    assert [len(family) for _, _, family in observed] == [16, 36, 56]
    assert all(type(candidate) is MiddleCurveCandidate for candidate in observed[0][2])
    assert all(type(candidate) is ZeroGuideLinkCandidate for candidate in observed[1][2])
    assert all(type(candidate) is MiddleCurveCandidate for candidate in observed[2][2])

    route = task13f.seeded_traversal.authority.route
    route_two_physical, route_two_traversal, route_two_family = observed[2]
    route_two_cursor = route_two_traversal.active_cursor
    assert route_two_cursor.route_step == route[2]
    assert route[1].exit_node_id != route[2].entry_node_id
    assert route[2].entry_node_id == route[0].entry_node_id
    assert route_two_physical.phase_point == Point2[WorldXY].build(5.0, 1.0)
    assert len(route_two_physical.operations) == 6
    assert route_two_cursor.cursor.cursor_identity.hex() == ("def1bf1471e2df355ba488378ffad7b9e20116ac81909d9f9822a5a62b6abbb0")
    assert hashlib.sha256(
        b"".join(bytes(candidate.identity) for candidate in route_two_family),
    ).hexdigest() == ("56d92fcf2089c1e771a9add79bd356d172904a69e6c7acb18ab02522f4ed093b")
    assert str(raised.value) == (
        "finite candidate family exhausted at cursor=def1bf1471e2df355ba488378ffad7b9e20116ac81909d9f9822a5a62b6abbb0; attempts=56; cap=0; gouge=56; degenerate-link=0."
    )

    containment = GougeContainment.build(
        task13f.identity.reachable_domain,
    )
    route_one_parent = observed[1][0]
    reverse_motion = ExactSegmentMotion.build(
        route_two_physical.phase_point,
        route_one_parent.phase_point,
    )
    containment.certify_segment(
        reverse_motion,
        task13f.identity.tool_radius,
    )

    direct_link_gouges = 0
    restored_links = 0
    contained_circles = 0
    restored_pairs = 0
    for candidate in route_two_family:
        assert type(candidate) is MiddleCurveCandidate
        phase_point = Point2[WorldXY].build(
            candidate.motion.center.x + candidate.motion.phase_vector.x,
            candidate.motion.center.y + candidate.motion.phase_vector.y,
        )
        try:
            containment.certify_segment(
                ExactSegmentMotion.build(
                    route_two_physical.phase_point,
                    phase_point,
                ),
                task13f.identity.tool_radius,
            )
        except GougeContainmentError:
            direct_link_gouges += 1
        try:
            containment.certify_segment(
                ExactSegmentMotion.build(
                    route_one_parent.phase_point,
                    phase_point,
                ),
                task13f.identity.tool_radius,
            )
        except GougeContainmentError:
            link_contained = False
        else:
            link_contained = True
        try:
            containment.certify_full_circle(
                candidate.motion,
                task13f.identity.tool_radius,
            )
        except GougeContainmentError:
            circle_contained = False
        else:
            circle_contained = True
        restored_links += int(link_contained)
        contained_circles += int(circle_contained)
        restored_pairs += int(link_contained and circle_contained)

    assert (
        direct_link_gouges,
        restored_links,
        contained_circles,
        restored_pairs,
    ) == (56, 10, 30, 6)
