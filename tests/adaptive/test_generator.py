"""Contracts for exact global adaptive-generation orchestration."""

from unittest.mock import Mock

import pytest

import compas_cgal.adaptive.generator as generator_module
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.generator import GenerationContinuation
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.generator import generate_exact_adaptive_continuation
from compas_cgal.adaptive.generator import materialize_active_candidate_family
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.traversal import MatTraversalState
from tests.adaptive.test_acceptance import _BranchFixture
from tests.adaptive.test_acceptance import _build_branch_fixture


@pytest.fixture(scope="module")
def branch() -> _BranchFixture:
    """Build the exact terminal-launch branch-switch fixture once."""
    return _build_branch_fixture()


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
