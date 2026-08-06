"""Fresh-replay contracts for exact zero-guide segment advances."""

from copy import copy

import pytest
from compas.geometry import Polygon

import compas_cgal.adaptive.replay as replay_module
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.errors import ReplayGrammarError
from compas_cgal.adaptive.errors import ReplayPairingError
from compas_cgal.adaptive.errors import ReplayTraversalError
from compas_cgal.adaptive.errors import ReplayZeroGuideCandidateError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.medial_axis import MatZeroGuideInventory
from compas_cgal.adaptive.medial_axis import MatZeroGuideRun
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import EdgeId
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.replay import replay_certificate
from compas_cgal.adaptive.replay_trace import FreshReplayTrace
from compas_cgal.adaptive.transaction import ZeroGuideLinkTransaction
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from tests.adaptive.task13f_fixture import TASK13F_OUTER
from tests.adaptive.task13f_fixture import Task13FFixture


def _pocket() -> Polygon:
    return Polygon(
        [[x, y, 0.0] for x, y in TASK13F_OUTER],
    )


def _replay(
    identity: InputIdentity,
    operations: tuple[CanonicalOperation, ...],
) -> None:
    replay_certificate(
        input_identity=identity,
        pocket=_pocket(),
        holes=(),
        cut_plane=identity.cut_plane,
        tool_radius=identity.tool_radius,
        user_cap=identity.user_cap,
        entry=identity.entry,
        operations=operations,
        candidate_policy=identity.candidate_policy,
        neck_policy=identity.neck_policy,
        depletion_policy=identity.depletion_policy,
        traversal_policy=identity.traversal_policy,
        cut_direction_policy=identity.cut_direction_policy,
    )


def _task13f_route_one_child() -> tuple[
    Task13FFixture,
    GenerationState,
    ZeroGuideLinkTransaction,
]:
    """Build the committed circle route and its first zero-guide advance.

    Returns:
        Authenticated fixture, route-1 physical child, and its exact accepted
        advancing-segment transaction.

    Raises:
        AssertionError: If global dispatch selects a foreign transaction
            variant for the proved zero-guide edge.
    """
    fixture = Task13FFixture.build()
    physical_zero, traversal_zero, _ = advance_active_candidate_family(
        evaluator=fixture.evaluator,
        physical=fixture.physical,
        traversal=fixture.traversal,
    )
    physical_one, _, commit = advance_active_candidate_family(
        evaluator=fixture.evaluator,
        physical=physical_zero,
        traversal=traversal_zero.activate_next(),
    )
    transaction = commit.transaction
    assert type(transaction) is ZeroGuideLinkTransaction
    return fixture, physical_one, transaction


def test_replay_reconstructs_task13f_zero_guide_before_terminal_gate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Rebuild the committed advancing segment from fresh MAT proof bytes.

    Task 13F remains globally nonterminal after route 1. Before that expected
    rejection, replay must select the unique spatial candidate owned by the
    rebuilt MAT and reproduce the live segment witness byte for byte.
    """
    fixture, physical_one, transaction = _task13f_route_one_child()
    selected: list[ZeroGuideLinkCandidate] = []
    traces: list[FreshReplayTrace] = []
    real_match = replay_module._match_zero_guide_candidate
    real_terminal = replay_module._require_terminal_replay

    def tracked_match(**kwargs: object) -> ZeroGuideLinkCandidate:
        candidate = real_match(**kwargs)  # type: ignore[arg-type]
        selected.append(candidate)
        return candidate

    def tracked_terminal(**kwargs: object) -> None:
        trace = kwargs["trace"]
        assert type(trace) is FreshReplayTrace
        traces.append(trace)
        real_terminal(**kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(
        replay_module,
        "_match_zero_guide_candidate",
        tracked_match,
    )
    monkeypatch.setattr(
        replay_module,
        "_require_terminal_replay",
        tracked_terminal,
    )

    with pytest.raises(
        ReplayTraversalError,
        match="fresh MAT traversal remains nonterminal",
    ):
        _replay(fixture.identity, physical_one.operations)

    assert selected == [transaction.candidate]
    assert len(traces) == 1
    assert traces[0].lateral_witnesses[-1].canonical_bytes == (transaction.segment_witness.canonical_bytes)


def test_replay_rejects_zero_guide_relabelling_and_proof_mutations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject every ambiguity between hold transport and proved advancement.

    These mutations preserve enough local geometry to look plausible while
    breaking operation grammar, fresh-edge ownership, endpoint identity,
    native proof presence, native proof bytes, or unique lattice lineage. None
    may fall through to ordinary circle reconstruction or physical replay.
    """
    fixture, physical_one, transaction = _task13f_route_one_child()
    operations = physical_one.operations
    route_link = operations[-3]
    route_circle = operations[-2]
    advance = operations[-1]
    assert type(route_link) is LinkSegmentOperation
    assert type(route_circle) is CutFullCircleOperation
    assert type(advance) is AdvanceSegmentOperation

    hold_relabelled_as_advance = AdvanceSegmentOperation.build(
        motion=route_link.motion,
        cut_z=route_link.cut_z,
        neck_scope=route_link.neck_scope,
        effective_cap_decision=route_link.effective_cap_decision,
        traversal_decision=route_circle.traversal_decision,
    )
    traversal = advance.traversal_decision
    advance_relabelled_as_hold = LinkSegmentOperation.build(
        motion=advance.motion,
        cut_z=advance.cut_z,
        neck_scope=advance.neck_scope,
        effective_cap_decision=advance.effective_cap_decision,
        traversal_decision=HoldTraversalDecision.build(
            component_id=traversal.component_id,
            edge_id=traversal.edge_id,
            branch_id=traversal.branch_id,
            cursor=traversal.cursor_before,
        ),
    )
    endpoint_step = float(
        fixture.identity.candidate_policy.spatial_resolution.value,
    )
    changed_endpoint = AdvanceSegmentOperation.build(
        motion=ExactSegmentMotion.build(
            advance.motion.start,
            Point2[WorldXY].build(
                advance.motion.end.x + endpoint_step,
                advance.motion.end.y,
            ),
        ),
        cut_z=advance.cut_z,
        neck_scope=advance.neck_scope,
        effective_cap_decision=advance.effective_cap_decision,
        traversal_decision=traversal,
    )
    axis = fixture.traversal.authority.axis
    foreign_edge = next(edge for edge in axis.edges if edge.identity not in axis.zero_guide_run_by_edge_id)
    changed_edge = AdvanceSegmentOperation.build(
        motion=advance.motion,
        cut_z=advance.cut_z,
        neck_scope=advance.neck_scope,
        effective_cap_decision=advance.effective_cap_decision,
        traversal_decision=AdvanceTraversalDecision.build(
            component_id=axis.component_by_edge_id[foreign_edge.identity],
            edge_id=EdgeId(bytes(foreign_edge.identity)),
            branch_id=foreign_edge.branch_id,
            cursor_before=traversal.cursor_before,
            cursor_after=traversal.cursor_after,
            makes_cursor_terminal=traversal.makes_cursor_terminal,
        ),
    )
    operation_mutations = (
        (*operations[:-3], hold_relabelled_as_advance, *operations[-2:]),
        (*operations[:-1], changed_endpoint),
        (*operations[:-1], changed_edge),
    )
    for mutated in operation_mutations:
        with pytest.raises(
            ReplayZeroGuideCandidateError,
            match="zero-guide",
        ):
            _replay(fixture.identity, mutated)
    with pytest.raises(
        ReplayPairingError,
        match="following circle",
    ):
        _replay(
            fixture.identity,
            (*operations[:-1], advance_relabelled_as_hold),
        )

    real_rebuild = replay_module._rebuild_medial_axis

    def rebuild_without_zero_guide(**kwargs: object) -> MedialAxis:
        rebuilt = real_rebuild(**kwargs)  # type: ignore[arg-type]
        mutated_axis = copy(rebuilt)
        mutated_proof = copy(rebuilt.proof)
        object.__setattr__(
            mutated_proof,
            "zero_guide_inventory",
            MatZeroGuideInventory.build(
                mat_certificate=rebuilt.mat_certificate,
                records=(),
            ),
        )
        object.__setattr__(mutated_axis, "proof", mutated_proof)
        return mutated_axis

    with monkeypatch.context() as scoped:
        scoped.setattr(
            replay_module,
            "_rebuild_medial_axis",
            rebuild_without_zero_guide,
        )
        with pytest.raises(
            ReplayZeroGuideCandidateError,
            match="zero-guide",
        ):
            _replay(fixture.identity, operations)
    real_enumerate = replay_module.enumerate_zero_guide_link_candidates

    def enumerate_with_foreign_run(
        **kwargs: object,
    ) -> tuple[ZeroGuideLinkCandidate, ...]:
        candidates = real_enumerate(**kwargs)  # type: ignore[arg-type]
        owned = next(
            (candidate for candidate in candidates if candidate.target == advance.motion.end and candidate.traversal_decision == traversal),
            None,
        )
        if owned is None:
            return candidates
        other_run = next(run for run in axis.zero_guide_inventory.runs if run.native_certificate != owned.zero_guide_run.native_certificate)
        foreign_run = MatZeroGuideRun.build(
            edge_id=owned.zero_guide_run.edge_id,
            mat_certificate=axis.mat_certificate,
            native_certificate=other_run.native_certificate,
        )
        foreign = ZeroGuideLinkCandidate.build(
            zero_guide_run=foreign_run,
            policy=owned.policy,
            spatial_progress=owned.spatial_progress,
            spatial_levels=owned.spatial_levels,
            target=owned.target,
            cursor_limit_identity=owned.cursor_limit_identity,
            neck_scope=owned.neck_scope,
            effective_cap_decision=owned.effective_cap_decision,
            traversal_decision=owned.traversal_decision,
        )
        return tuple(foreign if candidate == owned else candidate for candidate in candidates)

    with monkeypatch.context() as scoped:
        scoped.setattr(
            replay_module,
            "enumerate_zero_guide_link_candidates",
            enumerate_with_foreign_run,
        )
        with pytest.raises(
            ReplayZeroGuideCandidateError,
            match="zero-guide",
        ):
            _replay(fixture.identity, operations)

    def enumerate_with_ambiguous_match(
        **kwargs: object,
    ) -> tuple[ZeroGuideLinkCandidate, ...]:
        candidates = real_enumerate(**kwargs)  # type: ignore[arg-type]
        owned = next(
            (candidate for candidate in candidates if candidate == transaction.candidate),
            None,
        )
        if owned is None:
            return candidates
        duplicate_match = ZeroGuideLinkCandidate.build(
            zero_guide_run=owned.zero_guide_run,
            policy=owned.policy,
            spatial_progress=owned.spatial_progress,
            spatial_levels=(
                *owned.spatial_levels,
                owned.spatial_levels[-1] + 1,
            ),
            target=owned.target,
            cursor_limit_identity=owned.cursor_limit_identity,
            neck_scope=owned.neck_scope,
            effective_cap_decision=owned.effective_cap_decision,
            traversal_decision=owned.traversal_decision,
        )
        return (duplicate_match, *candidates)

    with monkeypatch.context() as scoped:
        scoped.setattr(
            replay_module,
            "enumerate_zero_guide_link_candidates",
            enumerate_with_ambiguous_match,
        )
        with pytest.raises(
            ReplayZeroGuideCandidateError,
            match="zero-guide",
        ):
            _replay(fixture.identity, operations)


def test_replay_rejects_complete_lateral_grammar_before_candidate_search(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject misplaced circles and advances before fresh MAT enumeration.

    Legal lateral types are insufficient: the first lateral operation must be
    the entry circle, and every later circle must be immediately preceded by
    its hold link. These failures must precede candidate reconstruction.
    """
    fixture, physical_one, _ = _task13f_route_one_child()
    operations = physical_one.operations
    entry_circle = operations[2]
    hold_link = operations[3]
    later_circle = operations[4]
    advance = operations[5]
    assert type(entry_circle) is CutFullCircleOperation
    assert type(hold_link) is LinkSegmentOperation
    assert type(later_circle) is CutFullCircleOperation
    assert type(advance) is AdvanceSegmentOperation

    def reject_candidate_search(**kwargs: object) -> None:
        raise AssertionError(
            f"malformed grammar reached candidate search ({kwargs!r})",
        )

    monkeypatch.setattr(
        replay_module,
        "_replay_candidate_stream",
        reject_candidate_search,
    )

    with pytest.raises(
        ReplayGrammarError,
        match="first lateral operation",
    ):
        _replay(
            fixture.identity,
            (
                *operations[:2],
                advance,
                entry_circle,
                hold_link,
                later_circle,
            ),
        )
    with pytest.raises(
        ReplayGrammarError,
        match="later circle",
    ):
        _replay(
            fixture.identity,
            (
                *operations[:3],
                later_circle,
                advance,
            ),
        )
