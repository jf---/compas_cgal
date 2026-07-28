import math
from dataclasses import dataclass
from dataclasses import replace
from fractions import Fraction
from itertools import permutations

import pytest
from compas.geometry import Polygon

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import CandidateSelectionError
from compas_cgal.adaptive.errors import CandidateStateMismatchError
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidCandidateTransactionError
from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.errors import StaleCandidateTransactionError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import MotionWitness
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import EffectiveCapDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MatSamplingPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.transaction import select_candidate_transaction
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType

OUTER = (
    (0.0, 0.0),
    (6.0, 0.0),
    (6.0, 2.0),
    (2.0, 2.0),
    (2.0, 6.0),
    (0.0, 6.0),
)


@dataclass(frozen=True)
class _StateFixture:
    identity: InputIdentity
    first: MiddleCurveCandidate
    second: MiddleCurveCandidate
    terminal_candidates: tuple[MiddleCurveCandidate, ...]
    state: GenerationState
    source_stock: Stock2Area
    source_coverage: CoverageLedger


def _pocket() -> Polygon:
    return Polygon([[x, y, 0.0] for x, y in OUTER])


def _ring() -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(tuple(Point2[WorldXY].build(x, y) for x, y in OUTER))


def _candidate_policy() -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.5),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=2,
        phase_count=4,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(0.25),
    )


def _axis() -> MedialAxis:
    return MedialAxis.build(
        design_boundary=_ring(),
        holes=(),
        tool_radius=ToolRadius.build(0.5),
        station_spacing=Spacing.build(0.75),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )


def _first_candidate_and_span(
    *,
    user_cap: EngagementCap,
    candidate_policy: CandidatePolicy,
) -> tuple[MiddleCurveCandidate, MiddleCurveSpan]:
    axis = _axis()
    edge = next(
        edge
        for edge in axis.edges
        if edge.curve_kind == "parabola" and tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)[-1].point == Point2[WorldXY].build(1.0, 2.0)
    )
    samples = tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[1],
        cursor_limit=samples[2],
    )
    full_cap = FullCapDecision.build(
        user_cap=user_cap,
        effective_cap=user_cap,
    )
    candidate = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    return candidate, span


def _first_candidate(
    *,
    user_cap: EngagementCap,
    candidate_policy: CandidatePolicy,
) -> MiddleCurveCandidate:
    candidate, _span = _first_candidate_and_span(
        user_cap=user_cap,
        candidate_policy=candidate_policy,
    )
    return candidate


def _input_identity(
    *,
    neck_cap_degrees: tuple[float, float, float] = (90.0, 80.0, 70.0),
) -> InputIdentity:
    design_boundary = _ring()
    tool_radius = ToolRadius.build(0.5)
    candidate_policy = _candidate_policy()
    user_cap = EngagementCap.build(math.radians(120.0))
    first = _first_candidate(
        user_cap=user_cap,
        candidate_policy=candidate_policy,
    )
    phase = _phase_point(first)
    reachable_domain = ReachableDomain.build(
        design_boundary=design_boundary,
        holes=(),
        tool_radius=tool_radius,
    )
    cut_plane = CutPlane.build(
        CutZ.build(-2.0),
        ClearanceZ.build(5.0),
    )
    entry = PreclearedEntry.build(
        reachable_domain=reachable_domain,
        center=phase,
        radius=EntryRadius.build(0.75),
        tool_radius=tool_radius,
        cut_plane=cut_plane,
        qualified_bore=QualifiedBore.build(
            cut_plane=cut_plane,
            process_identity=BoreProcessIdentity(b"predrill-cycle-revision-7"),
            evidence_bytes=b"qualified-bore-metrology-v1",
        ),
    )
    neck_policy = NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
        effective_caps={
            (neck_class, passage_state): EngagementCap.build(
                math.radians(neck_cap_degrees[passage_state.rank]),
            )
            for neck_class in range(2)
            for passage_state in ACTIVE_PASSAGE_STATES
        },
    )
    return InputIdentity.build(
        design_boundary=design_boundary,
        holes=(),
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        reachable_domain=reachable_domain,
        entry=entry,
        user_cap=user_cap,
        mat_sampling_policy=MatSamplingPolicy.build(
            station_spacing=Spacing.build(0.75),
            max_sagitta=ChordBound.build(0.02),
            max_refinement_depth=32,
        ),
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=DepletionPolicy.build(
            chord_bound=ChordBound.build(0.125),
            center_count_limit=4096,
        ),
        traversal_policy=TraversalPolicy.build(forward_window=4),
        cut_direction_policy=CutDirectionPolicy.build(CutIntent.CLIMB),
    )


def _candidate_pair(
    identity: InputIdentity,
) -> tuple[
    MiddleCurveCandidate,
    MiddleCurveCandidate,
    tuple[MiddleCurveCandidate, ...],
]:
    first, first_span = _first_candidate_and_span(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    full_cap = FullCapDecision.build(
        user_cap=identity.user_cap,
        effective_cap=identity.user_cap,
    )
    terminal_candidates = _terminal_candidates(
        identity=identity,
        first=first,
        first_span=first_span,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=full_cap,
    )
    second = next(
        candidate
        for candidate in terminal_candidates
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 3
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    return first, second, terminal_candidates


def _terminal_candidates(
    *,
    identity: InputIdentity,
    first: MiddleCurveCandidate,
    neck_scope: NoNeckScope | OrientedNeckScope,
    effective_cap_decision: EffectiveCapDecision,
    first_span: MiddleCurveSpan,
) -> tuple[MiddleCurveCandidate, ...]:
    axis = first_span.axis
    edge = next(edge for edge in axis.edges if bytes(edge.identity) == bytes(first.traversal_decision.edge_id))
    samples = tuple(
        sorted(
            (sample for sample in axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    terminal_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=DerivedCandidateCursor.build(
            span=first_span,
            candidate=first,
        ),
        cursor_limit=samples[2],
    )
    return enumerate_middle_curve_candidates(
        span=terminal_span,
        policy=identity.candidate_policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=neck_scope,
        effective_cap_decision=effective_cap_decision,
        makes_cursor_terminal_at_limit=True,
    )


def _oriented_candidate_pair(
    identity: InputIdentity,
) -> tuple[
    MiddleCurveCandidate,
    MiddleCurveCandidate,
    NeckPassage,
    tuple[MiddleCurveCandidate, ...],
]:
    axis = _axis()
    edge = next(
        edge
        for edge in axis.edges
        if edge.curve_kind == "parabola" and tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)[-1].point == Point2[WorldXY].build(1.0, 2.0)
    )
    samples = tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)
    first_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[1],
        cursor_limit=samples[2],
    )
    first_passage = (
        NeckInventory.build(
            axis=axis,
            policy=identity.neck_policy,
        )
        .necks[0]
        .forward
    )
    first_decision = first_passage.propose_cap_decision(identity.neck_policy)
    first = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=first_span,
            policy=identity.candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=first_passage.scope,
            effective_cap_decision=first_decision,
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    second_passage = first_passage.advance(first_decision)
    second_decision = second_passage.propose_cap_decision(identity.neck_policy)
    terminal_candidates = _terminal_candidates(
        identity=identity,
        first=first,
        first_span=first_span,
        neck_scope=second_passage.scope,
        effective_cap_decision=second_decision,
    )
    second = next(
        candidate
        for candidate in terminal_candidates
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 3
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    return first, second, second_passage, terminal_candidates


def _phase_point(candidate: MiddleCurveCandidate) -> Point2[WorldXY]:
    return Point2[WorldXY].build(
        candidate.motion.center.x + candidate.motion.phase_vector.x,
        candidate.motion.center.y + candidate.motion.phase_vector.y,
    )


def _circle_operation(
    identity: InputIdentity,
    candidate: MiddleCurveCandidate,
) -> CutFullCircleOperation:
    return CutFullCircleOperation.build(
        motion=candidate.motion,
        cut_z=identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
        neck_scope=candidate.neck_scope,
        effective_cap_decision=candidate.effective_cap_decision,
        traversal_decision=candidate.traversal_decision,
    )


def _passages(identity: InputIdentity) -> tuple[NeckPassage, ...]:
    inventory = NeckInventory.build(
        axis=_axis(),
        policy=identity.neck_policy,
    )
    return tuple(passage for neck in inventory.necks for passage in (neck.forward, neck.reverse))


def _state_fixture(
    *,
    oriented: bool = False,
    neck_cap_degrees: tuple[float, float, float] = (
        120.0,
        120.0,
        120.0,
    ),
) -> _StateFixture:
    identity = _input_identity(neck_cap_degrees=neck_cap_degrees)
    passages = _passages(identity)
    if oriented:
        first, second, second_passage, terminal_candidates = _oriented_candidate_pair(identity)
        passages = tuple(second_passage if passage.scope == second_passage.scope else passage for passage in passages)
    else:
        first, second, terminal_candidates = _candidate_pair(identity)
    stock = Stock2Area.build(Stock(_pocket(), []))
    stock.deplete(identity.entry)
    stock.deplete(
        first.motion,
        identity.tool_radius,
        identity.depletion_policy,
    )
    coverage = CoverageLedger.build(
        reachable_domain=identity.reachable_domain,
        precleared_center=identity.entry.center,
        precleared_radius=identity.entry.radius,
    )
    coverage.add_sweep(first.motion, identity.tool_radius)
    traversal = TraversalCursorState.before(first.traversal_decision).advance(first.traversal_decision)
    state = GenerationState.build(
        stock=stock,
        coverage=coverage,
        tool_radius=identity.tool_radius,
        phase_point=_phase_point(first),
        traversal=traversal,
        passages=passages,
        operations=(
            identity.entry.approach,
            identity.entry.plunge,
            _circle_operation(identity, first),
        ),
    )
    return _StateFixture(
        identity,
        first,
        second,
        terminal_candidates,
        state,
        stock,
        coverage,
    )


def _evaluator(identity: InputIdentity) -> CandidateEvaluator:
    return CandidateEvaluator.build(
        reachable_domain=identity.reachable_domain,
        tool_radius=identity.tool_radius,
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
        neck_policy=identity.neck_policy,
        depletion_policy=identity.depletion_policy,
        cut_direction_policy=identity.cut_direction_policy,
        cut_z=identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
    )


def _candidate_with_policy(
    candidate: MiddleCurveCandidate,
    policy: CandidatePolicy,
) -> MiddleCurveCandidate:
    return MiddleCurveCandidate.build(
        proposal=candidate.proposal,
        policy=policy,
        spatial_progress=candidate.spatial_progress,
        spatial_levels=candidate.spatial_levels,
        radius_levels=candidate.radius_levels,
        cursor_limit_identity=candidate.cursor_limit_identity,
        neck_scope=candidate.neck_scope,
        effective_cap_decision=candidate.effective_cap_decision,
        traversal_decision=candidate.traversal_decision,
    )


def _alternate_accepted_candidate(
    fixture: _StateFixture,
) -> MiddleCurveCandidate:
    return next(
        candidate
        for candidate in fixture.terminal_candidates
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )


def _real_gouging_candidate(
    fixture: _StateFixture,
) -> MiddleCurveCandidate:
    return next(
        candidate
        for candidate in fixture.terminal_candidates
        if candidate.guide_radius == Fraction(1, 4)
        and candidate.phase_index == 3
        and candidate.proposal.generator_site.kind == "point"
        and candidate.spatial_progress > Fraction(1, 4)
    )


def test_generation_state_clones_native_owners_and_binds_every_lineage() -> None:
    """Keep caller aliases and every accepted lineage outside authority.

    Candidate evaluation relies on one stable parent digest. Mutating the
    `Stock2Area` used to construct that parent must not mutate the authoritative
    snapshot or silently change which material the transaction certifies.
    """
    fixture = _state_fixture()
    digest = fixture.state.digest
    authoritative_lineage = fixture.state.fork_stock().lineage

    fixture.source_stock.deplete(
        fixture.second.motion,
        fixture.identity.tool_radius,
        fixture.identity.depletion_policy,
    )

    assert fixture.state.digest == digest
    assert fixture.state.fork_stock().lineage == authoritative_lineage
    assert fixture.state.fork_stock().lineage != fixture.source_stock.lineage


def test_generation_state_rejects_foreign_tool_radius() -> None:
    """Reject a state whose declared cutter differs from its proof lineage.

    Stock and coverage witnesses already bind the physical tool radius. A
    second radius on the state would let later containment and TEA queries
    certify a different cutter against the same parent digest.
    """
    fixture = _state_fixture()

    with pytest.raises(InvalidGenerationStateError, match="tool radius"):
        GenerationState.build(
            stock=fixture.source_stock,
            coverage=fixture.source_coverage,
            tool_radius=ToolRadius.build(0.25),
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=fixture.state.passages,
            operations=fixture.state.operations,
        )


def test_generation_state_rejects_stock_coverage_chronology_divergence() -> None:
    """Require stock depletion and coverage sweeps to name the same motions.

    A state with an extra stock depletion but no matching coverage sweep would
    under-report the region claimed by the eventual coverage certificate.
    """
    fixture = _state_fixture()
    divergent_stock = fixture.source_stock.fork()
    divergent_stock.deplete(
        fixture.second.motion,
        fixture.identity.tool_radius,
        fixture.identity.depletion_policy,
    )

    with pytest.raises(InvalidGenerationStateError, match="chronology"):
        GenerationState.build(
            stock=divergent_stock,
            coverage=fixture.source_coverage,
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=fixture.state.passages,
            operations=fixture.state.operations,
        )


def test_generation_state_rejects_duplicate_oriented_passage_scope() -> None:
    """Give each oriented separator passage exactly one authoritative state.

    Two values for the same owner and orientation would make effective-cap
    derivation dependent on lookup order rather than exact identity.
    """
    fixture = _state_fixture()
    duplicate = fixture.state.passages[0]

    with pytest.raises(InvalidGenerationStateError, match="passage"):
        GenerationState.build(
            stock=fixture.source_stock,
            coverage=fixture.source_coverage,
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=fixture.state.passages + (duplicate,),
            operations=fixture.state.operations,
        )


def test_generation_state_rejects_phase_or_cursor_outside_operation_prefix() -> None:
    """Bind the next link start and MAT cursor to the accepted circle prefix.

    Geometry-only phase agreement cannot repair a stale cursor, and a correct
    cursor cannot authorize a link starting anywhere except the last circle's
    exact phase point.
    """
    fixture = _state_fixture()

    with pytest.raises(InvalidGenerationStateError, match="phase"):
        GenerationState.build(
            stock=fixture.source_stock,
            coverage=fixture.source_coverage,
            tool_radius=fixture.identity.tool_radius,
            phase_point=Point2[WorldXY].build(1.0, 1.0),
            traversal=fixture.state.traversal,
            passages=fixture.state.passages,
            operations=fixture.state.operations,
        )

    with pytest.raises(InvalidGenerationStateError, match="traversal"):
        GenerationState.build(
            stock=fixture.source_stock,
            coverage=fixture.source_coverage,
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=TraversalCursorState.before(fixture.first.traversal_decision),
            passages=fixture.state.passages,
            operations=fixture.state.operations,
        )


def test_generation_state_rejects_unearned_passage_state() -> None:
    """Require untouched oriented passages to remain exactly unvisited.

    A future neck cannot begin at its second cap merely because its immutable
    endpoint is structurally valid. The accepted operation prefix must earn
    every passage transition in order.
    """
    fixture = _state_fixture()
    untouched = fixture.state.passages[0]
    forged = replace(
        untouched,
        state=PassageState.FIRST_PASS_COMPLETE,
    )

    with pytest.raises(InvalidGenerationStateError, match="passage"):
        GenerationState.build(
            stock=fixture.state.fork_stock(),
            coverage=fixture.state.clone_coverage(),
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=tuple(forged if passage.scope == forged.scope else passage for passage in fixture.state.passages),
            operations=fixture.state.operations,
        )


def test_generation_state_rejects_discontinuous_passage_history() -> None:
    """Replay neck transitions rather than trusting the latest endpoint.

    Replacing the first accepted cap decision with the second locally valid
    transition creates a history beginning at `FIRST_PASS_COMPLETE`. Matching
    the final passage to that forged operation must not authenticate the gap.
    """
    fixture = _state_fixture(
        oriented=True,
        neck_cap_degrees=(120.0, 120.0, 120.0),
    )
    scope = fixture.second.neck_scope
    second_decision = fixture.second.effective_cap_decision
    first_circle = fixture.state.operations[-1]
    assert type(scope) is OrientedNeckScope
    assert type(second_decision) is NeckCapDecision
    assert type(first_circle) is CutFullCircleOperation
    forged_circle = replace(
        first_circle,
        effective_cap_decision=second_decision,
    )
    forged_passage = fixture.state.passage(scope).advance(second_decision)

    with pytest.raises(InvalidGenerationStateError, match="passage"):
        GenerationState.build(
            stock=fixture.state.fork_stock(),
            coverage=fixture.state.clone_coverage(),
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=tuple(forged_passage if passage.scope == scope else passage for passage in fixture.state.passages),
            operations=fixture.state.operations[:-1] + (forged_circle,),
        )


def test_generation_state_rejects_lateral_cut_depth_drift() -> None:
    """Keep every accepted lateral motion on the qualified entry cut plane.

    Stock and coverage are planar, so their motion lineages cannot reveal a
    foreign Z value. The operation chronology must bind that typed value back
    to the plunge that established the machining plane.
    """
    fixture = _state_fixture()
    first_circle = fixture.state.operations[-1]
    assert type(first_circle) is CutFullCircleOperation
    foreign_circle = replace(
        first_circle,
        cut_z=CutZ.build(-1.0),
    )

    with pytest.raises(InvalidGenerationStateError, match="cut depth"):
        GenerationState.build(
            stock=fixture.state.fork_stock(),
            coverage=fixture.state.clone_coverage(),
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=fixture.state.passages,
            operations=fixture.state.operations[:-1] + (foreign_circle,),
        )


def test_evaluator_rejects_cut_depth_outside_parent_entry() -> None:
    """Reject a typed evaluator depth not authorized by the parent state.

    Without this check, a valid 2D candidate and all of its exact witnesses
    could be emitted at an unrelated Z coordinate because native stock and TEA
    authority operate in the cut plane.
    """
    fixture = _state_fixture()
    evaluator = CandidateEvaluator.build(
        reachable_domain=fixture.identity.reachable_domain,
        tool_radius=fixture.identity.tool_radius,
        user_cap=fixture.identity.user_cap,
        candidate_policy=fixture.identity.candidate_policy,
        neck_policy=fixture.identity.neck_policy,
        depletion_policy=fixture.identity.depletion_policy,
        cut_direction_policy=fixture.identity.cut_direction_policy,
        cut_z=CutZ.build(-1.0),
        material_side=MaterialSide.OUTSIDE,
    )
    before = fixture.state.digest

    with pytest.raises(CandidateStateMismatchError, match="cut depth"):
        evaluator.evaluate(fixture.state, fixture.second)

    assert fixture.state.digest == before


def test_generation_state_rejects_cross_wired_link_circle_history() -> None:
    """Keep accepted link metadata paired to its following circle.

    Stock and coverage bind only the exact motions. Replacing the link's hold
    cursor with the circle's result cursor preserves local operation validity
    but breaks the causal link/circle transaction recorded by the state.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    transaction = evaluator.evaluate(fixture.state, fixture.second)
    committed = evaluator.commit(fixture.state, transaction)
    link = committed.operations[-2]
    circle = committed.operations[-1]
    assert type(link) is LinkSegmentOperation
    assert type(circle) is CutFullCircleOperation
    traversal = circle.traversal_decision
    cross_wired_link = replace(
        link,
        traversal_decision=HoldTraversalDecision.build(
            component_id=traversal.component_id,
            edge_id=traversal.edge_id,
            branch_id=traversal.branch_id,
            cursor=traversal.cursor_after,
        ),
    )

    with pytest.raises(InvalidGenerationStateError, match="link.*circle"):
        GenerationState.build(
            stock=committed.fork_stock(),
            coverage=committed.clone_coverage(),
            tool_radius=fixture.identity.tool_radius,
            phase_point=committed.phase_point,
            traversal=committed.traversal,
            passages=committed.passages,
            operations=committed.operations[:-2] + (cross_wired_link, circle),
        )


def test_evaluator_rejects_foreign_candidate_policy() -> None:
    """Bind every candidate transaction to the declared input lattice.

    The same geometric proposal can be re-identified under another policy.
    Exact containment and TEA would still pass, so policy equality must fail
    before either native query.
    """
    fixture = _state_fixture()
    foreign_policy = CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.25),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=2,
        phase_count=4,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(0.25),
    )
    foreign = _candidate_with_policy(
        fixture.second,
        foreign_policy,
    )

    with pytest.raises(CandidateStateMismatchError, match="candidate policy"):
        _evaluator(fixture.identity).evaluate(
            fixture.state,
            foreign,
        )


def test_evaluator_rejects_cut_direction_mismatch() -> None:
    """Require material side, cut intent, and circle orientation to agree.

    The L-pocket candidates are counterclockwise for climb cutting with
    material outside. Conventional intent requires the opposite orientation,
    even though the circle geometry and engagement proof are unchanged.
    """
    fixture = _state_fixture()
    evaluator = CandidateEvaluator.build(
        reachable_domain=fixture.identity.reachable_domain,
        tool_radius=fixture.identity.tool_radius,
        user_cap=fixture.identity.user_cap,
        candidate_policy=fixture.identity.candidate_policy,
        neck_policy=fixture.identity.neck_policy,
        depletion_policy=fixture.identity.depletion_policy,
        cut_direction_policy=CutDirectionPolicy.build(CutIntent.CONVENTIONAL),
        cut_z=fixture.identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
    )

    with pytest.raises(CandidateStateMismatchError, match="cut direction"):
        evaluator.evaluate(
            fixture.state,
            fixture.second,
        )


def test_evaluator_rejects_foreign_parent_depletion_policy() -> None:
    """Keep historical stock construction under the declared policy.

    A state digest can faithfully bind depletions produced with another chord
    bound. The evaluator must reject that parent before extending it under a
    different policy and creating a mixed-strategy lineage.
    """
    fixture = _state_fixture()
    foreign_policy = DepletionPolicy.build(
        chord_bound=ChordBound.build(0.25),
        center_count_limit=4096,
    )
    stock = Stock2Area.build(Stock(_pocket(), []))
    stock.deplete(fixture.identity.entry)
    stock.deplete(
        fixture.first.motion,
        fixture.identity.tool_radius,
        foreign_policy,
    )
    foreign_state = GenerationState.build(
        stock=stock,
        coverage=fixture.state.clone_coverage(),
        tool_radius=fixture.identity.tool_radius,
        phase_point=fixture.state.phase_point,
        traversal=fixture.state.traversal,
        passages=fixture.state.passages,
        operations=fixture.state.operations,
    )

    with pytest.raises(CandidateStateMismatchError, match="depletion policy"):
        _evaluator(fixture.identity).evaluate(
            foreign_state,
            fixture.second,
        )


def test_evaluator_rejects_foreign_parent_cap_policy() -> None:
    """Revalidate historical cap decisions before extending a parent.

    Planar stock and coverage cannot reveal that a prior no-neck circle was
    recorded under another user cap. The evaluator must prevent that policy
    drift from entering a new transaction lineage.
    """
    fixture = _state_fixture()
    first_circle = fixture.state.operations[-1]
    assert type(first_circle) is CutFullCircleOperation
    foreign_cap = EngagementCap.build(math.radians(110.0))
    foreign_circle = replace(
        first_circle,
        effective_cap_decision=FullCapDecision.build(
            user_cap=foreign_cap,
            effective_cap=foreign_cap,
        ),
    )
    foreign_state = GenerationState.build(
        stock=fixture.state.fork_stock(),
        coverage=fixture.state.clone_coverage(),
        tool_radius=fixture.identity.tool_radius,
        phase_point=fixture.state.phase_point,
        traversal=fixture.state.traversal,
        passages=fixture.state.passages,
        operations=fixture.state.operations[:-1] + (foreign_circle,),
    )

    with pytest.raises(CandidateStateMismatchError, match="cap decision"):
        _evaluator(fixture.identity).evaluate(
            foreign_state,
            fixture.second,
        )


def test_successful_evaluation_is_deterministic_and_non_mutating() -> None:
    """Return byte-identical proof without changing authoritative owners.

    The candidate is a real accepted L-pocket link/circle pair. Repeating its
    evaluation must observe the same parent stock, produce the same exact
    witnesses, and leave both state lineage counts unchanged.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    parent_digest = fixture.state.digest
    stock_lineage = fixture.state.fork_stock().lineage
    coverage_lineage = fixture.state.clone_coverage().lineage

    first = evaluator.evaluate(fixture.state, fixture.second)
    second = evaluator.evaluate(fixture.state, fixture.second)

    assert type(first) is CandidateTransaction
    assert first.canonical_bytes == second.canonical_bytes
    assert first.digest == second.digest
    assert first.parent_state_digest == parent_digest
    assert fixture.state.digest == parent_digest
    assert fixture.state.fork_stock().lineage == stock_lineage
    assert fixture.state.clone_coverage().lineage == coverage_lineage


def test_circle_rejection_after_certified_link_is_atomic() -> None:
    """Discard a passing link when its joint circle exceeds the neck cap.

    The reproduced second passage uses the production `80 degrees` cap. Its
    direct link certifies, but its circle is exactly over cap at operation 4;
    the link's trial depletion and coverage sweep must not escape.
    """
    fixture = _state_fixture(
        oriented=True,
        neck_cap_degrees=(90.0, 80.0, 70.0),
    )
    evaluator = _evaluator(fixture.identity)
    before = fixture.state.digest

    with pytest.raises(
        EngagementCapExceededError,
        match="operation 4",
    ):
        evaluator.evaluate(fixture.state, fixture.second)

    assert fixture.state.digest == before
    assert len(fixture.state.fork_stock().lineage) == 2
    assert len(fixture.state.clone_coverage().lineage) == 1


def test_link_rejection_stops_before_circle_and_preserves_parent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Propagate the named link proof failure before any circle dispatch.

    Real link cap exceedance is gated by `MotionCertifier` tests. This consumer
    boundary injects that same named failure at operation 3 and proves the
    candidate evaluator neither dispatches the circle nor mutates its parent.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    calls: list[OperationType] = []
    original = MotionCertifier.certify

    def reject_link(
        self: MotionCertifier,
        *,
        operation_index: int,
        operation_kind: OperationType,
        motion: ExactSegmentMotion | ExactCircleMotion,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> MotionWitness:
        calls.append(operation_kind)
        if operation_kind is OperationType.LINK:
            raise EngagementCapExceededError(f"operation {operation_index} exceeds its exact effective cap.")
        return original(
            self,
            operation_index=operation_index,
            operation_kind=operation_kind,
            motion=motion,
            user_cap=user_cap,
            effective_cap=effective_cap,
        )

    monkeypatch.setattr(MotionCertifier, "certify", reject_link)
    before = fixture.state.digest

    with pytest.raises(EngagementCapExceededError, match="operation 3"):
        evaluator.evaluate(fixture.state, fixture.second)

    assert calls == [OperationType.LINK]
    assert fixture.state.digest == before


def test_unresolved_circle_discards_the_completed_link_trial(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep unresolved exact event authority distinct from rejection.

    The link follows the real native certified path. Injecting the already
    gated named unresolved verdict at the circle boundary proves that neither
    the link depletion nor its coverage sweep becomes authoritative.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    calls: list[OperationType] = []
    original = MotionCertifier.certify

    def reject_circle(
        self: MotionCertifier,
        *,
        operation_index: int,
        operation_kind: OperationType,
        motion: ExactSegmentMotion | ExactCircleMotion,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> MotionWitness:
        calls.append(operation_kind)
        if operation_kind is OperationType.CUT:
            raise UnresolvedMotionEventError(f"operation {operation_index} has an unresolved exact event partition.")
        return original(
            self,
            operation_index=operation_index,
            operation_kind=operation_kind,
            motion=motion,
            user_cap=user_cap,
            effective_cap=effective_cap,
        )

    monkeypatch.setattr(MotionCertifier, "certify", reject_circle)
    before = fixture.state.digest

    with pytest.raises(UnresolvedMotionEventError, match="operation 4"):
        evaluator.evaluate(fixture.state, fixture.second)

    assert calls == [OperationType.LINK, OperationType.CUT]
    assert fixture.state.digest == before


def test_real_circle_gouge_discards_its_certified_link_trial() -> None:
    """Reject a real finite-lattice circle outside the exact design domain.

    The larger point-site candidate's direct link is legal, but its full
    cutter sweep leaves the concave L pocket. Containment must fail before its
    circle TEA query and preserve the authoritative parent.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    gouging = _real_gouging_candidate(fixture)
    before = fixture.state.digest

    with pytest.raises(GougeContainmentError, match="full-circle"):
        evaluator.evaluate(fixture.state, gouging)

    assert fixture.state.digest == before


def test_commit_replays_winner_and_advances_every_owner_together() -> None:
    """Commit only independently reproduced link-and-circle evidence.

    Functional commit returns a new state with two stock depletions, two
    coverage sweeps, two operations, the candidate cursor, and circle phase;
    the parent remains a stable snapshot.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    transaction = evaluator.evaluate(fixture.state, fixture.second)
    parent_digest = fixture.state.digest

    committed = evaluator.commit(fixture.state, transaction)

    assert fixture.state.digest == parent_digest
    assert committed.digest == transaction.result_state_digest
    assert len(committed.fork_stock().lineage) == 4
    assert len(committed.clone_coverage().lineage) == 3
    assert len(committed.operations) == len(fixture.state.operations) + 2
    assert committed.operations[-2] == transaction.link_witness.operation
    assert committed.operations[-1] == transaction.circle_witness.operation
    assert committed.phase_point == _phase_point(fixture.second)
    assert committed.traversal == transaction.traversal_after


def test_commit_rejects_stale_parent_before_native_replay() -> None:
    """Prevent an accepted candidate from committing after state advances.

    The transaction names its complete parent digest. Reusing it against the
    committed child must fail before evaluating a now-terminal cursor.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    transaction = evaluator.evaluate(fixture.state, fixture.second)
    committed = evaluator.commit(fixture.state, transaction)

    with pytest.raises(StaleCandidateTransactionError, match="parent"):
        evaluator.commit(committed, transaction)


def test_transaction_rejects_deplete_before_certify_lineage() -> None:
    """Detect a witness captured from post-depletion instead of pre-state.

    Replacing the link motion witness lineage with the circle's post-link
    lineage simulates depletion before certification. The immutable transaction
    constructor must reject that reordered proof even though each witness is
    individually well formed.
    """
    fixture = _state_fixture()
    transaction = _evaluator(fixture.identity).evaluate(
        fixture.state,
        fixture.second,
    )
    post_link_motion = replace(
        transaction.link_witness.motion_witness,
        stock_lineage_digest=transaction.circle_witness.motion_witness.stock_lineage_digest,
    )
    reordered_link = replace(
        transaction.link_witness,
        motion_witness=post_link_motion,
    )

    with pytest.raises(
        InvalidCandidateTransactionError,
        match="certify-before-deplete",
    ):
        replace(transaction, link_witness=reordered_link)


def test_transaction_rejects_swapped_link_circle_bundles() -> None:
    """Keep the direct link and accepted circle in canonical causal order.

    Both proof bundles are valid in isolation. Swapping them must still fail
    because operation kinds, indices, and candidate ownership are relational
    transaction facts.
    """
    fixture = _state_fixture()
    transaction = _evaluator(fixture.identity).evaluate(
        fixture.state,
        fixture.second,
    )

    with pytest.raises(
        InvalidCandidateTransactionError,
        match="link.*circle",
    ):
        replace(
            transaction,
            link_witness=transaction.circle_witness,
            circle_witness=transaction.link_witness,
        )


def test_lateral_witness_rejects_cross_wired_containment() -> None:
    """Keep each containment certificate with its exact owned motion.

    Segment and circle certificates are independently valid, but replacing a
    link certificate with the following circle certificate must fail at the
    smallest proof-bundle boundary.
    """
    fixture = _state_fixture()
    transaction = _evaluator(fixture.identity).evaluate(
        fixture.state,
        fixture.second,
    )

    with pytest.raises(InvalidReplayTraceError, match="containment"):
        replace(
            transaction.link_witness,
            containment_certificate=transaction.circle_witness.containment_certificate,
        )


def test_transaction_rejects_reindexed_proof_bundle() -> None:
    """Bind internal witness validity to consecutive transaction position.

    Reindexing both the link bundle and its motion witness preserves their
    local agreement. It must still fail because the link can no longer precede
    its circle at the accepted operation boundary.
    """
    fixture = _state_fixture()
    transaction = _evaluator(fixture.identity).evaluate(
        fixture.state,
        fixture.second,
    )
    reindexed_motion = replace(
        transaction.link_witness.motion_witness,
        operation_index=transaction.circle_witness.operation_index,
    )
    reindexed_link = replace(
        transaction.link_witness,
        operation_index=transaction.circle_witness.operation_index,
        motion_witness=reindexed_motion,
    )

    with pytest.raises(
        InvalidCandidateTransactionError,
        match="indices",
    ):
        replace(
            transaction,
            link_witness=reindexed_link,
        )


def test_oriented_passage_advances_only_in_committed_state() -> None:
    """Keep neck policy state immutable throughout candidate evaluation.

    The accepted equal-cap control proposes the second oriented transition.
    Evaluation records its `passage_after`, but only the committed child may
    expose that state; pre-advancing the parent contradicts its operation
    chronology.
    """
    fixture = _state_fixture(
        oriented=True,
        neck_cap_degrees=(120.0, 120.0, 120.0),
    )
    evaluator = _evaluator(fixture.identity)
    transaction = evaluator.evaluate(fixture.state, fixture.second)
    scope = fixture.second.neck_scope
    decision = fixture.second.effective_cap_decision
    assert type(scope) is OrientedNeckScope
    assert type(decision) is NeckCapDecision
    parent_passage = fixture.state.passage(scope)
    assert transaction.passage_after is not None
    assert parent_passage.state is decision.passage_before

    committed = evaluator.commit(fixture.state, transaction)

    assert fixture.state.passage(scope) == parent_passage
    assert committed.passage(scope) == transaction.passage_after
    with pytest.raises(InvalidGenerationStateError, match="passage"):
        GenerationState.build(
            stock=fixture.state.fork_stock(),
            coverage=fixture.state.clone_coverage(),
            tool_radius=fixture.identity.tool_radius,
            phase_point=fixture.state.phase_point,
            traversal=fixture.state.traversal,
            passages=tuple(transaction.passage_after if passage.scope == scope else passage for passage in fixture.state.passages),
            operations=fixture.state.operations,
        )


def test_candidate_must_begin_at_authoritative_cursor() -> None:
    """Reject a previously accepted lattice cell at the next cursor state.

    Candidate identity carries cursor-before. Geometry and cap may still be
    valid, but replaying the first circle after it already advanced would break
    traversal lineage.
    """
    fixture = _state_fixture()

    with pytest.raises(CandidateStateMismatchError, match="traversal"):
        _evaluator(fixture.identity).evaluate(
            fixture.state,
            fixture.first,
        )


def test_winner_selection_is_permutation_invariant() -> None:
    """Select by exact policy order rather than evaluation completion order.

    Two real terminal candidates have equal progress and radius but different
    phase and canonical identity. Every input permutation must select the same
    lower canonical identity.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    alternate = _alternate_accepted_candidate(fixture)
    transactions = (
        evaluator.evaluate(fixture.state, fixture.second),
        evaluator.evaluate(fixture.state, alternate),
    )

    winners = {select_candidate_transaction(tuple(order)).digest for order in permutations(transactions)}

    assert winners == {transactions[1].digest}


def test_winner_selection_rejects_malformed_acceptance_sets() -> None:
    """Require one unique-policy acceptance set rooted at one state digest.

    Empty, duplicated, cross-parent, and cross-policy sets cannot represent
    one finite-lattice competition observed against one stock snapshot.
    """
    fixture = _state_fixture()
    evaluator = _evaluator(fixture.identity)
    first = evaluator.evaluate(fixture.state, fixture.second)
    second = evaluator.evaluate(
        fixture.state,
        _alternate_accepted_candidate(fixture),
    )
    foreign = replace(
        second,
        parent_state_digest=b"\xff" * 32,
    )
    foreign_policy = CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.25),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=2,
        phase_count=4,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(0.25),
    )
    foreign_candidate = _candidate_with_policy(
        second.candidate,
        foreign_policy,
    )
    cross_policy = replace(
        second,
        candidate=foreign_candidate,
    )

    with pytest.raises(CandidateSelectionError, match="nonempty"):
        select_candidate_transaction(())
    with pytest.raises(CandidateSelectionError, match="parent"):
        select_candidate_transaction((first, foreign))
    with pytest.raises(CandidateSelectionError, match="unique"):
        select_candidate_transaction((first, first))
    with pytest.raises(CandidateSelectionError, match="policy"):
        select_candidate_transaction((first, cross_policy))
