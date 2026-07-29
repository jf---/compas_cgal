"""Acceptance contracts for the exact entry-circle bootstrap."""

import math
from dataclasses import dataclass
from dataclasses import replace
from fractions import Fraction

import pytest
from compas.geometry import Polygon

from compas_cgal.adaptive.bootstrap import InitialCandidateEvaluator
from compas_cgal.adaptive.bootstrap import InitialCandidateTransaction
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import DegenerateSegmentMotionError
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import EngagementCapInfeasibleError
from compas_cgal.adaptive.errors import InitialCandidatePhaseError
from compas_cgal.adaptive.errors import InitialCandidateStateMismatchError
from compas_cgal.adaptive.errors import InvalidCandidateFamilyError
from compas_cgal.adaptive.errors import InvalidInitialCandidateTransactionError
from compas_cgal.adaptive.errors import InvalidTraversalCommitError
from compas_cgal.adaptive.errors import NeckTooTightError
from compas_cgal.adaptive.errors import NoFeasibleCandidateError
from compas_cgal.adaptive.errors import StaleInitialCandidateTransactionError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.generator import commit_traversal_candidate
from compas_cgal.adaptive.generator import evaluate_first_feasible_candidate
from compas_cgal.adaptive.generator import evaluate_traversal_candidate
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import MatSamplingPolicy
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.traversal import MatTraversalState
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
ENTRY_RADIUS = Fraction(9, 16)
LAUNCH_GUIDE_RADIUS = Fraction(1, 32)
LAUNCH_PROGRESS = Fraction(1, 2)
BRANCH_TOOL_RADIUS = Fraction(1, 4)
BRANCH_ENTRY_RADIUS = Fraction(35, 32)
BRANCH_LAUNCH_GUIDE_RADIUS = Fraction(27, 64)
BRANCH_GUIDE_RADIUS = Fraction(13, 32)
BRANCH_PROGRESS = Fraction(1, 4)


@dataclass(frozen=True)
class _LaunchFixture:
    identity: InputIdentity
    axis: MedialAxis
    traversal: MatTraversalState
    candidates: tuple[MiddleCurveCandidate, ...]
    candidate: MiddleCurveCandidate
    evaluator: InitialCandidateEvaluator


@dataclass(frozen=True)
class _BranchFixture:
    identity: InputIdentity
    seeded_traversal: MatTraversalState
    initial_evaluator: InitialCandidateEvaluator
    launch_transaction: InitialCandidateTransaction
    physical: GenerationState
    traversal: MatTraversalState
    candidate: MiddleCurveCandidate
    alternatives: tuple[MiddleCurveCandidate, ...]
    evaluator: CandidateEvaluator


def _ring() -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y in OUTER),
    )


def _pocket() -> Polygon:
    return Polygon([[x, y, 0.0] for x, y in OUTER])


def _sampling_policy(
    *,
    station_spacing: float = 0.75,
) -> MatSamplingPolicy:
    return MatSamplingPolicy.build(
        station_spacing=Spacing.build(station_spacing),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )


def _candidate_policy() -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.5),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(float(LAUNCH_GUIDE_RADIUS)),
        radius_refinement_levels=2,
        phase_count=4,
        minimum_guide_radius=GuideRadius.build(float(LAUNCH_GUIDE_RADIUS)),
        minimum_progress=Spacing.build(0.25),
    )


def _user_cap() -> EngagementCap:
    return EngagementCap.build(math.radians(120.0))


def _neck_policy(user_cap: EngagementCap) -> NeckPolicy:
    return NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
        effective_caps={
            (neck_class, passage_state): EngagementCap.build(
                math.radians(90.0 - 10.0 * passage_state.rank),
            )
            for neck_class in range(2)
            for passage_state in ACTIVE_PASSAGE_STATES
        },
    )


def _depletion_policy() -> DepletionPolicy:
    return DepletionPolicy.build(
        chord_bound=ChordBound.build(0.125),
        center_count_limit=4096,
    )


def _cut_plane(
    *,
    cut_z: float = -2.0,
) -> CutPlane:
    return CutPlane.build(
        CutZ.build(cut_z),
        ClearanceZ.build(5.0),
    )


def _axis(
    *,
    sampling_policy: MatSamplingPolicy,
    tool_radius: ToolRadius,
) -> MedialAxis:
    return MedialAxis.build(
        design_boundary=_ring(),
        holes=(),
        tool_radius=tool_radius,
        station_spacing=sampling_policy.station_spacing,
        max_sagitta=sampling_policy.max_sagitta,
        max_refinement_depth=sampling_policy.max_refinement_depth,
    )


def _seed_traversal(
    *,
    axis: MedialAxis,
    neck_policy: NeckPolicy,
    traversal_policy: TraversalPolicy,
) -> MatTraversalState:
    inventory = NeckInventory.build(
        axis=axis,
        policy=neck_policy,
    )
    entry_edge = next(edge for edge in axis.edges if edge.curve_kind == "parabola" and edge.target.point == Point2[WorldXY].build(1.0, 2.0))
    return MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=traversal_policy,
        entry_edge_id=entry_edge.identity,
        entry_node_id=entry_edge.target.identity,
    )


def _candidates(
    *,
    traversal: MatTraversalState,
    candidate_policy: CandidatePolicy,
    user_cap: EngagementCap,
    circle_orientation: CircleOrientation,
    terminal: bool = False,
) -> tuple[MiddleCurveCandidate, ...]:
    active = traversal.active_cursor
    span = MiddleCurveSpan.build(
        axis=traversal.authority.axis,
        cursor_before=active.cursor,
        cursor_limit=active.terminal_cursor,
    )
    full_cap = FullCapDecision.build(
        user_cap=user_cap,
        effective_cap=user_cap,
    )
    return enumerate_middle_curve_candidates(
        span=span,
        policy=candidate_policy,
        circle_orientation=circle_orientation,
        neck_scope=traversal.neck_scope,
        effective_cap_decision=full_cap,
        makes_cursor_terminal_at_limit=terminal,
    )


def _select_launch(
    candidates: tuple[MiddleCurveCandidate, ...],
    *,
    phase_index: int = 2,
) -> MiddleCurveCandidate:
    return next(
        candidate
        for candidate in candidates
        if candidate.spatial_progress == LAUNCH_PROGRESS
        and candidate.guide_radius == LAUNCH_GUIDE_RADIUS
        and candidate.phase_index == phase_index
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )


def _phase_point(
    candidate: MiddleCurveCandidate,
) -> Point2[WorldXY]:
    return Point2[WorldXY].build(
        candidate.motion.center.x + candidate.motion.phase_vector.x,
        candidate.motion.center.y + candidate.motion.phase_vector.y,
    )


def _identity(
    *,
    candidate: MiddleCurveCandidate,
    reachable_domain: ReachableDomain,
    user_cap: EngagementCap,
    sampling_policy: MatSamplingPolicy,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
    cut_plane: CutPlane,
    entry_radius: Fraction = ENTRY_RADIUS,
    process_identity: bytes = b"task-13d-qualified-bore-process-v1",
    evidence_bytes: bytes = b"task-13d-qualified-bore-evidence-v1",
) -> InputIdentity:
    tool_radius = reachable_domain.certificate.tool_radius
    entry = PreclearedEntry.build(
        reachable_domain=reachable_domain,
        center=_phase_point(candidate),
        radius=EntryRadius.build(float(entry_radius)),
        tool_radius=tool_radius,
        cut_plane=cut_plane,
        qualified_bore=QualifiedBore.build(
            cut_plane=cut_plane,
            process_identity=BoreProcessIdentity(process_identity),
            evidence_bytes=evidence_bytes,
        ),
    )
    return InputIdentity.build(
        design_boundary=_ring(),
        holes=(),
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        reachable_domain=reachable_domain,
        entry=entry,
        user_cap=user_cap,
        mat_sampling_policy=sampling_policy,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
    )


def _build_fixture(
    *,
    entry_radius: Fraction = ENTRY_RADIUS,
) -> _LaunchFixture:
    tool_radius = ToolRadius.build(0.5)
    user_cap = _user_cap()
    sampling_policy = _sampling_policy()
    candidate_policy = _candidate_policy()
    neck_policy = _neck_policy(user_cap)
    depletion_policy = _depletion_policy()
    traversal_policy = TraversalPolicy.build(forward_window=4)
    cut_direction_policy = CutDirectionPolicy.build(CutIntent.CLIMB)
    material_side = MaterialSide.OUTSIDE
    axis = _axis(
        sampling_policy=sampling_policy,
        tool_radius=tool_radius,
    )
    traversal = _seed_traversal(
        axis=axis,
        neck_policy=neck_policy,
        traversal_policy=traversal_policy,
    )
    candidates = _candidates(
        traversal=traversal,
        candidate_policy=candidate_policy,
        user_cap=user_cap,
        circle_orientation=cut_direction_policy.circle_orientation(
            material_side,
        ),
    )
    candidate = _select_launch(candidates)
    reachable_domain = ReachableDomain.build(
        design_boundary=_ring(),
        holes=(),
        tool_radius=tool_radius,
    )
    identity = _identity(
        candidate=candidate,
        reachable_domain=reachable_domain,
        user_cap=user_cap,
        sampling_policy=sampling_policy,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
        cut_plane=_cut_plane(),
        entry_radius=entry_radius,
    )
    evaluator = InitialCandidateEvaluator.build(
        input_identity=identity,
        material_side=material_side,
    )
    return _LaunchFixture(
        identity,
        axis,
        traversal,
        candidates,
        candidate,
        evaluator,
    )


def _build_branch_fixture() -> _BranchFixture:
    """Build one exact terminal-launch-to-adjacent-edge contract fixture."""
    tool_radius = ToolRadius.build(float(BRANCH_TOOL_RADIUS))
    user_cap = EngagementCap.build(math.pi)
    sampling_policy = _sampling_policy()
    candidate_policy = _candidate_policy()
    neck_policy = _neck_policy(user_cap)
    depletion_policy = _depletion_policy()
    traversal_policy = TraversalPolicy.build(forward_window=4)
    cut_direction_policy = CutDirectionPolicy.build(CutIntent.CLIMB)
    material_side = MaterialSide.OUTSIDE
    axis = _axis(
        sampling_policy=sampling_policy,
        tool_radius=tool_radius,
    )
    inventory = NeckInventory.build(
        axis=axis,
        policy=neck_policy,
    )
    entry_edge = axis.edges[0]
    traversal = MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=traversal_policy,
        entry_edge_id=entry_edge.identity,
        entry_node_id=entry_edge.source.identity,
    )
    launch_candidates = _candidates(
        traversal=traversal,
        candidate_policy=candidate_policy,
        user_cap=user_cap,
        circle_orientation=cut_direction_policy.circle_orientation(
            material_side,
        ),
        terminal=True,
    )
    launch = next(
        candidate
        for candidate in launch_candidates
        if candidate.traversal_decision.makes_cursor_terminal
        and candidate.guide_radius == BRANCH_LAUNCH_GUIDE_RADIUS
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    identity = _identity(
        candidate=launch,
        reachable_domain=ReachableDomain.build(
            design_boundary=_ring(),
            holes=(),
            tool_radius=tool_radius,
        ),
        user_cap=user_cap,
        sampling_policy=sampling_policy,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
        cut_plane=_cut_plane(),
        entry_radius=BRANCH_ENTRY_RADIUS,
        process_identity=b"task-13e-branch-bore-process-v1",
        evidence_bytes=b"task-13e-branch-bore-evidence-v1",
    )
    initial_evaluator = InitialCandidateEvaluator.build(
        input_identity=identity,
        material_side=material_side,
    )
    launch_transaction = initial_evaluator.evaluate(
        traversal,
        launch,
    )
    physical, traversal_after = initial_evaluator.commit(
        traversal,
        launch_transaction,
    )
    branch_parent = traversal_after.activate_next()
    active = branch_parent.active_cursor
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=active.cursor,
        cursor_limit=active.terminal_cursor,
    )
    full_cap = FullCapDecision.build(
        user_cap=user_cap,
        effective_cap=user_cap,
    )
    branch_candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=candidate_policy,
        circle_orientation=cut_direction_policy.circle_orientation(
            material_side,
        ),
        neck_scope=branch_parent.neck_scope,
        effective_cap_decision=full_cap,
        makes_cursor_terminal_at_limit=True,
    )
    candidate = next(
        candidate
        for candidate in branch_candidates
        if candidate.spatial_progress == BRANCH_PROGRESS
        and candidate.guide_radius == BRANCH_GUIDE_RADIUS
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    alternatives = tuple(
        alternative
        for alternative in branch_candidates
        if alternative.spatial_progress == candidate.spatial_progress and alternative.guide_radius == candidate.guide_radius and alternative.identity != candidate.identity
    )
    evaluator = CandidateEvaluator.build(
        reachable_domain=identity.reachable_domain,
        tool_radius=identity.tool_radius,
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
        neck_policy=identity.neck_policy,
        depletion_policy=identity.depletion_policy,
        cut_direction_policy=identity.cut_direction_policy,
        cut_z=identity.cut_plane.cut_z,
        material_side=material_side,
    )
    return _BranchFixture(
        identity=identity,
        seeded_traversal=traversal,
        initial_evaluator=initial_evaluator,
        launch_transaction=launch_transaction,
        physical=physical,
        traversal=branch_parent,
        candidate=candidate,
        alternatives=alternatives,
        evaluator=evaluator,
    )


@pytest.fixture(scope="module")
def launch() -> _LaunchFixture:
    return _build_fixture()


@pytest.fixture(scope="module")
def branch() -> _BranchFixture:
    return _build_branch_fixture()


def _terminal_candidate(
    state: MatTraversalState,
    fixture: _LaunchFixture,
) -> MiddleCurveCandidate:
    candidates = _candidates(
        traversal=state,
        candidate_policy=fixture.identity.candidate_policy,
        user_cap=fixture.identity.user_cap,
        circle_orientation=fixture.identity.cut_direction_policy.circle_orientation(
            MaterialSide.OUTSIDE,
        ),
        terminal=True,
    )
    active = state.active_cursor
    span = MiddleCurveSpan.build(
        axis=state.authority.axis,
        cursor_before=active.cursor,
        cursor_limit=active.terminal_cursor,
    )
    return next(candidate for candidate in candidates if candidate.spatial_progress == span.reported_length)


def test_entry_launch_evaluates_circle_only_and_commits_independently(
    launch: _LaunchFixture,
) -> None:
    """Start physical cutting with one proved circle, never a fake link.

    A launch differs from continuation because its phase is already the
    qualified bore center. The accepted prefix must therefore be approach,
    plunge, and a real circle whose stock/coverage witnesses begin immediately
    after the one-and-only entry depletion.
    """
    parent_bytes = launch.traversal.canonical_bytes

    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    state, traversal_after = launch.evaluator.commit(
        launch.traversal,
        transaction,
    )

    assert type(transaction) is InitialCandidateTransaction
    assert type(state) is GenerationState
    assert launch.traversal.canonical_bytes == parent_bytes
    assert transaction.traversal_before == launch.traversal
    assert transaction.traversal_after == traversal_after
    assert transaction.traversal_after == launch.traversal.advance(
        launch.candidate,
    )
    assert transaction.entry_depletion_witness.parent_lineage == ()
    assert transaction.initial_coverage_certificate.ordered_sweep_records == ()
    assert transaction.entry_circle_certificate.motion == launch.candidate.motion
    assert transaction.circle_witness.operation_index == 2
    assert transaction.circle_witness.depletion_witness.parent_lineage == (transaction.entry_depletion_witness.digest,)
    assert transaction.circle_witness.sweep_witness.parent_lineage == ()
    assert transaction.result_state_digest == state.digest
    assert state.operations[:2] == (
        launch.identity.entry.approach,
        launch.identity.entry.plunge,
    )
    assert type(state.operations[2]) is CutFullCircleOperation
    assert not any(type(operation) is LinkSegmentOperation for operation in state.operations)
    assert state.phase_point == launch.identity.entry.center
    assert all(passage.state is PassageState.UNVISITED for passage in state.passages)


def test_entry_launch_repeats_byte_identical_without_mutating_authority(
    launch: _LaunchFixture,
) -> None:
    """Make trial evidence deterministic and alias-safe.

    Candidate ranking may evaluate several launch cells. Repeating one trial
    must reproduce the complete transaction while leaving the global graph
    ledger at its exact parent digest until an explicit commit is selected.
    """
    parent_digest = launch.traversal.digest

    first = launch.evaluator.evaluate(launch.traversal, launch.candidate)
    second = launch.evaluator.evaluate(launch.traversal, launch.candidate)

    assert first.canonical_bytes == second.canonical_bytes
    assert first.digest == second.digest
    assert launch.traversal.digest == parent_digest


def test_entry_launch_requires_exact_phase_equality(
    launch: _LaunchFixture,
) -> None:
    """Reject a geometrically nearby circle whose phase is not the bore center.

    Phase is a causal operation-boundary identity, not a proximity test.
    Another finite phase of the same proposal therefore cannot inherit the
    qualified entry merely because both coordinates are representable.
    """
    wrong_phase = _select_launch(
        launch.candidates,
        phase_index=1,
    )

    with pytest.raises(
        InitialCandidatePhaseError,
        match="qualified entry center",
    ):
        launch.evaluator.evaluate(
            launch.traversal,
            wrong_phase,
        )


def test_entry_launch_rejects_circle_outside_qualified_void() -> None:
    """Keep entry containment separate from design-domain containment.

    The same circle remains safely inside the L pocket when the declared bore
    radius shrinks. It must still fail because design containment cannot prove
    that the cutter swept only already-cleared material.
    """
    undersized = _build_fixture(
        entry_radius=Fraction(11, 20),
    )

    with pytest.raises(GougeContainmentError):
        undersized.evaluator.evaluate(
            undersized.traversal,
            undersized.candidate,
        )


def test_entry_launch_has_no_zero_length_synthetic_link(
    launch: _LaunchFixture,
) -> None:
    """Make the missing launch link a closed grammar fact, not a convention.

    The entry center is already the candidate phase, so a link would have zero
    progress. The exact motion type rejects that segment before it can acquire
    operation, stock, coverage, or transaction identity.
    """
    with pytest.raises(
        DegenerateSegmentMotionError,
        match="nonzero exact progress",
    ):
        ExactSegmentMotion.build(
            launch.identity.entry.center,
            launch.identity.entry.center,
        )


def test_entry_launch_rejects_foreign_mat_sampling_authority(
    launch: _LaunchFixture,
) -> None:
    """Bind the launch to the MAT resolution named by `InputIdentity`.

    A newly built axis may report the same exact topology under another sample
    spacing. Its cursor identities and candidate grammar are nevertheless a
    different planning input and cannot be substituted at acceptance time.
    """
    foreign_axis = _axis(
        sampling_policy=_sampling_policy(station_spacing=0.5),
        tool_radius=launch.identity.tool_radius,
    )
    foreign_traversal = _seed_traversal(
        axis=foreign_axis,
        neck_policy=launch.identity.neck_policy,
        traversal_policy=launch.identity.traversal_policy,
    )

    with pytest.raises(
        InitialCandidateStateMismatchError,
        match="sampling",
    ):
        launch.evaluator.evaluate(
            foreign_traversal,
            launch.candidate,
        )


def test_entry_launch_rejects_foreign_tool_authority(
    launch: _LaunchFixture,
) -> None:
    """Reject a valid MAT graph built for another cutter radius.

    Exact edge identities include the native MAT build, but the launch also
    checks the framed typed tool owner directly so a changed-radius graph fails
    before any candidate or stock proof is attempted.
    """
    foreign_tool = ToolRadius.build(1.0)
    foreign_axis = _axis(
        sampling_policy=launch.identity.mat_sampling_policy,
        tool_radius=foreign_tool,
    )
    foreign_inventory = NeckInventory.build(
        axis=foreign_axis,
        policy=launch.identity.neck_policy,
    )
    entry_edge = foreign_axis.edges[0]
    foreign_traversal = MatTraversalState.seed(
        axis=foreign_axis,
        inventory=foreign_inventory,
        policy=launch.identity.traversal_policy,
        entry_edge_id=entry_edge.identity,
        entry_node_id=entry_edge.source.identity,
    )

    with pytest.raises(
        InitialCandidateStateMismatchError,
        match="tool",
    ):
        launch.evaluator.evaluate(
            foreign_traversal,
            launch.candidate,
        )


def test_entry_launch_rejects_stale_global_parent(
    launch: _LaunchFixture,
) -> None:
    """Prevent a valid trial from committing after its graph parent changes.

    Stock is reconstructed independently during commit, so the global parent
    digest is the stale-state guard that prevents two accepted candidates from
    consuming the same cursor boundary.
    """
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    changed_parent = launch.traversal.advance(launch.candidate)

    with pytest.raises(
        StaleInitialCandidateTransactionError,
        match="parent",
    ):
        launch.evaluator.commit(
            changed_parent,
            transaction,
        )


def test_entry_launch_rejects_foreign_qualified_entry_witness(
    launch: _LaunchFixture,
) -> None:
    """Keep bore process evidence in the launch transaction identity.

    Two entries can have identical XY geometry and cut interval while naming
    different metrology evidence. Their depletion records are not
    interchangeable because only one entry digest owns the circle-in-void
    certificate.
    """
    foreign_identity = _identity(
        candidate=launch.candidate,
        reachable_domain=launch.identity.reachable_domain,
        user_cap=launch.identity.user_cap,
        sampling_policy=launch.identity.mat_sampling_policy,
        candidate_policy=launch.identity.candidate_policy,
        neck_policy=launch.identity.neck_policy,
        depletion_policy=launch.identity.depletion_policy,
        traversal_policy=launch.identity.traversal_policy,
        cut_direction_policy=launch.identity.cut_direction_policy,
        cut_plane=launch.identity.cut_plane,
        process_identity=b"task-13d-foreign-bore-process-v1",
        evidence_bytes=b"task-13d-foreign-bore-evidence-v1",
    )
    foreign_evaluator = InitialCandidateEvaluator.build(
        input_identity=foreign_identity,
        material_side=MaterialSide.OUTSIDE,
    )
    foreign = foreign_evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )

    with pytest.raises(
        InvalidInitialCandidateTransactionError,
        match="entry-circle proof",
    ):
        replace(
            transaction,
            entry_depletion_witness=foreign.entry_depletion_witness,
        )


def test_entry_launch_rejects_foreign_design_domain_certificate(
    launch: _LaunchFixture,
) -> None:
    """Reject a valid containment proof issued by another exact domain.

    The candidate sweep fits inside both the L pocket and a large rectangle.
    Geometry equality alone is therefore insufficient: the containment
    certificate must name the reachable-domain digest authenticated by entry.
    """
    foreign_boundary = CanonicalRingV1.build_outer(
        (
            Point2[WorldXY].build(0.0, 0.0),
            Point2[WorldXY].build(8.0, 0.0),
            Point2[WorldXY].build(8.0, 8.0),
            Point2[WorldXY].build(0.0, 8.0),
        )
    )
    foreign_domain = ReachableDomain.build(
        design_boundary=foreign_boundary,
        holes=(),
        tool_radius=launch.identity.tool_radius,
    )
    foreign_containment = GougeContainment.build(
        foreign_domain,
    ).certify_full_circle(
        launch.candidate.motion,
        launch.identity.tool_radius,
    )
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    foreign_circle = replace(
        transaction.circle_witness,
        containment_certificate=foreign_containment,
    )

    with pytest.raises(
        InvalidInitialCandidateTransactionError,
        match="foreign domain",
    ):
        replace(
            transaction,
            circle_witness=foreign_circle,
        )


def test_entry_launch_rejects_certification_after_circle_depletion(
    launch: _LaunchFixture,
) -> None:
    """Reject a proof that observes stock only after removing its own sweep.

    Certifying after depletion can manufacture zero engagement for an unsafe
    motion. The transaction therefore binds the motion witness to the lineage
    containing the entry depletion and excludes the candidate depletion.
    """
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    stock = Stock2Area.build(Stock(_pocket(), []))
    stock.deplete(launch.identity.entry)
    stock.deplete(
        launch.candidate.motion,
        launch.identity.tool_radius,
        launch.identity.depletion_policy,
    )
    late_witness = MotionCertifier.build(
        stock=stock,
        tool_radius=launch.identity.tool_radius,
    ).certify(
        operation_index=2,
        operation_kind=OperationType.CUT,
        motion=launch.candidate.motion,
        user_cap=launch.identity.user_cap,
        effective_cap=launch.identity.user_cap,
    )
    reordered_circle = replace(
        transaction.circle_witness,
        motion_witness=late_witness,
    )

    with pytest.raises(
        InvalidInitialCandidateTransactionError,
        match="certify-before-deplete",
    ):
        replace(
            transaction,
            circle_witness=reordered_circle,
        )


def test_entry_launch_rejects_cross_wired_circle_witness(
    launch: _LaunchFixture,
) -> None:
    """Prevent an internally valid witness from authenticating another circle.

    Opposite material side changes cutter rotation while retaining the same
    geometric phase. Its proof bundle is valid on its own, which makes this a
    stronger lineage test than corrupting a field inside one witness.
    """
    inside_candidates = _candidates(
        traversal=launch.traversal,
        candidate_policy=launch.identity.candidate_policy,
        user_cap=launch.identity.user_cap,
        circle_orientation=launch.identity.cut_direction_policy.circle_orientation(
            MaterialSide.INSIDE,
        ),
    )
    inside_candidate = _select_launch(inside_candidates)
    inside_evaluator = InitialCandidateEvaluator.build(
        input_identity=launch.identity,
        material_side=MaterialSide.INSIDE,
    )
    inside_transaction = inside_evaluator.evaluate(
        launch.traversal,
        inside_candidate,
    )
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )

    with pytest.raises(
        InvalidInitialCandidateTransactionError,
        match="candidate",
    ):
        replace(
            transaction,
            circle_witness=inside_transaction.circle_witness,
        )


def test_entry_launch_rejects_foreign_cut_depth_witness(
    launch: _LaunchFixture,
) -> None:
    """Bind the first lateral operation to the qualified plunge depth.

    XY geometry and engagement evidence are unchanged by cut-Z, so a witness
    from another qualified plane can look geometrically valid. Transaction
    closure must still reject the cross-wired manufacturing interval.
    """
    foreign_identity = _identity(
        candidate=launch.candidate,
        reachable_domain=launch.identity.reachable_domain,
        user_cap=launch.identity.user_cap,
        sampling_policy=launch.identity.mat_sampling_policy,
        candidate_policy=launch.identity.candidate_policy,
        neck_policy=launch.identity.neck_policy,
        depletion_policy=launch.identity.depletion_policy,
        traversal_policy=launch.identity.traversal_policy,
        cut_direction_policy=launch.identity.cut_direction_policy,
        cut_plane=_cut_plane(cut_z=-3.0),
    )
    foreign_evaluator = InitialCandidateEvaluator.build(
        input_identity=foreign_identity,
        material_side=MaterialSide.OUTSIDE,
    )
    foreign = foreign_evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )

    with pytest.raises(
        InvalidInitialCandidateTransactionError,
        match="cut depth",
    ):
        replace(
            transaction,
            circle_witness=foreign.circle_witness,
        )


def test_entry_launch_rejects_child_with_two_advanced_cursors(
    launch: _LaunchFixture,
) -> None:
    """Require the launch to consume exactly one global edge cursor.

    A later valid traversal snapshot can contain two advanced cursors and pass
    its own graph invariants. It still cannot be substituted as the child of
    one launch candidate.
    """
    first_terminal = launch.traversal.advance(
        _terminal_candidate(launch.traversal, launch),
    )
    second_parent = first_terminal.activate_next()
    two_advanced = second_parent.advance(
        _terminal_candidate(second_parent, launch),
    )
    assert sum(cursor.accepted_candidate_count > 0 for cursor in two_advanced.cursors) == 2
    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )

    with pytest.raises(
        InvalidInitialCandidateTransactionError,
        match="exactly one",
    ):
        replace(
            transaction,
            traversal_after=two_advanced,
        )


def test_branch_switch_commit_binds_both_state_axes(
    branch: _BranchFixture,
) -> None:
    """Certify one real adjacent-edge switch through the Task 12 proof path.

    The terminal launch leaves the physical cursor on its completed edge,
    while deterministic route activation authorizes another exact edge. The
    continuation must certify one direct link then one circle and atomically
    bind both physical and global parent/child identities.
    """
    physical_edge_before = branch.physical.traversal.edge_id
    global_edge_before = branch.traversal.active_cursor.route_step.edge_id

    transaction = evaluate_traversal_candidate(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
        candidate=branch.candidate,
    )
    physical_after, traversal_after, commit = commit_traversal_candidate(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
        transaction=transaction,
    )

    changed_cursors = tuple(
        index
        for index, (before, after) in enumerate(
            zip(
                branch.traversal.cursors,
                traversal_after.cursors,
                strict=True,
            )
        )
        if before.canonical_bytes != after.canonical_bytes
    )
    assert physical_edge_before != global_edge_before
    assert type(transaction) is CandidateTransaction
    assert type(commit) is TraversalCommit
    assert transaction.link_witness.operation_index + 1 == transaction.circle_witness.operation_index
    assert transaction.parent_state_digest == branch.physical.digest
    assert transaction.result_state_digest == physical_after.digest
    assert physical_after.traversal == transaction.traversal_after
    assert commit.physical_parent_digest == branch.physical.digest
    assert commit.physical_child_digest == physical_after.digest
    assert commit.traversal_before == branch.traversal
    assert commit.traversal_after == traversal_after
    assert commit.transaction == transaction
    assert changed_cursors == (branch.traversal.active_route_index,)
    assert traversal_after.active_cursor.accepted_candidate_count == 1
    assert branch.traversal.pending_transit is None


def test_branch_switch_commit_rejects_stale_and_cross_wired_axes(
    branch: _BranchFixture,
) -> None:
    """Reject independently valid physical or global children from other steps.

    Content-addressing is useful only when the commit closes the square. A
    stale global parent, the unchanged physical parent, or a global child
    advanced by another exact candidate must fail before it can acquire one
    traversal-commit identity.
    """
    transaction = evaluate_traversal_candidate(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
        candidate=branch.candidate,
    )
    physical_after, traversal_after, _ = commit_traversal_candidate(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
        transaction=transaction,
    )
    stale_global = branch.traversal.advance(branch.candidate)

    with pytest.raises(
        InvalidTraversalCommitError,
        match="physical child",
    ):
        TraversalCommit.build(
            physical_before=branch.physical,
            traversal_before=branch.traversal,
            transaction=transaction,
            physical_after=branch.physical,
            traversal_after=traversal_after,
        )
    with pytest.raises(
        InvalidTraversalCommitError,
        match="global child",
    ):
        TraversalCommit.build(
            physical_before=branch.physical,
            traversal_before=branch.traversal,
            transaction=transaction,
            physical_after=physical_after,
            traversal_after=branch.traversal.advance(
                branch.alternatives[0],
            ),
        )
    with pytest.raises(
        InvalidTraversalCommitError,
        match="global parent",
    ):
        commit_traversal_candidate(
            evaluator=branch.evaluator,
            physical=branch.physical,
            traversal=stale_global,
            transaction=transaction,
        )


def _ordered_branch_family(
    branch: _BranchFixture,
) -> tuple[MiddleCurveCandidate, ...]:
    return branch.candidate.policy.order_candidates(
        (branch.candidate, *branch.alternatives),
        key=lambda candidate: candidate.order_key,
    )


def _first_causal_neck_family(
    branch: _BranchFixture,
) -> tuple[MatTraversalState, tuple[MiddleCurveCandidate, ...]]:
    traversal = branch.traversal
    while traversal.pending_transit is None:
        terminal = next(
            candidate
            for candidate in _candidates(
                traversal=traversal,
                candidate_policy=branch.identity.candidate_policy,
                user_cap=branch.identity.user_cap,
                circle_orientation=branch.identity.cut_direction_policy.circle_orientation(
                    branch.evaluator.material_side,
                ),
                terminal=True,
            )
            if candidate.traversal_decision.makes_cursor_terminal
        )
        traversal = traversal.advance(terminal).activate_next()
    passage = traversal.pending_transit.passage(branch.physical.passages)
    cap = passage.propose_cap_decision(branch.identity.neck_policy)
    active = traversal.active_cursor
    return traversal, enumerate_middle_curve_candidates(
        span=MiddleCurveSpan.build(
            axis=traversal.authority.axis,
            cursor_before=active.cursor,
            cursor_limit=active.terminal_cursor,
        ),
        policy=branch.identity.candidate_policy,
        circle_orientation=branch.identity.cut_direction_policy.circle_orientation(
            branch.evaluator.material_side,
        ),
        neck_scope=traversal.neck_scope,
        effective_cap_decision=cap,
        makes_cursor_terminal_at_limit=True,
    )


def test_finite_search_uses_order_and_stops_after_first_feasible(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Continue only proved local rejection, then stop at the dominant winner."""
    family = _ordered_branch_family(branch)
    winner = family[4]
    accepted = branch.evaluator.evaluate_from_cursor(
        branch.physical,
        TraversalCursorState.before(winner.traversal_decision),
        winner,
    )
    calls: list[bytes] = []
    failures = (
        GougeContainmentError("proved gouge"),
        EngagementCapExceededError("proved cap"),
        DegenerateSegmentMotionError("proved zero link"),
        EngagementCapExceededError("proved cap"),
    )

    def controlled_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        calls.append(candidate.identity)
        index = family.index(candidate)
        if index < len(failures):
            raise failures[index]
        if candidate == winner:
            return accepted
        raise AssertionError("lower-ranked candidate evaluated after winner")

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_from_cursor",
        controlled_trial,
    )

    selected = evaluate_first_feasible_candidate(
        evaluator=branch.evaluator,
        physical=branch.physical,
        traversal=branch.traversal,
        candidates=family,
    )

    assert selected == accepted
    assert calls == [candidate.identity for candidate in family[:5]]


def test_finite_search_rejects_noncanonical_family_order(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject reordered materialization before dispatching any exact proof."""
    calls: list[bytes] = []

    def forbidden_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        calls.append(candidate.identity)
        raise AssertionError("invalid family reached the proof engine")

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_from_cursor",
        forbidden_trial,
    )

    with pytest.raises(InvalidCandidateFamilyError, match="invariant order"):
        evaluate_first_feasible_candidate(
            evaluator=branch.evaluator,
            physical=branch.physical,
            traversal=branch.traversal,
            candidates=tuple(reversed(_ordered_branch_family(branch))),
        )

    assert calls == []


def test_finite_search_propagates_unresolved_authority_immediately(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Never convert an incomplete exact event partition into infeasibility."""
    family = _ordered_branch_family(branch)
    calls: list[bytes] = []

    def unresolved_second_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        calls.append(candidate.identity)
        if candidate == family[0]:
            raise GougeContainmentError("proved gouge")
        raise UnresolvedMotionEventError("unresolved exact event partition")

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_from_cursor",
        unresolved_second_trial,
    )

    with pytest.raises(
        UnresolvedMotionEventError,
        match="unresolved exact event",
    ):
        evaluate_first_feasible_candidate(
            evaluator=branch.evaluator,
            physical=branch.physical,
            traversal=branch.traversal,
            candidates=family,
        )

    assert calls == [family[0].identity, family[1].identity]


@pytest.mark.parametrize(
    ("mixed", "error_type"),
    (
        (False, EngagementCapInfeasibleError),
        (True, NoFeasibleCandidateError),
    ),
)
def test_finite_search_names_non_neck_exhaustion(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
    mixed: bool,
    error_type: type[RuntimeError],
) -> None:
    """Distinguish cap-only process infeasibility from mixed exact rejection."""
    family = _ordered_branch_family(branch)

    def reject_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        if mixed and candidate == family[0]:
            raise GougeContainmentError("proved gouge")
        raise EngagementCapExceededError("proved cap")

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_from_cursor",
        reject_trial,
    )

    with pytest.raises(
        error_type,
        match=rf"attempts={len(family)}; cap={len(family) - int(mixed)}; "
        rf"gouge={int(mixed)}; degenerate-link=0",
    ):
        evaluate_first_feasible_candidate(
            evaluator=branch.evaluator,
            physical=branch.physical,
            traversal=branch.traversal,
            candidates=family,
        )


def test_finite_search_names_causal_neck_cap_exhaustion(
    branch: _BranchFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Report cap-only failure under authenticated neck scope as too tight."""
    traversal, family = _first_causal_neck_family(branch)

    def reject_cap(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        raise EngagementCapExceededError("proved neck cap")

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_from_cursor",
        reject_cap,
    )

    with pytest.raises(
        NeckTooTightError,
        match=rf"attempts={len(family)}; cap={len(family)}; gouge=0; "
        r"degenerate-link=0",
    ):
        evaluate_first_feasible_candidate(
            evaluator=branch.evaluator,
            physical=branch.physical,
            traversal=traversal,
            candidates=family,
        )


def test_seeded_entry_launch_has_no_causal_passage_to_advance(
    launch: _LaunchFixture,
) -> None:
    """Keep first-side establishment distinct from a later neck crossing.

    A seeded route has no previously occupied neck side and therefore no
    pending transit. Launch commit advances one MAT cursor, while every
    independent oriented passage remains `UNVISITED` for later causal branch
    transitions.
    """
    assert launch.traversal.pending_transit is None
    assert type(launch.traversal.neck_scope) is NoNeckScope

    transaction = launch.evaluator.evaluate(
        launch.traversal,
        launch.candidate,
    )
    state, _ = launch.evaluator.commit(
        launch.traversal,
        transaction,
    )

    assert type(transaction.circle_witness.operation.neck_scope) is NoNeckScope
    assert all(passage.state is PassageState.UNVISITED for passage in state.passages)
