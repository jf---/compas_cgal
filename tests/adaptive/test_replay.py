import math
from fractions import Fraction

import pytest
from compas.geometry import Polygon

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import ReplayCandidateError
from compas_cgal.adaptive.errors import ReplayCutDirectionError
from compas_cgal.adaptive.errors import ReplayGrammarError
from compas_cgal.adaptive.errors import ReplayInputMismatchError
from compas_cgal.adaptive.errors import ReplayTraversalError
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MatSamplingPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.replay import replay_certificate
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

OUTER = (
    (0.0, 0.0),
    (6.0, 0.0),
    (6.0, 2.0),
    (2.0, 2.0),
    (2.0, 6.0),
    (0.0, 6.0),
)


def _pocket() -> Polygon:
    return Polygon([[x, y, 0.0] for x, y in OUTER])


def _ring() -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(tuple(Point2[WorldXY].build(x, y) for x, y in OUTER))


def _candidate_policy(*, spatial_resolution: float = 0.5) -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(spatial_resolution),
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


def _lattice_candidate(
    *,
    user_cap: EngagementCap,
    candidate_policy: CandidatePolicy,
) -> MiddleCurveCandidate:
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
    cap = FullCapDecision.build(
        user_cap=user_cap,
        effective_cap=user_cap,
    )
    return next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=cap,
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )


def _input_identity() -> InputIdentity:
    design_boundary = _ring()
    tool_radius = ToolRadius.build(0.5)
    candidate_policy = _candidate_policy()
    user_cap = EngagementCap.build(math.radians(120.0))
    candidate = _lattice_candidate(
        user_cap=user_cap,
        candidate_policy=candidate_policy,
    )
    phase = Point2[WorldXY].build(
        candidate.motion.center.x + candidate.motion.phase_vector.x,
        candidate.motion.center.y + candidate.motion.phase_vector.y,
    )
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
                math.radians(90.0 - 10.0 * passage_state.rank),
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


def _circle_operation(
    identity: InputIdentity,
    *,
    clockwise: bool,
) -> CutFullCircleOperation:
    candidate = _lattice_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    motion = candidate.motion
    if clockwise != motion.clockwise:
        motion = ExactCircleMotion.build(
            motion.center,
            motion.phase_vector,
            clockwise,
        )
    return CutFullCircleOperation.build(
        motion=motion,
        cut_z=identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
        neck_scope=candidate.neck_scope,
        effective_cap_decision=candidate.effective_cap_decision,
        traversal_decision=candidate.traversal_decision,
    )


def _replay(
    identity: InputIdentity,
    operations: tuple[
        object,
        ...,
    ],
    *,
    candidate_policy: CandidatePolicy | None = None,
) -> None:
    replay_certificate(
        input_identity=identity,
        pocket=_pocket(),
        holes=(),
        cut_plane=identity.cut_plane,
        tool_radius=identity.tool_radius,
        user_cap=identity.user_cap,
        entry=identity.entry,
        operations=operations,  # type: ignore[arg-type]
        candidate_policy=(identity.candidate_policy if candidate_policy is None else candidate_policy),
        neck_policy=identity.neck_policy,
        depletion_policy=identity.depletion_policy,
        traversal_policy=identity.traversal_policy,
        cut_direction_policy=identity.cut_direction_policy,
    )


def test_replay_rejects_grammar_before_rebuilding_exact_state() -> None:
    identity = _input_identity()

    with pytest.raises(ReplayGrammarError, match="approach, plunge"):
        _replay(
            identity,
            (
                identity.entry.approach,
                _circle_operation(identity, clockwise=False),
            ),
        )


def test_replay_rejects_policy_argument_that_differs_from_input_root() -> None:
    identity = _input_identity()

    with pytest.raises(ReplayInputMismatchError, match="candidate policy"):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _circle_operation(identity, clockwise=False),
            ),
            candidate_policy=_candidate_policy(spatial_resolution=0.25),
        )


def test_replay_rejects_orientation_not_caused_by_bound_material_side() -> None:
    identity = _input_identity()

    with pytest.raises(ReplayCutDirectionError, match="material side"):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _circle_operation(identity, clockwise=True),
            ),
        )


def test_replay_rejects_circle_without_one_fresh_lattice_candidate() -> None:
    identity = _input_identity()
    circle = _circle_operation(identity, clockwise=False)
    traversal = circle.traversal_decision
    mutated = CutFullCircleOperation.build(
        motion=circle.motion,
        cut_z=circle.cut_z,
        material_side=circle.material_side,
        neck_scope=circle.neck_scope,
        effective_cap_decision=circle.effective_cap_decision,
        traversal_decision=AdvanceTraversalDecision.build(
            component_id=traversal.component_id,
            edge_id=traversal.edge_id,
            branch_id=traversal.branch_id,
            cursor_before=traversal.cursor_before,
            cursor_after=CursorIdentity(b"\xff" * 32),
            makes_cursor_terminal=traversal.makes_cursor_terminal,
        ),
    )

    with pytest.raises(
        ReplayCandidateError,
        match="unique finite-lattice candidate",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                mutated,
            ),
        )


def test_replay_reconstructs_unique_candidate_before_stock_replay() -> None:
    identity = _input_identity()

    with pytest.raises(
        ReplayTraversalError,
        match="stock and coverage replay",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _circle_operation(identity, clockwise=False),
            ),
        )
