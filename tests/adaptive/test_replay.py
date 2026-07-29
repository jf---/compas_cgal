import hashlib
import math
from dataclasses import replace
from fractions import Fraction

import pytest
from compas.geometry import Polygon

import compas_cgal.adaptive.replay as replay_module
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import ReplayCandidateError
from compas_cgal.adaptive.errors import ReplayCutDirectionError
from compas_cgal.adaptive.errors import ReplayEffectiveCapError
from compas_cgal.adaptive.errors import ReplayGrammarError
from compas_cgal.adaptive.errors import ReplayInputMismatchError
from compas_cgal.adaptive.errors import ReplayPairingError
from compas_cgal.adaptive.errors import ReplayTraversalError
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckOwnerId
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import WidthClassId
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
from compas_cgal.adaptive.replay_trace import FreshReplayTrace
from compas_cgal.adaptive.stock_area import Stock2Area
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
from compas_cgal.toolpath import OperationType

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


def _terminal_lattice_candidate(
    *,
    identity: InputIdentity,
) -> MiddleCurveCandidate:
    axis = _axis()
    first = _lattice_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    edge = next(edge for edge in axis.edges if bytes(edge.identity) == bytes(first.traversal_decision.edge_id))
    samples = tuple(
        sorted(
            (sample for sample in axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    first_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[1],
        cursor_limit=samples[2],
    )
    cursor = DerivedCandidateCursor.build(
        span=first_span,
        candidate=first,
    )
    terminal_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=cursor,
        cursor_limit=samples[2],
    )
    cap = FullCapDecision.build(
        user_cap=identity.user_cap,
        effective_cap=identity.user_cap,
    )
    return next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=terminal_span,
            policy=identity.candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=cap,
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 3
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )


def _neck_lattice_candidates(
    *,
    identity: InputIdentity,
) -> tuple[MiddleCurveCandidate, MiddleCurveCandidate]:
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
    neck = NeckInventory.build(
        axis=axis,
        policy=identity.neck_policy,
    ).necks[0]
    first_passage = neck.forward
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
    terminal_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=DerivedCandidateCursor.build(
            span=first_span,
            candidate=first,
        ),
        cursor_limit=samples[2],
    )
    second = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=terminal_span,
            policy=identity.candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=second_passage.scope,
            effective_cap_decision=second_decision,
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 3
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    return first, second


def _input_identity(
    *,
    neck_cap_degrees: tuple[float, float, float] = (90.0, 80.0, 70.0),
) -> InputIdentity:
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


def _candidate_circle_operation(
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


def _phase_point(
    candidate: MiddleCurveCandidate,
) -> Point2[WorldXY]:
    return Point2[WorldXY].build(
        candidate.motion.center.x + candidate.motion.phase_vector.x,
        candidate.motion.center.y + candidate.motion.phase_vector.y,
    )


def _link_to_candidate(
    identity: InputIdentity,
    *,
    start: Point2[WorldXY],
    candidate: MiddleCurveCandidate,
) -> LinkSegmentOperation:
    traversal = candidate.traversal_decision
    return LinkSegmentOperation.build(
        motion=ExactSegmentMotion.build(
            start,
            _phase_point(candidate),
        ),
        cut_z=identity.cut_plane.cut_z,
        neck_scope=candidate.neck_scope,
        effective_cap_decision=candidate.effective_cap_decision,
        traversal_decision=HoldTraversalDecision.build(
            component_id=traversal.component_id,
            edge_id=traversal.edge_id,
            branch_id=traversal.branch_id,
            cursor=traversal.cursor_before,
        ),
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
    """Reject a missing plunge before constructing any fresh proof owner.

    Approach and plunge authenticate the qualified entry event.  Accepting a
    lateral circle without that exact pair would let replay infer an initial
    void that the submitted operation stream never established.
    """
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
    """Keep the finite candidate lattice inside `InputIdentity`.

    Spatial resolution changes native and derived cursor ownership even when
    the exact MAT topology is unchanged.  Replay must therefore reject an
    alternate policy instead of rebuilding a convenient lattice around the
    submitted circle.
    """
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
    """Require circle orientation to follow its recorded machining cause.

    A clockwise bit alone cannot prove climb or conventional intent.  The
    replayed `MaterialSide` and input-bound `CutDirectionPolicy` must derive
    the submitted orientation independently.
    """
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
    """Reject a circle whose cursor digest is not produced by fresh enumeration.

    The mutation preserves circle geometry, cap, scope, and cursor-before but
    invents cursor-after.  Coordinate agreement must not authorize that new
    traversal lineage.
    """
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
    """Reach the explicit nonterminal boundary through one genuine candidate.

    The submitted first circle is reconstructed from the fresh L-pocket MAT
    and passes entry, containment, engagement, depletion, and coverage.  The
    expected failure is solely that the remaining MAT edges were not traversed;
    replay must not mislabel this partial programme as a certificate.
    """
    identity = _input_identity()

    with pytest.raises(
        ReplayTraversalError,
        match="fresh MAT traversal remains nonterminal",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _circle_operation(identity, clockwise=False),
            ),
        )


def test_replay_reconstructs_unique_reverse_candidate_window() -> None:
    """Recover a reverse P-S candidate from exact finite-lattice lineage.

    Circle orientation describes cutting intent, not MAT graph direction.
    Fresh replay must search the native window on both sides of an initial
    cursor, select the unique cryptographic match, then retain the selected
    reverse direction in its derived continuation cursor.
    """
    identity = _input_identity()
    axis = _axis()
    edge = next(edge for edge in axis.edges if edge.curve_kind == "parabola" and len(tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)) >= 3)
    samples = tuple(
        sorted(
            (sample for sample in axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[2],
        cursor_limit=samples[1],
    )
    full_cap = FullCapDecision.build(
        user_cap=identity.user_cap,
        effective_cap=identity.user_cap,
    )
    candidate = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=identity.candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
    )
    operation = _candidate_circle_operation(identity, candidate)
    current_by_edge = {}

    reconstructed = replay_module._match_circle_candidate(
        axis=axis,
        operation=operation,
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
        neck_policy=identity.neck_policy,
        traversal_policy=identity.traversal_policy,
        cut_direction_policy=identity.cut_direction_policy,
        passages={},
        current_by_edge=current_by_edge,
    )

    assert reconstructed == candidate
    cursor = current_by_edge[edge.identity]
    assert type(cursor) is DerivedCandidateCursor
    assert cursor.ordinal_step == -1

    terminal_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=cursor,
        cursor_limit=samples[0],
    )
    terminal_candidate = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=terminal_span,
            policy=identity.candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.traversal_decision.makes_cursor_terminal and candidate.spatial_progress == terminal_span.reported_length
    )
    terminal_reconstructed = replay_module._match_circle_candidate(
        axis=axis,
        operation=_candidate_circle_operation(
            identity,
            terminal_candidate,
        ),
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
        neck_policy=identity.neck_policy,
        traversal_policy=identity.traversal_policy,
        cut_direction_policy=identity.cut_direction_policy,
        passages={},
        current_by_edge=current_by_edge,
    )

    assert terminal_reconstructed == terminal_candidate
    assert current_by_edge[edge.identity] is None


def test_replay_certifies_before_motion_depletion(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Lock the proof-before-mutation chronology of fresh replay.

    Entry depletion establishes the declared initial void.  Every lateral
    motion must then be certified against its frozen pre-cut stock, depleted
    only after certification, and recorded in coverage last.  Reordering these
    calls can make a cutter appear safe by testing it against material it has
    already removed.
    """
    identity = _input_identity()
    events: list[str] = []
    deplete = Stock2Area.deplete
    certify = MotionCertifier.certify
    add_sweep = CoverageLedger.add_sweep

    def tracked_deplete(
        stock: Stock2Area,
        depletion: object,
        *args: object,
    ) -> object:
        events.append("entry-deplete" if type(depletion) is PreclearedEntry else "motion-deplete")
        return deplete(stock, depletion, *args)  # type: ignore[misc]

    def tracked_certify(
        certifier: MotionCertifier,
        **kwargs: object,
    ) -> object:
        events.append("certify")
        return certify(certifier, **kwargs)  # type: ignore[arg-type]

    def tracked_add_sweep(
        ledger: CoverageLedger,
        motion: object,
        tool_radius: object,
    ) -> object:
        events.append("coverage-sweep")
        return add_sweep(ledger, motion, tool_radius)  # type: ignore[arg-type]

    monkeypatch.setattr(Stock2Area, "deplete", tracked_deplete)
    monkeypatch.setattr(MotionCertifier, "certify", tracked_certify)
    monkeypatch.setattr(CoverageLedger, "add_sweep", tracked_add_sweep)

    with pytest.raises(
        ReplayTraversalError,
        match="fresh MAT traversal remains nonterminal",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _circle_operation(identity, clockwise=False),
            ),
        )

    assert events == [
        "entry-deplete",
        "certify",
        "motion-deplete",
        "coverage-sweep",
    ]


def test_replay_rejects_dangling_link_before_fresh_state_mutation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject transport that has no following circle to own its MAT advance.

    The link is geometrically continuous from the first circle to a real
    terminal lattice candidate and carries that candidate's exact hold cursor.
    Omitting the candidate circle leaves the link without an advancing owner,
    so replay must reject the complete stream before applying even the
    authenticated entry depletion.
    """
    identity = _input_identity()
    first = _lattice_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    terminal = _terminal_lattice_candidate(identity=identity)
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=terminal,
    )

    def unexpected_deplete(
        stock: Stock2Area,
        depletion: object,
        *args: object,
    ) -> object:
        pytest.fail("pairing validation mutated fresh stock")

    monkeypatch.setattr(Stock2Area, "deplete", unexpected_deplete)

    with pytest.raises(
        ReplayPairingError,
        match="following circle",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
                link,
            ),
        )


def test_replay_rejects_link_hold_not_owned_by_following_candidate() -> None:
    """Reject a continuous link whose hold cursor is foreign to its circle.

    Geometry alone cannot pair transport with a candidate: the link must hold
    the exact component, edge, branch, and cursor-before that the immediately
    following circle advances.  This mutation changes only that relational
    proof while leaving both motions untouched.
    """
    identity = _input_identity()
    first = _lattice_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    terminal = _terminal_lattice_candidate(identity=identity)
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=terminal,
    )
    traversal = link.traversal_decision
    foreign_hold = HoldTraversalDecision.build(
        component_id=traversal.component_id,
        edge_id=traversal.edge_id,
        branch_id=traversal.branch_id,
        cursor=CursorIdentity(b"\xff" * 32),
    )
    mutated = LinkSegmentOperation.build(
        motion=link.motion,
        cut_z=link.cut_z,
        neck_scope=link.neck_scope,
        effective_cap_decision=link.effective_cap_decision,
        traversal_decision=foreign_hold,
    )

    with pytest.raises(
        ReplayPairingError,
        match="hold traversal",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
                mutated,
                _candidate_circle_operation(identity, terminal),
            ),
        )


def test_replay_certifies_paired_link_and_circle_before_terminal_check(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Lock exact proof chronology through a derived-cursor link/circle pair.

    The first radius-1/8 circle leaves an authenticated derived cursor.  A
    direct segment transports its phase to a terminal radius-1/8, phase-index-3
    circle generated by the concave point site `(2, 2)`.  Replay must certify,
    deplete, and cover the link and then independently certify the nonuniform
    circle against the post-link stock before mutating it.  Other MAT edges
    remain untouched, so the programme must still fail at the explicit
    nonterminal-traversal boundary rather than emit a partial certificate.
    """
    identity = _input_identity()
    first = _lattice_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    terminal = _terminal_lattice_candidate(identity=identity)
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=terminal,
    )
    events: list[str] = []
    deplete = Stock2Area.deplete
    certify = MotionCertifier.certify
    add_sweep = CoverageLedger.add_sweep

    def tracked_deplete(
        stock: Stock2Area,
        depletion: object,
        *args: object,
    ) -> object:
        if type(depletion) is PreclearedEntry:
            events.append("entry-deplete")
        elif type(depletion) is ExactSegmentMotion:
            events.append("link-deplete")
        else:
            events.append("circle-deplete")
        return deplete(stock, depletion, *args)  # type: ignore[misc]

    def tracked_certify(
        certifier: MotionCertifier,
        **kwargs: object,
    ) -> object:
        operation_kind = kwargs["operation_kind"]
        events.append("link-certify" if operation_kind is OperationType.LINK else "circle-certify")
        return certify(certifier, **kwargs)  # type: ignore[arg-type]

    def tracked_add_sweep(
        ledger: CoverageLedger,
        motion: object,
        tool_radius: object,
    ) -> object:
        events.append("link-coverage" if type(motion) is ExactSegmentMotion else "circle-coverage")
        return add_sweep(ledger, motion, tool_radius)  # type: ignore[arg-type]

    monkeypatch.setattr(Stock2Area, "deplete", tracked_deplete)
    monkeypatch.setattr(MotionCertifier, "certify", tracked_certify)
    monkeypatch.setattr(CoverageLedger, "add_sweep", tracked_add_sweep)

    with pytest.raises(
        ReplayTraversalError,
        match="fresh MAT traversal remains nonterminal",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
                link,
                _candidate_circle_operation(identity, terminal),
            ),
        )

    assert events == [
        "entry-deplete",
        "circle-certify",
        "circle-deplete",
        "circle-coverage",
        "link-certify",
        "link-deplete",
        "link-coverage",
        "circle-certify",
        "circle-deplete",
        "circle-coverage",
    ]


def test_replay_captures_one_ordered_fresh_proof_trace_before_terminal_rejection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Preserve every accepted proof in canonical operation chronology.

    The real three-motion prefix reaches containment, event-exact engagement,
    depletion, and coverage for a circle, its direct link, and the following
    circle.  Only untouched MAT edges make it nonterminal.  Intercepting the
    terminal gate exposes the already-computed trace without weakening that
    rejection or constructing substitute witnesses in the test.
    """
    identity = _input_identity()
    first = _lattice_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    terminal = _terminal_lattice_candidate(identity=identity)
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=terminal,
    )
    operations = (
        identity.entry.approach,
        identity.entry.plunge,
        _candidate_circle_operation(identity, first),
        link,
        _candidate_circle_operation(identity, terminal),
    )
    traces: list[FreshReplayTrace] = []
    require_terminal = replay_module._require_terminal_replay

    def tracked_require_terminal_replay(**kwargs: object) -> None:
        trace = kwargs["trace"]
        assert type(trace) is FreshReplayTrace
        traces.append(trace)
        require_terminal(**kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(
        replay_module,
        "_require_terminal_replay",
        tracked_require_terminal_replay,
    )

    for _ in range(2):
        with pytest.raises(
            ReplayTraversalError,
            match="fresh MAT traversal remains nonterminal",
        ):
            _replay(identity, operations)

    assert len(traces) == 2
    trace = traces[0]
    assert traces[1].canonical_bytes == trace.canonical_bytes
    assert require_canonical_record(trace.canonical_bytes) == trace.canonical_bytes
    assert trace.digest == hashlib.sha256(trace.canonical_bytes).digest()
    assert trace.initial_coverage_certificate.ordered_sweep_records == ()
    assert trace.entry_depletion_witness.motion.canonical_bytes == identity.entry.canonical_bytes
    assert trace.entry_depletion_witness.parent_lineage == ()
    assert [witness.operation_index for witness in trace.lateral_witnesses] == [
        2,
        3,
        4,
    ]
    assert [witness.motion_witness.operation_kind for witness in trace.lateral_witnesses] == [
        OperationType.CUT,
        OperationType.LINK,
        OperationType.CUT,
    ]
    assert trace.terminal_coverage_certificate.ordered_sweep_records == tuple(witness.sweep_witness.canonical_bytes for witness in trace.lateral_witnesses)
    assert trace.first_circle_entry_certificate.motion == trace.lateral_witnesses[0].operation.motion
    assert all(
        witness.operation.motion == witness.containment_certificate.motion == witness.motion_witness.motion == witness.depletion_witness.motion == witness.sweep_witness.motion
        for witness in trace.lateral_witnesses
    )

    with pytest.raises(InvalidReplayTraceError, match="cross-wired exact motions"):
        replace(
            trace.lateral_witnesses[1],
            sweep_witness=trace.lateral_witnesses[0].sweep_witness,
        )
    with pytest.raises(InvalidReplayTraceError, match="operation chronology"):
        replace(
            trace,
            lateral_witnesses=(
                trace.lateral_witnesses[1],
                trace.lateral_witnesses[0],
                trace.lateral_witnesses[2],
            ),
        )


def test_replay_certifies_first_neck_circle_at_reconstructed_cap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pass the freshly reconstructed neck cap to the real motion certifier.

    The first circle lies in the authenticated entry void, so it is physically
    valid at the narrower 90-degree first-passage cap.  Wrapping the real
    certifier makes the consumer boundary observable without replacing its
    engagement proof.
    """
    identity = _input_identity()
    first, _ = _neck_lattice_candidates(identity=identity)
    observed_caps: list[EngagementCap] = []
    certify = MotionCertifier.certify

    def tracked_certify(
        certifier: MotionCertifier,
        **kwargs: object,
    ) -> object:
        effective_cap = kwargs["effective_cap"]
        assert isinstance(effective_cap, EngagementCap)
        observed_caps.append(effective_cap)
        return certify(certifier, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(MotionCertifier, "certify", tracked_certify)

    with pytest.raises(
        ReplayTraversalError,
        match="fresh MAT traversal remains nonterminal",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
            ),
        )

    assert observed_caps == [EngagementCap.build(math.radians(90.0))]


def test_replay_reconstructs_two_oriented_neck_passages_before_motion() -> None:
    """Replay two accepted circles against one freshly owned neck passage.

    The first circle consumes the forward passage's unvisited transition.  Its
    paired link carries the second transition without advancing it, and the
    second circle alone advances first-pass-complete state.  Reaching the known
    nonterminal MAT boundary proves both state-bearing cap decisions survived
    fresh inventory, lattice, pairing, and motion replay.  Equal cap angles
    isolate state ownership from the separate candidate-feasibility problem.
    """
    identity = _input_identity(
        neck_cap_degrees=(120.0, 120.0, 120.0),
    )
    first, second = _neck_lattice_candidates(identity=identity)
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=second,
    )

    with pytest.raises(
        ReplayTraversalError,
        match="fresh MAT traversal remains nonterminal",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
                link,
                _candidate_circle_operation(identity, second),
            ),
        )


def test_replay_preserves_real_second_passage_cap_violation() -> None:
    """Keep passage replay distinct from geometric candidate feasibility.

    The same genuine two-circle lineage is legal in neck state, but its second
    circle exceeds the production policy's 80-degree cap against post-link
    stock.  Replay must preserve that exact negative instead of treating a valid
    state transition as motion acceptance or silently certifying at user cap.
    """
    identity = _input_identity()
    first, second = _neck_lattice_candidates(identity=identity)
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=second,
    )

    with pytest.raises(
        EngagementCapExceededError,
        match="operation 4",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
                link,
                _candidate_circle_operation(identity, second),
            ),
        )


def test_replay_rejects_reused_neck_transition_before_fresh_state_mutation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject replay that reuses an already consumed passage transition.

    The second circle keeps its genuine geometry and traversal but repeats the
    first circle's cap decision.  Accepting it would let a recorded operation
    stream choose passage state instead of deriving state from prior accepted
    circles.  The entire stream must fail before authenticated entry depletion.
    """
    identity = _input_identity(
        neck_cap_degrees=(120.0, 120.0, 120.0),
    )
    first, second = _neck_lattice_candidates(identity=identity)
    mutated_second = CutFullCircleOperation.build(
        motion=second.motion,
        cut_z=identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
        neck_scope=second.neck_scope,
        effective_cap_decision=first.effective_cap_decision,
        traversal_decision=second.traversal_decision,
    )
    link = _link_to_candidate(
        identity,
        start=_phase_point(first),
        candidate=second,
    )

    def unexpected_deplete(
        stock: Stock2Area,
        depletion: object,
        *args: object,
    ) -> object:
        pytest.fail("neck replay mutated fresh stock")

    monkeypatch.setattr(Stock2Area, "deplete", unexpected_deplete)

    with pytest.raises(
        ReplayEffectiveCapError,
        match="fresh oriented-neck passage",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                _candidate_circle_operation(identity, first),
                link,
                mutated_second,
            ),
        )


def test_replay_rejects_neck_owner_absent_from_fresh_inventory() -> None:
    """Reject a structurally valid scope that names no rebuilt exact neck.

    The cap decision, circle, and traversal all come from a genuine first
    passage.  Replacing only the scope owner must fail before candidate matching;
    otherwise arbitrary operation bytes could select policy state without an
    authenticated MAT-neck owner.
    """
    identity = _input_identity()
    first, _ = _neck_lattice_candidates(identity=identity)
    scope = first.neck_scope
    assert isinstance(scope, OrientedNeckScope)
    foreign_scope = OrientedNeckScope.build(
        neck_owner_id=NeckOwnerId(b"\xff" * 32),
        orientation=scope.orientation,
    )
    mutated = CutFullCircleOperation.build(
        motion=first.motion,
        cut_z=identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
        neck_scope=foreign_scope,
        effective_cap_decision=first.effective_cap_decision,
        traversal_decision=first.traversal_decision,
    )

    with pytest.raises(
        ReplayEffectiveCapError,
        match="no owner in the fresh exact inventory",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                mutated,
            ),
        )


@pytest.mark.parametrize(
    "mutated_field",
    ("evidence", "width-class", "effective-cap"),
)
def test_replay_rejects_each_mutated_neck_cap_field(
    mutated_field: str,
) -> None:
    """Bind every policy-bearing field of the fresh neck-cap decision.

    Each operation retains genuine motion, scope, traversal, and passage state.
    A replay implementation that compared only the state transition could still
    accept foreign evidence, a broader width class, or altered cap bytes.
    """
    identity = _input_identity()
    first, _ = _neck_lattice_candidates(identity=identity)
    decision = first.effective_cap_decision
    assert isinstance(decision, NeckCapDecision)
    evidence_digest = decision.neck_evidence_digest
    width_class_id = decision.width_class_id
    effective_cap = EngagementCap.build(math.radians(90.0))
    if mutated_field == "evidence":
        evidence_digest = IdentityDigest(b"\xff" * 32)
    elif mutated_field == "width-class":
        width_class_id = WidthClassId.build(1)
    else:
        effective_cap = EngagementCap.build(math.radians(85.0))
    mutated_decision = NeckCapDecision.build(
        neck_evidence_digest=evidence_digest,
        width_class_id=width_class_id,
        passage_before=decision.passage_before,
        passage_after=decision.passage_after,
        user_cap=identity.user_cap,
        effective_cap=effective_cap,
    )
    mutated = CutFullCircleOperation.build(
        motion=first.motion,
        cut_z=identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
        neck_scope=first.neck_scope,
        effective_cap_decision=mutated_decision,
        traversal_decision=first.traversal_decision,
    )

    with pytest.raises(
        ReplayEffectiveCapError,
        match="fresh oriented-neck passage",
    ):
        _replay(
            identity,
            (
                identity.entry.approach,
                identity.entry.plunge,
                mutated,
            ),
        )
