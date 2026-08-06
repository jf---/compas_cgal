import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast

from compas.geometry import Polygon  # type: ignore[import-untyped]

from compas_cgal.adaptive.advancing_segment_trial import evaluate_advancing_segment_trial
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import TraversalCandidate
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.candidates import enumerate_zero_guide_link_candidates
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import CircleContainmentCertificate
from compas_cgal.adaptive.containment import CircleInEntryCertificate
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import InvalidReplayCertificateError
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.errors import ReplayCandidateError
from compas_cgal.adaptive.errors import ReplayContinuityError
from compas_cgal.adaptive.errors import ReplayCutDirectionError
from compas_cgal.adaptive.errors import ReplayEffectiveCapError
from compas_cgal.adaptive.errors import ReplayGrammarError
from compas_cgal.adaptive.errors import ReplayInputMismatchError
from compas_cgal.adaptive.errors import ReplayPairingError
from compas_cgal.adaptive.errors import ReplayTraversalError
from compas_cgal.adaptive.errors import ReplayZeroGuideCandidateError
from compas_cgal.adaptive.errors import TerminalNeckPassageError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MatZeroGuideRun
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import EffectiveCapDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckScope
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.replay_trace import FreshReplayTrace
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType


def _canonical_ring(
    polygon: Polygon,
    *,
    outer: bool,
) -> CanonicalRingV1:
    if type(polygon) is not Polygon:
        raise ReplayInputMismatchError("replay pocket geometry must use exact compas Polygon inputs.")
    if any(point.z != 0.0 for point in polygon.points):
        raise ReplayInputMismatchError("replay pocket geometry must lie exactly in world XY.")
    points = tuple(Point2[WorldXY].build(point.x, point.y) for point in polygon.points)
    try:
        if outer:
            return CanonicalRingV1.build_outer(points)
        return CanonicalRingV1.build_hole(points)
    except ValueError as error:
        raise ReplayInputMismatchError(str(error)) from None


def _canonical_design(
    pocket: Polygon,
    holes: Sequence[Polygon],
) -> tuple[
    CanonicalRingV1,
    tuple[CanonicalRingV1, ...],
    tuple[Polygon, ...],
]:
    if isinstance(holes, (str, bytes)):
        raise ReplayInputMismatchError("replay holes must be a finite polygon sequence.")
    try:
        hole_polygons = tuple(holes)
    except TypeError:
        raise ReplayInputMismatchError("replay holes must be a finite polygon sequence.") from None
    design_boundary = _canonical_ring(pocket, outer=True)
    canonical_holes = tuple(
        sorted(
            (_canonical_ring(hole, outer=False) for hole in hole_polygons),
            key=lambda ring: ring.canonical_bytes,
        )
    )
    if len(canonical_holes) != len({hole.canonical_bytes for hole in canonical_holes}):
        raise ReplayInputMismatchError("replay holes contain a canonical duplicate.")
    return design_boundary, canonical_holes, hole_polygons


def _require_bound_inputs(
    *,
    input_identity: InputIdentity,
    design_boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    entry: PreclearedEntry,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> None:
    if type(input_identity) is not InputIdentity:
        raise ReplayInputMismatchError("replay requires one exact InputIdentity.")
    for submitted, expected, name in (
        (design_boundary, input_identity.design_boundary, "design boundary"),
        (holes, input_identity.holes, "hole set"),
        (cut_plane, input_identity.cut_plane, "cut plane"),
        (tool_radius, input_identity.tool_radius, "tool radius"),
        (user_cap, input_identity.user_cap, "user cap"),
        (
            candidate_policy,
            input_identity.candidate_policy,
            "candidate policy",
        ),
        (neck_policy, input_identity.neck_policy, "neck policy"),
        (
            depletion_policy,
            input_identity.depletion_policy,
            "depletion policy",
        ),
        (
            traversal_policy,
            input_identity.traversal_policy,
            "traversal policy",
        ),
        (
            cut_direction_policy,
            input_identity.cut_direction_policy,
            "cut-direction policy",
        ),
    ):
        if submitted != expected:
            raise ReplayInputMismatchError(f"submitted {name} differs from the authenticated input root.")
    if type(entry) is not PreclearedEntry:
        raise ReplayInputMismatchError("replay entry must use exact PreclearedEntry.")
    if entry.canonical_bytes != input_identity.entry.canonical_bytes:
        raise ReplayInputMismatchError("submitted entry differs from the authenticated input root.")


def _rebuild_input_identity(
    *,
    input_identity: InputIdentity,
    design_boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    entry: PreclearedEntry,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> tuple[InputIdentity, ReachableDomain, PreclearedEntry]:
    reachable_domain = ReachableDomain.build(
        design_boundary=design_boundary,
        holes=holes,
        tool_radius=tool_radius,
    )
    qualified_bore = QualifiedBore.build(
        cut_plane=cut_plane,
        process_identity=entry.qualified_bore.process_identity,
        evidence_bytes=entry.qualified_bore.evidence_bytes,
    )
    fresh_entry = PreclearedEntry.build(
        reachable_domain=reachable_domain,
        center=entry.center,
        radius=entry.radius,
        tool_radius=tool_radius,
        cut_plane=cut_plane,
        qualified_bore=qualified_bore,
    )
    rebuilt = InputIdentity.build(
        design_boundary=design_boundary,
        holes=holes,
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        reachable_domain=reachable_domain,
        entry=fresh_entry,
        user_cap=user_cap,
        mat_sampling_policy=input_identity.mat_sampling_policy,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
    )
    if rebuilt.canonical_bytes != input_identity.canonical_bytes or rebuilt.digest != input_identity.digest:
        raise ReplayInputMismatchError("fresh exact-state rebuild does not reproduce InputIdentity.")
    return rebuilt, reachable_domain, fresh_entry


def _validate_grammar(
    operations: tuple[CanonicalOperation, ...],
    entry: PreclearedEntry,
    cut_plane: CutPlane,
) -> None:
    if type(operations) is not tuple or len(operations) < 3:
        raise ReplayGrammarError("operation grammar requires approach, plunge, and lateral motion.")
    if type(operations[0]) is not ApproachOperation or type(operations[1]) is not PlungeOperation:
        raise ReplayGrammarError("operation grammar requires approach, plunge, then lateral motion.")
    if operations[0] != entry.approach or operations[1] != entry.plunge:
        raise ReplayGrammarError("approach and plunge must exactly match the qualified entry.")
    lateral = operations[2:]
    for operation in lateral:
        if type(operation) is LinkSegmentOperation:
            if operation.cut_z != cut_plane.cut_z:
                raise ReplayContinuityError(
                    "every lateral operation must remain on the authenticated cut plane.",
                )
            continue
        if type(operation) is AdvanceSegmentOperation:
            if operation.cut_z != cut_plane.cut_z:
                raise ReplayContinuityError("every lateral operation must remain on the authenticated cut plane.")
            continue
        if type(operation) is CutFullCircleOperation:
            if operation.cut_z != cut_plane.cut_z:
                raise ReplayContinuityError("every lateral operation must remain on the authenticated cut plane.")
            continue
        raise ReplayGrammarError(
            "operation grammar permits only hold links, advancing segments, and full circles after plunge.",
        )
    if type(lateral[0]) is not CutFullCircleOperation:
        raise ReplayGrammarError(
            "first lateral operation must be the authenticated entry full circle.",
        )

    lateral_index = 1
    while lateral_index < len(lateral):
        operation = lateral[lateral_index]
        if type(operation) is AdvanceSegmentOperation:
            if lateral_index + 1 < len(lateral) and type(lateral[lateral_index + 1]) is CutFullCircleOperation:
                raise ReplayZeroGuideCandidateError(
                    "zero-guide advance cannot replace the hold link of a later circle.",
                )
            lateral_index += 1
            continue
        if type(operation) is CutFullCircleOperation:
            raise ReplayGrammarError(
                "every later circle must immediately follow one hold link.",
            )
        if type(operation) is LinkSegmentOperation:
            if lateral_index + 1 >= len(lateral):
                raise ReplayPairingError(
                    "link has no immediately following circle.",
                )
            following = lateral[lateral_index + 1]
            if type(following) is LinkSegmentOperation:
                raise ReplayPairingError(
                    "consecutive links do not have one immediately following circle.",
                )
            if type(following) is AdvanceSegmentOperation:
                raise ReplayZeroGuideCandidateError(
                    "zero-guide advance cannot consume a preceding hold link.",
                )
            if type(following) is not CutFullCircleOperation:
                raise ReplayGrammarError(
                    "hold link must immediately precede one full circle.",
                )
            lateral_index += 2
            continue
        raise ReplayGrammarError(
            "operation grammar contains an unknown lateral type.",
        )


def _phase_point(
    operation: CutFullCircleOperation,
) -> Point2[WorldXY]:
    return Point2[WorldXY].build(
        operation.motion.center.x + operation.motion.phase_vector.x,
        operation.motion.center.y + operation.motion.phase_vector.y,
    )


def _validate_continuity(
    operations: tuple[CanonicalOperation, ...],
    entry: PreclearedEntry,
) -> None:
    current = entry.center
    for operation in operations[2:]:
        if type(operation) is LinkSegmentOperation:
            if operation.motion.start != current:
                raise ReplayContinuityError(
                    "link start does not equal the preceding exact phase endpoint.",
                )
            current = operation.motion.end
            continue
        if type(operation) is AdvanceSegmentOperation:
            if operation.motion.start != current:
                raise ReplayContinuityError(
                    "advancing segment start does not equal the preceding exact phase endpoint.",
                )
            current = operation.motion.end
            continue
        if type(operation) is CutFullCircleOperation:
            phase = _phase_point(operation)
            if phase != current:
                raise ReplayContinuityError("circle phase does not equal the preceding exact endpoint.")
            current = phase
            continue
        raise ReplayGrammarError("operation grammar contains an unknown lateral type.")


def _validate_cut_direction(
    operations: tuple[CanonicalOperation, ...],
    policy: CutDirectionPolicy,
) -> None:
    for operation in operations[2:]:
        if type(operation) is not CutFullCircleOperation:
            continue
        expected = policy.circle_orientation(operation.material_side)
        if operation.motion.clockwise is not expected.clockwise:
            raise ReplayCutDirectionError("circle orientation contradicts its bound material side and cut-direction policy.")


def _operation_chain(
    input_digest: IdentityDigest,
    operations: tuple[CanonicalOperation, ...],
) -> tuple[bytes, tuple[IdentityDigest, ...]]:
    operation_records = tuple(operation.canonical_bytes for operation in operations)
    ordered_digest = hashlib.sha256(
        encode_tagged_union(
            b"replay-ordered-operations-v1",
            encode_sequence(operation_records),
        )
    ).digest()
    parent = bytes(input_digest)
    chain: list[IdentityDigest] = []
    for index, record in enumerate(operation_records):
        link = IdentityDigest(
            hashlib.sha256(
                encode_tagged_union(
                    b"replay-operation-index-v1",
                    encode_component_map(
                        {
                            b"index": encode_integer(index),
                            b"operation": record,
                            b"parent": encode_bytes(parent),
                        }
                    ),
                )
            ).digest()
        )
        chain.append(link)
        parent = bytes(link)
    return ordered_digest, tuple(chain)


def _rebuild_medial_axis(
    *,
    input_identity: InputIdentity,
    design_boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
    tool_radius: ToolRadius,
    reachable_domain: ReachableDomain,
) -> MedialAxis:
    sampling = input_identity.mat_sampling_policy
    axis = MedialAxis.build(
        design_boundary=design_boundary,
        holes=holes,
        tool_radius=tool_radius,
        station_spacing=sampling.station_spacing,
        max_sagitta=sampling.max_sagitta,
        max_refinement_depth=sampling.max_refinement_depth,
    )
    if axis.center_domain_digest != reachable_domain.certificate.center_domain_digest:
        raise ReplayInputMismatchError("fresh MAT center-domain identity contradicts the authenticated reachable domain.")
    return axis


def _edge_samples(
    axis: MedialAxis,
    edge_id: MatEdgeId,
) -> tuple[MatSample, ...]:
    samples = tuple(
        sorted(
            (sample for sample in axis.samples if sample.edge_id == edge_id),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    if not samples:
        raise ReplayTraversalError("fresh MAT edge has no proposal cursor sequence.")
    return samples


def _initial_cursor(
    *,
    samples: tuple[MatSample, ...],
    operation: CutFullCircleOperation | AdvanceSegmentOperation,
) -> MatSample:
    matches = tuple(sample for sample in samples if sample.cursor_identity == operation.traversal_decision.cursor_before)
    if len(matches) != 1:
        raise ReplayTraversalError("recorded cursor-before is not one unique fresh MAT sample.")
    return matches[0]


def _directed_limits(
    *,
    samples: tuple[MatSample, ...],
    cursor_before: MatSample | DerivedCandidateCursor,
    window: int,
) -> tuple[MatSample, ...]:
    """Enumerate exact native limits without inferring from circle orientation."""
    if type(cursor_before) is MatSample:
        ordinal = cursor_before.ordinal_on_edge
        increasing = tuple(sample for sample in samples if sample.ordinal_on_edge > ordinal)[:window]
        decreasing = tuple(sample for sample in reversed(samples) if sample.ordinal_on_edge < ordinal)[:window]
        return (*increasing, *decreasing)

    derived_cursor = cast(DerivedCandidateCursor, cursor_before)
    next_ordinal = derived_cursor.next_limit_ordinal
    if derived_cursor.ordinal_step == 1:
        return tuple(sample for sample in samples if sample.ordinal_on_edge >= next_ordinal)[:window]
    return tuple(sample for sample in reversed(samples) if sample.ordinal_on_edge <= next_ordinal)[:window]


def _fresh_neck_passages(
    inventory: NeckInventory,
) -> dict[OrientedNeckScope, NeckPassage]:
    passages: dict[OrientedNeckScope, NeckPassage] = {}
    for neck in inventory.necks:
        for passage in (neck.forward, neck.reverse):
            passages[passage.scope] = passage
    return passages


def _reconstruct_cap_decision(
    *,
    operation: CutFullCircleOperation | AdvanceSegmentOperation,
    user_cap: EngagementCap,
    neck_policy: NeckPolicy,
    passages: dict[OrientedNeckScope, NeckPassage],
) -> tuple[
    NeckScope,
    EffectiveCapDecision,
    EngagementCap,
    NeckPassage | None,
]:
    if type(operation.neck_scope) is NoNeckScope:
        expected_scope = NoNeckScope.build()
        expected_full_cap = FullCapDecision.build(
            user_cap=user_cap,
            effective_cap=user_cap,
        )
        if operation.effective_cap_decision != expected_full_cap:
            raise ReplayEffectiveCapError("recorded no-neck cap decision differs from the fresh full-cap policy result.")
        return expected_scope, expected_full_cap, user_cap, None

    if type(operation.neck_scope) is not OrientedNeckScope:
        raise ReplayEffectiveCapError("recorded neck scope is outside the closed replay grammar.")
    passage = passages.get(operation.neck_scope)
    if passage is None:
        raise ReplayEffectiveCapError("recorded oriented-neck scope has no owner in the fresh exact inventory.")
    try:
        expected_neck_cap = passage.propose_cap_decision(neck_policy)
    except TerminalNeckPassageError:
        raise ReplayEffectiveCapError("recorded oriented-neck passage is already terminal.") from None
    if operation.effective_cap_decision != expected_neck_cap:
        raise ReplayEffectiveCapError("recorded cap decision differs from the fresh oriented-neck passage.")
    effective_cap = neck_policy.effective_cap(
        passage.neck.width_class_id.value,
        passage.state,
    )
    return (
        passage.scope,
        expected_neck_cap,
        effective_cap,
        passage.advance(expected_neck_cap),
    )


def _candidate_effective_cap(
    *,
    candidate: TraversalCandidate,
    user_cap: EngagementCap,
    neck_policy: NeckPolicy,
) -> EngagementCap:
    decision = candidate.effective_cap_decision
    if type(decision) is FullCapDecision:
        return user_cap
    if type(decision) is NeckCapDecision:
        return neck_policy.effective_cap(
            decision.width_class_id.value,
            decision.passage_before,
        )
    raise ReplayEffectiveCapError("reconstructed candidate cap is outside the closed replay grammar.")


def _candidate_matches_operation(
    candidate: MiddleCurveCandidate,
    operation: CutFullCircleOperation,
) -> bool:
    return (
        candidate.motion == operation.motion
        and candidate.neck_scope == operation.neck_scope
        and candidate.effective_cap_decision == operation.effective_cap_decision
        and candidate.traversal_decision == operation.traversal_decision
    )


def _zero_guide_candidate_matches_operation(
    candidate: ZeroGuideLinkCandidate,
    operation: AdvanceSegmentOperation,
) -> bool:
    expected = AdvanceSegmentOperation.build(
        motion=ExactSegmentMotion.build(
            operation.motion.start,
            candidate.target,
        ),
        cut_z=operation.cut_z,
        neck_scope=candidate.neck_scope,
        effective_cap_decision=candidate.effective_cap_decision,
        traversal_decision=candidate.traversal_decision,
    )
    return expected.canonical_bytes == operation.canonical_bytes


def _fresh_zero_guide_run(
    *,
    axis: MedialAxis,
    edge_id: MatEdgeId,
) -> MatZeroGuideRun:
    inventory_matches = tuple(run for run in axis.zero_guide_inventory.runs if bytes(run.edge_id) == bytes(edge_id))
    native_matches = tuple(record for record in axis.native_owner.zero_guide_records if record[0] == bytes(edge_id))
    if len(inventory_matches) != 1 or len(native_matches) != 1:
        raise ReplayZeroGuideCandidateError(
            "zero-guide replay requires one fresh typed and native proof record.",
        )
    run = inventory_matches[0]
    if run.mat_certificate_digest != hashlib.sha256(axis.mat_certificate).digest() or native_matches[0] != (bytes(run.edge_id), run.native_certificate):
        raise ReplayZeroGuideCandidateError(
            "zero-guide replay proof bytes contradict the fresh MAT owner.",
        )
    return run


def _match_circle_candidate(
    *,
    axis: MedialAxis,
    operation: CutFullCircleOperation,
    user_cap: EngagementCap,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
    passages: dict[OrientedNeckScope, NeckPassage],
    current_by_edge: dict[
        MatEdgeId,
        MatSample | DerivedCandidateCursor | None,
    ],
) -> MiddleCurveCandidate:
    edge_id = MatEdgeId(bytes(operation.traversal_decision.edge_id))
    edge = axis.edge_by_id.get(edge_id)
    if edge is None:
        raise ReplayTraversalError("recorded traversal edge is absent from the fresh MAT graph.")
    samples = _edge_samples(axis, edge_id)
    if edge_id in current_by_edge:
        cursor_before = current_by_edge[edge_id]
        if cursor_before is None:
            raise ReplayTraversalError("recorded traversal advances an already terminal MAT cursor.")
        if cursor_before.cursor_identity != operation.traversal_decision.cursor_before:
            raise ReplayTraversalError("recorded cursor-before does not continue the fresh edge state.")
    else:
        cursor_before = _initial_cursor(
            samples=samples,
            operation=operation,
        )

    expected_scope, expected_cap, _, next_passage = _reconstruct_cap_decision(
        operation=operation,
        user_cap=user_cap,
        neck_policy=neck_policy,
        passages=passages,
    )
    orientation = cut_direction_policy.circle_orientation(operation.material_side)
    limits = _directed_limits(
        samples=samples,
        cursor_before=cursor_before,
        window=traversal_policy.forward_window,
    )
    if not limits:
        raise ReplayTraversalError("recorded cursor has no fresh native limit inside its directed window.")

    matches: dict[
        bytes,
        tuple[MiddleCurveSpan, MiddleCurveCandidate],
    ] = {}
    for limit in limits:
        span = MiddleCurveSpan.build(
            axis=axis,
            cursor_before=cursor_before,
            cursor_limit=limit,
        )
        terminal_sample = samples[-1] if span.ordinal_step == 1 else samples[0]
        terminal_limit = limit.ordinal_on_edge == terminal_sample.ordinal_on_edge
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=candidate_policy,
            circle_orientation=orientation,
            neck_scope=expected_scope,
            effective_cap_decision=expected_cap,
            makes_cursor_terminal_at_limit=terminal_limit,
        ):
            if _candidate_matches_operation(candidate, operation):
                matches[bytes(candidate.identity)] = (
                    span,
                    candidate,
                )
    if len(matches) != 1:
        raise ReplayCandidateError("recorded circle does not identify one unique finite-lattice candidate.")
    span, selected = next(iter(matches.values()))
    if selected.traversal_decision.makes_cursor_terminal:
        current_by_edge[edge_id] = None
    elif selected.spatial_progress == span.reported_length:
        current_by_edge[edge_id] = span.cursor_limit
    else:
        current_by_edge[edge_id] = DerivedCandidateCursor.build(
            span=span,
            candidate=selected,
        )
    if next_passage is not None:
        passages[next_passage.scope] = next_passage
    return selected


def _match_zero_guide_candidate(
    *,
    axis: MedialAxis,
    operation: AdvanceSegmentOperation,
    user_cap: EngagementCap,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    traversal_policy: TraversalPolicy,
    passages: dict[OrientedNeckScope, NeckPassage],
    current_by_edge: dict[
        MatEdgeId,
        MatSample | DerivedCandidateCursor | None,
    ],
) -> ZeroGuideLinkCandidate:
    """Reconstruct one advancing segment from the fresh zero-guide proof.

    Args:
        axis: Independently rebuilt exact MAT owner.
        operation: Recorded advancing-segment operation.
        user_cap: Authenticated user engagement limit.
        candidate_policy: Authenticated finite-lattice policy.
        neck_policy: Authenticated causal neck policy.
        traversal_policy: Authenticated directed-window policy.
        passages: Fresh mutable replay view of oriented passage state.
        current_by_edge: Fresh replay cursor state by exact MAT edge.

    Returns:
        The unique fresh spatial candidate owning `operation`.

    Raises:
        ReplayZeroGuideCandidateError: If proof ownership, cursor lineage, or
            unique finite reconstruction fails.
    """
    if type(operation) is not AdvanceSegmentOperation:
        raise ReplayZeroGuideCandidateError(
            "zero-guide reconstruction requires one exact advancing segment.",
        )
    edge_id = MatEdgeId(bytes(operation.traversal_decision.edge_id))
    edge = axis.edge_by_id.get(edge_id)
    if edge is None:
        raise ReplayZeroGuideCandidateError(
            "zero-guide traversal edge is absent from the fresh MAT graph.",
        )
    owned_run = _fresh_zero_guide_run(
        axis=axis,
        edge_id=edge_id,
    )
    samples = _edge_samples(axis, edge_id)
    if edge_id in current_by_edge:
        cursor_before = current_by_edge[edge_id]
        if cursor_before is None:
            raise ReplayZeroGuideCandidateError(
                "zero-guide operation advances an already terminal fresh cursor.",
            )
        if cursor_before.cursor_identity != operation.traversal_decision.cursor_before:
            raise ReplayZeroGuideCandidateError(
                "zero-guide cursor-before does not continue fresh edge state.",
            )
    else:
        try:
            cursor_before = _initial_cursor(
                samples=samples,
                operation=operation,
            )
        except ReplayTraversalError as error:
            raise ReplayZeroGuideCandidateError(
                f"zero-guide initial cursor is not fresh ({error}).",
            ) from None

    expected_scope, expected_cap, _, next_passage = _reconstruct_cap_decision(
        operation=operation,
        user_cap=user_cap,
        neck_policy=neck_policy,
        passages=passages,
    )
    limits = _directed_limits(
        samples=samples,
        cursor_before=cursor_before,
        window=traversal_policy.forward_window,
    )
    if not limits:
        raise ReplayZeroGuideCandidateError(
            "zero-guide cursor has no fresh native limit in its directed window.",
        )

    matches: dict[
        bytes,
        tuple[MiddleCurveSpan, ZeroGuideLinkCandidate],
    ] = {}
    for limit in limits:
        span = MiddleCurveSpan.build(
            axis=axis,
            cursor_before=cursor_before,
            cursor_limit=limit,
        )
        terminal_sample = samples[-1] if span.ordinal_step == 1 else samples[0]
        terminal_limit = limit.ordinal_on_edge == terminal_sample.ordinal_on_edge
        candidates = enumerate_zero_guide_link_candidates(
            span=span,
            policy=candidate_policy,
            neck_scope=expected_scope,
            effective_cap_decision=expected_cap,
            makes_cursor_terminal_at_limit=terminal_limit,
        )
        for candidate in candidates:
            if (
                type(candidate) is not ZeroGuideLinkCandidate
                or candidate.zero_guide_run != owned_run
                or candidate.zero_guide_run.native_certificate != owned_run.native_certificate
            ):
                raise ReplayZeroGuideCandidateError(
                    "zero-guide candidate contradicts fresh native proof bytes.",
                )
            if _zero_guide_candidate_matches_operation(
                candidate,
                operation,
            ):
                matches[bytes(candidate.identity)] = (
                    span,
                    candidate,
                )
    if len(matches) != 1:
        raise ReplayZeroGuideCandidateError(
            "recorded zero-guide segment does not identify one unique fresh candidate.",
        )
    span, selected = next(iter(matches.values()))
    if selected.traversal_decision.makes_cursor_terminal:
        current_by_edge[edge_id] = None
    elif selected.spatial_progress == span.reported_length:
        current_by_edge[edge_id] = span.cursor_limit
    else:
        current_by_edge[edge_id] = DerivedCandidateCursor.build(
            span=span,
            candidate=selected,
        )
    if next_passage is not None:
        passages[next_passage.scope] = next_passage
    return selected


def _replay_candidate_stream(
    *,
    axis: MedialAxis,
    operations: tuple[CanonicalOperation, ...],
    user_cap: EngagementCap,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> tuple[
    tuple[TraversalCandidate, ...],
    dict[MatEdgeId, MatSample | DerivedCandidateCursor | None],
]:
    current_by_edge: dict[
        MatEdgeId,
        MatSample | DerivedCandidateCursor | None,
    ] = {}
    passages = _fresh_neck_passages(
        NeckInventory.build(
            axis=axis,
            policy=neck_policy,
        )
    )
    candidates: list[TraversalCandidate] = []
    for operation in operations[2:]:
        if type(operation) is LinkSegmentOperation:
            edge_id = MatEdgeId(
                bytes(operation.traversal_decision.edge_id),
            )
            if edge_id in axis.zero_guide_run_by_edge_id:
                raise ReplayZeroGuideCandidateError(
                    "zero-guide advancement cannot be relabelled as a hold link.",
                )
            continue
        if type(operation) is CutFullCircleOperation:
            candidates.append(
                _match_circle_candidate(
                    axis=axis,
                    operation=operation,
                    user_cap=user_cap,
                    candidate_policy=candidate_policy,
                    neck_policy=neck_policy,
                    traversal_policy=traversal_policy,
                    cut_direction_policy=cut_direction_policy,
                    passages=passages,
                    current_by_edge=current_by_edge,
                )
            )
            continue
        if type(operation) is AdvanceSegmentOperation:
            candidates.append(
                _match_zero_guide_candidate(
                    axis=axis,
                    operation=operation,
                    user_cap=user_cap,
                    candidate_policy=candidate_policy,
                    neck_policy=neck_policy,
                    traversal_policy=traversal_policy,
                    passages=passages,
                    current_by_edge=current_by_edge,
                )
            )
    return tuple(candidates), current_by_edge


def _pair_lateral_operations(
    *,
    operations: tuple[CanonicalOperation, ...],
    candidates: tuple[TraversalCandidate, ...],
) -> dict[int, TraversalCandidate]:
    paired: dict[int, TraversalCandidate] = {}
    candidate_index = 0
    pending_link: tuple[int, LinkSegmentOperation] | None = None
    for operation_index, operation in enumerate(operations[2:], start=2):
        if type(operation) is LinkSegmentOperation:
            if pending_link is not None:
                raise ReplayPairingError("consecutive links do not have one immediately following circle.")
            pending_link = operation_index, operation
            continue
        if type(operation) is AdvanceSegmentOperation:
            if pending_link is not None:
                raise ReplayZeroGuideCandidateError(
                    "zero-guide advance cannot consume a preceding hold link.",
                )
            if candidate_index >= len(candidates):
                raise ReplayZeroGuideCandidateError(
                    "fresh candidate lineage omits an advancing segment.",
                )
            candidate = candidates[candidate_index]
            candidate_index += 1
            if type(candidate) is not ZeroGuideLinkCandidate or not _zero_guide_candidate_matches_operation(
                candidate,
                operation,
            ):
                raise ReplayZeroGuideCandidateError(
                    "fresh zero-guide lineage diverged before state replay.",
                )
            paired[operation_index] = candidate
            continue
        if type(operation) is not CutFullCircleOperation:
            raise ReplayGrammarError("operation grammar contains an unknown lateral type.")
        if candidate_index >= len(candidates):
            raise ReplayCandidateError("fresh candidate lineage is shorter than the recorded circle stream.")
        candidate = candidates[candidate_index]
        candidate_index += 1
        if type(candidate) is not MiddleCurveCandidate or not _candidate_matches_operation(
            candidate,
            operation,
        ):
            raise ReplayCandidateError("fresh candidate lineage diverged before state replay.")
        paired[operation_index] = candidate
        if pending_link is None:
            continue

        link_index, link = pending_link
        traversal = candidate.traversal_decision
        expected_hold = HoldTraversalDecision.build(
            component_id=traversal.component_id,
            edge_id=traversal.edge_id,
            branch_id=traversal.branch_id,
            cursor=traversal.cursor_before,
        )
        if link.traversal_decision != expected_hold:
            raise ReplayPairingError("link hold traversal does not belong to its immediately following circle.")
        if link.neck_scope != candidate.neck_scope or link.effective_cap_decision != candidate.effective_cap_decision:
            raise ReplayPairingError("link scope and effective cap do not belong to their immediately following circle.")
        if link.motion.end != _phase_point(operation):
            raise ReplayPairingError("link endpoint does not equal its immediately following circle phase.")
        paired[link_index] = candidate
        pending_link = None

    if pending_link is not None:
        raise ReplayPairingError("link has no immediately following circle.")
    if candidate_index != len(candidates):
        raise ReplayCandidateError(
            "fresh candidate lineage is longer than the recorded lateral stream.",
        )
    return paired


def _replay_fresh_state(
    *,
    pocket: Polygon,
    holes: tuple[Polygon, ...],
    reachable_domain: ReachableDomain,
    entry: PreclearedEntry,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    operations: tuple[CanonicalOperation, ...],
    paired_candidates: dict[int, TraversalCandidate],
    axis: MedialAxis,
    current_by_edge: dict[
        MatEdgeId,
        MatSample | DerivedCandidateCursor | None,
    ],
) -> FreshReplayTrace:
    stock = Stock2Area.build(
        Stock(
            pocket,
            list(holes),
        )
    )
    pristine_state = MotionCertifier.build(
        stock=stock,
        tool_radius=tool_radius,
    )
    entry_depletion_witness = stock.deplete(entry)
    current_stock_state = MotionCertifier.build(
        stock=stock,
        tool_radius=tool_radius,
    )
    post_entry_stock_boundary_digest = current_stock_state.canonical_boundary_digest
    post_entry_stock_lineage_digest = current_stock_state.stock_lineage_digest
    coverage = CoverageLedger.build(
        reachable_domain=reachable_domain,
        precleared_center=entry.center,
        precleared_radius=entry.radius,
    )
    initial_coverage_certificate = coverage.certificate
    containment = GougeContainment.build(reachable_domain)
    lateral_witnesses: list[ReplayLateralWitness] = []
    first_circle_entry_certificate: CircleInEntryCertificate | None = None
    first_circle = True
    containment_certificate: SegmentContainmentCertificate | CircleContainmentCertificate
    for operation_index, operation in enumerate(operations):
        if type(operation) is LinkSegmentOperation:
            if operation_index not in paired_candidates:
                raise ReplayPairingError("fresh link has no validated circle pairing.")
            candidate = paired_candidates[operation_index]
            if type(candidate) is not MiddleCurveCandidate:
                raise ReplayPairingError(
                    "fresh hold link is paired with a foreign candidate variant.",
                )
            effective_cap = _candidate_effective_cap(
                candidate=candidate,
                user_cap=user_cap,
                neck_policy=neck_policy,
            )
            containment_certificate = containment.certify_segment(
                operation.motion,
                tool_radius,
            )
            motion_witness = current_stock_state.certify(
                operation_index=operation_index,
                operation_kind=OperationType.LINK,
                motion=operation.motion,
                user_cap=user_cap,
                effective_cap=effective_cap,
            )
            depletion_witness = stock.deplete(
                operation.motion,
                tool_radius,
                depletion_policy,
            )
            sweep_witness = coverage.add_sweep(
                operation.motion,
                tool_radius,
            )
            lateral_witnesses.append(
                ReplayLateralWitness(
                    operation_index=operation_index,
                    operation=operation,
                    effective_cap_decision=candidate.effective_cap_decision,
                    stock_boundary_digest=current_stock_state.canonical_boundary_digest,
                    containment_certificate=containment_certificate,
                    motion_witness=motion_witness,
                    depletion_witness=depletion_witness,
                    sweep_witness=sweep_witness,
                )
            )
            current_stock_state = MotionCertifier.build(
                stock=stock,
                tool_radius=tool_radius,
            )
            continue
        if type(operation) is AdvanceSegmentOperation:
            if operation_index not in paired_candidates:
                raise ReplayZeroGuideCandidateError(
                    "fresh advancing segment has no reconstructed zero-guide candidate.",
                )
            candidate = paired_candidates[operation_index]
            if type(candidate) is not ZeroGuideLinkCandidate:
                raise ReplayZeroGuideCandidateError(
                    "fresh advancing segment is paired with a foreign candidate variant.",
                )
            effective_cap = _candidate_effective_cap(
                candidate=candidate,
                user_cap=user_cap,
                neck_policy=neck_policy,
            )
            segment_witness = evaluate_advancing_segment_trial(
                containment_authority=containment,
                stock=stock,
                coverage=coverage,
                operation_index=operation_index,
                operation=operation,
                tool_radius=tool_radius,
                user_cap=user_cap,
                effective_cap=effective_cap,
                depletion_policy=depletion_policy,
            )
            lateral_witnesses.append(segment_witness)
            current_stock_state = MotionCertifier.build(
                stock=stock,
                tool_radius=tool_radius,
            )
            continue
        if type(operation) is not CutFullCircleOperation:
            continue
        if operation_index not in paired_candidates:
            raise ReplayCandidateError("fresh circle has no reconstructed candidate.")
        candidate = paired_candidates[operation_index]
        if type(candidate) is not MiddleCurveCandidate:
            raise ReplayCandidateError(
                "fresh circle is paired with a foreign candidate variant.",
            )
        effective_cap = _candidate_effective_cap(
            candidate=candidate,
            user_cap=user_cap,
            neck_policy=neck_policy,
        )
        if first_circle:
            first_circle_entry_certificate = entry.certify_first_circle(operation.motion)
            first_circle = False
        containment_certificate = containment.certify_full_circle(
            operation.motion,
            tool_radius,
        )
        motion_witness = current_stock_state.certify(
            operation_index=operation_index,
            operation_kind=OperationType.CUT,
            motion=operation.motion,
            user_cap=user_cap,
            effective_cap=effective_cap,
        )
        depletion_witness = stock.deplete(
            operation.motion,
            tool_radius,
            depletion_policy,
        )
        sweep_witness = coverage.add_sweep(
            operation.motion,
            tool_radius,
        )
        lateral_witnesses.append(
            ReplayLateralWitness(
                operation_index=operation_index,
                operation=operation,
                effective_cap_decision=candidate.effective_cap_decision,
                stock_boundary_digest=current_stock_state.canonical_boundary_digest,
                containment_certificate=containment_certificate,
                motion_witness=motion_witness,
                depletion_witness=depletion_witness,
                sweep_witness=sweep_witness,
            )
        )
        current_stock_state = MotionCertifier.build(
            stock=stock,
            tool_radius=tool_radius,
        )
    if first_circle_entry_certificate is None:
        raise InvalidReplayTraceError("fresh replay did not retain its required first-circle entry proof.")
    trace = FreshReplayTrace(
        pristine_stock_boundary_digest=pristine_state.canonical_boundary_digest,
        post_entry_stock_boundary_digest=post_entry_stock_boundary_digest,
        post_entry_stock_lineage_digest=post_entry_stock_lineage_digest,
        entry_depletion_witness=entry_depletion_witness,
        first_circle_entry_certificate=first_circle_entry_certificate,
        initial_coverage_certificate=initial_coverage_certificate,
        lateral_witnesses=tuple(lateral_witnesses),
        terminal_stock_boundary_digest=current_stock_state.canonical_boundary_digest,
        terminal_stock_lineage_digest=current_stock_state.stock_lineage_digest,
        terminal_coverage_certificate=coverage.certificate,
    )
    _require_terminal_replay(
        axis=axis,
        current_by_edge=current_by_edge,
        coverage=coverage,
        trace=trace,
    )
    return trace


def _require_terminal_replay(
    *,
    axis: MedialAxis,
    current_by_edge: dict[
        MatEdgeId,
        MatSample | DerivedCandidateCursor | None,
    ],
    coverage: CoverageLedger,
    trace: FreshReplayTrace,
) -> None:
    FreshReplayTrace.validate(trace)
    if set(current_by_edge) != set(axis.edge_by_id) or any(cursor is not None for cursor in current_by_edge.values()):
        raise ReplayTraversalError("fresh MAT traversal remains nonterminal.")
    coverage.require_complete()
    if coverage.certificate.canonical_bytes != trace.terminal_coverage_certificate.canonical_bytes:
        raise InvalidReplayTraceError("terminal replay trace does not match the live fresh coverage owner.")


@dataclass(frozen=True)
class ReplayCertificate:
    input_digest: IdentityDigest
    ordered_operation_digest: bytes
    operation_index_chain: tuple[IdentityDigest, ...]

    def __post_init__(self) -> None:
        if type(self.input_digest) is not bytes or len(self.input_digest) != 32:
            raise InvalidReplayCertificateError("replay certificate input identity must be one SHA-256 digest.")
        if type(self.ordered_operation_digest) is not bytes or len(self.ordered_operation_digest) != 32:
            raise InvalidReplayCertificateError("replay certificate operation digest must be one SHA-256 digest.")
        if (
            type(self.operation_index_chain) is not tuple
            or not self.operation_index_chain
            or any(type(digest) is not bytes or len(digest) != 32 for digest in self.operation_index_chain)
        ):
            raise InvalidReplayCertificateError("replay certificate operation chain must contain exact SHA-256 digests.")


def replay_certificate(
    *,
    input_identity: InputIdentity,
    pocket: Polygon,
    holes: Sequence[Polygon],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    entry: PreclearedEntry,
    operations: tuple[CanonicalOperation, ...],
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> ReplayCertificate:
    design_boundary, canonical_holes, hole_polygons = _canonical_design(
        pocket,
        holes,
    )
    _require_bound_inputs(
        input_identity=input_identity,
        design_boundary=design_boundary,
        holes=canonical_holes,
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        user_cap=user_cap,
        entry=entry,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
    )
    _validate_grammar(operations, entry, cut_plane)
    _validate_continuity(operations, entry)
    _validate_cut_direction(operations, cut_direction_policy)
    rebuilt, reachable_domain, fresh_entry = _rebuild_input_identity(
        input_identity=input_identity,
        design_boundary=design_boundary,
        holes=canonical_holes,
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        user_cap=user_cap,
        entry=entry,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
    )
    if fresh_entry.canonical_bytes != entry.canonical_bytes:
        raise ReplayInputMismatchError("fresh qualified entry does not reproduce the submitted entry.")
    _operation_chain(rebuilt.digest, operations)
    axis = _rebuild_medial_axis(
        input_identity=rebuilt,
        design_boundary=design_boundary,
        holes=canonical_holes,
        tool_radius=tool_radius,
        reachable_domain=reachable_domain,
    )
    candidates, current_by_edge = _replay_candidate_stream(
        axis=axis,
        operations=operations,
        user_cap=user_cap,
        candidate_policy=candidate_policy,
        neck_policy=neck_policy,
        traversal_policy=traversal_policy,
        cut_direction_policy=cut_direction_policy,
    )
    paired_candidates = _pair_lateral_operations(
        operations=operations,
        candidates=candidates,
    )
    _replay_fresh_state(
        pocket=pocket,
        holes=hole_polygons,
        reachable_domain=reachable_domain,
        entry=fresh_entry,
        tool_radius=tool_radius,
        user_cap=user_cap,
        neck_policy=neck_policy,
        depletion_policy=depletion_policy,
        operations=operations,
        paired_candidates=paired_candidates,
        axis=axis,
        current_by_edge=current_by_edge,
    )
    raise InvalidReplayCertificateError("complete fresh replay must bind the immutable certificate before return.")
