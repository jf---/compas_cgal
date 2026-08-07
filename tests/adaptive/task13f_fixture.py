"""Exact radius-1 fixture for Task 13F continuation contracts."""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Self

from compas_cgal.adaptive.bootstrap import InitialCandidateEvaluator
from compas_cgal.adaptive.bootstrap import InitialCandidateTransaction
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generator import _active_forward_limits
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import RouteNodeId
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import MatSamplingPolicy
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.traversal import MatTraversalState
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

TASK13F_OUTER = (
    (0.0, 0.0),
    (6.0, 0.0),
    (6.0, 2.0),
    (2.0, 2.0),
    (2.0, 6.0),
    (0.0, 6.0),
)
TASK13F_LAUNCH_PROGRESS = Fraction(
    6_234_160_216_152_465,
    18_014_398_509_481_984,
)
TASK13F_LAUNCH_GUIDE_RADIUS = Fraction(1, 32)
TASK13F_ENTRY_RADIUS = Fraction(17, 16)


def _ring() -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y in TASK13F_OUTER),
    )


def _phase_point(
    candidate: MiddleCurveCandidate,
) -> Point2[WorldXY]:
    motion = candidate.motion
    return Point2[WorldXY].build(
        motion.center.x + motion.phase_vector.x,
        motion.center.y + motion.phase_vector.y,
    )


@dataclass(frozen=True)
class Task13FFixture:
    """Authenticated launch child and continuation authority."""

    identity: InputIdentity
    seeded_traversal: MatTraversalState
    initial_evaluator: InitialCandidateEvaluator
    launch_transaction: InitialCandidateTransaction
    physical: GenerationState
    traversal: MatTraversalState
    evaluator: CandidateEvaluator

    @classmethod
    def build(cls) -> Self:
        ring = _ring()
        tool_radius = ToolRadius.build(1.0)
        user_cap = EngagementCap.build(math.pi)
        sampling_policy = MatSamplingPolicy.build(
            station_spacing=Spacing.build(2.0),
            max_sagitta=ChordBound.build(0.02),
            max_refinement_depth=32,
        )
        candidate_policy = CandidatePolicy.build(
            spatial_resolution=Spacing.build(0.125),
            spatial_refinement_levels=1,
            radius_resolution=Spacing.build(1.0),
            radius_refinement_levels=1,
            phase_count=2,
            minimum_guide_radius=GuideRadius.build(
                1.0 / 128.0,
            ),
            minimum_progress=Spacing.build(0.125),
        )
        neck_policy = NeckPolicy.build(
            user_cap=user_cap,
            squared_width_boundaries=(),
            effective_caps={(0, state): user_cap for state in ACTIVE_PASSAGE_STATES},
        )
        depletion_policy = DepletionPolicy.build(
            chord_bound=ChordBound.build(0.125),
            center_count_limit=4096,
        )
        traversal_policy = TraversalPolicy.build(
            forward_window=4,
        )
        cut_direction_policy = CutDirectionPolicy.build(
            CutIntent.CLIMB,
        )
        material_side = MaterialSide.OUTSIDE
        axis = MedialAxis.build(
            design_boundary=ring,
            holes=(),
            tool_radius=tool_radius,
            station_spacing=sampling_policy.station_spacing,
            max_sagitta=sampling_policy.max_sagitta,
            max_refinement_depth=(sampling_policy.max_refinement_depth),
        )
        inventory = NeckInventory.build(
            axis=axis,
            policy=neck_policy,
        )
        endpoint = Point2[WorldXY].build(2.0, 1.0)
        entry_edge = next(edge for edge in axis.edges if endpoint in (edge.source.point, edge.target.point))
        entry_node = entry_edge.target if entry_edge.source.point == endpoint else entry_edge.source
        seeded_traversal = MatTraversalState.seed(
            axis=axis,
            inventory=inventory,
            policy=traversal_policy,
            entry_edge_id=entry_edge.identity,
            entry_node_id=entry_node.identity,
        )
        active = seeded_traversal.active_cursor
        limit = _active_forward_limits(
            seeded_traversal,
        )[0]
        span = MiddleCurveSpan.build(
            axis=axis,
            cursor_before=active.cursor,
            cursor_limit=limit,
        )
        full_cap = FullCapDecision.build(
            user_cap=user_cap,
            effective_cap=user_cap,
        )
        launch = next(
            candidate
            for candidate in enumerate_middle_curve_candidates(
                span=span,
                policy=candidate_policy,
                circle_orientation=(
                    cut_direction_policy.circle_orientation(
                        material_side,
                    )
                ),
                neck_scope=seeded_traversal.neck_scope,
                effective_cap_decision=full_cap,
                makes_cursor_terminal_at_limit=(limit == active.terminal_cursor),
            )
            if candidate.spatial_progress == TASK13F_LAUNCH_PROGRESS
            and candidate.guide_radius == TASK13F_LAUNCH_GUIDE_RADIUS
            and candidate.phase_index == 1
            and candidate.proposal.generator_site.kind == "point"
        )
        reachable_domain = ReachableDomain.build(
            design_boundary=ring,
            holes=(),
            tool_radius=tool_radius,
        )
        cut_plane = CutPlane.build(
            CutZ.build(-2.0),
            ClearanceZ.build(5.0),
        )
        entry = PreclearedEntry.build(
            reachable_domain=reachable_domain,
            center=_phase_point(launch),
            radius=EntryRadius.build(
                float(TASK13F_ENTRY_RADIUS),
            ),
            tool_radius=tool_radius,
            cut_plane=cut_plane,
            qualified_bore=QualifiedBore.build(
                cut_plane=cut_plane,
                process_identity=BoreProcessIdentity(
                    b"task-13f-probe",
                ),
                evidence_bytes=(b"task-13f-probe-evidence"),
            ),
        )
        identity = InputIdentity.build(
            design_boundary=ring,
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
        initial_evaluator = InitialCandidateEvaluator.build(
            input_identity=identity,
            material_side=material_side,
        )
        launch_transaction = initial_evaluator.evaluate(
            seeded_traversal,
            launch,
        )
        physical, traversal = initial_evaluator.commit(
            seeded_traversal,
            launch_transaction,
        )
        evaluator = CandidateEvaluator.build(
            reachable_domain=identity.reachable_domain,
            tool_radius=identity.tool_radius,
            user_cap=identity.user_cap,
            candidate_policy=identity.candidate_policy,
            neck_policy=identity.neck_policy,
            depletion_policy=identity.depletion_policy,
            cut_direction_policy=(identity.cut_direction_policy),
            cut_z=identity.cut_plane.cut_z,
            material_side=material_side,
        )
        return cls(
            identity=identity,
            seeded_traversal=seeded_traversal,
            initial_evaluator=initial_evaluator,
            launch_transaction=launch_transaction,
            physical=physical,
            traversal=traversal,
            evaluator=evaluator,
        )


def task13f_route_one_terminal(
    fixture: Task13FFixture,
) -> tuple[
    GenerationState,
    MatTraversalState,
    tuple[TraversalCommit, TraversalCommit],
]:
    """Advance both established Task 13F routes to the retrace boundary.

    Args:
        fixture: Authenticated launch child and continuation authority.

    Returns:
        Physical child, route-one terminal traversal, and both commits.
    """
    physical, traversal, route_zero = advance_active_candidate_family(
        evaluator=fixture.evaluator,
        physical=fixture.physical,
        traversal=fixture.traversal,
    )
    assert traversal.active_cursor.terminal
    physical, traversal, route_one = advance_active_candidate_family(
        evaluator=fixture.evaluator,
        physical=physical,
        traversal=traversal.activate_next(),
    )
    assert traversal.active_route_index == 1
    assert traversal.active_cursor.terminal
    return physical, traversal, (route_zero, route_one)


def task13f_retrace_decision(
    *,
    physical: GenerationState,
    terminal: MatTraversalState,
    source_commit: TraversalCommit,
) -> RouteRetraceDecision:
    """Seal fixture preimages for physical-unit tests before generator wiring."""
    activated = terminal.activate_next()
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    return RouteRetraceDecision.build(
        completed_route_index=terminal.active_route_index,
        activated_route_index=activated.active_route_index,
        completed_exit_node_id=RouteNodeId(
            bytes(terminal.active_cursor.route_step.exit_node_id),
        ),
        activated_entry_node_id=RouteNodeId(
            bytes(activated.active_cursor.route_step.entry_node_id),
        ),
        terminal_traversal_digest=terminal.digest,
        activated_traversal_digest=activated.digest,
        source_commit_digest=source_commit.digest,
        source_transaction_digest=source_commit.transaction.digest,
        source_operation_index=len(physical.operations) - 1,
        source_operation_digest=IdentityDigest(
            hashlib.sha256(source.canonical_bytes).digest(),
        ),
    )
