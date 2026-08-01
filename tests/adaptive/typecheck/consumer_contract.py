from collections.abc import Mapping
from typing import Literal
from typing import assert_never
from typing import assert_type

import numpy as np

from compas_cgal import _medial_axis_2
from compas_cgal._medial_axis_2 import MedialAxisResult
from compas_cgal.adaptive.bootstrap import InitialCandidateEvaluator
from compas_cgal.adaptive.bootstrap import InitialCandidateTransaction
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import TraversalCandidate
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.candidates import enumerate_zero_guide_link_candidates
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.containment import CircleContainmentCertificate
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatZeroGuideInventory
from compas_cgal.adaptive.medial_axis import MatZeroGuideRun
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckScope
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.operation import TraversalDecision
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MatSamplingPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.traversal import DirectedEdgeCursor
from compas_cgal.adaptive.traversal import MatTraversalState
from compas_cgal.adaptive.traversal_graph import DirectedRouteStep
from compas_cgal.adaptive.traversal_graph import TraversalGraph
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock


class MachineXY:
    """Distinct machine-coordinate frame used only by this compile contract."""


mat_result = _medial_axis_2.segment_site_medial_axis(
    np.asarray(
        (
            (0.0, 0.0, 0.0),
            (6.0, 0.0, 0.0),
            (6.0, 2.0, 0.0),
            (2.0, 2.0, 0.0),
            (2.0, 6.0, 0.0),
            (0.0, 6.0, 0.0),
        ),
        dtype=np.float64,
    ),
    (),
    0.5,
    0.75,
    0.02,
    32,
)
assert_type(mat_result, MedialAxisResult)
mat_owner = _medial_axis_2.SegmentSiteMedialAxis.build(
    np.asarray(
        (
            (0.0, 0.0, 0.0),
            (6.0, 0.0, 0.0),
            (6.0, 2.0, 0.0),
            (2.0, 2.0, 0.0),
            (2.0, 6.0, 0.0),
            (0.0, 6.0, 0.0),
        ),
        dtype=np.float64,
    ),
    (),
    0.5,
    0.75,
    0.02,
    32,
)
assert_type(mat_owner, _medial_axis_2.SegmentSiteMedialAxis)

world_point_scalar = Point2[WorldXY].build(1.0, 2.0)
world_point_sequence = Point2[WorldXY].build([1.0, 2.0])
world_vector_scalar = Vector2[WorldXY].build(1.0, 2.0)
world_vector_sequence = Vector2[WorldXY].build((1.0, 2.0))

assert_type(world_point_scalar, Point2[WorldXY])
assert_type(world_point_sequence, Point2[WorldXY])
assert_type(world_vector_scalar, Vector2[WorldXY])
assert_type(world_vector_sequence, Vector2[WorldXY])

typed_boundary = CanonicalRingV1.build_outer(
    (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(6.0, 0.0),
        Point2[WorldXY].build(6.0, 2.0),
        Point2[WorldXY].build(2.0, 2.0),
        Point2[WorldXY].build(2.0, 6.0),
        Point2[WorldXY].build(0.0, 6.0),
    )
)
typed_axis = MedialAxis.build(
    design_boundary=typed_boundary,
    holes=(),
    tool_radius=ToolRadius.build(0.5),
    station_spacing=Spacing.build(0.75),
    max_sagitta=ChordBound.build(0.02),
    max_refinement_depth=32,
)
assert_type(typed_axis, MedialAxis)
assert_type(typed_axis.zero_guide_inventory, MatZeroGuideInventory)
assert_type(
    typed_axis.zero_guide_run_by_edge_id,
    Mapping[MatEdgeId, MatZeroGuideRun],
)
for zero_guide_run in typed_axis.zero_guide_inventory.runs:
    assert_type(zero_guide_run, MatZeroGuideRun)
    assert_type(zero_guide_run.edge_id, MatEdgeId)
    assert_type(zero_guide_run.native_certificate, bytes)
typed_traversal = MatTraversalState.seed(
    axis=typed_axis,
    inventory=NeckInventory.__new__(NeckInventory),
    policy=TraversalPolicy.build(forward_window=4),
    entry_edge_id=typed_axis.edges[0].identity,
    entry_node_id=typed_axis.edges[0].source.identity,
)
assert_type(typed_traversal, MatTraversalState)
assert_type(typed_traversal.authority.graph, TraversalGraph)
assert_type(typed_traversal.authority.route[0], DirectedRouteStep)
assert_type(typed_traversal.active_cursor, DirectedEdgeCursor)
assert_type(typed_traversal.neck_scope, NeckScope)
candidate_span = MiddleCurveSpan.build(
    axis=typed_axis,
    cursor_before=typed_axis.samples[0],
    cursor_limit=typed_axis.samples[1],
)
assert_type(candidate_span, MiddleCurveSpan)
assert_type(candidate_span.ordinal_step, Literal[-1, 1])
candidate_policy = CandidatePolicy.build(
    spatial_resolution=Spacing.build(0.5),
    spatial_refinement_levels=2,
    radius_resolution=Spacing.build(0.125),
    radius_refinement_levels=2,
    phase_count=4,
    minimum_guide_radius=GuideRadius.build(0.125),
    minimum_progress=Spacing.build(0.25),
)
candidate_cap = EngagementCap.build(1.0)
typed_candidates = enumerate_middle_curve_candidates(
    span=candidate_span,
    policy=candidate_policy,
    circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
    neck_scope=NoNeckScope.build(),
    effective_cap_decision=FullCapDecision.build(
        user_cap=candidate_cap,
        effective_cap=candidate_cap,
    ),
    makes_cursor_terminal_at_limit=False,
)
assert_type(typed_candidates, tuple[MiddleCurveCandidate, ...])
intermediate_candidate = next(candidate for candidate in typed_candidates if candidate.spatial_progress < candidate_span.reported_length)
derived_cursor = DerivedCandidateCursor.build(
    span=candidate_span,
    candidate=intermediate_candidate,
)
assert_type(derived_cursor, DerivedCandidateCursor)
assert_type(derived_cursor.ordinal_step, Literal[-1, 1])
assert_type(derived_cursor.next_limit_ordinal, int)
continued_span = MiddleCurveSpan.build(
    axis=typed_axis,
    cursor_before=derived_cursor,
    cursor_limit=typed_axis.samples[1],
)
assert_type(continued_span, MiddleCurveSpan)

typed_zero_axis = MedialAxis.build(
    design_boundary=typed_boundary,
    holes=(),
    tool_radius=ToolRadius.build(1.0),
    station_spacing=Spacing.build(0.75),
    max_sagitta=ChordBound.build(0.02),
    max_refinement_depth=32,
)
typed_zero_run = typed_zero_axis.zero_guide_inventory.runs[0]
typed_zero_samples = tuple(sample for sample in typed_zero_axis.samples if sample.edge_id == typed_zero_run.edge_id)
typed_zero_span = MiddleCurveSpan.build(
    axis=typed_zero_axis,
    cursor_before=typed_zero_samples[0],
    cursor_limit=typed_zero_samples[-1],
)
typed_zero_candidates = enumerate_zero_guide_link_candidates(
    span=typed_zero_span,
    policy=candidate_policy,
    neck_scope=NoNeckScope.build(),
    effective_cap_decision=FullCapDecision.build(
        user_cap=candidate_cap,
        effective_cap=candidate_cap,
    ),
    makes_cursor_terminal_at_limit=False,
)
assert_type(typed_zero_candidates, tuple[ZeroGuideLinkCandidate, ...])
typed_zero_candidate = typed_zero_candidates[-1]
assert_type(typed_zero_candidate, ZeroGuideLinkCandidate)
assert_type(typed_zero_candidate.zero_guide_run, MatZeroGuideRun)
zero_derived_cursor = DerivedCandidateCursor.build(
    span=typed_zero_span,
    candidate=typed_zero_candidate,
)
assert_type(zero_derived_cursor.candidate, TraversalCandidate)

approach_operation = ApproachOperation.build(
    position=world_point_scalar,
    clearance_z=ClearanceZ.build(5.0),
)
plunge_operation = PlungeOperation.build(
    position=world_point_scalar,
    clearance_z=ClearanceZ.build(5.0),
    cut_z=CutZ.build(-2.0),
)
assert_type(approach_operation, ApproachOperation)
assert_type(plunge_operation, PlungeOperation)

stock_area = Stock2Area.build(Stock.__new__(Stock))
depletion_policy = DepletionPolicy.build(
    chord_bound=ChordBound.build(0.25),
    center_count_limit=256,
)
segment_motion = ExactSegmentMotion.build(
    Point2[WorldXY].build(0.0, 0.0),
    Point2[WorldXY].build(1.0, 0.0),
)
circle_motion = ExactCircleMotion.build(
    Point2[WorldXY].build(0.0, 0.0),
    Vector2[WorldXY].build(1.0, 0.0),
    False,
)
segment_witness = stock_area.deplete(segment_motion, ToolRadius.build(0.5), depletion_policy)
circle_witness = stock_area.deplete(circle_motion, ToolRadius.build(0.5), depletion_policy)
assert_type(segment_witness, DepletionWitness)
assert_type(circle_witness, DepletionWitness)

reachable_domain = ReachableDomain.__new__(ReachableDomain)
containment = GougeContainment.build(reachable_domain)
segment_containment = containment.certify_segment(
    segment_motion,
    ToolRadius.build(0.5),
)
circle_containment = containment.certify_full_circle(
    circle_motion,
    ToolRadius.build(0.5),
)
assert_type(segment_containment, SegmentContainmentCertificate)
assert_type(circle_containment, CircleContainmentCertificate)

cut_plane = CutPlane.build(
    CutZ.build(-2.0),
    ClearanceZ.build(5.0),
)
qualified_bore = QualifiedBore.build(
    cut_plane=cut_plane,
    process_identity=BoreProcessIdentity(b"process"),
    evidence_bytes=b"evidence",
)
precleared_entry = PreclearedEntry.build(
    reachable_domain=reachable_domain,
    center=world_point_scalar,
    radius=EntryRadius.build(1.0),
    tool_radius=ToolRadius.build(0.5),
    cut_plane=cut_plane,
    qualified_bore=qualified_bore,
)
entry_witness = stock_area.deplete(precleared_entry)
assert_type(precleared_entry, PreclearedEntry)
assert_type(entry_witness, DepletionWitness)

input_identity = InputIdentity.build(
    design_boundary=typed_boundary,
    holes=(),
    cut_plane=cut_plane,
    tool_radius=ToolRadius.build(0.5),
    reachable_domain=reachable_domain,
    entry=precleared_entry,
    user_cap=candidate_cap,
    mat_sampling_policy=MatSamplingPolicy.build(
        station_spacing=Spacing.build(0.75),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    ),
    candidate_policy=candidate_policy,
    neck_policy=NeckPolicy.__new__(NeckPolicy),
    depletion_policy=depletion_policy,
    traversal_policy=TraversalPolicy.__new__(TraversalPolicy),
    cut_direction_policy=CutDirectionPolicy.__new__(CutDirectionPolicy),
)
assert_type(input_identity, InputIdentity)
initial_evaluator = InitialCandidateEvaluator.build(
    input_identity=input_identity,
    material_side=MaterialSide.OUTSIDE,
)
assert_type(initial_evaluator, InitialCandidateEvaluator)
initial_transaction = initial_evaluator.evaluate(
    typed_traversal,
    intermediate_candidate,
)
assert_type(initial_transaction, InitialCandidateTransaction)
initial_generation, initial_traversal = initial_evaluator.commit(
    typed_traversal,
    initial_transaction,
)
assert_type(initial_generation, GenerationState)
assert_type(initial_traversal, MatTraversalState)


def accepts_world_point(point: Point2[WorldXY]) -> None:
    pass


machine_point = Point2[MachineXY].build(1.0, 2.0)
accepts_world_point(machine_point)  # type: ignore[arg-type]
accepts_world_point(world_vector_scalar)  # type: ignore[arg-type]


def operation_kind(operation: CanonicalOperation) -> str:
    if isinstance(operation, ApproachOperation):
        return "approach"
    if isinstance(operation, PlungeOperation):
        return "plunge"
    if isinstance(operation, LinkSegmentOperation):
        return "link"
    if isinstance(operation, CutFullCircleOperation):
        return "circle"
    assert_never(operation)


def cap_kind(decision: FullCapDecision | NeckCapDecision) -> str:
    if isinstance(decision, FullCapDecision):
        return "full"
    if isinstance(decision, NeckCapDecision):
        return "neck"
    assert_never(decision)


def traversal_kind(decision: TraversalDecision) -> str:
    if isinstance(decision, HoldTraversalDecision):
        return "hold"
    if isinstance(decision, AdvanceTraversalDecision):
        return "advance"
    assert_never(decision)
