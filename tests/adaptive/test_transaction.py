import math
from dataclasses import dataclass
from fractions import Fraction

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
from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import CutFullCircleOperation
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
from compas_cgal.stock import Stock

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


def _first_candidate(
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
    full_cap = FullCapDecision.build(
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
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 2
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )


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
) -> tuple[MiddleCurveCandidate, MiddleCurveCandidate]:
    first = _first_candidate(
        user_cap=identity.user_cap,
        candidate_policy=identity.candidate_policy,
    )
    axis = _axis()
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
    terminal_span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=DerivedCandidateCursor.build(
            span=first_span,
            candidate=first,
        ),
        cursor_limit=samples[2],
    )
    full_cap = FullCapDecision.build(
        user_cap=identity.user_cap,
        effective_cap=identity.user_cap,
    )
    second = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=terminal_span,
            policy=identity.candidate_policy,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.spatial_progress == Fraction(1, 4)
        and candidate.guide_radius == Fraction(1, 8)
        and candidate.phase_index == 3
        and candidate.proposal.generator_site.kind == "point"
        and candidate.proposal.generator_site.source == Point2[WorldXY].build(2.0, 2.0)
    )
    return first, second


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


def _state_fixture() -> _StateFixture:
    identity = _input_identity(neck_cap_degrees=(120.0, 120.0, 120.0))
    first, second = _candidate_pair(identity)
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
        passages=_passages(identity),
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
        state,
        stock,
        coverage,
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
