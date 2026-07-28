import math
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import MathsmCircleProposal
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.errors import InvalidMathsmProposalError
from compas_cgal.adaptive.medial_axis import MatEdge
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

L_SHAPE = np.asarray(
    (
        (0.0, 0.0, 0.0),
        (6.0, 0.0, 0.0),
        (6.0, 2.0, 0.0),
        (2.0, 2.0, 0.0),
        (2.0, 6.0, 0.0),
        (0.0, 6.0, 0.0),
    ),
    dtype=np.float64,
)


def _axis() -> MedialAxis:
    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y, _ in L_SHAPE),
    )
    return MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(0.5),
        station_spacing=Spacing.build(0.75),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )


def _constant_clearance_span() -> tuple[MiddleCurveSpan, MatEdge, MatSample]:
    axis = _axis()
    edge_samples = {edge.identity: tuple(sample for sample in axis.samples if sample.edge_id == edge.identity) for edge in axis.edges}
    edge = next(
        edge
        for edge in axis.edges
        if edge.curve_kind == "line" and len(edge_samples[edge.identity]) == 5 and edge_samples[edge.identity][0].point == Point2[WorldXY].build(2.0, 1.0)
    )
    samples = edge_samples[edge.identity]
    return (
        MiddleCurveSpan.build(
            axis=axis,
            cursor_before=samples[0],
            cursor_limit=samples[1],
        ),
        edge,
        samples[1],
    )


def _parabolic_span() -> MiddleCurveSpan:
    axis = _axis()
    edge = next(edge for edge in axis.edges if edge.curve_kind == "parabola")
    samples = tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)
    return MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[0],
        cursor_limit=samples[1],
    )


def _policy() -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.5),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=2,
        phase_count=4,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(0.25),
    )


def _full_cap() -> FullCapDecision:
    cap = EngagementCap.build(math.radians(90.0))
    return FullCapDecision.build(user_cap=cap, effective_cap=cap)


def _neck_policy() -> NeckPolicy:
    user_cap = EngagementCap.build(math.radians(120.0))
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


def test_complete_candidate_lattice_matches_exhaustive_oracle_and_mathsm() -> None:
    span, edge, cursor_limit = _constant_clearance_span()
    policy = _policy()
    candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )

    expected_cells = {
        (progress, radius, phase, site_id)
        for progress in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
        for radius in (Fraction(1, 8), Fraction(3, 16), Fraction(1, 4))
        for phase in range(4)
        for site_id in edge.generator_site_ids
    }
    observed_cells = {
        (
            candidate.spatial_progress,
            candidate.guide_radius,
            candidate.phase_index,
            candidate.generator_site_id,
        )
        for candidate in candidates
    }

    assert observed_cells == expected_cells
    assert len(candidates) == len(expected_cells) == 72
    assert {level for candidate in candidates for level in candidate.spatial_levels} == {
        0,
        1,
    }
    assert {level for candidate in candidates for level in candidate.radius_levels} == {
        0,
        1,
    }
    assert all(
        type(candidate.traversal_decision) is AdvanceTraversalDecision and candidate.traversal_decision.cursor_before == span.cursor_before.cursor_identity
        for candidate in candidates
    )
    assert {candidate.traversal_decision.cursor_after for candidate in candidates if candidate.spatial_progress == Fraction(3, 4)} == {cursor_limit.cursor_identity}

    base = next(
        candidate
        for candidate in candidates
        if candidate.spatial_progress == Fraction(1, 2) and candidate.guide_radius == Fraction(1, 4) and candidate.phase_index == 0 and candidate.boundary_footpoint.y == 0.0
    )
    assert base.middle_point == Point2[WorldXY].build(2.5, 1.0)
    assert base.boundary_footpoint == Point2[WorldXY].build(2.5, 0.0)
    assert base.mathsm_contact_point == Point2[WorldXY].build(2.5, 0.5)
    assert base.motion.center == Point2[WorldXY].build(2.5, 0.75)
    assert base.motion.phase_vector.x == 0.0
    assert base.motion.phase_vector.y == -0.25
    assert base.motion.squared_radius == Fraction(1, 16)

    repeated = enumerate_middle_curve_candidates(
        span=span,
        policy=policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )
    assert tuple(candidate.identity for candidate in candidates) == tuple(candidate.identity for candidate in repeated)
    assert candidates == policy.order_candidates(
        candidates,
        key=lambda candidate: candidate.order_key,
    )
    assert candidates[0].spatial_progress == Fraction(3, 4)
    assert candidates[0].guide_radius == Fraction(1, 4)


def test_candidate_identity_binds_neck_evidence_and_passage_cap_transition() -> None:
    span, _, _ = _constant_clearance_span()
    policy = _policy()
    neck_policy = _neck_policy()
    neck = NeckInventory.build(
        axis=span.axis,
        policy=neck_policy,
    ).necks[0]
    first_passage = neck.forward
    first_decision = first_passage.propose_cap_decision(neck_policy)
    second_passage = first_passage.advance(first_decision)
    second_decision = second_passage.propose_cap_decision(neck_policy)

    first_candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=first_passage.scope,
        effective_cap_decision=first_decision,
        makes_cursor_terminal_at_limit=False,
    )
    second_candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=second_passage.scope,
        effective_cap_decision=second_decision,
        makes_cursor_terminal_at_limit=False,
    )

    assert {
        (
            candidate.spatial_progress,
            candidate.guide_radius,
            candidate.phase_index,
            candidate.generator_site_id,
        )
        for candidate in first_candidates
    } == {
        (
            candidate.spatial_progress,
            candidate.guide_radius,
            candidate.phase_index,
            candidate.generator_site_id,
        )
        for candidate in second_candidates
    }
    assert {candidate.identity for candidate in first_candidates}.isdisjoint(candidate.identity for candidate in second_candidates)
    assert all(
        candidate.neck_scope == first_passage.scope and candidate.effective_cap_decision == first_decision and neck.evidence_digest in candidate.canonical_bytes
        for candidate in first_candidates
    )


def test_spatial_lattice_never_backsteps_below_non_aligned_minimum() -> None:
    span, _, _ = _constant_clearance_span()
    policy = CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.5),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=1,
        phase_count=1,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(0.3),
    )

    candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )

    assert {candidate.spatial_progress for candidate in candidates} == {
        Fraction(1, 2),
        Fraction(3, 4),
    }
    assert all(candidate.spatial_progress >= Fraction.from_float(policy.minimum_progress.value) for candidate in candidates)


def test_terminal_clearance_leaf_ends_at_last_feasible_candidate_cell() -> None:
    axis = _axis()
    edge = next(
        edge
        for edge in axis.edges
        if (
            (samples := tuple(sample for sample in axis.samples if sample.edge_id == edge.identity))[-1].clearance.value == axis.tool_radius.value
            and samples[0].clearance.value > axis.tool_radius.value
        )
    )
    samples = tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[0],
        cursor_limit=samples[-1],
    )
    candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=True,
    )
    terminal = tuple(candidate for candidate in candidates if candidate.traversal_decision.makes_cursor_terminal)
    maximum_feasible_progress = max(candidate.spatial_progress for candidate in candidates)

    assert terminal
    assert {candidate.spatial_progress for candidate in terminal} == {maximum_feasible_progress}
    assert maximum_feasible_progress < span.reported_length
    assert all(candidate.traversal_decision.makes_cursor_terminal == (candidate.spatial_progress == maximum_feasible_progress) for candidate in candidates)

    nonterminal_candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )

    assert not any(candidate.traversal_decision.makes_cursor_terminal for candidate in nonterminal_candidates)


def test_mathsm_factory_fails_named_before_reading_untyped_geometry() -> None:
    span, edge, _ = _constant_clearance_span()
    site = span.axis.site_by_id[edge.generator_site_ids[0]]

    with pytest.raises(InvalidMathsmProposalError):
        MathsmCircleProposal.build(
            generator_site=site,
            middle_point=object(),  # type: ignore[arg-type]
            tool_radius=span.axis.tool_radius,
            guide_radius=Fraction(1, 8),
            phase_index=0,
            phase_count=4,
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        )


def test_spatial_refinement_follows_point_segment_parabola_not_sample_chord() -> None:
    span = _parabolic_span()
    half_span = float(span.reported_length) / 2.0
    policy = CandidatePolicy.build(
        spatial_resolution=Spacing.build(half_span),
        spatial_refinement_levels=1,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=1,
        phase_count=1,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(half_span),
    )

    candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=policy,
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )
    refined_points = {candidate.middle_point for candidate in candidates if candidate.spatial_progress < span.reported_length}

    assert refined_points == {Point2[WorldXY].build(1.75, 1.015625)}
    assert Point2[WorldXY].build(1.75, 1.03125) not in refined_points
