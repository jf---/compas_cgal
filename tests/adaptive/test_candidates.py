import math
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import ExhaustedCandidateCursor
from compas_cgal.adaptive.candidates import MathsmCircleProposal
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.candidates import enumerate_zero_guide_link_candidates
from compas_cgal.adaptive.errors import InvalidZeroGuideCandidateError
from compas_cgal.adaptive.errors import InvalidMathsmProposalError
from compas_cgal.adaptive.errors import InvalidMiddleCurveCursorError
from compas_cgal.adaptive.errors import InvalidMiddleCurveSpanError
from compas_cgal.adaptive.errors import UncertifiedZeroGuideEdgeError
from compas_cgal.adaptive.medial_axis import MatEdge
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import EdgeId
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


def _axis(*, tool_radius: float = 0.5) -> MedialAxis:
    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y, _ in L_SHAPE),
    )
    return MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(tool_radius),
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


def _clearance_leaf_span() -> MiddleCurveSpan:
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
    return MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[0],
        cursor_limit=samples[-1],
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


def _zero_guide_spans() -> tuple[tuple[MiddleCurveSpan, MiddleCurveSpan], ...]:
    axis = _axis(tool_radius=1.0)
    spans = []
    for run in axis.zero_guide_inventory.runs:
        samples = tuple(sample for sample in axis.samples if sample.edge_id == run.edge_id)
        spans.append(
            (
                MiddleCurveSpan.build(
                    axis=axis,
                    cursor_before=samples[0],
                    cursor_limit=samples[-1],
                ),
                MiddleCurveSpan.build(
                    axis=axis,
                    cursor_before=samples[-1],
                    cursor_limit=samples[0],
                ),
            )
        )
    return tuple(spans)


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


def test_derived_candidate_cursor_continues_the_same_exact_span_grammar() -> None:
    span, _, cursor_limit = _constant_clearance_span()
    candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )
    selected = next(candidate for candidate in candidates if candidate.spatial_progress == Fraction(1, 4))

    cursor = DerivedCandidateCursor.build(
        span=span,
        candidate=selected,
    )
    continuation = MiddleCurveSpan.build(
        axis=span.axis,
        cursor_before=cursor,
        cursor_limit=cursor_limit,
    )
    continued_candidates = enumerate_middle_curve_candidates(
        span=continuation,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )

    assert cursor.cursor_identity == selected.traversal_decision.cursor_after
    assert cursor.point == selected.middle_point
    assert continuation.cursor_before is cursor
    assert all(candidate.traversal_decision.cursor_before == selected.traversal_decision.cursor_after for candidate in continued_candidates)
    cross_wired = next(candidate for candidate in continued_candidates if candidate.spatial_progress < continuation.reported_length)
    with pytest.raises(
        InvalidMiddleCurveCursorError,
        match="claimed exact span",
    ):
        DerivedCandidateCursor.build(
            span=span,
            candidate=cross_wired,
        )


def test_derived_candidate_cursor_rejects_native_or_terminal_endpoint() -> None:
    span, _, _ = _constant_clearance_span()
    nonterminal_candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=False,
    )
    native_endpoint = next(candidate for candidate in nonterminal_candidates if candidate.spatial_progress == span.reported_length)
    with pytest.raises(
        InvalidMiddleCurveCursorError,
        match="native span endpoint",
    ):
        DerivedCandidateCursor.build(
            span=span,
            candidate=native_endpoint,
        )

    terminal_candidates = enumerate_middle_curve_candidates(
        span=span,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=True,
    )
    terminal = next(candidate for candidate in terminal_candidates if candidate.traversal_decision.makes_cursor_terminal)

    with pytest.raises(
        InvalidMiddleCurveCursorError,
        match="terminal candidate",
    ):
        DerivedCandidateCursor.build(
            span=span,
            candidate=terminal,
        )


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
    span = _clearance_leaf_span()
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


def test_exhausted_cursor_rejects_nonterminal_and_native_endpoint_candidates() -> None:
    """Close only a terminal positive-radius station before its native bound.

    `DerivedCandidateCursor` owns continuation and `MatSample` owns a reached
    native endpoint. The exhausted cursor must reject both states so terminal
    identity cannot silently change representation.
    """
    clearance_span = _clearance_leaf_span()
    clearance_candidates = enumerate_middle_curve_candidates(
        span=clearance_span,
        policy=_policy(),
        circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=True,
    )
    exhausted_candidate = next(candidate for candidate in clearance_candidates if candidate.traversal_decision.makes_cursor_terminal)
    exhausted = ExhaustedCandidateCursor.build(
        span=clearance_span,
        candidate=exhausted_candidate,
    )
    nonterminal = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=clearance_span,
            policy=_policy(),
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == exhausted_candidate.spatial_progress
    )
    native_span, _, _ = _constant_clearance_span()
    native_terminal = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=native_span,
            policy=_policy(),
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.traversal_decision.makes_cursor_terminal
    )

    assert exhausted.cursor_identity == exhausted_candidate.traversal_decision.cursor_after
    with pytest.raises(
        InvalidMiddleCurveCursorError,
        match="terminal candidate",
    ):
        ExhaustedCandidateCursor.build(
            span=clearance_span,
            candidate=nonterminal,
        )
    with pytest.raises(
        InvalidMiddleCurveCursorError,
        match="native terminal endpoint",
    ):
        ExhaustedCandidateCursor.build(
            span=native_span,
            candidate=native_terminal,
        )


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


def test_reverse_line_span_preserves_positive_progress_and_separate_identity() -> None:
    """Traverse a real constant-clearance MAT line from its higher ordinal.

    Connected graph discovery may enter this edge at either endpoint. The
    reverse lattice must keep positive machining progress and deterministic
    ordering while remaining cryptographically distinct from the forward
    traversal of the same native interval.
    """
    base, edge, _ = _constant_clearance_span()
    samples = tuple(
        sorted(
            (sample for sample in base.axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    forward = MiddleCurveSpan.build(
        axis=base.axis,
        cursor_before=samples[1],
        cursor_limit=samples[2],
    )
    reverse = MiddleCurveSpan.build(
        axis=base.axis,
        cursor_before=samples[2],
        cursor_limit=samples[1],
    )
    arguments = {
        "policy": _policy(),
        "circle_orientation": CircleOrientation.COUNTERCLOCKWISE,
        "neck_scope": NoNeckScope.build(),
        "effective_cap_decision": _full_cap(),
        "makes_cursor_terminal_at_limit": True,
    }
    forward_candidates = enumerate_middle_curve_candidates(
        span=forward,
        **arguments,
    )
    reverse_candidates = enumerate_middle_curve_candidates(
        span=reverse,
        **arguments,
    )
    repeated = enumerate_middle_curve_candidates(
        span=reverse,
        **arguments,
    )

    assert reverse.ordinal_step == -1
    assert {candidate.spatial_progress for candidate in reverse_candidates} == {candidate.spatial_progress for candidate in forward_candidates}
    assert all(candidate.spatial_progress > 0 for candidate in reverse_candidates)
    assert {candidate.identity for candidate in reverse_candidates}.isdisjoint(candidate.identity for candidate in forward_candidates)
    assert tuple(candidate.identity for candidate in reverse_candidates) == tuple(candidate.identity for candidate in repeated)
    terminal = tuple(candidate for candidate in reverse_candidates if candidate.traversal_decision.makes_cursor_terminal)
    assert terminal
    assert {candidate.spatial_progress for candidate in terminal} == {reverse.reported_length}


def test_zero_guide_lattice_covers_both_proved_arms_in_both_directions() -> None:
    """Enumerate only spatial progress from each exact constant-profile proof."""
    arguments = {
        "policy": _policy(),
        "neck_scope": NoNeckScope.build(),
        "effective_cap_decision": _full_cap(),
        "makes_cursor_terminal_at_limit": True,
    }
    span_pairs = _zero_guide_spans()

    assert len(span_pairs) == 2
    for forward, reverse in span_pairs:
        forward_candidates = enumerate_zero_guide_link_candidates(
            span=forward,
            **arguments,
        )
        reverse_candidates = enumerate_zero_guide_link_candidates(
            span=reverse,
            **arguments,
        )
        repeated = enumerate_zero_guide_link_candidates(
            span=reverse,
            **arguments,
        )

        assert len(forward_candidates) == len(reverse_candidates) == 12
        assert all(type(candidate) is ZeroGuideLinkCandidate for candidate in forward_candidates)
        assert all(candidate.zero_guide_run.edge_id == forward.edge.identity for candidate in forward_candidates)
        assert forward_candidates[0].spatial_progress == max(candidate.spatial_progress for candidate in forward_candidates)
        assert forward_candidates[0].traversal_decision.makes_cursor_terminal
        assert not hasattr(forward_candidates[0], "guide_radius")
        assert not hasattr(forward_candidates[0], "phase_index")
        assert not hasattr(forward_candidates[0], "motion")
        assert b"zero-guide-link-candidate-v1" in forward_candidates[0].canonical_bytes
        assert {candidate.identity for candidate in forward_candidates}.isdisjoint(candidate.identity for candidate in reverse_candidates)
        assert tuple(candidate.identity for candidate in reverse_candidates) == tuple(candidate.identity for candidate in repeated)


def test_zero_guide_candidate_produces_owned_derived_cursor() -> None:
    """An interior link station retains the same span and cursor identity grammar."""
    span = _zero_guide_spans()[0][0]
    candidate = next(
        candidate
        for candidate in enumerate_zero_guide_link_candidates(
            span=span,
            policy=_policy(),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
    )
    cursor = DerivedCandidateCursor.build(
        span=span,
        candidate=candidate,
    )

    assert cursor.candidate is candidate
    assert cursor.point == candidate.target
    assert cursor.cursor_identity == candidate.traversal_decision.cursor_after


def test_zero_guide_enumerator_rejects_unproved_positive_guide_edge() -> None:
    """A zero reporting value or endpoint touch cannot replace native proof."""
    axis = _axis(tool_radius=1.0)
    edge = next(edge for edge in axis.edges if edge.identity not in axis.zero_guide_run_by_edge_id)
    samples = tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=samples[0],
        cursor_limit=samples[-1],
    )

    with pytest.raises(UncertifiedZeroGuideEdgeError, match="native proof"):
        enumerate_zero_guide_link_candidates(
            span=span,
            policy=_policy(),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            makes_cursor_terminal_at_limit=True,
        )


def test_zero_guide_candidate_rejects_every_cross_wired_decision_field() -> None:
    """All candidate identity inputs remain closed under raw dataclass mutation."""
    span_pairs = _zero_guide_spans()
    span = span_pairs[0][0]
    candidate = enumerate_zero_guide_link_candidates(
        span=span,
        policy=_policy(),
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        makes_cursor_terminal_at_limit=True,
    )[0]
    foreign_run = span_pairs[1][0].axis.zero_guide_run_by_edge_id[span_pairs[1][0].edge.identity]
    foreign_traversal = AdvanceTraversalDecision.build(
        component_id=candidate.traversal_decision.component_id,
        edge_id=EdgeId(bytes(foreign_run.edge_id)),
        branch_id=candidate.traversal_decision.branch_id,
        cursor_before=candidate.traversal_decision.cursor_before,
        cursor_after=candidate.traversal_decision.cursor_after,
        makes_cursor_terminal=candidate.traversal_decision.makes_cursor_terminal,
    )
    mutations = (
        {"zero_guide_run": foreign_run},
        {"target": Point2[WorldXY].build(candidate.target.x, candidate.target.y + 0.125)},
        {"spatial_progress": Fraction(0)},
        {"spatial_levels": (1, 0)},
        {"cursor_limit_identity": b"foreign-cursor"},
        {"neck_scope": object()},
        {"effective_cap_decision": object()},
        {"traversal_decision": foreign_traversal},
    )

    for changes in mutations:
        with pytest.raises(InvalidZeroGuideCandidateError):
            replace(candidate, **changes)  # type: ignore[arg-type]


def test_reverse_parabola_uses_the_same_point_segment_geometry() -> None:
    """Mirror one real point-segment parabola without mirroring its formula.

    Reversing the ordered endpoints must reach the identical half-parameter
    point on the analytic parabola. Candidate lineage still differs because
    the exact cursor-before and cursor-limit records are reversed.
    """
    forward = _parabolic_span()
    reverse = MiddleCurveSpan.build(
        axis=forward.axis,
        cursor_before=forward.cursor_limit,
        cursor_limit=forward.cursor_before,
    )
    half_span = float(forward.reported_length) / 2.0
    policy = CandidatePolicy.build(
        spatial_resolution=Spacing.build(half_span),
        spatial_refinement_levels=1,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=1,
        phase_count=1,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(half_span),
    )
    arguments = {
        "policy": policy,
        "circle_orientation": CircleOrientation.COUNTERCLOCKWISE,
        "neck_scope": NoNeckScope.build(),
        "effective_cap_decision": _full_cap(),
        "makes_cursor_terminal_at_limit": False,
    }
    forward_candidates = enumerate_middle_curve_candidates(
        span=forward,
        **arguments,
    )
    reverse_candidates = enumerate_middle_curve_candidates(
        span=reverse,
        **arguments,
    )
    forward_midpoints = {candidate.middle_point for candidate in forward_candidates if candidate.spatial_progress < forward.reported_length}
    reverse_midpoints = {candidate.middle_point for candidate in reverse_candidates if candidate.spatial_progress < reverse.reported_length}

    assert reverse.ordinal_step == -1
    assert (
        forward_midpoints
        == reverse_midpoints
        == {
            Point2[WorldXY].build(1.75, 1.015625),
        }
    )
    assert {candidate.identity for candidate in reverse_candidates}.isdisjoint(candidate.identity for candidate in forward_candidates)


def test_derived_reverse_cursor_rejects_direction_reversal() -> None:
    """Retain reverse ordinal causality after stopping between native samples.

    A derived cursor has no native ordinal of its own. Its producing span must
    therefore carry the direction and nearest legal native limit; otherwise a
    continuation could backtrack while presenting one uninterrupted lineage.
    """
    forward, edge, _ = _constant_clearance_span()
    samples = tuple(
        sorted(
            (sample for sample in forward.axis.samples if sample.edge_id == edge.identity),
            key=lambda sample: sample.ordinal_on_edge,
        )
    )
    reverse = MiddleCurveSpan.build(
        axis=forward.axis,
        cursor_before=samples[2],
        cursor_limit=samples[1],
    )
    candidate = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=reverse,
            policy=_policy(),
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == Fraction(1, 4)
    )
    cursor = DerivedCandidateCursor.build(
        span=reverse,
        candidate=candidate,
    )
    continuation = MiddleCurveSpan.build(
        axis=forward.axis,
        cursor_before=cursor,
        cursor_limit=samples[0],
    )

    assert cursor.ordinal_step == -1
    assert cursor.next_limit_ordinal == samples[1].ordinal_on_edge
    assert continuation.ordinal_step == -1
    with pytest.raises(
        InvalidMiddleCurveSpanError,
        match="retained sample direction",
    ):
        MiddleCurveSpan.build(
            axis=forward.axis,
            cursor_before=cursor,
            cursor_limit=samples[2],
        )
