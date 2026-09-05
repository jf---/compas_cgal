import math
from dataclasses import replace
from fractions import Fraction

import pytest
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_path import build_hypothesis_figure5_path
from benchmarks.held_figure5_path import _canonicalize
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import ProjectionBoundaryVertexId
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


SQUARE = (_point(0.0, 0.0), _point(4.0, 0.0), _point(4.0, 4.0), _point(0.0, 4.0))
TOOL_RADIUS = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80.0))


def _hypothesis(
    run: int,
    x: float,
    *,
    ordinal: int = 0,
    center_y: float = 1.0,
    raw_shift: float = 0.0,
) -> ProjectionAdmissibleBoundaryHypothesis:
    site = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(0),
        Fraction.from_float(x / 4.0),
        len(SQUARE),
    )
    contact = _point(x + raw_shift, 0.0)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(run),
        GuideRunStationOrdinal(ordinal),
        site,
        (site,),
        Millimetre(x),
        _point(x + raw_shift, center_y + 1.0),
        contact,
        contact,
        _point(x + raw_shift, center_y),
        GuideRadius.build(center_y),
    )


def test_exact_circle_canonicalization_is_permutation_invariant_and_keeps_different_q() -> None:
    base = _hypothesis(2, 1.0)
    same_geometry = replace(base, run_id=GuideRunId(1))
    ulp_shifted = replace(
        base,
        run_id=GuideRunId(0),
        center=_point(math.nextafter(float(base.center.x), math.inf), float(base.center.y)),
    )
    different_q = replace(base, run_id=GuideRunId(3), contact_point=_point(1.5, 0.0))

    forward = _canonicalize((base, same_geometry, ulp_shifted, different_q))
    reverse = _canonicalize((different_q, ulp_shifted, same_geometry, base))

    assert forward == reverse
    assert len(forward) == 3
    assert forward[0].source_run_ids == (GuideRunId(1), GuideRunId(2))


def _vertex_hypothesis(
    ordinal: int,
    contact: tuple[float, float],
    center: tuple[float, float],
) -> ProjectionAdmissibleBoundaryHypothesis:
    site = ProjectionBoundarySite(
        ProjectionBoundarySegmentId(0),
        Fraction(0),
        ProjectionBoundaryVertexId(0),
    )
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(9),
        GuideRunStationOrdinal(ordinal),
        site,
        (site,),
        Millimetre(0.0),
        _point(*center),
        _point(*contact),
        _point(*contact),
        _point(*center),
        GuideRadius.build(1.0),
    )


def _side_hypothesis(
    *,
    run: int,
    ordinal: int,
    side: int,
    parameter: Fraction,
) -> ProjectionAdmissibleBoundaryHypothesis:
    start = SQUARE[side]
    end = SQUARE[(side + 1) % len(SQUARE)]
    contact = _point(
        float(start.x) + float(parameter) * (float(end.x) - float(start.x)),
        float(start.y) + float(parameter) * (float(end.y) - float(start.y)),
    )
    inward = ((0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0))[side]
    center = _point(float(contact.x) + inward[0], float(contact.y) + inward[1])
    site = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(side),
        parameter,
        len(SQUARE),
    )
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(run),
        GuideRunStationOrdinal(ordinal),
        site,
        (site,),
        Millimetre(0.0),
        center,
        contact,
        contact,
        center,
        GuideRadius.build(1.0),
    )


def _ownership_hypothesis(
    *,
    side: int,
    boundary_distance: float,
) -> ProjectionAdmissibleBoundaryHypothesis:
    parameter = Fraction(1, 2)
    start = SQUARE[side]
    end = SQUARE[(side + 1) % len(SQUARE)]
    contact = _point(
        (float(start.x) + float(end.x)) / 2.0,
        (float(start.y) + float(end.y)) / 2.0,
    )
    inward = ((0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0))[side]
    center = _point(float(contact.x) + inward[0], float(contact.y) + inward[1])
    site = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(side),
        parameter,
        len(SQUARE),
    )
    middle = _point(2.0, 2.0)
    footpoint = _point(2.0 + boundary_distance, 2.0)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(0),
        GuideRunStationOrdinal(0),
        site,
        (site,),
        Millimetre(0.0),
        middle,
        footpoint,
        contact,
        center,
        GuideRadius.build(1.0),
    )


def test_places_actual_circles_at_cap_and_preserves_reached_runs() -> None:
    target_spacing = 2.0 - math.sqrt(2.0 + 2.0 * math.cos(math.radians(80.0)))
    candidates = tuple(
        _hypothesis(
            0,
            0.5 + index * target_spacing,
            ordinal=index,
            raw_shift=index * 0.001,
        )
        for index in range(4)
    )

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=candidates,
        reached_run_ids=frozenset(candidate.run_id for candidate in candidates),
        published_start_center=candidates[0].contact_point,
        published_start_radius=Millimetre(float(candidates[0].guide_radius.value)),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )

    assert result.candidate_count == 4
    assert result.canonical_candidate_count == 4
    assert result.reached_run_ids == frozenset((GuideRunId(0),))
    assert len(result.path.circles) == 4
    assert all(circle.clockwise is False for circle in result.path.circles)
    assert len(result.placements) == 3
    assert all(float(record.placement.maximum_engagement) == pytest.approx(math.radians(80.0)) for record in result.placements)
    assert all(float(record.placement.center_spacing) == pytest.approx(target_spacing) for record in result.placements)
    assert float(candidates[1].center.x) - float(candidates[0].center.x) != pytest.approx(target_spacing)
    assert all(transition.boundary_progress > 0.0 for transition in result.path.transitions)


def test_placement_never_compares_interleaved_unrelated_guide_runs() -> None:
    candidates = (
        _hypothesis(0, 0.5, ordinal=0),
        _hypothesis(1, 0.6, ordinal=0, center_y=3.0),
        _hypothesis(0, 1.0, ordinal=1),
        _hypothesis(1, 1.1, ordinal=1, center_y=3.0),
    )

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=candidates,
        reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1))),
        published_start_center=candidates[0].contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=EngagementCap.build(math.radians(90.0)),
    )

    assert len(result.path.circles) == 4
    assert len(result.placements) == 2
    assert tuple(record.source_run_ids for record in result.placements) == (
        (GuideRunId(0),),
        (GuideRunId(1),),
    )


def test_placement_joins_compatible_continuous_side_fragments_across_runs() -> None:
    candidates = tuple(_hypothesis(run, x, ordinal=ordinal) for run, xs in ((0, (0.1, 0.2, 0.3)), (1, (0.4, 0.5, 0.6))) for ordinal, x in enumerate(xs))

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=candidates,
        reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1))),
        published_start_center=candidates[0].contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=EngagementCap.build(math.radians(80.0)),
    )

    assert len(result.path.circles) == 3
    assert len(result.placements) == 2
    assert result.placement_lane_count == 1
    assert result.reached_run_ids == frozenset((GuideRunId(0), GuideRunId(1)))


def test_placement_joins_unique_nearest_continuation_among_compatible_tails() -> None:
    candidates = tuple(_hypothesis(run, x, ordinal=ordinal) for run, xs in ((0, (0.1, 0.2)), (1, (0.05, 0.15)), (2, (0.3, 0.4))) for ordinal, x in enumerate(xs))

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=candidates,
        reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1), GuideRunId(2))),
        published_start_center=candidates[2].contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=EngagementCap.build(math.radians(80.0)),
    )

    assert result.placement_lane_count == 2


def test_placement_keeps_exactly_tied_continuations_as_distinct_fan_branches() -> None:
    candidates = tuple(
        _hypothesis(run, x, ordinal=ordinal, center_y=radius)
        for run, radius, xs in (
            (0, 0.5, (0.1, 0.2)),
            (1, 1.5, (0.1, 0.2)),
            (2, 1.0, (0.3, 0.4)),
        )
        for ordinal, x in enumerate(xs)
    )

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=candidates,
        reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1), GuideRunId(2))),
        published_start_center=candidates[0].contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=EngagementCap.build(math.radians(80.0)),
    )

    assert result.placement_lane_count == 3


def test_placement_keeps_one_sided_hypothesis_lanes_separate_within_run() -> None:
    candidates = tuple(
        _side_hypothesis(
            run=0,
            ordinal=ordinal,
            side=side,
            parameter=parameter,
        )
        for ordinal, parameter in enumerate((Fraction(1, 8), Fraction(1, 4)))
        for side in (0, 1)
    )

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=candidates,
        reached_run_ids=frozenset((GuideRunId(0),)),
        published_start_center=candidates[0].contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=EngagementCap.build(math.radians(90.0)),
    )

    assert len(result.path.circles) == 4
    assert len(result.placements) == 2


def test_station_owns_nearest_hypothesis_in_one_contiguous_boundary_neighborhood() -> None:
    nearest = _ownership_hypothesis(side=0, boundary_distance=1.0)
    adjacent_alternative = _ownership_hypothesis(side=1, boundary_distance=2.0)

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=(nearest, adjacent_alternative),
        reached_run_ids=frozenset((GuideRunId(0),)),
        published_start_center=nearest.contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )

    assert result.candidate_count == 2
    assert result.canonical_candidate_count == 1
    assert result.path.circles[0].source_boundary_site == nearest.projected_boundary_site


def test_station_preserves_disconnected_boundary_neighborhoods_and_exact_ties() -> None:
    first = _ownership_hypothesis(side=0, boundary_distance=1.0)
    disconnected = _ownership_hypothesis(side=2, boundary_distance=2.0)

    disconnected_result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=(first, disconnected),
        reached_run_ids=frozenset((GuideRunId(0),)),
        published_start_center=first.contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )
    tied_adjacent = _ownership_hypothesis(side=1, boundary_distance=1.0)
    tied_result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=(first, tied_adjacent),
        reached_run_ids=frozenset((GuideRunId(0),)),
        published_start_center=first.contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )

    assert disconnected_result.canonical_candidate_count == 2
    assert tied_result.canonical_candidate_count == 2


def test_publisher_start_witness_survives_neighborhood_minimum() -> None:
    nearest = _ownership_hypothesis(side=0, boundary_distance=1.0)
    published_start_witness = _ownership_hypothesis(side=1, boundary_distance=2.0)

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=(nearest, published_start_witness),
        reached_run_ids=frozenset((GuideRunId(0),)),
        published_start_center=published_start_witness.contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )

    assert result.canonical_candidate_count == 2
    assert result.path.circles[0].source_boundary_site == published_start_witness.projected_boundary_site


def test_merges_only_equivalent_same_site_hypotheses_and_preserves_contributors() -> None:
    first = _hypothesis(3, 1.0)
    equivalent = _hypothesis(7, 1.0)

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=(first, equivalent),
        reached_run_ids=frozenset((first.run_id, equivalent.run_id)),
        published_start_center=first.contact_point,
        published_start_radius=Millimetre(float(first.guide_radius.value)),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )

    assert result.candidate_count == 2
    assert result.canonical_candidate_count == 1
    assert result.circle_sources == ((first.run_id, equivalent.run_id),)

    fan_result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=(first, _hypothesis(7, 1.0, center_y=1.25)),
        reached_run_ids=frozenset((GuideRunId(3), GuideRunId(7))),
        published_start_center=first.contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=CAP,
    )

    assert fan_result.canonical_candidate_count == 2
    assert all(transition.within_fan for transition in fan_result.path.transitions)


def test_retains_non_equivalent_vertex_fan_in_run_station_order() -> None:
    fan = (
        _vertex_hypothesis(0, (0.1, 0.0), (0.1, 1.0)),
        _vertex_hypothesis(1, (0.0, 0.1), (1.0, 0.1)),
    )

    result = build_hypothesis_figure5_path(
        inward_components=(SQUARE,),
        hypotheses=fan,
        reached_run_ids=frozenset((GuideRunId(9),)),
        published_start_center=fan[0].contact_point,
        published_start_radius=Millimetre(1.0),
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
        tool_radius=TOOL_RADIUS,
        cap=EngagementCap.build(math.radians(160.0)),
    )

    assert result.canonical_candidate_count == 2
    assert tuple(circle.station_ordinal for circle in result.path.circles) == (
        GuideRunStationOrdinal(0),
        GuideRunStationOrdinal(1),
    )


def test_builds_canonical_figure5_path_from_every_raw_guide_run() -> None:
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)

    result = build_figure5_approximate_path(
        guide,
        figure7_inward_offset(case).components,
    )

    assert result.candidate_count == 31_742
    assert result.canonical_candidate_count == 9_370
    assert result.placement_lane_count == 187
    assert len(result.path.circles) == 1_532
    assert len(result.placements) == 1_345
    assert sum(record.placement.forced for record in result.placements) == 120
    assert result.reached_run_ids == frozenset(run.run_id for run in guide.runs)
    assert frozenset(run_id for sources in result.circle_sources for run_id in sources) == result.reached_run_ids
    assert all(circle.clockwise is False for circle in result.path.circles)
