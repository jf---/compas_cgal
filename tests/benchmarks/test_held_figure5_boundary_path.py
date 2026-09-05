from fractions import Fraction

import pytest

from benchmarks.held_figure5_boundary_path import AmbiguousFigure5BoundaryPathError
from benchmarks.held_figure5_boundary_path import DisconnectedFigure5BoundaryPathError
from benchmarks.held_figure5_boundary_path import MissingFigure5GuideRunError
from benchmarks.held_figure5_boundary_path import build_figure5_boundary_path
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


SQUARE = (_point(0.0, 0.0), _point(4.0, 0.0), _point(4.0, 4.0), _point(0.0, 4.0))
TOOL_RADIUS = ToolRadius.build(1.0)


def _candidate(
    run: int,
    side: int,
    fraction: Fraction,
    centre: tuple[float, float],
    radius: float = 1.0,
    *,
    station: int = 0,
) -> ProjectionAdmissibleBoundaryHypothesis:
    start = SQUARE[side]
    end = SQUARE[(side + 1) % len(SQUARE)]
    q = _point(
        float(start.x) + float(fraction) * (float(end.x) - float(start.x)),
        float(start.y) + float(fraction) * (float(end.y) - float(start.y)),
    )
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), fraction, len(SQUARE))
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(run),
        GuideRunStationOrdinal(station),
        site,
        (site,),
        Millimetre(side * 4.0 + float(fraction) * 4.0),
        _point(*centre),
        q,
        q,
        _point(*centre),
        GuideRadius.build(radius),
    )


def _build(candidates: tuple[ProjectionAdmissibleBoundaryHypothesis, ...], *, start_index: int = 0):
    start = candidates[start_index]
    return build_figure5_boundary_path(
        inward_components=(SQUARE,),
        side_candidates=candidates,
        reached_run_ids=frozenset(candidate.run_id for candidate in candidates),
        published_start_center=start.contact_point,
        published_start_radius=Millimetre(TOOL_RADIUS.value),
        tool_radius=TOOL_RADIUS,
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
    )


def test_builds_one_continuous_ccw_circle_and_transition_traversal() -> None:
    candidates = (
        _candidate(2, 2, Fraction(1, 4), (3.0, 3.0), 0.8),
        _candidate(0, 0, Fraction(1, 4), (1.0, 1.0), 0.7),
        _candidate(3, 3, Fraction(3, 4), (1.0, 1.1), 0.9),
        _candidate(1, 1, Fraction(1, 2), (3.0, 2.0), 0.6),
    )

    path = _build(candidates)

    assert tuple(circle.run_id for circle in path.circles) == (GuideRunId(2), GuideRunId(3), GuideRunId(0), GuideRunId(1))
    assert all(circle.clockwise is False for circle in path.circles)
    assert tuple(circle.radius for circle in path.circles) == tuple(candidate.guide_radius for candidate in (candidates[0], candidates[2], candidates[1], candidates[3]))
    for index, transition in enumerate(path.transitions):
        assert transition.boundary_progress > 0.0
        assert transition.samples[0] == path.circles[index].contact_point
        assert transition.samples[-1] == path.circles[(index + 1) % len(path.circles)].contact_point
        assert transition.samples[-1] == path.transitions[(index + 1) % len(path.transitions)].samples[0]


def test_mid_side_seed_keeps_the_single_ccw_wrap_at_its_true_position() -> None:
    candidates = (
        _candidate(0, 0, Fraction(1, 4), (1.0, 1.0)),
        _candidate(1, 0, Fraction(1, 2), (2.0, 1.0)),
        _candidate(2, 0, Fraction(3, 4), (3.0, 1.0)),
    )

    path = _build(candidates, start_index=1)

    assert tuple(circle.run_id for circle in path.circles) == (GuideRunId(1), GuideRunId(2), GuideRunId(0))
    assert len(path.transitions[0].samples) == 2
    assert len(path.transitions[1].samples) == 6
    assert len(path.transitions[2].samples) == 2
    assert all(transition.boundary_progress > 0.0 for transition in path.transitions)


def test_same_q_fan_preserves_run_station_order_and_explicit_zero_connectors() -> None:
    candidates = (
        _candidate(2, 0, Fraction(1, 4), (1.0, 0.8), 0.8, station=3),
        _candidate(1, 1, Fraction(1, 2), (3.0, 2.0), 0.6),
        _candidate(1, 0, Fraction(1, 4), (1.0, 1.0), 1.0, station=5),
        _candidate(1, 0, Fraction(1, 4), (1.0, 1.2), 1.2, station=6),
    )

    path = _build(candidates)

    assert tuple((circle.run_id, circle.station_ordinal) for circle in path.circles[:3]) == (
        (GuideRunId(1), GuideRunStationOrdinal(5)),
        (GuideRunId(1), GuideRunStationOrdinal(6)),
        (GuideRunId(2), GuideRunStationOrdinal(3)),
    )
    assert tuple(float(circle.radius.value) for circle in path.circles[:3]) == (1.0, 1.2, 0.8)
    assert tuple(transition.within_fan for transition in path.transitions) == (True, True, False, False)
    assert tuple(float(transition.boundary_progress) for transition in path.transitions[:2]) == (0.0, 0.0)
    assert all(transition.boundary_progress > 0.0 for transition in path.transitions[2:])


def test_rejects_endpoint_aliases_and_missing_reached_runs() -> None:
    endpoint_aliases = (
        _candidate(0, 0, Fraction(1), (4.0, 1.0)),
        _candidate(1, 1, Fraction(0), (3.0, 1.0)),
    )
    with pytest.raises(AmbiguousFigure5BoundaryPathError, match="endpoint"):
        _build(endpoint_aliases)

    candidates = (_candidate(0, 0, Fraction(1, 4), (1.0, 1.0)),)
    with pytest.raises(MissingFigure5GuideRunError, match="reached guide runs"):
        build_figure5_boundary_path(
            inward_components=(SQUARE,),
            side_candidates=candidates,
            reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1))),
            published_start_center=candidates[0].center,
            published_start_radius=Millimetre(TOOL_RADIUS.value),
            tool_radius=TOOL_RADIUS,
            boundary_evidence_bound=Millimetre(0.01),
            start_evidence_bound=Millimetre(0.01),
        )


def test_published_start_matches_tool_center_q_and_radius_one_evidence() -> None:
    matching_circle = _candidate(0, 0, Fraction(1, 4), (3.0, 3.0), 0.4)
    center_only_match = _candidate(1, 1, Fraction(1, 4), (1.0, 0.0), 2.0)

    path = build_figure5_boundary_path(
        inward_components=(SQUARE,),
        side_candidates=(matching_circle, center_only_match),
        reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1))),
        published_start_center=_point(1.0, 0.0),
        published_start_radius=Millimetre(1.0),
        tool_radius=TOOL_RADIUS,
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
    )

    assert path.circles[0].run_id == GuideRunId(0)

    with pytest.raises(AmbiguousFigure5BoundaryPathError, match="published start radius"):
        build_figure5_boundary_path(
            inward_components=(SQUARE,),
            side_candidates=(matching_circle,),
            reached_run_ids=frozenset((GuideRunId(0),)),
            published_start_center=_point(1.0, 0.0),
            published_start_radius=Millimetre(1.1),
            tool_radius=TOOL_RADIUS,
            boundary_evidence_bound=Millimetre(0.01),
            start_evidence_bound=Millimetre(0.01),
        )


def test_published_start_selects_exact_nearest_site_among_several_inside_bound() -> None:
    nearest = _candidate(0, 0, Fraction(1, 4), (1.0, 1.0))
    farther = _candidate(1, 0, Fraction(63, 250), (1.008, 1.0))

    path = build_figure5_boundary_path(
        inward_components=(SQUARE,),
        side_candidates=(farther, nearest),
        reached_run_ids=frozenset((GuideRunId(0), GuideRunId(1))),
        published_start_center=_point(1.001, 0.0),
        published_start_radius=Millimetre(1.0),
        tool_radius=TOOL_RADIUS,
        boundary_evidence_bound=Millimetre(0.01),
        start_evidence_bound=Millimetre(0.01),
    )

    assert path.circles[0].run_id == GuideRunId(0)

    equidistant = (
        _candidate(2, 0, Fraction(3, 16), (0.75, 1.0)),
        _candidate(3, 0, Fraction(5, 16), (1.25, 1.0)),
    )
    with pytest.raises(AmbiguousFigure5BoundaryPathError, match="unique nearest"):
        build_figure5_boundary_path(
            inward_components=(SQUARE,),
            side_candidates=equidistant,
            reached_run_ids=frozenset((GuideRunId(2), GuideRunId(3))),
            published_start_center=_point(1.0, 0.0),
            published_start_radius=Millimetre(1.0),
            tool_radius=TOOL_RADIUS,
            boundary_evidence_bound=Millimetre(0.01),
            start_evidence_bound=Millimetre(1.0),
        )


def test_rejects_more_than_one_inward_offset_component() -> None:
    candidates = (_candidate(0, 0, Fraction(1, 2), (2.0, 1.0)),)
    with pytest.raises(DisconnectedFigure5BoundaryPathError, match="one connected"):
        build_figure5_boundary_path(
            inward_components=(SQUARE, tuple(reversed(SQUARE))),
            side_candidates=candidates,
            reached_run_ids=frozenset((GuideRunId(0),)),
            published_start_center=candidates[0].center,
            published_start_radius=Millimetre(TOOL_RADIUS.value),
            tool_radius=TOOL_RADIUS,
            boundary_evidence_bound=Millimetre(0.01),
            start_evidence_bound=Millimetre(0.01),
        )


def test_actual_figure5_endpoint_hypothesis_resolves_locally_into_bounded_offset_motion() -> None:
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)
    hypothesis = next(
        candidate
        for candidate in build_distance_admissible_hypotheses(guide, guide.runs[2].stations[14])
        if candidate.projected_boundary_site.segment_id == ProjectionBoundarySegmentId(63)
    )
    assert hypothesis.projected_boundary_site.vertex_id is not None
    assert case.start_marker is not None
    assert case.start_marker_radius is not None

    path = build_figure5_boundary_path(
        inward_components=figure7_inward_offset(case).components,
        side_candidates=(hypothesis,),
        reached_run_ids=frozenset((GuideRunId(2),)),
        published_start_center=case.start_marker,
        published_start_radius=case.start_marker_radius,
        tool_radius=case.tool_radius,
        boundary_evidence_bound=guide.site_budget.admissible_distance_gap,
        start_evidence_bound=guide.site_budget.admissible_distance_gap,
    )

    assert path.circles[0].boundary_site.vertex_id is None
    assert path.circles[0].clockwise is False
    assert path.transitions[0].boundary_progress > 0.0
    assert path.transitions[0].samples[0] == path.transitions[0].samples[-1]
