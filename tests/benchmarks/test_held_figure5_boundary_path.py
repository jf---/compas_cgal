from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_boundary_path import AmbiguousFigure5BoundaryPathError
from benchmarks.held_figure5_boundary_path import DisconnectedFigure5BoundaryPathError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_boundary_path import MissingFigure5GuideRunError
from benchmarks.held_figure5_boundary_path import _project_to_offset_side
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
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


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


def _build(candidates: tuple[ProjectionAdmissibleBoundaryHypothesis, ...], *, start_index: int = 0) -> Figure5BoundaryPath:
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


# --- Closing seam connector -------------------------------------------------

SQUARE_PERIMETER_MM = 16.0


def test_closing_connector_runs_ccw_from_last_contact_to_first_along_the_loop() -> None:
    candidates = (
        _candidate(0, 0, Fraction(1, 4), (1.0, 1.0)),
        _candidate(1, 1, Fraction(1, 2), (3.0, 2.0)),
        _candidate(2, 3, Fraction(3, 4), (1.0, 1.0)),
    )
    circles = _build(candidates).circles

    connectors = Figure5BoundaryConnectors.build(SQUARE, circles)

    assert connectors.seam is Figure5SeamClosure.CONNECTOR
    assert len(connectors.transitions) == len(circles)
    closing = connectors.transitions[-1]
    assert closing.samples == (circles[-1].contact_point, _point(0.0, 0.0), circles[0].contact_point)
    assert closing.within_fan is False
    assert float(closing.boundary_progress) == 2.0
    assert connectors.seam_progress == closing.boundary_progress
    assert sum(float(transition.boundary_progress) for transition in connectors.transitions) == SQUARE_PERIMETER_MM


def test_coincident_seam_emits_no_zero_length_motion() -> None:
    candidates = (
        _candidate(0, 0, Fraction(1, 4), (1.0, 0.8), 0.8, station=0),
        _candidate(0, 0, Fraction(1, 4), (1.0, 1.2), 1.2, station=1),
        _candidate(1, 1, Fraction(1, 2), (3.0, 2.0)),
    )
    fan_first, fan_second, far = _build(candidates).circles
    assert fan_first.contact_point == fan_second.contact_point

    connectors = Figure5BoundaryConnectors.build(SQUARE, (fan_second, far, fan_first))

    assert connectors.seam is Figure5SeamClosure.COINCIDENT
    assert len(connectors.transitions) == 2
    assert float(connectors.seam_progress) == 0.0
    assert connectors.transitions[-1].samples[-1] == fan_first.contact_point


def test_single_circle_closes_with_one_full_lap() -> None:
    circles = _build((_candidate(0, 0, Fraction(1, 2), (2.0, 1.0)),)).circles

    connectors = Figure5BoundaryConnectors.build(SQUARE, circles)

    assert connectors.seam is Figure5SeamClosure.CONNECTOR
    assert len(connectors.transitions) == 1
    assert float(connectors.seam_progress) == SQUARE_PERIMETER_MM


def test_connectors_reject_an_empty_circle_sequence() -> None:
    with pytest.raises(DisconnectedFigure5BoundaryPathError, match="selected candidates"):
        Figure5BoundaryConnectors.build(SQUARE, ())


def test_boundary_path_transitions_are_the_closed_connectors() -> None:
    candidates = (
        _candidate(0, 0, Fraction(1, 4), (1.0, 1.0)),
        _candidate(1, 2, Fraction(1, 2), (2.0, 3.0)),
    )
    path = _build(candidates)

    assert path.transitions == Figure5BoundaryConnectors.build(SQUARE, path.circles).transitions


_SITES = st.sets(st.tuples(st.integers(0, 3), st.integers(1, 7)), min_size=1, max_size=12)


@given(sites=_SITES, data=st.data())
def test_closed_connectors_chain_contacts_and_sum_to_one_loop(sites: set[tuple[int, int]], data: st.DataObject) -> None:
    ordered = sorted(sites)
    candidates = tuple(_candidate(index, side, Fraction(eighth, 8), (2.0, 2.0)) for index, (side, eighth) in enumerate(ordered))
    start = data.draw(st.integers(0, len(candidates) - 1))
    circles = _build(candidates, start_index=start).circles

    connectors = Figure5BoundaryConnectors.build(SQUARE, circles)

    assert connectors.seam is Figure5SeamClosure.CONNECTOR
    assert len(connectors.transitions) == len(circles)
    for index, transition in enumerate(connectors.transitions):
        assert transition.samples[0] == circles[index].contact_point
        assert transition.samples[-1] == circles[(index + 1) % len(circles)].contact_point
        assert float(transition.boundary_progress) > 0.0
    # Dyadic sites on an integer square: every progress is exact in binary64.
    assert sum(float(transition.boundary_progress) for transition in connectors.transitions) == SQUARE_PERIMETER_MM


def _oracle_project_to_offset_side(
    point: tuple[Fraction, Fraction],
    component: tuple[Point2[WorldXY], ...],
    side_index: int,
) -> tuple[Fraction, Fraction]:
    """Verbatim pre-L10 ``_project_to_offset_side``: the Fraction clamp-and-project."""
    start = (Fraction.from_float(float(component[side_index].x)), Fraction.from_float(float(component[side_index].y)))
    end_point = component[(side_index + 1) % len(component)]
    end = (Fraction.from_float(float(end_point.x)), Fraction.from_float(float(end_point.y)))
    edge_x = end[0] - start[0]
    edge_y = end[1] - start[1]
    squared_length = edge_x * edge_x + edge_y * edge_y
    parameter = min(Fraction(1), max(Fraction(0), ((point[0] - start[0]) * edge_x + (point[1] - start[1]) * edge_y) / squared_length))
    delta_x = point[0] - (start[0] + parameter * edge_x)
    delta_y = point[1] - (start[1] + parameter * edge_y)
    return parameter, delta_x * delta_x + delta_y * delta_y


_PROJECTION_LATTICE = st.integers(min_value=-16, max_value=16).map(lambda value: value / 8)
_PROJECTION_COORDINATE = st.one_of(_PROJECTION_LATTICE, st.floats(min_value=-1e6, max_value=1e6, allow_nan=False, allow_infinity=False))


@settings(max_examples=1500, deadline=None)
@given(
    point=st.tuples(_PROJECTION_COORDINATE, _PROJECTION_COORDINATE),
    corners=st.lists(st.tuples(_PROJECTION_COORDINATE, _PROJECTION_COORDINATE), min_size=2, max_size=4, unique=True),
    side=st.integers(min_value=0, max_value=3),
)
def test_offset_side_projection_is_the_exact_fraction_projection(
    point: tuple[float, float],
    corners: list[tuple[float, float]],
    side: int,
) -> None:
    component = tuple(Point2[WorldXY].build(x, y) for x, y in corners)
    side_index = side % len(component)
    start, end = component[side_index], component[(side_index + 1) % len(component)]
    if start == end:
        return
    exact_point = (Fraction.from_float(point[0]), Fraction.from_float(point[1]))
    parameter, squared_distance = _project_to_offset_side(exact_point, component, side_index)
    expected_parameter, expected_squared_distance = _oracle_project_to_offset_side(exact_point, component, side_index)
    assert type(parameter) is Fraction and type(squared_distance) is Fraction
    assert (parameter, squared_distance) == (expected_parameter, expected_squared_distance)


def test_offset_side_projection_refuses_a_zero_length_side() -> None:
    corner = Point2[WorldXY].build(1.0, 2.0)
    with pytest.raises(DisconnectedFigure5BoundaryPathError):
        _project_to_offset_side((Fraction(0), Fraction(0)), (corner, corner), 0)
