import math
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from benchmarks.held_boundary_arclength import BoundaryArclength as _BoundaryArclength
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks import held_figure5_engagement_refinement as engagement_refinement
from benchmarks.held_connector_engagement import ConnectorEngagementLemma
from benchmarks.held_figure5_engagement_refinement import ConvexCornerEngagementError
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementResolutionError
from benchmarks.held_figure5_engagement_refinement import exterior_angle_exceeds_cap
from benchmarks.held_figure5_engagement_refinement import midpoint
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal import _coverage_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.benchmarks.held_engagement_reference_oracle import reference_engagement


BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
TOOL = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80))


def _circle(x: int, run: int) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x, 10), len(BOUNDARY))
    return Figure5CounterclockwiseCircle(
        GuideRunId(run),
        GuideRunStationOrdinal(0),
        site,
        site,
        Point2[WorldXY].build(x, 1),
        Point2[WorldXY].build(x, 0),
        GuideRadius.build(1),
    )


def test_refinement_preserves_sources_and_caps_every_final_pair() -> None:
    original = (_circle(2, 0), _circle(8, 1))
    result = refine_figure5_engagement(original, BOUNDARY, TOOL, CAP)
    assert len(result.circles) > len(original)
    assert tuple(result.circles[index] for index in result.original_indices) == original
    domain = _coverage_2.CutterCentreDomain2.build(np.array([(float(p.x), float(p.y), 0) for p in BOUNDARY]), [], 1)
    for circle in result.circles:
        assert domain.contains(float(circle.center.x), float(circle.center.y))
    for previous, circle in zip(result.circles, result.circles[1:]):
        assert reference_engagement(_paper_candidate(previous), _paper_candidate(circle), TOOL).within(CAP)


def test_unresolved_gap_fails_instead_of_forcing_successor() -> None:
    with pytest.raises(Figure5EngagementResolutionError, match="source interval 0"):
        refine_figure5_engagement((_circle(2, 0), _circle(8, 1)), BOUNDARY, TOOL, CAP, max_depth=0)


def _rectangle(width: Fraction, height: Fraction) -> tuple[Point2[WorldXY], ...]:
    return tuple(Point2[WorldXY].build(float(x), float(y)) for x, y in ((0, 0), (width, 0), (width, height), (0, height)))


def _station(boundary: tuple[Point2[WorldXY], ...], side: int, parameter: Fraction, ordinal: int, radius: Fraction) -> Figure5CounterclockwiseCircle:
    """A source circle of ``radius`` tangent to ``side`` at the exact side ``parameter``."""
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, len(boundary))
    contact = _canonical_offset_point(boundary, site)
    start, end = boundary[side], boundary[(side + 1) % len(boundary)]
    edge_x, edge_y = float(end.x) - float(start.x), float(end.y) - float(start.y)
    length = math.hypot(edge_x, edge_y)
    center = Point2[WorldXY].build(float(contact.x) - float(radius) * edge_y / length, float(contact.y) + float(radius) * edge_x / length)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(ordinal), site, site, center, contact, GuideRadius.build(float(radius)))


# Convex vertex (width, 0) with the source contacts the same exact arclength u
# before and after it: the arclength half is EXACTLY the vertex, so the third
# point (2u/3 along side 0) is the midpoint. The binary64 bookkeeping misses
# the junction and lands one rounding off the apex: (10, 4.4e-16) absorbed by
# the centre (10x10), or a 1.1e-16 radius the centre guard lets through (1x3,
# either side of the vertex).
EXACT_JUNCTIONS = (
    (Fraction(10), Fraction(10), Fraction(2)),
    (Fraction(1), Fraction(3), Fraction(21, 40)),
    (Fraction(1), Fraction(3), Fraction(3, 5)),
)


@pytest.mark.parametrize(("width", "height", "u"), EXACT_JUNCTIONS)
def test_midpoint_of_an_exact_vertex_half_takes_the_exact_third_point(width: Fraction, height: Fraction, u: Fraction) -> None:
    boundary = _rectangle(width, height)
    previous = _station(boundary, 0, 1 - u / width, 0, Fraction(1))
    successor = _station(boundary, 1, u / height, 1, Fraction(1))
    halfway = midpoint(previous, successor, boundary)
    third = 1 - u / width + (2 * u / 3) / width
    assert int(halfway.boundary_site.segment_id) == 0
    assert halfway.boundary_site.parameter == Fraction(float(third))
    assert halfway.contact_point == _canonical_offset_point(boundary, halfway.boundary_site)
    # The corner limit is the distance u/3 to side 1 (right angle, n·m = 0).
    assert float(halfway.radius.value) == pytest.approx(float(u / 3))


def test_exact_vertex_half_refines_across_the_convex_corner() -> None:
    """User ruling 2026-10-22: a right-angle convex corner under an 80-degree cap is refused by name.

    Lemma (docs/held_connector_engagement.md, "Convex corners"): the first
    arrival at a convex vertex of interior angle alpha engages at least
    pi - alpha = 90 > 80 degrees whatever was cut before, so no refinement
    across it can be certified.
    """
    boundary = _rectangle(Fraction(10), Fraction(10))
    first = _station(boundary, 0, Fraction(4, 5), 0, Fraction(1))
    second = _station(boundary, 1, Fraction(1, 5), 1, Fraction(1))
    for corner_approaches in (False, True):
        with pytest.raises(ConvexCornerEngagementError) as caught:
            refine_figure5_engagement((first, second), boundary, TOOL, CAP, corner_approaches=corner_approaches)
        assert (caught.value.vertex, caught.value.connector, caught.value.seam) == (1, 0, False)


def test_convex_corner_over_the_cap_is_refused_before_any_midpoint(monkeypatch: pytest.MonkeyPatch) -> None:
    """The refusal is decided from the vertex and the cuts so far, not by halving circles toward zero radius."""

    def no_midpoint(*args: object, **kwargs: object) -> Figure5CounterclockwiseCircle:
        raise AssertionError("midpoint called before the convex-corner refusal")

    monkeypatch.setattr(engagement_refinement, "midpoint", no_midpoint)
    boundary = _rectangle(Fraction(10), Fraction(10))
    first = _station(boundary, 0, Fraction(4, 5), 0, Fraction(1))
    second = _station(boundary, 1, Fraction(1, 5), 1, Fraction(1))
    with pytest.raises(ConvexCornerEngagementError) as caught:
        refine_figure5_engagement((first, second), boundary, TOOL, CAP)
    assert caught.value.vertex == 1


def test_connector_over_two_sharp_corners_is_refused_by_name_at_the_first() -> None:
    """Side 0 to side 2: midpoints split the connector until one sub-connector is the first arrival at (10, 0)."""
    boundary = _rectangle(Fraction(10), Fraction(10))
    first = _station(boundary, 0, Fraction(4, 5), 0, Fraction(1))
    second = _station(boundary, 2, Fraction(1, 5), 1, Fraction(1))
    with pytest.raises(ConvexCornerEngagementError) as caught:
        refine_figure5_engagement((first, second), boundary, TOOL, CAP)
    assert caught.value.vertex == 1 and not caught.value.seam


def test_seam_over_a_sharp_convex_corner_is_refused_by_name() -> None:
    boundary = _rectangle(Fraction(10), Fraction(10))
    first, second = _station(boundary, 0, Fraction(1, 5), 0, Fraction(1)), _station(boundary, 0, Fraction(4, 5), 1, Fraction(1))
    with pytest.raises(ConvexCornerEngagementError) as caught:
        refine_figure5_engagement((first, second), boundary, TOOL, CAP, close_seam=True)
    assert caught.value.seam and caught.value.vertex == 1


@pytest.mark.parametrize(("degrees", "sharp"), ((80, True), (89, True), (91, False), (100, False)))
def test_convex_corner_decision_compares_the_exterior_angle_with_the_cap_exactly(degrees: int, sharp: bool) -> None:
    boundary = _rectangle(Fraction(10), Fraction(10))
    cap = EngagementCap.build(math.radians(degrees))
    lemma = ConnectorEngagementLemma.build(TOOL, cap)
    assert all(exterior_angle_exceeds_cap(boundary, vertex, lemma.cos_sq_half_cap) is sharp for vertex in range(4))


def test_reflex_and_obtuse_vertices_are_never_sharp() -> None:
    # L-shape: one reflex vertex (3, 3) and right angles; a 135-degree octagon corner has exterior 45 < 80.
    l_shape = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (6, 0), (6, 3), (3, 3), (3, 6), (0, 6)))
    octagon = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (8, 0), (12, 4), (12, 12), (8, 16), (0, 16), (-4, 12), (-4, 4)))
    c = ConnectorEngagementLemma.build(TOOL, CAP).cos_sq_half_cap
    assert not exterior_angle_exceeds_cap(l_shape, 3, c)
    assert not any(exterior_angle_exceeds_cap(octagon, vertex, c) for vertex in range(len(octagon)))


def test_midpoint_exactly_at_a_convex_apex_refuses() -> None:
    # 4x1: arclength 2 | 1 | 3. The half (3) is exactly vertex (4, 1), the third
    # (2) exactly vertex (4, 0): the midpoint contact is a convex apex, where no
    # positive circle inside the wedge touches the boundary. Decided exactly.
    boundary = _rectangle(Fraction(4), Fraction(1))
    previous = _station(boundary, 0, Fraction(1, 2), 0, Fraction(1, 4))
    successor = _station(boundary, 2, Fraction(3, 4), 1, Fraction(1, 4))
    with pytest.raises(Figure5EngagementResolutionError, match="convex corner apex"):
        midpoint(previous, successor, boundary)


def test_genuinely_sub_ulp_corner_circle_is_refused_as_unrepresentable() -> None:
    # The exact half lies 5e-18 past vertex (10, 0) on side 1: a positive
    # circle exists (radius 5e-18), but its centre is not representable apart
    # from its contact in binary64. This keeps the representability guard
    # reachable after the exact-arclength fix.
    boundary = _rectangle(Fraction(10), Fraction(10))
    previous = _station(boundary, 0, Fraction(4, 5), 0, Fraction(1))
    successor = _station(boundary, 1, Fraction(1, 5) + Fraction(1, 10**18), 1, Fraction(1))
    with pytest.raises(Figure5EngagementResolutionError, match="representable"):
        midpoint(previous, successor, boundary)


def test_exact_junction_with_irrational_sides_in_rational_ratio() -> None:
    # Sides of length sqrt(2) and 3*sqrt(2) meeting at convex vertex (1, 1):
    # (1 - p0)*sqrt(2) == p1*3*sqrt(2) is an exact junction only the folded
    # radical arithmetic sees; the midpoint is the third point on side 0.
    # Binary64 lands past the apex and refuses "no positive circle".
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (1, 1), (-2, 4), (-3, 3)))
    p0, p1 = Fraction(1, 5), Fraction(4, 15)
    previous = _station(boundary, 0, p0, 0, Fraction(1, 10))
    successor = _station(boundary, 1, p1, 1, Fraction(1, 10))
    halfway = midpoint(previous, successor, boundary)
    assert int(halfway.boundary_site.segment_id) == 0
    assert halfway.boundary_site.parameter == Fraction(float(p0 + (2 * (1 - p0) / 3)))


@pytest.mark.parametrize(
    ("terms", "sign"),
    [
        (((Fraction(2), Fraction(3)), (Fraction(18), Fraction(-1))), 0),
        (((Fraction(2), Fraction(1)), (Fraction(3), Fraction(1)), (Fraction(49, 5), Fraction(-1))), 1),
        # sqrt(2) - 1.414213562373095049 = -1.98e-19: binary64 evaluates it as 0.
        (((Fraction(2), Fraction(1)), (Fraction(1), Fraction(-1414213562373095049, 10**18))), -1),
        (((Fraction(9, 4), Fraction(2)), (Fraction(1), Fraction(-3))), 0),
    ],
)
def test_boundary_arclength_sign_is_exact(terms: tuple[tuple[Fraction, Fraction], ...], sign: int) -> None:
    total = _BoundaryArclength.zero()
    for radicand, coefficient in terms:
        total = total + _BoundaryArclength.root(radicand, coefficient)
    assert total.sign() == sign


def test_concave_corner_approach_rotates_only_at_shared_contact() -> None:
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (8, 0), (8, 8), (4, 8), (4, 4), (0, 4)))
    incoming = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(3), Fraction(1, 2), len(boundary))
    outgoing = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(4), Fraction(1, 2), len(boundary))
    first = Figure5CounterclockwiseCircle(
        GuideRunId(0), GuideRunStationOrdinal(0), incoming, incoming, Point2[WorldXY].build(5, 6), Point2[WorldXY].build(4, 6), GuideRadius.build(1)
    )
    last = Figure5CounterclockwiseCircle(
        GuideRunId(1), GuideRunStationOrdinal(0), outgoing, outgoing, Point2[WorldXY].build(2, 3), Point2[WorldXY].build(2, 4), GuideRadius.build(1)
    )
    result = refine_figure5_engagement((first, last), boundary, TOOL, CAP, corner_approaches=True)
    assert result.circles[0] == first and result.circles[-1] == last
    corner_circles = tuple(circle for circle in result.circles if circle.contact_point == boundary[4])
    assert any(circle.center == Point2[WorldXY].build(5, 4) for circle in corner_circles)
    assert any(circle.center == Point2[WorldXY].build(4, 3) for circle in corner_circles)
    for previous, current in zip(result.circles, result.circles[1:]):
        assert reference_engagement(_paper_candidate(previous), _paper_candidate(current), TOOL).within(CAP)


def test_source_circle_is_limited_by_incident_corner_without_losing_contact() -> None:
    source = replace(_circle(9, 0), center=Point2[WorldXY].build(9, 2), radius=GuideRadius.build(2))
    result = refine_figure5_engagement((source,), BOUNDARY, TOOL, CAP)
    circle = result.circles[0]
    assert float(circle.radius.value) == 1
    assert circle.contact_point == source.contact_point
    assert circle.run_id == source.run_id
    assert circle.center == Point2[WorldXY].build(9, 1)


def test_concave_corner_uses_shared_contact_fan() -> None:
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)))
    circles = []
    for run, side, parameter, center, contact in (
        (0, 2, Fraction(119, 120), (4.05, 3.0), (4.05, 4.0)),
        (1, 3, Fraction(1, 120), (3.0, 4.05), (4.0, 4.05)),
    ):
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, len(boundary))
        circles.append(
            Figure5CounterclockwiseCircle(
                GuideRunId(run), GuideRunStationOrdinal(0), site, site, Point2[WorldXY].build(*center), Point2[WorldXY].build(*contact), GuideRadius.build(1)
            )
        )
    result = refine_figure5_engagement(tuple(circles), boundary, TOOL, CAP)
    assert any(circle.contact_point == boundary[3] for circle in result.circles)
    for previous, circle in zip(result.circles, result.circles[1:]):
        assert reference_engagement(_paper_candidate(previous), _paper_candidate(circle), TOOL).within(CAP)


def test_full_figure5_refinement_preserves_families_and_connects_capped_motion() -> None:
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)
    boundary = figure7_inward_offset(case).components[0]
    source = build_figure5_approximate_path(guide, (boundary,)).path.circles
    result = refine_figure5_engagement(source, boundary, case.tool_radius, CAP)
    assert len(result.original_indices) == len(source)
    for original, index in zip(source, result.original_indices, strict=True):
        repaired = result.circles[index]
        assert (repaired.run_id, repaired.station_ordinal, repaired.contact_point) == (original.run_id, original.station_ordinal, original.contact_point)
    # Held paths close: a CCW seam connector from the last contact back to the
    # first, unless those contacts coincide exactly (user ruling, e05674d3).
    connectors = Figure5BoundaryConnectors.build(boundary, result.circles)
    closed = len(result.circles) - (connectors.seam is Figure5SeamClosure.COINCIDENT)
    assert len(connectors.transitions) == closed
    successors = (*result.circles[1:], result.circles[0])
    for previous, circle, transition in zip(result.circles, successors, connectors.transitions):
        assert transition.samples[0] == previous.contact_point
        assert transition.samples[-1] == circle.contact_point
    for previous, circle in zip(result.circles, result.circles[1:]):
        assert reference_engagement(_paper_candidate(previous), _paper_candidate(circle), case.tool_radius).within(CAP)
