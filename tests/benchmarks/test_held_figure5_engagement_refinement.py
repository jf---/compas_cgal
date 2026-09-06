import math
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementResolutionError
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
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from compas_cgal import _coverage_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


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
        assert float(maximum_predecessor_engagement(_paper_candidate(previous), _paper_candidate(circle), TOOL)) <= float(CAP.theta)


def test_unresolved_gap_fails_instead_of_forcing_successor() -> None:
    with pytest.raises(Figure5EngagementResolutionError, match="source interval 0"):
        refine_figure5_engagement((_circle(2, 0), _circle(8, 1)), BOUNDARY, TOOL, CAP, max_depth=0)


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
        assert float(maximum_predecessor_engagement(_paper_candidate(previous), _paper_candidate(current), TOOL)) <= float(CAP.theta)


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
        assert float(maximum_predecessor_engagement(_paper_candidate(previous), _paper_candidate(circle), TOOL)) <= float(CAP.theta)


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
    lengths = _side_lengths(boundary)
    transitions = tuple(_transition(boundary, lengths, a, b) for a, b in zip(result.circles, result.circles[1:]))
    assert len(transitions) == len(result.circles) - 1
    for previous, circle, transition in zip(result.circles[:-1], result.circles[1:], transitions, strict=True):
        assert float(maximum_predecessor_engagement(_paper_candidate(previous), _paper_candidate(circle), case.tool_radius)) <= float(CAP.theta)
        assert transition.samples[0] == previous.contact_point
        assert transition.samples[-1] == circle.contact_point
