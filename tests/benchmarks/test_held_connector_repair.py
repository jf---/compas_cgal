"""Connector acceptance in the engagement refinement: a transition is emitted only with a certified connector.

User ruling 2026-10-21: connectors that tiers 1 + 1b leave without a
certificate are repaired by refinement (midpoint circles), and the seam
connector by new last circles (circles[0] is never re-cut). The post-gate
``ConnectorEngagementCertificate`` must then pass on every refined path.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import ConnectorEngagementCertificate
from benchmarks.held_connector_engagement import ConnectorEngagementLemma
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_engagement_refinement import ConnectorEngagementResolutionError
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementResolutionError
from benchmarks.held_figure5_engagement_refinement import circle_on_side
from benchmarks.held_figure5_engagement_refinement import midpoint
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.benchmarks.held_part_fixture import NO_MATERIAL
from tests.benchmarks.held_part_fixture import part_around

TOOL = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80))
SQUARE = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
# Obtuse (135 degree) convex corners: the corner limit leaves room for rho = 2
# near a vertex. At the square's right angles it shrinks midpoint circles
# towards the vertex, where the single-annulus lemma cannot certify a move
# away from a vanishing circle, so repair there exhausts its depth.
OCTAGON = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (8, 0), (12, 4), (12, 12), (8, 16), (0, 16), (-4, 12), (-4, 4)))
# A tooth: side 2 dives through the 1.8 guide circle tangent at (0, 0), past
# its centre (0, 1.8). Side 0's convex corner is gentle (14 degrees), so the
# corner limit (0.25 cot 7 degrees, about 2.0) keeps rho = 1.8.
TOOTH = tuple(Point2[WorldXY].build(x, y) for x, y in ((-3.0, 0.0), (0.25, 0.0), (0.5, 0.0625), (-0.5, 3.375), (-3.0, 3.375)))


def _on(boundary: tuple[Point2[WorldXY], ...], side: int, parameter: Fraction, rho: float, run: int) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, len(boundary))
    contact = _canonical_offset_point(boundary, site)
    source = Figure5CounterclockwiseCircle(GuideRunId(run), GuideRunStationOrdinal(0), site, site, contact, contact, GuideRadius.build(rho))
    return circle_on_side(source, site, contact, GuideRadius.build(rho), boundary)


def _open_assessment(boundary: tuple[Point2[WorldXY], ...], circles: tuple[Figure5CounterclockwiseCircle, ...]) -> ConnectorEngagementAssessment:
    """Tiers 1 + 1b on the open connectors only (the seam dropped)."""
    built = Figure5BoundaryConnectors.build(boundary, circles)
    opened = Figure5BoundaryConnectors(built.transitions[: len(circles) - 1], Figure5SeamClosure.COINCIDENT, built.seam_progress)
    return ConnectorEngagementAssessment.build(circles, opened, TOOL, CAP)


def _pair_within(previous: Figure5CounterclockwiseCircle, current: Figure5CounterclockwiseCircle) -> bool:
    return exact_predecessor_engagement(_paper_candidate(previous), _paper_candidate(current), TOOL).within(CAP)


def _convex_corner_pair() -> tuple[Figure5CounterclockwiseCircle, Figure5CounterclockwiseCircle]:
    """Circles either side of the octagon's convex corner (8, 0): the pair is within the cap, the connector is not proved."""
    return _on(OCTAGON, 0, Fraction(7, 8), 2.0, 0), _on(OCTAGON, 1, Fraction(1, 32), 2.0, 1)


def _dstar_pair() -> tuple[Figure5CounterclockwiseCircle, Figure5CounterclockwiseCircle]:
    """The tooth: the pair is within the cap; the connector is refused only at the interior hole minimum d*."""
    return _on(TOOTH, 0, Fraction(12, 13), 1.8, 0), _on(TOOTH, 2, Fraction(7, 8), 0.5, 1)


def test_convex_corner_fixture_is_a_connector_only_refusal() -> None:
    first, second = _convex_corner_pair()
    assert _pair_within(first, second)
    assert _open_assessment(OCTAGON, (first, second)).refused == (0,)


def test_refused_connector_forces_a_midpoint_and_the_post_gate_passes() -> None:
    first, second = _convex_corner_pair()
    result = refine_figure5_engagement((first, second), OCTAGON, TOOL, CAP, part=NO_MATERIAL)
    assert len(result.circles) == 3, "one midpoint circle repairs the connector"
    assert (result.circles[0].contact_point, result.circles[-1].contact_point) == (first.contact_point, second.contact_point)
    assert result.original_indices == (0, 2)
    assert not _open_assessment(OCTAGON, result.circles).refused
    for previous, current in zip(result.circles, result.circles[1:]):
        assert _pair_within(previous, current)


def test_dstar_fixture_is_refused_only_at_the_interior_hole_minimum(monkeypatch: pytest.MonkeyPatch) -> None:
    first, second = _dstar_pair()
    assert _pair_within(first, second)
    assert _open_assessment(TOOTH, (first, second)).refused == (0,)

    def endpoints_only(self: ConnectorEngagementLemma, d_min_sq: Fraction, d_max_sq: Fraction, inner: Fraction) -> bool:
        return self.hole_within(d_min_sq, inner) and self.hole_within(d_max_sq, inner)

    monkeypatch.setattr(ConnectorEngagementLemma, "hole_interval_within", endpoints_only)
    assert not _open_assessment(TOOTH, (first, second)).refused, "only the d* branch refuses the fixture"


def test_dstar_refusal_forces_a_midpoint() -> None:
    first, second = _dstar_pair()
    result = refine_figure5_engagement((first, second), TOOTH, TOOL, CAP, part=part_around(TOOTH, TOOL))
    assert len(result.circles) > 2
    assert (result.circles[0].contact_point, result.circles[-1].contact_point) == (first.contact_point, second.contact_point)
    assert not _open_assessment(TOOTH, result.circles).refused


def test_exhausted_connector_repair_raises_the_named_error() -> None:
    first, second = _convex_corner_pair()
    with pytest.raises(ConnectorEngagementResolutionError) as caught:
        refine_figure5_engagement((first, second), OCTAGON, TOOL, CAP, part=NO_MATERIAL, max_depth=0)
    assert caught.value.connector == 0
    assert caught.value.seam is False
    assert isinstance(caught.value, Figure5EngagementResolutionError), "callers catching the refinement's refusals still see it"


def test_right_angle_corner_with_vanishing_circles_exhausts_by_name() -> None:
    """Pinned: a right angle under an 80-degree cap is refused by name (ConvexCornerEngagementError, a subclass), never a pass."""
    first, second = _on(SQUARE, 0, Fraction(9, 10), 1.0, 0), _on(SQUARE, 1, Fraction(1, 20), 1.0, 1)
    assert _pair_within(first, second)
    with pytest.raises(ConnectorEngagementResolutionError):
        refine_figure5_engagement((first, second), SQUARE, TOOL, CAP, part=part_around(SQUARE, TOOL))


# Chords per quarter arc of the oracle's rounded square. At (10, 0) the sagitta
# 1 - cos(pi / 4096) = 2.9e-7 mm stays below every station distance the
# assertion reads (>= 2.4e-4 mm before V).
CORNER_CHORDS_AT_V = 1024
CORNER_CHORDS_ELSEWHERE = 16
# Corner-limited circles at (10 - 2^-k, 0), rho = 2^-k (tangent to both sides), k = 0..4.
CORNER_REFINEMENT_LEVELS = 4


def _inscribed_rounded_square() -> np.ndarray:
    """A polygon inside SQUARE (+) D_r: each rounded corner by chords with vertices ON its arc."""
    vertices = []
    for (vx, vy), start in (((10, 0), -90), ((10, 10), 0), ((0, 10), 90), ((0, 0), 180)):
        count = CORNER_CHORDS_AT_V if (vx, vy) == (10, 0) else CORNER_CHORDS_ELSEWHERE
        for j in range(count + 1):
            angle = math.radians(start + 90 * j / count)
            vertices.append((vx + math.cos(angle), vy + math.sin(angle), 0.0))
    return np.array(vertices, dtype=np.float64)


def test_right_angle_corner_refusal_is_a_true_engagement_above_the_cap() -> None:
    """The refusal above is not lemma looseness: the actual stock exceeds the cap on the approach to the vertex.

    Lemma (docs/held_connector_engagement.md, "Convex corners"): at a convex
    offset vertex of interior angle alpha, the first arrival of the tool centre
    engages, in the limit, at least pi - alpha of rim, whatever guide circles
    were cut before (all lie inside the centre domain, so none reaches the rim
    points in the vertex's normal cone). Here pi - alpha = 90 > 80 degrees.
    Oracle: Stock_OVER (contained in the true stock) on a polygon contained in
    SQUARE (+) D_r, so every reading is a lower bound on the true engagement.
    Cuts in emission order: (8, 0) rho 1, then the corner-limited circles down
    to rho 1/16 with their connectors; then the final approach to V = (10, 0).
    """
    r = float(TOOL.value)
    stock = _stock_2.Stock2(_inscribed_rounded_square(), [])
    previous = None
    for x, rho in [(8.0, 1.0)] + [(10 - 2.0**-k, 2.0**-k) for k in range(CORNER_REFINEMENT_LEVELS + 1)]:
        if previous is not None:
            stock.subtract_line_sweep_over(previous, 0.0, x, 0.0, r)
        stock.subtract_circle_sweep(x, rho, rho, r)
        previous = x
    assert previous is not None
    # Stations 7.8e-3 .. 2.4e-4 mm before V (piece length 1/16).
    parameters = [1 - 2.0**-j for j in range(3, 9)]
    rows = _stock_2.engagement_line_prefixes(stock, previous, 0.0, 10.0, 0.0, parameters, r, CAP.chord_ratio, 0.0)
    assert all(row[2] for row in rows), [math.degrees(row[1]) for row in rows]
    assert max(row[1] for row in rows) > math.pi / 2, "the reading passes the pi - alpha floor"


def test_pair_exhaustion_keeps_its_own_error() -> None:
    """A pair over the cap at depth 0 is not a connector refusal."""
    far = _on(SQUARE, 0, Fraction(2, 10), 1.0, 0), _on(SQUARE, 0, Fraction(8, 10), 1.0, 1)
    with pytest.raises(Figure5EngagementResolutionError) as caught:
        refine_figure5_engagement(far, SQUARE, TOOL, CAP, part=part_around(SQUARE, TOOL), max_depth=0)
    assert not isinstance(caught.value, ConnectorEngagementResolutionError)


def _seam_refused_path() -> tuple[Figure5CounterclockwiseCircle, ...]:
    """Two circles 0.25 apart on the octagon's bottom side: the open connector passes, the seam is the rest of the loop."""
    return _on(OCTAGON, 0, Fraction(1, 2), 1.0, 0), _on(OCTAGON, 0, Fraction(17, 32), 1.0, 1)


def test_seam_is_repaired_with_new_last_circles() -> None:
    sources = _seam_refused_path()
    unclosed = refine_figure5_engagement(sources, OCTAGON, TOOL, CAP, part=part_around(OCTAGON, TOOL))
    connectors = Figure5BoundaryConnectors.build(OCTAGON, unclosed.circles)
    assert connectors.seam is Figure5SeamClosure.CONNECTOR
    assert ConnectorEngagementAssessment.build(unclosed.circles, connectors, TOOL, CAP).refused == (len(unclosed.circles) - 1,)

    closed = refine_figure5_engagement(sources, OCTAGON, TOOL, CAP, part=part_around(OCTAGON, TOOL), close_seam=True)
    assert closed.circles[: len(unclosed.circles)] == unclosed.circles, "the seam repair only appends"
    assert closed.original_indices == unclosed.original_indices
    assert len(closed.circles) > len(unclosed.circles)
    closed_connectors = Figure5BoundaryConnectors.build(OCTAGON, closed.circles)
    certificate = ConnectorEngagementCertificate.build(closed.circles, closed_connectors, TOOL, CAP)
    assert certificate.verdicts[-1].seam
    for previous, current in zip(closed.circles, closed.circles[1:]):
        assert _pair_within(previous, current)


def test_seam_exhaustion_names_the_seam() -> None:
    with pytest.raises(ConnectorEngagementResolutionError) as caught:
        refine_figure5_engagement(_seam_refused_path(), OCTAGON, TOOL, CAP, part=part_around(OCTAGON, TOOL), max_depth=0, close_seam=True)
    assert caught.value.seam is True


def test_coincident_seam_needs_no_repair() -> None:
    first = _on(SQUARE, 0, Fraction(1, 2), 1.0, 0)
    fan = Figure5CounterclockwiseCircle(GuideRunId(1), GuideRunStationOrdinal(0), first.boundary_site, first.boundary_site, first.center, first.contact_point, first.radius)
    result = refine_figure5_engagement((first, fan), SQUARE, TOOL, CAP, part=part_around(SQUARE, TOOL), close_seam=True)
    assert result.circles == (first, fan)
    assert Figure5BoundaryConnectors.build(SQUARE, result.circles).seam is Figure5SeamClosure.COINCIDENT


def test_midpoint_across_the_seam_on_one_side_wraps_the_loop() -> None:
    """The seam's midpoint from a later to an earlier site on the same side goes once round, CCW."""
    last, first = _on(SQUARE, 0, Fraction(3, 4), 1.0, 1), _on(SQUARE, 0, Fraction(1, 4), 1.0, 0)
    halfway = midpoint(last, first, SQUARE)
    assert int(halfway.boundary_site.segment_id) == 2
    assert halfway.boundary_site.parameter == Fraction(1, 2)
