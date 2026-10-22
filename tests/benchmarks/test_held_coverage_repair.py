"""Contract tests for select, certify, repair (docs/held_coverage_preserving_selection.md).

Fixture: centre domain Omega = square (0,0)-(6,6), tool radius r = 1. A fan of
three circles at contact (3,0) (segment 0, parameter 1/2): R = 1, 2, 3 with
centres (3,1), (3,2), (3,3). One circle at contact (3,6) (segment 2, parameter
1/2): R = 1, centre (3,5). Connectors run CCW from (3,0) to (3,6) and close
the seam back to (3,0). All values are dyadic, so binary64-exact.

The connector band covers {dist(x, boundary of Omega) <= 1}; the fan covers
disk((3,3), 4). Without R = 2, (4.5, 3) and (1.5, 3) are uncovered: 1.5 from
(3,3), 2.5 from (3,1) and (3,5), 1.5 from the boundary. R = 3 is the only
member of run 1; the others are run 0. The fan's engagements are all 180
degrees (each circle surrounds its predecessor), so selection drops R = 2
only under the 180-degree cap.

Every certificate here is the real ``CoverageCertificate.build``, except the
one unit test of the greedy set cover, which is named as such.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import replace
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks import held_figure5_engagement_refinement as engagement_refinement
from benchmarks.held_boundary_order_spacing import InvalidBoundarySpacingInputError
from benchmarks.held_boundary_order_spacing import select_boundary_ordered_sources
from benchmarks.held_coverage_repair import certify_boundary_ordered_selection
from benchmarks.held_coverage_repair import path_sweeps
from benchmarks.held_coverage_repair import refine_boundary_ordered_path_certified
from benchmarks.held_coverage_repair import repair_coverage
from benchmarks.held_coverage_repair_cover import candidates_for
from benchmarks.held_coverage_repair_cover import contains_exactly
from benchmarks.held_coverage_repair_cover import outward_boxes
from benchmarks.held_coverage_repair_cover import witness_bracket
from benchmarks.held_coverage_repair_errors import CoverageNotRepairableBySelectionError
from benchmarks.held_coverage_repair_errors import CoverageRepairRefusedError
from benchmarks.held_coverage_repair_errors import InvalidCoverageRepairInputError
from benchmarks.held_coverage_repair_errors import RepairOrderViolationError
from benchmarks.held_coverage_repair_errors import UncoveredPocketWithoutWitnessesError
from benchmarks.held_coverage_repair_invariant_errors import CertificateRunCountError
from benchmarks.held_coverage_repair_invariant_errors import CertifiedLedgerMismatchError
from benchmarks.held_coverage_repair_invariant_errors import CertifiedSelectionRangeError
from benchmarks.held_coverage_repair_invariant_errors import CoverageRepairError
from benchmarks.held_coverage_repair_invariant_errors import CoverageRepairInvariantError
from benchmarks.held_coverage_repair_invariant_errors import DuplicateInsertionError
from benchmarks.held_coverage_repair_invariant_errors import EmissionIndexCountError
from benchmarks.held_coverage_repair_invariant_errors import EmissionIndexOrderError
from benchmarks.held_coverage_repair_invariant_errors import EmittedCountError
from benchmarks.held_coverage_repair_invariant_errors import InvalidRepairTimingError
from benchmarks.held_coverage_repair_invariant_errors import LedgerRoundError
from benchmarks.held_coverage_repair_invariant_errors import NegativeRefinementDepthError
from benchmarks.held_coverage_repair_invariant_errors import NegativeSourceIndexError
from benchmarks.held_coverage_repair_invariant_errors import NonPositiveCircleBudgetError
from benchmarks.held_coverage_repair_invariant_errors import NonPositiveRepairRoundError
from benchmarks.held_coverage_repair_invariant_errors import SpliceSlotError
from benchmarks.held_coverage_repair_invariant_errors import UncertifiedPathError
from benchmarks.held_coverage_repair_invariant_errors import UnorderedSelectionError
from benchmarks.held_coverage_repair_invariant_errors import UntracedRefinementError
from benchmarks.held_coverage_repair_splice import ContinuedSplice
from benchmarks.held_coverage_repair_splice import continue_with_early_stop
from benchmarks.held_coverage_repair_splice import cyclically_between
from benchmarks.held_coverage_repair_splice import emitted_form
from benchmarks.held_coverage_repair_splice import splice_source
from benchmarks.held_coverage_repair_types import BoundaryOrderedEmission
from benchmarks.held_coverage_repair_types import CertifiedBoundaryOrderedPath
from benchmarks.held_coverage_repair_types import CoverageOracle
from benchmarks.held_coverage_repair_types import CoverageRepairLedger
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_coverage_repair_types import RepairInsertion
from benchmarks.held_coverage_repair_types import RepairRefusal
from benchmarks.held_coverage_repair_types import Splice
from benchmarks.held_connector_engagement import ConnectorEngagementCertificate
from benchmarks.held_connector_engagement import UncertifiedConnectorEngagementError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementCosineRangeError
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementFragmentationError
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementRefinement
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementResolutionError
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementTraceError
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_figure5_engagement_trace import EngagementTrace
from benchmarks.held_figure5_engagement_refinement import circle_on_side
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import HypothesisFigure5Path
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import ExactEngagementCosineRangeError
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import PredecessorEngagementCosine
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

SQUARE = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (6, 0), (6, 6), (0, 6)))
TOOL = ToolRadius.build(1)
FAN_CAP = EngagementCap.build(math.pi)  # the fan's engagements are exactly 180 degrees
SETTINGS = EngagementRefinementSettings.build(TOOL, FAN_CAP)
OPEN_SETTINGS = EngagementRefinementSettings.build(TOOL, FAN_CAP, close_seam=False)
GAP_POINTS = ((Fraction(9, 2), Fraction(3)), (Fraction(3, 2), Fraction(3)))
Point = tuple[Fraction, Fraction]
# Connector-repair circles the refinement inserts on (3, 0) -> (3, 6) under the
# 180-degree cap (tier-1 bound O(d) <= pi iff d^2 <= R^2 - r^2), derived:
# leaving R = 3 about (3, 3), the corner (6, 0) has d^2 = 18 > 15 and no prior
# circle certifies it, so the midpoint (6, 3) (arclength 6 of 12, rho 2) goes
# in; that connector still crosses (6, 0), and its half is exactly the vertex,
# so the third point (5, 0) (rho corner-limited to 1) goes in and certifies
# both halves; leaving (6, 3) the corner (6, 6) has d^2 = 13 > 8 and C2 gives
# 18 > 15, so the third point (6, 5) (rho 1) goes in.
SELECTED_PATH_REPAIR_CONTACTS = ((5.0, 0.0), (6.0, 3.0), (6.0, 5.0))
# The closed path's seam (3, 6) -> (0, 6) -> (0, 0) -> (3, 0) is refined too (settings.close_seam):
# leaving R = 1 about (3, 5) the corner (0, 6) has d^2 = 10 > 3, so midpoint circles are appended
# as new last circles: (1, 6), (0, 3), and the corner-limited (0, 1), each R = 1, coordinates as the
# binary64 midpoint construction emits them (build/rs/probe_fan.py). The full post-gate certificate
# on the result is asserted in test_repair_flow_closes_and_gates_the_seam_connector.
SEAM_REPAIR_CONTACTS = ((0.9999999999999998, 6.0), (0.0, 3.0), (0.0, 0.9999999999999998))


def _circle(
    ordinal: int, segment: int, parameter: Fraction, contact: tuple[float, float], centre: tuple[float, float], radius: float, run: int = 0
) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(segment), parameter, 4)
    return Figure5CounterclockwiseCircle(
        GuideRunId(run), GuideRunStationOrdinal(ordinal), site, site, Point2[WorldXY].build(*centre), Point2[WorldXY].build(*contact), GuideRadius.build(radius)
    )


HALF = Fraction(1, 2)
C0 = _circle(0, 0, HALF, (3, 0), (3, 1), 1)
C1 = _circle(1, 0, HALF, (3, 0), (3, 2), 2)
C2 = _circle(2, 0, HALF, (3, 0), (3, 3), 3, run=1)
C3 = _circle(3, 2, HALF, (3, 6), (3, 5), 1)
FIXTURE = (C0, C1, C2, C3)
RUNS = ((GuideRunId(0),), (GuideRunId(0),), (GuideRunId(1),), (GuideRunId(0),))


def _squared_distance_to_segment(p: Point, a: Point, b: Point) -> Fraction:
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = min(max(((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy), Fraction(0)), Fraction(1))
    ex, ey = p[0] - a[0] - t * dx, p[1] - a[1] - t * dy
    return ex * ex + ey * ey


def _square_edges() -> list[tuple[Point, Point]]:
    corners = [(Fraction(float(p.x)), Fraction(float(p.y))) for p in SQUARE]
    return [(corners[k], corners[(k + 1) % 4]) for k in range(4)]


def _in_target(p: Point) -> bool:
    inside = Fraction(0) <= p[0] <= Fraction(6) and Fraction(0) <= p[1] <= Fraction(6)
    return inside or min(_squared_distance_to_segment(p, a, b) for a, b in _square_edges()) <= 1


def _in_boundary_band(p: Point) -> bool:
    """The full connector chain of the closed fixture sweeps every point within r of the square's boundary."""
    return min(_squared_distance_to_segment(p, a, b) for a, b in _square_edges()) <= 1


def _sources(circles: Sequence[Figure5CounterclockwiseCircle]) -> tuple[Figure5CounterclockwiseCircle, ...]:
    return tuple(emitted_form(circle, SQUARE) for circle in circles)


def test_boundary_selection_drops_the_closing_fan_circle_and_repair_reinserts_it() -> None:
    sources = _sources(FIXTURE)
    assert exact_predecessor_engagement(_paper_candidate(sources[0]), _paper_candidate(sources[2]), TOOL).within(FAN_CAP)
    assert select_boundary_ordered_sources(sources, RUNS, TOOL, FAN_CAP) == (0, 2, 3), "the defect: selection drops the covering R = 2"
    certified = certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, SETTINGS, CoverageCertificate.build)
    assert certified.certificate.exposed_edge_count == 0
    assert certified.ledger.inserted_sources == (1,)
    assert certified.emission.selected == (0, 1, 2, 3), "R = 2 goes back into order slot 1"
    assert certified.ledger.rounds == 1 and certified.ledger.certificate_runs == 2
    assert certified.ledger.emitted_after == certified.ledger.emitted_before + 1


def test_repair_witness_lies_in_the_reinserted_circle_exactly() -> None:
    certified = certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, SETTINGS, CoverageCertificate.build)
    (insertion,) = certified.ledger.insertions
    point = insertion.witness.uncovered_point
    assert isinstance(point[0], Fraction) and isinstance(point[1], Fraction)
    assert _in_target(point), "the witness lies in T"
    assert not _in_boundary_band(point), "the witness lies outside every connector capsule"
    kept = _sources((C0, C2, C3))
    assert not any(contains_exactly(emitted_form(circle, SQUARE), TOOL, point) for circle in kept), "the witness lies in T \\ U(0, 2, 3)"
    # Independent Fraction re-check against R = 2 about (3, 2): 1 <= |w - c| <= 3.
    distance2 = (point[0] - 3) ** 2 + (point[1] - 2) ** 2
    assert 1 <= distance2 <= 9, f"{insertion.witness} is not in the re-inserted annulus"


def test_repair_raises_when_the_circle_budget_refuses_the_only_closing_circle() -> None:
    """The refusal path: the budget admits the selected path's 9 circles and no more.

    The selected path (0, 2, 3) emits its 3 sources plus 3 connector-repair
    circles (``SELECTED_PATH_REPAIR_CONTACTS``) and 3 seam-repair circles
    (``SEAM_REPAIR_CONTACTS``), so the budget 9 is exactly its emission and
    re-inserting R = 2 needs a tenth.

    The design's cap variant (cap 1 degree, depth 0) is incoherent on this
    fixture: every fan engagement is 180 degrees, so a cap that refuses R = 2
    also refuses the selected path (0, 2, 3) itself. The circle budget is the
    other ``Figure5EngagementResolutionError`` the repair records as a refusal.
    """
    budget = EngagementRefinementSettings.build(TOOL, FAN_CAP, max_circles=3 + len(SELECTED_PATH_REPAIR_CONTACTS) + len(SEAM_REPAIR_CONTACTS))
    with pytest.raises(CoverageRepairRefusedError) as caught:
        certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, budget, CoverageCertificate.build)
    assert [refusal.candidate for refusal in caught.value.refusals] == [1]
    assert caught.value.refusals[0].interval == (0, 2)
    assert isinstance(caught.value.__cause__, Figure5EngagementResolutionError)
    assert contains_exactly(emitted_form(C1, SQUARE), TOOL, caught.value.witness.uncovered_point)


def test_repair_raises_when_no_dropped_circle_covers_the_witness() -> None:
    with pytest.raises(CoverageNotRepairableBySelectionError) as caught:
        certify_boundary_ordered_selection((C0, C2, C3), (RUNS[0], RUNS[2], RUNS[3]), SQUARE, SETTINGS, CoverageCertificate.build)
    points = [witness.uncovered_point for witness in caught.value.witnesses]
    assert all(_in_target(p) and not _in_boundary_band(p) for p in points)
    assert any((p[0] - g[0]) ** 2 + (p[1] - g[1]) ** 2 < 1 for p in points for g in GAP_POINTS), f"no witness near (4.5, 3) or (1.5, 3): {points}"


def test_repair_rejects_a_candidate_outside_its_boundary_slot() -> None:
    moved = _circle(1, 0, Fraction(1, 4), (3, 0), (3, 2), 2)  # site (0, 1/4) precedes its kept predecessor (0, 1/2)
    with pytest.raises(RepairOrderViolationError) as caught:
        certify_boundary_ordered_selection((C0, moved, C2, C3), RUNS, SQUARE, SETTINGS, CoverageCertificate.build)
    assert (caught.value.candidate, caught.value.predecessor, caught.value.successor) == (1, 0, 2)


def test_open_seam_is_not_repairable_by_selection() -> None:
    with pytest.raises(CoverageNotRepairableBySelectionError) as caught:
        certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, OPEN_SETTINGS, CoverageCertificate.build)
    kinds = {witness.piece_kind for witness in caught.value.witnesses}
    assert ExposedPieceKind.TARGET_BOUNDARY in kinds, f"no TARGET_BOUNDARY witness: {caught.value}"
    assert all(_in_target(w.uncovered_point) for w in caught.value.witnesses)


def test_repair_flow_closes_and_gates_the_seam_connector() -> None:
    """The seam connector (3, 6) -> (3, 0) round (0, 6) and (0, 0) is cut, so it is refined and post-gated too."""
    certified = certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, SETTINGS, CoverageCertificate.build)
    circles = certified.emission.refinement.circles
    assert len(circles) - 1 > certified.emission.refinement.original_indices[-1], "seam repair appends circles after the last source"
    certificate = ConnectorEngagementCertificate.build(circles, Figure5BoundaryConnectors.build(SQUARE, circles), TOOL, FAN_CAP)
    assert certificate.verdicts[-1].seam and certificate.verdicts[-1].certified


def test_post_gate_refuses_an_emission_whose_seam_was_never_refined() -> None:
    """The coverage verdict passes, but the seam of an emission refined open is not certified: the flow raises by name."""
    sources = _sources(FIXTURE)
    emission = BoundaryOrderedEmission.build((0, 1, 2, 3), refine_figure5_engagement(sources, SQUARE, TOOL, FAN_CAP))
    seam = len(emission.refinement.circles) - 1
    with pytest.raises(UncertifiedConnectorEngagementError) as caught:
        repair_coverage(sources, SQUARE, SETTINGS, emission, _scripted_oracle([_passing()], []))
    assert caught.value.connectors == (seam,)


def test_post_repair_gate_is_a_fresh_full_certificate() -> None:
    """User ruling: the post-repair verdict is a fresh full certificate over the whole repaired path.

    Replaces the design's windowed re-certification test: there is no local
    re-certification. The recorder delegates to the real certificate.
    """
    calls: list[tuple[int, int]] = []

    def recorder(domain: CentreDomain, tool_radius: ToolRadius, cuts: Sequence[FullCircleSweep], connectors: Sequence[BoundaryConnectorSweep], /) -> CoverageCertificate:
        calls.append((len(cuts), len(connectors)))
        return CoverageCertificate.build(domain, tool_radius, cuts, connectors)

    certified = certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, SETTINGS, recorder)
    circles = certified.emission.refinement.circles
    repairs = tuple((float(c.contact_point.x), float(c.contact_point.y)) for k, c in enumerate(circles) if k not in certified.emission.refinement.original_indices)
    assert repairs == SELECTED_PATH_REPAIR_CONTACTS + SEAM_REPAIR_CONTACTS
    domain = CentreDomain.build(list(SQUARE))
    cuts, connectors = path_sweeps(domain, SQUARE, circles)
    # 4 sources + 3 connector-repair + 3 seam-repair circles; the 2 moving connectors (fan -> (3, 6),
    # seam) are each split once per repair circle on them: 2 + 6 = 8.
    assert calls[-1] == (len(cuts), len(connectors)) == (len(FIXTURE) + len(repairs), 2 + len(repairs)) == (10, 8)
    assert len(calls) == certified.ledger.certificate_runs == certified.ledger.rounds + 1
    assert CoverageCertificate.build(domain, TOOL, cuts, connectors).exposed_edge_count == 0


def test_unit_fake_oracle_set_cover_prefers_most_witnesses_over_lowest_index() -> None:
    """Unit test of the greedy set cover with a FAKE ORACLE (hand-made witnesses, then a hand-made pass).

    Candidate #1 (R = 3/2 about (3, 3/2)) contains only the first witness;
    candidate #2 (R = 2 about (3, 2)) contains both, so it wins despite its
    higher index, and one insertion closes the round.
    """
    wide = _circle(1, 0, HALF, (3, 0), (3, 1.5), 1.5)
    circles = (C0, wide, C1, C2, C3)
    runs = (RUNS[0], RUNS[0], RUNS[1], RUNS[2], RUNS[3])
    witnesses = (
        UncoveredWitness(ExposedPieceKind.ANNULUS_OUTER, 0, (ExactMillimetre(Fraction(9, 2)), ExactMillimetre(Fraction(3)))),
        UncoveredWitness(ExposedPieceKind.ANNULUS_OUTER, 0, (ExactMillimetre(Fraction(3)), ExactMillimetre(Fraction(22, 5)))),
    )
    script: list[UncoveredPocketError | CoverageCertificate] = [UncoveredPocketError(witnesses, 0), CoverageCertificate(0, 0, 0, 0, 0, 0.0)]

    def fake_oracle(domain: CentreDomain, tool_radius: ToolRadius, cuts: Sequence[FullCircleSweep], connectors: Sequence[BoundaryConnectorSweep], /) -> CoverageCertificate:
        step = script.pop(0)
        if isinstance(step, UncoveredPocketError):
            raise step
        return step

    sources = _sources(circles)
    assert contains_exactly(emitted_form(sources[1], SQUARE), TOOL, witnesses[0].uncovered_point)
    assert not contains_exactly(emitted_form(sources[1], SQUARE), TOOL, witnesses[1].uncovered_point)
    certified = certify_boundary_ordered_selection(circles, runs, SQUARE, SETTINGS, fake_oracle)
    assert certified.ledger.inserted_sources == (2,)
    assert certified.ledger.certificate_runs == 2 and not script


# Splice == full re-run on the 10 x 10 square of test_held_boundary_order_spacing.py: contacts on side 0
# at x = 2 + index / 100. Contacts stay >= 2 mm from the corners, so every radius below fits; refinement
# across a convex corner fails on its own (a midpoint at the corner gets radius 0), which is not the splice.
SPLICE_BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
SPLICE_CAP = EngagementCap.build(math.radians(80))


def _side_circle(ordinal: int, index: int, radius: float) -> Figure5CounterclockwiseCircle:
    x = 2 + index / 100
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(200 + index, 1000), 4)
    return Figure5CounterclockwiseCircle(
        GuideRunId(0), GuideRunStationOrdinal(ordinal), site, site, Point2[WorldXY].build(x, radius), Point2[WorldXY].build(x, 0), GuideRadius.build(radius)
    )


@st.composite
def _splice_cases(draw: st.DrawFn) -> tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[int, ...], int, bool, int]:
    indices = draw(st.lists(st.integers(0, 600), min_size=3, max_size=14, unique=True))
    radii = draw(st.lists(st.sampled_from((1.0, 1.5, 2.0)), min_size=len(indices), max_size=len(indices)))
    circles = tuple(_side_circle(k, index, radius) for k, (index, radius) in enumerate(zip(sorted(indices), radii)))
    candidate = draw(st.integers(1, len(circles) - 2))
    keep = draw(st.lists(st.booleans(), min_size=len(circles), max_size=len(circles)))
    selected = (0, *(k for k in range(1, len(circles) - 1) if k != candidate and keep[k]), len(circles) - 1)
    return circles, selected, candidate, draw(st.booleans()), draw(st.integers(0, 4))


@settings(max_examples=150, deadline=None)
@given(_splice_cases())
def test_splice_equals_full_rerefinement(case: tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[int, ...], int, bool, int]) -> None:
    circles, selected, candidate, corners, slack = case
    sources = tuple(emitted_form(circle, SPLICE_BOUNDARY) for circle in circles)
    base = refine_figure5_engagement(tuple(sources[k] for k in selected), SPLICE_BOUNDARY, TOOL, SPLICE_CAP, corner_approaches=corners)
    budget = EngagementRefinementSettings.build(TOOL, SPLICE_CAP, max_circles=len(base.circles) + slack, corner_approaches=corners, close_seam=False)
    augmented = tuple(sorted((*selected, candidate)))
    try:
        full = refine_figure5_engagement(tuple(sources[k] for k in augmented), SPLICE_BOUNDARY, TOOL, SPLICE_CAP, max_circles=budget.max_circles, corner_approaches=corners)
    except Figure5EngagementResolutionError:
        with pytest.raises(Figure5EngagementResolutionError):
            splice_source(sources, SPLICE_BOUNDARY, budget, BoundaryOrderedEmission.build(selected, base), candidate)
        return
    spliced = splice_source(sources, SPLICE_BOUNDARY, budget, BoundaryOrderedEmission.build(selected, base), candidate)
    assert spliced.emission.selected == augmented
    assert spliced.emission.refinement == full
    assert spliced.added == len(full.circles) - len(base.circles)


def test_splice_sees_the_circles_emitted_before_its_interval() -> None:
    """Regression (Hypothesis find, 2026-10-22): the connector leaving x = 2.02 (R 1.5) for x = 2.71 (R 1)
    is certified by tier 1b from the earlier circles (x = 2.01, R 2), which a three-source local
    re-refinement never saw, so it inserted a midpoint the full re-refinement does not."""
    circles = tuple(_side_circle(k, index, radius) for k, (index, radius) in enumerate(((0, 1.0), (1, 2.0), (2, 1.5), (71, 1.0), (84, 1.0))))
    sources = tuple(emitted_form(circle, SPLICE_BOUNDARY) for circle in circles)
    base = refine_figure5_engagement(tuple(sources[k] for k in (0, 1, 2, 4)), SPLICE_BOUNDARY, TOOL, SPLICE_CAP)
    full = refine_figure5_engagement(sources, SPLICE_BOUNDARY, TOOL, SPLICE_CAP)
    spliced = splice_source(sources, SPLICE_BOUNDARY, SIDE_SETTINGS, BoundaryOrderedEmission.build((0, 1, 2, 4), base), 3)
    assert full.original_indices == (0, 4, 5, 6, 7)
    assert spliced.emission.refinement == full


# Exact early stop of the splice continuation (docs/held_coverage_preserving_selection.md, step 3).
# Fixtures on the 10 x 10 square, side 0, cap 80 degrees, open path; found by search
# (build/es/explore3.py, explore4.py) and asserted as preconditions here.
# EARLY_STOP: inserting x = 2.05 between 2.00 and 2.10 re-emits the next interval identically; the
# only later decision that reads other circles (the tier-1b connector 7.11 -> 7.85, certified by
# the R = 2 circle at 7.03) is more than rho + 2 r + H from every changed circle.
EARLY_STOP_SPEC = ((0, 1.0), (5, 1.0), (10, 1.0), (501, 1.5), (503, 2.0), (506, 1.5), (508, 1.5), (510, 1.5), (511, 1.5), (585, 1.0))
# NO_EARLY_STOP: inserting the R = 2 circle at x = 2.03 re-emits the next interval (-> 2.06)
# identically, but that circle certifies the later connector 2.11 -> 2.85 by tier 1b, which the
# reference refused (a midpoint), so adopting the reference's tail there would be wrong.
NO_EARLY_STOP_SPEC = ((1, 1.5), (3, 2.0), (6, 1.5), (8, 1.5), (10, 1.5), (11, 1.5), (85, 1.0))


def _early_stop_case(
    spec: Sequence[tuple[int, float]], candidate: int
) -> tuple[Figure5EngagementRefinement, Figure5EngagementRefinement, ContinuedSplice, tuple[Figure5CounterclockwiseCircle, ...]]:
    sources = tuple(emitted_form(_side_circle(k, index, radius), SPLICE_BOUNDARY) for k, (index, radius) in enumerate(spec))
    selected = tuple(k for k in range(len(spec)) if k != candidate)
    reference = refine_figure5_engagement(tuple(sources[k] for k in selected), SPLICE_BOUNDARY, TOOL, SPLICE_CAP)
    full = refine_figure5_engagement(sources, SPLICE_BOUNDARY, TOOL, SPLICE_CAP)
    later = tuple(sources[k] for k in selected[candidate:])
    return reference, full, continue_with_early_stop(reference, candidate, sources[candidate], later, SPLICE_BOUNDARY, SIDE_SETTINGS), later


def _footprints(refinement: Figure5EngagementRefinement) -> list[tuple[int, int | None, tuple[float, ...] | None]]:
    assert refinement.trace is not None
    return [
        (f.position, f.vertex, None if f.pieces is None else tuple(float(v) for v in (*f.pieces.mid_lo.ravel(), *f.pieces.mid_hi.ravel(), *f.pieces.half, *f.pieces.extent)))
        for f in refinement.trace.footprints
    ]


def test_splice_stops_early_when_no_changed_circle_reaches_a_later_decision() -> None:
    reference, full, spliced, later = _early_stop_case(EARLY_STOP_SPEC, 1)
    assert spliced.refinement == full
    assert spliced.rerefined == 1 < len(later), "only the interval after the inserted source is re-refined"
    assert reference.trace is not None
    stop = reference.original_indices[1] + 1
    assert any(f.position >= stop and f.pieces is not None for f in reference.trace.footprints), "a later tier-1b decision exists, so the stop was decided, not vacuous"
    assert _footprints(spliced.refinement) == _footprints(full), "the stitched trace is the full refinement's trace"


def test_splice_does_not_stop_early_when_a_changed_circle_certifies_a_later_connector() -> None:
    reference, full, spliced, later = _early_stop_case(NO_EARLY_STOP_SPEC, 1)
    first_interval = full.circles[full.original_indices[1] + 1 : full.original_indices[2] + 1]
    assert first_interval == reference.circles[reference.original_indices[0] + 1 : reference.original_indices[1] + 1], "precondition: re-emitted identically"
    assert reference.circles[reference.original_indices[1] + 1 :] != full.circles[full.original_indices[2] + 1 :], "precondition: the reference's tail is wrong here"
    assert reference.trace is not None
    assert reference.trace.last_influenced((full.circles[full.original_indices[1]],), reference.original_indices[1] + 1) > reference.original_indices[1]
    assert spliced.refinement == full
    assert spliced.rerefined == len(later)


def test_a_mismatch_resets_the_window_and_the_stop_reads_both_sides(monkeypatch: pytest.MonkeyPatch) -> None:
    """Match, mismatch, match: the influence test is re-run after the mismatch, on X and Y alike.

    NO_EARLY_STOP_SPEC plus x = 2.86 and 2.95 (R 1). The interval after the slot re-emits identically
    (first call, on X = {2.02, 2.03}); 2.11 -> 2.85 does not (the inserted circle certifies it, so the
    reference's midpoint 2.48 R 1.25 is gone); 2.85 -> 2.86 re-emits again. The second call must see the
    replaced reference circles (Y, the midpoint among them) besides X. No fixture among about 110k
    searched (build/es/search_inj2.py) makes a replaced circle flip a later decision that no new circle
    reaches, so this pins the window's state and both sides directly; the proof is in the docstring.
    """
    calls: list[tuple[set[tuple[float, float]], int]] = []
    original = EngagementTrace.last_influenced

    def spy(self: EngagementTrace, circles: Sequence[Figure5CounterclockwiseCircle], since: int) -> int:
        calls.append(({(float(c.contact_point.x), float(c.radius.value)) for c in circles}, since))
        return original(self, circles, since)

    monkeypatch.setattr(EngagementTrace, "last_influenced", spy)
    reference, full, spliced, later = _early_stop_case((*NO_EARLY_STOP_SPEC, (86, 1.0), (95, 1.0)), 1)
    assert spliced.refinement == full
    assert spliced.rerefined == 6 < len(later)
    assert [since for _, since in calls] == [reference.original_indices[1] + 1, reference.original_indices[6] + 1], "re-run after the mismatch"
    assert {(2.02, 1.75), (2.03, 2.0)} <= calls[0][0]
    assert {(2.02, 1.75), (2.03, 2.0), (2.48, 1.25)} <= calls[1][0], "X and the replaced midpoint (Y)"


def test_splice_without_a_reference_trace_is_a_contract_breach() -> None:
    reference, _, _, later = _early_stop_case(EARLY_STOP_SPEC, 1)
    sources = tuple(emitted_form(_side_circle(k, index, radius), SPLICE_BOUNDARY) for k, (index, radius) in enumerate(EARLY_STOP_SPEC))
    bare = Figure5EngagementRefinement(reference.circles, reference.original_indices)
    with pytest.raises(Figure5EngagementTraceError):
        continue_with_early_stop(bare, 1, sources[1], later, SPLICE_BOUNDARY, SIDE_SETTINGS)
    assert not issubclass(Figure5EngagementTraceError, Figure5EngagementResolutionError), "never recorded as a refusal"
    with pytest.raises(InvalidCoverageRepairInputError):
        continue_with_early_stop(reference, 1, sources[1], later[1:], SPLICE_BOUNDARY, SIDE_SETTINGS)
    with pytest.raises(InvalidCoverageRepairInputError):
        continue_with_early_stop(reference, 0, sources[1], (sources[0], *later), SPLICE_BOUNDARY, SIDE_SETTINGS)
    with pytest.raises(UntracedRefinementError):
        BoundaryOrderedEmission.build(tuple(k for k in range(len(EARLY_STOP_SPEC)) if k != 1), bare)


# Shipped configuration (close_seam=True, corner_approaches=True) on the octagon of
# test_held_connector_repair.py (obtuse convex corners, so the closed seam certifies): the
# reference appends 28 seam circles after its last source, and the splice of source 4 stops after
# 2 of 4 later sources, adopting them (found by build/es/seam_probe.py, asserted below).
OCTAGON = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (8, 0), (12, 4), (12, 12), (8, 16), (0, 16), (-4, 12), (-4, 4)))
SHIPPED_SPEC = (
    (0, Fraction(77, 100), 1.5),
    (2, Fraction(9, 20), 1.5),
    (2, Fraction(63, 100), 2.0),
    (3, Fraction(1, 4), 1.5),
    (3, Fraction(1, 2), 2.0),
    (3, Fraction(41, 50), 1.5),
    (5, Fraction(1, 5), 2.0),
    (5, Fraction(17, 20), 2.0),
    (7, Fraction(11, 25), 0.5),
)
# The tooth of test_held_connector_repair.py: the R = 1.8 circle tangent at (0, 0) leaves the wedge of
# the sharp vertex 2, so the arrival side 1 -> side 2 records floor footprints after the slot.
TOOTH = tuple(Point2[WorldXY].build(x, y) for x, y in ((-3.0, 0.0), (0.25, 0.0), (0.5, 0.0625), (-0.5, 3.375), (-3.0, 3.375)))
FLOOR_SPEC = ((0, Fraction(3, 10), 1.0), (0, Fraction(3, 5), 1.0), (0, Fraction(12, 13), 1.8), (1, Fraction(1, 2), 0.25), (2, Fraction(7, 8), 0.5))


def _polygon_sources(boundary: tuple[Point2[WorldXY], ...], spec: Sequence[tuple[int, Fraction, float]]) -> tuple[Figure5CounterclockwiseCircle, ...]:
    circles = []
    for ordinal, (side, parameter, radius) in enumerate(spec):
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, len(boundary))
        contact = _canonical_offset_point(boundary, site)
        raw = Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(ordinal), site, site, contact, contact, GuideRadius.build(radius))
        circles.append(circle_on_side(raw, site, contact, GuideRadius.build(radius), boundary))
    return tuple(circles)


def _polygon_splice(
    boundary: tuple[Point2[WorldXY], ...], spec: Sequence[tuple[int, Fraction, float]], candidate: int, settings_: EngagementRefinementSettings
) -> tuple[TracedEngagementRefinement, TracedEngagementRefinement, ContinuedSplice, int]:
    sources = _polygon_sources(boundary, spec)
    selected = tuple(k for k in range(len(sources)) if k != candidate)
    kw = {"corner_approaches": settings_.corner_approaches, "close_seam": settings_.close_seam}
    reference = refine_figure5_engagement(tuple(sources[k] for k in selected), boundary, TOOL, SPLICE_CAP, **kw)
    full = refine_figure5_engagement(sources, boundary, TOOL, SPLICE_CAP, **kw)
    later = tuple(sources[k] for k in selected[candidate:])
    return reference, full, continue_with_early_stop(reference, candidate, sources[candidate], later, boundary, settings_), len(later)


def test_splice_in_the_shipped_configuration_adopts_the_seam_circles() -> None:
    shipped = EngagementRefinementSettings.build(TOOL, SPLICE_CAP, corner_approaches=True, close_seam=True)
    reference, full, spliced, later = _polygon_splice(OCTAGON, SHIPPED_SPEC, 4, shipped)
    assert len(reference.circles) - 1 - reference.original_indices[-1] == 28, "precondition: the reference ends in seam circles"
    assert spliced.rerefined == 2 < later, "the seam circles are adopted, not re-refined"
    assert spliced.refinement == full
    assert spliced.refinement.trace == full.trace


def test_splice_with_floor_footprints_after_the_slot_equals_the_full_refinement() -> None:
    reference, full, spliced, later = _polygon_splice(TOOTH, FLOOR_SPEC, 1, SIDE_SETTINGS)
    start = reference.original_indices[0]
    assert any(f.vertex == 2 and f.position > start for f in reference.trace.footprints), "precondition: a floor decision after the slot"
    assert spliced.refinement == full
    assert spliced.refinement.trace == full.trace


@settings(max_examples=100, deadline=None)
@given(_splice_cases(), st.integers(0, 1000))
def test_chained_splices_equal_the_full_refinement_and_its_trace(case: tuple[tuple[Figure5CounterclockwiseCircle, ...], tuple[int, ...], int, bool, int], pick: int) -> None:
    """Two greedy steps: the second splice adopts from the first splice's stitched trace."""
    circles, selected, candidate, corners, _ = case
    sources = tuple(emitted_form(circle, SPLICE_BOUNDARY) for circle in circles)
    settings_ = EngagementRefinementSettings.build(TOOL, SPLICE_CAP, corner_approaches=corners, close_seam=False)
    try:
        base = refine_figure5_engagement(tuple(sources[k] for k in selected), SPLICE_BOUNDARY, TOOL, SPLICE_CAP, corner_approaches=corners)
        first = splice_source(sources, SPLICE_BOUNDARY, settings_, BoundaryOrderedEmission.build(selected, base), candidate).emission
    except Figure5EngagementResolutionError:
        return
    dropped = [k for k in range(1, len(circles) - 1) if k not in first.selected]
    if not dropped:
        return
    second = dropped[pick % len(dropped)]
    augmented = tuple(sorted((*first.selected, second)))
    try:
        full = refine_figure5_engagement(tuple(sources[k] for k in augmented), SPLICE_BOUNDARY, TOOL, SPLICE_CAP, corner_approaches=corners)
    except Figure5EngagementResolutionError:
        with pytest.raises(Figure5EngagementResolutionError):
            splice_source(sources, SPLICE_BOUNDARY, settings_, first, second)
        return
    spliced = splice_source(sources, SPLICE_BOUNDARY, settings_, first, second).emission.refinement
    assert spliced == full
    assert _footprints(spliced) == _footprints(full)


# ---------------------------------------------------------------------------
# Refusal ruling and greedy tie-break (review 2026-10-21, I1/I2, minors 1-3, 6, 7).
# Fixtures on the 10 x 10 square, side 0, cap 80 degrees: every source is a side circle
# (_side_circle), the kept selection is given by hand and repair_coverage runs directly.
# Witnesses come from a FAKE ORACLE (hand-made exact points, then a hand-made pass);
# each witness's membership is re-checked exactly in the test. The ``added`` values are
# asserted as preconditions with splice_source, so a fixture drift fails loudly.

# Open paths: the 10 x 10 square's seam turns three right-angle convex corners, which no refinement
# certifies under an 80-degree cap (docs/held_connector_engagement.md, convex-corner lemma).
SIDE_SETTINGS = EngagementRefinementSettings.build(TOOL, SPLICE_CAP, close_seam=False)


def _scripted_oracle(script: list[UncoveredPocketError | CoverageCertificate], calls: list[int]) -> CoverageOracle:
    def oracle(domain: CentreDomain, tool_radius: ToolRadius, cuts: Sequence[FullCircleSweep], connectors: Sequence[BoundaryConnectorSweep], /) -> CoverageCertificate:
        calls.append(len(cuts))
        step = script.pop(0)
        if isinstance(step, UncoveredPocketError):
            raise step
        return step

    return oracle


def _witness(x: Fraction, y: Fraction) -> UncoveredWitness:
    return UncoveredWitness(ExposedPieceKind.ANNULUS_OUTER, 0, (ExactMillimetre(x), ExactMillimetre(y)))


def _passing() -> CoverageCertificate:
    return CoverageCertificate(0, 0, 0, 0, 0, 0.0)


def _side_case(
    spec: Sequence[tuple[int, float]], selected: tuple[int, ...], settings: EngagementRefinementSettings = SIDE_SETTINGS
) -> tuple[tuple[Figure5CounterclockwiseCircle, ...], BoundaryOrderedEmission]:
    sources = tuple(emitted_form(_side_circle(k, index, radius), SPLICE_BOUNDARY) for k, (index, radius) in enumerate(spec))
    refinement = refine_figure5_engagement(tuple(sources[k] for k in selected), SPLICE_BOUNDARY, TOOL, SPLICE_CAP, max_circles=settings.max_circles)
    return sources, BoundaryOrderedEmission.build(selected, refinement)


def _added(sources: tuple[Figure5CounterclockwiseCircle, ...], emission: BoundaryOrderedEmission, candidate: int) -> int:
    return splice_source(sources, SPLICE_BOUNDARY, SIDE_SETTINGS, emission, candidate).added


def _members(sources: tuple[Figure5CounterclockwiseCircle, ...], point: Point) -> set[int]:
    return {k for k, source in enumerate(sources) if contains_exactly(source, TOOL, point)}


def _repair(
    sources: tuple[Figure5CounterclockwiseCircle, ...],
    emission: BoundaryOrderedEmission,
    witnesses: Sequence[UncoveredWitness],
    settings: EngagementRefinementSettings = SIDE_SETTINGS,
) -> CertifiedBoundaryOrderedPath:
    calls: list[int] = []
    script: list[UncoveredPocketError | CoverageCertificate] = [UncoveredPocketError(witnesses, 0), _passing()]
    certified = repair_coverage(sources, SPLICE_BOUNDARY, settings, emission, _scripted_oracle(script, calls))
    assert not script and len(calls) == 2
    return certified


# Tie on count and added: candidates #1 (x = 4.9) and #2 (x = 5.1), R = 1, both splice with added 0.
TIE_SPEC = ((0, 1.0), (290, 1.0), (310, 1.0), (600, 1.0))
TIE_POINT = (Fraction(5), Fraction(1))


def _tie_case() -> tuple[tuple[Figure5CounterclockwiseCircle, ...], BoundaryOrderedEmission]:
    sources, emission = _side_case(TIE_SPEC, (0, 3))
    assert {1, 2} <= _members(sources, TIE_POINT)
    assert _added(sources, emission, 1) == _added(sources, emission, 2) == 0
    return sources, emission


def test_tie_on_count_and_added_inserts_the_lowest_source_index() -> None:
    sources, emission = _tie_case()
    certified = _repair(sources, emission, (_witness(*TIE_POINT),))
    assert certified.ledger.inserted_sources == (1,)
    assert certified.ledger.refusals == ()


def test_tie_on_count_inserts_the_candidate_with_fewest_added_circles_over_a_lower_index() -> None:
    # #1: R = 2 at x = 3 (added 1); #2: R = 1 at x = 5 (added 0). (4, 2) lies in both annuli.
    sources, emission = _side_case(((0, 1.0), (100, 2.0), (300, 1.0), (600, 1.0)), (0, 3))
    point = (Fraction(4), Fraction(2))
    assert {1, 2} <= _members(sources, point)
    assert _added(sources, emission, 1) == 1 and _added(sources, emission, 2) == 0
    assert _repair(sources, emission, (_witness(*point),)).ledger.inserted_sources == (2,)


def test_more_witnesses_beat_fewer_added_circles() -> None:
    # #1 (x = 3, added 4) contains A only; #2 (x = 4, added 8) contains A and B: one insertion of #2.
    sources, emission = _side_case(((0, 1.0), (100, 1.0), (200, 1.0), (600, 1.0)), (0, 3))
    a, b = (Fraction(7, 2), Fraction(1)), (Fraction(11, 2), Fraction(1))
    assert {1, 2} <= _members(sources, a) and 2 in _members(sources, b) and 1 not in _members(sources, b)
    assert _added(sources, emission, 1) == 4 and _added(sources, emission, 2) == 8
    assert _repair(sources, emission, (_witness(*a), _witness(*b))).ledger.inserted_sources == (2,)


def test_budget_refusal_of_the_first_candidate_inserts_the_second() -> None:
    """The ruling: a refused candidate is recorded and the next candidate is tried.

    #1 (R = 2 at x = 3) needs one more circle, which the budget (the selected
    path's own count) refuses; #2 (R = 1 at x = 5) adds none and is inserted.
    """
    spec = ((0, 1.0), (100, 2.0), (300, 1.0), (600, 1.0))
    sources, emission = _side_case(spec, (0, 3))
    budget = EngagementRefinementSettings.build(TOOL, SPLICE_CAP, max_circles=len(emission.refinement.circles), close_seam=False)
    point = (Fraction(4), Fraction(2))
    certified = _repair(sources, emission, (_witness(*point),), budget)
    assert certified.ledger.inserted_sources == (2,)
    assert [(r.candidate, r.interval, r.round) for r in certified.ledger.refusals] == [(1, (0, 3), 1)]
    assert "budget" in certified.ledger.refusals[0].reason


def _refuse_target(monkeypatch: pytest.MonkeyPatch, refused: Figure5CounterclockwiseCircle, error: Exception) -> None:
    """Inject ``error`` from the exact engagement whenever ``refused`` is the transition's target.

    No natural fixture is known: 800,000 random quarter-unit circle pairs never
    raised either error (build/coverage-repair-fixes/frag_search*.py), so the
    engagement is replaced in the refinement module only (selection is not
    touched).
    """
    real = exact_predecessor_engagement

    def engagement(predecessor: PaperCircleCandidate, candidate: PaperCircleCandidate, tool_radius: ToolRadius) -> PredecessorEngagementCosine:
        if candidate.center == refused.center and candidate.guide_radius == refused.radius:
            raise error
        return real(predecessor, candidate, tool_radius)

    monkeypatch.setattr(engagement_refinement, "exact_predecessor_engagement", engagement)


PLACEMENT_ERRORS = (
    (StandardPlacementFragmentationError("Paper overlap correction has no real swept-disk intersection."), Figure5EngagementFragmentationError),
    (ExactEngagementCosineRangeError("Exact standard-model engagement cosine is outside [-1, 1]."), Figure5EngagementCosineRangeError),
)


@pytest.mark.parametrize(("raised", "wrapped"), PLACEMENT_ERRORS)
def test_refinement_wraps_placement_errors_as_named_resolution_errors(monkeypatch: pytest.MonkeyPatch, raised: Exception, wrapped: type[Exception]) -> None:
    sources, _ = _tie_case()
    _refuse_target(monkeypatch, sources[1], raised)
    with pytest.raises(wrapped, match="source interval 0") as caught:
        refine_figure5_engagement(sources[:2], SPLICE_BOUNDARY, TOOL, SPLICE_CAP)
    assert isinstance(caught.value, Figure5EngagementResolutionError)
    assert caught.value.__cause__ is raised


@pytest.mark.parametrize(("raised", "wrapped"), PLACEMENT_ERRORS)
def test_placement_error_refuses_the_candidate_and_the_next_is_inserted(monkeypatch: pytest.MonkeyPatch, raised: Exception, wrapped: type[Exception]) -> None:
    """User ruling 2026-10-21: fragmentation / cosine-range while re-refining a candidate is a refusal of it."""
    sources, emission = _tie_case()
    _refuse_target(monkeypatch, sources[1], raised)
    certified = _repair(sources, emission, (_witness(*TIE_POINT),))
    assert certified.ledger.inserted_sources == (2,), "the lowest index is refused, so the tie goes to #2"
    (refusal,) = certified.ledger.refusals
    assert (refusal.candidate, refusal.interval, refusal.round) == (1, (0, 3), 1)
    assert refusal.reason.startswith(wrapped.__name__)


@pytest.mark.parametrize(("raised", "wrapped"), PLACEMENT_ERRORS)
def test_placement_error_on_every_candidate_raises_refused_with_the_wrapped_cause(monkeypatch: pytest.MonkeyPatch, raised: Exception, wrapped: type[Exception]) -> None:
    sources, emission = _side_case(TIE_SPEC[:2] + TIE_SPEC[3:], (0, 2))
    _refuse_target(monkeypatch, sources[1], raised)
    with pytest.raises(CoverageRepairRefusedError) as caught:
        _repair(sources, emission, (_witness(*TIE_POINT),))
    assert [r.candidate for r in caught.value.refusals] == [1]
    assert isinstance(caught.value.__cause__, wrapped)


# Sources at x = 2, 3.5, 4, 5, 6, 7, 8 (R = 1 except #1, R = 3); kept (0, 2, 4, 6).
# #1 (refused) contains A and B; #3 contains A only; #5 contains B only. Two greedy steps
# in one round; #1's interval (0, 2) is unchanged by either insertion.
SPREAD_SPEC = ((0, 1.0), (150, 3.0), (200, 1.0), (300, 1.0), (400, 1.0), (500, 1.0), (600, 1.0))
SPREAD_A = (Fraction(9, 2), Fraction(1, 2))
SPREAD_B = (Fraction(13, 2), Fraction(5, 2))


def test_a_refused_candidate_is_recorded_once_per_interval_and_round(monkeypatch: pytest.MonkeyPatch) -> None:
    sources, emission = _side_case(SPREAD_SPEC, (0, 2, 4, 6))
    assert {1, 3} <= _members(sources, SPREAD_A) and 5 not in _members(sources, SPREAD_A)
    assert {1, 5} <= _members(sources, SPREAD_B) and 3 not in _members(sources, SPREAD_B)
    _refuse_target(monkeypatch, sources[1], PLACEMENT_ERRORS[0][0])
    certified = _repair(sources, emission, (_witness(*SPREAD_A), _witness(*SPREAD_B)))
    assert sorted(certified.ledger.inserted_sources) == [3, 5]
    assert [(r.candidate, r.interval, r.round) for r in certified.ledger.refusals] == [(1, (0, 2), 1)]


def test_refused_error_chains_a_refusal_of_the_named_witness(monkeypatch: pytest.MonkeyPatch) -> None:
    """Witness A's only candidates (#1, #3) are refused; #5's later refusal must not be the cause."""
    sources, emission = _side_case(SPREAD_SPEC, (0, 2, 4, 6))
    causes = {1: RuntimeError("one"), 3: RuntimeError("three"), 5: RuntimeError("five")}
    real = exact_predecessor_engagement

    def engagement(predecessor: PaperCircleCandidate, candidate: PaperCircleCandidate, tool_radius: ToolRadius) -> PredecessorEngagementCosine:
        for index, cause in causes.items():
            if candidate.center == sources[index].center and candidate.guide_radius == sources[index].radius:
                raise StandardPlacementFragmentationError(str(cause))
        return real(predecessor, candidate, tool_radius)

    monkeypatch.setattr(engagement_refinement, "exact_predecessor_engagement", engagement)
    with pytest.raises(CoverageRepairRefusedError) as caught:
        _repair(sources, emission, (_witness(*SPREAD_A), _witness(*SPREAD_B)))
    assert caught.value.witness.uncovered_point == SPREAD_A
    cause = caught.value.__cause__
    assert isinstance(cause, Figure5EngagementFragmentationError) and isinstance(cause.__cause__, StandardPlacementFragmentationError)
    assert str(cause.__cause__) in {"one", "three"}, f"cause {cause.__cause__} is not a refusal of the named witness"


def test_an_uncovered_pocket_without_witnesses_is_a_named_error_not_a_livelock() -> None:
    sources, emission = _tie_case()
    calls: list[int] = []
    script: list[UncoveredPocketError | CoverageCertificate] = [UncoveredPocketError((), 3), _passing()]
    with pytest.raises(UncoveredPocketWithoutWitnessesError) as caught:
        repair_coverage(sources, SPLICE_BOUNDARY, SIDE_SETTINGS, emission, _scripted_oracle(script, calls))
    assert caught.value.undecided_edge_count == 3
    assert isinstance(caught.value.__cause__, UncoveredPocketError)
    assert calls == [len(emission.refinement.circles)]


def test_run_ledger_is_checked_before_certification() -> None:
    """Minor 6: as in ``refine_boundary_ordered_path``, a ledger mismatch fails before any certificate runs."""
    calls: list[int] = []
    path = Figure5BoundaryPath(FIXTURE, (), frozenset({GuideRunId(0)}))
    source = HypothesisFigure5Path(path, (), RUNS, frozenset({GuideRunId(0)}), 0, 0, 0)
    with pytest.raises(InvalidBoundarySpacingInputError):
        refine_boundary_ordered_path_certified(source, SQUARE, TOOL, FAN_CAP, oracle=_scripted_oracle([_passing()], calls))
    assert calls == []


def _linear_position(key: tuple[int, Fraction], sides: int) -> Fraction:
    return Fraction(key[0]) + key[1] if key[0] < sides else Fraction(-1)


_SITE_KEYS = st.tuples(st.integers(0, 3), st.sampled_from((Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))))


@given(_SITE_KEYS, _SITE_KEYS, _SITE_KEYS)
def test_cyclically_between_matches_a_linearised_loop(first: tuple[int, Fraction], middle: tuple[int, Fraction], last: tuple[int, Fraction]) -> None:
    """Reference: on a loop of 4 sides, middle is between iff (middle - first) mod 4 <= (last - first) mod 4."""
    f, m, la = (_linear_position(key, 4) for key in (first, middle, last))
    assert cyclically_between(first, middle, last) == ((m - f) % 4 <= (la - f) % 4)


def test_cyclically_between_wraps_through_the_seam() -> None:
    assert cyclically_between((3, Fraction(1, 2)), (0, Fraction(1, 4)), (1, Fraction(0)))
    assert cyclically_between((3, Fraction(1, 2)), (3, Fraction(3, 4)), (1, Fraction(0)))
    assert not cyclically_between((3, Fraction(1, 2)), (2, Fraction(0)), (1, Fraction(0)))


_BOX_COORDINATE = st.floats(-1e6, 1e6, allow_nan=False, allow_infinity=False)
_BOX_RADIUS = st.floats(1e-6, 1e3, allow_nan=False, allow_infinity=False)


@settings(max_examples=300, deadline=None)
@given(_BOX_COORDINATE, _BOX_COORDINATE, _BOX_RADIUS, _BOX_RADIUS)
def test_outward_box_contains_the_exact_swept_disk(x: float, y: float, radius: float, tool: float) -> None:
    """The filter box brackets [c - (R + r), c + (R + r)] exactly (Fraction), per axis."""
    circle = Figure5CounterclockwiseCircle(
        GuideRunId(0),
        GuideRunStationOrdinal(0),
        C0.source_boundary_site,
        C0.boundary_site,
        Point2[WorldXY].build(x, y + radius),
        Point2[WorldXY].build(x, y),
        GuideRadius.build(radius),
    )
    tool_radius = ToolRadius.build(tool)
    ((xlo, xhi, ylo, yhi),) = outward_boxes((circle,), tool_radius).tolist()
    reach = Fraction(float(circle.radius.value)) + Fraction(tool)
    cx, cy = Fraction(float(circle.center.x)), Fraction(float(circle.center.y))
    assert Fraction(xlo) <= cx - reach and cx + reach <= Fraction(xhi)
    assert Fraction(ylo) <= cy - reach and cy + reach <= Fraction(yhi)


@given(st.fractions(min_value=-1000, max_value=1000))
def test_witness_bracket_contains_the_exact_coordinate(value: Fraction) -> None:
    low, high = witness_bracket(value)
    assert Fraction(low) <= value <= Fraction(high)


def test_a_witness_exactly_on_the_outer_rim_is_a_candidate() -> None:
    """Centre x = 0.1 (not dyadic), R = r = 1: the witness c + (R + r) on the rim is found, not filtered out."""
    circle = Figure5CounterclockwiseCircle(
        GuideRunId(0),
        GuideRunStationOrdinal(0),
        C0.source_boundary_site,
        C0.boundary_site,
        Point2[WorldXY].build(0.1, 1.0),
        Point2[WorldXY].build(0.1, 0.0),
        GuideRadius.build(1.0),
    )
    rim = (Fraction(0.1) + 2, Fraction(1))
    assert contains_exactly(circle, TOOL, rim)
    found = candidates_for(_witness(*rim), (7,), (circle,), outward_boxes((circle,), TOOL), TOOL)
    assert found == (7,)


# ---------------------------------------------------------------------------
# Value-type invariants: each ``build`` factory raises one named error per violated invariant;
# the raw constructors are unvalidated. Baseline: the real repair of the fan fixture
# (inserted (1,), one round, two certificate runs), rebuilt field by field.


@functools.cache
def _certified() -> CertifiedBoundaryOrderedPath:
    return certify_boundary_ordered_selection(FIXTURE, RUNS, SQUARE, SETTINGS, CoverageCertificate.build)


def _rebuild_ledger(ledger: CoverageRepairLedger) -> CoverageRepairLedger:
    return CoverageRepairLedger.build(
        ledger.insertions,
        ledger.refusals,
        rounds=ledger.rounds,
        emitted_before=ledger.emitted_before,
        emitted_after=ledger.emitted_after,
        repair_s=ledger.repair_s,
        certificate_s=ledger.certificate_s,
        certificate_runs=ledger.certificate_runs,
    )


def _rebuild_path(path: CertifiedBoundaryOrderedPath) -> CertifiedBoundaryOrderedPath:
    return CertifiedBoundaryOrderedPath.build(path.sources, path.emission, path.certificate, path.ledger)


def _raises_invariant(error: type[CoverageRepairInvariantError], build: Callable[[], object]) -> None:
    with pytest.raises(error) as caught:
        build()
    assert isinstance(caught.value, CoverageRepairError)


def test_build_round_trips_every_value_type_of_a_real_repair() -> None:
    certified = _certified()
    emission, ledger = certified.emission, certified.ledger
    (insertion,) = ledger.insertions
    assert EngagementRefinementSettings.build(TOOL, FAN_CAP, max_depth=0, max_circles=1) == EngagementRefinementSettings(TOOL, FAN_CAP, 0, 1)
    assert BoundaryOrderedEmission.build(list(emission.selected), emission.refinement) == emission
    assert RepairInsertion.build(insertion.source, insertion.witness, insertion.round) == insertion
    assert _rebuild_ledger(ledger) == ledger
    assert _rebuild_path(certified) == certified
    base = BoundaryOrderedEmission.build((0, 2, 3), refine_figure5_engagement(tuple(certified.sources[k] for k in (0, 2, 3)), SQUARE, TOOL, FAN_CAP))
    splice = splice_source(certified.sources, SQUARE, SETTINGS, base, 1)
    assert Splice.build(splice.emission, splice.candidate, splice.interval, splice.added) == splice
    assert RepairRefusal.build(1, (0, 2), 1, "budget") == RepairRefusal(1, (0, 2), 1, "budget")


@pytest.mark.parametrize(
    ("error", "max_depth", "max_circles"),
    ((NegativeRefinementDepthError, -1, 1), (NonPositiveCircleBudgetError, 0, 0)),
)
def test_settings_build_rejects_an_invalid_depth_or_budget(error: type[CoverageRepairInvariantError], max_depth: int, max_circles: int) -> None:
    _raises_invariant(error, lambda: EngagementRefinementSettings.build(TOOL, FAN_CAP, max_depth=max_depth, max_circles=max_circles))


def test_emission_build_rejects_each_broken_invariant() -> None:
    refinement = _certified().emission.refinement
    circles, indices = refinement.circles, refinement.original_indices
    _raises_invariant(NegativeSourceIndexError, lambda: BoundaryOrderedEmission.build((-1, 1, 2, 3), refinement))
    _raises_invariant(UnorderedSelectionError, lambda: BoundaryOrderedEmission.build((0, 2, 1, 3), refinement))
    _raises_invariant(UnorderedSelectionError, lambda: BoundaryOrderedEmission.build((0, 1, 1, 3), refinement))
    _raises_invariant(EmissionIndexCountError, lambda: BoundaryOrderedEmission.build((0, 1, 2), refinement))
    _raises_invariant(EmissionIndexOrderError, lambda: BoundaryOrderedEmission.build((0, 1, 2, 3), Figure5EngagementRefinement(circles, indices[::-1])))
    _raises_invariant(EmissionIndexOrderError, lambda: BoundaryOrderedEmission.build((0, 1, 2, 3), Figure5EngagementRefinement(circles, (*indices[:-1], len(circles)))))
    _raises_invariant(EmissionIndexOrderError, lambda: BoundaryOrderedEmission.build((0, 1, 2, 3), Figure5EngagementRefinement(circles, (-1, *indices[1:]))))


def test_splice_build_rejects_a_candidate_outside_its_interval() -> None:
    emission = _certified().emission  # selected (0, 1, 2, 3)
    assert Splice.build(emission, 1, (0, 2), 1).candidate == 1
    _raises_invariant(SpliceSlotError, lambda: Splice.build(emission, 1, (2, 0), 1))
    _raises_invariant(SpliceSlotError, lambda: Splice.build(emission, 1, (0, 3), 1))
    _raises_invariant(SpliceSlotError, lambda: Splice.build(emission, 7, (0, 2), 1))
    _raises_invariant(SpliceSlotError, lambda: Splice.build(emission, 0, (3, 1), 1))
    _raises_invariant(SpliceSlotError, lambda: Splice.build(emission, 3, (2, 0), 1))


def test_insertion_and_refusal_build_reject_negative_indices_and_rounds() -> None:
    witness = _certified().ledger.insertions[0].witness
    _raises_invariant(NegativeSourceIndexError, lambda: RepairInsertion.build(-1, witness, 1))
    _raises_invariant(NonPositiveRepairRoundError, lambda: RepairInsertion.build(1, witness, 0))
    _raises_invariant(NegativeSourceIndexError, lambda: RepairRefusal.build(-1, (0, 2), 1, "budget"))
    _raises_invariant(NegativeSourceIndexError, lambda: RepairRefusal.build(1, (-1, 2), 1, "budget"))
    _raises_invariant(NonPositiveRepairRoundError, lambda: RepairRefusal.build(1, (0, 2), 0, "budget"))


def _ledger_variants() -> list[tuple[type[CoverageRepairInvariantError], Callable[[CoverageRepairLedger], CoverageRepairLedger]]]:
    def second(ledger: CoverageRepairLedger, round_index: int) -> RepairInsertion:
        return replace(ledger.insertions[0], source=0, round=round_index)

    return [
        (CertificateRunCountError, lambda ledger: replace(ledger, certificate_runs=1)),
        (CertificateRunCountError, lambda ledger: replace(ledger, rounds=-1, certificate_runs=0)),
        (LedgerRoundError, lambda ledger: replace(ledger, rounds=2, certificate_runs=3)),
        (LedgerRoundError, lambda ledger: replace(ledger, insertions=(second(ledger, 2), ledger.insertions[0]), rounds=2, certificate_runs=3)),
        (LedgerRoundError, lambda ledger: replace(ledger, refusals=(RepairRefusal(1, (0, 2), 2, "budget"),))),
        (LedgerRoundError, lambda ledger: replace(ledger, refusals=(RepairRefusal(1, (0, 2), 0, "budget"),))),
        (LedgerRoundError, lambda ledger: replace(ledger, refusals=(RepairRefusal(1, (0, 2), 2, "a"), RepairRefusal(1, (0, 2), 1, "b")), rounds=2, certificate_runs=3)),
        (DuplicateInsertionError, lambda ledger: replace(ledger, insertions=ledger.insertions * 2)),
        (EmittedCountError, lambda ledger: replace(ledger, emitted_before=-1)),
        (EmittedCountError, lambda ledger: replace(ledger, emitted_after=-1)),
        (EmittedCountError, lambda ledger: replace(ledger, insertions=(), rounds=0, certificate_runs=1)),
        (InvalidRepairTimingError, lambda ledger: replace(ledger, repair_s=-1e-3)),
        (InvalidRepairTimingError, lambda ledger: replace(ledger, certificate_s=math.nan)),
        (InvalidRepairTimingError, lambda ledger: replace(ledger, repair_s=math.inf)),
    ]


@pytest.mark.parametrize(("error", "variant"), _ledger_variants())
def test_ledger_build_rejects_each_broken_invariant(error: type[CoverageRepairInvariantError], variant: Callable[[CoverageRepairLedger], CoverageRepairLedger]) -> None:
    ledger = _certified().ledger
    assert ledger.rounds == 1 and ledger.emitted_after == ledger.emitted_before + 1
    _raises_invariant(error, lambda: _rebuild_ledger(variant(ledger)))


def test_certified_path_build_rejects_each_broken_invariant() -> None:
    certified = _certified()
    ledger = certified.ledger
    _raises_invariant(CertifiedSelectionRangeError, lambda: _rebuild_path(replace(certified, sources=certified.sources[:3])))
    _raises_invariant(CertifiedLedgerMismatchError, lambda: _rebuild_path(replace(certified, ledger=replace(ledger, emitted_after=ledger.emitted_after + 1))))
    stranger = replace(ledger.insertions[0], source=9)
    _raises_invariant(CertifiedLedgerMismatchError, lambda: _rebuild_path(replace(certified, ledger=replace(ledger, insertions=(stranger,)))))
    exposed = replace(certified.certificate, exposed_edge_count=1)
    _raises_invariant(UncertifiedPathError, lambda: _rebuild_path(replace(certified, certificate=exposed)))


def test_raw_constructors_stay_unvalidated() -> None:
    """Bypass-safe: the dataclass constructors accept what ``build`` rejects."""
    assert EngagementRefinementSettings(TOOL, FAN_CAP, -1, 0).max_circles == 0
    assert BoundaryOrderedEmission((2, 1), _certified().emission.refinement).selected == (2, 1)
