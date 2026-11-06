"""Refinement decides the actual-stock engagement of every target (unit U4 of the 2026-11-04 design).

Fixtures: a 40 mm centre-domain square, its part polygon one tool radius outside it, and
circles on side 0 with their contact on the centre-domain wall. ``(rho0, rho1, dx)`` below
are found by ``build/u4/find_fixtures.py``; each states its own preconditions.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from fractions import Fraction

import numpy as np
import pytest
from hypothesis import HealthCheck
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_actual_stock_decision import ActualStockDecision
from benchmarks.held_actual_stock_decision import ActualStockViewMismatchError
from benchmarks.held_actual_stock_decision import EntryCircle
from benchmarks.held_actual_stock_decision import EntryCircleDecisionError
from benchmarks.held_actual_stock_decision import read_disk
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_coverage_repair_splice import ContinuedSplice
from benchmarks.held_coverage_repair_splice import continue_with_early_stop
from benchmarks.held_coverage_repair_splice import emitted_form
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_connector_gate import _ConnectorGate
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_engagement_circles import circle_on_side
from benchmarks.held_figure5_engagement_errors import ActualStockEngagementResolutionError
from benchmarks.held_figure5_engagement_errors import ConnectorEngagementResolutionError
from benchmarks.held_figure5_engagement_errors import Figure5EngagementResolutionError
from benchmarks.held_figure5_engagement_refinement import EngagementContinuation
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementTraceError
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_figure5_engagement_refinement import _Emission
from benchmarks.held_figure5_engagement_refinement import _transition_engagement
from benchmarks.held_figure5_engagement_refinement import continue_figure5_engagement
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_engagement_trace import DecisionSite
from benchmarks.held_figure5_engagement_trace import EngagementFootprint
from benchmarks.held_figure5_engagement_trace import EngagementRefinementKey
from benchmarks.held_figure5_engagement_trace import EngagementTrace
from benchmarks.held_figure5_engagement_trace import ReadDisk
from benchmarks.held_figure5_engagement_trace import _annulus_boxes
from benchmarks.held_figure5_engagement_trace import _read_boxes
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_nothing_circles import remove_nothing_circles
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.pipeline.admissibility import loop_parameters
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tools.actual_stock_profile import CircleSweep
from tests.benchmarks.held_part_fixture import part_around
from tools.actual_stock_profile import ReplayWorkload

SIDE = 40.0
BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((0.0, 0.0), (SIDE, 0.0), (SIDE, SIDE), (0.0, SIDE)))
PART = PartPolygon.build(tuple(Point2[WorldXY].build(x, y) for x, y in ((-1.0, -1.0), (SIDE + 1.0, -1.0), (SIDE + 1.0, SIDE + 1.0), (-1.0, SIDE + 1.0))))
OFFSET_AS_MATERIAL = PartPolygon.build(BOUNDARY)
TOOL = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80.0))
# The pair model and the connector gate accept these; the hole of the first circle is uncut, so C_1 refuses.
UNSEEN_HOLE = ((10.0, 2.0), (10.5, 1.5))
# The pair model accepts these and the connector gate refuses; only the part's wall band refuses on the stock.
WALL_BAND = ((10.0, 1.5), (10.75, 1.0))


def _side_circle(ordinal: int, x: float, rho: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x) / Fraction(SIDE), len(BOUNDARY))
    return Figure5CounterclockwiseCircle(
        GuideRunId(0), GuideRunStationOrdinal(ordinal), site, site, Point2[WorldXY].build(x, rho), Point2[WorldXY].build(x, 0.0), GuideRadius.build(rho)
    )


def _circles(spec: tuple[tuple[float, float], ...]) -> tuple[Figure5CounterclockwiseCircle, ...]:
    return tuple(_side_circle(k, x, rho) for k, (x, rho) in enumerate(spec))


def _pair_and_gate_accept(first: Figure5CounterclockwiseCircle, second: Figure5CounterclockwiseCircle) -> bool:
    gate = _ConnectorGate(BOUNDARY, TOOL, CAP)
    return _transition_engagement(first, second, TOOL, 0).within(CAP) and gate.certifies([first], gate.transition(first, second))


def _oracle(circles: tuple[Figure5CounterclockwiseCircle, ...], part: PartPolygon = PART) -> Last16Oracle:
    sweeps = tuple(CircleSweep(float(c.center.x), float(c.center.y), float(c.radius.value), float(c.contact_point.x), float(c.contact_point.y)) for c in circles)
    return Last16Oracle(ReplayWorkload("fixture", part.stock_rows(), float(TOOL.value), CAP.chord_ratio, sweeps, (), 0.0), sweeps=SweepView.CIRCLES)


def _decisions(refinement: TracedEngagementRefinement) -> list[EngagementFootprint]:
    return [footprint for footprint in refinement.trace.footprints if footprint.read is not None]


def test_unseen_hole_is_refused_so_a_midpoint_is_emitted_and_the_certificate_passes() -> None:
    first, second = _circles(UNSEEN_HOLE)
    assert _pair_and_gate_accept(first, second), "precondition: only the stock decision can refuse"
    assert _oracle((first, second)).certify(1).verdict is Verdict.REFUSED, "precondition: without the decision this emission violates the certificate"
    refined = refine_figure5_engagement((first, second), BOUNDARY, TOOL, CAP, part=PART)
    assert len(refined.circles) > 2
    assert (refined.circles[0], refined.circles[-1]) == (first, second)
    oracle = _oracle(refined.circles)
    assert oracle.certify(0).verdict is Verdict.REFUSED
    for t in range(1, len(refined.circles)):
        assert oracle.certify(t).verdict in (Verdict.CERTIFIED_LAST16, Verdict.CERTIFIED_FULL_WINDOW), t


def test_exhausted_stock_decision_raises_the_named_error() -> None:
    with pytest.raises(ActualStockEngagementResolutionError) as caught:
        refine_figure5_engagement(_circles(UNSEEN_HOLE), BOUNDARY, TOOL, CAP, part=PART, max_depth=0)
    assert caught.value.position == 1
    assert caught.value.seam is False
    assert isinstance(caught.value, Figure5EngagementResolutionError), "callers catching the refinement's refusals still see it"


def test_key_carries_the_material_domain() -> None:
    first, second = _circles(((10.0, 1.0), (10.25, 1.0)))
    refined = refine_figure5_engagement((first, second), BOUNDARY, TOOL, CAP, part=PART)
    assert refined.trace.key.part == PART
    with pytest.raises(Figure5EngagementTraceError):
        continue_figure5_engagement(refined, (_side_circle(2, 10.5, 1.0),), BOUNDARY, TOOL, CAP, part=OFFSET_AS_MATERIAL)


def test_every_emitted_circle_after_the_entry_has_its_own_decision() -> None:
    circles = _circles(((10.0, 1.0), (10.25, 1.0), (10.5, 1.0), (10.75, 1.0)))
    refined = refine_figure5_engagement(circles, BOUNDARY, TOOL, CAP, part=PART)
    assert refined.trace.filter_counters.decisions == len(_decisions(refined)), "the engine ran each decision, at every position"
    decided = {(f.position, f.read.cx, f.read.cy, f.read.rho, f.read.ex, f.read.ey) for f in _decisions(refined)}
    for t, circle in enumerate(refined.circles[1:], start=1):
        assert (t, float(circle.center.x), float(circle.center.y), float(circle.radius.value), float(circle.contact_point.x), float(circle.contact_point.y)) in decided
    assert min(f.position for f in _decisions(refined)) == 1, "the entry circle is never decided"


def test_refused_decisions_leave_a_footprint_and_the_native_count_agrees() -> None:
    refined = refine_figure5_engagement(_circles(UNSEEN_HOLE), BOUNDARY, TOOL, CAP, part=PART)
    footprints = _decisions(refined)
    assert len(footprints) > len(refined.circles) - 1, "some decision refused"
    assert refined.trace.filter_counters.decisions == len(footprints)
    assert refined.trace.filter_counters.refused == len(footprints) - (len(refined.circles) - 1)
    # The refused target's retry reads the same disk at a later position: both are recorded.
    target = refined.circles[-1]
    at_target = [f.position for f in footprints if f.read.cx == float(target.center.x) and f.read.rho == float(target.radius.value)]
    assert len(at_target) >= 2 and at_target[0] == 1


def test_stock_decisions_read_by_reads_the_read_disk_of_a_decision() -> None:
    key = EngagementRefinementKey(BOUNDARY, TOOL, CAP, 16, False, False, PART)
    decision = EngagementFootprint(7, read=ReadDisk(20.0, 5.0, 5.0, 20.0, 0.0), site=DecisionSite.ADVANCE, refused=False)
    trace = EngagementTrace(key, (decision,))
    # The read set is B((20, 5), 6) and B((20, 0), 1); its box is [14, 26] x [-1, 11].
    assert trace.stock_decisions_read_by((_side_circle(0, 30.0, 1.0),), 0) == (), "outer box [28, 32] misses it"
    assert trace.stock_decisions_read_by((_side_circle(0, 27.5, 1.0),), 0) == (decision,), "outer box [25.5, 29.5] meets it"
    assert trace.stock_decisions_read_by((_side_circle(0, 28.0, 1.0),), 0) == (decision,), "closed boxes: a shared edge at x = 26 reads"
    assert trace.stock_decisions_read_by((_side_circle(0, 27.5, 1.0),), 8) == (), "a decision before `since` is not asked"
    assert trace.last_influenced((_side_circle(0, 27.5, 1.0),), 0) == -1, "connector decisions only"


def test_entry_circle_is_exempt_by_name_and_only_the_entry() -> None:
    first = _side_circle(0, 10.0, 1.0)
    decision = ActualStockDecision.open(PART, TOOL, CAP, EntryCircle(first), ())
    assert decision.size == 1
    assert decision.exceeded(ReadDisk(*UNSEEN_HOLE[0], UNSEEN_HOLE[0][1], UNSEEN_HOLE[0][0], 0.0), 1) in (True, False)
    with pytest.raises(EntryCircleDecisionError):
        decision.exceeded(ReadDisk(10.0, 1.0, 1.0, 10.0, 0.0), 0)


def test_a_view_out_of_step_with_the_emission_is_refused() -> None:
    decision = ActualStockDecision.open(PART, TOOL, CAP, EntryCircle(_side_circle(0, 10.0, 1.0)), ())
    with pytest.raises(ActualStockViewMismatchError):
        decision.exceeded(ReadDisk(10.5, 1.0, 1.0, 10.5, 0.0), 2)


def test_every_wired_site_reaches_the_engine() -> None:
    octagon = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (8, 0), (12, 4), (12, 12), (8, 16), (0, 16), (-4, 12), (-4, 4)))
    part = PartPolygon.build(tuple(Point2[WorldXY].build(x, y) for x, y in ((-5, -1), (9, -1), (13, 3), (13, 13), (9, 17), (-1, 17), (-5, 13), (-5, 3))))

    def on_side_0(run: int, parameter: Fraction, x: float) -> Figure5CounterclockwiseCircle:
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), parameter, len(octagon))
        centre, contact = Point2[WorldXY].build(x, 1.0), Point2[WorldXY].build(x, 0.0)
        return Figure5CounterclockwiseCircle(GuideRunId(run), GuideRunStationOrdinal(0), site, site, centre, contact, GuideRadius.build(1.0))

    first, second = on_side_0(0, Fraction(1, 2), 4.0), on_side_0(1, Fraction(17, 32), 4.25)
    closed = refine_figure5_engagement((first, second), octagon, TOOL, CAP, part=part, close_seam=True)
    assert closed.trace.filter_counters.decisions == len(_decisions(closed)), "the engine ran the seam's decisions too"
    assert closed.trace.decisions_at(DecisionSite.ADVANCE) >= 1
    assert closed.trace.decisions_at(DecisionSite.CLOSE_SEAM) >= 1, "the seam repair's circles are decided"
    assert closed.trace.decisions_at(DecisionSite.CONTINUATION) == 0
    open_path = refine_figure5_engagement((first, second), octagon, TOOL, CAP, part=part)
    continued = continue_figure5_engagement(open_path, (on_side_0(2, Fraction(9, 16), 4.5),), octagon, TOOL, CAP, part=part)
    assert continued.trace.decisions_at(DecisionSite.CONTINUATION) >= 1


def test_the_part_rows_reach_the_engine_and_not_the_boundary() -> None:
    prefix = refine_figure5_engagement((_side_circle(0, 10.0, 1.0),), BOUNDARY, TOOL, CAP, part=PART)
    continuation = EngagementContinuation(prefix, 0, BOUNDARY, TOOL, CAP, part=PART, max_depth=16, max_circles=10, corner_approaches=False, close_seam=False)
    assert np.array_equal(continuation.stock.part.stock_rows(), PART.stock_rows())


def _spliced_dropping_a_midpoint() -> tuple[TracedEngagementRefinement, TracedEngagementRefinement, ContinuedSplice, list[set[tuple[float, float]]]]:
    """A splice whose reference emits a midpoint the augmented selection does not: Y is not empty."""
    spec = ((10.01, 1.5), (10.03, 2.0), (10.06, 1.5), (10.08, 1.5), (10.10, 1.5), (10.11, 1.5), (10.85, 1.0), (10.86, 1.0), (10.95, 1.0), (18.0, 1.0), (18.1, 1.0))
    sources = tuple(emitted_form(circle, BOUNDARY) for circle in _circles(spec))
    selected = tuple(k for k in range(len(sources)) if k != 1)
    reference = refine_figure5_engagement(tuple(sources[k] for k in selected), BOUNDARY, TOOL, CAP, part=PART)
    full = refine_figure5_engagement(sources, BOUNDARY, TOOL, CAP, part=PART)
    settings_ = EngagementRefinementSettings.build(TOOL, CAP, part=PART, close_seam=False)
    calls: list[set[tuple[float, float]]] = []
    original = EngagementTrace.stock_decisions_read_by

    def spy(self: EngagementTrace, circles: Sequence[Figure5CounterclockwiseCircle], since: int) -> tuple[EngagementFootprint, ...]:
        calls.append({(float(c.contact_point.x), float(c.radius.value)) for c in circles})
        return original(self, circles, since)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(EngagementTrace, "stock_decisions_read_by", spy)
        spliced = continue_with_early_stop(reference, 1, sources[1], tuple(sources[k] for k in selected[1:]), BOUNDARY, settings_)
    return reference, full, spliced, calls


def test_h4_a_dropped_midpoint_inside_a_later_read_disk_stops_the_adoption() -> None:
    reference, full, spliced, calls = _spliced_dropping_a_midpoint()
    dropped = [circle for circle in reference.circles if circle not in full.circles]
    assert dropped, "precondition: the augmented selection drops a reference circle (Y is not empty)"
    assert spliced.refinement == full
    assert spliced.refinement.trace == full.trace
    assert 0 < spliced.rerefined < len(reference.original_indices) - 1, "the stop was decided, not vacuous"
    after = reference.circles.index(dropped[-1]) + 1
    assert reference.trace.stock_decisions_read_by(tuple(dropped), after), "precondition: a dropped circle's annulus reaches a later decision's read disk"
    y = {(float(c.contact_point.x), float(c.radius.value)) for c in dropped}
    assert any(y <= seen for seen in calls), "the influence test reads the dropped circles (Y), not only the inserted ones (X)"
    assert spliced.refinement.trace.filter_counters.decisions > 0


def test_a_stock_verdict_that_flips_on_the_new_view_stops_the_adoption() -> None:
    spec = ((6.65, 2.5), (6.8875, 1.5), (6.9625, 2.0), (7.5625, 2.5), (7.65, 1.0), (8.7, 2.5))
    sources = tuple(emitted_form(circle, BOUNDARY) for circle in _circles(spec))
    selected = tuple(k for k in range(len(sources)) if k != 1)
    reference = refine_figure5_engagement(tuple(sources[k] for k in selected), BOUNDARY, TOOL, CAP, part=PART)
    full = refine_figure5_engagement(sources, BOUNDARY, TOOL, CAP, part=PART)
    settings_ = EngagementRefinementSettings.build(TOOL, CAP, part=PART, close_seam=False)
    later = tuple(sources[k] for k in selected[1:])
    assert continue_with_early_stop(reference, 1, sources[1], later, BOUNDARY, settings_).refinement == full
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(EngagementContinuation, "last_stock_flip", lambda *_: -1)
        assert continue_with_early_stop(reference, 1, sources[1], later, BOUNDARY, settings_).refinement != full, "precondition: adopting without the re-decision is wrong"


def test_last_stock_flip_is_the_last_verdict_that_differs_and_the_view_is_cut_back() -> None:
    reference = refine_figure5_engagement(_circles(((10.0, 1.0), (10.25, 1.0), (10.5, 1.0), (10.75, 1.0), (11.0, 1.0))), BOUNDARY, TOOL, CAP, part=PART)
    continuation = EngagementContinuation(
        TracedEngagementRefinement(reference.circles[:2], (0, 1), reference.trace),
        3,
        BOUNDARY,
        TOOL,
        CAP,
        part=PART,
        max_depth=16,
        max_circles=10,
        corner_approaches=False,
        close_seam=False,
    )
    size = continuation.stock.size
    # Every decision here passes on its view, so a footprint recorded as refused is a flip.
    decisions = tuple(EngagementFootprint(t, read=read_disk(reference.circles[t]), refused=t in (2, 4)) for t in (2, 3, 4))
    assert continuation.last_stock_flip(reference, 2, decisions) == 4
    assert continuation.last_stock_flip(reference, 2, decisions[:2]) == 2
    assert continuation.last_stock_flip(reference, 2, ()) == -1
    assert continuation.stock.size == size


def test_a_decision_that_a_splice_re_decides_is_at_the_continuation_site() -> None:
    _, _, spliced, _ = _spliced_dropping_a_midpoint()
    assert spliced.refinement.trace.decisions_at(DecisionSite.CONTINUATION) > 0


@st.composite
def _paths_with_nothing_circles(draw: st.DrawFn) -> tuple[Figure5CounterclockwiseCircle, ...]:
    steps = draw(st.lists(st.tuples(st.integers(1, 24), st.sampled_from((1.0, 1.25, 1.5))), min_size=3, max_size=9))
    x = 10.0
    spec: list[tuple[float, float]] = []
    for step, rho in steps:
        x += step / 64.0
        spec.append((x, rho))
    # An exact repeat of an earlier circle is a NOTHING circle by construction.
    repeat = draw(st.integers(1, len(spec) - 1))
    spec.insert(repeat + 1, spec[repeat])
    return _circles(tuple(spec))


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.filter_too_much, HealthCheck.too_slow])
@given(_paths_with_nothing_circles())
def test_h3_a_nothing_removal_leaves_every_later_view_unchanged(circles: tuple[Figure5CounterclockwiseCircle, ...]) -> None:
    pruned = remove_nothing_circles(circles, _paper_candidate, TOOL, CAP, closed=False)
    assume(pruned.removed)
    removed = set(pruned.removed)
    parameters = loop_parameters()
    for t in range(pruned.removed[0] + 1, len(circles)):
        with_all = _stock_2.AnnulusView2.open(PART.stock_rows(), float(TOOL.value))
        without = _stock_2.Stock2(PART.stock_rows(), [])
        for k in range(t):
            read = read_disk(circles[k])
            with_all.append(read.cx, read.cy, read.rho)
            if k not in removed:
                without.subtract_circle_sweep_local(read.cx, read.cy, read.rho, float(TOOL.value))
        read = read_disk(circles[t])
        rows = with_all.circle_matches_stock(without, 0, t, read.cx, read.cy, read.rho, read.ex, read.ey, False, parameters, CAP.chord_ratio, 0.0)
        assert all(disposition and runs for disposition, runs in rows), t


FINITE = st.floats(-1e3, 1e3, allow_nan=False, allow_infinity=False)
RADIUS = st.floats(1e-3, 1e3, allow_nan=False, allow_infinity=False)


@settings(max_examples=200, deadline=None)
@given(FINITE, FINITE, RADIUS, FINITE, FINITE, RADIUS)
def test_read_boxes_contain_the_exact_read_set_and_outer_boxes_the_exact_outer_disk(cx: float, cy: float, rho: float, ex: float, ey: float, tool: float) -> None:
    """A box miss is the proof of ``last_influenced``: the rounded boxes hold the exact ones, over the rationals the doubles are."""
    (min_x, min_y, max_x, max_y), = _read_boxes((ReadDisk(cx, cy, rho, ex, ey),), tool).tolist()
    exact_reach, exact_tool = Fraction(rho) + Fraction(tool), Fraction(tool)
    assert Fraction(min_x) <= min(Fraction(cx) - exact_reach, Fraction(ex) - exact_tool)
    assert Fraction(min_y) <= min(Fraction(cy) - exact_reach, Fraction(ey) - exact_tool)
    assert Fraction(max_x) >= max(Fraction(cx) + exact_reach, Fraction(ex) + exact_tool)
    assert Fraction(max_y) >= max(Fraction(cy) + exact_reach, Fraction(ey) + exact_tool)
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(1, 2), len(BOUNDARY))
    centre, contact = Point2[WorldXY].build(cx, cy), Point2[WorldXY].build(ex, ey)
    circle = Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(0), site, site, centre, contact, GuideRadius.build(rho))
    (low_x, low_y, high_x, high_y), = _annulus_boxes((circle,), tool).tolist()
    assert Fraction(low_x) <= Fraction(cx) - exact_reach and Fraction(low_y) <= Fraction(cy) - exact_reach
    assert Fraction(high_x) >= Fraction(cx) + exact_reach and Fraction(high_y) >= Fraction(cy) + exact_reach


def test_engine_decides_on_part_not_offset() -> None:
    first, second = _circles(WALL_BAND)
    assert _transition_engagement(first, second, TOOL, 0).within(CAP), "precondition: the pair model does not refuse"
    with pytest.raises(ActualStockEngagementResolutionError):
        refine_figure5_engagement((first, second), BOUNDARY, TOOL, CAP, part=PART, max_depth=0)
    with pytest.raises(ConnectorEngagementResolutionError):
        refine_figure5_engagement((first, second), BOUNDARY, TOOL, CAP, part=OFFSET_AS_MATERIAL, max_depth=0)
    target = read_disk(second)
    on_part = ActualStockDecision.open(PART, TOOL, CAP, EntryCircle(first), ())
    on_offset = ActualStockDecision.open(OFFSET_AS_MATERIAL, TOOL, CAP, EntryCircle(first), ())
    assert on_part.exceeded(target, 1)
    assert not on_offset.exceeded(target, 1)


# On the octagon of the connector tests: the seam repair's own circles meet the unseen holes of R = 3 circles.
OCTAGON = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (8, 0), (12, 4), (12, 12), (8, 16), (0, 16), (-4, 12), (-4, 4)))


def _on_octagon(parameter: Fraction, rho: float, run: int) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), parameter, len(OCTAGON))
    contact = _canonical_offset_point(OCTAGON, site)
    raw = Figure5CounterclockwiseCircle(GuideRunId(run), GuideRunStationOrdinal(0), site, site, contact, contact, GuideRadius.build(rho))
    return circle_on_side(raw, site, contact, GuideRadius.build(rho), OCTAGON)


def test_a_seam_circle_the_stock_refuses_gets_a_midpoint() -> None:
    sources = (_on_octagon(Fraction(1, 2), 3.0, 0), _on_octagon(Fraction(17, 32), 3.0, 1))
    part = part_around(OCTAGON, TOOL)
    closed = refine_figure5_engagement(sources, OCTAGON, TOOL, CAP, part=part, close_seam=True)
    original = _Emission.stock_refuses

    def seam_skipped(self: _Emission, target: Figure5CounterclockwiseCircle, *, seam: bool) -> bool:
        return False if seam else original(self, target, seam=seam)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(_Emission, "stock_refuses", seam_skipped)
        without = refine_figure5_engagement(sources, OCTAGON, TOOL, CAP, part=part, close_seam=True)
    seam_from = closed.original_indices[-1] + 1
    assert closed.circles[:seam_from] == without.circles[:seam_from]
    assert len(closed.circles) > len(without.circles), "the conjunct is live at the seam"
    oracle = _oracle(closed.circles, part)
    for t in range(seam_from, len(closed.circles)):
        assert oracle.certify(t).verdict in (Verdict.CERTIFIED_LAST16, Verdict.CERTIFIED_FULL_WINDOW), t
