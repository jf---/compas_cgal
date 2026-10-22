"""Selection and refinement decide the engagement cap exactly, on the chord-ratio surrogate.

Each fixture sits within a few ulps of the cap, where the float ``acos`` angle
and the exact surrogate ``cos theta >= 1 - chord_ratio / 2`` disagree. The
float verdict is asserted as the fixture's precondition; the exact verdict is
what selection and refinement must follow.
"""

import math
from collections.abc import Callable
from fractions import Fraction
from typing import TypeVar

import pytest
from hypothesis import event
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_boundary_order_spacing import select_boundary_ordered_sources
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import ConnectorEngagementResolutionError
from benchmarks.held_figure5_engagement_refinement import midpoint
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import PAPER_BISECTION_ANGLE_TOLERANCE
from benchmarks.held_standard_placement import ExactEngagementCosineRangeError
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from benchmarks.held_standard_placement import float_predecessor_engagement
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from benchmarks.held_standard_placement import select_next_standard_candidate
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.benchmarks.held_engagement_reference_oracle import ReferenceCosineRangeError
from tests.benchmarks.held_engagement_reference_oracle import ReferenceFragmentationError
from tests.benchmarks.held_engagement_reference_oracle import reference_engagement

T = TypeVar("T")

BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (10, 0), (10, 10), (0, 10)))
TOOL = ToolRadius.build(1.0)
PREDECESSOR = (4.0, 2.0)
# Float acos says theta <= cap (equal); the exact surrogate says it exceeds.
FLOAT_WITHIN_EXACT_EXCEEDS = ((5.3206062968093875, 2.0), EngagementCap.build(float.fromhex("0x1.12d36c84bb556p+1")))
# Float acos says theta > cap; the exact surrogate says it is within.
FLOAT_EXCEEDS_EXACT_WITHIN = ((4.378758879224697, 1.25), EngagementCap.build(float.fromhex("0x1.24d0b5d32b17fp-1")))


def _circle(ordinal: int, x: float, radius: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(x) / 10, 4)
    return Figure5CounterclockwiseCircle(
        GuideRunId(0), GuideRunStationOrdinal(ordinal), site, site, Point2[WorldXY].build(x, radius), Point2[WorldXY].build(x, 0), GuideRadius.build(radius)
    )


def _float_within(a: Figure5CounterclockwiseCircle, b: Figure5CounterclockwiseCircle, cap: EngagementCap) -> bool:
    return float(float_predecessor_engagement(_paper_candidate(a), _paper_candidate(b), TOOL)) <= float(cap.theta)


def _exact_within(a: Figure5CounterclockwiseCircle, b: Figure5CounterclockwiseCircle, cap: EngagementCap) -> bool:
    return exact_predecessor_engagement(_paper_candidate(a), _paper_candidate(b), TOOL).within(cap)


def test_fixtures_straddle_the_float_and_exact_verdicts() -> None:
    first = _circle(0, *PREDECESSOR)
    (successor, cap) = FLOAT_WITHIN_EXACT_EXCEEDS
    assert _float_within(first, _circle(1, *successor), cap) and not _exact_within(first, _circle(1, *successor), cap)
    (successor, cap) = FLOAT_EXCEEDS_EXACT_WITHIN
    assert not _float_within(first, _circle(1, *successor), cap) and _exact_within(first, _circle(1, *successor), cap)


def test_selection_retains_the_middle_circle_when_the_exact_engagement_exceeds() -> None:
    (successor, cap) = FLOAT_WITHIN_EXACT_EXCEEDS
    circles = (_circle(0, *PREDECESSOR), _circle(1, (PREDECESSOR[0] + successor[0]) / 2, 2.0), _circle(2, *successor))
    assert select_boundary_ordered_sources(circles, ((GuideRunId(0),),) * 3, TOOL, cap) == (0, 1, 2)


def test_selection_drops_the_middle_circle_when_the_exact_engagement_is_within() -> None:
    (successor, cap) = FLOAT_EXCEEDS_EXACT_WITHIN
    circles = (_circle(0, *PREDECESSOR), _circle(1, (PREDECESSOR[0] + successor[0]) / 2, 1.5), _circle(2, *successor))
    assert select_boundary_ordered_sources(circles, ((GuideRunId(0),),) * 3, TOOL, cap) == (0, 2)


def test_refinement_subdivides_when_the_exact_engagement_exceeds() -> None:
    (successor, cap) = FLOAT_WITHIN_EXACT_EXCEEDS
    result = refine_figure5_engagement((_circle(0, *PREDECESSOR), _circle(1, *successor)), BOUNDARY, TOOL, cap)
    assert len(result.circles) > 2
    assert all(reference_engagement(_paper_candidate(a), _paper_candidate(b), TOOL).within(cap) for a, b in zip(result.circles, result.circles[1:]))


def test_refinement_accepts_the_transition_when_the_exact_engagement_is_within() -> None:
    """The circle pair is accepted on the exact verdict; only its connector asks for a circle.

    Acceptance needs the pair within the cap AND the connector certified. At
    this 32.9 degree cap the connector (4, 0) -> (4.379, 0) leaving the
    rho = 2 circle has tier-1 bound O(d) = 2 pi - 2 acos((d^2 + r^2 - R^2) / (2 d r))
    = 37.3 degrees at its far end (d^2 = 4 + 0.379^2, R = 3), so the gate
    refuses it. With no depth the refusal is therefore the connector's
    (``ConnectorEngagementResolutionError``); a pair refused on the float
    verdict would raise the pair's own depth error instead. The one midpoint
    (x = 4.189, rho = 1.625) halves both connectors to 18.8 and 21.6 degrees.
    """
    (successor, cap) = FLOAT_EXCEEDS_EXACT_WITHIN
    pair = (_circle(0, *PREDECESSOR), _circle(1, *successor))
    with pytest.raises(ConnectorEngagementResolutionError):
        refine_figure5_engagement(pair, BOUNDARY, TOOL, cap, max_depth=0)
    result = refine_figure5_engagement(pair, BOUNDARY, TOOL, cap)
    assert len(result.circles) == 3 and result.original_indices == (0, 2)
    assert result.circles[1] == midpoint(*pair, BOUNDARY)


def _candidate(x: float, y: float, rho: float) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(center=Point2[WorldXY].build(x, y), guide_radius=GuideRadius.build(rho), contact_point=Point2[WorldXY].build(x, y - rho))


def _outcome(call: Callable[[], T]) -> T | str:
    try:
        return call()
    except (StandardPlacementFragmentationError, ReferenceFragmentationError):
        return "no real swept-disk intersection"
    except (ExactEngagementCosineRangeError, ReferenceCosineRangeError):
        return "cosine outside [-1, 1]"


# Quarter-unit lattice draws make exact ties (tangencies, concentric pairs,
# branch boundaries) common; generic doubles cover the rest.
def _length(low: float, high: float) -> st.SearchStrategy[float]:
    return st.one_of(st.integers(int(low * 4), int(high * 4)).map(lambda value: value / 4), st.floats(low, high))


# No float comparison here: the binary64 witness has no proven error bound
# across the model's branch discontinuities (the tangency test below pins one
# where it reports pi), so a float-vs-exact property would need a filter that
# can hide a systematic error. Exact agrees with the independent reference on
# every draw, with no assume.
# Axis-aligned offsets make the centre spacing rational, so exact branch ties
# (b = 0, q_y^2 = 0, W on the predecessor circle) occur; the binary64 pi cap's
# surrogate is exactly 4, the full-slot tie.
@settings(max_examples=1000, deadline=None)
@given(
    st.one_of(st.just((0.0, 0.0)), st.tuples(_length(-2.0, 2.0), st.just(0.0)), st.tuples(_length(-2.0, 2.0), _length(-2.0, 2.0))),
    _length(0.25, 4.0),
    _length(0.25, 4.0),
    _length(0.25, 2.0),
    st.one_of(st.floats(0.05, 3.1), st.just(math.pi)),
)
def test_exact_engagement_equals_the_independent_reference_on_every_draw(offset: tuple[float, float], previous_radius: float, radius: float, tool: float, theta: float) -> None:
    first, second, tool_radius, cap = _candidate(0.0, 0.0, previous_radius), _candidate(*offset, radius), ToolRadius.build(tool), EngagementCap.build(theta)
    exact = _outcome(lambda: exact_predecessor_engagement(first, second, tool_radius))
    reference = _outcome(lambda: reference_engagement(first, second, tool_radius))
    if isinstance(reference, str) or isinstance(exact, str):
        event(f"error: {reference}")
        assert exact == reference
        return
    event(reference.branch)
    assert reference.equals(exact)
    assert exact.within(cap) == reference.within(cap)
    assert exact.reaches(cap) == reference.reaches(cap)


def test_exact_tangency_takes_the_model_branch_the_float_model_rounds_away() -> None:
    # d = rho_prev + rho = 0.5 makes b + rho = r exactly, so q_y^2 = 0: not the
    # model's q_y^2 < 0 full-slot branch. Doubles round q_y^2 below zero.
    tool = ToolRadius.build(0.9005620300088643)
    first, second, cap = _candidate(0.0, 0.0, 0.25), _candidate(0.0, 0.5, 0.25), EngagementCap.build(3.0)
    assert float(float_predecessor_engagement(first, second, tool)) == math.pi
    exact = exact_predecessor_engagement(first, second, tool)
    assert exact.within(cap)
    assert float(exact.reported_angle()) < 3.0


# Section 2.4 bisection fixtures: predecessor (0, rho 1), tool 1, candidates on
# the x axis. The reported angle (exact branch, binary64 acos) sits within an
# ulp of the cap or of the paper's stopping band, where the old float tests
# ``theta <= cap`` and ``cap - theta <= 0.001`` disagree with the exact
# surrogates the bisection must follow.
BISECTION_PREDECESSOR = 0.0
# Reported angle <= cap; the exact surrogate says it exceeds: the lone candidate is forced.
BISECTION_FLOAT_WITHIN_EXACT_EXCEEDS = (0.5315897334612483, EngagementCap.build(float.fromhex("0x1.7e1b2f6432080p+0")))
# Reported angle > cap; the exact surrogate says within: the lone candidate is accepted.
BISECTION_FLOAT_EXCEEDS_EXACT_WITHIN = (0.3220174628065054, EngagementCap.build(float.fromhex("0x1.2696d69646a45p+0")))
# Within the cap; ``cap - reported <= 0.001`` says stop, the exact band floor says
# keep bisecting. The further candidate is within the cap and inside the band.
BISECTION_FLOAT_STOPS_EXACT_CONTINUES = ((0.37467351734612164, 0.375), EngagementCap.build(float.fromhex("0x1.3ec3e89127576p+0")))


def _axis_candidate(x: float) -> PaperCircleCandidate:
    return _candidate(x, 0.0, 1.0)


def _reported(x: float) -> float:
    return float(maximum_predecessor_engagement(_axis_candidate(BISECTION_PREDECESSOR), _axis_candidate(x), ToolRadius.build(1.0)))


def test_bisection_fixtures_straddle_the_float_and_exact_tests() -> None:
    tool = ToolRadius.build(1.0)
    first = _axis_candidate(BISECTION_PREDECESSOR)
    (x, cap) = BISECTION_FLOAT_WITHIN_EXACT_EXCEEDS
    assert _reported(x) <= float(cap.theta) and not exact_predecessor_engagement(first, _axis_candidate(x), tool).within(cap)
    (x, cap) = BISECTION_FLOAT_EXCEEDS_EXACT_WITHIN
    assert _reported(x) > float(cap.theta) and exact_predecessor_engagement(first, _axis_candidate(x), tool).within(cap)
    ((near, far), cap) = BISECTION_FLOAT_STOPS_EXACT_CONTINUES
    floor = EngagementCap.build(float(cap.theta) - float(PAPER_BISECTION_ANGLE_TOLERANCE))
    near_cosine = exact_predecessor_engagement(first, _axis_candidate(near), tool)
    assert near_cosine.within(cap) and float(cap.theta) - _reported(near) <= float(PAPER_BISECTION_ANGLE_TOLERANCE)
    assert not near_cosine.reaches(floor)
    far_cosine = exact_predecessor_engagement(first, _axis_candidate(far), tool)
    assert far_cosine.within(cap) and far_cosine.reaches(floor)


def test_bisection_forces_the_candidate_the_exact_surrogate_rejects() -> None:
    (x, cap) = BISECTION_FLOAT_WITHIN_EXACT_EXCEEDS
    placement = select_next_standard_candidate(_axis_candidate(BISECTION_PREDECESSOR), (_axis_candidate(x),), ToolRadius.build(1.0), cap)
    assert placement.forced


def test_bisection_accepts_the_candidate_the_exact_surrogate_admits() -> None:
    (x, cap) = BISECTION_FLOAT_EXCEEDS_EXACT_WITHIN
    placement = select_next_standard_candidate(_axis_candidate(BISECTION_PREDECESSOR), (_axis_candidate(x),), ToolRadius.build(1.0), cap)
    assert not placement.forced


def test_bisection_stops_only_inside_the_exact_band() -> None:
    ((near, far), cap) = BISECTION_FLOAT_STOPS_EXACT_CONTINUES
    candidates = (_axis_candidate(near), _axis_candidate(far))
    placement = select_next_standard_candidate(_axis_candidate(BISECTION_PREDECESSOR), candidates, ToolRadius.build(1.0), cap)
    assert placement.candidate == candidates[1]
    assert not placement.forced
