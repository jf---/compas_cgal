import math

import pytest

from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import StandardPlacementResolutionError
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from benchmarks.held_standard_placement import select_next_standard_candidate
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


TOOL_RADIUS = ToolRadius.build(1.0)


def _candidate(x: float, radius: float = 1.0) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=Point2[WorldXY].build(x, 0.0),
        guide_radius=GuideRadius.build(radius),
        contact_point=Point2[WorldXY].build(x, -radius),
    )


def test_predecessor_engagement_matches_paper_equation_7_on_actual_circles() -> None:
    engagement = maximum_predecessor_engagement(
        _candidate(0.0),
        _candidate(1.0),
        TOOL_RADIUS,
    )

    assert float(engagement) == pytest.approx(2.0 * math.pi / 3.0)


def test_predecessor_engagement_uses_paper_overlap_correction_branch() -> None:
    engagement = maximum_predecessor_engagement(
        _candidate(0.0, radius=0.1),
        _candidate(0.15, radius=0.1),
        TOOL_RADIUS,
    )

    assert float(engagement) == pytest.approx(2.060928921304565)


def test_bisection_selects_cap_target_on_transformed_candidate_circle() -> None:
    predecessor = _candidate(0.0)
    target_spacing = 0.5857864376269049
    candidates = (
        _candidate(0.25),
        _candidate(target_spacing),
        _candidate(0.75),
        _candidate(1.0),
    )

    placement = select_next_standard_candidate(
        predecessor,
        candidates,
        TOOL_RADIUS,
        EngagementCap.build(math.pi / 2.0),
    )

    assert placement.candidate == candidates[1]
    assert placement.candidate_advance == 2
    assert float(placement.center_spacing) == pytest.approx(target_spacing)
    assert float(placement.guide_progress) == pytest.approx(target_spacing)
    assert float(placement.maximum_engagement) == pytest.approx(math.pi / 2.0)
    assert 0.0 <= float(placement.cap_shortfall) <= 0.001
    assert float(placement.overlap_margin) == pytest.approx(2.0 - target_spacing)


def test_terminal_candidate_is_kept_when_whole_remainder_is_below_cap() -> None:
    candidates = (_candidate(0.1), _candidate(0.2), _candidate(0.3))

    placement = select_next_standard_candidate(
        _candidate(0.0),
        candidates,
        TOOL_RADIUS,
        EngagementCap.build(math.radians(80.0)),
    )

    assert placement.candidate == candidates[-1]
    assert placement.candidate_advance == 3
    assert placement.reached_terminal


def test_discrete_grid_keeps_furthest_under_cap_and_reports_shortfall() -> None:
    candidates = (_candidate(0.25), _candidate(0.5), _candidate(0.75))

    placement = select_next_standard_candidate(
        _candidate(0.0),
        candidates,
        TOOL_RADIUS,
        EngagementCap.build(math.radians(80.0)),
    )

    assert placement.candidate == candidates[0]
    assert not placement.reached_terminal
    assert float(placement.cap_shortfall) > 0.001


def test_coarse_grid_keeps_first_over_cap_successor_and_reports_excess() -> None:
    candidate = _candidate(0.75)

    placement = select_next_standard_candidate(
        _candidate(0.0),
        (candidate,),
        TOOL_RADIUS,
        EngagementCap.build(math.radians(80.0)),
    )

    assert placement.candidate == candidate
    assert placement.candidate_advance == 1
    assert float(placement.cap_shortfall) == 0.0
    assert float(placement.cap_excess) > 0.0
    assert float(placement.overlap_margin) > 0.0
    assert placement.forced


def test_radius_change_controls_overlap_without_forcing_uniform_center_spacing() -> None:
    placement = select_next_standard_candidate(
        _candidate(0.0, radius=2.0),
        (_candidate(2.5, radius=1.0),),
        TOOL_RADIUS,
        EngagementCap.build(math.radians(160.0)),
    )

    assert float(placement.center_spacing) == pytest.approx(2.5)
    assert float(placement.guide_progress) == pytest.approx(2.5)
    assert float(placement.overlap_margin) == pytest.approx(0.5)
    assert float(placement.maximum_engagement) == pytest.approx(math.radians(151.04497562814015))


def test_non_resolvable_candidate_restart_fails_loudly() -> None:
    predecessor = _candidate(1.0)
    next_x = math.nextafter(1.0, math.inf)

    with pytest.raises(StandardPlacementResolutionError):
        select_next_standard_candidate(
            predecessor,
            (_candidate(next_x),),
            TOOL_RADIUS,
            EngagementCap.build(math.radians(80.0)),
        )


def test_gap_beyond_paper_overlap_bound_fails_loudly() -> None:
    with pytest.raises(StandardPlacementFragmentationError):
        select_next_standard_candidate(
            _candidate(0.0),
            (_candidate(3.0),),
            TOOL_RADIUS,
            EngagementCap.build(math.radians(80.0)),
        )
