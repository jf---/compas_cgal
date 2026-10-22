"""The exact overlap margin is compared with 0 only, equality and truthiness included; the selected-successor guard is reachable."""

import math

import pytest

from benchmarks.held_standard_placement import OverlapMarginComparisonError
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import predecessor_overlap_margin
from benchmarks.held_standard_placement import select_next_standard_candidate
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)


def _paper(x: float, radius: float) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=Point2[WorldXY].build(x, 0.0),
        guide_radius=GuideRadius.build(radius),
        contact_point=Point2[WorldXY].build(x, -radius),
    )


# 2r - d - rho + rho_prev = 2 - 2 - 1 + 1 = 0 exactly.
ZERO = predecessor_overlap_margin(_paper(0.0, 1.0), _paper(2.0, 1.0), TOOL)
POSITIVE = predecessor_overlap_margin(_paper(0.0, 1.0), _paper(1.0, 1.0), TOOL)


def test_exact_zero_margin_equals_zero() -> None:
    assert ZERO.sign() == 0
    assert ZERO == 0
    assert not (ZERO != 0)
    assert POSITIVE != 0
    assert not (POSITIVE == 0)


def test_margin_equality_with_anything_but_zero_raises() -> None:
    with pytest.raises(OverlapMarginComparisonError):
        _ = ZERO == 0.0
    with pytest.raises(OverlapMarginComparisonError):
        _ = POSITIVE != 1


def test_margin_truthiness_raises() -> None:
    with pytest.raises(OverlapMarginComparisonError):
        bool(ZERO)
    with pytest.raises(OverlapMarginComparisonError):
        bool(POSITIVE)


def test_selected_successor_guard_refuses_an_exact_zero_margin_at_a_half_turn_cap() -> None:
    # At cap = pi the full-slot engagement outside the overlap domain is within
    # the cap, so the candidate is selected (not forced) and the selected guard decides.
    with pytest.raises(StandardPlacementFragmentationError, match="Selected successor"):
        select_next_standard_candidate(_paper(0.0, 1.0), (_paper(2.0, 1.0),), TOOL, EngagementCap.build(math.pi))
