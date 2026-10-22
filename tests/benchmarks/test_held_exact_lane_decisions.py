"""Lane split, lane join, fold and placement guards decide exactly, never on binary64.

Every fixture is chosen where the binary64 evaluation of the paper's Eq. 4
overlap margin ``2r - d - rho + rho_prev`` and its exact sign disagree, or
where two centres are distinct but closer than one ulp of their inputs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

import pytest

from benchmarks import held_omega_path
from benchmarks import held_standard_placement
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_path import _continuous_placement_lanes
from benchmarks.held_figure5_path import _PlacementLane
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import StandardPlacementResolutionError
from benchmarks.held_standard_placement import predecessor_overlap_margin
from benchmarks.held_standard_placement import select_next_in_span
from benchmarks.held_standard_placement import select_next_standard_candidate
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
CAP = EngagementCap.build(math.radians(80.0))

# (candidate x, candidate y, predecessor radius, candidate radius); predecessor at the origin, tool radius 1.
# Found by a search over circles at spacing ~ 2r - rho + rho_prev: binary64 and exact signs disagree.
EXACT_POSITIVE_FLOAT_ZERO = (-0.0636884872332801, 1.9989856869407385, 1.0, 1.0)
EXACT_POSITIVE_FLOAT_NEGATIVE = (-2.462649829701943, 1.1069579107938075, 1.0, 0.3)
EXACT_NEGATIVE_FLOAT_ZERO = (1.6285124375410025, -1.161011300881375, 1.0, 1.0)
EXACT_NEGATIVE_FLOAT_POSITIVE = (-1.3417931672865866, 1.9898721306212193, 1.3, 0.9)
DISAGREEMENTS = (
    (EXACT_POSITIVE_FLOAT_ZERO, 1),
    (EXACT_POSITIVE_FLOAT_NEGATIVE, 1),
    (EXACT_NEGATIVE_FLOAT_ZERO, -1),
    (EXACT_NEGATIVE_FLOAT_POSITIVE, -1),
)


def _paper(x: float, y: float, radius: float) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=Point2[WorldXY].build(x, y),
        guide_radius=GuideRadius.build(radius),
        contact_point=Point2[WorldXY].build(x + radius, y),
    )


def _pair(fixture: tuple[float, float, float, float]) -> tuple[PaperCircleCandidate, PaperCircleCandidate]:
    x, y, previous_radius, radius = fixture
    return _paper(0.0, 0.0, previous_radius), _paper(x, y, radius)


def _reference_sign(fixture: tuple[float, float, float, float]) -> int:
    """Independent exact sign of 2r - sqrt(x^2 + y^2) - rho + rho_prev over Fractions."""
    x, y, previous_radius, radius = fixture
    rational = 2 * Fraction(float(TOOL.value)) - Fraction(radius) + Fraction(previous_radius)
    squared_spacing = Fraction(x) ** 2 + Fraction(y) ** 2
    if rational <= 0:
        return -1
    difference = rational * rational - squared_spacing
    return (difference > 0) - (difference < 0)


def _binary64_sign(fixture: tuple[float, float, float, float]) -> int:
    x, y, previous_radius, radius = fixture
    margin = 2.0 * float(TOOL.value) - math.dist((0.0, 0.0), (x, y)) - radius + previous_radius
    return (margin > 0.0) - (margin < 0.0)


@pytest.mark.parametrize(("fixture", "expected"), DISAGREEMENTS)
def test_fixture_is_a_binary64_disagreement(fixture: tuple[float, float, float, float], expected: int) -> None:
    assert _reference_sign(fixture) == expected
    assert _binary64_sign(fixture) != expected


@pytest.mark.parametrize(("fixture", "expected"), DISAGREEMENTS)
def test_overlap_margin_compares_with_zero_exactly(fixture: tuple[float, float, float, float], expected: int) -> None:
    margin = predecessor_overlap_margin(*_pair(fixture), TOOL)
    assert (margin > 0) == (expected > 0)
    assert (margin <= 0) == (expected <= 0)
    assert margin.sign() == expected


@dataclass(frozen=True)
class _Circle:
    center: Point2[WorldXY]
    radius: GuideRadius

    def paper_candidate(self) -> PaperCircleCandidate:
        return PaperCircleCandidate.build(
            center=self.center,
            guide_radius=self.radius,
            contact_point=Point2[WorldXY].build(float(self.center.x) + float(self.radius.value), float(self.center.y)),
        )


def _circle(x: float, y: float, radius: float) -> _Circle:
    return _Circle(Point2[WorldXY].build(x, y), GuideRadius.build(radius))


def _lane_of(fixture: tuple[float, float, float, float]) -> _PlacementLane[_Circle]:
    x, y, previous_radius, radius = fixture
    return _PlacementLane(0, 1, (_circle(0.0, 0.0, previous_radius), _circle(x, y, radius)))


@pytest.mark.parametrize(("fixture", "expected"), DISAGREEMENTS)
def test_omega_lane_split_decides_the_overlap_margin_exactly(fixture: tuple[float, float, float, float], expected: int) -> None:
    lanes = held_omega_path._split_at_gaps((_lane_of(fixture),), TOOL)  # type: ignore[arg-type]
    assert len(lanes) == (1 if expected > 0 else 2)


def _site(parameter: Fraction) -> ProjectionBoundarySite:
    return ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), parameter, 1)


def _polygon_circle(run: int, station: int, parameter: Fraction, x: float, y: float, radius: float) -> Figure5CounterclockwiseCircle:
    site = _site(parameter)
    return Figure5CounterclockwiseCircle(
        GuideRunId(run),
        GuideRunStationOrdinal(station),
        site,
        site,
        Point2[WorldXY].build(x, y),
        Point2[WorldXY].build(x + radius, y),
        GuideRadius.build(radius),
    )


@pytest.mark.parametrize(("fixture", "expected"), DISAGREEMENTS)
def test_polygon_lane_joiner_decides_the_overlap_margin_exactly(fixture: tuple[float, float, float, float], expected: int) -> None:
    x, y, previous_radius, radius = fixture
    # Run 0 ends at the origin, run 1 starts at the fixture centre; both advance along segment 0.
    circles = (
        _polygon_circle(0, 0, Fraction(0), -0.5, 0.0, previous_radius),
        _polygon_circle(0, 1, Fraction(1, 8), 0.0, 0.0, previous_radius),
        _polygon_circle(1, 0, Fraction(2, 8), x, y, radius),
        _polygon_circle(1, 1, Fraction(3, 8), x, y + 0.5, radius),
    )
    lanes = _continuous_placement_lanes(circles, TOOL)
    assert len(lanes) == (1 if expected > 0 else 2)


def test_forced_successor_beyond_the_exact_overlap_bound_fails_loudly() -> None:
    # Binary64 reports a positive margin; the exact margin is negative, so the forced fallback must refuse it.
    predecessor, candidate = _pair(EXACT_NEGATIVE_FLOAT_POSITIVE)
    with pytest.raises(StandardPlacementFragmentationError):
        select_next_standard_candidate(predecessor, (candidate,), TOOL, CAP)


def test_exact_span_places_across_a_sub_ulp_restart() -> None:
    first = _paper(1.0, 0.0, 1.0)
    second = _paper(math.nextafter(1.0, math.inf), 0.0, 1.0)
    third = _paper(1.5, 0.0, 1.0)
    span = held_standard_placement.OrderedPlacementSpan.build((first, second, third), TOOL, CAP)
    placement = select_next_in_span(span, 0)
    assert placement.candidate in (second, third)


def test_exact_span_refuses_coincident_centres() -> None:
    first = _paper(1.0, 0.0, 1.0)
    with pytest.raises(StandardPlacementResolutionError):
        held_standard_placement.OrderedPlacementSpan.build((first, _paper(1.0, 0.0, 0.5)), TOOL, CAP)


def test_fold_keeps_sub_ulp_distinct_centres_and_folds_coincident_ones() -> None:
    origin = _circle(1.0, 0.0, 1.0)
    sub_ulp = _circle(math.nextafter(1.0, math.inf), 0.0, 1.0)
    coincident = _circle(math.nextafter(1.0, math.inf), 0.0, 0.5)
    runs_of = {id(origin): (GuideRunId(0),), id(sub_ulp): (GuideRunId(1),), id(coincident): (GuideRunId(2),)}
    lanes = held_omega_path._fold_coincident_centres((_PlacementLane(0, 1, (origin, sub_ulp, coincident)),), runs_of)  # type: ignore[arg-type]
    assert lanes[0].circles == (origin, sub_ulp)
    assert runs_of[id(sub_ulp)] == (GuideRunId(1), GuideRunId(2))
    assert runs_of[id(origin)] == (GuideRunId(0),)
