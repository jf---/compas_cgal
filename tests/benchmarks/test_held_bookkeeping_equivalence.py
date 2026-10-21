"""The Held bookkeeping fast paths decide and report exactly what the per-call code did.

Each oracle below is the pre-change implementation, kept verbatim, so these
properties pin bit-identical outputs (floats compared with ``==``, Fractions
exactly), not approximate agreement. The placement oracle keeps the pre-change
per-step structure with the decision rule injected (``_DecisionRule``): both
sides decide with the same live exact rule, so only the re-slicing is tested.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from types import SimpleNamespace
from typing import Sequence
from typing import TypeVar

from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks import held_standard_placement as placement_module
from benchmarks.held_figure5_path import _boundary_neighborhood_owners
from benchmarks.held_figure5_path import _cyclic_boundary_neighborhoods
from benchmarks.held_figure5_path import _exact_squared_distance
from benchmarks.held_figure5_raw_guide import Figure5BoundarySiteBudget
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import RationalSegment
from benchmarks.held_figure5_raw_guide import _boundary_arclength
from benchmarks.held_figure5_raw_guide import _boundary_segments
from benchmarks.held_figure5_raw_guide import _rational_point
from benchmarks.held_standard_placement import OrderedPlacementSpan
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import PredecessorEngagementCosine
from benchmarks.held_standard_placement import StandardPlacement
from benchmarks.held_standard_placement import select_next_in_span
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL_RADIUS = ToolRadius.build(1.0)
T = TypeVar("T")

# Quarter-unit lattice values make exact ties common; generic floats cover the rest.
_LATTICE = st.integers(min_value=-16, max_value=16).map(lambda value: value / 4)
_GENERIC = st.floats(min_value=-4.0, max_value=4.0, allow_nan=False, allow_infinity=False)
_COORDINATE = st.one_of(_LATTICE, _GENERIC)


# --- oracles: the pre-change implementations, verbatim ----------------------


def _old_boundary_segments(guide: Figure5RawGuide) -> tuple[RationalSegment, ...]:
    points = tuple(_rational_point(point) for point in guide.case.projection.points)
    return tuple(zip(points, (*points[1:], points[0])))


def _old_boundary_arclength(guide: Figure5RawGuide, site: ProjectionBoundarySite) -> Millimetre:
    points = guide.case.projection.points
    lengths = tuple(
        math.hypot(
            float(end.x) - float(start.x),
            float(end.y) - float(start.y),
        )
        for start, end in zip(points, (*points[1:], points[0]))
    )
    ordinal = int(site.segment_id)
    arclength = sum(lengths[:ordinal]) + float(site.parameter) * lengths[ordinal]
    return Millimetre(arclength % float(guide.boundary_length))


def _old_boundary_neighborhood_owners(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    segment_count: int,
    published_start_center: Point2[WorldXY],
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    by_station: dict[
        tuple[GuideRunId, GuideRunStationOrdinal],
        list[ProjectionAdmissibleBoundaryHypothesis],
    ] = {}
    for hypothesis in hypotheses:
        by_station.setdefault(
            (hypothesis.run_id, hypothesis.station_ordinal),
            [],
        ).append(hypothesis)

    owners: list[ProjectionAdmissibleBoundaryHypothesis] = []
    for station_hypotheses in by_station.values():
        neighborhoods = _cyclic_boundary_neighborhoods(
            {int(hypothesis.projected_boundary_site.segment_id) for hypothesis in station_hypotheses},
            segment_count,
        )
        for neighborhood in neighborhoods:
            candidates = tuple(hypothesis for hypothesis in station_hypotheses if int(hypothesis.projected_boundary_site.segment_id) in neighborhood)
            distances = {
                hypothesis: _exact_squared_distance(
                    hypothesis.middle_point,
                    hypothesis.boundary_footpoint,
                )
                for hypothesis in candidates
            }
            minimum = min(distances.values())
            owners.extend(hypothesis for hypothesis in candidates if distances[hypothesis] == minimum)

    start_distances = {
        hypothesis: _exact_squared_distance(
            published_start_center,
            hypothesis.contact_point,
        )
        for hypothesis in hypotheses
    }
    start_minimum = min(start_distances.values())
    owners.extend(hypothesis for hypothesis in hypotheses if start_distances[hypothesis] == start_minimum and hypothesis not in owners)
    return tuple(owners)


def _old_validate_ordered_progress(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    tool_radius: ToolRadius,
) -> tuple[float, ...]:
    if not candidates:
        raise placement_module.InvalidStandardPlacementInputError("Standard placement requires at least one successor candidate.")
    progress = []
    previous = predecessor
    total = 0.0
    for candidate in candidates:
        if type(candidate) is not PaperCircleCandidate:
            raise placement_module.InvalidStandardPlacementInputError("Standard placement candidates must be typed paper circles.")
        spacing = placement_module._center_spacing(previous, candidate)
        if spacing <= placement_module._resolvability_floor(previous, candidate, tool_radius):
            raise placement_module.StandardPlacementResolutionError("Ordered paper candidates contain a zero or binary64-unresolvable restart.")
        total += spacing
        progress.append(total)
        previous = candidate
    return tuple(progress)


@dataclass(frozen=True)
class _DecisionRule:
    """The bisection's decision rule, injected into the per-step oracle.

    The bookkeeping change under test is the re-slicing (one validation per
    span, cursor offsets, progress slices), not the rule. Both sides run the
    same live exact rule, so a rule change can never make this property a
    hybrid that never existed; the rule itself is pinned against an
    independent reference in ``test_held_exact_cap_decision.py``.
    """

    engagement: Callable[[PaperCircleCandidate, PaperCircleCandidate, ToolRadius], PredecessorEngagementCosine]
    within: Callable[[PredecessorEngagementCosine, EngagementCap], bool]
    stops: Callable[[PredecessorEngagementCosine, EngagementCap], bool]
    report: Callable[[PredecessorEngagementCosine], Radian]


LIVE_RULE = _DecisionRule(
    engagement=placement_module.exact_predecessor_engagement,
    within=lambda cosine, cap: cosine.within(cap),
    stops=lambda cosine, cap: placement_module._BisectionStopBand.build(cap).stops(cosine),
    report=lambda cosine: cosine.reported_angle(),
)


def _old_select_next_standard_candidate(
    predecessor: PaperCircleCandidate,
    candidates: Sequence[PaperCircleCandidate],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    rule: _DecisionRule,
) -> StandardPlacement:
    """Pre-bookkeeping ``select_next_standard_candidate`` verbatim; only the four rule calls are injected."""
    progress = _old_validate_ordered_progress(predecessor, candidates, tool_radius)
    cap_angle = float(cap.theta)
    low = 0
    high = len(candidates) - 1
    best: tuple[int, PredecessorEngagementCosine] | None = None
    while low <= high:
        middle = (low + high) // 2
        engagement = rule.engagement(
            predecessor,
            candidates[middle],
            tool_radius,
        )
        if rule.within(engagement, cap):
            best = middle, engagement
            if rule.stops(engagement, cap):
                break
            low = middle + 1
        else:
            high = middle - 1

    if best is None:
        forced = True
        first_margin = placement_module._overlap_margin(predecessor, candidates[0], tool_radius)
        if first_margin <= 0.0:
            raise placement_module.StandardPlacementFragmentationError("The first successor exceeds the paper swept-disk overlap bound.")
        selected_index = 0
        engagement = rule.engagement(
            predecessor,
            candidates[selected_index],
            tool_radius,
        )
    else:
        forced = False
        selected_index, engagement = best

    reported = rule.report(engagement)
    reached_terminal = selected_index == len(candidates) - 1
    cap_shortfall = max(0.0, cap_angle - float(reported))
    cap_excess = max(0.0, float(reported) - cap_angle)
    candidate = candidates[selected_index]
    spacing = placement_module._center_spacing(predecessor, candidate)
    margin = placement_module._overlap_margin(predecessor, candidate, tool_radius)
    if margin <= 0.0:
        raise placement_module.StandardPlacementFragmentationError("Selected successor violates the paper swept-disk overlap bound.")
    return StandardPlacement(
        candidate=candidate,
        candidate_advance=selected_index + 1,
        center_spacing=Millimetre(spacing),
        guide_progress=Millimetre(progress[selected_index]),
        maximum_engagement=reported,
        cap_shortfall=Radian(cap_shortfall),
        cap_excess=Radian(cap_excess),
        overlap_margin=Millimetre(margin),
        reached_terminal=reached_terminal,
        forced=forced,
    )


# --- raw guide: cached segments and arclength offsets ------------------------


def _ring_without_zero_length_segments(points: list[tuple[float, float]]) -> bool:
    return all(start != end for start, end in zip(points, (*points[1:], points[0])))


def _synthetic_guide(ring: list[tuple[float, float]]) -> Figure5RawGuide:
    """Bypass the Held-case factory: segments and arclength read only these fields."""
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ring)
    # The same expression Figure5RawGuide.build_reference uses.
    length = sum(math.hypot(float(end.x) - float(start.x), float(end.y) - float(start.y)) for start, end in zip(boundary, (*boundary[1:], boundary[0])))
    case = SimpleNamespace(projection=SimpleNamespace(points=boundary))
    budget = Figure5BoundarySiteBudget(Millimetre(0.0), Millimetre(0.0), Millimetre(0.0))
    return Figure5RawGuide(case, budget, Millimetre(length), ())  # type: ignore[arg-type]


@settings(max_examples=300, deadline=None)
@given(
    ring=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=3, max_size=40).filter(_ring_without_zero_length_segments),
    sites=st.lists(st.tuples(st.integers(min_value=0, max_value=39), st.fractions(min_value=0, max_value=1, max_denominator=1 << 20)), min_size=1, max_size=12),
)
def test_cached_boundary_bookkeeping_is_bit_identical(ring: list[tuple[float, float]], sites: list[tuple[int, Fraction]]) -> None:
    guide = _synthetic_guide(ring)
    assert _boundary_segments(guide) == _old_boundary_segments(guide)
    for segment, parameter in sites:
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(segment % len(ring)), parameter, len(ring))
        new = _boundary_arclength(guide, site)
        old = _old_boundary_arclength(guide, site)
        assert new == old, (site, new, old)


# --- path: neighbourhood ownership without hypothesis-keyed dicts ------------

_SQUARE_SIDES = 4
_DISTANCES = st.one_of(st.sampled_from((0.5, 1.0, 1.5, 2.0)), st.floats(min_value=0.25, max_value=2.0, allow_nan=False))


@st.composite
def _ownership_hypothesis(draw: st.DrawFn) -> ProjectionAdmissibleBoundaryHypothesis:
    side = draw(st.integers(min_value=0, max_value=_SQUARE_SIDES - 1))
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), draw(st.sampled_from((Fraction(1, 4), Fraction(1, 2)))), _SQUARE_SIDES)
    contact = Point2[WorldXY].build(draw(_LATTICE), draw(_LATTICE))
    middle = Point2[WorldXY].build(2.0, 2.0)
    footpoint = Point2[WorldXY].build(2.0 + draw(_DISTANCES), 2.0)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(draw(st.integers(min_value=0, max_value=2))),
        GuideRunStationOrdinal(draw(st.integers(min_value=0, max_value=2))),
        site,
        (site,),
        Millimetre(0.0),
        middle,
        footpoint,
        contact,
        contact,
        GuideRadius.build(1.0),
    )


@settings(max_examples=400, deadline=None)
@given(
    pool=st.lists(_ownership_hypothesis(), min_size=1, max_size=8),
    picks=st.lists(st.integers(min_value=0, max_value=7), min_size=1, max_size=24),
    start=st.tuples(_LATTICE, _LATTICE),
)
def test_neighbourhood_owners_are_identical_objects_in_identical_order(
    pool: list[ProjectionAdmissibleBoundaryHypothesis],
    picks: list[int],
    start: tuple[float, float],
) -> None:
    # Re-picking from a small pool makes equal and identical duplicates, the
    # cases a hypothesis-keyed dict and list membership treat specially.
    hypotheses = tuple(pool[pick % len(pool)] for pick in picks)
    center = Point2[WorldXY].build(*start)
    new = _boundary_neighborhood_owners(hypotheses, _SQUARE_SIDES, center)
    old = _old_boundary_neighborhood_owners(hypotheses, _SQUARE_SIDES, center)
    assert len(new) == len(old)
    assert all(first is second for first, second in zip(new, old, strict=True))


# --- placement: one validation per stop span ---------------------------------


def _candidate(x: float, y: float, radius: float) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=Point2[WorldXY].build(x, y),
        guide_radius=GuideRadius.build(radius),
        contact_point=Point2[WorldXY].build(x, y - radius),
    )


@st.composite
def _lane(draw: st.DrawFn) -> tuple[PaperCircleCandidate, ...]:
    count = draw(st.integers(min_value=2, max_value=14))
    # Steps include zero and one-ulp restarts so the resolvability check fires.
    steps = draw(st.lists(st.one_of(st.floats(min_value=0.02, max_value=1.6, allow_nan=False), st.just(0.0), st.just(5e-324)), min_size=count, max_size=count))
    radii = draw(st.lists(st.one_of(st.just(1.0), st.floats(min_value=0.2, max_value=2.0, allow_nan=False)), min_size=count, max_size=count))
    wobble = draw(st.lists(st.sampled_from((0.0, 0.0, 0.1, -0.05)), min_size=count, max_size=count))
    x = 0.0
    lane = []
    for step, radius, dy in zip(steps, radii, wobble, strict=True):
        x += step
        lane.append(_candidate(x, dy, radius))
    return tuple(lane)


def _outcome(call: Callable[[], T]) -> T | tuple[type[Exception], str]:
    try:
        return call()
    except (placement_module.StandardPlacementResolutionError, placement_module.StandardPlacementFragmentationError) as error:
        return type(error), str(error)


def _drive_old(candidates: tuple[PaperCircleCandidate, ...], stops: list[int], cap: EngagementCap) -> list[object]:
    """The caller's pre-change loop: re-slice and re-validate on every step."""
    events: list[object] = []
    cursor = 0
    for stop in stops:
        while cursor < stop:
            outcome = _outcome(lambda: _old_select_next_standard_candidate(candidates[cursor], candidates[cursor + 1 : stop + 1], TOOL_RADIUS, cap, LIVE_RULE))
            events.append(outcome)
            if not isinstance(outcome, StandardPlacement):
                return events
            cursor += outcome.candidate_advance
    return events


def _drive_new(candidates: tuple[PaperCircleCandidate, ...], stops: list[int], cap: EngagementCap) -> list[object]:
    """The caller's current loop in ``build_hypothesis_figure5_path``."""
    events: list[object] = []
    cursor = 0
    for stop in stops:
        if cursor >= stop:
            continue
        span_start = cursor
        span = _outcome(lambda: OrderedPlacementSpan.build(candidates[span_start : stop + 1], TOOL_RADIUS, cap))
        if not isinstance(span, OrderedPlacementSpan):
            events.append(span)
            return events
        while cursor < stop:
            outcome = _outcome(lambda: select_next_in_span(span, cursor - span_start))
            events.append(outcome)
            if not isinstance(outcome, StandardPlacement):
                return events
            cursor += outcome.candidate_advance
    return events


# Any cap in (0, pi]: the paper's caps, plus generic doubles that land a lane
# pair within an ulp of the cap or the band floor.
_CAP_THETA = st.one_of(st.sampled_from((30.0, 60.0, 80.0, 120.0, 160.0)).map(math.radians), st.floats(min_value=0.01, max_value=math.pi, exclude_min=True))


@settings(max_examples=500, deadline=None)
@given(
    lane=_lane(),
    stop_marks=st.lists(st.integers(min_value=1, max_value=13), max_size=4),
    cap_theta=_CAP_THETA,
)
# Review C1.1: a float-rule oracle and the exact rule straddle this cap.
@example(lane=(_candidate(0.0, 0.0, 1.0), _candidate(0.5315897334612483, 0.0, 1.0)), stop_marks=[1], cap_theta=float.fromhex("0x1.7e1b2f6432080p+0"))
def test_span_placement_matches_per_step_placement_bit_for_bit(
    lane: tuple[PaperCircleCandidate, ...],
    stop_marks: list[int],
    cap_theta: float,
) -> None:
    last = len(lane) - 1
    stops = sorted({mark for mark in stop_marks if mark < last}) + [last]
    cap = EngagementCap.build(cap_theta)
    assert _drive_new(lane, stops, cap) == _drive_old(lane, stops, cap)
