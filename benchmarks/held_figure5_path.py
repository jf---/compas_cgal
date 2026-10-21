"""Hypothesis-based Figure 5(a) placement and boundary-path integration."""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import replace
from fractions import Fraction
from typing import Any
from typing import Generic
from typing import NamedTuple
from typing import Protocol
from typing import TypeVar

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import MissingFigure5GuideRunError
from benchmarks.held_figure5_boundary_path import build_figure5_boundary_path
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_offset_contact_projection import project_bounded_offset_contacts
from benchmarks.held_reference_contacts import qualify_reference_contact_sites
from benchmarks.held_standard_placement import OrderedPlacementSpan
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacement
from benchmarks.held_standard_placement import predecessor_overlap_margin
from benchmarks.held_standard_placement import select_next_in_span
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


class Figure5PathIntegrationError(BenchmarkError):
    """Base failure while integrating Figure 5 path hypotheses."""


class AmbiguousFigure5PathHypothesisError(Figure5PathIntegrationError):
    """Distinct circle geometry claims one canonical boundary position."""


@dataclass(frozen=True)
class Figure5PlacementRecord:
    """One predecessor placement and every guide run it represents."""

    boundary_site: ProjectionBoundarySite
    source_run_ids: tuple[GuideRunId, ...]
    placement: StandardPlacement


@dataclass(frozen=True)
class HypothesisFigure5Path:
    """Approximate, hypothesis-based Figure 5 path with provenance ledger."""

    path: Figure5BoundaryPath
    placements: tuple[Figure5PlacementRecord, ...]
    circle_sources: tuple[tuple[GuideRunId, ...], ...]
    reached_run_ids: frozenset[GuideRunId]
    candidate_count: int
    canonical_candidate_count: int
    placement_lane_count: int
    rejected_contact_hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...] = ()


@dataclass(frozen=True)
class _CanonicalHypothesis:
    hypothesis: ProjectionAdmissibleBoundaryHypothesis
    source_run_ids: tuple[GuideRunId, ...]


class LaneCircle(Protocol):
    """What lane grouping reads of a placed circle (polygon-Ω or Ω* circles)."""

    @property
    def run_id(self) -> GuideRunId: ...

    @property
    def station_ordinal(self) -> GuideRunStationOrdinal: ...

    @property
    def source_boundary_site(self) -> ProjectionBoundarySite: ...

    @property
    def center(self) -> Point2[WorldXY]: ...

    @property
    def radius(self) -> GuideRadius: ...


LaneCircleT = TypeVar("LaneCircleT", bound=LaneCircle)


@dataclass(frozen=True)
class _PlacementLane(Generic[LaneCircleT]):
    """One monotone source-side sequence safe for predecessor placement."""

    source_segment: int
    direction: int
    circles: tuple[LaneCircleT, ...]


def _site_key(site: ProjectionBoundarySite) -> tuple[int, Fraction]:
    return int(site.segment_id), site.parameter


_CircleGeometryKey = tuple[float, float, float, float, float]


def _circle_geometry_key(
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
) -> _CircleGeometryKey:
    """Exact circle identity: centre, radius and contact as raw doubles.

    Finite doubles compare, order and hash exactly as their ``Fraction``
    values do (``-0.0 == 0.0`` in both), so the doubles are the exact key
    without ``Fraction.from_float`` normalisation.
    """
    return (
        float(hypothesis.center.x),
        float(hypothesis.center.y),
        float(hypothesis.guide_radius.value),
        float(hypothesis.contact_point.x),
        float(hypothesis.contact_point.y),
    )


def _hypothesis_order_key(
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
) -> tuple[_CircleGeometryKey, int, Fraction, int, int]:
    return (
        _circle_geometry_key(hypothesis),
        int(hypothesis.projected_boundary_site.segment_id),
        hypothesis.projected_boundary_site.parameter,
        int(hypothesis.run_id),
        int(hypothesis.station_ordinal),
    )


def _exact_squared_distance(
    first: Point2[WorldXY],
    second: Point2[WorldXY],
) -> Fraction:
    delta_x = Fraction.from_float(float(second.x)) - Fraction.from_float(float(first.x))
    delta_y = Fraction.from_float(float(second.y)) - Fraction.from_float(float(first.y))
    return delta_x * delta_x + delta_y * delta_y


class _DyadicSquaredDistance(NamedTuple):
    """Exact squared distance ``numerator / 4**exponent`` in mm².

    Every finite double is a dyadic rational ``n / 2**k``, so the squared
    distance between two double points is an integer over ``4**k`` for the
    largest coordinate ``k``. Comparison is integer shift and compare: the
    same exact decisions as the ``Fraction`` value, without gcd normalisation.
    """

    numerator: int
    exponent: int


def _dyadic_squared_distance(
    first: Point2[WorldXY],
    second: Point2[WorldXY],
) -> _DyadicSquaredDistance:
    first_x, first_x_scale = float(first.x).as_integer_ratio()
    first_y, first_y_scale = float(first.y).as_integer_ratio()
    second_x, second_x_scale = float(second.x).as_integer_ratio()
    second_y, second_y_scale = float(second.y).as_integer_ratio()
    # All four denominators are powers of two, so the largest is a common one.
    scale = max(first_x_scale, first_y_scale, second_x_scale, second_y_scale)
    delta_x = second_x * (scale // second_x_scale) - first_x * (scale // first_x_scale)
    delta_y = second_y * (scale // second_y_scale) - first_y * (scale // first_y_scale)
    return _DyadicSquaredDistance(delta_x * delta_x + delta_y * delta_y, scale.bit_length() - 1)


def _dyadic_compare(first: _DyadicSquaredDistance, second: _DyadicSquaredDistance) -> int:
    """Exact sign of ``first - second``."""
    shift = 2 * (first.exponent - second.exponent)
    if shift >= 0:
        left, right = first.numerator, second.numerator << shift
    else:
        left, right = first.numerator << -shift, second.numerator
    return (left > right) - (left < right)


def _dyadic_minimum(distances: tuple[_DyadicSquaredDistance, ...]) -> _DyadicSquaredDistance:
    minimum = distances[0]
    for distance in distances[1:]:
        if _dyadic_compare(distance, minimum) < 0:
            minimum = distance
    return minimum


def _cyclic_boundary_neighborhoods(
    segment_ids: set[int],
    segment_count: int,
) -> tuple[tuple[int, ...], ...]:
    starts = sorted(segment_id for segment_id in segment_ids if (segment_id - 1) % segment_count not in segment_ids)
    if not starts:
        return (tuple(sorted(segment_ids)),)
    neighborhoods = []
    for start in starts:
        neighborhood = []
        segment_id = start
        while segment_id in segment_ids:
            neighborhood.append(segment_id)
            segment_id = (segment_id + 1) % segment_count
        neighborhoods.append(tuple(neighborhood))
    return tuple(neighborhoods)


def _boundary_neighborhood_owner_indices(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    segment_count: int,
    published_start_center: Point2[WorldXY],
) -> tuple[int, ...]:
    """Indices of the owning hypotheses, in emission order.

    Indices, not hypotheses, so callers can carry parallel data (offset sites)
    through ownership without hashing any hypothesis.
    """
    by_station: dict[tuple[GuideRunId, GuideRunStationOrdinal], list[int]] = {}
    for index, hypothesis in enumerate(hypotheses):
        by_station.setdefault(
            (hypothesis.run_id, hypothesis.station_ordinal),
            [],
        ).append(index)

    owners: list[int] = []
    for station_indices in by_station.values():
        neighborhoods = _cyclic_boundary_neighborhoods(
            {int(hypotheses[index].projected_boundary_site.segment_id) for index in station_indices},
            segment_count,
        )
        for neighborhood in neighborhoods:
            candidates = tuple(index for index in station_indices if int(hypotheses[index].projected_boundary_site.segment_id) in neighborhood)
            distances = tuple(
                _dyadic_squared_distance(
                    hypotheses[index].middle_point,
                    hypotheses[index].boundary_footpoint,
                )
                for index in candidates
            )
            minimum = _dyadic_minimum(distances)
            owners.extend(index for index, distance in zip(candidates, distances, strict=True) if _dyadic_compare(distance, minimum) == 0)

    start_distances = tuple(
        _dyadic_squared_distance(
            published_start_center,
            hypothesis.contact_point,
        )
        for hypothesis in hypotheses
    )
    start_minimum = _dyadic_minimum(start_distances)
    # Only exact start-distance ties reach the membership scan (one or two per
    # workload), so a list scan beats hashing every owner into a set. The scan
    # is by hypothesis equality, as before: an equal hypothesis already owned
    # suppresses the start owner.
    for index, distance in enumerate(start_distances):
        if _dyadic_compare(distance, start_minimum) == 0 and hypotheses[index] not in (hypotheses[owner] for owner in owners):
            owners.append(index)
    return tuple(owners)


def _boundary_neighborhood_owners(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    segment_count: int,
    published_start_center: Point2[WorldXY],
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    return tuple(hypotheses[index] for index in _boundary_neighborhood_owner_indices(hypotheses, segment_count, published_start_center))


def _canonical_indices(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
) -> tuple[tuple[int, tuple[GuideRunId, ...]], ...]:
    """Representative index and source runs of each distinct circle, in key order."""
    grouped: dict[_CircleGeometryKey, list[int]] = {}
    for index, hypothesis in enumerate(hypotheses):
        grouped.setdefault(_circle_geometry_key(hypothesis), []).append(index)
    return tuple(
        (
            min(group, key=lambda index: _hypothesis_order_key(hypotheses[index])),
            tuple(sorted({hypotheses[index].run_id for index in group}, key=int)),
        )
        for _, group in sorted(grouped.items())
    )


def _canonicalize(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
) -> tuple[_CanonicalHypothesis, ...]:
    return tuple(_CanonicalHypothesis(hypotheses[index], runs) for index, runs in _canonical_indices(hypotheses))


def _paper_candidate(circle: Figure5CounterclockwiseCircle) -> PaperCircleCandidate:
    return PaperCircleCandidate.build(
        center=circle.center,
        guide_radius=circle.radius,
        contact_point=circle.contact_point,
    )


def _origin_key(
    run_id: GuideRunId,
    station_ordinal: GuideRunStationOrdinal,
    source_boundary_site: ProjectionBoundarySite,
) -> tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]:
    return run_id, station_ordinal, source_boundary_site


def _circle_station_key(
    circle: Figure5CounterclockwiseCircle,
) -> tuple[int, Fraction, int, Fraction, float]:
    return (
        int(circle.station_ordinal),
        circle.source_boundary_site.parameter,
        *_site_key(circle.boundary_site),
        float(circle.radius.value),
    )


def _parameter_direction(first: Fraction, second: Fraction) -> int:
    return (second > first) - (second < first)


def _continuation_score(
    tail: LaneCircle,
    head: LaneCircle,
) -> tuple[Fraction, float, float]:
    return (
        abs(head.source_boundary_site.parameter - tail.source_boundary_site.parameter),
        math.dist(
            (float(tail.center.x), float(tail.center.y)),
            (float(head.center.x), float(head.center.y)),
        ),
        abs(float(head.radius.value) - float(tail.radius.value)),
    )


def _run_side_fragment(
    circles: tuple[LaneCircleT, ...],
) -> _PlacementLane[LaneCircleT]:
    if not circles:
        raise Figure5PathIntegrationError("A continuous placement lane cannot be empty.")
    direction = next(
        (
            step
            for first, second in zip(circles, circles[1:], strict=False)
            if (
                step := _parameter_direction(
                    first.source_boundary_site.parameter,
                    second.source_boundary_site.parameter,
                )
            )
            != 0
        ),
        0,
    )
    return _PlacementLane(
        int(circles[0].source_boundary_site.segment_id),
        direction,
        circles,
    )


def _continuous_placement_lanes(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    tool_radius: ToolRadius,
) -> tuple[_PlacementLane[Figure5CounterclockwiseCircle], ...]:
    return continuous_placement_lanes(circles, tool_radius, _circle_station_key, _paper_candidate)


def continuous_placement_lanes(
    circles: tuple[LaneCircleT, ...],
    tool_radius: ToolRadius,
    station_key: Callable[[LaneCircleT], Any],
    paper_candidate: Callable[[LaneCircleT], PaperCircleCandidate],
) -> tuple[_PlacementLane[LaneCircleT], ...]:
    """Group circles into monotone per-run, per-source-side lanes and join compatible lanes.

    ``station_key`` orders the circles of one (run, source side); it must be a
    total order (the polygon path's key ends in the boundary site, the Ω* path's
    in the exact ∂Ω* site key).
    """
    by_run_side: dict[
        tuple[GuideRunId, int],
        list[LaneCircleT],
    ] = {}
    for circle in circles:
        key = circle.run_id, int(circle.source_boundary_site.segment_id)
        by_run_side.setdefault(key, []).append(circle)

    emitted_fragments = tuple(_run_side_fragment(tuple(sorted(by_run_side[key], key=station_key))) for key in sorted(by_run_side, key=lambda item: (int(item[0]), item[1])))
    fragments = tuple(
        sorted(
            emitted_fragments,
            key=lambda fragment: (
                fragment.source_segment,
                fragment.direction,
                fragment.circles[0].source_boundary_site.parameter if fragment.direction >= 0 else -fragment.circles[0].source_boundary_site.parameter,
                float(fragment.circles[0].center.x),
                float(fragment.circles[0].center.y),
                float(fragment.circles[0].radius.value),
            ),
        )
    )
    lanes: list[_PlacementLane[LaneCircleT]] = []
    # Lane indices per (source_segment, direction), in creation order; a merge
    # replaces lanes[index] with the same key, so the index never goes stale.
    lanes_by_key: dict[tuple[int, int], list[int]] = {}
    for fragment in fragments:
        lane_key = fragment.source_segment, fragment.direction
        if fragment.direction == 0:
            lanes_by_key.setdefault(lane_key, []).append(len(lanes))
            lanes.append(fragment)
            continue
        compatible = []
        head = fragment.circles[0]
        for index in lanes_by_key.get(lane_key, ()):
            tail = lanes[index].circles[-1]
            advances = _parameter_direction(
                tail.source_boundary_site.parameter,
                head.source_boundary_site.parameter,
            )
            if (
                tail.run_id != head.run_id
                and advances == fragment.direction
                and float(
                    predecessor_overlap_margin(
                        paper_candidate(tail),
                        paper_candidate(head),
                        tool_radius,
                    )
                )
                > 0.0
            ):
                compatible.append(index)
        if not compatible:
            lanes_by_key.setdefault(lane_key, []).append(len(lanes))
            lanes.append(fragment)
            continue
        scores = {index: _continuation_score(lanes[index].circles[-1], head) for index in compatible}
        best_score = min(scores.values())
        nearest = [index for index, score in scores.items() if score == best_score]
        if len(nearest) != 1:
            lanes_by_key.setdefault(lane_key, []).append(len(lanes))
            lanes.append(fragment)
            continue
        index = nearest[0]
        lane = lanes[index]
        lanes[index] = _PlacementLane(
            lane.source_segment,
            lane.direction,
            lane.circles + fragment.circles,
        )
    return tuple(lanes)


def build_hypothesis_figure5_path(
    *,
    inward_components: tuple[tuple[Point2[WorldXY], ...], ...],
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    reached_run_ids: frozenset[GuideRunId],
    published_start_center: Point2[WorldXY],
    published_start_radius: Millimetre,
    boundary_evidence_bound: Millimetre,
    start_evidence_bound: Millimetre,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    offset_sites: tuple[ProjectionBoundarySite, ...] | None = None,
    preserve_source_runs: bool = False,
) -> HypothesisFigure5Path:
    """Place distinct admissible hypotheses and emit one CCW boundary path.

    ``offset_sites``, when given, is parallel to ``hypotheses``: entry ``i`` is
    the resolved offset contact of hypothesis ``i``.
    """
    if not hypotheses or not reached_run_ids:
        raise Figure5PathIntegrationError("Figure 5 integration requires hypotheses and a reached-run ledger.")
    if len(inward_components) != 1 or not inward_components[0]:
        raise Figure5PathIntegrationError("Figure 5 neighborhood ownership requires one non-empty boundary cycle.")
    if offset_sites is not None and len(offset_sites) != len(hypotheses):
        raise Figure5PathIntegrationError("Offset sites must be parallel to the hypotheses.")
    owner_indices = _boundary_neighborhood_owner_indices(
        hypotheses,
        len(inward_components[0]),
        published_start_center,
    )
    owned_hypotheses = tuple(hypotheses[index] for index in owner_indices)
    canonical_indices = _canonical_indices(owned_hypotheses)
    canonical = tuple(_CanonicalHypothesis(owned_hypotheses[index], runs) for index, runs in canonical_indices)
    # Sites ride along by index; only the hypotheses handed to the boundary
    # path are ever hashed into its lookup.
    canonical_sites = None if offset_sites is None else tuple(offset_sites[owner_indices[index]] for index, _ in canonical_indices)
    by_origin = {
        _origin_key(
            item.hypothesis.run_id,
            item.hypothesis.station_ordinal,
            item.hypothesis.projected_boundary_site,
        ): position
        for position, item in enumerate(canonical)
    }
    representatives = tuple(item.hypothesis for item in canonical)
    representative_runs = frozenset(item.hypothesis.run_id for item in canonical)
    ordered_path = build_figure5_boundary_path(
        inward_components=inward_components,
        side_candidates=representatives,
        reached_run_ids=representative_runs,
        published_start_center=published_start_center,
        published_start_radius=published_start_radius,
        tool_radius=tool_radius,
        boundary_evidence_bound=boundary_evidence_bound,
        start_evidence_bound=start_evidence_bound,
        offset_sites=None if canonical_sites is None else dict(zip(representatives, canonical_sites, strict=True)),
    )
    selected_positions: list[int] = []
    placements = []
    placement_lanes = _continuous_placement_lanes(ordered_path.circles, tool_radius)
    for lane in placement_lanes:
        lane_circles = lane.circles
        lane_positions = tuple(
            by_origin[
                _origin_key(
                    circle.run_id,
                    circle.station_ordinal,
                    circle.source_boundary_site,
                )
            ]
            for circle in lane_circles
        )
        lane_items = tuple(canonical[position] for position in lane_positions)
        paper_candidates = tuple(_paper_candidate(circle) for circle in lane_circles)
        selected_positions.append(lane_positions[0])
        # Joining lanes must not let a spacing jump erase a whole source run.
        # The corpus draft anchors each run before engagement refinement.
        stops = []
        represented: set[GuideRunId] = set(lane_items[0].source_run_ids)
        if preserve_source_runs:
            for index, item in enumerate(lane_items[1:], start=1):
                if any(run not in represented for run in item.source_run_ids):
                    stops.append(index)
                    represented.update(item.source_run_ids)
        stops.append(len(paper_candidates) - 1)
        cursor = 0
        for stop in stops:
            if cursor >= stop:
                continue
            span_start = cursor
            span = OrderedPlacementSpan.build(paper_candidates[span_start : stop + 1], tool_radius, cap)
            while cursor < stop:
                placement = select_next_in_span(span, cursor - span_start)
                cursor += placement.candidate_advance
                selected_item = lane_items[cursor]
                selected_positions.append(lane_positions[cursor])
                placements.append(
                    Figure5PlacementRecord(
                        selected_item.hypothesis.projected_boundary_site,
                        selected_item.source_run_ids,
                        placement,
                    )
                )

    selected_items = tuple(canonical[position] for position in selected_positions)
    selected_runs = frozenset(run_id for item in selected_items for run_id in item.source_run_ids)
    if selected_runs != reached_run_ids:
        raise MissingFigure5GuideRunError("Engagement placement did not preserve every reached Figure 5 guide run.")
    selected_hypotheses = tuple(item.hypothesis for item in selected_items)
    selected_representatives = frozenset(item.hypothesis.run_id for item in selected_items)
    path = build_figure5_boundary_path(
        inward_components=inward_components,
        side_candidates=selected_hypotheses,
        reached_run_ids=selected_representatives,
        published_start_center=published_start_center,
        published_start_radius=published_start_radius,
        tool_radius=tool_radius,
        boundary_evidence_bound=boundary_evidence_bound,
        start_evidence_bound=start_evidence_bound,
        offset_sites=None if canonical_sites is None else dict(zip(selected_hypotheses, (canonical_sites[position] for position in selected_positions), strict=True)),
    )
    sources_by_origin = {
        _origin_key(
            item.hypothesis.run_id,
            item.hypothesis.station_ordinal,
            item.hypothesis.projected_boundary_site,
        ): item.source_run_ids
        for item in selected_items
    }
    return HypothesisFigure5Path(
        path,
        tuple(placements),
        tuple(
            sources_by_origin[
                _origin_key(
                    circle.run_id,
                    circle.station_ordinal,
                    circle.source_boundary_site,
                )
            ]
            for circle in path.circles
        ),
        reached_run_ids,
        len(hypotheses),
        len(canonical),
        len(placement_lanes),
    )


def build_figure5_approximate_path(
    guide: Figure5RawGuide,
    inward_components: tuple[tuple[Point2[WorldXY], ...], ...],
) -> HypothesisFigure5Path:
    """Build Figure 5(a) from every pre-thinning raw-guide hypothesis."""
    if type(guide) is not Figure5RawGuide or guide.case.name != "figure5":
        raise Figure5PathIntegrationError("Approximate path construction requires the canonical Figure 5 raw guide.")
    return build_held_reference_path(guide, inward_components, project_contacts=False)


def build_held_reference_path(
    guide: Figure5RawGuide,
    inward_components: tuple[tuple[Point2[WorldXY], ...], ...],
    *,
    project_contacts: bool = True,
    qualify_contacts: bool = False,
) -> HypothesisFigure5Path:
    """Generate the standard polygon draft on a prepared Held reference guide.

    Uses the case's publisher-derived start; does not claim contour awareness.
    """
    if type(guide) is not Figure5RawGuide:
        raise Figure5PathIntegrationError("Held path construction requires a reference raw guide.")
    if qualify_contacts and not project_contacts:
        raise Figure5PathIntegrationError("Contact qualification requires offset contact projection.")
    start = guide.case.start_marker
    start_radius = guide.case.start_marker_radius
    if start is None or start_radius is None:
        raise Figure5PathIntegrationError("Canonical Figure 5 requires its published start-circle evidence.")
    hypotheses = tuple(hypothesis for run in guide.runs for station in run.stations for hypothesis in build_distance_admissible_hypotheses(guide, station))
    offset_sites: tuple[ProjectionBoundarySite, ...] | None = None
    rejected_contacts: tuple[ProjectionAdmissibleBoundaryHypothesis, ...] = ()
    if project_contacts:
        if len(inward_components) != 1:
            raise Figure5PathIntegrationError("Reference draft requires one connected inward boundary.")
        component = inward_components[0]
        if qualify_contacts:
            qualified = qualify_reference_contact_sites(hypotheses, component, guide.site_budget.admissible_distance_gap)
            # Hypotheses are pairwise distinct by construction (one per native
            # site of one (run, station)), so the parallel form drops nothing
            # the hypothesis-keyed form would have merged.
            hypotheses, offset_sites, rejected_contacts = qualified.hypotheses, qualified.sites, qualified.rejected
        else:
            offset_sites = project_bounded_offset_contacts(
                tuple(hypothesis.contact_point for hypothesis in hypotheses),
                component,
                guide.site_budget.admissible_distance_gap,
            )
    reached = frozenset(run.run_id for run in guide.runs)
    start_evidence_bound = Millimetre(float(guide.site_budget.reconstruction_bound) + float(guide.site_budget.projection_bound))
    result = build_hypothesis_figure5_path(
        inward_components=inward_components,
        hypotheses=hypotheses,
        reached_run_ids=reached,
        published_start_center=start,
        published_start_radius=start_radius,
        boundary_evidence_bound=guide.site_budget.admissible_distance_gap,
        start_evidence_bound=start_evidence_bound,
        tool_radius=guide.case.tool_radius,
        cap=EngagementCap.build(math.radians(float(guide.case.tea_cap))),
        offset_sites=offset_sites,
        preserve_source_runs=project_contacts,
    )
    return replace(result, rejected_contact_hypotheses=rejected_contacts)
