"""Hypothesis-based Figure 5(a) placement and boundary-path integration."""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import MissingFigure5GuideRunError
from benchmarks.held_figure5_boundary_path import build_figure5_boundary_path
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacement
from benchmarks.held_standard_placement import predecessor_overlap_margin
from benchmarks.held_standard_placement import select_next_standard_candidate
from compas_cgal import _circle_geometry_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


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


@dataclass(frozen=True)
class _CanonicalHypothesis:
    hypothesis: ProjectionAdmissibleBoundaryHypothesis
    source_run_ids: tuple[GuideRunId, ...]


@dataclass(frozen=True)
class _PlacementLane:
    """One monotone source-side sequence safe for predecessor placement."""

    source_segment: int
    direction: int
    circles: tuple[Figure5CounterclockwiseCircle, ...]


def _site_key(site: ProjectionBoundarySite) -> tuple[int, Fraction]:
    return int(site.segment_id), site.parameter


def _circle_geometry_key(
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
    return (
        Fraction.from_float(float(hypothesis.center.x)),
        Fraction.from_float(float(hypothesis.center.y)),
        Fraction.from_float(float(hypothesis.guide_radius.value)),
        Fraction.from_float(float(hypothesis.contact_point.x)),
        Fraction.from_float(float(hypothesis.contact_point.y)),
    )


def _hypothesis_order_key(
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
) -> tuple[tuple[Fraction, Fraction, Fraction, Fraction, Fraction], int, Fraction, int, int]:
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


def _boundary_neighborhood_owners(
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


def _canonicalize(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
) -> tuple[_CanonicalHypothesis, ...]:
    grouped: dict[
        tuple[Fraction, Fraction, Fraction, Fraction, Fraction],
        list[ProjectionAdmissibleBoundaryHypothesis],
    ] = {}
    for hypothesis in hypotheses:
        grouped.setdefault(_circle_geometry_key(hypothesis), []).append(hypothesis)
    return tuple(
        _CanonicalHypothesis(
            min(group, key=_hypothesis_order_key),
            tuple(sorted({candidate.run_id for candidate in group}, key=int)),
        )
        for _, group in sorted(grouped.items())
    )


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
    tail: Figure5CounterclockwiseCircle,
    head: Figure5CounterclockwiseCircle,
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
    circles: tuple[Figure5CounterclockwiseCircle, ...],
) -> _PlacementLane:
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
) -> tuple[_PlacementLane, ...]:
    by_run_side: dict[
        tuple[GuideRunId, int],
        list[Figure5CounterclockwiseCircle],
    ] = {}
    for circle in circles:
        key = circle.run_id, int(circle.source_boundary_site.segment_id)
        by_run_side.setdefault(key, []).append(circle)

    emitted_fragments = tuple(
        _run_side_fragment(tuple(sorted(by_run_side[key], key=_circle_station_key))) for key in sorted(by_run_side, key=lambda item: (int(item[0]), item[1]))
    )
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
    lanes: list[_PlacementLane] = []
    for fragment in fragments:
        if fragment.direction == 0:
            lanes.append(fragment)
            continue
        compatible = []
        head = fragment.circles[0]
        for index, lane in enumerate(lanes):
            tail = lane.circles[-1]
            advances = _parameter_direction(
                tail.source_boundary_site.parameter,
                head.source_boundary_site.parameter,
            )
            if (
                lane.source_segment == fragment.source_segment
                and lane.direction == fragment.direction
                and tail.run_id != head.run_id
                and advances == fragment.direction
                and float(
                    predecessor_overlap_margin(
                        _paper_candidate(tail),
                        _paper_candidate(head),
                        tool_radius,
                    )
                )
                > 0.0
            ):
                compatible.append(index)
        if not compatible:
            lanes.append(fragment)
            continue
        scores = {index: _continuation_score(lanes[index].circles[-1], head) for index in compatible}
        best_score = min(scores.values())
        nearest = [index for index, score in scores.items() if score == best_score]
        if len(nearest) != 1:
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
    offset_sites: dict[ProjectionAdmissibleBoundaryHypothesis, ProjectionBoundarySite] | None = None,
    preserve_source_runs: bool = False,
) -> HypothesisFigure5Path:
    """Place distinct admissible hypotheses and emit one CCW boundary path."""
    if not hypotheses or not reached_run_ids:
        raise Figure5PathIntegrationError("Figure 5 integration requires hypotheses and a reached-run ledger.")
    if len(inward_components) != 1 or not inward_components[0]:
        raise Figure5PathIntegrationError("Figure 5 neighborhood ownership requires one non-empty boundary cycle.")
    owned_hypotheses = _boundary_neighborhood_owners(
        hypotheses,
        len(inward_components[0]),
        published_start_center,
    )
    canonical = _canonicalize(owned_hypotheses)
    by_origin = {
        _origin_key(
            item.hypothesis.run_id,
            item.hypothesis.station_ordinal,
            item.hypothesis.projected_boundary_site,
        ): item
        for item in canonical
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
        offset_sites=offset_sites,
    )
    selected_items_list = []
    placements = []
    placement_lanes = _continuous_placement_lanes(ordered_path.circles, tool_radius)
    for lane in placement_lanes:
        lane_circles = lane.circles
        lane_items = tuple(
            by_origin[
                _origin_key(
                    circle.run_id,
                    circle.station_ordinal,
                    circle.source_boundary_site,
                )
            ]
            for circle in lane_circles
        )
        paper_candidates = tuple(_paper_candidate(circle) for circle in lane_circles)
        selected_items_list.append(lane_items[0])
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
            while cursor < stop:
                placement = select_next_standard_candidate(
                    paper_candidates[cursor],
                    paper_candidates[cursor + 1 : stop + 1],
                    tool_radius,
                    cap,
                )
                cursor += placement.candidate_advance
                selected_item = lane_items[cursor]
                selected_items_list.append(selected_item)
                placements.append(
                    Figure5PlacementRecord(
                        selected_item.hypothesis.projected_boundary_site,
                        selected_item.source_run_ids,
                        placement,
                    )
                )

    selected_items = tuple(selected_items_list)
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
        offset_sites=offset_sites,
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
) -> HypothesisFigure5Path:
    """Generate the standard polygon draft on a prepared Held reference guide.

    Uses the case's publisher-derived start; does not claim contour awareness.
    """
    if type(guide) is not Figure5RawGuide:
        raise Figure5PathIntegrationError("Held path construction requires a reference raw guide.")
    start = guide.case.start_marker
    start_radius = guide.case.start_marker_radius
    if start is None or start_radius is None:
        raise Figure5PathIntegrationError("Canonical Figure 5 requires its published start-circle evidence.")
    hypotheses = tuple(hypothesis for run in guide.runs for station in run.stations for hypothesis in build_distance_admissible_hypotheses(guide, station))
    offset_sites = None
    if project_contacts:
        if len(inward_components) != 1:
            raise Figure5PathIntegrationError("Reference draft requires one connected inward boundary.")
        component = inward_components[0]
        boundary = [(float(p.x), float(p.y)) for p in component]
        offset_sites = {}
        for hypothesis in hypotheses:
            side, parameter = _circle_geometry_2.project_boundary_contact(
                boundary,
                (float(hypothesis.contact_point.x), float(hypothesis.contact_point.y)),
                float(guide.site_budget.admissible_distance_gap),
            )
            # Existing site representation; all projection decisions belong to CGAL.
            offset_sites[hypothesis] = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), Fraction(parameter), len(component))
    reached = frozenset(run.run_id for run in guide.runs)
    start_evidence_bound = Millimetre(float(guide.site_budget.reconstruction_bound) + float(guide.site_budget.projection_bound))
    return build_hypothesis_figure5_path(
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
