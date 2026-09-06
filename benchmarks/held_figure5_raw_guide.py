"""Projection hypotheses over Figure 5 pre-engagement-thinning guide runs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import NewType

from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import _guide_chains
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION

ProjectionBoundarySegmentId = NewType("ProjectionBoundarySegmentId", int)
ProjectionBoundaryVertexId = NewType("ProjectionBoundaryVertexId", int)
GuideRunId = NewType("GuideRunId", int)
GuideRunStationOrdinal = NewType("GuideRunStationOrdinal", int)

# The fixed Figure 5 configuration currently emits 126 runs. This bound is only
# a loud truncation guard for the rendering seam, not a machining policy.
FIGURE5_MAX_GUIDE_RUNS = 1000


class InvalidFigure5RawGuideError(ValueError):
    """The requested input cannot define the canonical Figure 5 raw guide."""


class InadmissibleFigure5BoundarySideError(ValueError):
    """A requested projected site is outside the distance-admissible set."""


class InvalidProjectionAdmissibleHypothesisError(ValueError):
    """A distance-admissible hypothesis cannot support positive geometry."""


@dataclass(frozen=True)
class ProjectionBoundarySite:
    """Canonical point identity on the closed polygonal projection."""

    segment_id: ProjectionBoundarySegmentId
    parameter: Fraction
    vertex_id: ProjectionBoundaryVertexId | None

    @classmethod
    def build(
        cls,
        segment_id: ProjectionBoundarySegmentId,
        parameter: Fraction,
        segment_count: int,
    ) -> "ProjectionBoundarySite":
        if segment_count <= 0 or not 0 <= int(segment_id) < segment_count:
            raise InvalidFigure5RawGuideError("Projected boundary site has no segment in the closed ring.")
        if not Fraction(0) <= parameter <= Fraction(1):
            raise InvalidFigure5RawGuideError("Projected boundary side parameter must lie in the closed unit interval.")
        canonical_segment = int(segment_id)
        canonical_parameter = parameter
        if parameter == 1:
            canonical_segment = (canonical_segment + 1) % segment_count
            canonical_parameter = Fraction(0)
        vertex_id = ProjectionBoundaryVertexId(canonical_segment) if canonical_parameter == 0 else None
        return cls(
            ProjectionBoundarySegmentId(canonical_segment),
            canonical_parameter,
            vertex_id,
        )


@dataclass(frozen=True)
class Figure5BoundarySiteBudget:
    """Distance-gap budget induced by both recorded boundary approximations."""

    reconstruction_bound: Millimetre
    projection_bound: Millimetre
    admissible_distance_gap: Millimetre

    @classmethod
    def build(cls, case: HeldReferenceCase) -> "Figure5BoundarySiteBudget":
        reconstruction = float(case.reconstruction.deviation_upper_bound)
        projection = float(case.projection.deviation_limit)
        # Each of two candidate distances may move by the full combined boundary
        # displacement, hence the factor two on their admissible difference.
        gap = 2.0 * (reconstruction + projection)
        return cls(
            Millimetre(reconstruction),
            Millimetre(projection),
            Millimetre(gap),
        )


@dataclass(frozen=True)
class Figure5RawGuideStation:
    """One pre-engagement-thinning emitted guide-run station in world XY."""

    run_id: GuideRunId
    ordinal_on_run: GuideRunStationOrdinal
    middle_point: Point2[WorldXY]

    @classmethod
    def build(
        cls,
        run_id: int,
        ordinal_on_run: int,
        station: _GuideStation,
    ) -> "Figure5RawGuideStation":
        if run_id < 0 or ordinal_on_run < 0:
            raise InvalidFigure5RawGuideError("Figure 5 guide identities must be non-negative.")
        return cls(
            GuideRunId(run_id),
            GuideRunStationOrdinal(ordinal_on_run),
            Point2[WorldXY].build(station.cx, station.cy),
        )


@dataclass(frozen=True)
class Figure5RawGuideRun:
    """One ordered run emitted by the current baseline guide generator."""

    run_id: GuideRunId
    stations: tuple[Figure5RawGuideStation, ...]

    @classmethod
    def build(
        cls,
        run_index: int,
        stations: tuple[Figure5RawGuideStation, ...],
    ) -> "Figure5RawGuideRun":
        if not stations or any(station.run_id != run_index or station.ordinal_on_run != ordinal for ordinal, station in enumerate(stations)):
            raise InvalidFigure5RawGuideError("Figure 5 raw guide run must be non-empty and ordered.")
        return cls(GuideRunId(run_index), stations)


@dataclass(frozen=True)
class Figure5RawGuide:
    """Boundary evidence plus every emitted pre-engagement-thinning station."""

    case: HeldReferenceCase
    site_budget: Figure5BoundarySiteBudget
    boundary_length: Millimetre
    runs: tuple[Figure5RawGuideRun, ...]

    @classmethod
    def build(
        cls,
        case: HeldReferenceCase,
        runs: tuple[Figure5RawGuideRun, ...],
    ) -> "Figure5RawGuide":
        if type(case) is not HeldReferenceCase or case.name != "figure5":
            raise InvalidFigure5RawGuideError("Raw guide extraction is restricted to canonical Figure 5.")
        return cls.build_reference(case, runs)

    @classmethod
    def build_reference(
        cls,
        case: HeldReferenceCase,
        runs: tuple[Figure5RawGuideRun, ...],
    ) -> "Figure5RawGuide":
        """Build a guide over any of the four prepared Held reference pockets."""
        if type(case) is not HeldReferenceCase or case.name not in CANONICAL_CASE_NAMES:
            raise InvalidFigure5RawGuideError("Raw guide requires a canonical Held reference case.")
        if not runs:
            raise InvalidFigure5RawGuideError("Held reference raw guide is empty.")
        points = case.projection.points
        boundary_length = sum(
            math.hypot(
                float(end.x) - float(start.x),
                float(end.y) - float(start.y),
            )
            for start, end in zip(points, (*points[1:], points[0]))
        )
        return cls(
            case,
            Figure5BoundarySiteBudget.build(case),
            Millimetre(boundary_length),
            runs,
        )

    @property
    def station_count(self) -> int:
        return sum(len(run.stations) for run in self.runs)

    @property
    def boundary_segment_count(self) -> int:
        return len(self.case.projection.points)


@dataclass(frozen=True)
class ProjectionAdmissibleBoundaryHypothesis:
    """Paper geometry for one distance-admissible projected-boundary hypothesis."""

    run_id: GuideRunId
    station_ordinal: GuideRunStationOrdinal
    projected_boundary_site: ProjectionBoundarySite
    distance_admissible_sites: tuple[ProjectionBoundarySite, ...]
    boundary_arclength: Millimetre
    middle_point: Point2[WorldXY]
    boundary_footpoint: Point2[WorldXY]
    contact_point: Point2[WorldXY]
    center: Point2[WorldXY]
    guide_radius: GuideRadius


RationalPoint = tuple[Fraction, Fraction]
RationalSegment = tuple[RationalPoint, RationalPoint]


def _rational_point(point: Point2[WorldXY]) -> RationalPoint:
    return Fraction.from_float(float(point.x)), Fraction.from_float(float(point.y))


def _boundary_segments(guide: Figure5RawGuide) -> tuple[RationalSegment, ...]:
    points = tuple(_rational_point(point) for point in guide.case.projection.points)
    return tuple(zip(points, (*points[1:], points[0])))


def _project_exact(
    point: RationalPoint,
    segment: RationalSegment,
) -> tuple[RationalPoint, Fraction]:
    start, end = segment
    edge_x = end[0] - start[0]
    edge_y = end[1] - start[1]
    squared_length = edge_x * edge_x + edge_y * edge_y
    if squared_length == 0:
        raise InvalidFigure5RawGuideError("Figure 5 projected boundary contains a zero-length segment.")
    numerator = (point[0] - start[0]) * edge_x + (point[1] - start[1]) * edge_y
    parameter = min(Fraction(1), max(Fraction(0), numerator / squared_length))
    return (
        (start[0] + parameter * edge_x, start[1] + parameter * edge_y),
        parameter,
    )


def _squared_distance(first: RationalPoint, second: RationalPoint) -> Fraction:
    delta_x = second[0] - first[0]
    delta_y = second[1] - first[1]
    return delta_x * delta_x + delta_y * delta_y


def _within_distance_gap(
    candidate_squared: Fraction,
    minimum_squared: Fraction,
    gap: Fraction,
) -> bool:
    """Decide ``sqrt(candidate)-sqrt(minimum) <= gap`` exactly."""
    if candidate_squared <= minimum_squared:
        return True
    gap_squared = gap * gap
    remainder = candidate_squared - minimum_squared - gap_squared
    if remainder <= 0:
        return True
    return remainder * remainder <= 4 * minimum_squared * gap_squared


def _station_records(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]:
    run_index = int(station.run_id)
    station_index = int(station.ordinal_on_run)
    if run_index < 0 or run_index >= len(guide.runs):
        raise InvalidFigure5RawGuideError("Boundary-side query requires a station owned by its raw guide.")
    run = guide.runs[run_index]
    if station_index < 0 or station_index >= len(run.stations) or run.stations[station_index] != station:
        raise InvalidFigure5RawGuideError("Boundary-side query requires a station owned by its raw guide.")
    middle = _rational_point(station.middle_point)
    records: dict[ProjectionBoundarySite, tuple[RationalPoint, Fraction]] = {}
    for ordinal, segment in enumerate(_boundary_segments(guide)):
        footpoint, parameter = _project_exact(middle, segment)
        site = ProjectionBoundarySite.build(
            ProjectionBoundarySegmentId(ordinal),
            parameter,
            guide.boundary_segment_count,
        )
        records[site] = (footpoint, _squared_distance(middle, footpoint))
    return tuple(
        (site, footpoint, squared_distance)
        for site, (footpoint, squared_distance) in sorted(
            records.items(),
            key=lambda item: (int(item[0].segment_id), item[0].parameter),
        )
    )


@dataclass(frozen=True)
class _StationBoundaryAnalysis:
    records: tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]
    distance_admissible_sites: tuple[ProjectionBoundarySite, ...]


def _analyze_station_boundary(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> _StationBoundaryAnalysis:
    records = _station_records(guide, station)
    minimum_squared = min(record[2] for record in records)
    gap = Fraction.from_float(float(guide.site_budget.admissible_distance_gap))
    admissible = tuple(record[0] for record in records if _within_distance_gap(record[2], minimum_squared, gap))
    return _StationBoundaryAnalysis(records, admissible)


def distance_admissible_boundary_sites(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> tuple[ProjectionBoundarySite, ...]:
    """Report every projected site inside the recorded distance-gap budget."""
    return _analyze_station_boundary(guide, station).distance_admissible_sites


def _boundary_arclength(
    guide: Figure5RawGuide,
    site: ProjectionBoundarySite,
) -> Millimetre:
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


def _build_projection_admissible_hypothesis(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
    analysis: _StationBoundaryAnalysis,
) -> ProjectionAdmissibleBoundaryHypothesis:
    admissible = analysis.distance_admissible_sites
    if site not in admissible:
        raise InadmissibleFigure5BoundarySideError(f"Projected boundary site {int(site.segment_id)}:{site.parameter} is outside the Figure 5 evidence budget.")
    footpoint_exact, squared_clearance = next((footpoint, distance) for record_site, footpoint, distance in analysis.records if record_site == site)

    tool_radius = float(guide.case.tool_radius.value)
    exact_tool_radius = Fraction.from_float(tool_radius)
    if squared_clearance <= exact_tool_radius * exact_tool_radius:
        raise InvalidProjectionAdmissibleHypothesisError("Figure 5 projected-boundary hypothesis has no positive paper-derived radius.")
    middle = station.middle_point
    footpoint = Point2[WorldXY].build(
        float(footpoint_exact[0]),
        float(footpoint_exact[1]),
    )
    clearance = math.sqrt(float(squared_clearance))
    direction_x = (float(middle.x) - float(footpoint.x)) / clearance
    direction_y = (float(middle.y) - float(footpoint.y)) / clearance
    radius = (clearance - tool_radius) / 2.0
    contact = Point2[WorldXY].build(
        float(footpoint.x) + tool_radius * direction_x,
        float(footpoint.y) + tool_radius * direction_y,
    )
    center = Point2[WorldXY].build(
        float(contact.x) + radius * direction_x,
        float(contact.y) + radius * direction_y,
    )
    return ProjectionAdmissibleBoundaryHypothesis(
        station.run_id,
        station.ordinal_on_run,
        site,
        admissible,
        _boundary_arclength(guide, site),
        middle,
        footpoint,
        contact,
        center,
        GuideRadius.build(radius),
    )


def build_projection_admissible_hypothesis(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
) -> ProjectionAdmissibleBoundaryHypothesis:
    """Construct paper geometry for one explicit distance-admissible hypothesis."""
    return _build_projection_admissible_hypothesis(
        guide,
        station,
        site,
        _analyze_station_boundary(guide, station),
    )


def build_distance_admissible_hypotheses(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    """Build every distance-admissible hypothesis without claiming ownership."""
    analysis = _analyze_station_boundary(guide, station)
    return tuple(
        _build_projection_admissible_hypothesis(
            guide,
            station,
            site,
            analysis,
        )
        for site in analysis.distance_admissible_sites
    )


def build_figure5_raw_guide(case: HeldReferenceCase) -> Figure5RawGuide:
    """Extract every emitted station before engagement-controlled thinning."""
    if type(case) is not HeldReferenceCase or case.name != "figure5":
        raise InvalidFigure5RawGuideError("Raw guide extraction is restricted to canonical Figure 5.")
    return build_held_reference_raw_guide(case)


def build_held_reference_raw_guide(case: HeldReferenceCase) -> Figure5RawGuide:
    """Extract unthinned stations for a prepared Figure 5 or Figure 8 pocket."""
    if type(case) is not HeldReferenceCase or case.name not in CANONICAL_CASE_NAMES:
        raise InvalidFigure5RawGuideError("Raw guide requires a canonical Held reference case.")
    spec = case.pocket_spec()
    source_runs = _guide_chains(
        spec.polygon,
        spec.tool_diameter,
        GUIDE_STEP_TOOL_DIAMETERS * spec.tool_diameter,
        RADIAL_CLEARANCE_FRACTION * spec.tool_diameter,
        False,
        FIGURE5_MAX_GUIDE_RUNS,
        list(spec.holes),
    )
    runs = tuple(
        Figure5RawGuideRun.build(
            run_index,
            tuple(Figure5RawGuideStation.build(run_index, ordinal, station) for ordinal, station in enumerate(source_run)),
        )
        for run_index, source_run in enumerate(source_runs)
    )
    return Figure5RawGuide.build_reference(case, runs)
