"""Projection hypotheses over Figure 5 pre-engagement-thinning guide runs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from functools import cached_property
from typing import NewType

import numpy as np

from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.spec import PocketSpec
from compas_cgal import _toolpath  # type: ignore[attr-defined]
from compas_cgal.engagement_toolpath import GUIDE_SAMPLES_PER_RADIAN
from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import GUIDE_STEPOVER_ADVANCE_MULTIPLE
from compas_cgal.engagement_toolpath import EmptyGuideError
from compas_cgal.held_boundary_clearance import BoundaryClearanceSite
from compas_cgal.held_boundary_clearance import HeldBoundaryClearanceIndex
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.toolpath import _holes_to_vertex_arrays
from compas_cgal.toolpath import _polygon_to_ccw_vertices
from compas_cgal.toolpath import _warn_if_truncated
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import SquaredMillimetre
from compas_cgal.units import WorldXY

ProjectionBoundarySegmentId = NewType("ProjectionBoundarySegmentId", int)
ProjectionBoundaryVertexId = NewType("ProjectionBoundaryVertexId", int)
GuideRunId = NewType("GuideRunId", int)
GuideRunStationOrdinal = NewType("GuideRunStationOrdinal", int)

# The fixed Figure 5 configuration currently emits 126 runs. This bound is only
# a loud truncation guard for the rendering seam, not a machining policy.
FIGURE5_MAX_GUIDE_RUNS = 1000

# The native guide request of ``compas_cgal.engagement_toolpath._guide_chains``,
# spelled out because the stations are read straight from the native rows: the
# arguments ``trochoidal_mat_toolpath_circular`` fills in by default (no radius
# clamp, unscaled clearance, 20 samples per arc, cut plane z = 0, no clearance
# plane) and the ones ``_guide_chains`` sets (no leads, links or final retract,
# nearest-entry order, conventional turn).
_GUIDE_MIN_TROCHOID_RADIUS = 0.0
_GUIDE_MAX_TROCHOID_RADIUS = math.inf
_GUIDE_MAT_SCALE = 1.0
_GUIDE_SAMPLES_PER_CYCLE = 20
_GUIDE_LEAD = 0.0
_GUIDE_CUT_Z = 0.0
_GUIDE_CLIMB = False

# Native primitive-type codes of a meta row (``toolpath._matrices_to_operations``):
# 0 line and 1 arc are bridges, 2 is a machining circle, the station.
_NATIVE_BRIDGE_PRIMITIVES = (0, 1)
_NATIVE_CIRCLE_PRIMITIVE = 2


RationalPoint = tuple[Fraction, Fraction]
# The canonical parameter of a site at a ring vertex, built once.
_ZERO_PARAMETER = Fraction(0)
RationalSegment = tuple[RationalPoint, RationalPoint]


def _rational_point(point: Point2[WorldXY]) -> RationalPoint:
    return Fraction.from_float(float(point.x)), Fraction.from_float(float(point.y))


class InvalidFigure5RawGuideError(ValueError):
    """The requested input cannot define the canonical Figure 5 raw guide."""


class UnknownNativeGuidePrimitiveError(ValueError):
    """A native guide row carries a primitive-type code outside line, arc and circle."""


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
        # Exact integer tests on the normalised rational n/d (d > 0): 0 <= n/d <= 1
        # is 0 <= n <= d, n/d == 1 is n == d, n/d == 0 is n == 0. No Fraction is
        # built or compared (484k sites on Monstera).
        numerator = parameter.numerator
        denominator = parameter.denominator
        if not 0 <= numerator <= denominator:
            raise InvalidFigure5RawGuideError("Projected boundary side parameter must lie in the closed unit interval.")
        canonical_segment = int(segment_id)
        canonical_parameter = parameter
        if numerator == denominator:
            canonical_segment = (canonical_segment + 1) % segment_count
            canonical_parameter = _ZERO_PARAMETER
            numerator = 0
        vertex_id = ProjectionBoundaryVertexId(canonical_segment) if numerator == 0 else None
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
        middle_point: Point2[WorldXY],
    ) -> "Figure5RawGuideStation":
        if run_id < 0 or ordinal_on_run < 0:
            raise InvalidFigure5RawGuideError("Figure 5 guide identities must be non-negative.")
        return cls(
            GuideRunId(run_id),
            GuideRunStationOrdinal(ordinal_on_run),
            middle_point,
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


@dataclass
class Figure5RawGuideDerived:
    """Memo of a guide's pure derived values: the only mutable part of a guide.

    Each slot is a pure function of the owning guide's frozen fields, filled
    on first read. The guide owns exactly one, created in its ``__init__``
    (``init=False``), so ``dataclasses.replace`` never shares it between guides
    and the frozen guide's own attributes are never written after construction.
    """

    rational_boundary_segments: tuple[RationalSegment, ...] | None = None
    boundary_side_lengths: tuple[float, ...] | None = None
    boundary_side_offsets: tuple[float, ...] | None = None
    station_boundary_analyses: tuple[tuple[_NativeStationBoundaryAnalysis, ...], ...] | None = None


@dataclass(frozen=True)
class Figure5RawGuide:
    """Boundary evidence plus every emitted pre-engagement-thinning station."""

    case: HeldReferenceCase
    site_budget: Figure5BoundarySiteBudget
    boundary_length: Millimetre
    runs: tuple[Figure5RawGuideRun, ...]
    derived: Figure5RawGuideDerived = field(default_factory=Figure5RawGuideDerived, init=False, repr=False, compare=False)

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

    @property
    def rational_boundary_segments(self) -> tuple[RationalSegment, ...]:
        """Closed projected ring as exact segments, injected once per guide."""
        if self.derived.rational_boundary_segments is None:
            points = tuple(_rational_point(point) for point in self.case.projection.points)
            self.derived.rational_boundary_segments = tuple(zip(points, (*points[1:], points[0])))
        return self.derived.rational_boundary_segments

    @property
    def boundary_side_lengths(self) -> tuple[float, ...]:
        """Reporting double length of every projected side, in ring order."""
        if self.derived.boundary_side_lengths is None:
            points = self.case.projection.points
            self.derived.boundary_side_lengths = tuple(
                math.hypot(
                    float(end.x) - float(start.x),
                    float(end.y) - float(start.y),
                )
                for start, end in zip(points, (*points[1:], points[0]))
            )
        return self.derived.boundary_side_lengths

    @property
    def boundary_side_offsets(self) -> tuple[float, ...]:
        """Reporting arclength at the start of every side.

        Each entry is ``sum(lengths[:k])`` itself, not a running total: Python's
        float ``sum`` is compensated, so a running prefix would differ in the
        last bits from the per-site sum this cache replaces.
        """
        if self.derived.boundary_side_offsets is None:
            lengths = self.boundary_side_lengths
            self.derived.boundary_side_offsets = tuple(sum(lengths[:ordinal]) for ordinal in range(len(lengths)))
        return self.derived.boundary_side_offsets

    @property
    def station_boundary_analyses(self) -> tuple[tuple[_NativeStationBoundaryAnalysis, ...], ...]:
        """Every station's exact boundary analysis, one native batch per guide, indexed [run][ordinal]."""
        if self.derived.station_boundary_analyses is None:
            self.derived.station_boundary_analyses = _native_station_boundary_analyses(self)
        return self.derived.station_boundary_analyses


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


def _boundary_segments(guide: Figure5RawGuide) -> tuple[RationalSegment, ...]:
    return guide.rational_boundary_segments


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


def _require_owned_station(guide: Figure5RawGuide, station: Figure5RawGuideStation) -> None:
    run_index = int(station.run_id)
    station_index = int(station.ordinal_on_run)
    if run_index < 0 or run_index >= len(guide.runs):
        raise InvalidFigure5RawGuideError("Boundary-side query requires a station owned by its raw guide.")
    run = guide.runs[run_index]
    if station_index < 0 or station_index >= len(run.stations) or run.stations[station_index] != station:
        raise InvalidFigure5RawGuideError("Boundary-side query requires a station owned by its raw guide.")


def _station_records(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]:
    _require_owned_station(guide, station)
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
    """Exact (site, footpoint, squared distance) records plus the admissible sites.

    The Fraction path records every ring segment; the native path records only
    the admissible sites. Hypotheses read records of admissible sites only.
    """

    records: tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]
    distance_admissible_sites: tuple[ProjectionBoundarySite, ...]


@dataclass(frozen=True)
class _NativeStationBoundaryAnalysis:
    """The native batch entry: the admissible sites and, aligned, their clearance sites.

    The clearance sites keep the native integers; ``records`` materialises the
    exact Fraction records of ``_StationBoundaryAnalysis`` only when read (the
    oracle comparisons and the single-hypothesis path), never on the bulk path.
    """

    distance_admissible_sites: tuple[ProjectionBoundarySite, ...]
    clearance_sites: tuple[BoundaryClearanceSite, ...]

    @cached_property
    def records(self) -> tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]:
        return tuple(
            (site, (clearance.footpoint[0], clearance.footpoint[1]), clearance.squared_distance)
            for site, clearance in zip(self.distance_admissible_sites, self.clearance_sites, strict=True)
        )


def _native_station_boundary_analyses(
    guide: Figure5RawGuide,
) -> tuple[tuple[_NativeStationBoundaryAnalysis, ...], ...]:
    """Analyse every station in one spatially indexed exact native batch.

    Decides the same admissible set as ``_analyze_station_boundary``: the ring
    and stations are injected exactly, and the gap is the same double the
    Fraction path injects with ``Fraction.from_float``.
    """
    index = HeldBoundaryClearanceIndex.build(guide.case.projection.points)
    stations = tuple(station.middle_point for run in guide.runs for station in run.stations)
    gap = Millimetre(float(guide.site_budget.admissible_distance_gap))
    rows = iter(index.distance_admissible_sites(stations, gap))
    segment_count = guide.boundary_segment_count
    analyses: list[tuple[_NativeStationBoundaryAnalysis, ...]] = []
    for run in guide.runs:
        run_analyses: list[_NativeStationBoundaryAnalysis] = []
        for _ in run.stations:
            clearance_sites = next(rows)
            sites = tuple(ProjectionBoundarySite.build(ProjectionBoundarySegmentId(site.segment), site.parameter, segment_count) for site in clearance_sites)
            run_analyses.append(_NativeStationBoundaryAnalysis(sites, clearance_sites))
        analyses.append(tuple(run_analyses))
    return tuple(analyses)


def _native_station_boundary(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> _NativeStationBoundaryAnalysis:
    """The generator's station analysis: the guide's cached native batch entry."""
    _require_owned_station(guide, station)
    return guide.station_boundary_analyses[int(station.run_id)][int(station.ordinal_on_run)]


def _analyze_station_boundary(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> _StationBoundaryAnalysis:
    """Brute-force Fraction analysis over every ring segment; the equivalence oracle."""
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
    return _native_station_boundary(guide, station).distance_admissible_sites


def _boundary_arclength(
    guide: Figure5RawGuide,
    site: ProjectionBoundarySite,
) -> Millimetre:
    ordinal = int(site.segment_id)
    # ``n / d`` in ints is ``float(Fraction)`` (one correctly rounded division),
    # without the ``numbers.Rational.__float__`` dispatch.
    parameter = site.parameter
    arclength = guide.boundary_side_offsets[ordinal] + parameter.numerator / parameter.denominator * guide.boundary_side_lengths[ordinal]
    return Millimetre(arclength % float(guide.boundary_length))


def _exact_square(value: float) -> Fraction:
    """Exact square of a double, injected once with ``Fraction.from_float``."""
    exact = Fraction.from_float(value)
    return exact * exact


def _hypothesis_from_record(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
    admissible: tuple[ProjectionBoundarySite, ...],
    footpoint_exact: RationalPoint,
    squared_clearance: Fraction,
    tool_radius: float,
    squared_tool_radius: Fraction,
) -> ProjectionAdmissibleBoundaryHypothesis:
    """Paper geometry from one admissible Fraction record.

    Decides ``clearance**2 > r**2`` exactly (``squared_tool_radius`` is the
    exact square of the double ``tool_radius``); everything after it reports.
    """
    if squared_clearance <= squared_tool_radius:
        raise InvalidProjectionAdmissibleHypothesisError("Figure 5 projected-boundary hypothesis has no positive paper-derived radius.")
    return _hypothesis_geometry(guide, station, site, admissible, float(footpoint_exact[0]), float(footpoint_exact[1]), float(squared_clearance), tool_radius)


def _hypothesis_from_clearance(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
    admissible: tuple[ProjectionBoundarySite, ...],
    clearance_site: BoundaryClearanceSite,
    tool_radius: float,
    squared_tool_radius: SquaredMillimetre,
) -> ProjectionAdmissibleBoundaryHypothesis:
    """``_hypothesis_from_record`` on the native integers: no Fraction is built.

    The verdict ``clearance**2 > r**2`` is decided exactly in integers; the
    doubles are the correctly rounded quotients, equal to ``float`` of the
    record's Fractions.
    """
    if not clearance_site.squared_distance_exceeds(squared_tool_radius):
        raise InvalidProjectionAdmissibleHypothesisError("Figure 5 projected-boundary hypothesis has no positive paper-derived radius.")
    footpoint_x, footpoint_y = clearance_site.footpoint_doubles
    return _hypothesis_geometry(guide, station, site, admissible, footpoint_x, footpoint_y, clearance_site.squared_distance_double, tool_radius)


def _hypothesis_geometry(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
    admissible: tuple[ProjectionBoundarySite, ...],
    footpoint_x: float,
    footpoint_y: float,
    squared_clearance: float,
    tool_radius: float,
) -> ProjectionAdmissibleBoundaryHypothesis:
    """The reported paper geometry of one admitted record; decides nothing."""
    middle = station.middle_point
    # ``Point2.build``, not the ``Point2[WorldXY]`` alias: the alias forwards every
    # attribute through typing (hot path, 95k calls on Figure 5); the frame is in
    # the annotation and the built value is the same.
    footpoint: Point2[WorldXY] = Point2.build(footpoint_x, footpoint_y)
    clearance = math.sqrt(squared_clearance)
    # The local doubles are the values the points hold (``Point2.build`` keeps a
    # finite float as is), so the arithmetic reads them, not the attributes.
    direction_x = (float(middle.x) - footpoint_x) / clearance
    direction_y = (float(middle.y) - footpoint_y) / clearance
    radius = (clearance - tool_radius) / 2.0
    contact_x = footpoint_x + tool_radius * direction_x
    contact_y = footpoint_y + tool_radius * direction_y
    contact: Point2[WorldXY] = Point2.build(contact_x, contact_y)
    center: Point2[WorldXY] = Point2.build(
        contact_x + radius * direction_x,
        contact_y + radius * direction_y,
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


def _build_projection_admissible_hypothesis(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
    analysis: _StationBoundaryAnalysis | _NativeStationBoundaryAnalysis,
) -> ProjectionAdmissibleBoundaryHypothesis:
    admissible = analysis.distance_admissible_sites
    if site not in admissible:
        raise InadmissibleFigure5BoundarySideError(f"Projected boundary site {int(site.segment_id)}:{site.parameter} is outside the Figure 5 evidence budget.")
    footpoint_exact, squared_clearance = next((footpoint, distance) for record_site, footpoint, distance in analysis.records if record_site == site)
    tool_radius = float(guide.case.tool_radius.value)
    return _hypothesis_from_record(guide, station, site, admissible, footpoint_exact, squared_clearance, tool_radius, _exact_square(tool_radius))


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
        _native_station_boundary(guide, station),
    )


def build_distance_admissible_hypotheses(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    """Build every distance-admissible hypothesis without claiming ownership.

    The native analysis holds exactly the admissible sites, aligned with
    their clearance sites, so each is used directly: no membership test, no
    record search, and no Fraction (the clearance verdict is decided in
    integers).
    """
    tool_radius = float(guide.case.tool_radius.value)
    squared_tool_radius = SquaredMillimetre(_exact_square(tool_radius))
    return _hypotheses_from_analysis(guide, station, _native_station_boundary(guide, station), tool_radius, squared_tool_radius)


def build_guide_distance_admissible_hypotheses(guide: Figure5RawGuide) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    """Every station's distance-admissible hypotheses, in run and station order, from one native batch.

    Equal to concatenating :func:`build_distance_admissible_hypotheses` over
    ``guide.runs`` and their stations. The guide's memoised station analyses
    are read when already filled; otherwise the batch runs once and is not
    memoised, so its analyses (two objects per admissible site) are released
    when this returns instead of living as long as the guide.
    """
    analyses = guide.derived.station_boundary_analyses
    if analyses is None:
        analyses = _native_station_boundary_analyses(guide)
    tool_radius = float(guide.case.tool_radius.value)
    squared_tool_radius = SquaredMillimetre(_exact_square(tool_radius))
    return tuple(
        hypothesis
        for run, run_analyses in zip(guide.runs, analyses, strict=True)
        for station, analysis in zip(run.stations, run_analyses, strict=True)
        for hypothesis in _hypotheses_from_analysis(guide, station, analysis, tool_radius, squared_tool_radius)
    )


def _hypotheses_from_analysis(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    analysis: _NativeStationBoundaryAnalysis,
    tool_radius: float,
    squared_tool_radius: SquaredMillimetre,
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    """One station's hypotheses from its native analysis, in admissible-site order."""
    admissible = analysis.distance_admissible_sites
    return tuple(
        _hypothesis_from_clearance(guide, station, site, admissible, clearance_site, tool_radius, squared_tool_radius)
        for site, clearance_site in zip(admissible, analysis.clearance_sites, strict=True)
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
    runs = tuple(
        Figure5RawGuideRun.build(
            run_index,
            tuple(Figure5RawGuideStation.build(run_index, ordinal, Point2.build(cx, cy)) for ordinal, (cx, cy) in enumerate(centres)),
        )
        for run_index, centres in enumerate(_guide_station_centres(case))
    )
    return Figure5RawGuide.build_reference(case, runs)


def _guide_station_centres(case: HeldReferenceCase) -> list[list[tuple[float, float]]]:
    """Every station centre of ``_guide_chains``' guide, read straight from the native rows.

    The same native request as ``compas_cgal.engagement_toolpath._guide_chains``;
    the circle rows, grouped by runs of equal ``path_index`` in emission order,
    are its chains. The raw guide reads only the centre doubles, so no compas
    operation is built for the rows (they were discarded after two doubles each).

    Raises:
        UnknownNativeGuidePrimitiveError: A row's primitive type is not a line, arc or circle.
        EmptyGuideError: The guide emitted no machining circle.
    """
    spec = case.guide_pocket_spec()
    meta, centres = _native_guide_rows(spec)
    runs: list[list[tuple[float, float]]] = []
    current_index: int | None = None
    meta_rows: list[list[float]] = np.asarray(meta, dtype=np.float64).tolist()
    centre_rows: list[list[float]] = np.asarray(centres, dtype=np.float64).tolist()
    for meta_row, centre in zip(meta_rows, centre_rows, strict=True):
        primitive = int(meta_row[1])
        if primitive != _NATIVE_CIRCLE_PRIMITIVE:
            if primitive not in _NATIVE_BRIDGE_PRIMITIVES:
                raise UnknownNativeGuidePrimitiveError(f"native guide row has primitive type {primitive}; expected 0 (line), 1 (arc) or 2 (circle).")
            continue
        path_index = int(meta_row[0])
        if path_index != current_index:
            runs.append([])
            current_index = path_index
        runs[-1].append((centre[0], centre[1]))
    if not runs:
        raise EmptyGuideError(
            f"The straight-skeleton guide produced no machining circles for tool_diameter={spec.tool_diameter!r}; the pocket admits no gouge-free trochoid at this tool size."
        )
    return runs


def _native_guide_rows(spec: PocketSpec) -> tuple[np.ndarray, np.ndarray]:
    """The native guide request of ``_guide_station_centres``: its meta and centre rows.

    ``trochoidal_mat_toolpath_circular`` decides ring simplicity exactly
    (``data_to_polygon``: CGAL ``Polygon_2::is_simple`` on the injected doubles)
    and raises ``ValueError`` otherwise; this is the simplicity check of the
    guide path.
    """
    guide_step = GUIDE_STEP_TOOL_DIAMETERS * spec.tool_diameter
    meta, _starts, _ends, centres, _radii, _polyline, _start_tangents, _end_tangents, skipped = _toolpath.trochoidal_mat_toolpath_circular(
        _polygon_to_ccw_vertices(spec.polygon),
        _holes_to_vertex_arrays(list(spec.holes)),
        float(spec.tool_diameter),
        float(GUIDE_STEPOVER_ADVANCE_MULTIPLE * guide_step),
        float(guide_step),
        _GUIDE_MIN_TROCHOID_RADIUS,
        _GUIDE_MAX_TROCHOID_RADIUS,
        _GUIDE_MAT_SCALE,
        float(RADIAL_CLEARANCE_FRACTION * spec.tool_diameter),
        _GUIDE_SAMPLES_PER_CYCLE,
        FIGURE5_MAX_GUIDE_RUNS,
        _GUIDE_CLIMB,
        _GUIDE_LEAD,
        _GUIDE_LEAD,
        False,
        True,
        _GUIDE_CUT_Z,
        0.0,
        False,
        False,
        float(GUIDE_SAMPLES_PER_RADIAN),
    )
    _warn_if_truncated(int(skipped), FIGURE5_MAX_GUIDE_RUNS)
    return meta, centres
