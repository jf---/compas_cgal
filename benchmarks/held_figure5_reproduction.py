"""Figure 5 comparison geometry on explicitly approximate guide stations.

This module is a rendering adapter, not a generator or a geometric authority.
It reads the circular stations emitted by the existing straight-skeleton
generator and applies Held and Pfeiffer's paper-derived ``q``/``c``/``rho``
construction to them. The algebra is preserved; the station locus is explicitly
approximate and is not a segment-site medial axis.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from compas.geometry import Circle

from benchmarks.errors import InvalidFigure5ApproximateStationError
from benchmarks.errors import InvalidFigure5ReproductionInputError
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.runner import CLEARANCE_Z_TOOL_DIAMETERS
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement_toolpath import engagement_controlled_toolpath
from compas_cgal.toolpath import ToolpathResult

APPROXIMATION_PROVENANCE: Literal["straight-skeleton station approximation"] = "straight-skeleton station approximation"


@dataclass(frozen=True)
class Figure5ApproximateStation:
    """Paper-derived geometry conditional on one approximate guide station."""

    path_index: int
    middle_point: Point2[WorldXY]
    boundary_footpoint: Point2[WorldXY]
    contact_point: Point2[WorldXY]
    center: Point2[WorldXY]
    guide_radius: Millimetre


@dataclass(frozen=True)
class Figure5StationApproximation:
    """Ordered Figure 5 rendering stations with truthful provenance."""

    provenance: Literal["straight-skeleton station approximation"]
    source_path: ToolpathResult
    stations: tuple[Figure5ApproximateStation, ...]

    def __post_init__(self) -> None:
        if self.provenance != APPROXIMATION_PROVENANCE or type(self.source_path) is not ToolpathResult or not self.stations:
            raise InvalidFigure5ReproductionInputError("Figure 5 approximation requires labeled, non-empty stations.")

    @classmethod
    def build(
        cls,
        source_path: ToolpathResult,
        stations: tuple[Figure5ApproximateStation, ...],
    ) -> "Figure5StationApproximation":
        """Validate and label one ordered approximation."""
        return cls(APPROXIMATION_PROVENANCE, source_path, stations)


RationalPoint = tuple[Fraction, Fraction]
RationalSegment = tuple[RationalPoint, RationalPoint]


def _rational_point(x: float, y: float) -> RationalPoint:
    return Fraction.from_float(x), Fraction.from_float(y)


def _boundary_segments(
    case: HeldReferenceCase,
) -> tuple[RationalSegment, ...]:
    boundary = tuple(_rational_point(float(point.x), float(point.y)) for point in case.projection.points)
    return tuple(zip(boundary, (*boundary[1:], boundary[0])))


def _closest_point_exact(point: RationalPoint, segment: RationalSegment) -> RationalPoint:
    start, end = segment
    edge_x = end[0] - start[0]
    edge_y = end[1] - start[1]
    squared_length = edge_x * edge_x + edge_y * edge_y
    if squared_length == 0:
        raise InvalidFigure5ReproductionInputError("Figure 5 projected boundary contains a zero-length segment.")
    parameter_numerator = (point[0] - start[0]) * edge_x + (point[1] - start[1]) * edge_y
    if parameter_numerator <= 0:
        return start
    if parameter_numerator >= squared_length:
        return end
    parameter = parameter_numerator / squared_length
    return start[0] + parameter * edge_x, start[1] + parameter * edge_y


def _squared_distance(first: RationalPoint, second: RationalPoint) -> Fraction:
    delta_x = second[0] - first[0]
    delta_y = second[1] - first[1]
    return delta_x * delta_x + delta_y * delta_y


def _boundary_footpoint(
    *,
    middle_point: RationalPoint,
    phase_point: RationalPoint,
    boundary_segments: tuple[RationalSegment, ...],
) -> tuple[RationalPoint, Fraction]:
    candidates = []
    for ordinal, segment in enumerate(boundary_segments):
        footpoint = _closest_point_exact(middle_point, segment)
        squared_distance = _squared_distance(middle_point, footpoint)
        phase_alignment = (footpoint[0] - middle_point[0]) * (phase_point[0] - middle_point[0]) + (footpoint[1] - middle_point[1]) * (phase_point[1] - middle_point[1])
        candidates.append((squared_distance, -phase_alignment, ordinal, footpoint))
    squared_distance, _negated_alignment, _ordinal, footpoint = min(candidates)
    return footpoint, squared_distance


def _adapt_station(
    *,
    path_index: int,
    circle: Circle,
    boundary_segments: tuple[RationalSegment, ...],
    tool_radius: Millimetre,
) -> Figure5ApproximateStation:
    source_center = circle.frame.point
    middle_exact = _rational_point(
        float(source_center.x),
        float(source_center.y),
    )
    phase_exact = _rational_point(
        float(source_center.x) + float(circle.radius) * float(circle.frame.xaxis.x),
        float(source_center.y) + float(circle.radius) * float(circle.frame.xaxis.y),
    )
    footpoint_exact, squared_clearance = _boundary_footpoint(
        middle_point=middle_exact,
        phase_point=phase_exact,
        boundary_segments=boundary_segments,
    )
    tool_radius_exact = Fraction.from_float(float(tool_radius))
    if squared_clearance <= tool_radius_exact * tool_radius_exact:
        raise InvalidFigure5ApproximateStationError("Figure 5 approximate station has no positive paper-derived radius.")

    middle_point = Point2[WorldXY].build(
        float(middle_exact[0]),
        float(middle_exact[1]),
    )
    boundary_footpoint = Point2[WorldXY].build(
        float(footpoint_exact[0]),
        float(footpoint_exact[1]),
    )
    delta_x = float(middle_point.x) - float(boundary_footpoint.x)
    delta_y = float(middle_point.y) - float(boundary_footpoint.y)
    clearance = math.sqrt(float(squared_clearance))
    radius = (clearance - float(tool_radius)) / 2.0
    direction_x = delta_x / clearance
    direction_y = delta_y / clearance
    contact_point = Point2[WorldXY].build(
        float(boundary_footpoint.x) + float(tool_radius) * direction_x,
        float(boundary_footpoint.y) + float(tool_radius) * direction_y,
    )
    center = Point2[WorldXY].build(
        float(contact_point.x) + radius * direction_x,
        float(contact_point.y) + radius * direction_y,
    )
    return Figure5ApproximateStation(
        path_index=path_index,
        middle_point=middle_point,
        boundary_footpoint=boundary_footpoint,
        contact_point=contact_point,
        center=center,
        guide_radius=Millimetre(radius),
    )


def _adapt_figure5_stations(
    case: HeldReferenceCase,
    result: ToolpathResult,
) -> Figure5StationApproximation:
    """Apply paper equations to canonical Figure 5 straight-skeleton stations.

    Args:
        case: The canonical reconstructed Figure 5 reference case.
        result: Existing generator output; only its circular guide stations are
            projected into comparison geometry.

    Returns:
        Ordered, explicitly labeled rendering geometry.

    Raises:
        InvalidFigure5ReproductionInputError: The case/result is not the
            canonical Figure 5 rendering input or has no circular stations.
        InvalidFigure5ApproximateStationError: A station has insufficient
            boundary clearance for a positive ``rho``.
    """
    if type(case) is not HeldReferenceCase or case.name != "figure5":
        raise InvalidFigure5ReproductionInputError("The station approximation is restricted to canonical Figure 5.")
    if type(result) is not ToolpathResult:
        raise InvalidFigure5ReproductionInputError("Figure 5 station approximation requires one ToolpathResult.")

    segments = _boundary_segments(case)
    stations = tuple(
        _adapt_station(
            path_index=operation.path_index,
            circle=operation.geometry,
            boundary_segments=segments,
            tool_radius=Millimetre(float(case.tool_radius.value)),
        )
        for operation in result.operations
        if isinstance(operation.geometry, Circle)
    )
    if not stations:
        raise InvalidFigure5ReproductionInputError("Figure 5 station approximation requires circular guide stations.")
    return Figure5StationApproximation.build(result, stations)


def build_figure5_station_approximation(
    case: HeldReferenceCase,
) -> Figure5StationApproximation:
    """Build the comparison view from the sampled engagement-controlled path.

    The declared 80-degree setting controls only the source generator's sampled
    station selection. This adapter makes no continuous cap or certification
    claim.
    """
    if type(case) is not HeldReferenceCase or case.name != "figure5":
        raise InvalidFigure5ReproductionInputError("The station approximation is restricted to canonical Figure 5.")
    spec = case.pocket_spec()
    result = engagement_controlled_toolpath(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=float(case.tea_cap),
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )
    return _adapt_figure5_stations(case, result)
