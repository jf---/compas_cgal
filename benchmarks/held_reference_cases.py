"""Typed loader for the committed Held-Pfeiffer reference cases."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import NewType
from typing import Self
from typing import TypeAlias
from typing import cast

import jsonschema
from compas.geometry import Point
from compas.geometry import Polygon

from benchmarks.errors import BenchmarkError
from benchmarks.errors import MalformedHeldReferenceCaseError
from benchmarks.errors import UnknownHeldReferenceCaseError
from benchmarks.errors import UnsupportedHeldReferenceVersionError
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import PolygonProjection
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from benchmarks.held_reference_geometry import ReferenceReconstruction
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourcePrimitive
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import project_boundary
from benchmarks.held_reference_geometry import reconstruct_source_path
from benchmarks.spec import PocketSpec
from compas_cgal.adaptive.errors import InvalidUnitValueError
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

Degree = NewType("Degree", float)
JsonValue: TypeAlias = None | bool | int | float | str | list["JsonValue"] | dict[str, "JsonValue"]

CANONICAL_CASE_NAMES = (
    "figure5",
    "figure8_upper",
    "figure8_crossed_skis",
    "figure8_monstera",
)
DATA_DIRECTORY = Path(__file__).with_name("data") / "held_pfeiffer_2025"
_CITATION_AUTHORS = ("Martin Held", "Josef Pfeiffer")
_CITATION_TITLE = "Trochoidal Tool Paths for Pocket Machining with Full Control of the Tool Engagement Angle"
_CITATION_DOI = "10.14733/cadaps.2025.731-747"
_CASE_METADATA = {
    "figure5": (742, 12, "5", None, (40.0, 70.0, 265.0, 223.0), (5, 10)),
    "figure8_upper": (746, 16, "8", "upper", (132.0, 70.0, 410.0, 259.0), (4, 22)),
    "figure8_crossed_skis": (746, 16, "8", "crossed_skis", (132.0, 253.0, 407.0, 322.0), (2, 28)),
    "figure8_monstera": (746, 16, "8", "monstera_deliciosa", (133.0, 329.0, 407.0, 617.0), (103, 107)),
}

_POINT = {
    "type": "array",
    "prefixItems": [{"type": "number"}, {"type": "number"}],
    "items": False,
    "minItems": 2,
    "maxItems": 2,
}


def _closed(properties: Mapping[str, object], required: tuple[str, ...]) -> dict[str, object]:
    return {
        "type": "object",
        "properties": dict(properties),
        "required": list(required),
        "additionalProperties": False,
    }


_SOURCE_LINE = _closed(
    {"kind": {"const": "line"}, "start": _POINT, "end": _POINT},
    ("kind", "start", "end"),
)
_SOURCE_CUBIC = _closed(
    {
        "kind": {"const": "cubic"},
        "start": _POINT,
        "control1": _POINT,
        "control2": _POINT,
        "end": _POINT,
    },
    ("kind", "start", "control1", "control2", "end"),
)
_REFERENCE_LINE = _closed(
    {"kind": {"const": "line"}, "start": _POINT, "end": _POINT},
    ("kind", "start", "end"),
)
_REFERENCE_ARC = _closed(
    {"kind": {"const": "arc"}, "start": _POINT, "end": _POINT, "centre": _POINT, "sweep_rad": {"type": "number"}},
    ("kind", "start", "end", "centre", "sweep_rad"),
)

CASE_SCHEMA = _closed(
    {
        "schema": {"const": "held-pfeiffer-reference-case"},
        "version": {"type": "integer"},
        "name": {"enum": list(CANONICAL_CASE_NAMES)},
        "citation": _closed(
            {
                "authors": {"type": "array", "items": {"type": "string"}, "minItems": 1},
                "title": {"type": "string", "minLength": 1},
                "doi": {"type": "string", "minLength": 1},
            },
            ("authors", "title", "doi"),
        ),
        "publication": _closed(
            {
                "publication_page": {"type": "integer", "minimum": 1},
                "pdf_page": {"type": "integer", "minimum": 1},
                "figure": {"type": "string", "minLength": 1},
                "subfigure": {"type": ["string", "null"]},
            },
            ("publication_page", "pdf_page", "figure", "subfigure"),
        ),
        "source": _closed(
            {
                "frame": {"const": "publisher_pdf_page"},
                "unit": {"const": "pdf_point"},
                "crop": _closed({"minimum": _POINT, "maximum": _POINT}, ("minimum", "maximum")),
                "bounds": _closed({"minimum": _POINT, "maximum": _POINT}, ("minimum", "maximum")),
                "primitives": {"type": "array", "items": {"oneOf": [_SOURCE_LINE, _SOURCE_CUBIC]}, "minItems": 1},
            },
            ("frame", "unit", "crop", "bounds", "primitives"),
        ),
        "normalization": _closed(
            {
                "statement": {"const": "depicted_tool_radius_equals_1_mm"},
                "physical_scale_published": {"const": False},
                "source_tool": _closed(
                    {"centre": _POINT, "radius_pdf_point": {"type": "number", "exclusiveMinimum": 0}},
                    ("centre", "radius_pdf_point"),
                ),
                "millimetres_per_pdf_point": {"type": "number", "exclusiveMinimum": 0},
                "source_origin": _POINT,
                "world_origin": _POINT,
                "reflect_source_y": {"const": True},
            },
            (
                "statement",
                "physical_scale_published",
                "source_tool",
                "millimetres_per_pdf_point",
                "source_origin",
                "world_origin",
                "reflect_source_y",
            ),
        ),
        "machining": _closed(
            {
                "tool_radius_mm": {"const": 1.0},
                "tea_cap_deg": {"const": 80.0},
                "start_marker": {
                    "oneOf": [
                        {"type": "null"},
                        _closed(
                            {"centre": _POINT, "radius_mm": {"type": "number", "exclusiveMinimum": 0}},
                            ("centre", "radius_mm"),
                        ),
                    ]
                },
            },
            ("tool_radius_mm", "tea_cap_deg", "start_marker"),
        ),
        "analytic_boundary": _closed(
            {
                "frame": {"const": "world_xy"},
                "unit": {"const": "mm"},
                "orientation": {"const": "CCW"},
                "normalized_stroke_width_mm": {"type": "number", "exclusiveMinimum": 0},
                "reconstruction_limit_mm": {"type": "number", "exclusiveMinimum": 0},
                "certified_deviation_upper_bound_mm": {"type": "number", "minimum": 0},
                "line_count": {"type": "integer", "minimum": 0},
                "arc_count": {"type": "integer", "minimum": 0},
                "primitive_count": {"type": "integer", "minimum": 1},
                "primitives": {"type": "array", "items": {"oneOf": [_REFERENCE_LINE, _REFERENCE_ARC]}, "minItems": 1},
            },
            (
                "frame",
                "unit",
                "orientation",
                "normalized_stroke_width_mm",
                "reconstruction_limit_mm",
                "certified_deviation_upper_bound_mm",
                "line_count",
                "arc_count",
                "primitive_count",
                "primitives",
            ),
        ),
        "polygon_projection": _closed(
            {
                "frame": {"const": "world_xy"},
                "unit": {"const": "mm"},
                "deviation_limit_mm": {"type": "number", "exclusiveMinimum": 0},
                "observed_deviation_mm": {"type": "number", "minimum": 0},
                "vertex_count": {"type": "integer", "minimum": 3},
                "points": {"type": "array", "items": _POINT, "minItems": 3},
            },
            ("frame", "unit", "deviation_limit_mm", "observed_deviation_mm", "vertex_count", "points"),
        ),
        "figure7_observation": {
            "oneOf": [
                {"type": "null"},
                _closed(
                    {
                        "publication_page": {"type": "integer", "minimum": 1},
                        "pdf_page": {"type": "integer", "minimum": 1},
                        "panels": {"type": "array", "items": {"type": "string"}, "minItems": 1},
                        "role": {"const": "shape_only_tool_centre_observation"},
                        "boundary_authority": {"const": False},
                        "numeric_fidelity_gate": {"const": False},
                    },
                    ("publication_page", "pdf_page", "panels", "role", "boundary_authority", "numeric_fidelity_gate"),
                ),
            ]
        },
    },
    (
        "schema",
        "version",
        "name",
        "citation",
        "publication",
        "source",
        "normalization",
        "machining",
        "analytic_boundary",
        "polygon_projection",
        "figure7_observation",
    ),
)


@dataclass(frozen=True)
class Figure7Observation:
    publication_page: int
    pdf_page: int
    panels: tuple[str, ...]
    role: str
    boundary_authority: bool
    numeric_fidelity_gate: bool

    @classmethod
    def build(
        cls,
        *,
        publication_page: int,
        pdf_page: int,
        panels: tuple[str, ...],
        role: str,
        boundary_authority: bool,
        numeric_fidelity_gate: bool,
    ) -> Self:
        if publication_page < 1 or pdf_page < 1 or not panels:
            raise MalformedHeldReferenceCaseError("Figure 7 observation requires pages and named panels.")
        if role != "shape_only_tool_centre_observation" or boundary_authority or numeric_fidelity_gate:
            raise MalformedHeldReferenceCaseError("Figure 7 is shape-only non-authoritative evidence.")
        return cls(publication_page, pdf_page, panels, role, boundary_authority, numeric_fidelity_gate)


@dataclass(frozen=True)
class HeldReferenceCase:
    name: str
    authors: tuple[str, ...]
    title: str
    doi: str
    publication_page: int
    pdf_page: int
    figure: str
    subfigure: str | None
    source_crop: tuple[PdfPoint2, PdfPoint2]
    boundary: ReferenceBoundary
    reconstruction: ReferenceReconstruction
    projection: PolygonProjection
    tool_radius: ToolRadius
    tea_cap: Degree
    start_marker: Point2[WorldXY] | None
    figure7_observation: Figure7Observation | None

    @classmethod
    def build(
        cls,
        *,
        name: str,
        authors: tuple[str, ...],
        title: str,
        doi: str,
        publication_page: int,
        pdf_page: int,
        figure: str,
        subfigure: str | None,
        source_crop: tuple[PdfPoint2, PdfPoint2],
        boundary: ReferenceBoundary,
        reconstruction: ReferenceReconstruction,
        projection: PolygonProjection,
        tool_radius: ToolRadius,
        tea_cap: Degree,
        start_marker: Point2[WorldXY] | None,
        figure7_observation: Figure7Observation | None,
    ) -> Self:
        if name not in CANONICAL_CASE_NAMES or not authors or not title or not doi:
            raise MalformedHeldReferenceCaseError("Reference case publication identity is incomplete.")
        if publication_page < 1 or pdf_page < 1 or not figure:
            raise MalformedHeldReferenceCaseError("Reference case publication location is invalid.")
        case = cls(
            name,
            authors,
            title,
            doi,
            publication_page,
            pdf_page,
            figure,
            subfigure,
            source_crop,
            boundary,
            reconstruction,
            projection,
            tool_radius,
            tea_cap,
            start_marker,
            figure7_observation,
        )
        if case.analytic_signed_area <= 0.0 or case.polygon_signed_area <= 0.0:
            raise MalformedHeldReferenceCaseError(f"{name}: analytic and polygon boundaries must both be CCW.")
        if (figure7_observation is not None) != (name == "figure5"):
            raise MalformedHeldReferenceCaseError("Figure 7 metadata belongs only to Figure 5.")
        return case

    @property
    def projection_vertex_count(self) -> int:
        return len(self.projection.points)

    @property
    def analytic_signed_area(self) -> float:
        twice_area = 0.0
        for primitive in self.boundary.primitives:
            if isinstance(primitive, ReferenceLine):
                twice_area += float(primitive.start.x) * float(primitive.end.y) - float(primitive.start.y) * float(primitive.end.x)
            else:
                twice_area += (
                    float(primitive.centre.x) * (float(primitive.end.y) - float(primitive.start.y))
                    - float(primitive.centre.y) * (float(primitive.end.x) - float(primitive.start.x))
                    + _arc_radius_squared(primitive) * float(primitive.sweep)
                )
        return twice_area / 2.0

    @property
    def polygon_signed_area(self) -> float:
        points = self.projection.points
        return sum(float(start.x) * float(end.y) - float(start.y) * float(end.x) for start, end in zip(points, (*points[1:], points[0]))) / 2.0

    def pocket_spec(self) -> PocketSpec:
        if self.analytic_signed_area <= 0.0 or self.polygon_signed_area <= 0.0:
            raise MalformedHeldReferenceCaseError(f"{self.name}: PocketSpec requires independently verified CCW rings.")
        polygon = Polygon([Point(float(point.x), float(point.y), 0.0) for point in self.projection.points])
        return PocketSpec.build(
            self.name,
            "held_pfeiffer_2025",
            polygon,
            2.0,
            80.0,
            holes=(),
        )


def _arc_radius_squared(arc: ReferenceArc) -> float:
    return (float(arc.start.x) - float(arc.centre.x)) ** 2 + (float(arc.start.y) - float(arc.centre.y)) ** 2


def _reject_duplicate_pairs(pairs: list[tuple[str, JsonValue]]) -> dict[str, JsonValue]:
    result: dict[str, JsonValue] = {}
    for key, value in pairs:
        if key in result:
            raise MalformedHeldReferenceCaseError(f"Duplicate JSON key {key!r}.")
        result[key] = value
    return result


def _reject_constant(value: str) -> float:
    raise MalformedHeldReferenceCaseError(f"Non-finite JSON constant {value!r}.")


def _mapping(value: object) -> dict[str, JsonValue]:
    if not isinstance(value, dict):
        raise MalformedHeldReferenceCaseError("Expected a JSON object.")
    return cast(dict[str, JsonValue], value)


def _sequence(value: object) -> list[JsonValue]:
    if not isinstance(value, list):
        raise MalformedHeldReferenceCaseError("Expected a JSON array.")
    return cast(list[JsonValue], value)


def _number(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MalformedHeldReferenceCaseError("Expected a JSON number.")
    return float(value)


def _integer(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise MalformedHeldReferenceCaseError("Expected a JSON integer.")
    return value


def _require_finite_json(value: JsonValue) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise MalformedHeldReferenceCaseError("JSON numbers must be finite.")
    if isinstance(value, list):
        for item in value:
            _require_finite_json(item)
    elif isinstance(value, dict):
        for item in value.values():
            _require_finite_json(item)


def _point_pdf(value: object) -> PdfPoint2:
    point = _sequence(value)
    return PdfPoint2.build(_number(point[0]), _number(point[1]))


def _point_world(value: object) -> Point2[WorldXY]:
    point = _sequence(value)
    return Point2[WorldXY].build(_number(point[0]), _number(point[1]))


def _parse_sources(values: object) -> tuple[SourcePrimitive, ...]:
    parsed: list[SourcePrimitive] = []
    for value in _sequence(values):
        primitive = _mapping(value)
        if primitive["kind"] == "line":
            parsed.append(SourceLine.build(_point_pdf(primitive["start"]), _point_pdf(primitive["end"])))
        elif primitive["kind"] == "cubic":
            parsed.append(
                SourceCubic.build(
                    _point_pdf(primitive["start"]),
                    _point_pdf(primitive["control1"]),
                    _point_pdf(primitive["control2"]),
                    _point_pdf(primitive["end"]),
                )
            )
        else:
            raise MalformedHeldReferenceCaseError(f"Unknown source primitive {primitive['kind']!r}.")
    return tuple(parsed)


def _source_coordinate(source: SourceCubic, coordinate: str, parameter: float) -> float:
    values = tuple(float(getattr(point, coordinate)) for point in (source.start, source.control1, source.control2, source.end))
    complement = 1.0 - parameter
    return complement**3 * values[0] + 3.0 * complement**2 * parameter * values[1] + 3.0 * complement * parameter**2 * values[2] + parameter**3 * values[3]


def _source_derivative_roots(source: SourceCubic, coordinate: str) -> tuple[float, ...]:
    p0, p1, p2, p3 = (float(getattr(point, coordinate)) for point in (source.start, source.control1, source.control2, source.end))
    a = -p0 + 3.0 * p1 - 3.0 * p2 + p3
    b = 2.0 * (p0 - 2.0 * p1 + p2)
    c = p1 - p0
    if a == 0.0:
        if b == 0.0:
            return ()
        root = -c / b
        return (root,) if 0.0 < root < 1.0 else ()
    discriminant = b * b - 4.0 * a * c
    if discriminant <= 0.0:
        return ()
    square_root = math.sqrt(discriminant)
    numerator = -0.5 * (b + math.copysign(square_root, b))
    roots = (numerator / a, c / numerator) if numerator != 0.0 else (-b / (2.0 * a),)
    return tuple(sorted({root for root in roots if 0.0 < root < 1.0}))


def _source_bounds(sources: tuple[SourcePrimitive, ...]) -> tuple[PdfPoint2, PdfPoint2]:
    xs: list[float] = []
    ys: list[float] = []
    for source in sources:
        xs.extend((float(source.start.x), float(source.end.x)))
        ys.extend((float(source.start.y), float(source.end.y)))
        if isinstance(source, SourceCubic):
            xs.extend(_source_coordinate(source, "x", root) for root in _source_derivative_roots(source, "x"))
            ys.extend(_source_coordinate(source, "y", root) for root in _source_derivative_roots(source, "y"))
    return PdfPoint2.build(min(xs), min(ys)), PdfPoint2.build(max(xs), max(ys))


def _source_signed_area_twice(sources: tuple[SourcePrimitive, ...]) -> Fraction:
    total = Fraction(0)
    for source in sources:
        if isinstance(source, SourceLine):
            total += Fraction.from_float(float(source.start.x)) * Fraction.from_float(float(source.end.y))
            total -= Fraction.from_float(float(source.start.y)) * Fraction.from_float(float(source.end.x))
            continue
        x = _source_power(source.start.x, source.control1.x, source.control2.x, source.end.x)
        y = _source_power(source.start.y, source.control1.y, source.control2.y, source.end.y)
        dx = tuple(Fraction(index) * x[index] for index in range(1, 4))
        dy = tuple(Fraction(index) * y[index] for index in range(1, 4))
        coefficients = [Fraction(0) for _ in range(6)]
        for left_index, left in enumerate(x):
            for right_index, right in enumerate(dy):
                coefficients[left_index + right_index] += left * right
        for left_index, left in enumerate(y):
            for right_index, right in enumerate(dx):
                coefficients[left_index + right_index] -= left * right
        total += sum((coefficient / (index + 1) for index, coefficient in enumerate(coefficients)), Fraction(0))
    return total


def _source_power(
    start: PdfPointUnit,
    control1: PdfPointUnit,
    control2: PdfPointUnit,
    end: PdfPointUnit,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    p0, p1, p2, p3 = (Fraction.from_float(float(value)) for value in (start, control1, control2, end))
    return p0, 3 * (p1 - p0), 3 * (p2 - 2 * p1 + p0), p3 - 3 * p2 + 3 * p1 - p0


def _parse_reference(values: object) -> tuple[ReferencePrimitive, ...]:
    parsed: list[ReferencePrimitive] = []
    for value in _sequence(values):
        primitive = _mapping(value)
        if primitive["kind"] == "line":
            parsed.append(ReferenceLine.build(_point_world(primitive["start"]), _point_world(primitive["end"])))
        elif primitive["kind"] == "arc":
            parsed.append(
                ReferenceArc.build(
                    _point_world(primitive["start"]),
                    _point_world(primitive["end"]),
                    _point_world(primitive["centre"]),
                    Radian(_number(primitive["sweep_rad"])),
                )
            )
        else:
            raise MalformedHeldReferenceCaseError(f"Unknown reference primitive {primitive['kind']!r}.")
    return tuple(parsed)


def _same_float(recorded: object, computed: float) -> bool:
    return _number(recorded) == computed


def validate_case_payload(payload: object, *, expected_name: str) -> HeldReferenceCase:
    document = _mapping(payload)
    _require_finite_json(document)
    if "version" in document and document["version"] != 1:
        raise UnsupportedHeldReferenceVersionError(f"Unsupported Held-Pfeiffer version {document.get('version')!r}.")
    try:
        jsonschema.validate(document, CASE_SCHEMA)
        if document["name"] != expected_name:
            raise MalformedHeldReferenceCaseError(f"Filename {expected_name!r} contains case {document['name']!r}.")
        source = _mapping(document["source"])
        normalization = _mapping(document["normalization"])
        machining = _mapping(document["machining"])
        analytic = _mapping(document["analytic_boundary"])
        projection_record = _mapping(document["polygon_projection"])
        citation = _mapping(document["citation"])
        publication = _mapping(document["publication"])
        sources = _parse_sources(source["primitives"])
        name = str(document["name"])
        publication_page, pdf_page, figure, subfigure, crop_coordinates, source_census = _CASE_METADATA[name]
        if tuple(str(author) for author in _sequence(citation["authors"])) != _CITATION_AUTHORS or citation["title"] != _CITATION_TITLE or citation["doi"] != _CITATION_DOI:
            raise MalformedHeldReferenceCaseError("Citation differs from the approved publication.")
        if (
            publication["publication_page"],
            publication["pdf_page"],
            publication["figure"],
            publication["subfigure"],
        ) != (publication_page, pdf_page, figure, subfigure):
            raise MalformedHeldReferenceCaseError("Publication location differs from the approved figure.")
        if (
            sum(isinstance(source_primitive, SourceLine) for source_primitive in sources),
            sum(isinstance(source_primitive, SourceCubic) for source_primitive in sources),
        ) != source_census:
            raise MalformedHeldReferenceCaseError("Publisher source primitive census changed.")
        if _source_signed_area_twice(sources) >= 0:
            raise MalformedHeldReferenceCaseError("Publisher source cycle must be clockwise before Y reflection.")
        source_bounds = _mapping(source["bounds"])
        minimum, maximum = _source_bounds(sources)
        if (_point_pdf(source_bounds["minimum"]), _point_pdf(source_bounds["maximum"])) != (minimum, maximum):
            raise MalformedHeldReferenceCaseError("Recorded source bounds differ from the true cubic extrema.")
        required_origin = PdfPoint2.build(minimum.x, maximum.y)
        if _point_pdf(normalization["source_origin"]) != required_origin or _point_world(normalization["world_origin"]) != Point2[WorldXY].build(0.0, 0.0):
            raise MalformedHeldReferenceCaseError("Normalization must place the reflected source bounds at the world origin.")
        transform = SourceToWorld.build(
            source_origin=_point_pdf(normalization["source_origin"]),
            world_origin=_point_world(normalization["world_origin"]),
            scale=MillimetresPerPdfPoint(_number(normalization["millimetres_per_pdf_point"])),
            reflect_source_y=bool(normalization["reflect_source_y"]),
        )
        source_tool = _mapping(normalization["source_tool"])
        if not _same_float(normalization["millimetres_per_pdf_point"], 1.0 / _number(source_tool["radius_pdf_point"])):
            raise MalformedHeldReferenceCaseError("Normalization scale does not make the depicted radius one millimetre.")
        if not _same_float(analytic["normalized_stroke_width_mm"], 4.0 * _number(analytic["reconstruction_limit_mm"])):
            raise MalformedHeldReferenceCaseError("Reconstruction limit must equal one quarter of the normalized stroke width.")
        if not _same_float(projection_record["deviation_limit_mm"], _number(analytic["reconstruction_limit_mm"]) / 2.0):
            raise MalformedHeldReferenceCaseError("Projection limit must be half the reconstruction limit.")
        reconstruction = reconstruct_source_path(sources, transform, Millimetre(_number(analytic["reconstruction_limit_mm"])))
        recorded_primitives = _parse_reference(analytic["primitives"])
        recorded_reconstruction = ReferenceReconstruction.build(
            recorded_primitives,
            Millimetre(_number(analytic["certified_deviation_upper_bound_mm"])),
        )
        if reconstruction != recorded_reconstruction:
            raise MalformedHeldReferenceCaseError("Recorded analytic reconstruction differs from deterministic source evidence.")
        line_count = sum(isinstance(primitive, ReferenceLine) for primitive in recorded_primitives)
        arc_count = sum(isinstance(primitive, ReferenceArc) for primitive in recorded_primitives)
        if (line_count, arc_count, len(recorded_primitives)) != (
            analytic["line_count"],
            analytic["arc_count"],
            analytic["primitive_count"],
        ):
            raise MalformedHeldReferenceCaseError("Recorded analytic primitive census is inconsistent.")
        tool_radius = ToolRadius.build(_number(machining["tool_radius_mm"]))
        boundary = ReferenceBoundary.build(
            recorded_primitives,
            tool_radius,
            Millimetre(_number(analytic["normalized_stroke_width_mm"])),
        )
        computed_projection = project_boundary(boundary, Millimetre(_number(projection_record["deviation_limit_mm"])))
        recorded_points = tuple(_point_world(point) for point in _sequence(projection_record["points"]))
        recorded_projection = PolygonProjection.build(
            recorded_points,
            Millimetre(_number(projection_record["deviation_limit_mm"])),
            Millimetre(_number(projection_record["observed_deviation_mm"])),
        )
        if computed_projection != recorded_projection or len(recorded_points) != projection_record["vertex_count"]:
            raise MalformedHeldReferenceCaseError("Recorded polygon differs from deterministic projection evidence.")
        figure7_record = document["figure7_observation"]
        figure7 = None
        if figure7_record is not None:
            observed = _mapping(figure7_record)
            figure7 = Figure7Observation.build(
                publication_page=_integer(observed["publication_page"]),
                pdf_page=_integer(observed["pdf_page"]),
                panels=tuple(str(panel) for panel in _sequence(observed["panels"])),
                role=str(observed["role"]),
                boundary_authority=bool(observed["boundary_authority"]),
                numeric_fidelity_gate=bool(observed["numeric_fidelity_gate"]),
            )
            if (figure7.publication_page, figure7.pdf_page, figure7.panels) != (744, 14, ("a", "b", "c")):
                raise MalformedHeldReferenceCaseError("Figure 7 location or panels differ from the approved observation.")
        marker_record = machining["start_marker"]
        start_marker = None if marker_record is None else _point_world(_mapping(marker_record)["centre"])
        crop = _mapping(source["crop"])
        crop_minimum = _point_pdf(crop["minimum"])
        crop_maximum = _point_pdf(crop["maximum"])
        if (float(crop_minimum.x), float(crop_minimum.y), float(crop_maximum.x), float(crop_maximum.y)) != crop_coordinates:
            raise MalformedHeldReferenceCaseError("Source crop differs from the approved publisher-page crop.")
        return HeldReferenceCase.build(
            name=str(document["name"]),
            authors=tuple(str(author) for author in _sequence(citation["authors"])),
            title=str(citation["title"]),
            doi=str(citation["doi"]),
            publication_page=_integer(publication["publication_page"]),
            pdf_page=_integer(publication["pdf_page"]),
            figure=str(publication["figure"]),
            subfigure=None if publication["subfigure"] is None else str(publication["subfigure"]),
            source_crop=(crop_minimum, crop_maximum),
            boundary=boundary,
            reconstruction=recorded_reconstruction,
            projection=recorded_projection,
            tool_radius=tool_radius,
            tea_cap=Degree(_number(machining["tea_cap_deg"])),
            start_marker=start_marker,
            figure7_observation=figure7,
        )
    except (jsonschema.ValidationError, BenchmarkError, InvalidUnitValueError, KeyError, TypeError, ValueError, OverflowError) as error:
        if isinstance(error, (MalformedHeldReferenceCaseError, UnsupportedHeldReferenceVersionError)):
            raise
        raise MalformedHeldReferenceCaseError(f"{expected_name}: malformed Held-Pfeiffer case: {error}.") from error


def decode_case_document(data: bytes, *, expected_name: str) -> HeldReferenceCase:
    try:
        payload = json.loads(data, object_pairs_hook=_reject_duplicate_pairs, parse_constant=_reject_constant)
    except (UnicodeError, json.JSONDecodeError) as error:
        raise MalformedHeldReferenceCaseError(f"{expected_name}: invalid JSON: {error}.") from error
    return validate_case_payload(payload, expected_name=expected_name)


def load_held_reference_case(name: str) -> HeldReferenceCase:
    if name not in CANONICAL_CASE_NAMES:
        raise UnknownHeldReferenceCaseError(f"Unknown Held-Pfeiffer reference case {name!r}.")
    try:
        data = (DATA_DIRECTORY / f"{name}.json").read_bytes()
    except OSError as error:
        raise MalformedHeldReferenceCaseError(f"Cannot read committed case {name!r}: {error}.") from error
    return decode_case_document(data, expected_name=name)


def load_all_held_reference_cases() -> tuple[HeldReferenceCase, ...]:
    return tuple(load_held_reference_case(name) for name in CANONICAL_CASE_NAMES)
