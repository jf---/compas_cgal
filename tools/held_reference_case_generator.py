"""Generate the committed Held-Pfeiffer reference-case corpus."""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Mapping
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path
from typing import TypeAlias
from typing import cast

from benchmarks.errors import MalformedHeldReferenceCaseError
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourcePrimitive
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import project_boundary
from benchmarks.held_reference_geometry import reconstruct_source_path
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from tools.held_reference_extractor import FIGURE_CROPS
from tools.held_reference_extractor import ExtractedCase
from tools.held_reference_extractor import FigureCrop
from tools.held_reference_extractor import extract_reference_sources

JsonValue: TypeAlias = None | bool | int | float | str | list["JsonValue"] | dict[str, "JsonValue"]

CANONICAL_FILENAMES = {
    "figure5": "figure5.json",
    "figure8_upper": "figure8_upper.json",
    "figure8_crossed_skis": "figure8_crossed_skis.json",
    "figure8_monstera": "figure8_monstera.json",
}
_SOURCE_TO_CANONICAL = {
    "figure-5": "figure5",
    "figure-8-upper": "figure8_upper",
    "figure-8-skis": "figure8_crossed_skis",
    "figure-8-monstera": "figure8_monstera",
}
_PUBLICATION = {
    "figure5": (742, "5", None),
    "figure8_upper": (746, "8", "upper"),
    "figure8_crossed_skis": (746, "8", "crossed_skis"),
    "figure8_monstera": (746, "8", "monstera_deliciosa"),
}
DATA_DIRECTORY = Path(__file__).parents[1] / "benchmarks" / "data" / "held_pfeiffer_2025"


def _fraction(value: float) -> Fraction:
    return Fraction.from_float(value)


def _source_signed_area_twice(sources: Sequence[SourcePrimitive]) -> Fraction:
    """Return the exact line integral of the stored binary64 source cycle."""
    total = Fraction(0)
    for source in sources:
        if isinstance(source, SourceLine):
            total += _fraction(float(source.start.x)) * _fraction(float(source.end.y))
            total -= _fraction(float(source.start.y)) * _fraction(float(source.end.x))
            continue
        x = _power_coefficients(source.start.x, source.control1.x, source.control2.x, source.end.x)
        y = _power_coefficients(source.start.y, source.control1.y, source.control2.y, source.end.y)
        dx = tuple(Fraction(index) * x[index] for index in range(1, 4))
        dy = tuple(Fraction(index) * y[index] for index in range(1, 4))
        integrand = [Fraction(0) for _ in range(6)]
        for left_index, left in enumerate(x):
            for right_index, right in enumerate(dy):
                integrand[left_index + right_index] += left * right
        for left_index, left in enumerate(y):
            for right_index, right in enumerate(dx):
                integrand[left_index + right_index] -= left * right
        total += sum((coefficient / (index + 1) for index, coefficient in enumerate(integrand)), Fraction(0))
    return total


def _power_coefficients(
    start: PdfPointUnit,
    control1: PdfPointUnit,
    control2: PdfPointUnit,
    end: PdfPointUnit,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    p0, p1, p2, p3 = map(_fraction, map(float, (start, control1, control2, end)))
    return p0, 3 * (p1 - p0), 3 * (p2 - 2 * p1 + p0), p3 - 3 * p2 + 3 * p1 - p0


def _cubic_coordinate(source: SourceCubic, coordinate: str, parameter: float) -> float:
    values = tuple(float(getattr(point, coordinate)) for point in (source.start, source.control1, source.control2, source.end))
    one_minus = 1.0 - parameter
    return one_minus**3 * values[0] + 3.0 * one_minus**2 * parameter * values[1] + 3.0 * one_minus * parameter**2 * values[2] + parameter**3 * values[3]


def _derivative_roots(source: SourceCubic, coordinate: str) -> tuple[float, ...]:
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


def _source_bounds(sources: Sequence[SourcePrimitive]) -> tuple[PdfPoint2, PdfPoint2]:
    xs: list[float] = []
    ys: list[float] = []
    for source in sources:
        xs.extend((float(source.start.x), float(source.end.x)))
        ys.extend((float(source.start.y), float(source.end.y)))
        if isinstance(source, SourceCubic):
            xs.extend(_cubic_coordinate(source, "x", root) for root in _derivative_roots(source, "x"))
            ys.extend(_cubic_coordinate(source, "y", root) for root in _derivative_roots(source, "y"))
    return PdfPoint2.build(min(xs), min(ys)), PdfPoint2.build(max(xs), max(ys))


def _normalization_transform(sources: Sequence[SourcePrimitive], source_tool_radius: PdfPointUnit) -> SourceToWorld:
    minimum, maximum = _source_bounds(sources)
    return SourceToWorld.build(
        source_origin=PdfPoint2.build(minimum.x, maximum.y),
        world_origin=Point2[WorldXY].build(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0 / float(source_tool_radius)),
        reflect_source_y=True,
    )


def _pdf_point(point: PdfPoint2) -> list[JsonValue]:
    return [float(point.x), float(point.y)]


def _world_point(point: Point2[WorldXY]) -> list[JsonValue]:
    return [float(point.x), float(point.y)]


def _source_primitive(source: SourcePrimitive) -> dict[str, JsonValue]:
    if isinstance(source, SourceLine):
        return {"kind": "line", "start": _pdf_point(source.start), "end": _pdf_point(source.end)}
    return {
        "kind": "cubic",
        "start": _pdf_point(source.start),
        "control1": _pdf_point(source.control1),
        "control2": _pdf_point(source.control2),
        "end": _pdf_point(source.end),
    }


def _reference_primitive(primitive: ReferenceLine | ReferenceArc) -> dict[str, JsonValue]:
    if isinstance(primitive, ReferenceLine):
        return {"kind": "line", "start": _world_point(primitive.start), "end": _world_point(primitive.end)}
    return {
        "kind": "arc",
        "start": _world_point(primitive.start),
        "end": _world_point(primitive.end),
        "centre": _world_point(primitive.centre),
        "sweep_rad": float(primitive.sweep),
    }


def _case_payload(extracted: ExtractedCase, crop: FigureCrop) -> dict[str, JsonValue]:
    canonical_name = _SOURCE_TO_CANONICAL.get(extracted.name)
    if canonical_name is None or extracted.page != crop.page or extracted.name != crop.name:
        raise MalformedHeldReferenceCaseError(f"Unexpected extracted case {extracted.name!r}.")
    line_count = sum(isinstance(source, SourceLine) for source in extracted.sources)
    cubic_count = sum(isinstance(source, SourceCubic) for source in extracted.sources)
    if (line_count, cubic_count) != (crop.expected_line_count, crop.expected_cubic_count):
        raise MalformedHeldReferenceCaseError(f"{canonical_name}: extracted primitive census changed.")
    if _source_signed_area_twice(extracted.sources) >= 0:
        raise MalformedHeldReferenceCaseError(f"{canonical_name}: source cycle is not clockwise before reflection.")
    if len(extracted.start_markers) > 1:
        raise MalformedHeldReferenceCaseError(f"{canonical_name}: ambiguous start-marker observation.")

    minimum, maximum = _source_bounds(extracted.sources)
    transform = _normalization_transform(extracted.sources, extracted.tool_circle.radius)
    reconstruction_limit = Millimetre(float(extracted.boundary_stroke_width) * float(transform.scale) / 4.0)
    reconstruction = reconstruct_source_path(extracted.sources, transform, reconstruction_limit)
    normalized_stroke_width = Millimetre(float(extracted.boundary_stroke_width) * float(transform.scale))
    boundary = ReferenceBoundary.build(reconstruction.primitives, ToolRadius.build(1.0), normalized_stroke_width)
    projection_limit = Millimetre(float(reconstruction_limit) / 2.0)
    projection = project_boundary(boundary, projection_limit)
    reconstructed_lines = sum(isinstance(primitive, ReferenceLine) for primitive in reconstruction.primitives)
    reconstructed_arcs = sum(isinstance(primitive, ReferenceArc) for primitive in reconstruction.primitives)
    publication_page, figure, subfigure = _PUBLICATION[canonical_name]
    source_origin = transform.source_origin
    marker: JsonValue = None
    if extracted.start_markers:
        observed = extracted.start_markers[0]
        marker = {
            "centre": _world_point(transform.point(observed.centre)),
            "radius_mm": float(observed.radius) * float(transform.scale),
        }

    return {
        "schema": "held-pfeiffer-reference-case",
        "version": 1,
        "name": canonical_name,
        "citation": {
            "authors": ["Martin Held", "Josef Pfeiffer"],
            "title": "Trochoidal Tool Paths for Pocket Machining with Full Control of the Tool Engagement Angle",
            "doi": "10.14733/cadaps.2025.731-747",
        },
        "publication": {
            "publication_page": publication_page,
            "pdf_page": crop.page,
            "figure": figure,
            "subfigure": subfigure,
        },
        "source": {
            "frame": "publisher_pdf_page",
            "unit": "pdf_point",
            "crop": {"minimum": _pdf_point(crop.crop.minimum), "maximum": _pdf_point(crop.crop.maximum)},
            "bounds": {"minimum": _pdf_point(minimum), "maximum": _pdf_point(maximum)},
            "primitives": [_source_primitive(source) for source in extracted.sources],
        },
        "normalization": {
            "statement": "depicted_tool_radius_equals_1_mm",
            "physical_scale_published": False,
            "source_tool": {
                "centre": _pdf_point(extracted.tool_circle.centre),
                "radius_pdf_point": float(extracted.tool_circle.radius),
            },
            "millimetres_per_pdf_point": float(transform.scale),
            "source_origin": _pdf_point(source_origin),
            "world_origin": _world_point(transform.world_origin),
            "reflect_source_y": True,
        },
        "machining": {"tool_radius_mm": 1.0, "tea_cap_deg": 80.0, "start_marker": marker},
        "analytic_boundary": {
            "frame": "world_xy",
            "unit": "mm",
            "orientation": "CCW",
            "normalized_stroke_width_mm": float(normalized_stroke_width),
            "reconstruction_limit_mm": float(reconstruction_limit),
            "certified_deviation_upper_bound_mm": float(reconstruction.deviation_upper_bound),
            "line_count": reconstructed_lines,
            "arc_count": reconstructed_arcs,
            "primitive_count": len(reconstruction.primitives),
            "primitives": [_reference_primitive(primitive) for primitive in reconstruction.primitives],
        },
        "polygon_projection": {
            "frame": "world_xy",
            "unit": "mm",
            "deviation_limit_mm": float(projection.deviation_limit),
            "observed_deviation_mm": float(projection.observed_deviation),
            "vertex_count": len(projection.points),
            "points": [_world_point(point) for point in projection.points],
        },
        "figure7_observation": (
            {
                "publication_page": 744,
                "pdf_page": 14,
                "panels": ["a", "b", "c"],
                "role": "shape_only_tool_centre_observation",
                "boundary_authority": False,
                "numeric_fidelity_gate": False,
            }
            if canonical_name == "figure5"
            else None
        ),
    }


def _generate_payloads(pdf_path: Path) -> dict[str, JsonValue]:
    extracted_cases = extract_reference_sources(pdf_path)
    crop_by_name = {crop.name: crop for crop in FIGURE_CROPS}
    if tuple(case.name for case in extracted_cases) != tuple(_SOURCE_TO_CANONICAL):
        raise MalformedHeldReferenceCaseError("Publisher extraction did not return the four approved cases in order.")
    return {_SOURCE_TO_CANONICAL[case.name]: _case_payload(case, crop_by_name[case.name]) for case in extracted_cases}


def _reject_duplicate_pairs(pairs: list[tuple[str, JsonValue]]) -> dict[str, JsonValue]:
    result: dict[str, JsonValue] = {}
    for key, value in pairs:
        if key in result:
            raise MalformedHeldReferenceCaseError(f"Duplicate JSON key {key!r}.")
        result[key] = value
    return result


def _decode_semantic(path: Path) -> JsonValue:
    try:
        return cast(
            JsonValue,
            json.loads(
                path.read_text(),
                object_pairs_hook=_reject_duplicate_pairs,
                parse_constant=lambda value: (_ for _ in ()).throw(MalformedHeldReferenceCaseError(f"Non-finite constant {value}.")),
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MalformedHeldReferenceCaseError(f"Cannot decode {path.name}: {error}.") from error


def _write_or_check_payloads(
    payloads: Mapping[str, JsonValue],
    *,
    directory: Path,
    check: bool,
) -> None:
    if tuple(payloads) != tuple(CANONICAL_FILENAMES):
        raise MalformedHeldReferenceCaseError("Generator payload names or order changed.")
    directory.mkdir(parents=True, exist_ok=True)
    for name, filename in CANONICAL_FILENAMES.items():
        payload = payloads[name]
        destination = directory / filename
        if check:
            if _decode_semantic(destination) != payload:
                raise MalformedHeldReferenceCaseError(f"Committed {filename} differs semantically from the publisher PDF.")
            continue
        destination.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("pdf", type=Path)
    arguments = parser.parse_args(argv)
    payloads = _generate_payloads(arguments.pdf)
    _write_or_check_payloads(payloads, directory=DATA_DIRECTORY, check=arguments.check)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
