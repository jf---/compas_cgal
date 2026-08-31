"""Extract Held-Pfeiffer source vectors from the publisher PDF."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from typing import Self
from typing import TypeAlias

from benchmarks.errors import AmbiguousPublishedBoundaryError
from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import MissingPublishedToolCircleError
from benchmarks.errors import UnsupportedPdfBoundaryOperatorError
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine

SourcePrimitive: TypeAlias = SourceLine | SourceCubic

_NUMBER_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_TOKEN_PATTERN = re.compile(rf"[A-Za-z]|{_NUMBER_PATTERN}")
_MATRIX_PATTERN = re.compile(
    rf"matrix\(\s*({_NUMBER_PATTERN})\s*,\s*({_NUMBER_PATTERN})\s*,\s*({_NUMBER_PATTERN})\s*,\s*({_NUMBER_PATTERN})\s*,\s*({_NUMBER_PATTERN})\s*,\s*({_NUMBER_PATTERN})\s*\)"
)

BOUNDARY_GREEN = "rgb(17.999268%, 54.499817%, 34.098816%)"
TOOL_RED = "rgb(100%, 0%, 0%)"
BOUNDARY_STROKE_WIDTH = 2.0
TOOL_CIRCLE_STROKE_WIDTH = 0.8
# Poppler emits page coordinates on the PDF's 1/256-point coordinate grid.
PDF_COORDINATE_QUANTUM_PT = 1.0 / 256.0
# Six-decimal Poppler affine coefficients contribute at most half a printed
# decimal unit when a one-quantum local seam enters the page frame.
PDF_AFFINE_EXPORT_ROUNDOFF_PT = 0.5 / 1_000_000.0
# A diameter compares two independently quantized extrema on each axis.
PDF_DIAMETER_EXTREMA_BOUND_PT = 2.0 * PDF_COORDINATE_QUANTUM_PT


@dataclass(frozen=True)
class AffineTransform:
    """One SVG six-value affine transform in PDF-page coordinates."""

    a: float
    b: float
    c: float
    d: float
    e: float
    f: float

    def __post_init__(self) -> None:
        coefficients = (self.a, self.b, self.c, self.d, self.e, self.f)
        if not all(math.isfinite(value) for value in coefficients):
            raise InvalidPublishedPrimitiveError("An SVG affine transform requires six finite coefficients.")
        if self.a * self.d - self.b * self.c == 0.0:
            raise InvalidPublishedPrimitiveError("An SVG affine transform must be invertible.")

    @classmethod
    def build(
        cls,
        a: float,
        b: float,
        c: float,
        d: float,
        e: float,
        f: float,
        /,
    ) -> Self:
        return cls(a, b, c, d, e, f)

    @classmethod
    def identity(cls) -> Self:
        return cls.build(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)

    def point(self, x: float, y: float) -> PdfPoint2:
        return PdfPoint2.build(
            self.a * x + self.c * y + self.e,
            self.b * x + self.d * y + self.f,
        )


@dataclass(frozen=True)
class PdfCrop:
    """Axis-aligned crop in transformed PDF-page coordinates."""

    minimum: PdfPoint2
    maximum: PdfPoint2

    def __post_init__(self) -> None:
        if float(self.minimum.x) >= float(self.maximum.x) or float(self.minimum.y) >= float(self.maximum.y):
            raise InvalidPublishedPrimitiveError("A PDF crop requires strictly increasing page-coordinate bounds.")

    @classmethod
    def build(cls, minimum_x: float, minimum_y: float, maximum_x: float, maximum_y: float, /) -> Self:
        return cls(PdfPoint2.build(minimum_x, minimum_y), PdfPoint2.build(maximum_x, maximum_y))

    def contains(self, point: PdfPoint2) -> bool:
        return float(self.minimum.x) <= float(point.x) <= float(self.maximum.x) and float(self.minimum.y) <= float(point.y) <= float(self.maximum.y)

    @property
    def bounds(self) -> tuple[PdfPointUnit, PdfPointUnit, PdfPointUnit, PdfPointUnit]:
        return self.minimum.x, self.minimum.y, self.maximum.x, self.maximum.y


@dataclass(frozen=True)
class FigureCrop:
    """Approved page region and vector-family census for one reference case."""

    name: str
    page: int
    crop: PdfCrop
    transform: AffineTransform
    expected_line_count: int
    expected_cubic_count: int
    expected_tool_circle_count: int
    normalization_radius_local: PdfPointUnit

    def __post_init__(self) -> None:
        if not self.name or self.page < 1:
            raise InvalidPublishedPrimitiveError("A figure crop requires a name and positive PDF page.")
        if self.expected_line_count < 0 or self.expected_cubic_count < 0:
            raise InvalidPublishedPrimitiveError("Expected primitive counts cannot be negative.")
        if self.expected_line_count + self.expected_cubic_count == 0:
            raise InvalidPublishedPrimitiveError("A figure crop must describe a non-empty boundary family.")
        if self.expected_tool_circle_count < 1:
            raise InvalidPublishedPrimitiveError("A figure crop must describe at least one normalization circle.")
        if not math.isfinite(float(self.normalization_radius_local)) or float(self.normalization_radius_local) <= 0.0:
            raise InvalidPublishedPrimitiveError("A figure crop requires a finite positive local normalization radius.")

    @classmethod
    def build(
        cls,
        *,
        name: str,
        page: int,
        crop: PdfCrop,
        transform: AffineTransform,
        expected_line_count: int,
        expected_cubic_count: int,
        expected_tool_circle_count: int,
        normalization_radius_local: float,
    ) -> Self:
        return cls(
            name,
            page,
            crop,
            transform,
            expected_line_count,
            expected_cubic_count,
            expected_tool_circle_count,
            PdfPointUnit(normalization_radius_local),
        )


@dataclass(frozen=True)
class SvgPath:
    """One selected SVG path with transforms already applied."""

    stroke_width: float
    primitives: tuple[SourcePrimitive, ...]

    @classmethod
    def build(cls, stroke_width: float, primitives: Sequence[SourcePrimitive]) -> Self:
        numeric_width = float(stroke_width)
        ordered = tuple(primitives)
        if not math.isfinite(numeric_width) or numeric_width <= 0.0 or not ordered:
            raise InvalidPublishedPrimitiveError("A selected SVG path requires positive width and source geometry.")
        return cls(numeric_width, ordered)


@dataclass(frozen=True)
class PublishedToolCircle:
    """The isolated depicted tool circle that establishes normalized scale."""

    centre: PdfPoint2
    radius: PdfPointUnit

    @classmethod
    def build(cls, centre: PdfPoint2, radius: float) -> Self:
        if not math.isfinite(radius) or radius <= 0.0:
            raise InvalidPublishedPrimitiveError("A depicted tool circle requires a finite positive radius.")
        return cls(centre, PdfPointUnit(radius))


@dataclass(frozen=True)
class ExtractedCase:
    """One selected source cycle and its scale/marker observations."""

    name: str
    page: int
    sources: tuple[SourcePrimitive, ...]
    tool_circle: PublishedToolCircle
    start_markers: tuple[PublishedToolCircle, ...]
    boundary_stroke_width: PdfPointUnit

    @classmethod
    def build(
        cls,
        *,
        name: str,
        page: int,
        sources: Sequence[SourcePrimitive],
        tool_circle: PublishedToolCircle,
        start_markers: Sequence[PublishedToolCircle],
        boundary_stroke_width: float,
    ) -> Self:
        ordered = tuple(sources)
        if not ordered or not all(current.end == following.start for current, following in zip(ordered, (*ordered[1:], ordered[0]))):
            raise DisconnectedPublishedBoundaryError("An extracted case requires one exactly closed ordered source cycle.")
        if not math.isfinite(boundary_stroke_width) or boundary_stroke_width <= 0.0:
            raise InvalidPublishedPrimitiveError("Extracted boundary stroke width must be finite and positive.")
        return cls(
            name,
            page,
            ordered,
            tool_circle,
            tuple(start_markers),
            PdfPointUnit(boundary_stroke_width),
        )


# Publisher-page coordinates and primitive censuses measured from Figures 5/8.
FIGURE_CROPS = (
    FigureCrop.build(
        name="figure-5",
        page=12,
        crop=PdfCrop.build(40.0, 70.0, 265.0, 223.0),
        transform=AffineTransform.build(0.476643, 0.0, 0.0, -0.476643, 40.067609, 229.645005),
        expected_line_count=5,
        expected_cubic_count=10,
        expected_tool_circle_count=1,
        normalization_radius_local=6.769349,
    ),
    FigureCrop.build(
        name="figure-8-upper",
        page=16,
        crop=PdfCrop.build(132.0, 70.0, 410.0, 259.0),
        transform=AffineTransform.build(0.604974, 0.0, 0.0, -0.604974, 92.657287, 274.411007),
        expected_line_count=4,
        expected_cubic_count=22,
        expected_tool_circle_count=1,
        normalization_radius_local=5.977675,
    ),
    FigureCrop.build(
        name="figure-8-skis",
        page=16,
        crop=PdfCrop.build(132.0, 253.0, 407.0, 322.0),
        transform=AffineTransform.build(0.485758, 0.0, 0.0, -0.485758, 123.390035, 354.2883),
        expected_line_count=2,
        expected_cubic_count=28,
        expected_tool_circle_count=1,
        normalization_radius_local=5.118423,
    ),
    FigureCrop.build(
        name="figure-8-monstera",
        page=16,
        crop=PdfCrop.build(133.0, 329.0, 407.0, 617.0),
        transform=AffineTransform.build(0.581077, 0.0, 0.0, -0.581077, 120.141587, 620.884297),
        expected_line_count=103,
        expected_cubic_count=107,
        expected_tool_circle_count=1,
        normalization_radius_local=4.863683,
    ),
)


def parse_pdf_svg_path(
    path_data: str,
    transform: AffineTransform,
) -> tuple[SourcePrimitive, ...]:
    """Parse the absolute line/cubic subset emitted by the publisher PDF."""
    tokens = _tokens(path_data)
    if len(tokens) < 4 or tokens[0] != "M":
        raise UnsupportedPdfBoundaryOperatorError("A publisher boundary path must begin with one absolute move.")
    start = transform.point(_number(tokens[1]), _number(tokens[2]))
    operator = tokens[3]
    if operator == "L":
        if len(tokens) != 6:
            raise UnsupportedPdfBoundaryOperatorError("A publisher line path must contain exactly one line operator.")
        return (SourceLine.build(start, transform.point(_number(tokens[4]), _number(tokens[5]))),)
    if operator != "C":
        raise UnsupportedPdfBoundaryOperatorError("Only absolute line and cubic publisher paths are supported.")

    cubics: list[SourceCubic] = []
    cursor = start
    index = 3
    while index < len(tokens):
        if tokens[index] != "C" or index + 6 >= len(tokens):
            raise UnsupportedPdfBoundaryOperatorError("A publisher cubic path has malformed or mixed operators.")
        control1 = transform.point(_number(tokens[index + 1]), _number(tokens[index + 2]))
        control2 = transform.point(_number(tokens[index + 3]), _number(tokens[index + 4]))
        end = transform.point(_number(tokens[index + 5]), _number(tokens[index + 6]))
        cubics.append(SourceCubic.build(cursor, control1, control2, end))
        cursor = end
        index += 7
    if not 1 <= len(cubics) <= 4:
        raise UnsupportedPdfBoundaryOperatorError("A publisher cubic path must contain one through four cubic operators.")
    return tuple(cubics)


def select_boundary_paths(svg_path: Path, crop: FigureCrop) -> tuple[SvgPath, ...]:
    """Select the unique closed boundary component matching an approved crop."""
    candidates = _styled_paths(
        svg_path,
        crop,
        stroke=BOUNDARY_GREEN,
        stroke_width=BOUNDARY_STROKE_WIDTH,
    )
    canonical = _canonicalized_paths(candidates)
    components = _path_components(canonical)
    matches = tuple(component for component in components if _matches_boundary_family(component, crop))
    if len(matches) > 1:
        raise AmbiguousPublishedBoundaryError(f"Crop {crop.name!r} contains more than one matching closed boundary.")
    if not matches:
        raise DisconnectedPublishedBoundaryError(f"Crop {crop.name!r} contains no matching closed boundary.")
    return matches[0]


def select_tool_circle(svg_path: Path, crop: FigureCrop) -> PublishedToolCircle:
    """Select the isolated red normalization circle in an approved crop."""
    candidates = _select_circles(svg_path, crop, TOOL_CIRCLE_STROKE_WIDTH)
    if len(candidates) != crop.expected_tool_circle_count or len(candidates) != 1:
        raise MissingPublishedToolCircleError(f"Crop {crop.name!r} requires one unambiguous normalization circle; found {len(candidates)}.")
    observed = candidates[0]
    scale = math.sqrt(abs(crop.transform.a * crop.transform.d - crop.transform.b * crop.transform.c))
    measured_radius = float(crop.normalization_radius_local) * scale
    if abs(float(observed.radius) - measured_radius) > PDF_DIAMETER_EXTREMA_BOUND_PT:
        raise MissingPublishedToolCircleError(f"Crop {crop.name!r} normalization circle disagrees with its measured radius.")
    return PublishedToolCircle.build(observed.centre, measured_radius)


def extract_case_from_svg(svg_path: Path, crop: FigureCrop) -> ExtractedCase:
    """Extract one approved case from a Poppler SVG page."""
    selected_paths = select_boundary_paths(svg_path, crop)
    sources = _order_cycle(tuple(primitive for path in selected_paths for primitive in path.primitives))
    scale = math.sqrt(abs(crop.transform.a * crop.transform.d - crop.transform.b * crop.transform.c))
    return ExtractedCase.build(
        name=crop.name,
        page=crop.page,
        sources=sources,
        tool_circle=select_tool_circle(svg_path, crop),
        start_markers=_select_circles(svg_path, crop, BOUNDARY_STROKE_WIDTH),
        boundary_stroke_width=BOUNDARY_STROKE_WIDTH * scale,
    )


def extract_reference_sources(pdf_path: Path) -> tuple[ExtractedCase, ...]:
    """Convert the publisher PDF pages and extract all four approved cases."""
    source_pdf = pdf_path.resolve(strict=True)
    build_directory = Path("build")
    build_directory.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="held-reference-", dir=build_directory) as temporary:
        temporary_directory = Path(temporary)
        page_svgs: dict[int, Path] = {}
        for page in sorted({crop.page for crop in FIGURE_CROPS}):
            svg_path = temporary_directory / f"page-{page}.svg"
            subprocess.run(
                [
                    "pdftocairo",
                    "-svg",
                    "-f",
                    str(page),
                    "-l",
                    str(page),
                    str(source_pdf),
                    str(svg_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            page_svgs[page] = svg_path
        return tuple(extract_case_from_svg(page_svgs[crop.page], crop) for crop in FIGURE_CROPS)


def _select_circles(
    svg_path: Path,
    crop: FigureCrop,
    stroke_width: float,
) -> tuple[PublishedToolCircle, ...]:
    candidates: list[PublishedToolCircle] = []
    for path in _styled_paths(
        svg_path,
        crop,
        stroke=TOOL_RED,
        stroke_width=stroke_width,
    ):
        if len(path.primitives) != 4 or not all(isinstance(primitive, SourceCubic) for primitive in path.primitives):
            continue
        endpoints = tuple(primitive.start for primitive in path.primitives)
        if path.primitives[-1].end != path.primitives[0].start:
            continue
        centre_x = sum(float(point.x) for point in endpoints) / len(endpoints)
        centre_y = sum(float(point.y) for point in endpoints) / len(endpoints)
        radii = tuple(math.hypot(float(point.x) - centre_x, float(point.y) - centre_y) for point in endpoints)
        if max(radii) - min(radii) > PDF_DIAMETER_EXTREMA_BOUND_PT:
            continue
        candidates.append(
            PublishedToolCircle.build(
                PdfPoint2.build(centre_x, centre_y),
                min(radii),
            )
        )
    return tuple(candidates)


def canonicalize_degree_one_endpoints(
    primitives: Sequence[SourcePrimitive],
) -> tuple[SourcePrimitive, ...]:
    """Join only mutually unique one-quantum gaps between degree-one endpoints."""
    ordered = tuple(primitives)
    degree = Counter(point for primitive in ordered for point in (primitive.start, primitive.end))
    loose = tuple(point for point, count in degree.items() if count == 1)
    near: dict[PdfPoint2, tuple[PdfPoint2, ...]] = {point: tuple(other for other in loose if other != point and _within_coordinate_quantum(point, other)) for point in loose}
    replacements: dict[PdfPoint2, PdfPoint2] = {}
    for point, neighbours in near.items():
        if len(neighbours) != 1:
            continue
        other = neighbours[0]
        if near.get(other) != (point,):
            continue
        canonical = PdfPoint2.build(
            (float(point.x) + float(other.x)) / 2.0,
            (float(point.y) + float(other.y)) / 2.0,
        )
        replacements[point] = canonical
        replacements[other] = canonical
    return tuple(_replace_primitive_endpoints(primitive, replacements) for primitive in ordered)


def _styled_paths(
    svg_path: Path,
    crop: FigureCrop,
    *,
    stroke: str,
    stroke_width: float,
) -> tuple[SvgPath, ...]:
    root = ET.parse(svg_path).getroot()
    selected: list[SvgPath] = []
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] != "path":
            continue
        if element.get("fill") != "none" or element.get("stroke") != stroke:
            continue
        if element.get("stroke-linejoin") != "round":
            continue
        try:
            width = float(element.get("stroke-width", ""))
        except ValueError:
            continue
        if width != stroke_width:
            continue
        transform_text = element.get("transform")
        path_data = element.get("d")
        if transform_text is None or path_data is None:
            continue
        transform = _parse_transform(transform_text)
        if transform != crop.transform:
            continue
        primitives = parse_pdf_svg_path(path_data, transform)
        if not all(crop.crop.contains(point) for primitive in primitives for point in _primitive_points(primitive)):
            continue
        selected.append(SvgPath.build(width, primitives))
    return tuple(selected)


def _parse_transform(value: str) -> AffineTransform:
    match = _MATRIX_PATTERN.fullmatch(value)
    if match is None:
        raise UnsupportedPdfBoundaryOperatorError("A selected publisher path has an unsupported SVG transform.")
    return AffineTransform.build(*(float(coefficient) for coefficient in match.groups()))


def _primitive_points(primitive: SourcePrimitive) -> tuple[PdfPoint2, ...]:
    if isinstance(primitive, SourceLine):
        return primitive.start, primitive.end
    return primitive.start, primitive.control1, primitive.control2, primitive.end


def _canonicalized_paths(paths: Sequence[SvgPath]) -> tuple[SvgPath, ...]:
    flat = tuple(primitive for path in paths for primitive in path.primitives)
    canonical = canonicalize_degree_one_endpoints(flat)
    result: list[SvgPath] = []
    offset = 0
    for path in paths:
        count = len(path.primitives)
        result.append(SvgPath.build(path.stroke_width, canonical[offset : offset + count]))
        offset += count
    return tuple(result)


def _path_components(paths: Sequence[SvgPath]) -> tuple[tuple[SvgPath, ...], ...]:
    endpoint_paths: dict[PdfPoint2, set[int]] = defaultdict(set)
    for index, path in enumerate(paths):
        for primitive in path.primitives:
            endpoint_paths[primitive.start].add(index)
            endpoint_paths[primitive.end].add(index)
    adjacency: dict[int, set[int]] = defaultdict(set)
    for incident in endpoint_paths.values():
        for index in incident:
            adjacency[index].update(incident - {index})
    unseen = set(range(len(paths)))
    components: list[tuple[SvgPath, ...]] = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        indexes: set[int] = set()
        while stack:
            index = stack.pop()
            if index in indexes:
                continue
            indexes.add(index)
            stack.extend(adjacency[index] - indexes)
        unseen -= indexes
        components.append(tuple(paths[index] for index in sorted(indexes)))
    return tuple(components)


def _matches_boundary_family(paths: Sequence[SvgPath], crop: FigureCrop) -> bool:
    primitives = tuple(primitive for path in paths for primitive in path.primitives)
    endpoint_degree = Counter(point for primitive in primitives for point in (primitive.start, primitive.end))
    if not endpoint_degree or any(degree != 2 for degree in endpoint_degree.values()):
        return False
    line_count = sum(isinstance(primitive, SourceLine) for primitive in primitives)
    cubic_count = sum(isinstance(primitive, SourceCubic) for primitive in primitives)
    return line_count == crop.expected_line_count and cubic_count == crop.expected_cubic_count


def _within_coordinate_quantum(first: PdfPoint2, second: PdfPoint2) -> bool:
    page_quantum_bound = PDF_COORDINATE_QUANTUM_PT + PDF_AFFINE_EXPORT_ROUNDOFF_PT
    return abs(float(first.x) - float(second.x)) <= page_quantum_bound and abs(float(first.y) - float(second.y)) <= page_quantum_bound


def _replace_primitive_endpoints(
    primitive: SourcePrimitive,
    replacements: dict[PdfPoint2, PdfPoint2],
) -> SourcePrimitive:
    start = replacements.get(primitive.start, primitive.start)
    end = replacements.get(primitive.end, primitive.end)
    if isinstance(primitive, SourceLine):
        return SourceLine.build(start, end)
    return SourceCubic.build(start, primitive.control1, primitive.control2, end)


def _order_cycle(primitives: Sequence[SourcePrimitive]) -> tuple[SourcePrimitive, ...]:
    if not primitives:
        raise DisconnectedPublishedBoundaryError("A published source cycle cannot be empty.")
    unused = list(primitives)
    ordered = [unused.pop(0)]
    while unused:
        endpoint = ordered[-1].end
        incident = tuple((index, primitive) for index, primitive in enumerate(unused) if primitive.start == endpoint or primitive.end == endpoint)
        if len(incident) != 1:
            raise DisconnectedPublishedBoundaryError("A published source cycle branches or is disconnected.")
        index, primitive = incident[0]
        unused.pop(index)
        ordered.append(primitive if primitive.start == endpoint else _reverse_primitive(primitive))
    if ordered[-1].end != ordered[0].start:
        raise DisconnectedPublishedBoundaryError("A published source cycle does not close exactly.")
    return tuple(ordered)


def _reverse_primitive(primitive: SourcePrimitive) -> SourcePrimitive:
    if isinstance(primitive, SourceLine):
        return SourceLine.build(primitive.end, primitive.start)
    return SourceCubic.build(primitive.end, primitive.control2, primitive.control1, primitive.start)


def _case_payload(case: ExtractedCase) -> dict[str, object]:
    return {
        "name": case.name,
        "page": case.page,
        "source_count": len(case.sources),
        "line_count": sum(isinstance(source, SourceLine) for source in case.sources),
        "cubic_count": sum(isinstance(source, SourceCubic) for source in case.sources),
        "tool_radius_pdf_points": float(case.tool_circle.radius),
        "boundary_stroke_width_pdf_points": float(case.boundary_stroke_width),
        "start_marker_count": len(case.start_markers),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("publisher_pdf", type=Path)
    arguments = parser.parse_args(argv)
    print(json.dumps({"cases": [_case_payload(case) for case in extract_reference_sources(arguments.publisher_pdf)]}, indent=2))
    return 0


def _tokens(path_data: str) -> tuple[str, ...]:
    matches = tuple(_TOKEN_PATTERN.finditer(path_data))
    cursor = 0
    tokens: list[str] = []
    for match in matches:
        if path_data[cursor : match.start()].strip(" ,\t\r\n"):
            raise UnsupportedPdfBoundaryOperatorError("A publisher path contains an unsupported token.")
        tokens.append(match.group())
        cursor = match.end()
    if path_data[cursor:].strip(" ,\t\r\n"):
        raise UnsupportedPdfBoundaryOperatorError("A publisher path contains an unsupported token.")
    return tuple(tokens)


def _number(token: str) -> float:
    if re.fullmatch(_NUMBER_PATTERN, token) is None:
        raise UnsupportedPdfBoundaryOperatorError("A publisher path coordinate is malformed.")
    value = float(token)
    if not math.isfinite(value):
        raise UnsupportedPdfBoundaryOperatorError("A publisher path coordinate must be finite.")
    return value


if __name__ == "__main__":
    sys.exit(main())
