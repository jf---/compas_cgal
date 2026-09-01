"""Semantic overlays for the Held-Pfeiffer reference corpus."""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Literal
from typing import NewType
from typing import Self
from typing import TypeAlias

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MatplotlibPath
from numpy.typing import NDArray
from PIL import Image

from benchmarks.errors import AmbiguousFigurePanelRegistrationError
from benchmarks.errors import EmptyFigureColourSamplesError
from benchmarks.errors import InvalidFigureInwardOffsetError
from benchmarks.errors import InvalidReferenceOverlayError
from benchmarks.errors import MissingFigureAxesError
from benchmarks.held_reference_cases import Figure7Observation
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_all_held_reference_cases
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourcePrimitive
from benchmarks.held_reference_geometry import SourceToWorld
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.straight_skeleton_2 import offset_polygon
from tools.held_reference_extractor import ExtractedCase
from tools.held_reference_extractor import extract_reference_sources

RasterPixel = NewType("RasterPixel", int)
DotsPerInch = NewType("DotsPerInch", int)
MillimetresPerRasterPixel = NewType("MillimetresPerRasterPixel", float)
OverlayRole: TypeAlias = Literal["source", "analytic", "junction", "projection", "start", "tool"]
Figure7PanelName: TypeAlias = Literal["a", "b", "c"]

OVERLAY_LONG_SIDE = RasterPixel(2400)
OVERLAY_RENDER_DPI = DotsPerInch(300)
FIGURE7_POPPLER_DPI = DotsPerInch(600)
FIGURE7_PDF_PAGE = 14
FIGURE7_PAGE_WIDTH = RasterPixel(4500)
FIGURE7_PAGE_HEIGHT = RasterPixel(6180)
MINIMUM_SHORT_SIDE = RasterPixel(900)
ARC_RENDER_SEGMENTS = 256
# Black/grey axes are well below this RGB8 ceiling on the named Poppler render.
AXIS_DARK_CHANNEL_CEILING = 96
# 25 pixels is three PDF points at 600 DPI and contains the complete tick stroke.
AXIS_BORDER_BAND_PX = 25
# Published ticks extend by at least one PDF point (8.33 pixels at 600 DPI).
AXIS_TICK_MINIMUM_DARK_PIXELS = 8
# Panel c's approved half-open crop includes one white raster column outside
# each vertical axis; the border itself must still span every interior column.
AXIS_CROP_EXTERIOR_COLUMNS = 2
EXPECTED_X_TICKS = 9
EXPECTED_Y_TICKS = 7
# The analyst-stable mask excludes near-white anti-aliasing and low-chroma grey axes.
COLOUR_SAMPLE_MINIMUM_CHANNEL = 64
COLOUR_SAMPLE_MINIMUM_CHROMA = 48
PANEL_INTERIOR_BORDER_PX = 5
DEFAULT_OUTPUT_DIRECTORY = Path("docs/assets/images")
_EXTRACTED_TO_CANONICAL_NAME = {
    "figure-5": "figure5",
    "figure-8-upper": "figure8_upper",
    "figure-8-skis": "figure8_crossed_skis",
    "figure-8-monstera": "figure8_monstera",
}


@dataclass(frozen=True)
class RasterPoint2:
    """A pixel location in the full 600-DPI Poppler page frame."""

    x: RasterPixel
    y: RasterPixel

    @classmethod
    def build(cls, x: int, y: int) -> Self:
        if x < 0 or y < 0:
            raise AmbiguousFigurePanelRegistrationError("Raster coordinates must be non-negative pixel indices.")
        return cls(RasterPixel(x), RasterPixel(y))


@dataclass(frozen=True)
class RasterCrop:
    """A half-open RGB page crop in raster-pixel coordinates."""

    minimum: RasterPoint2
    maximum: RasterPoint2

    @classmethod
    def build(cls, minimum_x: int, minimum_y: int, maximum_x: int, maximum_y: int) -> Self:
        minimum = RasterPoint2.build(minimum_x, minimum_y)
        maximum = RasterPoint2.build(maximum_x, maximum_y)
        if minimum_x >= maximum_x or minimum_y >= maximum_y:
            raise AmbiguousFigurePanelRegistrationError("A raster crop requires strictly increasing half-open bounds.")
        return cls(minimum, maximum)


@dataclass(frozen=True)
class RasterBounds:
    """Inclusive extrema of detected raster support."""

    minimum: RasterPoint2
    maximum: RasterPoint2

    @classmethod
    def build(cls, minimum: RasterPoint2, maximum: RasterPoint2) -> Self:
        if int(minimum.x) >= int(maximum.x) or int(minimum.y) >= int(maximum.y):
            raise AmbiguousFigurePanelRegistrationError("Colour support requires two-dimensional raster extent.")
        return cls(minimum, maximum)


@dataclass(frozen=True)
class WorldBounds:
    """Inclusive normalized-world support extrema."""

    minimum: Point2[WorldXY]
    maximum: Point2[WorldXY]

    @classmethod
    def build(cls, minimum: Point2[WorldXY], maximum: Point2[WorldXY]) -> Self:
        values = (float(minimum.x), float(minimum.y), float(maximum.x), float(maximum.y))
        if not all(math.isfinite(value) for value in values) or values[0] >= values[2] or values[1] >= values[3]:
            raise InvalidFigureInwardOffsetError("An inward offset requires finite, strictly increasing world support bounds.")
        return cls(minimum, maximum)


@dataclass(frozen=True)
class AxisAlignedDisplayAffine:
    """An anisotropic raster-to-world display map with mandatory Y reflection."""

    raster_minimum: RasterPoint2
    raster_maximum: RasterPoint2
    world_bounds: WorldBounds
    x_scale: MillimetresPerRasterPixel
    y_scale: MillimetresPerRasterPixel

    @classmethod
    def build(
        cls,
        *,
        raster_minimum: RasterPoint2,
        raster_maximum: RasterPoint2,
        world_bounds: WorldBounds,
    ) -> Self:
        x_pixels = int(raster_maximum.x) - int(raster_minimum.x)
        y_pixels = int(raster_maximum.y) - int(raster_minimum.y)
        if x_pixels <= 0 or y_pixels <= 0:
            raise AmbiguousFigurePanelRegistrationError("Panel colour support cannot determine two display scales.")
        x_scale = (float(world_bounds.maximum.x) - float(world_bounds.minimum.x)) / x_pixels
        y_scale = (float(world_bounds.maximum.y) - float(world_bounds.minimum.y)) / y_pixels
        if not math.isfinite(x_scale) or not math.isfinite(y_scale) or x_scale <= 0.0 or y_scale <= 0.0:
            raise AmbiguousFigurePanelRegistrationError("Panel support produced an invalid axis-aligned display affine.")
        return cls(
            raster_minimum,
            raster_maximum,
            world_bounds,
            MillimetresPerRasterPixel(x_scale),
            MillimetresPerRasterPixel(y_scale),
        )

    def point(self, raster: RasterPoint2) -> Point2[WorldXY]:
        return Point2[WorldXY].build(
            float(self.world_bounds.minimum.x) + (int(raster.x) - int(self.raster_minimum.x)) * float(self.x_scale),
            float(self.world_bounds.maximum.y) - (int(raster.y) - int(self.raster_minimum.y)) * float(self.y_scale),
        )

    def raster_point(self, world: Point2[WorldXY]) -> tuple[float, float]:
        return (
            int(self.raster_minimum.x) + (float(world.x) - float(self.world_bounds.minimum.x)) / float(self.x_scale),
            int(self.raster_minimum.y) + (float(self.world_bounds.maximum.y) - float(world.y)) / float(self.y_scale),
        )


@dataclass(frozen=True)
class Figure7PanelGrammar:
    """Approved page-14 plot and colourbar regions for one named panel."""

    panel: Figure7PanelName
    plot_crop: RasterCrop
    colourbar_crop: RasterCrop

    @classmethod
    def build(cls, *, panel: Figure7PanelName, plot_crop: RasterCrop, colourbar_crop: RasterCrop) -> Self:
        return cls(panel, plot_crop, colourbar_crop)


FIGURE7_PANEL_GRAMMARS = (
    Figure7PanelGrammar.build(
        panel="a",
        plot_crop=RasterCrop.build(442, 667, 2009, 1759),
        colourbar_crop=RasterCrop.build(2044, 667, 2126, 1759),
    ),
    Figure7PanelGrammar.build(
        panel="b",
        plot_crop=RasterCrop.build(2387, 667, 3954, 1759),
        colourbar_crop=RasterCrop.build(3990, 667, 4071, 1759),
    ),
    Figure7PanelGrammar.build(
        panel="c",
        plot_crop=RasterCrop.build(1414, 2014, 2982, 3106),
        colourbar_crop=RasterCrop.build(3017, 2014, 3099, 3106),
    ),
)


@dataclass(frozen=True)
class Figure7PanelEvidence:
    """Detected colour support and its shape-only raster-to-world registration."""

    panel: Figure7PanelName
    plot_crop: RasterCrop
    rgb: NDArray[np.uint8]
    colour_mask: NDArray[np.bool_]
    colour_support: RasterBounds
    registration: AxisAlignedDisplayAffine


@dataclass(frozen=True)
class InwardOffsetEvidence:
    """A projection-derived comparator plus analytic registration support."""

    components: tuple[tuple[Point2[WorldXY], ...], ...]
    bounds: WorldBounds
    component_provenance: Literal["certified_polygon_projection"] = field(default="certified_polygon_projection", init=False)
    support_provenance: Literal["analytic_boundary_extrema"] = field(default="analytic_boundary_extrema", init=False)

    @classmethod
    def build(cls, components: Sequence[Sequence[Point2[WorldXY]]], bounds: WorldBounds) -> Self:
        ordered = tuple(tuple(component) for component in components)
        if len(ordered) != 1 or len(ordered[0]) < 3:
            raise InvalidFigureInwardOffsetError("Figure 7 requires one unambiguous inward-offset component.")
        return cls(ordered, bounds)


@dataclass(frozen=True)
class SourceCrop:
    """Publisher centreline plus its explicit PDF-to-world transform."""

    name: str
    primitives: tuple[SourcePrimitive, ...]
    transform: SourceToWorld

    @classmethod
    def build(
        cls,
        *,
        name: str,
        primitives: Sequence[SourcePrimitive],
        transform: SourceToWorld,
    ) -> Self:
        ordered = tuple(primitives)
        if not name or not ordered:
            raise InvalidReferenceOverlayError("A source crop requires a named, non-empty publisher centreline.")
        return cls(name, ordered, transform)

    @classmethod
    def normalized(
        cls,
        *,
        name: str,
        primitives: Sequence[SourcePrimitive],
        source_tool_radius: PdfPointUnit,
    ) -> Self:
        radius = float(source_tool_radius)
        if not math.isfinite(radius) or radius <= 0.0:
            raise InvalidReferenceOverlayError("Source normalization requires a finite positive depicted tool radius.")
        ordered = tuple(primitives)
        minimum, maximum = _source_bounds(ordered)
        transform = SourceToWorld.build(
            source_origin=PdfPoint2.build(minimum.x, maximum.y),
            world_origin=Point2[WorldXY].build(0.0, 0.0),
            scale=MillimetresPerPdfPoint(1.0 / radius),
            reflect_source_y=True,
        )
        return cls.build(name=name, primitives=ordered, transform=transform)


@dataclass(frozen=True)
class OverlayCircle:
    centre: Point2[WorldXY]
    radius: Millimetre

    @classmethod
    def build(cls, centre: Point2[WorldXY], radius: Millimetre) -> Self:
        numeric_radius = float(radius)
        if not math.isfinite(numeric_radius) or numeric_radius <= 0.0:
            raise InvalidReferenceOverlayError("An overlay circle requires a finite positive radius.")
        return cls(centre, radius)


@dataclass(frozen=True)
class SourceOverlayMark:
    source: SourceCrop
    role: Literal["source"] = field(default="source", init=False)


@dataclass(frozen=True)
class AnalyticOverlayMark:
    primitives: tuple[ReferencePrimitive, ...]
    role: Literal["analytic"] = field(default="analytic", init=False)


@dataclass(frozen=True)
class JunctionOverlayMark:
    points: tuple[Point2[WorldXY], ...]
    role: Literal["junction"] = field(default="junction", init=False)


@dataclass(frozen=True)
class ProjectionOverlayMark:
    points: tuple[Point2[WorldXY], ...]
    role: Literal["projection"] = field(default="projection", init=False)


@dataclass(frozen=True)
class StartOverlayMark:
    circle: OverlayCircle
    role: Literal["start"] = field(default="start", init=False)


@dataclass(frozen=True)
class ToolOverlayMark:
    circle: OverlayCircle
    role: Literal["tool"] = field(default="tool", init=False)


OverlayMark: TypeAlias = SourceOverlayMark | AnalyticOverlayMark | JunctionOverlayMark | ProjectionOverlayMark | StartOverlayMark | ToolOverlayMark


def reference_overlay_marks(case: HeldReferenceCase, source: SourceCrop) -> tuple[OverlayMark, ...]:
    """Expose every semantic layer before raster rendering."""
    if source.name != case.name:
        raise InvalidReferenceOverlayError(f"Source crop {source.name!r} does not belong to case {case.name!r}.")
    if case.start_marker is None or case.start_marker_radius is None:
        raise InvalidReferenceOverlayError(f"Case {case.name!r} has no published start-marker observation.")
    tool_centre = source.transform.point(case.source_tool_centre)
    return (
        SourceOverlayMark(source),
        AnalyticOverlayMark(case.boundary.primitives),
        JunctionOverlayMark(tuple(primitive.start for primitive in case.boundary.primitives)),
        ProjectionOverlayMark(case.projection.points),
        StartOverlayMark(OverlayCircle.build(case.start_marker, case.start_marker_radius)),
        ToolOverlayMark(OverlayCircle.build(tool_centre, Millimetre(float(case.tool_radius.value)))),
    )


def reference_overlay_caption(case: HeldReferenceCase) -> str:
    """Return the evidence-bearing caption rendered below one overlay."""
    figure = f"Figure {case.figure}"
    if case.subfigure is not None:
        figure = f"{figure} {case.subfigure.replace('_', ' ')}"
    return (
        f"{figure} | {case.name} | r = 1 mm | {len(case.boundary.primitives)} primitives | "
        f"{case.projection_vertex_count} projection vertices | "
        f"reconstruction upper bound = {float(case.reconstruction.deviation_upper_bound):.6g} mm | "
        f"observed projection deviation = {float(case.projection.observed_deviation):.6g} mm"
    )


def reference_overlay_caption_lines(case: HeldReferenceCase) -> tuple[str, str]:
    """Split a long overlay caption at its semantic evidence boundary."""
    caption = reference_overlay_caption(case)
    prefix, separator, evidence = caption.partition(" | reconstruction upper bound")
    if not separator:
        raise InvalidReferenceOverlayError(f"Case {case.name!r} produced an incomplete overlay caption.")
    return prefix, f"reconstruction upper bound{evidence}"


def render_reference_overlay(
    case: HeldReferenceCase,
    source: SourceCrop,
    output: Path,
    *,
    figure7_panels: Sequence[Figure7PanelEvidence] = (),
) -> None:
    """Render one semantic reference overlay to an exact-size PNG."""
    if output.suffix.lower() != ".png":
        raise InvalidReferenceOverlayError("A reference overlay output must be a PNG path.")
    marks = reference_overlay_marks(case, source)
    world_points = _overlay_extent_points(marks)
    panels = tuple(figure7_panels)
    if panels:
        if case.name != "figure5" or tuple(panel.panel for panel in panels) != ("a", "b", "c"):
            raise AmbiguousFigurePanelRegistrationError("Only Figure 5 may embed the measured Figure 7 panels a, b, and c.")
        width = height = int(OVERLAY_LONG_SIDE)
    else:
        width, height = _raster_dimensions(world_points)
    figure = Figure(
        figsize=(width / int(OVERLAY_RENDER_DPI), height / int(OVERLAY_RENDER_DPI)),
        dpi=int(OVERLAY_RENDER_DPI),
        facecolor="white",
    )
    FigureCanvasAgg(figure)
    axes = figure.add_axes((0.055, 0.39, 0.89, 0.56) if panels else (0.055, 0.13, 0.89, 0.82))
    _draw_marks(axes, marks)
    _configure_axes(axes, world_points)
    if panels:
        _draw_figure7_panels(figure, panels, figure7_inward_offset(case))
    caption = reference_overlay_caption(case) if panels else "\n".join(reference_overlay_caption_lines(case))
    figure.text(0.5, 0.025 if panels else 0.045, caption, ha="center", va="center", fontsize=6.5, color="#202428")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, format="png", dpi=int(OVERLAY_RENDER_DPI), facecolor="white")


def render_reference_overlays(
    cases: Sequence[HeldReferenceCase],
    sources: Sequence[SourceCrop],
    output_directory: Path,
    *,
    figure7_panels: Sequence[Figure7PanelEvidence] = (),
) -> tuple[Path, ...]:
    """Render one deterministically named PNG for every case/source pair."""
    if len(cases) != len(sources):
        raise InvalidReferenceOverlayError("Reference cases and publisher source crops must have equal cardinality.")
    outputs: list[Path] = []
    for case, source in zip(cases, sources):
        output = output_directory / f"held_reference_{case.name}.png"
        embedded_panels = figure7_panels if case.name == "figure5" else ()
        render_reference_overlay(case, source, output, figure7_panels=embedded_panels)
        outputs.append(output)
    return tuple(outputs)


def source_crops_from_extracted(extracted_cases: Sequence[ExtractedCase]) -> tuple[SourceCrop, ...]:
    """Convert publisher extraction records into normalized overlay sources."""
    source_crops: list[SourceCrop] = []
    for extracted in extracted_cases:
        try:
            canonical_name = _EXTRACTED_TO_CANONICAL_NAME[extracted.name]
        except KeyError as error:
            raise InvalidReferenceOverlayError(f"Unknown extracted publisher figure {extracted.name!r}.") from error
        source_crops.append(
            SourceCrop.normalized(
                name=canonical_name,
                primitives=extracted.sources,
                source_tool_radius=extracted.tool_circle.radius,
            )
        )
    return tuple(source_crops)


def measure_figure7_panel(
    page_rgb: NDArray[np.uint8],
    grammar: Figure7PanelGrammar,
    inward_bounds: WorldBounds,
) -> Figure7PanelEvidence:
    """Detect axes and coloured support, then derive one shape-only display affine."""
    _validate_rgb_page(page_rgb)
    panel_rgb = _crop_rgb(page_rgb, grammar.plot_crop)
    colourbar_rgb = _crop_rgb(page_rgb, grammar.colourbar_crop)
    _validate_panel_axes(panel_rgb)
    colour_mask = _colour_mask(panel_rgb, exclude_border=True)
    colourbar_mask = _colour_mask(colourbar_rgb, exclude_border=False)
    if not np.any(colour_mask) or not np.any(colourbar_mask):
        raise EmptyFigureColourSamplesError(f"Figure 7 panel {grammar.panel!r} contains no stable RGB8 colour samples.")
    y_indexes, x_indexes = np.nonzero(colour_mask)
    support = RasterBounds.build(
        RasterPoint2.build(int(grammar.plot_crop.minimum.x) + int(x_indexes.min()), int(grammar.plot_crop.minimum.y) + int(y_indexes.min())),
        RasterPoint2.build(int(grammar.plot_crop.minimum.x) + int(x_indexes.max()), int(grammar.plot_crop.minimum.y) + int(y_indexes.max())),
    )
    registration = AxisAlignedDisplayAffine.build(
        raster_minimum=support.minimum,
        raster_maximum=support.maximum,
        world_bounds=inward_bounds,
    )
    return Figure7PanelEvidence(grammar.panel, grammar.plot_crop, panel_rgb.copy(), colour_mask, support, registration)


def figure7_inward_offset(case: HeldReferenceCase) -> InwardOffsetEvidence:
    """Offset the certified projection while retaining analytic support authority."""
    radius = float(case.tool_radius.value)
    boundary_bounds = _analytic_boundary_bounds(case.boundary.primitives)
    inward_bounds = WorldBounds.build(
        Point2[WorldXY].build(float(boundary_bounds.minimum.x) + radius, float(boundary_bounds.minimum.y) + radius),
        Point2[WorldXY].build(float(boundary_bounds.maximum.x) - radius, float(boundary_bounds.maximum.y) - radius),
    )
    polygon_points = [(float(point.x), float(point.y), 0.0) for point in case.projection.points]
    try:
        polygons = offset_polygon(polygon_points, radius)
    except (RuntimeError, ValueError) as error:
        raise InvalidFigureInwardOffsetError(f"Case {case.name!r} cannot construct its one-radius inward component.") from error
    components = tuple(tuple(Point2[WorldXY].build(float(point.x), float(point.y)) for point in polygon.points) for polygon in polygons if len(polygon.points) >= 3)
    try:
        return InwardOffsetEvidence.build(components, inward_bounds)
    except InvalidFigureInwardOffsetError as error:
        raise InvalidFigureInwardOffsetError(f"Case {case.name!r} has no unique valid one-radius inward component.") from error


def measure_figure7_panels(pdf_path: Path, case: HeldReferenceCase) -> tuple[Figure7PanelEvidence, ...]:
    """Render and validate all three Figure 7 shape-only panels at 600 DPI."""
    if case.name != "figure5" or case.figure7_observation is None:
        raise AmbiguousFigurePanelRegistrationError("Figure 7 evidence belongs only to the Figure 5 reference case.")
    page_rgb = _render_figure7_page(pdf_path)
    inward = figure7_inward_offset(case)
    panels = tuple(measure_figure7_panel(page_rgb, grammar, inward.bounds) for grammar in FIGURE7_PANEL_GRAMMARS)
    _validate_panel_translation_agreement(panels)
    return panels


def measure_figure7_observation(pdf_path: Path, case: HeldReferenceCase) -> Figure7Observation:
    """Validate the live panels and return their non-numeric corpus observation."""
    measure_figure7_panels(pdf_path, case)
    observation = case.figure7_observation
    if observation is None:
        raise AmbiguousFigurePanelRegistrationError("Figure 5 has no declared Figure 7 observation.")
    return observation


def _render_figure7_page(pdf_path: Path) -> NDArray[np.uint8]:
    source_pdf = pdf_path.resolve(strict=True)
    build_directory = Path("build")
    build_directory.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="held-figure7-", dir=build_directory) as temporary:
        output_prefix = Path(temporary) / "page-14"
        try:
            subprocess.run(
                [
                    "pdftoppm",
                    "-f",
                    str(FIGURE7_PDF_PAGE),
                    "-l",
                    str(FIGURE7_PDF_PAGE),
                    "-r",
                    str(int(FIGURE7_POPPLER_DPI)),
                    "-png",
                    "-singlefile",
                    str(source_pdf),
                    str(output_prefix),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError) as error:
            raise InvalidReferenceOverlayError("Poppler could not render publisher PDF page 14 at 600 DPI.") from error
        rendered_path = output_prefix.with_suffix(".png")
        with Image.open(rendered_path) as image:
            if image.mode != "RGB" or image.size != (int(FIGURE7_PAGE_WIDTH), int(FIGURE7_PAGE_HEIGHT)):
                raise AmbiguousFigurePanelRegistrationError("Figure 7 Poppler render must be exactly 4500x6180 RGB at 600 DPI.")
            return np.asarray(image, dtype=np.uint8).copy()


def _validate_panel_translation_agreement(panels: Sequence[Figure7PanelEvidence]) -> None:
    if tuple(panel.panel for panel in panels) != ("a", "b", "c"):
        raise AmbiguousFigurePanelRegistrationError("Figure 7 requires panels a, b, and c exactly once and in publication order.")
    reference = panels[0].colour_support
    for panel in panels[1:]:
        support = panel.colour_support
        minimum_translation = (
            int(support.minimum.x) - int(reference.minimum.x),
            int(support.minimum.y) - int(reference.minimum.y),
        )
        maximum_translation = (
            int(support.maximum.x) - int(reference.maximum.x),
            int(support.maximum.y) - int(reference.maximum.y),
        )
        if minimum_translation != maximum_translation:
            raise AmbiguousFigurePanelRegistrationError("Figure 7 panel colour supports disagree after page-frame translation.")


def _analytic_boundary_bounds(primitives: Sequence[ReferencePrimitive]) -> WorldBounds:
    points: list[Point2[WorldXY]] = []
    for primitive in primitives:
        points.extend((primitive.start, primitive.end))
        if isinstance(primitive, ReferenceArc):
            centre_x = float(primitive.centre.x)
            centre_y = float(primitive.centre.y)
            radius = math.hypot(float(primitive.start.x) - centre_x, float(primitive.start.y) - centre_y)
            start_angle = math.atan2(float(primitive.start.y) - centre_y, float(primitive.start.x) - centre_x)
            sweep = float(primitive.sweep)
            for quadrant, direction in enumerate(((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0))):
                angle = quadrant * math.pi / 2.0
                travel = (angle - start_angle) % math.tau if sweep > 0.0 else -((start_angle - angle) % math.tau)
                if (sweep > 0.0 and 0.0 <= travel <= sweep) or (sweep < 0.0 and sweep <= travel <= 0.0):
                    points.append(Point2[WorldXY].build(centre_x + direction[0] * radius, centre_y + direction[1] * radius))
    if not points:
        raise InvalidFigureInwardOffsetError("An analytic boundary is required to derive inward support bounds.")
    x_values = tuple(float(point.x) for point in points)
    y_values = tuple(float(point.y) for point in points)
    return WorldBounds.build(
        Point2[WorldXY].build(min(x_values), min(y_values)),
        Point2[WorldXY].build(max(x_values), max(y_values)),
    )


def _validate_rgb_page(page_rgb: NDArray[np.uint8]) -> None:
    if page_rgb.dtype != np.uint8 or page_rgb.ndim != 3 or page_rgb.shape[2] != 3:
        raise AmbiguousFigurePanelRegistrationError("Figure 7 measurement requires one RGB8 Poppler page.")


def _crop_rgb(page_rgb: NDArray[np.uint8], crop: RasterCrop) -> NDArray[np.uint8]:
    minimum_x = int(crop.minimum.x)
    minimum_y = int(crop.minimum.y)
    maximum_x = int(crop.maximum.x)
    maximum_y = int(crop.maximum.y)
    if maximum_x > page_rgb.shape[1] or maximum_y > page_rgb.shape[0]:
        raise AmbiguousFigurePanelRegistrationError("An approved Figure 7 crop lies outside the rendered page.")
    return page_rgb[minimum_y:maximum_y, minimum_x:maximum_x]


def _validate_panel_axes(panel_rgb: NDArray[np.uint8]) -> None:
    dark = np.max(panel_rgb, axis=2) <= AXIS_DARK_CHANNEL_CEILING
    if panel_rgb.shape[0] < 2 * AXIS_BORDER_BAND_PX or panel_rgb.shape[1] < 2 * AXIS_BORDER_BAND_PX:
        raise MissingFigureAxesError("A Figure 7 panel is too small to contain the published axes.")
    horizontal_border = max(int(dark[row].sum()) for row in (*range(4), *range(panel_rgb.shape[0] - 4, panel_rgb.shape[0])))
    vertical_border = max(int(dark[:, column].sum()) for column in (*range(4), *range(panel_rgb.shape[1] - 4, panel_rgb.shape[1])))
    if horizontal_border < panel_rgb.shape[1] - AXIS_CROP_EXTERIOR_COLUMNS or vertical_border < panel_rgb.shape[0] - AXIS_CROP_EXTERIOR_COLUMNS:
        raise MissingFigureAxesError("A Figure 7 panel lacks its complete rectangular plot axes.")
    top = _dark_run_clusters(np.sum(dark[:AXIS_BORDER_BAND_PX], axis=0) > AXIS_TICK_MINIMUM_DARK_PIXELS)
    bottom = _dark_run_clusters(np.sum(dark[-AXIS_BORDER_BAND_PX:], axis=0) > AXIS_TICK_MINIMUM_DARK_PIXELS)
    left = _dark_run_clusters(np.sum(dark[:, :AXIS_BORDER_BAND_PX], axis=1) > AXIS_TICK_MINIMUM_DARK_PIXELS)
    right = _dark_run_clusters(np.sum(dark[:, -AXIS_BORDER_BAND_PX:], axis=1) > AXIS_TICK_MINIMUM_DARK_PIXELS)
    if top != bottom or left != right or len(top) != EXPECTED_X_TICKS or len(left) != EXPECTED_Y_TICKS:
        raise AmbiguousFigurePanelRegistrationError("Figure 7 panel axes do not contain the unique published 9-by-7 tick lattice.")


def _dark_run_clusters(signal: NDArray[np.bool_]) -> tuple[tuple[int, int], ...]:
    indexes = tuple(int(index) for index in np.flatnonzero(signal))
    if not indexes:
        return ()
    clusters: list[tuple[int, int]] = []
    start = indexes[0]
    previous = start
    for index in indexes[1:]:
        if index != previous + 1:
            clusters.append((start, previous))
            start = index
        previous = index
    clusters.append((start, previous))
    return tuple(clusters)


def _colour_mask(rgb: NDArray[np.uint8], *, exclude_border: bool) -> NDArray[np.bool_]:
    channels = rgb.astype(np.int16)
    maximum = np.max(channels, axis=2)
    chroma = maximum - np.min(channels, axis=2)
    mask = (maximum >= COLOUR_SAMPLE_MINIMUM_CHANNEL) & (chroma >= COLOUR_SAMPLE_MINIMUM_CHROMA)
    if exclude_border:
        mask[:PANEL_INTERIOR_BORDER_PX] = False
        mask[-PANEL_INTERIOR_BORDER_PX:] = False
        mask[:, :PANEL_INTERIOR_BORDER_PX] = False
        mask[:, -PANEL_INTERIOR_BORDER_PX:] = False
    return mask


def _draw_marks(axes: Axes, marks: Sequence[OverlayMark]) -> None:
    for mark in marks:
        if isinstance(mark, SourceOverlayMark):
            _draw_source(axes, mark.source)
        elif isinstance(mark, AnalyticOverlayMark):
            _draw_analytic(axes, mark.primitives)
        elif isinstance(mark, JunctionOverlayMark):
            axes.scatter(
                [float(point.x) for point in mark.points],
                [float(point.y) for point in mark.points],
                s=7.0,
                c="#111111",
                zorder=5,
            )
        elif isinstance(mark, ProjectionOverlayMark):
            points = (*mark.points, mark.points[0])
            axes.plot(
                [float(point.x) for point in points],
                [float(point.y) for point in points],
                color="#2767b1",
                linewidth=1.15,
                linestyle=(0, (4, 3)),
                zorder=4.5,
            )
        elif isinstance(mark, StartOverlayMark):
            _draw_circle(axes, mark.circle, edge="#c62828", linewidth=1.0, linestyle="--", zorder=6)
        else:
            _draw_circle(axes, mark.circle, edge="#d7191c", linewidth=1.25, linestyle="-", zorder=6)


def _draw_figure7_panels(
    figure: Figure,
    panels: Sequence[Figure7PanelEvidence],
    inward: InwardOffsetEvidence,
) -> None:
    figure.text(
        0.5,
        0.345,
        "Figure 7 tool-centre distributions - projection-derived comparator registered to analytic support\nshape-only falsification evidence; no numeric fidelity gate",
        ha="center",
        va="center",
        fontsize=6.5,
        color="#202428",
    )
    for index, panel in enumerate(panels):
        axes = figure.add_axes((0.035 + index * 0.325, 0.075, 0.29, 0.245))
        axes.imshow(panel.rgb, origin="upper")
        for component in inward.components:
            closed = (*component, component[0])
            raster_points = tuple(panel.registration.raster_point(point) for point in closed)
            x_values = [point[0] - int(panel.plot_crop.minimum.x) for point in raster_points]
            y_values = [point[1] - int(panel.plot_crop.minimum.y) for point in raster_points]
            axes.plot(x_values, y_values, color="white", linewidth=1.8, zorder=3)
            axes.plot(x_values, y_values, color="#159447", linewidth=0.85, zorder=4)
        axes.set_title(f"({panel.panel}) projection-derived one-radius comparator", fontsize=6.0, pad=2.0)
        axes.set_xlim(0.0, float(panel.rgb.shape[1]))
        axes.set_ylim(float(panel.rgb.shape[0]), 0.0)
        axes.set_axis_off()


def _draw_source(axes: Axes, source: SourceCrop) -> None:
    for primitive in source.primitives:
        if isinstance(primitive, SourceLine):
            start = source.transform.point(primitive.start)
            end = source.transform.point(primitive.end)
            axes.plot(
                (float(start.x), float(end.x)),
                (float(start.y), float(end.y)),
                color="#9a9a96",
                linewidth=2.0,
                zorder=1,
            )
            continue
        points = tuple(source.transform.point(point) for point in (primitive.start, primitive.control1, primitive.control2, primitive.end))
        path = MatplotlibPath(
            [(float(point.x), float(point.y)) for point in points],
            (MatplotlibPath.MOVETO, MatplotlibPath.CURVE4, MatplotlibPath.CURVE4, MatplotlibPath.CURVE4),
        )
        axes.add_patch(PathPatch(path, facecolor="none", edgecolor="#9a9a96", linewidth=2.0, zorder=1))


def _draw_analytic(axes: Axes, primitives: Sequence[ReferencePrimitive]) -> None:
    for primitive in primitives:
        points = _reference_primitive_points(primitive)
        axes.plot(
            [float(point.x) for point in points],
            [float(point.y) for point in points],
            color="#159447",
            linewidth=1.15,
            zorder=4,
        )


def _draw_circle(
    axes: Axes,
    circle: OverlayCircle,
    *,
    edge: str,
    linewidth: float,
    linestyle: str,
    zorder: int,
) -> None:
    axes.add_patch(
        Circle(
            (float(circle.centre.x), float(circle.centre.y)),
            float(circle.radius),
            facecolor="none",
            edgecolor=edge,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=zorder,
        )
    )


def _reference_primitive_points(primitive: ReferencePrimitive) -> tuple[Point2[WorldXY], ...]:
    if isinstance(primitive, ReferenceLine):
        return primitive.start, primitive.end
    radius_x = float(primitive.start.x) - float(primitive.centre.x)
    radius_y = float(primitive.start.y) - float(primitive.centre.y)
    points: list[Point2[WorldXY]] = []
    for index in range(ARC_RENDER_SEGMENTS + 1):
        angle = float(primitive.sweep) * index / ARC_RENDER_SEGMENTS
        cosine = math.cos(angle)
        sine = math.sin(angle)
        points.append(
            Point2[WorldXY].build(
                float(primitive.centre.x) + cosine * radius_x - sine * radius_y,
                float(primitive.centre.y) + sine * radius_x + cosine * radius_y,
            )
        )
    return tuple(points)


def _source_bounds(primitives: Sequence[SourcePrimitive]) -> tuple[PdfPoint2, PdfPoint2]:
    if not primitives:
        raise InvalidReferenceOverlayError("Source normalization requires non-empty publisher geometry.")
    x_values: list[float] = []
    y_values: list[float] = []
    for primitive in primitives:
        x_values.extend((float(primitive.start.x), float(primitive.end.x)))
        y_values.extend((float(primitive.start.y), float(primitive.end.y)))
        if isinstance(primitive, SourceCubic):
            x_values.extend(_source_coordinate(primitive, "x", parameter) for parameter in _source_derivative_roots(primitive, "x"))
            y_values.extend(_source_coordinate(primitive, "y", parameter) for parameter in _source_derivative_roots(primitive, "y"))
    return PdfPoint2.build(min(x_values), min(y_values)), PdfPoint2.build(max(x_values), max(y_values))


def _source_coordinate(source: SourceCubic, axis: Literal["x", "y"], parameter: float) -> float:
    values = tuple(float(getattr(point, axis)) for point in (source.start, source.control1, source.control2, source.end))
    complement = 1.0 - parameter
    return complement**3 * values[0] + 3.0 * complement**2 * parameter * values[1] + 3.0 * complement * parameter**2 * values[2] + parameter**3 * values[3]


def _source_derivative_roots(source: SourceCubic, axis: Literal["x", "y"]) -> tuple[float, ...]:
    p0, p1, p2, p3 = (float(getattr(point, axis)) for point in (source.start, source.control1, source.control2, source.end))
    quadratic = -p0 + 3.0 * p1 - 3.0 * p2 + p3
    linear = 2.0 * (p0 - 2.0 * p1 + p2)
    constant = p1 - p0
    if quadratic == 0.0:
        if linear == 0.0:
            return ()
        root = -constant / linear
        return (root,) if 0.0 < root < 1.0 else ()
    discriminant = linear * linear - 4.0 * quadratic * constant
    if discriminant <= 0.0:
        return ()
    square_root = math.sqrt(discriminant)
    numerator = -0.5 * (linear + math.copysign(square_root, linear))
    roots = (numerator / quadratic, constant / numerator) if numerator != 0.0 else (-linear / (2.0 * quadratic),)
    return tuple(sorted({root for root in roots if 0.0 < root < 1.0}))


def _overlay_extent_points(marks: Sequence[OverlayMark]) -> tuple[Point2[WorldXY], ...]:
    points: list[Point2[WorldXY]] = []
    for mark in marks:
        if isinstance(mark, SourceOverlayMark):
            for source_primitive in mark.source.primitives:
                source_points: Sequence[PdfPoint2] = (source_primitive.start, source_primitive.end)
                if isinstance(source_primitive, SourceCubic):
                    source_points = (
                        source_primitive.start,
                        source_primitive.control1,
                        source_primitive.control2,
                        source_primitive.end,
                    )
                points.extend(mark.source.transform.point(point) for point in source_points)
        elif isinstance(mark, AnalyticOverlayMark):
            for reference_primitive in mark.primitives:
                points.extend(_reference_primitive_points(reference_primitive))
        elif isinstance(mark, (JunctionOverlayMark, ProjectionOverlayMark)):
            points.extend(mark.points)
        else:
            radius = float(mark.circle.radius)
            points.extend(
                (
                    Point2[WorldXY].build(float(mark.circle.centre.x) - radius, float(mark.circle.centre.y) - radius),
                    Point2[WorldXY].build(float(mark.circle.centre.x) + radius, float(mark.circle.centre.y) + radius),
                )
            )
    if not points:
        raise InvalidReferenceOverlayError("A reference overlay contains no drawable world geometry.")
    return tuple(points)


def _raster_dimensions(points: Sequence[Point2[WorldXY]]) -> tuple[int, int]:
    x_values = tuple(float(point.x) for point in points)
    y_values = tuple(float(point.y) for point in points)
    x_span = max(x_values) - min(x_values)
    y_span = max(y_values) - min(y_values)
    if not math.isfinite(x_span) or not math.isfinite(y_span) or x_span <= 0.0 or y_span <= 0.0:
        raise InvalidReferenceOverlayError("A reference overlay requires a finite two-dimensional extent.")
    if x_span >= y_span:
        return int(OVERLAY_LONG_SIDE), max(int(MINIMUM_SHORT_SIDE), round(int(OVERLAY_LONG_SIDE) * y_span / x_span))
    return max(int(MINIMUM_SHORT_SIDE), round(int(OVERLAY_LONG_SIDE) * x_span / y_span)), int(OVERLAY_LONG_SIDE)


def _configure_axes(axes: Axes, points: Sequence[Point2[WorldXY]]) -> None:
    x_values = tuple(float(point.x) for point in points)
    y_values = tuple(float(point.y) for point in points)
    x_span = max(x_values) - min(x_values)
    y_span = max(y_values) - min(y_values)
    axes.set_xlim(min(x_values) - 0.04 * x_span, max(x_values) + 0.04 * x_span)
    axes.set_ylim(min(y_values) - 0.04 * y_span, max(y_values) + 0.04 * y_span)
    axes.set_aspect("equal", adjustable="box")
    axes.set_axis_off()


def _parse_arguments(arguments: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("publisher_pdf", type=Path)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    return parser.parse_args(arguments)


def main(arguments: Sequence[str] | None = None) -> int:
    parsed = _parse_arguments(arguments if arguments is not None else sys.argv[1:])
    cases = load_all_held_reference_cases()
    sources = source_crops_from_extracted(extract_reference_sources(parsed.publisher_pdf))
    figure7_panels = measure_figure7_panels(parsed.publisher_pdf, cases[0])
    render_reference_overlays(cases, sources, parsed.output_directory, figure7_panels=figure7_panels)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
