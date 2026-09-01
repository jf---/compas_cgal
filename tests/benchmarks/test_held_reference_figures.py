from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest
import numpy as np
from PIL import Image
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

import benchmarks.held_reference_figures as held_figures
from benchmarks.errors import AmbiguousFigurePanelRegistrationError
from benchmarks.errors import EmptyFigureColourSamplesError
from benchmarks.errors import InvalidFigureInwardOffsetError
from benchmarks.errors import InvalidReferenceOverlayError
from benchmarks.errors import MissingFigureAxesError
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import AxisAlignedDisplayAffine
from benchmarks.held_reference_figures import Figure7PanelGrammar
from benchmarks.held_reference_figures import InwardOffsetEvidence
from benchmarks.held_reference_figures import RasterCrop
from benchmarks.held_reference_figures import RasterPoint2
from benchmarks.held_reference_figures import SourceCrop
from benchmarks.held_reference_figures import WorldBounds
from benchmarks.held_reference_figures import measure_figure7_panel
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_reference_figures import reference_overlay_caption
from benchmarks.held_reference_figures import reference_overlay_caption_lines
from benchmarks.held_reference_figures import reference_overlay_marks
from benchmarks.held_reference_figures import render_reference_overlay
from benchmarks.held_reference_figures import render_reference_overlays
from benchmarks.held_reference_figures import source_crops_from_extracted
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourceToWorld
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from tools.held_reference_extractor import ExtractedCase
from tools.held_reference_extractor import PublishedToolCircle


def _figure5_case() -> HeldReferenceCase:
    return load_held_reference_case("figure5")


def _figure5_source() -> SourceCrop:
    corners = (
        PdfPoint2.build(0.0, 0.0),
        PdfPoint2.build(2.0, 0.0),
        PdfPoint2.build(2.0, 1.0),
        PdfPoint2.build(0.0, 1.0),
    )
    transform = SourceToWorld.build(
        source_origin=PdfPoint2.build(0.0, 1.0),
        world_origin=Point2[WorldXY].build(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0),
        reflect_source_y=True,
    )
    return SourceCrop.build(
        name="figure5",
        primitives=tuple(SourceLine.build(start, end) for start, end in zip(corners, (*corners[1:], corners[0]))),
        transform=transform,
    )


def test_overlay_contains_every_evidence_layer() -> None:
    marks = reference_overlay_marks(_figure5_case(), _figure5_source())

    assert {mark.role for mark in marks} == {
        "source",
        "analytic",
        "junction",
        "projection",
        "start",
        "tool",
    }


def test_projection_dashes_render_above_the_analytic_boundary() -> None:
    figure = Figure()
    FigureCanvasAgg(figure)
    axes = figure.add_subplot()

    held_figures._draw_marks(axes, reference_overlay_marks(_figure5_case(), _figure5_source()))

    analytic = next(line for line in axes.lines if line.get_color() == "#159447")
    projection = next(line for line in axes.lines if line.get_color() == "#2767b1")
    assert projection.get_zorder() > analytic.get_zorder()
    assert projection.get_linestyle() != "-"
    assert projection.get_linewidth() >= analytic.get_linewidth()


def test_overlay_caption_reports_the_reconstruction_evidence() -> None:
    case = _figure5_case()

    assert reference_overlay_caption(case) == (
        "Figure 5 | figure5 | r = 1 mm | 31 primitives | 65 projection vertices | reconstruction upper bound = 0.0738004 mm | observed projection deviation = 0.0321165 mm"
    )


def test_long_overlay_caption_is_split_at_the_evidence_boundary() -> None:
    case = load_held_reference_case("figure8_upper")

    lines = reference_overlay_caption_lines(case)

    assert lines == (
        "Figure 8 upper | figure8_upper | r = 1 mm | 70 primitives | 144 projection vertices",
        "reconstruction upper bound = 0.0813111 mm | observed projection deviation = 0.0302284 mm",
    )


def test_rendered_overlay_has_exact_long_side_and_rgb_pixels(tmp_path: Path) -> None:
    output = tmp_path / "figure5.png"

    render_reference_overlay(_figure5_case(), _figure5_source(), output)

    with Image.open(output) as rendered:
        assert max(rendered.size) == 2400
        assert rendered.mode in {"RGB", "RGBA"}


def test_render_reference_overlays_names_each_case_output(tmp_path: Path) -> None:
    outputs = render_reference_overlays((_figure5_case(),), (_figure5_source(),), tmp_path)

    assert outputs == (tmp_path / "held_reference_figure5.png",)
    assert outputs[0].is_file()


def test_rendered_figure5_embeds_three_shape_only_figure7_panels(tmp_path: Path) -> None:
    panel = measure_figure7_panel(_synthetic_figure7_page(), _synthetic_panel_grammar(), _synthetic_inward_bounds())
    panels = (panel, replace(panel, panel="b"), replace(panel, panel="c"))
    output = tmp_path / "figure5-with-figure7.png"

    render_reference_overlay(_figure5_case(), _figure5_source(), output, figure7_panels=panels)

    with Image.open(output) as rendered:
        assert rendered.size == (2400, 2400)


def test_overlay_rejects_source_for_another_case() -> None:
    source = _figure5_source()
    mismatched = SourceCrop.build(name="figure8_upper", primitives=source.primitives, transform=source.transform)

    with pytest.raises(InvalidReferenceOverlayError):
        reference_overlay_marks(_figure5_case(), mismatched)


def test_source_crop_rejects_empty_publisher_centreline() -> None:
    source = _figure5_source()

    with pytest.raises(InvalidReferenceOverlayError):
        SourceCrop.build(name="figure5", primitives=(), transform=source.transform)


def test_normalized_source_crop_maps_reflected_lower_left_to_world_origin() -> None:
    source = _figure5_source()

    normalized = SourceCrop.normalized(
        name="figure5",
        primitives=source.primitives,
        source_tool_radius=PdfPointUnit(2.0),
    )

    assert normalized.transform.point(PdfPoint2.build(0.0, 1.0)) == Point2[WorldXY].build(0.0, 0.0)
    assert float(normalized.transform.scale) == 0.5


def test_extracted_figure_name_maps_to_canonical_case_name() -> None:
    source = _figure5_source()
    extracted = ExtractedCase.build(
        name="figure-5",
        page=12,
        sources=source.primitives,
        tool_circle=PublishedToolCircle.build(PdfPoint2.build(4.0, 4.0), PdfPointUnit(2.0)),
        boundary_markers=(),
        start_markers=(),
        boundary_stroke_width=PdfPointUnit(1.0),
    )

    mapped = source_crops_from_extracted((extracted,))

    assert tuple(crop.name for crop in mapped) == ("figure5",)


def _synthetic_figure7_page(*, x_tick_count: int = 9, colour_samples: bool = True, side_padding: int = 0) -> np.ndarray:
    page = np.full((120, 210, 3), 255, dtype=np.uint8)
    panel = page[10:110, 10:170]
    left = side_padding
    right = panel.shape[1] - 1 - side_padding
    panel[0:2, left : right + 1] = 0
    panel[-2:, left : right + 1] = 0
    panel[:, left : left + 2] = 0
    panel[:, right - 1 : right + 1] = 0
    for x in np.linspace(left, right, x_tick_count, dtype=int):
        panel[:12, x : x + 2] = 0
        panel[-12:, x : x + 2] = 0
    for y in np.linspace(0, panel.shape[0] - 1, 7, dtype=int):
        panel[y : y + 2, left : left + 12] = 0
        panel[y : y + 2, right - 11 : right + 1] = 0
    if colour_samples:
        panel[30:80, 35:145] = (128, 16, 224)
        page[10:110, 180:190] = (224, 32, 96)
    return page


def _synthetic_panel_grammar() -> Figure7PanelGrammar:
    return Figure7PanelGrammar.build(
        panel="a",
        plot_crop=RasterCrop.build(10, 10, 170, 110),
        colourbar_crop=RasterCrop.build(180, 10, 190, 110),
    )


def _synthetic_inward_bounds() -> WorldBounds:
    return WorldBounds.build(
        Point2[WorldXY].build(1.0, 2.0),
        Point2[WorldXY].build(65.0, 46.0),
    )


def test_axis_aligned_display_affine_preserves_anisotropic_axis_scales() -> None:
    transform = AxisAlignedDisplayAffine.build(
        raster_minimum=RasterPoint2.build(10, 20),
        raster_maximum=RasterPoint2.build(110, 70),
        world_bounds=_synthetic_inward_bounds(),
    )

    assert float(transform.x_scale) == 0.64
    assert float(transform.y_scale) == 0.88
    assert not hasattr(transform, "reflect_raster_y")
    assert transform.point(RasterPoint2.build(10, 20)) == Point2[WorldXY].build(1.0, 46.0)
    assert transform.point(RasterPoint2.build(110, 70)) == Point2[WorldXY].build(65.0, 2.0)


def test_figure7_panel_measurement_detects_axes_and_colour_support() -> None:
    evidence = measure_figure7_panel(_synthetic_figure7_page(), _synthetic_panel_grammar(), _synthetic_inward_bounds())

    assert evidence.panel == "a"
    assert evidence.colour_support.minimum == RasterPoint2.build(45, 40)
    assert evidence.colour_support.maximum == RasterPoint2.build(154, 89)


def test_figure7_axes_accept_the_approved_one_pixel_side_padding() -> None:
    evidence = measure_figure7_panel(
        _synthetic_figure7_page(side_padding=1),
        _synthetic_panel_grammar(),
        _synthetic_inward_bounds(),
    )

    assert evidence.panel == "a"


def test_figure7_panel_rejects_missing_axes() -> None:
    with pytest.raises(MissingFigureAxesError):
        measure_figure7_panel(np.full((120, 210, 3), 255, dtype=np.uint8), _synthetic_panel_grammar(), _synthetic_inward_bounds())


def test_figure7_panel_rejects_empty_colour_samples() -> None:
    with pytest.raises(EmptyFigureColourSamplesError):
        measure_figure7_panel(_synthetic_figure7_page(colour_samples=False), _synthetic_panel_grammar(), _synthetic_inward_bounds())


def test_figure7_panel_rejects_ambiguous_tick_lattice() -> None:
    with pytest.raises(AmbiguousFigurePanelRegistrationError):
        measure_figure7_panel(_synthetic_figure7_page(x_tick_count=8), _synthetic_panel_grammar(), _synthetic_inward_bounds())


def test_inward_support_requires_strictly_increasing_world_bounds() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)

    with pytest.raises(InvalidFigureInwardOffsetError):
        WorldBounds.build(point, point)


def test_figure5_inward_offset_has_one_component_and_analytic_support_bounds() -> None:
    inward = figure7_inward_offset(_figure5_case())

    assert len(inward.components) == 1
    assert inward.component_provenance == "certified_polygon_projection"
    assert inward.support_provenance == "analytic_boundary_extrema"
    assert inward.bounds == WorldBounds.build(
        Point2[WorldXY].build(1.0002419667658535, 1.0),
        Point2[WorldXY].build(65.7718836534693, 44.784909922871265),
    )


def test_figure7_comparator_moves_with_a_perturbed_polygon_projection() -> None:
    case = _figure5_case()
    original = figure7_inward_offset(case)
    shifted_projection = replace(
        case.projection,
        points=tuple(Point2[WorldXY].build(float(point.x) + 5.0, float(point.y) - 3.0) for point in case.projection.points),
    )

    perturbed = figure7_inward_offset(replace(case, projection=shifted_projection))

    assert perturbed.components != original.components
    assert perturbed.bounds == original.bounds


def test_figure7_analytic_support_moves_independently_of_the_comparator() -> None:
    case = _figure5_case()
    original = figure7_inward_offset(case)
    translated_primitives = []
    for primitive in case.boundary.primitives:
        translated = {
            "start": Point2[WorldXY].build(float(primitive.start.x) + 2.0, float(primitive.start.y) + 4.0),
            "end": Point2[WorldXY].build(float(primitive.end.x) + 2.0, float(primitive.end.y) + 4.0),
        }
        if hasattr(primitive, "centre"):
            translated["centre"] = Point2[WorldXY].build(float(primitive.centre.x) + 2.0, float(primitive.centre.y) + 4.0)
        translated_primitives.append(replace(primitive, **translated))
    translated_case = replace(case, boundary=replace(case.boundary, primitives=tuple(translated_primitives)))

    moved = figure7_inward_offset(translated_case)

    assert moved.components == original.components
    expected_minimum = Point2[WorldXY].build(float(original.bounds.minimum.x) + 2.0, float(original.bounds.minimum.y) + 4.0)
    expected_maximum = Point2[WorldXY].build(float(original.bounds.maximum.x) + 2.0, float(original.bounds.maximum.y) + 4.0)
    for actual, expected in (
        (float(moved.bounds.minimum.x), float(expected_minimum.x)),
        (float(moved.bounds.minimum.y), float(expected_minimum.y)),
        (float(moved.bounds.maximum.x), float(expected_maximum.x)),
        (float(moved.bounds.maximum.y), float(expected_maximum.y)),
    ):
        assert actual in {expected, math.nextafter(expected, -math.inf), math.nextafter(expected, math.inf)}


def test_figure7_panel_titles_distinguish_comparator_from_analytic_support() -> None:
    panel = measure_figure7_panel(_synthetic_figure7_page(), _synthetic_panel_grammar(), _synthetic_inward_bounds())
    panels = (panel, replace(panel, panel="b"), replace(panel, panel="c"))
    figure = Figure()
    FigureCanvasAgg(figure)

    held_figures._draw_figure7_panels(figure, panels, figure7_inward_offset(_figure5_case()))

    assert "analytic support" in figure.texts[0].get_text()
    assert all("projection-derived" in axes.get_title() for axes in figure.axes)


def test_inward_offset_rejects_missing_or_duplicate_components() -> None:
    bounds = _synthetic_inward_bounds()
    triangle = (
        Point2[WorldXY].build(1.0, 2.0),
        Point2[WorldXY].build(2.0, 2.0),
        Point2[WorldXY].build(1.0, 3.0),
    )

    with pytest.raises(InvalidFigureInwardOffsetError):
        InwardOffsetEvidence.build((), bounds)
    with pytest.raises(InvalidFigureInwardOffsetError):
        InwardOffsetEvidence.build((triangle, triangle), bounds)


def test_figure7_rejects_panel_supports_without_one_page_translation() -> None:
    panel = measure_figure7_panel(_synthetic_figure7_page(), _synthetic_panel_grammar(), _synthetic_inward_bounds())
    shifted_support = replace(
        panel.colour_support,
        minimum=RasterPoint2.build(int(panel.colour_support.minimum.x) + 3, int(panel.colour_support.minimum.y)),
        maximum=RasterPoint2.build(int(panel.colour_support.maximum.x) + 4, int(panel.colour_support.maximum.y)),
    )

    with pytest.raises(AmbiguousFigurePanelRegistrationError):
        held_figures._validate_panel_translation_agreement((panel, replace(panel, panel="b"), replace(panel, panel="c", colour_support=shifted_support)))
