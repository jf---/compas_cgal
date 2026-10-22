"""Emission-order sequence plot: colour order is monotone in emission index and a tiny path renders."""

from __future__ import annotations

import io
import math
from pathlib import Path

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.image import imread

from benchmarks.tools.held_sequence_plot import ColourMode
from benchmarks.tools.held_sequence_plot import CoverageVerdict
from benchmarks.tools.held_sequence_plot import EmptySequenceError
from benchmarks.tools.held_sequence_plot import InvalidCoverageVerdictError
from benchmarks.tools.held_sequence_plot import SequenceCircle
from benchmarks.tools.held_sequence_plot import PAGE_COLOR
from benchmarks.tools.held_sequence_plot import SequencePlotData
from benchmarks.tools.held_sequence_plot import draw_sweep
from benchmarks.tools.held_sequence_plot import order_colors
from benchmarks.tools.held_sequence_plot import order_positions
from benchmarks.tools.held_sequence_plot import render_sequence


def _relative_luminance(rgba: np.ndarray) -> np.ndarray:
    linear = np.where(rgba[:, :3] <= 0.04045, rgba[:, :3] / 12.92, ((rgba[:, :3] + 0.055) / 1.055) ** 2.4)
    return np.asarray(linear @ np.array([0.2126, 0.7152, 0.0722]))


def _fixture() -> SequencePlotData:
    # Two runs of a small fan sweeping left to right inside a 10 x 6 mm box.
    circles = tuple(SequenceCircle(2.0 + 0.6 * index, 3.0, 1.0 + 0.1 * (index % 4), 0.0, 0 if index < 5 else 1) for index in range(9))
    boundary = ((0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0))
    contacts = [(c.x + c.radius, c.y) for c in circles]
    connectors = tuple((contacts[index], contacts[index + 1]) for index in range(len(circles) - 1))
    return SequencePlotData.build(
        workload="fixture",
        title="fixture",
        boundary=boundary,
        circles=circles,
        connectors=connectors,
        seam_connector=(contacts[-1], contacts[0]),
        engagement_deg=(40.0, 55.0, math.nan, 70.0, 20.0, 60.0, 79.0, 30.0),
        cap_deg=80.0,
        tool_radius=0.5,
        coverage=CoverageVerdict.build(covered=False, witnesses=((9.5, 5.5),), piece_kinds=("annulus_outer",)),
    )


def test_order_colour_is_monotone_in_emission_index() -> None:
    positions = order_positions(_fixture().circles, ColourMode.ORDER)
    assert np.all(np.diff(positions) > 0)
    assert positions[0] == 0.0 and positions[-1] == 1.0
    luminance = _relative_luminance(order_colors(_fixture().circles, ColourMode.ORDER))
    assert np.all(np.diff(luminance) > 0)


def test_run_colour_restarts_per_run_and_is_monotone_within() -> None:
    circles = _fixture().circles
    positions = order_positions(circles, ColourMode.RUN)
    luminance = _relative_luminance(order_colors(circles, ColourMode.RUN))
    for block in (slice(0, 5), slice(5, 9)):
        assert positions[block][0] == 0.0 and positions[block][-1] == 1.0
        assert np.all(np.diff(luminance[block]) > 0)


def test_tiny_fixture_renders(tmp_path: Path) -> None:
    for mode in ColourMode:
        output = tmp_path / f"sequence_{mode.value}.png"
        render_sequence(_fixture(), output, mode)
        image = imread(str(output))
        assert image.shape[0] > 100 and image.shape[1] > 100


def test_json_round_trip_keeps_order_and_unresolved() -> None:
    data = _fixture()
    again = SequencePlotData.from_json(data.to_json())
    assert again.circles == data.circles
    assert again.seam_connector == data.seam_connector
    assert math.isnan(again.engagement_deg[2])
    assert again.engagement_deg[3] == data.engagement_deg[3]
    assert again.tool_radius == data.tool_radius
    assert again.coverage == data.coverage


def test_empty_sequence_is_rejected() -> None:
    with pytest.raises(EmptySequenceError):
        SequencePlotData.build(
            workload="x", title="x", boundary=(), circles=(), connectors=(), seam_connector=None, engagement_deg=(), cap_deg=80.0, tool_radius=1.0, coverage=None
        )


def test_coverage_verdict_rejects_contradicting_witnesses() -> None:
    with pytest.raises(InvalidCoverageVerdictError):
        CoverageVerdict.build(covered=True, witnesses=((1.0, 2.0),), piece_kinds=("disk",))
    with pytest.raises(InvalidCoverageVerdictError):
        CoverageVerdict.build(covered=False, witnesses=(), piece_kinds=())
    assert "CERTIFIED COVERED" in CoverageVerdict.build(covered=True, witnesses=(), piece_kinds=()).caption(1.0)


def test_witnesses_group_into_sites_ordered_by_y() -> None:
    verdict = CoverageVerdict.build(
        covered=False,
        witnesses=((34.02, 71.41), (34.03, 71.35), (45.57, 61.18), (45.43, 61.15)),
        piece_kinds=("annulus_outer",) * 4,
    )
    sites = verdict.sites(1.0)
    assert [count for _, _, count in sites] == [2, 2]
    assert sites[0][1] < sites[1][1]
    assert "2 site(s): A (45.50, 61.16) \u00d72, B (34.03, 71.38) \u00d72" in verdict.caption(1.0)


def _gap_page_pixels(tool_radius: float) -> int:
    """Exact page-colour pixels on the segment between two centre circles whose gap is 0.6 mm."""
    circles = (SequenceCircle(3.0, 3.0, 1.0, 0.0, 0), SequenceCircle(5.6, 3.0, 1.0, math.pi, 1))
    data = SequencePlotData.build(
        workload="gap",
        title="gap",
        boundary=((0.0, 0.0), (8.6, 0.0), (8.6, 6.0), (0.0, 6.0)),
        circles=circles,
        connectors=(),
        seam_connector=None,
        engagement_deg=(10.0,),
        cap_deg=80.0,
        tool_radius=tool_radius,
        coverage=None,
    )
    figure = Figure(figsize=(6, 4.2), dpi=180)
    FigureCanvasAgg(figure)
    axis = figure.add_subplot()
    draw_sweep(axis, data, np.array([to_rgba("#3b0f70"), to_rgba("#fe9f6d")]), legend=False)
    axis.grid(False)
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=figure.dpi)
    buffer.seek(0)
    image = np.round(imread(buffer)[:, :, :3] * 255).astype(int)
    page = np.round(np.array(to_rgba(PAGE_COLOR)[:3]) * 255).astype(int)
    # Sample the gap between the centre circles' facing points (x in [4.0, 4.6], y = 3.2).
    xs = np.linspace(4.02, 4.58, 60)
    pixels = axis.transData.transform(np.column_stack((xs, np.full(len(xs), 3.2))))
    rows = image.shape[0] - 1 - np.round(pixels[:, 1]).astype(int)
    cols = np.round(pixels[:, 0]).astype(int)
    return int(np.sum(np.all(image[rows, cols] == page, axis=1)))


def test_gap_narrower_than_tool_radius_renders_no_white() -> None:
    # The rim lies inside its band in a shade of the band colour, so a page
    # pixel is exactly uncovered material, rim included.
    assert _gap_page_pixels(tool_radius=0.5) == 0
    # The same probe sees the gap when the bands are too thin to close it.
    assert _gap_page_pixels(tool_radius=0.05) > 30


# Page-colour match in 8-bit channels: the old near-white rim (#fdfbf4) is 3
# away from the page (#faf9f5) and must count as page; any band colour of the
# ramp is far more than this away.
PAGE_MATCH_LEVELS = 6
# Pixels around the start marker (a white disk drawn on top) are not swept-band paint.
START_MARKER_EXCLUSION_PX = 16


def _render_sweep(data: SequencePlotData) -> tuple[np.ndarray, Figure, object]:
    figure = Figure(figsize=(4.4, 4.4), dpi=180)
    canvas = FigureCanvasAgg(figure)
    axis = figure.add_subplot()
    draw_sweep(axis, data, order_colors(data.circles, ColourMode.ORDER), legend=False)
    axis.grid(False)
    canvas.draw()
    image = np.asarray(canvas.buffer_rgba())[:, :, :3].astype(int)
    return image, figure, axis


def _page_mask(image: np.ndarray) -> np.ndarray:
    page = np.round(np.array(to_rgba(PAGE_COLOR)[:3]) * 255).astype(int)
    return np.all(np.abs(image - page) <= PAGE_MATCH_LEVELS, axis=2)


def _pocket_mask(image: np.ndarray, axis: object, boundary: tuple[tuple[float, float], ...]) -> np.ndarray:
    """Pixels whose centre lies in the pocket, eroded by one pixel (4-neighbourhood)."""
    from matplotlib.path import Path as MplPath

    height, width = image.shape[:2]
    outline = MplPath(axis.transData.transform(np.array(boundary)))  # type: ignore[attr-defined]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    centres = np.column_stack(((cols + 0.5).ravel(), (height - rows - 0.5).ravel()))
    inside = outline.contains_points(centres).reshape(height, width)
    eroded = inside.copy()
    eroded[1:, :] &= inside[:-1, :]
    eroded[:-1, :] &= inside[1:, :]
    eroded[:, 1:] &= inside[:, :-1]
    eroded[:, :-1] &= inside[:, 1:]
    return eroded


def _pixel(axis: object, x: float, y: float, height: int) -> tuple[int, int]:
    px, py = axis.transData.transform((x, y))  # type: ignore[attr-defined]
    return int(height - py), int(px)


def _covered_corner_fan() -> SequencePlotData:
    # Square pocket [0, 4]^2, r = 0.5. A grid of circles R = 0.3 (each sweeps the
    # full disk of radius 0.8) at spacing 0.8 covers the square exactly: every
    # point is within 0.4 * sqrt(2) < 0.8 of a centre. Then a tight fan of 60
    # small circles turns in the convex corner (4, 4), each step about 0.3 px,
    # so every rim drawn outside its own band would stack into a page-white crescent.
    grid = [SequenceCircle(0.4 + 0.8 * i, 0.4 + 0.8 * j, 0.3, 0.0, 0) for j in range(5) for i in range(5)]
    angles = np.linspace(math.pi, 1.5 * math.pi, 60)
    fan = [SequenceCircle(3.4 + 0.1 * math.cos(a), 3.4 + 0.1 * math.sin(a), 0.25, 0.0, 1) for a in angles]
    circles = tuple(grid + fan)
    return SequencePlotData.build(
        workload="corner-fan",
        title="corner-fan",
        boundary=((0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)),
        circles=circles,
        connectors=(),
        seam_connector=None,
        engagement_deg=(10.0,) * (len(circles) - 1),
        cap_deg=80.0,
        tool_radius=0.5,
        coverage=None,
    )


def test_covered_pocket_with_a_tight_corner_fan_renders_no_page_pixel() -> None:
    # Invariant: nothing is painted outside path (+) D_r and every pixel inside it
    # is band paint, so a covered pocket has no page-colour pixel (raster: 1 px erosion).
    data = _covered_corner_fan()
    image, _, axis = _render_sweep(data)
    pocket = _pocket_mask(image, axis, data.boundary)
    first = data.circles[0]
    row, col = _pixel(axis, first.x + first.radius, first.y, image.shape[0])
    rows, cols = np.ogrid[: image.shape[0], : image.shape[1]]
    pocket &= (rows - row) ** 2 + (cols - col) ** 2 > START_MARKER_EXCLUSION_PX**2
    assert pocket.sum() > 100_000
    assert int(np.sum(pocket & _page_mask(image))) == 0


def _wall_sliver(sliver_mm: float) -> SequencePlotData:
    # Pocket [0.13, 4.13] x [0.13, 4.13], r = 0.5. Vertical connector capsules sweep
    # exactly the strips [x_k - r, x_k + r]; the first starts at 0.13 + sliver_mm,
    # so the strip [0.13, 0.13 + sliver_mm] along the left wall is uncut, exactly.
    left = 0.13
    start = left + sliver_mm
    xs = [start + 0.5 + 0.9 * k for k in range(5)]
    # The circles (each sweeping a disk of radius 0.55) sit on the last strip, far from the sliver.
    circles = tuple(SequenceCircle(xs[-1], 0.6 + 0.6 * k, 0.05, 0.0, 0) for k in range(len(xs) + 1))
    connectors = tuple(((x, -1.0), (x, 5.3)) for x in xs)
    return SequencePlotData.build(
        workload="wall-sliver",
        title="wall-sliver",
        boundary=((left, left), (4.13, left), (4.13, 4.13), (left, 4.13)),
        circles=circles,
        connectors=connectors,
        seam_connector=None,
        engagement_deg=(10.0,) * (len(circles) - 1),
        cap_deg=80.0,
        tool_radius=0.5,
        coverage=None,
    )


def test_uncut_sliver_along_the_wall_shows_page_colour() -> None:
    # The converse: an uncut strip wider than 2 px shows page colour, including
    # the pixels next to the wall (the boundary line must not cover the interior).
    sliver_mm = 0.05
    data = _wall_sliver(sliver_mm)
    image, _, axis = _render_sweep(data)
    pocket = _pocket_mask(image, axis, data.boundary)
    _, sliver_col = _pixel(axis, 0.13 + sliver_mm, 2.13, image.shape[0])
    _, wall_col = _pixel(axis, 0.13, 2.13, image.shape[0])
    assert sliver_col - wall_col > 4
    cols = np.arange(image.shape[1])[None, :]
    core = pocket & (cols < sliver_col - 1)
    assert core.sum() > 0
    assert bool(np.all(_page_mask(image)[core]))


def _write_summary(directory: Path, stem: str, *, circles: int, covered: bool) -> None:
    from benchmarks.tools.held_sequence_summary import CoverageSummary
    from benchmarks.tools.held_sequence_summary import FigureSummary
    from benchmarks.tools.held_sequence_summary import summary_path

    coverage = CoverageSummary.build(covered=covered, witnesses=0 if covered else 8, sites=0 if covered else 2)
    summary = FigureSummary.build(figure=f"{stem}.png", circles=circles, connectors=circles, over_cap=0, unresolved=0, max_known_deg=79.991, seam_deg=180.0, coverage=coverage)
    summary_path(directory / f"{stem}.png").write_text(summary.to_json())


def test_figure_summary_round_trips_and_groups_sites_like_the_plot() -> None:
    from benchmarks.tools.held_sequence_summary import CoverageSummary
    from benchmarks.tools.held_sequence_summary import FigureSummary

    verdict = CoverageVerdict.build(covered=False, witnesses=((1.0, 1.0), (1.1, 1.0), (5.0, 5.0)), piece_kinds=("a", "a", "a"))
    coverage = CoverageSummary.from_verdict(verdict, join_mm=0.5)
    assert (coverage.covered, coverage.witnesses, coverage.sites) == (False, 3, 2)
    summary = FigureSummary.build(figure="x.png", circles=3, connectors=3, over_cap=0, unresolved=0, max_known_deg=12.5, seam_deg=None, coverage=coverage)
    assert FigureSummary.from_json(summary.to_json()) == summary


def test_coverage_summary_rejects_witness_counts_that_contradict_the_verdict() -> None:
    from benchmarks.tools.held_sequence_summary import CoverageSummary
    from benchmarks.tools.held_sequence_summary import InvalidFigureSummaryError

    with pytest.raises(InvalidFigureSummaryError):
        CoverageSummary.build(covered=True, witnesses=2, sites=1)
    with pytest.raises(InvalidFigureSummaryError):
        CoverageSummary.build(covered=False, witnesses=0, sites=0)


def test_gallery_values_rewrite_only_marked_spans_from_the_summaries(tmp_path: Path) -> None:
    from benchmarks.tools.held_sequence_gallery_values import gallery_values
    from benchmarks.tools.held_sequence_gallery_values import rewrite_marked_values

    _write_summary(tmp_path, "held_figure8_monstera_standard_refined_boundary_qualified_contacts_corner_approaches", circles=4440, covered=False)
    values = gallery_values(tmp_path)
    text = "Monstera <!-- gv:monstera.circles -->4,366<!-- /gv --> circles, <!-- gv:monstera.verdict -->x<!-- /gv --> (8 kept)."
    once = rewrite_marked_values(text, values)
    assert once == "Monstera <!-- gv:monstera.circles -->4,440<!-- /gv --> circles, <!-- gv:monstera.verdict -->NOT COVERED<!-- /gv --> (8 kept)."
    assert rewrite_marked_values(once, values) == once


def test_gallery_values_fail_loud_on_a_key_without_a_generated_value(tmp_path: Path) -> None:
    from benchmarks.tools.held_sequence_gallery_values import UnknownGalleryValueError
    from benchmarks.tools.held_sequence_gallery_values import gallery_values
    from benchmarks.tools.held_sequence_gallery_values import rewrite_marked_values

    with pytest.raises(UnknownGalleryValueError, match="skis.circles"):
        rewrite_marked_values("<!-- gv:skis.circles -->1<!-- /gv -->", gallery_values(tmp_path))


def test_gallery_values_check_names_a_stale_page(tmp_path: Path) -> None:
    from benchmarks.tools.held_sequence_gallery_values import StaleGalleryValuesError
    from benchmarks.tools.held_sequence_gallery_values import update_pages

    _write_summary(tmp_path, "held_figure8_crossed_skis_standard_refined_boundary", circles=1621, covered=True)
    page = tmp_path / "page.md"
    page.write_text("<!-- gv:skis.circles -->1,616<!-- /gv -->\n")
    with pytest.raises(StaleGalleryValuesError, match="page.md"):
        update_pages((page,), tmp_path, check=True)
    update_pages((page,), tmp_path, check=False)
    assert page.read_text() == "<!-- gv:skis.circles -->1,621<!-- /gv -->\n"
    update_pages((page,), tmp_path, check=True)


def test_repaired_figure_summary_names_its_reinserted_sources(tmp_path: Path) -> None:
    from benchmarks.tools.held_sequence_gallery_values import gallery_values
    from benchmarks.tools.held_sequence_summary import CoverageSummary
    from benchmarks.tools.held_sequence_summary import FigureSummary
    from benchmarks.tools.held_sequence_summary import summary_path

    covered = CoverageSummary.build(covered=True, witnesses=0, sites=0)
    summary = FigureSummary.build(
        figure="r.png", circles=4430, connectors=4430, over_cap=0, unresolved=0, max_known_deg=79.9, seam_deg=70.0, coverage=covered, inserted_sources=(744, 1273)
    )
    assert FigureSummary.from_json(summary.to_json()) == summary
    name = "held_figure8_monstera_standard_refined_boundary_qualified_contacts_corner_approaches_certified_repair.png"
    summary_path(tmp_path / name).write_text(summary.to_json())
    values = gallery_values(tmp_path)
    assert values["monstera_repaired.inserted_sources"] == "744 and 1273"
    assert values["monstera_repaired.inserted_count"] == "2"
    assert values["monstera_repaired.verdict"] == "CERTIFIED COVERED"


def test_inserted_sources_must_be_distinct_and_non_negative() -> None:
    from benchmarks.tools.held_sequence_summary import FigureSummary
    from benchmarks.tools.held_sequence_summary import InvalidFigureSummaryError

    for bad in ((3, 3), (-1,)):
        with pytest.raises(InvalidFigureSummaryError):
            FigureSummary.build(figure="r.png", circles=2, connectors=2, over_cap=0, unresolved=0, max_known_deg=1.0, seam_deg=None, coverage=None, inserted_sources=bad)


def test_a_zoom_whose_bands_agg_would_clip_at_the_canvas_fails_loud() -> None:
    # Agg cuts unfilled stroke paths at the canvas edge (+1 px) before stroking, so a
    # band wider than the gap between the axes and the canvas edge loses paint inside
    # the axes: page colour where the band cuts. Such a zoom must raise, not mislead.
    from benchmarks.tools.held_sequence_plot import StrokeClippedAtCanvasError

    data = _fixture()
    figure = Figure(figsize=(2.0, 2.0), dpi=100)
    canvas = FigureCanvasAgg(figure)
    axis = figure.add_axes((0.0, 0.0, 1.0, 1.0))
    draw_sweep(axis, data, order_colors(data.circles, ColourMode.ORDER), legend=False)
    axis.set_xlim(4.9, 5.1)
    axis.set_ylim(2.9, 3.1)
    with pytest.raises(StrokeClippedAtCanvasError):
        canvas.draw()
