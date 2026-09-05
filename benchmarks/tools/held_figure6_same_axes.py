"""Render publisher Figure 6 observations and bounded repository lengths."""

from __future__ import annotations

import math
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.image import imread

from benchmarks.held_figure6_comparison import RepositoryFigure6Comparison
from benchmarks.held_figure6_comparison import measure_repository_figure6
from benchmarks.held_figure6_publisher import PublisherFigure6Evidence
from benchmarks.held_figure6_publisher import load_figure6_publisher_evidence

OUTPUT_ROOT = Path(__file__).parents[2] / "docs" / "assets" / "images"
OUTPUT_PNG = OUTPUT_ROOT / "held_figure6_same_axes_comparison.png"
OUTPUT_SVG = OUTPUT_ROOT / "held_figure6_same_axes_comparison.svg"
PUBLISHER_IMAGE = OUTPUT_ROOT / "held_reference_figure6.png"
LITERAL_OVERLAY_PNG = OUTPUT_ROOT / "held_figure6_literal_overlay.png"
PUBLISHER_COLOURS = ("#9400d3", "#009e73", "#56b4e9")
REPOSITORY_COLOUR = "#d55e00"

# The tracked publisher crop is page 13 pixels [155:1175, 175:825]. Subtracting
# this origin from the recorded page-space plot frame registers the overlay.
PUBLISHER_CROP_LEFT_PX = 155
PUBLISHER_CROP_TOP_PX = 175


def build_figure(publisher: PublisherFigure6Evidence, repository: RepositoryFigure6Comparison) -> Figure:
    """Build the comparison figure without generating paths or writing files."""
    figure = Figure(figsize=(10.8, 6.4), facecolor="white")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    for series, colour in zip(publisher.series, PUBLISHER_COLOURS, strict=True):
        axis.plot(
            [float(point.engagement_deg) for point in series.points],
            [float(point.path_length) for point in series.points],
            color=colour,
            linewidth=1.7,
            linestyle=(0, (4, 3)),
            label=f"publisher {series.label} (digitized)",
        )
    axis.plot(
        [float(point.requested_cap_deg) for point in repository.points],
        [float(point.path_length_mm) for point in repository.points],
        color=REPOSITORY_COLOUR,
        linewidth=2.0,
        marker="o",
        markersize=4.5,
        label="repository measured requested-cap; compliance unaudited",
    )
    axes = publisher.axes
    axis.set_xlim(float(axes.x_min_deg), float(axes.x_max_deg))
    axis.set_ylim(float(axes.y_min), float(axes.y_max))
    axis.set_yscale("log")
    tick_count = round((float(axes.x_max_deg) - float(axes.x_min_deg)) / float(axes.x_tick_deg))
    axis.set_xticks([float(axes.x_min_deg) + index * float(axes.x_tick_deg) for index in range(tick_count + 1)])
    axis.set_xlabel("engagement angle / requested cap (degrees)")
    axis.set_ylabel("path length (publisher graphical unit; repository mm)")
    axis.grid(True, which="major", color="#c7c7c7", linewidth=0.65, linestyle=(0, (1, 4)))
    axis.set_title("Held--Pfeiffer Figure 6 / repository requested-cap lengths on identical axes")
    axis.legend(loc="upper right", frameon=False, fontsize=8.5)
    figure.text(
        0.5,
        0.02,
        "Dashed publisher series: pixel-digitized graphical observations. Solid points: repository path lengths in mm.\n"
        "Cap compliance is unaudited; constant-spacing repository curve unavailable because it requires the deferred exact replay.\n"
        "Same axes support visual comparison only: no unpublished numeric parity.",
        ha="center",
        fontsize=8.2,
    )
    figure.subplots_adjust(left=0.105, right=0.98, top=0.91, bottom=0.19)
    return figure


def _publisher_crop_pixel(
    publisher: PublisherFigure6Evidence,
    requested_cap_deg: float,
    path_length_mm: float,
) -> tuple[float, float]:
    axes = publisher.axes
    plot_left = axes.pixel_left - PUBLISHER_CROP_LEFT_PX
    plot_right = axes.pixel_right - PUBLISHER_CROP_LEFT_PX
    plot_top = axes.pixel_top - PUBLISHER_CROP_TOP_PX
    plot_bottom = axes.pixel_bottom - PUBLISHER_CROP_TOP_PX
    x_fraction = (requested_cap_deg - float(axes.x_min_deg)) / (float(axes.x_max_deg) - float(axes.x_min_deg))
    log_min = math.log10(float(axes.y_min))
    log_max = math.log10(float(axes.y_max))
    y_fraction = (log_max - math.log10(path_length_mm)) / (log_max - log_min)
    return plot_left + x_fraction * (plot_right - plot_left), plot_top + y_fraction * (plot_bottom - plot_top)


def build_literal_overlay(publisher: PublisherFigure6Evidence, repository: RepositoryFigure6Comparison) -> Figure:
    """Place repository measurements directly over the unchanged publisher crop."""
    publisher_image = imread(PUBLISHER_IMAGE)
    image_height, image_width = publisher_image.shape[:2]
    figure = Figure(figsize=(image_width / 100.0, (image_height + 80) / 100.0), facecolor="white")
    FigureCanvasAgg(figure)
    axis = figure.add_axes((0.0, 80.0 / (image_height + 80), 1.0, image_height / (image_height + 80)))
    axis.imshow(publisher_image)
    pixels = [_publisher_crop_pixel(publisher, float(point.requested_cap_deg), float(point.path_length_mm)) for point in repository.points]
    axis.plot(
        [pixel[0] for pixel in pixels],
        [pixel[1] for pixel in pixels],
        color=REPOSITORY_COLOUR,
        linewidth=3.0,
        marker="o",
        markersize=7.0,
        markeredgecolor="white",
        markeredgewidth=1.2,
        label="repository requested cap / path length in mm; compliance unaudited",
    )
    axis.set_xlim(-0.5, image_width - 0.5)
    axis.set_ylim(image_height - 0.5, -0.5)
    axis.set_axis_off()
    figure.text(
        0.5,
        0.045,
        "Orange: repository requested-cap path length in mm; cap compliance unaudited. Publisher pixels unchanged.\n"
        "Literal same-numeric-axis inspection only: publisher ordinate is graphical, not common-unit or numeric-parity evidence.",
        ha="center",
        va="center",
        fontsize=8.2,
    )
    return figure


def render(png_path: Path = OUTPUT_PNG, svg_path: Path = OUTPUT_SVG) -> RepositoryFigure6Comparison:
    """Measure once and write the tracked PNG and SVG comparison."""
    publisher = load_figure6_publisher_evidence()
    repository = measure_repository_figure6()
    figure = build_figure(publisher, repository)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_path, dpi=220, facecolor="white")
    figure.savefig(svg_path, facecolor="white", metadata={"Date": None})
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text("\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n", encoding="utf-8")
    literal_overlay = build_literal_overlay(publisher, repository)
    literal_overlay.savefig(LITERAL_OVERLAY_PNG, dpi=100, facecolor="white")
    return repository


if __name__ == "__main__":
    result = render()
    for point in result.points:
        print(f"requested_cap_deg={float(point.requested_cap_deg):g} path_length_mm={float(point.path_length_mm):.9f} compliance={point.compliance}")
