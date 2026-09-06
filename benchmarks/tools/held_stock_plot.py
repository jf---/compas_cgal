"""Render circle-only native stock at one generated toolpath prefix."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_stock_replay import replay_circle_stock

DISPLAY_SAMPLES_PER_TOOL_RADIUS = 4  # Raster display pitch only; native stock is not rasterized.


def render_circle_stock_prefix(case: HeldReferenceCase, circles: tuple[Figure5CounterclockwiseCircle, ...], output: Path) -> None:
    """Show supplied machining circles and remaining material after their sweeps.

    Pixel membership comes from CGAL. Entry and connectors are absent from
    the replay, so the panel is an intermediate consumer for Figure 5(c), not
    a completed reproduction of its contour-aware stock state.
    """
    stock = replay_circle_stock(case.projection.points, circles, case.tool_radius)
    uncut = replay_circle_stock(case.projection.points, (), case.tool_radius)
    boundary = np.array([(float(point.x), float(point.y)) for point in case.projection.points])
    pitch = float(case.tool_radius.value) / DISPLAY_SAMPLES_PER_TOOL_RADIUS
    xs = np.arange(boundary[:, 0].min(), boundary[:, 0].max(), pitch)
    ys = np.arange(boundary[:, 1].min(), boundary[:, 1].max(), pitch)
    pixels = np.zeros((len(ys), len(xs)), dtype=np.uint8)
    for row, y in enumerate(ys):
        for column, x in enumerate(xs):
            if stock.contains(float(x), float(y)):
                pixels[row, column] = 2
            elif uncut.contains(float(x), float(y)):
                pixels[row, column] = 1
    figure = Figure(figsize=(12, 6), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    motion_axis, stock_axis = figure.subplots(1, 2)
    closed = np.vstack((boundary, boundary[0]))
    for axis in (motion_axis, stock_axis):
        axis.plot(*closed.T, color="#34403f", linewidth=1)
        axis.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    for circle in circles:
        motion_axis.add_patch(Circle((float(circle.center.x), float(circle.center.y)), float(circle.radius.value), fill=False, edgecolor="#087e8b", linewidth=0.4))
    motion_axis.set_title(f"First {len(circles):,} machining circles")
    # Pixel centres coincide with the coordinates queried from the exact stock.
    stock_axis.imshow(
        pixels,
        origin="lower",
        extent=(xs[0] - pitch / 2, xs[-1] + pitch / 2, ys[0] - pitch / 2, ys[-1] + pitch / 2),
        interpolation="nearest",
        cmap=ListedColormap(("#faf9f5", "#91d0ce", "#575c60")),
        vmin=0,
        vmax=2,
    )
    stock_axis.set_title("Remaining material (gray) · circle sweeps (teal)")
    padding = float(case.tool_radius.value)
    for axis in (motion_axis, stock_axis):
        axis.set_xlim(boundary[:, 0].min() - padding, boundary[:, 0].max() + padding)
        axis.set_ylim(boundary[:, 1].min() - padding, boundary[:, 1].max() + padding)
    figure.suptitle(f"{case.name.replace('_', ' ')} · accumulated circle-stock prefix", x=0.07, ha="left", fontsize=17)
    figure.text(0.07, 0.055, "Exact CGAL annulus unions on stored circles · uncut centre islands preserved · entry and connector clearing excluded.", fontsize=9)
    figure.text(0.07, 0.025, f"Display membership sampled at {pitch:g} mm pitch. Stock state is native geometry; no contour-aware placement or engagement claim.", fontsize=9)
    figure.subplots_adjust(left=0.07, right=0.98, top=0.85, bottom=0.15, wspace=0.24)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
