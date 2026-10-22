"""Visual comparison for native coverage-preserving circle deletion."""

from __future__ import annotations

import json
import math
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_motion_coverage import InvalidCoverageMotionError
from benchmarks.held_motion_coverage import connector_pairs
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.tools.held_contour_bound_plot import _length
from benchmarks.tools.held_sequence_plot import DPI
from benchmarks.tools.held_sequence_plot import ENGAGEMENT_PANEL_IN
from benchmarks.tools.held_sequence_plot import FIGURE_WIDTH_IN
from benchmarks.tools.held_sequence_plot import FOOTER_IN
from benchmarks.tools.held_sequence_plot import HEADER_IN
from benchmarks.tools.held_sequence_plot import NOTE_LINE_IN
from benchmarks.tools.held_sequence_plot import PAGE_COLOR
from benchmarks.tools.held_sequence_plot import PANEL_GAP_IN
from benchmarks.tools.held_sequence_plot import ColourMode
from benchmarks.tools.held_sequence_plot import CoverageVerdict
from benchmarks.tools.held_sequence_plot import add_order_colorbar
from benchmarks.tools.held_sequence_plot import colour_legend_text
from benchmarks.tools.held_sequence_plot import draw_engagement
from benchmarks.tools.held_sequence_plot import draw_sweep
from benchmarks.tools.held_sequence_plot import order_colors
from benchmarks.tools.held_sequence_plot import plot_data
from benchmarks.tools.held_sequence_summary import CoverageSummary
from compas_cgal.units import Radian

PAIR_PANEL_IN = 3.9  # Height of the before/after pocket panels, inches.


class InvalidCoveragePreservingPlotError(ValueError):
    """A supplied path lacks its connectors or successor bounds."""


def render_coverage_preserving_comparison(
    case: HeldReferenceCase,
    before: Figure5BoundaryPath,
    after: Figure5BoundaryPath,
    bounds: tuple[Radian, ...],
    output: Path,
    *,
    before_coverage: CoverageVerdict | None = None,
    after_coverage: CoverageVerdict | None = None,
    dense_circles: int | None = None,
    mode: ColourMode = ColourMode.RUN,
) -> None:
    """Render actual supplied paths; no absolute coverage claim is inferred.

    Both paths are open or closed connector sequences (``connector_pairs``);
    the bounds cover the successor pairs only, never the seam pair. Both
    paths are drawn as sequence plots (swept bands in emission order,
    ``held_sequence_plot``); the coverage verdicts, if given, are the exact
    certificate's on each path. ``dense_circles`` is the size of the dense
    source the deletion started from, recorded in the JSON report.

    Raises:
        InvalidCoveragePreservingPlotError: a path is empty, its connector count
            is neither open nor closed, a connector misses its circle contacts,
            or ``bounds`` is not one per emitted successor.
    """
    try:
        connector_pairs(before.circles, before.transitions)
        connector_pairs(after.circles, after.transitions)
    except InvalidCoverageMotionError as error:
        raise InvalidCoveragePreservingPlotError(f"Comparison requires connected circle sequences: {error}") from error
    if len(bounds) != len(after.circles) - 1:
        raise InvalidCoveragePreservingPlotError("Comparison requires one bound per emitted successor.")
    angles = tuple(map(math.degrees, bounds))
    before_length = _length(before.circles, before.transitions)
    after_length = _length(after.circles, after.transitions)
    nan = (math.nan,) * (len(before.circles) - 1)  # The before panel plots no bound.
    before_data = plot_data(workload="before", title="before", case=case, circles=before.circles, transitions=before.transitions, engagement_deg=nan, coverage=before_coverage)
    after_data = plot_data(workload="after", title="after", case=case, circles=after.circles, transitions=after.transitions, engagement_deg=angles, coverage=after_coverage)
    after_colors = order_colors(after_data.circles, mode)
    lines = (
        "Local native deletion preserves dense repaired source coverage and actual connector waypoints.",
        "Source completeness, containment and entry remain unqualified. Preserved coverage is not proof of a cleared pocket.",
        *(() if before_coverage is None else (f"Before · {before_coverage.caption(before_data.tool_radius)}",)),
        *(() if after_coverage is None else (f"After · {after_coverage.caption(after_data.tool_radius)}",)),
        colour_legend_text(mode),
    )
    footer = FOOTER_IN + NOTE_LINE_IN * (len(lines) - 1)
    height = HEADER_IN + PAIR_PANEL_IN + PANEL_GAP_IN + ENGAGEMENT_PANEL_IN + footer
    figure = Figure(figsize=(FIGURE_WIDTH_IN, height), facecolor=PAGE_COLOR, dpi=DPI)
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 3, height_ratios=(PAIR_PANEL_IN, ENGAGEMENT_PANEL_IN), width_ratios=(1, 1, 0.036), wspace=0.18)
    before_axes, after_axes = figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])
    engagement = figure.add_subplot(grid[1, :2])
    draw_sweep(before_axes, before_data, order_colors(before_data.circles, mode), legend=False)
    draw_sweep(after_axes, after_data, after_colors, legend=False)
    add_order_colorbar(figure, figure.add_subplot(grid[:, 2]), mode, len(after.circles))
    for axes, label, path, length, coverage in (
        (before_axes, "Before: engagement-only", before, before_length, before_coverage),
        (after_axes, "After: coverage-preserving deletion", after, after_length, after_coverage),
    ):
        verdict = "" if coverage is None else (" · CERTIFIED COVERED" if coverage.covered else f" · NOT COVERED ({len(coverage.sites(after_data.tool_radius))} site(s))")
        axes.set_title(f"{label}\n{len(path.circles):,} circles · {length:,.1f} mm{verdict}", fontsize=12)
    handles, labels = after_axes.get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.93, 1 - 0.62 / height), ncol=len(labels), frameon=False, fontsize=9)
    cap = float(case.tea_cap)
    draw_engagement(engagement, after_data, after_colors, ylabel="Whole-orbit bound / degrees", ylim=(0.0, max(90.0, cap + 10, max(angles, default=0) + 5)))
    engagement.set_xlabel("Successor in emitted order (after)")
    change = 100 * (after_length / before_length - 1)
    figure.suptitle(f"{case.name.replace('_', ' ')} · preserve swept coverage", x=0.06, y=1 - 0.12 / height, ha="left", va="top", fontsize=19, fontweight="bold")
    figure.text(0.06, 1 - 0.55 / height, f"Length change {change:+.1f}% · maximum reported bound {max(angles, default=0):.3f}° · native deletion decisions", fontsize=11)
    for row, line in enumerate(reversed(lines)):
        figure.text(0.06, (0.15 + NOTE_LINE_IN * row) / height, line, fontsize=9)
    figure.subplots_adjust(left=0.06, right=0.93, top=1 - (HEADER_IN + 0.35) / height, bottom=footer / height, hspace=PANEL_GAP_IN / ((PAIR_PANEL_IN + ENGAGEMENT_PANEL_IN) / 2))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=figure.dpi, facecolor=figure.get_facecolor())
    report = {
        "case": case.name,
        "algorithm": "coverage_preserving_thinning",
        "before_circles": len(before.circles),
        "dense_circles": dense_circles,
        "circles": len(after.circles),
        "connectors": len(after.transitions),
        "before_length_mm": before_length,
        "length_mm": after_length,
        "length_change_percent": change,
        "max_reported_bound_degrees": max(angles, default=0),
        "before_coverage": None if before_coverage is None else CoverageSummary.from_verdict(before_coverage, join_mm=before_data.tool_radius).to_dict(),
        "after_coverage": None if after_coverage is None else CoverageSummary.from_verdict(after_coverage, join_mm=after_data.tool_radius).to_dict(),
        "successor_bounds_degrees": angles,
        "scope": "local native deletion preserves dense repaired source coverage and actual connector waypoints; source completeness, containment and entry unqualified",
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
