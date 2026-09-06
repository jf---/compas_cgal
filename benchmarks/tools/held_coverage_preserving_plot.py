"""Visual comparison for native coverage-preserving circle deletion."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.tools.held_contour_bound_plot import _length
from benchmarks.tools.held_toolpath_progress import _draw_motion
from compas_cgal.adaptive.units import Radian


class InvalidCoveragePreservingPlotError(ValueError):
    """A supplied path lacks its connectors or successor bounds."""


def render_coverage_preserving_comparison(
    case: HeldReferenceCase,
    before: Figure5BoundaryPath,
    after: Figure5BoundaryPath,
    bounds: tuple[Radian, ...],
    output: Path,
) -> None:
    """Render actual supplied paths; no absolute coverage claim is inferred."""
    if any(not path.circles or len(path.transitions) != len(path.circles) - 1 for path in (before, after)) or len(bounds) != len(after.circles) - 1:
        raise InvalidCoveragePreservingPlotError("Comparison requires connected paths and one bound per emitted successor.")
    angles = tuple(map(math.degrees, bounds))
    before_length = _length(before.circles, before.transitions)
    after_length = _length(after.circles, after.transitions)
    figure = Figure(figsize=(13, 9), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 2, height_ratios=(3, 1))
    before_axes, after_axes = figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])
    engagement = figure.add_subplot(grid[1, :])
    cap = float(case.tea_cap)
    _draw_motion(before_axes, case, before.circles, before.transitions, (0.0,) * len(before.transitions), cap)
    _draw_motion(after_axes, case, after.circles, after.transitions, angles, cap)
    before_axes.set_title(f"Before: engagement-only\n{len(before.circles):,} circles · {before_length:,.1f} mm", fontsize=12)
    after_axes.set_title(f"After: coverage-preserving deletion\n{len(after.circles):,} circles · {after_length:,.1f} mm", fontsize=12)
    engagement.plot(np.arange(1, len(after.circles)), angles, color="#087e8b", linewidth=0.8)
    engagement.axhline(cap, color="#c34c32", linestyle="--", label=f"{cap:g}° cap")
    engagement.set(xlabel="Successor in emitted order", ylabel="Whole-orbit bound / degrees", ylim=(0, max(90, cap + 10, max(angles, default=0) + 5)))
    engagement.legend(loc="upper right", frameon=False)
    engagement.grid(alpha=0.2)
    change = 100 * (after_length / before_length - 1)
    figure.suptitle(f"{case.name.replace('_', ' ')} · preserve swept coverage", x=0.075, ha="left", fontsize=19, fontweight="bold")
    figure.text(0.075, 0.925, f"Length change {change:+.1f}% · maximum reported bound {max(angles, default=0):.3f}° · native deletion decisions", fontsize=11)
    figure.text(0.075, 0.045, "Local native deletion preserves dense repaired source coverage and actual connector waypoints.", fontsize=9)
    figure.text(0.075, 0.02, "Source completeness, containment and entry remain unqualified. Preserved coverage is not proof of a cleared pocket.", fontsize=9)
    figure.subplots_adjust(left=0.075, right=0.97, top=0.85, bottom=0.13, hspace=0.3, wspace=0.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    report = {
        "case": case.name,
        "algorithm": "coverage_preserving_thinning",
        "before_circles": len(before.circles),
        "circles": len(after.circles),
        "connectors": len(after.transitions),
        "before_length_mm": before_length,
        "length_mm": after_length,
        "length_change_percent": change,
        "max_reported_bound_degrees": max(angles, default=0),
        "successor_bounds_degrees": angles,
        "scope": "local native deletion preserves dense repaired source coverage and actual connector waypoints; source completeness, containment and entry unqualified",
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
