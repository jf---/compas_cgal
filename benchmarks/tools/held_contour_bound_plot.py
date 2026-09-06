"""Compare an existing standard draft with a generated contour-bound draft."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.tools.held_toolpath_progress import _draw_motion
from compas_cgal.adaptive.units import Radian


class InvalidContourBoundPlotError(ValueError):
    """Motion and successor-bound arrays do not correspond."""


def _length(circles: tuple[Figure5CounterclockwiseCircle, ...], transitions: tuple[Figure5BoundaryTransition, ...]) -> float:
    return math.tau * math.fsum(float(circle.radius.value) for circle in circles) + math.fsum(float(transition.boundary_progress) for transition in transitions)


def render_contour_bound_comparison(
    case: HeldReferenceCase,
    baseline: tuple[Figure5CounterclockwiseCircle, ...],
    baseline_transitions: tuple[Figure5BoundaryTransition, ...],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    bounds: tuple[Radian, ...],
    output: Path,
) -> None:
    """Render supplied motion and native bound reports without reselection."""
    if not baseline or not circles or len(bounds) != len(circles) - 1 or len(transitions) != len(bounds) or len(baseline_transitions) != len(baseline) - 1:
        raise InvalidContourBoundPlotError("Comparison requires connected circle sequences and one bound per successor.")
    angles = tuple(map(math.degrees, bounds))
    baseline_length, current_length = _length(baseline, baseline_transitions), _length(circles, transitions)
    figure = Figure(figsize=(13, 9), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 2, height_ratios=(3, 1))
    before, after = figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])
    engagement = figure.add_subplot(grid[1, :])
    cap = float(case.tea_cap)
    _draw_motion(before, case, baseline, baseline_transitions, (0.0,) * (len(baseline) - 1), cap)
    _draw_motion(after, case, circles, transitions, angles, cap)
    before.set_title(f"Standard draft · {len(baseline):,} circles\n{baseline_length:,.1f} mm", fontsize=12)
    after.set_title(f"Contour-bound draft · {len(circles):,} circles\n{current_length:,.1f} mm", fontsize=12)
    x = np.arange(1, len(circles))
    engagement.plot(x, angles, color="#087e8b", linewidth=0.8)
    engagement.axhline(cap, color="#c34c32", linestyle="--", label=f"{cap:g}° cap")
    engagement.set(xlabel="Successor in emitted order", ylabel="Whole-orbit bound / degrees", ylim=(0, max(90, cap + 10)))
    engagement.legend(loc="upper right", frameon=False)
    engagement.grid(alpha=0.2)
    reduction = 100 * (1 - current_length / baseline_length)
    figure.suptitle(f"{case.name.replace('_', ' ')} · contour-bound spacing", x=0.075, ha="left", fontsize=19, fontweight="bold")
    figure.text(0.075, 0.925, f"{reduction:.1f}% shorter · maximum reported bound {max(angles, default=0):.3f}° · native squared-chord cap decisions", fontsize=11)
    figure.text(0.075, 0.045, "All earlier outer disks determine the bound. This is a conservative contour-aware variant, not Held's ordered-arc shortcut.", fontsize=9)
    figure.text(0.075, 0.02, "First outer disk assumed cleared. Entry, connector engagement, containment and continuous pocket coverage remain unqualified.", fontsize=9)
    figure.subplots_adjust(left=0.075, right=0.97, top=0.85, bottom=0.13, hspace=0.3, wspace=0.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    report = {
        "case": case.name,
        "algorithm": "contour_clearance_bound",
        "baseline_circles": len(baseline),
        "circles": len(circles),
        "connectors": len(transitions),
        "baseline_length_mm": baseline_length,
        "length_mm": current_length,
        "length_reduction_percent": reduction,
        "max_reported_bound_degrees": max(angles, default=0),
        "successor_bounds_degrees": angles,
        "scope": "leading-semicircle bound under filled-disk history and cleared initial outer disk; entry/connectors/containment/coverage unqualified",
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
