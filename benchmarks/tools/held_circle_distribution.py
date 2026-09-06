"""Reporting-only circle distribution diagnostics for generated Held drafts."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.held_figure5_engagement_refinement import Figure5EngagementRefinement
from benchmarks.held_figure5_path import HypothesisFigure5Path
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_standard_placement import maximum_predecessor_engagement


def render_distribution(case: HeldReferenceCase, initial: HypothesisFigure5Path, refined: Figure5EngagementRefinement, output: Path) -> None:
    """Plot source/inserted centres and final pairwise engagement, without selection.

    The JSON sibling records reporting measurements. Angles are predecessor-only;
    circle circumference excludes connectors and is not complete path length.
    """
    circles = refined.circles
    original = frozenset(refined.original_indices)
    angles = np.array([math.degrees(float(maximum_predecessor_engagement(_paper_candidate(a), _paper_candidate(b), case.tool_radius))) for a, b in zip(circles, circles[1:])])
    centers = np.array([(float(circle.center.x), float(circle.center.y)) for circle in circles])
    inserted = np.array([index not in original for index in range(len(circles))])
    cap = float(case.tea_cap)
    report = {
        "case": case.name,
        "model": "standard predecessor only; approximate polygon draft",
        "initial_circles": len(initial.path.circles),
        "refined_circles": len(circles),
        "inserted_circles": int(inserted.sum()),
        "placement_lanes": initial.placement_lane_count,
        "source_runs": len(initial.reached_run_ids),
        "initial_terminal_selections": sum(record.placement.reached_terminal for record in initial.placements),
        "initial_forced_selections": sum(record.placement.forced for record in initial.placements),
        "final_pairs_over_cap": int((angles > cap).sum()),
        "final_pairs_below_half_cap": int((angles < cap / 2).sum()),
        "final_angle_quantiles_degrees": dict(zip(("min", "q25", "median", "q75", "max"), np.quantile(angles, [0, 0.25, 0.5, 0.75, 1]).tolist())),
        "initial_circle_circumference_mm": sum(math.tau * float(circle.radius.value) for circle in initial.path.circles),
        "refined_circle_circumference_mm": sum(math.tau * float(circle.radius.value) for circle in circles),
        "inserted_circle_circumference_mm": sum(math.tau * float(circle.radius.value) for index, circle in enumerate(circles) if index not in original),
    }
    figure = Figure(figsize=(12, 9), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    source_axis, angle_axis, count_axis, histogram_axis = figure.subplots(2, 2).flat
    boundary = [(float(point.x), float(point.y)) for point in case.projection.points]
    for axis in (source_axis, angle_axis):
        axis.plot(*zip(*(boundary + boundary[:1])), color="#34403f", linewidth=0.9)
        axis.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    source_axis.scatter(*centers[~inserted].T, s=3, color="#087e8b", label="Repaired source circles")
    source_axis.scatter(*centers[inserted].T, s=4, color="#c34c32", label="Inserted circles")
    source_axis.set_title("Where refinement adds circle centres")
    source_axis.legend(frameon=False, fontsize=8)
    colored = angle_axis.scatter(*centers[1:].T, c=angles, s=5, vmin=0, vmax=cap, cmap="viridis")
    angle_axis.set_title("Final predecessor engagement at each centre")
    figure.colorbar(colored, ax=angle_axis, label="Degrees")
    bars = count_axis.bar(("Initial", "Refined"), (len(initial.path.circles), len(circles)), color=("#087e8b", "#c34c32"))
    count_axis.bar_label(bars, padding=3)
    count_axis.set(ylabel="Machining circles", ylim=(0, len(circles) * 1.15), title=f"{len(original):,} source circles + {int(inserted.sum()):,} insertions")
    histogram_axis.hist(angles, bins=np.linspace(0, max(cap, float(angles.max())), 21), color="#087e8b")
    histogram_axis.axvline(cap, color="#c34c32", linestyle="--", label=f"{cap:g}° cap")
    histogram_axis.set(xlabel="Final predecessor engagement / degrees", ylabel="Successor pairs", title="Cap compliance does not imply efficient spacing")
    histogram_axis.legend(frameon=False)
    figure.suptitle(f"{case.name.replace('_', ' ')} · circle distribution", fontsize=18, x=0.08, ha="left")
    figure.text(0.08, 0.035, "Centre maps, not toolpaths. Publisher-derived start; no contour-aware, coverage or machining qualification.", fontsize=9)
    figure.subplots_adjust(left=0.08, right=0.96, top=0.91, bottom=0.11, hspace=0.35, wspace=0.32)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
