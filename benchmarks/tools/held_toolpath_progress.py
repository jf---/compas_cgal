"""Draw complete emitted Figure 5 motion and its predecessor engagement."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import maximum_predecessor_engagement

TURN_DISPLAY_SAMPLES = 96  # Display resolution only; no geometric decisions.
PATH_COLOR = "#087e8b"
FAILURE_COLOR = "#c34c32"


def render_path(
    case: HeldReferenceCase,
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    output: Path,
    *,
    title: str,
    scope: str,
) -> None:
    """Render the supplied execution order without sorting or selecting moves.

    The lower panel recomputes the standard predecessor model from that order;
    it does not certify depleted-stock engagement or entry/connector safety.
    """
    measured = []
    unresolved = []
    for index, (first, second) in enumerate(zip(circles, circles[1:]), start=1):
        try:
            angle = math.degrees(float(maximum_predecessor_engagement(_paper_candidate(first), _paper_candidate(second), case.tool_radius)))
        except StandardPlacementFragmentationError:
            # Missing model evidence is visibly unresolved, never a passing
            # angle. NaN breaks the plot at this adjacency.
            angle = math.nan
            unresolved.append(index)
        measured.append(angle)
    angles = tuple(measured)
    cap = float(case.tea_cap)
    figure = Figure(figsize=(13, 10), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 1, height_ratios=(3, 1))
    path_axis = figure.add_subplot(grid[0])
    engagement_axis = figure.add_subplot(grid[1])
    _draw_motion(path_axis, case, circles, transitions, angles, cap)
    indices = np.arange(1, len(circles))
    engagement_axis.plot(indices, angles, color=PATH_COLOR, linewidth=0.8)
    failed = [(index, angle) for index, angle in zip(indices, angles) if angle > cap]
    if failed:
        engagement_axis.scatter([item[0] for item in failed], [item[1] for item in failed], color=FAILURE_COLOR, s=9, zorder=3)
    if unresolved:
        engagement_axis.scatter(unresolved, [185] * len(unresolved), marker="x", color="#6e39a6", s=24, label="Unresolved model (marker only)")
    engagement_axis.axhline(cap, color=FAILURE_COLOR, linestyle="--", linewidth=1, label=f"{cap:g}° cap")
    engagement_axis.set(xlabel="Successor index in emitted order", ylabel="Engagement / degrees", ylim=(0, 190))
    engagement_axis.legend(loc="upper right", frameon=False)
    engagement_axis.grid(alpha=0.2)
    maximum = max((angle for angle in angles if math.isfinite(angle)), default=0.0)
    figure.suptitle(title, x=0.075, ha="left", fontsize=19, fontweight="bold")
    figure.text(0.075, 0.935, f"{len(circles):,} circles · {len(transitions):,} connectors · {len(failed):,} over cap · {len(unresolved)} unresolved · max known {maximum:.2f}°")
    figure.text(0.075, 0.035, scope, fontsize=9)
    figure.text(0.075, 0.015, "Predecessor-circle model only; entry, connectors, depleted-stock engagement and coverage require separate checks.", fontsize=9)
    figure.subplots_adjust(left=0.075, right=0.97, top=0.90, bottom=0.10, hspace=0.22)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, facecolor=figure.get_facecolor())


def _draw_motion(
    axis: Axes,
    case: HeldReferenceCase,
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    angles: tuple[float, ...],
    cap: float,
) -> None:
    boundary = [(float(point.x), float(point.y)) for point in case.projection.points]
    axis.plot(*zip(*(boundary + boundary[:1])), color="#34403f", linewidth=1.2)
    turns = []
    colors = []
    for index, circle in enumerate(circles):
        center = np.array((float(circle.center.x), float(circle.center.y)))
        phase = math.atan2(float(circle.contact_point.y) - center[1], float(circle.contact_point.x) - center[0])
        parameters = np.linspace(phase, phase + math.tau, TURN_DISPLAY_SAMPLES + 1)
        points = center + float(circle.radius.value) * np.column_stack((np.cos(parameters), np.sin(parameters)))
        # Preserve the emitted contact in the display, avoiding trigonometric
        # roundoff gaps at the start/end of the plotted full circle.
        points[0] = points[-1] = (float(circle.contact_point.x), float(circle.contact_point.y))
        turns.append(points)
        if index > 0 and not math.isfinite(angles[index - 1]):
            colors.append("#6e39a6")
        else:
            colors.append(FAILURE_COLOR if index > 0 and angles[index - 1] > cap else PATH_COLOR)
    axis.add_collection(LineCollection(turns, colors=colors, linewidths=0.38, alpha=0.7))
    axis.add_collection(
        LineCollection(
            [[(float(point.x), float(point.y)) for point in transition.samples] for transition in transitions],
            colors="#b77c29",
            linewidths=0.8,
            alpha=0.9,
        )
    )
    if circles:
        terminal = transitions[-1].samples[-1] if len(transitions) == len(circles) else circles[-1].contact_point
        for point, marker, label in ((circles[0].contact_point, "o", "Start"), (terminal, "s", "Terminal")):
            axis.scatter(float(point.x), float(point.y), marker=marker, s=45, edgecolor="white", label=label, zorder=5)
    axis.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    axis.legend(loc="upper right", frameon=False)
    axis.set_facecolor("#faf9f5")


def main() -> None:
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)
    result = build_figure5_approximate_path(guide, figure7_inward_offset(case).components)
    output = Path("docs/assets/images/held_figure5_toolpath_current.png")
    render_path(
        case,
        result.path.circles,
        result.path.transitions,
        output,
        title="Figure 5 · full-path failure map",
        scope="Polygon guide, publisher-derived start, original lane placement. Over-cap moves are defects, not accepted motion.",
    )
    print(output)


if __name__ == "__main__":
    main()
