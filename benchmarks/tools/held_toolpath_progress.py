"""Draw complete emitted Figure 5 motion and its predecessor engagement."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from benchmarks.tools.held_contour_contact_plot import render_contour_contacts
from benchmarks.tools.held_stock_plot import render_circle_stock_prefix
from compas_cgal import _coverage_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

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


def _containment_report(circles: tuple[Figure5CounterclockwiseCircle, ...], boundary: tuple[Point2[WorldXY], ...]) -> tuple[int, Millimetre]:
    polygon = np.array([(float(point.x), float(point.y), 0.0) for point in boundary])
    starts = polygon[:, :2]
    edges = np.roll(starts, -1, axis=0) - starts
    rejected = 0
    maximum_excess = 0.0
    for circle in circles:
        radius = float(circle.radius.value)
        center = np.array((float(circle.center.x), float(circle.center.y)))
        if _coverage_2.CutterCentreDomain2.build(polygon, [], radius).contains(float(center[0]), float(center[1])):
            continue
        rejected += 1
        # Reporting only: the exact native predicate above owns the rejection.
        parameters = np.clip(np.sum((center - starts) * edges, axis=1) / np.sum(edges * edges, axis=1), 0.0, 1.0)
        clearance = float(np.linalg.norm(center - starts - parameters[:, None] * edges, axis=1).min())
        maximum_excess = max(maximum_excess, radius - clearance)
    return rejected, Millimetre(maximum_excess)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refine", action="store_true", help="Repair corner circles and subdivide over-cap gaps; report remaining containment failures.")
    parser.add_argument("--stock-prefix", type=int, help="Also render native circle-only stock after this many emitted circles.")
    parser.add_argument("--contour-prefix", type=int, help="Inspect native exposed-contour contacts on this many existing circles.")
    args = parser.parse_args()
    if args.stock_prefix is not None and args.stock_prefix <= 0:
        parser.error("Stock prefix must be positive.")
    if args.contour_prefix is not None and args.contour_prefix < 2:
        parser.error("Contour prefix requires at least two circles.")
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)
    components = figure7_inward_offset(case).components
    result = build_figure5_approximate_path(guide, components)
    circles, transitions = result.path.circles, result.path.transitions
    output = Path("docs/assets/images/held_figure5_toolpath_current.png")
    title = "Figure 5 · full-path failure map"
    scope = "Polygon guide, publisher-derived start, original lane placement. Over-cap moves are defects, not accepted motion."
    if args.refine:
        refined = refine_figure5_engagement(circles, components[0], case.tool_radius, EngagementCap.build(math.radians(float(case.tea_cap))))
        circles = refined.circles
        lengths = _side_lengths(components[0])
        transitions = tuple(_transition(components[0], lengths, a, b) for a, b in zip(circles, circles[1:]))
        rejected, excess = _containment_report(circles, components[0])
        output = Path("docs/assets/images/held_figure5_engagement_refined.png")
        title = "Figure 5 · corner-aware engagement refinement"
        scope = f"Publisher-start polygon proposal · {rejected:,} exact containment rejections · max reported protrusion {float(excess):.2g} mm · not machining-qualified."
        print(scope)
    if args.stock_prefix is not None:
        if args.stock_prefix > len(circles):
            parser.error("Stock prefix exceeds emitted circle count.")
        render_circle_stock_prefix(case, circles[: args.stock_prefix], Path("docs/assets/images") / f"held_figure5_circle_stock_{args.stock_prefix}.png")
    if args.contour_prefix is not None:
        if args.contour_prefix > len(circles):
            parser.error("Contour prefix exceeds emitted circle count.")
        render_contour_contacts(case, circles[: args.contour_prefix], Path("docs/assets/images") / f"held_figure5_contour_contacts_{args.contour_prefix}.png")
    render_path(
        case,
        circles,
        transitions,
        output,
        title=title,
        scope=scope,
    )
    print(output)


if __name__ == "__main__":
    main()
