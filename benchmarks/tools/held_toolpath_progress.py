"""Draw complete emitted Figure 5 motion and its predecessor engagement."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection

from benchmarks.held_containment_report import containment_report
from benchmarks.held_exact_motion_coverage import require_full_motion_coverage
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_toolpath_benchmark import generate
from benchmarks.tools.held_contour_contact_plot import render_contour_contacts
from benchmarks.tools.held_sequence_coverage import certify_path_coverage
from benchmarks.tools.held_sequence_plot import ColourMode
from benchmarks.tools.held_sequence_plot import CoverageVerdict
from benchmarks.tools.held_sequence_plot import emitted_engagement
from benchmarks.tools.held_sequence_plot import plot_data
from benchmarks.tools.held_sequence_plot import render_sequence
from benchmarks.tools.held_sequence_summary import CoverageSummary
from benchmarks.tools.held_sequence_summary import FigureSummary
from benchmarks.tools.held_sequence_summary import summary_path
from benchmarks.tools.held_stock_plot import render_circle_stock_prefix

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
    coverage: CoverageVerdict | None = None,
    mode: ColourMode = ColourMode.RUN,
    inserted_sources: tuple[int, ...] | None = None,
) -> None:
    """Render the supplied execution order without sorting or selecting moves, as a sequence plot.

    Each motion's swept band is painted in emission order (``held_sequence_plot``),
    coloured per guide run by default. The lower panel recomputes the standard
    predecessor model from that order; it does not certify depleted-stock
    engagement or entry/connector safety. A closed path (one connector per
    circle) also reports the seam pair, last circle back to the first; it is a
    reported value only. ``coverage`` is the certificate verdict drawn over
    the path, if one was run. The printed metrics and verdict are also
    written beside the figure (``held_sequence_summary``) for the docs;
    ``inserted_sources`` (a certified repair's re-inserted sources) is
    printed on the figure and recorded there.
    """
    angles, seam = emitted_engagement(case, circles, transitions)
    cap = float(case.tea_cap)
    # Successor pairs only; the seam pair is reported on its own.
    failed = sum(1 for angle in angles if angle > cap)
    unresolved = sum(1 for angle in angles if not math.isfinite(angle))
    maximum = max((angle for angle in angles if math.isfinite(angle)), default=0.0)
    seam_text = "" if seam is None else (" · seam unresolved" if not math.isfinite(seam) else f" · seam {seam:.2f}°{' (over cap)' if seam > cap else ''}")
    data = plot_data(workload=output.stem, title=title, case=case, circles=circles, transitions=transitions, engagement_deg=angles, coverage=coverage)
    render_sequence(
        data,
        output,
        mode,
        subtitle=f"{len(circles):,} circles · {len(transitions):,} connectors · {failed:,} over cap · {unresolved} unresolved · max known {maximum:.2f}°{seam_text}",
        notes=(
            scope,
            "Predecessor-circle model only; entry, connectors and depleted-stock engagement require separate checks.",
            *(() if inserted_sources is None else (f"Certificate-gated coverage repair re-inserted source(s) {', '.join(map(str, inserted_sources)) or 'none'}.",)),
        ),
    )
    summary = FigureSummary.build(
        figure=output.name,
        circles=len(circles),
        connectors=len(transitions),
        over_cap=failed,
        unresolved=unresolved,
        max_known_deg=maximum,
        seam_deg=seam,
        coverage=None if coverage is None else CoverageSummary.from_verdict(coverage, join_mm=data.tool_radius),
        inserted_sources=inserted_sources,
    )
    summary_path(output).write_text(summary.to_json())


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refine", action="store_true", help="The generator's refined path (figure5_refined: seam closed, connector post-gate); report containment failures.")
    parser.add_argument("--stock-prefix", type=int, help="Also render native circle-only stock after this many emitted circles.")
    parser.add_argument("--contour-prefix", type=int, help="Inspect native exposed-contour contacts on this many existing circles.")
    parser.add_argument("--check-coverage", action="store_true", help="Also run the exact full-motion coverage replay; unworkable at this size, see docs/held_motion_coverage.md.")
    parser.add_argument("--output-dir", type=Path, default=Path("docs/assets/images"), help="Directory every figure is written to.")
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
    output = args.output_dir / "held_figure5_toolpath_current.png"
    title = "Figure 5 · full-path failure map"
    scope = "Polygon guide, publisher-derived start, original lane placement. Over-cap moves are defects, not accepted motion."
    if args.refine:
        # The published path is the benchmark generator's own (seam closed, connector post-gate), not a re-derivation.
        generated = generate("figure5_refined")
        circles, transitions = generated.circles, generated.transitions
        rejected, excess = containment_report(circles, components[0])
        output = args.output_dir / "held_figure5_engagement_refined.png"
        title = "Figure 5 · corner-aware engagement refinement"
        scope = f"Publisher-start polygon proposal · {rejected:,} exact containment rejections · max reported protrusion {float(excess):.2g} mm · not machining-qualified."
        print(scope)
    if args.stock_prefix is not None:
        if args.stock_prefix > len(circles):
            parser.error("Stock prefix exceeds emitted circle count.")
        render_circle_stock_prefix(case, circles[: args.stock_prefix], args.output_dir / f"held_figure5_circle_stock_{args.stock_prefix}.png")
    if args.contour_prefix is not None:
        if args.contour_prefix > len(circles):
            parser.error("Contour prefix exceeds emitted circle count.")
        render_contour_contacts(case, circles[: args.contour_prefix], args.output_dir / f"held_figure5_contour_contacts_{args.contour_prefix}.png")
    render_path(
        case,
        circles,
        transitions,
        output,
        title=title,
        scope=scope,
        coverage=certify_path_coverage(components[0], circles, case.tool_radius),
    )
    if args.check_coverage:
        require_full_motion_coverage(case, circles, transitions, output.with_suffix(".exact-coverage.json"))
    else:
        print("exact full-motion coverage NOT checked (set aside 2026-09-08; pass --check-coverage to run the replay)", flush=True)
    print(output)


if __name__ == "__main__":
    main()
