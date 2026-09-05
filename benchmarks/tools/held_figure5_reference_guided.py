"""Render the tracked Figure 5 reference-guided circle correspondence."""

from __future__ import annotations

import math
import statistics
from pathlib import Path

from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from benchmarks.held_figure5_publisher import Figure5PublisherCubicEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherLineEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherPathEvidence
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_reference_guided import ReferenceGuidedFigure5Correspondence
from benchmarks.held_figure5_reference_guided import build_reference_guided_figure5_correspondence
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset

OUTPUT = Path(__file__).parents[2] / "docs/assets/images/held_figure5_reference_guided_overlay.png"
PUBLISHER_COLOUR = "#8f35f4"
REPOSITORY_COLOUR = "#e63946"
BOUNDARY_COLOUR = "#238b57"
TURN_SAMPLES = 96


def _publisher_segments(evidence: Figure5PublisherPathEvidence) -> tuple[tuple[tuple[float, float], ...], ...]:
    segments: list[tuple[tuple[float, float], ...]] = []
    for turn in evidence.turns:
        phase = math.atan2(float(turn.entry.y) - float(turn.center.y), float(turn.entry.x) - float(turn.center.x))
        segments.append(
            tuple(
                (
                    float(turn.center.x) + float(turn.radius) * math.cos(phase + float(turn.signed_sweep) * index / TURN_SAMPLES),
                    float(turn.center.y) + float(turn.radius) * math.sin(phase + float(turn.signed_sweep) * index / TURN_SAMPLES),
                )
                for index in range(TURN_SAMPLES + 1)
            )
        )
        for primitive in turn.connector_after.primitives:
            if isinstance(primitive, Figure5PublisherLineEvidence):
                segments.append(((float(primitive.start.x), float(primitive.start.y)), (float(primitive.end.x), float(primitive.end.y))))
            elif isinstance(primitive, Figure5PublisherCubicEvidence):
                points = (primitive.start, primitive.control1, primitive.control2, primitive.end)
                samples = []
                for index in range(17):
                    parameter = index / 16
                    weights = ((1 - parameter) ** 3, 3 * (1 - parameter) ** 2 * parameter, 3 * (1 - parameter) * parameter**2, parameter**3)
                    samples.append(
                        (
                            sum(weight * float(point.x) for weight, point in zip(weights, points, strict=True)),
                            sum(weight * float(point.y) for weight, point in zip(weights, points, strict=True)),
                        )
                    )
                segments.append(tuple(samples))
    return tuple(segments)


def _repository_segments(result: ReferenceGuidedFigure5Correspondence) -> tuple[tuple[tuple[float, float], ...], ...]:
    segments = []
    for circle in result.circles:
        item = circle.hypothesis
        radius = float(item.guide_radius.value)
        phase = math.atan2(float(item.contact_point.y) - float(item.center.y), float(item.contact_point.x) - float(item.center.x))
        segments.append(
            tuple(
                (
                    float(item.center.x) + radius * math.cos(phase + math.tau * index / TURN_SAMPLES),
                    float(item.center.y) + radius * math.sin(phase + math.tau * index / TURN_SAMPLES),
                )
                for index in range(TURN_SAMPLES + 1)
            )
        )
    return tuple(segments)


def _draw_boundary(axis: Axes, points: tuple[tuple[float, float], ...]) -> None:
    closed = points + (points[0],)
    axis.plot(*(zip(*closed, strict=True)), color=BOUNDARY_COLOUR, linewidth=1.4, zorder=3)


def render(output: Path = OUTPUT) -> ReferenceGuidedFigure5Correspondence:
    """Render from tracked evidence and unchanged repository hypotheses."""
    case = load_held_reference_case("figure5")
    evidence = load_figure5_publisher_path()
    guide = build_figure5_raw_guide(case)
    result = build_reference_guided_figure5_correspondence(guide, figure7_inward_offset(case).components)
    figure = Figure(figsize=(12, 8), facecolor="white")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    axis.add_collection(LineCollection(_publisher_segments(evidence), colors=PUBLISHER_COLOUR, linewidths=0.45, alpha=0.8, zorder=2))
    axis.add_collection(LineCollection(_repository_segments(result), colors=REPOSITORY_COLOUR, linewidths=0.68, linestyles=(0, (4, 3)), alpha=0.9, zorder=4))
    boundary = tuple((float(point.x), float(point.y)) for point in case.projection.points)
    _draw_boundary(axis, boundary)
    axis.set_title(f"Figure 5(a) reference-guided circle correspondence ({len(result.circles)} circles)", fontsize=14)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-2, 69)
    axis.set_ylim(-2, 48)
    axis.set_axis_off()
    residuals = tuple(map(float, result.correspondence.circle_locus_residuals))
    figure.text(
        0.5,
        0.018,
        "Solid violet: fitted publisher turns + tracked connector primitives   Dashed red: unchanged repository circles\n"
        f"circle-locus median/max {statistics.median(residuals):.3f}/{max(residuals):.3f} mm; "
        f"{result.correspondence.outlier_count} above 0.150243 mm; 264 unresolved transitions; no independent path/topology parity",
        ha="center",
        fontsize=9.5,
    )
    figure.subplots_adjust(left=0.015, right=0.985, top=0.94, bottom=0.075)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=240, facecolor="white")
    return result


if __name__ == "__main__":
    render()
