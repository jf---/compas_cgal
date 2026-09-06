"""Plot exact Figure 5 boundary contacts and the remaining toolpath scope."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import _coverage_2

OUTPUT = Path("docs/assets/images/held_figure5_exact_connector.png")
ARC_DISPLAY_SAMPLES = 81  # Rendering density; never used by native geometry.
LINE_COLOR = "#087e8b"
ARC_COLOR = "#bd4c21"


def _display_points(primitive: _coverage_2.ReachableBoundaryPrimitive2) -> npt.NDArray[np.float64]:
    if primitive.kind == "line":
        return np.asarray((primitive.start_mm, primitive.end_mm))
    center = np.asarray(primitive.arc_center_mm)
    start = np.asarray(primitive.start_mm) - center
    end = np.asarray(primitive.end_mm) - center
    first = math.atan2(start[1], start[0])
    last = math.atan2(end[1], end[0])
    sweep = (last - first) % math.tau if primitive.arc_counterclockwise else -((first - last) % math.tau)
    angles = np.linspace(first, first + sweep, ARC_DISPLAY_SAMPLES)
    points = center + primitive.arc_radius_mm * np.column_stack((np.cos(angles), np.sin(angles)))
    points[0], points[-1] = primitive.start_mm, primitive.end_mm
    return points


def main() -> None:
    case = load_held_reference_case("figure5")
    boundary = np.asarray([(float(point.x), float(point.y), 0.0) for point in case.projection.points])
    cycle = _coverage_2.build_center_boundary_cycle(boundary, [], float(case.tool_radius.value))
    primitives = cycle.primitives
    # Choose a visible bend for display only; this never drives generation.
    arc_index = max(
        (
            index
            for index, primitive in enumerate(primitives)
            if primitive.kind == "arc" and primitives[(index - 1) % len(primitives)].kind == "line" and primitives[(index + 1) % len(primitives)].kind == "line"
        ),
        key=lambda index: float(np.linalg.norm(np.diff(_display_points(primitives[index]), axis=0), axis=1).sum()),
    )
    expected = tuple(primitives[(arc_index + offset) % len(primitives)] for offset in (-1, 0, 1))
    connector = cycle.ccw_transition(expected[0].start, expected[-1].end)
    if len(connector) != len(expected):
        raise _coverage_2.ReachableArrangementTopologyError("Display connector lost its source interval.")
    for actual, source in zip(connector, expected):
        if actual.start != source.start or actual.end != source.end or actual.source_piece_records != source.source_piece_records:
            raise _coverage_2.ReachableArrangementTopologyError("Display connector changed exact contacts or source lineage.")
    if any(left.end != right.start for left, right in zip(connector, connector[1:])):
        raise _coverage_2.ReachableArrangementTopologyError("Display connector has disconnected native contacts.")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(13, 7), gridspec_kw={"width_ratios": [1, 1.1]})
    fig.patch.set_facecolor("#faf9f5")
    ring = np.vstack((boundary[:, :2], boundary[0, :2]))
    for ax in axes:
        ax.set_facecolor("#faf9f5")
        ax.fill(ring[:, 0], ring[:, 1], color="#eceae3", zorder=0)
        ax.plot(ring[:, 0], ring[:, 1], color="#8d918e", linewidth=1.2)
        for primitive in primitives:
            points = _display_points(primitive)
            ax.plot(points[:, 0], points[:, 1], color="#babfba", linewidth=1)
        for primitive in connector:
            points = _display_points(primitive)
            ax.plot(points[:, 0], points[:, 1], color=LINE_COLOR if primitive.kind == "line" else ARC_COLOR, linewidth=3)
        ax.set_aspect("equal")
        ax.set_xlabel("World X / mm")
        ax.set_ylabel("World Y / mm")
        ax.grid(color="#dfdfd8", linewidth=0.5, zorder=-1)

    axes[0].set_title("Figure 5 · exact offset of the polygon input", loc="left", fontsize=12, pad=16)
    axes[1].set_title("Native line → arc → line junction", loc="left", fontsize=12, pad=16)
    arc_points = _display_points(connector[1])
    low, high = arc_points.min(axis=0), arc_points.max(axis=0)
    margin = max(float(max(high - low)) * 0.65, float(case.tool_radius.value) * 0.05)
    axes[1].set_xlim(low[0] - margin, high[0] + margin)
    axes[1].set_ylim(low[1] - margin, high[1] + margin)
    for index, point in enumerate((connector[1].start_mm, connector[1].end_mm)):
        axes[1].scatter(*point, s=48, color="#222b35", edgecolor="white", zorder=5)
        axes[1].annotate(
            "line → arc" if index == 0 else "arc → line",
            point,
            xytext=(15, 22 if index == 0 else -30),
            textcoords="offset points",
            fontsize=10,
            arrowprops={"arrowstyle": "-", "color": "#222b35"},
        )
    axes[1].plot([], [], color=LINE_COLOR, linewidth=3, label="Native line")
    axes[1].plot([], [], color=ARC_COLOR, linewidth=3, label="Native circular arc")
    axes[1].legend(loc="best", frameon=False)
    fig.suptitle("Task 5 progress: source-owned offset connectors", x=0.06, y=0.97, ha="left", fontsize=19, fontweight="bold")
    fig.text(0.06, 0.91, f"{len(primitives)} boundary primitives · selected connector: {len(connector)} primitives · every native junction equal", color="#365848")
    fig.text(0.06, 0.055, "Scope: 65-edge polygon projection. Native endpoint and source-lineage checks pass; curves sampled only for display.", fontsize=10)
    fig.text(0.06, 0.025, "Still pending: analytic segment/arc MAT traversal, machining circles, and the full 80° standard toolpath.", fontsize=10, color="#8b432b")
    fig.subplots_adjust(left=0.06, right=0.97, top=0.83, bottom=0.16, wspace=0.25)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=190, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
