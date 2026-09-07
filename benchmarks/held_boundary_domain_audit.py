"""Approximate source-curve audit; no machinability or coverage acceptance claim.

Run ``pixi run --as-is python -m benchmarks.held_boundary_domain_audit``.
All numeric classifications here are reporting diagnostics, never generator inputs.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import _reference_primitive_points
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

OUTPUT = Path("docs/assets/images/held_boundary_domain_audit.png")
# Display-only cutoff separates visible tangent defects from binary64 residuals.
# Every unrounded join angle is retained in JSON; this is not an acceptance tolerance.
DISPLAY_JOIN_ANGLE_DEG = 0.001


def _tangent(primitive: ReferencePrimitive, *, end: bool) -> tuple[float, float]:
    if isinstance(primitive, ReferenceLine):
        return float(primitive.end.x - primitive.start.x), float(primitive.end.y - primitive.start.y)
    point = primitive.end if end else primitive.start
    x, y = float(point.x - primitive.centre.x), float(point.y - primitive.centre.y)
    direction = 1.0 if float(primitive.sweep) > 0.0 else -1.0
    return -direction * y, direction * x


def _turn(first: tuple[float, float], second: tuple[float, float]) -> float:
    return math.degrees(math.atan2(first[0] * second[1] - first[1] * second[0], first[0] * second[0] + first[1] * second[1]))


def _xy(point: Point2[WorldXY]) -> tuple[float, float]:
    return float(point.x), float(point.y)


def main() -> None:
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 2, figsize=(14, 20))
    reports: list[dict[str, object]] = []
    for row, name in enumerate(CANONICAL_CASE_NAMES):
        case = load_held_reference_case(name)
        primitives = case.boundary.primitives
        tool = float(case.tool_radius.value)
        joins: list[dict[str, object]] = []
        arcs: list[dict[str, object]] = []
        convex_points: list[tuple[float, float]] = []
        convex_radii: list[float] = []
        join_angles: list[float] = []
        for index, primitive in enumerate(primitives):
            previous = primitives[index - 1]
            angle = _turn(_tangent(previous, end=True), _tangent(primitive, end=False))
            join_angles.append(angle)
            joins.append({"index": index, "xy_mm": _xy(primitive.start), "turn_deg": angle, "endpoint_gap_mm": math.dist(_xy(previous.end), _xy(primitive.start))})
            if angle > DISPLAY_JOIN_ANGLE_DEG:
                convex_points.append(_xy(primitive.start))
            points = [_xy(point) for point in _reference_primitive_points(primitive)]
            color = "#087e8b"
            if isinstance(primitive, ReferenceArc):
                radius = math.dist(_xy(primitive.start), _xy(primitive.centre))
                convex = float(primitive.sweep) > 0.0
                if convex:
                    convex_radii.append(radius)
                arcs.append({"index": index, "radius_mm": radius, "convex": convex, "convex_radius_le_tool_reporting": convex and radius <= tool})
                if convex and radius <= tool:
                    color = "#ca3433"
            axes[row, 0].plot(*zip(*points), color=color, linewidth=2.0 if color == "#ca3433" else 0.9)
        if convex_points:
            axes[row, 0].scatter(*zip(*convex_points), color="#d28b12", s=18, zorder=5)
        polygon = [_xy(point) for point in case.projection.points]
        turns = [
            _turn(
                (point[0] - polygon[index - 1][0], point[1] - polygon[index - 1][1]),
                (polygon[(index + 1) % len(polygon)][0] - point[0], polygon[(index + 1) % len(polygon)][1] - point[1]),
            )
            for index, point in enumerate(polygon)
        ]
        projected_convex = [point for point, turn in zip(polygon, turns) if turn > DISPLAY_JOIN_ANGLE_DEG]
        axes[row, 1].plot(*zip(*(polygon + polygon[:1])), color="#59636b", linewidth=0.8)
        if projected_convex:
            axes[row, 1].scatter(*zip(*projected_convex), color="#d28b12", s=5, zorder=5)
        subtool = sum(bool(arc["convex_radius_le_tool_reporting"]) for arc in arcs)
        max_join = max(map(abs, join_angles))
        axes[row, 0].set_title(f"{name}: preserved curves · tool r={tool:g} mm\n{len(convex_points)} convex joins > {DISPLAY_JOIN_ANGLE_DEG:g}° · {subtool} convex arcs ≤ r")
        axes[row, 1].set_title(f"Polygon projection: {len(polygon)} vertices\n{len(projected_convex)} convex turns > {DISPLAY_JOIN_ANGLE_DEG:g}°")
        for ax in axes[row]:
            ax.set_aspect("equal")
            ax.set_xlabel("World X [mm]")
            ax.set_ylabel("World Y [mm]")
        reports.append(
            {
                "case": name,
                "tool_radius_mm": tool,
                "primitive_count": len(primitives),
                "minimum_convex_arc_radius_mm": min(convex_radii) if convex_radii else None,
                "convex_arcs_le_tool_count_reporting": subtool,
                "maximum_absolute_join_turn_deg": max_join,
                "display_convex_join_count": len(convex_points),
                "projection_vertex_count": len(polygon),
                "display_projection_convex_count": len(projected_convex),
                "joins": joins,
                "arcs": arcs,
                "projection_turns_deg": turns,
            }
        )
        print(
            f"{name}: min convex R={min(convex_radii):.9g}, r={tool:g}, subtool={subtool}, "
            f"convex joins={len(convex_points)}, max |join|={max_join:.9g}°, projection convex={len(projected_convex)}",
            flush=True,
        )
    fig.suptitle(
        "Preserved boundary vs polygonization — approximate domain audit\nTeal: source curves · amber: convex joins/vertices · red: convex arc radius ≤ cutter radius", fontsize=15
    )
    fig.text(
        0.5,
        0.012,
        "Reporting only: marker cutoff 0.001°; all raw join angles retained in JSON. No exact smoothness, machinability, or coverage acceptance.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 0.96))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=160)
    plt.close(fig)
    OUTPUT.with_suffix(".json").write_text(
        json.dumps({"status": "approximate_reporting_only", "display_join_angle_deg": DISPLAY_JOIN_ANGLE_DEG, "cases": reports}, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
