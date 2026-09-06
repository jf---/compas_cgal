"""Continuous native coverage of complete Held figure motions, initially uncut."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_motion_coverage import InvalidCoverageMotionError
from benchmarks.held_reference_cases import HeldReferenceCase
from compas_cgal import _coverage_2
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

SAMPLES_PER_TOOL_RADIUS = 4  # Plot inspection only; native residual emptiness owns acceptance.


class IncompleteMotionCoverageError(BenchmarkError):
    """The complete emitted motion leaves material in the declared target."""


def replay_exact_motion_coverage(
    boundary: tuple[Point2[WorldXY], ...],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    tool_radius: ToolRadius,
) -> _coverage_2.ExactRegion2:
    """Return exact target-minus-sweeps; no precleared disk is assumed.

    Full circles use their stored radius, injected separately from the center.
    Connector samples define the emitted straight segments. All constructions,
    unions, and residual predicates belong to native CGAL.
    """
    if not circles or len(transitions) != len(circles) - 1:
        raise InvalidCoverageMotionError("A complete circle path requires exactly one connector per adjacency.")
    for first, second, transition in zip(circles, circles[1:], transitions):
        if len(transition.samples) < 2 or transition.samples[0] != first.contact_point or transition.samples[-1] != second.contact_point:
            raise InvalidCoverageMotionError("Connector endpoints must match the emitted circle contacts.")
    polygon = np.array([(float(p.x), float(p.y), 0.0) for p in boundary], dtype=np.float64)
    target = _coverage_2.ExactRegion2.from_polygon(polygon, [])
    circle_rows = np.array([(float(circle.center.x), float(circle.center.y), float(circle.radius.value)) for circle in circles], dtype=np.float64)
    segments: list[tuple[float, float, float, float]] = []
    disks: list[tuple[float, float]] = []
    for transition in transitions:
        for start, end in zip(transition.samples, transition.samples[1:]):
            if start != end:
                segments.append((float(start.x), float(start.y), float(end.x), float(end.y)))
            else:
                disks.append((float(start.x), float(start.y)))
    return _coverage_2.remaining_material(
        target,
        circle_rows,
        np.array(segments, dtype=np.float64).reshape((-1, 4)),
        np.array(disks, dtype=np.float64).reshape((-1, 2)),
        float(tool_radius.value),
    )


def report_full_motion_coverage(
    case: HeldReferenceCase,
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    output: Path,
) -> _coverage_2.ExactRegion2:
    """Write exact acceptance and a sampled residual plot; return exact residual.

    The declared design polygon is the target. Cutter-inaccessible material is
    not silently excused; reachable-target qualification is a separate claim.
    """
    residual = replay_exact_motion_coverage(case.projection.points, circles, transitions, case.tool_radius)
    complete = residual.is_empty()
    component_count = residual.component_count()
    polygon = np.array([(float(p.x), float(p.y)) for p in case.projection.points])
    pitch = float(case.tool_radius.value) / SAMPLES_PER_TOOL_RADIUS
    remaining = []
    if not complete:
        for y in np.arange(polygon[:, 1].min(), polygon[:, 1].max(), pitch):
            for x in np.arange(polygon[:, 0].min(), polygon[:, 0].max(), pitch):
                if residual.contains(float(x), float(y)):
                    remaining.append((float(x), float(y)))
    report = {
        "case": case.name,
        "target": "declared design polygon",
        "initial_stock": "uncut",
        "coverage_complete": complete,
        "residual_component_count": component_count,
        "circle_count": len(circles),
        "connector_count": len(transitions),
        "grid_pitch_mm": pitch,
        "sampled_remaining_xy_mm": remaining,
        "scope": (
            "Exact native circle annuli and segment capsules. Continuous residual emptiness owns acceptance; "
            "grid only visualizes residuals. No entry clearing assumed; containment and engagement are separate."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    figure = Figure(figsize=(8, 8), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    closed = np.vstack((polygon, polygon[0]))
    axis.fill(closed[:, 0], closed[:, 1], color="#d6ecec")
    axis.plot(closed[:, 0], closed[:, 1], color="#34403f", linewidth=0.8)
    if remaining:
        points = np.array(remaining)
        axis.scatter(points[:, 0], points[:, 1], s=9, marker="s", color="#c0392b", label="Sampled residual stock")
        axis.legend(loc="lower left")
    axis.set(aspect="equal", xlabel="x / mm", ylabel="y / mm")
    figure.suptitle(f"{case.name} · {'COMPLETE' if complete else 'INCOMPLETE'}\n{component_count} exact residual components", fontsize=15)
    figure.text(0.08, 0.025, f"Residual display pitch {pitch:g} mm; sub-grid gaps still fail the exact gate.", fontsize=9)
    figure.subplots_adjust(top=0.88, bottom=0.12)
    figure.savefig(output.with_suffix(".png"), dpi=180)
    return residual


def require_full_motion_coverage(
    case: HeldReferenceCase,
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    output: Path,
) -> None:
    """Report a complete workload and reject any nonempty exact residual."""
    residual = report_full_motion_coverage(case, circles, transitions, output)
    if not residual.is_empty():
        raise IncompleteMotionCoverageError(f"{case.name}: {residual.component_count()} exact residual components; see {output.with_suffix('.json')}")
