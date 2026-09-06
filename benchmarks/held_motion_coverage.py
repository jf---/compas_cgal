"""Sample native remaining stock after full circles and boundary connectors."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_stock_replay import replay_circle_stock
from compas_cgal import _stock_2
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

SAMPLES_PER_TOOL_RADIUS = 4  # Quarter-radius inspection grid; not a continuous proof.


class InvalidCoverageMotionError(ValueError):
    """A circle sequence lacks its connecting motions."""


def replay_circle_connector_stock(
    boundary: tuple[Point2[WorldXY], ...],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    tool_radius: ToolRadius,
) -> _stock_2.Stock2:
    """Return final stock; exact annuli plus certified under-cover capsules.

    Final subtraction commutes, so circles precede connectors in this replay.
    This function does not claim to measure engagement during the motion.
    """
    if not circles or len(transitions) != len(circles) - 1:
        raise InvalidCoverageMotionError("Coverage replay requires a nonempty circle sequence and its connectors.")
    stock = replay_circle_stock(boundary, circles, tool_radius)
    for transition in transitions:
        if len(transition.samples) < 2:
            raise InvalidCoverageMotionError("A connector requires at least its two endpoints.")
        for first, second in zip(transition.samples, transition.samples[1:]):
            if first == second:
                continue
            stock.subtract_capsule_quad(float(first.x), float(first.y), float(second.x), float(second.y), float(tool_radius.value))
    return stock


def compare_sampled_motion_coverage(
    case: HeldReferenceCase,
    baseline: tuple[Figure5CounterclockwiseCircle, ...],
    baseline_transitions: tuple[Figure5BoundaryTransition, ...],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    output: Path,
) -> tuple[int, int]:
    """Write grid comparison and return lost/newly removed sample counts.

    Initially uncut stock is explicit: no seed clearing or immersion is invented.
    Counts cannot establish continuous pocket coverage or contour containment.
    """
    before = replay_circle_connector_stock(case.projection.points, baseline, baseline_transitions, case.tool_radius)
    after = replay_circle_connector_stock(case.projection.points, circles, transitions, case.tool_radius)
    uncut = replay_circle_stock(case.projection.points, (), case.tool_radius)
    boundary = np.array([(float(point.x), float(point.y)) for point in case.projection.points])
    pitch = float(case.tool_radius.value) / SAMPLES_PER_TOOL_RADIUS
    xs = np.arange(boundary[:, 0].min(), boundary[:, 0].max(), pitch)
    ys = np.arange(boundary[:, 1].min(), boundary[:, 1].max(), pitch)
    inside = baseline_cleared = draft_cleared = 0
    lost: list[tuple[float, float]] = []
    gained: list[tuple[float, float]] = []
    for y in ys:
        for x in xs:
            if not uncut.contains(float(x), float(y)):
                continue
            inside += 1
            old = not before.contains(float(x), float(y))
            new = not after.contains(float(x), float(y))
            baseline_cleared += int(old)
            draft_cleared += int(new)
            if old and not new:
                lost.append((float(x), float(y)))
            elif new and not old:
                gained.append((float(x), float(y)))
    report = {
        "case": case.name,
        "grid_pitch_mm": pitch,
        "grid_origin_xy_mm": [float(xs[0]), float(ys[0])],
        "inside_sample_count": inside,
        "baseline_removed_count": baseline_cleared,
        "draft_removed_count": draft_cleared,
        "lost_count": len(lost),
        "newly_removed_count": len(gained),
        "lost_xy_mm": lost,
        "newly_removed_xy_mm": gained,
        "scope": "Sampled initially uncut stock; exact annular sweeps and certified under-cover connector capsules; no seed clearing or continuous proof.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return len(lost), len(gained)
