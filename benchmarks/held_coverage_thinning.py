"""Coverage-preserving deletion from an already refined contour-bound path."""

from __future__ import annotations

import math
from collections import Counter

import numpy as np

from benchmarks.errors import BenchmarkError
from benchmarks.held_contour_bound_spacing import InvalidContourSpacingInputError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_motion_coverage import InvalidCoverageMotionError
from compas_cgal import _stock_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


class UnqualifiedDenseSourceError(BenchmarkError):
    """The source must be refined before coverage-preserving circle deletion."""


def _row(circle: Figure5CounterclockwiseCircle) -> tuple[float, float, float]:
    return float(circle.center.x), float(circle.center.y), float(circle.radius.value)


def _merge(first: Figure5BoundaryTransition, second: Figure5BoundaryTransition) -> Figure5BoundaryTransition:
    if first.samples[-1] != second.samples[0]:
        raise InvalidCoverageMotionError("Adjacent connector chains must meet at the same stored contact.")
    # Keep the actual route, including kinks and complete boundary wraps.
    # Reconstructing a direct connector would invalidate the deletion proof.
    return Figure5BoundaryTransition(
        first.samples + second.samples[1:],
        Millimetre(math.fsum((first.boundary_progress, second.boundary_progress))),
        first.within_fan and second.within_fan,
    )


def thin_covered_contour_path(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    design_boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
) -> tuple[tuple[int, ...], Figure5BoundaryPath, tuple[Radian, ...], _stock_2.Stock2]:
    """Delete circles only when native predicates preserve required swept stock.

    Start after source geometry repair and engagement refinement. No circle is
    moved or inserted here. Original connector chains survive every deletion.
    A conservative native remaining-stock model authorizes deletion only when
    the removed circle contributes no required material beyond the replacement.
    Returns selected indices, motion, bound reports, and conservative remaining
    stock. This preserves source coverage; it does not certify source completeness.
    """
    if not circles or len(circles) != len(circle_sources) or any(not runs or len(set(runs)) != len(runs) for runs in circle_sources):
        raise InvalidContourSpacingInputError("Each dense source circle requires distinct nonempty source attribution.")
    if len(transitions) != len(circles) - 1:
        raise InvalidCoverageMotionError("The dense path requires one connector per circle adjacency.")
    for first, second, transition in zip(circles, circles[1:], transitions):
        if len(transition.samples) < 2 or transition.samples[0] != first.contact_point or transition.samples[-1] != second.contact_point:
            raise InvalidCoverageMotionError("Dense connector endpoints must match their circle contacts.")
    radius = float(tool_radius.value)
    stock = _stock_2.Stock2(np.array([(float(p.x), float(p.y), 0.0) for p in design_boundary], dtype=np.float64), [])
    first_row = _row(circles[0])
    stock.subtract_circle_sweep(*first_row, radius)
    contour = _stock_2.HeldDiskContour2((first_row[0], first_row[1]), first_row[2], radius)
    selected = [0]
    connectors: list[Figure5BoundaryTransition] = []
    bounds: list[Radian] = []
    remaining = Counter(run for runs in circle_sources for run in runs)
    remaining.subtract(circle_sources[0])
    represented = set(circle_sources[0])
    pending = transitions[0] if transitions else None
    for index in range(1, len(circles)):
        if pending is None:
            raise InvalidCoverageMotionError("A retained successor lacks its original connector chain.")
        remaining.subtract(circle_sources[index])
        required = index == len(circles) - 1 or any(run not in represented and remaining[run] == 0 for run in circle_sources[index])
        if not required:
            following = _row(circles[index + 1])
            _, exceeded = contour.engagement_bound((following[0], following[1]), following[2], cap.chord_ratio)
            if not exceeded:
                merged = _merge(pending, transitions[index])
                samples = np.array([(float(p.x), float(p.y)) for p in merged.samples], dtype=np.float64)
                if stock.can_remove_circle(_row(circles[selected[-1]]), _row(circles[index]), following, samples, radius):
                    pending = merged
                    continue
        current = _row(circles[index])
        angle, exceeded = contour.engagement_bound((current[0], current[1]), current[2], cap.chord_ratio)
        if exceeded:
            raise UnqualifiedDenseSourceError(f"Dense source {index} exceeds the contour bound; refine before thinning.")
        selected.append(index)
        connectors.append(pending)
        bounds.append(Radian(angle))
        represented.update(circle_sources[index])
        stock.subtract_circle_sweep(*current, radius)
        for start, end in zip(pending.samples, pending.samples[1:]):
            stock.subtract_capsule_quad(float(start.x), float(start.y), float(end.x), float(end.y), radius)
        contour.append((current[0], current[1]), current[2])
        if index < len(transitions):
            pending = transitions[index]
    expected = {run for runs in circle_sources for run in runs}
    if represented != expected:
        raise InvalidContourSpacingInputError("Coverage thinning lost source attribution.")
    path = Figure5BoundaryPath(tuple(circles[index] for index in selected), tuple(connectors), frozenset(represented))
    return tuple(selected), path, tuple(bounds), stock
