"""Conservative contour-bound drafts on the existing ordered circle proposals.

CGAL's whole-orbit leading-rim bound owns every cap decision. This algorithm
does not assume the ordered-MATHSM contour shortcut or certify physical stock,
initial immersion, connector engagement, containment, or continuous coverage.
"""

from __future__ import annotations

import math
from collections import Counter

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _validate_boundary
from benchmarks.held_figure5_engagement_refinement import MAX_REFINED_CIRCLES
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementRefinement
from benchmarks.held_figure5_engagement_refinement import _circle_on_side
from benchmarks.held_figure5_engagement_refinement import _concave_corner_bridges
from benchmarks.held_figure5_engagement_refinement import _midpoint
from benchmarks.held_figure5_raw_guide import GuideRunId
from compas_cgal import _stock_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


class InvalidContourSpacingInputError(BenchmarkError):
    """Circle order, boundary, source attribution, or resource budget is invalid."""


class ContourSpacingResolutionError(BenchmarkError):
    """A required target cannot be reached within the declared refinement budget."""


def _xy(circle: Figure5CounterclockwiseCircle) -> tuple[float, float]:
    return float(circle.center.x), float(circle.center.y)


def select_contour_bound_path(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    corner_approaches: bool = False,
    required_sources: frozenset[int] = frozenset(),
    max_depth: int = 16,
    max_circles: int = MAX_REFINED_CIRCLES,
) -> tuple[tuple[int, ...], Figure5EngagementRefinement, tuple[Radian, ...]]:
    """Select and refine in emission order against all accepted outer disks.

    Returns selected source indices, refined circles/source positions, and one
    reported bound per successor. The first circle is the declared cleared
    seed; there is no fabricated engagement measurement for its immersion.
    """
    if not circles or len(circles) != len(circle_sources) or any(not runs or len(set(runs)) != len(runs) for runs in circle_sources):
        raise InvalidContourSpacingInputError("Each circle requires distinct nonempty source-run attribution.")
    if max_depth < 0 or max_circles < 1:
        raise InvalidContourSpacingInputError("Contour spacing requires nonnegative depth and a positive circle budget.")
    if any(index < 0 or index >= len(circles) for index in required_sources):
        raise InvalidContourSpacingInputError("Required source indices must belong to the supplied circle sequence.")
    _validate_boundary(boundary)
    targets = tuple(_circle_on_side(circle, circle.boundary_site, circle.contact_point, circle.radius, boundary) for circle in circles)
    emitted = [targets[0]]
    selected = [0]
    original_indices = [0]
    bounds: list[Radian] = []
    contour = _stock_2.HeldDiskContour2(_xy(targets[0]), float(targets[0].radius.value), float(tool_radius.value))
    remaining = Counter(run for runs in circle_sources for run in runs)
    represented = set(circle_sources[0])
    remaining.subtract(circle_sources[0])

    def bound(target: Figure5CounterclockwiseCircle) -> tuple[float, bool]:
        return contour.engagement_bound(_xy(target), float(target.radius.value), cap.chord_ratio)

    for index in range(1, len(targets)):
        remaining.subtract(circle_sources[index])
        required = index in required_sources or index == len(targets) - 1 or any(run not in represented and remaining[run] == 0 for run in circle_sources[index])
        if not required and not bound(targets[index + 1])[1]:
            continue

        successor = targets[index]
        bridges = _concave_corner_bridges(emitted[-1], successor, boundary) if corner_approaches else ()
        pending = [(target, 0) for target in reversed((*bridges, successor))]
        while pending:
            target, depth = pending.pop()
            angle, exceeded = bound(target)
            if not exceeded:
                if len(emitted) >= max_circles:
                    raise ContourSpacingResolutionError(f"Circle budget exhausted at source {index}.")
                contour.append(_xy(target), float(target.radius.value))
                emitted.append(target)
                bounds.append(Radian(angle))
                continue
            if depth >= max_depth:
                raise ContourSpacingResolutionError(f"Depth exhausted at source {index}: bound {math.degrees(angle):.6f} degrees.")
            midpoint = _midpoint(emitted[-1], target, boundary, corner_approaches=corner_approaches)
            pending.append((target, depth + 1))
            pending.append((midpoint, depth + 1))
        selected.append(index)
        original_indices.append(len(emitted) - 1)
        represented.update(circle_sources[index])
    expected = {run for runs in circle_sources for run in runs}
    if represented != expected:
        raise InvalidContourSpacingInputError("Contour selection lost a required source run.")
    return tuple(selected), Figure5EngagementRefinement(tuple(emitted), tuple(original_indices)), tuple(bounds)
