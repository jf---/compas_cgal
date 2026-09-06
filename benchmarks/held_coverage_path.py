"""Refine the dense source before any coverage-preserving circle deletion."""

from __future__ import annotations

from benchmarks.held_contour_bound_spacing import select_contour_bound_path
from benchmarks.held_coverage_thinning import thin_covered_contour_path
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_raw_guide import GuideRunId
from compas_cgal import _stock_2
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def build_coverage_preserving_path(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    center_boundary: tuple[Point2[WorldXY], ...],
    design_boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    corner_approaches: bool = False,
) -> tuple[Figure5BoundaryPath, Figure5BoundaryPath, tuple[Radian, ...], _stock_2.Stock2]:
    """Return dense source, thinned motion, bounds, and conservative final stock.

    Geometry reconstruction and bound refinement precede the coverage invariant.
    Subsequent deletion preserves actual connector samples and never changes a
    circle. Source completeness and containment remain separate qualifications.
    """
    _, refined, _ = select_contour_bound_path(
        circles,
        circle_sources,
        center_boundary,
        tool_radius,
        cap,
        corner_approaches=corner_approaches,
        required_sources=frozenset(range(len(circles))),
    )
    dense_sources: list[tuple[GuideRunId, ...]] = [(circle.run_id,) for circle in refined.circles]
    for original, emitted in enumerate(refined.original_indices):
        dense_sources[emitted] = circle_sources[original]
    lengths = _side_lengths(center_boundary)
    transitions = tuple(_transition(center_boundary, lengths, first, second) for first, second in zip(refined.circles, refined.circles[1:]))
    reached = frozenset(run for sources in dense_sources for run in sources)
    dense = Figure5BoundaryPath(refined.circles, transitions, reached)
    _, path, bounds, stock = thin_covered_contour_path(
        dense.circles,
        tuple(dense_sources),
        dense.transitions,
        design_boundary,
        tool_radius,
        cap,
    )
    return dense, path, bounds, stock
