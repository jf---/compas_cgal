"""Experimental standard-model spacing in final boundary order.

Retains source-run provenance, not a continuous swept-coverage certificate.
The established lane generator remains a separate, unchanged draft algorithm.
"""

from __future__ import annotations

from collections import Counter

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementRefinement
from benchmarks.held_figure5_engagement_refinement import circle_on_side
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import HypothesisFigure5Path
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


class InvalidBoundarySpacingInputError(BenchmarkError):
    """Circle order and nonempty source-run attribution must correspond."""


def select_boundary_ordered_sources(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
) -> tuple[int, ...]:
    """Select source indices with one predecessor across the supplied order.

    Keep endpoints and the last available representative of an unseen run.
    Otherwise retain the current circle when advancing past it would exceed
    the predecessor cap or lose the model's overlap domain. This is a linear
    lookahead, not bisection: engagement need not be monotone across families.
    Remaining over-cap gaps must pass refinement before motion is emitted.
    """
    if not circles or len(circles) != len(circle_sources) or any(not runs or len(set(runs)) != len(runs) for runs in circle_sources):
        raise InvalidBoundarySpacingInputError("Each ordered circle needs distinct, nonempty source-run attribution.")
    remaining = Counter(run for runs in circle_sources for run in runs)
    represented = set(circle_sources[0])
    selected = [0]
    remaining.subtract(circle_sources[0])
    for index in range(1, len(circles) - 1):
        remaining.subtract(circle_sources[index])
        required = any(run not in represented and remaining[run] == 0 for run in circle_sources[index])
        if not required:
            try:
                engagement = exact_predecessor_engagement(_paper_candidate(circles[selected[-1]]), _paper_candidate(circles[index + 1]), tool_radius)
            except StandardPlacementFragmentationError:
                # An inadmissible jump retains its source for bounded repair;
                # it never becomes an accepted unresolved transition.
                required = True
            else:
                required = not engagement.within(cap)
        if required:
            selected.append(index)
            represented.update(circle_sources[index])
    if len(circles) > 1:
        selected.append(len(circles) - 1)
    return tuple(selected)


def refine_boundary_ordered_path(
    source: HypothesisFigure5Path,
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    part: PartPolygon,
    corner_approaches: bool = False,
    close_seam: bool = False,
) -> tuple[tuple[int, ...], Figure5EngagementRefinement]:
    """Repair source geometry, select in emitted order, then cap every gap.

    Returns original source indices and refined motion. The refinement's
    original_indices refer to the selected subsequence, not the full input.
    Every emitted connector is certified (tiers 1 + 1b) and every emitted
    circle after the entry decided on its circle view of ``part`` by
    ``refine_figure5_engagement``; ``close_seam`` certifies the seam connector
    too, appending circles after the last selected source.
    """
    circles = tuple(circle_on_side(circle, circle.boundary_site, circle.contact_point, circle.radius, boundary) for circle in source.path.circles)
    selected = select_boundary_ordered_sources(circles, source.circle_sources, tool_radius, cap)
    represented = frozenset(run for index in selected for run in source.circle_sources[index])
    if represented != source.reached_run_ids:
        raise InvalidBoundarySpacingInputError("Boundary spacing source attribution differs from the reached-run ledger.")
    picked = tuple(circles[index] for index in selected)
    refined = refine_figure5_engagement(picked, boundary, tool_radius, cap, part=part, corner_approaches=corner_approaches, close_seam=close_seam)
    return selected, refined
