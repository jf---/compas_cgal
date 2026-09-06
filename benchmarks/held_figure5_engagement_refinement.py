"""Bounded engagement refinement proposals on the existing polygon path.

This is an approximate spacing experiment, not machining qualification or MAT
construction. Source circles keep their contacts and family attribution while
their normals and corner radii are repaired. Added circles interpolate source
intervals; original_indices locates repaired source stations. Containment must
be assessed independently.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _validate_boundary
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

MAX_REFINED_CIRCLES = 20000  # Explicit resource ceiling for the Figure 5 experiment.


class Figure5EngagementResolutionError(BenchmarkError):
    """A source interval cannot be refined within the declared budget."""


@dataclass(frozen=True)
class Figure5EngagementRefinement:
    """Refined motion proposals with original guide stations distinguished."""

    circles: tuple[Figure5CounterclockwiseCircle, ...]
    original_indices: tuple[int, ...]


def _corner_radius_limit(boundary: tuple[Point2[WorldXY], ...], side: int, contact: Point2[WorldXY], proposed: GuideRadius) -> GuideRadius:
    radius = float(proposed.value)
    start, end = boundary[side], boundary[(side + 1) % len(boundary)]
    ex, ey = Fraction(float(end.x)) - Fraction(float(start.x)), Fraction(float(end.y)) - Fraction(float(start.y))
    length = math.hypot(float(ex), float(ey))
    for neighbor in ((side - 1) % len(boundary), (side + 1) % len(boundary)):
        a, b = boundary[neighbor], boundary[(neighbor + 1) % len(boundary)]
        fx, fy = Fraction(float(b.x)) - Fraction(float(a.x)), Fraction(float(b.y)) - Fraction(float(a.y))
        turn = fx * ey - fy * ex if neighbor == (side - 1) % len(boundary) else ex * fy - ey * fx
        if turn <= 0:
            continue
        neighbor_length = math.hypot(float(fx), float(fy))
        # A circle centered at q+r*n has signed distance d+r*(n·m)
        # to its neighboring side. Tangency requires r <= d/(1-n·m).
        cosine = float(ex * fx + ey * fy) / (length * neighbor_length)
        denominator = 1 - cosine
        if denominator <= 0:
            raise Figure5EngagementResolutionError("Convex corner normal could not be resolved.")
        qx, qy = Fraction(float(contact.x)) - Fraction(float(a.x)), Fraction(float(contact.y)) - Fraction(float(a.y))
        distance = float(fx * qy - fy * qx) / neighbor_length
        radius = min(radius, distance / denominator)
    if radius <= 0:
        raise Figure5EngagementResolutionError("No positive circle fits at the proposed convex corner contact.")
    return GuideRadius.build(radius)


def _circle_on_side(
    source: Figure5CounterclockwiseCircle,
    site: ProjectionBoundarySite,
    contact: Point2[WorldXY],
    proposed_radius: GuideRadius,
    boundary: tuple[Point2[WorldXY], ...],
) -> Figure5CounterclockwiseCircle:
    side = int(site.segment_id)
    # A bridge must follow the active side's inward normal. Interpolating
    # phase vectors from different sides tilts it through the boundary.
    start, end = boundary[side % len(boundary)], boundary[(side + 1) % len(boundary)]
    edge_x = Fraction(float(end.x)) - Fraction(float(start.x))
    edge_y = Fraction(float(end.y)) - Fraction(float(start.y))
    length = math.hypot(float(edge_x), float(edge_y))
    if length == 0:
        raise Figure5EngagementResolutionError("Interpolation encountered a zero-length boundary side.")
    limited_radius = float(_corner_radius_limit(boundary, side, contact, proposed_radius).value)
    center = Point2[WorldXY].build(float(contact.x) - limited_radius * float(edge_y) / length, float(contact.y) + limited_radius * float(edge_x) / length)
    radius = GuideRadius.build(math.hypot(float(center.x) - float(contact.x), float(center.y) - float(contact.y)))
    # These identifiers retain the left source attribution; original_indices
    # is authoritative for distinguishing actual stations from interpolation.
    return Figure5CounterclockwiseCircle(source.run_id, source.station_ordinal, source.source_boundary_site, site, center, contact, radius)


def _midpoint(
    previous: Figure5CounterclockwiseCircle,
    successor: Figure5CounterclockwiseCircle,
    boundary: tuple[Point2[WorldXY], ...],
) -> Figure5CounterclockwiseCircle:
    previous_side, successor_side = int(previous.boundary_site.segment_id), int(successor.boundary_site.segment_id)
    corner = None
    if successor_side == (previous_side + 1) % len(boundary):
        corner = successor_side
    elif previous.boundary_site.vertex_id is not None:
        corner = previous_side
    elif successor.boundary_site.vertex_id is not None:
        corner = successor_side
    if corner is not None:
        before, vertex, after = boundary[(corner - 1) % len(boundary)], boundary[corner], boundary[(corner + 1) % len(boundary)]
        ex, ey = Fraction(float(vertex.x)) - Fraction(float(before.x)), Fraction(float(vertex.y)) - Fraction(float(before.y))
        fx, fy = Fraction(float(after.x)) - Fraction(float(vertex.x)), Fraction(float(after.y)) - Fraction(float(vertex.y))
        if ex * fy - ey * fx < 0:
            site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(corner), Fraction(0), len(boundary))
            dx = (
                (Fraction(float(previous.center.x)) - Fraction(float(previous.contact_point.x)))
                + (Fraction(float(successor.center.x)) - Fraction(float(successor.contact_point.x)))
            ) / 2
            dy = (
                (Fraction(float(previous.center.y)) - Fraction(float(previous.contact_point.y)))
                + (Fraction(float(successor.center.y)) - Fraction(float(successor.contact_point.y)))
            ) / 2
            center = Point2[WorldXY].build(float(Fraction(float(vertex.x)) + dx), float(Fraction(float(vertex.y)) + dy))
            radius = GuideRadius.build(math.hypot(float(center.x) - float(vertex.x), float(center.y) - float(vertex.y)))
            return Figure5CounterclockwiseCircle(previous.run_id, previous.station_ordinal, previous.source_boundary_site, site, center, vertex, radius)
    lengths = tuple(Fraction(length) for length in _side_lengths(boundary))
    ranges = []
    side = previous_side
    begin = previous.boundary_site.parameter
    while side != successor_side:
        ranges.append((side, begin, Fraction(1)))
        side = (side + 1) % len(boundary)
        begin = Fraction(0)
    ranges.append((side, begin, successor.boundary_site.parameter))
    # Subdivide physical boundary progress: equal side-index progress is not
    # equal distance when a short curved-outline chord meets a long straight.
    total = sum(((end - start) * lengths[index] for index, start, end in ranges), Fraction(0))
    remaining = total / 2
    cumulative = Fraction(0)
    for index, start, end in ranges[:-1]:
        cumulative += (end - start) * lengths[index]
        if cumulative == remaining:
            remaining = total / 3
            break
    for side, begin, end in ranges:
        span = (end - begin) * lengths[side]
        if remaining <= span:
            parameter = begin + remaining / lengths[side]
            break
        remaining -= span
    else:
        raise Figure5EngagementResolutionError("Boundary midpoint exceeded its source interval.")
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, len(boundary))
    contact = _canonical_offset_point(boundary, site)
    return _circle_on_side(previous, site, contact, GuideRadius.build((float(previous.radius.value) + float(successor.radius.value)) / 2), boundary)


def refine_figure5_engagement(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    max_depth: int = 16,
    max_circles: int = MAX_REFINED_CIRCLES,
) -> Figure5EngagementRefinement:
    """Repair corner radii and resolve gaps, preserving source order and q.

    Raises:
        Figure5EngagementResolutionError: Empty input or refinement exhaustion.

    Circle and connector containment are not implied by this model-only step.
    """
    if not circles or len(boundary) < 3 or max_depth < 0 or max_circles < len(circles):
        raise Figure5EngagementResolutionError("Refinement requires circles, a polygon, and sufficient non-negative budgets.")
    _validate_boundary(boundary)
    targets = tuple(_circle_on_side(circle, circle.boundary_site, circle.contact_point, circle.radius, boundary) for circle in circles)
    emitted = [targets[0]]
    original_indices = [0]
    for source_interval, successor in enumerate(targets[1:]):
        pending = [(successor, 0)]
        while pending:
            target, depth = pending.pop()
            predecessor = emitted[-1]
            engagement = maximum_predecessor_engagement(_paper_candidate(predecessor), _paper_candidate(target), tool_radius)
            if float(engagement) <= float(cap.theta):
                if len(emitted) >= max_circles:
                    raise Figure5EngagementResolutionError(f"Circle budget exhausted at source interval {source_interval}.")
                emitted.append(target)
                continue
            if depth >= max_depth:
                raise Figure5EngagementResolutionError(f"Depth exhausted at source interval {source_interval}: {math.degrees(float(engagement)):.6f} degrees.")
            midpoint = _midpoint(predecessor, target, boundary)
            pending.append((target, depth + 1))
            pending.append((midpoint, depth + 1))
        original_indices.append(len(emitted) - 1)
    return Figure5EngagementRefinement(tuple(emitted), tuple(original_indices))
