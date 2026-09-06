"""Continuous boundary-ordered motion for the approximate Figure 5 guide."""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Literal

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


class Figure5BoundaryPathError(BenchmarkError):
    """Base failure while constructing a boundary-ordered Figure 5 path."""


class DisconnectedFigure5BoundaryPathError(Figure5BoundaryPathError):
    """The radius-one inward offset is not one connected closed cycle."""


class AmbiguousFigure5BoundaryPathError(Figure5BoundaryPathError):
    """Boundary ownership or the published-start circle is not unique."""


class MissingFigure5GuideRunError(Figure5BoundaryPathError):
    """Selected circles do not preserve every reached approximate guide run."""


RationalPoint = tuple[Fraction, Fraction]


def _rational(point: Point2[WorldXY]) -> RationalPoint:
    return Fraction.from_float(float(point.x)), Fraction.from_float(float(point.y))


def _squared_distance(first: Point2[WorldXY], second: Point2[WorldXY]) -> Fraction:
    first_exact = _rational(first)
    second_exact = _rational(second)
    delta_x = second_exact[0] - first_exact[0]
    delta_y = second_exact[1] - first_exact[1]
    return delta_x * delta_x + delta_y * delta_y


@dataclass(frozen=True)
class Figure5CounterclockwiseCircle:
    """One selected approximate machining circle based at canonical offset q."""

    run_id: GuideRunId
    station_ordinal: GuideRunStationOrdinal
    source_boundary_site: ProjectionBoundarySite
    boundary_site: ProjectionBoundarySite
    center: Point2[WorldXY]
    contact_point: Point2[WorldXY]
    radius: GuideRadius
    clockwise: Literal[False] = field(default=False, init=False)


@dataclass(frozen=True)
class Figure5BoundaryTransition:
    """CCW offset motion, or an explicit zero-progress motion inside a fan."""

    samples: tuple[Point2[WorldXY], ...]
    boundary_progress: Millimetre
    within_fan: bool


@dataclass(frozen=True)
class Figure5BoundaryPath:
    """Published-start phased circles and their continuous closed transitions."""

    circles: tuple[Figure5CounterclockwiseCircle, ...]
    transitions: tuple[Figure5BoundaryTransition, ...]
    reached_run_ids: frozenset[GuideRunId]


def _validate_boundary(component: tuple[Point2[WorldXY], ...]) -> None:
    if len(component) < 3:
        raise DisconnectedFigure5BoundaryPathError("Figure 5 requires one connected radius-one inward-offset cycle.")
    exact = tuple(_rational(point) for point in component)
    if any(exact[index] == exact[(index + 1) % len(exact)] for index in range(len(exact))):
        raise DisconnectedFigure5BoundaryPathError("Figure 5 inward-offset cycle contains a zero-length side.")
    twice_area = sum(point[0] * exact[(index + 1) % len(exact)][1] - exact[(index + 1) % len(exact)][0] * point[1] for index, point in enumerate(exact))
    if twice_area <= 0:
        raise AmbiguousFigure5BoundaryPathError("Figure 5 inward-offset cycle must have an unambiguous CCW orientation.")


def _site_key(site: ProjectionBoundarySite) -> tuple[int, Fraction]:
    return int(site.segment_id), site.parameter


def _canonical_offset_point(
    component: tuple[Point2[WorldXY], ...],
    site: ProjectionBoundarySite,
) -> Point2[WorldXY]:
    side_index = int(site.segment_id)
    if side_index >= len(component):
        raise AmbiguousFigure5BoundaryPathError("A Figure 5 hypothesis names no radius-one offset side.")
    start = _rational(component[side_index])
    end = _rational(component[(side_index + 1) % len(component)])
    point = (
        start[0] + site.parameter * (end[0] - start[0]),
        start[1] + site.parameter * (end[1] - start[1]),
    )
    return Point2[WorldXY].build(float(point[0]), float(point[1]))


def _project_to_offset_side(
    point: RationalPoint,
    component: tuple[Point2[WorldXY], ...],
    side_index: int,
) -> tuple[Fraction, Fraction]:
    start = _rational(component[side_index])
    end = _rational(component[(side_index + 1) % len(component)])
    edge_x = end[0] - start[0]
    edge_y = end[1] - start[1]
    squared_length = edge_x * edge_x + edge_y * edge_y
    parameter = min(
        Fraction(1),
        max(
            Fraction(0),
            ((point[0] - start[0]) * edge_x + (point[1] - start[1]) * edge_y) / squared_length,
        ),
    )
    projected = (
        start[0] + parameter * edge_x,
        start[1] + parameter * edge_y,
    )
    delta_x = point[0] - projected[0]
    delta_y = point[1] - projected[1]
    return parameter, delta_x * delta_x + delta_y * delta_y


def _resolved_site(
    component: tuple[Point2[WorldXY], ...],
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
    evidence_bound: Millimetre,
) -> ProjectionBoundarySite:
    source = hypothesis.projected_boundary_site
    if source.vertex_id is None:
        return source
    current = int(source.segment_id)
    adjacent = ((current - 1) % len(component), current)
    projections = tuple(
        (distance, side, parameter) for side in adjacent for parameter, distance in (_project_to_offset_side(_rational(hypothesis.contact_point), component, side),)
    )
    ordered = tuple(sorted(projections))
    if ordered[0][0] == ordered[1][0]:
        raise AmbiguousFigure5BoundaryPathError("Figure 5 boundary endpoint has no unique adjacent offset side.")
    distance, side, parameter = ordered[0]
    bound = Fraction.from_float(float(evidence_bound))
    if distance > bound * bound:
        raise AmbiguousFigure5BoundaryPathError("Figure 5 endpoint exceeds the recorded boundary evidence bound.")
    resolved = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(side),
        parameter,
        len(component),
    )
    if resolved.vertex_id is not None:
        raise AmbiguousFigure5BoundaryPathError("Figure 5 boundary endpoint ownership remains ambiguous after local resolution.")
    return resolved


def _circle(
    component: tuple[Point2[WorldXY], ...],
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
    evidence_bound: Millimetre,
    resolved_site: ProjectionBoundarySite | None = None,
) -> Figure5CounterclockwiseCircle:
    site = _resolved_site(component, hypothesis, evidence_bound) if resolved_site is None else resolved_site
    contact = _canonical_offset_point(component, site)
    shift_x = float(contact.x) - float(hypothesis.contact_point.x)
    shift_y = float(contact.y) - float(hypothesis.contact_point.y)
    center = Point2[WorldXY].build(
        float(hypothesis.center.x) + shift_x,
        float(hypothesis.center.y) + shift_y,
    )
    return Figure5CounterclockwiseCircle(
        hypothesis.run_id,
        hypothesis.station_ordinal,
        hypothesis.projected_boundary_site,
        site,
        center,
        contact,
        hypothesis.guide_radius,
    )


def _side_lengths(component: tuple[Point2[WorldXY], ...]) -> tuple[float, ...]:
    return tuple(
        math.dist(
            (float(start.x), float(start.y)),
            (float(end.x), float(end.y)),
        )
        for start, end in zip(component, (*component[1:], component[0]))
    )


def _transition(
    component: tuple[Point2[WorldXY], ...],
    lengths: tuple[float, ...],
    start: Figure5CounterclockwiseCircle,
    end: Figure5CounterclockwiseCircle,
) -> Figure5BoundaryTransition:
    if start.boundary_site == end.boundary_site and start is not end:
        if start.contact_point != end.contact_point:
            raise AmbiguousFigure5BoundaryPathError("A Figure 5 fan must share one resolved contact q.")
        return Figure5BoundaryTransition(
            (start.contact_point, end.contact_point),
            Millimetre(0.0),
            True,
        )
    start_side, start_parameter = _site_key(start.boundary_site)
    end_side, end_parameter = _site_key(end.boundary_site)
    wraps = start is end or _site_key(end.boundary_site) < _site_key(start.boundary_site)
    samples = [start.contact_point]
    progress = 0.0
    side = start_side
    if wraps or side != end_side:
        progress += (1.0 - float(start_parameter)) * lengths[side]
        while True:
            vertex = component[(side + 1) % len(component)]
            if vertex != samples[-1]:
                samples.append(vertex)
            side = (side + 1) % len(component)
            if side == end_side:
                break
            progress += lengths[side]
        progress += float(end_parameter) * lengths[end_side]
    else:
        progress = float(end_parameter - start_parameter) * lengths[start_side]
    if end.contact_point != samples[-1] or wraps:
        samples.append(end.contact_point)
    if progress <= 0.0:
        raise AmbiguousFigure5BoundaryPathError("Figure 5 transition requires strictly positive CCW boundary progress.")
    return Figure5BoundaryTransition(tuple(samples), Millimetre(progress), False)


def build_figure5_boundary_path(
    *,
    inward_components: tuple[tuple[Point2[WorldXY], ...], ...],
    side_candidates: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    reached_run_ids: frozenset[GuideRunId],
    published_start_center: Point2[WorldXY],
    published_start_radius: Millimetre,
    tool_radius: ToolRadius,
    boundary_evidence_bound: Millimetre,
    start_evidence_bound: Millimetre,
    offset_sites: dict[ProjectionAdmissibleBoundaryHypothesis, ProjectionBoundarySite] | None = None,
) -> Figure5BoundaryPath:
    """Emit one bounded, continuous CCW traversal over selected hypotheses."""
    if len(inward_components) != 1:
        raise DisconnectedFigure5BoundaryPathError("Figure 5 requires one connected radius-one inward-offset cycle.")
    if not side_candidates:
        raise DisconnectedFigure5BoundaryPathError("Figure 5 requires selected candidates on the inward-offset cycle.")
    if float(published_start_radius) <= 0.0 or float(boundary_evidence_bound) <= 0.0 or float(start_evidence_bound) <= 0.0 or type(tool_radius) is not ToolRadius:
        raise AmbiguousFigure5BoundaryPathError("Published start-circle evidence requires positive radius and bound.")
    radius_gap = abs(Fraction.from_float(float(published_start_radius)) - Fraction.from_float(float(tool_radius.value)))
    if radius_gap > Fraction.from_float(float(start_evidence_bound)):
        raise AmbiguousFigure5BoundaryPathError("Figure 5 published start radius does not match the radius-one tool evidence.")
    component = inward_components[0]
    _validate_boundary(component)
    circles = tuple(_circle(component, candidate, boundary_evidence_bound, None if offset_sites is None else offset_sites[candidate]) for candidate in side_candidates)
    selected_runs = frozenset(circle.run_id for circle in circles)
    if selected_runs != reached_run_ids:
        raise MissingFigure5GuideRunError("Selected circles must cover exactly the reached guide runs.")
    ordered = tuple(
        sorted(
            circles,
            key=lambda circle: (
                *_site_key(circle.boundary_site),
                int(circle.run_id),
                int(circle.station_ordinal),
                float(circle.radius.value),
            ),
        )
    )
    distances = tuple(_squared_distance(circle.contact_point, published_start_center) for circle in ordered)
    minimum_distance = min(distances)
    start_bound = Fraction.from_float(float(start_evidence_bound))
    if minimum_distance > start_bound * start_bound:
        raise AmbiguousFigure5BoundaryPathError("Nearest Figure 5 start q exceeds the published-start evidence bound.")
    nearest = tuple(index for index, distance in enumerate(distances) if distance == minimum_distance)
    seed_sites = {_site_key(ordered[index].boundary_site) for index in nearest}
    if len(seed_sites) != 1:
        raise AmbiguousFigure5BoundaryPathError("Figure 5 published start has no unique nearest boundary site.")
    seed_site = next(iter(seed_sites))
    seed = next(index for index, circle in enumerate(ordered) if _site_key(circle.boundary_site) == seed_site)
    phased = ordered[seed:] + ordered[:seed]
    lengths = _side_lengths(component)
    transitions = tuple(_transition(component, lengths, circle, phased[(index + 1) % len(phased)]) for index, circle in enumerate(phased))
    return Figure5BoundaryPath(phased, transitions, reached_run_ids)
