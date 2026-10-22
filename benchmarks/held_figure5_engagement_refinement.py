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
from typing import Self

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_boundary_path import _validate_boundary
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal import _circle_geometry_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

MAX_REFINED_CIRCLES = 20000  # Explicit resource ceiling for the Figure 5 experiment.


class Figure5EngagementResolutionError(BenchmarkError):
    """A source interval cannot be refined within the declared budget."""


@dataclass(frozen=True)
class Figure5EngagementRefinement:
    """Refined motion proposals with original guide stations distinguished."""

    circles: tuple[Figure5CounterclockwiseCircle, ...]
    original_indices: tuple[int, ...]


FIRST_ENCLOSURE_BITS = 64  # Initial binary digits of each certified square-root enclosure; doubled until a decision separates.


def _rational_square_root(value: Fraction) -> Fraction | None:
    """The exact rational square root of ``value >= 0``, or None when it is irrational."""
    numerator, denominator = math.isqrt(value.numerator), math.isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Fraction(numerator, denominator)
    return None


def _root_enclosure(radicand: Fraction, bits: int) -> tuple[Fraction, Fraction]:
    """Rational bounds ``lo <= sqrt(radicand) <= hi`` with ``hi - lo <= 2**-bits / denominator``."""
    scaled = radicand.numerator * radicand.denominator << (2 * bits)
    root = math.isqrt(scaled)
    scale = radicand.denominator << bits
    return Fraction(root, scale), Fraction(root + (root * root != scaled), scale)


@dataclass(frozen=True)
class _BoundaryArclength:
    """Exact arclength ``sum(c * sqrt(r))`` along the polygonal Ω, ``r > 0`` rational.

    Folded so no two radicands differ by a rational-square factor; rational
    terms carry radicand 1. Square roots of positive rationals pairwise in
    distinct square classes are linearly independent over Q (Besicovitch),
    so the sum is zero iff it has no terms, and a nonzero sum's sign is
    settled by certified ``isqrt`` enclosures refined until they exclude 0.
    """

    terms: tuple[tuple[Fraction, Fraction], ...]

    @classmethod
    def zero(cls) -> Self:
        return cls(())

    @classmethod
    def root(cls, radicand: Fraction, coefficient: Fraction) -> Self:
        """``coefficient * sqrt(radicand)`` for ``radicand >= 0``."""
        return cls.zero()._folded(((radicand, coefficient),))

    def _folded(self, extra: tuple[tuple[Fraction, Fraction], ...]) -> Self:
        classes = dict(self.terms)
        for radicand, coefficient in extra:
            if radicand == 0 or coefficient == 0:
                continue
            square_root = _rational_square_root(radicand)
            if square_root is not None:
                radicand, coefficient = Fraction(1), coefficient * square_root
            for representative in classes:
                ratio = _rational_square_root(radicand * representative)
                if ratio is not None:
                    # sqrt(radicand) == sqrt(radicand * representative) / representative * sqrt(representative)
                    classes[representative] += coefficient * ratio / representative
                    break
            else:
                classes[radicand] = coefficient
        return type(self)(tuple((radicand, coefficient) for radicand, coefficient in classes.items() if coefficient != 0))

    def __add__(self, other: _BoundaryArclength) -> _BoundaryArclength:
        return self._folded(other.terms)

    def __sub__(self, other: _BoundaryArclength) -> _BoundaryArclength:
        return self._folded(tuple((radicand, -coefficient) for radicand, coefficient in other.terms))

    def scaled(self, factor: Fraction) -> _BoundaryArclength:
        return self.zero()._folded(tuple((radicand, coefficient * factor) for radicand, coefficient in self.terms))

    def per_root(self, radicand: Fraction) -> _BoundaryArclength:
        """``self / sqrt(radicand)`` for ``radicand > 0``: ``c*sqrt(r)/sqrt(s) == (c/s)*sqrt(r*s)``."""
        return self.zero()._folded(tuple((term * radicand, coefficient / radicand) for term, coefficient in self.terms))

    def rational(self) -> Fraction | None:
        """The exact value when it is rational (only the radicand-1 class), else None."""
        if not self.terms:
            return Fraction(0)
        if len(self.terms) == 1 and self.terms[0][0] == 1:
            return self.terms[0][1]
        return None

    def enclosure(self, bits: int) -> tuple[Fraction, Fraction]:
        low = high = Fraction(0)
        for radicand, coefficient in self.terms:
            root_low, root_high = _root_enclosure(radicand, bits)
            low += coefficient * (root_low if coefficient > 0 else root_high)
            high += coefficient * (root_high if coefficient > 0 else root_low)
        return low, high

    def sign(self) -> int:
        exact = self.rational()
        if exact is not None:
            return (exact > 0) - (exact < 0)
        bits = FIRST_ENCLOSURE_BITS
        while True:
            low, high = self.enclosure(bits)
            if low > 0:
                return 1
            if high < 0:
                return -1
            bits *= 2

    def nearest_binary64(self) -> float:
        """The correctly rounded binary64 value: enclosures refine until both ends round alike."""
        exact = self.rational()
        if exact is not None:
            return float(exact)
        bits = FIRST_ENCLOSURE_BITS
        while True:
            low, high = self.enclosure(bits)
            # Rounding boundaries are rational, so an irrational value is never on one.
            if float(low) == float(high):
                return float(low)
            bits *= 2


def _side_parameter(position: _BoundaryArclength, begin: Fraction, end: Fraction) -> Fraction:
    """Representation of an exact side parameter in ``[begin, end]``.

    The binary64 nearest ``position`` (the representation the float bookkeeping
    stored), unless that rounding leaves the exactly decided range; then the
    exact rational, or for an irrational position an enclosure end strictly
    inside the range (it is then strictly inside, since ``begin``/``end`` are
    rational).
    """
    nearest = Fraction(position.nearest_binary64())
    if begin <= nearest <= end:
        return nearest
    exact = position.rational()
    if exact is not None:
        return exact
    bits = FIRST_ENCLOSURE_BITS
    while True:
        low, high = position.enclosure(bits)
        if begin < low and high < end:
            return low
        bits *= 2


def _turn(a: Point2[WorldXY], b: Point2[WorldXY], c: Point2[WorldXY]) -> int:
    return int(_circle_geometry_2.orientation((float(a.x), float(a.y)), (float(b.x), float(b.y)), (float(c.x), float(c.y))))


def _squared_side_length(boundary: tuple[Point2[WorldXY], ...], side: int) -> Fraction:
    start, end = boundary[side], boundary[(side + 1) % len(boundary)]
    edge_x = Fraction(float(end.x)) - Fraction(float(start.x))
    edge_y = Fraction(float(end.y)) - Fraction(float(start.y))
    return edge_x * edge_x + edge_y * edge_y


def _constructed_radius(center: Point2[WorldXY], contact: Point2[WorldXY], cause: str) -> GuideRadius:
    """Radius of the circle as constructed, refused when its centre collapsed onto the contact.

    A representability guard, not a geometric decision: whether a positive
    circle exists is decided exactly upstream (``_corner_radius_limit``).
    Centres are binary64 points, so a positive radius below the contact
    coordinate's ulp is absorbed by the construction; exact equality of the
    constructed points detects that. Reachable only when the exact midpoint
    lies within about an ulp of a convex vertex (``cause`` names the caller's
    construction).
    """
    if center == contact:
        raise Figure5EngagementResolutionError(f"No representable positive circle touches contact ({float(contact.x)!r}, {float(contact.y)!r}): {cause}.")
    return GuideRadius.build(math.hypot(float(center.x) - float(contact.x), float(center.y) - float(contact.y)))


def _corner_radius_limit(boundary: tuple[Point2[WorldXY], ...], side: int, contact: Point2[WorldXY], proposed: GuideRadius) -> GuideRadius:
    radius = float(proposed.value)
    start, end = boundary[side], boundary[(side + 1) % len(boundary)]
    ex, ey = float(end.x) - float(start.x), float(end.y) - float(start.y)
    length = math.hypot(float(ex), float(ey))
    for neighbor in ((side - 1) % len(boundary), (side + 1) % len(boundary)):
        a, b = boundary[neighbor], boundary[(neighbor + 1) % len(boundary)]
        fx, fy = float(b.x) - float(a.x), float(b.y) - float(a.y)
        turn = _turn(a, start, end) if neighbor == (side - 1) % len(boundary) else _turn(start, end, b)
        if turn <= 0:
            continue
        # Exact predicate: inside a convex wedge a positive circle tangent at
        # the contact exists iff the contact is strictly left of the
        # neighbouring side, i.e. it is not the apex.
        if _turn(a, b, contact) <= 0:
            raise Figure5EngagementResolutionError("No positive circle fits at the convex corner apex contact.")
        neighbor_length = math.hypot(float(fx), float(fy))
        # A circle centered at q+r*n has signed distance d+r*(n·m)
        # to its neighboring side. Tangency requires r <= d/(1-n·m).
        cosine = float(ex * fx + ey * fy) / (length * neighbor_length)
        denominator = 1 - cosine
        if denominator <= 0:
            raise Figure5EngagementResolutionError("Convex corner normal could not be resolved.")
        qx, qy = float(contact.x) - float(a.x), float(contact.y) - float(a.y)
        distance = float(fx * qy - fy * qx) / neighbor_length
        radius = min(radius, distance / denominator)
    if radius <= 0:
        # Construction guard: the exact predicate above proved the radius
        # positive; binary64 evaluation lost it below representability.
        raise Figure5EngagementResolutionError("Corner-limited radius is positive but not representable at this contact.")
    return GuideRadius.build(radius)


def circle_on_side(
    source: Figure5CounterclockwiseCircle,
    site: ProjectionBoundarySite,
    contact: Point2[WorldXY],
    proposed_radius: GuideRadius,
    boundary: tuple[Point2[WorldXY], ...],
) -> Figure5CounterclockwiseCircle:
    return _circle_with_side_normal(source, site, contact, proposed_radius, boundary, int(site.segment_id))


def _circle_with_side_normal(
    source: Figure5CounterclockwiseCircle,
    site: ProjectionBoundarySite,
    contact: Point2[WorldXY],
    proposed_radius: GuideRadius,
    boundary: tuple[Point2[WorldXY], ...],
    side: int,
) -> Figure5CounterclockwiseCircle:
    """Construct a circle using one incident side at a canonical contact.

    A concave vertex has incoming and outgoing normals at the same site;
    contact identity alone does not select the phase of its machining circle.
    """
    # A bridge must follow the active side's inward normal. Interpolating
    # phase vectors from different sides tilts it through the boundary.
    start, end = boundary[side % len(boundary)], boundary[(side + 1) % len(boundary)]
    edge_x = float(end.x) - float(start.x)
    edge_y = float(end.y) - float(start.y)
    length = math.hypot(float(edge_x), float(edge_y))
    if length == 0:
        raise Figure5EngagementResolutionError("Interpolation encountered a zero-length boundary side.")
    limited_radius = float(_corner_radius_limit(boundary, side, contact, proposed_radius).value)
    center = Point2[WorldXY].build(float(contact.x) - limited_radius * float(edge_y) / length, float(contact.y) + limited_radius * float(edge_x) / length)
    radius = _constructed_radius(center, contact, "corner-limited radius below binary64 resolution")
    # These identifiers retain the left source attribution; original_indices
    # is authoritative for distinguishing actual stations from interpolation.
    return Figure5CounterclockwiseCircle(source.run_id, source.station_ordinal, source.source_boundary_site, site, center, contact, radius)


def _midpoint(
    previous: Figure5CounterclockwiseCircle,
    successor: Figure5CounterclockwiseCircle,
    boundary: tuple[Point2[WorldXY], ...],
    *,
    corner_approaches: bool = False,
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
        shared_contact = previous.contact_point == vertex and successor.contact_point == vertex
        if _turn(before, vertex, after) < 0 and (not corner_approaches or shared_contact):
            site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(corner), Fraction(0), len(boundary))
            dx = ((float(previous.center.x) - float(previous.contact_point.x)) + (float(successor.center.x) - float(successor.contact_point.x))) / 2
            dy = ((float(previous.center.y) - float(previous.contact_point.y)) + (float(successor.center.y) - float(successor.contact_point.y))) / 2
            center = Point2[WorldXY].build(float(vertex.x) + dx, float(vertex.y) + dy)
            radius = _constructed_radius(center, vertex, "concave corner normals cancel")
            return Figure5CounterclockwiseCircle(previous.run_id, previous.station_ordinal, previous.source_boundary_site, site, center, vertex, radius)
    # Subdivide physical boundary progress: equal side-index progress is not
    # equal distance when a short curved-outline chord meets a long straight.
    # Every decision (half on a junction? which side?) is exact on the true
    # arclength sum(w_i * |side_i|) over the doubles' rationals; binary64
    # bookkeeping put halves that are exactly a convex vertex one rounding
    # past it, where the corner-limited radius collapses.
    ranges: list[tuple[int, Fraction, Fraction]] = []
    side = previous_side
    begin = previous.boundary_site.parameter
    while side != successor_side:
        ranges.append((side, begin, Fraction(1)))
        side = (side + 1) % len(boundary)
        begin = Fraction(0)
    ranges.append((side, begin, successor.boundary_site.parameter))
    squared_lengths = {index: _squared_side_length(boundary, index) for index, _, _ in ranges}
    cumulative = [_BoundaryArclength.zero()]
    for index, start, end in ranges:
        cumulative.append(cumulative[-1] + _BoundaryArclength.root(squared_lengths[index], end - start))
    total = cumulative[-1]
    target = total.scaled(Fraction(1, 2))
    if any((reached - target).sign() == 0 for reached in cumulative[1:-1]):
        target = total.scaled(Fraction(1, 3))
    for (side, begin, end), entered, reached in zip(ranges, cumulative, cumulative[1:]):
        if (reached - target).sign() >= 0:
            break
    else:
        raise Figure5EngagementResolutionError("Boundary midpoint exceeded its source interval.")
    position = (target - entered).per_root(squared_lengths[side]) + _BoundaryArclength.root(Fraction(1), begin)
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), _side_parameter(position, begin, end), len(boundary))
    contact = _canonical_offset_point(boundary, site)
    return circle_on_side(previous, site, contact, GuideRadius.build((float(previous.radius.value) + float(successor.radius.value)) / 2), boundary)


def _concave_corner_bridges(
    previous: Figure5CounterclockwiseCircle,
    successor: Figure5CounterclockwiseCircle,
    boundary: tuple[Point2[WorldXY], ...],
) -> tuple[Figure5CounterclockwiseCircle, ...]:
    """Partition crossed concave corners into incoming and outgoing phases."""
    side = int(previous.boundary_site.segment_id)
    target_side = int(successor.boundary_site.segment_id)
    bridges = []
    radius = GuideRadius.build((float(previous.radius.value) + float(successor.radius.value)) / 2)
    while side != target_side:
        corner = (side + 1) % len(boundary)
        vertex = boundary[corner]
        if _turn(boundary[side], vertex, boundary[(corner + 1) % len(boundary)]) < 0:
            site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(corner), Fraction(0), len(boundary))
            bridges.append(_circle_with_side_normal(previous, site, vertex, radius, boundary, side))
            bridges.append(_circle_with_side_normal(previous, site, vertex, radius, boundary, corner))
        side = corner
    return tuple(bridges)


def refine_figure5_engagement(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    max_depth: int = 16,
    max_circles: int = MAX_REFINED_CIRCLES,
    corner_approaches: bool = False,
) -> Figure5EngagementRefinement:
    """Repair corner radii and resolve gaps, preserving source order and q.

    Raises:
        Figure5EngagementResolutionError: Empty input, refinement exhaustion,
            or a midpoint whose circle collapses onto its contact (a convex
            corner apex admits no positive circle).

    Circle and connector containment are not implied by this model-only step.
    """
    if not circles or len(boundary) < 3 or max_depth < 0 or max_circles < len(circles):
        raise Figure5EngagementResolutionError("Refinement requires circles, a polygon, and sufficient non-negative budgets.")
    _validate_boundary(boundary)
    targets = tuple(circle_on_side(circle, circle.boundary_site, circle.contact_point, circle.radius, boundary) for circle in circles)
    emitted = [targets[0]]
    original_indices = [0]
    for source_interval, successor in enumerate(targets[1:]):
        bridges = _concave_corner_bridges(emitted[-1], successor, boundary) if corner_approaches else ()
        pending = [(target, 0) for target in reversed((*bridges, successor))]
        while pending:
            target, depth = pending.pop()
            predecessor = emitted[-1]
            engagement = exact_predecessor_engagement(_paper_candidate(predecessor), _paper_candidate(target), tool_radius)
            if engagement.within(cap):
                if len(emitted) >= max_circles:
                    raise Figure5EngagementResolutionError(f"Circle budget exhausted at source interval {source_interval}.")
                emitted.append(target)
                continue
            if depth >= max_depth:
                raise Figure5EngagementResolutionError(f"Depth exhausted at source interval {source_interval}: {math.degrees(float(engagement.reported_angle())):.6f} degrees.")
            midpoint = _midpoint(predecessor, target, boundary, corner_approaches=corner_approaches)
            pending.append((target, depth + 1))
            pending.append((midpoint, depth + 1))
        original_indices.append(len(emitted) - 1)
    return Figure5EngagementRefinement(tuple(emitted), tuple(original_indices))
