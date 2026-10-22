"""Circle construction of the engagement refinement: side placement, corner limits, midpoints, concave-corner bridges."""

from __future__ import annotations

import math
from fractions import Fraction

from benchmarks.held_boundary_arclength import FIRST_ENCLOSURE_BITS
from benchmarks.held_boundary_arclength import BoundaryArclength
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_engagement_errors import Figure5EngagementResolutionError
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal import _circle_geometry_2
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


def _side_parameter(position: BoundaryArclength, begin: Fraction, end: Fraction) -> Fraction:
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
    """Return the circle `refine_figure5_engagement` emits for `source` at `site` on its own side.

    The circle touches the boundary side `site.segment_id` at `contact`, with its centre on that
    side's inward normal and `proposed_radius` capped by the convex-corner limit of the adjacent
    sides. Callers that must reproduce the refinement's emitted form call this, never a copy.

    Args:
        source: Circle whose run and station attribution the result keeps.
        site: Boundary site; its side selects the inward normal.
        contact: Tangency point on that side.
        proposed_radius: Radius before the corner limit.
        boundary: Counterclockwise boundary polygon.

    Returns:
        The emitted counterclockwise circle tangent at `contact`.

    Raises:
        Figure5EngagementResolutionError: The side has zero length, `contact` is a convex corner
            apex, or the corner-limited radius is not representable at `contact`.
    """
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


def midpoint(
    previous: Figure5CounterclockwiseCircle,
    successor: Figure5CounterclockwiseCircle,
    boundary: tuple[Point2[WorldXY], ...],
    *,
    corner_approaches: bool = False,
) -> Figure5CounterclockwiseCircle:
    """Construct the refinement circle halfway between two circles along the boundary.

    When the interval meets a concave corner, the circle touches that vertex and its centre is
    offset by the mean of the two centre-to-contact vectors. Otherwise its contact bisects the
    exact boundary arclength from `previous` to `successor`, or trisects it when the half lands
    exactly on a side junction. That circle is `circle_on_side` with the mean of the two radii.

    Args:
        previous: Circle the interval starts at; the result keeps its run and station attribution.
        successor: Circle the interval ends at, counterclockwise from `previous`.
        boundary: Counterclockwise boundary polygon.
        corner_approaches: Take the concave-corner branch only when both circles touch the vertex.

    Returns:
        The interpolated counterclockwise circle.

    Raises:
        Figure5EngagementResolutionError: The concave-corner normals cancel, the arclength target
            leaves the interval, or `circle_on_side` refuses the contact (zero-length side, convex
            corner apex, or a corner-limited radius that is not representable).
    """
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
    # The seam runs from the last circle back to circles[0]: on one side, a
    # successor at or before ``previous`` (or ``previous`` itself, a lone
    # circle's lap) is reached once round the loop, as ``_transition`` walks it.
    if side == successor_side and (previous is successor or successor.boundary_site.parameter < begin):
        ranges.append((side, begin, Fraction(1)))
        side = (side + 1) % len(boundary)
        begin = Fraction(0)
    while side != successor_side:
        ranges.append((side, begin, Fraction(1)))
        side = (side + 1) % len(boundary)
        begin = Fraction(0)
    ranges.append((side, begin, successor.boundary_site.parameter))
    squared_lengths = {index: _squared_side_length(boundary, index) for index, _, _ in ranges}
    cumulative = [BoundaryArclength.zero()]
    for index, start, end in ranges:
        cumulative.append(cumulative[-1] + BoundaryArclength.root(squared_lengths[index], end - start))
    total = cumulative[-1]
    target = total.scaled(Fraction(1, 2))
    if any((reached - target).sign() == 0 for reached in cumulative[1:-1]):
        target = total.scaled(Fraction(1, 3))
    for (side, begin, end), entered, reached in zip(ranges, cumulative, cumulative[1:]):
        if (reached - target).sign() >= 0:
            break
    else:
        raise Figure5EngagementResolutionError("Boundary midpoint exceeded its source interval.")
    position = (target - entered).per_root(squared_lengths[side]) + BoundaryArclength.root(Fraction(1), begin)
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), _side_parameter(position, begin, end), len(boundary))
    contact = _canonical_offset_point(boundary, site)
    return circle_on_side(previous, site, contact, GuideRadius.build((float(previous.radius.value) + float(successor.radius.value)) / 2), boundary)


def concave_corner_bridges(
    previous: Figure5CounterclockwiseCircle,
    successor: Figure5CounterclockwiseCircle,
    boundary: tuple[Point2[WorldXY], ...],
) -> tuple[Figure5CounterclockwiseCircle, ...]:
    """Partition crossed concave corners into incoming and outgoing phases.

    Walking the sides counterclockwise from `previous` to `successor`, each concave vertex crossed
    gets two circles touching it. The first follows the incoming side's inward normal and the
    second follows the outgoing side's. Both use the mean of the two radii, capped by the
    convex-corner limit.

    Args:
        previous: Circle the interval starts at; the bridges keep its run and station attribution.
        successor: Circle the interval ends at, counterclockwise from `previous`.
        boundary: Counterclockwise boundary polygon.

    Returns:
        The bridge circles in emission order, two per concave corner crossed. The tuple is empty
        when the interval crosses no concave corner.

    Raises:
        Figure5EngagementResolutionError: A bridge circle is refused at its vertex (zero-length
            side, convex corner apex, or a corner-limited radius that is not representable).
    """
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
