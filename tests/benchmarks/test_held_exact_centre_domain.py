"""Contract tests for the exact Held centre domain Ω* (docs/held_exact_centre_domain.md).

Every fixture is Pythagorean by construction (integers, 3-4-5 triples, r = 1),
so expected values are exact rationals or one-root numbers. The tests named in
the design's "Failing tests to write first" are red until the native
arrangement (``src/held_centre_domain_2``) exists: each one first checks the
pure-rational parts that are implemented in Python (supports, local vertices),
then calls :meth:`CentreDomain.build`, which raises
:class:`NativeCentreDomainUnavailableError`. The "old" assertions document the
defect of the float straight-skeleton domain and are green.
"""

from __future__ import annotations

import hashlib
import math
import time
from fractions import Fraction
from functools import cache

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_centre_domain_supports import CentreDomainSupports
from benchmarks.held_centre_domain_supports import InvalidCentreDomainToolRadiusError
from benchmarks.held_centre_domain_supports import JunctionVertexKind
from benchmarks.held_centre_domain_supports import NoSupportError
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_declared_boundary_certificate import DeclaredBoundaryCertificate
from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import CentreDomainPiece
from benchmarks.held_exact_centre_domain import CentreDomainSite
from benchmarks.held_exact_centre_domain import NeckSplitsCentreDomainError
from benchmarks.held_exact_centre_domain import PointLocation
from benchmarks.held_exact_centre_domain import UnreachableCause
from benchmarks.held_one_root import CrossRootArithmeticError
from benchmarks.held_one_root import MixedRootPointError
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from benchmarks.held_reference_geometry import project_boundary
from compas_cgal.straight_skeleton_2 import offset_polygon
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

F = Fraction
R = F(1)
# Chord deviation limit (mm) of the old domain's polygonal projection in these
# fixtures: the published cases use the same order (0.02-0.05 mm sagittae).
PROJECTION_LIMIT = F(1, 20)
# Boundary stroke width (mm) required by ReferenceBoundary; unused by the offset.
FIXTURE_STROKE_WIDTH = 0.1
# Ω* construction budget on Monstera (docs/held_exact_centre_domain.md, "Time
# budget"); every small pocket must build well within it.
CENTRE_DOMAIN_BUILD_BUDGET_SECONDS = 2.0
# Neck widths: exactly 2r, and 2r plus one binary64-exact 2^-20 mm.
NECK_EXACT_TWO_R = F(2)
NECK_JUST_WIDER = F(2) + F(1, 2**20)

_Vector = tuple[Fraction, Fraction]
_Polygon = tuple[_Vector, ...]


# ---------------------------------------------------------------------------
# Fixture construction (declared boundaries, exact by construction)
# ---------------------------------------------------------------------------


def _point(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(F(x), F(y))


class _Pocket:
    """Sequential builder of a declared CCW cycle starting at ``start``."""

    def __init__(self, start: ExactPoint2[WorldXY]) -> None:
        self.cursor = start
        self.primitives: list[DeclaredLine | DeclaredArc] = []

    def line(self, quarter: int, length: Fraction | int | float) -> _Pocket:
        line = DeclaredLine.build(self.cursor, PythagoreanDirection.build(F(0), quarter), F(length), len(self.primitives))
        self.primitives.append(line)
        self.cursor = line.end
        return self

    def arc(self, end: ExactPoint2[WorldXY], t: Fraction, quarter: int, sense: int) -> _Pocket:
        assert sense in (1, -1)
        arc = DeclaredArc.build(self.cursor, end, PythagoreanDirection.build(t, quarter), 1 if sense == 1 else -1, len(self.primitives))
        self.primitives.append(arc)
        self.cursor = arc.end
        return self

    def close(self, name: str) -> DeclaredBoundary:
        return DeclaredBoundary.build(name, tuple(self.primitives), hashlib.sha256(name.encode()).hexdigest())


def _square_with_bottom_arc(name: str, t: Fraction, sense: int) -> DeclaredBoundary:
    return _Pocket(_point(0, 0)).arc(_point(20, 0), t, 0, sense).line(1, 20).line(2, 20).line(3, 20).close(name)


def concave_arc_pocket() -> DeclaredBoundary:
    """Arc about (10, -7.5), ρ = 12.5, clockwise from (0,0) to (20,0): bulges into P to (10, 5)."""
    return _square_with_bottom_arc("concave_arc", F(1, 2), -1)


def convex_arc_pocket() -> DeclaredBoundary:
    """Arc about (10, 7.5), ρ = 12.5, counter-clockwise from (0,0) to (20,0): bulges out to (10, -5)."""
    return _square_with_bottom_arc("convex_arc", F(-1, 2), 1)


def l_shape_pocket() -> DeclaredBoundary:
    """(0,0), (20,0), (20,10), (10,10), (10,20), (0,20): reflex corner at junction 3 = (10,10)."""
    return _Pocket(_point(0, 0)).line(0, 20).line(1, 10).line(2, 10).line(1, 10).line(2, 10).line(3, 20).close("l_shape")


def neck_pocket(width: Fraction) -> DeclaredBoundary:
    """Two 10 x 10 squares joined by a 4 mm channel of ``width`` centred on y = 5."""
    low, high = 5 - width / 2, 5 + width / 2
    pocket = _Pocket(_point(0, 0)).line(0, 10).line(1, low).line(0, 4).line(3, low).line(0, 10).line(1, 10)
    return pocket.line(2, 10).line(3, 10 - high).line(2, 4).line(1, 10 - high).line(2, 10).line(3, 10).close(f"neck_{width}")


def tight_convex_arc_pocket() -> DeclaredBoundary:
    """20 x 20 square whose (20,20) corner is a convex arc of ρ = 0.5 < r about (19.5, 19.5) (primitive 2)."""
    return _Pocket(_point(0, 0)).line(0, 20).line(1, F(39, 2)).arc(_point(F(39, 2), 20), F(0), 1, 1).line(2, F(39, 2)).line(3, 20).close("tight_convex_arc")


# ---------------------------------------------------------------------------
# The old domain: float straight-skeleton offset of the polygonal projection
# ---------------------------------------------------------------------------


def _reference_boundary(boundary: DeclaredBoundary) -> ReferenceBoundary:
    primitives: list[ReferencePrimitive] = []
    for primitive in boundary.primitives:
        start = Point2[WorldXY].build(*primitive.start.as_floats())
        end = Point2[WorldXY].build(*primitive.end.as_floats())
        if isinstance(primitive, DeclaredLine):
            primitives.append(ReferenceLine.build(start, end))
            continue
        centre = primitive.centre.as_floats()
        start_angle = math.atan2(float(primitive.start.y) - centre[1], float(primitive.start.x) - centre[0])
        end_angle = math.atan2(float(primitive.end.y) - centre[1], float(primitive.end.x) - centre[0])
        sweep = primitive.sense * ((primitive.sense * (end_angle - start_angle)) % math.tau)
        primitives.append(ReferenceArc.build(start, end, Point2[WorldXY].build(*centre), Radian(sweep)))
    return ReferenceBoundary.build(primitives, ToolRadius.build(float(R)), Millimetre(FIXTURE_STROKE_WIDTH))


def old_domain(boundary: DeclaredBoundary) -> tuple[_Polygon, ...]:
    """The old Ω exactly as ``figure7_inward_offset`` builds it, with exact vertex coordinates."""
    projection = project_boundary(_reference_boundary(boundary), Millimetre(float(PROJECTION_LIMIT)))
    polygons = offset_polygon([(float(point.x), float(point.y), 0.0) for point in projection.points], float(R))
    return tuple(tuple((F(point[0]), F(point[1])) for point in polygon.points) for polygon in polygons)


def _edges(polygons: tuple[_Polygon, ...]) -> list[tuple[_Vector, _Vector]]:
    return [(polygon[index - 1], polygon[index]) for polygon in polygons for index in range(len(polygon))]


def _nearest_on_segment(point: _Vector, start: _Vector, end: _Vector) -> tuple[Fraction, _Vector]:
    segment = (end[0] - start[0], end[1] - start[1])
    relative = (point[0] - start[0], point[1] - start[1])
    parameter = min(F(1), max(F(0), (relative[0] * segment[0] + relative[1] * segment[1]) / (segment[0] ** 2 + segment[1] ** 2)))
    foot = (start[0] + parameter * segment[0], start[1] + parameter * segment[1])
    return (point[0] - foot[0]) ** 2 + (point[1] - foot[1]) ** 2, foot


def _strictly_inside(polygons: tuple[_Polygon, ...], point: _Vector) -> bool:
    """Exact even-odd crossing test (the point is off every edge in these fixtures)."""
    inside = False
    for start, end in _edges(polygons):
        if (start[1] > point[1]) != (end[1] > point[1]):
            crossing = start[0] + (point[1] - start[1]) * (end[0] - start[0]) / (end[1] - start[1])
            inside ^= point[0] < crossing
    return inside


def _cross(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def _owner(kind: BoundaryFeatureKind, index: int) -> BoundaryFeature:
    return BoundaryFeature(kind, index)


def _supports_of(domain: CentreDomain, piece: CentreDomainPiece) -> None:
    """Every piece lies on the exact r-offset of each of its owners (clearance r by construction)."""
    for owner in piece.owners:
        assert domain.supports.support_of(owner) == piece.support
    assert piece.support.contains(piece.source) and piece.support.contains(piece.target)


# ---------------------------------------------------------------------------
# Design tests 1-5: the fixtures
# ---------------------------------------------------------------------------


def test_concave_arc_pocket_old_domain_gouges_new_does_not() -> None:
    boundary = concave_arc_pocket()
    centre, rho = (F(10), F(-15, 2)), F(25, 2)

    # Old: some old-Ω centre is strictly closer than ρ + r to the arc centre, in
    # the arc's wedge, so its tool disk crosses the arc (a gouge); the overshoot
    # is bounded by the projection's chord deviation (the sagitta).
    nearest_squared, foot = min(_nearest_on_segment(centre, start, end) for start, end in _edges(old_domain(boundary)))
    radial = (foot[0] - centre[0], foot[1] - centre[1])
    assert _cross((F(-10), F(15, 2)), radial) <= 0 <= _cross((F(10), F(15, 2)), radial)
    assert nearest_squared < (rho + R) ** 2
    assert nearest_squared >= (rho + R - PROJECTION_LIMIT) ** 2

    # New, pure-rational part: the offset arc and the one-root corner vertex.
    supports = CentreDomainSupports.build(boundary, R)
    arc_support = supports.support_of(_owner(BoundaryFeatureKind.ARC, 0))
    assert isinstance(arc_support, OffsetCircleSupport)
    assert arc_support.kind is PieceKind.OFFSET_ARC
    assert (arc_support.centre.x, arc_support.centre.y, arc_support.radius) == (F(10), F(-15, 2), F(27, 2))
    corner = supports.junction_vertex(1)
    assert corner.kind is JunctionVertexKind.CONVEX
    expected = OneRootPoint2.build(OneRootNumber.of(19), OneRootNumber.build(F(-15, 2), 1, F(405, 4)))
    assert corner.point == expected
    assert corner.point == OneRootPoint2.build(OneRootNumber.of(19), OneRootNumber.build(F(-15, 2), F(9, 2), 5))
    assert arc_support.contains(corner.point)
    assert supports.support_of(_owner(BoundaryFeatureKind.LINE, 1)).contains(corner.point)

    # New: Ω*.
    domain = CentreDomain.build(boundary, R)
    arc_pieces = [piece for piece in domain.pieces if piece.kind is PieceKind.OFFSET_ARC]
    assert len(arc_pieces) == 1
    assert arc_pieces[0].support == arc_support
    assert any(piece.source == corner.point or piece.target == corner.point for piece in domain.pieces)
    for piece in domain.pieces:
        _supports_of(domain, piece)
    apex = _point(10, 6)
    assert domain.locate(apex) is PointLocation.BOUNDARY
    assert (apex.y - centre[1]) - rho == R


def test_convex_arc_pocket_old_domain_leaves_a_sliver_new_reaches_it() -> None:
    boundary = convex_arc_pocket()
    midpoint = (F(10), F(-5))

    # Old: the arc midpoint is farther than r from the old Ω (unreachable sliver).
    nearest_squared = min(_nearest_on_segment(midpoint, start, end)[0] for start, end in _edges(old_domain(boundary)))
    assert nearest_squared > R**2

    supports = CentreDomainSupports.build(boundary, R)
    arc_support = supports.support_of(_owner(BoundaryFeatureKind.ARC, 0))
    assert isinstance(arc_support, OffsetCircleSupport)
    assert (arc_support.centre.x, arc_support.centre.y, arc_support.radius) == (F(10), F(15, 2), F(23, 2))
    assert arc_support.contains(OneRootPoint2.rational(10, -4))

    domain = CentreDomain.build(boundary, R)
    assert [piece.support for piece in domain.pieces if piece.kind is PieceKind.OFFSET_ARC] == [arc_support]
    assert domain.squared_distance_to(_point(*midpoint)) == R**2
    assert domain.locate(_point(10, -4)) is PointLocation.BOUNDARY


def test_reflex_corner_has_a_junction_arc() -> None:
    boundary = l_shape_pocket()
    witness = (F(47, 5), F(46, 5))

    # Old: the mitre vertex (9, 9) is √2 from the corner (lost reach), and the
    # witness on the true ∂Ω* is outside the old Ω.
    old = old_domain(boundary)
    assert (F(9), F(9)) in {vertex for polygon in old for vertex in polygon}
    assert (F(9) - 10) ** 2 + (F(9) - 10) ** 2 == 2
    assert not _strictly_inside(old, witness)

    supports = CentreDomainSupports.build(boundary, R)
    junction = supports.support_of(_owner(BoundaryFeatureKind.JUNCTION, 3))
    assert isinstance(junction, OffsetCircleSupport)
    assert junction.kind is PieceKind.JUNCTION_ARC
    assert (junction.centre, junction.radius, junction.sense) == (_point(10, 10), R, -1)
    assert (junction.raw_source, junction.raw_target) == (_point(10, 9), _point(9, 10))
    assert junction.contains(OneRootPoint2.rational(*witness))

    domain = CentreDomain.build(boundary, R)
    (piece,) = [piece for piece in domain.pieces if piece.kind is PieceKind.JUNCTION_ARC]
    assert piece.owners == frozenset({_owner(BoundaryFeatureKind.JUNCTION, 3)})
    assert (piece.source, piece.target) == (OneRootPoint2.rational(10, 9), OneRootPoint2.rational(9, 10))
    assert domain.locate(_point(*witness)) is PointLocation.BOUNDARY
    site = CentreDomainSite.build(piece, OneRootPoint2.rational(*witness))
    assert CentreDomainSite.build(piece, piece.source).compare(site) < 0 < CentreDomainSite.build(piece, piece.target).compare(site)


def test_neck_of_width_exactly_2r_is_a_hard_error() -> None:
    boundary = neck_pocket(NECK_EXACT_TWO_R)
    supports = CentreDomainSupports.build(boundary, R)
    channel_floor = supports.support_of(_owner(BoundaryFeatureKind.LINE, 2))
    channel_ceiling = supports.support_of(_owner(BoundaryFeatureKind.LINE, 8))
    assert isinstance(channel_floor, OffsetLineSupport) and isinstance(channel_ceiling, OffsetLineSupport)
    # The two channel offsets coincide on y = 5: Ω* there is a measure-zero segment.
    assert channel_floor.raw_source.y == channel_ceiling.raw_source.y == 5

    with pytest.raises(NeckSplitsCentreDomainError) as raised:
        CentreDomain.build(boundary, R)
    assert raised.value.component_count == 2


def test_neck_wider_than_2r_keeps_one_component() -> None:
    boundary = neck_pocket(NECK_JUST_WIDER)
    assert all(Fraction(float(primitive.length)) == primitive.length for primitive in boundary.primitives if isinstance(primitive, DeclaredLine))
    domain = CentreDomain.build(boundary, R)
    assert len(domain.components) == 1
    assert domain.locate(_point(12, 5)) is PointLocation.INTERIOR


def test_convex_arc_tighter_than_tool_has_no_offset_piece() -> None:
    boundary = tight_convex_arc_pocket()
    arc = boundary.primitives[2]
    assert isinstance(arc, DeclaredArc) and (arc.centre, arc.radius) == (_point(F(39, 2), F(39, 2)), F(1, 2))

    supports = CentreDomainSupports.build(boundary, R)
    assert supports.tight_convex_arcs == (2,)
    assert supports.offset_arcs == ()
    with pytest.raises(NoSupportError):
        supports.support_of(_owner(BoundaryFeatureKind.ARC, 2))
    assert supports.deferred_junctions == (2, 3)

    domain = CentreDomain.build(boundary, R)
    assert not [piece for piece in domain.pieces if piece.kind is PieceKind.OFFSET_ARC]
    assert any(vertex.centre == OneRootPoint2.rational(19, 19) for vertex in domain.vertex_disks)
    tight = [part for part in domain.unreachable if part.cause is UnreachableCause.TIGHT_CONVEX_ARC]
    assert [part.owners for part in tight] == [frozenset({_owner(BoundaryFeatureKind.ARC, 2)})]


# ---------------------------------------------------------------------------
# Design tests 6-7: the four reference cases
# ---------------------------------------------------------------------------


@cache
def _case(name: str) -> HeldReferenceCase:
    return load_held_reference_case(name)


@cache
def _declared(name: str) -> DeclaredBoundary:
    return declared_boundary_from_reference(name, _case(name).boundary)


def _tool_radius(name: str) -> Fraction:
    return F(float(_case(name).tool_radius.value))


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_centre_domain_boundary_closes_exactly(name: str) -> None:
    supports = CentreDomainSupports.build(_declared(name), _tool_radius(name))
    for line in supports.offset_lines:
        assert line.normal[0] ** 2 + line.normal[1] ** 2 == 1
        assert line.contains(OneRootPoint2.from_exact(line.raw_source)) and line.contains(OneRootPoint2.from_exact(line.raw_target))
    for circle in (*supports.offset_arcs, *supports.junction_arcs):
        assert circle.radius > 0
        assert circle.contains(OneRootPoint2.from_exact(circle.raw_source)) and circle.contains(OneRootPoint2.from_exact(circle.raw_target))
    for vertex in supports.junction_vertices:
        assert supports.support_of(vertex.incoming).contains(vertex.point)
        assert supports.support_of(vertex.outgoing).contains(vertex.point)

    domain = CentreDomain.build(_declared(name), _tool_radius(name))
    assert len(domain.components) == 1
    for component in domain.components:
        for ccb in (component.outer, *component.holes):
            for current, following in zip(ccb, (*ccb[1:], ccb[0])):
                assert current.target == following.source
                _supports_of(domain, current)


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_centre_domain_is_built_on_the_recertified_boundary(name: str) -> None:
    declared = _declared(name)
    certificate = DeclaredBoundaryCertificate.build(_case(name), declared)
    assert certificate.declared_bound <= certificate.reconstruction_limit
    assert CentreDomainSupports.build(declared, _tool_radius(name)).boundary == certificate.boundary

    domain = CentreDomain.build(declared, _tool_radius(name))
    assert domain.identity.boundary == certificate.boundary
    assert domain.identity.tool_radius == _tool_radius(name)
    assert domain.identity.sha256 == CentreDomain.build(declared, _tool_radius(name)).identity.sha256


# ---------------------------------------------------------------------------
# Design test 8: generic-double witness
# ---------------------------------------------------------------------------


@st.composite
def generic_concave_pockets(draw: st.DrawFn) -> DeclaredBoundary:
    """Rectangle with a concave bottom arc, every coordinate a generic binary64.

    ``x0, y0 ∈ [100, 200]`` and ``w, h <= 50`` make ``x0 + w - x0`` exact
    (Sterbenz), so every line length is a binary64.
    """
    x0 = draw(st.floats(100.0, 200.0))
    y0 = draw(st.floats(100.0, 200.0))
    x1 = x0 + draw(st.floats(8.0, 40.0))
    y1 = y0 + draw(st.floats(48.0, 50.0))
    width, height = F(x1) - F(x0), F(y1) - F(y0)
    t = F(draw(st.integers(2**40, 2**54)), 2**56)
    start, end = ExactPoint2[WorldXY].from_binary64(x0, y0), ExactPoint2[WorldXY].from_binary64(x1, y0)
    return _Pocket(start).arc(end, t, 0, -1).line(1, height).line(2, width).line(3, height).close(f"generic_{x0!r}_{y0!r}_{x1!r}_{y1!r}_{t}")


@settings(max_examples=20, deadline=None)
@given(generic_concave_pockets())
def test_generic_double_witness(boundary: DeclaredBoundary) -> None:
    supports = CentreDomainSupports.build(boundary, R)
    for vertex in supports.junction_vertices:
        assert vertex.kind is JunctionVertexKind.CONVEX
        assert supports.support_of(vertex.incoming).contains(vertex.point)
        assert supports.support_of(vertex.outgoing).contains(vertex.point)
    assert any(not vertex.point.x.is_rational or not vertex.point.y.is_rational for vertex in supports.junction_vertices)

    started = time.perf_counter()
    domain = CentreDomain.build(boundary, R)
    assert time.perf_counter() - started <= CENTRE_DOMAIN_BUILD_BUDGET_SECONDS
    arc = boundary.primitives[0]
    assert isinstance(arc, DeclaredArc)
    apex = arc.centre.offset((F(0), F(1)), arc.radius + R)
    assert domain.locate(apex) is PointLocation.BOUNDARY
    for piece in domain.pieces:
        _supports_of(domain, piece)


# ---------------------------------------------------------------------------
# Pure-rational parts (green): one-root numbers, order, contact map
# ---------------------------------------------------------------------------


def test_one_root_numbers_fold_and_compare_exactly() -> None:
    assert OneRootNumber.build(1, 2, F(9, 4)).is_rational and OneRootNumber.build(1, 2, F(9, 4)) == 4
    assert OneRootNumber.build(0, 1, 2) > F(14142, 10000) and OneRootNumber.build(0, 1, 2) < F(14143, 10000)
    assert OneRootNumber.build(3, 1, 2) > OneRootNumber.build(2, 1, 5)  # 4.414 > 4.236, across roots
    assert OneRootNumber.build(F(-15, 2), 1, F(405, 4)) == OneRootNumber.build(F(-15, 2), F(9, 2), 5)
    assert OneRootNumber.build(1, -1, 2).sign() == -1 and OneRootNumber.build(-1, 1, 2).sign() == 1
    product = OneRootNumber.build(1, 1, 2) * OneRootNumber.build(1, -1, 2)
    assert product.is_rational and product == -1


def test_cross_root_arithmetic_and_mixed_points_are_refused() -> None:
    with pytest.raises(CrossRootArithmeticError):
        OneRootNumber.build(0, 1, 2) + OneRootNumber.build(0, 1, 3)
    with pytest.raises(MixedRootPointError):
        OneRootPoint2.build(OneRootNumber.build(0, 1, 2), OneRootNumber.build(0, 1, 3))


def test_tool_radius_must_be_a_positive_binary64() -> None:
    with pytest.raises(InvalidCentreDomainToolRadiusError):
        CentreDomainSupports.build(l_shape_pocket(), F(1, 3))
    with pytest.raises(InvalidCentreDomainToolRadiusError):
        CentreDomainSupports.build(l_shape_pocket(), F(0))


def test_junction_arc_sites_are_ordered_exactly_across_roots() -> None:
    junction = CentreDomainSupports.build(l_shape_pocket(), R).support_of(_owner(BoundaryFeatureKind.JUNCTION, 3))
    assert isinstance(junction, OffsetCircleSupport)
    source, target = OneRootPoint2.from_exact(junction.raw_source), OneRootPoint2.from_exact(junction.raw_target)
    piece = CentreDomainPiece(0, 0, 0, frozenset({junction.owner}), junction, source, target)
    root_half = OneRootNumber.build(10, -1, F(1, 2))
    ordered = [
        source,
        OneRootPoint2.rational(F(47, 5), F(46, 5)),
        OneRootPoint2.build(root_half, root_half),
        OneRootPoint2.rational(F(46, 5), F(47, 5)),
        target,
    ]
    sites = [CentreDomainSite.build(piece, point) for point in ordered]
    for index, site in enumerate(sites):
        for other_index, other in enumerate(sites):
            assert site.compare(other) == (index > other_index) - (index < other_index)
    assert compare_on_arc(junction.centre, 1, source, ordered[1], ordered[3]) == 1


def test_contact_map_adds_no_root() -> None:
    supports = CentreDomainSupports.build(concave_arc_pocket(), R)
    arc_support = supports.support_of(_owner(BoundaryFeatureKind.ARC, 0))
    assert isinstance(arc_support, OffsetCircleSupport)
    corner = supports.junction_vertex(1).point
    contact = arc_support.contact(corner)
    assert contact.radicand == corner.radicand
    dx, dy = contact.minus(arc_support.centre)
    assert dx * dx + dy * dy == F(25, 2) ** 2
    line_support = supports.support_of(_owner(BoundaryFeatureKind.LINE, 1))
    assert line_support.contact(corner).x == 20
