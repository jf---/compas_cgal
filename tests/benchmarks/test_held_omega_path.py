"""Held paths on Ω*: contacts on ∂Ω*, exact order, arc connectors, exact containment, zero gouges.

The concave-arc fixture (docs/held_exact_centre_domain.md, design test 1) is
where the polygon Ω gouges: its chords cut inside the offset arc. The same
binary64 contacts placed by the polygon-Ω path and by the Ω* path are audited
by one exact oracle against the declared boundary.
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_exact_centre_domain import AmbiguousContactProjectionError
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import CentreDomainSite
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_path import _PlacementLane
from benchmarks.held_figure5_engagement_refinement import circle_on_side
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_gouge_audit import AuditCircle
from benchmarks.held_gouge_audit import AuditConnector
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_audit import HeldGougeAudit
from benchmarks.held_offset_contact_projection import project_bounded_offset_contacts
from benchmarks.held_omega_circles import NoContainedCircleError
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_circles import contain
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_gouge_audit import audit_polygon_path
from benchmarks.held_omega_path import _split_at_gaps
from benchmarks.held_omega_refinement import omega_midpoint_site
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_omega_sites import single_loop_pieces
from benchmarks.held_omega_sites import sort_sites
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import _reference_boundary
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import convex_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_exact_centre_domain import old_domain

F = Fraction
R = F(1)
# Hypothesis contacts: points of the true offset arc (centre (10, -7.5), radius
# 13.5) at these x, rounded to binary64 as a guide would deliver them.
ARC_CONTACT_XS = (4.0, 7.0, 10.0, 13.0, 16.0)
# Guide radius (mm) of the fixture circles: smaller than every chord.
GUIDE_RADIUS = 0.5
# Admissible gap (mm) between a hypothesis contact and ∂Ω*: covers the
# fixtures' rounding and the polygon's sagitta (<= 1/20 mm).
GAP = ExactMillimetre(F(1, 4))


def _arc_contacts() -> list[tuple[float, float]]:
    return [(x, -7.5 + math.sqrt(13.5**2 - (x - 10.0) ** 2)) for x in ARC_CONTACT_XS]


def _provenance(index: int) -> tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]:
    return GuideRunId(0), GuideRunStationOrdinal(index), ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), F(0), 1)


def _omega_circles(domain: CentreDomain, clearance: DeclaredClearance, contacts: list[tuple[float, float]]) -> tuple[OmegaCircle, ...]:
    decisions = domain.decide_contacts(contacts, GAP).contacts
    loop = OmegaLoop.build(domain)
    sites = []
    for decision in decisions:
        assert decision is not None
        sites.append(OmegaSite.decided(loop, decision))
    order = sort_sites(sites)
    return tuple(contain(place_local(sites[index], GUIDE_RADIUS, _provenance(index)), clearance, domain) for index in order)


def _polygon_path(boundary: DeclaredBoundary, contacts: list[tuple[float, float]]) -> tuple[list[tuple[Point2[WorldXY], float]], list[tuple[Point2[WorldXY], ...]]]:
    """The polygon-Ω path's own placement: nearest polygon site, side-normal circle, polyline connectors."""
    (polygon,) = old_domain(boundary)
    component = tuple(Point2[WorldXY].build(float(x), float(y)) for x, y in polygon)
    points = tuple(Point2[WorldXY].build(x, y) for x, y in contacts)
    sites = project_bounded_offset_contacts(points, component, Millimetre(float(GAP)))
    circles = []
    for index, site in enumerate(sites):
        contact = _canonical_offset_point(component, site)
        source = Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(index), site, site, contact, contact, GuideRadius.build(GUIDE_RADIUS))
        circles.append(circle_on_side(source, site, contact, GuideRadius.build(GUIDE_RADIUS), component))
    ordered = tuple(sorted(circles, key=lambda circle: (int(circle.boundary_site.segment_id), circle.boundary_site.parameter)))
    connectors = Figure5BoundaryConnectors.build(component, ordered)
    return [(circle.center, float(circle.radius.value)) for circle in ordered], [transition.samples for transition in connectors.transitions]


def test_concave_arc_polygon_path_gouges_and_the_omega_path_does_not() -> None:
    boundary = concave_arc_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    contacts = _arc_contacts()

    old_circles, old_polylines = _polygon_path(boundary, contacts)
    old = audit_polygon_path(clearance, old_circles, old_polylines)
    assert old.gouging_circles and old.gouging_connectors
    assert all(item.primitive_kind == "concave-arc" for item in old.gouges)

    domain = CentreDomain.build(boundary, R)
    circles = _omega_circles(domain, clearance, contacts)
    connectors = OmegaConnectors.build(single_loop_pieces(domain), circles)
    new = audit_omega_path(clearance, circles, connectors.connectors)
    assert not new.gouges
    assert new.touching_connectors


def test_omega_connectors_follow_the_offset_arc_and_touch_it_exactly() -> None:
    boundary = concave_arc_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    circles = _omega_circles(domain, clearance, _arc_contacts())
    connectors = OmegaConnectors.build(single_loop_pieces(domain), circles)
    assert connectors.seam is Figure5SeamClosure.CONNECTOR
    assert len(connectors.connectors) == len(circles)
    for connector in connectors.connectors[:-1]:
        (piece,) = connector.pieces
        assert piece.piece.kind is PieceKind.OFFSET_ARC
        element = piece.element()
        assert isinstance(element, ArcElement)
        assert clearance.element(element).clearance is Clearance.TOUCH
    seam_kinds = {piece.piece.kind for piece in connectors.connectors[-1].pieces}
    assert seam_kinds == {PieceKind.OFFSET_ARC, PieceKind.OFFSET_LINE}


def test_lone_circle_seam_is_one_full_lap() -> None:
    boundary = l_shape_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    (circle,) = _omega_circles(domain, clearance, [(5.0, 1.0)])
    pieces = single_loop_pieces(domain)
    connectors = OmegaConnectors.build(pieces, (circle,))
    (lap,) = connectors.connectors
    assert connectors.seam is Figure5SeamClosure.CONNECTOR
    assert len(lap.pieces) == len(pieces) + 1
    assert lap.start == lap.end == OneRootPoint2.rational(5, 1)


def test_circles_on_omega_lines_touch_exactly_and_are_never_shrunk_there() -> None:
    boundary = l_shape_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    circles = _omega_circles(domain, clearance, [(5.0, 0.75), (5.0, 1.25)])
    for circle in circles:
        assert (float(circle.center.x), float(circle.center.y), float(circle.radius.value)) == (5.0, 1.5, 0.5)
        assert clearance.circle(circle.element()).clearance is Clearance.TOUCH
        assert not circle.limited


def test_circle_at_a_convex_corner_is_shrunk_until_contained() -> None:
    boundary = l_shape_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    # Site on the bottom offset line one guide radius short of the convex
    # corner (19, 1): a radius-2 circle there crosses the right offset line.
    decision = domain.decide_contacts([(18.0, 1.0)], GAP).contacts[0]
    assert decision is not None
    site = OmegaSite.decided(OmegaLoop.build(domain), decision)
    local = place_local(site, 2.0, _provenance(0))
    assert clearance.circle(local.element()).clearance is Clearance.GOUGE
    contained = contain(local, clearance, domain)
    assert contained.limited and float(contained.radius.value) <= 1.0
    assert clearance.circle(contained.element()).clearance is not Clearance.GOUGE


def test_site_at_a_convex_vertex_admits_no_circle() -> None:
    boundary = l_shape_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    decision = domain.decide_contacts([(19.125, 0.875)], GAP).contacts[0]
    assert decision is not None and decision.is_vertex
    site = OmegaSite.decided(OmegaLoop.build(domain), decision)
    with pytest.raises(NoContainedCircleError):
        contain(place_local(site, 1.0, _provenance(0)), clearance, domain)


def test_clearance_decides_touch_and_gouge_exactly() -> None:
    boundary = convex_arc_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    # Convex arc about (10, 7.5), ρ = 12.5, top side y = 20: the circle about
    # (10, 15) of radius 4 reaches y = 19 and 11.5 from the arc centre, so it
    # touches both at distance r = 1; 2^-40 mm larger gouges both.
    centre = ExactPoint2[WorldXY].build(F(10), F(15))
    assert clearance.circle(CircleElement(centre, F(4))).clearance is Clearance.TOUCH
    assert clearance.circle(CircleElement(centre, F(4) + F(1, 2**40))).clearance is Clearance.GOUGE
    left, right = ExactPoint2[WorldXY].build(F(5), F(19)), ExactPoint2[WorldXY].build(F(15), F(19))
    assert clearance.element(SegmentElement.between(left, right)).clearance is Clearance.TOUCH
    lower = ExactPoint2[WorldXY].build(F(15), F(19) + F(1, 2**40))
    assert clearance.element(SegmentElement.between(left, lower)).clearance is Clearance.GOUGE


@settings(max_examples=60, deadline=None)
@given(st.lists(st.tuples(st.floats(min_value=-0.5, max_value=20.5), st.floats(min_value=-0.5, max_value=20.5)), min_size=2, max_size=8))
def test_site_keys_order_as_the_exact_site_comparison(points: list[tuple[float, float]]) -> None:
    domain = _domain_cache("concave")
    try:
        decisions = domain.decide_contacts(points, ExactMillimetre(F(40))).contacts
    except AmbiguousContactProjectionError:
        assume(False)
        return
    loop = OmegaLoop.build(domain)
    sites = [OmegaSite.decided(loop, decision) for decision in decisions if decision is not None]
    exact = [CentreDomainSite.build(site.piece, site.point()) for site in sites]
    for first, first_exact in zip(sites, exact):
        for second, second_exact in zip(sites, exact):
            expected = first_exact.compare(second_exact)
            assert (first.key < second.key) == (expected < 0)
            assert first.same_position(second) == (expected == 0)


_DOMAINS: dict[str, CentreDomain] = {}


def _domain_cache(name: str) -> CentreDomain:
    if name not in _DOMAINS:
        _DOMAINS[name] = CentreDomain.build(concave_arc_pocket(), R)
    return _DOMAINS[name]


# ---------------------------------------------------------------------------
# Contract: the declared-boundary oracle against the stored-boundary audit
# ---------------------------------------------------------------------------

_FIXTURES = {"concave": concave_arc_pocket, "convex": convex_arc_pocket, "l_shape": l_shape_pocket}


@settings(max_examples=150, deadline=None)
@given(
    st.sampled_from(tuple(_FIXTURES)),
    st.floats(min_value=-2, max_value=22),
    st.floats(min_value=-2, max_value=22),
    st.floats(min_value=2**-10, max_value=12),
    st.floats(min_value=-2, max_value=22),
    st.floats(min_value=-2, max_value=22),
)
def test_declared_oracle_agrees_with_the_stored_dyadic_audit(fixture: str, x: float, y: float, radius: float, u: float, v: float) -> None:
    """Fixtures whose stored floats are exact: two independent exact audits must decide alike (filters included)."""
    boundary = _FIXTURES[fixture]()
    clearance = DeclaredClearance.build(boundary, R)
    stored = HeldGougeAudit.build(_reference_boundary(boundary).primitives, Millimetre(float(R)))
    centre, other = Point2[WorldXY].build(x, y), Point2[WorldXY].build(u, v)
    report = stored.audit([AuditCircle(centre, Millimetre(radius))], [AuditConnector((centre, other))])
    expected_circle = next((item.clearance for item in report.elements if item.element == "circle"), Clearance.CLEAR)
    expected_segment = next((item.clearance for item in report.elements if item.element == "connector"), Clearance.CLEAR)
    exact_centre = ExactPoint2[WorldXY].from_binary64(x, y)
    assert clearance.circle(CircleElement(exact_centre, F(radius))).clearance is expected_circle
    if (x, y) != (u, v):
        assert clearance.element(SegmentElement.between(exact_centre, ExactPoint2[WorldXY].from_binary64(u, v))).clearance is expected_segment


@pytest.mark.parametrize("fixture", sorted(_FIXTURES))
def test_declared_oracle_decides_exact_touches_like_the_stored_audit(fixture: str) -> None:
    """Every ∂Ω* piece of the fixture touches the boundary exactly; full-circle pieces and lines included."""
    boundary = _FIXTURES[fixture]()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    for piece in single_loop_pieces(domain):
        support = piece.support
        element: ArcElement | SegmentElement
        if isinstance(support, OffsetLineSupport):
            element = SegmentElement(piece.source, piece.target, support.direction, support.normal, support.level)
        else:
            element = ArcElement(support.centre, support.radius, support.sense, piece.source, piece.target)
        assert clearance.element(element).clearance is Clearance.TOUCH


def test_guide_radius_above_a_convex_offset_arc_is_capped_at_its_radius() -> None:
    """Convex arc about (10, 7.5), ρ = 12.5: its offset arc has ρ' = 11.5; a 20 mm guide circle is capped, not refused."""
    boundary = convex_arc_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    decision = domain.decide_contacts([(10.0, -4.0)], GAP).contacts[0]
    assert decision is not None and decision.piece.kind is PieceKind.OFFSET_ARC
    circle = place_local(OmegaSite.decided(OmegaLoop.build(domain), decision), 20.0, _provenance(0))
    assert circle.limited and float(circle.radius.value) <= 11.5
    assert clearance.circle(contain(circle, clearance, domain).element()).clearance is not Clearance.GOUGE


# ---------------------------------------------------------------------------
# Midpoints along ∂Ω*
# ---------------------------------------------------------------------------


def _sites(boundary: DeclaredBoundary, contacts: list[tuple[float, float]]) -> tuple[OmegaLoop, list[OmegaSite]]:
    domain = CentreDomain.build(boundary, R)
    loop = OmegaLoop.build(domain)
    decisions = domain.decide_contacts(contacts, GAP).contacts
    return loop, [OmegaSite.decided(loop, decision) for decision in decisions if decision is not None]


def test_midpoint_on_one_arc_piece_lies_strictly_between() -> None:
    loop, (first, second) = _sites(concave_arc_pocket(), _arc_contacts()[:2])
    middle = omega_midpoint_site(loop, first, second)
    assert middle.index == first.index == second.index
    assert first.key < middle.key < second.key


def test_midpoint_across_a_convex_corner_lies_on_the_chain() -> None:
    loop, (first, second) = _sites(l_shape_pocket(), [(15.0, 1.0), (19.0, 5.0)])
    middle = omega_midpoint_site(loop, first, second)
    assert first.key < middle.key < second.key
    assert not middle.at_vertex


def test_midpoint_of_a_fan_is_its_shared_site() -> None:
    loop, (first, second) = _sites(l_shape_pocket(), [(5.0, 0.75), (5.0, 1.25)])
    assert first.same_position(second)
    assert omega_midpoint_site(loop, first, second) is first


def test_lanes_split_where_consecutive_circles_leave_the_overlap_domain() -> None:
    """A lane that lost its middle circles (refused contacts) splits instead of failing placement."""
    boundary = l_shape_pocket()
    clearance = DeclaredClearance.build(boundary, R)
    domain = CentreDomain.build(boundary, R)
    circles = _omega_circles(domain, clearance, [(3.0, 1.0), (4.0, 1.0), (15.0, 1.0), (16.0, 1.0)])
    lanes = _split_at_gaps((_PlacementLane(0, 1, circles),), ToolRadius.build(1.0))
    assert [len(lane.circles) for lane in lanes] == [2, 2]
