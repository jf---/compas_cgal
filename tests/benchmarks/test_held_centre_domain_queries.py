"""Queries on the exact Held centre domain Ω*: locate, squared distance, contact projection.

Fixtures are Pythagorean by construction (integers, r = 1), reused from the
contract tests, so every expected value is an exact rational or one-root
number. Design: docs/held_exact_centre_domain.md, "Projection of a hypothesis
contact".
"""

from __future__ import annotations

import functools
import math
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_centre_domain_distance import NearestOnPiece
from benchmarks.held_centre_domain_distance import nearest_on_piece
from benchmarks.held_centre_domain_native import encode_boundary
from benchmarks.held_centre_domain_native import native_nearest
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_exact_centre_domain import AmbiguousContactProjectionError
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import ContactBeyondGapError
from benchmarks.held_exact_centre_domain import DecidedContact
from benchmarks.held_exact_centre_domain import PointLocation
from benchmarks.held_exact_centre_domain import UnreachableCause
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from compas_cgal import _held_centre_domain_2
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import _declared
from tests.benchmarks.test_held_exact_centre_domain import _Pocket
from tests.benchmarks.test_held_exact_centre_domain import _tool_radius
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import convex_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket

F = Fraction
R = F(1)
# Admissible projection gaps (mm) for the fixtures: one generous, one below the
# exact distance 1 - √½ ≈ 0.293 of the reflex-corner query.
WIDE_GAP = ExactMillimetre(F(10))
NARROW_GAP = ExactMillimetre(F(1, 4))
# Rational brackets of irrational vertex coordinates, 2^-70 mm wide: far below
# a binary64 enclosure of the same coordinate (about 2^-48 mm at 20 mm), so a
# query built from them sits off a vertex direction by less than any
# floating-point margin.
SUB_ULP_BRACKET_BITS = 70
# Radial offsets (mm) of the arc-extreme and vertex-direction queries.
RADIAL_OFFSETS = (F(-1, 8), F(1, 8))
# Reference cases the per-piece oracle sweeps run on: the pure-Fraction oracle
# costs O(pieces) per query on ~1000-bit numbers, so Monstera (487 pieces) and
# upper are left to the Hypothesis and contract tests.
ORACLE_REFERENCE_CASES = ("figure5", "figure8_crossed_skis")


def _point(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(F(x), F(y))


def _equal_radius_corner_pocket() -> DeclaredBoundary:
    """20 x 20 square whose (20,20) corner is a convex arc of ρ = r = 1 about (19, 19)."""
    return _Pocket(_point(0, 0)).line(0, 20).line(1, 19).arc(_point(19, 20), F(0), 1, 1).line(2, 19).line(3, 20).close("equal_radius_corner")


def test_reflex_corner_projects_onto_the_junction_arc_at_a_one_root_point() -> None:
    domain = CentreDomain.build(l_shape_pocket(), R)
    query = _point(F(19, 2), F(19, 2))
    assert domain.locate(query) is PointLocation.EXTERIOR
    root_half = OneRootNumber.build(10, -1, F(1, 2))
    assert domain.squared_distance_to(query) == OneRootNumber.build(F(3, 2), -2, F(1, 2))  # (1 - √½)²
    site = domain.project_contact(query, WIDE_GAP)
    assert site.piece.kind is PieceKind.JUNCTION_ARC
    assert site.point == OneRootPoint2.build(root_half, root_half)
    assert site.contact() == OneRootPoint2.rational(10, 10)


def test_projection_beyond_the_gap_raises() -> None:
    domain = CentreDomain.build(l_shape_pocket(), R)
    with pytest.raises(ContactBeyondGapError):
        domain.project_contact(_point(F(19, 2), F(19, 2)), NARROW_GAP)


def test_a_vertex_site_is_reported_on_the_following_piece() -> None:
    domain = CentreDomain.build(l_shape_pocket(), R)
    site = domain.project_contact(_point(F(1, 2), F(1, 2)), WIDE_GAP)
    assert site.point == OneRootPoint2.rational(1, 1)
    assert site.piece.source == site.point


def test_equidistant_distinct_contacts_are_ambiguous() -> None:
    domain = CentreDomain.build(l_shape_pocket(), R)
    assert domain.locate(_point(5, 5)) is PointLocation.INTERIOR
    assert domain.squared_distance_to(_point(5, 5)) == 0
    with pytest.raises(AmbiguousContactProjectionError):
        domain.project_contact(_point(5, 5), WIDE_GAP)


def test_concave_offset_arc_projection_lands_on_the_offset_circle() -> None:
    domain = CentreDomain.build(concave_arc_pocket(), R)
    site = domain.project_contact(_point(10, 7), WIDE_GAP)
    assert site.piece.kind is PieceKind.OFFSET_ARC
    assert site.point == OneRootPoint2.rational(10, 6)
    assert site.contact() == OneRootPoint2.rational(10, 5)


def test_arc_with_radius_equal_to_the_tool_is_reached_not_unreachable() -> None:
    domain = CentreDomain.build(_equal_radius_corner_pocket(), R)
    assert domain.supports.tight_convex_arcs == (2,)
    assert not [part for part in domain.unreachable if part.cause is UnreachableCause.TIGHT_CONVEX_ARC]
    assert any(vertex.centre == OneRootPoint2.rational(19, 19) for vertex in domain.vertex_disks)
    assert domain.locate(_point(19, 19)) is PointLocation.BOUNDARY
    corners = {part.owners for part in domain.unreachable if part.cause is UnreachableCause.CONVEX_KINK_CUSP}
    assert corners == {frozenset({BoundaryFeature(BoundaryFeatureKind.JUNCTION, junction)}) for junction in (0, 1, 4)}


def test_native_construction_rejects_a_non_unit_line_normal() -> None:
    wire = encode_boundary(l_shape_pocket())
    kind, sense, values = wire[0]
    wire[0] = (kind, sense, [*values[:4], ("0", "1"), ("2", "1")])
    with pytest.raises(_held_centre_domain_2.InvalidCentreDomainInputError):
        _held_centre_domain_2.CentreDomain2.build(wire, ("1", "1"))


def _circle_pocket() -> DeclaredBoundary:
    """Two semicircles of ρ = 5 about the origin: Ω* is the full circle of radius 4 (a ``piece.full`` piece)."""
    return _Pocket(_point(5, 0)).arc(_point(-5, 0), F(0), 1, 1).arc(_point(5, 0), F(0), 3, 1).close("circle")


def _stadium_pocket() -> DeclaredBoundary:
    """A 20 x 10 stadium, caps of ρ = 5 about (20, 5) and (0, 5): each cap's offset spans exactly π."""
    return _Pocket(_point(0, 0)).line(0, 20).arc(_point(20, 10), F(0), 0, 1).line(2, 20).arc(_point(0, 0), F(0), 2, 1).close("stadium")


def _rich_pocket() -> DeclaredBoundary:
    """30 x 20 box: a concave semicircular bump (ρ = 5 about (10, 0)) on the floor, a convex semicircular cap (ρ = 10) on the right."""
    pocket = _Pocket(_point(0, 0)).line(0, 5).arc(_point(15, 0), F(0), 1, -1).line(0, 15).arc(_point(30, 20), F(0), 0, 1)
    return pocket.line(2, 30).line(3, 20).close("rich")


FIXTURES = {
    "l_shape": l_shape_pocket,
    "concave": concave_arc_pocket,
    "convex": convex_arc_pocket,
    "circle": _circle_pocket,
    "stadium": _stadium_pocket,
    "rich": _rich_pocket,
}


@functools.cache
def _domain(fixture: str) -> CentreDomain:
    if fixture in FIXTURES:
        return CentreDomain.build(FIXTURES[fixture](), R)
    return CentreDomain.build(_declared(fixture), _tool_radius(fixture))


def _by_distance(first: NearestOnPiece, second: NearestOnPiece) -> int:
    return first.squared_distance.compare(second.squared_distance)


def _assert_native_matches_oracle(domain: CentreDomain, query: ExactPoint2[WorldXY]) -> None:
    """The native nearest equals the pure-Fraction oracle in distance, point and ambiguity (both directions)."""
    candidates = [nearest_on_piece(piece.support, piece.source, piece.target, query) for piece in domain.pieces]
    oracle = min(candidates, key=functools.cmp_to_key(_by_distance))
    tied = [candidate for candidate in candidates if candidate.squared_distance == oracle.squared_distance]
    distinct: list[OneRootPoint2] = []
    for candidate in tied:
        if not any(candidate.point == seen for seen in distinct):
            distinct.append(candidate.point)
    oracle_ambiguous = len(distinct) > 1 or any(candidate.ambiguous for candidate in tied)
    native = native_nearest(domain.native, query)
    assert native.squared_distance == oracle.squared_distance
    assert native.ambiguous == oracle_ambiguous
    if not oracle_ambiguous:
        assert native.point == oracle.point


@settings(max_examples=150, deadline=None)
@given(
    st.sampled_from(tuple(FIXTURES)),
    st.fractions(min_value=-6, max_value=41, max_denominator=64),
    st.fractions(min_value=-6, max_value=23, max_denominator=64),
)
def test_native_nearest_matches_the_python_oracle(fixture: str, x: Fraction, y: Fraction) -> None:
    _assert_native_matches_oracle(_domain(fixture), _point(x, y))


def test_batch_projection_equals_single_projection() -> None:
    domain = _domain("concave")
    points = [(10.0, 7.0), (0.5, 0.5), (5.25, 3.0), (30.0, 30.0)]
    contacts = domain.project_contacts(points, ExactMillimetre(F(3)))
    assert contacts[3] is None
    for point, contact in zip(points[:3], contacts[:3]):
        single = domain.project_contact(ExactPoint2[WorldXY].from_binary64(*point), ExactMillimetre(F(3)))
        assert isinstance(contact, DecidedContact)
        site = contact.site()
        assert contact.piece == single.piece and site.piece == single.piece and site.point == single.point


def test_batch_decisions_return_their_exact_fallback_count() -> None:
    """A gap tie whose double enclosures overlap (0.3 is not dyadic) is decided by the exact stage, and counted."""
    domain = _domain("l_shape")
    tie_gap = ExactMillimetre(1 - F(0.3))
    tie = domain.decide_contacts([(0.3, 5.0)], tie_gap)
    assert tie.contacts[0] is not None
    assert tie.exact_fallbacks >= 1
    clear = domain.decide_contacts([(0.3, 5.0)], ExactMillimetre(F(10)))
    assert clear.exact_fallbacks < tie.exact_fallbacks


# ---------------------------------------------------------------------------
# M3: tests that kill the defects the review found surviving
# ---------------------------------------------------------------------------


def test_contact_exactly_at_the_gap_is_projected_not_beyond() -> None:
    domain = _domain("l_shape")
    gap = ExactMillimetre(F(1, 2))
    query = _point(F(1, 2), 5)
    assert domain.squared_distance_to(query) == F(1, 4)
    site = domain.project_contact(query, gap)
    assert site.point == OneRootPoint2.rational(1, 5)
    (decided,) = domain.decide_contacts([(0.5, 5.0)], gap).contacts
    assert decided is not None and decided.piece == site.piece
    (beyond,) = domain.decide_contacts([(0.5, 5.0)], ExactMillimetre(F(1, 2) - F(1, 2**52))).contacts
    assert beyond is None
    # The same tie where the double enclosures overlap, so the exact stage decides it.
    tie_gap = ExactMillimetre(1 - F(0.3))
    (inexact,) = domain.decide_contacts([(0.3, 5.0)], tie_gap).contacts
    assert inexact is not None and inexact.site().point == OneRootPoint2.rational(1, 5)
    (inexact_beyond,) = domain.decide_contacts([(0.3, 5.0)], ExactMillimetre(tie_gap - F(1, 2**60))).contacts
    assert inexact_beyond is None


def test_query_at_a_full_circle_centre_is_ambiguous() -> None:
    domain = _domain("circle")
    (piece,) = domain.pieces
    assert piece.source == piece.target
    nearest = native_nearest(domain.native, _point(0, 0))
    assert nearest.ambiguous and nearest.squared_distance == 16
    with pytest.raises(AmbiguousContactProjectionError):
        domain.project_contact(_point(0, 0), WIDE_GAP)


@pytest.mark.parametrize(
    ("query", "squared", "point"),
    [((F(9, 2), 0), F(1, 4), (4, 0)), ((-2, 0), F(4), (-4, 0)), ((3, 0), F(1), (4, 0)), ((0, F(-9, 2)), F(1, 4), (0, -4)), ((F(12, 5), F(16, 5)), F(0), (F(12, 5), F(16, 5)))],
)
def test_full_circle_piece_projects_radially(query: tuple[Fraction | int, Fraction | int], squared: Fraction, point: tuple[Fraction | int, Fraction | int]) -> None:
    domain = _domain("circle")
    nearest = native_nearest(domain.native, _point(*query))
    assert not nearest.ambiguous
    assert nearest.squared_distance == squared
    assert nearest.point == OneRootPoint2.rational(*point)


def test_half_turn_piece_decides_both_sides_of_its_diameter() -> None:
    domain = _domain("stadium")
    caps = [piece for piece in domain.pieces if piece.kind is PieceKind.OFFSET_ARC]
    assert len(caps) == 2
    outside = native_nearest(domain.native, _point(F(49, 2), 5))
    assert not outside.ambiguous and outside.squared_distance == F(1, 4) and outside.point == OneRootPoint2.rational(24, 5)
    behind = native_nearest(domain.native, _point(17, 5))
    assert behind.ambiguous and behind.squared_distance == 16
    _assert_native_matches_oracle(domain, _point(F(49, 2), 5))
    _assert_native_matches_oracle(domain, _point(17, 5))


def _circle_extremes(support: OffsetCircleSupport) -> list[tuple[Fraction, Fraction]]:
    centre, radius = support.centre, support.radius
    return [(centre.x + radius, centre.y), (centre.x, centre.y + radius), (centre.x - radius, centre.y), (centre.x, centre.y - radius)]


@pytest.mark.parametrize("fixture", [*FIXTURES, *ORACLE_REFERENCE_CASES])
def test_queries_beyond_every_arc_extreme_match_the_oracle(fixture: str) -> None:
    """The certified piece box must contain each arc's axis extremes, or its pruning drops the true nearest piece."""
    domain = _domain(fixture)
    checked = 0
    for piece in domain.pieces:
        support = piece.support
        if not isinstance(support, OffsetCircleSupport):
            continue
        for x, y in _circle_extremes(support):
            extreme = OneRootPoint2.rational(x, y)
            full = piece.source == piece.target
            if not full and compare_on_arc(support.centre, support.sense, piece.source, extreme, piece.target) > 0:
                continue
            for offset in RADIAL_OFFSETS:
                scale = 1 + offset / support.radius
                _assert_native_matches_oracle(domain, _point(support.centre.x + scale * (x - support.centre.x), support.centre.y + scale * (y - support.centre.y)))
                checked += 1
    assert checked > 0


def _bracket(value: OneRootNumber) -> tuple[Fraction, Fraction]:
    """Rational ``[low, high]`` around ``value``, 2^-SUB_ULP_BRACKET_BITS wide (exact: integer square roots)."""
    if value.is_rational:
        return value.rational, value.rational
    radicand = value.radicand
    scale = 1 << SUB_ULP_BRACKET_BITS
    root_low = Fraction(math.isqrt(radicand.numerator * radicand.denominator * scale * scale), radicand.denominator * scale)
    root_high = root_low + Fraction(1, radicand.denominator * scale)
    ends = (value.rational + value.coefficient * root_low, value.rational + value.coefficient * root_high)
    return min(ends), max(ends)


@pytest.mark.parametrize("fixture", ["concave", "rich", "figure5"])
def test_queries_off_irrational_vertex_directions_by_less_than_an_ulp_match_the_oracle(fixture: str) -> None:
    """Interval enclosures of one-root vertices must be sound: a point enclosure flips these sub-ulp side decisions."""
    domain = _domain(fixture)
    checked = 0
    for piece in domain.pieces:
        support = piece.support
        vertex = piece.source
        if not isinstance(support, OffsetCircleSupport) or vertex.radicand == 0:
            continue
        x_low, x_high = _bracket(vertex.x)
        y_low, y_high = _bracket(vertex.y)
        for x, y in ((x_low, y_low), (x_high, y_high)):
            for offset in RADIAL_OFFSETS:
                scale = 1 + offset / support.radius
                _assert_native_matches_oracle(domain, _point(support.centre.x + scale * (x - support.centre.x), support.centre.y + scale * (y - support.centre.y)))
                checked += 1
    assert checked > 0
