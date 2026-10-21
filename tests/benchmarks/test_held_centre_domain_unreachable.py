"""Unreachable attribution, vertex disks and input guards of the exact Held centre domain Ω*.

Attribution rule (docs/held_exact_centre_domain.md, "Unreachable attribution"):
a boundary point b is reached iff its inward centre b + r·n̂ lies in Ω*, so
the reached part of a primitive is the contact-map image of its owned ∂Ω*
pieces plus the points the vertex disks touch. The complement intervals are
the unreachable boundary, grouped into stretches along ∂P and attributed.

Fixtures are Pythagorean by construction (integers, r = 1), so every expected
interval end is an exact rational point.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_centre_domain_supports import ClockwiseDeclaredBoundaryError
from benchmarks.held_centre_domain_supports import NonSimpleDeclaredBoundaryError
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import EmptyCentreDomainError
from benchmarks.held_exact_centre_domain import PointLocation
from benchmarks.held_exact_centre_domain import UnmatchedCentreDomainPieceError
from benchmarks.held_exact_centre_domain import UnreachableAttribution
from benchmarks.held_exact_centre_domain import UnreachableCause
from benchmarks.held_exact_centre_domain import UnreachableInterval
from benchmarks.held_exact_centre_domain import _matched_piece
from benchmarks.held_centre_domain_native import NativeCircle
from benchmarks.held_centre_domain_native import NativePiece
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import _Pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_exact_centre_domain import tight_convex_arc_pocket

F = Fraction
R = F(1)


def _point(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(F(x), F(y))


def _line(index: int) -> BoundaryFeature:
    return BoundaryFeature(BoundaryFeatureKind.LINE, index)


def _arc(index: int) -> BoundaryFeature:
    return BoundaryFeature(BoundaryFeatureKind.ARC, index)


def _junction(index: int) -> BoundaryFeature:
    return BoundaryFeature(BoundaryFeatureKind.JUNCTION, index)


def _interval(
    primitive: BoundaryFeature, start: tuple[Fraction | int, Fraction | int], end: tuple[Fraction | int, Fraction | int], closed: tuple[bool, bool]
) -> UnreachableInterval:
    return UnreachableInterval(primitive, OneRootPoint2.rational(*start), OneRootPoint2.rational(*end), closed[0], closed[1])


def _same_intervals(actual: tuple[UnreachableInterval, ...], expected: list[UnreachableInterval]) -> bool:
    return len(actual) == len(expected) and all(
        (got.primitive, got.start_closed, got.end_closed) == (want.primitive, want.start_closed, want.end_closed) and got.start == want.start and got.end == want.end
        for got, want in zip(actual, expected)
    )


def _parts(domain: CentreDomain, cause: UnreachableCause) -> list[UnreachableAttribution]:
    return [part for part in domain.unreachable if part.cause is cause]


def dead_end_slot_pocket(width: Fraction) -> DeclaredBoundary:
    """20 x 20 square with a 10 mm dead-end slot of ``width`` centred on x = 10, capped by a semicircle of ρ = width / 2.

    Primitives: 0 bottom, 1 right, 2 top (right part), 3 slot right wall, 4 cap,
    5 slot left wall, 6 top (left part), 7 left. Junctions 3 and 6 are the
    reflex slot mouth.
    """
    half = width / 2
    pocket = _Pocket(_point(0, 0)).line(0, 20).line(1, 20).line(2, 10 - half).line(1, 10)
    return pocket.arc(_point(10 - half, 30), F(0), 1, 1).line(3, 10).line(2, 10 - half).line(3, 20).close(f"dead_end_slot_{width}")


def cusp_pocket() -> DeclaredBoundary:
    """A passage of exactly 2r between two concave bumps, closed on the left by a wall r beyond the pinch.

    Primitives: 0 left wall x = 9 (down), 1 bottom bump (CW about (10, 11/4),
    ρ = 5/4, apex (10, 4)), 2 line down, 3 bottom, 4 right, 5 top, 6 line down,
    7 top bump (CW about (10, 29/4), apex (10, 6)). The two bump offsets
    (radius 9/4) touch at w = (10, 5), which the wall offset x = 10 also passes
    through: ∂Ω* turns back on itself there (a π-turn cusp), and the wall is
    reached at the single point (9, 5).
    """
    pocket = _Pocket(_point(9, F(13, 2))).line(3, 3).arc(_point(F(45, 4), F(11, 4)), F(1, 2), 0, -1).line(3, F(11, 4))
    pocket = pocket.line(0, F(75, 4)).line(1, 10).line(2, F(75, 4)).line(3, F(11, 4))
    return pocket.arc(_point(9, F(13, 2)), F(0), 3, -1).close("cusp")


def square_pocket(side: int) -> DeclaredBoundary:
    return _Pocket(_point(0, 0)).line(0, side).line(1, side).line(2, side).line(3, side).close(f"square_{side}")


# ---------------------------------------------------------------------------
# I1: attribution by unreached boundary
# ---------------------------------------------------------------------------


def test_dead_end_slot_exactly_2r_wide_is_attributed_to_a_neck() -> None:
    domain = CentreDomain.build(dead_end_slot_pocket(F(2)), R)
    assert domain.locate(_point(10, 25)) is PointLocation.EXTERIOR
    assert domain.locate(_point(10, 30)) is PointLocation.EXTERIOR
    (neck,) = _parts(domain, UnreachableCause.NECK)
    assert neck.owners == frozenset({_line(3), _arc(4), _line(5)})
    expected = [
        _interval(_line(3), (11, 20), (11, 30), (False, True)),
        _interval(_arc(4), (11, 30), (9, 30), (True, True)),
        _interval(_line(5), (9, 30), (9, 20), (True, False)),
    ]
    assert _same_intervals(neck.intervals, expected)
    corners = {part.owners for part in _parts(domain, UnreachableCause.CONVEX_KINK_CUSP)}
    assert corners == {frozenset({_junction(junction)}) for junction in (0, 1, 2, 7)}
    assert not _parts(domain, UnreachableCause.TIGHT_CONVEX_ARC)


def test_dead_end_slot_wider_than_2r_is_reached() -> None:
    domain = CentreDomain.build(dead_end_slot_pocket(F(5, 2)), R)
    assert domain.locate(_point(10, 25)) is PointLocation.INTERIOR
    assert not _parts(domain, UnreachableCause.NECK)
    corners = {part.owners for part in _parts(domain, UnreachableCause.CONVEX_KINK_CUSP)}
    assert corners == {frozenset({_junction(junction)}) for junction in (0, 1, 2, 7)}


def test_square_corner_cusp_intervals_are_exact() -> None:
    domain = CentreDomain.build(l_shape_pocket(), R)
    (corner,) = [part for part in domain.unreachable if part.owners == frozenset({_junction(0)})]
    assert corner.cause is UnreachableCause.CONVEX_KINK_CUSP
    expected = [_interval(_line(0), (0, 0), (1, 0), (True, False)), _interval(_line(5), (0, 1), (0, 0), (False, True))]
    # Stretches are reported in boundary order from the first primitive; this one wraps 5 -> 0.
    assert _same_intervals(tuple(sorted(corner.intervals, key=lambda interval: interval.primitive.index)), expected)
    assert {part.owners for part in domain.unreachable} == {frozenset({_junction(junction)}) for junction in (0, 1, 2, 4, 5)}


def test_tight_convex_arc_part_spans_the_arc_and_the_line_tails() -> None:
    domain = CentreDomain.build(tight_convex_arc_pocket(), R)
    (tight,) = _parts(domain, UnreachableCause.TIGHT_CONVEX_ARC)
    expected = [
        _interval(_line(1), (20, 19), (20, F(39, 2)), (False, True)),
        _interval(_arc(2), (20, F(39, 2)), (F(39, 2), 20), (True, True)),
        _interval(_line(3), (F(39, 2), 20), (19, 20), (True, False)),
    ]
    assert _same_intervals(tight.intervals, expected)


def test_equal_radius_cap_is_exempt_only_when_its_centre_is_in_the_domain() -> None:
    slot = CentreDomain.build(dead_end_slot_pocket(F(2)), R)
    assert slot.locate(_point(10, 30)) is PointLocation.EXTERIOR
    (neck,) = _parts(slot, UnreachableCause.NECK)
    assert _arc(4) in {interval.primitive for interval in neck.intervals}


def test_every_unreached_interval_is_disjoint_from_the_reached_contacts() -> None:
    """The inward centre of an interior point of every unreached interval is outside Ω* (the lemma, exactly)."""
    domain = CentreDomain.build(dead_end_slot_pocket(F(2)), R)
    for part in domain.unreachable:
        for interval in part.intervals:
            if interval.primitive.kind is not BoundaryFeatureKind.LINE or not interval.start.x.is_rational or not interval.end.x.is_rational:
                continue
            middle = _point((interval.start.x.rational + interval.end.x.rational) / 2, (interval.start.y.rational + interval.end.y.rational) / 2)
            support = domain.supports.support_of(interval.primitive)
            assert isinstance(support, OffsetLineSupport)
            inward = _point(middle.x + R * support.normal[0], middle.y + R * support.normal[1])
            assert domain.locate(inward) is PointLocation.EXTERIOR


# ---------------------------------------------------------------------------
# M1: π-turn cusps of Ω* get a vertex disk; a vertex disk reaches a point
# ---------------------------------------------------------------------------


def test_pi_turn_cusp_of_the_domain_gets_a_vertex_disk() -> None:
    domain = CentreDomain.build(dead_end_slot_pocket(F(2)), R)
    assert any(disk.centre == OneRootPoint2.rational(10, 20) for disk in domain.vertex_disks)


def test_vertex_disk_reaches_an_isolated_boundary_point() -> None:
    domain = CentreDomain.build(cusp_pocket(), R)
    assert domain.locate(_point(10, 5)) is PointLocation.BOUNDARY
    assert any(disk.centre == OneRootPoint2.rational(10, 5) for disk in domain.vertex_disks)
    wall = sorted(
        (interval for part in domain.unreachable for interval in part.intervals if interval.primitive == _line(0)),
        key=lambda interval: -interval.start.y.rational,
    )
    expected = [_interval(_line(0), (9, F(13, 2)), (9, 5), (True, False)), _interval(_line(0), (9, 5), (9, F(7, 2)), (False, True))]
    assert _same_intervals(tuple(wall), expected)


# ---------------------------------------------------------------------------
# M5: the multi-owner sense check runs over every candidate
# ---------------------------------------------------------------------------


def test_coincident_support_of_opposite_sense_does_not_block_the_match() -> None:
    centre = _point(10, 10)
    offset_arc = OffsetCircleSupport(_arc(7), centre, ExactMillimetre(R), ExactMillimetre(2 * R), 1, _point(11, 10), _point(10, 11))
    junction = OffsetCircleSupport(_junction(3), centre, ExactMillimetre(R), ExactMillimetre(F(0)), -1, _point(10, 9), _point(9, 10))
    native = NativePiece(NativeCircle(F(10), F(10), F(1)), -1, OneRootPoint2.rational(10, 9), OneRootPoint2.rational(9, 10))
    piece = _matched_piece(native, (0, 0, 0), ({}, {(F(10), F(10), F(1)): [offset_arc, junction]}))
    assert piece.owners == frozenset({_junction(3)})
    assert piece.support == junction
    with pytest.raises(UnmatchedCentreDomainPieceError):
        _matched_piece(native, (0, 0, 0), ({}, {(F(10), F(10), F(1)): [offset_arc]}))


# ---------------------------------------------------------------------------
# I2 / I4 stopgaps: named errors for input the native construction mishandles
# ---------------------------------------------------------------------------


def clockwise_l_shape_pocket() -> DeclaredBoundary:
    """The L-shape of the contract tests traversed clockwise."""
    return _Pocket(_point(0, 0)).line(1, 20).line(0, 10).line(3, 10).line(0, 10).line(3, 10).line(2, 20).close("l_shape_clockwise")


def bow_tie_pocket() -> DeclaredBoundary:
    """(0,0) → (20,0) → (20,10) → (10,10) → (10,-10) → (0,-10): the fourth edge crosses the first at (10, 0)."""
    return _Pocket(_point(0, 0)).line(0, 20).line(1, 10).line(2, 10).line(3, 20).line(2, 10).line(1, 10).close("bow_tie")


def touching_squares_pocket() -> DeclaredBoundary:
    """Two 10 x 10 squares that share only the corner (10, 10), visited twice by the cycle (a pinch)."""
    pocket = _Pocket(_point(0, 0)).line(0, 10).line(1, 10).line(0, 10).line(1, 10)
    return pocket.line(2, 10).line(3, 10).line(2, 10).line(3, 10).close("touching_squares")


def test_clockwise_cycle_is_a_named_error() -> None:
    with pytest.raises(ClockwiseDeclaredBoundaryError):
        CentreDomain.build(clockwise_l_shape_pocket(), R)


def test_crossing_cycle_is_a_named_error() -> None:
    with pytest.raises(NonSimpleDeclaredBoundaryError):
        CentreDomain.build(bow_tie_pocket(), R)


def test_self_touching_cycle_is_a_named_error() -> None:
    with pytest.raises(NonSimpleDeclaredBoundaryError):
        CentreDomain.build(touching_squares_pocket(), R)


def test_tool_that_fits_nowhere_is_an_empty_domain_error() -> None:
    with pytest.raises(EmptyCentreDomainError):
        CentreDomain.build(square_pocket(4), F(3))


def split_cap_stadium_pocket() -> DeclaredBoundary:
    """A 20 x 10 stadium whose right cap (ρ = 5 about (20, 5)) is declared as two co-circular quarter arcs (primitives 1, 2)."""
    pocket = _Pocket(_point(0, 0)).line(0, 20).arc(_point(25, 5), F(0), 0, 1).arc(_point(20, 10), F(0), 1, 1)
    return pocket.line(2, 20).arc(_point(0, 0), F(0), 2, 1).close("split_cap_stadium")


def test_co_circular_arcs_share_one_piece_and_are_both_reached() -> None:
    domain = CentreDomain.build(split_cap_stadium_pocket(), R)
    (cap,) = [piece for piece in domain.pieces if piece.owners & {_arc(1), _arc(2)}]
    assert cap.owners == frozenset({_arc(1), _arc(2)})
    assert not domain.unreachable
