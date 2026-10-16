"""Inter-loop links (Rulings 1.1-R12, 1.1-D1, 1.1-C1-A, 1.1-C2)."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from compas.geometry import Polygon
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.quality import TANGENT_CONTINUITY_SLACK

from compas_cgal.engagement import _operation_removal
from compas_cgal.engagement import _subtract_operation
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import arc_operation
from compas_cgal.pipeline.linking import BeltLink
from compas_cgal.pipeline.linking import LinePiece
from compas_cgal.pipeline.linking import NoLink
from compas_cgal.pipeline.linking import SameCircleArc
from compas_cgal.pipeline.linking import link_between
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

TOOL_RADIUS = 1.0
CUT = CutZ.build(0.0)
# Test-side floor on the belt length below which its direction is not measured:
# a sub-micron segment's direction is rounding noise (measurement, not a decision).
MIN_BELT_FOR_DIRECTION_MM = 1e-3

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])


def _loop(cx: float, cy: float, rho: float, phase: float = 0.0) -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(cx, cy), GuideRadius.build(rho), Radian(phase))


@pytest.mark.parametrize(("phase", "turn"), [(0.3, 1.6), (1.1, 2.9), (2.2, 3.7), (4.4, 0.7)])
@pytest.mark.parametrize("clockwise", [True, False])
def test_a_same_circle_link_is_an_air_arc_and_removes_nothing(phase: float, turn: float, clockwise: bool) -> None:
    """C2-B-minimal: after its loop, an arc on the loop's own circle is decided NOTHING and depletes nothing.

    The arc leaves the loop's exact entry and ends at a double construction on
    the circle (as a belt's overrun arc does), ulps off the authored radius.
    """
    stock = Stock(RECT_12X8)
    loop = _loop(0.4, -0.3, 2.0, phase)
    stock.subtract_circle_sweep_local(loop.centre.x, loop.centre.y, loop.radius.value, TOOL_RADIUS)
    assert stock.circle_sweep_removal(loop.centre.x, loop.centre.y, loop.radius.value, TOOL_RADIUS) is SweepRemoval.NOTHING
    angle = loop.phase + turn
    end = Point2[WorldXY].build(loop.centre.x + loop.radius.value * math.cos(angle), loop.centre.y + loop.radius.value * math.sin(angle))
    arc = arc_operation(loop.centre, loop.entry, end, clockwise, CUT, 0, OperationType.LINK)
    assert _operation_removal(stock, arc, TOOL_RADIUS) is SweepRemoval.NOTHING
    before = stock.clone()
    _subtract_operation(stock, arc, TOOL_RADIUS)
    assert stock.exactly_equals(before)


def _unit(dx: float, dy: float) -> tuple[float, float]:
    length = math.hypot(dx, dy)
    return dx / length, dy / length


def _travel(centre: Point2[WorldXY], point: Point2[WorldXY], clockwise: bool) -> tuple[float, float]:
    nx, ny = _unit(point.x - centre.x, point.y - centre.y)
    return (ny, -nx) if clockwise else (-ny, nx)


@settings(max_examples=200, deadline=None)
@given(
    ax=st.floats(-5.0, 5.0),
    ay=st.floats(-5.0, 5.0),
    ra=st.floats(1.05, 4.0),
    pa=st.floats(-math.pi, math.pi),
    dx=st.floats(-6.0, 6.0),
    dy=st.floats(-6.0, 6.0),
    rb=st.floats(1.05, 4.0),
    pb=st.floats(-math.pi, math.pi),
    clockwise=st.booleans(),
)
def test_the_belt_link_leaves_and_meets_both_loops_tangentially(ax: float, ay: float, ra: float, pa: float, dx: float, dy: float, rb: float, pb: float, clockwise: bool) -> None:
    a, b = _loop(ax, ay, ra, pa), _loop(ax + dx, ay + dy, rb, pb)
    link = link_between(a, b, clockwise, clockwise)
    assume(isinstance(link, BeltLink))
    assert isinstance(link, BeltLink)
    assume(link.tangent_length > MIN_BELT_FOR_DIRECTION_MM)
    u = _unit(link.tangent_b.x - link.tangent_a.x, link.tangent_b.y - link.tangent_a.y)
    for centre, point in ((a.centre, link.tangent_a), (b.centre, link.tangent_b)):
        tx, ty = _travel(centre, point, clockwise)
        assert u[0] * tx + u[1] * ty >= 1.0 - TANGENT_CONTINUITY_SLACK
    pieces = link.pieces
    assert pieces[0].start == a.entry and pieces[-1].end == b.entry
    assert all(earlier.end == later.start for earlier, later in zip(pieces, pieces[1:]))
    assert link.next_loop is b and link.next_loop.certifies_incidence()


def test_a_wall_parallel_belt_leaves_and_meets_exactly_at_the_entries() -> None:
    """Both loops touch one wall offset at their axis-phase entries: decided exactly, the belt is the chord.

    The double construction of the tangent points lands ulps off both entries
    here (a zero-length-arc hazard); the exact decision does not.
    """
    a = _loop(0.0, 0.0, 3.3, math.pi / 2.0)
    b = _loop(2.2, 1.1, a.radius.value - 1.1, math.pi / 2.0)
    assert a.entry.y == b.entry.y, "fixture precondition: both entries on one wall offset"
    link = link_between(a, b, True, True)
    assert isinstance(link, BeltLink)
    assert link.pieces == (LinePiece(start=a.entry, end=b.entry),)


def test_contained_circles_have_no_link() -> None:
    """No external common tangent for nested or internally tangent circles, decided exactly."""
    outer = _loop(0.0, 0.0, 3.0)
    assert isinstance(link_between(outer, _loop(0.5, 0.0, 1.5), True, True), NoLink)
    assert isinstance(link_between(outer, _loop(1.0, 0.0, 2.0), True, True), NoLink), "internally tangent: (3-2)^2 == 1^2"
    assert isinstance(link_between(outer, _loop(1.0, 0.0, 2.0078125), True, True), BeltLink), "one dyadic step past internal tangency"


def test_a_zero_length_same_circle_link_is_not_emitted() -> None:
    """The stage-0 zero-length row: a link from a loop's entry to itself has no pieces at all."""
    loop = _loop(0.25, -0.5, 2.0, 1.3)
    link = link_between(loop, loop, True, True)
    assert isinstance(link, SameCircleArc)
    assert link.pieces == ()


def test_a_belt_touching_both_entries_off_axis_is_the_chord_decided_exactly() -> None:
    """Both entries lie on one common tangent at the (3, 4, 5) chart direction: no arc, not even an ulp one."""
    a = _loop(0.0, 0.0, 5.0, math.atan2(4.0, 3.0))
    b = _loop(5.5, -1.0, 2.5, math.atan2(4.0, 3.0))
    assert (a.entry, b.entry) == (Point2[WorldXY].build(3.0, 4.0), Point2[WorldXY].build(7.0, 1.0))
    link = link_between(a, b, True, True)
    assert isinstance(link, BeltLink)
    assert link.pieces == (LinePiece(start=a.entry, end=b.entry),)


def test_concentric_circles_of_different_radius_are_not_a_same_circle_link() -> None:
    """I5: a shared centre is not a shared circle; with different radii no arc reaches the next entry, and no belt exists."""
    a = _loop(0.25, -0.5, 2.0, 0.0)
    b = _loop(0.25, -0.5, 2.5, 1.3)
    assert a.centre == b.centre and a.radius != b.radius
    assert isinstance(link_between(a, b, True, True), NoLink)


def test_a_same_circle_link_needs_one_turn_direction() -> None:
    """I5: counter-rotating loops share no co-rotating link (as the chain-ordered prototype required)."""
    a = _loop(0.25, -0.5, 2.0, 0.0)
    b = _loop(0.25, -0.5, 2.0, math.pi / 2.0)
    assert (a.centre, a.radius) == (b.centre, b.radius)
    assert isinstance(link_between(a, b, True, True), SameCircleArc)
    assert isinstance(link_between(a, b, True, False), NoLink)
    assert isinstance(link_between(a, _loop(5.0, 0.0, 1.5), True, False), NoLink)


@pytest.mark.parametrize("clockwise", [True, False])
def test_the_belt_side_is_the_exact_co_rotating_one_near_internal_tangency(clockwise: bool) -> None:
    """R1-exact: the side is fixed by the turn direction (lateral > 0 is decided exactly), never by a double dot.

    One dyadic step past internal tangency the belt is short; its travel at the
    departure still runs from A towards B, read on the exact emitted points.
    """
    a, b = _loop(0.0, 0.0, 3.0), _loop(1.0, 0.0, 2.0078125)
    link = link_between(a, b, clockwise, clockwise)
    assert isinstance(link, BeltLink)
    nx, ny = Fraction(link.tangent_a.x) - Fraction(a.centre.x), Fraction(link.tangent_a.y) - Fraction(a.centre.y)
    tx, ty = (ny, -nx) if clockwise else (-ny, nx)
    assert (Fraction(link.tangent_b.x) - Fraction(link.tangent_a.x)) * tx + (Fraction(link.tangent_b.y) - Fraction(link.tangent_a.y)) * ty > 0
