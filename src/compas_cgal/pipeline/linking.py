"""The one link constructor between two loops (Rulings 1.1-R12, 1.1-D1, 1.1-C1-A).

`link_between(loop_k, next_loop, clockwise, next_clockwise)` is the only place a cut-plane
motion between two loops is shaped; the ladder and the next-chain choice both
call it. It returns one of:

* `SameCircleArc` -- the two loops share their circle exactly (centre and
  authored radius): an arc about its centre
  from ``loop_k``'s entry to ``next_loop``'s, in the turn direction;
* `StraightLink` -- the two entries were placed on the pair's external common
  tangent by construction (consecutive narrow stations of one chain,
  ``tangent_entries=True``): the straight chord between them is the link;
* `BeltLink` -- otherwise, the co-rotating external common tangent ("belt"):
  an overrun arc on ``loop_k`` from its entry to the tangent point, the
  tangent segment, and an approach arc on ``next_loop`` to its entry. Every
  junction is G1 at any radius difference;
* `NoLink` -- the loops turn opposite ways, or the circles are nested,
  concentric or internally tangent, so no co-rotating link exists.

The belt's EXISTENCE is decided exactly over the rationals of the doubles,
``(rho_a - rho_b)^2 < |c_b - c_a|^2``. Whether the belt touches a loop EXACTLY
at its authored entry (a wall-parallel belt between two loops tangent to the
same wall offset, at an axis phase) is decided exactly too, and then that end
needs no arc. Otherwise the tangent points are a double CONSTRUCTION after the
existence decision; they place geometry and decide nothing (the pieces'
removal and admissibility are decided on what is emitted).

`StraightLink` is ruled by the construction that placed its entries, not
re-derived: a chord between two exactly incident entries cannot be decided
tangent exactly (Ruling 1.1-C17).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


@dataclass(frozen=True)
class LinePiece:
    """A straight piece of a link.

    Attributes:
        start: Start.
        end: End.
    """

    start: Point2[WorldXY]
    end: Point2[WorldXY]


@dataclass(frozen=True)
class ArcPiece:
    """A circular piece of a link, travelled in the link's turn direction.

    Attributes:
        centre: Arc centre.
        start: Start, on the arc's circle.
        end: End.
    """

    centre: Point2[WorldXY]
    start: Point2[WorldXY]
    end: Point2[WorldXY]


LinkPiece = LinePiece | ArcPiece


@dataclass(frozen=True)
class SameCircleArc:
    """An arc about the shared centre of two loops.

    Attributes:
        next_loop: The loop the link arrives at.
        departure: Where it leaves (the previous loop's entry).
    """

    next_loop: IncidentLoop
    departure: Point2[WorldXY]

    @property
    def pieces(self) -> tuple[LinkPiece, ...]:
        """The arc, or nothing when the two entries coincide exactly."""
        if self.departure == self.next_loop.entry:
            return ()
        return (ArcPiece(centre=self.next_loop.centre, start=self.departure, end=self.next_loop.entry),)


@dataclass(frozen=True)
class StraightLink:
    """The chord between two entries placed on their external tangent by construction.

    Attributes:
        next_loop: The loop the link arrives at.
        departure: Where it leaves (the previous loop's entry).
    """

    next_loop: IncidentLoop
    departure: Point2[WorldXY]

    @property
    def pieces(self) -> tuple[LinkPiece, ...]:
        """The chord, or nothing when the two entries coincide exactly."""
        if self.departure == self.next_loop.entry:
            return ()
        return (LinePiece(start=self.departure, end=self.next_loop.entry),)


@dataclass(frozen=True)
class BeltLink:
    """Overrun arc, external common tangent, approach arc.

    Attributes:
        next_loop: The loop the link arrives at.
        departure: Where it leaves (the previous loop's entry).
        departure_centre: The previous loop's centre.
        tangent_a: The tangent point on the previous loop's circle.
        tangent_b: The tangent point on the next loop's circle.
    """

    next_loop: IncidentLoop
    departure: Point2[WorldXY]
    departure_centre: Point2[WorldXY]
    tangent_a: Point2[WorldXY]
    tangent_b: Point2[WorldXY]

    @property
    def pieces(self) -> tuple[LinkPiece, ...]:
        """The non-empty pieces, in travel order."""
        pieces: list[LinkPiece] = []
        if self.departure != self.tangent_a:
            pieces.append(ArcPiece(centre=self.departure_centre, start=self.departure, end=self.tangent_a))
        if self.tangent_a != self.tangent_b:
            pieces.append(LinePiece(start=self.tangent_a, end=self.tangent_b))
        if self.tangent_b != self.next_loop.entry:
            pieces.append(ArcPiece(centre=self.next_loop.centre, start=self.tangent_b, end=self.next_loop.entry))
        return tuple(pieces)

    @property
    def tangent_length(self) -> float:
        """Length of the tangent segment (reporting and next-chain ordering only)."""
        return math.hypot(self.tangent_b.x - self.tangent_a.x, self.tangent_b.y - self.tangent_a.y)


@dataclass(frozen=True)
class NoLink:
    """No co-rotating link exists: opposite turn directions, or nested, concentric or internally tangent circles."""


Link = SameCircleArc | BeltLink | StraightLink


def belt_exists(loop_a: IncidentLoop, loop_b: IncidentLoop) -> bool:
    """Whether the two circles have an external common tangent: ``(rho_a - rho_b)^2 < |c_b - c_a|^2``, exactly.

    Args:
        loop_a: First loop.
        loop_b: Second loop.

    Returns:
        True iff neither circle lies inside the other (nor touches it from inside).
    """
    rise = Fraction(loop_b.radius.value) - Fraction(loop_a.radius.value)
    dx = Fraction(loop_b.centre.x) - Fraction(loop_a.centre.x)
    dy = Fraction(loop_b.centre.y) - Fraction(loop_a.centre.y)
    return rise * rise < dx * dx + dy * dy


def _belt_tangent_points(loop_a: IncidentLoop, loop_b: IncidentLoop, clockwise: bool) -> tuple[Point2[WorldXY], Point2[WorldXY]]:
    """The co-rotating external tangent's touch points (a construction; existence decided by `belt_exists`).

    The common normal is ``n = along * u + side * lateral * v`` with ``u`` the
    unit centre difference, ``v`` its counter-clockwise perpendicular,
    ``along = -(rho_b - rho_a) / d`` and ``lateral = sqrt(1 - along^2)``. WHICH
    side is the co-rotating belt is not decided on doubles (Ruling
    1.1-R1-exact): the travel direction at the touch points is ``perp(n)``, and
    ``(t_b - t_a) . perp(n) = (c_b - c_a) . perp(n) = d * side * lateral``
    for clockwise travel (negated for counter-clockwise), because
    ``n . perp(n) = 0`` and ``cross(u, v) = 1``. `belt_exists` decides
    ``lateral > 0`` exactly (``(rho_a - rho_b)^2 < d^2``), so the belt is
    ``side = +1`` for clockwise travel and ``side = -1`` for counter-clockwise,
    with no comparison at all.
    """
    ax, ay, ra = loop_a.centre.x, loop_a.centre.y, loop_a.radius.value
    bx, by, rb = loop_b.centre.x, loop_b.centre.y, loop_b.radius.value
    dx, dy = bx - ax, by - ay
    distance = math.hypot(dx, dy)
    ux, uy = dx / distance, dy / distance
    vx, vy = -uy, ux
    along = -(rb - ra) / distance
    lateral = math.sqrt(max(0.0, 1.0 - along * along))
    side = 1.0 if clockwise else -1.0
    nx, ny = along * ux + side * lateral * vx, along * uy + side * lateral * vy
    return Point2[WorldXY].build(ax + ra * nx, ay + ra * ny), Point2[WorldXY].build(bx + rb * nx, by + rb * ny)


def _exact(point: Point2[WorldXY]) -> tuple[Fraction, Fraction]:
    return Fraction(point.x), Fraction(point.y)


def _touches_at_entry(on: IncidentLoop, other: IncidentLoop, clockwise: bool, leaving: bool) -> tuple[Fraction, Fraction] | None:
    """The other touch point when the belt touches *on* exactly at its entry, else None (exact).

    The tangent line at *on*'s entry ``e`` (normal ``e - c``) is the co-rotating
    external common tangent iff *other*'s centre lies on *on*'s side of it at
    distance exactly ``rho_other``,
    ``(c_other - e) . (e - c) == -rho_other * rho``, and the travel direction at
    ``e`` runs toward *other*'s touch point when *leaving* (away from it when
    arriving). Every quantity is a rational of the doubles, so the predicate is
    exact; the other touch point ``c_other + rho_other (e - c) / rho`` is then
    rational too.

    Args:
        on: The loop whose entry is tested.
        other: The loop at the other end of the belt.
        clockwise: The common turn direction.
        leaving: True when the belt leaves *on* (else it arrives at *on*).

    Returns:
        The other touch point, exactly, or None.
    """
    (ex, ey), (cx, cy) = _exact(on.entry), _exact(on.centre)
    (ox, oy), rho, rho_other = _exact(other.centre), Fraction(on.radius.value), Fraction(other.radius.value)
    nx, ny = ex - cx, ey - cy
    if (ox - ex) * nx + (oy - ey) * ny != -rho_other * rho:
        return None
    tx, ty = (ny, -nx) if clockwise else (-ny, nx)
    touch = (ox + rho_other * nx / rho, oy + rho_other * ny / rho)
    along = (touch[0] - ex) * tx + (touch[1] - ey) * ty
    if (along > 0) != leaving:
        return None
    return touch


def _as_point(exact: tuple[Fraction, Fraction], snap_to: Point2[WorldXY]) -> Point2[WorldXY]:
    """*snap_to* when it IS the exact point, else the double rounding of the exact point (placement)."""
    if _exact(snap_to) == exact:
        return snap_to
    return Point2[WorldXY].build(float(exact[0]), float(exact[1]))


def link_between(loop_k: IncidentLoop, next_loop: IncidentLoop, clockwise: bool, next_clockwise: bool, *, tangent_entries: bool = False) -> Link | NoLink:
    """The cut-plane link from *loop_k*'s entry to *next_loop*'s.

    Every link here is co-rotating: loops turning opposite ways get `NoLink`.
    A same-circle arc needs the SAME circle, exactly: equal centre and equal
    authored radius (Ruling 1.1-I5; two loops on one centre can differ in
    authored radius through their phase-dependent chart quanta, and then no
    arc reaches the next entry). Concentric circles of different radius have
    no external tangent either, so they get `NoLink`.

    Args:
        loop_k: The loop the tool is on (it has completed its turn and is back
            at its entry).
        next_loop: The loop to reach.
        clockwise: *loop_k*'s turn direction.
        next_clockwise: *next_loop*'s turn direction.
        tangent_entries: The two entries were placed on the pair's external
            common tangent by construction (consecutive narrow stations of one
            chain); the chord is then the link.

    Returns:
        The link, or `NoLink`.
    """
    if clockwise != next_clockwise:
        return NoLink()
    if loop_k.centre == next_loop.centre and loop_k.radius == next_loop.radius:
        return SameCircleArc(next_loop=next_loop, departure=loop_k.entry)
    if tangent_entries:
        return StraightLink(next_loop=next_loop, departure=loop_k.entry)
    if not belt_exists(loop_k, next_loop):
        return NoLink()
    leaving = _touches_at_entry(loop_k, next_loop, clockwise, leaving=True)
    arriving = _touches_at_entry(next_loop, loop_k, clockwise, leaving=False)
    if leaving is not None:
        tangent_a, tangent_b = loop_k.entry, _as_point(leaving, next_loop.entry)
    elif arriving is not None:
        tangent_a, tangent_b = _as_point(arriving, loop_k.entry), next_loop.entry
    else:
        tangent_a, tangent_b = _belt_tangent_points(loop_k, next_loop, clockwise)
    return BeltLink(next_loop=next_loop, departure=loop_k.entry, departure_centre=loop_k.centre, tangent_a=tangent_a, tangent_b=tangent_b)
