"""Exact detection of a convex ∂Ω* corner whose engagement floor exceeds the cap (no refinement can certify it).

The polygon path's convex-corner lemma (``docs/held_connector_engagement.md``,
"Convex corners") on Ω*: at a convex vertex V of the centre domain with
interior angle α, the motion that first reaches V engages at least π − α of rim
on the actual stock, whatever was cut before, provided every tool-centre point
cut before within 2r of V lies in V's tangent wedge T. With exterior angle
β = π − α above the cap, no tier-1 / 1b certificate exists and no midpoint
circle can create one, so the connector repair refuses by name instead of
halving circles toward the vertex (user ruling 2026-10-22 for the polygon path).

Scope: V joins two offset LINE pieces (the polygon lemma's straight wedge,
verbatim). :func:`sharp_vertex` answers only when everything is decided
exactly; otherwise it returns ``None`` and the repair's depth budget names
the failure (the conservative direction).

The hypothesis is checked on what has been cut, exactly:

- every emitted circle either lies with its junction line in T (its disk:
  ``n̂_k·c - level_k >= ρ`` for both lines; its site ``q`` in T; the junction
  ``J -> q`` joins a disk point to ``q``, so it is in the convex T), or the
  disk of radius ``max(ρ, |q - c|)`` about ``c``, which holds the circle and its
  junction, is provably farther than ``2r`` from V;
- every other ∂Ω* piece lies in T (a line: both ends; an arc: its whole
  circle) or is provably farther than ``2r`` from V (a rational box).

The second gives Ω* ∩ B(V, 2r) ⊂ T: B(V, 2r) minus T is connected, ∂Ω* does
not meet it, and Ω* is locally T at a line-line convex vertex, so B minus T
lies outside Ω*. Earlier ∂Ω* connectors are therefore in T near V. "Provably
farther" uses rational interval boxes and exact comparisons (the hull radius
``max(ρ, √F)`` is compared without a square root) and only ever claims "far"
when it is true.

Precondition, checked exactly: the connector is the FIRST arrival at V. The
lemma's "every earlier cut keeps δ > 0 from V" covers circles (a disk in a
wedge misses the apex) and junction lines (``q ≠ V``), but an earlier ∂Ω*
connector in T could pass through V itself. The generator rules that out: it
emits circles in exact ∂Ω* site-key order, one counter-clockwise lap from
circle 0, and repair midpoints lie strictly between their neighbours. So
:func:`sharp_vertex` requires the emitted sites to be non-decreasing along the
lap that starts at circle 0's site, with V strictly after the predecessor's
site on that lap (then no earlier connector, each of which spans the lap
between two consecutive emitted sites, reaches V). A breach is a contract
violation of the caller and raises :class:`OmegaCornerFirstArrivalError`.
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_dyadic_interval import FRACTION_BITS
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connectors import _chain
from benchmarks.held_omega_refinement import OmegaRefinementError
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_omega_sites import SiteKey
from benchmarks.held_one_root import OneRootPoint2

_Rational2 = tuple[Fraction, Fraction]


class OmegaCornerFirstArrivalError(OmegaRefinementError):
    """The connector is not the first arrival at the convex vertex: the emitted sites do not lap ∂Ω* once in order before it.

    ``piece`` is the outgoing piece (the vertex is its source).
    """

    def __init__(self, piece: int, detail: str) -> None:
        super().__init__(f"Convex vertex at the source of piece {piece}: {detail}")
        self.piece = piece


class UnknownOmegaSupportError(OmegaRefinementError):
    """A ∂Ω* piece has a support that is neither an offset line nor an offset circle."""


def exterior_angle_exceeds(incoming: OffsetLineSupport, outgoing: OffsetLineSupport, cos_sq_half_cap: Fraction) -> bool:
    """Exact: the turn from ``incoming`` to ``outgoing`` is a convex (left) turn whose exterior angle β exceeds the cap.

    Both travel directions are exactly unit, so ``cos β = u0·u1``; the cap enters
    as its surrogate ``C = cos²(cap/2)``, ``cos(cap) = 2C - 1``, and β > cap iff
    ``u0·u1 < 2C - 1`` (both angles in [0, π]).
    """
    (ax, ay), (bx, by) = incoming.direction, outgoing.direction
    if ax * by - ay * bx <= 0:
        return False
    return ax * bx + ay * by < 2 * cos_sq_half_cap - 1


def _bounds(box: Interval) -> tuple[Fraction, Fraction]:
    return Fraction(box.lo, 1 << FRACTION_BITS), Fraction(box.hi, 1 << FRACTION_BITS)


class _Wedge:
    """The tangent wedge T at a line-line vertex V and the exact tests against it."""

    def __init__(self, incoming: OffsetLineSupport, outgoing: OffsetLineSupport, vertex: _Rational2, tool_radius: Fraction) -> None:
        self.lines = (incoming, outgoing)
        self.vertex = vertex
        self.reach = 2 * tool_radius

    def contains(self, point: OneRootPoint2) -> bool:
        return all((line.normal[0] * point.x + line.normal[1] * point.y).compare(line.level) >= 0 for line in self.lines)

    def holds_disk(self, centre: _Rational2, radius: Fraction) -> bool:
        return all(line.normal[0] * centre[0] + line.normal[1] * centre[1] - line.level >= radius for line in self.lines)

    def box_far(self, x: tuple[Fraction, Fraction], y: tuple[Fraction, Fraction]) -> bool:
        """The axis box ``x × y`` is farther than ``2r`` from V (exact, rational)."""
        vx, vy = self.vertex
        dx = max(x[0] - vx, Fraction(0), vx - x[1])
        dy = max(y[0] - vy, Fraction(0), vy - y[1])
        return dx * dx + dy * dy > self.reach * self.reach

    def hull_far(self, centre: _Rational2, radius: Fraction, squared_reach: Fraction) -> bool:
        """The disk of radius ``max(radius, √squared_reach)`` about ``centre`` is farther than ``2r`` from V, exactly.

        With ``D² = |centre - V|²`` and ``R = 2r``: if ``radius² >= F`` it is
        ``D² > (R + radius)²``; otherwise ``D² > (R + √F)²`` iff
        ``X = D² - R² - F > 0`` and ``X² > 4R²F`` (both sides positive).
        """
        dx, dy = centre[0] - self.vertex[0], centre[1] - self.vertex[1]
        squared = dx * dx + dy * dy
        if radius * radius >= squared_reach:
            return squared > (radius + self.reach) ** 2
        excess = squared - self.reach * self.reach - squared_reach
        return excess > 0 and excess * excess > 4 * self.reach * self.reach * squared_reach


def _circle_clear(wedge: _Wedge, circle: OmegaCircle) -> bool:
    centre = (Fraction(float(circle.center.x)), Fraction(float(circle.center.y)))
    radius = Fraction(float(circle.radius.value))
    qx, qy = (_bounds(box) for box in circle.site.point_box())
    far_x = max(abs(qx[0] - centre[0]), abs(qx[1] - centre[0]))
    far_y = max(abs(qy[0] - centre[1]), abs(qy[1] - centre[1]))
    if wedge.hull_far(centre, radius, far_x * far_x + far_y * far_y):
        return True
    return wedge.holds_disk(centre, radius) and wedge.contains(circle.site.point())


def _piece_clear(wedge: _Wedge, loop: OmegaLoop, index: int) -> bool:
    piece = loop.pieces[index]
    support = piece.support
    if isinstance(support, OffsetLineSupport):
        if wedge.contains(piece.source) and wedge.contains(piece.target):
            return True
        frame = loop.frames[index]
        (sx, sy), (tx, ty) = frame.source_box, frame.target_box
        xs, ys = (_bounds(sx), _bounds(tx)), (_bounds(sy), _bounds(ty))
        return wedge.box_far((min(xs[0][0], xs[1][0]), max(xs[0][1], xs[1][1])), (min(ys[0][0], ys[1][0]), max(ys[0][1], ys[1][1])))
    if not isinstance(support, OffsetCircleSupport):
        raise UnknownOmegaSupportError(f"Piece {index} has support {type(support).__name__}.")
    centre, radius = (support.centre.x, support.centre.y), support.radius
    if wedge.holds_disk(centre, radius):
        return True
    return wedge.box_far((centre[0] - radius, centre[0] + radius), (centre[1] - radius, centre[1] + radius))


def _lap_compare(first: SiteKey, second: SiteKey, origin: SiteKey) -> int:
    """Exact three-way order along the counter-clockwise lap of ∂Ω* that starts at ``origin``."""
    first_wraps, second_wraps = first.compare(origin) < 0, second.compare(origin) < 0
    if first_wraps != second_wraps:
        return 1 if first_wraps else -1
    return first.compare(second)


def require_first_arrival(loop: OmegaLoop, emitted: Sequence[OmegaCircle], piece: int) -> None:
    """Exact check that no connector between ``emitted`` circles reaches the source vertex of ``piece``.

    Raises:
        OmegaCornerFirstArrivalError: The emitted sites are not non-decreasing
            along the lap from ``emitted[0]``, or the vertex is not strictly
            after the last emitted site on that lap.
    """
    origin = emitted[0].site.key
    keys = [circle.site.key for circle in emitted]
    for k, (first, second) in enumerate(zip(keys, keys[1:])):
        if _lap_compare(first, second, origin) > 0:
            raise OmegaCornerFirstArrivalError(piece, f"emitted circle {k + 1} lies before circle {k} on the lap from circle 0")
    vertex = OmegaSite.at_source(loop.frames[piece]).key
    if _lap_compare(vertex, keys[-1], origin) <= 0:
        raise OmegaCornerFirstArrivalError(piece, "an earlier connector already reaches the vertex")


def sharp_vertex(loop: OmegaLoop, emitted: Sequence[OmegaCircle], target: OmegaCircle, cos_sq_half_cap: Fraction, tool_radius: Fraction) -> int | None:
    """The piece whose source is the convex vertex the connector ``emitted[-1] -> target`` cannot pass under the cap, if proved.

    Returns the index of the outgoing piece when, all exactly: the connector
    crosses exactly one ∂Ω* vertex V, between two offset line pieces, with the
    two sites on them (neither at V); V is convex with exterior angle above the
    cap (:func:`exterior_angle_exceeds`); and the lemma's hypothesis holds on
    the emitted circles and on ∂Ω* (module docstring). Otherwise ``None``.

    Precondition (module docstring): the connector is the first arrival at V,
    the emitted sites lapping ∂Ω* once in order from circle 0.

    Raises:
        OmegaCornerFirstArrivalError: The precondition fails at a candidate vertex.
        UnknownOmegaSupportError: A ∂Ω* piece has an unknown support kind.
    """
    predecessor = emitted[-1]
    start, end = predecessor.site, target.site
    if start.same_position(end) or end.at_vertex:
        return None
    chain = _chain(loop.pieces, start, end)
    if len(chain) != 2 or chain[0].index == chain[1].index:
        return None
    incoming, outgoing = chain[0].piece.support, chain[1].piece.support
    if not isinstance(incoming, OffsetLineSupport) or not isinstance(outgoing, OffsetLineSupport):
        return None
    if not exterior_angle_exceeds(incoming, outgoing, cos_sq_half_cap):
        return None
    corner = chain[1].piece.source
    if not (corner.x.is_rational and corner.y.is_rational) or start.point() == corner:
        return None
    require_first_arrival(loop, emitted, chain[1].index)
    wedge = _Wedge(incoming, outgoing, (corner.x.rational, corner.y.rational), tool_radius)
    adjacent = {chain[0].index, chain[1].index}
    if not all(_piece_clear(wedge, loop, index) for index in range(len(loop.pieces)) if index not in adjacent):
        return None
    if not all(_circle_clear(wedge, circle) for circle in emitted):
        return None
    return chain[1].index
