"""Exact analytic Held pocket boundary: lines and arcs with closed-form distance candidates.

Every published double enters once, exactly, as an integer on a common dyadic
grid ``2^-m`` (:class:`DyadicGrid`). A line is the segment between its two
endpoints. An arc is the circle about ``centre`` of radius
``rho = |start - centre|`` (``rho²`` exact) over the closed angular span from
``dir(start - centre)`` to ``dir(end - centre)`` in the sense of
``sign(sweep)``; its true end is ``centre + rho·unit(end - centre)``, i.e.
``centre + (√K / |w|²)·w`` with ``w = end - centre`` and ``K = rho²·|w|²``.

Each query returns a finite set of exact SQUARED distances (grid units²),
every one realised by an actual pair of points, that provably contains the
extreme. Each is a :class:`~benchmarks.held_gouge_surd.Surd` of one quadratic
extension, so a comparison against an integer threshold is one exact sign.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from benchmarks.errors import BenchmarkError
from benchmarks.held_gouge_surd import Surd
from benchmarks.held_gouge_surd import sign_surd
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine

Vec = tuple[int, int]
PrimitiveKind = Literal["line", "convex-arc", "concave-arc"]
ReferencePrimitive = ReferenceLine | ReferenceArc


class GougeBoundaryError(BenchmarkError):
    """Base error of the exact gouge-audit boundary."""


class UnknownBoundaryPrimitiveError(GougeBoundaryError):
    """A boundary primitive is neither a ReferenceLine nor a ReferenceArc."""


class OffGridValueError(GougeBoundaryError):
    """A double is finer than the dyadic grid it is being injected into."""


@dataclass(frozen=True, slots=True)
class DyadicGrid:
    """The grid ``2^-exponent`` mm on which every injected double is an integer."""

    exponent: int

    @staticmethod
    def exponent_of(value: float) -> int:
        """Smallest ``m`` with ``value·2^m`` an integer."""
        return Fraction(value).denominator.bit_length() - 1

    def integer(self, value: float) -> int:
        """Exact ``value·2^exponent``.

        Raises:
            OffGridValueError: the double is not on this grid.
        """
        exact = Fraction(value)
        shift = self.exponent - (exact.denominator.bit_length() - 1)
        if shift < 0:
            raise OffGridValueError(f"{value!r} needs grid 2^-{self.exponent - shift}, finer than 2^-{self.exponent}")
        return exact.numerator << shift

    def point(self, x: float, y: float) -> Vec:
        return (self.integer(x), self.integer(y))


@dataclass(frozen=True, slots=True)
class Box:
    """Closed axis-aligned box in doubles (mm) that CONTAINS its exact geometry (outward rounded)."""

    xmin: float
    ymin: float
    xmax: float
    ymax: float


def round_down(value: Fraction) -> float:
    """Largest double that is ``<= value``."""
    nearest = float(value)
    return nearest if Fraction(nearest) <= value else math.nextafter(nearest, -math.inf)


def round_up(value: Fraction) -> float:
    """Smallest double that is ``>= value``."""
    nearest = float(value)
    return nearest if Fraction(nearest) >= value else math.nextafter(nearest, math.inf)


def primitive_exponent(primitive: ReferencePrimitive) -> int:
    """Grid exponent needed to inject a primitive exactly."""
    points = [primitive.start, primitive.end] + ([primitive.centre] if isinstance(primitive, ReferenceArc) else [])
    return max(DyadicGrid.exponent_of(v) for p in points for v in (p.x, p.y))


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1])


def dot(a: Vec, b: Vec) -> int:
    return a[0] * b[0] + a[1] * b[1]


def cross(a: Vec, b: Vec) -> int:
    return a[0] * b[1] - a[1] * b[0]


def _sign(value: int) -> int:
    return (value > 0) - (value < 0)


def point_segment_squared(x: Vec, a: Vec, b: Vec) -> Surd:
    """Exact squared distance from ``x`` to the closed segment ``ab`` (``a != b``)."""
    v = sub(b, a)
    xa = sub(x, a)
    along = dot(xa, v)
    if along <= 0:
        return Surd.of(dot(xa, xa))
    length2 = dot(v, v)
    if along >= length2:
        xb = sub(x, b)
        return Surd.of(dot(xb, xb))
    c = cross(v, xa)
    return Surd.of(c * c, length2)


def _orientation(a: Vec, b: Vec, c: Vec) -> int:
    return _sign(cross(sub(b, a), sub(c, a)))


def _wedge_from_signs(start_side: int, start_dot: int, end_in_upper: bool, to_end: int) -> bool:
    """Closed angular-span membership from exact signs, angles measured from the arc start in its sense.

    ``start_side`` / ``start_dot`` are the oriented cross and the dot of the
    start radius with the query direction, ``to_end`` the oriented cross of the
    query direction with the end radius, ``end_in_upper`` whether the end
    radius lies in the second half-turn. Within one half-turn the angular gap
    is below a half turn, so the oriented cross alone orders two directions.
    """
    query_in_upper = not (start_side > 0 or (start_side == 0 and start_dot > 0))
    if query_in_upper != end_in_upper:
        return not query_in_upper
    return to_end >= 0


@dataclass(frozen=True)
class ExactLine:
    """A boundary line segment on the grid."""

    index: int
    start: Vec
    end: Vec

    @property
    def kind(self) -> PrimitiveKind:
        return "line"

    def nearest(self, c: Vec) -> tuple[Surd, ...]:
        """Squared-distance candidates whose minimum is ``min |y - c|`` over the line."""
        return (point_segment_squared(c, self.start, self.end),)

    def farthest(self, c: Vec) -> tuple[Surd, ...]:
        """Squared-distance candidates whose maximum is ``max |y - c|`` over the line."""
        ca, cb = sub(c, self.start), sub(c, self.end)
        return (Surd.of(max(dot(ca, ca), dot(cb, cb))),)

    def segment_candidates(self, p: Vec, q: Vec) -> tuple[Surd, ...] | None:
        """Squared-distance candidates to segment ``pq``; ``None`` when the segments cross."""
        if p == q:
            return self.nearest(p)
        a, b = self.start, self.end
        if _orientation(p, q, a) * _orientation(p, q, b) < 0 and _orientation(a, b, p) * _orientation(a, b, q) < 0:
            return None
        # Touching / collinear contact makes one of these four exactly zero.
        return (point_segment_squared(p, a, b), point_segment_squared(q, a, b), point_segment_squared(a, p, q), point_segment_squared(b, p, q))


@dataclass(frozen=True)
class ExactArc:
    """A boundary circular arc on the grid; see the module docstring for its semantics."""

    index: int
    centre: Vec
    start: Vec
    start_radius: Vec
    end_radius: Vec
    orientation: int
    radius_squared: int
    end_radius_squared: int
    end_radicand: int
    full_turn: bool
    end_in_upper: bool

    @property
    def kind(self) -> PrimitiveKind:
        # The pocket is counter-clockwise: a CCW arc bulges outward (convex).
        return "convex-arc" if self.orientation > 0 else "concave-arc"

    def in_wedge(self, d: Vec) -> bool:
        """Whether the integer direction ``d`` lies in the closed angular span."""
        if self.full_turn or d == (0, 0):
            return True
        s = self.orientation
        return _wedge_from_signs(s * _sign(cross(self.start_radius, d)), _sign(dot(self.start_radius, d)), self.end_in_upper, s * _sign(cross(d, self.end_radius)))

    def in_wedge_surd(self, dx: tuple[int, int], dy: tuple[int, int], k: int) -> bool:
        """Whether the direction ``(dx[0] + dx[1]√k, dy[0] + dy[1]√k)`` lies in the span."""
        if self.full_turn or (sign_surd(dx[0], dx[1], k) == 0 and sign_surd(dy[0], dy[1], k) == 0):
            return True
        (ux, uy), (wx, wy), s = self.start_radius, self.end_radius, self.orientation
        side = s * sign_surd(ux * dy[0] - uy * dx[0], ux * dy[1] - uy * dx[1], k)
        along = sign_surd(ux * dx[0] + uy * dy[0], ux * dx[1] + uy * dy[1], k)
        to_end = s * sign_surd(dx[0] * wy - dy[0] * wx, dx[1] * wy - dy[1] * wx, k)
        return _wedge_from_signs(side, along, self.end_in_upper, to_end)

    def _end_squared(self, c: Vec) -> Surd:
        """Exact ``|c - true_end|²`` = ``((|c-o|² + rho²)|w|² - 2 (c-o)·w √K) / |w|²``."""
        if self.full_turn:
            ca = sub(c, self.start)
            return Surd.of(dot(ca, ca))
        co = sub(c, self.centre)
        w2 = self.end_radius_squared
        return Surd((dot(co, co) + self.radius_squared) * w2, -2 * dot(co, self.end_radius), self.end_radicand, w2)

    def _extreme(self, c: Vec, sign: int) -> tuple[Surd, ...]:
        """Candidates for the min (``sign = -1``) or max (``sign = 1``) of ``|y - c|`` over the arc.

        The extreme point of the full circle is ``o - sign·rho·dir(c - o)``
        (nearest for -1, farthest for 1); when it lies in the span it IS the extreme
        over the arc, otherwise ``|y - c|`` is monotone along the arc away from
        it and the extreme is at an endpoint.
        """
        co = sub(c, self.centre)
        d2 = dot(co, co)
        rho2 = self.radius_squared
        if d2 == 0:
            return (Surd.of(rho2),)
        if self.in_wedge(co if sign < 0 else (-co[0], -co[1])):
            # (d ∓ rho)² = d² + rho² ∓ 2√(d²rho²)
            return (Surd(d2 + rho2, 2 * sign, d2 * rho2, 1),)
        ca = sub(c, self.start)
        return (Surd.of(dot(ca, ca)), self._end_squared(c))

    def nearest(self, c: Vec) -> tuple[Surd, ...]:
        """Squared-distance candidates whose minimum is ``min |y - c|`` over the arc."""
        return self._extreme(c, -1)

    def farthest(self, c: Vec) -> tuple[Surd, ...]:
        """Squared-distance candidates whose maximum is ``max |y - c|`` over the arc."""
        return self._extreme(c, 1)

    def _end_to_segment_squared(self, p: Vec, q: Vec) -> Surd:
        """Exact squared distance from the arc's true end to segment ``pq`` (``p != q``)."""
        if self.full_turn:
            return point_segment_squared(self.start, p, q)
        w, w2, k, rho2 = self.end_radius, self.end_radius_squared, self.end_radicand, self.radius_squared
        v = sub(q, p)
        length2 = dot(v, v)
        op = sub(self.centre, p)
        # (true_end - p)·v = op·v + (w·v)·√K / |w|²
        if sign_surd(dot(op, v) * w2, dot(w, v), k) <= 0:
            return Surd((dot(op, op) + rho2) * w2, 2 * dot(op, w), k, w2)
        if sign_surd((dot(op, v) - length2) * w2, dot(w, v), k) >= 0:
            oq = sub(self.centre, q)
            return Surd((dot(oq, oq) + rho2) * w2, 2 * dot(oq, w), k, w2)
        c0, c1 = cross(v, op), cross(v, w)
        return Surd(c0 * c0 * w2 + c1 * c1 * rho2, 2 * c0 * c1, k, length2 * w2)

    def _meets_segment(self, p: Vec, q: Vec) -> bool:
        v = sub(q, p)
        po = sub(p, self.centre)
        a = dot(v, v)
        b = 2 * dot(v, po)
        disc = b * b - 4 * a * (dot(po, po) - self.radius_squared)
        if disc < 0:
            return False
        for s in (1, -1):
            # 2a·t = -b + s·√disc must lie in [0, 2a]
            if sign_surd(-b, s, disc) < 0 or sign_surd(-b - 2 * a, s, disc) > 0:
                continue
            # 2a·(point - o) = 2a·po + v·(-b + s√disc), a positive multiple of the direction.
            dx = (2 * a * po[0] - b * v[0], s * v[0])
            dy = (2 * a * po[1] - b * v[1], s * v[1])
            if self.in_wedge_surd(dx, dy, disc):
                return True
        return False

    def segment_candidates(self, p: Vec, q: Vec) -> tuple[Surd, ...] | None:
        """Squared-distance candidates to segment ``pq``; ``None`` when they meet.

        A closest pair is: an intersection; a segment endpoint against the arc;
        an arc endpoint against the segment; or two interior points, where the
        arc point is ``o ± rho·n`` for the segment normal ``n``.
        """
        if p == q:
            return self.nearest(p)
        if self._meets_segment(p, q):
            return None
        out = [*self.nearest(p), *self.nearest(q), point_segment_squared(self.start, p, q), self._end_to_segment_squared(p, q)]
        v = sub(q, p)
        length2 = dot(v, v)
        op = sub(self.centre, p)
        if 0 < dot(op, v) < length2:
            normal = (-v[1], v[0])
            height = dot(op, normal)
            product = self.radius_squared * length2
            for s in (1, -1):
                if self.in_wedge((s * normal[0], s * normal[1])):
                    # (height/√A + s·rho)² = (height² + rho²A + 2·s·height·√(rho²A)) / A
                    out.append(Surd(height * height + product, 2 * s * height, product, length2))
        return tuple(out)


ExactPrimitive = ExactLine | ExactArc


def exact_primitive(index: int, primitive: ReferencePrimitive, grid: DyadicGrid) -> ExactPrimitive:
    """Inject one published primitive exactly onto ``grid``.

    Raises:
        UnknownBoundaryPrimitiveError: not a ReferenceLine or ReferenceArc.
    """
    if isinstance(primitive, ReferenceLine):
        return ExactLine(index, grid.point(primitive.start.x, primitive.start.y), grid.point(primitive.end.x, primitive.end.y))
    if isinstance(primitive, ReferenceArc):
        o = grid.point(primitive.centre.x, primitive.centre.y)
        a = grid.point(primitive.start.x, primitive.start.y)
        b = grid.point(primitive.end.x, primitive.end.y)
        u, w = sub(a, o), sub(b, o)
        s = 1 if float(primitive.sweep) > 0 else -1
        side = s * cross(u, w)
        end_in_upper = not (side > 0 or (side == 0 and dot(u, w) > 0))
        rho2, w2 = dot(u, u), dot(w, w)
        return ExactArc(index, o, a, u, w, s, rho2, w2, rho2 * w2, a == b, end_in_upper)
    raise UnknownBoundaryPrimitiveError(f"primitive {index} is {type(primitive).__name__}")


def _root_bounds(square: int) -> tuple[int, int]:
    """Integers ``lo <= √square <= hi`` with ``hi - lo <= 1``."""
    lo = math.isqrt(square)
    return lo, lo + (lo * lo != square)


def exact_box(primitive: ExactPrimitive, grid: DyadicGrid) -> Box:
    """Outward-rounded box (mm doubles) that contains the exact primitive.

    An arc's box spans its start, an exact integer bracket of its true end, and
    the axis-extreme points ``o ± rho·e`` of every axis direction ``e`` in its span.
    """
    scale = Fraction(1, 1 << grid.exponent)
    if isinstance(primitive, ExactLine):
        xs: list[Fraction] = [Fraction(primitive.start[0]), Fraction(primitive.end[0])]
        ys: list[Fraction] = [Fraction(primitive.start[1]), Fraction(primitive.end[1])]
    else:
        (ox, oy), (wx, wy), w2 = primitive.centre, primitive.end_radius, primitive.end_radius_squared
        xs, ys = [Fraction(primitive.start[0])], [Fraction(primitive.start[1])]
        if not primitive.full_turn:
            for root in _root_bounds(primitive.end_radicand):
                xs.append(ox + Fraction(wx * root, w2))
                ys.append(oy + Fraction(wy * root, w2))
        rho_hi = _root_bounds(primitive.radius_squared)[1]
        for direction, values, base, step in (((1, 0), xs, ox, rho_hi), ((-1, 0), xs, ox, -rho_hi), ((0, 1), ys, oy, rho_hi), ((0, -1), ys, oy, -rho_hi)):
            if primitive.in_wedge(direction):
                values.append(Fraction(base + step))
    return Box(round_down(min(xs) * scale), round_down(min(ys) * scale), round_up(max(xs) * scale), round_up(max(ys) * scale))
