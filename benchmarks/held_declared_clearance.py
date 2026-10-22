"""Exact clearance of Held path elements against the declared (Pythagorean) boundary.

A path element is a locus of the tool centre: a binary64 machining circle, a
segment or an arc of an exact ∂Ω* connector (one-root endpoints, rational
support). The swept tool disk of radius ``r`` leaves the pocket exactly when
``dist(element, ∂P) < r``; ``== r`` is a touch (allowed) and ``> r`` is clear.

Every decision here is exact over the declared boundary's rationals and the
elements' one-root numbers (:mod:`benchmarks.held_one_root`). A distance is
one closed-form candidate (:class:`_Distance`) and every comparison with a
rational threshold reduces to the sign of a one-root number. Reported depths
are binary64 views and never decide.

Distances between two curves are the minimum over the closed-form critical
pairs: each endpoint against the other curve, and the interior-interior pairs
with a common normal (the arc point on the normal of a line through its
centre; the two arc points on the line of centres). Every candidate is the
distance of a pair of real points of the two curves, and the minimum is among
them, unless the curves cross. Precondition (stated, as in
:mod:`benchmarks.held_gouge_audit`): every element lies in the closed pocket, so
a curve meets ∂P only tangentially or at an endpoint, which the candidates
cover. Circles are checked against the precondition exactly: the caller locates
one point of the circle in Ω* (:meth:`DeclaredClearance.circle`).

Pruning is exact: a primitive is skipped only when the outward-rounded boxes
are more than ``r`` apart, decided on the boxes' rationals.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Protocol
from typing import Self
from typing import TypeAlias

import numpy as np
import numpy.typing as npt

from benchmarks.errors import BenchmarkError
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import Sense
from benchmarks.held_dyadic_interval import FRACTION_BITS
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_dyadic_interval import maximum
from benchmarks.held_dyadic_interval import minimum
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_boundary import PrimitiveKind
from benchmarks.held_gouge_boundary import round_down
from benchmarks.held_gouge_boundary import round_up
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from compas_cgal.units import WorldXY

_Vector: TypeAlias = tuple[Fraction, Fraction]
# Reporting conversion of an interval end to millimetres (never a decision).
_INTERVAL_SCALE = 2.0**-FRACTION_BITS
# Bits of the integer square-root bracket that encloses a one-root coordinate
# for box pruning; pruning precision only, never a decision.
ENCLOSURE_BITS = 64


class DeclaredClearanceError(BenchmarkError):
    """Base error of the declared-boundary clearance oracle."""


class NonPositiveClearanceRadiusError(DeclaredClearanceError):
    """The tool radius or a circle radius is not positive."""


class DegenerateClearanceElementError(DeclaredClearanceError):
    """A segment element has no direction or an arc element no radius."""


# ---------------------------------------------------------------------------
# Distances as closed forms
# ---------------------------------------------------------------------------


def _sqrt_sign(square: OneRootNumber, threshold: Fraction) -> int:
    """Exact ``sign(√square - threshold)`` for ``square >= 0``."""
    if threshold < 0:
        return 1
    return (square - threshold * threshold).sign()


@dataclass(frozen=True)
class _Root:
    """A distance ``√square``."""

    square: OneRootNumber

    def sign_minus(self, threshold: Fraction) -> int:
        return _sqrt_sign(self.square, threshold)

    def approx(self) -> float:
        return math.sqrt(max(float(self.square), 0.0))


@dataclass(frozen=True)
class _Absolute:
    """A distance ``|value|``."""

    value: OneRootNumber

    def sign_minus(self, threshold: Fraction) -> int:
        return max((self.value - threshold).sign(), (-self.value - threshold).sign())

    def approx(self) -> float:
        return abs(float(self.value))


@dataclass(frozen=True)
class _ShiftedRoot:
    """A distance ``|√square + shift|`` (a point and a circle, or two circles on their centre line)."""

    square: OneRootNumber
    shift: Fraction

    def sign_minus(self, threshold: Fraction) -> int:
        # |√D + k| - m = max(√D + k - m, -√D - k - m).
        return max(_sqrt_sign(self.square, threshold - self.shift), -_sqrt_sign(self.square, -self.shift - threshold))

    def approx(self) -> float:
        return abs(math.sqrt(max(float(self.square), 0.0)) + float(self.shift))


_Distance: TypeAlias = "_Root | _Absolute | _ShiftedRoot"


# ---------------------------------------------------------------------------
# Elements and primitives
# ---------------------------------------------------------------------------


def _dot(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _one_root_dot(direction: _Vector, point: OneRootPoint2) -> OneRootNumber:
    return direction[0] * point.x + direction[1] * point.y


def _squared_distance(point: OneRootPoint2, origin: ExactPoint2[WorldXY]) -> OneRootNumber:
    dx, dy = point.minus(origin)
    return dx * dx + dy * dy


@dataclass(frozen=True)
class SegmentElement:
    """A straight piece of tool-centre motion on the rational line ``normal·x = level``.

    ``direction`` is any rational vector along the line (``normal·direction == 0``);
    the endpoints may lie in different quadratic extensions.
    """

    start: OneRootPoint2
    end: OneRootPoint2
    direction: _Vector
    normal: _Vector
    level: Fraction

    @classmethod
    def between(cls, start: ExactPoint2[WorldXY], end: ExactPoint2[WorldXY]) -> Self:
        """The segment between two rational points (a binary64 polyline piece).

        Raises:
            DegenerateClearanceElementError: If the points coincide.
        """
        direction = end.minus(start)
        if direction[0] == 0 and direction[1] == 0:
            raise DegenerateClearanceElementError("A segment element needs two distinct points.")
        normal = (-direction[1], direction[0])
        return cls(OneRootPoint2.from_exact(start), OneRootPoint2.from_exact(end), direction, normal, _dot(normal, (start.x, start.y)))

    def contains_parameter_of(self, point: ExactPoint2[WorldXY]) -> bool:
        """Whether the foot of ``point`` on the line lies in the closed segment (exact, across extensions)."""
        position = OneRootNumber.of(_dot(self.direction, (point.x, point.y)))
        return _one_root_dot(self.direction, self.start) <= position <= _one_root_dot(self.direction, self.end)

    def distance_to_point(self, point: ExactPoint2[WorldXY]) -> _Distance:
        """Exact distance from a rational point."""
        if self.contains_parameter_of(point):
            offset = _dot(self.normal, (point.x, point.y)) - self.level
            return _Root(OneRootNumber.of(offset * offset / _dot(self.normal, self.normal)))
        return _Root(min(_squared_distance(self.start, point), _squared_distance(self.end, point)))


@dataclass(frozen=True)
class ArcElement:
    """A piece of tool-centre motion on the rational circle (``centre``, ``radius``) from ``start`` to ``end`` in ``sense``.

    ``start == end`` is the full circle.
    """

    centre: ExactPoint2[WorldXY]
    radius: Fraction
    sense: Sense
    start: OneRootPoint2
    end: OneRootPoint2

    @property
    def full(self) -> bool:
        return self.start == self.end

    def spans(self, probe: OneRootPoint2) -> bool:
        """Whether the ray from the centre through ``probe`` meets this arc (exact)."""
        return self.full or compare_on_arc(self.centre, self.sense, self.start, probe, self.end) <= 0

    def distance_to_point(self, point: ExactPoint2[WorldXY]) -> _Distance:
        """Exact distance from a rational point."""
        offset = point.minus(self.centre)
        square = OneRootNumber.of(_dot(offset, offset))
        if square.sign() == 0:
            return _Absolute(OneRootNumber.of(self.radius))
        if self.spans(OneRootPoint2.from_exact(point)):
            return _ShiftedRoot(square, -self.radius)
        return _Root(min(_squared_distance(self.start, point), _squared_distance(self.end, point)))


@dataclass(frozen=True)
class CircleElement:
    """A binary64 machining circle (rational centre and radius)."""

    centre: ExactPoint2[WorldXY]
    radius: Fraction


Element: TypeAlias = "SegmentElement | ArcElement"


def _cross(left: _Vector, right: _Vector) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def _filtered_sign(enclosure: Interval, exact: Fraction | Callable[[], Fraction]) -> int:
    """Sign decided by the enclosure when it excludes zero, else by the exact rational."""
    decided = enclosure.sign()
    if decided is not None:
        return decided
    value = exact() if callable(exact) else exact
    return (value > 0) - (value < 0)


def _box(vector: _Vector) -> tuple[Interval, Interval]:
    return Interval.of(vector[0]), Interval.of(vector[1])


def _rational_half(start: _Vector, start_box: tuple[Interval, Interval], sense: Sense, vector: _Vector, vector_box: tuple[Interval, Interval]) -> int:
    side = sense * _filtered_sign(start_box[0] * vector_box[1] - start_box[1] * vector_box[0], lambda: _cross(start, vector))
    if side != 0:
        return 0 if side > 0 else 1
    return 0 if _filtered_sign(start_box[0] * vector_box[0] + start_box[1] * vector_box[1], lambda: _dot(start, vector)) > 0 else 1


def _rational_spans(start: _Vector, end: _Vector, sense: Sense, vector: _Vector, boxes: tuple[tuple[Interval, Interval], tuple[Interval, Interval], int]) -> bool:
    """``compare_on_arc(c, σ, S, c + vector, E) <= 0`` on rational vectors from the centre (filtered, no roots)."""
    start_box, end_box, end_half = boxes
    vector_box = _box(vector)
    probe_half = _rational_half(start, start_box, sense, vector, vector_box)
    if probe_half != end_half:
        return probe_half < end_half
    return sense * _filtered_sign(vector_box[0] * end_box[1] - vector_box[1] * end_box[0], lambda: _cross(vector, end)) >= 0


def _rational_squared(point: ExactPoint2[WorldXY], origin: ExactPoint2[WorldXY]) -> OneRootNumber:
    dx, dy = point.x - origin.x, point.y - origin.y
    return OneRootNumber.of(dx * dx + dy * dy)


@dataclass(frozen=True)
class _Enclosed:
    """Interval enclosures of a primitive's exact data (filter stage)."""

    start: tuple[Interval, Interval]
    end: tuple[Interval, Interval]
    centre: tuple[Interval, Interval]
    radius: Interval
    unit: tuple[Interval, Interval]
    normal: tuple[Interval, Interval]
    length: Interval


def _pair(x: Fraction, y: Fraction) -> tuple[Interval, Interval]:
    return Interval.of(x), Interval.of(y)


def _iv_distance(point: tuple[Interval, Interval], other: tuple[Interval, Interval]) -> Interval:
    return ((point[0] - other[0]).square() + (point[1] - other[1]).square()).sqrt()


@dataclass(frozen=True)
class _Line:
    index: int
    start: ExactPoint2[WorldXY]
    end: ExactPoint2[WorldXY]
    unit: _Vector
    normal: _Vector
    length: Fraction
    enclosed: _Enclosed

    def foot_in(self, point: OneRootPoint2) -> bool:
        position = _one_root_dot(self.unit, point) - _dot(self.unit, (self.start.x, self.start.y))
        return position.sign() >= 0 and (position - self.length).sign() <= 0

    def distance_to(self, point: OneRootPoint2) -> _Distance:
        if self.foot_in(point):
            return _Absolute(_one_root_dot(self.normal, point) - _dot(self.normal, (self.start.x, self.start.y)))
        return _Root(min(_squared_distance(point, self.start), _squared_distance(point, self.end)))

    def nearest(self, centre: ExactPoint2[WorldXY]) -> tuple[_Distance, ...]:
        offset = centre.minus(self.start)
        position = _dot(self.unit, offset)
        if 0 <= position <= self.length:
            return (_Absolute(OneRootNumber.of(_dot(self.normal, offset))),)
        return (_Root(_rational_squared(centre, self.start if position < 0 else self.end)),)

    def farthest(self, centre: ExactPoint2[WorldXY]) -> tuple[_Distance, ...]:
        return (_Root(_rational_squared(centre, self.start)), _Root(_rational_squared(centre, self.end)))

    def circle_bounds(self, centre: tuple[Interval, Interval]) -> tuple[Interval, Interval, Interval]:
        """Enclosures of a lower bound of the nearest, and a lower and an upper bound of the farthest distance from ``centre``."""
        box = self.enclosed
        dx, dy = centre[0] - box.start[0], centre[1] - box.start[1]
        position = box.unit[0] * dx + box.unit[1] * dy
        to_start, to_end = _iv_distance(centre, box.start), _iv_distance(centre, box.end)
        if position.hi < 0:
            near = to_start
        elif (position - box.length).lo > 0:
            near = to_end
        else:
            # The distance to the line never exceeds the distance to the segment.
            near = (box.normal[0] * dx + box.normal[1] * dy).absolute()
        far = maximum([to_start, to_end])
        return near, far, far


@dataclass(frozen=True)
class _Arc:
    index: int
    start: ExactPoint2[WorldXY]
    end: ExactPoint2[WorldXY]
    centre: ExactPoint2[WorldXY]
    radius: Fraction
    sense: Sense
    enclosed: _Enclosed

    def spans(self, probe: OneRootPoint2) -> bool:
        if probe.radicand == 0:
            return self.spans_rational(probe.x.rational - self.centre.x, probe.y.rational - self.centre.y)
        return compare_on_arc(self.centre, self.sense, OneRootPoint2.from_exact(self.start), probe, OneRootPoint2.from_exact(self.end)) <= 0

    def spans_rational(self, x: Fraction, y: Fraction) -> bool:
        start, end = self.start.minus(self.centre), self.end.minus(self.centre)
        start_box, end_box = _box(start), _box(end)
        end_half = _rational_half(start, start_box, self.sense, end, end_box)
        return _rational_spans(start, end, self.sense, (x, y), (start_box, end_box, end_half))

    def distance_to(self, point: OneRootPoint2) -> _Distance:
        square = _squared_distance(point, self.centre)
        if square.sign() == 0:
            return _Absolute(OneRootNumber.of(self.radius))
        if self.spans(point):
            return _ShiftedRoot(square, -self.radius)
        return _Root(min(_squared_distance(point, self.start), _squared_distance(point, self.end)))

    def nearest(self, centre: ExactPoint2[WorldXY]) -> tuple[_Distance, ...]:
        dx, dy = centre.minus(self.centre)
        square = OneRootNumber.of(dx * dx + dy * dy)
        if square.sign() == 0:
            return (_Absolute(OneRootNumber.of(self.radius)),)
        if self.spans_rational(dx, dy):
            return (_ShiftedRoot(square, -self.radius),)
        return (_Root(min(_rational_squared(centre, self.start), _rational_squared(centre, self.end))),)

    def farthest(self, centre: ExactPoint2[WorldXY]) -> tuple[_Distance, ...]:
        dx, dy = centre.minus(self.centre)
        square = OneRootNumber.of(dx * dx + dy * dy)
        if square.sign() == 0:
            return (_Absolute(OneRootNumber.of(self.radius)),)
        ends = (_Root(_rational_squared(centre, self.start)), _Root(_rational_squared(centre, self.end)))
        return (*ends, _ShiftedRoot(square, self.radius)) if self.spans_rational(-dx, -dy) else ends

    def circle_bounds(self, centre: tuple[Interval, Interval]) -> tuple[Interval, Interval, Interval]:
        """Enclosures of a lower bound of the nearest, and a lower and an upper bound of the farthest distance from ``centre``."""
        box = self.enclosed
        radial = _iv_distance(centre, box.centre)
        # The full circle is never farther than the arc; its far side never nearer
        # than the arc's farthest point; the arc's ends never farther than it.
        near = (radial - box.radius).absolute()
        ends = maximum([_iv_distance(centre, box.start), _iv_distance(centre, box.end)])
        return near, ends, maximum([ends, radial + box.radius])


_Primitive: TypeAlias = "_Line | _Arc"


def _primitive(index: int, primitive: DeclaredLine | DeclaredArc) -> _Primitive:
    zero = Fraction(0)
    start, end = _pair(primitive.start.x, primitive.start.y), _pair(primitive.end.x, primitive.end.y)
    if isinstance(primitive, DeclaredLine):
        unit = primitive.direction.vector
        normal = (-unit[1], unit[0])
        enclosed = _Enclosed(start, end, _pair(zero, zero), Interval.of(0), _pair(*unit), _pair(*normal), Interval.of(primitive.length))
        return _Line(index, primitive.start, primitive.end, unit, normal, primitive.length, enclosed)
    centre = primitive.centre
    enclosed = _Enclosed(start, end, _pair(centre.x, centre.y), Interval.of(primitive.radius), _pair(zero, zero), _pair(zero, zero), Interval.of(0))
    return _Arc(index, primitive.start, primitive.end, centre, primitive.radius, primitive.sense, enclosed)


# ---------------------------------------------------------------------------
# Curve-to-primitive candidates
# ---------------------------------------------------------------------------


def _segment_candidates(element: SegmentElement, primitive: _Primitive) -> list[_Distance]:
    candidates = [primitive.distance_to(element.start), primitive.distance_to(element.end)]
    candidates += [element.distance_to_point(primitive.start), element.distance_to_point(primitive.end)]
    if isinstance(primitive, _Arc):
        # Arc points on the line's normal through the arc centre: C ∓ ρ·n̂.
        norm_squared = _dot(element.normal, element.normal)
        height = _dot(element.normal, (primitive.centre.x, primitive.centre.y)) - element.level
        if element.contains_parameter_of(primitive.centre):
            for side in (1, -1):
                probe = OneRootPoint2.rational(primitive.centre.x - side * element.normal[0], primitive.centre.y - side * element.normal[1])
                if primitive.spans(probe):
                    # Signed distance of C - side·ρ·n̂ to the line: h/|n| - side·ρ.
                    candidates.append(_Absolute(OneRootNumber.build(-side * primitive.radius, height / norm_squared, norm_squared)))
    return candidates


def _arc_candidates(element: ArcElement, primitive: _Primitive) -> list[_Distance]:
    candidates: list[_Distance] = [element.distance_to_point(primitive.start), element.distance_to_point(primitive.end)]
    if not element.full:
        candidates += [primitive.distance_to(element.start), primitive.distance_to(element.end)]
    if isinstance(primitive, _Line):
        position = _dot(primitive.unit, element.centre.minus(primitive.start))
        if 0 <= position <= primitive.length:
            height = _dot(primitive.normal, element.centre.minus(primitive.start))
            for side in (1, -1):
                probe = OneRootPoint2.rational(element.centre.x - side * primitive.normal[0], element.centre.y - side * primitive.normal[1])
                if element.spans(probe):
                    candidates.append(_Absolute(OneRootNumber.of(height - side * element.radius)))
        return candidates
    between = primitive.centre.minus(element.centre)
    square = _dot(between, between)
    if square == 0:
        return candidates
    for own in (1, -1):
        probe = OneRootPoint2.rational(element.centre.x + own * between[0], element.centre.y + own * between[1])
        if not element.spans(probe):
            continue
        for other in (1, -1):
            # x_E = C' + own·ρ'·v̂, x_X = C + other·ρ·v̂ (v = C - C'): |x_X - x_E| = |√V + other·ρ - own·ρ'|.
            target = OneRootPoint2.rational(primitive.centre.x + other * between[0], primitive.centre.y + other * between[1])
            if primitive.spans(target):
                candidates.append(_ShiftedRoot(OneRootNumber.of(square), other * primitive.radius - own * element.radius))
    return candidates


# ---------------------------------------------------------------------------
# Crossings (distance zero): decided, not assumed away
# ---------------------------------------------------------------------------


def _line_circle_points(point: _Vector, direction: _Vector, centre: ExactPoint2[WorldXY], radius: Fraction) -> list[tuple[OneRootNumber, OneRootPoint2]]:
    """Points ``point + λ·direction`` on the circle, with their ``λ`` (one root each)."""
    offset = (point[0] - centre.x, point[1] - centre.y)
    a = _dot(direction, direction)
    b = 2 * _dot(direction, offset)
    c = _dot(offset, offset) - radius * radius
    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return []
    found = []
    for side in (1, -1) if discriminant > 0 else (1,):
        along = OneRootNumber.build(-b / (2 * a), Fraction(side) / (2 * a), discriminant)
        found.append((along, OneRootPoint2.build(point[0] + along * direction[0], point[1] + along * direction[1])))
    return found


def _circle_circle_points(first: ExactPoint2[WorldXY], first_radius: Fraction, second: ExactPoint2[WorldXY], second_radius: Fraction) -> list[OneRootPoint2]:
    """Intersection points of two non-concentric circles (one root)."""
    dx, dy = second.x - first.x, second.y - first.y
    distance = dx * dx + dy * dy
    if distance == 0:
        return []
    along = (first_radius * first_radius - second_radius * second_radius + distance) / (2 * distance)
    height = first_radius * first_radius / distance - along * along
    if height < 0:
        return []
    base = (first.x + along * dx, first.y + along * dy)
    return [OneRootPoint2.build(OneRootNumber.build(base[0], -side * dy, height), OneRootNumber.build(base[1], side * dx, height)) for side in ((1, -1) if height > 0 else (1,))]


def _segment_range(element: SegmentElement) -> tuple[OneRootNumber, OneRootNumber]:
    return _one_root_dot(element.direction, element.start), _one_root_dot(element.direction, element.end)


def _line_anchor(element: SegmentElement) -> _Vector:
    scale = element.level / _dot(element.normal, element.normal)
    return element.normal[0] * scale, element.normal[1] * scale


def _segment_crosses(element: SegmentElement, primitive: _Primitive) -> bool:
    low, high = _segment_range(element)
    if isinstance(primitive, _Line):
        start_side = _dot(element.normal, (primitive.start.x, primitive.start.y)) - element.level
        end_side = _dot(element.normal, (primitive.end.x, primitive.end.y)) - element.level
        if start_side * end_side > 0:
            return False
        if start_side == 0 and end_side == 0:
            ends = sorted((_dot(element.direction, (primitive.start.x, primitive.start.y)), _dot(element.direction, (primitive.end.x, primitive.end.y))))
            return not (low > ends[1] or high < ends[0])
        fraction = start_side / (start_side - end_side)
        crossing = (primitive.start.x + fraction * (primitive.end.x - primitive.start.x), primitive.start.y + fraction * (primitive.end.y - primitive.start.y))
        return low <= OneRootNumber.of(_dot(element.direction, crossing)) <= high
    anchor = _line_anchor(element)
    base = _dot(element.direction, anchor)
    scale = _dot(element.direction, element.direction)
    return any(
        low <= base + along * scale <= high and primitive.spans(point) for along, point in _line_circle_points(anchor, element.direction, primitive.centre, primitive.radius)
    )


def _arc_crosses(element: ArcElement, primitive: _Primitive) -> bool:
    if isinstance(primitive, _Line):
        anchor = (primitive.start.x, primitive.start.y)
        points = _line_circle_points(anchor, primitive.unit, element.centre, element.radius)
        return any(along.sign() >= 0 and along <= primitive.length and element.spans(point) for along, point in points)
    return any(element.spans(point) and primitive.spans(point) for point in _circle_circle_points(element.centre, element.radius, primitive.centre, primitive.radius))


# ---------------------------------------------------------------------------
# Segment filter: an enclosure of a lower bound of dist(segment, primitive)
# ---------------------------------------------------------------------------

_IntervalPair: TypeAlias = "tuple[Interval, Interval]"


@dataclass(frozen=True)
class _SegmentBox:
    """Enclosures of a segment element: ends, direction, unit normal and parameter range."""

    start: _IntervalPair
    end: _IntervalPair
    direction: _IntervalPair
    unit_normal: _IntervalPair
    length_squared: Interval

    @classmethod
    def of(cls, element: SegmentElement) -> _SegmentBox | None:
        """The enclosures, or ``None`` when the normal is below the enclosure resolution (the filter then declines)."""
        start = (Interval.of_one_root(element.start.x), Interval.of_one_root(element.start.y))
        end = (Interval.of_one_root(element.end.x), Interval.of_one_root(element.end.y))
        direction = (end[0] - start[0], end[1] - start[1])
        norm = Interval.sqrt_of(_dot(element.normal, element.normal))
        if norm.lo <= 0:
            return None
        unit = (Interval.of(element.normal[0]) / norm, Interval.of(element.normal[1]) / norm)
        return cls(start, end, direction, unit, direction[0].square() + direction[1].square())

    def point_lower(self, point: _IntervalPair) -> Interval:
        """Enclosure of a lower bound of the distance from ``point`` to the segment."""
        dx, dy = point[0] - self.start[0], point[1] - self.start[1]
        along = self.direction[0] * dx + self.direction[1] * dy
        if along.hi < 0:
            return _iv_distance(point, self.start)
        if (along - self.length_squared).lo > 0:
            return _iv_distance(point, self.end)
        # The distance to the line never exceeds the distance to the segment.
        return (self.unit_normal[0] * dx + self.unit_normal[1] * dy).absolute()


def _side(point: _IntervalPair, origin: _IntervalPair, normal: _IntervalPair) -> int | None:
    return (normal[0] * (point[0] - origin[0]) + normal[1] * (point[1] - origin[1])).sign()


def _segment_lower(box: _SegmentBox, primitive: _Primitive) -> Interval | None:
    """Enclosure of a lower bound of ``dist(segment, primitive)``; ``None`` when a crossing is not excluded."""
    enclosed = primitive.enclosed
    if isinstance(primitive, _Arc):
        # Every arc point is on its circle: dist >= max(lo_C - ρ, ρ - hi_C) for the segment's distances to C.
        near = box.point_lower(enclosed.centre)
        far = maximum([_iv_distance(enclosed.centre, box.start), _iv_distance(enclosed.centre, box.end)])
        return maximum([near - enclosed.radius, enclosed.radius - far])
    first, second = _side(box.start, enclosed.start, enclosed.normal), _side(box.end, enclosed.start, enclosed.normal)
    if first is None or first != second:
        return None
    to_line = minimum([_line_point_lower(enclosed, box.start), _line_point_lower(enclosed, box.end)])
    return minimum([to_line, box.point_lower(enclosed.start), box.point_lower(enclosed.end)])


def _line_point_lower(enclosed: _Enclosed, point: _IntervalPair) -> Interval:
    dx, dy = point[0] - enclosed.start[0], point[1] - enclosed.start[1]
    along = enclosed.unit[0] * dx + enclosed.unit[1] * dy
    if along.hi < 0:
        return _iv_distance(point, enclosed.start)
    if (along - enclosed.length).lo > 0:
        return _iv_distance(point, enclosed.end)
    return (enclosed.normal[0] * dx + enclosed.normal[1] * dy).absolute()


# ---------------------------------------------------------------------------
# Boxes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Box:
    """Axis-aligned box with binary64 ends rounded outward from exact bounds."""

    xmin: float
    ymin: float
    xmax: float
    ymax: float


def _enclose(value: OneRootNumber) -> tuple[Fraction, Fraction]:
    """Rational bracket of a one-root number from an integer square root."""
    if value.is_rational:
        return value.rational, value.rational
    scale = 1 << ENCLOSURE_BITS
    numerator, denominator = value.radicand.numerator, value.radicand.denominator
    root = math.isqrt(numerator * denominator * scale * scale)
    low, high = Fraction(root, denominator * scale), Fraction(root + 1, denominator * scale)
    first, second = value.rational + value.coefficient * low, value.rational + value.coefficient * high
    return min(first, second), max(first, second)


def _outward(low: Fraction, high: Fraction) -> tuple[float, float]:
    return round_down(low), round_up(high)


def _point_box(points: Sequence[OneRootPoint2]) -> Box:
    xs = [bound for point in points for bound in _enclose(point.x)]
    ys = [bound for point in points for bound in _enclose(point.y)]
    xmin, xmax = _outward(min(xs), max(xs))
    ymin, ymax = _outward(min(ys), max(ys))
    return Box(xmin, ymin, xmax, ymax)


def _disk_box(centre: ExactPoint2[WorldXY], radius: Fraction) -> Box:
    xmin, xmax = _outward(centre.x - radius, centre.x + radius)
    ymin, ymax = _outward(centre.y - radius, centre.y + radius)
    return Box(xmin, ymin, xmax, ymax)


def _primitive_box(primitive: _Primitive) -> Box:
    if isinstance(primitive, _Arc):
        return _disk_box(primitive.centre, primitive.radius)
    return _point_box((OneRootPoint2.from_exact(primitive.start), OneRootPoint2.from_exact(primitive.end)))


def element_box(element: Element | CircleElement) -> Box:
    """Outward box of an element (an arc's full circle)."""
    if isinstance(element, (ArcElement, CircleElement)):
        return _disk_box(element.centre, element.radius)
    return _point_box((element.start, element.end))


# ---------------------------------------------------------------------------
# The oracle
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ClearanceVerdict:
    """Exact clearance of one element and the primitive that decides it.

    ``primitive`` is the declared primitive of the worst (then nearest) pair,
    ``None`` when every primitive was pruned (CLEAR). ``depth_mm`` reports
    ``r - dist`` of that pair (binary64; ``0.0`` unless GOUGE).
    """

    clearance: Clearance
    primitive: int | None
    depth_mm: float


_CLEAR = ClearanceVerdict(Clearance.CLEAR, None, 0.0)


class CircleClearanceOracle(Protocol):
    """Exact clearance of binary64 machining circles, asked one batch at a time.

    Implemented by :class:`DeclaredClearance` (the Python reference) and by
    :class:`benchmarks.held_declared_clearance_native.NativeDeclaredClearance`
    (one native call per batch); both decide the same ``clearance`` per circle.
    """

    def circles(self, elements: Sequence[CircleElement]) -> list[ClearanceVerdict]:
        """Verdicts of the circles, in order."""
        ...


@dataclass(frozen=True)
class DeclaredClearance:
    """Exact clearance oracle over one declared boundary and one tool radius."""

    boundary: DeclaredBoundary
    tool_radius: Fraction
    primitives: tuple[_Primitive, ...]
    boxes: npt.NDArray[np.float64]

    @classmethod
    def build(cls, boundary: DeclaredBoundary, tool_radius: Fraction) -> Self:
        """Index the declared primitives.

        Raises:
            NonPositiveClearanceRadiusError: If ``tool_radius <= 0``.
        """
        if tool_radius <= 0 or Fraction(float(tool_radius)) != tool_radius:
            raise NonPositiveClearanceRadiusError(f"Tool radius {tool_radius} is not a positive binary64.")
        primitives = tuple(_primitive(index, primitive) for index, primitive in enumerate(boundary.primitives))
        boxes = [_primitive_box(primitive) for primitive in primitives]
        array = np.array([(box.xmin, box.ymin, box.xmax, box.ymax) for box in boxes], dtype=np.float64)
        return cls(boundary, Fraction(tool_radius), primitives, array)

    def _near(self, box: Box) -> tuple[_Primitive, ...]:
        """Primitives not provably farther than ``r`` from the box.

        The box ends and ``r`` are doubles and IEEE subtraction is correctly
        rounded, hence monotone: ``fl(x - y) > r`` implies ``x - y > r``, so a
        skipped primitive is provably farther than ``r`` (the gouge audit's
        argument).
        """
        b = self.boxes
        gap = np.maximum.reduce([b[:, 0] - box.xmax, box.xmin - b[:, 2], b[:, 1] - box.ymax, box.ymin - b[:, 3]])
        return tuple(self.primitives[int(index)] for index in np.nonzero(~(gap > float(self.tool_radius)))[0])

    def _verdict(self, decided: list[tuple[int, float, int]]) -> ClearanceVerdict:
        if not decided:
            return _CLEAR
        worst = min(sign for sign, _, _ in decided)
        distance, index = min((distance, index) for sign, distance, index in decided if sign == worst)
        depth = float(self.tool_radius) - distance if worst < 0 else 0.0
        return ClearanceVerdict(Clearance(worst), index, depth)

    def circle(self, element: CircleElement) -> ClearanceVerdict:
        """Exact ``sign(dist(circle, ∂P) - r)``: ``max(lo - R, R - hi, 0)`` against ``r`` per primitive.

        Raises:
            NonPositiveClearanceRadiusError: If the circle radius is not positive.
        """
        if element.radius <= 0:
            raise NonPositiveClearanceRadiusError(f"Circle radius {element.radius} is not positive.")
        outer, inner = element.radius + self.tool_radius, element.radius - self.tool_radius
        outer_iv, inner_iv = Interval.of(outer), Interval.of(inner)
        centre_iv = _pair(element.centre.x, element.centre.y)
        decided: list[tuple[int, float, int]] = []
        for primitive in self._near(element_box(element)):
            near, far_low, far = primitive.circle_bounds(centre_iv)
            # Filter: provably CLEAR when lo > R + r or hi < R - r on the enclosures.
            if near.lo > outer_iv.hi or far.hi < inner_iv.lo:
                continue
            lows = primitive.nearest(element.centre)
            low_sign = min(low.sign_minus(outer) for low in lows)
            if low_sign > 0:
                continue
            low = min(item.approx() for item in lows)
            if low_sign == 0:
                # lo = R + r forces hi >= lo > R - r: the far side cannot clear it, a touch.
                decided.append((Clearance.TOUCH, max(low - float(element.radius), 0.0), primitive.index))
                continue
            # lo < R + r: clear only when the whole primitive sits inside the hole (hi < R - r).
            if far_low.lo > inner_iv.hi:
                high_sign = 1
                high = far_low.lo * _INTERVAL_SCALE
            else:
                highs = primitive.farthest(element.centre)
                high_sign = max(high.sign_minus(inner) for high in highs)
                high = max(item.approx() for item in highs)
            sign = max(low_sign, -high_sign, -1)
            if sign > 0:
                continue
            decided.append((sign, max(low - float(element.radius), float(element.radius) - high, 0.0), primitive.index))
        return self._verdict(decided)

    def circles(self, elements: Sequence[CircleElement]) -> list[ClearanceVerdict]:
        """:meth:`circle` of each element, in order (the :class:`CircleClearanceOracle` batch).

        Raises:
            NonPositiveClearanceRadiusError: If a circle radius is not positive.
        """
        return [self.circle(element) for element in elements]

    def element(self, element: Element) -> ClearanceVerdict:
        """Exact ``sign(dist(element, ∂P) - r)`` of a segment or arc element (precondition: in the closed pocket)."""
        decided: list[tuple[int, float, int]] = []
        threshold = Interval.of(self.tool_radius)
        circle = None
        segment = None
        if isinstance(element, ArcElement):
            circle = (_pair(element.centre.x, element.centre.y), Interval.of(element.radius))
        else:
            segment = _SegmentBox.of(element)
        for primitive in self._near(element_box(element)):
            if circle is not None:
                # Filter: the arc's full circle is never farther from a primitive than the arc.
                near, _, far = primitive.circle_bounds(circle[0])
                lower: Interval | None = maximum([near - circle[1], circle[1] - far])
            else:
                lower = None if segment is None else _segment_lower(segment, primitive)
            if lower is not None and lower.lo > threshold.hi:
                continue
            crosses = _segment_crosses(element, primitive) if isinstance(element, SegmentElement) else _arc_crosses(element, primitive)
            if crosses:
                decided.append((Clearance.GOUGE, 0.0, primitive.index))
                continue
            candidates = _segment_candidates(element, primitive) if isinstance(element, SegmentElement) else _arc_candidates(element, primitive)
            sign = min(candidate.sign_minus(self.tool_radius) for candidate in candidates)
            if sign > 0:
                continue
            decided.append((sign, min(candidate.approx() for candidate in candidates), primitive.index))
        return self._verdict(decided)

    def primitive_kind(self, index: int) -> PrimitiveKind:
        """``"line"``, ``"convex-arc"`` or ``"concave-arc"`` (the gouge audit's vocabulary)."""
        primitive = self.primitives[index]
        if isinstance(primitive, _Line):
            return "line"
        return "convex-arc" if primitive.sense == 1 else "concave-arc"
