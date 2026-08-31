from __future__ import annotations

import math
import sys
from collections.abc import Sequence
from fractions import Fraction
from typing import TypeAlias

from benchmarks.errors import InvalidReferenceProjectionError

BINARY64_UNIT_ROUNDOFF = sys.float_info.epsilon / 2.0
# Coefficient dot/sum paths use at most 3 roundings, term paths 4, and the two
# final accumulations make 6; gamma(7) reserves one outward-sum conversion.
BIARC_POLYNOMIAL_OPERATION_COUNT = 7
MAX_ROOT_ENCLOSURE_EXPANSIONS = 64

_XY: TypeAlias = tuple[float, float]
_FractionXY: TypeAlias = tuple[Fraction, Fraction]


def gamma(operation_count: int) -> float:
    accumulated = operation_count * BINARY64_UNIT_ROUNDOFF
    return accumulated / (1.0 - accumulated)


def outward_product(left: float, right: float) -> float:
    return math.nextafter(left * right, math.inf)


def outward_sum(left: float, right: float) -> float:
    return math.nextafter(left + right, math.inf)


def fraction_xy(point: _XY) -> _FractionXY:
    return Fraction.from_float(point[0]), Fraction.from_float(point[1])


def outward_sqrt_fraction(value: Fraction) -> float:
    represented = math.sqrt(float(value))
    while Fraction.from_float(represented) ** 2 < value:
        represented = math.nextafter(represented, math.inf)
    return represented


def certify_biarc_root(
    quadratic: float,
    linear_half: float,
    chord_squared: float,
    distance: float,
) -> float | None:
    values = (quadratic, linear_half, chord_squared, distance)
    if any(not _finite_normal_or_zero(value) for value in values) or distance < sys.float_info.min:
        return None
    coefficient_a = Fraction.from_float(quadratic)
    coefficient_b = Fraction.from_float(2.0 * linear_half)
    coefficient_c = -Fraction.from_float(chord_squared)
    represented_distance = Fraction.from_float(distance)
    residual = coefficient_a * represented_distance**2 + coefficient_b * represented_distance + coefficient_c
    absolute_sum = abs(coefficient_a) * represented_distance**2 + abs(coefficient_b) * represented_distance + abs(coefficient_c)
    unit_roundoff = Fraction.from_float(BINARY64_UNIT_ROUNDOFF)
    accumulated = BIARC_POLYNOMIAL_OPERATION_COUNT * unit_roundoff
    backward_bound = accumulated / (1 - accumulated) * absolute_sum
    if abs(residual) > backward_bound:
        return None

    derivative = 2 * coefficient_a * represented_distance + coefficient_b
    if derivative == 0:
        return None
    initial_uncertainty = float(backward_bound / abs(derivative))
    if initial_uncertainty == 0.0:
        initial_uncertainty = math.ulp(distance)
    uncertainty = math.nextafter(initial_uncertainty, math.inf)
    for _ in range(MAX_ROOT_ENCLOSURE_EXPANSIONS):
        exact_uncertainty = Fraction.from_float(uncertainty)
        lower = represented_distance - exact_uncertainty
        upper = represented_distance + exact_uncertainty
        if lower > 0:
            lower_value = coefficient_a * lower**2 + coefficient_b * lower + coefficient_c
            upper_value = coefficient_a * upper**2 + coefficient_b * upper + coefficient_c
            if lower_value * upper_value <= 0:
                return uncertainty
        uncertainty *= 2.0
        if not math.isfinite(uncertainty):
            return None
    return None


def exact_reflection(tangent: _XY, chord: _XY) -> _FractionXY:
    tangent_exact = fraction_xy(tangent)
    chord_exact = fraction_xy(chord)
    chord_squared = chord_exact[0] ** 2 + chord_exact[1] ** 2
    projection = tangent_exact[0] * chord_exact[0] + tangent_exact[1] * chord_exact[1]
    return (
        2 * projection * chord_exact[0] / chord_squared - tangent_exact[0],
        2 * projection * chord_exact[1] / chord_squared - tangent_exact[1],
    )


def exact_vector_error(represented: _XY, exact: _FractionXY) -> float:
    differences = (
        Fraction.from_float(represented[0]) - exact[0],
        Fraction.from_float(represented[1]) - exact[1],
    )
    return outward_sqrt_fraction(differences[0] ** 2 + differences[1] ** 2)


def exact_segment_distance_squared(point: _XY, start: _XY, end: _XY) -> Fraction:
    point_exact = fraction_xy(point)
    start_exact = fraction_xy(start)
    end_exact = fraction_xy(end)
    segment = end_exact[0] - start_exact[0], end_exact[1] - start_exact[1]
    relative = point_exact[0] - start_exact[0], point_exact[1] - start_exact[1]
    length_squared = segment[0] ** 2 + segment[1] ** 2
    projection_numerator = relative[0] * segment[0] + relative[1] * segment[1]
    if projection_numerator <= 0:
        return relative[0] ** 2 + relative[1] ** 2
    if projection_numerator >= length_squared:
        end_relative = point_exact[0] - end_exact[0], point_exact[1] - end_exact[1]
        return end_relative[0] ** 2 + end_relative[1] ** 2
    cross = segment[0] * relative[1] - segment[1] * relative[0]
    return cross**2 / length_squared


def control_hull_radius_bounds(
    control_points: tuple[_XY, _XY, _XY, _XY],
    centre: _XY,
) -> tuple[float, float]:
    hull = _convex_hull(control_points)
    maximum = max(math.dist(point, centre) for point in hull)
    if _point_in_convex_polygon(centre, hull):
        return 0.0, maximum
    minimum = min(point_segment_distance(centre, start, end) for start, end in zip(hull, (*hull[1:], hull[0])))
    return minimum, maximum


def _convex_hull(points: Sequence[_XY]) -> tuple[_XY, ...]:
    ordered = sorted(set(points))
    if len(ordered) <= 1:
        return tuple(ordered)

    def half(sequence: Sequence[_XY]) -> list[_XY]:
        result: list[_XY] = []
        for point in sequence:
            while len(result) >= 2 and _cross(_subtract(result[-1], result[-2]), _subtract(point, result[-1])) <= 0.0:
                result.pop()
            result.append(point)
        return result

    return tuple(half(ordered)[:-1] + half(tuple(reversed(ordered)))[:-1])


def _point_in_convex_polygon(point: _XY, polygon: Sequence[_XY]) -> bool:
    if len(polygon) < 3:
        return False
    signs = [_cross(_subtract(end, start), _subtract(point, start)) for start, end in zip(polygon, (*polygon[1:], polygon[0]))]
    return all(value >= 0.0 for value in signs) or all(value <= 0.0 for value in signs)


def point_segment_distance(point: _XY, start: _XY, end: _XY) -> float:
    segment = _subtract(end, start)
    length_squared = _dot(segment, segment)
    if length_squared == 0.0:
        return math.dist(point, start)
    parameter = max(0.0, min(1.0, _dot(_subtract(point, start), segment) / length_squared))
    projected = start[0] + segment[0] * parameter, start[1] + segment[1] * parameter
    return math.dist(point, projected)


def maximum_arc_segment_distance(
    start_radius: _XY,
    start_offset: float,
    end_offset: float,
    segment_start: _XY,
    segment_end: _XY,
) -> float:
    segment = _subtract(segment_end, segment_start)
    if _dot(segment, segment) == 0.0:
        raise InvalidReferenceProjectionError("An emitted projection chord collapsed to one point.")

    segment_quarter_turn = -segment[1], segment[0]
    start_quarter_turn = -segment_start[1], segment_start[0]
    end_quarter_turn = -segment_end[1], segment_end[0]
    candidates = {start_offset, end_offset}
    equations = (
        (segment, _dot(segment, segment_start)),
        (segment, _dot(segment, segment_end)),
        (start_quarter_turn, 0.0),
        (end_quarter_turn, 0.0),
        (segment_quarter_turn, _dot(segment_quarter_turn, segment_start)),
        (segment, 0.0),
    )
    for normal, value in equations:
        candidates.update(
            _linear_circle_root_offsets(
                start_radius,
                normal,
                value,
                start_offset,
                end_offset,
            )
        )
    return max(
        point_segment_distance(
            _rotate_vector(start_radius, offset),
            segment_start,
            segment_end,
        )
        for offset in candidates
    )


def _linear_circle_root_offsets(
    start_radius: _XY,
    normal: _XY,
    value: float,
    start_offset: float,
    end_offset: float,
) -> tuple[float, ...]:
    quarter_turn = -start_radius[1], start_radius[0]
    cosine_coefficient = _dot(normal, start_radius)
    sine_coefficient = _dot(normal, quarter_turn)
    amplitude = math.hypot(cosine_coefficient, sine_coefficient)
    if amplitude == 0.0 or value < -amplitude or value > amplitude:
        return ()

    phase = math.atan2(sine_coefficient, cosine_coefficient)
    root_offset = math.acos(max(-1.0, min(1.0, value / amplitude)))
    lower = min(start_offset, end_offset)
    upper = max(start_offset, end_offset)
    roots: set[float] = set()
    for base in (phase - root_offset, phase + root_offset):
        minimum_turn = math.ceil((lower - base) / math.tau)
        maximum_turn = math.floor((upper - base) / math.tau)
        roots.update(base + turn * math.tau for turn in range(minimum_turn, maximum_turn + 1))
    return tuple(sorted(roots))


def _rotate_vector(vector: _XY, angle: float) -> _XY:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return cosine * vector[0] - sine * vector[1], sine * vector[0] + cosine * vector[1]


def _subtract(left: _XY, right: _XY) -> _XY:
    return left[0] - right[0], left[1] - right[1]


def _dot(left: _XY, right: _XY) -> float:
    return left[0] * right[0] + left[1] * right[1]


def _cross(left: _XY, right: _XY) -> float:
    return left[0] * right[1] - left[1] * right[0]


def _finite_normal_or_zero(value: float) -> bool:
    return math.isfinite(value) and (value == 0.0 or abs(value) >= sys.float_info.min)
