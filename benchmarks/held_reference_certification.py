from __future__ import annotations

import math
import sys
from collections.abc import Sequence
from fractions import Fraction
from functools import cache
from typing import TypeAlias

from benchmarks.errors import InvalidReferenceProjectionError

BINARY64_UNIT_ROUNDOFF = sys.float_info.epsilon / 2.0
# Coefficient dot/sum paths use at most 3 roundings, term paths 4, and the two
# final accumulations make 6; gamma(7) reserves one outward-sum conversion.
BIARC_POLYNOMIAL_OPERATION_COUNT = 7
MAX_ROOT_ENCLOSURE_EXPANSIONS = 64

_XY: TypeAlias = tuple[float, float]
_FractionXY: TypeAlias = tuple[Fraction, Fraction]
_ArcProof: TypeAlias = tuple[_XY, _XY, _XY, float]


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
    if value < 0:
        raise ArithmeticError("Cannot enclose a negative square root.")
    represented = math.sqrt(float(value))
    if not math.isfinite(represented):
        raise ArithmeticError("Square-root enclosure is not finite.")
    while Fraction.from_float(represented) ** 2 < value:
        represented = math.nextafter(represented, math.inf)
    return represented


def inward_sqrt_fraction(value: Fraction) -> float:
    if value < 0:
        raise ArithmeticError("Cannot enclose a negative square root.")
    represented = math.sqrt(float(value))
    if not math.isfinite(represented):
        raise ArithmeticError("Square-root enclosure is not finite.")
    while Fraction.from_float(represented) ** 2 > value:
        represented = math.nextafter(represented, -math.inf)
    return represented


def biarc_solver_coefficients(
    start_tangent: _XY,
    end_tangent: _XY,
    chord: _XY,
) -> tuple[float, float, float] | None:
    coefficients = _biarc_polynomial(start_tangent, end_tangent, chord)
    if coefficients is None:
        return None
    coefficient_a, coefficient_b, coefficient_c = coefficients
    represented = float(coefficient_a), float(coefficient_b / 2), float(-coefficient_c)
    if any(not _finite_normal_or_zero(value) for value in represented):
        return None
    return represented


def certify_biarc_root(
    start_tangent: _XY,
    end_tangent: _XY,
    chord: _XY,
    distance: float,
) -> float | None:
    coefficients = _biarc_polynomial(start_tangent, end_tangent, chord)
    if coefficients is None or not _finite_normal_or_zero(distance) or distance <= 0.0:
        return None
    coefficient_a, coefficient_b, coefficient_c = coefficients
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
    initial_uncertainty = max(float(backward_bound / abs(derivative)), math.ulp(distance))
    if not math.isfinite(initial_uncertainty):
        return None
    uncertainty = math.nextafter(initial_uncertainty, math.inf)
    for _ in range(MAX_ROOT_ENCLOSURE_EXPANSIONS):
        exact_uncertainty = Fraction.from_float(uncertainty)
        lower = represented_distance - exact_uncertainty
        upper = represented_distance + exact_uncertainty
        if lower > 0:
            lower_value = coefficient_a * lower**2 + coefficient_b * lower + coefficient_c
            upper_value = coefficient_a * upper**2 + coefficient_b * upper + coefficient_c
            if lower_value <= 0 <= upper_value:
                return uncertainty
        uncertainty *= 2.0
        if not math.isfinite(uncertainty):
            return None
    return None


def _biarc_polynomial(
    start_tangent: _XY,
    end_tangent: _XY,
    chord: _XY,
) -> tuple[Fraction, Fraction, Fraction] | None:
    components = (*start_tangent, *end_tangent, *chord)
    if any(not _finite_normal_or_zero(value) for value in components):
        return None
    start_exact = fraction_xy(start_tangent)
    end_exact = fraction_xy(end_tangent)
    chord_exact = fraction_xy(chord)
    tangent_dot = _fraction_dot(start_exact, end_exact)
    coefficient_a = 2 * (1 - tangent_dot)
    coefficient_b = 2 * _fraction_dot(
        chord_exact,
        (start_exact[0] + end_exact[0], start_exact[1] + end_exact[1]),
    )
    coefficient_c = -_fraction_dot(chord_exact, chord_exact)
    if coefficient_a <= 0 or coefficient_c >= 0:
        return None
    return coefficient_a, coefficient_b, coefficient_c


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


def affine_join_error(
    stored_join: _XY,
    chord: _XY,
    tangent_difference: _XY,
    distance: float,
    uncertainty: float,
) -> float | None:
    values = (*stored_join, *chord, *tangent_difference, distance, uncertainty)
    if any(not _finite_normal_or_zero(value) for value in values) or uncertainty < 0.0:
        return None
    stored_exact = fraction_xy(stored_join)
    chord_exact = fraction_xy(chord)
    tangent_exact = fraction_xy(tangent_difference)
    distance_exact = Fraction.from_float(distance)
    uncertainty_exact = Fraction.from_float(uncertainty)
    component_errors: list[Fraction] = []
    for index in range(2):
        endpoints = tuple((chord_exact[index] + root * tangent_exact[index]) / 2 for root in (distance_exact - uncertainty_exact, distance_exact + uncertainty_exact))
        component_errors.append(max(abs(stored_exact[index] - endpoint) for endpoint in endpoints))
    try:
        return outward_sqrt_fraction(component_errors[0] ** 2 + component_errors[1] ** 2)
    except (ArithmeticError, OverflowError, ValueError):
        return None


def reflection_perturbation_bound(join_error: float, chord: _XY) -> float | None:
    if not _finite_normal_or_zero(join_error) or join_error < 0.0:
        return None
    chord_squared = _fraction_dot(fraction_xy(chord), fraction_xy(chord))
    try:
        chord_lower = inward_sqrt_fraction(chord_squared)
    except (ArithmeticError, OverflowError, ValueError):
        return None
    denominator = Fraction.from_float(chord_lower) - Fraction.from_float(join_error)
    if denominator <= 0:
        return None
    return _outward_float_fraction(2 * Fraction.from_float(join_error) / denominator)


def mapped_radius_tangent_bound(
    local_point: _XY,
    local_centre: _XY,
    mapped_point: _XY,
    mapped_centre: _XY,
    scale: float,
) -> float | None:
    values = (*local_point, *local_centre, *mapped_point, *mapped_centre, scale)
    if any(not _finite_normal_or_zero(value) for value in values) or scale <= 0.0:
        return None
    scale_exact = Fraction.from_float(scale)
    local_radius = _fraction_subtract(fraction_xy(local_point), fraction_xy(local_centre))
    ideal_radius = scale_exact * local_radius[0], scale_exact * local_radius[1]
    stored_radius = _fraction_subtract(fraction_xy(mapped_point), fraction_xy(mapped_centre))
    difference = _fraction_subtract(stored_radius, ideal_radius)
    try:
        error_upper = outward_sqrt_fraction(_fraction_dot(difference, difference))
        radius_lower = inward_sqrt_fraction(_fraction_dot(ideal_radius, ideal_radius))
    except (ArithmeticError, OverflowError, ValueError):
        return None
    denominator = Fraction.from_float(radius_lower) - Fraction.from_float(error_upper)
    if denominator <= 0:
        return None
    return _outward_float_fraction(2 * Fraction.from_float(error_upper) / denominator)


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


def exact_control_hull_radius_squared_bounds(
    control_points: Sequence[_XY],
    centre: _XY,
) -> tuple[Fraction, Fraction]:
    hull = tuple(fraction_xy(point) for point in exact_convex_hull(control_points))
    return _exact_control_hull_radius_squared_bounds(hull, fraction_xy(centre))


def _exact_control_hull_radius_squared_bounds(
    control_points: Sequence[_FractionXY],
    centre: _FractionXY,
) -> tuple[Fraction, Fraction]:
    hull_exact = _exact_convex_hull_fraction(control_points)
    maximum = max(_fraction_distance_squared(point, centre) for point in hull_exact)
    if _exact_point_in_convex_polygon(centre, hull_exact):
        return Fraction(0), maximum
    minimum = min(_exact_segment_distance_squared(centre, start, end) for start, end in zip(hull_exact, (*hull_exact[1:], hull_exact[0])))
    return minimum, maximum


def circle_deviation_upper_bound(
    control_points: Sequence[_XY],
    centre: _XY,
    radius_point: _XY,
) -> float | None:
    try:
        minimum_squared, maximum_squared = exact_control_hull_radius_squared_bounds(control_points, centre)
        radius_squared = _fraction_distance_squared(fraction_xy(radius_point), fraction_xy(centre))
        radius_lower = inward_sqrt_fraction(radius_squared)
        radius_upper = outward_sqrt_fraction(radius_squared)
        minimum_lower = inward_sqrt_fraction(minimum_squared)
        maximum_upper = outward_sqrt_fraction(maximum_squared)
    except (ArithmeticError, OverflowError, ValueError):
        return None
    bound = max(radius_upper - minimum_lower, maximum_upper - radius_lower)
    if not math.isfinite(bound):
        return None
    return math.nextafter(bound, math.inf)


def certify_cubic_interval_circle(
    control_points: Sequence[_XY],
    start_parameter: Fraction,
    end_parameter: Fraction,
    centre: _XY,
    radius_point: _XY,
    direction: float,
    deviation_limit: float,
    maximum_depth: int,
) -> float | None:
    if len(control_points) != 4 or not 0 <= start_parameter < end_parameter <= 1 or direction not in (-1.0, 1.0) or not math.isfinite(deviation_limit) or deviation_limit < 0:
        return None
    exact_controls = tuple(fraction_xy(point) for point in control_points)
    interval_controls = _exact_cubic_interval(
        (exact_controls[0], exact_controls[1], exact_controls[2], exact_controls[3]),
        start_parameter,
        end_parameter,
    )
    centre_exact = fraction_xy(centre)
    if not _certify_monotone_polar_angle(
        interval_controls,
        centre_exact,
        1 if direction > 0 else -1,
        maximum_depth,
        depth=0,
    ):
        return None
    radius_squared = _fraction_distance_squared(fraction_xy(radius_point), centre_exact)
    return _exact_circle_interval_bound(
        interval_controls,
        centre_exact,
        radius_squared,
        Fraction.from_float(deviation_limit),
        maximum_depth,
        depth=0,
    )


def _exact_circle_interval_bound(
    control_points: tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    centre: _FractionXY,
    radius_squared: Fraction,
    deviation_limit: Fraction,
    maximum_depth: int,
    *,
    depth: int,
) -> float | None:
    minimum_squared, maximum_squared = _exact_control_hull_radius_squared_bounds(control_points, centre)
    radius_lower = inward_sqrt_fraction(radius_squared)
    radius_upper = outward_sqrt_fraction(radius_squared)
    bound = math.nextafter(
        max(
            radius_upper - inward_sqrt_fraction(minimum_squared),
            outward_sqrt_fraction(maximum_squared) - radius_lower,
        ),
        math.inf,
    )
    if Fraction.from_float(bound) <= deviation_limit:
        return bound
    if depth >= maximum_depth:
        return None
    left, right = _exact_split_cubic(control_points)
    left_bound = _exact_circle_interval_bound(
        left,
        centre,
        radius_squared,
        deviation_limit,
        maximum_depth,
        depth=depth + 1,
    )
    if left_bound is None:
        return None
    right_bound = _exact_circle_interval_bound(
        right,
        centre,
        radius_squared,
        deviation_limit,
        maximum_depth,
        depth=depth + 1,
    )
    return None if right_bound is None else max(left_bound, right_bound)


def exact_convex_hull(points: Sequence[_XY]) -> tuple[_XY, ...]:
    ordered = sorted(set(points))
    if len(ordered) <= 1:
        return tuple(ordered)

    def half(sequence: Sequence[_XY]) -> list[_XY]:
        result: list[_XY] = []
        for point in sequence:
            while len(result) >= 2 and _exact_turn(result[-2], result[-1], point) <= 0:
                result.pop()
            result.append(point)
        return result

    return tuple(half(ordered)[:-1] + half(tuple(reversed(ordered)))[:-1])


def exact_cubic_point(
    control_points: Sequence[_XY],
    parameter: Fraction,
) -> _FractionXY:
    if len(control_points) != 4 or not 0 <= parameter <= 1:
        raise ArithmeticError("Cubic evaluation requires four controls and a unit parameter.")
    level = list(map(fraction_xy, control_points))
    complement = 1 - parameter
    while len(level) > 1:
        level = [
            (
                complement * left[0] + parameter * right[0],
                complement * left[1] + parameter * right[1],
            )
            for left, right in zip(level, level[1:])
        ]
    return level[0]


def exact_polar_bernstein_coefficients(
    control_points: Sequence[_XY],
    centre: _XY,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    if len(control_points) != 4:
        raise ArithmeticError("Polar certification requires four cubic controls.")
    exact_controls = tuple(fraction_xy(control_points[index]) for index in range(4))
    return _exact_polar_bernstein_coefficients(
        (exact_controls[0], exact_controls[1], exact_controls[2], exact_controls[3]),
        fraction_xy(centre),
    )


def certify_monotone_polar_angle(
    control_points: Sequence[_XY],
    centre: _XY,
    direction: float,
    maximum_depth: int,
) -> bool:
    if len(control_points) != 4 or direction not in (-1.0, 1.0):
        return False
    exact_controls = tuple(fraction_xy(control_points[index]) for index in range(4))
    return _certify_monotone_polar_angle(
        (exact_controls[0], exact_controls[1], exact_controls[2], exact_controls[3]),
        fraction_xy(centre),
        1 if direction > 0 else -1,
        maximum_depth,
        depth=0,
    )


def _certify_monotone_polar_angle(
    control_points: tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    centre: _FractionXY,
    direction: int,
    maximum_depth: int,
    *,
    depth: int,
) -> bool:
    coefficients = _exact_polar_bernstein_coefficients(control_points, centre)
    directed = tuple(direction * coefficient for coefficient in coefficients)
    minimum_squared, _ = _exact_control_hull_radius_squared_bounds(control_points, centre)
    if min(directed) >= 0 and minimum_squared > 0:
        return True
    midpoint_value = sum(Fraction(math.comb(5, index), 32) * coefficient for index, coefficient in enumerate(directed))
    if midpoint_value < 0 or depth >= maximum_depth:
        return False
    left, right = _exact_split_cubic(control_points)
    return _certify_monotone_polar_angle(
        left,
        centre,
        direction,
        maximum_depth,
        depth=depth + 1,
    ) and _certify_monotone_polar_angle(
        right,
        centre,
        direction,
        maximum_depth,
        depth=depth + 1,
    )


def _exact_polar_bernstein_coefficients(
    exact_controls: tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    centre_exact: _FractionXY,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    relative = tuple(_fraction_subtract(point, centre_exact) for point in exact_controls)
    derivatives = tuple(
        (
            3 * (exact_controls[index + 1][0] - exact_controls[index][0]),
            3 * (exact_controls[index + 1][1] - exact_controls[index][1]),
        )
        for index in range(3)
    )
    coefficients: list[Fraction] = []
    for degree in range(6):
        coefficient = Fraction(0)
        for cubic_index in range(max(0, degree - 2), min(3, degree) + 1):
            derivative_index = degree - cubic_index
            weight = Fraction(
                math.comb(3, cubic_index) * math.comb(2, derivative_index),
                math.comb(5, degree),
            )
            coefficient += weight * _fraction_cross(relative[cubic_index], derivatives[derivative_index])
        coefficients.append(coefficient)
    return (
        coefficients[0],
        coefficients[1],
        coefficients[2],
        coefficients[3],
        coefficients[4],
        coefficients[5],
    )


def _exact_split_cubic(
    control_points: tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
) -> tuple[
    tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
]:
    start, control1, control2, end = control_points
    first = _fraction_midpoint(start, control1)
    second = _fraction_midpoint(control1, control2)
    third = _fraction_midpoint(control2, end)
    fourth = _fraction_midpoint(first, second)
    fifth = _fraction_midpoint(second, third)
    midpoint = _fraction_midpoint(fourth, fifth)
    return (start, first, fourth, midpoint), (midpoint, fifth, third, end)


def _exact_split_cubic_at(
    control_points: tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    parameter: Fraction,
) -> tuple[
    tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
]:
    start, control1, control2, end = control_points
    complement = 1 - parameter

    def interpolate(left: _FractionXY, right: _FractionXY) -> _FractionXY:
        return (
            complement * left[0] + parameter * right[0],
            complement * left[1] + parameter * right[1],
        )

    first = interpolate(start, control1)
    second = interpolate(control1, control2)
    third = interpolate(control2, end)
    fourth = interpolate(first, second)
    fifth = interpolate(second, third)
    point = interpolate(fourth, fifth)
    return (start, first, fourth, point), (point, fifth, third, end)


def _exact_cubic_interval(
    control_points: tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY],
    start_parameter: Fraction,
    end_parameter: Fraction,
) -> tuple[_FractionXY, _FractionXY, _FractionXY, _FractionXY]:
    if start_parameter == 0 and end_parameter == 1:
        return control_points
    left, _ = _exact_split_cubic_at(control_points, end_parameter)
    if start_parameter == 0:
        return left
    _, interval = _exact_split_cubic_at(left, start_parameter / end_parameter)
    return interval


def biarc_branch_at_parameter(
    first_radius_squared: Fraction,
    first_sweep_squared: Fraction,
    second_radius_squared: Fraction,
    second_sweep_squared: Fraction,
    parameter: Fraction,
) -> int:
    if min(first_radius_squared, first_sweep_squared, second_radius_squared, second_sweep_squared) <= 0:
        raise ArithmeticError("Biarc branch lengths must be strictly positive.")
    if not 0 <= parameter <= 1:
        raise ArithmeticError("Biarc correspondence parameter must lie in the unit interval.")
    first_length_squared = first_radius_squared * first_sweep_squared
    second_length_squared = second_radius_squared * second_sweep_squared
    return 1 if parameter**2 * second_length_squared <= (1 - parameter) ** 2 * first_length_squared else 2


def audited_sin_cos(angle: Fraction) -> tuple[float, float, float, float]:
    represented_angle = float(angle)
    if not math.isfinite(represented_angle):
        raise ArithmeticError("Trigonometric argument is not finite.")
    exact_argument = Fraction.from_float(represented_angle)
    sine_polynomial = sum(
        ((-1) ** index * exact_argument ** (2 * index + 1) / math.factorial(2 * index + 1) for index in range(20)),
        Fraction(0),
    )
    cosine_polynomial = sum(
        ((-1) ** index * exact_argument ** (2 * index) / math.factorial(2 * index) for index in range(21)),
        Fraction(0),
    )
    remainder = abs(exact_argument) ** 41 / math.factorial(41)
    angle_radius = abs(angle - exact_argument)
    represented_sine = math.sin(represented_angle)
    represented_cosine = math.cos(represented_angle)
    if not math.isfinite(represented_sine) or not math.isfinite(represented_cosine):
        raise ArithmeticError("Trigonometric result is not finite.")
    sine_error = abs(Fraction.from_float(represented_sine) - sine_polynomial) + remainder + angle_radius
    cosine_error = abs(Fraction.from_float(represented_cosine) - cosine_polynomial) + remainder + angle_radius
    return (
        represented_sine,
        outward_sqrt_fraction(sine_error**2),
        represented_cosine,
        outward_sqrt_fraction(cosine_error**2),
    )


def closed_arc_point_box(
    arc: _ArcProof,
    parameter: Fraction,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    start, end, centre, sweep = arc
    if parameter == 0:
        exact = fraction_xy(start)
        return (exact[0], exact[0]), (exact[1], exact[1])
    if parameter == 1:
        exact = fraction_xy(end)
        return (exact[0], exact[0]), (exact[1], exact[1])
    if not 0 < parameter < 1:
        raise ArithmeticError("Closed arc parameter must lie in the unit interval.")
    ideal_box = _ideal_arc_point_box(arc, parameter)
    ideal_end_box = _ideal_arc_point_box(arc, Fraction(1))
    stored_end = fraction_xy(end)
    result: list[tuple[Fraction, Fraction]] = []
    for coordinate in range(2):
        ideal_lower, ideal_upper = ideal_box[coordinate]
        endpoint_lower, endpoint_upper = ideal_end_box[coordinate]
        correction_lower = stored_end[coordinate] - endpoint_upper
        correction_upper = stored_end[coordinate] - endpoint_lower
        result.append(
            (
                ideal_lower + parameter * correction_lower,
                ideal_upper + parameter * correction_upper,
            )
        )
    return result[0], result[1]


@cache
def _ideal_arc_point_box(
    arc: _ArcProof,
    parameter: Fraction,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    start, _, centre, sweep = arc
    start_exact = fraction_xy(start)
    centre_exact = fraction_xy(centre)
    radius = _fraction_subtract(start_exact, centre_exact)
    sine, sine_error, cosine, cosine_error = audited_sin_cos(Fraction.from_float(sweep) * parameter)
    sine_exact = Fraction.from_float(sine)
    cosine_exact = Fraction.from_float(cosine)
    proxy = (
        centre_exact[0] + cosine_exact * radius[0] - sine_exact * radius[1],
        centre_exact[1] + sine_exact * radius[0] + cosine_exact * radius[1],
    )
    error = Fraction.from_float(max(sine_error, cosine_error)) * (abs(radius[0]) + abs(radius[1]))
    return (proxy[0] - error, proxy[0] + error), (proxy[1] - error, proxy[1] + error)


def biarc_correspondence_node_bound(
    control_points: Sequence[_XY],
    biarc: tuple[_ArcProof, _ArcProof],
    start_parameter: Fraction,
    end_parameter: Fraction,
) -> float | None:
    if len(control_points) != 4 or not 0 <= start_parameter < end_parameter <= 1:
        return None
    try:
        length_intervals = (_arc_length_interval(biarc[0]), _arc_length_interval(biarc[1]))
        if any(lower <= 0 for lower, _ in length_intervals):
            return None
        total_interval = (
            length_intervals[0][0] + length_intervals[1][0],
            length_intervals[0][1] + length_intervals[1][1],
        )
        endpoint_corrections = (
            _arc_endpoint_correction_upper(biarc[0]),
            _arc_endpoint_correction_upper(biarc[1]),
        )
        locus_correction = max(endpoint_corrections)
        cubic_speed_upper = _cubic_speed_upper(control_points)
        start_branch = _biarc_branch(biarc, start_parameter)
        end_branch = _biarc_branch(biarc, end_parameter)
        if start_branch != end_branch:
            source_to_closed = _first_order_biarc_node_bound(
                control_points,
                biarc,
                length_intervals,
                total_interval,
                endpoint_corrections,
                cubic_speed_upper,
                start_parameter,
                end_parameter,
            )
        else:
            source_to_closed = _second_order_biarc_node_bound(
                control_points,
                biarc,
                length_intervals,
                total_interval,
                endpoint_corrections,
                cubic_speed_upper,
                start_parameter,
                end_parameter,
                start_branch,
            )
        if source_to_closed is None:
            return None
        bound = source_to_closed + locus_correction
        return _outward_float_fraction(bound)
    except (ArithmeticError, OverflowError, ValueError, ZeroDivisionError):
        return None


def biarc_correspondence_sample_exceeds(
    control_points: Sequence[_XY],
    biarc: tuple[_ArcProof, _ArcProof],
    parameter: Fraction,
    deviation_limit: float,
) -> bool:
    """Prove that one paired source/biarc sample violates the limit."""
    if len(control_points) != 4 or not 0 <= parameter <= 1 or not math.isfinite(deviation_limit) or deviation_limit < 0:
        return False
    try:
        length_intervals = (_arc_length_interval(biarc[0]), _arc_length_interval(biarc[1]))
        if any(lower <= 0 for lower, _ in length_intervals):
            return False
        total_interval = (
            length_intervals[0][0] + length_intervals[1][0],
            length_intervals[0][1] + length_intervals[1][1],
        )
        point_box = _closed_biarc_point_box(biarc, length_intervals, total_interval, parameter)
        if point_box is None:
            return False
        minimum_squared = _point_to_box_minimum_distance_squared(
            exact_cubic_point(control_points, parameter),
            point_box,
        )
        locus_correction = max(_arc_endpoint_correction_upper(arc) for arc in biarc)
        allowed = Fraction.from_float(deviation_limit) + locus_correction
        return minimum_squared > allowed**2
    except (ArithmeticError, OverflowError, ValueError, ZeroDivisionError):
        return False


def _first_order_biarc_node_bound(
    control_points: Sequence[_XY],
    biarc: tuple[_ArcProof, _ArcProof],
    length_intervals: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    total_interval: tuple[Fraction, Fraction],
    endpoint_corrections: tuple[Fraction, Fraction],
    cubic_speed_upper: Fraction,
    start_parameter: Fraction,
    end_parameter: Fraction,
) -> Fraction | None:
    midpoint = (start_parameter + end_parameter) / 2
    point_box = _closed_biarc_point_box(biarc, length_intervals, total_interval, midpoint)
    if point_box is None:
        return None
    cubic_point = exact_cubic_point(control_points, midpoint)
    sample_upper = Fraction.from_float(outward_sqrt_fraction(_point_to_box_distance_squared(cubic_point, point_box)))
    correction_speed = max(endpoint_corrections[index] * total_interval[1] / length_intervals[index][0] for index in range(2))
    closed_speed_upper = total_interval[1] + correction_speed
    return sample_upper + (end_parameter - start_parameter) / 2 * (cubic_speed_upper + closed_speed_upper)


def _second_order_biarc_node_bound(
    control_points: Sequence[_XY],
    biarc: tuple[_ArcProof, _ArcProof],
    length_intervals: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    total_interval: tuple[Fraction, Fraction],
    endpoint_corrections: tuple[Fraction, Fraction],
    cubic_speed_upper: Fraction,
    start_parameter: Fraction,
    end_parameter: Fraction,
    branch: int,
) -> Fraction | None:
    endpoint_squared: list[Fraction] = []
    for parameter in (start_parameter, end_parameter):
        point_box = _closed_biarc_point_box(biarc, length_intervals, total_interval, parameter)
        if point_box is None:
            return None
        endpoint_squared.append(_point_to_box_distance_squared(exact_cubic_point(control_points, parameter), point_box))

    arc_index = branch - 1
    correction_speed = endpoint_corrections[arc_index] * total_interval[1] / length_intervals[arc_index][0]
    first_derivative = cubic_speed_upper + total_interval[1] + correction_speed
    cubic_second_derivative = _cubic_second_derivative_upper(control_points)
    radius_lower = Fraction.from_float(inward_sqrt_fraction(_arc_radius_squared(biarc[arc_index])))
    if radius_lower <= 0:
        return None
    second_derivative = cubic_second_derivative + total_interval[1] ** 2 / radius_lower
    position_bound = _position_difference_upper(control_points, biarc[arc_index], endpoint_corrections[arc_index])
    second_derivative_squared_distance = 2 * (first_derivative**2 + position_bound * second_derivative)
    width = end_parameter - start_parameter
    squared_bound = max(endpoint_squared) + width**2 / 8 * second_derivative_squared_distance
    return Fraction.from_float(outward_sqrt_fraction(squared_bound))


def _biarc_branch(biarc: tuple[_ArcProof, _ArcProof], parameter: Fraction) -> int:
    return biarc_branch_at_parameter(
        _arc_radius_squared(biarc[0]),
        Fraction.from_float(abs(biarc[0][3])) ** 2,
        _arc_radius_squared(biarc[1]),
        Fraction.from_float(abs(biarc[1][3])) ** 2,
        parameter,
    )


def _cubic_speed_upper(control_points: Sequence[_XY]) -> Fraction:
    return 3 * max(
        Fraction.from_float(
            outward_sqrt_fraction(
                _fraction_distance_squared(
                    fraction_xy(control_points[index]),
                    fraction_xy(control_points[index + 1]),
                )
            )
        )
        for index in range(3)
    )


def _cubic_second_derivative_upper(control_points: Sequence[_XY]) -> Fraction:
    exact = tuple(fraction_xy(point) for point in control_points)
    second_differences = tuple(
        (
            6 * (exact[index + 2][0] - 2 * exact[index + 1][0] + exact[index][0]),
            6 * (exact[index + 2][1] - 2 * exact[index + 1][1] + exact[index][1]),
        )
        for index in range(2)
    )
    return max(Fraction.from_float(outward_sqrt_fraction(_fraction_dot(vector, vector))) for vector in second_differences)


def _position_difference_upper(
    control_points: Sequence[_XY],
    arc: _ArcProof,
    endpoint_correction: Fraction,
) -> Fraction:
    centre = fraction_xy(arc[2])
    control_radius = max(Fraction.from_float(outward_sqrt_fraction(_fraction_distance_squared(fraction_xy(point), centre))) for point in control_points)
    arc_radius = Fraction.from_float(outward_sqrt_fraction(_arc_radius_squared(arc)))
    return control_radius + arc_radius + endpoint_correction


def _point_to_box_distance_squared(
    point: _FractionXY,
    box: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> Fraction:
    return sum(
        (_interval_component_error(point[index], box[index]) ** 2 for index in range(2)),
        Fraction(0),
    )


def _point_to_box_minimum_distance_squared(
    point: _FractionXY,
    box: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> Fraction:
    squared = Fraction(0)
    for value, (lower, upper) in zip(point, box):
        if value < lower:
            squared += (lower - value) ** 2
        elif value > upper:
            squared += (value - upper) ** 2
    return squared


def _closed_biarc_point_box(
    biarc: tuple[_ArcProof, _ArcProof],
    length_intervals: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    total_interval: tuple[Fraction, Fraction],
    parameter: Fraction,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]] | None:
    branch = biarc_branch_at_parameter(
        _arc_radius_squared(biarc[0]),
        Fraction.from_float(abs(biarc[0][3])) ** 2,
        _arc_radius_squared(biarc[1]),
        Fraction.from_float(abs(biarc[1][3])) ** 2,
        parameter,
    )
    if branch == 1:
        lower = parameter * total_interval[0] / length_intervals[0][1]
        upper = parameter * total_interval[1] / length_intervals[0][0]
    else:
        lower = (parameter * total_interval[0] - length_intervals[0][1]) / length_intervals[1][1]
        upper = (parameter * total_interval[1] - length_intervals[0][0]) / length_intervals[1][0]
    lower = max(Fraction(0), lower)
    upper = min(Fraction(1), upper)
    if lower > upper:
        return None
    return _closed_arc_point_interval_box(biarc[branch - 1], lower, upper)


def _closed_arc_point_interval_box(
    arc: _ArcProof,
    lower: Fraction,
    upper: Fraction,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    if lower == upper:
        return closed_arc_point_box(arc, lower)
    midpoint = (lower + upper) / 2
    radius = _fraction_subtract(fraction_xy(arc[0]), fraction_xy(arc[2]))
    angle_radius = abs(Fraction.from_float(arc[3])) * (upper - lower) / 2
    ideal_box = _ideal_arc_point_box(arc, midpoint)
    rotation_radius = (abs(radius[0]) + abs(radius[1])) * angle_radius
    ideal_box = (
        (ideal_box[0][0] - rotation_radius, ideal_box[0][1] + rotation_radius),
        (ideal_box[1][0] - rotation_radius, ideal_box[1][1] + rotation_radius),
    )
    ideal_end_box = _ideal_arc_point_box(arc, Fraction(1))
    stored_end = fraction_xy(arc[1])
    result: list[tuple[Fraction, Fraction]] = []
    for coordinate in range(2):
        correction = (
            stored_end[coordinate] - ideal_end_box[coordinate][1],
            stored_end[coordinate] - ideal_end_box[coordinate][0],
        )
        correction_interval = _interval_product((lower, upper), correction)
        result.append(
            (
                ideal_box[coordinate][0] + correction_interval[0],
                ideal_box[coordinate][1] + correction_interval[1],
            )
        )
    return result[0], result[1]


@cache
def _arc_radius_squared(arc: _ArcProof) -> Fraction:
    return _fraction_distance_squared(fraction_xy(arc[0]), fraction_xy(arc[2]))


@cache
def _arc_length_squared(arc: _ArcProof) -> Fraction:
    return _arc_radius_squared(arc) * Fraction.from_float(abs(arc[3])) ** 2


@cache
def _arc_length_interval(arc: _ArcProof) -> tuple[Fraction, Fraction]:
    squared = _arc_length_squared(arc)
    return Fraction.from_float(inward_sqrt_fraction(squared)), Fraction.from_float(outward_sqrt_fraction(squared))


@cache
def _arc_endpoint_correction_upper(arc: _ArcProof) -> Fraction:
    ideal_end = _ideal_arc_point_box(arc, Fraction(1))
    stored_end = fraction_xy(arc[1])
    squared = sum(
        (_interval_component_error(stored_end[index], ideal_end[index]) ** 2 for index in range(2)),
        Fraction(0),
    )
    return Fraction.from_float(outward_sqrt_fraction(squared))


def _interval_component_error(value: Fraction, interval: tuple[Fraction, Fraction]) -> Fraction:
    return max(abs(value - interval[0]), abs(value - interval[1]))


def _interval_product(
    left: tuple[Fraction, Fraction],
    right: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    products = tuple(left_value * right_value for left_value in left for right_value in right)
    return min(products), max(products)


def _outward_float_fraction(value: Fraction) -> float:
    represented = float(value)
    if not math.isfinite(represented):
        raise ArithmeticError("Outward conversion is not finite.")
    while Fraction.from_float(represented) < value:
        represented = math.nextafter(represented, math.inf)
    return represented


def _exact_point_in_convex_polygon(point: _FractionXY, polygon: Sequence[_FractionXY]) -> bool:
    if len(polygon) < 3:
        return False
    signs = [_fraction_cross(_fraction_subtract(end, start), _fraction_subtract(point, start)) for start, end in zip(polygon, (*polygon[1:], polygon[0]))]
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


def _fraction_subtract(left: _FractionXY, right: _FractionXY) -> _FractionXY:
    return left[0] - right[0], left[1] - right[1]


def _fraction_dot(left: _FractionXY, right: _FractionXY) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _fraction_cross(left: _FractionXY, right: _FractionXY) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def _fraction_distance_squared(left: _FractionXY, right: _FractionXY) -> Fraction:
    difference = _fraction_subtract(left, right)
    return _fraction_dot(difference, difference)


def _fraction_midpoint(left: _FractionXY, right: _FractionXY) -> _FractionXY:
    return (left[0] + right[0]) / 2, (left[1] + right[1]) / 2


def _exact_convex_hull_fraction(points: Sequence[_FractionXY]) -> tuple[_FractionXY, ...]:
    ordered = sorted(set(points))
    if len(ordered) <= 1:
        return tuple(ordered)

    def half(sequence: Sequence[_FractionXY]) -> list[_FractionXY]:
        result: list[_FractionXY] = []
        for point in sequence:
            while (
                len(result) >= 2
                and _fraction_cross(
                    _fraction_subtract(result[-1], result[-2]),
                    _fraction_subtract(point, result[-1]),
                )
                <= 0
            ):
                result.pop()
            result.append(point)
        return result

    return tuple(half(ordered)[:-1] + half(tuple(reversed(ordered)))[:-1])


def _exact_segment_distance_squared(
    point: _FractionXY,
    start: _FractionXY,
    end: _FractionXY,
) -> Fraction:
    segment = _fraction_subtract(end, start)
    relative = _fraction_subtract(point, start)
    length_squared = _fraction_dot(segment, segment)
    if length_squared == 0:
        return _fraction_dot(relative, relative)
    projection = _fraction_dot(relative, segment)
    if projection <= 0:
        return _fraction_dot(relative, relative)
    if projection >= length_squared:
        end_relative = _fraction_subtract(point, end)
        return _fraction_dot(end_relative, end_relative)
    return _fraction_cross(segment, relative) ** 2 / length_squared


def _exact_turn(start: _XY, middle: _XY, end: _XY) -> Fraction:
    start_exact = fraction_xy(start)
    return _fraction_cross(
        _fraction_subtract(fraction_xy(middle), start_exact),
        _fraction_subtract(fraction_xy(end), start_exact),
    )


def _finite_normal_or_zero(value: float) -> bool:
    return math.isfinite(value) and (value == 0.0 or abs(value) >= sys.float_info.min)
