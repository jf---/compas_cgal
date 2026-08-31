from __future__ import annotations

import math
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from typing import NewType
from typing import Self
from typing import TypeAlias
from typing import overload

from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import InvalidReferenceProjectionError
from benchmarks.errors import UnresolvedPublishedCurveError
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

PdfPointUnit = NewType("PdfPointUnit", float)
MillimetresPerPdfPoint = NewType("MillimetresPerPdfPoint", float)

# Binary64 roundoff allowance for normalized tangent and reconstructed-radius
# comparisons. These guard arithmetic singularities, not geometric fidelity.
ROUNDOFF_FACTOR = 128.0
NORMALIZED_ROUNDOFF = ROUNDOFF_FACTOR * sys.float_info.epsilon
MAX_RECONSTRUCTION_SUBDIVISIONS = 24

_XY: TypeAlias = tuple[float, float]


def _pdf_coordinate(value: object, name: str) -> PdfPointUnit:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidPublishedPrimitiveError(f"{name} must be a finite real number.")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise InvalidPublishedPrimitiveError(f"{name} must be finite.")
    return PdfPointUnit(numeric)


@dataclass(frozen=True)
class PdfPoint2:
    x: PdfPointUnit
    y: PdfPointUnit

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", _pdf_coordinate(self.x, "PdfPoint2.x"))
        object.__setattr__(self, "y", _pdf_coordinate(self.y, "PdfPoint2.y"))

    @classmethod
    @overload
    def build(cls, x: float, y: float, /) -> Self: ...

    @classmethod
    @overload
    def build(cls, components: Sequence[float], /) -> Self: ...

    @classmethod
    def build(cls, *args: object) -> Self:  # pyright: ignore[reportInconsistentOverload]
        if len(args) == 1 and isinstance(args[0], Sequence) and not isinstance(args[0], (str, bytes)):
            components = args[0]
            if len(components) != 2:
                raise InvalidPublishedPrimitiveError("PdfPoint2 requires exactly two coordinates.")
            return cls(
                _pdf_coordinate(components[0], "PdfPoint2.x"),
                _pdf_coordinate(components[1], "PdfPoint2.y"),
            )
        if len(args) == 2:
            return cls(
                _pdf_coordinate(args[0], "PdfPoint2.x"),
                _pdf_coordinate(args[1], "PdfPoint2.y"),
            )
        raise InvalidPublishedPrimitiveError("PdfPoint2 requires two scalars or one two-coordinate sequence.")


@dataclass(frozen=True)
class SourceLine:
    start: PdfPoint2
    end: PdfPoint2

    def __post_init__(self) -> None:
        if self.start == self.end:
            raise InvalidPublishedPrimitiveError("A published line must have distinct endpoints.")

    @classmethod
    def build(cls, start: PdfPoint2, end: PdfPoint2) -> Self:
        return cls(start, end)


@dataclass(frozen=True)
class SourceCubic:
    start: PdfPoint2
    control1: PdfPoint2
    control2: PdfPoint2
    end: PdfPoint2

    def __post_init__(self) -> None:
        if self.start == self.end:
            raise InvalidPublishedPrimitiveError("A published cubic must have distinct endpoints.")
        if self.start == self.control1 or self.control2 == self.end:
            raise InvalidPublishedPrimitiveError("A published cubic must have non-zero endpoint tangents.")

    @classmethod
    def build(
        cls,
        start: PdfPoint2,
        control1: PdfPoint2,
        control2: PdfPoint2,
        end: PdfPoint2,
    ) -> Self:
        return cls(start, control1, control2, end)


@dataclass(frozen=True)
class SourceToWorld:
    source_origin: PdfPoint2
    world_origin: Point2[WorldXY]
    scale: MillimetresPerPdfPoint

    def __post_init__(self) -> None:
        scale = float(self.scale)
        if not math.isfinite(scale) or scale <= 0.0:
            raise InvalidPublishedPrimitiveError("Source-to-world scale must be finite and positive.")

    @classmethod
    def build(
        cls,
        *,
        source_origin: PdfPoint2,
        world_origin: Point2[WorldXY],
        scale: MillimetresPerPdfPoint,
    ) -> Self:
        return cls(source_origin, world_origin, scale)

    def point(self, source: PdfPoint2) -> Point2[WorldXY]:
        return Point2[WorldXY].build(
            float(self.world_origin.x) + (float(source.x) - float(self.source_origin.x)) * float(self.scale),
            float(self.world_origin.y) + (float(source.y) - float(self.source_origin.y)) * float(self.scale),
        )


@dataclass(frozen=True)
class ReferenceLine:
    start: Point2[WorldXY]
    end: Point2[WorldXY]

    def __post_init__(self) -> None:
        if self.start == self.end:
            raise InvalidPublishedPrimitiveError("A reference line must have distinct endpoints.")

    @classmethod
    def build(cls, start: Point2[WorldXY], end: Point2[WorldXY]) -> Self:
        return cls(start, end)


@dataclass(frozen=True)
class ReferenceArc:
    start: Point2[WorldXY]
    end: Point2[WorldXY]
    centre: Point2[WorldXY]
    sweep: Radian

    def __post_init__(self) -> None:
        sweep = float(self.sweep)
        if not math.isfinite(sweep) or sweep == 0.0 or abs(sweep) > math.tau:
            raise InvalidPublishedPrimitiveError("A reference arc sweep must be finite, non-zero, and at most one turn.")
        start_radius = _subtract(_point_xy(self.start), _point_xy(self.centre))
        end_radius = _subtract(_point_xy(self.end), _point_xy(self.centre))
        if _length(start_radius) == 0.0:
            raise InvalidPublishedPrimitiveError("A reference arc requires distinct endpoints on one circle.")
        predicted_end_radius = _rotate_vector(start_radius, sweep)
        if not _vectors_align(predicted_end_radius, end_radius):
            raise InvalidPublishedPrimitiveError("A reference arc sweep does not terminate at its declared endpoint.")

    @classmethod
    def build(
        cls,
        start: Point2[WorldXY],
        end: Point2[WorldXY],
        centre: Point2[WorldXY],
        sweep: Radian,
    ) -> Self:
        return cls(start, end, centre, sweep)


ReferencePrimitive: TypeAlias = ReferenceLine | ReferenceArc


@dataclass(frozen=True)
class ReferenceBoundary:
    primitives: tuple[ReferencePrimitive, ...]
    tool_radius: ToolRadius
    boundary_stroke_width: Millimetre

    def __post_init__(self) -> None:
        _validate_closed_cycle(self.primitives)
        stroke_width = float(self.boundary_stroke_width)
        if not math.isfinite(stroke_width) or stroke_width <= 0.0:
            raise InvalidPublishedPrimitiveError("Boundary stroke width must be finite and positive.")

    @classmethod
    def build(
        cls,
        primitives: Sequence[ReferencePrimitive],
        tool_radius: ToolRadius,
        boundary_stroke_width: Millimetre,
    ) -> Self:
        return cls(tuple(primitives), tool_radius, boundary_stroke_width)


@dataclass(frozen=True)
class PolygonProjection:
    points: tuple[Point2[WorldXY], ...]
    deviation_limit: Millimetre
    observed_deviation: Millimetre

    def __post_init__(self) -> None:
        limit = float(self.deviation_limit)
        observed = float(self.observed_deviation)
        if not math.isfinite(limit) or limit <= 0.0:
            raise InvalidReferenceProjectionError("Projection deviation limit must be finite and positive.")
        if not math.isfinite(observed) or observed < 0.0 or observed > limit:
            raise InvalidReferenceProjectionError("Observed projection deviation exceeds its declared bound.")
        if len(self.points) < 3 or self.points[0] == self.points[-1]:
            raise InvalidReferenceProjectionError("A projection requires an open-storage ring with at least three points.")

    @classmethod
    def build(
        cls,
        points: Sequence[Point2[WorldXY]],
        deviation_limit: Millimetre,
        observed_deviation: Millimetre,
    ) -> Self:
        return cls(tuple(points), deviation_limit, observed_deviation)


def _validate_closed_cycle(primitives: Sequence[ReferencePrimitive]) -> None:
    if not primitives:
        raise DisconnectedPublishedBoundaryError("A published boundary cannot be empty.")
    for current, following in zip(primitives, (*primitives[1:], primitives[0])):
        if current.end != following.start:
            raise DisconnectedPublishedBoundaryError("Published boundary primitive endpoints are disconnected.")


def project_boundary(
    boundary: ReferenceBoundary,
    deviation_limit: Millimetre,
) -> PolygonProjection:
    limit = float(deviation_limit)
    if not math.isfinite(limit) or limit <= 0.0:
        raise InvalidReferenceProjectionError("Projection deviation limit must be finite and positive.")

    points = [boundary.primitives[0].start]
    maximum_deviation = 0.0
    for primitive in boundary.primitives:
        if isinstance(primitive, ReferenceLine):
            points.append(primitive.end)
            continue

        radius = _distance(_point_xy(primitive.start), _point_xy(primitive.centre))
        if limit > radius:
            raise InvalidReferenceProjectionError("Projection deviation limit cannot exceed an arc radius.")
        maximum_angle = 2.0 * math.acos(1.0 - limit / radius)
        if maximum_angle == 0.0:
            raise InvalidReferenceProjectionError("Projection deviation limit is below binary64 angular resolution.")
        segment_count = math.ceil(abs(float(primitive.sweep)) / maximum_angle)
        segment_angle = abs(float(primitive.sweep)) / segment_count
        sagitta = radius * (1.0 - math.cos(segment_angle / 2.0))
        while sagitta > limit:
            segment_count += 1
            segment_angle = abs(float(primitive.sweep)) / segment_count
            sagitta = radius * (1.0 - math.cos(segment_angle / 2.0))
        centre = _point_xy(primitive.centre)
        start_radius = _subtract(_point_xy(primitive.start), centre)
        for index in range(1, segment_count + 1):
            if index == segment_count:
                emitted_end = primitive.end
            else:
                emitted_end = _world_point(
                    _add(
                        centre,
                        _rotate_vector(
                            start_radius,
                            float(primitive.sweep) * index / segment_count,
                        ),
                    )
                )
            actual_start_radius = _subtract(_point_xy(points[-1]), centre)
            actual_end_radius = _subtract(_point_xy(emitted_end), centre)
            ideal_end_radius = _rotate_vector(
                start_radius,
                float(primitive.sweep) * index / segment_count,
            )
            ideal_midpoint_radius = _rotate_vector(
                start_radius,
                float(primitive.sweep) * (index - 0.5) / segment_count,
            )
            emitted_coordinate_error = _distance(ideal_end_radius, actual_end_radius)
            emitted_chord_deviation = _point_segment_distance(
                ideal_midpoint_radius,
                actual_start_radius,
                actual_end_radius,
            )
            maximum_deviation = max(
                maximum_deviation,
                emitted_coordinate_error,
                emitted_chord_deviation,
            )
            if maximum_deviation > limit:
                raise InvalidReferenceProjectionError("Emitted polygon coordinates exceed the projection deviation limit.")
            points.append(emitted_end)

    if points[-1] != points[0]:
        raise InvalidReferenceProjectionError("Projected primitives did not preserve the analytic cycle closure.")
    points.pop()
    return PolygonProjection.build(
        points,
        deviation_limit,
        Millimetre(maximum_deviation),
    )


def reconstruct_cubic(
    source: SourceCubic,
    transform: SourceToWorld,
    deviation_limit: Millimetre,
) -> tuple[ReferenceArc, ...]:
    limit = float(deviation_limit)
    if not math.isfinite(limit) or limit <= 0.0:
        raise InvalidPublishedPrimitiveError("Published-curve deviation limit must be finite and positive.")
    control_points = (
        _point_xy(transform.point(source.start)),
        _point_xy(transform.point(source.control1)),
        _point_xy(transform.point(source.control2)),
        _point_xy(transform.point(source.end)),
    )
    return _reconstruct_control_points(control_points, limit, depth=0)


def _reconstruct_control_points(
    control_points: tuple[_XY, _XY, _XY, _XY],
    deviation_limit: float,
    *,
    depth: int,
) -> tuple[ReferenceArc, ...]:
    circle = _single_circle_candidate(control_points)
    if (
        circle is not None
        and _cubic_has_monotone_polar_angle(
            control_points,
            circle[0],
            1.0 if circle[1] > 0.0 else -1.0,
            depth=0,
        )
        and _cubic_within_circle(
            control_points,
            circle,
            deviation_limit,
            depth=0,
        )
    ):
        centre, sweep = circle
        return (
            ReferenceArc.build(
                _world_point(control_points[0]),
                _world_point(control_points[3]),
                _world_point(centre),
                Radian(sweep),
            ),
        )
    biarc = _equal_distance_biarc(control_points)
    if biarc is not None and _cubic_within_biarc(
        control_points,
        biarc,
        deviation_limit,
        start_parameter=0.0,
        end_parameter=1.0,
        depth=0,
    ):
        return biarc
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        raise UnresolvedPublishedCurveError("Published cubic did not admit a bounded circular reconstruction.")
    if _control_points_are_collinear(control_points):
        raise UnresolvedPublishedCurveError("A collinear published cubic must be represented as a source line.")
    left, right = _split_cubic(control_points)
    return (
        *_reconstruct_control_points(left, deviation_limit, depth=depth + 1),
        *_reconstruct_control_points(right, deviation_limit, depth=depth + 1),
    )


def _single_circle_candidate(
    control_points: tuple[_XY, _XY, _XY, _XY],
) -> tuple[_XY, float] | None:
    start, control1, control2, end = control_points
    start_tangent = _unit(_subtract(control1, start))
    end_tangent = _unit(_subtract(end, control2))
    chord = _subtract(end, start)
    start_normal = (-start_tangent[1], start_tangent[0])
    divisor = 2.0 * _dot(chord, start_normal)
    if abs(divisor) <= NORMALIZED_ROUNDOFF * _length(chord):
        return None

    centre = _add(start, _scale(start_normal, _dot(chord, chord) / divisor))

    start_radius = _subtract(start, centre)
    end_radius = _subtract(end, centre)
    end_radius_unit = _unit(end_radius)
    if abs(_dot(end_tangent, end_radius_unit)) > NORMALIZED_ROUNDOFF:
        return None
    ccw_start_tangent = (-start_radius[1], start_radius[0])
    ccw_end_tangent = (-end_radius[1], end_radius[0])
    start_direction = _dot(start_tangent, ccw_start_tangent)
    end_direction = _dot(end_tangent, ccw_end_tangent)
    if start_direction * end_direction <= 0.0:
        return None
    direction = 1.0 if start_direction > 0.0 else -1.0
    sweep = _sweep_between(start, end, centre, direction)
    if sweep == 0.0:
        return None
    return centre, sweep


def _equal_distance_biarc(
    control_points: tuple[_XY, _XY, _XY, _XY],
) -> tuple[ReferenceArc, ReferenceArc] | None:
    start, control1, control2, end = control_points
    start_tangent = _unit(_subtract(control1, start))
    end_tangent = _unit(_subtract(end, control2))
    chord = _subtract(end, start)
    tangent_dot = _dot(start_tangent, end_tangent)
    denominator = 2.0 * (1.0 - tangent_dot)
    tangent_sum = _add(start_tangent, end_tangent)
    chord_dot_tangents = _dot(chord, tangent_sum)
    chord_squared = _dot(chord, chord)

    if denominator <= NORMALIZED_ROUNDOFF:
        chord_dot_end_tangent = _dot(chord, end_tangent)
        if abs(chord_dot_end_tangent) <= NORMALIZED_ROUNDOFF * _length(chord) and _tangents_align(start_tangent, end_tangent):
            return _opposed_semicircle_biarc(start, end, end_tangent)
    if denominator == 0.0:
        return None
    discriminant = chord_dot_tangents**2 + denominator * chord_squared
    root = math.sqrt(discriminant)
    stable_divisor = root + chord_dot_tangents
    if stable_divisor > NORMALIZED_ROUNDOFF * root:
        distance = chord_squared / stable_divisor
    else:
        distance = (-chord_dot_tangents + root) / denominator
    if not math.isfinite(distance) or distance <= 0.0:
        return None

    join = _scale(
        _add(_add(start, end), _scale(_subtract(start_tangent, end_tangent), distance)),
        0.5,
    )
    first = _arc_from_start_tangent(start, join, start_tangent)
    second = _arc_from_end_tangent(join, end, end_tangent)
    if first is None or second is None:
        return None
    if not (
        _tangents_align(_arc_tangent(first, at_end=False), start_tangent)
        and _tangents_align(_arc_tangent(first, at_end=True), _arc_tangent(second, at_end=False))
        and _tangents_align(_arc_tangent(second, at_end=True), end_tangent)
    ):
        return None
    return first, second


def _opposed_semicircle_biarc(start: _XY, end: _XY, tangent: _XY) -> tuple[ReferenceArc, ReferenceArc] | None:
    chord = _subtract(end, start)
    if _cross(chord, tangent) == 0.0:
        return None
    join = _midpoint(start, end)
    first_centre = _add(start, _scale(chord, 0.25))
    second_centre = _add(start, _scale(chord, 0.75))
    first_sweep = math.pi if _cross(chord, tangent) < 0.0 else -math.pi
    second_sweep = -first_sweep
    return (
        ReferenceArc.build(
            _world_point(start),
            _world_point(join),
            _world_point(first_centre),
            Radian(first_sweep),
        ),
        ReferenceArc.build(
            _world_point(join),
            _world_point(end),
            _world_point(second_centre),
            Radian(second_sweep),
        ),
    )


def _arc_from_start_tangent(start: _XY, end: _XY, start_tangent: _XY) -> ReferenceArc | None:
    chord = _subtract(end, start)
    normal = (-start_tangent[1], start_tangent[0])
    divisor = 2.0 * _dot(normal, chord)
    if abs(divisor) <= NORMALIZED_ROUNDOFF * _length(chord):
        return None
    centre = _add(start, _scale(normal, _dot(chord, chord) / divisor))
    sweep = _directed_sweep(start, end, centre, start_tangent)
    return ReferenceArc.build(_world_point(start), _world_point(end), _world_point(centre), Radian(sweep))


def _arc_from_end_tangent(start: _XY, end: _XY, end_tangent: _XY) -> ReferenceArc | None:
    chord = _subtract(start, end)
    normal = (-end_tangent[1], end_tangent[0])
    divisor = 2.0 * _dot(normal, chord)
    if abs(divisor) <= NORMALIZED_ROUNDOFF * _length(chord):
        return None
    centre = _add(end, _scale(normal, _dot(chord, chord) / divisor))
    end_radius = _subtract(end, centre)
    direction = 1.0 if _dot(end_tangent, (-end_radius[1], end_radius[0])) > 0.0 else -1.0
    sweep = _sweep_between(start, end, centre, direction)
    return ReferenceArc.build(_world_point(start), _world_point(end), _world_point(centre), Radian(sweep))


def _directed_sweep(start: _XY, end: _XY, centre: _XY, tangent: _XY) -> float:
    start_radius = _subtract(start, centre)
    direction = 1.0 if _dot(tangent, (-start_radius[1], start_radius[0])) > 0.0 else -1.0
    return _sweep_between(start, end, centre, direction)


def _sweep_between(start: _XY, end: _XY, centre: _XY, direction: float) -> float:
    start_radius = _subtract(start, centre)
    end_radius = _subtract(end, centre)
    start_angle = math.atan2(start_radius[1], start_radius[0])
    end_angle = math.atan2(end_radius[1], end_radius[0])
    if direction > 0.0:
        return (end_angle - start_angle) % math.tau
    return -((start_angle - end_angle) % math.tau)


def _cubic_within_biarc(
    control_points: tuple[_XY, _XY, _XY, _XY],
    biarc: tuple[ReferenceArc, ReferenceArc],
    deviation_limit: float,
    *,
    start_parameter: float,
    end_parameter: float,
    depth: int,
) -> bool:
    midpoint_parameter = (start_parameter + end_parameter) / 2.0
    cubic_point = _evaluate_cubic(control_points, midpoint_parameter)
    biarc_point = _evaluate_biarc(biarc, midpoint_parameter)
    observed = _distance(cubic_point, biarc_point)
    if observed > deviation_limit:
        return False

    half_width = (end_parameter - start_parameter) / 2.0
    cubic_speed_bound = 3.0 * max(_distance(control_points[index], control_points[index + 1]) for index in range(3))
    biarc_speed = sum(_arc_length(arc) for arc in biarc)
    if observed + half_width * (cubic_speed_bound + biarc_speed) <= deviation_limit:
        return True
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        return False
    return _cubic_within_biarc(
        control_points,
        biarc,
        deviation_limit,
        start_parameter=start_parameter,
        end_parameter=midpoint_parameter,
        depth=depth + 1,
    ) and _cubic_within_biarc(
        control_points,
        biarc,
        deviation_limit,
        start_parameter=midpoint_parameter,
        end_parameter=end_parameter,
        depth=depth + 1,
    )


def _evaluate_biarc(biarc: tuple[ReferenceArc, ReferenceArc], parameter: float) -> _XY:
    first_length = _arc_length(biarc[0])
    second_length = _arc_length(biarc[1])
    total_length = first_length + second_length
    travelled = parameter * total_length
    if travelled <= first_length:
        return _evaluate_arc(biarc[0], travelled / first_length)
    return _evaluate_arc(biarc[1], (travelled - first_length) / second_length)


def _evaluate_arc(arc: ReferenceArc, parameter: float) -> _XY:
    return _rotate_about(
        _point_xy(arc.start),
        _point_xy(arc.centre),
        float(arc.sweep) * parameter,
    )


def _arc_length(arc: ReferenceArc) -> float:
    return _distance(_point_xy(arc.start), _point_xy(arc.centre)) * abs(float(arc.sweep))


def _arc_tangent(arc: ReferenceArc, *, at_end: bool) -> _XY:
    point = arc.end if at_end else arc.start
    radius = _subtract(_point_xy(point), _point_xy(arc.centre))
    direction = 1.0 if float(arc.sweep) > 0.0 else -1.0
    return _unit(_scale((-radius[1], radius[0]), direction))


def _tangents_align(left: _XY, right: _XY) -> bool:
    return _distance(left, right) <= NORMALIZED_ROUNDOFF


def _cubic_within_circle(
    control_points: tuple[_XY, _XY, _XY, _XY],
    candidate: tuple[_XY, float],
    deviation_limit: float,
    *,
    depth: int,
) -> bool:
    centre, _ = candidate
    radius = _distance(control_points[0], centre)
    midpoint = _evaluate_cubic(control_points, 0.5)
    if abs(_distance(midpoint, centre) - radius) > deviation_limit:
        return False

    minimum_radius, maximum_radius = _control_hull_radius_bounds(control_points, centre)
    if max(radius - minimum_radius, maximum_radius - radius) <= deviation_limit:
        return True
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        return False
    left, right = _split_cubic(control_points)
    return _cubic_within_circle(left, candidate, deviation_limit, depth=depth + 1) and _cubic_within_circle(right, candidate, deviation_limit, depth=depth + 1)


def _cubic_has_monotone_polar_angle(
    control_points: tuple[_XY, _XY, _XY, _XY],
    centre: _XY,
    direction: float,
    *,
    depth: int,
) -> bool:
    coefficients = _polar_derivative_bernstein_coefficients(control_points, centre)
    directed = tuple(direction * coefficient for coefficient in coefficients)
    if min(directed) >= 0.0:
        return True

    midpoint = _evaluate_cubic(control_points, 0.5)
    midpoint_derivative = _evaluate_cubic_derivative(control_points, 0.5)
    if direction * _cross(_subtract(midpoint, centre), midpoint_derivative) < 0.0:
        return False
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        return False
    left, right = _split_cubic(control_points)
    return _cubic_has_monotone_polar_angle(left, centre, direction, depth=depth + 1) and _cubic_has_monotone_polar_angle(right, centre, direction, depth=depth + 1)


def _polar_derivative_bernstein_coefficients(control_points: tuple[_XY, _XY, _XY, _XY], centre: _XY) -> tuple[float, float, float, float, float, float]:
    relative = tuple(_subtract(point, centre) for point in control_points)
    derivatives = tuple(_scale(_subtract(control_points[index + 1], control_points[index]), 3.0) for index in range(3))
    coefficients: list[float] = []
    for degree in range(6):
        coefficient = 0.0
        for cubic_index in range(max(0, degree - 2), min(3, degree) + 1):
            derivative_index = degree - cubic_index
            weight = math.comb(3, cubic_index) * math.comb(2, derivative_index) / math.comb(5, degree)
            coefficient += weight * _cross(relative[cubic_index], derivatives[derivative_index])
        coefficients.append(coefficient)
    return (
        coefficients[0],
        coefficients[1],
        coefficients[2],
        coefficients[3],
        coefficients[4],
        coefficients[5],
    )


def _evaluate_cubic_derivative(control_points: tuple[_XY, _XY, _XY, _XY], parameter: float) -> _XY:
    complement = 1.0 - parameter
    edges = tuple(_subtract(control_points[index + 1], control_points[index]) for index in range(3))
    return _scale(
        _add(
            _add(
                _scale(edges[0], complement**2),
                _scale(edges[1], 2.0 * complement * parameter),
            ),
            _scale(edges[2], parameter**2),
        ),
        3.0,
    )


def _control_hull_radius_bounds(control_points: tuple[_XY, _XY, _XY, _XY], centre: _XY) -> tuple[float, float]:
    hull = _convex_hull(control_points)
    maximum = max(_distance(point, centre) for point in hull)
    if _point_in_convex_polygon(centre, hull):
        return 0.0, maximum
    minimum = min(_point_segment_distance(centre, start, end) for start, end in zip(hull, (*hull[1:], hull[0])))
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


def _point_segment_distance(point: _XY, start: _XY, end: _XY) -> float:
    segment = _subtract(end, start)
    length_squared = _dot(segment, segment)
    if length_squared == 0.0:
        return _distance(point, start)
    parameter = max(0.0, min(1.0, _dot(_subtract(point, start), segment) / length_squared))
    return _distance(point, _add(start, _scale(segment, parameter)))


def _split_cubic(
    control_points: tuple[_XY, _XY, _XY, _XY],
) -> tuple[tuple[_XY, _XY, _XY, _XY], tuple[_XY, _XY, _XY, _XY]]:
    start, control1, control2, end = control_points
    first = _midpoint(start, control1)
    second = _midpoint(control1, control2)
    third = _midpoint(control2, end)
    fourth = _midpoint(first, second)
    fifth = _midpoint(second, third)
    midpoint = _midpoint(fourth, fifth)
    return (start, first, fourth, midpoint), (midpoint, fifth, third, end)


def _evaluate_cubic(control_points: tuple[_XY, _XY, _XY, _XY], parameter: float) -> _XY:
    complement = 1.0 - parameter
    start, control1, control2, end = control_points
    return (
        complement**3 * start[0] + 3.0 * complement**2 * parameter * control1[0] + 3.0 * complement * parameter**2 * control2[0] + parameter**3 * end[0],
        complement**3 * start[1] + 3.0 * complement**2 * parameter * control1[1] + 3.0 * complement * parameter**2 * control2[1] + parameter**3 * end[1],
    )


def _control_points_are_collinear(
    control_points: tuple[_XY, _XY, _XY, _XY],
) -> bool:
    chord = _subtract(control_points[3], control_points[0])
    chord_length = _length(chord)
    return all(
        abs(_cross(chord, _subtract(point, control_points[0]))) <= NORMALIZED_ROUNDOFF * chord_length * _distance(point, control_points[0]) for point in control_points[1:3]
    )


def _vectors_align(left: _XY, right: _XY) -> bool:
    local_scale = max(_length(left), _length(right))
    return _distance(left, right) <= NORMALIZED_ROUNDOFF * local_scale


def _rotate_vector(vector: _XY, angle: float) -> _XY:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return cosine * vector[0] - sine * vector[1], sine * vector[0] + cosine * vector[1]


def _rotate_about(point: _XY, centre: _XY, angle: float) -> _XY:
    radius = _subtract(point, centre)
    rotated = _rotate_vector(radius, angle)
    return centre[0] + rotated[0], centre[1] + rotated[1]


def _point_xy(point: Point2[WorldXY]) -> _XY:
    return float(point.x), float(point.y)


def _world_point(point: _XY) -> Point2[WorldXY]:
    return Point2[WorldXY].build(*point)


def _midpoint(left: _XY, right: _XY) -> _XY:
    return (left[0] + right[0]) / 2.0, (left[1] + right[1]) / 2.0


def _add(left: _XY, right: _XY) -> _XY:
    return left[0] + right[0], left[1] + right[1]


def _subtract(left: _XY, right: _XY) -> _XY:
    return left[0] - right[0], left[1] - right[1]


def _scale(vector: _XY, factor: float) -> _XY:
    return vector[0] * factor, vector[1] * factor


def _dot(left: _XY, right: _XY) -> float:
    return left[0] * right[0] + left[1] * right[1]


def _cross(left: _XY, right: _XY) -> float:
    return left[0] * right[1] - left[1] * right[0]


def _length(vector: _XY) -> float:
    return math.hypot(*vector)


def _unit(vector: _XY) -> _XY:
    length = _length(vector)
    if length == 0.0:
        raise InvalidPublishedPrimitiveError("A published tangent must be non-zero.")
    return vector[0] / length, vector[1] / length


def _distance(left: _XY, right: _XY) -> float:
    return _length(_subtract(left, right))
