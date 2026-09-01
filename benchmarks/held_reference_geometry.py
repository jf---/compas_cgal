from __future__ import annotations

import math
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import NewType
from typing import Self
from typing import TypeAlias
from typing import overload

from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import InvalidReferenceProjectionError
from benchmarks.errors import InvalidReferenceReconstructionError
from benchmarks.errors import UnresolvedPublishedCurveError
from benchmarks.held_reference_certification import affine_join_error
from benchmarks.held_reference_certification import biarc_correspondence_node_bound
from benchmarks.held_reference_certification import biarc_correspondence_sample_exceeds
from benchmarks.held_reference_certification import biarc_solver_coefficients
from benchmarks.held_reference_certification import certify_biarc_root
from benchmarks.held_reference_certification import certify_cubic_interval_circle
from benchmarks.held_reference_certification import certify_monotone_polar_angle
from benchmarks.held_reference_certification import circle_deviation_upper_bound
from benchmarks.held_reference_certification import exact_reflection
from benchmarks.held_reference_certification import exact_segment_distance_squared
from benchmarks.held_reference_certification import exact_vector_error
from benchmarks.held_reference_certification import fraction_xy
from benchmarks.held_reference_certification import mapped_radius_tangent_bound
from benchmarks.held_reference_certification import maximum_arc_segment_distance
from benchmarks.held_reference_certification import outward_sqrt_fraction
from benchmarks.held_reference_certification import outward_sum
from benchmarks.held_reference_certification import reflection_perturbation_bound
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
_CubicWitness: TypeAlias = tuple[tuple[_XY, _XY, _XY, _XY], Fraction, Fraction]


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
    reflect_source_y: bool = False

    def __post_init__(self) -> None:
        scale = float(self.scale)
        if not math.isfinite(scale) or scale <= 0.0:
            raise InvalidPublishedPrimitiveError("Source-to-world scale must be finite and positive.")
        if not isinstance(self.reflect_source_y, bool):
            raise InvalidPublishedPrimitiveError("Source-to-world Y reflection must be a bool.")

    @classmethod
    def build(
        cls,
        *,
        source_origin: PdfPoint2,
        world_origin: Point2[WorldXY],
        scale: MillimetresPerPdfPoint,
        reflect_source_y: bool = False,
    ) -> Self:
        return cls(source_origin, world_origin, scale, reflect_source_y)

    def point(self, source: PdfPoint2) -> Point2[WorldXY]:
        y_direction = -1.0 if self.reflect_source_y else 1.0
        return Point2[WorldXY].build(
            float(self.world_origin.x) + (float(source.x) - float(self.source_origin.x)) * float(self.scale),
            float(self.world_origin.y) + y_direction * (float(source.y) - float(self.source_origin.y)) * float(self.scale),
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
SourcePrimitive: TypeAlias = SourceLine | SourceCubic


@dataclass(frozen=True)
class ReferenceReconstruction:
    primitives: tuple[ReferencePrimitive, ...]
    deviation_upper_bound: Millimetre

    def __post_init__(self) -> None:
        if not self.primitives:
            raise InvalidReferenceReconstructionError("A reference reconstruction cannot be empty.")
        for current, following in zip(self.primitives, self.primitives[1:]):
            if current.end != following.start:
                raise InvalidReferenceReconstructionError("Reference reconstruction primitives are disconnected.")
        upper_bound = float(self.deviation_upper_bound)
        if not math.isfinite(upper_bound) or upper_bound < 0.0:
            raise InvalidReferenceReconstructionError("Reference reconstruction bound must be finite and non-negative.")

    @classmethod
    def build(
        cls,
        primitives: Sequence[ReferencePrimitive],
        deviation_upper_bound: Millimetre,
    ) -> Self:
        return cls(tuple(primitives), deviation_upper_bound)


@dataclass(frozen=True)
class _CertifiedSpan:
    control_points: tuple[_XY, _XY, _XY, _XY]
    primitives: tuple[ReferencePrimitive, ...]
    deviation_upper_bound: float


@dataclass(frozen=True)
class _MergeEntry:
    primitive: ReferencePrimitive
    source_witnesses: tuple[_CubicWitness, ...] | None
    deviation_upper_bound: float


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
            emitted_chord_deviation = maximum_arc_segment_distance(
                start_radius,
                float(primitive.sweep) * (index - 1) / segment_count,
                float(primitive.sweep) * index / segment_count,
                actual_start_radius,
                actual_end_radius,
            )
            maximum_deviation = max(
                maximum_deviation,
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
) -> tuple[ReferencePrimitive, ...]:
    return reconstruct_cubic_certified(source, transform, deviation_limit).primitives


def reconstruct_cubic_certified(
    source: SourceCubic,
    transform: SourceToWorld,
    deviation_limit: Millimetre,
) -> ReferenceReconstruction:
    limit = float(deviation_limit)
    if not math.isfinite(limit) or limit <= 0.0:
        raise InvalidPublishedPrimitiveError("Published-curve deviation limit must be finite and positive.")
    control_points = (
        _point_xy(transform.point(source.start)),
        _point_xy(transform.point(source.control1)),
        _point_xy(transform.point(source.control2)),
        _point_xy(transform.point(source.end)),
    )
    spans = _reconstruct_control_points(control_points, limit, depth=0)
    primitives = tuple(primitive for span in spans for primitive in span.primitives)
    return ReferenceReconstruction.build(
        primitives,
        Millimetre(max(span.deviation_upper_bound for span in spans)),
    )


def reconstruct_source_path(
    sources: Sequence[SourcePrimitive],
    transform: SourceToWorld,
    deviation_limit: Millimetre,
) -> ReferenceReconstruction:
    limit = float(deviation_limit)
    if not math.isfinite(limit) or limit <= 0.0:
        raise InvalidPublishedPrimitiveError("Published-path deviation limit must be finite and positive.")
    entries: list[_MergeEntry] = []
    for source in sources:
        if isinstance(source, SourceLine):
            primitive = ReferenceLine.build(transform.point(source.start), transform.point(source.end))
            entries.append(_MergeEntry(primitive, None, 0.0))
            continue
        control_points = (
            _point_xy(transform.point(source.start)),
            _point_xy(transform.point(source.control1)),
            _point_xy(transform.point(source.control2)),
            _point_xy(transform.point(source.end)),
        )
        for span in _reconstruct_control_points(control_points, limit, depth=0):
            if len(span.primitives) == 1 and isinstance(span.primitives[0], ReferenceArc):
                entries.append(
                    _MergeEntry(
                        span.primitives[0],
                        ((span.control_points, Fraction(0), Fraction(1)),),
                        span.deviation_upper_bound,
                    )
                )
            elif len(span.primitives) == 2 and all(isinstance(primitive, ReferenceArc) for primitive in span.primitives):
                first_arc = span.primitives[0]
                second_arc = span.primitives[1]
                assert isinstance(first_arc, ReferenceArc) and isinstance(second_arc, ReferenceArc)
                first_length = _distance(
                    _point_xy(first_arc.start),
                    _point_xy(first_arc.centre),
                ) * abs(float(first_arc.sweep))
                second_length = _distance(
                    _point_xy(second_arc.start),
                    _point_xy(second_arc.centre),
                ) * abs(float(second_arc.sweep))
                breakpoint = Fraction.from_float(first_length / (first_length + second_length))
                intervals = ((Fraction(0), breakpoint), (breakpoint, Fraction(1)))
                entries.extend(
                    _MergeEntry(
                        primitive,
                        ((span.control_points, *interval),),
                        span.deviation_upper_bound,
                    )
                    for primitive, interval in zip(span.primitives, intervals)
                )
            else:
                entries.extend(_MergeEntry(primitive, None, span.deviation_upper_bound) for primitive in span.primitives)
    if not entries:
        raise InvalidReferenceReconstructionError("A source path cannot be empty.")

    merged: list[_MergeEntry] = []
    for entry in entries:
        candidate = _merge_arc_entries(merged[-1], entry, limit) if merged else None
        if candidate is None:
            merged.append(entry)
        else:
            merged[-1] = candidate
    return ReferenceReconstruction.build(
        tuple(entry.primitive for entry in merged),
        Millimetre(max(entry.deviation_upper_bound for entry in merged)),
    )


def _merge_arc_entries(
    first: _MergeEntry,
    second: _MergeEntry,
    deviation_limit: float,
) -> _MergeEntry | None:
    if not isinstance(first.primitive, ReferenceArc) or not isinstance(second.primitive, ReferenceArc):
        return None
    if first.source_witnesses is None or second.source_witnesses is None:
        return None
    if first.primitive.end != second.primitive.start:
        return None
    first_sweep = float(first.primitive.sweep)
    second_sweep = float(second.primitive.sweep)
    if first_sweep * second_sweep <= 0.0:
        return None
    combined_sweep = first_sweep + second_sweep
    if abs(combined_sweep) > math.tau:
        return None
    try:
        combined = ReferenceArc.build(
            first.primitive.start,
            second.primitive.end,
            first.primitive.centre,
            Radian(combined_sweep),
        )
    except InvalidPublishedPrimitiveError:
        return None

    centre = _point_xy(combined.centre)
    bounds: list[float] = []
    witnesses = (*first.source_witnesses, *second.source_witnesses)
    for control_points, start_parameter, end_parameter in witnesses:
        bound = certify_cubic_interval_circle(
            control_points,
            start_parameter,
            end_parameter,
            centre,
            _point_xy(combined.start),
            1.0 if combined_sweep > 0.0 else -1.0,
            deviation_limit,
            MAX_RECONSTRUCTION_SUBDIVISIONS,
        )
        if bound is None:
            return None
        bounds.append(bound)
    return _MergeEntry(combined, witnesses, max(bounds))


def _reconstruct_control_points(
    control_points: tuple[_XY, _XY, _XY, _XY],
    deviation_limit: float,
    *,
    depth: int,
) -> tuple[_CertifiedSpan, ...]:
    circle = _single_circle_candidate(control_points)
    circle_bound = None
    if circle is not None and _cubic_has_monotone_polar_angle(
        control_points,
        circle[0],
        1.0 if circle[1] > 0.0 else -1.0,
        depth=0,
    ):
        circle_bound = _cubic_within_circle(
            control_points,
            circle,
            deviation_limit,
            depth=0,
        )
    if circle is not None and circle_bound is not None:
        centre, sweep = circle
        return (
            _CertifiedSpan(
                control_points,
                (
                    ReferenceArc.build(
                        _world_point(control_points[0]),
                        _world_point(control_points[3]),
                        _world_point(centre),
                        Radian(sweep),
                    ),
                ),
                circle_bound,
            ),
        )
    biarc = _equal_distance_biarc(control_points)
    if biarc is not None:
        biarc_bound = _cubic_within_biarc(
            control_points,
            biarc,
            deviation_limit,
            start_parameter=Fraction(0),
            end_parameter=Fraction(1),
            depth=0,
        )
        if biarc_bound is not None:
            return (_CertifiedSpan(control_points, biarc, biarc_bound),)
    if _control_points_are_exactly_collinear(control_points):
        line = _certified_line(control_points, deviation_limit)
        if line is not None:
            primitive, line_bound = line
            return (_CertifiedSpan(control_points, (primitive,), line_bound),)
        raise UnresolvedPublishedCurveError("A collinear published cubic must advance along its authored chord.")
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        raise UnresolvedPublishedCurveError("Published cubic did not admit a bounded circular reconstruction.")
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
    chord_length = _length(chord)
    if not math.isfinite(chord_length) or chord_length < sys.float_info.min:
        return None
    local_end = _scale(chord, 1.0 / chord_length)
    solver_coefficients = biarc_solver_coefficients(start_tangent, end_tangent, local_end)
    if solver_coefficients is None:
        return None
    denominator, chord_dot_tangents, chord_squared = solver_coefficients

    if denominator <= NORMALIZED_ROUNDOFF:
        chord_dot_end_tangent = _dot(local_end, end_tangent)
        if abs(chord_dot_end_tangent) <= NORMALIZED_ROUNDOFF and _tangents_align(start_tangent, end_tangent):
            return _opposed_semicircle_biarc(start, end, end_tangent)
    if denominator == 0.0:
        return None
    discriminant = chord_dot_tangents**2 + denominator * chord_squared
    if not _finite_normal_or_zero(discriminant) or discriminant < 0.0:
        return None
    root = math.sqrt(discriminant)
    if not _finite_normal_or_zero(root):
        return None
    stable_divisor = root + chord_dot_tangents
    if not _finite_normal_or_zero(stable_divisor):
        return None
    if stable_divisor > NORMALIZED_ROUNDOFF * root:
        distance = chord_squared / stable_divisor
    else:
        distance = (-chord_dot_tangents + root) / denominator
    if not _finite_normal_or_zero(distance) or distance <= 0.0:
        return None
    root_uncertainty = certify_biarc_root(
        start_tangent,
        end_tangent,
        local_end,
        distance,
    )
    if root_uncertainty is None:
        return None

    local_start = (0.0, 0.0)
    local_join = _scale(
        _add(local_end, _scale(_subtract(start_tangent, end_tangent), distance)),
        0.5,
    )
    first_local = _arc_from_start_tangent(local_start, local_join, start_tangent)
    second_local = _arc_from_end_tangent(local_join, local_end, end_tangent)
    if first_local is None or second_local is None:
        return None
    local_tangent_bounds = _biarc_tangent_bounds(
        local_end,
        local_join,
        start_tangent,
        end_tangent,
        distance,
        root_uncertainty,
        first_local,
        second_local,
    )
    if local_tangent_bounds is None:
        return None
    first, second = _map_local_biarc(
        first_local,
        second_local,
        start,
        end,
        chord_length,
    )
    storage_bounds = _mapped_biarc_tangent_bounds(
        (first_local, second_local),
        (first, second),
        chord_length,
    )
    if storage_bounds is None:
        return None
    if not (
        _distance(_arc_tangent(first, at_end=False), start_tangent) <= local_tangent_bounds[0] + storage_bounds[0]
        and _distance(_arc_tangent(first, at_end=True), _arc_tangent(second, at_end=False)) <= local_tangent_bounds[1] + storage_bounds[1]
        and _distance(_arc_tangent(second, at_end=True), end_tangent) <= local_tangent_bounds[2] + storage_bounds[2]
    ):
        return None
    return first, second


def _biarc_tangent_bounds(
    local_end: _XY,
    local_join: _XY,
    start_tangent: _XY,
    end_tangent: _XY,
    distance: float,
    root_uncertainty: float,
    first_arc: ReferenceArc,
    second_arc: ReferenceArc,
) -> tuple[float, float, float] | None:
    tangent_difference = _subtract(start_tangent, end_tangent)
    join_error = affine_join_error(
        local_join,
        local_end,
        tangent_difference,
        distance,
        root_uncertainty,
    )
    if join_error is None:
        return None

    first_chord = local_join
    second_chord = _subtract(local_end, local_join)
    first_geometric = reflection_perturbation_bound(join_error, first_chord)
    second_geometric = reflection_perturbation_bound(join_error, second_chord)
    if first_geometric is None or second_geometric is None:
        return None
    first_exact_join = exact_reflection(start_tangent, first_chord)
    second_exact_join = exact_reflection(end_tangent, second_chord)
    start_evaluation = exact_vector_error(_arc_tangent(first_arc, at_end=False), fraction_xy(start_tangent))
    first_join_evaluation = exact_vector_error(_arc_tangent(first_arc, at_end=True), first_exact_join)
    second_join_evaluation = exact_vector_error(_arc_tangent(second_arc, at_end=False), second_exact_join)
    end_evaluation = exact_vector_error(_arc_tangent(second_arc, at_end=True), fraction_xy(end_tangent))
    return (
        start_evaluation,
        outward_sum(
            outward_sum(first_geometric, second_geometric),
            outward_sum(first_join_evaluation, second_join_evaluation),
        ),
        end_evaluation,
    )


def _map_local_point(point: Point2[WorldXY], origin: _XY, scale: float) -> Point2[WorldXY]:
    return _world_point(_add(origin, _scale(_point_xy(point), scale)))


def _map_local_biarc(
    first: ReferenceArc,
    second: ReferenceArc,
    start: _XY,
    end: _XY,
    scale: float,
) -> tuple[ReferenceArc, ReferenceArc]:
    join = _map_local_point(first.end, start, scale)
    return (
        ReferenceArc.build(
            _world_point(start),
            join,
            _map_local_point(first.centre, start, scale),
            first.sweep,
        ),
        ReferenceArc.build(
            join,
            _world_point(end),
            _map_local_point(second.centre, start, scale),
            second.sweep,
        ),
    )


def _mapped_arc_tangent_bound(
    local_arc: ReferenceArc,
    mapped_arc: ReferenceArc,
    scale: float,
    *,
    at_end: bool,
) -> float | None:
    local_point = local_arc.end if at_end else local_arc.start
    mapped_point = mapped_arc.end if at_end else mapped_arc.start
    return mapped_radius_tangent_bound(
        _point_xy(local_point),
        _point_xy(local_arc.centre),
        _point_xy(mapped_point),
        _point_xy(mapped_arc.centre),
        scale,
    )


def _mapped_biarc_tangent_bounds(
    local_biarc: tuple[ReferenceArc, ReferenceArc],
    mapped_biarc: tuple[ReferenceArc, ReferenceArc],
    scale: float,
) -> tuple[float, float, float] | None:
    bounds = (
        _mapped_arc_tangent_bound(local_biarc[0], mapped_biarc[0], scale, at_end=False),
        _mapped_arc_tangent_bound(local_biarc[0], mapped_biarc[0], scale, at_end=True),
        _mapped_arc_tangent_bound(local_biarc[1], mapped_biarc[1], scale, at_end=False),
        _mapped_arc_tangent_bound(local_biarc[1], mapped_biarc[1], scale, at_end=True),
    )
    if any(bound is None for bound in bounds):
        return None
    finite = tuple(bound for bound in bounds if bound is not None)
    combined = finite[0], outward_sum(finite[1], finite[2]), finite[3]
    if any(not math.isfinite(bound) for bound in combined):
        return None
    return combined


def _certified_line(
    control_points: tuple[_XY, _XY, _XY, _XY],
    deviation_limit: float,
) -> tuple[ReferenceLine, float] | None:
    start, control1, control2, end = control_points
    start_exact = fraction_xy(start)
    end_exact = fraction_xy(end)
    control1_exact = fraction_xy(control1)
    control2_exact = fraction_xy(control2)
    chord = end_exact[0] - start_exact[0], end_exact[1] - start_exact[1]
    start_direction = control1_exact[0] - start_exact[0], control1_exact[1] - start_exact[1]
    end_direction = end_exact[0] - control2_exact[0], end_exact[1] - control2_exact[1]
    if start_direction[0] * chord[0] + start_direction[1] * chord[1] <= 0:
        return None
    if end_direction[0] * chord[0] + end_direction[1] * chord[1] <= 0:
        return None
    squared_bound = max(exact_segment_distance_squared(point, start, end) for point in control_points)
    exact_limit = fraction_xy((deviation_limit, deviation_limit))[0]
    if squared_bound > exact_limit**2:
        return None
    upper_bound = outward_sqrt_fraction(squared_bound)
    return ReferenceLine.build(_world_point(start), _world_point(end)), upper_bound


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
    start_parameter: Fraction,
    end_parameter: Fraction,
    depth: int,
) -> float | None:
    proof_biarc = tuple(
        (
            _point_xy(arc.start),
            _point_xy(arc.end),
            _point_xy(arc.centre),
            float(arc.sweep),
        )
        for arc in biarc
    )
    midpoint_parameter = (start_parameter + end_parameter) / 2
    if biarc_correspondence_sample_exceeds(
        control_points,
        (proof_biarc[0], proof_biarc[1]),
        midpoint_parameter,
        deviation_limit,
    ):
        return None
    enclosure = biarc_correspondence_node_bound(
        control_points,
        (proof_biarc[0], proof_biarc[1]),
        start_parameter,
        end_parameter,
    )
    if enclosure is None:
        return None
    if enclosure <= deviation_limit:
        return enclosure
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        return None
    left_bound = _cubic_within_biarc(
        control_points,
        biarc,
        deviation_limit,
        start_parameter=start_parameter,
        end_parameter=midpoint_parameter,
        depth=depth + 1,
    )
    if left_bound is None:
        return None
    right_bound = _cubic_within_biarc(
        control_points,
        biarc,
        deviation_limit,
        start_parameter=midpoint_parameter,
        end_parameter=end_parameter,
        depth=depth + 1,
    )
    if right_bound is None:
        return None
    return max(left_bound, right_bound)


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
) -> float | None:
    return _cubic_within_stored_circle(
        control_points,
        candidate,
        control_points[0],
        deviation_limit,
        depth=depth,
    )


def _cubic_within_stored_circle(
    control_points: tuple[_XY, _XY, _XY, _XY],
    candidate: tuple[_XY, float],
    radius_point: _XY,
    deviation_limit: float,
    *,
    depth: int,
) -> float | None:
    centre, _ = candidate
    enclosure = circle_deviation_upper_bound(control_points, centre, radius_point)
    if enclosure is None:
        return None
    if enclosure <= deviation_limit:
        return enclosure
    if depth >= MAX_RECONSTRUCTION_SUBDIVISIONS:
        return None
    left, right = _split_cubic(control_points)
    left_bound = _cubic_within_stored_circle(
        left,
        candidate,
        radius_point,
        deviation_limit,
        depth=depth + 1,
    )
    right_bound = _cubic_within_stored_circle(
        right,
        candidate,
        radius_point,
        deviation_limit,
        depth=depth + 1,
    )
    if left_bound is None or right_bound is None:
        return None
    return max(left_bound, right_bound)


def _cubic_has_monotone_polar_angle(
    control_points: tuple[_XY, _XY, _XY, _XY],
    centre: _XY,
    direction: float,
    *,
    depth: int,
) -> bool:
    return certify_monotone_polar_angle(
        control_points,
        centre,
        direction,
        MAX_RECONSTRUCTION_SUBDIVISIONS - depth,
    )


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


def _control_points_are_exactly_collinear(
    control_points: tuple[_XY, _XY, _XY, _XY],
) -> bool:
    exact = tuple(fraction_xy(point) for point in control_points)
    chord = exact[3][0] - exact[0][0], exact[3][1] - exact[0][1]
    return all(chord[0] * (point[1] - exact[0][1]) - chord[1] * (point[0] - exact[0][0]) == 0 for point in exact[1:3])


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


def _finite_normal_or_zero(value: float) -> bool:
    return math.isfinite(value) and (value == 0.0 or abs(value) >= sys.float_info.min)
