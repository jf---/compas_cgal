"""Immutable structural snapshots of generated Held toolpath operations."""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional
from typing import Type
from typing import TypeVar
from typing import Union
from typing import cast

import numpy as np
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import angle_vectors
from compas.tolerance import TOL
from typing_extensions import Self
from typing_extensions import TypeAlias

from benchmarks.errors import InvalidHeldOperationSnapshotError
from benchmarks.errors import MutatedHeldToolpathError
from benchmarks.units import OperationIndex
from benchmarks.units import operation_index
from compas_cgal.helix import Helix
from compas_cgal.helix import HelixDescentError
from compas_cgal.helix import HelixTurnsError
from compas_cgal.incident_loop import INCIDENCE_CHART_DENOMINATOR
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.incident_loop import LoopEntryIncidenceError
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import Direction3
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Point3
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY
from compas_cgal.units import WorldXYZ

TangentSnapshot: TypeAlias = Optional[Direction3[WorldXYZ]]
SnapshotT = TypeVar("SnapshotT")


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidHeldOperationSnapshotError(f"{name} must be a finite real number.")
    try:
        numeric = float(value)
    except OverflowError:
        raise InvalidHeldOperationSnapshotError(f"{name} exceeds the binary64 range.") from None
    if not math.isfinite(numeric):
        raise InvalidHeldOperationSnapshotError(f"{name} must be finite.")
    return numeric


def _components3(value: object, *, name: str) -> tuple[float, float, float]:
    if not isinstance(value, Iterable) or isinstance(value, (str, bytes)):
        raise InvalidHeldOperationSnapshotError(f"{name} must expose exactly three finite coordinates.")
    components = tuple(value)
    if len(components) != 3:
        raise InvalidHeldOperationSnapshotError(f"{name} must expose exactly three finite coordinates.")
    return (
        _finite(components[0], name=f"{name}[0]"),
        _finite(components[1], name=f"{name}[1]"),
        _finite(components[2], name=f"{name}[2]"),
    )


def _point3(value: object, *, name: str) -> Point3[WorldXYZ]:
    return Point3[WorldXYZ].build(_components3(value, name=name))


def _direction3(value: object, *, name: str) -> Direction3[WorldXYZ]:
    return Direction3[WorldXYZ].build(_components3(value, name=name))


def _validate_point(value: object, *, name: str) -> None:
    if type(value) is not Point3:
        raise InvalidHeldOperationSnapshotError(f"{name} must be one typed world-XYZ point.")
    _finite(getattr(value, "x", None), name=f"{name}.x")
    _finite(getattr(value, "y", None), name=f"{name}.y")
    _finite(getattr(value, "z", None), name=f"{name}.z")


def _validate_direction(value: object, *, name: str) -> Direction3[WorldXYZ]:
    if type(value) is not Direction3:
        raise InvalidHeldOperationSnapshotError(f"{name} must be one typed world-XYZ direction.")
    x = _finite(getattr(value, "x", None), name=f"{name}.x")
    y = _finite(getattr(value, "y", None), name=f"{name}.y")
    z = _finite(getattr(value, "z", None), name=f"{name}.z")
    squared_norm = x * x + y * y + z * z
    if not TOL.is_between(squared_norm, 1.0, 1.0, atol=TOL.absolute):
        raise InvalidHeldOperationSnapshotError(f"{name} must be a unit direction.")
    return cast(Direction3[WorldXYZ], value)


def _validate_metadata(
    ordinal: object,
    operation: object,
    path_index: object,
    clockwise: object,
    start_tangent: object,
    end_tangent: object,
) -> None:
    if type(ordinal) is not int or ordinal < 0:
        raise InvalidHeldOperationSnapshotError("snapshot ordinal must be a non-negative operation index.")
    if type(operation) is not OperationType:
        raise InvalidHeldOperationSnapshotError("snapshot operation role must be an exact OperationType.")
    if type(path_index) is not int or path_index < 0:
        raise InvalidHeldOperationSnapshotError("snapshot path index must be an exact non-negative integer.")
    if type(clockwise) is not bool:
        raise InvalidHeldOperationSnapshotError("snapshot clockwise flag must be an exact bool.")
    if start_tangent is not None:
        _validate_direction(start_tangent, name="snapshot start tangent")
    if end_tangent is not None:
        _validate_direction(end_tangent, name="snapshot end tangent")


def _validate_curve(
    centre: object,
    xaxis: object,
    yaxis: object,
    radius: object,
) -> None:
    _validate_point(centre, name="snapshot centre")
    try:
        typed_xaxis = _validate_direction(xaxis, name="snapshot frame axes X")
        typed_yaxis = _validate_direction(yaxis, name="snapshot frame axes Y")
    except InvalidHeldOperationSnapshotError as error:
        raise InvalidHeldOperationSnapshotError(f"snapshot frame axes are invalid: {error}") from None
    dot = float(typed_xaxis.x) * float(typed_yaxis.x) + float(typed_xaxis.y) * float(typed_yaxis.y) + float(typed_xaxis.z) * float(typed_yaxis.z)
    if not TOL.is_between(dot, 0.0, 0.0, atol=TOL.absolute):
        raise InvalidHeldOperationSnapshotError("snapshot frame axes must be orthogonal.")
    numeric_radius = _finite(radius, name="snapshot radius")
    if numeric_radius <= 0.0:
        raise InvalidHeldOperationSnapshotError("snapshot radius must be positive.")


def _validate_angles(start_angle: object, end_angle: object) -> None:
    angles: list[float] = []
    for value, name in (
        (start_angle, "snapshot start angle"),
        (end_angle, "snapshot end angle"),
    ):
        angle = _finite(value, name=name)
        if not 0.0 <= angle <= math.tau:
            raise InvalidHeldOperationSnapshotError(f"{name} must lie in [0, math.tau].")
        angles.append(angle)
    if angles[1] < angles[0]:
        raise InvalidHeldOperationSnapshotError("snapshot end angle must not precede its start angle.")


def _direction_components(direction: Direction3[WorldXYZ]) -> tuple[float, float, float]:
    return float(direction.x), float(direction.y), float(direction.z)


def _validate_tangent_agreement(
    tangent: TangentSnapshot,
    expected: tuple[float, float, float] | None,
    *,
    name: str,
) -> None:
    if tangent is None:
        return
    if expected is None or not TOL.is_angle_zero(angle_vectors(_direction_components(tangent), expected)):
        raise InvalidHeldOperationSnapshotError(f"{name} disagrees with the primitive's travel direction.")


def _line_tangent(
    start: Point3[WorldXYZ],
    end: Point3[WorldXYZ],
) -> tuple[float, float, float] | None:
    delta = (
        float(end.x) - float(start.x),
        float(end.y) - float(start.y),
        float(end.z) - float(start.z),
    )
    if delta == (0.0, 0.0, 0.0):
        return None
    length = math.hypot(*delta)
    return delta[0] / length, delta[1] / length, delta[2] / length


def _curve_tangent(
    xaxis: Direction3[WorldXYZ],
    yaxis: Direction3[WorldXYZ],
    angle: float,
    *,
    clockwise: bool,
) -> tuple[float, float, float]:
    turn = -1.0 if clockwise else 1.0
    sine = math.sin(angle)
    cosine = math.cos(angle)
    x_components = _direction_components(xaxis)
    y_components = _direction_components(yaxis)
    return (
        turn * (-sine * x_components[0] + cosine * y_components[0]),
        turn * (-sine * x_components[1] + cosine * y_components[1]),
        turn * (-sine * x_components[2] + cosine * y_components[2]),
    )


def _build_record(record_type: Type[SnapshotT], values: dict[str, object]) -> SnapshotT:
    record = object.__new__(record_type)
    for name, value in values.items():
        object.__setattr__(record, name, value)
    return record


@dataclass(frozen=True, init=False)
class HeldLineSnapshot:
    """One immutable, world-XYZ line-operation observation."""

    ordinal: OperationIndex
    operation: OperationType
    path_index: int
    clockwise: bool
    start: Point3[WorldXYZ]
    end: Point3[WorldXYZ]
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __init__(self) -> None:
        raise TypeError("HeldLineSnapshot must be created with HeldLineSnapshot.build().")

    @classmethod
    def build(
        cls,
        *,
        ordinal: OperationIndex,
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        start: Point3[WorldXYZ],
        end: Point3[WorldXYZ],
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        _validate_metadata(ordinal, operation, path_index, clockwise, start_tangent, end_tangent)
        _validate_point(start, name="line start")
        _validate_point(end, name="line end")
        expected_tangent = _line_tangent(start, end)
        _validate_tangent_agreement(start_tangent, expected_tangent, name="snapshot start tangent")
        _validate_tangent_agreement(end_tangent, expected_tangent, name="snapshot end tangent")
        return _build_record(
            cls,
            {
                "ordinal": ordinal,
                "operation": operation,
                "path_index": path_index,
                "clockwise": clockwise,
                "start": start,
                "end": end,
                "start_tangent": start_tangent,
                "end_tangent": end_tangent,
            },
        )


@dataclass(frozen=True, init=False)
class HeldArcSnapshot:
    """One immutable, world-XYZ circular-arc operation observation."""

    ordinal: OperationIndex
    operation: OperationType
    path_index: int
    clockwise: bool
    centre: Point3[WorldXYZ]
    xaxis: Direction3[WorldXYZ]
    yaxis: Direction3[WorldXYZ]
    radius: Millimetre
    start_angle: Radian
    end_angle: Radian
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __init__(self) -> None:
        raise TypeError("HeldArcSnapshot must be created with HeldArcSnapshot.build().")

    @classmethod
    def build(
        cls,
        *,
        ordinal: OperationIndex,
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        centre: Point3[WorldXYZ],
        xaxis: Direction3[WorldXYZ],
        yaxis: Direction3[WorldXYZ],
        radius: Millimetre,
        start_angle: Radian,
        end_angle: Radian,
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        _validate_metadata(ordinal, operation, path_index, clockwise, start_tangent, end_tangent)
        _validate_curve(centre, xaxis, yaxis, radius)
        _validate_angles(start_angle, end_angle)
        _validate_tangent_agreement(
            start_tangent,
            _curve_tangent(xaxis, yaxis, float(start_angle), clockwise=clockwise),
            name="snapshot start tangent",
        )
        _validate_tangent_agreement(
            end_tangent,
            _curve_tangent(xaxis, yaxis, float(end_angle), clockwise=clockwise),
            name="snapshot end tangent",
        )
        return _build_record(
            cls,
            {
                "ordinal": ordinal,
                "operation": operation,
                "path_index": path_index,
                "clockwise": clockwise,
                "centre": centre,
                "xaxis": xaxis,
                "yaxis": yaxis,
                "radius": radius,
                "start_angle": start_angle,
                "end_angle": end_angle,
                "start_tangent": start_tangent,
                "end_tangent": end_tangent,
            },
        )


@dataclass(frozen=True, init=False)
class HeldCircleSnapshot:
    """One immutable, world-XYZ circle-operation observation."""

    ordinal: OperationIndex
    operation: OperationType
    path_index: int
    clockwise: bool
    centre: Point3[WorldXYZ]
    xaxis: Direction3[WorldXYZ]
    yaxis: Direction3[WorldXYZ]
    radius: Millimetre
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __init__(self) -> None:
        raise TypeError("HeldCircleSnapshot must be created with HeldCircleSnapshot.build().")

    @classmethod
    def build(
        cls,
        *,
        ordinal: OperationIndex,
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        centre: Point3[WorldXYZ],
        xaxis: Direction3[WorldXYZ],
        yaxis: Direction3[WorldXYZ],
        radius: Millimetre,
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        _validate_metadata(ordinal, operation, path_index, clockwise, start_tangent, end_tangent)
        _validate_curve(centre, xaxis, yaxis, radius)
        expected_tangent = _curve_tangent(xaxis, yaxis, 0.0, clockwise=clockwise)
        _validate_tangent_agreement(start_tangent, expected_tangent, name="snapshot start tangent")
        _validate_tangent_agreement(end_tangent, expected_tangent, name="snapshot end tangent")
        return _build_record(
            cls,
            {
                "ordinal": ordinal,
                "operation": operation,
                "path_index": path_index,
                "clockwise": clockwise,
                "centre": centre,
                "xaxis": xaxis,
                "yaxis": yaxis,
                "radius": radius,
                "start_tangent": start_tangent,
                "end_tangent": end_tangent,
            },
        )


class InvalidHeldHelixSnapshotError(InvalidHeldOperationSnapshotError):
    """A HELIX operation or helix record is malformed or contradicts itself."""


HELIX_CHART_COUNT = 4  # Quarter charts of `IncidentLoop.chart`, [0, 3]; incident_loop exposes no public count (`_QUARTERS` is private).


def _helix_point(value: object, *, name: str) -> Point2[WorldXY]:
    if type(value) is not Point2:
        raise InvalidHeldHelixSnapshotError(f"{name} must be one typed world-XY point.")
    _finite(getattr(value, "x", None), name=f"{name}.x")
    _finite(getattr(value, "y", None), name=f"{name}.y")
    return cast(Point2[WorldXY], value)


def _footprint_chart_agrees(footprint: IncidentLoop) -> bool:
    """The exact check that chart, numerator and phase describe the entry direction.

    Mirrors the public atlas definition documented in `compas_cgal.incident_loop`
    (quarter 0 is ``((Q^2 - j^2), 2 j Q) / (Q^2 + j^2)`` reduced, each further
    quarter a quarter-turn rotation; ``Q = INCIDENCE_CHART_DENOMINATOR``):
    ``(P - c) * m == (a, b) * rho`` over the rationals, and ``phase`` is the
    authored ``atan2(b, a)``.
    """
    q = INCIDENCE_CHART_DENOMINATOR
    j = footprint.chart_numerator
    a, b, m = q * q - j * j, 2 * j * q, q * q + j * j
    divisor = math.gcd(a, b)
    a, b, m = a // divisor, b // divisor, m // divisor
    for _ in range(footprint.chart):
        a, b = -b, a
    rho = Fraction(footprint.radius.value)
    dx = Fraction(footprint.entry.x) - Fraction(footprint.centre.x)
    dy = Fraction(footprint.entry.y) - Fraction(footprint.centre.y)
    return dx * m == a * rho and dy * m == b * rho and float(footprint.phase) == math.atan2(b, a)


def _helix_tangent(helix: Helix) -> Direction3[WorldXYZ]:
    """The unit XY travel direction at the helix entry (= exit), z = 0.

    The emitter convention (`compas_cgal.pipeline.emit.helix_operation`, like
    every emitter): the tangent is horizontal; the descent slope is not in it,
    it is exact in ``Helix.pitch``. Doubles, used only for the TOL-gated
    tangent agreement this module applies to every primitive.
    """
    footprint = helix.footprint
    dx = float(footprint.entry.x) - float(footprint.centre.x)
    dy = float(footprint.entry.y) - float(footprint.centre.y)
    length = math.hypot(dx, dy)
    ux, uy = dx / length, dy / length
    tx, ty = (uy, -ux) if helix.clockwise else (-uy, ux)
    return Direction3[WorldXYZ].build(tx, ty, 0.0)


@dataclass(frozen=True, init=False)
class HeldHelixSnapshot:
    """One immutable helical-entry observation.

    The record carries every defining value of the core `Helix` as emitted:
    the footprint loop's centre, exact entry ``P``, radius, phase and chart
    parameter (doubles, each an exact rational), the top and bottom Z, the
    integer turn count and the exact pitch ``(top_z - bottom_z) / turns`` as a
    `Fraction`. `helix` rebuilds the core value; the round trip is equality.
    """

    ordinal: OperationIndex
    operation: OperationType
    path_index: int
    clockwise: bool
    centre: Point2[WorldXY]
    entry: Point2[WorldXY]
    radius: Millimetre
    phase: Radian
    chart: int
    chart_numerator: int
    top_z: Millimetre
    bottom_z: Millimetre
    turns: int
    pitch: ExactMillimetre
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __init__(self) -> None:
        raise TypeError("HeldHelixSnapshot must be created with HeldHelixSnapshot.build().")

    @classmethod
    def build(
        cls,
        *,
        ordinal: OperationIndex,
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        centre: Point2[WorldXY],
        entry: Point2[WorldXY],
        radius: Millimetre,
        phase: Radian,
        chart: int,
        chart_numerator: int,
        top_z: Millimetre,
        bottom_z: Millimetre,
        turns: int,
        pitch: ExactMillimetre,
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        """Validate one helix record.

        Raises:
            InvalidHeldHelixSnapshotError: the role is not HELIX, a value is
                untyped or non-finite, the turn count is not a positive integer,
                the descent is not positive, the entry is not exactly on the
                footprint circle, the chart numerator is outside ``[0, Q]``, the
                chart, numerator or phase contradicts the exact entry direction,
                the pitch is not the exact ``descent / turns``, or a tangent is
                not the unit XY travel direction (z = 0) at the entry.
        """
        try:
            _validate_metadata(ordinal, operation, path_index, clockwise, start_tangent, end_tangent)
        except InvalidHeldOperationSnapshotError as error:
            raise InvalidHeldHelixSnapshotError(str(error)) from None
        if operation is not OperationType.HELIX:
            raise InvalidHeldHelixSnapshotError(f"a helix record must carry the HELIX role, not {operation!r}.")
        if type(chart) is not int or not 0 <= chart < HELIX_CHART_COUNT or type(chart_numerator) is not int or not 0 <= chart_numerator <= INCIDENCE_CHART_DENOMINATOR:
            raise InvalidHeldHelixSnapshotError("helix footprint chart parameter must be exact integers in range.")
        if type(turns) is not int:
            raise InvalidHeldHelixSnapshotError(f"helix turns must be an exact int; got {turns!r}.")
        if type(pitch) is not Fraction:
            raise InvalidHeldHelixSnapshotError("helix pitch must be one exact Fraction.")
        try:
            numeric_radius = _finite(radius, name="helix footprint radius")
            numeric_phase = _finite(phase, name="helix footprint phase")
            numeric_top = _finite(top_z, name="helix top Z")
            numeric_bottom = _finite(bottom_z, name="helix bottom Z")
        except InvalidHeldOperationSnapshotError as error:
            raise InvalidHeldHelixSnapshotError(str(error)) from None
        if numeric_radius <= 0.0:
            raise InvalidHeldHelixSnapshotError("helix footprint radius must be positive.")
        footprint = IncidentLoop(
            centre=_helix_point(centre, name="helix footprint centre"),
            radius=GuideRadius(Millimetre(numeric_radius)),
            entry=_helix_point(entry, name="helix footprint entry"),
            phase=Radian(numeric_phase),
            chart=chart,
            chart_numerator=chart_numerator,
        )
        if not footprint.certifies_incidence():
            raise InvalidHeldHelixSnapshotError(f"helix footprint entry {footprint.entry!r} is not exactly on its circle.")
        if not _footprint_chart_agrees(footprint):
            raise InvalidHeldHelixSnapshotError("helix footprint chart, numerator or phase contradicts its exact entry direction.")
        try:
            helix = Helix(
                footprint=footprint,
                top_z=ClearanceZ(Millimetre(numeric_top)),
                bottom_z=CutZ(Millimetre(numeric_bottom)),
                turns=turns,
                clockwise=clockwise,
            )
        except (HelixTurnsError, HelixDescentError, LoopEntryIncidenceError) as error:
            raise InvalidHeldHelixSnapshotError(f"helix record is not a valid helix: {error}") from None
        if pitch != helix.pitch:
            raise InvalidHeldHelixSnapshotError(f"helix pitch {pitch!r} is not the exact descent / turns {helix.pitch!r}.")
        expected_tangent = _direction_components(_helix_tangent(helix))
        try:
            _validate_tangent_agreement(start_tangent, expected_tangent, name="snapshot start tangent")
            _validate_tangent_agreement(end_tangent, expected_tangent, name="snapshot end tangent")
        except InvalidHeldOperationSnapshotError as error:
            raise InvalidHeldHelixSnapshotError(str(error)) from None
        return _build_record(
            cls,
            {
                "ordinal": ordinal,
                "operation": operation,
                "path_index": path_index,
                "clockwise": clockwise,
                "centre": footprint.centre,
                "entry": footprint.entry,
                "radius": footprint.radius.value,
                "phase": footprint.phase,
                "chart": chart,
                "chart_numerator": chart_numerator,
                "top_z": helix.top_z.value,
                "bottom_z": helix.bottom_z.value,
                "turns": turns,
                "pitch": helix.pitch,
                "start_tangent": start_tangent,
                "end_tangent": end_tangent,
            },
        )

    def helix(self) -> Helix:
        """The core `Helix` this record observed; equal to the emitted one."""
        footprint = IncidentLoop(
            centre=self.centre,
            radius=GuideRadius(self.radius),
            entry=self.entry,
            phase=self.phase,
            chart=self.chart,
            chart_numerator=self.chart_numerator,
        )
        return Helix(
            footprint=footprint,
            top_z=ClearanceZ(self.top_z),
            bottom_z=CutZ(self.bottom_z),
            turns=self.turns,
            clockwise=self.clockwise,
        )

    def analytic_tangent(self) -> Direction3[WorldXYZ]:
        """The unit travel direction at the entry, as the tangent check uses it."""
        return _helix_tangent(self.helix())


HeldOperationSnapshot: TypeAlias = Union[HeldLineSnapshot, HeldArcSnapshot, HeldCircleSnapshot, HeldHelixSnapshot]


def _tangent(value: object, *, name: str, operation: object) -> TangentSnapshot:
    if value is None:
        return None
    if type(value) is not np.ndarray or value.shape != (3,):
        raise InvalidHeldOperationSnapshotError(f"{name} must be None or an exact three-component ndarray.")
    tangent = _direction3(value, name=name)
    if operation is OperationType.PLUNGE or operation is OperationType.RETRACT:
        if tangent.x == 0.0 and tangent.y == 0.0 and tangent.z == 0.0:
            return None
    _validate_direction(tangent, name=name)
    return tangent


def _metadata(
    source: ToolpathOperation,
    *,
    ordinal: OperationIndex,
) -> tuple[OperationType, int, bool, TangentSnapshot, TangentSnapshot]:
    operation = source.operation
    path_index = source.path_index
    clockwise = source.clockwise
    start_tangent = _tangent(source.start_tangent, name="start tangent", operation=operation)
    end_tangent = _tangent(source.end_tangent, name="end tangent", operation=operation)
    _validate_metadata(ordinal, operation, path_index, clockwise, start_tangent, end_tangent)
    return operation, path_index, clockwise, start_tangent, end_tangent


def _snapshot_helix(
    source: ToolpathOperation,
    geometry: Helix,
    *,
    ordinal: OperationIndex,
    start_tangent: TangentSnapshot,
    end_tangent: TangentSnapshot,
) -> HeldHelixSnapshot:
    if type(geometry) is not Helix or type(geometry.footprint) is not IncidentLoop:
        raise InvalidHeldHelixSnapshotError("helix geometry must be an exact Helix on an exact IncidentLoop, not a subclass.")
    if source.operation is not OperationType.HELIX:
        raise InvalidHeldHelixSnapshotError(f"Helix geometry must carry the HELIX role, not {source.operation!r}.")
    if source.clockwise is not geometry.clockwise:
        raise InvalidHeldHelixSnapshotError("operation clockwise flag contradicts the helix turn sense.")
    entry = geometry.entry_xy
    if source.entry_xy is None or tuple(source.entry_xy) != (entry.x, entry.y):
        raise InvalidHeldHelixSnapshotError(f"a HELIX operation must carry its exact entry {(entry.x, entry.y)!r} as entry_xy; got {source.entry_xy!r}.")
    footprint = geometry.footprint
    return HeldHelixSnapshot.build(
        ordinal=ordinal,
        operation=source.operation,
        path_index=source.path_index,
        clockwise=geometry.clockwise,
        centre=footprint.centre,
        entry=footprint.entry,
        radius=footprint.radius.value,
        phase=footprint.phase,
        chart=footprint.chart,
        chart_numerator=footprint.chart_numerator,
        top_z=geometry.top_z.value,
        bottom_z=geometry.bottom_z.value,
        turns=geometry.turns,
        pitch=geometry.pitch,
        start_tangent=start_tangent,
        end_tangent=end_tangent,
    )


def _snapshot_operation(
    source: ToolpathOperation,
    *,
    ordinal: OperationIndex,
) -> HeldOperationSnapshot:
    if type(source) is not ToolpathOperation:
        raise InvalidHeldOperationSnapshotError("operation must be an exact ToolpathOperation, not a subclass.")
    operation, path_index, clockwise, start_tangent, end_tangent = _metadata(source, ordinal=ordinal)
    geometry = source.geometry
    if isinstance(geometry, Helix):
        return _snapshot_helix(source, geometry, ordinal=ordinal, start_tangent=start_tangent, end_tangent=end_tangent)
    if operation is OperationType.HELIX:
        raise InvalidHeldHelixSnapshotError(f"a HELIX operation must carry Helix geometry, not {type(geometry).__name__!r}.")
    if type(geometry) is Line:
        return HeldLineSnapshot.build(
            ordinal=ordinal,
            operation=operation,
            path_index=path_index,
            clockwise=clockwise,
            start=_point3(geometry.start, name="line start"),
            end=_point3(geometry.end, name="line end"),
            start_tangent=start_tangent,
            end_tangent=end_tangent,
        )
    if type(geometry) is Arc:
        frame = geometry.frame
        return HeldArcSnapshot.build(
            ordinal=ordinal,
            operation=operation,
            path_index=path_index,
            clockwise=clockwise,
            centre=_point3(frame.point, name="arc centre"),
            xaxis=_direction3(frame.xaxis, name="arc frame X axis"),
            yaxis=_direction3(frame.yaxis, name="arc frame Y axis"),
            radius=Millimetre(_finite(geometry.radius, name="arc radius")),
            start_angle=Radian(_finite(geometry.start_angle, name="arc start angle")),
            end_angle=Radian(_finite(geometry.end_angle, name="arc end angle")),
            start_tangent=start_tangent,
            end_tangent=end_tangent,
        )
    if type(geometry) is Circle:
        frame = geometry.frame
        return HeldCircleSnapshot.build(
            ordinal=ordinal,
            operation=operation,
            path_index=path_index,
            clockwise=clockwise,
            centre=_point3(frame.point, name="circle centre"),
            xaxis=_direction3(frame.xaxis, name="circle frame X axis"),
            yaxis=_direction3(frame.yaxis, name="circle frame Y axis"),
            radius=Millimetre(_finite(geometry.radius, name="circle radius")),
            start_tangent=start_tangent,
            end_tangent=end_tangent,
        )
    raise InvalidHeldOperationSnapshotError(f"geometry type {type(geometry).__name__!r} is not a supported exact primitive.")


def snapshot_toolpath(result: ToolpathResult) -> tuple[HeldOperationSnapshot, ...]:
    """Read a mutable toolpath once into behavior-complete immutable records."""
    if not isinstance(result, ToolpathResult):
        raise InvalidHeldOperationSnapshotError("result must be a ToolpathResult.")
    if type(result.operations) is not list:
        raise InvalidHeldOperationSnapshotError("result operations must be one exact list.")
    operations = tuple(result.operations)
    operation_count = len(operations)
    if operation_count == 0:
        return ()
    return tuple(_snapshot_operation(source, ordinal=operation_index(index, operation_count=operation_count)) for index, source in enumerate(operations))


def assert_toolpath_matches_snapshot(
    result: ToolpathResult,
    snapshot: tuple[HeldOperationSnapshot, ...],
) -> None:
    """Reject any structural change to the characterized operation stream."""
    try:
        observed = snapshot_toolpath(result)
    except InvalidHeldOperationSnapshotError as error:
        raise MutatedHeldToolpathError("Generated toolpath differs from the characterized operation snapshot.") from error
    if observed != snapshot:
        raise MutatedHeldToolpathError("Generated toolpath differs from the characterized operation snapshot.")
