"""Immutable scalar snapshots and canonical identity for COMPAS ingress."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Final
from typing import NewType
from typing import Self
from typing import TypeAlias

import numpy as np
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.engagement_audit.errors import InvalidAuditOperationError
from compas_cgal.engagement_audit.errors import NonFiniteAuditGeometryError
from compas_cgal.engagement_audit.records import OperationDigest
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

OPERATION_STREAM_VERSION: Final[bytes] = b"toolpath-operation-stream-v1"
OperationStreamDigest = NewType("OperationStreamDigest", bytes)
TangentSnapshot: TypeAlias = Direction3[WorldXYZ] | None


@dataclass(frozen=True)
class LineOperationSnapshot:
    """One coherent scalar observation of a mutable COMPAS line operation."""

    start: Point3[WorldXYZ]
    end: Point3[WorldXYZ]
    operation: OperationType
    path_index: int
    clockwise: bool
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __post_init__(self) -> None:
        if type(self.start) is not Point3 or type(self.end) is not Point3:
            raise InvalidAuditOperationError("line snapshot requires typed world-XYZ endpoints.")
        _validate_metadata_snapshot(self.operation, self.path_index, self.clockwise, self.start_tangent, self.end_tangent)

    @classmethod
    def build(
        cls,
        *,
        start: Point3[WorldXYZ],
        end: Point3[WorldXYZ],
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        return cls(start, end, operation, path_index, clockwise, start_tangent, end_tangent)


@dataclass(frozen=True)
class CircleOperationSnapshot:
    """One coherent scalar observation of a mutable COMPAS circle operation."""

    center: Point3[WorldXYZ]
    xaxis: Direction3[WorldXYZ]
    yaxis: Direction3[WorldXYZ]
    radius: GuideRadius
    operation: OperationType
    path_index: int
    clockwise: bool
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __post_init__(self) -> None:
        _validate_curve_snapshot(self.center, self.xaxis, self.yaxis, self.radius)
        _validate_metadata_snapshot(self.operation, self.path_index, self.clockwise, self.start_tangent, self.end_tangent)

    @classmethod
    def build(
        cls,
        *,
        center: Point3[WorldXYZ],
        xaxis: Direction3[WorldXYZ],
        yaxis: Direction3[WorldXYZ],
        radius: GuideRadius,
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        return cls(center, xaxis, yaxis, radius, operation, path_index, clockwise, start_tangent, end_tangent)


@dataclass(frozen=True)
class ArcOperationSnapshot:
    """One coherent scalar observation of a mutable COMPAS arc operation."""

    center: Point3[WorldXYZ]
    xaxis: Direction3[WorldXYZ]
    yaxis: Direction3[WorldXYZ]
    radius: GuideRadius
    start_angle: Radian
    end_angle: Radian
    operation: OperationType
    path_index: int
    clockwise: bool
    start_tangent: TangentSnapshot
    end_tangent: TangentSnapshot

    def __post_init__(self) -> None:
        _validate_curve_snapshot(self.center, self.xaxis, self.yaxis, self.radius)
        if type(self.start_angle) is not float or not math.isfinite(self.start_angle):
            raise InvalidAuditOperationError("arc snapshot requires one finite typed start angle.")
        if type(self.end_angle) is not float or not math.isfinite(self.end_angle):
            raise InvalidAuditOperationError("arc snapshot requires one finite typed end angle.")
        _validate_metadata_snapshot(self.operation, self.path_index, self.clockwise, self.start_tangent, self.end_tangent)

    @classmethod
    def build(
        cls,
        *,
        center: Point3[WorldXYZ],
        xaxis: Direction3[WorldXYZ],
        yaxis: Direction3[WorldXYZ],
        radius: GuideRadius,
        start_angle: Radian,
        end_angle: Radian,
        operation: OperationType,
        path_index: int,
        clockwise: bool,
        start_tangent: TangentSnapshot,
        end_tangent: TangentSnapshot,
    ) -> Self:
        return cls(
            center,
            xaxis,
            yaxis,
            radius,
            start_angle,
            end_angle,
            operation,
            path_index,
            clockwise,
            start_tangent,
            end_tangent,
        )


OperationSnapshot: TypeAlias = LineOperationSnapshot | CircleOperationSnapshot | ArcOperationSnapshot


def _validate_metadata_snapshot(
    operation: object,
    path_index: object,
    clockwise: object,
    start_tangent: object,
    end_tangent: object,
) -> None:
    if type(operation) is not OperationType:
        raise InvalidAuditOperationError("snapshot role must be exact OperationType.")
    if type(path_index) is not int or path_index < 0:
        raise InvalidAuditOperationError("snapshot path index must be an exact non-negative integer.")
    if type(clockwise) is not bool:
        raise InvalidAuditOperationError("snapshot orientation must be an exact bool.")
    if start_tangent is not None and type(start_tangent) is not Direction3:
        raise InvalidAuditOperationError("snapshot start tangent must be typed world-XYZ direction or absent.")
    if end_tangent is not None and type(end_tangent) is not Direction3:
        raise InvalidAuditOperationError("snapshot end tangent must be typed world-XYZ direction or absent.")


def _validate_curve_snapshot(
    center: object,
    xaxis: object,
    yaxis: object,
    radius: object,
) -> None:
    if type(center) is not Point3:
        raise InvalidAuditOperationError("curve snapshot center must be one typed world-XYZ point.")
    if type(xaxis) is not Direction3 or type(yaxis) is not Direction3:
        raise InvalidAuditOperationError("curve snapshot axes must be typed world-XYZ directions.")
    if type(radius) is not GuideRadius:
        raise InvalidAuditOperationError("curve snapshot radius must be one typed guide radius.")


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidAuditOperationError(f"{name} must be a real binary64 value.")
    try:
        numeric = float(value)
    except OverflowError:
        raise NonFiniteAuditGeometryError(f"{name} exceeds the finite binary64 range.") from None
    if not math.isfinite(numeric):
        raise NonFiniteAuditGeometryError(f"{name} must be finite.")
    return numeric


def _components3(point: object, name: str) -> tuple[float, float, float]:
    if not isinstance(point, Iterable) or isinstance(point, (str, bytes)):
        raise InvalidAuditOperationError(f"{name} must expose exactly three coordinates.")
    coordinates = tuple(point)
    if len(coordinates) != 3:
        raise InvalidAuditOperationError(f"{name} must expose exactly three coordinates.")
    return (
        _finite(coordinates[0], f"{name}[0]"),
        _finite(coordinates[1], f"{name}[1]"),
        _finite(coordinates[2], f"{name}[2]"),
    )


def _point3(point: object, name: str) -> Point3[WorldXYZ]:
    return Point3[WorldXYZ].build(_components3(point, name))


def _direction3(direction: object, name: str) -> Direction3[WorldXYZ]:
    return Direction3[WorldXYZ].build(_components3(direction, name))


def _tangent(tangent: object, name: str) -> TangentSnapshot:
    if tangent is None:
        return None
    if type(tangent) is not np.ndarray or tangent.shape != (3,):
        raise InvalidAuditOperationError(f"{name} must be None or an exact three-component ndarray.")
    return _direction3(tangent, name)


def _operation_metadata(
    operation: ToolpathOperation,
) -> tuple[OperationType, int, bool, TangentSnapshot, TangentSnapshot]:
    role = operation.operation
    path_index = operation.path_index
    clockwise = operation.clockwise
    start_tangent = _tangent(operation.start_tangent, "start tangent")
    end_tangent = _tangent(operation.end_tangent, "end tangent")
    if type(role) is not OperationType:
        raise InvalidAuditOperationError("operation role must be an exact OperationType.")
    if type(path_index) is not int or path_index < 0:
        raise InvalidAuditOperationError("path index must be an exact non-negative integer.")
    if type(clockwise) is not bool:
        raise InvalidAuditOperationError("operation orientation must be an exact bool.")
    return role, path_index, clockwise, start_tangent, end_tangent


def snapshot_toolpath_operation(operation: ToolpathOperation) -> OperationSnapshot:
    """Read one mutable COMPAS operation exactly once into immutable scalars."""
    if type(operation) is not ToolpathOperation:
        raise InvalidAuditOperationError("operation must be an exact ToolpathOperation, not a subclass.")
    role, path_index, clockwise, start_tangent, end_tangent = _operation_metadata(operation)
    geometry = operation.geometry
    if isinstance(geometry, Line):
        if type(geometry) is not Line:
            raise InvalidAuditOperationError("line geometry must be exact Line, not a subclass.")
        return LineOperationSnapshot.build(
            start=_point3(geometry.start, "line start"),
            end=_point3(geometry.end, "line end"),
            operation=role,
            path_index=path_index,
            clockwise=clockwise,
            start_tangent=start_tangent,
            end_tangent=end_tangent,
        )
    if isinstance(geometry, Arc):
        if type(geometry) is not Arc:
            raise InvalidAuditOperationError("arc geometry must be exact Arc, not a subclass.")
        frame = geometry.frame
        radius = _finite(geometry.radius, "arc radius")
        if radius <= 0.0:
            raise InvalidAuditOperationError("arc radius must be positive.")
        return ArcOperationSnapshot.build(
            center=_point3(frame.point, "arc center"),
            xaxis=_direction3(frame.xaxis, "arc frame X axis"),
            yaxis=_direction3(frame.yaxis, "arc frame Y axis"),
            radius=GuideRadius.build(radius),
            start_angle=Radian(_finite(geometry.start_angle, "arc start angle")),
            end_angle=Radian(_finite(geometry.end_angle, "arc end angle")),
            operation=role,
            path_index=path_index,
            clockwise=clockwise,
            start_tangent=start_tangent,
            end_tangent=end_tangent,
        )
    if isinstance(geometry, Circle):
        if type(geometry) is not Circle:
            raise InvalidAuditOperationError("circle geometry must be exact Circle, not a subclass.")
        frame = geometry.frame
        radius = _finite(geometry.radius, "circle radius")
        if radius <= 0.0:
            raise InvalidAuditOperationError("circle radius must be positive.")
        return CircleOperationSnapshot.build(
            center=_point3(frame.point, "circle center"),
            xaxis=_direction3(frame.xaxis, "circle frame X axis"),
            yaxis=_direction3(frame.yaxis, "circle frame Y axis"),
            radius=GuideRadius.build(radius),
            operation=role,
            path_index=path_index,
            clockwise=clockwise,
            start_tangent=start_tangent,
            end_tangent=end_tangent,
        )
    raise InvalidAuditOperationError(f"geometry type {type(geometry).__name__!r} is not a supported exact primitive.")


def _coordinate_bytes(x: float, y: float, z: float) -> bytes:
    return encode_sequence((encode_binary64(x), encode_binary64(y), encode_binary64(z)))


def _point3_bytes(point: Point3[WorldXYZ]) -> bytes:
    return _coordinate_bytes(point.x, point.y, point.z)


def _direction3_bytes(direction: Direction3[WorldXYZ]) -> bytes:
    return _coordinate_bytes(direction.x, direction.y, direction.z)


def _tangent_bytes(tangent: TangentSnapshot) -> bytes:
    if tangent is None:
        return encode_tagged_union(b"absent-tangent-v1", b"")
    return encode_tagged_union(b"tangent-v1", _direction3_bytes(tangent))


def _frame_bytes(
    center: Point3[WorldXYZ],
    xaxis: Direction3[WorldXYZ],
    yaxis: Direction3[WorldXYZ],
) -> bytes:
    return encode_component_map(
        {
            b"center": _point3_bytes(center),
            b"x-axis": _direction3_bytes(xaxis),
            b"y-axis": _direction3_bytes(yaxis),
        }
    )


def _geometry_bytes(snapshot: OperationSnapshot) -> bytes:
    if isinstance(snapshot, LineOperationSnapshot):
        return encode_tagged_union(
            b"compas-line-v1",
            encode_component_map(
                {
                    b"end": _point3_bytes(snapshot.end),
                    b"start": _point3_bytes(snapshot.start),
                }
            ),
        )
    if isinstance(snapshot, ArcOperationSnapshot):
        return encode_tagged_union(
            b"compas-arc-v1",
            encode_component_map(
                {
                    b"end-angle": encode_binary64(float(snapshot.end_angle)),
                    b"frame": _frame_bytes(snapshot.center, snapshot.xaxis, snapshot.yaxis),
                    b"radius": encode_binary64(float(snapshot.radius.value)),
                    b"start-angle": encode_binary64(float(snapshot.start_angle)),
                }
            ),
        )
    return encode_tagged_union(
        b"compas-circle-v1",
        encode_component_map(
            {
                b"frame": _frame_bytes(snapshot.center, snapshot.xaxis, snapshot.yaxis),
                b"radius": encode_binary64(float(snapshot.radius.value)),
            }
        ),
    )


def canonical_operation_snapshot_bytes(snapshot: OperationSnapshot) -> bytes:
    """Encode every audit-relevant scalar from one immutable observation."""
    return encode_tagged_union(
        b"toolpath-operation-v1",
        encode_component_map(
            {
                b"clockwise": encode_boolean(snapshot.clockwise),
                b"end-tangent": _tangent_bytes(snapshot.end_tangent),
                b"geometry": _geometry_bytes(snapshot),
                b"operation-role": snapshot.operation.value.encode("ascii"),
                b"path-index": encode_integer(snapshot.path_index),
                b"start-tangent": _tangent_bytes(snapshot.start_tangent),
            }
        ),
    )


def canonical_toolpath_operation_bytes(operation: ToolpathOperation) -> bytes:
    """Snapshot and encode one mutable legacy operation."""
    return canonical_operation_snapshot_bytes(snapshot_toolpath_operation(operation))


def operation_digest(source_bytes: bytes) -> OperationDigest:
    return OperationDigest(hashlib.sha256(source_bytes).digest())


def operation_stream_digest(source_operations: tuple[bytes, ...]) -> OperationStreamDigest:
    stream_bytes = encode_tagged_union(OPERATION_STREAM_VERSION, encode_sequence(source_operations))
    return OperationStreamDigest(hashlib.sha256(stream_bytes).digest())
