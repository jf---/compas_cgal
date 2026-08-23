"""One-shot typed ingress into native exact audit classification."""

from __future__ import annotations

from typing import assert_never

from compas_cgal import _stock_2
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.engagement_audit.errors import ContradictoryOperationOrientationError
from compas_cgal.engagement_audit.errors import ContradictoryOperationRoleError
from compas_cgal.engagement_audit.errors import InvalidAuditPlaneError
from compas_cgal.engagement_audit.errors import InvalidEngagementAuditInputError
from compas_cgal.engagement_audit.errors import MultipleCutPlaneError
from compas_cgal.engagement_audit.errors import NonFiniteAuditGeometryError
from compas_cgal.engagement_audit.errors import UnsupportedAuditGeometryError
from compas_cgal.engagement_audit.operation_identity import ArcOperationSnapshot
from compas_cgal.engagement_audit.operation_identity import CircleOperationSnapshot
from compas_cgal.engagement_audit.operation_identity import LineOperationSnapshot
from compas_cgal.engagement_audit.operation_identity import OperationSnapshot
from compas_cgal.engagement_audit.operation_identity import canonical_operation_snapshot_bytes
from compas_cgal.engagement_audit.operation_identity import operation_digest
from compas_cgal.engagement_audit.operation_identity import snapshot_toolpath_operation
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedNonEngagingOperation
from compas_cgal.engagement_audit.records import AuthenticatedOperation
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.toolpath import ToolpathOperation


def _point3_native(point: Point3[WorldXYZ]) -> tuple[float, float, float]:
    return float(point.x), float(point.y), float(point.z)


def _direction3_native(direction: Direction3[WorldXYZ]) -> tuple[float, float, float]:
    return float(direction.x), float(direction.y), float(direction.z)


def _translate_native_error(error: ValueError) -> InvalidEngagementAuditInputError:
    if isinstance(error, _stock_2.AuditInvalidPlaneError):
        return InvalidAuditPlaneError(str(error))
    if isinstance(error, _stock_2.AuditOffPlaneError):
        return MultipleCutPlaneError(str(error))
    if isinstance(error, _stock_2.AuditContradictoryOrientationError):
        return ContradictoryOperationOrientationError(str(error))
    if isinstance(error, _stock_2.AuditContradictoryRoleError):
        return ContradictoryOperationRoleError(str(error))
    if isinstance(error, _stock_2.AuditNonFiniteInputError):
        return NonFiniteAuditGeometryError(str(error))
    if isinstance(error, _stock_2.AuditUnsupportedGeometryError):
        return UnsupportedAuditGeometryError(str(error))
    return InvalidEngagementAuditInputError(str(error))


def _classify_line(
    snapshot: LineOperationSnapshot,
    cut_plane: CutPlane,
) -> _stock_2.AuditLineClassification2:
    return _stock_2.classify_audit_line(
        _point3_native(snapshot.start),
        _point3_native(snapshot.end),
        cut_plane.cut_z.value,
        cut_plane.clearance_z.value,
        snapshot.operation.value,
    )


def _classify_circle(
    snapshot: CircleOperationSnapshot,
    cut_plane: CutPlane,
) -> _stock_2.AuditCircleClassification2:
    return _stock_2.classify_audit_circle(
        _point3_native(snapshot.center),
        _direction3_native(snapshot.xaxis),
        _direction3_native(snapshot.yaxis),
        snapshot.radius.value,
        snapshot.clockwise,
        cut_plane.cut_z.value,
        cut_plane.clearance_z.value,
        snapshot.operation.value,
    )


def _classify_arc(
    snapshot: ArcOperationSnapshot,
    cut_plane: CutPlane,
) -> _stock_2.AuditArcClassification2:
    return _stock_2.classify_audit_arc(
        _point3_native(snapshot.center),
        _direction3_native(snapshot.xaxis),
        _direction3_native(snapshot.yaxis),
        snapshot.radius.value,
        snapshot.start_angle,
        snapshot.end_angle,
        snapshot.clockwise,
        cut_plane.cut_z.value,
        cut_plane.clearance_z.value,
        snapshot.operation.value,
    )


def _native_classification(
    snapshot: OperationSnapshot,
    cut_plane: CutPlane,
) -> _stock_2.AuditLineClassification2 | _stock_2.AuditCircleClassification2 | _stock_2.AuditArcClassification2:
    try:
        if isinstance(snapshot, LineOperationSnapshot):
            return _classify_line(snapshot, cut_plane)
        if isinstance(snapshot, CircleOperationSnapshot):
            return _classify_circle(snapshot, cut_plane)
        if isinstance(snapshot, ArcOperationSnapshot):
            return _classify_arc(snapshot, cut_plane)
        assert_never(snapshot)
    except (
        _stock_2.AuditInvalidPlaneError,
        _stock_2.AuditOffPlaneError,
        _stock_2.AuditContradictoryOrientationError,
        _stock_2.AuditContradictoryRoleError,
        _stock_2.AuditNonFiniteInputError,
        _stock_2.AuditUnsupportedGeometryError,
    ) as error:
        raise _translate_native_error(error) from None


def classify_operation_snapshot(
    snapshot: OperationSnapshot,
    cut_plane: CutPlane,
    *,
    operation_index: int,
) -> AuthenticatedOperation:
    """Classify one immutable ingress snapshot into an opaque native motion."""
    if type(cut_plane) is not CutPlane:
        raise InvalidEngagementAuditInputError("classification requires exact CutPlane.")
    if type(operation_index) is not int or operation_index < 0:
        raise InvalidEngagementAuditInputError("operation index must be an exact non-negative stream ordinal.")
    source_bytes = canonical_operation_snapshot_bytes(snapshot)
    digest = operation_digest(source_bytes)
    result = _native_classification(snapshot, cut_plane)
    if isinstance(
        result,
        (
            _stock_2.AuditSegmentMotion2,
            _stock_2.AuditCircleMotion2,
            _stock_2.AuditArcMotion2,
        ),
    ):
        return AuthenticatedLateralOperation.build(
            operation_index=operation_index,
            operation_digest=digest,
            motion=result,
        )
    if isinstance(result, _stock_2.AuditVerticalPlunge2):
        if not isinstance(snapshot, LineOperationSnapshot):
            raise InvalidEngagementAuditInputError("native plunge classification requires one line snapshot.")
        return AuthenticatedPlungeOperation.build(
            operation_index=operation_index,
            operation_digest=digest,
            motion=result,
        )
    if isinstance(
        result,
        (_stock_2.AuditVerticalRetract2, _stock_2.AuditClearanceTransport2),
    ):
        return AuthenticatedNonEngagingOperation.build(
            operation_index=operation_index,
            operation_digest=digest,
            motion=result,
        )
    assert_never(result)


def classify_operation(
    operation: ToolpathOperation,
    cut_plane: CutPlane,
    *,
    operation_index: int,
) -> AuthenticatedOperation:
    """Capture a mutable COMPAS operation once, then classify it natively."""
    return classify_operation_snapshot(
        snapshot_toolpath_operation(operation),
        cut_plane,
        operation_index=operation_index,
    )
