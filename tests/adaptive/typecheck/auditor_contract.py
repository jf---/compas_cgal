"""Strict consumer contract for authoritative audit identity and records."""

from typing import assert_never
from typing import assert_type

from compas_cgal import _stock_2
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion_oracle_cache import NativeMotionVerdict
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import ObservedRadius
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.operation_identity import ArcOperationSnapshot
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedNonEngagingOperation
from compas_cgal.engagement_audit.records import AuthenticatedOperation
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.engagement_audit.records import MeasuredOperationAudit
from compas_cgal.engagement_audit.records import NonEngagingOperationAudit
from compas_cgal.engagement_audit.records import NonEngagingReason
from compas_cgal.engagement_audit.records import OperationAudit
from compas_cgal.engagement_audit.records import SupportedLateralMotion


def consume_build_and_operation(
    build: BuildIdentity,
    operation: OperationAudit,
) -> IdentityDigest:
    if isinstance(operation, MeasuredOperationAudit):
        assert_type(operation.verdict, NativeMotionVerdict)
    elif isinstance(operation, NonEngagingOperationAudit):
        assert_type(operation.reason, NonEngagingReason)
    else:
        assert_never(operation)
    return build.digest


def consume_native_motion_request(motion: SupportedLateralMotion) -> None:
    if isinstance(motion, _stock_2.AuditSegmentMotion2):
        assert_type(motion, _stock_2.AuditSegmentMotion2)
    elif isinstance(motion, _stock_2.AuditArcMotion2):
        assert_type(motion, _stock_2.AuditArcMotion2)
    elif isinstance(motion, _stock_2.AuditCircleMotion2):
        assert_type(motion, _stock_2.AuditCircleMotion2)
    else:
        assert_never(motion)


def consume_authenticated_operation(operation: AuthenticatedOperation) -> None:
    assert_type(operation.canonical_bytes, bytes)
    assert_type(operation.digest, IdentityDigest)
    if isinstance(operation, AuthenticatedLateralOperation):
        consume_native_motion_request(operation.motion)
    elif isinstance(operation, AuthenticatedPlungeOperation):
        assert_type(operation.motion, _stock_2.AuditVerticalPlunge2)
    elif isinstance(operation, AuthenticatedNonEngagingOperation):
        assert_type(
            operation.motion,
            _stock_2.AuditVerticalRetract2 | _stock_2.AuditClearanceTransport2,
        )
    else:
        assert_never(operation)


def consume_arc_snapshot(snapshot: ArcOperationSnapshot) -> None:
    assert_type(snapshot.center, Point3[WorldXYZ])
    assert_type(snapshot.xaxis, Direction3[WorldXYZ])
    assert_type(snapshot.yaxis, Direction3[WorldXYZ])
    assert_type(snapshot.radius, ObservedRadius)
    assert_type(snapshot.start_angle, Radian)
    assert_type(snapshot.end_angle, Radian)
