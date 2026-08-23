"""Strict consumer contract for authoritative audit identity and records."""

from typing import assert_never
from typing import assert_type

from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion_oracle_cache import NativeMotionVerdict
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.records import MeasuredOperationAudit
from compas_cgal.engagement_audit.records import NonEngagingOperationAudit
from compas_cgal.engagement_audit.records import NonEngagingReason
from compas_cgal.engagement_audit.records import OperationAudit


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
