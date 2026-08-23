"""Immutable, disjoint records for measured and non-engaging operations."""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Final
from typing import Literal
from typing import NewType
from typing import Self
from typing import TypeAlias

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion_oracle_cache import NativeMotionVerdict
from compas_cgal.adaptive.units import Radian
from compas_cgal.engagement_audit.errors import InvalidAuthenticatedLateralOperationError
from compas_cgal.engagement_audit.errors import InvalidAuthenticatedNonEngagingOperationError
from compas_cgal.engagement_audit.errors import InvalidAuthenticatedPlungeOperationError
from compas_cgal.engagement_audit.errors import InvalidMeasuredOperationAuditError
from compas_cgal.engagement_audit.errors import InvalidMotionVerdictError
from compas_cgal.engagement_audit.errors import InvalidNonEngagingOperationAuditError

MEASURED_OPERATION_AUDIT_VERSION: Final[bytes] = b"measured-operation-audit-v1"
NON_ENGAGING_OPERATION_AUDIT_VERSION: Final[bytes] = b"non-engaging-operation-audit-v1"
AUTHENTICATED_LATERAL_OPERATION_VERSION: Final[bytes] = b"authenticated-lateral-operation-v1"
AUTHENTICATED_PLUNGE_OPERATION_VERSION: Final[bytes] = b"authenticated-plunge-operation-v1"
AUTHENTICATED_NON_ENGAGING_OPERATION_VERSION: Final[bytes] = b"authenticated-non-engaging-operation-v1"

OperationDigest = NewType("OperationDigest", bytes)
StockLineageDigest = NewType("StockLineageDigest", bytes)
MotionCertificateDigest = NewType("MotionCertificateDigest", bytes)
NonEngagingReason: TypeAlias = Literal[
    "vertical_plunge",
    "vertical_retract",
    "clearance_transport",
]

_NATIVE_VERDICTS: Final[frozenset[str]] = frozenset({"certified", "cap_exceeded", "unresolved"})
_NON_ENGAGING_REASONS: Final[frozenset[str]] = frozenset({"vertical_plunge", "vertical_retract", "clearance_transport"})
SupportedLateralMotion: TypeAlias = _stock_2.AuditSegmentMotion2 | _stock_2.AuditCircleMotion2 | _stock_2.AuditArcMotion2
SupportedNonEngagingMotion: TypeAlias = _stock_2.AuditVerticalRetract2 | _stock_2.AuditClearanceTransport2


@dataclass(frozen=True)
class AuthenticatedLateralOperation:
    """Immutable classified motion tied to its source operation identity."""

    operation_index: int
    operation_digest: OperationDigest
    motion: SupportedLateralMotion

    def __post_init__(self) -> None:
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidAuthenticatedLateralOperationError("operation index must be an exact non-negative integer.")
        _require_digest(
            self.operation_digest,
            "operation digest",
            InvalidAuthenticatedLateralOperationError,
        )
        if type(self.motion) not in (
            _stock_2.AuditSegmentMotion2,
            _stock_2.AuditCircleMotion2,
            _stock_2.AuditArcMotion2,
        ):
            raise InvalidAuthenticatedLateralOperationError("lateral motion must be one exact supported motion type.")

    @classmethod
    def build(
        cls,
        *,
        operation_index: int,
        operation_digest: OperationDigest,
        motion: SupportedLateralMotion,
    ) -> Self:
        return cls(operation_index=operation_index, operation_digest=operation_digest, motion=motion)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not AuthenticatedLateralOperation:
            raise InvalidAuthenticatedLateralOperationError("authenticated lateral operation must be exact, not a subclass.")
        return _authenticated_operation_bytes(
            AUTHENTICATED_LATERAL_OPERATION_VERSION,
            self.operation_index,
            self.operation_digest,
            self.motion,
            InvalidAuthenticatedLateralOperationError,
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class AuthenticatedPlungeOperation:
    """Native-proved plunge retained as opaque exact geometry."""

    operation_index: int
    operation_digest: OperationDigest
    motion: _stock_2.AuditVerticalPlunge2

    def __post_init__(self) -> None:
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidAuthenticatedPlungeOperationError("operation index must be an exact non-negative integer.")
        _require_digest(
            self.operation_digest,
            "operation digest",
            InvalidAuthenticatedPlungeOperationError,
        )
        if type(self.motion) is not _stock_2.AuditVerticalPlunge2:
            raise InvalidAuthenticatedPlungeOperationError("plunge motion must be one native exact plunge value.")

    @classmethod
    def build(
        cls,
        *,
        operation_index: int,
        operation_digest: OperationDigest,
        motion: _stock_2.AuditVerticalPlunge2,
    ) -> Self:
        return cls(operation_index, operation_digest, motion)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not AuthenticatedPlungeOperation:
            raise InvalidAuthenticatedPlungeOperationError("authenticated plunge operation must be exact, not a subclass.")
        return _authenticated_operation_bytes(
            AUTHENTICATED_PLUNGE_OPERATION_VERSION,
            self.operation_index,
            self.operation_digest,
            self.motion,
            InvalidAuthenticatedPlungeOperationError,
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class AuthenticatedNonEngagingOperation:
    """Native-proved retract or clearance transport retained for replay."""

    operation_index: int
    operation_digest: OperationDigest
    motion: SupportedNonEngagingMotion

    def __post_init__(self) -> None:
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidAuthenticatedNonEngagingOperationError("operation index must be an exact non-negative integer.")
        _require_digest(
            self.operation_digest,
            "operation digest",
            InvalidAuthenticatedNonEngagingOperationError,
        )
        if type(self.motion) not in (
            _stock_2.AuditVerticalRetract2,
            _stock_2.AuditClearanceTransport2,
        ):
            raise InvalidAuthenticatedNonEngagingOperationError("non-engaging motion must be one supported native exact value.")

    @classmethod
    def build(
        cls,
        *,
        operation_index: int,
        operation_digest: OperationDigest,
        motion: SupportedNonEngagingMotion,
    ) -> Self:
        return cls(operation_index, operation_digest, motion)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not AuthenticatedNonEngagingOperation:
            raise InvalidAuthenticatedNonEngagingOperationError("authenticated non-engaging operation must be exact, not a subclass.")
        return _authenticated_operation_bytes(
            AUTHENTICATED_NON_ENGAGING_OPERATION_VERSION,
            self.operation_index,
            self.operation_digest,
            self.motion,
            InvalidAuthenticatedNonEngagingOperationError,
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


AuthenticatedOperation: TypeAlias = AuthenticatedLateralOperation | AuthenticatedPlungeOperation | AuthenticatedNonEngagingOperation


def _require_digest(value: object, name: str, error: type[ValueError]) -> bytes:
    if type(value) is not bytes or len(value) != hashlib.sha256().digest_size:
        raise error(f"{name} must be exactly one 32-byte SHA-256 digest.")
    return value


def _native_classification_tag(motion: object, error: type[ValueError]) -> bytes:
    if type(motion) is _stock_2.AuditSegmentMotion2:
        return b"segment"
    if type(motion) is _stock_2.AuditCircleMotion2:
        return b"circle"
    if type(motion) is _stock_2.AuditArcMotion2:
        return b"arc"
    if type(motion) is _stock_2.AuditVerticalPlunge2:
        return b"vertical-plunge"
    if type(motion) is _stock_2.AuditVerticalRetract2:
        return b"vertical-retract"
    if type(motion) is _stock_2.AuditClearanceTransport2:
        return b"clearance-transport"
    raise error("native motion is outside the closed authenticated classification domain.")


def _authenticated_operation_bytes(
    version: bytes,
    operation_index: int,
    operation_digest: OperationDigest,
    motion: object,
    error: type[ValueError],
) -> bytes:
    return encode_tagged_union(
        version,
        encode_component_map(
            {
                b"native-classification": _native_classification_tag(motion, error),
                b"operation-digest": bytes(operation_digest),
                b"operation-index": encode_integer(operation_index),
            }
        ),
    )


@dataclass(frozen=True)
class MeasuredOperationAudit:
    """Native evidence for one cut-height lateral operation."""

    operation_index: int
    operation_digest: OperationDigest
    verdict: NativeMotionVerdict
    max_tea: Radian
    station_count: int
    pre_motion_stock_lineage: StockLineageDigest
    motion_certificate_digest: MotionCertificateDigest

    def __post_init__(self) -> None:
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidMeasuredOperationAuditError("operation index must be an exact non-negative integer.")
        _require_digest(
            self.operation_digest,
            "operation digest",
            InvalidMeasuredOperationAuditError,
        )
        if type(self.verdict) is not str or self.verdict not in _NATIVE_VERDICTS:
            raise InvalidMotionVerdictError(f"native motion verdict {self.verdict!r} is foreign to the closed domain.")
        if type(self.max_tea) is not float or not math.isfinite(self.max_tea) or self.max_tea < 0.0:
            raise InvalidMeasuredOperationAuditError("maximum TEA must be a finite non-negative Radian.")
        if type(self.station_count) is not int or self.station_count <= 0:
            raise InvalidMeasuredOperationAuditError("station count must be an exact positive integer.")
        _require_digest(
            self.pre_motion_stock_lineage,
            "pre-motion stock lineage",
            InvalidMeasuredOperationAuditError,
        )
        _require_digest(
            self.motion_certificate_digest,
            "motion certificate digest",
            InvalidMeasuredOperationAuditError,
        )

    @classmethod
    def build(
        cls,
        *,
        operation_index: int,
        operation_digest: OperationDigest,
        verdict: NativeMotionVerdict,
        max_tea: Radian,
        station_count: int,
        pre_motion_stock_lineage: StockLineageDigest,
        motion_certificate_digest: MotionCertificateDigest,
    ) -> Self:
        return cls(
            operation_index=operation_index,
            operation_digest=operation_digest,
            verdict=verdict,
            max_tea=max_tea,
            station_count=station_count,
            pre_motion_stock_lineage=pre_motion_stock_lineage,
            motion_certificate_digest=motion_certificate_digest,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not MeasuredOperationAudit:
            raise InvalidMeasuredOperationAuditError("measured operation audit must be exact, not a subclass.")
        return encode_tagged_union(
            MEASURED_OPERATION_AUDIT_VERSION,
            encode_component_map(
                {
                    b"max-tea-radian": encode_binary64(float(self.max_tea)),
                    b"motion-certificate-digest": bytes(self.motion_certificate_digest),
                    b"operation-digest": bytes(self.operation_digest),
                    b"operation-index": encode_integer(self.operation_index),
                    b"pre-motion-stock-lineage": bytes(self.pre_motion_stock_lineage),
                    b"station-count": encode_integer(self.station_count),
                    b"verdict": self.verdict.encode("ascii"),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class NonEngagingOperationAudit:
    """Geometric proof that one operation cannot engage at the cut plane."""

    operation_index: int
    operation_digest: OperationDigest
    reason: NonEngagingReason

    def __post_init__(self) -> None:
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidNonEngagingOperationAuditError("operation index must be an exact non-negative integer.")
        _require_digest(
            self.operation_digest,
            "operation digest",
            InvalidNonEngagingOperationAuditError,
        )
        if type(self.reason) is not str or self.reason not in _NON_ENGAGING_REASONS:
            raise InvalidNonEngagingOperationAuditError(f"non-engaging reason {self.reason!r} is foreign to the closed domain.")

    @classmethod
    def build(
        cls,
        *,
        operation_index: int,
        operation_digest: OperationDigest,
        reason: NonEngagingReason,
    ) -> Self:
        return cls(
            operation_index=operation_index,
            operation_digest=operation_digest,
            reason=reason,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not NonEngagingOperationAudit:
            raise InvalidNonEngagingOperationAuditError("non-engaging operation audit must be exact, not a subclass.")
        return encode_tagged_union(
            NON_ENGAGING_OPERATION_AUDIT_VERSION,
            encode_component_map(
                {
                    b"operation-digest": bytes(self.operation_digest),
                    b"operation-index": encode_integer(self.operation_index),
                    b"reason": self.reason.encode("ascii"),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


OperationAudit: TypeAlias = MeasuredOperationAudit | NonEngagingOperationAudit
