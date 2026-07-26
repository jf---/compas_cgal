import hashlib
from dataclasses import dataclass
from typing import Literal
from typing import cast

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import InvalidMotionCertificateError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.toolpath import OperationType

_DIGEST_SIZE = hashlib.sha256().digest_size


def _digest_sequence(domain: bytes, values: tuple[bytes, ...]) -> bytes:
    return hashlib.sha256(encode_tagged_union(domain, encode_sequence(values))).digest()


def _boundary_digest(stock: _stock_2.Stock2) -> bytes:
    feature_ids = tuple(record.feature_id for record in _continuous_tea_2.extract_boundary_records(stock))
    return _digest_sequence(b"exact-stock-boundary-v1", feature_ids)


def _lineage_digest(stock: Stock2Area) -> bytes:
    return _digest_sequence(
        b"exact-stock-lineage-v1",
        tuple(bytes(witness.digest) for witness in stock.lineage),
    )


@dataclass(frozen=True)
class MotionWitness:
    operation_index: int
    operation_kind: OperationType
    motion: ExactSegmentMotion | ExactCircleMotion
    user_cap_bytes: bytes
    effective_cap_bytes: bytes
    strategy_identity: bytes
    stock_lineage_digest: bytes
    event_trace_digest: bytes
    verdict: Literal["certified"]
    event_cell_count: int
    unresolved_count: int

    def __post_init__(self) -> None:
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidMotionCertificateError("motion witness operation index must be nonnegative.")
        if type(self.operation_kind) is not OperationType:
            raise InvalidMotionCertificateError("motion witness operation kind must be exact OperationType.")
        if type(self.motion) not in (ExactSegmentMotion, ExactCircleMotion):
            raise InvalidMotionCertificateError("motion witness requires one exact motion.")
        if any(type(value) is not bytes or len(value) != _DIGEST_SIZE for value in (self.stock_lineage_digest, self.event_trace_digest)):
            raise InvalidMotionCertificateError("motion witness digests must be exact SHA-256 bytes.")
        if type(self.strategy_identity) is not bytes or not self.strategy_identity:
            raise InvalidMotionCertificateError("motion witness requires a native strategy identity.")
        if self.verdict != "certified" or self.unresolved_count != 0:
            raise InvalidMotionCertificateError("motion witness may represent only resolved certification.")
        if type(self.event_cell_count) is not int or self.event_cell_count < 0:
            raise InvalidMotionCertificateError("motion witness event-cell count must be nonnegative.")


@dataclass(frozen=True)
class MotionCertifier:
    _stock: _stock_2.Stock2
    tool_radius: ToolRadius
    stock_lineage_digest: bytes
    canonical_boundary_digest: bytes

    def __post_init__(self) -> None:
        if type(self._stock) is not _stock_2.Stock2 or type(self.tool_radius) is not ToolRadius:
            raise InvalidMotionCertificateError("motion certifier requires exact stock and tool radius.")
        if any(type(value) is not bytes or len(value) != _DIGEST_SIZE for value in (self.stock_lineage_digest, self.canonical_boundary_digest)):
            raise InvalidMotionCertificateError("motion certifier state digests must be exact SHA-256 bytes.")

    @classmethod
    def build(
        cls,
        *,
        stock: Stock2Area,
        tool_radius: ToolRadius,
    ) -> "MotionCertifier":
        if type(stock) is not Stock2Area or type(tool_radius) is not ToolRadius:
            raise InvalidMotionCertificateError("motion certifier factory requires exact stock and tool radius.")
        snapshot = stock.raw
        return cls(
            snapshot,
            tool_radius,
            _lineage_digest(stock),
            _boundary_digest(snapshot),
        )

    def contains(self, x: float, y: float) -> bool:
        return self._stock.contains(x, y)

    def certify(
        self,
        *,
        operation_index: int,
        operation_kind: OperationType,
        motion: ExactSegmentMotion | ExactCircleMotion,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> MotionWitness:
        if type(operation_index) is not int or operation_index < 0:
            raise InvalidMotionCertificateError("operation index must be a nonnegative exact integer.")
        if type(operation_kind) is not OperationType:
            raise InvalidMotionCertificateError("operation kind must be exact OperationType.")
        if operation_kind not in (OperationType.CUT, OperationType.LINK):
            raise InvalidMotionCertificateError("motion certification is restricted to lateral operations.")
        if type(user_cap) is not EngagementCap or type(effective_cap) is not EngagementCap:
            raise InvalidMotionCertificateError("motion certification requires exact engagement caps.")
        if not _stock_2.cap_chord_ratio_le(effective_cap.chord_ratio, user_cap.chord_ratio):
            raise InvalidMotionCertificateError("effective cap exceeds the user cap.")
        if operation_kind is OperationType.LINK and type(motion) is not ExactSegmentMotion:
            raise InvalidMotionCertificateError("operation kind and motion are incompatible.")
        if operation_kind is OperationType.CUT and type(motion) is not ExactCircleMotion:
            raise InvalidMotionCertificateError("operation kind and motion are incompatible.")
        if type(motion) is ExactSegmentMotion:
            raise UnresolvedMotionEventError("event-exact segment oracle is absent from the current Task 6 substrate.")

        circle = cast(ExactCircleMotion, motion)
        verdict, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
            self._stock,
            circle.center.x,
            circle.center.y,
            circle.phase_vector.x,
            circle.phase_vector.y,
            circle.clockwise,
            self.tool_radius.value,
            effective_cap.chord_ratio,
        )
        if verdict == "cap_exceeded":
            raise EngagementCapExceededError(f"operation {operation_index} exceeds its exact effective cap.")
        if verdict != "certified":
            raise UnresolvedMotionEventError(f"operation {operation_index} has an unresolved exact event partition.")
        return MotionWitness(
            operation_index,
            operation_kind,
            motion,
            user_cap.chord_ratio_bytes,
            effective_cap.chord_ratio_bytes,
            trace.oracle_strategy_version.encode(),
            self.stock_lineage_digest,
            trace.canonical_digest,
            "certified",
            trace.event_cell_count,
            0,
        )
