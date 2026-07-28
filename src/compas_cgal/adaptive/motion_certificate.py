import hashlib
import struct
from dataclasses import dataclass
from typing import Final
from typing import Literal
from typing import cast

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import InvalidMotionCertificateError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.toolpath import OperationType

_DIGEST_SIZE = hashlib.sha256().digest_size
_BINARY64_SIZE = struct.calcsize(">d")
MOTION_CERTIFICATE_SCHEMA_VERSION: Final[bytes] = b"motion-certificate-schema-v1"


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


def _validate_native_trace(
    verdict: str,
    trace: _continuous_tea_2.EventTrace2,
) -> None:
    if type(verdict) is not str or verdict not in {"certified", "cap_exceeded", "unresolved"}:
        raise InvalidMotionCertificateError("native motion oracle returned an unknown verdict.")
    if type(trace) is not _continuous_tea_2.EventTrace2:
        raise InvalidMotionCertificateError("native motion oracle returned an invalid event trace.")
    if trace.exact_verdict != verdict:
        raise InvalidMotionCertificateError("native motion oracle verdict contradicts its event trace.")
    if type(trace.canonical_bytes) is not bytes or type(trace.canonical_digest) is not bytes or hashlib.sha256(trace.canonical_bytes).digest() != trace.canonical_digest:
        raise InvalidMotionCertificateError("native event-trace digest contradicts its canonical bytes.")
    if (
        type(trace.decision_authority_bytes) is not bytes
        or type(trace.decision_authority_digest) is not bytes
        or hashlib.sha256(trace.decision_authority_bytes).digest() != trace.decision_authority_digest
        or trace.decision_authority_digest not in trace.canonical_bytes
    ):
        raise InvalidMotionCertificateError("native event trace does not bind its deciding authority.")
    if type(trace.oracle_strategy_version) is not str or not trace.oracle_strategy_version:
        raise InvalidMotionCertificateError("native event trace omits its strategy identity.")


@dataclass(frozen=True)
class MotionWitness:
    """Immutable proof binding for one certified lateral motion.

    Attributes:
        operation_index: Canonical operation ordinal.
        operation_kind: Semantic lateral operation kind.
        motion: Exact segment or full-circle motion.
        user_cap_bytes: Canonical binary64 user-cap surrogate.
        effective_cap_bytes: Canonical binary64 effective-cap surrogate.
        strategy_identity: Native event-oracle strategy identifier.
        stock_lineage_digest: SHA-256 digest of the observed stock lineage.
        event_trace_digest: SHA-256 digest of the verified native event trace.
        verdict: Exact certified verdict.
        event_cell_count: Number of sign-invariant event cells.
        unresolved_count: Zero for every constructible witness.
    """

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
        if type(self) is not MotionWitness:
            raise InvalidMotionCertificateError("motion witness must use the exact owned type.")
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidMotionCertificateError("motion witness operation index must be nonnegative.")
        if type(self.operation_kind) is not OperationType:
            raise InvalidMotionCertificateError("motion witness operation kind must be exact OperationType.")
        if self.operation_kind not in (OperationType.CUT, OperationType.LINK):
            raise InvalidMotionCertificateError("motion witness is restricted to lateral operations.")
        if type(self.motion) not in (ExactSegmentMotion, ExactCircleMotion):
            raise InvalidMotionCertificateError("motion witness requires one exact motion.")
        if (self.operation_kind is OperationType.LINK and type(self.motion) is not ExactSegmentMotion) or (
            self.operation_kind is OperationType.CUT and type(self.motion) is not ExactCircleMotion
        ):
            raise InvalidMotionCertificateError("motion witness operation kind and motion are incompatible.")
        if any(type(value) is not bytes or len(value) != _BINARY64_SIZE for value in (self.user_cap_bytes, self.effective_cap_bytes)):
            raise InvalidMotionCertificateError("motion witness cap surrogates must be exact binary64 bytes.")
        user_cap = struct.unpack(">d", self.user_cap_bytes)[0]
        effective_cap = struct.unpack(">d", self.effective_cap_bytes)[0]
        try:
            effective_within_user_cap = _stock_2.cap_chord_ratio_le(
                effective_cap,
                user_cap,
            )
        except ValueError:
            raise InvalidMotionCertificateError("motion witness cap surrogates are outside the exact native domain.") from None
        if not effective_within_user_cap:
            raise InvalidMotionCertificateError("motion witness effective cap exceeds its user cap.")
        if any(type(value) is not bytes or len(value) != _DIGEST_SIZE for value in (self.stock_lineage_digest, self.event_trace_digest)):
            raise InvalidMotionCertificateError("motion witness digests must be exact SHA-256 bytes.")
        if type(self.strategy_identity) is not bytes or not self.strategy_identity:
            raise InvalidMotionCertificateError("motion witness requires a native strategy identity.")
        if self.verdict != "certified" or self.unresolved_count != 0:
            raise InvalidMotionCertificateError("motion witness may represent only resolved certification.")
        if type(self.event_cell_count) is not int or self.event_cell_count < 0:
            raise InvalidMotionCertificateError("motion witness event-cell count must be nonnegative.")

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete versioned proof record.

        Returns:
            Canonical CCAN bytes binding every witness field.
        """
        if self.operation_kind is OperationType.CUT:
            operation_tag = b"operation-cut-v1"
        elif self.operation_kind is OperationType.LINK:
            operation_tag = b"operation-link-v1"
        else:
            raise InvalidMotionCertificateError("motion witness operation kind has no canonical lateral tag.")
        return encode_tagged_union(
            b"motion-witness-v1",
            encode_component_map(
                {
                    b"effective-cap": encode_bytes(self.effective_cap_bytes),
                    b"event-cell-count": encode_integer(self.event_cell_count),
                    b"event-trace-digest": encode_bytes(self.event_trace_digest),
                    b"motion": canonical_task1_bytes(self.motion),
                    b"operation-index": encode_integer(self.operation_index),
                    b"operation-kind": encode_tagged_union(operation_tag, b""),
                    b"stock-lineage-digest": encode_bytes(self.stock_lineage_digest),
                    b"strategy-identity": encode_bytes(self.strategy_identity),
                    b"unresolved-count": encode_integer(self.unresolved_count),
                    b"user-cap": encode_bytes(self.user_cap_bytes),
                    b"verdict": encode_tagged_union(
                        b"motion-verdict-certified-v1",
                        b"",
                    ),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`.

        Returns:
            Immutable witness identity digest.
        """
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class MotionCertifier:
    """Read-only exact event certifier over one immutable stock snapshot.

    The certifier owns a native stock copy plus digests for the observed
    boundary and depletion lineage. Certification never mutates that state.
    """

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
        """Build a certifier from a typed stock owner.

        Args:
            stock: Exact stock-area owner whose current state is snapshotted.
            tool_radius: Typed cutter radius used by every motion audit.

        Returns:
            Read-only certifier with canonical boundary and lineage digests.

        Raises:
            InvalidMotionCertificateError: If either argument has the wrong
                exact type.
        """
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
        """Query the immutable native stock snapshot.

        Args:
            x: World-XY x coordinate in millimetres.
            y: World-XY y coordinate in millimetres.

        Returns:
            Whether the exact stock contains the queried point.
        """
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
        """Certify one exact lateral motion against its effective TEA cap.

        Args:
            operation_index: Canonical nonnegative operation ordinal.
            operation_kind: `CUT` for circles or `LINK` for segments.
            motion: Exact segment or full-circle motion.
            user_cap: User engagement limit.
            effective_cap: Exact policy-derived cap, no greater than
                `user_cap`.

        Returns:
            Immutable witness binding the native event trace and stock
            pre-state.

        Raises:
            InvalidMotionCertificateError: If local inputs or the native trace
                violate their structural contract.
            EngagementCapExceededError: If the exact oracle proves cap
                exceedance.
            UnresolvedMotionEventError: If the exact event proof is incomplete.
        """
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
            verdict, trace = _continuous_tea_2.audit_segment_tea_event_exact(
                self._stock,
                motion.start.x,
                motion.start.y,
                motion.end.x,
                motion.end.y,
                self.tool_radius.value,
                effective_cap.chord_ratio,
            )
        else:
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
        _validate_native_trace(verdict, trace)
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
