import hashlib
import struct
from dataclasses import dataclass
from typing import Final
from typing import Literal
from typing import cast

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
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
SWEPT_PREFIX_STRATEGY_VERSION: Final[bytes] = _continuous_tea_2.swept_prefix_strategy_version()
SWEPT_PREFIX_THEOREM_VERSION: Final[bytes] = _continuous_tea_2.swept_prefix_theorem_version()
SWEPT_PREFIX_MOTION_STRATA: Final[int] = 2


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


def _validate_cap_surrogates(
    user_cap_bytes: object,
    effective_cap_bytes: object,
    witness_kind: str,
) -> None:
    if type(user_cap_bytes) is not bytes or len(user_cap_bytes) != _BINARY64_SIZE:
        raise InvalidMotionCertificateError(f"{witness_kind} cap surrogates must be exact binary64 bytes.")
    if type(effective_cap_bytes) is not bytes or len(effective_cap_bytes) != _BINARY64_SIZE:
        raise InvalidMotionCertificateError(f"{witness_kind} cap surrogates must be exact binary64 bytes.")
    user_cap = struct.unpack(">d", user_cap_bytes)[0]
    effective_cap = struct.unpack(">d", effective_cap_bytes)[0]
    try:
        effective_within_user_cap = _stock_2.cap_chord_ratio_le(
            effective_cap,
            user_cap,
        )
    except ValueError:
        raise InvalidMotionCertificateError(f"{witness_kind} cap surrogates are outside the exact native domain.") from None
    if not effective_within_user_cap:
        raise InvalidMotionCertificateError(f"{witness_kind} effective cap exceeds its user cap.")


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


def _validate_swept_prefix_audit(
    *,
    audit: object,
    expected_stock_boundary_digest: bytes,
    expected_source_bytes: bytes,
) -> None:
    if type(audit) is not _continuous_tea_2.SweptPrefixSegmentTeaAudit2:
        raise InvalidMotionCertificateError("native swept-prefix oracle did not return one sealed audit record.")
    if audit.exact_verdict not in {"certified", "unresolved"}:
        raise InvalidMotionCertificateError("native swept-prefix oracle returned an unknown verdict.")
    if not audit.is_self_consistent:
        raise InvalidMotionCertificateError("native swept-prefix audit failed exact self-verification.")
    if type(audit.canonical_bytes) is not bytes or not audit.canonical_bytes:
        raise InvalidMotionCertificateError("native swept-prefix oracle returned invalid proof bytes.")
    if type(audit.canonical_digest) is not bytes or len(audit.canonical_digest) != _DIGEST_SIZE or hashlib.sha256(audit.canonical_bytes).digest() != audit.canonical_digest:
        raise InvalidMotionCertificateError("native swept-prefix proof digest contradicts its bytes.")
    if audit.strategy_version != SWEPT_PREFIX_STRATEGY_VERSION:
        raise InvalidMotionCertificateError("native swept-prefix audit uses a foreign strategy identity.")
    if audit.theorem_version != SWEPT_PREFIX_THEOREM_VERSION:
        raise InvalidMotionCertificateError("native swept-prefix audit uses a foreign theorem identity.")
    if audit.source_canonical_bytes != expected_source_bytes:
        raise InvalidMotionCertificateError("native swept-prefix proof contradicts its exact motion source.")
    if type(audit.stock_boundary_digest) is not bytes or len(audit.stock_boundary_digest) != _DIGEST_SIZE or audit.stock_boundary_digest != expected_stock_boundary_digest:
        raise InvalidMotionCertificateError("native swept-prefix proof contradicts its exact stock boundary.")
    if type(audit.motion_stratum_count) is not int or audit.motion_stratum_count != SWEPT_PREFIX_MOTION_STRATA:
        raise InvalidMotionCertificateError("native swept-prefix proof contradicts its exact motion strata.")


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
        _validate_cap_surrogates(
            self.user_cap_bytes,
            self.effective_cap_bytes,
            "motion witness",
        )
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
class SweptPrefixMotionWitness:
    """Certified advancing-cut witness carrying one sealed native theorem."""

    operation_index: int
    motion: ExactSegmentMotion
    tool_radius: ToolRadius
    user_cap_bytes: bytes
    effective_cap_bytes: bytes
    stock_lineage_digest: bytes
    stock_boundary_digest: bytes
    native_audit: _continuous_tea_2.SweptPrefixSegmentTeaAudit2

    def __post_init__(self) -> None:
        if type(self) is not SweptPrefixMotionWitness:
            raise InvalidMotionCertificateError("swept-prefix witness must use the exact owned type.")
        if type(self.operation_index) is not int or self.operation_index < 0:
            raise InvalidMotionCertificateError("swept-prefix witness operation index must be nonnegative.")
        if type(self.motion) is not ExactSegmentMotion or type(self.tool_radius) is not ToolRadius:
            raise InvalidMotionCertificateError("swept-prefix witness requires exact segment geometry and tool radius.")
        _validate_cap_surrogates(
            self.user_cap_bytes,
            self.effective_cap_bytes,
            "swept-prefix witness",
        )
        if any(type(value) is not bytes or len(value) != _DIGEST_SIZE for value in (self.stock_lineage_digest, self.stock_boundary_digest)):
            raise InvalidMotionCertificateError("swept-prefix witness state identities must be exact SHA-256 digests.")
        try:
            source = _continuous_tea_2.SegmentEventSource2.from_binary64(
                self.motion.start.x,
                self.motion.start.y,
                self.motion.end.x,
                self.motion.end.y,
                self.tool_radius.value,
                struct.unpack(">d", self.effective_cap_bytes)[0],
            )
        except (
            _continuous_tea_2.InvalidCapChordRatioError,
            _continuous_tea_2.NonFiniteSegmentInputError,
            _continuous_tea_2.NonPositiveToolRadiusError,
            _continuous_tea_2.ZeroLengthSegmentMotionError,
        ) as error:
            raise InvalidMotionCertificateError(
                f"swept-prefix witness violates its exact native source ({type(error).__name__}).",
            ) from error
        _validate_swept_prefix_audit(
            audit=self.native_audit,
            expected_stock_boundary_digest=self.stock_boundary_digest,
            expected_source_bytes=source.canonical_bytes,
        )
        if self.native_audit.exact_verdict != "certified":
            raise InvalidMotionCertificateError("swept-prefix witness requires one certified native verdict.")

    @property
    def operation_kind(self) -> OperationType:
        """Return the link operation kind owned by an advancing segment."""
        return OperationType.LINK

    @property
    def strategy_identity(self) -> bytes:
        """Return the exact native strategy identity."""
        return self.native_audit.strategy_version

    @property
    def theorem_identity(self) -> bytes:
        """Return the exact geometric theorem identity."""
        return self.native_audit.theorem_version

    @property
    def event_trace_digest(self) -> bytes:
        """Return the sealed native audit digest."""
        return self.native_audit.canonical_digest

    @property
    def event_cell_count(self) -> int:
        """Return the exact start/open-translation stratum count."""
        return self.native_audit.motion_stratum_count

    @property
    def verdict(self) -> Literal["certified"]:
        """Return the only constructible witness verdict."""
        return "certified"

    @property
    def unresolved_count(self) -> int:
        """Return zero because unresolved audits cannot construct witnesses."""
        return 0

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete advancing-cut proof record."""
        return encode_tagged_union(
            b"swept-prefix-motion-witness-v1",
            encode_component_map(
                {
                    b"effective-cap": encode_bytes(self.effective_cap_bytes),
                    b"motion": canonical_task1_bytes(self.motion),
                    b"native-audit": encode_bytes(self.native_audit.canonical_bytes),
                    b"operation-index": encode_integer(self.operation_index),
                    b"stock-boundary-digest": encode_bytes(self.stock_boundary_digest),
                    b"stock-lineage-digest": encode_bytes(self.stock_lineage_digest),
                    b"strategy-identity": encode_bytes(self.strategy_identity),
                    b"theorem-identity": encode_bytes(self.theorem_identity),
                    b"tool-radius": encode_tagged_union(
                        b"tool-radius-mm-v1",
                        encode_binary64(float(self.tool_radius.value)),
                    ),
                    b"user-cap": encode_bytes(self.user_cap_bytes),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
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
        try:
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
        except (
            _continuous_tea_2.EventPartitionVerificationError,
            _continuous_tea_2.EventTraceVerificationError,
            _continuous_tea_2.InvalidCapChordRatioError,
            _continuous_tea_2.NonFiniteSegmentInputError,
            _continuous_tea_2.NonPositiveToolRadiusError,
            _continuous_tea_2.ZeroLengthSegmentMotionError,
        ) as error:
            raise InvalidMotionCertificateError(
                f"operation {operation_index} violates the native event contract ({type(error).__name__}).",
            ) from error
        except (
            _continuous_tea_2.BoundaryExtractionError,
            _continuous_tea_2.EventSubstrateError,
        ) as error:
            raise UnresolvedMotionEventError(
                f"operation {operation_index} has an unresolved exact event substrate ({type(error).__name__}).",
            ) from error
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

    def certify_swept_prefix_segment(
        self,
        *,
        operation_index: int,
        motion: ExactSegmentMotion,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> SweptPrefixMotionWitness:
        """Certify an advancing segment from its already-cleared sweep prefix.

        This theorem-backed path is intentionally separate from `certify`.
        It proves that a clear start disk and the cutter's own translation
        prefix restrict possible engagement to the closed forward semicircle.

        Args:
            operation_index: Canonical nonnegative operation ordinal.
            motion: Exact nonzero advancing segment.
            user_cap: User engagement limit.
            effective_cap: Exact policy-derived cap, no greater than
                `user_cap`; the native theorem currently admits only pi.

        Returns:
            Immutable witness binding the native swept-prefix proof and stock
            pre-state.

        Raises:
            InvalidMotionCertificateError: If inputs, proof bytes, or stock
                identity violate their exact structural contract.
            UnresolvedMotionEventError: If start clearance or the exact pi-cap
                theorem is unavailable.
        """
        if type(operation_index) is not int or operation_index < 0:
            raise InvalidMotionCertificateError("operation index must be a nonnegative exact integer.")
        if type(motion) is not ExactSegmentMotion:
            raise InvalidMotionCertificateError("swept-prefix certification requires an exact segment motion.")
        if type(user_cap) is not EngagementCap or type(effective_cap) is not EngagementCap:
            raise InvalidMotionCertificateError("swept-prefix certification requires exact engagement caps.")
        if not _stock_2.cap_chord_ratio_le(effective_cap.chord_ratio, user_cap.chord_ratio):
            raise InvalidMotionCertificateError("effective cap exceeds the user cap.")
        try:
            source = _continuous_tea_2.SegmentEventSource2.from_binary64(
                motion.start.x,
                motion.start.y,
                motion.end.x,
                motion.end.y,
                self.tool_radius.value,
                effective_cap.chord_ratio,
            )
            audit = _continuous_tea_2.audit_swept_prefix_segment_tea_exact(
                self._stock,
                motion.start.x,
                motion.start.y,
                motion.end.x,
                motion.end.y,
                self.tool_radius.value,
                effective_cap.chord_ratio,
            )
        except (
            _continuous_tea_2.InvalidCapChordRatioError,
            _continuous_tea_2.NonFiniteSegmentInputError,
            _continuous_tea_2.NonPositiveToolRadiusError,
            _continuous_tea_2.ZeroLengthSegmentMotionError,
        ) as error:
            raise InvalidMotionCertificateError(
                f"operation {operation_index} violates the native swept-prefix contract ({type(error).__name__}).",
            ) from error
        except (
            _continuous_tea_2.BoundaryExtractionError,
            _continuous_tea_2.EventSubstrateError,
        ) as error:
            raise UnresolvedMotionEventError(
                f"operation {operation_index} has an unresolved swept-prefix theorem ({type(error).__name__}).",
            ) from error
        _validate_swept_prefix_audit(
            audit=audit,
            expected_stock_boundary_digest=self.canonical_boundary_digest,
            expected_source_bytes=source.canonical_bytes,
        )
        if audit.exact_verdict != "certified":
            raise UnresolvedMotionEventError(
                f"operation {operation_index} has unresolved swept-prefix preconditions.",
            )
        return SweptPrefixMotionWitness(
            operation_index,
            motion,
            self.tool_radius,
            user_cap.chord_ratio_bytes,
            effective_cap.chord_ratio_bytes,
            self.stock_lineage_digest,
            self.canonical_boundary_digest,
            audit,
        )
