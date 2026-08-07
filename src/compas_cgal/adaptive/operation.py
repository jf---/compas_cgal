import hashlib
import math
import struct
from dataclasses import dataclass
from enum import Enum
from typing import NewType
from typing import Self
from typing import TypeAlias

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import canonical_clearance_z_bytes
from compas_cgal.adaptive.canonical import canonical_cut_z_bytes
from compas_cgal.adaptive.canonical import canonical_point2_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_material_side
from compas_cgal.adaptive.canonical import encode_passage_state
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import ArtifactIdentityError
from compas_cgal.adaptive.errors import CanonicalEncodingError
from compas_cgal.adaptive.errors import DegenerateSegmentMotionError
from compas_cgal.adaptive.errors import InvalidAdvanceSegmentOperationError
from compas_cgal.adaptive.errors import InvalidEffectiveCapDecisionError
from compas_cgal.adaptive.errors import InvalidOperationIdentityError
from compas_cgal.adaptive.errors import InvalidRetraceSegmentOperationError
from compas_cgal.adaptive.errors import InvalidRouteRetraceDecisionError
from compas_cgal.adaptive.errors import InvalidTraversalDecisionError
from compas_cgal.adaptive.identity import OPERATION_SCHEMA_VERSION
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY

EdgeId = NewType("EdgeId", bytes)
CursorIdentity = NewType("CursorIdentity", bytes)
NeckOwnerId = NewType("NeckOwnerId", bytes)
RouteNodeId = NewType("RouteNodeId", bytes)

ROUTE_RETRACE_STRATEGY_VERSION = b"reverse-final-zero-guide-v1"

_BINARY64 = struct.Struct(">d")


def _identity_bytes(
    value: object,
    name: str,
    error_type: type[ArtifactIdentityError],
    *,
    digest: bool = False,
) -> bytes:
    if type(value) is not bytes or not value:
        raise error_type(f"{name} must be nonempty identity bytes.")
    if digest and len(value) != hashlib.sha256().digest_size:
        raise error_type(f"{name} must be exactly one 32-byte SHA-256 digest.")
    return value


def _cap_surrogate_bytes(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _BINARY64.size:
        raise InvalidEffectiveCapDecisionError(f"{name} must be exactly one big-endian binary64 surrogate.")
    surrogate = _BINARY64.unpack(value)[0]
    if not math.isfinite(surrogate) or not 0.0 < surrogate <= 4.0:
        raise InvalidEffectiveCapDecisionError(f"{name} must lie in the exact engagement surrogate domain (0, 4].")
    return value


def _cap_bytes(cap: EngagementCap, name: str) -> bytes:
    if type(cap) is not EngagementCap:
        raise InvalidEffectiveCapDecisionError(f"{name} must be exact EngagementCap from the native boundary.")
    return _cap_surrogate_bytes(cap.chord_ratio_bytes, name)


@dataclass(frozen=True)
class WidthClassId:
    value: int

    def __post_init__(self) -> None:
        if type(self.value) is not int or self.value < 0:
            raise InvalidEffectiveCapDecisionError("width-class ID must be an exact nonnegative integer.")

    @classmethod
    def build(cls, value: int) -> Self:
        return cls(value)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not WidthClassId:
            raise InvalidEffectiveCapDecisionError("width-class ID must be exact WidthClassId, not a subclass.")
        return encode_tagged_union(b"width-class-id-v1", encode_integer(self.value))


class NeckTraversalOrientation(Enum):
    FORWARD = b"neck-forward-v1"
    REVERSE = b"neck-reverse-v1"


@dataclass(frozen=True)
class NoNeckScope:
    @classmethod
    def build(cls) -> Self:
        return cls()

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not NoNeckScope:
            raise InvalidOperationIdentityError("neck scope must be exact NoNeckScope, not a subclass.")
        return encode_tagged_union(b"no-neck-scope-v1", b"")


@dataclass(frozen=True)
class OrientedNeckScope:
    neck_owner_id: NeckOwnerId
    orientation: NeckTraversalOrientation

    def __post_init__(self) -> None:
        _identity_bytes(self.neck_owner_id, "neck owner ID", InvalidOperationIdentityError)
        if type(self.orientation) is not NeckTraversalOrientation:
            raise InvalidOperationIdentityError("neck traversal orientation must be exact.")

    @classmethod
    def build(
        cls,
        *,
        neck_owner_id: NeckOwnerId,
        orientation: NeckTraversalOrientation,
    ) -> Self:
        return cls(neck_owner_id, orientation)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not OrientedNeckScope:
            raise InvalidOperationIdentityError("neck scope must be exact OrientedNeckScope, not a subclass.")
        return encode_tagged_union(
            b"oriented-neck-scope-v1",
            encode_component_map(
                {
                    b"neck-owner-id": bytes(self.neck_owner_id),
                    b"orientation": encode_tagged_union(bytes(self.orientation.value), b""),
                }
            ),
        )


NeckScope: TypeAlias = NoNeckScope | OrientedNeckScope


@dataclass(frozen=True)
class FullCapDecision:
    user_cap_bytes: bytes
    effective_cap_bytes: bytes

    def __post_init__(self) -> None:
        _cap_surrogate_bytes(self.user_cap_bytes, "user cap")
        _cap_surrogate_bytes(self.effective_cap_bytes, "effective cap")
        if self.user_cap_bytes != self.effective_cap_bytes:
            raise InvalidEffectiveCapDecisionError("full-cap decision requires equal user and effective surrogate bytes.")

    @classmethod
    def build(
        cls,
        *,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> Self:
        return cls(
            _cap_bytes(user_cap, "user cap"),
            _cap_bytes(effective_cap, "effective cap"),
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not FullCapDecision:
            raise InvalidEffectiveCapDecisionError("effective-cap decision must be exact FullCapDecision, not a subclass.")
        return encode_tagged_union(
            b"full-cap-decision-v1",
            encode_component_map(
                {
                    b"effective-cap": encode_bytes(self.effective_cap_bytes),
                    b"user-cap": encode_bytes(self.user_cap_bytes),
                }
            ),
        )


_PASSAGE_ADVANCES = {
    (PassageState.UNVISITED, PassageState.FIRST_PASS_COMPLETE),
    (PassageState.FIRST_PASS_COMPLETE, PassageState.SECOND_PASS_COMPLETE),
    (PassageState.SECOND_PASS_COMPLETE, PassageState.TERMINAL),
}


@dataclass(frozen=True)
class NeckCapDecision:
    neck_evidence_digest: IdentityDigest
    width_class_id: WidthClassId
    passage_before: PassageState
    passage_after: PassageState
    user_cap_bytes: bytes
    effective_cap_bytes: bytes

    def __post_init__(self) -> None:
        _identity_bytes(
            self.neck_evidence_digest,
            "neck evidence digest",
            InvalidEffectiveCapDecisionError,
            digest=True,
        )
        if type(self.width_class_id) is not WidthClassId:
            raise InvalidEffectiveCapDecisionError("neck-cap decision requires exact WidthClassId.")
        if type(self.passage_before) is not PassageState or type(self.passage_after) is not PassageState:
            raise InvalidEffectiveCapDecisionError("neck-cap decision requires typed passage states.")
        if (self.passage_before, self.passage_after) not in _PASSAGE_ADVANCES:
            raise InvalidEffectiveCapDecisionError("neck-cap passage state must advance exactly one canonical step.")
        _cap_surrogate_bytes(self.user_cap_bytes, "user cap")
        _cap_surrogate_bytes(self.effective_cap_bytes, "effective cap")
        user_cap = _BINARY64.unpack(self.user_cap_bytes)[0]
        effective_cap = _BINARY64.unpack(self.effective_cap_bytes)[0]
        if not _stock_2.cap_chord_ratio_le(effective_cap, user_cap):
            raise InvalidEffectiveCapDecisionError("effective neck cap exceeds the user cap.")

    @classmethod
    def build(
        cls,
        *,
        neck_evidence_digest: IdentityDigest,
        width_class_id: WidthClassId,
        passage_before: PassageState,
        passage_after: PassageState,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> Self:
        return cls(
            neck_evidence_digest,
            width_class_id,
            passage_before,
            passage_after,
            _cap_bytes(user_cap, "user cap"),
            _cap_bytes(effective_cap, "effective cap"),
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not NeckCapDecision:
            raise InvalidEffectiveCapDecisionError("effective-cap decision must be exact NeckCapDecision, not a subclass.")
        return encode_tagged_union(
            b"neck-cap-decision-v1",
            encode_component_map(
                {
                    b"effective-cap": encode_bytes(self.effective_cap_bytes),
                    b"neck-evidence-digest": bytes(self.neck_evidence_digest),
                    b"passage-after": encode_passage_state(self.passage_after),
                    b"passage-before": encode_passage_state(self.passage_before),
                    b"user-cap": encode_bytes(self.user_cap_bytes),
                    b"width-class-id": self.width_class_id.canonical_bytes,
                }
            ),
        )


EffectiveCapDecision: TypeAlias = FullCapDecision | NeckCapDecision


def _validate_scope_cap_decision(
    neck_scope: NeckScope,
    effective_cap_decision: EffectiveCapDecision,
    error_type: type[ArtifactIdentityError] = InvalidOperationIdentityError,
) -> None:
    if type(neck_scope) is NoNeckScope and type(effective_cap_decision) is not FullCapDecision:
        raise error_type("no-neck scope requires an exact full-cap decision.")
    if type(neck_scope) is OrientedNeckScope and type(effective_cap_decision) is not NeckCapDecision:
        raise error_type("oriented neck scope requires an exact neck-evidence cap decision.")


def _validate_traversal_ids(
    component_id: ComponentId,
    edge_id: EdgeId,
    branch_id: BranchId,
    cursor_before: CursorIdentity,
    cursor_after: CursorIdentity,
) -> None:
    _identity_bytes(component_id, "component ID", InvalidTraversalDecisionError)
    _identity_bytes(edge_id, "edge ID", InvalidTraversalDecisionError)
    _identity_bytes(branch_id, "branch ID", InvalidTraversalDecisionError)
    _identity_bytes(cursor_before, "cursor-before identity", InvalidTraversalDecisionError)
    _identity_bytes(cursor_after, "cursor-after identity", InvalidTraversalDecisionError)


def _traversal_bytes(
    *,
    tag: bytes,
    component_id: ComponentId,
    edge_id: EdgeId,
    branch_id: BranchId,
    cursor_before: CursorIdentity,
    cursor_after: CursorIdentity,
    makes_cursor_terminal: bool,
) -> bytes:
    return encode_tagged_union(
        tag,
        encode_component_map(
            {
                b"branch-id": bytes(branch_id),
                b"component-id": bytes(component_id),
                b"cursor-after": bytes(cursor_after),
                b"cursor-before": bytes(cursor_before),
                b"edge-id": bytes(edge_id),
                b"makes-cursor-terminal": encode_boolean(makes_cursor_terminal),
            }
        ),
    )


@dataclass(frozen=True)
class HoldTraversalDecision:
    component_id: ComponentId
    edge_id: EdgeId
    branch_id: BranchId
    cursor_before: CursorIdentity
    cursor_after: CursorIdentity
    makes_cursor_terminal: bool

    def __post_init__(self) -> None:
        _validate_traversal_ids(
            self.component_id,
            self.edge_id,
            self.branch_id,
            self.cursor_before,
            self.cursor_after,
        )
        if self.cursor_before != self.cursor_after:
            raise InvalidTraversalDecisionError("hold traversal must preserve the exact cursor identity.")
        if self.makes_cursor_terminal is not False:
            raise InvalidTraversalDecisionError("hold traversal cannot make the cursor terminal.")

    @classmethod
    def build(
        cls,
        *,
        component_id: ComponentId,
        edge_id: EdgeId,
        branch_id: BranchId,
        cursor: CursorIdentity,
    ) -> Self:
        return cls(component_id, edge_id, branch_id, cursor, cursor, False)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not HoldTraversalDecision:
            raise InvalidTraversalDecisionError("traversal decision must be exact HoldTraversalDecision, not a subclass.")
        return _traversal_bytes(
            tag=b"hold-traversal-v1",
            component_id=self.component_id,
            edge_id=self.edge_id,
            branch_id=self.branch_id,
            cursor_before=self.cursor_before,
            cursor_after=self.cursor_after,
            makes_cursor_terminal=self.makes_cursor_terminal,
        )


@dataclass(frozen=True)
class AdvanceTraversalDecision:
    component_id: ComponentId
    edge_id: EdgeId
    branch_id: BranchId
    cursor_before: CursorIdentity
    cursor_after: CursorIdentity
    makes_cursor_terminal: bool

    def __post_init__(self) -> None:
        _validate_traversal_ids(
            self.component_id,
            self.edge_id,
            self.branch_id,
            self.cursor_before,
            self.cursor_after,
        )
        if self.cursor_before == self.cursor_after:
            raise InvalidTraversalDecisionError("advance traversal requires different exact cursor identities.")
        if type(self.makes_cursor_terminal) is not bool:
            raise InvalidTraversalDecisionError("advance traversal terminal state must be explicitly boolean.")

    @classmethod
    def build(
        cls,
        *,
        component_id: ComponentId,
        edge_id: EdgeId,
        branch_id: BranchId,
        cursor_before: CursorIdentity,
        cursor_after: CursorIdentity,
        makes_cursor_terminal: bool,
    ) -> Self:
        return cls(
            component_id,
            edge_id,
            branch_id,
            cursor_before,
            cursor_after,
            makes_cursor_terminal,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not AdvanceTraversalDecision:
            raise InvalidTraversalDecisionError("traversal decision must be exact AdvanceTraversalDecision, not a subclass.")
        return _traversal_bytes(
            tag=b"advance-traversal-v1",
            component_id=self.component_id,
            edge_id=self.edge_id,
            branch_id=self.branch_id,
            cursor_before=self.cursor_before,
            cursor_after=self.cursor_after,
            makes_cursor_terminal=self.makes_cursor_terminal,
        )


TraversalDecision: TypeAlias = HoldTraversalDecision | AdvanceTraversalDecision


@dataclass(frozen=True)
class RouteRetraceDecision:
    """Authenticated causal boundary for one nonincident route switch."""

    completed_route_index: int
    activated_route_index: int
    completed_exit_node_id: RouteNodeId
    activated_entry_node_id: RouteNodeId
    terminal_traversal_digest: IdentityDigest
    activated_traversal_digest: IdentityDigest
    source_commit_digest: IdentityDigest
    source_transaction_digest: IdentityDigest
    source_operation_index: int
    source_operation_digest: IdentityDigest
    strategy_identity: bytes = ROUTE_RETRACE_STRATEGY_VERSION

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        if type(self) is not RouteRetraceDecision:
            raise InvalidRouteRetraceDecisionError(
                "route retrace decision must use the exact owned type.",
            )
        if type(self.completed_route_index) is not int or self.completed_route_index < 0:
            raise InvalidRouteRetraceDecisionError(
                "completed route index must be an exact nonnegative integer.",
            )
        if type(self.activated_route_index) is not int or self.activated_route_index < 0:
            raise InvalidRouteRetraceDecisionError(
                "activated route index must be an exact nonnegative integer.",
            )
        if self.activated_route_index != self.completed_route_index + 1:
            raise InvalidRouteRetraceDecisionError(
                "activated route index must immediately follow the completed route.",
            )
        completed_node = _identity_bytes(
            self.completed_exit_node_id,
            "completed exit node ID",
            InvalidRouteRetraceDecisionError,
        )
        activated_node = _identity_bytes(
            self.activated_entry_node_id,
            "activated entry node ID",
            InvalidRouteRetraceDecisionError,
        )
        if completed_node == activated_node:
            raise InvalidRouteRetraceDecisionError(
                "route retrace decision requires distinct nonincident boundary nodes.",
            )
        for digest, name in (
            (self.terminal_traversal_digest, "terminal traversal digest"),
            (self.activated_traversal_digest, "activated traversal digest"),
            (self.source_commit_digest, "source commit digest"),
            (self.source_transaction_digest, "source transaction digest"),
            (self.source_operation_digest, "source operation digest"),
        ):
            _identity_bytes(
                digest,
                name,
                InvalidRouteRetraceDecisionError,
                digest=True,
            )
        if type(self.source_operation_index) is not int or self.source_operation_index < 2:
            raise InvalidRouteRetraceDecisionError(
                "source operation index must be an exact integer of at least two.",
            )
        if type(self.strategy_identity) is not bytes or self.strategy_identity != ROUTE_RETRACE_STRATEGY_VERSION:
            raise InvalidRouteRetraceDecisionError(
                "route retrace decision requires the fixed exact strategy identity.",
            )

    @classmethod
    def build(
        cls,
        *,
        completed_route_index: int,
        activated_route_index: int,
        completed_exit_node_id: RouteNodeId,
        activated_entry_node_id: RouteNodeId,
        terminal_traversal_digest: IdentityDigest,
        activated_traversal_digest: IdentityDigest,
        source_commit_digest: IdentityDigest,
        source_transaction_digest: IdentityDigest,
        source_operation_index: int,
        source_operation_digest: IdentityDigest,
    ) -> Self:
        """Build one exact adjacent-route retrace decision.

        Args:
            completed_route_index: Terminal route ordinal.
            activated_route_index: Immediately following route ordinal.
            completed_exit_node_id: Exact completed-route endpoint identity.
            activated_entry_node_id: Exact activated-route entry identity.
            terminal_traversal_digest: Terminal global traversal identity.
            activated_traversal_digest: Activated global traversal identity.
            source_commit_digest: Accepted source commit identity.
            source_transaction_digest: Accepted source transaction identity.
            source_operation_index: Source ordinal in the physical programme.
            source_operation_digest: Canonical source operation identity.

        Returns:
            Validated content-addressed route-boundary decision.

        Raises:
            InvalidRouteRetraceDecisionError: If an identity is malformed or
                the route boundary and admitted source are cross-wired.
        """
        return cls(
            completed_route_index,
            activated_route_index,
            completed_exit_node_id,
            activated_entry_node_id,
            terminal_traversal_digest,
            activated_traversal_digest,
            source_commit_digest,
            source_transaction_digest,
            source_operation_index,
            source_operation_digest,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return CCAN bytes binding every causal identity field."""
        self._validate()
        return encode_tagged_union(
            b"route-retrace-decision-v1",
            encode_component_map(
                {
                    b"activated-entry-node-id": encode_bytes(
                        bytes(self.activated_entry_node_id),
                    ),
                    b"activated-route-index": encode_integer(
                        self.activated_route_index,
                    ),
                    b"activated-traversal-digest": encode_bytes(
                        bytes(self.activated_traversal_digest),
                    ),
                    b"completed-exit-node-id": encode_bytes(
                        bytes(self.completed_exit_node_id),
                    ),
                    b"completed-route-index": encode_integer(
                        self.completed_route_index,
                    ),
                    b"source-commit-digest": encode_bytes(
                        bytes(self.source_commit_digest),
                    ),
                    b"source-operation-digest": encode_bytes(
                        bytes(self.source_operation_digest),
                    ),
                    b"source-operation-index": encode_integer(
                        self.source_operation_index,
                    ),
                    b"source-transaction-digest": encode_bytes(
                        bytes(self.source_transaction_digest),
                    ),
                    b"strategy-identity": encode_bytes(self.strategy_identity),
                    b"terminal-traversal-digest": encode_bytes(
                        bytes(self.terminal_traversal_digest),
                    ),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
        return IdentityDigest(
            hashlib.sha256(self.canonical_bytes).digest(),
        )


def _operation_bytes(tag: bytes, components: dict[bytes, bytes]) -> bytes:
    return encode_tagged_union(
        OPERATION_SCHEMA_VERSION,
        encode_tagged_union(tag, encode_component_map(components)),
    )


@dataclass(frozen=True)
class ApproachOperation:
    position: Point2[WorldXY]
    clearance_z: ClearanceZ

    def __post_init__(self) -> None:
        if type(self.position) is not Point2 or type(self.clearance_z) is not ClearanceZ:
            raise InvalidOperationIdentityError("approach requires exact Point2[WorldXY] and ClearanceZ.")

    @classmethod
    def build(cls, *, position: Point2[WorldXY], clearance_z: ClearanceZ) -> Self:
        return cls(position, clearance_z)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not ApproachOperation:
            raise InvalidOperationIdentityError("operation must be exact ApproachOperation, not a subclass.")
        return _operation_bytes(
            b"approach-v1",
            {
                b"clearance-z": canonical_clearance_z_bytes(self.clearance_z),
                b"position": canonical_point2_bytes(self.position),
            },
        )


@dataclass(frozen=True)
class PlungeOperation:
    position: Point2[WorldXY]
    clearance_z: ClearanceZ
    cut_z: CutZ

    def __post_init__(self) -> None:
        if type(self.position) is not Point2 or type(self.clearance_z) is not ClearanceZ or type(self.cut_z) is not CutZ:
            raise InvalidOperationIdentityError("plunge requires exact Point2[WorldXY], ClearanceZ, and CutZ.")
        if self.clearance_z.value <= self.cut_z.value:
            raise InvalidOperationIdentityError("plunge clearance Z must be strictly above cut Z.")

    @classmethod
    def build(
        cls,
        *,
        position: Point2[WorldXY],
        clearance_z: ClearanceZ,
        cut_z: CutZ,
    ) -> Self:
        return cls(position, clearance_z, cut_z)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not PlungeOperation:
            raise InvalidOperationIdentityError("operation must be exact PlungeOperation, not a subclass.")
        return _operation_bytes(
            b"plunge-v1",
            {
                b"clearance-z": canonical_clearance_z_bytes(self.clearance_z),
                b"cut-z": canonical_cut_z_bytes(self.cut_z),
                b"position": canonical_point2_bytes(self.position),
            },
        )


@dataclass(frozen=True)
class LinkSegmentOperation:
    motion: ExactSegmentMotion
    cut_z: CutZ
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: HoldTraversalDecision

    def __post_init__(self) -> None:
        if type(self.motion) is not ExactSegmentMotion or type(self.cut_z) is not CutZ:
            raise InvalidOperationIdentityError("link segment requires exact ExactSegmentMotion and CutZ.")
        if type(self.motion.start) is not Point2 or type(self.motion.end) is not Point2:
            raise InvalidOperationIdentityError("link segment motion endpoints must be exact Point2[WorldXY].")
        if type(self.neck_scope) not in (NoNeckScope, OrientedNeckScope):
            raise InvalidOperationIdentityError("link segment requires an exact closed neck scope.")
        if type(self.effective_cap_decision) not in (FullCapDecision, NeckCapDecision):
            raise InvalidOperationIdentityError("link segment requires an exact effective-cap decision.")
        _validate_scope_cap_decision(self.neck_scope, self.effective_cap_decision)
        if type(self.traversal_decision) is not HoldTraversalDecision:
            raise InvalidOperationIdentityError("link segment must hold traversal using exact HoldTraversalDecision; the accepted circle owns the advance.")

    @classmethod
    def build(
        cls,
        *,
        motion: ExactSegmentMotion,
        cut_z: CutZ,
        neck_scope: NeckScope,
        effective_cap_decision: EffectiveCapDecision,
        traversal_decision: HoldTraversalDecision,
    ) -> Self:
        return cls(
            motion,
            cut_z,
            neck_scope,
            effective_cap_decision,
            traversal_decision,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not LinkSegmentOperation:
            raise InvalidOperationIdentityError("operation must be exact LinkSegmentOperation, not a subclass.")
        return _operation_bytes(
            b"link-segment-v1",
            {
                b"cut-z": canonical_cut_z_bytes(self.cut_z),
                b"effective-cap-decision": self.effective_cap_decision.canonical_bytes,
                b"motion": canonical_task1_bytes(self.motion),
                b"neck-scope": self.neck_scope.canonical_bytes,
                b"traversal-decision": self.traversal_decision.canonical_bytes,
            },
        )


@dataclass(frozen=True)
class AdvanceSegmentOperation:
    """One exact segment that owns a MAT traversal advance."""

    motion: ExactSegmentMotion
    cut_z: CutZ
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: AdvanceTraversalDecision

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        if type(self.motion) is not ExactSegmentMotion or type(self.cut_z) is not CutZ:
            raise InvalidAdvanceSegmentOperationError(
                "advancing segment requires exact ExactSegmentMotion and CutZ.",
            )
        if type(self.motion.start) is not Point2 or type(self.motion.end) is not Point2:
            raise InvalidAdvanceSegmentOperationError(
                "advancing segment motion endpoints must be exact Point2[WorldXY].",
            )
        if type(self.neck_scope) not in (NoNeckScope, OrientedNeckScope):
            raise InvalidAdvanceSegmentOperationError(
                "advancing segment requires an exact closed neck scope.",
            )
        if type(self.effective_cap_decision) not in (
            FullCapDecision,
            NeckCapDecision,
        ):
            raise InvalidAdvanceSegmentOperationError(
                "advancing segment requires an exact effective-cap decision.",
            )
        _validate_scope_cap_decision(
            self.neck_scope,
            self.effective_cap_decision,
            InvalidAdvanceSegmentOperationError,
        )
        if type(self.traversal_decision) is not AdvanceTraversalDecision:
            raise InvalidAdvanceSegmentOperationError(
                "advancing segment must own one exact advance traversal decision.",
            )
        try:
            self.motion.__post_init__()
            canonical_task1_bytes(self.motion)
            canonical_cut_z_bytes(self.cut_z)
            if type(self.neck_scope) is OrientedNeckScope:
                self.neck_scope.__post_init__()
            self.effective_cap_decision.__post_init__()
            if type(self.effective_cap_decision) is NeckCapDecision:
                self.effective_cap_decision.width_class_id.__post_init__()
            self.traversal_decision.__post_init__()
        except (
            CanonicalEncodingError,
            DegenerateSegmentMotionError,
            InvalidEffectiveCapDecisionError,
            InvalidOperationIdentityError,
            InvalidTraversalDecisionError,
        ) as error:
            raise InvalidAdvanceSegmentOperationError(
                "advancing segment contains malformed nested domain state.",
            ) from error

    @classmethod
    def build(
        cls,
        *,
        motion: ExactSegmentMotion,
        cut_z: CutZ,
        neck_scope: NeckScope,
        effective_cap_decision: EffectiveCapDecision,
        traversal_decision: AdvanceTraversalDecision,
    ) -> Self:
        """Build one exact link-only traversal operation.

        Args:
            motion: Nondegenerate exact centre-line segment.
            cut_z: Qualified cutting plane depth.
            neck_scope: Closed causal neck scope for this motion.
            effective_cap_decision: Cap decision owned by the scope.
            traversal_decision: The one cursor advance owned by the segment.

        Returns:
            Validated advancing segment operation.

        Raises:
            InvalidAdvanceSegmentOperationError: If any field is foreign or
                the scope, cap, and traversal contracts disagree.
        """
        return cls(
            motion,
            cut_z,
            neck_scope,
            effective_cap_decision,
            traversal_decision,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete canonical advancing-operation record.

        Returns:
            Versioned CCAN bytes binding motion and every owned decision.

        Raises:
            InvalidAdvanceSegmentOperationError: If invoked on a subclass.
        """
        if type(self) is not AdvanceSegmentOperation:
            raise InvalidAdvanceSegmentOperationError(
                "operation must be exact AdvanceSegmentOperation, not a subclass.",
            )
        return _operation_bytes(
            b"advance-segment-v1",
            {
                b"cut-z": canonical_cut_z_bytes(self.cut_z),
                b"effective-cap-decision": self.effective_cap_decision.canonical_bytes,
                b"motion": canonical_task1_bytes(self.motion),
                b"neck-scope": self.neck_scope.canonical_bytes,
                b"traversal-decision": self.traversal_decision.canonical_bytes,
            },
        )


@dataclass(frozen=True, init=False)
class RetraceSegmentOperation:
    """Exact source-derived reversal restoring an adjacent route entry."""

    motion: ExactSegmentMotion
    cut_z: CutZ
    effective_cap_decision: FullCapDecision
    decision: RouteRetraceDecision

    def __init__(self) -> None:
        raise InvalidRetraceSegmentOperationError(
            "retrace operation must be built from an admitted source.",
        )

    def _validate(self) -> None:
        if type(self) is not RetraceSegmentOperation:
            raise InvalidRetraceSegmentOperationError(
                "retrace operation must use the exact owned type.",
            )
        if type(self.motion) is not ExactSegmentMotion or type(self.cut_z) is not CutZ:
            raise InvalidRetraceSegmentOperationError(
                "retrace requires exact ExactSegmentMotion and CutZ.",
            )
        if type(self.motion.start) is not Point2 or type(self.motion.end) is not Point2:
            raise InvalidRetraceSegmentOperationError(
                "retrace motion endpoints must be exact Point2[WorldXY].",
            )
        if type(self.effective_cap_decision) is not FullCapDecision:
            raise InvalidRetraceSegmentOperationError(
                "retrace requires one exact full-cap decision.",
            )
        if type(self.decision) is not RouteRetraceDecision or self.decision.strategy_identity != ROUTE_RETRACE_STRATEGY_VERSION:
            raise InvalidRetraceSegmentOperationError(
                "retrace requires one exact fixed-strategy route decision.",
            )

    @classmethod
    def build(
        cls,
        *,
        source_operation: AdvanceSegmentOperation,
        decision: RouteRetraceDecision,
    ) -> Self:
        """Derive the only admitted reverse motion from its accepted source.

        Args:
            source_operation: Accepted route-terminal advancing segment.
            decision: Authenticated nonincident route boundary and source.

        Returns:
            Validated retrace with source depth and full-cap authority.

        Raises:
            InvalidRetraceSegmentOperationError: If the source is foreign,
                neck-owned, non-full-cap, nonterminal, or misidentified.
        """
        if type(source_operation) is not AdvanceSegmentOperation:
            raise InvalidRetraceSegmentOperationError(
                "retrace source must be one exact advancing segment.",
            )
        if type(decision) is not RouteRetraceDecision:
            raise InvalidRetraceSegmentOperationError(
                "retrace requires one exact route retrace decision.",
            )
        decision._validate()
        try:
            source_operation._validate()
        except InvalidAdvanceSegmentOperationError as error:
            raise InvalidRetraceSegmentOperationError(
                "retrace source must be one valid exact advancing segment.",
            ) from error
        if (
            type(source_operation.neck_scope) is not NoNeckScope
            or type(source_operation.effective_cap_decision) is not FullCapDecision
            or source_operation.traversal_decision.makes_cursor_terminal is not True
        ):
            raise InvalidRetraceSegmentOperationError(
                "retrace source must be no-neck, full-cap, and route-terminal.",
            )
        source_digest = IdentityDigest(
            hashlib.sha256(source_operation.canonical_bytes).digest(),
        )
        if decision.source_operation_digest != source_digest:
            raise InvalidRetraceSegmentOperationError(
                "retrace decision does not name its source operation.",
            )
        motion = ExactSegmentMotion.build(
            source_operation.motion.end,
            source_operation.motion.start,
        )
        operation = object.__new__(cls)
        object.__setattr__(operation, "motion", motion)
        object.__setattr__(operation, "cut_z", source_operation.cut_z)
        object.__setattr__(
            operation,
            "effective_cap_decision",
            source_operation.effective_cap_decision,
        )
        object.__setattr__(operation, "decision", decision)
        operation._validate()
        return operation

    @property
    def canonical_bytes(self) -> bytes:
        """Return the distinct source-derived retrace operation record."""
        self._validate()
        return _operation_bytes(
            b"retrace-segment-v1",
            {
                b"cut-z": canonical_cut_z_bytes(self.cut_z),
                b"decision": self.decision.canonical_bytes,
                b"effective-cap-decision": self.effective_cap_decision.canonical_bytes,
                b"motion": canonical_task1_bytes(self.motion),
            },
        )


@dataclass(frozen=True)
class CutFullCircleOperation:
    motion: ExactCircleMotion
    cut_z: CutZ
    material_side: MaterialSide
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: AdvanceTraversalDecision

    def __post_init__(self) -> None:
        if type(self.motion) is not ExactCircleMotion or type(self.cut_z) is not CutZ:
            raise InvalidOperationIdentityError("cut full circle requires exact ExactCircleMotion and CutZ.")
        if type(self.motion.center) is not Point2 or type(self.motion.phase_vector) is not Vector2:
            raise InvalidOperationIdentityError("cut full circle motion requires exact Point2[WorldXY] center and exact Vector2[WorldXY] phase.")
        if type(self.material_side) is not MaterialSide:
            raise InvalidOperationIdentityError("cut full circle requires the exact radial material side that caused its orientation.")
        if type(self.neck_scope) not in (NoNeckScope, OrientedNeckScope):
            raise InvalidOperationIdentityError("cut full circle requires an exact closed neck scope.")
        if type(self.effective_cap_decision) not in (FullCapDecision, NeckCapDecision):
            raise InvalidOperationIdentityError("cut full circle requires an exact effective-cap decision.")
        _validate_scope_cap_decision(self.neck_scope, self.effective_cap_decision)
        if type(self.traversal_decision) is not AdvanceTraversalDecision:
            raise InvalidOperationIdentityError("cut full circle must own the one accepted exact AdvanceTraversalDecision.")

    @classmethod
    def build(
        cls,
        *,
        motion: ExactCircleMotion,
        cut_z: CutZ,
        material_side: MaterialSide,
        neck_scope: NeckScope,
        effective_cap_decision: EffectiveCapDecision,
        traversal_decision: AdvanceTraversalDecision,
    ) -> Self:
        return cls(
            motion,
            cut_z,
            material_side,
            neck_scope,
            effective_cap_decision,
            traversal_decision,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not CutFullCircleOperation:
            raise InvalidOperationIdentityError("operation must be exact CutFullCircleOperation, not a subclass.")
        return _operation_bytes(
            b"cut-full-circle-v1",
            {
                b"cut-z": canonical_cut_z_bytes(self.cut_z),
                b"effective-cap-decision": self.effective_cap_decision.canonical_bytes,
                b"material-side": encode_material_side(self.material_side),
                b"motion": canonical_task1_bytes(self.motion),
                b"neck-scope": self.neck_scope.canonical_bytes,
                b"traversal-decision": self.traversal_decision.canonical_bytes,
            },
        )


CanonicalOperation: TypeAlias = ApproachOperation | PlungeOperation | LinkSegmentOperation | CutFullCircleOperation | AdvanceSegmentOperation | RetraceSegmentOperation
