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
from compas_cgal.adaptive.canonical import encode_passage_state
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import ArtifactIdentityError
from compas_cgal.adaptive.errors import InvalidEffectiveCapDecisionError
from compas_cgal.adaptive.errors import InvalidOperationIdentityError
from compas_cgal.adaptive.errors import InvalidTraversalDecisionError
from compas_cgal.adaptive.identity import OPERATION_SCHEMA_VERSION
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY

EdgeId = NewType("EdgeId", bytes)
CursorIdentity = NewType("CursorIdentity", bytes)
NeckOwnerId = NewType("NeckOwnerId", bytes)

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
) -> None:
    if type(neck_scope) is NoNeckScope and type(effective_cap_decision) is not FullCapDecision:
        raise InvalidOperationIdentityError("no-neck scope requires an exact full-cap decision.")
    if type(neck_scope) is OrientedNeckScope and type(effective_cap_decision) is not NeckCapDecision:
        raise InvalidOperationIdentityError("oriented neck scope requires an exact neck-evidence cap decision.")


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
class CutFullCircleOperation:
    motion: ExactCircleMotion
    cut_z: CutZ
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: AdvanceTraversalDecision

    def __post_init__(self) -> None:
        if type(self.motion) is not ExactCircleMotion or type(self.cut_z) is not CutZ:
            raise InvalidOperationIdentityError("cut full circle requires exact ExactCircleMotion and CutZ.")
        if type(self.motion.center) is not Point2 or type(self.motion.phase_vector) is not Vector2:
            raise InvalidOperationIdentityError("cut full circle motion requires exact Point2[WorldXY] center and exact Vector2[WorldXY] phase.")
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
        neck_scope: NeckScope,
        effective_cap_decision: EffectiveCapDecision,
        traversal_decision: AdvanceTraversalDecision,
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
        if type(self) is not CutFullCircleOperation:
            raise InvalidOperationIdentityError("operation must be exact CutFullCircleOperation, not a subclass.")
        return _operation_bytes(
            b"cut-full-circle-v1",
            {
                b"cut-z": canonical_cut_z_bytes(self.cut_z),
                b"effective-cap-decision": self.effective_cap_decision.canonical_bytes,
                b"motion": canonical_task1_bytes(self.motion),
                b"neck-scope": self.neck_scope.canonical_bytes,
                b"traversal-decision": self.traversal_decision.canonical_bytes,
            },
        )


CanonicalOperation: TypeAlias = ApproachOperation | PlungeOperation | LinkSegmentOperation | CutFullCircleOperation
