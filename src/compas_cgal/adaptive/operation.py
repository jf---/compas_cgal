import hashlib
import math
import struct
from dataclasses import dataclass
from enum import Enum
from typing import NewType
from typing import Self
from typing import TypeAlias

from compas_cgal import _stock_2
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
from compas_cgal.adaptive.identity import OPERATION_SCHEMA_VERSION
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

EdgeId = NewType("EdgeId", bytes)
CursorIdentity = NewType("CursorIdentity", bytes)
NeckOwnerId = NewType("NeckOwnerId", bytes)

_BINARY64 = struct.Struct(">d")


def _identity_bytes(value: object, name: str, *, digest: bool = False) -> bytes:
    if type(value) is not bytes or not value:
        raise ArtifactIdentityError(f"{name} must be nonempty identity bytes.")
    if digest and len(value) != hashlib.sha256().digest_size:
        raise ArtifactIdentityError(f"{name} must be exactly one 32-byte SHA-256 digest.")
    return value


def _cap_surrogate_bytes(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _BINARY64.size:
        raise ArtifactIdentityError(f"{name} must be exactly one big-endian binary64 surrogate.")
    surrogate = _BINARY64.unpack(value)[0]
    if not math.isfinite(surrogate) or not 0.0 < surrogate <= 4.0:
        raise ArtifactIdentityError(f"{name} must lie in the exact engagement surrogate domain (0, 4].")
    return value


def _cap_bytes(cap: EngagementCap, name: str) -> bytes:
    if not isinstance(cap, EngagementCap):
        raise ArtifactIdentityError(f"{name} must come from the native cap boundary.")
    return _cap_surrogate_bytes(cap.chord_ratio_bytes, name)


@dataclass(frozen=True)
class WidthClassId:
    value: int

    def __post_init__(self) -> None:
        if type(self.value) is not int or self.value < 0:
            raise ArtifactIdentityError("width-class ID must be an exact nonnegative integer.")

    @classmethod
    def build(cls, value: int) -> Self:
        return cls(value)

    @property
    def canonical_bytes(self) -> bytes:
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
        return encode_tagged_union(b"no-neck-scope-v1", b"")


@dataclass(frozen=True)
class OrientedNeckScope:
    neck_owner_id: NeckOwnerId
    orientation: NeckTraversalOrientation

    def __post_init__(self) -> None:
        _identity_bytes(self.neck_owner_id, "neck owner ID")
        if not isinstance(self.orientation, NeckTraversalOrientation):
            raise ArtifactIdentityError("neck traversal orientation must be explicit.")

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
            raise ArtifactIdentityError("full-cap decision requires equal user and effective surrogate bytes.")

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
        _identity_bytes(self.neck_evidence_digest, "neck evidence digest", digest=True)
        if not isinstance(self.width_class_id, WidthClassId):
            raise ArtifactIdentityError("neck-cap decision requires a typed width-class ID.")
        if not isinstance(self.passage_before, PassageState) or not isinstance(self.passage_after, PassageState):
            raise ArtifactIdentityError("neck-cap decision requires typed passage states.")
        if (self.passage_before, self.passage_after) not in _PASSAGE_ADVANCES:
            raise ArtifactIdentityError("neck-cap passage state must advance exactly one canonical step.")
        _cap_surrogate_bytes(self.user_cap_bytes, "user cap")
        _cap_surrogate_bytes(self.effective_cap_bytes, "effective cap")
        user_cap = _BINARY64.unpack(self.user_cap_bytes)[0]
        effective_cap = _BINARY64.unpack(self.effective_cap_bytes)[0]
        if not _stock_2.cap_chord_ratio_le(effective_cap, user_cap):
            raise ArtifactIdentityError("effective neck cap exceeds the user cap.")

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
    if isinstance(neck_scope, NoNeckScope) and not isinstance(effective_cap_decision, FullCapDecision):
        raise ArtifactIdentityError("no-neck scope requires a full-cap decision.")
    if isinstance(neck_scope, OrientedNeckScope) and not isinstance(effective_cap_decision, NeckCapDecision):
        raise ArtifactIdentityError("oriented neck scope requires a neck-evidence cap decision.")


def _validate_traversal_ids(
    component_id: ComponentId,
    edge_id: EdgeId,
    branch_id: BranchId,
    cursor_before: CursorIdentity,
    cursor_after: CursorIdentity,
) -> None:
    _identity_bytes(component_id, "component ID")
    _identity_bytes(edge_id, "edge ID")
    _identity_bytes(branch_id, "branch ID")
    _identity_bytes(cursor_before, "cursor-before identity")
    _identity_bytes(cursor_after, "cursor-after identity")


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
            raise ArtifactIdentityError("hold traversal must preserve the exact cursor identity.")
        if self.makes_cursor_terminal is not False:
            raise ArtifactIdentityError("hold traversal cannot make the cursor terminal.")

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
            raise ArtifactIdentityError("advance traversal requires different exact cursor identities.")
        if type(self.makes_cursor_terminal) is not bool:
            raise ArtifactIdentityError("advance traversal terminal state must be explicitly boolean.")

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
    cut_plane: CutPlane

    def __post_init__(self) -> None:
        if not isinstance(self.position, Point2) or not isinstance(self.cut_plane, CutPlane):
            raise ArtifactIdentityError("approach requires a typed world-XY position and complete cut plane.")

    @classmethod
    def build(cls, *, position: Point2[WorldXY], cut_plane: CutPlane) -> Self:
        return cls(position, cut_plane)

    @property
    def canonical_bytes(self) -> bytes:
        return _operation_bytes(
            b"approach-v1",
            {
                b"cut-plane": canonical_task1_bytes(self.cut_plane),
                b"position": canonical_point2_bytes(self.position),
            },
        )


@dataclass(frozen=True)
class PlungeOperation:
    position: Point2[WorldXY]
    cut_plane: CutPlane

    def __post_init__(self) -> None:
        if not isinstance(self.position, Point2) or not isinstance(self.cut_plane, CutPlane):
            raise ArtifactIdentityError("plunge requires a typed world-XY position and complete cut plane.")

    @classmethod
    def build(cls, *, position: Point2[WorldXY], cut_plane: CutPlane) -> Self:
        return cls(position, cut_plane)

    @property
    def canonical_bytes(self) -> bytes:
        return _operation_bytes(
            b"plunge-v1",
            {
                b"cut-plane": canonical_task1_bytes(self.cut_plane),
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
        if not isinstance(self.motion, ExactSegmentMotion) or not isinstance(self.cut_z, CutZ):
            raise ArtifactIdentityError("link segment requires typed exact motion and cut Z.")
        if not isinstance(self.neck_scope, (NoNeckScope, OrientedNeckScope)):
            raise ArtifactIdentityError("link segment requires a closed neck scope.")
        if not isinstance(self.effective_cap_decision, (FullCapDecision, NeckCapDecision)):
            raise ArtifactIdentityError("link segment requires a closed effective-cap decision.")
        _validate_scope_cap_decision(self.neck_scope, self.effective_cap_decision)
        if not isinstance(self.traversal_decision, HoldTraversalDecision):
            raise ArtifactIdentityError("link segment must hold traversal; the accepted circle owns the advance.")

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
        if not isinstance(self.motion, ExactCircleMotion) or not isinstance(self.cut_z, CutZ):
            raise ArtifactIdentityError("cut full circle requires typed exact motion and cut Z.")
        if not isinstance(self.neck_scope, (NoNeckScope, OrientedNeckScope)):
            raise ArtifactIdentityError("cut full circle requires a closed neck scope.")
        if not isinstance(self.effective_cap_decision, (FullCapDecision, NeckCapDecision)):
            raise ArtifactIdentityError("cut full circle requires a closed effective-cap decision.")
        _validate_scope_cap_decision(self.neck_scope, self.effective_cap_decision)
        if not isinstance(self.traversal_decision, AdvanceTraversalDecision):
            raise ArtifactIdentityError("cut full circle must own the one accepted cursor advance.")

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
