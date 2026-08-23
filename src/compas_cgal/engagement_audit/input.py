"""Immutable authoritative engagement-audit request."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Final
from typing import Self

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import canonical_clearance_z_bytes
from compas_cgal.adaptive.canonical import canonical_cut_z_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.engagement_audit.classification import classify_operation_snapshot
from compas_cgal.engagement_audit.errors import EmptyToolpathAuditError
from compas_cgal.engagement_audit.errors import InvalidAuditOperationError
from compas_cgal.engagement_audit.errors import InvalidEngagementAuditInputError
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.operation_identity import OperationStreamDigest
from compas_cgal.engagement_audit.operation_identity import canonical_operation_snapshot_bytes
from compas_cgal.engagement_audit.operation_identity import operation_stream_digest
from compas_cgal.engagement_audit.operation_identity import snapshot_toolpath_operation
from compas_cgal.engagement_audit.records import AuthenticatedOperation
from compas_cgal.toolpath import ToolpathOperation

AUDIT_INPUT_VERSION: Final[bytes] = b"engagement-audit-input-v1"


def _validated_holes(
    holes: object,
    *,
    require_canonical_order: bool,
) -> tuple[CanonicalRingV1, ...]:
    if type(holes) is not tuple or any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in holes):
        raise InvalidEngagementAuditInputError("holes must be an exact tuple of canonical hole rings.")
    ordered = tuple(sorted(holes, key=lambda hole: hole.canonical_bytes))
    if len({hole.canonical_bytes for hole in ordered}) != len(ordered):
        raise InvalidEngagementAuditInputError("canonical hole rings must be unique.")
    if require_canonical_order and holes != ordered:
        raise InvalidEngagementAuditInputError("raw hole rings must already use canonical order.")
    return ordered


def _validate_authoritative_fields(
    design_boundary: object,
    cut_plane: object,
    tool_radius: object,
    engagement_cap: object,
    build_identity: object,
) -> None:
    if type(design_boundary) is not CanonicalRingV1 or not design_boundary.is_outer:
        raise InvalidEngagementAuditInputError("design boundary must be one exact canonical outer ring.")
    if type(cut_plane) is not CutPlane:
        raise InvalidEngagementAuditInputError("cut plane must be exact CutPlane.")
    if type(tool_radius) is not ToolRadius:
        raise InvalidEngagementAuditInputError("tool radius must be exact ToolRadius.")
    if type(engagement_cap) is not EngagementCap:
        raise InvalidEngagementAuditInputError("engagement cap must be exact EngagementCap.")
    if type(build_identity) is not BuildIdentity:
        raise InvalidEngagementAuditInputError("build identity must be exact BuildIdentity.")


def _classify_operations(
    operations: object,
    cut_plane: CutPlane,
) -> tuple[tuple[AuthenticatedOperation, ...], OperationStreamDigest]:
    if type(operations) is not tuple:
        raise InvalidEngagementAuditInputError("operations must be an exact tuple.")
    if not operations:
        raise EmptyToolpathAuditError("engagement audit requires at least one operation.")
    if any(type(operation) is not ToolpathOperation for operation in operations):
        raise InvalidAuditOperationError("every operation must be an exact ToolpathOperation.")
    snapshots = tuple(snapshot_toolpath_operation(operation) for operation in operations)
    source_bytes = tuple(canonical_operation_snapshot_bytes(snapshot) for snapshot in snapshots)
    classified = tuple(classify_operation_snapshot(snapshot, cut_plane, operation_index=operation_index) for operation_index, snapshot in enumerate(snapshots))
    return classified, operation_stream_digest(source_bytes)


@dataclass(frozen=True, init=False)
class EngagementAuditInput:
    """Authenticated request containing no caller-owned mutable geometry."""

    design_boundary: CanonicalRingV1
    holes: tuple[CanonicalRingV1, ...]
    cut_plane: CutPlane
    tool_radius: ToolRadius
    engagement_cap: EngagementCap
    operations: tuple[AuthenticatedOperation, ...]
    operation_stream_digest: OperationStreamDigest
    build_identity: BuildIdentity

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise InvalidEngagementAuditInputError("EngagementAuditInput must be created by EngagementAuditInput.build().")

    @classmethod
    def build(
        cls,
        *,
        design_boundary: CanonicalRingV1,
        holes: tuple[CanonicalRingV1, ...],
        cut_plane: CutPlane,
        tool_radius: ToolRadius,
        engagement_cap: EngagementCap,
        operations: tuple[ToolpathOperation, ...],
        build_identity: BuildIdentity,
    ) -> Self:
        """Validate, classify natively, snapshot, and content-address one request."""
        _validate_authoritative_fields(design_boundary, cut_plane, tool_radius, engagement_cap, build_identity)
        validated_holes = _validated_holes(holes, require_canonical_order=False)
        classified_operations, stream_digest = _classify_operations(operations, cut_plane)
        instance = object.__new__(cls)
        object.__setattr__(instance, "design_boundary", design_boundary)
        object.__setattr__(instance, "holes", validated_holes)
        object.__setattr__(instance, "cut_plane", cut_plane)
        object.__setattr__(instance, "tool_radius", tool_radius)
        object.__setattr__(instance, "engagement_cap", engagement_cap)
        object.__setattr__(instance, "operations", classified_operations)
        object.__setattr__(instance, "operation_stream_digest", stream_digest)
        object.__setattr__(instance, "build_identity", build_identity)
        return instance

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not EngagementAuditInput:
            raise InvalidEngagementAuditInputError("audit input must be exact EngagementAuditInput, not a subclass.")
        return encode_tagged_union(
            AUDIT_INPUT_VERSION,
            encode_component_map(
                {
                    b"build-identity": self.build_identity.canonical_bytes,
                    b"arc-phase-strategy": _stock_2.audit_arc_phase_strategy_version(),
                    b"clearance-z": canonical_clearance_z_bytes(self.cut_plane.clearance_z),
                    b"cut-z": canonical_cut_z_bytes(self.cut_plane.cut_z),
                    b"design-boundary": self.design_boundary.canonical_bytes,
                    b"engagement-cap-chord-ratio": self.engagement_cap.chord_ratio_bytes,
                    b"holes": encode_sequence(tuple(hole.canonical_bytes for hole in self.holes)),
                    b"operation-stream-digest": bytes(self.operation_stream_digest),
                    b"tool-radius": encode_binary64(float(self.tool_radius.value)),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())
