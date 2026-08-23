from __future__ import annotations

from dataclasses import replace
import hashlib
import math
from typing import cast

import numpy as np
import pytest
from compas.geometry import Line

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.identity import ComponentDomainTag
from compas_cgal.adaptive.identity import ComponentIdentity
from compas_cgal.adaptive.identity import NativeSourceTreeDigest
from compas_cgal.adaptive.identity import SourceRevision
from compas_cgal.adaptive.identity import StrategyVersion
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement_audit.errors import ContradictoryOperationRoleError
from compas_cgal.engagement_audit.errors import EmptyToolpathAuditError
from compas_cgal.engagement_audit.errors import InvalidAuthenticatedLateralOperationError
from compas_cgal.engagement_audit.errors import InvalidAuditOperationError
from compas_cgal.engagement_audit.errors import InvalidEngagementAuditInputError
from compas_cgal.engagement_audit.errors import MultipleCutPlaneError
from compas_cgal.engagement_audit.errors import NonFiniteAuditGeometryError
from compas_cgal.engagement_audit.errors import UnsupportedAuditGeometryError
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.identity import PixiLockDigest
from compas_cgal.engagement_audit.identity import PythonSourceTreeDigest
from compas_cgal.engagement_audit.input import EngagementAuditInput
from compas_cgal.engagement_audit.input import OperationStreamDigest
from compas_cgal.engagement_audit.operation_identity import canonical_toolpath_operation_bytes
from compas_cgal.engagement_audit.operation_identity import OperationSnapshot
from compas_cgal.engagement_audit.operation_identity import operation_digest
from compas_cgal.engagement_audit.operation_identity import operation_stream_digest
from compas_cgal.engagement_audit.operation_identity import snapshot_toolpath_operation
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedOperation
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

CUT_Z = 0.0
CLEARANCE_Z = 5.0


def _digest(seed: bytes) -> bytes:
    return hashlib.sha256(seed).digest()


def _build_identity(*, revision: bytes = b"source-v1") -> BuildIdentity:
    component = ComponentIdentity.build(
        component_domain=ComponentDomainTag(b"engagement-auditor"),
        strategy_version=StrategyVersion(b"v1"),
        source_revision=SourceRevision(revision),
        native_source_tree_digest=NativeSourceTreeDigest(_digest(b"component-native")),
        canonical_parameter_bytes=encode_tagged_union(b"test-v1", encode_bytes(b"parameters")),
    )
    return BuildIdentity.build(
        components=(component,),
        native_source_tree_digest=NativeSourceTreeDigest(_digest(b"native")),
        python_source_tree_digest=PythonSourceTreeDigest(_digest(b"python")),
        pixi_lock_digest=PixiLockDigest(_digest(b"lock")),
    )


def _ring(*, extent: float = 10.0) -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(tuple(Point2[WorldXY].build(x, y) for x, y in ((0.0, 0.0), (extent, 0.0), (extent, extent), (0.0, extent))))


def _hole() -> CanonicalRingV1:
    return CanonicalRingV1.build_hole(tuple(Point2[WorldXY].build(x, y) for x, y in ((2.0, 2.0), (2.0, 3.0), (3.0, 3.0), (3.0, 2.0))))


def _line(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    operation: OperationType = OperationType.CUT,
    *,
    path_index: int = 0,
    clockwise: bool = False,
    start_tangent: np.ndarray | None = None,
    end_tangent: np.ndarray | None = None,
) -> ToolpathOperation:
    return ToolpathOperation(
        geometry=Line(start, end),
        operation=operation,
        path_index=path_index,
        clockwise=clockwise,
        start_tangent=start_tangent,
        end_tangent=end_tangent,
    )


def _audit_input(
    *,
    operations: tuple[ToolpathOperation, ...] | None = None,
    design_boundary: CanonicalRingV1 | None = None,
    holes: tuple[CanonicalRingV1, ...] = (),
    cut_plane: CutPlane | None = None,
    tool_radius: ToolRadius | None = None,
    engagement_cap: EngagementCap | None = None,
    build_identity: BuildIdentity | None = None,
) -> EngagementAuditInput:
    return EngagementAuditInput.build(
        design_boundary=design_boundary or _ring(),
        holes=holes,
        cut_plane=cut_plane or CutPlane.build(CutZ.build(CUT_Z), ClearanceZ.build(CLEARANCE_Z)),
        tool_radius=tool_radius or ToolRadius.build(2.0),
        engagement_cap=engagement_cap or EngagementCap.build(math.pi / 2.0),
        operations=operations if operations is not None else (_line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z)),),
        build_identity=build_identity or _build_identity(),
    )


def test_empty_toolpath_cannot_audit_clean() -> None:
    with pytest.raises(EmptyToolpathAuditError, match="at least one operation"):
        _audit_input(operations=())


def test_factory_rejects_untyped_cut_plane_with_named_input_error() -> None:
    with pytest.raises(InvalidEngagementAuditInputError, match="cut plane"):
        _audit_input(cut_plane=cast(CutPlane, object()))


def test_cut_plane_is_input_not_inferred_from_lower_motion() -> None:
    operations = (
        _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z)),
        _line((1.0, 2.0, -1.0), (4.0, 2.0, -1.0)),
    )

    with pytest.raises(MultipleCutPlaneError, match="exact declared plane"):
        _audit_input(operations=operations)


def test_ramped_lateral_motion_is_rejected_before_audit() -> None:
    operations = (_line((1.0, 1.0, CLEARANCE_Z), (4.0, 1.0, CUT_Z)),)

    with pytest.raises(UnsupportedAuditGeometryError, match="ramp"):
        _audit_input(operations=operations)


def test_nonfinite_geometry_is_rejected_before_canonical_encoding() -> None:
    operations = (_line((1.0, 1.0, CUT_Z), (math.nan, 1.0, CUT_Z)),)

    with pytest.raises(NonFiniteAuditGeometryError, match="finite"):
        _audit_input(operations=operations)


def test_mislabeled_motion_is_rejected_at_input_boundary() -> None:
    operations = (
        _line(
            (1.0, 1.0, CUT_Z),
            (4.0, 1.0, CUT_Z),
            OperationType.RETRACT,
        ),
    )

    with pytest.raises(ContradictoryOperationRoleError, match="exact lateral motion"):
        _audit_input(operations=operations)


def test_off_plane_motion_is_rejected_at_input_boundary() -> None:
    operations = (_line((1.0, 1.0, 2.0), (4.0, 1.0, 2.0), OperationType.LINK),)

    with pytest.raises(MultipleCutPlaneError, match="exact declared plane"):
        _audit_input(operations=operations)


def test_input_digest_binds_every_authoritative_field() -> None:
    first = _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z), path_index=0)
    second = _line((4.0, 1.0, CUT_Z), (4.0, 4.0, CUT_Z), path_index=1)
    baseline = _audit_input(operations=(first, second))
    mutations = (
        _audit_input(design_boundary=_ring(extent=11.0), operations=(first, second)),
        _audit_input(holes=(_hole(),), operations=(first, second)),
        _audit_input(
            cut_plane=CutPlane.build(CutZ.build(-1.0), ClearanceZ.build(CLEARANCE_Z)),
            operations=(
                _line((1.0, 1.0, -1.0), (4.0, 1.0, -1.0), path_index=0),
                _line((4.0, 1.0, -1.0), (4.0, 4.0, -1.0), path_index=1),
            ),
        ),
        _audit_input(
            cut_plane=CutPlane.build(CutZ.build(CUT_Z), ClearanceZ.build(6.0)),
            operations=(first, second),
        ),
        _audit_input(tool_radius=ToolRadius.build(2.5), operations=(first, second)),
        _audit_input(engagement_cap=EngagementCap.build(math.pi / 3.0), operations=(first, second)),
        _audit_input(
            operations=(
                _line((1.0, 1.0, CUT_Z), (5.0, 1.0, CUT_Z), path_index=0),
                second,
            )
        ),
        _audit_input(
            operations=(
                _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z), OperationType.LEAD_IN, path_index=0),
                second,
            )
        ),
        _audit_input(
            operations=(
                _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z), path_index=0, clockwise=True),
                second,
            )
        ),
        _audit_input(
            operations=(
                _line(
                    (1.0, 1.0, CUT_Z),
                    (4.0, 1.0, CUT_Z),
                    path_index=0,
                    start_tangent=np.array([1.0, 0.0, 0.0]),
                ),
                second,
            )
        ),
        _audit_input(
            operations=(
                _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z), path_index=2),
                second,
            )
        ),
        _audit_input(operations=(second, first)),
        _audit_input(build_identity=_build_identity(revision=b"source-v2"), operations=(first, second)),
    )

    assert len({baseline.digest, *(mutation.digest for mutation in mutations)}) == 1 + len(mutations)


def test_raw_construction_cannot_replace_authenticated_stream_digest() -> None:
    audit_input = _audit_input()

    with pytest.raises(InvalidEngagementAuditInputError, match="must be created"):
        replace(audit_input, operation_stream_digest=OperationStreamDigest(b"forged"))


def test_input_snapshots_legacy_operation_before_caller_mutation() -> None:
    operation = _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z))
    audit_input = _audit_input(operations=(operation,))
    original_digest = audit_input.digest
    original_operation = audit_input.operations[0]

    assert isinstance(original_operation, AuthenticatedLateralOperation)

    operation.operation = OperationType.RETRACT
    operation.geometry = Line((50.0, 50.0, CLEARANCE_Z), (60.0, 50.0, CLEARANCE_Z))
    operation.path_index = 99
    operation.clockwise = True

    assert audit_input.operations[0] == original_operation
    assert audit_input.digest == original_digest


def test_one_capture_owns_identity_and_native_classification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import compas_cgal.engagement_audit.input as input_module

    operation = _line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z))
    expected_source = canonical_toolpath_operation_bytes(operation)
    classify_snapshot = input_module.classify_operation_snapshot

    def mutate_before_native_call(
        snapshot: OperationSnapshot,
        cut_plane: CutPlane,
        *,
        operation_index: int,
    ) -> AuthenticatedOperation:
        operation.geometry = Line(
            (50.0, 50.0, CLEARANCE_Z),
            (60.0, 50.0, CLEARANCE_Z),
        )
        operation.operation = OperationType.LINK
        return classify_snapshot(
            snapshot,
            cut_plane,
            operation_index=operation_index,
        )

    monkeypatch.setattr(input_module, "classify_operation_snapshot", mutate_before_native_call)

    audit_input = _audit_input(operations=(operation,))
    classified = audit_input.operations[0]

    assert isinstance(classified, AuthenticatedLateralOperation)
    assert isinstance(classified.motion, _stock_2.AuditSegmentMotion2)
    assert classified.operation_digest == operation_digest(expected_source)
    assert audit_input.operation_stream_digest == operation_stream_digest((expected_source,))


def test_input_retains_opaque_plunge_after_caller_mutation() -> None:
    operation = _line(
        (2.0, 3.0, CLEARANCE_Z),
        (2.0, 3.0, CUT_Z),
        OperationType.PLUNGE,
    )
    audit_input = _audit_input(operations=(operation,))
    snapshot = audit_input.operations[0]

    assert isinstance(snapshot, AuthenticatedPlungeOperation)
    assert isinstance(snapshot.motion, _stock_2.AuditVerticalPlunge2)
    assert not hasattr(snapshot, "endpoint")

    operation.geometry = Line((8.0, 9.0, CLEARANCE_Z), (8.0, 9.0, CUT_Z))

    assert audit_input.operations[0] == snapshot


def test_input_identity_binds_exact_cap_surrogate() -> None:
    audit_input = _audit_input()

    assert audit_input.engagement_cap.chord_ratio_bytes in audit_input.canonical_bytes


def test_input_identity_binds_native_arc_phase_strategy() -> None:
    audit_input = _audit_input()

    assert _stock_2.audit_arc_phase_strategy_version() in audit_input.canonical_bytes


def test_authenticated_lateral_operation_has_its_own_failure_model() -> None:
    operation = _audit_input().operations[0]

    assert isinstance(operation, AuthenticatedLateralOperation)
    with pytest.raises(InvalidAuthenticatedLateralOperationError, match="operation index"):
        replace(operation, operation_index=-1)


def test_raw_operation_snapshot_cannot_bypass_typed_geometry() -> None:
    snapshot = snapshot_toolpath_operation(_line((1.0, 1.0, CUT_Z), (4.0, 1.0, CUT_Z)))

    with pytest.raises(InvalidAuditOperationError, match="world-XYZ"):
        replace(snapshot, start=object())
