from __future__ import annotations

import hashlib
import inspect
import math
from collections.abc import Callable

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
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
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement_audit.decision_limits import AuditDecisionLimits
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.identity import PixiLockDigest
from compas_cgal.engagement_audit.identity import PythonSourceTreeDigest
from compas_cgal.engagement_audit.input import EngagementAuditInput
from compas_cgal.engagement_audit.records import AuthenticatedNonEngagingOperation
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

CUT_Z = 0.0
CLEARANCE_Z = 5.0


def _digest(label: str) -> bytes:
    return hashlib.sha256(f"public-replay:{label}".encode()).digest()


def _build_identity() -> BuildIdentity:
    component = ComponentIdentity.build(
        component_domain=ComponentDomainTag(b"public-replay-test"),
        strategy_version=StrategyVersion(b"v1"),
        source_revision=SourceRevision(b"test-source"),
        native_source_tree_digest=NativeSourceTreeDigest(_digest("component-native")),
        canonical_parameter_bytes=encode_tagged_union(b"test-v1", encode_bytes(b"parameters")),
    )
    return BuildIdentity.build(
        components=(component,),
        native_source_tree_digest=NativeSourceTreeDigest(_digest("native")),
        python_source_tree_digest=PythonSourceTreeDigest(_digest("python")),
        pixi_lock_digest=PixiLockDigest(_digest("lock")),
    )


def _ring(
    x_min: float,
    y_min: float,
    x_max: float,
    y_max: float,
) -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(
        tuple(
            Point2[WorldXY].build(x, y)
            for x, y in (
                (x_min, y_min),
                (x_max, y_min),
                (x_max, y_max),
                (x_min, y_max),
            )
        )
    )


def _line(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    operation: OperationType,
    *,
    path_index: int,
) -> ToolpathOperation:
    return ToolpathOperation(
        geometry=Line(start, end),
        operation=operation,
        path_index=path_index,
    )


def _audit_input(
    operations: tuple[ToolpathOperation, ...],
    *,
    boundary: CanonicalRingV1 | None = None,
) -> EngagementAuditInput:
    return EngagementAuditInput.build(
        design_boundary=boundary or _ring(0.0, 0.0, 10.0, 10.0),
        holes=(),
        cut_plane=CutPlane.build(CutZ.build(CUT_Z), ClearanceZ.build(CLEARANCE_Z)),
        tool_radius=ToolRadius.build(0.5),
        engagement_cap=EngagementCap.build(math.pi / 2.0),
        depletion_policy=DepletionPolicy.build(
            chord_bound=ChordBound.build(0.02),
            center_count_limit=4096,
        ),
        decision_limits=AuditDecisionLimits.build(
            spatial_floor_mm=0.015625,
            max_depth=16,
            max_nodes=8192,
        ),
        operations=operations,
        build_identity=_build_identity(),
    )


def _closed_route_input() -> EngagementAuditInput:
    circle = Circle(radius=1.0, frame=Frame([30.0, 0.0, CUT_Z]))
    arc = Arc(
        radius=1.0,
        start_angle=0.0,
        end_angle=math.pi / 2.0,
        frame=Frame([40.0, 0.0, CUT_Z]),
    )
    return _audit_input(
        (
            _line((20.0, 0.0, CUT_Z), (21.0, 0.0, CUT_Z), OperationType.CUT, path_index=0),
            ToolpathOperation(
                geometry=circle,
                operation=OperationType.CUT,
                path_index=1,
                clockwise=False,
            ),
            ToolpathOperation(
                geometry=arc,
                operation=OperationType.CUT,
                path_index=2,
                clockwise=False,
            ),
            _line((3.0, 3.0, CLEARANCE_Z), (3.0, 3.0, CUT_Z), OperationType.PLUNGE, path_index=3),
            _line((3.0, 3.0, CUT_Z), (3.0, 3.0, CLEARANCE_Z), OperationType.RETRACT, path_index=4),
            _line((3.0, 3.0, CLEARANCE_Z), (4.0, 3.0, CLEARANCE_Z), OperationType.LINK, path_index=5),
        )
    )


def _ring_rows(ring: CanonicalRingV1) -> np.ndarray:
    return np.array(tuple((point.x, point.y) for point in ring.vertices), dtype=np.float64)


def _direct_native_run(
    audit_input: EngagementAuditInput,
) -> tuple[tuple[object, ...], _stock_2.AuditReplayCompletion2]:
    policy = _stock_2.build_audit_policy(
        audit_input.tool_radius.value,
        audit_input.engagement_cap.theta,
        audit_input.engagement_cap.chord_ratio,
        audit_input.depletion_policy.chord_bound.value,
        audit_input.depletion_policy.center_count_limit,
    )
    request = _stock_2.build_audit_native_request_identity(
        _ring_rows(audit_input.design_boundary),
        [_ring_rows(hole) for hole in audit_input.holes],
        policy,
        audit_input.decision_limits.native,
        tuple(operation.motion for operation in audit_input.operations),
    )
    replay = _stock_2.begin_audit_replay(
        _ring_rows(audit_input.design_boundary),
        [_ring_rows(hole) for hole in audit_input.holes],
        bytes(audit_input.digest),
        request,
        policy,
        audit_input.decision_limits.native,
        tuple(bytes(operation.digest) for operation in audit_input.operations),
    )
    results: list[object] = []
    for operation in audit_input.operations:
        if type(operation) is AuthenticatedLateralOperation:
            if type(operation.motion) is _stock_2.AuditSegmentMotion2:
                result = _stock_2.audit_deplete_segment(replay, operation.motion, bytes(operation.digest))
            elif type(operation.motion) is _stock_2.AuditCircleMotion2:
                result = _stock_2.audit_deplete_circle(replay, operation.motion, bytes(operation.digest))
            else:
                assert type(operation.motion) is _stock_2.AuditArcMotion2
                result = _stock_2.audit_deplete_arc(replay, operation.motion, bytes(operation.digest))
        elif type(operation) is AuthenticatedPlungeOperation:
            result = _stock_2.deplete_audit_plunge(replay, operation.motion, bytes(operation.digest))
        else:
            assert type(operation) is AuthenticatedNonEngagingOperation
            if type(operation.motion) is _stock_2.AuditVerticalRetract2:
                result = _stock_2.record_audit_retract(replay, operation.motion, bytes(operation.digest))
            else:
                assert type(operation.motion) is _stock_2.AuditClearanceTransport2
                result = _stock_2.record_audit_clearance(replay, operation.motion, bytes(operation.digest))
        results.append(result)
    return tuple(results), _stock_2.finish_audit_replay(replay)


def test_public_replay_accepts_only_authenticated_input() -> None:
    # Production mutation caught: public replay grows a second raw-geometry or
    # caller-owned stock argument.
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    assert tuple(inspect.signature(audit_toolpath_engagement).parameters) == ("audit_input",)


def test_package_root_exports_truthful_audit_surface() -> None:
    from compas_cgal.engagement_audit import EngagementAuditReport
    from compas_cgal.engagement_audit import MeasuredOperationAudit
    from compas_cgal.engagement_audit import NonEngagingOperationAudit
    from compas_cgal.engagement_audit import PlungeOperationAudit
    from compas_cgal.engagement_audit import audit_toolpath_engagement

    assert EngagementAuditReport.__name__ == "EngagementAuditReport"
    assert MeasuredOperationAudit.__name__ == "MeasuredOperationAudit"
    assert PlungeOperationAudit.__name__ == "PlungeOperationAudit"
    assert NonEngagingOperationAudit.__name__ == "NonEngagingOperationAudit"
    assert callable(audit_toolpath_engagement)


def test_public_replay_rejects_foreign_input_before_native_state_exists() -> None:
    # Production mutation caught: a structurally similar caller object enters
    # native replay without the factory-owned input identity.
    from compas_cgal.engagement_audit.errors import InvalidPublicAuditInputError
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    with pytest.raises(InvalidPublicAuditInputError):
        audit_toolpath_engagement(object())  # type: ignore[arg-type]


def test_public_replay_dispatches_each_authenticated_operation_once_in_order() -> None:
    # Production mutations caught: route dispatch skips, duplicates, or
    # reorders an authenticated stream member.
    from compas_cgal.engagement_audit.replay import _audit_toolpath_engagement

    audit_input = _closed_route_input()
    calls: list[tuple[str, int | None]] = []

    def observe(event: str, index: int | None, invoke: Callable[[], object]) -> object:
        calls.append((event, index))
        return invoke()

    report = _audit_toolpath_engagement(audit_input, observe)

    assert calls == [
        ("begin", None),
        ("segment", 0),
        ("circle", 1),
        ("arc", 2),
        ("plunge", 3),
        ("retract", 4),
        ("clearance", 5),
        ("finish", None),
    ]
    assert tuple(record.operation_index for record in report.operations) == tuple(range(6))
    assert tuple(record.authenticated_operation_digest for record in report.operations) == tuple(operation.digest for operation in audit_input.operations)


def test_middle_replay_exception_never_returns_a_partial_report() -> None:
    # Production mutation caught: replay catches a native failure and returns a
    # prefix report as if cardinality and finalization had succeeded.
    from compas_cgal.engagement_audit.replay import _audit_toolpath_engagement

    class InjectedNativeReplayFailure(RuntimeError):
        pass

    calls: list[tuple[str, int | None]] = []

    def fail_arc(event: str, index: int | None, invoke: Callable[[], object]) -> object:
        calls.append((event, index))
        if event == "arc":
            raise InjectedNativeReplayFailure("injected after committed prefix")
        return invoke()

    with pytest.raises(InjectedNativeReplayFailure, match="committed prefix"):
        _audit_toolpath_engagement(_closed_route_input(), fail_arc)
    assert calls == [
        ("begin", None),
        ("segment", 0),
        ("circle", 1),
        ("arc", 2),
    ]


def test_finish_exception_never_mints_a_report() -> None:
    from compas_cgal.engagement_audit.replay import _audit_toolpath_engagement

    class InjectedNativeFinalizationFailure(RuntimeError):
        pass

    calls: list[tuple[str, int | None]] = []

    def fail_finish(event: str, index: int | None, invoke: Callable[[], object]) -> object:
        calls.append((event, index))
        if event == "finish":
            raise InjectedNativeFinalizationFailure("injected finalization")
        return invoke()

    with pytest.raises(InjectedNativeFinalizationFailure, match="finalization"):
        _audit_toolpath_engagement(_closed_route_input(), fail_finish)
    assert calls[-1] == ("finish", None)


@pytest.mark.parametrize(
    ("event_to_fail", "native_error_name", "public_error_name"),
    (
        (
            "begin",
            "AuditReplayRequestIdentityError",
            "InvalidNativeAuditReplayRequestError",
        ),
        (
            "segment",
            "AuditReplayOperationIdentityError",
            "InvalidNativeAuditReplayOperationError",
        ),
        (
            "segment",
            "AuditReportingObservationError",
            "InvalidNativeReportingObservationError",
        ),
        (
            "finish",
            "AuditReplayIncompleteError",
            "InvalidNativeAuditReplayCompletionError",
        ),
    ),
)
def test_public_boundary_translates_native_failures_with_chained_cause(
    event_to_fail: str,
    native_error_name: str,
    public_error_name: str,
) -> None:
    import compas_cgal.engagement_audit.errors as public_errors
    from compas_cgal.engagement_audit.replay import _audit_toolpath_engagement

    native_error = getattr(_stock_2, native_error_name)
    public_error = getattr(public_errors, public_error_name)

    def fail_native(
        event: str,
        _index: int | None,
        invoke: Callable[[], object],
    ) -> object:
        if event == event_to_fail:
            raise native_error(f"injected {event_to_fail} failure")
        return invoke()

    with pytest.raises(public_error, match=event_to_fail) as caught:
        _audit_toolpath_engagement(_closed_route_input(), fail_native)
    assert isinstance(caught.value.__cause__, native_error)


def test_public_replay_preserves_native_seed_and_complete_lineage_chain() -> None:
    # Production mutations caught: Python recomputes a lineage, drops an
    # unresolved transition, or reports a terminal lineage not finalized by
    # native replay.
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    audit_input = _closed_route_input()
    _native_results, completion = _direct_native_run(audit_input)
    report = audit_toolpath_engagement(audit_input)

    assert report.operations[0].pre_motion_stock_lineage == completion.seed_lineage
    for previous, following in zip(report.operations, report.operations[1:], strict=False):
        assert previous.post_motion_stock_lineage == following.pre_motion_stock_lineage
    assert report.terminal_stock_lineage == completion.terminal_lineage
    assert report.native_completion_digest == completion.digest


def test_public_replay_has_no_legacy_or_reclassification_route(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    audit_input = _closed_route_input()

    def forbidden(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("legacy or reclassification route reached")

    monkeypatch.setattr(_stock_2, "Stock2", forbidden)
    monkeypatch.setattr(_stock_2, "engagement_at", forbidden)
    monkeypatch.setattr(_stock_2, "certify_segment_tea", forbidden)

    import compas_cgal.engagement_audit.classification as classification

    monkeypatch.setattr(classification, "classify_operation", forbidden)
    monkeypatch.setattr(classification, "classify_operation_snapshot", forbidden)

    import compas_cgal.engagement_audit.native as native
    import compas_cgal.engagement_audit.replay as replay

    source = inspect.getsource(replay) + inspect.getsource(native)
    for forbidden_symbol in (
        "Stock2",
        "engagement_at",
        "certify_segment_tea",
        "classify_operation",
        "ToolpathOperation",
        ".geometry",
    ):
        assert forbidden_symbol not in source

    report = replay.audit_toolpath_engagement(audit_input)
    assert report.operation_count == len(audit_input.operations)


def test_public_records_are_closed_factory_owned_projections_of_native_results() -> None:
    # Production mutations caught: output records expose mutable geometry,
    # admit raw construction, merge plunge with non-engaging transport, or
    # conflate exact evidence with reporting work.
    from compas_cgal.engagement_audit.records import MeasuredOperationAudit
    from compas_cgal.engagement_audit.records import NonEngagingOperationAudit
    from compas_cgal.engagement_audit.records import PlungeOperationAudit
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    measured, circle, arc, plunge, retract, clearance = audit_toolpath_engagement(_closed_route_input()).operations

    assert type(measured) is MeasuredOperationAudit
    assert type(circle) is MeasuredOperationAudit
    assert type(arc) is MeasuredOperationAudit
    assert type(plunge) is PlungeOperationAudit
    assert type(retract) is NonEngagingOperationAudit
    assert type(clearance) is NonEngagingOperationAudit
    assert retract.reason == "vertical_retract"
    assert clearance.reason == "clearance_transport"
    assert not hasattr(plunge, "reason")
    assert not hasattr(retract, "verdict")
    assert not hasattr(measured, "motion")
    assert measured.evidence_count > 0
    assert measured.reporting_station_count > 0
    assert isinstance(measured.native_decision_digest, bytes)
    assert isinstance(measured.depletion_witness_digest, bytes)
    assert isinstance(measured.reporting_observation_digest, bytes)
    assert isinstance(measured.native_result_digest, bytes)
    assert isinstance(plunge.native_result_digest, bytes)
    assert isinstance(retract.native_result_digest, bytes)
    from compas_cgal.engagement_audit.errors import InvalidMeasuredOperationAuditError
    from compas_cgal.engagement_audit.errors import InvalidNonEngagingOperationAuditError
    from compas_cgal.engagement_audit.errors import InvalidPlungeOperationAuditError

    with pytest.raises(InvalidMeasuredOperationAuditError):
        MeasuredOperationAudit()  # type: ignore[call-arg]
    with pytest.raises(InvalidPlungeOperationAuditError):
        PlungeOperationAudit()  # type: ignore[call-arg]
    with pytest.raises(InvalidNonEngagingOperationAuditError):
        NonEngagingOperationAudit()  # type: ignore[call-arg]
