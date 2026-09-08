from __future__ import annotations

import hashlib
import inspect
import math
from dataclasses import replace

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Frame
from compas.geometry import Line

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
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
from compas_cgal import _stock_2
from compas_cgal.engagement_audit.decision_limits import AuditDecisionLimits
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.identity import PixiLockDigest
from compas_cgal.engagement_audit.identity import PythonSourceTreeDigest
from compas_cgal.engagement_audit.input import EngagementAuditInput
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedNonEngagingOperation
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

CUT_Z = 0.0
CLEARANCE_Z = 5.0


def _digest(label: str) -> bytes:
    return hashlib.sha256(f"truthful-report:{label}".encode()).digest()


def _identity(source_revision: bytes = b"test-source") -> BuildIdentity:
    component = ComponentIdentity.build(
        component_domain=ComponentDomainTag(b"truthful-report-test"),
        strategy_version=StrategyVersion(b"v1"),
        source_revision=SourceRevision(source_revision),
        native_source_tree_digest=NativeSourceTreeDigest(_digest("component-native")),
        canonical_parameter_bytes=encode_tagged_union(b"test-v1", encode_bytes(b"parameters")),
    )
    return BuildIdentity.build(
        components=(component,),
        native_source_tree_digest=NativeSourceTreeDigest(_digest("native")),
        python_source_tree_digest=PythonSourceTreeDigest(_digest("python")),
        pixi_lock_digest=PixiLockDigest(_digest("lock")),
    )


def _boundary(
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


def _input(
    boundary: CanonicalRingV1,
    operations: tuple[ToolpathOperation, ...],
    *,
    build_identity: BuildIdentity | None = None,
) -> EngagementAuditInput:
    return EngagementAuditInput.build(
        design_boundary=boundary,
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
        build_identity=build_identity or _identity(),
    )


def _line_operation(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    role: OperationType = OperationType.CUT,
) -> ToolpathOperation:
    return ToolpathOperation(
        geometry=Line(start, end),
        operation=role,
        path_index=0,
    )


def _ring_rows(ring: CanonicalRingV1) -> np.ndarray:
    return np.array(
        tuple((point.x, point.y) for point in ring.vertices),
        dtype=np.float64,
    )


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


def _certified_input() -> EngagementAuditInput:
    return _input(
        _boundary(-10.0, -10.0, 10.0, -0.4),
        (_line_operation((-1.0, 0.0, CUT_Z), (1.0, 0.0, CUT_Z)),),
    )


def _cap_exceeded_input() -> EngagementAuditInput:
    return _input(
        _boundary(-1.6, -0.6, -0.4, 0.6),
        (_line_operation((-2.0, 0.0, CUT_Z), (2.0, 0.0, CUT_Z)),),
    )


def _unresolved_input() -> EngagementAuditInput:
    arc = Arc(
        radius=2.0,
        start_angle=0.0,
        end_angle=math.pi / 2.0,
        frame=Frame([0.0, 0.0, CUT_Z]),
    )
    return _input(
        _boundary(-2.6, -0.6, -1.4, 0.6),
        (
            ToolpathOperation(
                geometry=arc,
                operation=OperationType.CUT,
                path_index=0,
                clockwise=False,
            ),
        ),
    )


def _retract_only_input() -> EngagementAuditInput:
    return _input(
        _boundary(0.0, 0.0, 10.0, 10.0),
        (
            _line_operation(
                (2.0, 2.0, CUT_Z),
                (2.0, 2.0, CLEARANCE_Z),
                OperationType.RETRACT,
            ),
        ),
    )


def _plunge_only_input() -> EngagementAuditInput:
    return _input(
        _boundary(0.0, 0.0, 10.0, 10.0),
        (
            _line_operation(
                (2.0, 2.0, CLEARANCE_Z),
                (2.0, 2.0, CUT_Z),
                OperationType.PLUNGE,
            ),
        ),
    )


def _two_certified_input(
    *,
    build_identity: BuildIdentity | None = None,
) -> EngagementAuditInput:
    return _input(
        _boundary(-10.0, -10.0, 10.0, -0.4),
        (
            _line_operation((-1.0, 0.0, CUT_Z), (0.0, 0.0, CUT_Z)),
            ToolpathOperation(
                geometry=Line((0.0, 0.0, CUT_Z), (1.0, 0.0, CUT_Z)),
                operation=OperationType.CUT,
                path_index=1,
            ),
        ),
        build_identity=build_identity,
    )


def _certified_with_nonmeasured_input() -> EngagementAuditInput:
    return _input(
        _boundary(-10.0, -10.0, 10.0, -0.4),
        (
            _line_operation((-1.0, 0.0, CUT_Z), (1.0, 0.0, CUT_Z)),
            ToolpathOperation(
                geometry=Line((20.0, 0.0, CLEARANCE_Z), (20.0, 0.0, CUT_Z)),
                operation=OperationType.PLUNGE,
                path_index=1,
            ),
            ToolpathOperation(
                geometry=Line((20.0, 0.0, CUT_Z), (20.0, 0.0, CLEARANCE_Z)),
                operation=OperationType.RETRACT,
                path_index=2,
            ),
        ),
    )


def _mixed_exceeded_unresolved_input() -> EngagementAuditInput:
    boundary = CanonicalRingV1.build_outer(
        tuple(
            Point2[WorldXY].build(x, y)
            for x, y in (
                (-11.6, -0.6),
                (-10.4, -0.6),
                (-10.4, 0.5),
                (7.4, 0.5),
                (7.4, -0.6),
                (8.6, -0.6),
                (8.6, 0.6),
                (7.4, 0.6),
                (7.4, 0.55),
                (-10.4, 0.55),
                (-10.4, 0.6),
                (-11.6, 0.6),
            )
        )
    )
    arc = Arc(
        radius=2.0,
        start_angle=0.0,
        end_angle=math.pi / 2.0,
        frame=Frame([10.0, 0.0, CUT_Z]),
    )
    return _input(
        boundary,
        (
            _line_operation((-12.0, 0.0, CUT_Z), (-8.0, 0.0, CUT_Z)),
            ToolpathOperation(
                geometry=arc,
                operation=OperationType.CUT,
                path_index=1,
                clockwise=False,
            ),
        ),
    )


def _all_five_arms_input() -> EngagementAuditInput:
    boundary = CanonicalRingV1.build_outer(
        tuple(
            Point2[WorldXY].build(x, y)
            for x, y in (
                (-11.6, -0.6),
                (-10.4, -0.6),
                (-10.4, 0.5),
                (7.4, 0.5),
                (7.4, -0.6),
                (8.6, -0.6),
                (8.6, 0.6),
                (7.4, 0.6),
                (7.4, 0.55),
                (-10.4, 0.55),
                (-10.4, 0.6),
                (-11.6, 0.6),
            )
        )
    )
    arc = Arc(
        radius=2.0,
        start_angle=0.0,
        end_angle=math.pi / 2.0,
        frame=Frame([10.0, 0.0, CUT_Z]),
    )
    return _input(
        boundary,
        (
            _line_operation((-12.0, 0.0, CUT_Z), (-8.0, 0.0, CUT_Z)),
            ToolpathOperation(
                geometry=arc,
                operation=OperationType.CUT,
                path_index=1,
                clockwise=False,
            ),
            ToolpathOperation(
                geometry=Line((20.0, 0.0, CUT_Z), (21.0, 0.0, CUT_Z)),
                operation=OperationType.CUT,
                path_index=2,
            ),
            ToolpathOperation(
                geometry=Line((0.0, 0.0, CLEARANCE_Z), (0.0, 0.0, CUT_Z)),
                operation=OperationType.PLUNGE,
                path_index=3,
            ),
            ToolpathOperation(
                geometry=Line((0.0, 0.0, CUT_Z), (0.0, 0.0, CLEARANCE_Z)),
                operation=OperationType.RETRACT,
                path_index=4,
            ),
        ),
    )


def test_report_derives_closed_counts_and_content_identity() -> None:
    # Production mutations caught: report accepts caller aggregates, conflates
    # verdicts, or omits its ordered operations/input identity from the digest.
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    audit_input = _certified_input()
    report = audit_toolpath_engagement(audit_input)
    _results, completion = _direct_native_run(audit_input)

    assert report.operation_count == 1
    assert report.certified_count == 1
    assert report.cap_exceeded_count == 0
    assert report.unresolved_count == 0
    assert report.plunge_count == 0
    assert report.non_engaging_count == 0
    assert report.certified_count + report.cap_exceeded_count + report.unresolved_count + report.plunge_count + report.non_engaging_count == report.operation_count
    assert report.max_tea == report.operations[0].max_tea
    assert report.seed_stock_lineage == completion.seed_lineage
    assert report.terminal_stock_lineage == report.operations[-1].post_motion_stock_lineage
    assert report.native_completion_digest == completion.digest
    expected_canonical = encode_tagged_union(
        b"engagement-audit-report-v1",
        encode_component_map(
            {
                b"audit-input-digest": bytes(audit_input.digest),
                b"cap-exceeded-count": encode_integer(report.cap_exceeded_count),
                b"certified-count": encode_integer(report.certified_count),
                b"max-tea-radian": encode_binary64(float(report.max_tea)),
                b"native-completion-digest": report.native_completion_digest,
                b"native-request-digest": bytes(audit_input.native_request_digest),
                b"non-engaging-count": encode_integer(report.non_engaging_count),
                b"operation-count": encode_integer(report.operation_count),
                b"operation-digests": encode_sequence(
                    tuple(bytes(operation.digest) for operation in report.operations)
                ),
                b"plunge-count": encode_integer(report.plunge_count),
                b"seed-stock-lineage": report.seed_stock_lineage,
                b"terminal-stock-lineage": report.terminal_stock_lineage,
                b"unresolved-count": encode_integer(report.unresolved_count),
            }
        ),
    )
    assert report.canonical_bytes == expected_canonical
    assert report.digest == hashlib.sha256(report.canonical_bytes).digest()


def test_unresolved_is_not_a_cap_violation_or_certification() -> None:
    # Production mutation caught: UNRESOLVED collapses into either Boolean
    # success or a proved cap violation.
    from compas_cgal.engagement_audit.errors import UnresolvedEngagementAuditError
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_unresolved_input())

    assert report.certified_count == 0
    assert report.cap_exceeded_count == 0
    assert report.unresolved_count == 1
    with pytest.raises(UnresolvedEngagementAuditError, match="1 unresolved"):
        report.require_certified()


def test_proved_cap_exceedance_has_its_own_certification_failure() -> None:
    # Production mutation caught: a proved violation is reported as generic
    # unresolved work or silently accepted.
    from compas_cgal.engagement_audit.errors import CapExceededToolpathError
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_cap_exceeded_input())

    assert report.certified_count == 0
    assert report.cap_exceeded_count == 1
    assert report.unresolved_count == 0
    with pytest.raises(CapExceededToolpathError, match="1 cap-exceeded"):
        report.require_certified()


def test_report_requires_a_measured_lateral_motion() -> None:
    # Production mutation caught: a retract-only stream is advertised as a
    # certified engagement measurement.
    from compas_cgal.engagement_audit.errors import NoMeasuredLateralMotionError
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_retract_only_input())

    assert report.operation_count == 1
    assert report.non_engaging_count == 1
    assert report.certified_count == 0
    with pytest.raises(NoMeasuredLateralMotionError):
        report.require_certified()

    plunge_report = audit_toolpath_engagement(_plunge_only_input())
    assert plunge_report.plunge_count == 1
    with pytest.raises(NoMeasuredLateralMotionError):
        plunge_report.require_certified()


def test_certified_lateral_with_plunge_and_transport_is_certified() -> None:
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_certified_with_nonmeasured_input())

    assert report.certified_count == 1
    assert report.plunge_count == 1
    assert report.non_engaging_count == 1
    assert report.require_certified() is None


def test_report_factory_accepts_only_native_completion_not_caller_aggregates() -> None:
    # Production mutation caught: callers supply counts, maximum, ordered
    # digest, or terminal lineage instead of the report deriving them.
    from compas_cgal.engagement_audit.report import EngagementAuditReport

    assert tuple(inspect.signature(EngagementAuditReport.build).parameters) == (
        "audit_input",
        "operations",
        "completion",
    )


def test_reporting_observation_is_authenticated_but_cannot_change_native_truth() -> None:
    # Production mutations caught: reporting work overwrites the native
    # verdict/digest or is accepted without its own authenticated identity.
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_certified_input())
    record = report.operations[0]

    assert record.verdict == "certified"
    assert record.max_tea >= 0.0
    assert record.reporting_station_count > 0
    assert isinstance(record.reporting_observation_digest, bytes)
    assert len(record.reporting_observation_digest) == hashlib.sha256().digest_size
    assert isinstance(record.native_result_digest, bytes)
    assert record.native_decision_digest != record.reporting_observation_digest
    assert record.evidence_count > 0


def test_all_five_report_arms_are_disjoint_and_derived() -> None:
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_all_five_arms_input())

    assert (
        report.certified_count,
        report.cap_exceeded_count,
        report.unresolved_count,
        report.plunge_count,
        report.non_engaging_count,
    ) == (1, 1, 1, 1, 1)
    assert report.certified_count + report.cap_exceeded_count + report.unresolved_count + report.plunge_count + report.non_engaging_count == report.operation_count == 5


def test_cap_exceeded_precedes_unresolved_but_reports_both() -> None:
    from compas_cgal.engagement_audit.errors import CapExceededToolpathError
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    report = audit_toolpath_engagement(_mixed_exceeded_unresolved_input())

    assert report.cap_exceeded_count == 1
    assert report.unresolved_count == 1
    with pytest.raises(
        CapExceededToolpathError,
        match=r"1 cap-exceeded.*1 unresolved",
    ):
        report.require_certified()


def test_report_factory_rejects_foreign_completion_and_cross_replay_lineage() -> None:
    from compas_cgal.engagement_audit.errors import ReportCompletionError
    from compas_cgal.engagement_audit.errors import ReportLineageContinuityError
    from compas_cgal.engagement_audit.report import EngagementAuditReport
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    audit_input = _two_certified_input()
    report = audit_toolpath_engagement(audit_input)
    _native_results, completion = _direct_native_run(audit_input)
    rebuilt = EngagementAuditReport.build(audit_input, report.operations, completion)
    assert rebuilt.digest == report.digest

    foreign_input = _two_certified_input(build_identity=_identity(b"foreign-build-source"))
    foreign_report = audit_toolpath_engagement(foreign_input)
    _foreign_results, foreign_completion = _direct_native_run(foreign_input)

    with pytest.raises(ReportCompletionError, match="input|seed"):
        EngagementAuditReport.build(
            audit_input,
            report.operations,
            foreign_completion,
        )
    with pytest.raises(ReportLineageContinuityError, match="seed|first"):
        EngagementAuditReport.build(
            audit_input,
            foreign_report.operations,
            completion,
        )


def test_report_rejects_index_reorder_duplicate_and_hole() -> None:
    from compas_cgal.engagement_audit.errors import ReportOperationCoverageError
    from compas_cgal.engagement_audit.report import EngagementAuditReport
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    audit_input = _two_certified_input()
    report = audit_toolpath_engagement(audit_input)
    _results, completion = _direct_native_run(audit_input)
    first, second = report.operations
    ordered = encode_sequence((bytes(first.digest), bytes(second.digest)))
    reversed_order = encode_sequence((bytes(second.digest), bytes(first.digest)))
    assert ordered in report.canonical_bytes
    assert reversed_order not in report.canonical_bytes

    for invalid in ((second, first), (first, first), (second,)):
        with pytest.raises(ReportOperationCoverageError):
            EngagementAuditReport.build(audit_input, invalid, completion)


def test_report_is_factory_owned_sealed_and_revalidates_derived_state() -> None:
    from compas_cgal.engagement_audit.errors import InvalidEngagementAuditReportError
    from compas_cgal.engagement_audit.report import EngagementAuditReport
    from compas_cgal.engagement_audit.replay import audit_toolpath_engagement

    with pytest.raises(InvalidEngagementAuditReportError):
        EngagementAuditReport()  # type: ignore[call-arg]

    forged = object.__new__(EngagementAuditReport)
    with pytest.raises(InvalidEngagementAuditReportError):
        _ = forged.canonical_bytes

    with pytest.raises(TypeError):

        class ForgedReport(EngagementAuditReport):
            pass

    report = audit_toolpath_engagement(_certified_input())
    with pytest.raises(InvalidEngagementAuditReportError):
        replace(report, certified_count=99)

    object.__setattr__(report, "certified_count", 99)
    with pytest.raises(InvalidEngagementAuditReportError):
        _ = report.canonical_bytes
    with pytest.raises(InvalidEngagementAuditReportError):
        _ = report.digest
