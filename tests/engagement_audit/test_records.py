from __future__ import annotations

import hashlib
from copy import copy
from dataclasses import fields
from dataclasses import replace

import numpy as np
import pytest

from compas_cgal import _stock_2
from compas_cgal.engagement_audit.errors import InvalidMeasuredOperationAuditError
from compas_cgal.engagement_audit.errors import InvalidNonEngagingOperationAuditError
from compas_cgal.engagement_audit.errors import InvalidPlungeOperationAuditError
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedNonEngagingOperation
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.engagement_audit.records import MeasuredOperationAudit
from compas_cgal.engagement_audit.records import NonEngagingOperationAudit
from compas_cgal.engagement_audit.records import OperationDigest
from compas_cgal.engagement_audit.records import PlungeOperationAudit


def _digest(seed: bytes) -> bytes:
    return hashlib.sha256(seed).digest()


def _line_classification(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    role: str,
) -> _stock_2.AuditLineClassification2:
    return _stock_2.classify_audit_line(start, end, 0.0, 5.0, role)


def _native_record_inputs(
    *,
    source_suffix: bytes = b"baseline",
    stock_floor_y: float = -0.4,
) -> tuple[
    tuple[
        AuthenticatedLateralOperation,
        AuthenticatedPlungeOperation,
        AuthenticatedNonEngagingOperation,
    ],
    tuple[
        _stock_2.AuditLateralResult2,
        _stock_2.AuditPlungeResult2,
        _stock_2.AuditNonEngagingResult2,
    ],
]:
    segment = _line_classification((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), "cut")
    plunge = _line_classification((20.0, 0.0, 5.0), (20.0, 0.0, 0.0), "plunge")
    retract = _line_classification((20.0, 0.0, 0.0), (20.0, 0.0, 5.0), "retract")
    assert isinstance(segment, _stock_2.AuditSegmentMotion2)
    assert isinstance(plunge, _stock_2.AuditVerticalPlunge2)
    assert isinstance(retract, _stock_2.AuditVerticalRetract2)

    carriers = (
        AuthenticatedLateralOperation.build(
            operation_index=0,
            operation_digest=OperationDigest(_digest(b"segment:" + source_suffix)),
            motion=segment,
        ),
        AuthenticatedPlungeOperation.build(
            operation_index=1,
            operation_digest=OperationDigest(_digest(b"plunge:" + source_suffix)),
            motion=plunge,
        ),
        AuthenticatedNonEngagingOperation.build(
            operation_index=2,
            operation_digest=OperationDigest(_digest(b"retract:" + source_suffix)),
            motion=retract,
        ),
    )
    boundary = np.array(
        (
            (-5.0, -5.0),
            (5.0, -5.0),
            (5.0, stock_floor_y),
            (-5.0, stock_floor_y),
        ),
        dtype=np.float64,
    )
    cap = np.pi / 2.0
    policy = _stock_2.build_audit_policy(0.5, cap, 2.0, 0.02, 4096)
    limits = _stock_2.build_audit_decision_limits(0.015625, 8, 256)
    request = _stock_2.build_audit_native_request_identity(
        boundary,
        [],
        policy,
        limits,
        tuple(carrier.motion for carrier in carriers),
    )
    replay = _stock_2.begin_audit_replay(
        boundary,
        [],
        _digest(b"record-input"),
        request,
        policy,
        limits,
        tuple(bytes(carrier.digest) for carrier in carriers),
    )
    results = (
        _stock_2.audit_deplete_segment(
            replay,
            segment,
            bytes(carriers[0].digest),
        ),
        _stock_2.deplete_audit_plunge(
            replay,
            plunge,
            bytes(carriers[1].digest),
        ),
        _stock_2.record_audit_retract(
            replay,
            retract,
            bytes(carriers[2].digest),
        ),
    )
    _stock_2.finish_audit_replay(replay)
    return carriers, results


def test_authenticated_carrier_digest_binds_domain_index_and_source() -> None:
    operation_digest = OperationDigest(_digest(b"source"))
    segment = _line_classification((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "cut")
    plunge = _line_classification((0.0, 0.0, 5.0), (0.0, 0.0, 0.0), "plunge")
    retract = _line_classification((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), "retract")
    assert isinstance(segment, _stock_2.AuditSegmentMotion2)
    assert isinstance(plunge, _stock_2.AuditVerticalPlunge2)
    assert isinstance(retract, _stock_2.AuditVerticalRetract2)
    lateral = AuthenticatedLateralOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=segment,
    )
    authenticated_plunge = AuthenticatedPlungeOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=plunge,
    )
    non_engaging = AuthenticatedNonEngagingOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=retract,
    )

    variants = (
        lateral,
        replace(lateral, operation_index=1),
        replace(lateral, operation_digest=OperationDigest(_digest(b"other-source"))),
        authenticated_plunge,
        non_engaging,
    )

    assert len({variant.digest for variant in variants}) == len(variants)
    assert all(bytes(variant.digest) == hashlib.sha256(variant.canonical_bytes).digest() for variant in variants)


def test_authenticated_carrier_digest_binds_closed_native_tag() -> None:
    operation_digest = OperationDigest(_digest(b"same-source"))
    segment = _line_classification((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "cut")
    circle = _stock_2.classify_audit_circle(
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        1.0,
        False,
        0.0,
        5.0,
        "cut",
    )
    arc = _stock_2.classify_audit_arc(
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        1.0,
        0.0,
        1.0,
        False,
        0.0,
        5.0,
        "cut",
    )
    retract = _line_classification((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), "retract")
    clearance = _line_classification((0.0, 0.0, 5.0), (1.0, 0.0, 5.0), "link")
    assert isinstance(segment, _stock_2.AuditSegmentMotion2)
    assert isinstance(circle, _stock_2.AuditCircleMotion2)
    assert isinstance(arc, _stock_2.AuditArcMotion2)
    assert isinstance(retract, _stock_2.AuditVerticalRetract2)
    assert isinstance(clearance, _stock_2.AuditClearanceTransport2)

    variants = (
        AuthenticatedLateralOperation.build(operation_index=0, operation_digest=operation_digest, motion=segment),
        AuthenticatedLateralOperation.build(operation_index=0, operation_digest=operation_digest, motion=circle),
        AuthenticatedLateralOperation.build(operation_index=0, operation_digest=operation_digest, motion=arc),
        AuthenticatedNonEngagingOperation.build(operation_index=0, operation_digest=operation_digest, motion=retract),
        AuthenticatedNonEngagingOperation.build(operation_index=0, operation_digest=operation_digest, motion=clearance),
    )

    assert len({variant.digest for variant in variants}) == len(variants)


def test_authenticated_carrier_v2_binds_native_motion_digest() -> None:
    operation_digest = OperationDigest(_digest(b"same-source"))
    first = _line_classification((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "cut")
    second = _line_classification((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), "cut")
    assert isinstance(first, _stock_2.AuditSegmentMotion2)
    assert isinstance(second, _stock_2.AuditSegmentMotion2)

    first_carrier = AuthenticatedLateralOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=first,
    )
    second_carrier = AuthenticatedLateralOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=second,
    )

    assert b"authenticated-lateral-operation-v2" in first_carrier.canonical_bytes
    assert first.digest in first_carrier.canonical_bytes
    assert first_carrier.digest != second_carrier.digest


def test_all_carrier_v2_variants_bind_native_motion_digest() -> None:
    operation_digest = OperationDigest(_digest(b"same-source"))
    first_plunge = _line_classification((0.0, 0.0, 5.0), (0.0, 0.0, 0.0), "plunge")
    second_plunge = _line_classification((1.0, 0.0, 5.0), (1.0, 0.0, 0.0), "plunge")
    first_clearance = _line_classification((0.0, 0.0, 5.0), (1.0, 0.0, 5.0), "link")
    second_clearance = _line_classification((0.0, 0.0, 5.0), (2.0, 0.0, 5.0), "link")
    first_retract = _line_classification((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), "retract")
    second_retract = _line_classification((1.0, 0.0, 0.0), (1.0, 0.0, 5.0), "retract")
    assert isinstance(first_plunge, _stock_2.AuditVerticalPlunge2)
    assert isinstance(second_plunge, _stock_2.AuditVerticalPlunge2)
    assert isinstance(first_clearance, _stock_2.AuditClearanceTransport2)
    assert isinstance(second_clearance, _stock_2.AuditClearanceTransport2)
    assert isinstance(first_retract, _stock_2.AuditVerticalRetract2)
    assert isinstance(second_retract, _stock_2.AuditVerticalRetract2)

    plunge = AuthenticatedPlungeOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=first_plunge,
    )
    changed_plunge = replace(plunge, motion=second_plunge)
    clearance = AuthenticatedNonEngagingOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=first_clearance,
    )
    changed_clearance = replace(clearance, motion=second_clearance)
    retract = AuthenticatedNonEngagingOperation.build(
        operation_index=0,
        operation_digest=operation_digest,
        motion=first_retract,
    )
    changed_retract = replace(retract, motion=second_retract)

    assert b"authenticated-plunge-operation-v2" in plunge.canonical_bytes
    assert first_plunge.digest in plunge.canonical_bytes
    assert plunge.digest != changed_plunge.digest
    assert b"authenticated-non-engaging-operation-v2" in clearance.canonical_bytes
    assert first_clearance.digest in clearance.canonical_bytes
    assert clearance.digest != changed_clearance.digest
    assert first_retract.digest in retract.canonical_bytes
    assert retract.digest != changed_retract.digest


def test_output_records_derive_only_from_exact_native_result_pairs() -> None:
    carriers, results = _native_record_inputs()

    measured = MeasuredOperationAudit.from_native(
        carriers[0],
        results[0],
        expected_native_request_digest=results[0].request_digest,
    )
    plunge = PlungeOperationAudit.from_native(
        carriers[1],
        results[1],
        expected_native_request_digest=results[1].request_digest,
    )
    non_engaging = NonEngagingOperationAudit.from_native(
        carriers[2],
        results[2],
        expected_native_request_digest=results[2].request_digest,
    )

    assert measured.operation_index == 0
    assert measured.authenticated_operation_digest == carriers[0].digest
    assert measured.verdict == results[0].verdict
    assert measured.evidence_count == results[0].evidence_count
    assert measured.reporting_station_count == results[0].reporting_observation.station_count
    assert measured.max_tea == results[0].reporting_observation.max_tea
    assert measured.native_decision_digest == results[0].decision_digest
    assert measured.depletion_witness_digest == results[0].depletion_witness_digest
    assert measured.reporting_observation_digest == results[0].reporting_observation.digest
    assert measured.native_result_digest == results[0].digest
    assert plunge.depletion_witness_digest == results[1].depletion_witness_digest
    assert plunge.native_result_digest == results[1].digest
    assert non_engaging.reason == "vertical_retract"
    assert non_engaging.pre_motion_stock_lineage == non_engaging.post_motion_stock_lineage
    assert non_engaging.native_result_digest == results[2].digest
    assert bytes(measured.digest) == hashlib.sha256(measured.canonical_bytes).digest()
    assert bytes(plunge.digest) == hashlib.sha256(plunge.canonical_bytes).digest()
    assert bytes(non_engaging.digest) == hashlib.sha256(non_engaging.canonical_bytes).digest()


def test_output_record_union_has_disjoint_chronology_fields() -> None:
    carriers, results = _native_record_inputs()
    measured = MeasuredOperationAudit.from_native(
        carriers[0],
        results[0],
        expected_native_request_digest=results[0].request_digest,
    )
    plunge = PlungeOperationAudit.from_native(
        carriers[1],
        results[1],
        expected_native_request_digest=results[1].request_digest,
    )
    non_engaging = NonEngagingOperationAudit.from_native(
        carriers[2],
        results[2],
        expected_native_request_digest=results[2].request_digest,
    )

    for record in (measured, plunge, non_engaging):
        assert not hasattr(record, "motion")
        assert not hasattr(record, "operation_digest")
    for forbidden in ("reason",):
        assert not hasattr(measured, forbidden)
        assert not hasattr(plunge, forbidden)
    for forbidden in (
        "verdict",
        "max_tea",
        "evidence_count",
        "reporting_station_count",
        "native_decision_digest",
        "reporting_observation_digest",
    ):
        assert not hasattr(plunge, forbidden)
        assert not hasattr(non_engaging, forbidden)
    assert not hasattr(non_engaging, "depletion_witness_digest")


def test_output_record_factories_reject_wrong_kinds_and_foreign_pairs() -> None:
    carriers, results = _native_record_inputs()
    _foreign_carriers, foreign_results = _native_record_inputs(source_suffix=b"foreign")
    _same_carriers, foreign_request_results = _native_record_inputs(
        stock_floor_y=-0.3
    )

    with pytest.raises(InvalidMeasuredOperationAuditError):
        MeasuredOperationAudit.from_native(
            carriers[0],
            results[1],  # type: ignore[arg-type]
            expected_native_request_digest=results[0].request_digest,
        )
    with pytest.raises(InvalidMeasuredOperationAuditError, match="operation|digest"):
        MeasuredOperationAudit.from_native(
            carriers[0],
            foreign_results[0],
            expected_native_request_digest=results[0].request_digest,
        )
    with pytest.raises(InvalidMeasuredOperationAuditError, match="request"):
        MeasuredOperationAudit.from_native(
            carriers[0],
            foreign_request_results[0],
            expected_native_request_digest=results[0].request_digest,
        )
    with pytest.raises(InvalidPlungeOperationAuditError):
        PlungeOperationAudit.from_native(
            carriers[1],
            results[0],  # type: ignore[arg-type]
            expected_native_request_digest=results[1].request_digest,
        )
    with pytest.raises(InvalidNonEngagingOperationAuditError):
        NonEngagingOperationAudit.from_native(
            carriers[2],
            results[1],  # type: ignore[arg-type]
            expected_native_request_digest=results[2].request_digest,
        )


@pytest.mark.parametrize(
    ("record_type", "error"),
    (
        (MeasuredOperationAudit, InvalidMeasuredOperationAuditError),
        (PlungeOperationAudit, InvalidPlungeOperationAuditError),
        (NonEngagingOperationAudit, InvalidNonEngagingOperationAuditError),
    ),
)
def test_output_records_are_factory_owned_and_reject_raw_construction(
    record_type: type[object],
    error: type[ValueError],
) -> None:
    with pytest.raises(error):
        record_type()

    forged = object.__new__(record_type)
    with pytest.raises(error):
        _ = forged.canonical_bytes  # type: ignore[attr-defined]


def test_output_records_reject_replace_subclass_and_post_build_tamper() -> None:
    carriers, results = _native_record_inputs()
    measured = MeasuredOperationAudit.from_native(
        carriers[0],
        results[0],
        expected_native_request_digest=results[0].request_digest,
    )
    plunge = PlungeOperationAudit.from_native(
        carriers[1],
        results[1],
        expected_native_request_digest=results[1].request_digest,
    )
    non_engaging = NonEngagingOperationAudit.from_native(
        carriers[2],
        results[2],
        expected_native_request_digest=results[2].request_digest,
    )
    records_and_errors = (
        (measured, InvalidMeasuredOperationAuditError),
        (plunge, InvalidPlungeOperationAuditError),
        (non_engaging, InvalidNonEngagingOperationAuditError),
    )

    for record, error in records_and_errors:
        with pytest.raises(error):
            replace(record, operation_index=999)

        for field in fields(record):
            tampered = copy(record)
            value = getattr(tampered, field.name)
            if isinstance(value, bytes):
                changed: object = _digest(b"tampered-field")
            elif isinstance(value, str):
                changed = f"{value}-tampered"
            elif isinstance(value, int):
                changed = value + 1
            elif isinstance(value, float):
                changed = value + 0.125
            else:
                changed = object()
            object.__setattr__(tampered, field.name, changed)
            with pytest.raises(error):
                _ = tampered.canonical_bytes
            with pytest.raises(error):
                _ = tampered.digest

    for record_type in (
        MeasuredOperationAudit,
        PlungeOperationAudit,
        NonEngagingOperationAudit,
    ):
        with pytest.raises(TypeError):

            class ForgedRecord(record_type):  # type: ignore[valid-type,misc]
                pass
