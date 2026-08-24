from __future__ import annotations

import hashlib
import math

import numpy as np
import pytest

from compas_cgal import _stock_2


def _digest(label: str) -> bytes:
    return hashlib.sha256(f"native-replay:{label}".encode()).digest()


def _boundary() -> np.ndarray:
    return np.array(
        ((-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)),
        dtype=np.float64,
    )


def _policy() -> _stock_2.AuditPolicy2:
    cap = math.pi / 2.0
    return _stock_2.build_audit_policy(
        0.5,
        cap,
        _stock_2.cap_chord_ratio(cap),
        0.0625,
        4096,
    )


def _limits() -> _stock_2.AuditDecisionLimits2:
    return _stock_2.build_audit_decision_limits(0.015625, 8, 256)


def _motions() -> tuple[object, ...]:
    return (
        _stock_2.classify_audit_line((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.0, 5.0, "cut"),
        _stock_2.classify_audit_arc(
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            1.0,
            0.0,
            math.pi / 2.0,
            False,
            0.0,
            5.0,
            "cut",
        ),
        _stock_2.classify_audit_circle(
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            1.0,
            False,
            0.0,
            5.0,
            "cut",
        ),
        _stock_2.classify_audit_line((3.0, 3.0, 5.0), (3.0, 3.0, 0.0), 0.0, 5.0, "plunge"),
        _stock_2.classify_audit_line((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), 0.0, 5.0, "retract"),
        _stock_2.classify_audit_line((-1.0, 0.0, 5.0), (1.0, 0.0, 5.0), 0.0, 5.0, "link"),
    )


def _assert_read_only(instance: object, names: tuple[str, ...]) -> None:
    for name in names:
        value = getattr(instance, name)
        with pytest.raises(AttributeError):
            setattr(instance, name, value)


def test_native_replay_is_opaque_and_executes_the_closed_motion_domain() -> None:
    boundary = _boundary()
    policy = _policy()
    limits = _limits()
    motions = _motions()
    operation_digests = tuple(_digest(str(index)) for index in range(6))
    request = _stock_2.build_audit_native_request_identity(boundary, [], policy, limits, motions)
    replay = _stock_2.begin_audit_replay(
        boundary,
        [],
        _digest("input"),
        request,
        policy,
        limits,
        operation_digests,
    )

    segment = _stock_2.audit_deplete_segment(replay, motions[0], operation_digests[0])
    arc = _stock_2.audit_deplete_arc(replay, motions[1], operation_digests[1])
    circle = _stock_2.audit_deplete_circle(replay, motions[2], operation_digests[2])
    plunge = _stock_2.deplete_audit_plunge(replay, motions[3], operation_digests[3])
    retract = _stock_2.record_audit_retract(replay, motions[4], operation_digests[4])
    clearance = _stock_2.record_audit_clearance(replay, motions[5], operation_digests[5])
    completion = _stock_2.finish_audit_replay(replay)

    for result, operation_digest in zip(
        (segment, arc, circle, plunge, retract, clearance),
        operation_digests,
        strict=True,
    ):
        assert result.authenticated_operation_digest == operation_digest
        assert len(result.digest) == 32
        assert len(result.pre_lineage) == 32
        assert len(result.post_lineage) == 32
        assert all(
            isinstance(getattr(result, name), bytes)
            for name in (
                "authenticated_operation_digest",
                "digest",
                "pre_lineage",
                "post_lineage",
            )
        )
    for result in (segment, arc, circle):
        assert result.verdict in {"certified", "cap_exceeded", "unresolved"}
        assert result.evidence_count > 0
        assert len(result.decision_digest) == 32
        assert len(result.depletion_witness_digest) == 32
    assert len(plunge.depletion_witness_digest) == 32
    assert retract.reason == "vertical_retract"
    assert clearance.reason == "clearance_transport"
    assert retract.pre_lineage == retract.post_lineage
    assert clearance.pre_lineage == clearance.post_lineage
    assert completion.request_digest == request.digest
    assert completion.operation_count == 6
    assert len(completion.terminal_lineage) == 32
    assert len(completion.digest) == 32
    assert all(
        isinstance(value, bytes)
        for value in (
            request.canonical_bytes,
            request.digest,
            limits.canonical_bytes,
            segment.decision_digest,
            segment.depletion_witness_digest,
            plunge.depletion_witness_digest,
            completion.request_digest,
            completion.terminal_lineage,
            completion.digest,
        )
    )

    _assert_read_only(request, ("canonical_bytes", "digest"))
    _assert_read_only(limits, ("canonical_bytes",))
    _assert_read_only(
        segment,
        (
            "verdict",
            "evidence_count",
            "authenticated_operation_digest",
            "decision_digest",
            "depletion_witness_digest",
            "pre_lineage",
            "post_lineage",
            "digest",
        ),
    )
    _assert_read_only(
        plunge,
        (
            "authenticated_operation_digest",
            "depletion_witness_digest",
            "pre_lineage",
            "post_lineage",
            "digest",
        ),
    )
    _assert_read_only(
        retract,
        (
            "reason",
            "authenticated_operation_digest",
            "pre_lineage",
            "post_lineage",
            "digest",
        ),
    )
    _assert_read_only(
        completion,
        ("request_digest", "operation_count", "terminal_lineage", "digest"),
    )

    for opaque_type in (
        _stock_2.AuditReplay2,
        _stock_2.AuditLateralResult2,
        _stock_2.AuditPlungeResult2,
        _stock_2.AuditNonEngagingResult2,
        _stock_2.AuditReplayCompletion2,
    ):
        with pytest.raises(TypeError):
            opaque_type()
    assert not hasattr(_stock_2, "AuditDepletionWitness2")
    with pytest.raises(AttributeError):
        segment.digest = b"forged"
    assert not hasattr(replay, "stock")
    assert not hasattr(segment, "max_tea")
    assert {name for name in dir(_stock_2) if name.startswith(("audit_deplete_", "deplete_audit_", "record_audit_"))} == {
        "audit_deplete_segment",
        "audit_deplete_arc",
        "audit_deplete_circle",
        "deplete_audit_plunge",
        "record_audit_retract",
        "record_audit_clearance",
    }
    assert not hasattr(_stock_2, "audit_deplete_motion")
    assert not hasattr(_stock_2, "audit_replay_operation")


def test_native_replay_names_ingress_identity_and_lifecycle_failures() -> None:
    boundary = _boundary()
    policy = _policy()
    limits = _limits()
    motion = _motions()[0]
    operation_digest = _digest("operation")
    request = _stock_2.build_audit_native_request_identity(boundary, [], policy, limits, (motion,))

    for operation_digests in ((), (operation_digest, _digest("surplus"))):
        with pytest.raises(_stock_2.AuditReplayCardinalityError):
            _stock_2.begin_audit_replay(
                boundary,
                [],
                _digest("input"),
                request,
                policy,
                limits,
                operation_digests,
            )

    with pytest.raises(_stock_2.AuditDigestSizeError):
        _stock_2.begin_audit_replay(boundary, [], b"short", request, policy, limits, (operation_digest,))
    with pytest.raises(_stock_2.AuditReplayRequestIdentityError):
        _stock_2.begin_audit_replay(
            boundary,
            [],
            _digest("input"),
            request,
            policy,
            _stock_2.build_audit_decision_limits(0.0078125, 8, 256),
            (operation_digest,),
        )

    replay = _stock_2.begin_audit_replay(
        boundary,
        [],
        _digest("input"),
        request,
        policy,
        limits,
        (operation_digest,),
    )
    with pytest.raises(_stock_2.AuditReplayIncompleteError):
        _stock_2.finish_audit_replay(replay)
    with pytest.raises(_stock_2.AuditReplayOperationIdentityError):
        _stock_2.audit_deplete_segment(replay, motion, _digest("foreign"))
    with pytest.raises(TypeError):
        _stock_2.audit_deplete_segment(replay, _motions()[2], operation_digest)
    with pytest.raises(TypeError):
        _stock_2.audit_deplete_segment(replay, (0.0, 0.0), operation_digest)
    with pytest.raises(TypeError):
        _stock_2.audit_deplete_segment(_stock_2.Stock2(boundary, []), motion, operation_digest)
    _stock_2.audit_deplete_segment(replay, motion, operation_digest)
    with pytest.raises(_stock_2.AuditReplayOperationExhaustedError):
        _stock_2.audit_deplete_segment(replay, motion, operation_digest)
    _stock_2.finish_audit_replay(replay)
    with pytest.raises(_stock_2.AuditReplayFinalizedError):
        _stock_2.finish_audit_replay(replay)
    with pytest.raises(_stock_2.AuditReplayFinalizedError):
        _stock_2.audit_deplete_segment(replay, motion, operation_digest)


def test_native_replay_rejects_out_of_order_and_wrong_typed_motion() -> None:
    boundary = _boundary()
    policy = _policy()
    limits = _limits()
    first = _stock_2.classify_audit_line((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.0, 5.0, "cut")
    second = _stock_2.classify_audit_line((-1.0, 1.0, 0.0), (1.0, 1.0, 0.0), 0.0, 5.0, "cut")
    first_digest = _digest("first")
    second_digest = _digest("second")
    request = _stock_2.build_audit_native_request_identity(boundary, [], policy, limits, (first, second))
    replay = _stock_2.begin_audit_replay(
        boundary,
        [],
        _digest("ordered-input"),
        request,
        policy,
        limits,
        (first_digest, second_digest),
    )

    with pytest.raises(_stock_2.AuditReplayOperationIdentityError):
        _stock_2.audit_deplete_segment(replay, second, second_digest)
    with pytest.raises(_stock_2.AuditReplayMotionIdentityError):
        _stock_2.audit_deplete_segment(replay, second, first_digest)
    _stock_2.audit_deplete_segment(replay, first, first_digest)
    _stock_2.audit_deplete_segment(replay, second, second_digest)
