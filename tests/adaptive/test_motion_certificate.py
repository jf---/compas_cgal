import hashlib
import math
from dataclasses import replace

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import InvalidMotionCertificateError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import MotionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.toolpath import OperationType

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)


def _certifier() -> MotionCertifier:
    return MotionCertifier.build(
        stock=Stock2Area(_stock_2.Stock2(SQUARE, []), ()),
        tool_radius=ToolRadius.build(0.5),
    )


def _custom_certifier(
    stock: _stock_2.Stock2,
    tool_radius: float,
) -> MotionCertifier:
    return MotionCertifier.build(
        stock=Stock2Area(stock, ()),
        tool_radius=ToolRadius.build(tool_radius),
    )


def _circle() -> ExactCircleMotion:
    return ExactCircleMotion.build(
        Point2[WorldXY].build(5.0, 5.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )


def _slotting_segment() -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(4.0, 5.0),
        Point2[WorldXY].build(6.0, 5.0),
    )


def _clear_segment() -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(2.0, -2.0),
        Point2[WorldXY].build(8.0, -2.0),
    )


def test_certifier_rejects_invalid_local_contract_before_native_dispatch() -> None:
    certifier = _certifier()

    with pytest.raises(InvalidMotionCertificateError, match="nonnegative"):
        certifier.certify(
            operation_index=-1,
            operation_kind=OperationType.CUT,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="exceeds"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.CUT,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="incompatible"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="lateral"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.PLUNGE,
            motion=_slotting_segment(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )


def test_predeplete_virgin_slotting_is_named_cap_exceedance_and_const() -> None:
    certifier = _certifier()
    boundary_digest = certifier.canonical_boundary_digest
    lineage_digest = certifier.stock_lineage_digest
    contains = certifier.contains(5.0, 5.0)

    with pytest.raises(EngagementCapExceededError, match="exceeds"):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=_slotting_segment(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )

    assert certifier.canonical_boundary_digest == boundary_digest
    assert certifier.stock_lineage_digest == lineage_digest
    assert certifier.contains(5.0, 5.0) is contains


def test_circle_cap_exceedance_is_named() -> None:
    with pytest.raises(EngagementCapExceededError, match="exceeds"):
        _certifier().certify(
            operation_index=1,
            operation_kind=OperationType.CUT,
            motion=_circle(),
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )


def test_real_native_unresolved_is_named_and_const() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(0.0, 5.0, 2.0)
    certifier = _custom_certifier(stock, 2.0)
    motion = ExactCircleMotion.build(
        Point2[WorldXY].build(1.0, 5.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )
    boundary_digest = certifier.canonical_boundary_digest
    lineage_digest = certifier.stock_lineage_digest
    contains = certifier.contains(5.0, 5.0)

    with pytest.raises(UnresolvedMotionEventError, match="unresolved"):
        certifier.certify(
            operation_index=2,
            operation_kind=OperationType.CUT,
            motion=motion,
            user_cap=EngagementCap.build(math.pi),
            effective_cap=EngagementCap.build(math.pi),
        )

    assert certifier.canonical_boundary_digest == boundary_digest
    assert certifier.stock_lineage_digest == lineage_digest
    assert certifier.contains(5.0, 5.0) is contains


def test_segment_certification_binds_exact_native_trace() -> None:
    certifier = _certifier()
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)

    witness = certifier.certify(
        operation_index=3,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )
    verdict, trace = _continuous_tea_2.audit_segment_tea_event_exact(
        _stock_2.Stock2(SQUARE, []),
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        0.5,
        cap.chord_ratio,
    )

    assert verdict == "certified"
    assert witness.operation_index == 3
    assert witness.operation_kind is OperationType.LINK
    assert witness.motion is motion
    assert witness.user_cap_bytes == cap.chord_ratio_bytes
    assert witness.effective_cap_bytes == cap.chord_ratio_bytes
    assert witness.strategy_identity == trace.oracle_strategy_version.encode()
    assert witness.stock_lineage_digest == certifier.stock_lineage_digest
    assert witness.event_trace_digest == trace.canonical_digest
    assert witness.event_cell_count == trace.event_cell_count
    assert witness.verdict == "certified"
    assert witness.unresolved_count == 0


def test_circle_certification_binds_exact_native_trace() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 100.0)
    certifier = _custom_certifier(stock, 0.5)
    motion = _circle()
    cap = EngagementCap.build(math.pi)

    witness = certifier.certify(
        operation_index=4,
        operation_kind=OperationType.CUT,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )
    verdict, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        motion.center.x,
        motion.center.y,
        motion.phase_vector.x,
        motion.phase_vector.y,
        motion.clockwise,
        0.5,
        cap.chord_ratio,
    )

    assert verdict == "certified"
    assert witness.operation_index == 4
    assert witness.operation_kind is OperationType.CUT
    assert witness.motion is motion
    assert witness.strategy_identity == trace.oracle_strategy_version.encode()
    assert witness.event_trace_digest == trace.canonical_digest
    assert witness.event_cell_count == trace.event_cell_count


def test_certifier_rejects_native_verdict_trace_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 100.0)
    motion = _circle()
    cap = EngagementCap.build(math.pi)
    _, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        motion.center.x,
        motion.center.y,
        motion.phase_vector.x,
        motion.phase_vector.y,
        motion.clockwise,
        0.5,
        cap.chord_ratio,
    )
    monkeypatch.setattr(
        _continuous_tea_2,
        "audit_full_circle_tea_event_exact",
        lambda *_args: ("unresolved", trace),
    )

    with pytest.raises(InvalidMotionCertificateError, match="verdict"):
        _custom_certifier(stock, 0.5).certify(
            operation_index=5,
            operation_kind=OperationType.CUT,
            motion=motion,
            user_cap=cap,
            effective_cap=cap,
        )


def test_witness_type_is_frozen() -> None:
    assert MotionWitness.__dataclass_params__.frozen


def test_witness_rejects_cross_wired_kind_and_cap_bytes() -> None:
    cap = EngagementCap.build(math.pi)
    witness = _certifier().certify(
        operation_index=6,
        operation_kind=OperationType.LINK,
        motion=_clear_segment(),
        user_cap=cap,
        effective_cap=cap,
    )

    with pytest.raises(InvalidMotionCertificateError, match="incompatible"):
        replace(witness, operation_kind=OperationType.CUT)
    with pytest.raises(InvalidMotionCertificateError, match="effective cap"):
        replace(
            witness,
            user_cap_bytes=EngagementCap.build(math.pi / 2.0).chord_ratio_bytes,
        )


def test_motion_witness_is_one_canonical_content_addressed_record() -> None:
    """Make the accepted safety proof reusable by replay and artifact identity.

    Replay must order the exact witnesses returned by the real certifier.  A
    frozen Python object without canonical bytes would force each downstream
    consumer to invent its own serialization and could make the same proof hash
    differently across the replay and artifact boundaries.
    """
    cap = EngagementCap.build(math.pi)
    witness = _certifier().certify(
        operation_index=7,
        operation_kind=OperationType.LINK,
        motion=_clear_segment(),
        user_cap=cap,
        effective_cap=cap,
    )
    repeated = _certifier().certify(
        operation_index=7,
        operation_kind=OperationType.LINK,
        motion=_clear_segment(),
        user_cap=cap,
        effective_cap=cap,
    )

    assert require_canonical_record(witness.canonical_bytes) == witness.canonical_bytes
    assert repeated.canonical_bytes == witness.canonical_bytes
    assert repeated.digest == witness.digest
    assert witness.digest == hashlib.sha256(witness.canonical_bytes).digest()

    variants = (
        replace(witness, operation_index=8),
        replace(witness, effective_cap_bytes=EngagementCap.build(math.pi / 2.0).chord_ratio_bytes),
        replace(witness, strategy_identity=b"alternate-exact-oracle-v1"),
        replace(
            witness,
            stock_lineage_digest=hashlib.sha256(b"alternate-stock-lineage").digest(),
        ),
        replace(
            witness,
            event_trace_digest=hashlib.sha256(b"alternate-event-trace").digest(),
        ),
        replace(witness, event_cell_count=witness.event_cell_count + 1),
    )
    assert all(variant.canonical_bytes != witness.canonical_bytes for variant in variants)
    assert len({variant.digest for variant in variants}) == len(variants)
