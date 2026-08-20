from __future__ import annotations

import hashlib
import math
from collections.abc import Iterator

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive import motion_oracle_cache
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
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
TOOL_RADIUS_MM = 0.5


@pytest.fixture(autouse=True)
def _isolated_memo() -> Iterator[None]:
    """Measure hits and misses against an empty memo in both directions."""
    motion_oracle_cache.clear_native_motion_audit_cache()
    yield
    motion_oracle_cache.clear_native_motion_audit_cache()


def _certifier() -> MotionCertifier:
    return MotionCertifier.build(
        stock=Stock2Area(_stock_2.Stock2(SQUARE, []), ()),
        tool_radius=ToolRadius.build(TOOL_RADIUS_MM),
    )


def _clear_segment() -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(2.0, -2.0),
        Point2[WorldXY].build(8.0, -2.0),
    )


def _circle() -> ExactCircleMotion:
    return ExactCircleMotion.build(
        Point2[WorldXY].build(5.0, 5.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )


def _count_segment_oracle(monkeypatch: pytest.MonkeyPatch) -> list[tuple[object, ...]]:
    """Route every segment dispatch through one counting native delegate."""
    native = _continuous_tea_2.audit_segment_tea_event_exact
    dispatches: list[tuple[object, ...]] = []

    def counting(*arguments: object) -> object:
        dispatches.append(arguments)
        return native(*arguments)

    monkeypatch.setattr(
        _continuous_tea_2,
        "audit_segment_tea_event_exact",
        counting,
    )
    return dispatches


def test_repeat_certification_consults_the_native_oracle_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Prove the measured duplicate certification pays the oracle once.

    The two certifiers are distinct instances built from independently
    constructed stock. They share only content identity, which is exactly the
    duplication observed on the route-retrace path.
    """
    dispatches = _count_segment_oracle(monkeypatch)
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)

    first = _certifier().certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )
    second = _certifier().certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    assert len(dispatches) == 1
    assert first.verdict == second.verdict == "certified"
    assert first.strategy_identity == second.strategy_identity
    assert first.event_trace_digest == second.event_trace_digest
    assert first.event_cell_count == second.event_cell_count
    assert first.canonical_bytes == second.canonical_bytes
    assert first.digest == second.digest


def test_memoized_certification_equals_its_uncached_native_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Hold the memo to byte equality with a direct native audit."""
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)
    verdict, trace = _continuous_tea_2.audit_segment_tea_event_exact(
        _stock_2.Stock2(SQUARE, []),
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        TOOL_RADIUS_MM,
        cap.chord_ratio,
    )
    dispatches = _count_segment_oracle(monkeypatch)

    witnesses = tuple(
        _certifier().certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=motion,
            user_cap=cap,
            effective_cap=cap,
        )
        for _ in range(3)
    )

    assert verdict == "certified"
    assert len(dispatches) == 1
    for witness in witnesses:
        assert witness.verdict == "certified"
        assert witness.strategy_identity == trace.oracle_strategy_version.encode()
        assert witness.event_trace_digest == trace.canonical_digest
        assert witness.event_cell_count == trace.event_cell_count


def test_operation_ordinal_is_never_memoized() -> None:
    """Keep the caller's ordinal outside the memoized native result."""
    certifier = _certifier()
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)

    first = certifier.certify(
        operation_index=11,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )
    second = certifier.certify(
        operation_index=12,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    assert first.operation_index == 11
    assert second.operation_index == 12
    assert first.event_trace_digest == second.event_trace_digest
    assert first.canonical_bytes != second.canonical_bytes
    assert first.digest != second.digest


def test_foreign_stock_lineage_never_reuses_a_memoized_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bind every memo entry to the exact stock identity that produced it."""
    dispatches = _count_segment_oracle(monkeypatch)
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)
    owned = _certifier()
    owned.certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    foreign = MotionCertifier(
        _stock_2.Stock2(SQUARE, []),
        ToolRadius.build(TOOL_RADIUS_MM),
        hashlib.sha256(b"foreign-stock-lineage").digest(),
        owned.canonical_boundary_digest,
    )
    foreign.certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    assert foreign.stock_lineage_digest != owned.stock_lineage_digest
    assert len(dispatches) == 2


def test_foreign_stock_boundary_never_reuses_a_memoized_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Refuse a memo hit when the exact native boundary content differs."""
    dispatches = _count_segment_oracle(monkeypatch)
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)
    owned = _certifier()
    owned.certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    foreign = MotionCertifier(
        _stock_2.Stock2(SQUARE, []),
        ToolRadius.build(TOOL_RADIUS_MM),
        owned.stock_lineage_digest,
        hashlib.sha256(b"foreign-stock-boundary").digest(),
    )
    foreign.certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    assert len(dispatches) == 2


def test_tool_radius_and_effective_cap_each_discriminate_the_memo(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Complete the key over every native argument the oracle consumes."""
    dispatches = _count_segment_oracle(monkeypatch)
    motion = _clear_segment()
    user_cap = EngagementCap.build(math.pi)
    stock_area = Stock2Area(_stock_2.Stock2(SQUARE, []), ())

    for tool_radius, effective_cap in (
        (TOOL_RADIUS_MM, user_cap),
        (TOOL_RADIUS_MM, EngagementCap.build(math.pi / 2.0)),
        (TOOL_RADIUS_MM * 2.0, user_cap),
    ):
        MotionCertifier.build(
            stock=stock_area,
            tool_radius=ToolRadius.build(tool_radius),
        ).certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=motion,
            user_cap=user_cap,
            effective_cap=effective_cap,
        )

    assert len(dispatches) == 3


def test_segment_and_circle_motions_never_share_a_memo_entry() -> None:
    """Separate the two native entry points inside one memoized stock.

    A returned `cap_exceeded` verdict is a native result, not a native failure,
    so the circle audit is memoized exactly like the certified segment audit and
    the caller still observes its named cap exceedance on every call.
    """
    certifier = _certifier()
    cap = EngagementCap.build(math.pi)

    for ordinal in range(2):
        link = certifier.certify(
            operation_index=ordinal,
            operation_kind=OperationType.LINK,
            motion=_clear_segment(),
            user_cap=cap,
            effective_cap=cap,
        )
        assert link.verdict == "certified"
        with pytest.raises(EngagementCapExceededError, match="exceeds"):
            certifier.certify(
                operation_index=ordinal,
                operation_kind=OperationType.CUT,
                motion=_circle(),
                user_cap=cap,
                effective_cap=cap,
            )

    assert motion_oracle_cache.native_motion_audit_cache_size() == 2


def test_memo_never_grows_past_its_declared_capacity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bound retained native traces across an arbitrarily long run."""
    monkeypatch.setattr(motion_oracle_cache, "MOTION_ORACLE_CACHE_CAPACITY", 2)
    certifier = _certifier()
    cap = EngagementCap.build(math.pi)

    for offset in range(4):
        certifier.certify(
            operation_index=offset,
            operation_kind=OperationType.LINK,
            motion=ExactSegmentMotion.build(
                Point2[WorldXY].build(2.0, -2.0 - float(offset)),
                Point2[WorldXY].build(8.0, -2.0 - float(offset)),
            ),
            user_cap=cap,
            effective_cap=cap,
        )
        assert motion_oracle_cache.native_motion_audit_cache_size() <= 2


def test_a_replaced_native_oracle_is_never_answered_from_the_memo(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep a memo entry the property of the exact oracle that produced it.

    A memo that outlived its deciding authority would silently answer for a
    different oracle. Consumers that install a substitute oracle must observe
    that substitute, not a result the real native oracle returned earlier.
    """
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)
    _certifier().certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )
    substitute_calls = 0

    def substitute(*_arguments: object) -> tuple[str, object]:
        nonlocal substitute_calls
        substitute_calls += 1
        raise _continuous_tea_2.IncompleteSegmentPartitionError(
            "substitute oracle refuses this motion",
        )

    monkeypatch.setattr(
        _continuous_tea_2,
        "audit_segment_tea_event_exact",
        substitute,
    )

    with pytest.raises(
        UnresolvedMotionEventError,
        match="IncompleteSegmentPartitionError",
    ):
        _certifier().certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=motion,
            user_cap=cap,
            effective_cap=cap,
        )

    assert substitute_calls == 1
