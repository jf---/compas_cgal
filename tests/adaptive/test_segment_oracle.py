import hashlib
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import DegenerateSegmentMotionError
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


HALF_PLANE_WINDOW = np.array(
    [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 10.0, 0.0],
        [0.0, 10.0, 0.0],
    ],
    dtype=np.float64,
)

CONCAVE_WINDOW = np.array(
    [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 10.0, 0.0],
        [5.0, 10.0, 0.0],
        [5.0, 5.0, 0.0],
        [0.0, 5.0, 0.0],
    ],
    dtype=np.float64,
)

DYADIC_FALSIFIER_DEPTH = 5


def _motion(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
) -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(x0, y0),
        Point2[WorldXY].build(x1, y1),
    )


def _audit(
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion,
    tool_radius: float,
    cap_chord_ratio: float,
):
    return _continuous_tea_2.audit_segment_tea_event_exact(
        stock,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        tool_radius,
        cap_chord_ratio,
    )


def _events(trace, kind: str) -> tuple[object, ...]:
    return tuple(event for event in trace.events if event.kind == kind)


def _root(
    trace,
    coefficients: tuple[str, ...],
    ordinal: int = 0,
):
    return next(root for root in trace.partition.roots if tuple(root.factor_coefficients) == coefficients and root.root_ordinal == ordinal)


def _assert_trace_binds_verified_partition(
    verdict: str,
    trace,
) -> None:
    assert trace.exact_verdict == verdict
    assert trace.decision_authority_bytes != trace.partition.canonical_bytes
    assert trace.partition.canonical_bytes in trace.decision_authority_bytes
    assert hashlib.sha256(trace.decision_authority_bytes).digest() == trace.decision_authority_digest
    assert trace.decision_authority_digest in trace.canonical_bytes
    assert trace.canonical_digest
    assert trace.motion_chart_id
    assert trace.effective_cap_bytes
    assert trace.oracle_strategy_version == "segment-event-exact-v1"
    assert trace.event_cell_count == len(trace.partition.cells)


def _assert_no_confirmed_dyadic_violation(
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion,
    tool_radius: float,
    cap_chord_ratio: float,
) -> None:
    denominator = 2**DYADIC_FALSIFIER_DEPTH
    start_x = Fraction.from_float(motion.start.x)
    start_y = Fraction.from_float(motion.start.y)
    end_x = Fraction.from_float(motion.end.x)
    end_y = Fraction.from_float(motion.end.y)
    for numerator in range(denominator + 1):
        parameter = Fraction(numerator, denominator)
        exact_x = start_x + parameter * (end_x - start_x)
        exact_y = start_y + parameter * (end_y - start_y)
        probe_x = float(exact_x)
        probe_y = float(exact_y)
        probe_exceeded = _stock_2.engagement_at(
            stock,
            probe_x,
            probe_y,
            tool_radius,
            cap_chord_ratio,
            0.0,
        )[2]
        if not probe_exceeded:
            continue
        if Fraction.from_float(probe_x) == exact_x and Fraction.from_float(probe_y) == exact_y:
            pytest.fail(
                f"certified segment has an exact dyadic violation at {numerator}/{denominator}",
            )
        exact_exceeded = _continuous_tea_2.segment_station_cap_exceeded_exact(
            stock,
            motion.start.x,
            motion.start.y,
            motion.end.x,
            motion.end.y,
            numerator,
            denominator,
            tool_radius,
            cap_chord_ratio,
        )
        assert not exact_exceeded


def test_line_half_plane_pi_equality_certifies() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    motion = _motion(2.0, 0.0, 8.0, 0.0)

    verdict, trace = _audit(
        stock,
        motion,
        1.0,
        4.0,
    )

    assert verdict == "certified"
    _assert_trace_binds_verified_partition(verdict, trace)
    _assert_no_confirmed_dyadic_violation(stock, motion, 1.0, 4.0)


def test_line_half_plane_tighter_cap_exceeds() -> None:
    verdict, trace = _audit(
        _stock_2.Stock2(HALF_PLANE_WINDOW, []),
        _motion(2.0, 0.0, 8.0, 0.0),
        1.0,
        2.0,
    )

    assert verdict == "cap_exceeded"
    _assert_trace_binds_verified_partition(verdict, trace)


def test_fully_clear_and_fully_material_rims_are_not_conflated() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    clear_motion = _motion(2.0, -2.0, 8.0, -2.0)
    material_motion = _motion(2.0, 5.0, 8.0, 5.0)

    clear_verdict, clear_trace = _audit(stock, clear_motion, 0.5, 4.0)
    material_verdict, material_trace = _audit(stock, material_motion, 0.5, 4.0)

    assert clear_verdict == "certified"
    assert material_verdict == "cap_exceeded"
    assert clear_trace.whole_rim_disposition == "clear"
    assert material_trace.whole_rim_disposition == "material"
    _assert_trace_binds_verified_partition(clear_verdict, clear_trace)
    _assert_trace_binds_verified_partition(material_verdict, material_trace)
    _assert_no_confirmed_dyadic_violation(
        stock,
        clear_motion,
        0.5,
        4.0,
    )


def test_swept_prefix_segment_certifies_only_clear_start_and_pi_cap() -> None:
    """Prove a slot-opening translation from its exact moving sweep prefix.

    Frozen stock sees material around more than a semicircle after the cutter
    leaves its clear start disk. The physical motion has already removed every
    strictly trailing rim point, so the native theorem may certify the closed
    forward semicircle at an exact pi cap, while refusing an uncleared start or
    a tighter cap.
    """
    cleared = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    cleared.subtract_disk(2.0, 5.0, 1.0)
    motion = _motion(2.0, 5.0, 8.0, 5.0)

    audit = _continuous_tea_2.audit_swept_prefix_segment_tea_exact(
        cleared,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        1.0,
        4.0,
    )
    uncleared = _continuous_tea_2.audit_swept_prefix_segment_tea_exact(
        _stock_2.Stock2(HALF_PLANE_WINDOW, []),
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        1.0,
        4.0,
    )
    tighter = _continuous_tea_2.audit_swept_prefix_segment_tea_exact(
        cleared,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        1.0,
        2.0,
    )

    assert type(audit) is _continuous_tea_2.SweptPrefixSegmentTeaAudit2
    assert audit.exact_verdict == "certified"
    assert audit.is_self_consistent
    assert hashlib.sha256(audit.canonical_bytes).digest() == audit.canonical_digest
    assert audit.strategy_version == b"swept-prefix-segment-tea-exact-v1"
    assert audit.theorem_version == b"translation-swept-prefix-forward-semicircle-v1"
    assert audit.source_canonical_bytes
    assert len(audit.stock_boundary_digest) == hashlib.sha256().digest_size
    assert audit.motion_stratum_count == 2
    assert uncleared.exact_verdict == tighter.exact_verdict == "unresolved"
    with pytest.raises(TypeError):
        _continuous_tea_2.SweptPrefixSegmentTeaAudit2()


@pytest.mark.parametrize(
    ("reverse", "disposition"),
    [
        (False, "birth"),
        (True, "death"),
    ],
)
def test_line_tangency_root_owns_run_birth_and_death(
    reverse: bool,
    disposition: str,
) -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    endpoints = ((5.0, -2.0), (5.0, 0.0))
    if reverse:
        endpoints = tuple(reversed(endpoints))
    motion = _motion(*endpoints[0], *endpoints[1])

    verdict, trace = _audit(stock, motion, 1.0, 4.0)

    assert verdict == "certified"
    tangent = _events(trace, "tangent")
    assert len(tangent) == 1
    parameter_root = _root(trace, ("-1", "2"))
    assert parameter_root.multiplicity == 1
    assert tangent[0].root_id == parameter_root.root_id
    assert tangent[0].multiplicity == 2
    assert tangent[0].disposition == disposition
    _assert_no_confirmed_dyadic_violation(stock, motion, 1.0, 4.0)


def test_circle_circle_tangency_root_is_not_omitted() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    stock.subtract_disk(5.0, 1.5, 1.0)
    motion = _motion(3.0, -0.5, 7.0, -0.5)

    verdict, trace = _audit(stock, motion, 1.0, 4.0)

    assert verdict == "certified"
    tangent = _events(trace, "tangent")
    assert len(tangent) == 1
    parameter_root = _root(trace, ("-1", "2"))
    assert parameter_root.multiplicity == 2
    assert tangent[0].root_id == parameter_root.root_id
    assert tangent[0].multiplicity == 2
    _assert_no_confirmed_dyadic_violation(stock, motion, 1.0, 4.0)


def test_cutter_rim_vertex_fibre_is_not_omitted() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    motion = _motion(-2.0, 0.0, 0.0, 0.0)

    verdict, trace = _audit(stock, motion, 1.0, 4.0)

    assert verdict == "certified"
    vertex = _events(trace, "vertex")
    assert len(vertex) == 1
    assert vertex[0].root_id == _root(trace, ("-1", "2")).root_id
    assert vertex[0].feature_ids
    _assert_no_confirmed_dyadic_violation(stock, motion, 1.0, 4.0)


def test_concave_vertex_run_split_and_merge_roots_are_not_omitted() -> None:
    stock = _stock_2.Stock2(CONCAVE_WINDOW, [])
    motion = _motion(4.0, 4.0, 6.0, 6.0)

    verdict, trace = _audit(stock, motion, 1.0, 4.0)

    topology = _events(trace, "endpoint-order")
    assert verdict == "cap_exceeded"
    assert [event.disposition for event in topology] == ["split", "merge"]
    assert all(event.root_id and event.multiplicity == 2 and len(event.feature_ids) == 2 and len(event.branch_ids) == 2 for event in topology)
    assert _continuous_tea_2.segment_station_cap_exceeded_exact(
        stock,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        1,
        2,
        1.0,
        4.0,
    )


def test_coincident_tool_disk_boundary_is_never_empty_by_omission() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    stock.subtract_disk(5.0, 5.0, 1.0)

    verdict, trace = _audit(
        stock,
        _motion(5.0, 5.0, 5.0 + 2.0**-20, 5.0),
        1.0,
        4.0,
    )

    assert verdict in {"cap_exceeded", "unresolved"}
    assert _events(trace, "overlap")
    assert trace.whole_rim_disposition != "clear"
    _assert_trace_binds_verified_partition(verdict, trace)


def test_cap_equality_at_algebraic_root_is_not_omitted() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    motion = _motion(5.0, -1.0, 5.0, 0.0)

    verdict, trace = _audit(stock, motion, 1.0, 2.0)

    assert verdict == "cap_exceeded"
    cap_root = _root(trace, ("1", "-4", "2"), ordinal=0)
    cap_events = _events(trace, "cap")
    assert len(cap_events) == 1
    assert cap_events[0].root_id == cap_root.root_id
    assert cap_events[0].disposition == "equal"
    _assert_trace_binds_verified_partition(verdict, trace)
    partition = _continuous_tea_2.construct_segment_event_partition(
        stock,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        1.0,
        2.0,
    )
    mutated = _continuous_tea_2.mutate_segment_event_partition(
        partition,
        "delete-projection-cap-crossing",
    )
    assert (
        _continuous_tea_2.verify_segment_event_partition(
            stock,
            partition.source,
            mutated,
        ).verdict.name
        == "UNRESOLVED_DEGENERACY"
    )


def test_zero_length_motion_is_rejected_before_native_dispatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dispatched = False

    def record_dispatch(*_args: object) -> None:
        nonlocal dispatched
        dispatched = True

    monkeypatch.setattr(
        _continuous_tea_2,
        "audit_segment_tea_event_exact",
        record_dispatch,
        raising=False,
    )
    point = Point2[WorldXY].build(1.0, 2.0)

    with pytest.raises(DegenerateSegmentMotionError):
        _audit(
            _stock_2.Stock2(HALF_PLANE_WINDOW, []),
            ExactSegmentMotion.build(point, point),
            1.0,
            4.0,
        )

    assert not dispatched


def test_trace_identity_is_independent_of_boundary_insertion_order() -> None:
    forward = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    reverse = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    disks = (
        (4.5, 5.0, 1.0),
        (5.5, 5.0, 1.0),
    )
    for disk in disks:
        forward.subtract_disk(*disk)
    for disk in reversed(disks):
        reverse.subtract_disk(*disk)
    motion = _motion(5.0, 3.0, 5.0, 7.0)

    forward_verdict, forward_trace = _audit(forward, motion, 0.5, 4.0)
    reverse_verdict, reverse_trace = _audit(reverse, motion, 0.5, 4.0)

    assert forward_verdict == reverse_verdict
    assert forward_trace.canonical_bytes == reverse_trace.canonical_bytes
    assert forward_trace.canonical_digest == reverse_trace.canonical_digest
    assert [(event.root_id, event.kind, event.feature_ids, event.branch_ids) for event in forward_trace.events] == [
        (event.root_id, event.kind, event.feature_ids, event.branch_ids) for event in reverse_trace.events
    ]
