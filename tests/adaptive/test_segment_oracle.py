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


def _certificate_bytes(
    certificate: object,
    **replacements: object,
) -> bytes:
    fields = {
        "build_evidence": certificate.build_evidence,
        "charts": certificate.charts,
        "projections": certificate.projections,
        "roots": certificate.roots,
        "cells": certificate.cells,
        "fibres": certificate.fibres,
        "overlaps": certificate.overlaps,
        "seams": certificate.seams,
        "source_kind": certificate.source_kind,
        "source_payload": certificate.source_payload,
    }
    fields.update(replacements)
    return _continuous_tea_2.encode_event_partition_certificate(**fields)


def _assert_trace_binds_verified_partition(
    verdict: str,
    trace,
) -> None:
    assert trace.exact_verdict == verdict
    assert trace.partition.canonical_digest == hashlib.sha256(
        trace.partition.canonical_bytes,
    ).digest()
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


@pytest.mark.parametrize(
    ("reverse", "disposition"),
    [
        (False, "merge"),
        (True, "split"),
    ],
)
def test_slotted_run_merge_and_split_root_is_not_omitted(
    reverse: bool,
    disposition: str,
) -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    stock.subtract_capsule(5.0, 4.3, 5.0, 4.53, 0.05)
    endpoints = ((4.7, 5.0), (4.725, 5.0))
    if reverse:
        endpoints = tuple(reversed(endpoints))
    motion = _motion(*endpoints[0], *endpoints[1])

    verdict, trace = _audit(stock, motion, 0.5, 4.0)

    assert verdict == "certified"
    topology = _events(trace, "endpoint-order")
    assert len(topology) == 1
    assert topology[0].root_id
    assert topology[0].disposition == disposition
    _assert_no_confirmed_dyadic_violation(stock, motion, 0.5, 4.0)


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
    root_only_mutation = _certificate_bytes(
        trace.partition,
        roots=tuple(root for root in trace.partition.roots if root.root_id != cap_root.root_id),
    )
    with pytest.raises(_continuous_tea_2.EventPartitionVerificationError):
        _continuous_tea_2.verify_event_partition(root_only_mutation)


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


def _shared_trace_certificate(
    cap_chord_ratio: float,
) -> _continuous_tea_2.EventPartitionCertificate2:
    cap = Fraction.from_float(cap_chord_ratio)
    return _continuous_tea_2.partition_cap_crossings(
        ("0", "1"),
        ("1",),
        ("0", "3"),
        ("1",),
        str(cap.numerator),
        str(cap.denominator),
        _continuous_tea_2.PartitionEvent2(
            "cap",
            b"feature-cap",
            b"support-cap",
            b"trim-cap",
            b"vertex-cap",
            b"branch-cap",
            "equal",
        ),
    )


def test_shared_event_trace_requires_verification_and_canonicalizes_events() -> None:
    cap_chord_ratio = 64.0 / 65.0
    certificate = _shared_trace_certificate(cap_chord_ratio)
    first_root = certificate.roots[0].root_id
    event = _continuous_tea_2.EventTraceEvent2(
        first_root,
        first_root,
        "endpoint-order",
        (b"feature-b", b"feature-a"),
        (b"branch-b", b"branch-a"),
        1,
        "merge",
        0,
    )
    reversed_event = _continuous_tea_2.EventTraceEvent2(
        first_root,
        first_root,
        "endpoint-order",
        (b"feature-a", b"feature-b"),
        (b"branch-a", b"branch-b"),
        1,
        "merge",
        0,
    )

    trace = _continuous_tea_2.build_event_trace(
        certificate,
        "segment-linear-0-1-v1",
        certificate.source_payload,
        certificate.canonical_digest,
        _continuous_tea_2.ContinuousTeaVerdict.CERTIFIED,
        "partial",
        "segment-event-exact-v1",
        (event,),
    )
    reordered = _continuous_tea_2.build_event_trace(
        certificate,
        "segment-linear-0-1-v1",
        certificate.source_payload,
        certificate.canonical_digest,
        _continuous_tea_2.ContinuousTeaVerdict.CERTIFIED,
        "partial",
        "segment-event-exact-v1",
        (reversed_event,),
    )

    assert trace.exact_verdict == "certified"
    assert trace.partition.canonical_digest == certificate.canonical_digest
    assert trace.events[0].feature_ids == (b"feature-a", b"feature-b")
    assert trace.events[0].branch_ids == (b"branch-a", b"branch-b")
    assert trace.events[0].canonical_id == reordered.events[0].canonical_id
    assert trace.canonical_bytes == reordered.canonical_bytes
    assert trace.canonical_digest == hashlib.sha256(trace.canonical_bytes).digest()
    assert trace.event_cell_count == len(certificate.cells)


def test_shared_event_trace_rejects_unverified_partition() -> None:
    cap_chord_ratio = 64.0 / 65.0
    certificate = _shared_trace_certificate(cap_chord_ratio)
    mutation = _continuous_tea_2.mutate_certificate_record(
        certificate,
        "alter-multiplicity",
    )

    with pytest.raises(_continuous_tea_2.EventTraceVerificationError):
        _continuous_tea_2.build_event_trace(
            mutation,
            "segment-linear-0-1-v1",
            certificate.source_payload,
            certificate.canonical_digest,
            _continuous_tea_2.ContinuousTeaVerdict.CERTIFIED,
            "partial",
            "segment-event-exact-v1",
            (),
        )
