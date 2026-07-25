from collections.abc import Iterable

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

DYADIC_FALSIFIER_PARAMETERS = tuple(index / 32.0 for index in range(33))


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
    return next(root for root in trace.partition_certificate.roots if tuple(root.factor_coefficients) == coefficients and root.root_ordinal == ordinal)


def _assert_trace_binds_verified_partition(
    verdict: str,
    trace,
) -> None:
    assert trace.exact_verdict == verdict
    assert trace.partition_certificate_digest == trace.partition_certificate.canonical_digest
    assert trace.canonical_digest
    assert trace.motion_chart_id
    assert trace.effective_cap_bytes
    assert trace.oracle_strategy_version == "segment-event-exact-v1"
    assert trace.event_cell_count == len(trace.partition_certificate.cells)


def _assert_dyadic_falsifier_clean(
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion,
    tool_radius: float,
    cap_chord_ratio: float,
    parameters: Iterable[float] = DYADIC_FALSIFIER_PARAMETERS,
) -> None:
    for parameter in parameters:
        cx = motion.start.x + parameter * (motion.end.x - motion.start.x)
        cy = motion.start.y + parameter * (motion.end.y - motion.start.y)
        assert not _stock_2.engagement_at(
            stock,
            cx,
            cy,
            tool_radius,
            cap_chord_ratio,
            0.0,
        )[2]


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
    _assert_dyadic_falsifier_clean(stock, motion, 1.0, 4.0)


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
    _assert_dyadic_falsifier_clean(stock, clear_motion, 0.5, 4.0)


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
    assert tangent[0].root_id == _root(trace, ("-1", "2")).root_id
    assert tangent[0].multiplicity == 2
    assert tangent[0].disposition == disposition
    _assert_dyadic_falsifier_clean(stock, motion, 1.0, 4.0)


def test_circle_circle_tangency_root_is_not_omitted() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    stock.subtract_disk(5.0, 1.5, 1.0)
    motion = _motion(3.0, -0.5, 7.0, -0.5)

    verdict, trace = _audit(stock, motion, 1.0, 4.0)

    assert verdict == "certified"
    tangent = _events(trace, "tangent")
    assert len(tangent) == 1
    assert tangent[0].root_id == _root(trace, ("-1", "2")).root_id
    assert tangent[0].multiplicity == 2
    _assert_dyadic_falsifier_clean(stock, motion, 1.0, 4.0)


def test_cutter_rim_vertex_fibre_is_not_omitted() -> None:
    stock = _stock_2.Stock2(HALF_PLANE_WINDOW, [])
    motion = _motion(-2.0, 0.0, 0.0, 0.0)

    verdict, trace = _audit(stock, motion, 1.0, 4.0)

    assert verdict == "certified"
    vertex = _events(trace, "vertex")
    assert len(vertex) == 1
    assert vertex[0].root_id == _root(trace, ("-1", "2")).root_id
    assert vertex[0].feature_ids
    _assert_dyadic_falsifier_clean(stock, motion, 1.0, 4.0)


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
    topology = _events(trace, "merge")
    assert len(topology) == 1
    assert topology[0].root_id
    assert topology[0].disposition == disposition
    _assert_dyadic_falsifier_clean(stock, motion, 0.5, 4.0)


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
    reverse = _stock_2.Stock2(HALF_PLANE_WINDOW[::-1].copy(), [])
    motion = _motion(2.0, 0.0, 8.0, 0.0)

    forward_verdict, forward_trace = _audit(forward, motion, 1.0, 4.0)
    reverse_verdict, reverse_trace = _audit(reverse, motion, 1.0, 4.0)

    assert forward_verdict == reverse_verdict == "certified"
    assert forward_trace.canonical_bytes == reverse_trace.canonical_bytes
    assert forward_trace.canonical_digest == reverse_trace.canonical_digest
    assert [(event.root_id, event.kind, event.feature_ids, event.branch_ids) for event in forward_trace.events] == [
        (event.root_id, event.kind, event.feature_ids, event.branch_ids) for event in reverse_trace.events
    ]
