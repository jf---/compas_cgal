from __future__ import annotations

from collections.abc import Callable
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import DegenerateCircleMotionError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)
LOWER_HALF = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 5.0, 0.0], [0.0, 5.0, 0.0]],
    dtype=np.float64,
)
CENTER_CHART_IDS = (
    "center-quarter-0-v1",
    "center-quarter-1-v1",
    "center-quarter-2-v1",
    "center-quarter-3-v1",
)
CIRCLE_FALSIFIER_DEPTH = 6


def _audit(
    stock: _stock_2.Stock2,
    *,
    center: tuple[float, float] = (5.0, 5.0),
    phase: tuple[float, float] = (1.0, 0.0),
    clockwise: bool = False,
    tool_radius: float = 0.5,
    cap_chord_ratio: float = 4.0,
) -> tuple[str, object]:
    motion = ExactCircleMotion.build(
        Point2[WorldXY].build(*center),
        Vector2[WorldXY].build(*phase),
        clockwise,
    )
    return _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        motion.center.x,
        motion.center.y,
        motion.phase_vector.x,
        motion.phase_vector.y,
        motion.clockwise,
        tool_radius,
        cap_chord_ratio,
    )


def _partition_events(trace: object) -> tuple[object, ...]:
    return tuple(event for fibre in trace.partition_certificate.fibres for event in fibre.events)


def _oriented_fibre_blocks(trace: object) -> tuple[tuple[bytes, ...], ...]:
    blocks: list[list[bytes]] = []
    previous_fibre_id: bytes | None = None
    for event in trace.events:
        if event.global_fibre_id != previous_fibre_id:
            blocks.append([])
            previous_fibre_id = event.global_fibre_id
        blocks[-1].append(event.canonical_id)
    return tuple(tuple(block) for block in blocks)


def _stock_with_disks(
    *disks: tuple[float, float, float],
) -> _stock_2.Stock2:
    stock = _stock_2.Stock2(SQUARE, [])
    for disk in disks:
        stock.subtract_disk(*disk)
    return stock


def _slotted_lower_half() -> _stock_2.Stock2:
    stock = _stock_2.Stock2(LOWER_HALF, [])
    stock.subtract_capsule(5.0, 4.30, 5.0, 4.53, 0.05)
    return stock


def _two_disk_merge_split_trace() -> object:
    return _audit(
        _stock_with_disks(
            (4.5, 5.0, 1.0),
            (5.5, 5.0, 1.0),
        ),
        center=(5.0, 6.25),
        phase=(0.0, -0.125),
        tool_radius=0.5,
    )[1]


def _verified_circle_order_fixture() -> tuple[
    _continuous_tea_2.VerifiedEventPartition2,
    tuple[_continuous_tea_2.EventTraceEvent2, ...],
]:
    certificate = _continuous_tea_2.partition_cap_crossings(
        ("0", "1"),
        ("1",),
        ("0", "3"),
        ("1",),
        "64",
        "65",
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
    verified = _continuous_tea_2.verify_event_partition(certificate)
    center_charts = tuple(chart for chart in verified.partition.charts if chart.family == "center-circle")
    center_seams = tuple(seam for seam in verified.partition.seams if seam.owner_chart_id in CENTER_CHART_IDS)
    assert verified.verdict is _continuous_tea_2.ContinuousTeaVerdict.CERTIFIED
    assert tuple(chart.chart_id for chart in center_charts) == CENTER_CHART_IDS
    assert tuple(seam.owner_chart_id for seam in center_seams) == CENTER_CHART_IDS

    first_root = verified.partition.roots[0].root_id
    second_root = verified.partition.roots[1].root_id

    def event(
        label: bytes,
        *,
        root_id: bytes,
        global_fibre_id: bytes,
        kind: str,
        motion_order: int,
        feature_id: bytes | None = None,
    ) -> _continuous_tea_2.EventTraceEvent2:
        return _continuous_tea_2.EventTraceEvent2(
            root_id,
            global_fibre_id,
            kind,
            (feature_id or b"feature-" + label,),
            (b"branch-" + label,),
            1,
            "owned" if kind == "seam" else "transverse",
            motion_order,
        )

    return (
        verified,
        (
            event(
                b"quarter-2",
                root_id=second_root,
                global_fibre_id=b"circle-fibre-2",
                kind="endpoint-order",
                motion_order=2,
            ),
            event(
                b"quarter-1-b",
                root_id=first_root,
                global_fibre_id=b"circle-fibre-1",
                kind="endpoint-order",
                motion_order=1,
            ),
            event(
                b"anchor",
                root_id=first_root,
                global_fibre_id=b"circle-fibre-anchor",
                kind="seam",
                motion_order=0,
                feature_id=center_charts[0].start_seam_id,
            ),
            event(
                b"quarter-3",
                root_id=first_root,
                global_fibre_id=b"circle-fibre-3",
                kind="endpoint-order",
                motion_order=3,
            ),
            event(
                b"quarter-1-a",
                root_id=first_root,
                global_fibre_id=b"circle-fibre-1",
                kind="endpoint-order",
                motion_order=1,
            ),
        ),
    )


def _event_blocks(
    events: tuple[_continuous_tea_2.EventTraceEvent2, ...],
) -> tuple[tuple[bytes, ...], ...]:
    blocks: list[list[bytes]] = []
    previous_fibre_id: bytes | None = None
    for event in events:
        if event.global_fibre_id != previous_fibre_id:
            blocks.append([])
            previous_fibre_id = event.global_fibre_id
        blocks[-1].append(event.canonical_id)
    return tuple(tuple(block) for block in blocks)


def _quarter_phase(
    phase: tuple[float, float],
    chart: int,
    parameter: Fraction,
) -> tuple[Fraction, Fraction]:
    denominator = 1 + parameter * parameter
    cosine = (1 - parameter * parameter) / denominator
    sine = 2 * parameter / denominator
    px = Fraction.from_float(phase[0])
    py = Fraction.from_float(phase[1])
    quarter_turns = (
        (px, py),
        (-py, px),
        (-px, -py),
        (py, -px),
    )
    tangent_turns = (
        (-py, px),
        (-px, -py),
        (py, -px),
        (px, py),
    )
    start_x, start_y = quarter_turns[chart]
    tangent_x, tangent_y = tangent_turns[chart]
    return (
        cosine * start_x + sine * tangent_x,
        cosine * start_y + sine * tangent_y,
    )


def _assert_no_violating_dyadic_probe(
    stock: _stock_2.Stock2,
    *,
    center: tuple[float, float],
    phase: tuple[float, float],
    tool_radius: float,
    cap_chord_ratio: float,
) -> None:
    denominator = 2**CIRCLE_FALSIFIER_DEPTH
    center_x = Fraction.from_float(center[0])
    center_y = Fraction.from_float(center[1])
    for chart in range(4):
        for numerator in range(denominator):
            dx, dy = _quarter_phase(
                phase,
                chart,
                Fraction(numerator, denominator),
            )
            rounded_exceeded = _stock_2.engagement_at(
                stock,
                float(center_x + dx),
                float(center_y + dy),
                tool_radius,
                cap_chord_ratio,
                0.0,
            )[2]
            if not rounded_exceeded:
                continue
            # The binary64 station is only near the rational chart point. A
            # hit is a falsifier trigger, never proof, until the original
            # rational motion point is checked without the rounding detour.
            exact_exceeded = _continuous_tea_2.full_circle_rational_probe_exceeds_cap_exact(
                stock,
                center[0],
                center[1],
                phase[0],
                phase[1],
                chart,
                numerator,
                denominator,
                tool_radius,
                cap_chord_ratio,
            )
            assert not exact_exceeded


def test_zero_phase_is_rejected_by_the_typed_factory() -> None:
    with pytest.raises(DegenerateCircleMotionError):
        ExactCircleMotion.build(
            Point2[WorldXY].build(5.0, 5.0),
            Vector2[WorldXY].build(0.0, 0.0),
            False,
        )


def test_zero_phase_is_rejected_before_native_dispatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dispatched = False

    def record_dispatch(*_args: object) -> None:
        nonlocal dispatched
        dispatched = True

    monkeypatch.setattr(
        _continuous_tea_2,
        "audit_full_circle_tea_event_exact",
        record_dispatch,
        raising=False,
    )

    with pytest.raises(DegenerateCircleMotionError):
        _audit(
            _stock_2.Stock2(SQUARE, []),
            phase=(0.0, 0.0),
        )

    assert not dispatched


def test_full_circle_order_is_ccw_canonical_and_insertion_independent() -> None:
    verified, events = _verified_circle_order_fixture()

    ordered = tuple(
        _continuous_tea_2.order_full_circle_events(
            verified,
            False,
            events,
        )
    )
    permuted = tuple(
        _continuous_tea_2.order_full_circle_events(
            verified,
            False,
            tuple(reversed(events)),
        )
    )
    expected = tuple(
        event.canonical_id
        for event in sorted(
            events,
            key=lambda event: (
                event.motion_order,
                event.canonical_id,
            ),
        )
    )

    assert tuple(event.canonical_id for event in ordered) == expected
    assert tuple(event.canonical_id for event in permuted) == expected


def test_clockwise_full_circle_order_keeps_anchor_and_reverses_remaining_fibres() -> None:
    verified, events = _verified_circle_order_fixture()

    counterclockwise = tuple(
        _continuous_tea_2.order_full_circle_events(
            verified,
            False,
            events,
        )
    )
    clockwise = tuple(
        _continuous_tea_2.order_full_circle_events(
            verified,
            True,
            tuple(reversed(events)),
        )
    )
    counterclockwise_blocks = _event_blocks(counterclockwise)
    clockwise_blocks = _event_blocks(clockwise)

    assert counterclockwise[0].kind == clockwise[0].kind == "seam"
    assert counterclockwise_blocks[0] == clockwise_blocks[0]
    assert clockwise_blocks[1:] == tuple(reversed(counterclockwise_blocks[1:]))
    assert all(block == tuple(sorted(block)) for block in clockwise_blocks)


def test_full_circle_trace_covers_four_center_charts_and_owns_each_seam_once() -> None:
    _, trace = _audit(_stock_2.Stock2(SQUARE, []))

    center_charts = tuple(chart for chart in trace.partition_certificate.charts if chart.family == "center-circle")
    assert tuple(chart.chart_id for chart in center_charts) == CENTER_CHART_IDS
    assert all((chart.domain_low, chart.domain_high, chart.orientation) == ("0", "1", "ccw") for chart in center_charts)
    center_seams = tuple(seam for seam in trace.partition_certificate.seams if seam.owner_chart_id in CENTER_CHART_IDS)
    assert len(center_seams) == 4
    assert tuple(seam.owner_chart_id for seam in center_seams) == CENTER_CHART_IDS
    assert len({seam.seam_id for seam in center_seams}) == 4


def test_clockwise_and_counterclockwise_have_same_verdict_and_reverse_event_order() -> None:
    counterclockwise_verdict, counterclockwise = _audit(
        _stock_2.Stock2(LOWER_HALF, []),
        clockwise=False,
    )
    clockwise_verdict, clockwise = _audit(
        _stock_2.Stock2(LOWER_HALF, []),
        clockwise=True,
    )

    assert counterclockwise_verdict == clockwise_verdict == "cap_exceeded"
    counterclockwise_blocks = _oriented_fibre_blocks(counterclockwise)
    clockwise_blocks = _oriented_fibre_blocks(clockwise)
    assert counterclockwise.events[0].kind == clockwise.events[0].kind == "seam"
    assert counterclockwise_blocks[0] == clockwise_blocks[0]
    assert clockwise_blocks[1:] == tuple(reversed(counterclockwise_blocks[1:]))
    assert counterclockwise.canonical_digest != clockwise.canonical_digest


def test_coincident_disk_is_one_overlap_fibre_and_near_coincidence_stays_distinct() -> None:
    coincident_verdict, coincident = _audit(
        _stock_with_disks((6.0, 5.0, 1.0)),
        tool_radius=1.0,
    )
    near_verdict, near = _audit(
        _stock_with_disks((6.0 + 2.0**-20, 5.0, 1.0)),
        tool_radius=1.0,
    )

    coincident_events = tuple(event for event in _partition_events(coincident) if (event.kind, event.disposition) == ("overlap", "isolated-coincidence"))
    near_events = tuple(event for event in _partition_events(near) if event.kind == "tangency" and event.disposition == "external")
    assert coincident_verdict == near_verdict == "cap_exceeded"
    assert len(coincident_events) == 1
    assert not any(event.kind == "overlap" for event in _partition_events(near))
    assert len({event.global_fibre_id for event in near_events}) == 2


@pytest.mark.parametrize(
    ("disk", "tool_radius", "expected_disposition"),
    (
        ((7.0, 5.0, 0.5), 0.5, "external"),
        ((6.5, 5.0, 0.5), 1.0, "internal"),
    ),
)
def test_external_and_internal_tangency_are_exact_fibres(
    disk: tuple[float, float, float],
    tool_radius: float,
    expected_disposition: str,
) -> None:
    _, trace = _audit(
        _stock_with_disks(disk),
        tool_radius=tool_radius,
    )

    tangencies = tuple(event for event in _partition_events(trace) if event.kind == "tangency" and event.disposition == expected_disposition)
    assert len({event.global_fibre_id for event in tangencies}) == 1
    assert all(event.original_equations_rechecked for event in tangencies)
    assert all(event.trim_disposition == "accepted" for event in tangencies)


def test_two_disk_vertex_events_drive_one_merge_and_one_split() -> None:
    # The unit holes meet at V+=(5, 5+sqrt(3)/2) and V- below. Cutter-center
    # events occur when the radius-1/2 cutter rim contains V. For the guide
    # center O=(5,25/4) and guide radius 1/8,
    #
    #   3/8 < |O-V+| = 5/4-sqrt(3)/2 < 5/8,
    #
    # because 3 < 49/16 and 25/16 < 3. The guide therefore crosses the V+
    # event circle exactly twice. V- is farther than 5/8, so it contributes
    # none. At the start C=(5,49/8), Q=(5,45/8) lies on the cutter rim and
    # inside both holes since |Q-(4.5,5)|^2=|Q-(5.5,5)|^2=41/64<1: the void
    # arcs overlap. Opposite, the rim has y>=47/8>5+sqrt(3)/2, while each
    # hole still intersects it (center-distance squared 137/64<9/4), so the
    # void arcs are disjoint. The two transverse V+ events are thus one
    # material-run split and one merge.
    trace = _two_disk_merge_split_trace()

    dispositions = tuple(event.disposition for event in _partition_events(trace) if event.kind == "endpoint-order")
    assert set(dispositions) == {"merge", "split"}
    assert dispositions.count("merge") == 1
    assert dispositions.count("split") == 1


def test_slotted_full_circle_is_certified_by_events_not_gap_closure_sampling() -> None:
    stock = _slotted_lower_half()
    verdict, trace = _audit(
        stock,
        center=(5.0, 5.05),
        phase=(0.0, -0.05),
        tool_radius=0.5,
        cap_chord_ratio=4.0,
    )

    assert verdict == "certified"
    assert any(event.kind == "endpoint-order" for event in _partition_events(trace))
    _assert_no_violating_dyadic_probe(
        stock,
        center=(5.0, 5.05),
        phase=(0.0, -0.05),
        tool_radius=0.5,
        cap_chord_ratio=4.0,
    )


def test_contour_vertex_retains_all_simultaneous_event_incidences() -> None:
    _, trace = _audit(
        _stock_with_disks(
            (4.0, 5.0, 1.0),
            (6.0, 5.0, 1.0),
        ),
        center=(5.0, 3.0),
        phase=(0.0, 1.0),
        tool_radius=1.0,
    )

    simultaneous = tuple(fibre for fibre in trace.partition_certificate.fibres if len(fibre.events) >= 3 and any(event.kind == "vertex" for event in fibre.events))
    assert len(simultaneous) == 1
    assert len({event.canonical_id for event in simultaneous[0].events}) >= 3


def test_exact_pi_cap_equality_certifies_and_retains_the_cap_root() -> None:
    stock = _stock_2.Stock2(LOWER_HALF, [])
    verdict, trace = _audit(
        stock,
        center=(5.0, 6.0),
        phase=(0.0, -1.0),
        tool_radius=0.5,
        cap_chord_ratio=4.0,
    )

    cap_contacts = tuple(event for event in _partition_events(trace) if (event.kind, event.disposition) == ("cap-contact", "touch"))
    assert verdict == "certified"
    assert len(cap_contacts) == 1
    assert cap_contacts[0].multiplicity % 2 == 0
    assert cap_contacts[0].adjacent_cap_dispositions == ("below", "below")
    assert cap_contacts[0].adjacent_run_counts == (1, 1)
    _assert_no_violating_dyadic_probe(
        stock,
        center=(5.0, 6.0),
        phase=(0.0, -1.0),
        tool_radius=0.5,
        cap_chord_ratio=4.0,
    )


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


def _assert_rejects_mutation(
    trace: object,
    mutate: Callable[[object], bytes],
) -> None:
    mutant = mutate(trace.partition_certificate)
    with pytest.raises(_continuous_tea_2.EventPartitionVerificationError):
        _continuous_tea_2.verify_event_partition(mutant)


def test_certificate_rejects_removed_center_chart_and_seam() -> None:
    _, trace = _audit(_stock_2.Stock2(SQUARE, []))

    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            charts=tuple(chart for chart in certificate.charts if chart.chart_id != "center-quarter-2-v1"),
        ),
    )
    removed_seam = next(seam for seam in trace.partition_certificate.seams if seam.owner_chart_id == "center-quarter-2-v1")
    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            seams=tuple(seam for seam in certificate.seams if seam.seam_id != removed_seam.seam_id),
        ),
    )


@pytest.mark.parametrize(
    ("trace_factory", "event_match"),
    (
        (
            lambda: _audit(
                _stock_with_disks((6.0, 5.0, 1.0)),
                tool_radius=1.0,
            )[1],
            lambda event: (
                (
                    event.kind,
                    event.disposition,
                )
                == ("overlap", "isolated-coincidence")
            ),
        ),
        (
            _two_disk_merge_split_trace,
            lambda event: (
                (
                    event.kind,
                    event.disposition,
                )
                == ("endpoint-order", "merge")
            ),
        ),
        (
            lambda: _audit(
                _stock_2.Stock2(LOWER_HALF, []),
                center=(5.0, 6.0),
                phase=(0.0, -1.0),
                tool_radius=0.5,
                cap_chord_ratio=4.0,
            )[1],
            lambda event: (
                (
                    event.kind,
                    event.disposition,
                )
                == ("cap-contact", "touch")
            ),
        ),
    ),
    ids=("coincidence", "merge", "cap-root"),
)
def test_certificate_rejects_removed_required_root_with_fibre_retained(
    trace_factory: Callable[[], object],
    event_match: Callable[[object], bool],
) -> None:
    trace = trace_factory()
    fibre = next(fibre for fibre in trace.partition_certificate.fibres if any(event_match(event) for event in fibre.events))

    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            roots=tuple(root for root in certificate.roots if root.global_fibre_id != fibre.global_fibre_id),
        ),
    )


@pytest.mark.parametrize(
    ("trace_factory", "event_match"),
    (
        (
            lambda: _audit(
                _stock_with_disks((6.0, 5.0, 1.0)),
                tool_radius=1.0,
            )[1],
            lambda event: (
                (
                    event.kind,
                    event.disposition,
                )
                == ("overlap", "isolated-coincidence")
            ),
        ),
        (
            _two_disk_merge_split_trace,
            lambda event: (
                (
                    event.kind,
                    event.disposition,
                )
                == ("endpoint-order", "merge")
            ),
        ),
        (
            lambda: _audit(
                _stock_2.Stock2(LOWER_HALF, []),
                center=(5.0, 6.0),
                phase=(0.0, -1.0),
                tool_radius=0.5,
                cap_chord_ratio=4.0,
            )[1],
            lambda event: (
                (
                    event.kind,
                    event.disposition,
                )
                == ("cap-contact", "touch")
            ),
        ),
    ),
    ids=("coincidence", "merge", "cap-root"),
)
def test_certificate_rejects_removed_required_fibre_with_root_retained(
    trace_factory: Callable[[], object],
    event_match: Callable[[object], bool],
) -> None:
    trace = trace_factory()
    fibre = next(fibre for fibre in trace.partition_certificate.fibres if any(event_match(event) for event in fibre.events))

    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            fibres=tuple(candidate for candidate in certificate.fibres if candidate.global_fibre_id != fibre.global_fibre_id),
        ),
    )
