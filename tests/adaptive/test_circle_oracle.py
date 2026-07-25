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
    return _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        center[0],
        center[1],
        phase[0],
        phase[1],
        clockwise,
        tool_radius,
        cap_chord_ratio,
    )


def _partition_events(trace: object) -> tuple[object, ...]:
    return tuple(
        event
        for fibre in trace.partition.fibres
        for event in fibre.events
    )


def _event_ids(trace: object) -> tuple[bytes, ...]:
    return tuple(event.canonical_id for event in trace.events)


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
            assert not _stock_2.engagement_at(
                stock,
                float(center_x + dx),
                float(center_y + dy),
                tool_radius,
                cap_chord_ratio,
                0.0,
            )[2]


def test_zero_phase_is_rejected_by_the_typed_factory() -> None:
    with pytest.raises(DegenerateCircleMotionError):
        ExactCircleMotion.build(
            Point2[WorldXY].build(5.0, 5.0),
            Vector2[WorldXY].build(0.0, 0.0),
            False,
        )


def test_full_circle_trace_covers_four_center_charts_and_owns_each_seam_once() -> None:
    _, trace = _audit(_stock_2.Stock2(SQUARE, []))

    center_charts = tuple(
        chart
        for chart in trace.partition.charts
        if chart.family == "center-circle"
    )
    assert tuple(chart.chart_id for chart in center_charts) == CENTER_CHART_IDS
    assert all(
        (chart.domain_low, chart.domain_high, chart.orientation)
        == ("0", "1", "ccw")
        for chart in center_charts
    )
    center_seams = tuple(
        seam
        for seam in trace.partition.seams
        if seam.owner_chart_id in CENTER_CHART_IDS
    )
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
    assert _event_ids(counterclockwise)
    assert _event_ids(clockwise) == tuple(reversed(_event_ids(counterclockwise)))
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

    coincident_events = tuple(
        event
        for event in _partition_events(coincident)
        if (event.kind, event.disposition)
        == ("overlap", "isolated-coincidence")
    )
    near_events = tuple(
        event
        for event in _partition_events(near)
        if event.kind == "tangency"
        and event.disposition == "external"
    )
    assert coincident_verdict == near_verdict == "cap_exceeded"
    assert len(coincident_events) == 1
    assert not any(
        event.kind == "overlap"
        for event in _partition_events(near)
    )
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

    tangencies = tuple(
        event
        for event in _partition_events(trace)
        if event.kind == "tangency"
        and event.disposition == expected_disposition
    )
    assert len({event.global_fibre_id for event in tangencies}) == 1
    assert all(event.original_equations_rechecked for event in tangencies)
    assert all(event.trim_disposition == "accepted" for event in tangencies)


def test_two_disk_vertex_events_drive_one_merge_and_one_split() -> None:
    _, trace = _audit(
        _stock_with_disks(
            (4.5, 5.0, 1.0),
            (5.5, 5.0, 1.0),
        ),
        tool_radius=1.0,
    )

    dispositions = tuple(
        event.disposition
        for event in _partition_events(trace)
        if event.kind == "endpoint-order"
    )
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
    assert any(
        event.kind == "endpoint-order"
        for event in _partition_events(trace)
    )
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

    simultaneous = tuple(
        fibre
        for fibre in trace.partition.fibres
        if len(fibre.events) >= 3
        and any(event.kind == "vertex" for event in fibre.events)
    )
    assert len(simultaneous) == 1
    assert len(
        {
            event.canonical_id
            for event in simultaneous[0].events
        }
    ) >= 3


def test_exact_pi_cap_equality_certifies_and_retains_the_cap_root() -> None:
    stock = _stock_2.Stock2(LOWER_HALF, [])
    verdict, trace = _audit(
        stock,
        center=(5.0, 6.0),
        phase=(0.0, -1.0),
        tool_radius=0.5,
        cap_chord_ratio=4.0,
    )

    cap_equalities = tuple(
        event
        for event in _partition_events(trace)
        if (event.kind, event.disposition)
        == ("cap-crossing", "equal")
    )
    assert verdict == "certified"
    assert len(cap_equalities) == 1
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
    mutant = mutate(trace.partition)
    with pytest.raises(_continuous_tea_2.EventPartitionVerificationError):
        _continuous_tea_2.verify_event_partition(mutant)


def test_certificate_rejects_removed_center_chart_and_seam() -> None:
    _, trace = _audit(_stock_2.Stock2(SQUARE, []))

    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            charts=tuple(
                chart
                for chart in certificate.charts
                if chart.chart_id != "center-quarter-2-v1"
            ),
        ),
    )
    removed_seam = next(
        seam
        for seam in trace.partition.seams
        if seam.owner_chart_id == "center-quarter-2-v1"
    )
    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            seams=tuple(
                seam
                for seam in certificate.seams
                if seam.seam_id != removed_seam.seam_id
            ),
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
                event.kind,
                event.disposition,
            )
            == ("overlap", "isolated-coincidence"),
        ),
        (
            lambda: _audit(
                _stock_with_disks(
                    (4.5, 5.0, 1.0),
                    (5.5, 5.0, 1.0),
                ),
                tool_radius=1.0,
            )[1],
            lambda event: (
                event.kind,
                event.disposition,
            )
            == ("endpoint-order", "merge"),
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
                event.kind,
                event.disposition,
            )
            == ("cap-crossing", "equal"),
        ),
    ),
    ids=("coincidence", "merge", "cap-root"),
)
def test_certificate_rejects_removed_required_event_root(
    trace_factory: Callable[[], object],
    event_match: Callable[[object], bool],
) -> None:
    trace = trace_factory()
    fibre = next(
        fibre
        for fibre in trace.partition.fibres
        if any(event_match(event) for event in fibre.events)
    )

    _assert_rejects_mutation(
        trace,
        lambda certificate: _certificate_bytes(
            certificate,
            roots=tuple(
                root
                for root in certificate.roots
                if root.global_fibre_id != fibre.global_fibre_id
            ),
            fibres=tuple(
                candidate
                for candidate in certificate.fibres
                if candidate.global_fibre_id != fibre.global_fibre_id
            ),
        ),
    )
