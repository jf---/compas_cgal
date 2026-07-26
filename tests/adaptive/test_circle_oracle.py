from __future__ import annotations

import numpy as np

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)
CENTER_CHART_IDS = (
    "center-quarter-0-v1",
    "center-quarter-1-v1",
    "center-quarter-2-v1",
    "center-quarter-3-v1",
)


def _audit(stock: _stock_2.Stock2) -> tuple[str, _continuous_tea_2.EventTrace2]:
    return _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        5.0,
        5.0,
        1.0,
        0.0,
        False,
        0.5,
        4.0,
    )


def test_full_circle_trace_owns_all_four_exact_center_seams() -> None:
    verdict, trace = _audit(_stock_2.Stock2(SQUARE, []))

    charts = tuple(chart for chart in trace.partition.charts if chart.family == "center-circle")
    seams = tuple(seam for seam in trace.partition.seams if seam.owner_chart_id in CENTER_CHART_IDS)
    assert verdict == "cap_exceeded"
    assert tuple(chart.chart_id for chart in charts) == CENTER_CHART_IDS
    assert tuple(seam.owner_chart_id for seam in seams) == CENTER_CHART_IDS
    assert len({seam.seam_id for seam in seams}) == 4


def test_empty_exact_stock_certifies_every_full_circle() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 100.0)
    assert stock.is_empty()

    verdict, trace = _audit(stock)

    assert verdict == "certified"
    assert trace.exact_verdict == "certified"
    assert trace.whole_rim_disposition == "clear"


def test_virgin_slotting_proves_cap_exceeded_without_sampling() -> None:
    verdict, trace = _audit(_stock_2.Stock2(SQUARE, []))

    assert verdict == "cap_exceeded"
    assert trace.exact_verdict == "cap_exceeded"
    assert trace.whole_rim_disposition == "material"


def test_uniform_circle_certificate_replays_and_rejects_mutation() -> None:
    _, trace = _audit(_stock_2.Stock2(SQUARE, []))

    verified = _continuous_tea_2.verify_event_partition(trace.partition)
    mutated = _continuous_tea_2.mutate_certificate_record(
        trace.partition,
        "alter-disposition",
    )

    assert verified.verdict.name == "CERTIFIED"
    assert _continuous_tea_2.verify_event_partition(mutated).verdict.name == "UNRESOLVED_DEGENERACY"


def test_nonuniform_axis_circle_replays_every_exact_line_rim_pullback() -> None:
    verdict, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        _stock_2.Stock2(SQUARE, []),
        0.5,
        0.5,
        1.0,
        0.0,
        False,
        1.0,
        4.0,
    )

    assert verdict == "unresolved"
    assert trace.partition.source_kind == "full-circle-line-pullbacks-v1"
    pullbacks = tuple(projection for projection in trace.partition.projections if projection.degree_bound_id == "full-circle-line-(2,2)-v1")
    event_projections = tuple(projection for projection in trace.partition.projections if projection.degree_bound_id.startswith("exact-univariate-"))
    assert len(pullbacks) == 32
    assert len(event_projections) == 32
    assert trace.partition.roots
    assert len(trace.partition.fibres) == len(trace.partition.roots)
    assert {event.kind for fibre in trace.partition.fibres for event in fibre.events} >= {"tangent", "endpoint-order"}
    assert _continuous_tea_2.verify_event_partition(trace.partition).verdict.name == "CERTIFIED"


def test_exact_rational_probe_uses_the_circle_chart_not_binary64_sampling() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 100.0)

    assert not _continuous_tea_2.full_circle_rational_probe_exceeds_cap_exact(
        stock,
        5.0,
        5.0,
        1.0,
        0.0,
        0,
        1,
        3,
        0.5,
        4.0,
    )
