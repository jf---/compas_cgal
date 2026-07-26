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

    assert verdict == "cap_exceeded"
    assert trace.exact_verdict == "cap_exceeded"
    assert trace.whole_rim_disposition == "partial"
    assert trace.partition.source_kind == "full-circle-boundary-pullbacks-v2"
    pullbacks = tuple(projection for projection in trace.partition.projections if projection.degree_bound_id == "full-circle-line-(2,2)-v1")
    event_projections = tuple(projection for projection in trace.partition.projections if projection.degree_bound_id.startswith("exact-univariate-"))
    assert len(pullbacks) == 32
    assert len(event_projections) == 48
    assert trace.partition.roots
    assert len(trace.partition.fibres) == len(trace.partition.roots)
    assert {event.kind for fibre in trace.partition.fibres for event in fibre.events} >= {"tangent", "endpoint-order"}
    assert _continuous_tea_2.verify_event_partition(trace.partition).verdict.name == "CERTIFIED"


def test_slotted_stock_replays_circle_tangency_and_coincidence_fibres() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(0.0, 5.0, 2.0)
    line_support_ids = {record.support_id for record in _continuous_tea_2.extract_boundary_records(stock) if record.support_kind == "line"}

    verdict, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        1.0,
        5.0,
        1.0,
        0.0,
        False,
        2.0,
        4.0,
    )

    circle_pullbacks = tuple(projection for projection in trace.partition.projections if projection.degree_bound_id == "full-circle-circle-(4,4)-v1")
    event_kinds = {event.kind for fibre in trace.partition.fibres for event in fibre.events}
    directed_fibres = tuple(fibre for fibre in trace.partition.fibres if fibre.ccw_direction in {"merge", "split"})
    assert verdict == "unresolved"
    assert circle_pullbacks
    assert event_kinds >= {"tangent", "support-overlap", "cap-crossing"}
    assert any(event.kind == "cap-crossing" and event.support_id in line_support_ids for fibre in trace.partition.fibres for event in fibre.events)
    assert directed_fibres
    assert all({fibre.ccw_direction, fibre.cw_direction} == {"merge", "split"} for fibre in directed_fibres)
    assert all(branch.branch_id and branch.support_id for fibre in directed_fibres for branch in (*fibre.left_active_branches, *fibre.right_active_branches))
    assert all(
        {branch.support_id for branch in (*fibre.left_active_branches, *fibre.right_active_branches)}
        == {event.support_id for event in fibre.events if event.kind == "endpoint-order"}
        for fibre in directed_fibres
    )
    assert any(
        len({event.support_id for event in fibre.events if event.kind == "endpoint-order"}) >= 2
        and {event.disposition for event in fibre.events if event.kind == "endpoint-order"} == {"merge-split"}
        for fibre in trace.partition.fibres
    )
    assert _continuous_tea_2.verify_event_partition(trace.partition).verdict.name == "CERTIFIED"


def test_negative_vertex_radial_predicate_preserves_merge_split_sign() -> None:
    stock = _stock_2.Stock2(
        np.array(
            [
                [-1.5, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [5.0, 5.0, 0.0],
                [-1.5, 5.0, 0.0],
            ],
            dtype=np.float64,
        ),
        [],
    )
    records = _continuous_tea_2.extract_boundary_records(stock)
    bottom = next(record for record in records if tuple(record.primitive_coefficients) == ("0", "1", "0"))
    left = next(record for record in records if tuple(record.primitive_coefficients) == ("2", "0", "3"))
    target_supports = {bottom.support_id, left.support_id}

    _, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        0.0,
        0.0,
        1.0,
        0.0,
        False,
        2.0,
        4.0,
    )

    target_fibres = tuple(fibre for fibre in trace.partition.fibres if {event.support_id for event in fibre.events if event.kind == "endpoint-order"} == target_supports)
    directions = tuple(fibre.ccw_direction for fibre in target_fibres)
    assert set(directions) == {"split", "merge"}
    assert directions == tuple(sorted(directions, key={"split": 0, "merge": 1}.__getitem__))
    assert any(projection.signed_predicate_coefficients for projection in trace.partition.projections)
    assert _continuous_tea_2.verify_event_partition(trace.partition).verdict.name == "CERTIFIED"


def test_general_phase_replays_tangencies_and_oriented_trace() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(0.0, 5.0, 2.0)

    ccw_verdict, ccw = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        2.0,
        5.0,
        1.2,
        1.6,
        False,
        1.0,
        4.0,
    )
    cw_verdict, cw = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        2.0,
        5.0,
        1.2,
        1.6,
        True,
        1.0,
        4.0,
    )

    assert ccw_verdict == cw_verdict == "unresolved"
    assert ccw.events and cw.events
    assert ccw.canonical_digest != cw.canonical_digest
    assert tuple(event.canonical_id for event in ccw.events) != tuple(event.canonical_id for event in cw.events)
    assert {"external-contact", "internal-contact"} <= {event.disposition for event in ccw.events if event.kind == "tangent"}
    mutated = _continuous_tea_2.mutate_certificate_record(
        ccw.partition,
        "alter-fibre-direction",
    )
    assert _continuous_tea_2.verify_event_partition(mutated).verdict.name == "UNRESOLVED_DEGENERACY"
    assert _continuous_tea_2.verify_event_partition(ccw.partition).verdict.name == "CERTIFIED"


def test_same_support_endpoint_sheets_remain_distinct() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 2.0)

    _, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        5.0,
        5.0,
        2.0,
        1.0,
        False,
        1.0,
        4.0,
    )

    same_support = tuple(
        fibre
        for fibre in trace.partition.fibres
        if fibre.ccw_direction in {"merge", "split"} and len({event.support_id for event in fibre.events if event.kind == "endpoint-order"}) == 1
    )
    assert same_support
    assert all(
        len({branch.branch_id for branch in active}) == len(active) and len({branch.feature_id for branch in active}) == len(active)
        for fibre in same_support
        for active in (fibre.left_active_branches, fibre.right_active_branches)
        if active
    )


def test_phase_seam_is_one_owned_fibre_with_oriented_cyclic_order() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 2.0)
    audits = tuple(
        _continuous_tea_2.audit_full_circle_tea_event_exact(
            stock,
            5.0,
            5.0,
            2.0,
            1.0,
            clockwise,
            1.0,
            4.0,
        )[1]
        for clockwise in (False, True)
    )
    ccw, cw = audits

    seam_fibres = tuple(fibre for fibre in ccw.partition.fibres if len(fibre.local_root_ids) == 2)
    assert len(seam_fibres) == 1
    assert seam_fibres[0].events
    assert seam_fibres[0].seam_id == next(seam.seam_id for seam in ccw.partition.seams if seam.owner_chart_id == "center-quarter-0-v1")
    assert (seam_fibres[0].ccw_direction, seam_fibres[0].cw_direction) == (
        "merge",
        "split",
    )
    assert {witness.local_root_id for witness in seam_fibres[0].local_event_witnesses} == set(seam_fibres[0].local_root_ids)
    assert len(seam_fibres[0].local_event_witnesses) > len(
        {
            (
                witness.kind,
                witness.feature_id,
                witness.support_id,
                witness.trim_id,
                witness.vertex_id,
            )
            for witness in seam_fibres[0].local_event_witnesses
        }
    )

    def fibre_order(trace: _continuous_tea_2.EventTrace2) -> tuple[bytes, ...]:
        return tuple(dict.fromkeys(event.global_fibre_id for event in trace.events))

    ccw_order = fibre_order(ccw)
    cw_order = fibre_order(cw)
    assert cw_order == (ccw_order[0], *reversed(ccw_order[1:]))
    assert len({event.canonical_id for event in ccw.events}) == len(ccw.events)
    assert sum(global_id == ccw_order[0] for global_id in ccw_order) == 1
    seam_trace_events = tuple(event for event in ccw.events if event.global_fibre_id == ccw_order[0])
    assert len(seam_trace_events) == len(
        {
            (
                witness.kind,
                witness.feature_id,
                witness.support_id,
                witness.trim_id,
                witness.vertex_id,
            )
            for witness in seam_fibres[0].local_event_witnesses
        }
    )

    for mutation in (
        "delete-signed-predicate",
        "delete-active-sheet",
        "delete-seam",
    ):
        mutated = _continuous_tea_2.mutate_certificate_record(
            ccw.partition,
            mutation,
        )
        assert _continuous_tea_2.verify_event_partition(mutated).verdict.name == "UNRESOLVED_DEGENERACY"


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
