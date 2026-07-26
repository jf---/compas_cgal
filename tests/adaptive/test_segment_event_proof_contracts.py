from __future__ import annotations

import numpy as np

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2


SQUARE = np.array(
    [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 10.0, 0.0],
        [0.0, 10.0, 0.0],
    ],
    dtype=np.float64,
)


def test_cross_feature_cap_root_binds_the_ordered_pair_event() -> None:
    partition = _continuous_tea_2.construct_segment_event_partition(
        _stock_2.Stock2(SQUARE, []),
        0.6,
        0.3,
        0.6,
        0.6,
        1.0,
        3.75,
    )
    branches = {branch.branch_id: branch for branch in partition.branches}
    transitions = []
    for stratum in partition.strata:
        if stratum.kind != "fibre":
            continue
        left = {pair.pair_sheet_id: pair for pair in stratum.left_pair_dispositions}
        at = {pair.pair_sheet_id: pair for pair in stratum.pair_dispositions}
        right = {pair.pair_sheet_id: pair for pair in stratum.right_pair_dispositions}
        for pair_sheet_id in left.keys() & at.keys() & right.keys():
            pair = at[pair_sheet_id]
            if branches[pair.first_branch_id].feature_id != branches[pair.second_branch_id].feature_id and (
                left[pair_sheet_id].cap_disposition,
                pair.cap_disposition,
                right[pair_sheet_id].cap_disposition,
            ) == ("below-cap", "equal-cap", "above-cap"):
                transitions.append((stratum, pair))

    assert len(transitions) == 1
    stratum, pair = transitions[0]
    pair_events = [event for event in stratum.events if event.kind == "cap-crossing" and event.pair_sheet_id == pair.pair_sheet_id]
    assert len(pair_events) == 1
    assert (
        pair_events[0].first_branch_id,
        pair_events[0].second_branch_id,
    ) == (
        pair.first_branch_id,
        pair.second_branch_id,
    )


def test_identically_equal_cross_feature_cap_is_an_explicit_interval() -> None:
    partition = _continuous_tea_2.construct_segment_event_partition(
        _stock_2.Stock2(SQUARE, []),
        0.8,
        0.8,
        0.9,
        0.9,
        1.0,
        2.0,
    )
    intervals = [overlap for overlap in partition.certificate.overlaps if overlap.kind == "identically-equal-cap-interval"]

    assert intervals
    assert {
        (
            overlap.domain_low,
            overlap.domain_high,
            overlap.orientation_disposition,
        )
        for overlap in intervals
    } == {("0", "1", "ccw")}


def test_segment_vertex_projection_filters_circle_conjugates_by_incidence() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(4.5, 5.0, 1.0)
    stock.subtract_disk(5.5, 5.0, 1.0)
    partition = _continuous_tea_2.construct_segment_event_partition(
        stock,
        5.0,
        5.0,
        5.0,
        7.0,
        1.0,
        4.0,
    )
    circle_features = {branch.feature_id for branch in partition.branches if branch.support_kind == "circle"}
    vertex_fibres = []
    for stratum in partition.strata:
        endpoint_events = [event for event in stratum.events if event.kind == "endpoint-order" and event.feature_id in circle_features]
        support_ids = {event.support_id for event in endpoint_events}
        if len(support_ids) >= 2:
            vertex_fibres.append((stratum, endpoint_events))

    assert len(vertex_fibres) == 2
    assert len({event.vertex_id for _, events in vertex_fibres for event in events}) == 2
    assert all(len({event.vertex_id for event in events}) == 1 and len({event.support_id for event in events}) >= 2 for _, events in vertex_fibres)
    assert vertex_fibres[0][0].root_id != vertex_fibres[1][0].root_id
    branches = {branch.branch_id: branch for branch in partition.branches}
    assert all(
        event.branch_id in stratum.active_branch_ids and branches[event.branch_id].feature_id == event.feature_id and branches[event.branch_id].support_id == event.support_id
        for stratum, events in vertex_fibres
        for event in events
    )


def test_every_segment_event_projection_has_a_complete_degree_contract() -> None:
    partition = _continuous_tea_2.construct_segment_event_partition(
        _stock_2.Stock2(SQUARE, []),
        2.0,
        0.0,
        8.0,
        0.0,
        1.0,
        4.0,
    )

    pullback_count = len(partition.projections)
    event_projections = partition.certificate.projections[pullback_count:]
    assert event_projections
    for projection in event_projections:
        coefficients = projection.coefficient_rows[0]
        actual_motion_degree = max(index for index, coefficient in enumerate(coefficients) if int(coefficient) != 0)
        assert projection.actual_degree == (actual_motion_degree, 0)
        assert projection.degree_bound[0] >= actual_motion_degree and projection.degree_bound[1] == 0
        assert projection.degree_bound_id
        assert projection.normalized_coefficient_bytes
