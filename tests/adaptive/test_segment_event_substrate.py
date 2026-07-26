from __future__ import annotations

import hashlib
import math

import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2

import numpy as np


SQUARE = np.array(
    [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 10.0, 0.0],
        [0.0, 10.0, 0.0],
    ],
    dtype=np.float64,
)

LOWER_HALF_PLANE_WINDOW = np.array(
    [
        [-10.0, -10.0, 0.0],
        [10.0, -10.0, 0.0],
        [10.0, 3.0, 0.0],
        [-10.0, 3.0, 0.0],
    ],
    dtype=np.float64,
)

IDENTITY_ORDER_SQUARE = np.array(
    [
        [-4.0, -4.0, 0.0],
        [0.0, -4.0, 0.0],
        [0.0, 0.0, 0.0],
        [-4.0, 0.0, 0.0],
    ],
    dtype=np.float64,
)


def test_segment_source_exact_lifts_each_binary64_once() -> None:
    source = _continuous_tea_2.SegmentEventSource2.from_binary64(
        0.1,
        -0.5,
        1.25,
        2.0,
        0.25,
        4.0,
    )

    assert (source.x0.numerator, source.x0.denominator) == (
        "3602879701896397",
        "36028797018963968",
    )
    assert (source.y0.numerator, source.y0.denominator) == ("-1", "2")
    assert (source.x1.numerator, source.x1.denominator) == ("5", "4")
    assert (source.y1.numerator, source.y1.denominator) == ("2", "1")
    assert (source.tool_radius.numerator, source.tool_radius.denominator) == ("1", "4")
    assert (source.cap_chord_ratio.numerator, source.cap_chord_ratio.denominator) == ("4", "1")
    assert source.cap_chord_ratio.canonical_bytes
    assert source.motion_data == (
        "3602879701896397/36028797018963968",
        "-1/2",
        "5/4",
        "2",
    )
    assert source.canonical_digest == hashlib.sha256(source.canonical_bytes).digest()


@pytest.mark.parametrize(
    ("values", "error"),
    (
        ((math.nan, 0.0, 1.0, 0.0, 0.25, 1.0), "NonFiniteSegmentInputError"),
        ((0.0, 0.0, math.inf, 0.0, 0.25, 1.0), "NonFiniteSegmentInputError"),
        ((0.0, 0.0, 0.0, 0.0, 0.25, 1.0), "ZeroLengthSegmentMotionError"),
        ((0.0, 0.0, 1.0, 0.0, 0.0, 1.0), "NonPositiveToolRadiusError"),
        ((0.0, 0.0, 1.0, 0.0, -0.25, 1.0), "NonPositiveToolRadiusError"),
        ((0.0, 0.0, 1.0, 0.0, 0.25, 0.0), "InvalidCapChordRatioError"),
        ((0.0, 0.0, 1.0, 0.0, 0.25, 4.000000000000001), "InvalidCapChordRatioError"),
    ),
)
def test_segment_source_rejects_each_invalid_boundary_with_named_error(
    values: tuple[float, float, float, float, float, float],
    error: str,
) -> None:
    with pytest.raises(getattr(_continuous_tea_2, error)):
        _continuous_tea_2.SegmentEventSource2.from_binary64(*values)


def test_segment_partition_derives_line_pullbacks_and_complete_strata() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    partition = _continuous_tea_2.construct_segment_event_partition(
        stock,
        2.0,
        0.0,
        8.0,
        0.0,
        1.0,
        4.0,
    )

    assert partition.boundary_feature_ids == tuple(sorted(partition.boundary_feature_ids))
    assert {projection.degree_bound_id for projection in partition.projections} == {
        "segment-line-(1,2)-v1",
    }
    assert partition.strata
    assert all(stratum.active_branch_ids for stratum in partition.strata)
    assert all(
        {
            (pair.first_branch_id, pair.second_branch_id)
            for pair in stratum.pair_dispositions
        }
        == {
            (first, second)
            for first in stratum.active_branch_ids
            for second in stratum.active_branch_ids
            if first != second
        }
        for stratum in partition.strata
    )
    assert all(
        len({pair.pair_sheet_id for pair in stratum.pair_dispositions})
        == len(stratum.pair_dispositions)
        for stratum in partition.strata
    )
    assert {
        (
            pair.orientation_disposition,
            pair.cap_disposition,
        )
        for pair in partition.strata[0].pair_dispositions
    } == {("ccw-pi", "equal-cap")}
    assert partition.canonical_digest == hashlib.sha256(partition.canonical_bytes).digest()


def test_segment_partition_derives_circle_pullbacks_from_stock_records() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 1.0)

    partition = _continuous_tea_2.construct_segment_event_partition(
        stock,
        3.0,
        5.0,
        7.0,
        5.0,
        0.5,
        4.0,
    )

    assert {
        projection.degree_bound_id for projection in partition.projections
    } == {
        "segment-line-(1,2)-v1",
        "segment-circle-(2,4)-v1",
    }
    assert any(branch.support_kind == "circle" for branch in partition.branches)


def test_segment_partition_coalesces_simultaneous_vertex_events() -> None:
    partition = _continuous_tea_2.construct_segment_event_partition(
        _stock_2.Stock2(SQUARE, []),
        -2.0,
        0.0,
        0.0,
        0.0,
        1.0,
        4.0,
    )

    vertex_fibres = [
        stratum
        for stratum in partition.strata
        if stratum.kind == "fibre"
        and sum(event.kind == "endpoint-order" for event in stratum.events) >= 2
    ]
    assert len(vertex_fibres) == 1
    endpoint_events = tuple(
        event
        for event in vertex_fibres[0].events
        if event.kind == "endpoint-order"
    )
    assert len({event.vertex_id for event in endpoint_events}) == 1


def test_segment_verifier_reconstructs_source_and_rejects_mutation() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    partition = _continuous_tea_2.construct_segment_event_partition(
        stock,
        2.0,
        0.0,
        8.0,
        0.0,
        1.0,
        4.0,
    )

    verified = _continuous_tea_2.verify_segment_event_partition(
        stock,
        partition.source,
        partition,
    )
    mutated = _continuous_tea_2.mutate_segment_event_partition(
        partition,
        "delete-pair-disposition",
    )
    rejected = _continuous_tea_2.verify_segment_event_partition(
        stock,
        partition.source,
        mutated,
    )

    assert verified.verdict.name == "CERTIFIED"
    assert verified.partition.canonical_digest == partition.canonical_digest
    assert rejected.verdict.name == "UNRESOLVED_DEGENERACY"


def test_audit_red_cell_order_is_exact_ccw_not_feature_identity_order() -> None:
    source = _continuous_tea_2.SegmentEventSource2.from_binary64(
        -1.25,
        -1.0,
        -0.75,
        -1.0,
        1.25,
        4.0,
    )
    cell = _continuous_tea_2.construct_segment_cell_stratum(
        _stock_2.Stock2(IDENTITY_ORDER_SQUARE, []),
        source,
        "1",
        "2",
    )
    branches = {branch.branch_id: branch for branch in cell.branches}
    geometric_order = cell.stratum.active_branch_ids
    assert [
        (
            branches[branch_id].rim_chart_id,
            tuple(branches[branch_id].rim_factor_coefficients),
            branches[branch_id].rim_root_ordinal,
        )
        for branch_id in geometric_order
    ] == [
        ("rim-half-0-v1", ("1", "3"), 0),
        ("rim-half-0-v1", ("-1", "3"), 0),
        ("rim-half-0-v1", ("-1", "2"), 0),
        ("rim-half-1-v1", ("1", "2"), 0),
    ]
    assert geometric_order != tuple(sorted(geometric_order))
    ordered_pairs = {
        (pair.first_branch_id, pair.second_branch_id)
        for pair in cell.stratum.pair_dispositions
    }
    assert ordered_pairs == {
        (first, second)
        for first in geometric_order
        for second in geometric_order
        if first != second
    }


def test_audit_red_cell_pair_classifies_major_complement_exactly() -> None:
    assert _continuous_tea_2.segment_pair_literal_signs() == (
        "positive",
        "negative",
        "negative",
        "positive",
        "positive",
        "negative",
        "negative",
        "zero",
    )

    source = _continuous_tea_2.SegmentEventSource2.from_binary64(
        -1.0,
        0.0,
        1.0,
        0.0,
        5.0,
        4.0,
    )
    cell = _continuous_tea_2.construct_segment_cell_stratum(
        _stock_2.Stock2(LOWER_HALF_PLANE_WINDOW, []),
        source,
        "1",
        "2",
    )

    assert {
        (
            pair.orientation_disposition,
            pair.cap_disposition,
        )
        for pair in cell.stratum.pair_dispositions
    } == {
        ("ccw-major", "above-cap"),
        ("ccw-minor", "below-cap"),
    }


def test_segment_cell_branch_identity_is_station_independent() -> None:
    source = _continuous_tea_2.SegmentEventSource2.from_binary64(
        5.0,
        -1.0,
        5.0,
        0.0,
        1.0,
        4.0,
    )
    stock = _stock_2.Stock2(SQUARE, [])
    first = _continuous_tea_2.construct_segment_cell_stratum(
        stock,
        source,
        "1",
        "4",
    )
    second = _continuous_tea_2.construct_segment_cell_stratum(
        stock,
        source,
        "3",
        "4",
    )

    def by_sheet(
        cell: _continuous_tea_2.SegmentCellStratum2,
    ) -> dict[tuple[bytes, str, int], _continuous_tea_2.SegmentBoundaryBranch2]:
        return {
            (
                branch.feature_id,
                branch.rim_chart_id,
                branch.rim_sheet_ordinal,
            ): branch
            for branch in cell.branches
        }

    first_by_sheet = by_sheet(first)
    second_by_sheet = by_sheet(second)
    assert first_by_sheet.keys() == second_by_sheet.keys()
    assert all(
        first_by_sheet[sheet].branch_id == second_by_sheet[sheet].branch_id
        for sheet in first_by_sheet
    )
    assert any(
        (
            first_by_sheet[sheet].rim_factor_coefficients,
            first_by_sheet[sheet].rim_root_ordinal,
        )
        != (
            second_by_sheet[sheet].rim_factor_coefficients,
            second_by_sheet[sheet].rim_root_ordinal,
        )
        for sheet in first_by_sheet
    )


def test_segment_rational_square_root_normalization_cases() -> None:
    assert _continuous_tea_2.segment_rational_square_root_cases() == (
        "2/3",
        "non-square",
        "zero",
        "negative",
    )


def test_audit_red_fibre_preserves_pair_disposition_change() -> None:
    partition = _continuous_tea_2.construct_segment_event_partition(
        _stock_2.Stock2(SQUARE, []),
        5.0,
        -1.0,
        5.0,
        0.0,
        1.0,
        2.0,
    )
    changing_fibres = [
        stratum
        for stratum in partition.strata
        if stratum.kind == "fibre"
        and tuple(stratum.root_factor_coefficients) == ("1", "-4", "2")
        and stratum.root_ordinal == 0
    ]

    assert len(changing_fibres) == 1
    transition = changing_fibres[0]
    left = {pair.pair_sheet_id: pair for pair in transition.left_pair_dispositions}
    at = {pair.pair_sheet_id: pair for pair in transition.pair_dispositions}
    right = {pair.pair_sheet_id: pair for pair in transition.right_pair_dispositions}
    common_sheets = left.keys() & at.keys() & right.keys()
    assert any(
        (
            left[sheet].cap_disposition,
            at[sheet].cap_disposition,
            right[sheet].cap_disposition,
        )
        == ("below-cap", "equal-cap", "above-cap")
        for sheet in common_sheets
    )


def test_audit_red_equal_cardinality_reversal_is_merge_split_event() -> None:
    assert hasattr(_continuous_tea_2.PartitionEvent2, "left_active_count")
    assert hasattr(_continuous_tea_2.PartitionEvent2, "right_active_count")
    assert hasattr(
        _continuous_tea_2.PartitionEvent2,
        "incidence_permutation_rechecked",
    )

    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(5.0, 4.3, 5.0, 4.53, 0.05)
    forward = _continuous_tea_2.construct_segment_event_partition(
        stock,
        4.7,
        5.0,
        4.725,
        5.0,
        0.5,
        4.0,
    )
    reverse = _continuous_tea_2.construct_segment_event_partition(
        stock,
        4.725,
        5.0,
        4.7,
        5.0,
        0.5,
        4.0,
    )

    def equal_cardinality_endpoint_events(
        partition: _continuous_tea_2.SegmentEventPartition2,
    ) -> dict[bytes, _continuous_tea_2.PartitionEvent2]:
        return {
            event.vertex_id: event
            for stratum in partition.strata
            if stratum.kind == "fibre"
            for event in stratum.events
            if event.kind == "endpoint-order"
            and event.left_active_count == event.right_active_count
        }

    forward_events = equal_cardinality_endpoint_events(forward)
    reverse_events = equal_cardinality_endpoint_events(reverse)
    common_vertices = forward_events.keys() & reverse_events.keys()

    assert common_vertices
    assert all(
        {
            forward_events[vertex_id].disposition,
            reverse_events[vertex_id].disposition,
        }
        == {"merge", "split"}
        for vertex_id in common_vertices
    )
    assert all(
        forward_events[vertex_id].incidence_permutation_rechecked
        and reverse_events[vertex_id].incidence_permutation_rechecked
        for vertex_id in common_vertices
    )


def test_audit_red_fibre_is_evaluated_at_algebraic_root() -> None:
    partition = _continuous_tea_2.construct_segment_event_partition(
        _stock_2.Stock2(SQUARE, []),
        5.0,
        -1.0,
        5.0,
        0.0,
        1.0,
        2.0,
    )
    fibres = [
        stratum
        for stratum in partition.strata
        if stratum.kind == "fibre"
        and tuple(stratum.root_factor_coefficients) == ("1", "-4", "2")
    ]

    assert {stratum.root_ordinal for stratum in fibres} == {0}
    assert all(stratum.algebraic_root_evaluated for stratum in fibres)
    assert all(stratum.original_equations_rechecked for stratum in fibres)
    assert all(stratum.orientation_rechecked for stratum in fibres)
    assert all(stratum.trim_predicates_rechecked for stratum in fibres)


def test_audit_red_projection_completeness_rejects_missing_event_class() -> None:
    assert hasattr(
        _continuous_tea_2.SegmentEventStratum2,
        "root_factor_coefficients",
    )

    stock = _stock_2.Stock2(SQUARE, [])
    tangency = _continuous_tea_2.construct_segment_event_partition(
        stock,
        5.0,
        -2.0,
        5.0,
        0.0,
        1.0,
        4.0,
    )
    vertex = _continuous_tea_2.construct_segment_event_partition(
        stock,
        -2.0,
        0.0,
        0.0,
        0.0,
        1.0,
        4.0,
    )
    cap = _continuous_tea_2.construct_segment_event_partition(
        stock,
        5.0,
        -1.0,
        5.0,
        0.0,
        1.0,
        2.0,
    )
    overlap_stock = _stock_2.Stock2(SQUARE, [])
    overlap_stock.subtract_disk(5.0, 5.0, 1.0)
    overlap = _continuous_tea_2.construct_segment_event_partition(
        overlap_stock,
        5.0,
        5.0,
        5.0 + 2.0**-20,
        5.0,
        1.0,
        4.0,
    )

    def exact_fibre(
        partition: _continuous_tea_2.SegmentEventPartition2,
        factor: tuple[str, ...],
        ordinal: int,
    ) -> _continuous_tea_2.SegmentEventStratum2:
        return next(
            stratum
            for stratum in partition.strata
            if stratum.kind == "fibre"
            and tuple(stratum.root_factor_coefficients) == factor
            and stratum.root_ordinal == ordinal
        )

    tangency_fibre = exact_fibre(tangency, ("-1", "2"), 0)
    assert any(event.kind == "tangent" for event in tangency_fibre.events)
    assert any(
        pair.orientation_disposition == "ccw-zero"
        for pair in tangency_fibre.pair_dispositions
    )
    vertex_fibre = exact_fibre(vertex, ("-1", "2"), 0)
    assert any(event.kind == "endpoint-order" for event in vertex_fibre.events)
    cap_fibre = exact_fibre(cap, ("1", "-4", "2"), 0)
    assert any(
        pair.cap_disposition == "equal-cap"
        for pair in cap_fibre.pair_dispositions
    )
    overlap_fibre = exact_fibre(overlap, ("0", "1"), 0)
    assert any(event.kind == "support-overlap" for event in overlap_fibre.events)
    assert all(
        stratum.original_equations_rechecked
        and stratum.trim_predicates_rechecked
        for stratum in (
            tangency_fibre,
            vertex_fibre,
            cap_fibre,
            overlap_fibre,
        )
    )

    class_partitions = {
        "support-tangency": (stock, tangency),
        "vertex-passage": (stock, vertex),
        "support-overlap": (overlap_stock, overlap),
        "cap-crossing": (stock, cap),
        "endpoint-order-resultant": (stock, vertex),
        "pair-orientation-boundary": (stock, tangency),
    }
    for event_class, (event_stock, partition) in class_partitions.items():
        assert (
            _continuous_tea_2.verify_segment_event_partition(
                event_stock,
                partition.source,
                partition,
            ).verdict.name
            == "CERTIFIED"
        )
        mutated = _continuous_tea_2.mutate_segment_event_partition(
            partition,
            f"delete-projection-{event_class}",
        )
        assert mutated.canonical_bytes != partition.canonical_bytes
        rejected = _continuous_tea_2.verify_segment_event_partition(
            event_stock,
            partition.source,
            mutated,
        )
        assert rejected.verdict.name == "UNRESOLVED_DEGENERACY"
