import numpy as np

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)


def _stock_with_two_overlapping_disks() -> _stock_2.Stock2:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(4.5, 5.0, 1.0)
    stock.subtract_disk(5.5, 5.0, 1.0)
    return stock


def test_result_contract_has_no_boolean_compatibility() -> None:
    assert [verdict.name for verdict in _continuous_tea_2.ContinuousTeaVerdict] == [
        "CERTIFIED",
        "CAP_EXCEEDED",
        "UNRESOLVED_DEGENERACY",
    ]
    assert _continuous_tea_2.BOUNDARY_EVENT_KINDS == (
        "transverse",
        "tangent",
        "vertex",
        "overlap",
        "seam",
    )


def test_square_boundary_records_are_exact_sorted_and_stable() -> None:
    forward = _stock_2.Stock2(SQUARE, [])
    reversed_input = _stock_2.Stock2(SQUARE[::-1].copy(), [])
    records = _continuous_tea_2.extract_boundary_records(forward)
    reversed_records = _continuous_tea_2.extract_boundary_records(reversed_input)
    assert len(records) == 4
    assert [record.feature_id for record in records] == sorted(record.feature_id for record in records)
    assert [record.feature_id for record in records] == [record.feature_id for record in reversed_records]
    assert {record.support_kind for record in records} == {"line"}
    assert {record.material_side for record in records} == {"left"}
    assert all(record.source_exact and record.target_exact for record in records)
    assert all(record.source_vertex_id and record.target_vertex_id for record in records)
    assert all(record.trim_predicate for record in records)


def test_regularized_stock_exposes_line_and_circle_supports() -> None:
    records = _continuous_tea_2.extract_boundary_records(_stock_with_two_overlapping_disks())
    assert {record.support_kind for record in records} == {"line", "circle"}
    assert all(record.support_coefficients for record in records)


def test_trimmed_intersections_exclude_infinite_support_hits() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    events = _continuous_tea_2.extract_boundary_events(stock)
    assert len(events) == 4
    assert {event.kind for event in events} == {"vertex"}
    assert len({event.vertex_id for event in events}) == 4


def test_positive_length_overlap_is_not_flattened_to_tangent() -> None:
    stock = _stock_with_two_overlapping_disks()
    records = _continuous_tea_2.extract_boundary_records(stock)
    circle_index = next(index for index, record in enumerate(records) if record.support_kind == "circle")
    events = _continuous_tea_2.classify_boundary_pair(stock, circle_index, circle_index)
    assert [event.kind for event in events] == ["overlap"]


def test_quadratic_circle_vertex_identity_is_operand_order_invariant() -> None:
    stock = _stock_with_two_overlapping_disks()
    records = _continuous_tea_2.extract_boundary_records(stock)
    circle_indices = [index for index, record in enumerate(records) if record.support_kind == "circle"]
    found: tuple[int, int] | None = None
    for first in circle_indices:
        for second in circle_indices:
            if first >= second:
                continue
            events = _continuous_tea_2.classify_boundary_pair(stock, first, second)
            if any(event.kind == "vertex" for event in events):
                found = (first, second)
                break
        if found is not None:
            break
    assert found is not None
    first, second = found
    forward = _continuous_tea_2.classify_boundary_pair(stock, first, second)
    reverse = _continuous_tea_2.classify_boundary_pair(stock, second, first)
    assert [(event.kind, event.vertex_id) for event in forward] == [(event.kind, event.vertex_id) for event in reverse]
