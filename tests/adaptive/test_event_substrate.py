import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.identity import IncidentSupport
from compas_cgal.adaptive.identity import IncidentSupportIdV1
from compas_cgal.adaptive.identity import IntersectionBoundaryVertexIdV1
from compas_cgal.adaptive.identity import SupportKind
from compas_cgal.adaptive.identity import TrimIncidenceOrientation

CORPUS_PATH = Path(__file__).parent / "fixtures" / "event_corpus.json"
SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)


def _case(identifier: str) -> dict[str, object]:
    corpus = json.loads(CORPUS_PATH.read_text(encoding="utf-8"))
    return next(case for case in corpus["cases"] if case["id"] == identifier)


def _stock_with_disks(identifier: str, key: str = "disks") -> _stock_2.Stock2:
    geometry = _case(identifier)["geometry"]
    stock = _stock_2.Stock2(SQUARE, [])
    for disk in geometry[key]:
        stock.subtract_disk(*disk)
    return stock


def _stock_with_two_overlapping_disks() -> _stock_2.Stock2:
    return _stock_with_disks("circle-circle-regularization-vertex")


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
    assert issubclass(
        _continuous_tea_2.DegenerateBoundarySupportError,
        _continuous_tea_2.BoundaryExtractionError,
    )
    assert issubclass(
        _continuous_tea_2.MissingBoundaryEndpointError,
        _continuous_tea_2.BoundaryExtractionError,
    )
    assert issubclass(
        _continuous_tea_2.MissingBoundaryIntersectionError,
        _continuous_tea_2.BoundaryExtractionError,
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


def test_native_boundary_ids_byte_match_task1a_canonical_identity() -> None:
    records = _continuous_tea_2.extract_boundary_records(_stock_2.Stock2(SQUARE, []))
    support_by_id: dict[bytes, IncidentSupportIdV1] = {}
    incidents_by_vertex: dict[bytes, list[IncidentSupport]] = {}
    for record in records:
        support = IncidentSupportIdV1(
            SupportKind.LINE,
            tuple(int(value) for value in record.primitive_coefficients),
        )
        assert record.support_id == support.canonical_bytes
        support_by_id[record.support_id] = support
    for record in records:
        support = support_by_id[record.support_id]
        incidents_by_vertex.setdefault(record.source_vertex_id, []).append(IncidentSupport(support, TrimIncidenceOrientation.LEAVING))
        incidents_by_vertex.setdefault(record.target_vertex_id, []).append(IncidentSupport(support, TrimIncidenceOrientation.ENTERING))
    for vertex_id, incidents in incidents_by_vertex.items():
        expected = IntersectionBoundaryVertexIdV1.build(
            incident_supports=incidents,
            intersection_ordinal=0,
        )
        assert vertex_id == expected.canonical_bytes


def test_regularized_stock_exposes_line_and_circle_supports() -> None:
    records = _continuous_tea_2.extract_boundary_records(_stock_with_two_overlapping_disks())
    assert {record.support_kind for record in records} == {"line", "circle"}
    assert all(record.support_coefficients for record in records)


def test_trimmed_intersections_exclude_infinite_support_hits() -> None:
    pentagon = np.array(
        [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [6.0, 3.0, 0.0], [3.0, 6.0, 0.0], [0.0, 4.0, 0.0]],
        dtype=np.float64,
    )
    stock = _stock_2.Stock2(pentagon, [])
    records = _continuous_tea_2.extract_boundary_records(stock)
    events = _continuous_tea_2.extract_boundary_events(stock)
    assert len(events) == 5
    assert {event.kind for event in events} == {"vertex"}
    assert len({event.vertex_id for event in events}) == 5
    nonparallel_off_trim_pairs = []
    for first in range(len(records)):
        for second in range(first + 1, len(records)):
            first_vertices = {records[first].source_vertex_id, records[first].target_vertex_id}
            second_vertices = {records[second].source_vertex_id, records[second].target_vertex_id}
            a0, a1, _ = (int(value) for value in records[first].primitive_coefficients)
            b0, b1, _ = (int(value) for value in records[second].primitive_coefficients)
            if first_vertices.isdisjoint(second_vertices) and a0 * b1 != a1 * b0:
                nonparallel_off_trim_pairs.append((first, second))
    assert nonparallel_off_trim_pairs
    assert all(not _continuous_tea_2.classify_boundary_pair(stock, first, second) for first, second in nonparallel_off_trim_pairs)


def test_positive_length_overlap_is_not_flattened_to_tangent() -> None:
    stock = _stock_with_two_overlapping_disks()
    records = _continuous_tea_2.extract_boundary_records(stock)
    circle_index = next(index for index, record in enumerate(records) if record.support_kind == "circle")
    events = _continuous_tea_2.classify_boundary_pair(stock, circle_index, circle_index)
    assert [event.kind for event in events] == ["overlap"]
    assert events[0].exact_overlap_record
    assert require_canonical_record(events[0].overlap_id) == events[0].overlap_id
    assert events[0].overlap_id != events[0].exact_overlap_record
    assert events[0].multiplicity > 0
    assert records[circle_index].overlap_multiplicity > 0


def test_segment_and_full_circle_split_are_explicit_seams() -> None:
    segmented_square = np.array(
        [
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 0.0],
            [0.0, 10.0, 0.0],
        ],
        dtype=np.float64,
    )
    segmented_stock = _stock_2.Stock2(segmented_square, [])
    segment_records = _continuous_tea_2.extract_boundary_records(segmented_stock)
    segment_indices = [index for index, record in enumerate(segment_records) if tuple(record.primitive_coefficients) == ("0", "1", "0")]
    segment_events = _continuous_tea_2.classify_boundary_pair(
        segmented_stock,
        *segment_indices,
    )
    assert len(segment_events) == 1
    assert {event.kind for event in segment_events} == {"seam"}

    geometry = _case("segment-full-circle-parameter-seams")["geometry"]
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(*geometry["disk"])
    records = _continuous_tea_2.extract_boundary_records(stock)
    circle_indices = [index for index, record in enumerate(records) if record.support_kind == "circle"]
    assert len(circle_indices) == 2
    assert records[circle_indices[0]].support_id == records[circle_indices[1]].support_id
    events = _continuous_tea_2.classify_boundary_pair(stock, *circle_indices)
    assert len(events) == 2
    assert {event.kind for event in events} == {"seam"}


@pytest.mark.parametrize("key", ["external_disks", "internal_disk"])
def test_external_and_internal_tangencies_retain_multiplicity(key: str) -> None:
    geometry = _case("external-internal-tangencies")["geometry"]
    stock = _stock_2.Stock2(SQUARE, [])
    disks = geometry[key]
    if key == "internal_disk":
        disks = [disks]
    for disk in disks:
        stock.subtract_disk(*disk)
    tangent_events = [event for event in _continuous_tea_2.extract_boundary_events(stock) if event.multiplicity > 1]
    assert tangent_events
    assert {event.kind for event in tangent_events} == {"tangent"}


def test_regularization_vertices_cover_line_circle_and_two_disk_merge_split() -> None:
    line_circle_geometry = _case("line-circle-regularization-vertex")["geometry"]
    line_circle = _stock_2.Stock2(SQUARE, [])
    line_circle.subtract_disk(*line_circle_geometry["disk"])
    line_circle_records = _continuous_tea_2.extract_boundary_records(line_circle)
    support_kind = {record.feature_id: record.support_kind for record in line_circle_records}
    assert any(
        {support_kind[event.first_feature_id], support_kind[event.second_feature_id]} == {"line", "circle"}
        for event in _continuous_tea_2.extract_boundary_events(line_circle)
        if event.kind == "vertex"
    )

    two_disk = _stock_with_disks("two-disk-run-merge-split")
    two_disk_records = _continuous_tea_2.extract_boundary_records(two_disk)
    two_disk_kind = {record.feature_id: record.support_kind for record in two_disk_records}
    two_disk_support = {record.feature_id: record.support_id for record in two_disk_records}
    circle_vertices = [
        event
        for event in _continuous_tea_2.extract_boundary_events(two_disk)
        if event.kind == "vertex"
        and two_disk_kind[event.first_feature_id] == "circle"
        and two_disk_kind[event.second_feature_id] == "circle"
        and two_disk_support[event.first_feature_id] != two_disk_support[event.second_feature_id]
    ]
    assert len({event.vertex_id for event in circle_vertices}) == 2


def test_simultaneous_pair_incidences_are_not_collapsed() -> None:
    geometry = _case("three-simultaneous-intersections")["geometry"]
    stock = _stock_2.Stock2(SQUARE, [])
    for disk in geometry["tangent_disks"]:
        stock.subtract_disk(*disk)
    incidence_counts = Counter(event.vertex_id for event in _continuous_tea_2.extract_boundary_events(stock) if event.vertex_id)
    assert max(incidence_counts.values()) >= geometry["minimum_pair_events"]


def test_feature_index_failure_is_named() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(_continuous_tea_2.BoundaryFeatureIndexError):
        _continuous_tea_2.classify_boundary_pair(stock, 0, 99)


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

    geometry = _case("quadratic-circle-circle-stable-vertex")["geometry"]
    reversed_stock = _stock_2.Stock2(SQUARE, [])
    for disk in reversed(geometry["disks"]):
        reversed_stock.subtract_disk(*disk)
    forward_vertex_ids = sorted(event.vertex_id for event in _continuous_tea_2.extract_boundary_events(stock) if event.kind == "vertex")
    reversed_vertex_ids = sorted(event.vertex_id for event in _continuous_tea_2.extract_boundary_events(reversed_stock) if event.kind == "vertex")
    assert forward_vertex_ids == reversed_vertex_ids
