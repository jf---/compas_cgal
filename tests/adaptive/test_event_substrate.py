from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.identity import IncidentSupport
from compas_cgal.adaptive.identity import IncidentSupportIdV1
from compas_cgal.adaptive.identity import IntersectionBoundaryVertexIdV1
from compas_cgal.adaptive.identity import MultiIncidenceIntersectionBoundaryVertexIdV1
from compas_cgal.adaptive.identity import ParameterSeamBoundaryVertexIdV1
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
        assert require_canonical_record(record.source_incidence) == record.source_incidence
        assert require_canonical_record(record.target_incidence) == record.target_incidence
        assert record.source_incidence != record.target_incidence
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
    segment_support = IncidentSupportIdV1(
        SupportKind.LINE,
        tuple(int(value) for value in segment_records[segment_indices[0]].primitive_coefficients),
    )
    segment_incidents = (
        IncidentSupport(segment_support, TrimIncidenceOrientation.ENTERING),
        IncidentSupport(segment_support, TrimIncidenceOrientation.LEAVING),
    )
    assert (
        segment_events[0].vertex_id
        == ParameterSeamBoundaryVertexIdV1.build(
            incident_supports=segment_incidents,
            seam_ordinal=0,
        ).canonical_bytes
    )

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
    circle_record = records[circle_indices[0]]
    circle_support = IncidentSupportIdV1(
        SupportKind.CIRCLE,
        tuple(int(value) for value in circle_record.primitive_coefficients),
    )
    circle_incidents = (
        IncidentSupport(circle_support, TrimIncidenceOrientation.ENTERING),
        IncidentSupport(circle_support, TrimIncidenceOrientation.LEAVING),
    )
    assert {event.vertex_id for event in events} == {
        ParameterSeamBoundaryVertexIdV1.build(
            incident_supports=circle_incidents,
            seam_ordinal=ordinal,
        ).canonical_bytes
        for ordinal in (0, 1)
    }


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


def test_external_tangent_vertex_preserves_all_oriented_incidences() -> None:
    stock = _stock_with_disks("external-internal-tangencies", "external_disks")
    records = _continuous_tea_2.extract_boundary_records(stock)
    tangent = next(event for event in _continuous_tea_2.extract_boundary_events(stock) if event.kind == "tangent")
    incidents: list[IncidentSupport] = []
    for record in records:
        support = IncidentSupportIdV1(
            SupportKind.LINE if record.support_kind == "line" else SupportKind.CIRCLE,
            tuple(int(value) for value in record.primitive_coefficients),
        )
        if record.source_vertex_id == tangent.vertex_id:
            incident = IncidentSupport(support, TrimIncidenceOrientation.LEAVING)
            assert record.source_incidence == incident.canonical_bytes
            incidents.append(incident)
        if record.target_vertex_id == tangent.vertex_id:
            incident = IncidentSupport(support, TrimIncidenceOrientation.ENTERING)
            assert record.target_incidence == incident.canonical_bytes
            incidents.append(incident)

    assert len(incidents) == 4
    assert (
        tangent.vertex_id
        == MultiIncidenceIntersectionBoundaryVertexIdV1.build(
            incident_supports=incidents,
            intersection_ordinal=0,
        ).canonical_bytes
    )


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


def _root_id(coefficients: tuple[int, ...], ordinal: int) -> bytes:
    return encode_tagged_union(
        b"algebraic-root-id-v1",
        encode_component_map(
            {
                b"coefficients": encode_sequence(tuple(encode_integer(value) for value in coefficients)),
                b"root-ordinal": encode_integer(ordinal),
            }
        ),
    )


def _event(
    suffix: bytes,
    *,
    kind: str = "tangency",
    disposition: str = "accepted",
) -> _continuous_tea_2.PartitionEvent2:
    return _continuous_tea_2.PartitionEvent2(
        kind,
        b"feature-" + suffix,
        b"support-" + suffix,
        b"trim-" + suffix,
        b"vertex-" + suffix,
        b"branch-" + suffix,
        disposition,
    )


def _projection_input(
    identifier: str,
    coefficients: tuple[str, ...],
    events: tuple[_continuous_tea_2.PartitionEvent2, ...],
) -> _continuous_tea_2.ProjectionInput2:
    return _continuous_tea_2.ProjectionInput2(
        identifier,
        coefficients,
        events,
    )


def _cell_witnesses(
    certificate: _continuous_tea_2.EventPartitionCertificate2,
) -> tuple[Fraction, ...]:
    return tuple(
        Fraction(
            int(cell.witness_numerator),
            int(cell.witness_denominator),
        )
        for cell in certificate.cells
    )


def test_locked_algebraic_backend_evidence_is_immutable_and_exact() -> None:
    evidence = _continuous_tea_2.exact_algebraic_backend_evidence()

    assert evidence.cgal_version == "6.0.1"
    assert evidence.integer_backend == "CORE::BigInt(boost::multiprecision::cpp_int)"
    assert evidence.algebraic_kernel_1 == "CGAL::Algebraic_kernel_d_1<CORE::BigInt>"
    assert evidence.algebraic_kernel_2 == "CGAL::Algebraic_kernel_d_2<CORE::BigInt>"
    assert evidence.arrangement_traits == "CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>"
    assert evidence.compile_definitions == (
        "CGAL_DISABLE_GMP=1",
        "CGAL_USE_BOOST_MP=1",
        "CGAL_USE_CORE=1",
        "CGAL_CORE_USE_BOOST_BACKEND=1",
    )
    with pytest.raises(AttributeError):
        evidence.cgal_version = "changed"


def test_frozen_center_and_rim_charts_cover_each_seam_once() -> None:
    charts = _continuous_tea_2.parameter_charts()

    assert [chart.chart_id for chart in charts] == [
        "center-quarter-0-v1",
        "center-quarter-1-v1",
        "center-quarter-2-v1",
        "center-quarter-3-v1",
        "rim-half-0-v1",
        "rim-half-1-v1",
    ]
    assert [(chart.domain_low, chart.domain_high, chart.orientation) for chart in charts] == [
        ("0", "1", "ccw"),
        ("0", "1", "ccw"),
        ("0", "1", "ccw"),
        ("0", "1", "ccw"),
        ("-1", "1", "ccw"),
        ("-1", "1", "ccw"),
    ]
    seam_owners: Counter[bytes] = Counter()
    for chart in charts:
        assert chart.owns_start_seam
        assert not chart.owns_end_seam
        seam_owners[chart.start_seam_id] += 1
        assert chart.end_seam_id in {candidate.start_seam_id for candidate in charts}
    assert sorted(seam_owners.values()) == [1, 1, 1, 1, 1, 1]

    verified = _continuous_tea_2.verify_chart_coverage(charts)
    assert verified.verdict.name == "CERTIFIED"
    assert len(verified.partition.seams) == 6


@pytest.mark.parametrize(
    (
        "motion_kind",
        "motion_data",
        "support_kind",
        "support_data",
        "center_chart",
        "expected_rows",
        "expected_degree",
        "bound_id",
    ),
    (
        (
            "segment",
            ("0", "0", "1", "0"),
            "line",
            ("1", "0", "0"),
            "",
            (("1", "0", "-1"), ("1", "0", "1")),
            (1, 2),
            "segment-line-(1,2)-v1",
        ),
        (
            "segment",
            ("-1", "0", "1", "0"),
            "circle",
            ("0", "0", "1"),
            "",
            (
                ("-1", "0", "2", "0", "3"),
                ("0", "0", "-8", "0", "-8"),
                ("4", "0", "8", "0", "4"),
            ),
            (2, 4),
            "segment-circle-(2,4)-v1",
        ),
        (
            "full-circle",
            ("0", "0", "1"),
            "line",
            ("1", "0", "0"),
            "center-quarter-0-v1",
            (("1", "0", "0"), ("0", "0", "0"), ("0", "0", "-1")),
            (2, 2),
            "full-circle-line-(2,2)-v1",
        ),
        (
            "full-circle",
            ("0", "0", "1"),
            "circle",
            ("3", "0", "1"),
            "center-quarter-0-v1",
            (
                ("0", "0", "1", "0", "1"),
                ("0", "1", "0", "1", "0"),
                ("1", "0", "5", "0", "4"),
                ("0", "1", "0", "1", "0"),
                ("1", "0", "4", "0", "3"),
            ),
            (4, 4),
            "full-circle-circle-(4,4)-v1",
        ),
    ),
)
def test_exact_pullbacks_match_literal_coefficient_grids_and_bounds(
    motion_kind: str,
    motion_data: tuple[str, ...],
    support_kind: str,
    support_data: tuple[str, ...],
    center_chart: str,
    expected_rows: tuple[tuple[str, ...], ...],
    expected_degree: tuple[int, int],
    bound_id: str,
) -> None:
    projection = _continuous_tea_2.construct_pullback(
        motion_kind,
        motion_data,
        support_kind,
        support_data,
        "1",
        center_chart,
        "rim-half-0-v1",
    )

    assert projection.coefficient_rows == expected_rows
    assert projection.actual_degree == expected_degree
    assert projection.degree_bound == expected_degree
    assert projection.degree_bound_id == bound_id
    assert projection.normalized_coefficient_bytes


def test_projection_degree_overflow_fails_before_partitioning() -> None:
    with pytest.raises(
        _continuous_tea_2.ProjectionDegreeBoundError,
        match="degree bound",
    ):
        _continuous_tea_2.projection_from_grid(
            "overflow",
            (("1", "0", "0", "1"), ("0",), ("1",)),
            "segment-line-(1,2)-v1",
        )


def test_repeated_tangency_root_is_one_fibre_with_multiplicity_two() -> None:
    certificate = _continuous_tea_2.partition_projections(
        (
            _projection_input(
                "repeated",
                ("1", "-4", "4"),
                (_event(b"repeated"),),
            ),
        )
    )

    assert len(certificate.roots) == 1
    root = certificate.roots[0]
    assert root.root_id == _root_id((-1, 2), 0)
    assert root.factor_coefficients == ("-1", "2")
    assert root.multiplicity == 2
    assert Fraction(root.interval_low) < Fraction(1, 2) < Fraction(root.interval_high)
    assert _cell_witnesses(certificate) == (Fraction(1, 4), Fraction(3, 4))
    assert len(certificate.fibres) == 1
    assert len(certificate.fibres[0].events) == 1


def test_same_root_from_reducible_and_irreducible_projections_merges() -> None:
    certificate = _continuous_tea_2.partition_projections(
        (
            _projection_input(
                "irreducible",
                ("-1", "2"),
                (_event(b"irreducible"),),
            ),
            _projection_input(
                "reducible",
                ("-1", "1", "2"),
                (_event(b"reducible"),),
            ),
        )
    )

    assert [root.root_id for root in certificate.roots] == [
        _root_id((-1, 2), 0),
    ]
    assert [root.factor_coefficients for root in certificate.roots] == [
        ("-1", "2"),
    ]
    assert [event.feature_id for event in certificate.fibres[0].events] == [
        b"feature-irreducible",
        b"feature-reducible",
    ]


def test_irreducible_quadratic_owns_its_algebraic_root_identity() -> None:
    certificate = _continuous_tea_2.partition_projections(
        (
            _projection_input(
                "irreducible-quadratic",
                ("-1", "1", "1"),
                (_event(b"quadratic"),),
            ),
        )
    )

    assert [root.root_id for root in certificate.roots] == [
        _root_id((-1, 1, 1), 1),
    ]
    assert certificate.roots[0].factor_coefficients == (
        "-1",
        "1",
        "1",
    )


def test_reducible_quartic_without_rational_roots_splits_exactly() -> None:
    certificate = _continuous_tea_2.partition_projections(
        (
            _projection_input(
                "quadratic-pair",
                ("1", "0", "-5", "0", "6"),
                (_event(b"quadratic-pair"),),
            ),
        )
    )

    assert certificate.projections[0].factor_coefficients == (
        ("-1", "0", "2"),
        ("-1", "0", "3"),
    )
    assert [root.root_id for root in certificate.roots] == [
        _root_id((-1, 0, 3), 1),
        _root_id((-1, 0, 2), 1),
    ]


def test_degree_above_frozen_task5_contract_fails_before_solving() -> None:
    with pytest.raises(
        _continuous_tea_2.UnsupportedAlgebraicDegreeError,
        match="degree 4",
    ):
        _continuous_tea_2.partition_projections(
            (
                _projection_input(
                    "degree-five",
                    ("1", "0", "0", "0", "0", "1"),
                    (_event(b"degree-five"),),
                ),
            )
        )
