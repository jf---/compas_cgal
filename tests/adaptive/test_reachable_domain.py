import math
from dataclasses import FrozenInstanceError
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest
from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _coverage_2
from compas_cgal.adaptive import reachable_domain as reachable_domain_module
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.errors import InvalidReachableDomainInputError
from compas_cgal.adaptive.errors import InvalidReachableDomainCertificateError
from compas_cgal.adaptive.errors import PocketNotMachinableError
from compas_cgal.adaptive.errors import ReachableArrangementTopologyError
from compas_cgal.adaptive.errors import ReachableMaterialContainmentError
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.reachable_domain import ReachableDomainCertificate
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


def _array(points: tuple[tuple[float, float], ...]) -> np.ndarray:
    return np.asarray(tuple((x, y, 0.0) for x, y in points), dtype=np.float64)


RECTANGLE = ((0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (0.0, 8.0))
ACUTE = ((0.0, 0.0), (10.0, 0.0), (2.0, 7.0))
REFLEX = (
    (0.0, 0.0),
    (9.0, 0.0),
    (9.0, 3.0),
    (4.0, 3.0),
    (4.0, 9.0),
    (0.0, 9.0),
)
NARROW_NECK = (
    (0.0, 0.0),
    (4.0, 0.0),
    (4.0, 1.5),
    (8.0, 1.5),
    (8.0, 0.0),
    (12.0, 0.0),
    (12.0, 6.0),
    (8.0, 6.0),
    (8.0, 4.5),
    (4.0, 4.5),
    (4.0, 6.0),
    (0.0, 6.0),
)
DISCONNECTED_NECK = (
    (0.0, 0.0),
    (4.0, 0.0),
    (4.0, 2.25),
    (8.0, 2.25),
    (8.0, 0.0),
    (12.0, 0.0),
    (12.0, 6.0),
    (8.0, 6.0),
    (8.0, 3.75),
    (4.0, 3.75),
    (4.0, 6.0),
    (0.0, 6.0),
)
ISLAND = ((4.0, 3.0), (6.0, 3.0), (6.0, 5.0), (4.0, 5.0))
ACUTE_CERTIFICATE_DIGEST = bytes.fromhex("54f55334eb9a9a8d2b6b955f5ec8bf7fd57a2f41b09821d09c7c7e6c79876049")


def _ring(
    points: tuple[tuple[float, float], ...],
    *,
    outer: bool,
) -> CanonicalRingV1:
    framed = tuple(Point2[WorldXY].build(x, y) for x, y in points)
    return CanonicalRingV1.build_outer(framed) if outer else CanonicalRingV1.build_hole(framed)


def _domain(
    boundary: tuple[tuple[float, float], ...] = RECTANGLE,
    holes: tuple[tuple[tuple[float, float], ...], ...] = (),
    radius: float = 1.0,
) -> ReachableDomain:
    return ReachableDomain.build(
        design_boundary=_ring(boundary, outer=True),
        holes=tuple(_ring(hole, outer=False) for hole in holes),
        tool_radius=ToolRadius.build(radius),
    )


def _record_fields(record: bytes, tag: bytes) -> tuple[bytes, ...]:
    prefix = tag + b"\x00"
    assert record.startswith(prefix)
    offset = len(prefix)
    count = int.from_bytes(record[offset : offset + 8], "big")
    offset += 8
    fields = []
    for _ in range(count):
        size = int.from_bytes(record[offset : offset + 8], "big")
        offset += 8
        fields.append(record[offset : offset + size])
        offset += size
    assert offset == len(record)
    return tuple(fields)


def _tagged_record(tag: bytes, fields: tuple[bytes, ...]) -> bytes:
    return tag + b"\x00" + len(fields).to_bytes(8, "big") + b"".join(len(field).to_bytes(8, "big") + field for field in fields)


def _replace_sorted_record(
    records: tuple[bytes, ...],
    original: bytes,
    replacement: bytes,
) -> tuple[bytes, ...]:
    return tuple(sorted(replacement if record == original else record for record in records))


def _replace_first_cycle_element(
    cell: bytes,
    element_fields: tuple[bytes, ...],
) -> bytes:
    cell_fields = _record_fields(cell, b"selected-cell-v1")
    cycle_fields = _record_fields(
        cell_fields[0],
        b"selected-boundary-cycle-v1",
    )
    elements = list(_record_fields(cycle_fields[1], b"cycle-elements-v1"))
    elements[0] = _tagged_record(
        b"selected-cycle-element-v1",
        element_fields,
    )
    outer = _tagged_record(
        b"selected-boundary-cycle-v1",
        (
            cycle_fields[0],
            _tagged_record(b"cycle-elements-v1", tuple(elements)),
        ),
    )
    return _tagged_record(b"selected-cell-v1", (outer, cell_fields[1]))


def _malformed_source_certificate(
    certificate: ReachableDomainCertificate,
    field_index: int,
    replacement: bytes,
) -> ReachableDomainCertificate:
    source = certificate.source_curve_records[0]
    fields = list(_record_fields(source, b"source-curve-v2"))
    fields[field_index] = replacement
    malformed = _tagged_record(b"source-curve-v2", tuple(fields))
    return replace(
        certificate,
        source_curve_records=_replace_sorted_record(
            certificate.source_curve_records,
            source,
            malformed,
        ),
    )


def _globally_relabelled_source_certificate(
    certificate: ReachableDomainCertificate,
) -> ReachableDomainCertificate:
    source = next(record for record in certificate.source_curve_records if _record_fields(record, b"source-curve-v2")[3] == b"offset-minus")
    fields = list(_record_fields(source, b"source-curve-v2"))
    fields[3] = b"invalid-role"
    malformed = _tagged_record(b"source-curve-v2", tuple(fields))

    def relabel(records: tuple[bytes, ...]) -> tuple[bytes, ...]:
        return tuple(sorted(record.replace(source, malformed) for record in records))

    return replace(
        certificate,
        source_curve_records=relabel(certificate.source_curve_records),
        arrangement_vertex_records=relabel(certificate.arrangement_vertex_records),
        selected_cell_records=relabel(certificate.selected_cell_records),
        component_records=relabel(certificate.component_records),
    )


def _malformed_cell_certificate(
    certificate: ReachableDomainCertificate,
    field_index: int,
    replacement: bytes,
) -> ReachableDomainCertificate:
    cell = certificate.selected_cell_records[0]
    cell_fields = _record_fields(cell, b"selected-cell-v1")
    cycle_fields = _record_fields(
        cell_fields[0],
        b"selected-boundary-cycle-v1",
    )
    elements = _record_fields(cycle_fields[1], b"cycle-elements-v1")
    element_fields = list(_record_fields(elements[0], b"selected-cycle-element-v1"))
    element_fields[field_index] = replacement
    malformed = _replace_first_cycle_element(cell, tuple(element_fields))
    return replace(
        certificate,
        selected_cell_records=_replace_sorted_record(
            certificate.selected_cell_records,
            cell,
            malformed,
        ),
    )


def _malformed_component_certificate(
    certificate: ReachableDomainCertificate,
    *,
    cells: bytes | None = None,
    adjacencies: bytes | None = None,
    boundary_cycles: bytes | None = None,
) -> ReachableDomainCertificate:
    component = certificate.component_records[0]
    fields = list(_record_fields(component, b"selected-component-v1"))
    replacements = (cells, adjacencies, boundary_cycles)
    for index, replacement in enumerate(replacements):
        if replacement is not None:
            fields[index] = replacement
    malformed = _tagged_record(b"selected-component-v1", tuple(fields))
    return replace(certificate, component_records=(malformed,))


def _assert_certificate_rejected(
    certificate: ReachableDomainCertificate,
    match: str,
) -> None:
    with pytest.raises(InvalidReachableDomainCertificateError, match=match):
        ReachableDomainCertificate.validate(certificate)
    with pytest.raises(InvalidReachableDomainCertificateError, match=match):
        _ = certificate.digest


@pytest.mark.parametrize(
    ("boundary", "holes"),
    (
        (RECTANGLE, ()),
        (ACUTE, ()),
        (REFLEX, ()),
        (NARROW_NECK, ()),
        (RECTANGLE, (ISLAND,)),
    ),
)
def test_reachable_fixtures_have_exact_closed_certificate(
    boundary: tuple[tuple[float, float], ...],
    holes: tuple[tuple[tuple[float, float], ...], ...],
) -> None:
    domain = _domain(boundary, holes)
    certificate = domain.certificate

    assert certificate.exact_cell_selection
    assert certificate.complete_source_provenance
    assert certificate.reachable_subset_of_design
    assert certificate.source_curve_records == tuple(sorted(certificate.source_curve_records))
    assert certificate.arrangement_vertex_records == tuple(sorted(certificate.arrangement_vertex_records))
    assert certificate.arrangement_vertex_records
    assert certificate.selected_cell_records == tuple(sorted(certificate.selected_cell_records))
    assert certificate.component_records == tuple(sorted(certificate.component_records))
    assert len(certificate.source_curve_records) == 3 * (len(boundary) + sum(len(hole) for hole in holes))
    assert certificate.selected_cell_records
    assert len(certificate.component_records) == 1
    assert len(certificate.digest) == 32
    assert len(certificate.center_domain_digest) == 32
    assert len(certificate.reachable_material_digest) == 32
    assert len(certificate.unreachable_residual_digest) == 32


def test_disconnected_center_domain_fails_phase_one_entry_contract() -> None:
    with pytest.raises(PocketNotMachinableError, match="connected"):
        _domain(DISCONNECTED_NECK)


def test_rectangle_center_domain_accepts_exact_tangency_and_rejects_across_each_boundary() -> None:
    domain = _domain()
    center = domain.center_domain
    tangent_points = ((1.0, 4.0), (9.0, 4.0), (5.0, 1.0), (5.0, 7.0))
    across_points = (
        (np.nextafter(1.0, 0.0), 4.0),
        (np.nextafter(9.0, 10.0), 4.0),
        (5.0, np.nextafter(1.0, 0.0)),
        (5.0, np.nextafter(7.0, 8.0)),
    )

    assert all(center.contains(x, y) for x, y in tangent_points)
    assert not any(center.contains(x, y) for x, y in across_points)


def test_island_center_domain_uses_exact_disk_containment_on_hole_boundary() -> None:
    domain = _domain(RECTANGLE, (ISLAND,))
    center = domain.center_domain

    assert center.contains(3.0, 4.0)
    assert not center.contains(np.nextafter(3.0, 4.0), 4.0)
    assert center.contains(7.0, 4.0)
    assert not center.contains(np.nextafter(7.0, 6.0), 4.0)


@pytest.mark.parametrize(
    ("boundary", "holes", "queries"),
    (
        pytest.param(
            RECTANGLE,
            (),
            ((5.0, 4.0), (1.0, 4.0), (0.5, 4.0), (11.0, 4.0)),
            id="convex",
        ),
        pytest.param(
            REFLEX,
            (),
            ((2.0, 2.0), (8.0, 1.5), (4.5, 4.0), (3.0, 8.0)),
            id="concave",
        ),
        pytest.param(
            RECTANGLE,
            (ISLAND,),
            ((2.0, 4.0), (3.0, 4.0), (3.5, 4.0), (5.0, 4.0), (7.0, 4.0), (8.0, 4.0)),
            id="holed",
        ),
    ),
)
def test_lightweight_cutter_centre_domain_matches_exact_region(
    boundary: tuple[tuple[float, float], ...],
    holes: tuple[tuple[tuple[float, float], ...], ...],
    queries: tuple[tuple[float, float], ...],
) -> None:
    lightweight = _coverage_2.CutterCentreDomain2.build(
        _array(boundary),
        [_array(hole) for hole in holes],
        1.0,
    )
    exact_region = _coverage_2.ReachableDomain2(
        _array(boundary),
        [_array(hole) for hole in holes],
        1.0,
    ).center_domain()

    assert tuple(lightweight.contains(x, y) for x, y in queries) == tuple(exact_region.contains(x, y) for x, y in queries)


def test_lightweight_cutter_centre_domain_keeps_exact_closed_boundary() -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(RECTANGLE), [], 1.0)

    assert domain.contains(5.0, 4.0)
    assert domain.contains(1.0, 4.0)
    assert not domain.contains(np.nextafter(1.0, 0.0), 4.0)
    assert not domain.contains(0.5, 4.0)
    assert not domain.contains(-1.0, 4.0)


def test_lightweight_cutter_centre_domain_keeps_exact_hole_boundary() -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(RECTANGLE), [_array(ISLAND)], 1.0)

    assert domain.contains(3.0, 4.0)
    assert not domain.contains(np.nextafter(3.0, 4.0), 4.0)
    assert domain.contains(7.0, 4.0)
    assert not domain.contains(np.nextafter(7.0, 6.0), 4.0)
    assert not domain.contains(5.0, 4.0)


@pytest.mark.parametrize(
    ("boundary", "radius"),
    (
        pytest.param(RECTANGLE[:2], 1.0, id="short-ring"),
        pytest.param(RECTANGLE, 0.0, id="zero-radius"),
        pytest.param(RECTANGLE, float("nan"), id="nonfinite-radius"),
        pytest.param(((0.0, 0.0), (10.0, 0.0), (float("nan"), 8.0), (0.0, 8.0)), 1.0, id="nonfinite-vertex"),
    ),
)
def test_lightweight_cutter_centre_domain_preserves_named_input_errors(
    boundary: tuple[tuple[float, float], ...],
    radius: float,
) -> None:
    with pytest.raises(InvalidReachableDomainInputError):
        _coverage_2.CutterCentreDomain2.build(_array(boundary), [], radius)


def test_lightweight_cutter_centre_domain_rejects_invalid_hole_relationship() -> None:
    outside_hole = ((11.0, 3.0), (13.0, 3.0), (13.0, 5.0), (11.0, 5.0))

    with pytest.raises(InvalidReachableDomainInputError):
        _coverage_2.CutterCentreDomain2.build(_array(RECTANGLE), [_array(outside_hole)], 1.0)


def test_lightweight_predicate_does_not_claim_global_connectivity_certification() -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(DISCONNECTED_NECK), [], 1.0)

    assert domain.contains(2.0, 3.0)
    assert domain.contains(10.0, 3.0)


@pytest.mark.parametrize(
    ("x", "y"),
    (
        pytest.param(float("nan"), 4.0, id="nan-x"),
        pytest.param(5.0, float("nan"), id="nan-y"),
        pytest.param(float("inf"), 4.0, id="positive-infinity"),
        pytest.param(float("-inf"), 4.0, id="negative-infinity"),
    ),
)
def test_lightweight_cutter_centre_domain_rejects_nonfinite_queries(x: float, y: float) -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(RECTANGLE), [], 1.0)

    with pytest.raises(InvalidReachableDomainInputError):
        domain.contains(x, y)


def test_lightweight_cutter_centre_domain_disables_direct_construction() -> None:
    with pytest.raises(TypeError):
        _coverage_2.CutterCentreDomain2(_array(RECTANGLE), [], 1.0)


def test_lightweight_cutter_centre_domain_preserves_oblique_exact_tangency() -> None:
    boundary = ((0.0, 0.0), (4.0, 3.0), (0.0, 3.0))
    tangent = (1.625, 2.0)
    adjacent_illegal = (np.nextafter(tangent[0], 4.0), tangent[1])
    lightweight = _coverage_2.CutterCentreDomain2.build(_array(boundary), [], 0.625)
    exact_region = _coverage_2.ReachableDomain2(_array(boundary), [], 0.625).center_domain()

    assert lightweight.contains(*tangent)
    assert exact_region.contains(*tangent)
    assert not lightweight.contains(*adjacent_illegal)
    assert not exact_region.contains(*adjacent_illegal)


@pytest.mark.parametrize(
    ("boundary", "holes"),
    (
        pytest.param(RECTANGLE, (), id="convex"),
        pytest.param(REFLEX, (), id="concave"),
        pytest.param(RECTANGLE, (ISLAND,), id="holed"),
        pytest.param(NARROW_NECK, (), id="narrow-bay"),
    ),
)
def test_material_predicate_matches_legacy_exact_reachable_material_pointwise(
    boundary: tuple[tuple[float, float], ...],
    holes: tuple[tuple[tuple[float, float], ...], ...],
) -> None:
    native_holes = [_array(hole) for hole in holes]
    expected = _coverage_2.ReachableDomain2(_array(boundary), native_holes, 1.0).reachable_material()
    actual = _coverage_2.ReachableMaterialPredicate2.build(_array(boundary), native_holes, 1.0)
    xs = np.linspace(min(x for x, _ in boundary) - 1.0, max(x for x, _ in boundary) + 1.0, 17)
    ys = np.linspace(min(y for _, y in boundary) - 1.0, max(y for _, y in boundary) + 1.0, 17)

    assert tuple(actual.contains(x, y) for x in xs for y in ys) == tuple(expected.contains(x, y) for x in xs for y in ys)


def test_material_predicate_is_factory_only_and_preserves_named_errors() -> None:
    with pytest.raises(TypeError):
        _coverage_2.ReachableMaterialPredicate2(_array(RECTANGLE), [], 1.0)
    with pytest.raises(InvalidReachableDomainInputError):
        _coverage_2.ReachableMaterialPredicate2.build(_array(RECTANGLE[:2]), [], 1.0)
    outside_hole = ((11.0, 3.0), (13.0, 3.0), (13.0, 5.0), (11.0, 5.0))
    with pytest.raises(InvalidReachableDomainInputError):
        _coverage_2.ReachableMaterialPredicate2.build(_array(RECTANGLE), [_array(outside_hole)], 1.0)
    with pytest.raises(PocketNotMachinableError):
        _coverage_2.ReachableMaterialPredicate2.build(_array(DISCONNECTED_NECK), [], 1.0)

    predicate = _coverage_2.ReachableMaterialPredicate2.build(_array(RECTANGLE), [], 1.0)
    with pytest.raises(InvalidReachableDomainInputError):
        predicate.contains(float("nan"), 0.0)


def test_reachable_material_and_residual_are_exact_partition_of_design() -> None:
    native = _coverage_2.ReachableDomain2(_array(RECTANGLE), [], 1.0)
    center = native.center_domain()
    reachable = native.reachable_material()
    residual = native.unreachable_residual()
    certificate = native.certificate()

    assert center.component_count() == 1
    assert reachable.is_subset_of(native.design_region())
    assert not residual.is_empty()
    assert residual.component_count() == 4
    assert certificate.reachable_subset_of_design
    assert not reachable.contains(0.0, 0.0)
    assert residual.contains(0.0, 0.0)


def test_native_regions_are_opaque_immutable_clones() -> None:
    native = _coverage_2.ReachableDomain2(_array(RECTANGLE), [], 1.0)
    first = native.reachable_material()
    second = native.reachable_material()

    assert type(first) is _coverage_2.ExactRegion2
    assert first is not second
    assert first.exactly_equals(second)
    assert not hasattr(first, "difference")
    assert not hasattr(first, "join")
    with pytest.raises(TypeError):
        _coverage_2.ExactRegion2()


def test_structural_identity_is_invariant_to_ring_and_hole_insertion_order() -> None:
    outer_rotated = RECTANGLE[2:] + RECTANGLE[:2]
    outer_reversed = tuple(reversed(RECTANGLE))
    left_island = ((2.5, 3.0), (3.5, 3.0), (3.5, 4.0), (2.5, 4.0))
    right_island = ((6.5, 4.0), (7.5, 4.0), (7.5, 5.0), (6.5, 5.0))
    variants = (
        _domain(RECTANGLE, (left_island, right_island)),
        _domain(outer_rotated, (right_island, left_island)),
        _domain(outer_reversed, (tuple(reversed(left_island)), tuple(reversed(right_island)))),
    )

    assert all(domain.center_domain.component_count() == 1 for domain in variants)
    assert len({domain.certificate.digest for domain in variants}) == 1
    assert len({domain.certificate.source_curve_records for domain in variants}) == 1
    assert len({domain.certificate.arrangement_vertex_records for domain in variants}) == 1
    assert len({domain.certificate.selected_cell_records for domain in variants}) == 1
    assert len({domain.certificate.component_records for domain in variants}) == 1
    assert variants[0].center_domain.exactly_equals(variants[1].center_domain)
    assert variants[0].reachable_material.exactly_equals(variants[2].reachable_material)


@pytest.mark.parametrize(
    "boundary",
    (
        pytest.param(ACUTE, id="baseline"),
        pytest.param(ACUTE[1:] + ACUTE[:1], id="rotated"),
        pytest.param(tuple(reversed(ACUTE)), id="reversed"),
    ),
)
def test_circle_vertical_trim_extrema_have_distinct_stable_structural_ordinals(
    boundary: tuple[tuple[float, float], ...],
) -> None:
    certificate = _domain(boundary).certificate
    assert certificate.digest == ACUTE_CERTIFICATE_DIGEST

    repeated_source_members: set[tuple[bytes, bytes]] = set()
    for vertex in certificate.arrangement_vertex_records:
        incidence, ordinal = _record_fields(
            vertex,
            b"arrangement-vertex-v1",
        )
        pieces = _record_fields(
            incidence,
            b"arrangement-incidence-multiset-v1",
        )
        source_ids = tuple(
            _record_fields(
                piece,
                b"source-piece-v1",
            )[0]
            for piece in pieces
        )
        if len(source_ids) == 2 and source_ids[0] == source_ids[1]:
            repeated_source_members.add((incidence, ordinal))
    assert len(repeated_source_members) >= 2


def test_native_binding_preserves_task4_surface_and_private_audit() -> None:
    assert hasattr(_coverage_2, "ReachableDomainConstructionError")
    assert hasattr(_coverage_2, "ReachableDomainCertificate2")
    assert hasattr(_coverage_2, "ReachableDomain2")
    assert not hasattr(
        _coverage_2.ReachableDomain2,
        "build_audit_for_native_gate",
    )


def test_python_certificate_rejects_unknown_strategy() -> None:
    certificate = replace(
        _domain().certificate,
        strategy_version=b"exact-reachable-arrangement-v999",
    )

    _assert_certificate_rejected(certificate, "strategy")


@pytest.mark.parametrize(
    ("field_index", "replacement"),
    (
        pytest.param(0, b"unknown", id="ring-role"),
        pytest.param(1, (1).to_bytes(8, "big"), id="ring-ordinal"),
        pytest.param(2, (99).to_bytes(8, "big"), id="feature-ordinal"),
        pytest.param(3, b"unknown", id="construction-role"),
        pytest.param(4, bytes(8), id="radius"),
    ),
)
def test_python_certificate_rejects_malformed_source_schema(
    field_index: int,
    replacement: bytes,
) -> None:
    certificate = _malformed_source_certificate(
        _domain().certificate,
        field_index,
        replacement,
    )

    _assert_certificate_rejected(certificate, "source")


def test_python_certificate_rejects_consistently_relabelled_invalid_source() -> None:
    certificate = _globally_relabelled_source_certificate(_domain().certificate)

    _assert_certificate_rejected(certificate, "source")


@pytest.mark.parametrize(
    "field_index",
    (
        pytest.param(0, id="source-endpoint"),
        pytest.param(1, id="target-endpoint"),
        pytest.param(2, id="source-set"),
        pytest.param(3, id="direction"),
    ),
)
def test_python_certificate_rejects_malformed_selected_cell_schema(
    field_index: int,
) -> None:
    certificate = _domain().certificate
    cell = certificate.selected_cell_records[0]
    cell_fields = _record_fields(cell, b"selected-cell-v1")
    cycle_fields = _record_fields(
        cell_fields[0],
        b"selected-boundary-cycle-v1",
    )
    element = _record_fields(cycle_fields[1], b"cycle-elements-v1")[0]
    element_fields = _record_fields(element, b"selected-cycle-element-v1")
    replacement = _tagged_record(
        b"invalid-selected-cycle-field-v1",
        _record_fields(element_fields[field_index], element_fields[field_index].split(b"\x00", 1)[0]),
    )
    malformed = _malformed_cell_certificate(
        certificate,
        field_index,
        replacement,
    )

    _assert_certificate_rejected(malformed, "cycle")


def test_python_certificate_rejects_malformed_inner_ccb_schema() -> None:
    certificate = _domain(RECTANGLE, (ISLAND,)).certificate
    cell = certificate.selected_cell_records[0]
    cell_fields = _record_fields(cell, b"selected-cell-v1")
    holes = list(_record_fields(cell_fields[1], b"selected-cell-holes-v1"))
    hole_fields = _record_fields(holes[0], b"selected-boundary-cycle-v1")
    holes[0] = _tagged_record(
        b"selected-boundary-cycle-v1",
        (b"outer", hole_fields[1]),
    )
    malformed_cell = _tagged_record(
        b"selected-cell-v1",
        (
            cell_fields[0],
            _tagged_record(b"selected-cell-holes-v1", tuple(holes)),
        ),
    )
    malformed = replace(certificate, selected_cell_records=(malformed_cell,))

    _assert_certificate_rejected(malformed, "hole")


def test_python_certificate_rejects_noncanonical_cycle_rotation() -> None:
    certificate = _domain().certificate
    cell = certificate.selected_cell_records[0]
    cell_fields = _record_fields(cell, b"selected-cell-v1")
    cycle_fields = _record_fields(
        cell_fields[0],
        b"selected-boundary-cycle-v1",
    )
    elements = _record_fields(cycle_fields[1], b"cycle-elements-v1")
    rotated_elements = elements[1:] + elements[:1]
    malformed_cell = _tagged_record(
        b"selected-cell-v1",
        (
            _tagged_record(
                b"selected-boundary-cycle-v1",
                (
                    cycle_fields[0],
                    _tagged_record(
                        b"cycle-elements-v1",
                        rotated_elements,
                    ),
                ),
            ),
            cell_fields[1],
        ),
    )
    malformed = replace(
        certificate,
        selected_cell_records=(malformed_cell,),
    )

    _assert_certificate_rejected(malformed, "canonical rotation")


def test_minimal_rotation_index_handles_repeated_prefixes() -> None:
    elements = (b"a", b"a", b"a", b"d", b"a", b"a", b"a", b"c")

    assert reachable_domain_module._minimal_rotation_index(elements) == 4


def test_python_certificate_rejects_unknown_component_cell() -> None:
    certificate = _domain().certificate
    unknown_cell = _tagged_record(
        b"selected-cell-id-v1",
        ((1).to_bytes(8, "big"),),
    )
    malformed = _malformed_component_certificate(
        certificate,
        cells=_tagged_record(
            b"selected-component-cells-v1",
            (unknown_cell,),
        ),
    )

    _assert_certificate_rejected(malformed, "cell")


def test_python_certificate_rejects_self_component_adjacency() -> None:
    certificate = _domain().certificate
    component_fields = _record_fields(
        certificate.component_records[0],
        b"selected-component-v1",
    )
    cell_id = _record_fields(
        component_fields[0],
        b"selected-component-cells-v1",
    )[0]
    malformed = _malformed_component_certificate(
        certificate,
        adjacencies=_tagged_record(
            b"selected-component-adjacencies-v1",
            (
                _tagged_record(
                    b"selected-cell-adjacency-v1",
                    (cell_id, cell_id),
                ),
            ),
        ),
    )

    _assert_certificate_rejected(malformed, "adjacency")


def test_python_certificate_rejects_malformed_component_boundary_cycle() -> None:
    certificate = _domain().certificate
    component_fields = _record_fields(
        certificate.component_records[0],
        b"selected-component-v1",
    )
    boundary_cycle = _record_fields(
        component_fields[2],
        b"selected-component-boundary-cycles-v1",
    )[0]
    cycle_fields = _record_fields(
        boundary_cycle,
        b"selected-boundary-cycle-v1",
    )
    malformed = _malformed_component_certificate(
        certificate,
        boundary_cycles=_tagged_record(
            b"selected-component-boundary-cycles-v1",
            (
                _tagged_record(
                    b"selected-boundary-cycle-v1",
                    (b"unknown", cycle_fields[1]),
                ),
            ),
        ),
    )

    _assert_certificate_rejected(malformed, "boundary")


def test_certificate_is_immutable_and_rejects_invalid_structural_records() -> None:
    certificate = _domain().certificate
    changed_source = replace(
        certificate,
        source_curve_records=certificate.source_curve_records[:-1],
    )
    changed_vertex = replace(
        certificate,
        arrangement_vertex_records=(certificate.arrangement_vertex_records[0] + b"x",),
    )
    changed_cell = replace(
        certificate,
        selected_cell_records=(certificate.selected_cell_records[0] + b"x",),
    )
    changed_component = replace(
        certificate,
        component_records=(certificate.component_records[0] + b"x",),
    )

    assert len(certificate.digest) == 32
    _assert_certificate_rejected(changed_source, "source")
    _assert_certificate_rejected(changed_vertex, "trailing bytes")
    _assert_certificate_rejected(changed_cell, "trailing bytes")
    _assert_certificate_rejected(changed_component, "trailing bytes")
    with pytest.raises(FrozenInstanceError):
        certificate.strategy_version = b"changed"  # type: ignore[misc]


def test_python_owner_rejects_native_certificate_relation_mutation() -> None:
    certificate = _domain().certificate
    invalid_selection = replace(certificate, exact_cell_selection=False)
    invalid_provenance = replace(
        certificate,
        complete_source_provenance=False,
    )

    with pytest.raises(InvalidReachableDomainCertificateError, match="cell selection"):
        ReachableDomainCertificate.validate(invalid_selection)
    with pytest.raises(InvalidReachableDomainCertificateError, match="source provenance"):
        ReachableDomainCertificate.validate(invalid_provenance)


def test_native_reachable_errors_are_distinct_and_audit_stays_private() -> None:
    assert InvalidReachableDomainInputError is _coverage_2.InvalidReachableDomainInputError
    assert ReachableArrangementTopologyError is _coverage_2.ReachableArrangementTopologyError
    assert PocketNotMachinableError is _coverage_2.PocketNotMachinableError
    assert ReachableMaterialContainmentError is _coverage_2.ReachableMaterialContainmentError
    assert (
        len(
            {
                InvalidReachableDomainInputError,
                ReachableArrangementTopologyError,
                PocketNotMachinableError,
                ReachableMaterialContainmentError,
            }
        )
        == 4
    )

    with pytest.raises(InvalidReachableDomainInputError):
        _coverage_2.ReachableDomain2(_array(RECTANGLE[:2]), [], 1.0)
    native = _coverage_2.ReachableDomain2(_array(RECTANGLE), [], 1.0)
    assert not hasattr(native, "build_audit_for_native_gate")


# --- Exact segment containment (CutterCentreDomain2.contains_segment) --------
#
# A segment [a, b] lies in the closed cutter-centre domain iff a does and the
# segment keeps squared distance >= r^2 from every ring edge. The reference
# below decides the same statement over Fractions of the injected doubles; it
# shares no code with the native kernel.

GATE_L_SHAPE = (
    (0.0, 0.0),
    (12.0, 0.0),
    (12.0, 6.0),
    (6.0, 6.0),
    (6.0, 10.0),
    (0.0, 10.0),
)
# A chord of length 0.985 r whose ends lie 1.0001 r from the reflex vertex
# (6, 6) and whose middle passes 0.1296 r inside that vertex's r-disk.
REFLEX_CHORD = (
    (5.036273183746056, 5.732734900084135),
    (5.732734900084135, 5.036273183746056),
)
# Interior stations at k / SOUNDNESS_STATIONS of a segment, exact rationals.
SOUNDNESS_STATIONS = 34

_FPt = tuple[Fraction, Fraction]


def _f_sub(p: _FPt, q: _FPt) -> _FPt:
    return (p[0] - q[0], p[1] - q[1])


def _f_dot(u: _FPt, v: _FPt) -> Fraction:
    return u[0] * v[0] + u[1] * v[1]


def _f_cross(u: _FPt, v: _FPt) -> Fraction:
    return u[0] * v[1] - u[1] * v[0]


def _f_sign(x: Fraction) -> int:
    return (x > 0) - (x < 0)


def _f_point_segment_sqdist(p: _FPt, a: _FPt, b: _FPt) -> Fraction:
    d = _f_sub(b, a)
    w = _f_sub(p, a)
    dd = _f_dot(d, d)
    if dd == 0:
        return _f_dot(w, w)
    t = _f_dot(w, d)
    if t <= 0:
        return _f_dot(w, w)
    if t >= dd:
        v = _f_sub(p, b)
        return _f_dot(v, v)
    c = _f_cross(d, w)
    return c * c / dd


def _f_on_box(p: _FPt, a: _FPt, b: _FPt) -> bool:
    return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def _f_segments_intersect(p0: _FPt, p1: _FPt, q0: _FPt, q1: _FPt) -> bool:
    o1 = _f_sign(_f_cross(_f_sub(p1, p0), _f_sub(q0, p0)))
    o2 = _f_sign(_f_cross(_f_sub(p1, p0), _f_sub(q1, p0)))
    o3 = _f_sign(_f_cross(_f_sub(q1, q0), _f_sub(p0, q0)))
    o4 = _f_sign(_f_cross(_f_sub(q1, q0), _f_sub(p1, q0)))
    if o1 * o2 < 0 and o3 * o4 < 0:
        return True
    return (
        (o1 == 0 and _f_on_box(q0, p0, p1))
        or (o2 == 0 and _f_on_box(q1, p0, p1))
        or (o3 == 0 and _f_on_box(p0, q0, q1))
        or (o4 == 0 and _f_on_box(p1, q0, q1))
    )


def _f_segment_segment_sqdist(p0: _FPt, p1: _FPt, q0: _FPt, q1: _FPt) -> Fraction:
    if _f_segments_intersect(p0, p1, q0, q1):
        return Fraction(0)
    return min(
        _f_point_segment_sqdist(p0, q0, q1),
        _f_point_segment_sqdist(p1, q0, q1),
        _f_point_segment_sqdist(q0, p0, p1),
        _f_point_segment_sqdist(q1, p0, p1),
    )


def _f_edges(rings: tuple[tuple[tuple[float, float], ...], ...]) -> list[tuple[_FPt, _FPt]]:
    edges: list[tuple[_FPt, _FPt]] = []
    for ring in rings:
        points = [(Fraction(x), Fraction(y)) for x, y in ring]
        edges += [(points[i], points[(i + 1) % len(points)]) for i in range(len(points))]
    return edges


def _f_clears(
    edges: list[tuple[_FPt, _FPt]],
    radius: float,
    a: tuple[float, float],
    b: tuple[float, float],
) -> bool:
    fa = (Fraction(a[0]), Fraction(a[1]))
    fb = (Fraction(b[0]), Fraction(b[1]))
    squared_radius = Fraction(radius) ** 2
    return all(_f_segment_segment_sqdist(fa, fb, e0, e1) >= squared_radius for e0, e1 in edges)


def test_cutter_centre_domain_segment_refuses_a_reflex_corner_chord() -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(GATE_L_SHAPE), [], 1.0)
    a, b = REFLEX_CHORD

    assert domain.contains(*a)
    assert domain.contains(*b)
    assert not domain.contains_segment(*a, *b)
    assert not domain.contains_segment(*b, *a)


def test_cutter_centre_domain_segment_keeps_exact_closed_clearance() -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(RECTANGLE), [], 1.0)
    x0 = 0.1 + 0.2 + 1.0
    x1 = math.pi + 2.0
    below = float(np.nextafter(1.0, 0.0))

    # Clearance exactly r from the bottom wall: the domain is closed.
    assert domain.contains_segment(x0, 1.0, x1, 1.0)
    assert not domain.contains_segment(x0, below, x1, below)
    # Legal start, end one ulp past the offset wall: the clearance decides.
    assert domain.contains(x0, 2.0)
    assert not domain.contains(x1, below)
    assert not domain.contains_segment(x0, 2.0, x1, below)
    # Interior point to a point past the wall: only the in-domain end passes.
    inside, outside = (5.0, 4.0), (5.0, 0.5)
    assert domain.contains(*inside)
    assert not domain.contains(*outside)
    assert not domain.contains_segment(*inside, *outside)
    assert not domain.contains_segment(*outside, *inside)
    # A degenerate segment is its point.
    assert domain.contains_segment(*inside, *inside)
    assert not domain.contains_segment(*outside, *outside)


@pytest.mark.parametrize("coordinates", [(math.nan, 0.0), (0.0, math.inf), (-math.inf, 1.0)])
def test_cutter_centre_domain_segment_rejects_nonfinite_ends(coordinates: tuple[float, float]) -> None:
    domain = _coverage_2.CutterCentreDomain2.build(_array(RECTANGLE), [], 1.0)

    with pytest.raises(InvalidReachableDomainInputError):
        domain.contains_segment(5.0, 4.0, *coordinates)
    with pytest.raises(InvalidReachableDomainInputError):
        domain.contains_segment(*coordinates, 5.0, 4.0)


_SEGMENT_FIXTURES = {
    "holed": (RECTANGLE, (ISLAND,)),
    "l_shape": (GATE_L_SHAPE, ()),
}
_SEGMENT_DOMAINS = {
    name: _coverage_2.CutterCentreDomain2.build(_array(outer), [_array(hole) for hole in holes], 1.0)
    for name, (outer, holes) in _SEGMENT_FIXTURES.items()
}
_SEGMENT_EDGES = {name: _f_edges((outer, *holes)) for name, (outer, holes) in _SEGMENT_FIXTURES.items()}
# Bounding box of both fixtures, padded by one tool radius (mm).
_BOX_COORDINATE = st.floats(min_value=-1.0, max_value=13.0, allow_nan=False, allow_infinity=False)
# Offsets for short segments, up to two tool radii (mm).
_SHORT_OFFSET = st.floats(min_value=-2.0, max_value=2.0, allow_nan=False, allow_infinity=False)


@st.composite
def _segment_query(draw: st.DrawFn) -> tuple[str, tuple[float, float], tuple[float, float]]:
    name = draw(st.sampled_from(sorted(_SEGMENT_FIXTURES)))
    a = (draw(_BOX_COORDINATE), draw(_BOX_COORDINATE))
    if draw(st.booleans()):
        b = (a[0] + draw(_SHORT_OFFSET), a[1] + draw(_SHORT_OFFSET))
    else:
        b = (draw(_BOX_COORDINATE), draw(_BOX_COORDINATE))
    return name, a, b


@settings(max_examples=300, deadline=None)
@given(_segment_query())
@example(("l_shape", REFLEX_CHORD[0], REFLEX_CHORD[1]))
@example(("holed", (0.1 + 0.2 + 1.0, 1.0), (math.pi + 2.0, 1.0)))
@example(("holed", (1.5, 2.0), (8.5, 2.0)))
def test_cutter_centre_domain_segment_matches_rational_oracle(
    query: tuple[str, tuple[float, float], tuple[float, float]],
) -> None:
    name, a, b = query
    domain = _SEGMENT_DOMAINS[name]
    edges = _SEGMENT_EDGES[name]

    native = domain.contains_segment(*a, *b)

    assert native == (domain.contains(*a) and _f_clears(edges, 1.0, a, b))
    if native:
        assert domain.contains(*a)
        assert domain.contains(*b)
        fa = (Fraction(a[0]), Fraction(a[1]))
        fb = (Fraction(b[0]), Fraction(b[1]))
        for k in range(1, SOUNDNESS_STATIONS):
            t = Fraction(k, SOUNDNESS_STATIONS)
            station = (fa[0] + t * (fb[0] - fa[0]), fa[1] + t * (fb[1] - fa[1]))
            assert all(_f_point_segment_sqdist(station, e0, e1) >= 1 for e0, e1 in edges)
