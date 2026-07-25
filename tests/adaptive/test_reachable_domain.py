from dataclasses import FrozenInstanceError
from dataclasses import replace

import numpy as np
import pytest

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
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


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
