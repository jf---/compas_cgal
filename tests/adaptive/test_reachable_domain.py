from dataclasses import FrozenInstanceError
from dataclasses import replace

import numpy as np
import pytest

from compas_cgal import _coverage_2
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


def test_circle_vertical_trim_extrema_have_distinct_stable_structural_ordinals() -> None:
    variants = (
        ACUTE,
        ACUTE[1:] + ACUTE[:1],
        tuple(reversed(ACUTE)),
    )
    native = _coverage_2.ReachableDomain2(_array(ACUTE), [], 1.0)
    certificate = native.certificate()
    assert all(certificate.matches_exact_inputs(_array(variant), [], 1.0) for variant in variants)
    assert tuple(certificate.arrangement_vertex_records) == tuple(native.certificate().arrangement_vertex_records)

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


def test_certificate_is_immutable_and_binds_every_structural_record_and_recipe() -> None:
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

    assert (
        len(
            {
                certificate.digest,
                changed_source.digest,
                changed_vertex.digest,
                changed_cell.digest,
                changed_component.digest,
            }
        )
        == 5
    )
    with pytest.raises(InvalidReachableDomainCertificateError, match="trailing bytes"):
        ReachableDomainCertificate.validate(changed_vertex)
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
