import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Self

import numpy as np

from compas_cgal import _coverage_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import canonical_polygon_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidReachableDomainCertificateError
from compas_cgal.adaptive.units import ToolRadius


def _native_ring(ring: CanonicalRingV1) -> np.ndarray:
    return np.asarray(
        tuple((point.x, point.y, 0.0) for point in ring.vertices),
        dtype=np.float64,
    )


def _bytes_tuple(values: Sequence[bytes], name: str) -> tuple[bytes, ...]:
    try:
        frozen = tuple(values)
    except TypeError:
        raise InvalidReachableDomainCertificateError(f"{name} must be a finite native byte sequence.") from None
    if any(type(value) is not bytes or not value for value in frozen):
        raise InvalidReachableDomainCertificateError(f"{name} must contain nonempty exact bytes.")
    return frozen


def _record_fields(record: bytes, tag: bytes, name: str) -> tuple[bytes, ...]:
    prefix = tag + b"\x00"
    if not record.startswith(prefix):
        raise InvalidReachableDomainCertificateError(f"{name} has the wrong exact record tag.")
    offset = len(prefix)
    if len(record) - offset < 8:
        raise InvalidReachableDomainCertificateError(f"{name} is missing its exact field count.")
    count = int.from_bytes(record[offset : offset + 8], "big")
    offset += 8
    fields: list[bytes] = []
    for _ in range(count):
        if len(record) - offset < 8:
            raise InvalidReachableDomainCertificateError(f"{name} has a truncated exact field length.")
        size = int.from_bytes(record[offset : offset + 8], "big")
        offset += 8
        if size > len(record) - offset:
            raise InvalidReachableDomainCertificateError(f"{name} has a truncated exact field.")
        fields.append(record[offset : offset + size])
        offset += size
    if offset != len(record):
        raise InvalidReachableDomainCertificateError(f"{name} has trailing bytes.")
    return tuple(fields)


def _validate_arrangement_vertex_record(
    record: bytes,
    source_curve_records: tuple[bytes, ...],
) -> None:
    fields = _record_fields(
        record,
        b"arrangement-vertex-v1",
        "arrangement vertex record",
    )
    if len(fields) != 2:
        raise InvalidReachableDomainCertificateError("arrangement vertex record must bind incidence and structural ordinal.")
    incidence_fields = _record_fields(
        fields[0],
        b"arrangement-incidence-multiset-v1",
        "arrangement vertex incidence record",
    )
    for source_piece in incidence_fields:
        source_fields = _record_fields(
            source_piece,
            b"source-piece-v1",
            "arrangement vertex source-piece record",
        )
        if len(source_fields) != 2 or source_fields[0] not in source_curve_records or len(source_fields[1]) != 8:
            raise InvalidReachableDomainCertificateError("arrangement vertex source-piece record has invalid identity or ordinal.")
    if len(fields[1]) != 8:
        raise InvalidReachableDomainCertificateError("arrangement vertex structural ordinal must be one exact unsigned integer.")


def _region_digest(role: bytes, certificate_digest: bytes) -> bytes:
    recipe = encode_tagged_union(
        role,
        encode_component_map({b"reachable-domain-certificate-digest": encode_bytes(certificate_digest)}),
    )
    return hashlib.sha256(recipe).digest()


@dataclass(frozen=True)
class ReachableDomainCertificate:
    design_boundary: CanonicalRingV1
    holes: tuple[CanonicalRingV1, ...]
    tool_radius: ToolRadius
    strategy_version: bytes
    source_curve_records: tuple[bytes, ...]
    arrangement_vertex_records: tuple[bytes, ...]
    selected_cell_records: tuple[bytes, ...]
    component_records: tuple[bytes, ...]
    exact_cell_selection: bool
    complete_source_provenance: bool
    reachable_subset_of_design: bool

    @classmethod
    def validate(cls, certificate: "ReachableDomainCertificate") -> None:
        if type(certificate) is not cls:
            raise InvalidReachableDomainCertificateError("reachable-domain certificate must use the exact owned type.")
        if type(certificate.design_boundary) is not CanonicalRingV1:
            raise InvalidReachableDomainCertificateError("reachable-domain certificate requires one canonical outer ring.")
        if not certificate.design_boundary.is_outer:
            raise InvalidReachableDomainCertificateError("reachable-domain certificate design boundary must be outer.")
        if type(certificate.holes) is not tuple or any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in certificate.holes):
            raise InvalidReachableDomainCertificateError("reachable-domain certificate holes must be canonical hole rings.")
        hole_bytes = tuple(hole.canonical_bytes for hole in certificate.holes)
        if hole_bytes != tuple(sorted(hole_bytes)) or len(hole_bytes) != len(set(hole_bytes)):
            raise InvalidReachableDomainCertificateError("reachable-domain certificate holes must be sorted and unique.")
        if type(certificate.tool_radius) is not ToolRadius:
            raise InvalidReachableDomainCertificateError("reachable-domain certificate requires an exact typed tool radius.")
        if type(certificate.strategy_version) is not bytes or not certificate.strategy_version:
            raise InvalidReachableDomainCertificateError("reachable-domain strategy version must be nonempty bytes.")
        for name, records in (
            ("source curve records", certificate.source_curve_records),
            ("arrangement vertex records", certificate.arrangement_vertex_records),
            ("selected cell records", certificate.selected_cell_records),
            ("component records", certificate.component_records),
        ):
            if type(records) is not tuple or any(type(record) is not bytes or not record for record in records):
                raise InvalidReachableDomainCertificateError(f"{name} must be an exact nonempty byte tuple.")
            if records != tuple(sorted(records)) or len(records) != len(set(records)):
                raise InvalidReachableDomainCertificateError(f"{name} must be sorted and unique.")
        expected_sources = 3 * (certificate.design_boundary.vertex_count + sum(hole.vertex_count for hole in certificate.holes))
        if len(certificate.source_curve_records) != expected_sources:
            raise InvalidReachableDomainCertificateError("source curve reconstruction does not bind every boundary feature.")
        if not certificate.arrangement_vertex_records:
            raise InvalidReachableDomainCertificateError("reachable-domain reconstruction binds no arrangement vertices.")
        for record in certificate.arrangement_vertex_records:
            _validate_arrangement_vertex_record(
                record,
                certificate.source_curve_records,
            )
        if not certificate.selected_cell_records:
            raise InvalidReachableDomainCertificateError("reachable-domain reconstruction selected no exact cells.")
        if len(certificate.component_records) != 1:
            raise InvalidReachableDomainCertificateError("reachable-domain reconstruction requires one exact component.")
        for name, relation in (
            ("cell selection", certificate.exact_cell_selection),
            ("source provenance", certificate.complete_source_provenance),
            ("reachable subset", certificate.reachable_subset_of_design),
        ):
            if type(relation) is not bool or not relation:
                raise InvalidReachableDomainCertificateError(f"reachable-domain certificate {name} relation is not exact.")

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"reachable-domain-certificate-v2",
            encode_component_map(
                {
                    b"complete-source-provenance": encode_boolean(self.complete_source_provenance),
                    b"component-records": encode_sequence(self.component_records),
                    b"design": canonical_polygon_bytes(
                        self.design_boundary.vertices,
                        tuple(hole.vertices for hole in self.holes),
                    ),
                    b"exact-cell-selection": encode_boolean(self.exact_cell_selection),
                    b"arrangement-vertex-records": encode_sequence(self.arrangement_vertex_records),
                    b"reachable-subset-of-design": encode_boolean(self.reachable_subset_of_design),
                    b"selected-cell-records": encode_sequence(self.selected_cell_records),
                    b"source-curve-records": encode_sequence(self.source_curve_records),
                    b"strategy-version": encode_bytes(self.strategy_version),
                    b"tool-radius": encode_binary64(float(self.tool_radius.value)),
                }
            ),
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()

    @property
    def center_domain_digest(self) -> bytes:
        return _region_digest(b"exact-region-center-domain-v1", self.digest)

    @property
    def reachable_material_digest(self) -> bytes:
        return _region_digest(b"exact-region-reachable-material-v1", self.digest)

    @property
    def unreachable_residual_digest(self) -> bytes:
        return _region_digest(b"exact-region-unreachable-residual-v1", self.digest)


class ReachableDomain:
    def __init__(
        self,
        native: _coverage_2.ReachableDomain2,
        certificate: ReachableDomainCertificate,
    ) -> None:
        self._native = native
        self._certificate = certificate

    @classmethod
    def build(
        cls,
        *,
        design_boundary: CanonicalRingV1,
        holes: Sequence[CanonicalRingV1],
        tool_radius: ToolRadius,
    ) -> Self:
        if type(design_boundary) is not CanonicalRingV1 or not design_boundary.is_outer:
            raise InvalidReachableDomainCertificateError("reachable domain requires one canonical outer design boundary.")
        try:
            canonical_holes = tuple(holes)
        except TypeError:
            raise InvalidReachableDomainCertificateError("reachable domain holes must be a finite canonical sequence.") from None
        if any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in canonical_holes):
            raise InvalidReachableDomainCertificateError("reachable domain holes must be canonical hole rings.")
        canonical_holes = tuple(sorted(canonical_holes, key=lambda hole: hole.canonical_bytes))
        if len(canonical_holes) != len({hole.canonical_bytes for hole in canonical_holes}):
            raise InvalidReachableDomainCertificateError("reachable domain holes must be canonically unique.")
        if type(tool_radius) is not ToolRadius:
            raise InvalidReachableDomainCertificateError("reachable domain requires an exact typed tool radius.")

        boundary_array = _native_ring(design_boundary)
        hole_arrays = [_native_ring(hole) for hole in canonical_holes]
        native = _coverage_2.ReachableDomain2(
            boundary_array,
            hole_arrays,
            tool_radius.value,
        )
        native_certificate = native.certificate()
        if not native_certificate.matches_exact_inputs(
            boundary_array,
            hole_arrays,
            tool_radius.value,
        ):
            raise InvalidReachableDomainCertificateError("native reachable-domain certificate does not match exact inputs.")
        certificate = ReachableDomainCertificate(
            design_boundary=design_boundary,
            holes=canonical_holes,
            tool_radius=tool_radius,
            strategy_version=bytes(native_certificate.strategy_version),
            source_curve_records=_bytes_tuple(
                native_certificate.source_curve_records,
                "source curve records",
            ),
            arrangement_vertex_records=_bytes_tuple(
                native_certificate.arrangement_vertex_records,
                "arrangement vertex records",
            ),
            selected_cell_records=_bytes_tuple(
                native_certificate.selected_cell_records,
                "selected cell records",
            ),
            component_records=_bytes_tuple(
                native_certificate.component_records,
                "component records",
            ),
            exact_cell_selection=native_certificate.exact_cell_selection,
            complete_source_provenance=native_certificate.complete_source_provenance,
            reachable_subset_of_design=native_certificate.reachable_subset_of_design,
        )
        ReachableDomainCertificate.validate(certificate)
        return cls(native, certificate)

    @property
    def certificate(self) -> ReachableDomainCertificate:
        return self._certificate

    @property
    def center_domain(self) -> _coverage_2.ExactRegion2:
        return self._native.center_domain()

    @property
    def reachable_material(self) -> _coverage_2.ExactRegion2:
        return self._native.reachable_material()

    @property
    def unreachable_residual(self) -> _coverage_2.ExactRegion2:
        return self._native.unreachable_residual()
