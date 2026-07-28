import hashlib
import struct
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final
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

_U64_BYTES: Final[int] = 8
_BINARY64: Final[struct.Struct] = struct.Struct(">d")
_REACHABLE_STRATEGY_VERSION: Final[bytes] = b"exact-reachable-arrangement-v2"
_SOURCE_PIECE_COUNTS: Final[dict[bytes, int]] = {
    b"offset-minus": 1,
    b"offset-plus": 1,
    b"vertex-circle": 2,
}
_DIRECTION_ROLES: Final[frozenset[bytes]] = frozenset(
    {
        b"linear-right",
        b"linear-left",
        b"circular-counterclockwise",
        b"circular-clockwise",
    }
)


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
    if len(record) - offset < _U64_BYTES:
        raise InvalidReachableDomainCertificateError(f"{name} is missing its exact field count.")
    count = int.from_bytes(record[offset : offset + _U64_BYTES], "big")
    offset += _U64_BYTES
    if count > (len(record) - offset) // _U64_BYTES:
        raise InvalidReachableDomainCertificateError(f"{name} declares impossible exact fields.")
    fields: list[bytes] = []
    for _ in range(count):
        if len(record) - offset < _U64_BYTES:
            raise InvalidReachableDomainCertificateError(f"{name} has a truncated exact field length.")
        size = int.from_bytes(record[offset : offset + _U64_BYTES], "big")
        offset += _U64_BYTES
        if size > len(record) - offset:
            raise InvalidReachableDomainCertificateError(f"{name} has a truncated exact field.")
        fields.append(record[offset : offset + size])
        offset += size
    if offset != len(record):
        raise InvalidReachableDomainCertificateError(f"{name} has trailing bytes.")
    return tuple(fields)


def _decode_u64(field: bytes, name: str) -> int:
    if len(field) != _U64_BYTES:
        raise InvalidReachableDomainCertificateError(f"{name} must be one exact unsigned integer.")
    return int.from_bytes(field, "big")


def _require_field_count(
    fields: tuple[bytes, ...],
    expected: int,
    name: str,
) -> None:
    if len(fields) != expected:
        raise InvalidReachableDomainCertificateError(f"{name} has the wrong exact field count.")


def _validate_source_curve_records(
    certificate: "ReachableDomainCertificate",
) -> dict[bytes, int]:
    ring_sizes: dict[tuple[bytes, int], int] = {
        (b"outer", 0): certificate.design_boundary.vertex_count,
    }
    for ordinal, hole in enumerate(certificate.holes):
        ring_sizes[(b"hole", ordinal)] = hole.vertex_count

    expected = {
        (role, ring_ordinal, feature_ordinal, construction_role)
        for (role, ring_ordinal), vertex_count in ring_sizes.items()
        for feature_ordinal in range(vertex_count)
        for construction_role in _SOURCE_PIECE_COUNTS
    }
    observed: set[tuple[bytes, int, int, bytes]] = set()
    contracts: dict[bytes, int] = {}
    expected_radius = _BINARY64.pack(float(certificate.tool_radius.value))
    for record in certificate.source_curve_records:
        fields = _record_fields(record, b"source-curve-v2", "source curve record")
        _require_field_count(fields, 5, "source curve record")
        ring_ordinal = _decode_u64(fields[1], "source ring ordinal")
        feature_ordinal = _decode_u64(fields[2], "source feature ordinal")
        ring = (fields[0], ring_ordinal)
        if ring not in ring_sizes:
            raise InvalidReachableDomainCertificateError("source curve record references an unknown canonical ring.")
        if feature_ordinal >= ring_sizes[ring]:
            raise InvalidReachableDomainCertificateError("source curve record references an unknown ring feature.")
        piece_count = _SOURCE_PIECE_COUNTS.get(fields[3])
        if piece_count is None:
            raise InvalidReachableDomainCertificateError("source curve record has an unknown construction role.")
        if fields[4] != expected_radius:
            raise InvalidReachableDomainCertificateError("source curve record radius contradicts the exact tool radius.")
        observed.add((fields[0], ring_ordinal, feature_ordinal, fields[3]))
        contracts[record] = piece_count
    if observed != expected or len(contracts) != len(expected):
        raise InvalidReachableDomainCertificateError("source curve records do not cover the exact construction universe.")
    return contracts


def _validate_source_piece(
    record: bytes,
    source_contracts: dict[bytes, int],
    name: str,
) -> bytes:
    fields = _record_fields(record, b"source-piece-v1", name)
    _require_field_count(fields, 2, name)
    piece_count = source_contracts.get(fields[0])
    if piece_count is None:
        raise InvalidReachableDomainCertificateError(f"{name} references an unknown source curve.")
    if _decode_u64(fields[1], f"{name} ordinal") >= piece_count:
        raise InvalidReachableDomainCertificateError(f"{name} has an invalid source-piece ordinal.")
    return fields[0]


def _validate_arrangement_vertex_record(
    record: bytes,
    source_contracts: dict[bytes, int],
    *,
    require_source_incidence: bool,
) -> tuple[bytes, int, set[bytes]]:
    fields = _record_fields(
        record,
        b"arrangement-vertex-v1",
        "arrangement vertex record",
    )
    _require_field_count(fields, 2, "arrangement vertex record")
    incidence_fields = _record_fields(
        fields[0],
        b"arrangement-incidence-multiset-v1",
        "arrangement vertex incidence record",
    )
    if require_source_incidence and not incidence_fields:
        raise InvalidReachableDomainCertificateError("selected cycle vertex has no exact source incidence.")
    if incidence_fields != tuple(sorted(incidence_fields)):
        raise InvalidReachableDomainCertificateError("arrangement vertex incidences are not canonical.")
    referenced_sources: set[bytes] = set()
    for source_piece in incidence_fields:
        referenced_sources.add(
            _validate_source_piece(
                source_piece,
                source_contracts,
                "arrangement vertex source piece",
            )
        )
    return (
        fields[0],
        _decode_u64(
            fields[1],
            "arrangement vertex structural ordinal",
        ),
        referenced_sources,
    )


def _validate_arrangement_vertex_records(
    records: tuple[bytes, ...],
    source_contracts: dict[bytes, int],
) -> set[bytes]:
    ordinals_by_incidence: dict[bytes, list[int]] = {}
    referenced_sources: set[bytes] = set()
    for record in records:
        incidence, ordinal, sources = _validate_arrangement_vertex_record(
            record,
            source_contracts,
            require_source_incidence=False,
        )
        ordinals_by_incidence.setdefault(incidence, []).append(ordinal)
        referenced_sources.update(sources)
    for ordinals in ordinals_by_incidence.values():
        if sorted(ordinals) != list(range(len(ordinals))):
            raise InvalidReachableDomainCertificateError("arrangement vertex structural ordinals are not complete and canonical.")
    if referenced_sources != set(source_contracts):
        raise InvalidReachableDomainCertificateError("arrangement vertex records do not bind every exact source curve.")
    return set(records)


def _minimal_rotation_index(values: tuple[bytes, ...]) -> int:
    size = len(values)
    left = 0
    right = 1
    offset = 0
    while left < size and right < size and offset < size:
        left_value = values[(left + offset) % size]
        right_value = values[(right + offset) % size]
        if left_value == right_value:
            offset += 1
            continue
        if left_value > right_value:
            left += offset + 1
            if left == right:
                left += 1
        else:
            right += offset + 1
            if left == right:
                right += 1
        offset = 0
    return min(left, right)


def _validate_boundary_cycle(
    record: bytes,
    *,
    expected_role: bytes | None,
    vertices: set[bytes],
    source_contracts: dict[bytes, int],
    name: str,
) -> bytes:
    fields = _record_fields(record, b"selected-boundary-cycle-v1", name)
    _require_field_count(fields, 2, name)
    role = fields[0]
    if role not in (b"outer", b"hole") or (expected_role is not None and role != expected_role):
        raise InvalidReachableDomainCertificateError(f"{name} has an invalid boundary role.")
    elements = _record_fields(fields[1], b"cycle-elements-v1", f"{name} elements")
    if not elements:
        raise InvalidReachableDomainCertificateError(f"{name} must contain exact cycle elements.")
    if _minimal_rotation_index(elements) != 0:
        raise InvalidReachableDomainCertificateError(f"{name} does not use canonical rotation.")

    endpoints: list[tuple[bytes, bytes]] = []
    for element in elements:
        element_fields = _record_fields(
            element,
            b"selected-cycle-element-v1",
            f"{name} cycle element",
        )
        _require_field_count(element_fields, 4, f"{name} cycle element")
        source_endpoint = _record_fields(
            element_fields[0],
            b"source-endpoint-v1",
            f"{name} source endpoint",
        )
        target_endpoint = _record_fields(
            element_fields[1],
            b"target-endpoint-v1",
            f"{name} target endpoint",
        )
        _require_field_count(source_endpoint, 1, f"{name} source endpoint")
        _require_field_count(target_endpoint, 1, f"{name} target endpoint")
        source_vertex = source_endpoint[0]
        target_vertex = target_endpoint[0]
        if source_vertex == target_vertex or source_vertex not in vertices or target_vertex not in vertices:
            raise InvalidReachableDomainCertificateError(f"{name} cycle element references invalid vertices.")
        _validate_arrangement_vertex_record(
            source_vertex,
            source_contracts,
            require_source_incidence=True,
        )
        _validate_arrangement_vertex_record(
            target_vertex,
            source_contracts,
            require_source_incidence=True,
        )
        source_pieces = _record_fields(
            element_fields[2],
            b"curve-source-set-v1",
            f"{name} curve source set",
        )
        if not source_pieces or source_pieces != tuple(sorted(source_pieces)) or len(source_pieces) != len(set(source_pieces)):
            raise InvalidReachableDomainCertificateError(f"{name} curve source set is empty or noncanonical.")
        for source_piece in source_pieces:
            _validate_source_piece(
                source_piece,
                source_contracts,
                f"{name} curve source piece",
            )
        direction = _record_fields(
            element_fields[3],
            b"curve-direction-v1",
            f"{name} curve direction",
        )
        if len(direction) != 1 or direction[0] not in _DIRECTION_ROLES:
            raise InvalidReachableDomainCertificateError(f"{name} has an invalid exact curve direction.")
        endpoints.append((source_vertex, target_vertex))
    if any(target != endpoints[(index + 1) % len(endpoints)][0] for index, (_, target) in enumerate(endpoints)):
        raise InvalidReachableDomainCertificateError(f"{name} does not form a closed exact cycle.")
    return role


def _validate_selected_cell_records(
    records: tuple[bytes, ...],
    vertices: set[bytes],
    source_contracts: dict[bytes, int],
) -> None:
    for record in records:
        fields = _record_fields(record, b"selected-cell-v1", "selected cell record")
        _require_field_count(fields, 2, "selected cell record")
        _validate_boundary_cycle(
            fields[0],
            expected_role=b"outer",
            vertices=vertices,
            source_contracts=source_contracts,
            name="selected cell outer cycle",
        )
        holes = _record_fields(
            fields[1],
            b"selected-cell-holes-v1",
            "selected cell hole collection",
        )
        if holes != tuple(sorted(holes)) or len(holes) != len(set(holes)):
            raise InvalidReachableDomainCertificateError("selected cell hole cycles are not canonical.")
        for hole in holes:
            _validate_boundary_cycle(
                hole,
                expected_role=b"hole",
                vertices=vertices,
                source_contracts=source_contracts,
                name="selected cell hole cycle",
            )


def _validate_component_record(
    record: bytes,
    selected_cell_records: tuple[bytes, ...],
    vertices: set[bytes],
    source_contracts: dict[bytes, int],
) -> None:
    fields = _record_fields(record, b"selected-component-v1", "selected component record")
    _require_field_count(fields, 3, "selected component record")
    cells = _record_fields(
        fields[0],
        b"selected-component-cells-v1",
        "selected component cell collection",
    )
    if cells != tuple(sorted(cells)) or len(cells) != len(set(cells)) or len(cells) != len(selected_cell_records):
        raise InvalidReachableDomainCertificateError("selected component cell membership is incomplete or noncanonical.")
    ordinals: list[int] = []
    for cell in cells:
        cell_fields = _record_fields(cell, b"selected-cell-id-v1", "selected component cell identity")
        _require_field_count(cell_fields, 1, "selected component cell identity")
        ordinals.append(_decode_u64(cell_fields[0], "selected component cell ordinal"))
    if ordinals != list(range(len(selected_cell_records))):
        raise InvalidReachableDomainCertificateError("selected component cell identities do not bind selected records.")

    adjacencies = _record_fields(
        fields[1],
        b"selected-component-adjacencies-v1",
        "selected component adjacency collection",
    )
    if adjacencies != tuple(sorted(adjacencies)) or len(adjacencies) != len(set(adjacencies)):
        raise InvalidReachableDomainCertificateError("selected component adjacency records are not canonical.")
    cell_members = set(cells)
    for adjacency in adjacencies:
        pair = _record_fields(
            adjacency,
            b"selected-cell-adjacency-v1",
            "selected component adjacency",
        )
        if len(pair) != 2 or pair[0] >= pair[1] or pair[0] not in cell_members or pair[1] not in cell_members:
            raise InvalidReachableDomainCertificateError("selected component adjacency is self-referential or outside cell membership.")

    boundary_cycles = _record_fields(
        fields[2],
        b"selected-component-boundary-cycles-v1",
        "selected component boundary cycle collection",
    )
    if not boundary_cycles or boundary_cycles != tuple(sorted(boundary_cycles)) or len(boundary_cycles) != len(set(boundary_cycles)):
        raise InvalidReachableDomainCertificateError("selected component boundary cycles are empty or noncanonical.")
    roles = [
        _validate_boundary_cycle(
            cycle,
            expected_role=None,
            vertices=vertices,
            source_contracts=source_contracts,
            name="selected component boundary cycle",
        )
        for cycle in boundary_cycles
    ]
    if roles.count(b"outer") != 1:
        raise InvalidReachableDomainCertificateError("selected component boundary cycles require exactly one outer cycle.")


def _validate_structural_certificate(
    certificate: "ReachableDomainCertificate",
) -> None:
    source_contracts = _validate_source_curve_records(certificate)
    vertices = _validate_arrangement_vertex_records(
        certificate.arrangement_vertex_records,
        source_contracts,
    )
    _validate_selected_cell_records(
        certificate.selected_cell_records,
        vertices,
        source_contracts,
    )
    _validate_component_record(
        certificate.component_records[0],
        certificate.selected_cell_records,
        vertices,
        source_contracts,
    )


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
        if type(certificate.strategy_version) is not bytes or certificate.strategy_version != _REACHABLE_STRATEGY_VERSION:
            raise InvalidReachableDomainCertificateError("reachable-domain strategy version is unknown.")
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
        if not certificate.arrangement_vertex_records:
            raise InvalidReachableDomainCertificateError("reachable-domain reconstruction binds no arrangement vertices.")
        if not certificate.selected_cell_records:
            raise InvalidReachableDomainCertificateError("reachable-domain reconstruction selected no exact cells.")
        if len(certificate.component_records) != 1:
            raise InvalidReachableDomainCertificateError("reachable-domain reconstruction requires one exact component.")
        _validate_structural_certificate(certificate)
        for name, relation in (
            ("cell selection", certificate.exact_cell_selection),
            ("source provenance", certificate.complete_source_provenance),
            ("reachable subset", certificate.reachable_subset_of_design),
        ):
            if type(relation) is not bool or not relation:
                raise InvalidReachableDomainCertificateError(f"reachable-domain certificate {name} relation is not exact.")

    @property
    def canonical_bytes(self) -> bytes:
        type(self).validate(self)
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
    def design_region(self) -> _coverage_2.ExactRegion2:
        """Return an immutable shared-storage view of exact design pocket `D`."""
        return self._native.design_region()

    @property
    def center_domain(self) -> _coverage_2.ExactRegion2:
        return self._native.center_domain()

    @property
    def reachable_material(self) -> _coverage_2.ExactRegion2:
        return self._native.reachable_material()

    @property
    def unreachable_residual(self) -> _coverage_2.ExactRegion2:
        return self._native.unreachable_residual()
