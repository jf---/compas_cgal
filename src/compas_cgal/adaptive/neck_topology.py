"""Typed exact topology carried by classified MAT necks."""

import hashlib
from dataclasses import dataclass
from typing import NewType
from typing import Self
from typing import TypeAlias

from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidNeckEvidenceError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatNodeId

NeckSideId = NewType("NeckSideId", bytes)

STRICT_EDGE_NECK_TAG = b"mat-neck-strict-edge-v1"
CLEARANCE_ENDPOINT_NECK_TAG = b"mat-neck-clearance-endpoint-v1"
SHARED_VERTEX_NECK_TAG = b"mat-neck-shared-vertex-v1"
PLATEAU_NECK_TAG = b"mat-neck-plateau-v1"

_SHA256_BYTES = hashlib.sha256().digest_size


def _exact_bytes(value: object, name: str, *, digest: bool = False) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidNeckEvidenceError(f"{name} must be nonempty exact bytes.")
    if digest and len(value) != _SHA256_BYTES:
        raise InvalidNeckEvidenceError(f"{name} must be one SHA-256 digest.")
    return value


def _canonical_ids(
    values: object,
    name: str,
    *,
    allow_empty: bool = False,
) -> tuple[bytes, ...]:
    if type(values) is not tuple:
        raise InvalidNeckEvidenceError(f"{name} must be one immutable tuple.")
    identities = values
    if (not identities and not allow_empty) or any(type(identity) is not bytes or not identity for identity in identities):
        raise InvalidNeckEvidenceError(f"{name} must contain canonical exact identities.")
    if identities != tuple(sorted(set(identities))):
        raise InvalidNeckEvidenceError(f"{name} must be unique and canonical.")
    return identities


@dataclass(frozen=True)
class StrictEdgeNeckLocus:
    edge_id: MatEdgeId
    parameter_root_id: bytes

    def __post_init__(self) -> None:
        _exact_bytes(self.edge_id, "strict-edge neck locus edge")
        _exact_bytes(self.parameter_root_id, "strict-edge neck parameter root")

    @classmethod
    def build(
        cls,
        *,
        edge_id: MatEdgeId,
        parameter_root_id: bytes,
    ) -> Self:
        return cls(edge_id, parameter_root_id)

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            STRICT_EDGE_NECK_TAG,
            encode_component_map(
                {
                    b"edge-id": encode_bytes(bytes(self.edge_id)),
                    b"parameter-root-id": encode_bytes(self.parameter_root_id),
                }
            ),
        )


@dataclass(frozen=True)
class ClearanceEndpointNeckLocus:
    edge_id: MatEdgeId
    node_id: MatNodeId
    parameter_root_id: bytes

    def __post_init__(self) -> None:
        _exact_bytes(self.edge_id, "clearance-endpoint neck locus edge")
        _exact_bytes(self.node_id, "clearance-endpoint neck locus node")
        _exact_bytes(self.parameter_root_id, "clearance-endpoint neck parameter root")

    @classmethod
    def build(
        cls,
        *,
        edge_id: MatEdgeId,
        node_id: MatNodeId,
        parameter_root_id: bytes,
    ) -> Self:
        return cls(edge_id, node_id, parameter_root_id)

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            CLEARANCE_ENDPOINT_NECK_TAG,
            encode_component_map(
                {
                    b"edge-id": encode_bytes(bytes(self.edge_id)),
                    b"node-id": encode_bytes(bytes(self.node_id)),
                    b"parameter-root-id": encode_bytes(self.parameter_root_id),
                }
            ),
        )


@dataclass(frozen=True)
class SharedVertexNeckLocus:
    node_id: MatNodeId
    minimizing_edge_ids: tuple[MatEdgeId, ...]

    def __post_init__(self) -> None:
        _exact_bytes(self.node_id, "shared-vertex neck locus node")
        edges = _canonical_ids(self.minimizing_edge_ids, "shared-vertex minimizing edges")
        if len(edges) < 2:
            raise InvalidNeckEvidenceError("shared-vertex neck locus requires at least two minimizing edges.")

    @classmethod
    def build(
        cls,
        *,
        node_id: MatNodeId,
        minimizing_edge_ids: tuple[MatEdgeId, ...],
    ) -> Self:
        return cls(node_id, minimizing_edge_ids)

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            SHARED_VERTEX_NECK_TAG,
            encode_component_map(
                {
                    b"minimizing-edge-ids": encode_sequence(tuple(encode_bytes(bytes(edge_id)) for edge_id in self.minimizing_edge_ids)),
                    b"node-id": encode_bytes(bytes(self.node_id)),
                }
            ),
        )


@dataclass(frozen=True)
class PlateauNeckLocus:
    node_ids: tuple[MatNodeId, ...]
    edge_ids: tuple[MatEdgeId, ...]

    def __post_init__(self) -> None:
        _canonical_ids(self.node_ids, "plateau neck locus nodes")
        _canonical_ids(self.edge_ids, "plateau neck locus edges")

    @classmethod
    def build(
        cls,
        *,
        node_ids: tuple[MatNodeId, ...],
        edge_ids: tuple[MatEdgeId, ...],
    ) -> Self:
        return cls(node_ids, edge_ids)

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            PLATEAU_NECK_TAG,
            encode_component_map(
                {
                    b"edge-ids": encode_sequence(tuple(encode_bytes(bytes(edge_id)) for edge_id in self.edge_ids)),
                    b"node-ids": encode_sequence(tuple(encode_bytes(bytes(node_id)) for node_id in self.node_ids)),
                }
            ),
        )


NeckLocus: TypeAlias = StrictEdgeNeckLocus | ClearanceEndpointNeckLocus | SharedVertexNeckLocus | PlateauNeckLocus

NECK_LOCUS_TYPES = (
    StrictEdgeNeckLocus,
    ClearanceEndpointNeckLocus,
    SharedVertexNeckLocus,
    PlateauNeckLocus,
)


def build_neck_locus(
    *,
    tag: bytes,
    edge_ids: tuple[MatEdgeId, ...],
    node_ids: tuple[MatNodeId, ...],
    parameter_root_id: bytes | None,
) -> NeckLocus:
    """Build the closed typed neck location projected by the native owner.

    Args:
        tag: Versioned native location tag.
        edge_ids: Canonical locus or minimizing edge identities.
        node_ids: Canonical locus node identities.
        parameter_root_id: Exact algebraic parameter identity when applicable.

    Returns:
        One validated exact neck-locus variant.

    Raises:
        InvalidNeckEvidenceError: If the fields contradict the selected
            location variant.
    """
    _exact_bytes(tag, "neck location tag")
    if tag == STRICT_EDGE_NECK_TAG:
        if len(edge_ids) != 1 or node_ids or parameter_root_id is None:
            raise InvalidNeckEvidenceError("strict-edge neck location fields are incomplete.")
        return StrictEdgeNeckLocus.build(
            edge_id=edge_ids[0],
            parameter_root_id=parameter_root_id,
        )
    if tag == CLEARANCE_ENDPOINT_NECK_TAG:
        if len(edge_ids) != 1 or len(node_ids) != 1 or parameter_root_id is None:
            raise InvalidNeckEvidenceError("clearance-endpoint neck location fields are incomplete.")
        return ClearanceEndpointNeckLocus.build(
            edge_id=edge_ids[0],
            node_id=node_ids[0],
            parameter_root_id=parameter_root_id,
        )
    if tag == SHARED_VERTEX_NECK_TAG:
        if len(edge_ids) < 2 or len(node_ids) != 1 or parameter_root_id is not None:
            raise InvalidNeckEvidenceError("shared-vertex neck location fields are incomplete.")
        return SharedVertexNeckLocus.build(
            node_id=node_ids[0],
            minimizing_edge_ids=edge_ids,
        )
    if tag == PLATEAU_NECK_TAG:
        if not edge_ids or not node_ids or parameter_root_id is not None:
            raise InvalidNeckEvidenceError("plateau neck location fields are incomplete.")
        return PlateauNeckLocus.build(
            node_ids=node_ids,
            edge_ids=edge_ids,
        )
    raise InvalidNeckEvidenceError("native neck location uses an unknown versioned tag.")


def neck_locus_edge_ids(locus: NeckLocus) -> tuple[MatEdgeId, ...]:
    if type(locus) is StrictEdgeNeckLocus:
        return (locus.edge_id,)
    if type(locus) is ClearanceEndpointNeckLocus:
        return (locus.edge_id,)
    if type(locus) is SharedVertexNeckLocus:
        return locus.minimizing_edge_ids
    if type(locus) is PlateauNeckLocus:
        return locus.edge_ids
    raise InvalidNeckEvidenceError("neck locus is outside the closed exact type.")


def neck_locus_node_ids(locus: NeckLocus) -> tuple[MatNodeId, ...]:
    if type(locus) is StrictEdgeNeckLocus:
        return ()
    if type(locus) is ClearanceEndpointNeckLocus:
        return (locus.node_id,)
    if type(locus) is SharedVertexNeckLocus:
        return (locus.node_id,)
    if type(locus) is PlateauNeckLocus:
        return locus.node_ids
    raise InvalidNeckEvidenceError("neck locus is outside the closed exact type.")


def _neck_side_bytes(
    *,
    neck_evidence_digest: IdentityDigest,
    partition_ordinal: int,
    edge_ids: tuple[MatEdgeId, ...],
) -> bytes:
    return encode_tagged_union(
        b"mat-neck-side-v1",
        encode_component_map(
            {
                b"edge-ids": encode_sequence(tuple(encode_bytes(bytes(edge_id)) for edge_id in edge_ids)),
                b"neck-evidence-digest": encode_bytes(bytes(neck_evidence_digest)),
                b"partition-ordinal": encode_integer(partition_ordinal),
            }
        ),
    )


@dataclass(frozen=True)
class NeckSide:
    identity: NeckSideId
    neck_evidence_digest: IdentityDigest
    partition_ordinal: int
    edge_ids: tuple[MatEdgeId, ...]

    def __post_init__(self) -> None:
        _exact_bytes(self.identity, "neck side identity", digest=True)
        _exact_bytes(self.neck_evidence_digest, "neck side evidence digest", digest=True)
        if type(self.partition_ordinal) is not int or self.partition_ordinal < 0:
            raise InvalidNeckEvidenceError("neck side partition ordinal must be non-negative.")
        _canonical_ids(self.edge_ids, "neck side edges")
        if self.identity != hashlib.sha256(self.canonical_bytes).digest():
            raise InvalidNeckEvidenceError("neck side identity contradicts its exact partition.")

    @classmethod
    def build(
        cls,
        *,
        neck_evidence_digest: IdentityDigest,
        partition_ordinal: int,
        edge_ids: tuple[MatEdgeId, ...],
    ) -> Self:
        canonical = _neck_side_bytes(
            neck_evidence_digest=neck_evidence_digest,
            partition_ordinal=partition_ordinal,
            edge_ids=edge_ids,
        )
        return cls(
            NeckSideId(hashlib.sha256(canonical).digest()),
            neck_evidence_digest,
            partition_ordinal,
            edge_ids,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return _neck_side_bytes(
            neck_evidence_digest=self.neck_evidence_digest,
            partition_ordinal=self.partition_ordinal,
            edge_ids=self.edge_ids,
        )
