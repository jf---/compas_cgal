import hashlib
import math
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal
from typing import NewType
from typing import Self
from typing import TypeAlias
from typing import cast

import numpy as np
import numpy.typing as npt

from compas_cgal import _medial_axis_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidMedialAxisProjectionError
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import Clearance
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

MatNodeId = NewType("MatNodeId", bytes)
MatEdgeId = NewType("MatEdgeId", bytes)
MatSiteId = NewType("MatSiteId", bytes)
MatOriginalDualId = NewType("MatOriginalDualId", bytes)
MatParameterId = NewType("MatParameterId", bytes)
MatParameter = NewType("MatParameter", float)
MatCurveKind: TypeAlias = Literal["line", "parabola"]
MatOriginalDualKind: TypeAlias = Literal[
    "point-point",
    "point-segment",
    "segment-segment",
]
MatSiteKind: TypeAlias = Literal["point", "open-segment"]

_CURVE_KINDS: dict[int, MatCurveKind] = {
    0: "line",
    1: "parabola",
}
_DUAL_KINDS: dict[int, MatOriginalDualKind] = {
    0: "point-point",
    1: "point-segment",
    2: "segment-segment",
}
_SITE_KINDS: dict[int, MatSiteKind] = {
    0: "point",
    1: "open-segment",
}
_ENDPOINT_ORIGINAL_VERTEX = 1
_ENDPOINT_DOMAIN_BOUNDARY = 2
_ENDPOINT_CLEARANCE_ROOT = 4
_SHA256_BYTES = hashlib.sha256().digest_size
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))


def _identity(value: object, name: str) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidMedialAxisProjectionError(f"{name} must be nonempty exact bytes.")
    return value


def _int64_array(
    value: object,
    *,
    rank: int,
    name: str,
) -> npt.NDArray[np.int64]:
    if not isinstance(value, np.ndarray) or value.dtype != np.dtype(np.int64) or value.ndim != rank:
        raise InvalidMedialAxisProjectionError(f"{name} must be one rank-{rank} native int64 array.")
    return cast(npt.NDArray[np.int64], value)


def _float64_array(
    value: object,
    *,
    rank: int,
    name: str,
) -> npt.NDArray[np.float64]:
    if not isinstance(value, np.ndarray) or value.dtype != np.dtype(np.float64) or value.ndim != rank:
        raise InvalidMedialAxisProjectionError(f"{name} must be one rank-{rank} native float64 array.")
    if not np.all(np.isfinite(value)):
        raise InvalidMedialAxisProjectionError(f"{name} must contain finite reporting values.")
    return cast(npt.NDArray[np.float64], value)


def _csr_offsets(
    value: object,
    *,
    expected_size: int,
    terminal: int,
    name: str,
) -> npt.NDArray[np.int64]:
    offsets = _int64_array(value, rank=1, name=name)
    if offsets.shape != (expected_size,) or offsets[0] != 0 or offsets[-1] != terminal or np.any(offsets[1:] < offsets[:-1]):
        raise InvalidMedialAxisProjectionError(f"{name} is not one complete canonical CSR offset vector.")
    return offsets


def _component_identity(
    node_ids: tuple[MatNodeId, ...],
    edge_ids: tuple[MatEdgeId, ...],
) -> ComponentId:
    canonical = encode_tagged_union(
        b"mat-component-v1",
        encode_component_map(
            {
                b"edges": encode_sequence(tuple(encode_bytes(bytes(identity)) for identity in edge_ids)),
                b"nodes": encode_sequence(tuple(encode_bytes(bytes(identity)) for identity in node_ids)),
            }
        ),
    )
    return ComponentId(hashlib.sha256(canonical).digest())


def _branch_identity(edge_id: MatEdgeId) -> BranchId:
    canonical = encode_tagged_union(b"mat-edge-branch-v1", encode_bytes(bytes(edge_id)))
    return BranchId(hashlib.sha256(canonical).digest())


def _cursor_identity(edge_id: MatEdgeId, parameter_id: MatParameterId) -> CursorIdentity:
    canonical = encode_tagged_union(
        b"mat-native-sample-cursor-v1",
        encode_component_map(
            {
                b"edge-id": encode_bytes(bytes(edge_id)),
                b"parameter-id": encode_bytes(bytes(parameter_id)),
            }
        ),
    )
    return CursorIdentity(hashlib.sha256(canonical).digest())


def _native_ring(ring: CanonicalRingV1) -> npt.NDArray[np.float64]:
    return np.asarray(
        tuple((point.x, point.y, 0.0) for point in ring.vertices),
        dtype=np.float64,
    )


@dataclass(frozen=True)
class MatSite:
    identity: MatSiteId
    kind: MatSiteKind
    ring_ordinal: int
    feature_ordinal: int
    source: Point2[WorldXY]
    target: Point2[WorldXY] | None

    def __post_init__(self) -> None:
        _identity(self.identity, "MAT site identity")
        if self.kind not in ("point", "open-segment"):
            raise InvalidMedialAxisProjectionError("MAT site has an unknown closed kind.")
        if type(self.ring_ordinal) is not int or self.ring_ordinal < 0:
            raise InvalidMedialAxisProjectionError("MAT site ring ordinal must be non-negative.")
        if type(self.feature_ordinal) is not int or self.feature_ordinal < 0:
            raise InvalidMedialAxisProjectionError("MAT site feature ordinal must be non-negative.")
        if type(self.source) is not Point2:
            raise InvalidMedialAxisProjectionError("MAT site source must be one world-XY point.")
        if (self.kind == "point") != (self.target is None):
            raise InvalidMedialAxisProjectionError("MAT site kind contradicts its boundary geometry.")
        if self.target is not None and type(self.target) is not Point2:
            raise InvalidMedialAxisProjectionError("MAT segment site target must be one world-XY point.")


@dataclass(frozen=True)
class MatBoundaryFeature:
    domain_kind: int
    component: int
    curve_kind: int
    source_site_or_ring: int
    derived_feature_index: int

    def __post_init__(self) -> None:
        if any(
            type(value) is not int or value < 0
            for value in (
                self.domain_kind,
                self.component,
                self.curve_kind,
                self.source_site_or_ring,
                self.derived_feature_index,
            )
        ):
            raise InvalidMedialAxisProjectionError("MAT endpoint feature fields must be non-negative integers.")


@dataclass(frozen=True)
class MatEndpointState:
    flags: int
    features: tuple[MatBoundaryFeature, ...]

    def __post_init__(self) -> None:
        if type(self.flags) is not int or self.flags <= 0:
            raise InvalidMedialAxisProjectionError("MAT endpoint requires nonempty exact provenance flags.")
        if self.flags & ~(_ENDPOINT_ORIGINAL_VERTEX | _ENDPOINT_DOMAIN_BOUNDARY | _ENDPOINT_CLEARANCE_ROOT):
            raise InvalidMedialAxisProjectionError("MAT endpoint contains an unknown provenance flag.")
        if type(self.features) is not tuple or any(type(feature) is not MatBoundaryFeature for feature in self.features):
            raise InvalidMedialAxisProjectionError("MAT endpoint features must be one immutable typed tuple.")

    @property
    def original_voronoi_vertex(self) -> bool:
        return bool(self.flags & _ENDPOINT_ORIGINAL_VERTEX)

    @property
    def domain_boundary(self) -> bool:
        return bool(self.flags & _ENDPOINT_DOMAIN_BOUNDARY)

    @property
    def clearance_root(self) -> bool:
        return bool(self.flags & _ENDPOINT_CLEARANCE_ROOT)


@dataclass(frozen=True)
class MatNode:
    identity: MatNodeId
    point: Point2[WorldXY]
    generator_site_ids: tuple[MatSiteId, ...]

    def __post_init__(self) -> None:
        _identity(self.identity, "MAT node identity")
        if type(self.point) is not Point2:
            raise InvalidMedialAxisProjectionError("MAT node requires one world-XY reporting point.")
        if type(self.generator_site_ids) is not tuple or not self.generator_site_ids or self.generator_site_ids != tuple(sorted(set(self.generator_site_ids))):
            raise InvalidMedialAxisProjectionError("MAT node site identities must be nonempty, unique, and canonical.")


@dataclass(frozen=True)
class MatEdge:
    identity: MatEdgeId
    source: MatNode
    target: MatNode
    curve_kind: MatCurveKind
    generator_site_ids: tuple[MatSiteId, MatSiteId]
    original_dual_kind: MatOriginalDualKind
    original_dual_identity: MatOriginalDualId
    clip_component_index: int
    source_endpoint: MatEndpointState
    target_endpoint: MatEndpointState
    branch_id: BranchId

    def __post_init__(self) -> None:
        _identity(self.identity, "MAT edge identity")
        if type(self.source) is not MatNode or type(self.target) is not MatNode or self.source == self.target:
            raise InvalidMedialAxisProjectionError("MAT edge requires two distinct exact nodes.")
        if self.curve_kind not in ("line", "parabola"):
            raise InvalidMedialAxisProjectionError("MAT edge has an unknown curve kind.")
        if type(self.generator_site_ids) is not tuple or len(self.generator_site_ids) != 2 or self.generator_site_ids != tuple(sorted(set(self.generator_site_ids))):
            raise InvalidMedialAxisProjectionError("MAT edge requires two canonical generator-site identities.")
        if self.original_dual_kind not in ("point-point", "point-segment", "segment-segment"):
            raise InvalidMedialAxisProjectionError("MAT edge has an unknown original-dual kind.")
        _identity(self.original_dual_identity, "MAT original-dual identity")
        if type(self.clip_component_index) is not int or self.clip_component_index < 0:
            raise InvalidMedialAxisProjectionError("MAT edge clip-component index must be non-negative.")
        if type(self.source_endpoint) is not MatEndpointState or type(self.target_endpoint) is not MatEndpointState:
            raise InvalidMedialAxisProjectionError("MAT edge endpoints must retain exact clipped state.")
        if type(self.branch_id) is not bytes or len(self.branch_id) != _SHA256_BYTES:
            raise InvalidMedialAxisProjectionError("MAT edge branch ID must be one SHA-256 digest.")


@dataclass(frozen=True)
class MatToolFitRun:
    parent_edge_id: MatEdgeId
    source_node_id: MatNodeId
    target_node_id: MatNodeId
    source_endpoint: MatEndpointState
    target_endpoint: MatEndpointState

    def __post_init__(self) -> None:
        _identity(self.parent_edge_id, "tool-fit run parent edge")
        _identity(self.source_node_id, "tool-fit run source node")
        _identity(self.target_node_id, "tool-fit run target node")
        if self.source_node_id == self.target_node_id:
            raise InvalidMedialAxisProjectionError("tool-fit run endpoints must be distinct.")


@dataclass(frozen=True)
class MatSample:
    exact_parameter_id: MatParameterId
    cursor_identity: CursorIdentity
    edge_id: MatEdgeId
    point: Point2[WorldXY]
    clearance: Clearance
    reported_guide_radius: Clearance
    parameter: MatParameter
    ordinal_on_edge: int

    def __post_init__(self) -> None:
        _identity(self.exact_parameter_id, "MAT sample parameter identity")
        if type(self.cursor_identity) is not bytes or len(self.cursor_identity) != _SHA256_BYTES:
            raise InvalidMedialAxisProjectionError("MAT sample cursor must be one SHA-256 identity.")
        _identity(self.edge_id, "MAT sample edge identity")
        if type(self.point) is not Point2 or type(self.clearance) is not Clearance:
            raise InvalidMedialAxisProjectionError("MAT sample requires typed reporting geometry.")
        if type(self.reported_guide_radius) is not Clearance:
            raise InvalidMedialAxisProjectionError("MAT sample guide-radius report must be typed non-negative length.")
        if type(self.parameter) is not float or not math.isfinite(self.parameter):
            raise InvalidMedialAxisProjectionError("MAT sample parameter must be one finite reporting scalar.")
        if type(self.ordinal_on_edge) is not int or self.ordinal_on_edge < 0:
            raise InvalidMedialAxisProjectionError("MAT sample edge ordinal must be non-negative.")


@dataclass(frozen=True)
class MatComponent:
    identity: ComponentId
    node_ids: tuple[MatNodeId, ...]
    edge_ids: tuple[MatEdgeId, ...]

    def __post_init__(self) -> None:
        if type(self.identity) is not bytes or len(self.identity) != _SHA256_BYTES:
            raise InvalidMedialAxisProjectionError("MAT component ID must be one SHA-256 digest.")
        if type(self.node_ids) is not tuple or not self.node_ids or self.node_ids != tuple(sorted(set(self.node_ids))):
            raise InvalidMedialAxisProjectionError("MAT component nodes must be nonempty and canonical.")
        if type(self.edge_ids) is not tuple or not self.edge_ids or self.edge_ids != tuple(sorted(set(self.edge_ids))):
            raise InvalidMedialAxisProjectionError("MAT component edges must be nonempty and canonical.")


@dataclass(frozen=True)
class MatTopology:
    sites: tuple[MatSite, ...]
    nodes: tuple[MatNode, ...]
    edges: tuple[MatEdge, ...]
    components: tuple[MatComponent, ...]
    tool_fit_runs: tuple[MatToolFitRun, ...]
    site_by_id: Mapping[MatSiteId, MatSite]
    node_by_id: Mapping[MatNodeId, MatNode]
    edge_by_id: Mapping[MatEdgeId, MatEdge]
    component_by_edge_id: Mapping[MatEdgeId, ComponentId]

    def __post_init__(self) -> None:
        if type(self.sites) is not tuple or not self.sites or any(type(site) is not MatSite for site in self.sites):
            raise InvalidMedialAxisProjectionError("MAT topology sites must be one nonempty immutable typed tuple.")
        if type(self.nodes) is not tuple or not self.nodes or any(type(node) is not MatNode for node in self.nodes):
            raise InvalidMedialAxisProjectionError("MAT topology nodes must be one nonempty immutable typed tuple.")
        if type(self.edges) is not tuple or not self.edges or any(type(edge) is not MatEdge for edge in self.edges):
            raise InvalidMedialAxisProjectionError("MAT topology edges must be one nonempty immutable typed tuple.")
        if type(self.components) is not tuple or not self.components or any(type(component) is not MatComponent for component in self.components):
            raise InvalidMedialAxisProjectionError("MAT topology components must be one nonempty immutable typed tuple.")
        if type(self.tool_fit_runs) is not tuple or any(type(run) is not MatToolFitRun for run in self.tool_fit_runs):
            raise InvalidMedialAxisProjectionError("MAT tool-fit runs must be one immutable typed tuple.")
        if any(
            type(mapping) is not _MAPPING_PROXY_TYPE
            for mapping in (
                self.site_by_id,
                self.node_by_id,
                self.edge_by_id,
                self.component_by_edge_id,
            )
        ):
            raise InvalidMedialAxisProjectionError("MAT topology lookup views must be immutable mapping proxies.")

        expected_sites = {site.identity: site for site in self.sites}
        expected_nodes = {node.identity: node for node in self.nodes}
        expected_edges = {edge.identity: edge for edge in self.edges}
        if (
            len(expected_sites) != len(self.sites)
            or len(expected_nodes) != len(self.nodes)
            or len(expected_edges) != len(self.edges)
            or dict(self.site_by_id) != expected_sites
            or dict(self.node_by_id) != expected_nodes
            or dict(self.edge_by_id) != expected_edges
        ):
            raise InvalidMedialAxisProjectionError("MAT topology lookup views contradict their owned records.")
        if any(
            edge.source != expected_nodes.get(edge.source.identity)
            or edge.target != expected_nodes.get(edge.target.identity)
            or any(site_id not in expected_sites for site_id in edge.generator_site_ids)
            for edge in self.edges
        ):
            raise InvalidMedialAxisProjectionError("MAT edge ownership contradicts topology node or site records.")

        expected_component_by_edge: dict[MatEdgeId, ComponentId] = {}
        component_node_ids: set[MatNodeId] = set()
        for component in self.components:
            if any(node_id not in expected_nodes for node_id in component.node_ids):
                raise InvalidMedialAxisProjectionError("MAT component references an unknown node.")
            component_node_ids.update(component.node_ids)
            for edge_id in component.edge_ids:
                if edge_id not in expected_edges or edge_id in expected_component_by_edge:
                    raise InvalidMedialAxisProjectionError("MAT component edge ownership is incomplete or duplicated.")
                expected_component_by_edge[edge_id] = component.identity
        if (
            component_node_ids != set(expected_nodes)
            or expected_component_by_edge.keys() != expected_edges.keys()
            or dict(self.component_by_edge_id) != expected_component_by_edge
        ):
            raise InvalidMedialAxisProjectionError("MAT component lookup contradicts the complete graph partition.")

        runs_by_edge = {run.parent_edge_id: run for run in self.tool_fit_runs}
        if len(runs_by_edge) != len(self.tool_fit_runs) or runs_by_edge.keys() != expected_edges.keys():
            raise InvalidMedialAxisProjectionError("MAT tool-fit runs must cover every exact edge once.")
        if any(
            run.source_node_id != expected_edges[edge_id].source.identity
            or run.target_node_id != expected_edges[edge_id].target.identity
            or run.source_endpoint != expected_edges[edge_id].source_endpoint
            or run.target_endpoint != expected_edges[edge_id].target_endpoint
            for edge_id, run in runs_by_edge.items()
        ):
            raise InvalidMedialAxisProjectionError("MAT tool-fit run contradicts its exact parent edge.")


@dataclass(frozen=True)
class MatProposalSampling:
    samples: tuple[MatSample, ...]

    def __post_init__(self) -> None:
        if type(self.samples) is not tuple or not self.samples or any(type(sample) is not MatSample for sample in self.samples):
            raise InvalidMedialAxisProjectionError("MAT proposal samples must be one nonempty immutable typed tuple.")
        cursor_ids = tuple(sample.cursor_identity for sample in self.samples)
        if len(cursor_ids) != len(set(cursor_ids)):
            raise InvalidMedialAxisProjectionError("MAT proposal cursor identities must be unique.")


@dataclass(frozen=True)
class MatProof:
    center_domain_digest: bytes
    mat_certificate: bytes
    native_owner: _medial_axis_2.SegmentSiteMedialAxis

    def __post_init__(self) -> None:
        if type(self.center_domain_digest) is not bytes or len(self.center_domain_digest) != _SHA256_BYTES:
            raise InvalidMedialAxisProjectionError("MAT center-domain identity must be one SHA-256 digest.")
        _identity(self.mat_certificate, "MAT certificate")
        if type(self.native_owner) is not _medial_axis_2.SegmentSiteMedialAxis:
            raise InvalidMedialAxisProjectionError("MAT proof must retain its exact native owner.")


@dataclass(frozen=True)
class MedialAxis:
    design_boundary: CanonicalRingV1
    holes: tuple[CanonicalRingV1, ...]
    tool_radius: ToolRadius
    station_spacing: Spacing
    max_sagitta: ChordBound
    max_refinement_depth: int
    topology: MatTopology
    proposal_sampling: MatProposalSampling
    proof: MatProof

    def __post_init__(self) -> None:
        if type(self.design_boundary) is not CanonicalRingV1 or not self.design_boundary.is_outer:
            raise InvalidMedialAxisProjectionError("typed MAT requires one canonical outer design boundary.")
        if type(self.holes) is not tuple or any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in self.holes):
            raise InvalidMedialAxisProjectionError("typed MAT holes must be one immutable canonical tuple.")
        if tuple(sorted(self.holes, key=lambda hole: hole.canonical_bytes)) != self.holes:
            raise InvalidMedialAxisProjectionError("typed MAT holes must use canonical identity order.")
        if type(self.tool_radius) is not ToolRadius or type(self.station_spacing) is not Spacing or type(self.max_sagitta) is not ChordBound:
            raise InvalidMedialAxisProjectionError("typed MAT policies must retain their unit-bearing types.")
        if type(self.max_refinement_depth) is not int or self.max_refinement_depth <= 0:
            raise InvalidMedialAxisProjectionError("typed MAT refinement depth must be one positive integer.")
        if type(self.topology) is not MatTopology or type(self.proposal_sampling) is not MatProposalSampling or type(self.proof) is not MatProof:
            raise InvalidMedialAxisProjectionError("typed MAT aggregates must use their exact closed domain types.")
        edge_ids = set(self.topology.edge_by_id)
        sampled_edge_ids = {sample.edge_id for sample in self.proposal_sampling.samples}
        if sampled_edge_ids != edge_ids:
            raise InvalidMedialAxisProjectionError("typed MAT proposal samples must cover every exact edge.")

    @classmethod
    def build(
        cls,
        *,
        design_boundary: CanonicalRingV1,
        holes: Sequence[CanonicalRingV1],
        tool_radius: ToolRadius,
        station_spacing: Spacing,
        max_sagitta: ChordBound,
        max_refinement_depth: int,
    ) -> Self:
        if type(design_boundary) is not CanonicalRingV1 or not design_boundary.is_outer:
            raise InvalidMedialAxisProjectionError("MAT design boundary must be one canonical outer ring.")
        try:
            canonical_holes = tuple(holes)
        except TypeError:
            raise InvalidMedialAxisProjectionError("MAT holes must be one finite canonical sequence.") from None
        if any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in canonical_holes):
            raise InvalidMedialAxisProjectionError("MAT holes must contain only canonical hole rings.")
        if tuple(sorted(canonical_holes, key=lambda hole: hole.canonical_bytes)) != canonical_holes:
            raise InvalidMedialAxisProjectionError("MAT holes must use canonical identity order.")
        if type(tool_radius) is not ToolRadius:
            raise InvalidMedialAxisProjectionError("MAT tool radius must be typed.")
        if type(station_spacing) is not Spacing:
            raise InvalidMedialAxisProjectionError("MAT station spacing must be typed.")
        if type(max_sagitta) is not ChordBound:
            raise InvalidMedialAxisProjectionError("MAT sagitta bound must be typed.")
        if type(max_refinement_depth) is not int or max_refinement_depth <= 0:
            raise InvalidMedialAxisProjectionError("MAT refinement depth must be one positive integer.")

        native_owner = _medial_axis_2.SegmentSiteMedialAxis.build(
            _native_ring(design_boundary),
            tuple(_native_ring(hole) for hole in canonical_holes),
            float(tool_radius.value),
            float(station_spacing.value),
            float(max_sagitta.value),
            max_refinement_depth,
        )
        topology, proposal_sampling, proof = _project_native(
            native_owner,
            design_boundary,
            canonical_holes,
        )
        return cls(
            design_boundary,
            canonical_holes,
            tool_radius,
            station_spacing,
            max_sagitta,
            max_refinement_depth,
            topology,
            proposal_sampling,
            proof,
        )

    @property
    def sites(self) -> tuple[MatSite, ...]:
        return self.topology.sites

    @property
    def nodes(self) -> tuple[MatNode, ...]:
        return self.topology.nodes

    @property
    def edges(self) -> tuple[MatEdge, ...]:
        return self.topology.edges

    @property
    def components(self) -> tuple[MatComponent, ...]:
        return self.topology.components

    @property
    def tool_fit_runs(self) -> tuple[MatToolFitRun, ...]:
        return self.topology.tool_fit_runs

    @property
    def samples(self) -> tuple[MatSample, ...]:
        return self.proposal_sampling.samples

    @property
    def site_by_id(self) -> Mapping[MatSiteId, MatSite]:
        return self.topology.site_by_id

    @property
    def node_by_id(self) -> Mapping[MatNodeId, MatNode]:
        return self.topology.node_by_id

    @property
    def edge_by_id(self) -> Mapping[MatEdgeId, MatEdge]:
        return self.topology.edge_by_id

    @property
    def component_by_edge_id(self) -> Mapping[MatEdgeId, ComponentId]:
        return self.topology.component_by_edge_id

    @property
    def center_domain_digest(self) -> bytes:
        return self.proof.center_domain_digest

    @property
    def mat_certificate(self) -> bytes:
        return self.proof.mat_certificate

    @property
    def native_owner(self) -> _medial_axis_2.SegmentSiteMedialAxis:
        return self.proof.native_owner


def _sites(
    site_ids: tuple[bytes, ...],
    provenance: npt.NDArray[np.int64],
    boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
) -> tuple[MatSite, ...]:
    if provenance.shape != (len(site_ids), 3):
        raise InvalidMedialAxisProjectionError("MAT site provenance cardinality contradicts native site identities.")
    rings = (boundary, *holes)
    sites: list[MatSite] = []
    for index, row in enumerate(provenance):
        kind_value, ring_ordinal, feature_ordinal = (int(value) for value in row)
        kind = _SITE_KINDS.get(kind_value)
        if kind is None or ring_ordinal >= len(rings):
            raise InvalidMedialAxisProjectionError("MAT site provenance references an unknown kind or ring.")
        ring = rings[ring_ordinal]
        if feature_ordinal >= ring.vertex_count:
            raise InvalidMedialAxisProjectionError("MAT site provenance references an unknown boundary feature.")
        source = ring.vertices[feature_ordinal]
        target = None if kind == "point" else ring.vertices[(feature_ordinal + 1) % ring.vertex_count]
        sites.append(
            MatSite(
                MatSiteId(_identity(site_ids[index], "MAT site identity")),
                kind,
                ring_ordinal,
                feature_ordinal,
                source,
                target,
            )
        )
    return tuple(sites)


def _endpoint_features(
    features: npt.NDArray[np.int64],
    offsets: npt.NDArray[np.int64],
    endpoint_index: int,
) -> tuple[MatBoundaryFeature, ...]:
    begin = int(offsets[endpoint_index])
    end = int(offsets[endpoint_index + 1])
    return tuple(MatBoundaryFeature(*(int(value) for value in row)) for row in features[begin:end])


def _connected_components(
    nodes: tuple[MatNode, ...],
    edges: tuple[MatEdge, ...],
) -> tuple[tuple[MatComponent, ...], Mapping[MatEdgeId, ComponentId]]:
    incident: dict[MatNodeId, list[MatEdge]] = {node.identity: [] for node in nodes}
    for edge in edges:
        incident[edge.source.identity].append(edge)
        incident[edge.target.identity].append(edge)

    unseen = set(incident)
    components: list[MatComponent] = []
    edge_components: dict[MatEdgeId, ComponentId] = {}
    while unseen:
        seed = min(unseen)
        pending = [seed]
        component_nodes: set[MatNodeId] = set()
        component_edges: set[MatEdgeId] = set()
        while pending:
            node_id = pending.pop()
            if node_id in component_nodes:
                continue
            component_nodes.add(node_id)
            unseen.discard(node_id)
            for edge in incident[node_id]:
                component_edges.add(edge.identity)
                other = edge.target.identity if edge.source.identity == node_id else edge.source.identity
                if other not in component_nodes:
                    pending.append(other)
        ordered_nodes = tuple(sorted(component_nodes))
        ordered_edges = tuple(sorted(component_edges))
        component_id = _component_identity(ordered_nodes, ordered_edges)
        component = MatComponent(component_id, ordered_nodes, ordered_edges)
        components.append(component)
        for edge_id in ordered_edges:
            edge_components[edge_id] = component_id
    components.sort(key=lambda component: component.identity)
    return tuple(components), MappingProxyType(edge_components)


def _project_native(
    native_owner: _medial_axis_2.SegmentSiteMedialAxis,
    boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
) -> tuple[MatTopology, MatProposalSampling, MatProof]:
    projection = native_owner.projection
    if type(projection) is not tuple or len(projection) != 20:
        raise InvalidMedialAxisProjectionError("MAT native owner must project the fixed 20-field contract.")

    node_coordinates = _float64_array(projection[0], rank=2, name="MAT nodes")
    edge_rows = _int64_array(projection[1], rank=2, name="MAT edges")
    node_site_ids = _int64_array(projection[3], rank=1, name="MAT node site IDs")
    site_provenance = _int64_array(projection[4], rank=2, name="MAT site provenance")
    endpoint_flags = _int64_array(projection[5], rank=2, name="MAT endpoint flags")
    endpoint_features = _int64_array(projection[7], rank=2, name="MAT endpoint features")
    edge_exact_flags = _int64_array(projection[8], rank=2, name="MAT edge exact flags")
    sample_centers = _float64_array(projection[9], rank=2, name="MAT samples")
    sample_clearance = _float64_array(projection[10], rank=1, name="MAT sample clearance")
    sample_guide_radius = _float64_array(projection[11], rank=1, name="MAT sample guide radius")
    sample_flags = _int64_array(projection[12], rank=2, name="MAT sample exact flags")
    sample_parameters = _float64_array(projection[14], rank=1, name="MAT sample parameters")
    neck_cut_edge_ids = _int64_array(projection[17], rank=1, name="MAT neck cut edge IDs")

    node_count = node_coordinates.shape[0]
    edge_count = edge_rows.shape[0]
    site_count = site_provenance.shape[0]
    sample_count = sample_centers.shape[0]
    if (
        node_coordinates.shape[1:] != (3,)
        or edge_rows.shape[1:] != (8,)
        or site_provenance.shape[1:] != (3,)
        or endpoint_flags.shape != (edge_count, 2)
        or endpoint_features.shape[1:] != (5,)
        or edge_exact_flags.shape != (edge_count, 3)
        or sample_centers.shape[1:] != (3,)
        or sample_clearance.shape != (sample_count,)
        or sample_guide_radius.shape != (sample_count,)
        or sample_flags.shape != (sample_count, 2)
        or sample_parameters.shape != (sample_count,)
    ):
        raise InvalidMedialAxisProjectionError("MAT native projection has contradictory ranks or cardinalities.")
    if np.any(node_coordinates[:, 2] != 0.0) or np.any(sample_centers[:, 2] != 0.0):
        raise InvalidMedialAxisProjectionError("MAT reporting geometry must lie in world XY.")
    if np.any(edge_exact_flags != 1) or np.any(sample_flags != 1):
        raise InvalidMedialAxisProjectionError("MAT native projection contains unverified exact verdicts.")
    if np.any(sample_clearance < 0.0) or np.any(sample_guide_radius < 0.0):
        raise InvalidMedialAxisProjectionError("MAT sample reporting lengths must be non-negative.")

    node_site_offsets = _csr_offsets(
        projection[2],
        expected_size=node_count + 1,
        terminal=node_site_ids.shape[0],
        name="MAT node-site offsets",
    )
    endpoint_feature_offsets = _csr_offsets(
        projection[6],
        expected_size=2 * edge_count + 1,
        terminal=endpoint_features.shape[0],
        name="MAT endpoint-feature offsets",
    )
    edge_sample_offsets = _csr_offsets(
        projection[13],
        expected_size=edge_count + 1,
        terminal=sample_count,
        name="MAT edge-sample offsets",
    )
    neck_evidence = projection[15]
    if type(neck_evidence) is not tuple or any(type(record) is not bytes or not record for record in neck_evidence):
        raise InvalidMedialAxisProjectionError("MAT neck evidence must be one immutable exact byte tuple.")
    _csr_offsets(
        projection[16],
        expected_size=len(neck_evidence) + 1,
        terminal=neck_cut_edge_ids.shape[0],
        name="MAT neck-cut offsets",
    )
    if np.any(neck_cut_edge_ids < 0) or np.any(neck_cut_edge_ids >= edge_count):
        raise InvalidMedialAxisProjectionError("MAT neck cut references an unknown edge.")

    node_ids = tuple(MatNodeId(_identity(identity, "MAT node identity")) for identity in native_owner.node_ids)
    edge_ids = tuple(MatEdgeId(_identity(identity, "MAT edge identity")) for identity in native_owner.edge_ids)
    site_ids = tuple(_identity(identity, "MAT site identity") for identity in native_owner.site_ids)
    original_dual_ids = tuple(MatOriginalDualId(_identity(identity, "MAT original-dual identity")) for identity in native_owner.original_dual_ids)
    parameter_ids = tuple(MatParameterId(_identity(identity, "MAT sample parameter identity")) for identity in native_owner.sample_parameter_ids)
    if len(node_ids) != node_count or len(edge_ids) != edge_count or len(parameter_ids) != sample_count:
        raise InvalidMedialAxisProjectionError("MAT exact identity cardinality contradicts its numeric projection.")
    if node_ids != tuple(sorted(set(node_ids))) or edge_ids != tuple(sorted(set(edge_ids))):
        raise InvalidMedialAxisProjectionError("MAT node and edge identities must be unique and canonical.")

    sites = _sites(site_ids, site_provenance, boundary, holes)
    site_identity_by_index = tuple(site.identity for site in sites)
    if np.any(node_site_ids < 0) or np.any(node_site_ids >= site_count):
        raise InvalidMedialAxisProjectionError("MAT node-site CSR references an unknown site.")
    nodes = tuple(
        MatNode(
            node_ids[index],
            Point2[WorldXY].build(node_coordinates[index, 0], node_coordinates[index, 1]),
            tuple(site_identity_by_index[int(site_index)] for site_index in node_site_ids[int(node_site_offsets[index]) : int(node_site_offsets[index + 1])]),
        )
        for index in range(node_count)
    )

    edges: list[MatEdge] = []
    for index, row in enumerate(edge_rows):
        source_index, target_index = int(row[0]), int(row[1])
        if source_index < 0 or source_index >= node_count or target_index < 0 or target_index >= node_count:
            raise InvalidMedialAxisProjectionError("MAT edge references an unknown node.")
        curve_kind = _CURVE_KINDS.get(int(row[2]))
        dual_kind = _DUAL_KINDS.get(int(row[5]))
        first_site, second_site = int(row[3]), int(row[4])
        dual_index = int(row[6])
        if (
            curve_kind is None
            or dual_kind is None
            or first_site < 0
            or second_site >= site_count
            or first_site >= second_site
            or dual_index < 0
            or dual_index >= len(original_dual_ids)
        ):
            raise InvalidMedialAxisProjectionError("MAT edge row contains an unknown native identity index.")
        edges.append(
            MatEdge(
                edge_ids[index],
                nodes[source_index],
                nodes[target_index],
                curve_kind,
                (site_identity_by_index[first_site], site_identity_by_index[second_site]),
                dual_kind,
                original_dual_ids[dual_index],
                int(row[7]),
                MatEndpointState(
                    int(endpoint_flags[index, 0]),
                    _endpoint_features(endpoint_features, endpoint_feature_offsets, 2 * index),
                ),
                MatEndpointState(
                    int(endpoint_flags[index, 1]),
                    _endpoint_features(endpoint_features, endpoint_feature_offsets, 2 * index + 1),
                ),
                _branch_identity(edge_ids[index]),
            )
        )
    frozen_edges = tuple(edges)
    edge_by_id = MappingProxyType({edge.identity: edge for edge in frozen_edges})
    if len(edge_by_id) != edge_count:
        raise InvalidMedialAxisProjectionError("MAT edge identities are not unique.")

    samples: list[MatSample] = []
    for edge_index, edge in enumerate(frozen_edges):
        begin = int(edge_sample_offsets[edge_index])
        end = int(edge_sample_offsets[edge_index + 1])
        if begin == end:
            raise InvalidMedialAxisProjectionError("MAT edge has no proposal sample.")
        for sample_index in range(begin, end):
            parameter_id = parameter_ids[sample_index]
            samples.append(
                MatSample(
                    parameter_id,
                    _cursor_identity(edge.identity, parameter_id),
                    edge.identity,
                    Point2[WorldXY].build(
                        sample_centers[sample_index, 0],
                        sample_centers[sample_index, 1],
                    ),
                    Clearance.build(sample_clearance[sample_index]),
                    Clearance.build(sample_guide_radius[sample_index]),
                    MatParameter(float(sample_parameters[sample_index])),
                    sample_index - begin,
                )
            )

    components, component_by_edge_id = _connected_components(nodes, frozen_edges)
    tool_fit_runs = tuple(
        MatToolFitRun(
            edge.identity,
            edge.source.identity,
            edge.target.identity,
            edge.source_endpoint,
            edge.target_endpoint,
        )
        for edge in frozen_edges
    )
    topology = MatTopology(
        sites,
        nodes,
        frozen_edges,
        components,
        tool_fit_runs,
        MappingProxyType({site.identity: site for site in sites}),
        MappingProxyType({node.identity: node for node in nodes}),
        edge_by_id,
        component_by_edge_id,
    )
    proposal_sampling = MatProposalSampling(tuple(samples))
    proof = MatProof(
        _identity(projection[18], "MAT center-domain digest"),
        _identity(projection[19], "MAT certificate"),
        native_owner,
    )
    return topology, proposal_sampling, proof
