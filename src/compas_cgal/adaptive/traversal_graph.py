"""Immutable exact MAT graph projection and deterministic route discovery."""

import hashlib
from dataclasses import dataclass
from dataclasses import field
from types import MappingProxyType
from typing import Mapping
from typing import Self

from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidTraversalGraphError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatNodeId
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import TraversalPolicy

_SHA256_BYTES = hashlib.sha256().digest_size


def _identity(
    value: object,
    name: str,
    *,
    digest: bool = False,
) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidTraversalGraphError(f"{name} must be nonempty exact bytes.")
    if digest and len(value) != _SHA256_BYTES:
        raise InvalidTraversalGraphError(f"{name} must be one SHA-256 digest.")
    return value


@dataclass(frozen=True)
class TraversalSampleIndex:
    """One immutable native-sample index shared by graph-lifetime consumers."""

    axis: MedialAxis
    samples_by_edge: Mapping[MatEdgeId, tuple[MatSample, ...]] = field(
        init=False,
        repr=False,
        compare=False,
    )
    sample_by_cursor_id: Mapping[CursorIdentity, MatSample] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if type(self.axis) is not MedialAxis:
            raise InvalidTraversalGraphError("traversal sample index requires one exact typed MAT.")
        mutable_by_edge: dict[MatEdgeId, list[MatSample]] = {edge.identity: [] for edge in self.axis.edges}
        sample_by_cursor_id: dict[CursorIdentity, MatSample] = {}
        for sample in self.axis.samples:
            edge_samples = mutable_by_edge.get(sample.edge_id)
            if edge_samples is None:
                raise InvalidTraversalGraphError("native traversal sample references an unknown exact edge.")
            if sample.cursor_identity in sample_by_cursor_id:
                raise InvalidTraversalGraphError("native traversal cursor identities must be globally unique.")
            edge_samples.append(sample)
            sample_by_cursor_id[sample.cursor_identity] = sample

        samples_by_edge: dict[MatEdgeId, tuple[MatSample, ...]] = {}
        for edge_id, edge_samples in mutable_by_edge.items():
            ordered = tuple(
                sorted(
                    edge_samples,
                    key=lambda sample: sample.ordinal_on_edge,
                )
            )
            if len(ordered) < 2 or tuple(sample.ordinal_on_edge for sample in ordered) != tuple(range(len(ordered))):
                raise InvalidTraversalGraphError("every traversal edge requires complete canonical native samples.")
            samples_by_edge[edge_id] = ordered

        object.__setattr__(
            self,
            "samples_by_edge",
            MappingProxyType(samples_by_edge),
        )
        object.__setattr__(
            self,
            "sample_by_cursor_id",
            MappingProxyType(sample_by_cursor_id),
        )

    @classmethod
    def build(cls, axis: MedialAxis) -> Self:
        """Build both native sample lookup views in one complete pass.

        Args:
            axis: Exact MAT owner of every sample.

        Returns:
            Immutable edge and cursor indexes bound to `axis`.
        """
        return cls(axis)

    @property
    def authority_digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.axis.mat_certificate).digest())

    @property
    def canonical_bytes(self) -> bytes:
        samples = tuple(
            encode_tagged_union(
                b"mat-traversal-sample-v1",
                encode_component_map(
                    {
                        b"cursor-id": encode_bytes(bytes(sample.cursor_identity)),
                        b"edge-id": encode_bytes(bytes(sample.edge_id)),
                        b"exact-parameter-id": encode_bytes(bytes(sample.exact_parameter_id)),
                        b"ordinal": encode_integer(sample.ordinal_on_edge),
                    }
                ),
            )
            for edge_id in sorted(self.samples_by_edge)
            for sample in self.samples_by_edge[edge_id]
        )
        return encode_tagged_union(
            b"mat-traversal-sample-index-v1",
            encode_component_map(
                {
                    b"authority-digest": encode_bytes(bytes(self.authority_digest)),
                    b"samples": encode_sequence(samples),
                }
            ),
        )


@dataclass(frozen=True)
class TraversalGraphEdge:
    """Static exact graph and endpoint-cursor authority for one MAT edge."""

    component_id: ComponentId
    edge_id: MatEdgeId
    branch_id: BranchId
    source_node_id: MatNodeId
    target_node_id: MatNodeId
    source_cursor_id: CursorIdentity
    target_cursor_id: CursorIdentity

    def __post_init__(self) -> None:
        _identity(self.component_id, "traversal component identity", digest=True)
        _identity(self.edge_id, "traversal edge identity")
        _identity(self.branch_id, "traversal branch identity", digest=True)
        _identity(self.source_node_id, "traversal source-node identity")
        _identity(self.target_node_id, "traversal target-node identity")
        _identity(self.source_cursor_id, "traversal source cursor", digest=True)
        _identity(self.target_cursor_id, "traversal target cursor", digest=True)
        if self.source_node_id == self.target_node_id:
            raise InvalidTraversalGraphError("traversal edge endpoints must be distinct.")
        if self.source_cursor_id == self.target_cursor_id:
            raise InvalidTraversalGraphError("traversal edge endpoint cursors must be distinct.")

    @classmethod
    def build(
        cls,
        *,
        component_id: ComponentId,
        edge_id: MatEdgeId,
        branch_id: BranchId,
        source_node_id: MatNodeId,
        target_node_id: MatNodeId,
        source_cursor_id: CursorIdentity,
        target_cursor_id: CursorIdentity,
    ) -> Self:
        """Build one static traversal edge.

        Args:
            component_id: Exact connected-component identity.
            edge_id: Exact MAT edge identity.
            branch_id: Stable branch-order identity.
            source_node_id: Canonical MAT source-node identity.
            target_node_id: Canonical MAT target-node identity.
            source_cursor_id: Native cursor at the source endpoint.
            target_cursor_id: Native cursor at the target endpoint.

        Returns:
            Validated immutable graph edge.
        """
        return cls(
            component_id,
            edge_id,
            branch_id,
            source_node_id,
            target_node_id,
            source_cursor_id,
            target_cursor_id,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"mat-traversal-graph-edge-v1",
            encode_component_map(
                {
                    b"branch-id": encode_bytes(bytes(self.branch_id)),
                    b"component-id": encode_bytes(bytes(self.component_id)),
                    b"edge-id": encode_bytes(bytes(self.edge_id)),
                    b"source-cursor-id": encode_bytes(bytes(self.source_cursor_id)),
                    b"source-node-id": encode_bytes(bytes(self.source_node_id)),
                    b"target-cursor-id": encode_bytes(bytes(self.target_cursor_id)),
                    b"target-node-id": encode_bytes(bytes(self.target_node_id)),
                }
            ),
        )


@dataclass(frozen=True)
class DirectedRouteStep:
    """One deterministic edge direction in the complete global route."""

    component_id: ComponentId
    edge_id: MatEdgeId
    branch_id: BranchId
    entry_node_id: MatNodeId
    exit_node_id: MatNodeId
    initial_cursor_id: CursorIdentity
    terminal_cursor_id: CursorIdentity

    def __post_init__(self) -> None:
        _identity(self.component_id, "route component identity", digest=True)
        _identity(self.edge_id, "route edge identity")
        _identity(self.branch_id, "route branch identity", digest=True)
        _identity(self.entry_node_id, "route entry-node identity")
        _identity(self.exit_node_id, "route exit-node identity")
        _identity(self.initial_cursor_id, "route initial cursor", digest=True)
        _identity(self.terminal_cursor_id, "route terminal cursor", digest=True)
        if self.entry_node_id == self.exit_node_id:
            raise InvalidTraversalGraphError("directed route step requires distinct endpoint identities.")
        if self.initial_cursor_id == self.terminal_cursor_id:
            raise InvalidTraversalGraphError("directed route step requires distinct endpoint cursors.")

    @classmethod
    def build(
        cls,
        *,
        edge: TraversalGraphEdge,
        entry_node_id: MatNodeId,
    ) -> Self:
        """Orient one exact graph edge from its discovered node.

        Args:
            edge: Static exact traversal edge.
            entry_node_id: Exact endpoint from which the route discovered it.

        Returns:
            Directed immutable route step.

        Raises:
            InvalidTraversalGraphError: If the discovery node is not an edge
                endpoint.
        """
        if type(edge) is not TraversalGraphEdge:
            raise InvalidTraversalGraphError("directed route requires one exact graph edge.")
        if entry_node_id == edge.source_node_id:
            return cls(
                edge.component_id,
                edge.edge_id,
                edge.branch_id,
                edge.source_node_id,
                edge.target_node_id,
                edge.source_cursor_id,
                edge.target_cursor_id,
            )
        if entry_node_id == edge.target_node_id:
            return cls(
                edge.component_id,
                edge.edge_id,
                edge.branch_id,
                edge.target_node_id,
                edge.source_node_id,
                edge.target_cursor_id,
                edge.source_cursor_id,
            )
        raise InvalidTraversalGraphError("route discovery node is not incident to its exact edge.")

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"directed-mat-route-step-v1",
            encode_component_map(
                {
                    b"branch-id": encode_bytes(bytes(self.branch_id)),
                    b"component-id": encode_bytes(bytes(self.component_id)),
                    b"edge-id": encode_bytes(bytes(self.edge_id)),
                    b"entry-node-id": encode_bytes(bytes(self.entry_node_id)),
                    b"exit-node-id": encode_bytes(bytes(self.exit_node_id)),
                    b"initial-cursor-id": encode_bytes(bytes(self.initial_cursor_id)),
                    b"terminal-cursor-id": encode_bytes(bytes(self.terminal_cursor_id)),
                }
            ),
        )


@dataclass(frozen=True)
class TraversalGraph:
    """Static content-addressed graph authority used by route discovery."""

    authority_digest: IdentityDigest
    edges: tuple[TraversalGraphEdge, ...]

    def __post_init__(self) -> None:
        _identity(self.authority_digest, "traversal graph authority", digest=True)
        if type(self.edges) is not tuple or not self.edges or any(type(edge) is not TraversalGraphEdge for edge in self.edges):
            raise InvalidTraversalGraphError("traversal graph requires one nonempty immutable edge tuple.")
        canonical = tuple(
            sorted(
                self.edges,
                key=lambda edge: (
                    edge.component_id,
                    edge.branch_id,
                    edge.edge_id,
                ),
            )
        )
        if self.edges != canonical:
            raise InvalidTraversalGraphError("traversal graph edges must use canonical component/branch order.")
        if len({edge.edge_id for edge in self.edges}) != len(self.edges):
            raise InvalidTraversalGraphError("traversal graph edge identities must be unique.")
        if len({edge.branch_id for edge in self.edges}) != len(self.edges):
            raise InvalidTraversalGraphError("traversal graph branch identities must be unique.")
        for component_edges in _edges_by_component(self.edges).values():
            if not _component_is_connected(component_edges):
                raise InvalidTraversalGraphError("every declared traversal component must be connected by exact node identities.")

    @classmethod
    def build(
        cls,
        *,
        authority_digest: IdentityDigest,
        edges: tuple[TraversalGraphEdge, ...],
    ) -> Self:
        """Build a canonical static graph authority.

        Args:
            authority_digest: SHA-256 identity of the owning exact MAT proof.
            edges: Complete exact edge records.

        Returns:
            Canonically ordered traversal graph.
        """
        if type(edges) is not tuple:
            raise InvalidTraversalGraphError("traversal graph factory requires one immutable edge tuple.")
        return cls(
            authority_digest,
            tuple(
                sorted(
                    edges,
                    key=lambda edge: (
                        edge.component_id,
                        edge.branch_id,
                        edge.edge_id,
                    ),
                )
            ),
        )

    @classmethod
    def from_axis(
        cls,
        *,
        axis: MedialAxis,
        sample_index: TraversalSampleIndex,
    ) -> Self:
        """Project the static route authority from one exact MAT owner.

        Args:
            axis: Typed MAT retaining the native proof owner.
            sample_index: Complete graph-lifetime native sample index.

        Returns:
            Complete traversal graph with native endpoint cursors.

        Raises:
            InvalidTraversalGraphError: If samples do not cover each exact
                edge in complete canonical order.
        """
        if type(axis) is not MedialAxis:
            raise InvalidTraversalGraphError("traversal graph projection requires one exact typed MAT.")
        if type(sample_index) is not TraversalSampleIndex or sample_index.axis is not axis:
            raise InvalidTraversalGraphError("traversal graph sample index must retain the same exact MAT owner.")
        edges: list[TraversalGraphEdge] = []
        for edge in axis.edges:
            samples = sample_index.samples_by_edge[edge.identity]
            edges.append(
                TraversalGraphEdge.build(
                    component_id=axis.component_by_edge_id[edge.identity],
                    edge_id=edge.identity,
                    branch_id=edge.branch_id,
                    source_node_id=edge.source.identity,
                    target_node_id=edge.target.identity,
                    source_cursor_id=samples[0].cursor_identity,
                    target_cursor_id=samples[-1].cursor_identity,
                )
            )
        return cls.build(
            authority_digest=IdentityDigest(hashlib.sha256(axis.mat_certificate).digest()),
            edges=tuple(edges),
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"mat-traversal-graph-v1",
            encode_component_map(
                {
                    b"authority-digest": encode_bytes(bytes(self.authority_digest)),
                    b"edges": encode_sequence(tuple(edge.canonical_bytes for edge in self.edges)),
                }
            ),
        )

    @property
    def identity(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())

    def route(
        self,
        *,
        policy: TraversalPolicy,
        entry_edge_id: MatEdgeId,
        entry_node_id: MatNodeId,
    ) -> tuple[DirectedRouteStep, ...]:
        """Discover one complete deterministic directed edge route.

        Args:
            policy: Stable component and branch ordering.
            entry_edge_id: Exact edge containing the authenticated entry.
            entry_node_id: Exact endpoint that fixes the first direction.

        Returns:
            One directed step for every graph edge.

        Raises:
            InvalidTraversalGraphError: If entry authority is foreign or route
                discovery fails to cover the complete graph.
        """
        if type(policy) is not TraversalPolicy:
            raise InvalidTraversalGraphError("route discovery requires one exact traversal policy.")
        edges_by_id = {edge.edge_id: edge for edge in self.edges}
        entry_edge = edges_by_id.get(entry_edge_id)
        if entry_edge is None:
            raise InvalidTraversalGraphError("route entry edge is not owned by the exact graph.")
        if entry_node_id not in (
            entry_edge.source_node_id,
            entry_edge.target_node_id,
        ):
            raise InvalidTraversalGraphError("route entry node is not incident to the entry edge.")

        edges_by_component = _edges_by_component(self.edges)
        component_ids = policy.order_components(tuple(edges_by_component))
        ordered_components = (
            entry_edge.component_id,
            *(component_id for component_id in component_ids if component_id != entry_edge.component_id),
        )
        route: list[DirectedRouteStep] = []
        for component_id in ordered_components:
            component_edges = edges_by_component[component_id]
            if component_id == entry_edge.component_id:
                root_edge = entry_edge
                root_node_id = entry_node_id
            else:
                root_edge = _ordered_incident_edges(
                    component_edges,
                    policy,
                )[0]
                root_node_id = min(
                    root_edge.source_node_id,
                    root_edge.target_node_id,
                )
            route.extend(
                _route_component(
                    component_edges,
                    policy=policy,
                    root_edge=root_edge,
                    root_node_id=root_node_id,
                )
            )
        if len(route) != len(self.edges) or {step.edge_id for step in route} != set(edges_by_id):
            raise InvalidTraversalGraphError("directed route does not cover every exact MAT edge once.")
        return tuple(route)


def _edges_by_component(
    edges: tuple[TraversalGraphEdge, ...],
) -> dict[ComponentId, tuple[TraversalGraphEdge, ...]]:
    mutable: dict[ComponentId, list[TraversalGraphEdge]] = {}
    for edge in edges:
        mutable.setdefault(edge.component_id, []).append(edge)
    return {component_id: tuple(component_edges) for component_id, component_edges in mutable.items()}


def _component_is_connected(edges: tuple[TraversalGraphEdge, ...]) -> bool:
    nodes = {node_id for edge in edges for node_id in (edge.source_node_id, edge.target_node_id)}
    incident: dict[MatNodeId, list[MatNodeId]] = {node_id: [] for node_id in nodes}
    for edge in edges:
        incident[edge.source_node_id].append(edge.target_node_id)
        incident[edge.target_node_id].append(edge.source_node_id)
    pending = [min(nodes)]
    visited: set[MatNodeId] = set()
    while pending:
        node_id = pending.pop()
        if node_id in visited:
            continue
        visited.add(node_id)
        pending.extend(other for other in incident[node_id] if other not in visited)
    return visited == nodes


def _ordered_incident_edges(
    edges: tuple[TraversalGraphEdge, ...],
    policy: TraversalPolicy,
) -> tuple[TraversalGraphEdge, ...]:
    by_branch = {edge.branch_id: edge for edge in edges}
    return tuple(by_branch[branch_id] for branch_id in policy.order_branches(tuple(by_branch)))


def _route_component(
    edges: tuple[TraversalGraphEdge, ...],
    *,
    policy: TraversalPolicy,
    root_edge: TraversalGraphEdge,
    root_node_id: MatNodeId,
) -> tuple[DirectedRouteStep, ...]:
    mutable_incident: dict[MatNodeId, list[TraversalGraphEdge]] = {}
    for edge in edges:
        mutable_incident.setdefault(edge.source_node_id, []).append(edge)
        mutable_incident.setdefault(edge.target_node_id, []).append(edge)
    incident = {
        node_id: _ordered_incident_edges(
            tuple(node_edges),
            policy,
        )
        for node_id, node_edges in mutable_incident.items()
    }
    node_ids = incident.keys()

    route = [
        DirectedRouteStep.build(
            edge=root_edge,
            entry_node_id=root_node_id,
        )
    ]
    visited = {root_edge.edge_id}
    positions = {node_id: 0 for node_id in node_ids}

    def explore(start_node_id: MatNodeId) -> None:
        stack = [start_node_id]
        while stack:
            node_id = stack[-1]
            node_edges = incident[node_id]
            position = positions[node_id]
            while position < len(node_edges) and node_edges[position].edge_id in visited:
                position += 1
            positions[node_id] = position
            if position == len(node_edges):
                stack.pop()
                continue
            edge = node_edges[position]
            positions[node_id] = position + 1
            step = DirectedRouteStep.build(
                edge=edge,
                entry_node_id=node_id,
            )
            route.append(step)
            visited.add(edge.edge_id)
            stack.append(step.exit_node_id)

    explore(route[0].exit_node_id)
    explore(route[0].entry_node_id)
    if len(visited) != len(edges):
        raise InvalidTraversalGraphError("component route discovery omitted an exact edge.")
    return tuple(route)
