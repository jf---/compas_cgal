"""Immutable global MAT traversal ledger and causal neck state."""

import hashlib
from dataclasses import dataclass
from dataclasses import field
from types import MappingProxyType
from typing import Mapping
from typing import Self
from typing import TypeAlias

from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import AmbiguousNeckSideError
from compas_cgal.adaptive.errors import InvalidCausalNeckTransitError
from compas_cgal.adaptive.errors import InvalidMatTraversalStateError
from compas_cgal.adaptive.errors import NonterminalMatTraversalError
from compas_cgal.adaptive.errors import OverlappingNeckTransitError
from compas_cgal.adaptive.errors import StaleTraversalCursorError
from compas_cgal.adaptive.errors import TerminalTraversalCursorError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatNodeId
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.neck import ClassifiedNeck
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.neck_topology import NeckSide
from compas_cgal.adaptive.neck_topology import neck_locus_edge_ids
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import NeckScope
from compas_cgal.adaptive.operation import NeckTraversalOrientation
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.traversal_graph import DirectedRouteStep
from compas_cgal.adaptive.traversal_graph import TraversalGraph
from compas_cgal.adaptive.traversal_graph import TraversalSampleIndex


def _state_identity(
    value: object,
    name: str,
) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidMatTraversalStateError(f"{name} must be nonempty exact bytes.")
    return value


@dataclass(frozen=True)
class CausalNeckTransit:
    """One exact transition between two certified sides of one neck."""

    neck: ClassifiedNeck
    source_side: NeckSide
    target_side: NeckSide
    orientation: NeckTraversalOrientation

    def __post_init__(self) -> None:
        if type(self.neck) is not ClassifiedNeck:
            raise InvalidCausalNeckTransitError("causal transit requires one exact classified neck.")
        if type(self.source_side) is not NeckSide or type(self.target_side) is not NeckSide:
            raise InvalidCausalNeckTransitError("causal transit requires two exact neck sides.")
        if self.source_side not in self.neck.sides or self.target_side not in self.neck.sides:
            raise InvalidCausalNeckTransitError("causal transit references a non-owned side.")
        if self.source_side == self.target_side:
            raise InvalidCausalNeckTransitError("causal transit requires distinct source and target sides.")
        expected = NeckTraversalOrientation.FORWARD if self.source_side.identity < self.target_side.identity else NeckTraversalOrientation.REVERSE
        if self.orientation is not expected:
            raise InvalidCausalNeckTransitError("causal transit orientation contradicts exact side order.")

    @classmethod
    def build(
        cls,
        *,
        neck: ClassifiedNeck,
        source_side: NeckSide,
        target_side: NeckSide,
    ) -> Self:
        """Build an oriented transition between two owned sides.

        Args:
            neck: Exact classified neck owning both sides.
            source_side: Certified side occupied before the crossing.
            target_side: Certified side occupied after the crossing.

        Returns:
            Immutable causal transition.
        """
        if type(source_side) is not NeckSide or type(target_side) is not NeckSide:
            raise InvalidCausalNeckTransitError("causal transit factory requires two exact neck sides.")
        orientation = NeckTraversalOrientation.FORWARD if source_side.identity < target_side.identity else NeckTraversalOrientation.REVERSE
        return cls(neck, source_side, target_side, orientation)

    @property
    def scope(self) -> OrientedNeckScope:
        return OrientedNeckScope.build(
            neck_owner_id=self.neck.owner_id,
            orientation=self.orientation,
        )

    def passage(self, passages: tuple[NeckPassage, ...]) -> NeckPassage:
        """Select the unique current passage for this transit.

        Args:
            passages: Complete independent oriented-passage inventory.

        Returns:
            Passage matching this exact neck and orientation.

        Raises:
            InvalidCausalNeckTransitError: If passage state is missing,
                duplicated, foreign, or cross-wired.
        """
        if type(passages) is not tuple or any(type(passage) is not NeckPassage for passage in passages):
            raise InvalidCausalNeckTransitError("causal transit passage lookup requires one exact immutable tuple.")
        matches = tuple(passage for passage in passages if passage.neck == self.neck and passage.scope == self.scope)
        if len(matches) != 1:
            raise InvalidCausalNeckTransitError("causal transit requires one unique owned passage state.")
        return matches[0]

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"causal-neck-transit-v1",
            encode_component_map(
                {
                    b"neck-evidence-digest": encode_bytes(bytes(self.neck.evidence_digest)),
                    b"neck-owner-id": encode_bytes(bytes(self.neck.owner_id)),
                    b"orientation": encode_tagged_union(bytes(self.orientation.value), b""),
                    b"source-side": self.source_side.canonical_bytes,
                    b"target-side": self.target_side.canonical_bytes,
                }
            ),
        )


@dataclass(frozen=True)
class NeckSidePosition:
    """Current certified side, or an exact not-yet-entered component state."""

    neck_owner_id: bytes
    side: NeckSide | None

    def __post_init__(self) -> None:
        _state_identity(self.neck_owner_id, "neck-position owner identity")
        if self.side is not None and type(self.side) is not NeckSide:
            raise InvalidMatTraversalStateError("neck position must contain one exact side or remain unentered.")

    @classmethod
    def build(
        cls,
        *,
        neck: ClassifiedNeck,
        side: NeckSide | None,
    ) -> Self:
        if type(neck) is not ClassifiedNeck:
            raise InvalidMatTraversalStateError("neck position requires one exact classified neck.")
        if side is not None and side not in neck.sides:
            raise InvalidMatTraversalStateError("neck position references a side not owned by its neck.")
        return cls(neck.owner_id, side)

    @property
    def canonical_bytes(self) -> bytes:
        side = encode_tagged_union(b"unentered-neck-component-v1", b"") if self.side is None else self.side.canonical_bytes
        return encode_tagged_union(
            b"neck-side-position-v1",
            encode_component_map(
                {
                    b"neck-owner-id": encode_bytes(self.neck_owner_id),
                    b"side": side,
                }
            ),
        )


@dataclass(frozen=True)
class VisitedEdgeIncidence:
    """One explicitly visited endpoint incidence of a directed MAT edge."""

    edge_id: MatEdgeId
    node_id: MatNodeId

    def __post_init__(self) -> None:
        _state_identity(self.edge_id, "visited-incidence edge identity")
        _state_identity(self.node_id, "visited-incidence node identity")

    @classmethod
    def build(
        cls,
        *,
        edge_id: MatEdgeId,
        node_id: MatNodeId,
    ) -> Self:
        return cls(edge_id, node_id)

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"visited-mat-edge-incidence-v1",
            encode_component_map(
                {
                    b"edge-id": encode_bytes(bytes(self.edge_id)),
                    b"node-id": encode_bytes(bytes(self.node_id)),
                }
            ),
        )


TraversalCursor: TypeAlias = MatSample | DerivedCandidateCursor


@dataclass(frozen=True)
class DirectedEdgeCursor:
    """Evolving accepted cursor on one statically directed route edge."""

    route_step: DirectedRouteStep
    cursor: TraversalCursor
    terminal_cursor: MatSample
    accepted_candidate_count: int
    terminal: bool

    def __post_init__(self) -> None:
        if type(self.route_step) is not DirectedRouteStep:
            raise InvalidMatTraversalStateError("directed cursor requires one exact route step.")
        if type(self.cursor) not in (MatSample, DerivedCandidateCursor):
            raise InvalidMatTraversalStateError("directed cursor requires one native or proof-carrying cursor.")
        if type(self.terminal_cursor) is not MatSample:
            raise InvalidMatTraversalStateError("directed cursor terminal bound must be one native sample.")
        if self.cursor.edge_id != self.route_step.edge_id or self.terminal_cursor.edge_id != self.route_step.edge_id:
            raise InvalidMatTraversalStateError("directed cursor records must remain on their exact edge.")
        if self.terminal_cursor.cursor_identity != self.route_step.terminal_cursor_id:
            raise InvalidMatTraversalStateError("directed cursor terminal sample contradicts its route.")
        if type(self.accepted_candidate_count) is not int or self.accepted_candidate_count < 0:
            raise InvalidMatTraversalStateError("directed cursor accepted count must be non-negative.")
        if type(self.terminal) is not bool:
            raise InvalidMatTraversalStateError("directed cursor terminal state must be explicitly boolean.")
        if self.accepted_candidate_count == 0:
            if type(self.cursor) is not MatSample or self.cursor.cursor_identity != self.route_step.initial_cursor_id or self.terminal:
                raise InvalidMatTraversalStateError("untouched directed cursor must remain at its native route entry.")
        elif self.terminal:
            if type(self.cursor) is not MatSample or self.cursor.cursor_identity != self.route_step.terminal_cursor_id:
                raise InvalidMatTraversalStateError("terminal directed cursor must retain its native route limit.")
        elif type(self.cursor) is not DerivedCandidateCursor:
            raise InvalidMatTraversalStateError("advanced nonterminal cursor must retain proof-carrying candidate lineage.")

    @classmethod
    def build(
        cls,
        *,
        axis: MedialAxis,
        route_step: DirectedRouteStep,
    ) -> Self:
        """Seed one route cursor from native endpoint samples.

        Args:
            axis: Exact MAT owner of the route.
            route_step: Static directed edge authority.

        Returns:
            Untouched nonterminal cursor at the directed entry.
        """
        sample_index = TraversalSampleIndex.build(axis)
        return cls.build_all(
            sample_index=sample_index,
            route=(route_step,),
        )[0]

    @classmethod
    def build_all(
        cls,
        *,
        sample_index: TraversalSampleIndex,
        route: tuple[DirectedRouteStep, ...],
    ) -> tuple[Self, ...]:
        """Seed every route cursor through one immutable native-sample index.

        Args:
            sample_index: Graph-lifetime native sample lookup authority.
            route: Complete deterministic directed edge route.

        Returns:
            One untouched cursor per route edge.

        Raises:
            InvalidMatTraversalStateError: If a route endpoint cursor is
                absent, foreign, or cross-wired to another edge.
        """
        if type(sample_index) is not TraversalSampleIndex:
            raise InvalidMatTraversalStateError("cursor batch requires one exact traversal sample index.")
        if type(route) is not tuple or not route or any(type(step) is not DirectedRouteStep for step in route):
            raise InvalidMatTraversalStateError("cursor batch requires one nonempty immutable route.")
        cursors: list[Self] = []
        for route_step in route:
            initial = sample_index.sample_by_cursor_id.get(route_step.initial_cursor_id)
            terminal = sample_index.sample_by_cursor_id.get(route_step.terminal_cursor_id)
            if initial is None or terminal is None or initial.edge_id != route_step.edge_id or terminal.edge_id != route_step.edge_id:
                raise InvalidMatTraversalStateError("route endpoint cursors must resolve uniquely in their exact MAT owner.")
            cursors.append(cls(route_step, initial, terminal, 0, False))
        return tuple(cursors)

    @property
    def cursor_identity(self) -> CursorIdentity:
        return self.cursor.cursor_identity

    def _advance(
        self,
        *,
        axis: MedialAxis,
        sample_index: TraversalSampleIndex,
        candidate: MiddleCurveCandidate,
    ) -> Self:
        if self.terminal:
            raise TerminalTraversalCursorError("candidate cannot advance a terminal global MAT cursor.")
        if type(candidate) is not MiddleCurveCandidate:
            raise StaleTraversalCursorError("global cursor advance requires one exact middle-curve candidate.")
        decision = candidate.traversal_decision
        if (
            decision.component_id != self.route_step.component_id
            or bytes(decision.edge_id) != bytes(self.route_step.edge_id)
            or decision.branch_id != self.route_step.branch_id
            or decision.cursor_before != self.cursor_identity
        ):
            raise StaleTraversalCursorError("candidate does not begin at the active global MAT cursor.")
        if type(sample_index) is not TraversalSampleIndex or sample_index.axis is not axis:
            raise StaleTraversalCursorError("candidate cursor index is not owned by the active exact MAT.")
        limit = sample_index.sample_by_cursor_id.get(candidate.cursor_limit_identity)
        if limit is None or limit.edge_id != self.route_step.edge_id:
            raise StaleTraversalCursorError("candidate limit is not one owned native cursor.")
        span = MiddleCurveSpan.build(
            axis=axis,
            cursor_before=self.cursor,
            cursor_limit=limit,
        )
        if decision.makes_cursor_terminal:
            if decision.cursor_after != self.route_step.terminal_cursor_id or limit != self.terminal_cursor or candidate.spatial_progress != span.reported_length:
                raise StaleTraversalCursorError("terminal candidate does not reach the directed edge limit.")
            cursor: TraversalCursor = self.terminal_cursor
        else:
            cursor = DerivedCandidateCursor.build(
                span=span,
                candidate=candidate,
            )
        return type(self)(
            self.route_step,
            cursor,
            self.terminal_cursor,
            self.accepted_candidate_count + 1,
            decision.makes_cursor_terminal,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"directed-edge-cursor-v1",
            encode_component_map(
                {
                    b"accepted-candidate-count": encode_integer(self.accepted_candidate_count),
                    b"cursor-id": encode_bytes(bytes(self.cursor_identity)),
                    b"route-step": self.route_step.canonical_bytes,
                    b"terminal": encode_boolean(self.terminal),
                    b"terminal-cursor-id": encode_bytes(bytes(self.terminal_cursor.cursor_identity)),
                }
            ),
        )


@dataclass(frozen=True)
class RouteNeckFrontier:
    """Precomputed causal side state before and after one route edge begins."""

    positions_before: tuple[NeckSidePosition, ...]
    positions_after: tuple[NeckSidePosition, ...]
    transit: CausalNeckTransit | None

    def __post_init__(self) -> None:
        if (
            type(self.positions_before) is not tuple
            or type(self.positions_after) is not tuple
            or any(type(position) is not NeckSidePosition for position in (*self.positions_before, *self.positions_after))
            or len(self.positions_before) != len(self.positions_after)
        ):
            raise InvalidMatTraversalStateError("route neck frontier requires two complete immutable position tuples.")
        if self.transit is not None and type(self.transit) is not CausalNeckTransit:
            raise InvalidMatTraversalStateError("route neck frontier transit must use the closed exact type.")

    @property
    def canonical_bytes(self) -> bytes:
        transit = encode_tagged_union(b"no-route-neck-transit-v1", b"") if self.transit is None else self.transit.canonical_bytes
        return encode_tagged_union(
            b"route-neck-frontier-v1",
            encode_component_map(
                {
                    b"positions-after": encode_sequence(tuple(position.canonical_bytes for position in self.positions_after)),
                    b"positions-before": encode_sequence(tuple(position.canonical_bytes for position in self.positions_before)),
                    b"transit": transit,
                }
            ),
        )


@dataclass(frozen=True)
class MatTraversalAuthority:
    """Static exact MAT, neck, route, and policy authority built once."""

    axis: MedialAxis
    inventory: NeckInventory
    policy: TraversalPolicy
    entry_edge_id: MatEdgeId
    entry_node_id: MatNodeId
    sample_index: TraversalSampleIndex = field(init=False, repr=False)
    graph: TraversalGraph = field(init=False)
    route: tuple[DirectedRouteStep, ...] = field(init=False)
    neck_frontiers: tuple[RouteNeckFrontier, ...] = field(init=False)
    _digest: IdentityDigest = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if type(self.axis) is not MedialAxis:
            raise InvalidMatTraversalStateError("traversal authority requires one exact typed MAT.")
        if type(self.inventory) is not NeckInventory or self.inventory.axis is not self.axis:
            raise InvalidMatTraversalStateError("traversal neck inventory must retain the same exact MAT owner.")
        if type(self.policy) is not TraversalPolicy:
            raise InvalidMatTraversalStateError("traversal authority requires one exact traversal policy.")
        _state_identity(self.entry_edge_id, "traversal entry edge identity")
        _state_identity(self.entry_node_id, "traversal entry node identity")
        sample_index = TraversalSampleIndex.build(self.axis)
        graph = TraversalGraph.from_axis(
            axis=self.axis,
            sample_index=sample_index,
        )
        route = graph.route(
            policy=self.policy,
            entry_edge_id=self.entry_edge_id,
            entry_node_id=self.entry_node_id,
        )
        neck_frontiers = _build_neck_frontiers(
            axis=self.axis,
            inventory=self.inventory,
            route=route,
        )
        object.__setattr__(self, "sample_index", sample_index)
        object.__setattr__(self, "graph", graph)
        object.__setattr__(self, "route", route)
        object.__setattr__(self, "neck_frontiers", neck_frontiers)
        object.__setattr__(
            self,
            "_digest",
            IdentityDigest(hashlib.sha256(self.canonical_bytes).digest()),
        )

    @classmethod
    def build(
        cls,
        *,
        axis: MedialAxis,
        inventory: NeckInventory,
        policy: TraversalPolicy,
        entry_edge_id: MatEdgeId,
        entry_node_id: MatNodeId,
    ) -> Self:
        return cls(
            axis,
            inventory,
            policy,
            entry_edge_id,
            entry_node_id,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"mat-traversal-authority-v1",
            encode_component_map(
                {
                    b"entry-edge-id": encode_bytes(bytes(self.entry_edge_id)),
                    b"entry-node-id": encode_bytes(bytes(self.entry_node_id)),
                    b"graph": self.graph.canonical_bytes,
                    b"neck-inventory": encode_sequence(tuple(neck.canonical_bytes for neck in self.inventory.necks)),
                    b"neck-frontiers": encode_sequence(tuple(frontier.canonical_bytes for frontier in self.neck_frontiers)),
                    b"route": encode_sequence(tuple(step.canonical_bytes for step in self.route)),
                    b"sampling-policy": encode_tagged_union(
                        b"mat-traversal-sampling-policy-v1",
                        encode_component_map(
                            {
                                b"max-refinement-depth": encode_integer(self.axis.max_refinement_depth),
                                b"max-sagitta": encode_binary64(float(self.axis.max_sagitta.value)),
                                b"station-spacing": encode_binary64(float(self.axis.station_spacing.value)),
                            }
                        ),
                    ),
                    b"sample-index": self.sample_index.canonical_bytes,
                    b"traversal-forward-window": encode_integer(self.policy.forward_window),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return self._digest


@dataclass(frozen=True)
class MatTraversalState:
    """Immutable graph-lifetime state independent of stock and coverage."""

    authority: MatTraversalAuthority
    cursors: tuple[DirectedEdgeCursor, ...]
    active_route_index: int | None
    visited_incidences: tuple[VisitedEdgeIncidence, ...]
    neck_positions: tuple[NeckSidePosition, ...]
    pending_transit: CausalNeckTransit | None
    _digest: IdentityDigest = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if type(self.authority) is not MatTraversalAuthority:
            raise InvalidMatTraversalStateError("global traversal state requires one static exact authority.")
        if type(self.cursors) is not tuple or len(self.cursors) != len(self.authority.route) or any(type(cursor) is not DirectedEdgeCursor for cursor in self.cursors):
            raise InvalidMatTraversalStateError("global traversal requires one immutable cursor per route edge.")
        if tuple(cursor.route_step for cursor in self.cursors) != self.authority.route:
            raise InvalidMatTraversalStateError("global cursors contradict canonical route order.")
        if self.active_route_index is not None and (type(self.active_route_index) is not int or self.active_route_index < 0 or self.active_route_index >= len(self.cursors)):
            raise InvalidMatTraversalStateError("active route index is outside the complete cursor ledger.")
        if type(self.visited_incidences) is not tuple or any(type(incidence) is not VisitedEdgeIncidence for incidence in self.visited_incidences):
            raise InvalidMatTraversalStateError("visited incidences must use one immutable exact tuple.")
        if len(set(self.visited_incidences)) != len(self.visited_incidences):
            raise InvalidMatTraversalStateError("visited edge incidences must remain unique.")
        if type(self.neck_positions) is not tuple or any(type(position) is not NeckSidePosition for position in self.neck_positions):
            raise InvalidMatTraversalStateError("neck positions must use one immutable exact tuple.")
        if tuple(position.neck_owner_id for position in self.neck_positions) != tuple(neck.owner_id for neck in self.authority.inventory.necks):
            raise InvalidMatTraversalStateError("neck positions contradict exact inventory order.")
        for neck, position in zip(
            self.authority.inventory.necks,
            self.neck_positions,
            strict=True,
        ):
            if position.side is not None and position.side not in neck.sides:
                raise InvalidMatTraversalStateError("neck position is not owned by its classified neck.")
        if self.pending_transit is not None and type(self.pending_transit) is not CausalNeckTransit:
            raise InvalidMatTraversalStateError("pending neck transit must use the closed exact type.")
        _validate_cursor_frontier(self.cursors, self.active_route_index)
        if self.visited_incidences != _expected_visited_incidences(
            self.cursors,
            self.active_route_index,
        ):
            raise InvalidMatTraversalStateError("visited incidences contradict the directed cursor frontier.")
        if self.active_route_index is None:
            expected_positions = self.authority.neck_frontiers[-1].positions_after
            expected_transit = None
        else:
            frontier = self.authority.neck_frontiers[self.active_route_index]
            if self.active_cursor.accepted_candidate_count == 0:
                expected_positions = frontier.positions_before
                expected_transit = frontier.transit
            else:
                expected_positions = frontier.positions_after
                expected_transit = None
        if self.neck_positions != expected_positions or self.pending_transit != expected_transit:
            raise InvalidMatTraversalStateError("neck positions or transit contradict the precomputed causal route frontier.")
        object.__setattr__(
            self,
            "_digest",
            IdentityDigest(hashlib.sha256(self.canonical_bytes).digest()),
        )

    @classmethod
    def seed(
        cls,
        *,
        axis: MedialAxis,
        inventory: NeckInventory,
        policy: TraversalPolicy,
        entry_edge_id: MatEdgeId,
        entry_node_id: MatNodeId,
    ) -> Self:
        """Seed complete global traversal at one authenticated MAT endpoint.

        Args:
            axis: Exact MAT and native proof owner.
            inventory: Exact neck inventory retaining that same owner.
            policy: Stable component and branch order.
            entry_edge_id: Exact first edge.
            entry_node_id: Exact endpoint fixing first direction.

        Returns:
            Immutable state with one untouched cursor per edge.
        """
        authority = MatTraversalAuthority.build(
            axis=axis,
            inventory=inventory,
            policy=policy,
            entry_edge_id=entry_edge_id,
            entry_node_id=entry_node_id,
        )
        cursors = DirectedEdgeCursor.build_all(
            sample_index=authority.sample_index,
            route=authority.route,
        )
        positions = authority.neck_frontiers[0].positions_before
        entry = VisitedEdgeIncidence.build(
            edge_id=authority.route[0].edge_id,
            node_id=authority.route[0].entry_node_id,
        )
        return cls(
            authority,
            cursors,
            0,
            (entry,),
            positions,
            None,
        )

    @property
    def active_cursor(self) -> DirectedEdgeCursor:
        if self.active_route_index is None:
            raise NonterminalMatTraversalError("terminal traversal has no active cursor.")
        return self.cursors[self.active_route_index]

    @property
    def neck_scope(self) -> NeckScope:
        if self.pending_transit is None:
            return NoNeckScope.build()
        return self.pending_transit.scope

    @property
    def canonical_bytes(self) -> bytes:
        active_index = -1 if self.active_route_index is None else self.active_route_index
        pending = encode_tagged_union(b"no-pending-neck-transit-v1", b"") if self.pending_transit is None else self.pending_transit.canonical_bytes
        return encode_tagged_union(
            b"mat-traversal-state-v1",
            encode_component_map(
                {
                    b"active-route-index": encode_integer(active_index),
                    b"authority-digest": encode_bytes(bytes(self.authority.digest)),
                    b"cursors": encode_sequence(tuple(cursor.canonical_bytes for cursor in self.cursors)),
                    b"neck-positions": encode_sequence(tuple(position.canonical_bytes for position in self.neck_positions)),
                    b"pending-transit": pending,
                    b"visited-incidences": encode_sequence(tuple(incidence.canonical_bytes for incidence in self.visited_incidences)),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return self._digest

    def advance(self, candidate: MiddleCurveCandidate) -> Self:
        """Advance exactly the active cursor with one accepted candidate.

        Args:
            candidate: Exact candidate beginning at the active cursor.

        Returns:
            Child traversal state; every other cursor remains byte-identical.
        """
        active_index = self.active_route_index
        if active_index is None:
            raise TerminalTraversalCursorError("terminal global traversal cannot accept another candidate.")
        cursor_after = self.active_cursor._advance(
            axis=self.authority.axis,
            sample_index=self.authority.sample_index,
            candidate=candidate,
        )
        if candidate.neck_scope != self.neck_scope:
            raise StaleTraversalCursorError("candidate neck scope contradicts global causal side history.")
        cursors = list(self.cursors)
        cursors[active_index] = cursor_after
        visited = self.visited_incidences
        if cursor_after.terminal:
            visited = (
                *visited,
                VisitedEdgeIncidence.build(
                    edge_id=cursor_after.route_step.edge_id,
                    node_id=cursor_after.route_step.exit_node_id,
                ),
            )
        positions = self.authority.neck_frontiers[active_index].positions_after
        return type(self)(
            self.authority,
            tuple(cursors),
            active_index,
            visited,
            positions,
            None,
        )

    def activate_next(self) -> Self:
        """Activate the next deterministic route edge after terminalization.

        Returns:
            State at the next edge entry, or terminal state after the last
            edge.

        Raises:
            NonterminalMatTraversalError: If the active edge is unfinished.
        """
        active_index = self.active_route_index
        if active_index is None:
            raise NonterminalMatTraversalError("terminal traversal has no next route edge.")
        if not self.active_cursor.terminal:
            raise NonterminalMatTraversalError("active cursor must be terminal before branch activation.")
        if self.pending_transit is not None:
            raise InvalidMatTraversalStateError("consumed active edge cannot retain a pending neck transit.")
        next_index = active_index + 1
        if next_index == len(self.cursors):
            return type(self)(
                self.authority,
                self.cursors,
                None,
                self.visited_incidences,
                self.neck_positions,
                None,
            )
        frontier = self.authority.neck_frontiers[next_index]
        visited = (
            *self.visited_incidences,
            VisitedEdgeIncidence.build(
                edge_id=self.authority.route[next_index].edge_id,
                node_id=self.authority.route[next_index].entry_node_id,
            ),
        )
        return type(self)(
            self.authority,
            self.cursors,
            next_index,
            visited,
            frontier.positions_before,
            frontier.transit,
        )

    def require_terminal(self) -> None:
        """Require every edge and incidence to be terminal.

        Raises:
            NonterminalMatTraversalError: If any route state remains active or
                incomplete.
        """
        if (
            self.active_route_index is not None
            or self.pending_transit is not None
            or any(not cursor.terminal for cursor in self.cursors)
            or self.visited_incidences != _expected_visited_incidences(self.cursors, None)
        ):
            raise NonterminalMatTraversalError("global MAT traversal is not terminal.")


def _neck_component_id(
    axis: MedialAxis,
    neck: ClassifiedNeck,
) -> ComponentId:
    edge_ids = {edge_id for side in neck.sides for edge_id in side.edge_ids}
    edge_ids.update(neck_locus_edge_ids(neck.locus))
    if not edge_ids or any(edge_id not in axis.component_by_edge_id for edge_id in edge_ids):
        raise AmbiguousNeckSideError("neck topology references an edge outside the exact traversal graph.")
    component_ids = {axis.component_by_edge_id[edge_id] for edge_id in edge_ids}
    if len(component_ids) != 1:
        raise AmbiguousNeckSideError("one neck topology cannot span multiple MAT components.")
    return next(iter(component_ids))


NeckEdgeLookup: TypeAlias = Mapping[MatEdgeId, NeckSide | None]
NeckRouteLookup: TypeAlias = tuple[ComponentId, NeckEdgeLookup]


def _neck_edge_lookup(neck: ClassifiedNeck) -> NeckEdgeLookup:
    if type(neck) is not ClassifiedNeck:
        raise AmbiguousNeckSideError("neck lookup requires one exact classified neck.")
    by_edge: dict[MatEdgeId, NeckSide | None] = {}
    for side in neck.sides:
        for edge_id in side.edge_ids:
            if edge_id in by_edge:
                raise AmbiguousNeckSideError("every route edge must resolve one unique certified neck side.")
            by_edge[edge_id] = side
    for edge_id in neck_locus_edge_ids(neck.locus):
        if edge_id in by_edge:
            raise AmbiguousNeckSideError("an exact neck locus edge cannot also belong to a certified side.")
        by_edge[edge_id] = None
    if not by_edge:
        raise AmbiguousNeckSideError("neck lookup cannot be empty.")
    return MappingProxyType(by_edge)


def _build_neck_lookups(
    *,
    axis: MedialAxis,
    inventory: NeckInventory,
) -> tuple[NeckRouteLookup, ...]:
    lookups: list[NeckRouteLookup] = []
    for neck in inventory.necks:
        component_id = _neck_component_id(axis, neck)
        by_edge = _neck_edge_lookup(neck)
        component_edge_ids = {edge_id for edge_id, owned_component_id in axis.component_by_edge_id.items() if owned_component_id == component_id}
        if set(by_edge) != component_edge_ids:
            raise AmbiguousNeckSideError("neck sides plus locus must partition its complete exact MAT component.")
        lookups.append((component_id, by_edge))
    return tuple(lookups)


def resolve_neck_side(
    *,
    neck: ClassifiedNeck,
    edge_id: MatEdgeId,
) -> NeckSide | None:
    """Resolve one route edge against an exact neck partition.

    Args:
        neck: Exact classified neck retaining disjoint certified sides.
        edge_id: Exact MAT edge identity to resolve.

    Returns:
        The unique owning side, or `None` when the edge belongs to the exact
        neck locus.

    Raises:
        AmbiguousNeckSideError: If the inputs are invalid or the edge is
            absent from, or duplicated across, the neck topology.
    """
    if type(neck) is not ClassifiedNeck or type(edge_id) is not bytes or not edge_id:
        raise AmbiguousNeckSideError("side resolution requires one exact neck and nonempty edge identity.")
    by_edge = _neck_edge_lookup(neck)
    if edge_id not in by_edge:
        raise AmbiguousNeckSideError("route edge must resolve one unique certified neck side or exact locus.")
    return by_edge[edge_id]


def _initial_neck_positions(
    *,
    inventory: NeckInventory,
    lookups: tuple[NeckRouteLookup, ...],
    entry_step: DirectedRouteStep,
) -> tuple[NeckSidePosition, ...]:
    positions: list[NeckSidePosition] = []
    for neck, (component_id, by_edge) in zip(
        inventory.necks,
        lookups,
        strict=True,
    ):
        if component_id != entry_step.component_id:
            positions.append(NeckSidePosition.build(neck=neck, side=None))
            continue
        side = by_edge[entry_step.edge_id]
        if side is None:
            raise AmbiguousNeckSideError("entry edge must resolve one unique exact side of every neck in its component.")
        positions.append(
            NeckSidePosition.build(
                neck=neck,
                side=side,
            )
        )
    return tuple(positions)


def _resolve_next_neck_state(
    *,
    axis: MedialAxis,
    inventory: NeckInventory,
    positions: tuple[NeckSidePosition, ...],
    target_step: DirectedRouteStep,
) -> tuple[tuple[NeckSidePosition, ...], CausalNeckTransit | None]:
    return _resolve_next_neck_state_indexed(
        inventory=inventory,
        lookups=_build_neck_lookups(
            axis=axis,
            inventory=inventory,
        ),
        positions=positions,
        target_step=target_step,
    )


def _resolve_next_neck_state_indexed(
    *,
    inventory: NeckInventory,
    lookups: tuple[NeckRouteLookup, ...],
    positions: tuple[NeckSidePosition, ...],
    target_step: DirectedRouteStep,
) -> tuple[tuple[NeckSidePosition, ...], CausalNeckTransit | None]:
    target_component = target_step.component_id
    updated = list(positions)
    transits: list[CausalNeckTransit] = []
    for index, (neck, position, lookup) in enumerate(
        zip(
            inventory.necks,
            positions,
            lookups,
            strict=True,
        )
    ):
        component_id, by_edge = lookup
        if component_id != target_component:
            continue
        target_side = by_edge[target_step.edge_id]
        if position.side is None:
            if target_side is None:
                raise AmbiguousNeckSideError("component entry edge must resolve one unique exact neck side.")
            updated[index] = NeckSidePosition.build(
                neck=neck,
                side=target_side,
            )
            continue
        if target_side is not None:
            if target_side != position.side:
                transits.append(
                    CausalNeckTransit.build(
                        neck=neck,
                        source_side=position.side,
                        target_side=target_side,
                    )
                )
            continue
        # Locus edges retain the last causal side until a certified target
        # side is entered; the transition is attributed to that first edge.
    if len(transits) > 1:
        raise OverlappingNeckTransitError("one route activation would cross multiple independently active necks.")
    return tuple(updated), (transits[0] if transits else None)


def _build_neck_frontiers(
    *,
    axis: MedialAxis,
    inventory: NeckInventory,
    route: tuple[DirectedRouteStep, ...],
) -> tuple[RouteNeckFrontier, ...]:
    lookups = _build_neck_lookups(
        axis=axis,
        inventory=inventory,
    )
    positions = _initial_neck_positions(
        inventory=inventory,
        lookups=lookups,
        entry_step=route[0],
    )
    frontiers = [
        RouteNeckFrontier(
            positions,
            positions,
            None,
        )
    ]
    for target_step in route[1:]:
        positions_before, transit = _resolve_next_neck_state_indexed(
            inventory=inventory,
            lookups=lookups,
            positions=positions,
            target_step=target_step,
        )
        positions_after = tuple(
            NeckSidePosition.build(
                neck=neck,
                side=(transit.target_side if transit is not None and neck == transit.neck else position.side),
            )
            for neck, position in zip(
                inventory.necks,
                positions_before,
                strict=True,
            )
        )
        frontiers.append(
            RouteNeckFrontier(
                positions_before,
                positions_after,
                transit,
            )
        )
        positions = positions_after
    return tuple(frontiers)


def _validate_cursor_frontier(
    cursors: tuple[DirectedEdgeCursor, ...],
    active_route_index: int | None,
) -> None:
    terminal_prefix = len(cursors) if active_route_index is None else active_route_index
    if any(not cursor.terminal for cursor in cursors[:terminal_prefix]):
        raise InvalidMatTraversalStateError("every cursor before the active route edge must be terminal.")
    if active_route_index is None:
        return
    if any(cursor.accepted_candidate_count != 0 or cursor.terminal for cursor in cursors[active_route_index + 1 :]):
        raise InvalidMatTraversalStateError("future route cursors must remain untouched.")


def _expected_visited_incidences(
    cursors: tuple[DirectedEdgeCursor, ...],
    active_route_index: int | None,
) -> tuple[VisitedEdgeIncidence, ...]:
    last_index = len(cursors) if active_route_index is None else active_route_index
    visited: list[VisitedEdgeIncidence] = []
    for cursor in cursors[:last_index]:
        visited.extend(
            (
                VisitedEdgeIncidence.build(
                    edge_id=cursor.route_step.edge_id,
                    node_id=cursor.route_step.entry_node_id,
                ),
                VisitedEdgeIncidence.build(
                    edge_id=cursor.route_step.edge_id,
                    node_id=cursor.route_step.exit_node_id,
                ),
            )
        )
    if active_route_index is not None:
        active = cursors[active_route_index]
        visited.append(
            VisitedEdgeIncidence.build(
                edge_id=active.route_step.edge_id,
                node_id=active.route_step.entry_node_id,
            )
        )
        if active.terminal:
            visited.append(
                VisitedEdgeIncidence.build(
                    edge_id=active.route_step.edge_id,
                    node_id=active.route_step.exit_node_id,
                )
            )
    return tuple(visited)
