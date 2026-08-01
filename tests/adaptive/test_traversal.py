"""Contracts for exact identity-driven global MAT traversal."""

import hashlib
import math
import operator
from dataclasses import replace
from fractions import Fraction

import pytest

import compas_cgal.adaptive.traversal as traversal_module
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.candidates import ExhaustedCandidateCursor
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.candidates import enumerate_zero_guide_link_candidates
from compas_cgal.adaptive.errors import InvalidCausalNeckTransitError
from compas_cgal.adaptive.errors import InvalidMatTraversalStateError
from compas_cgal.adaptive.errors import InvalidTraversalGraphError
from compas_cgal.adaptive.errors import AmbiguousNeckSideError
from compas_cgal.adaptive.errors import NonterminalMatTraversalError
from compas_cgal.adaptive.errors import OverlappingNeckTransitError
from compas_cgal.adaptive.errors import TerminalTraversalCursorError
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.neck import NeckInventory
from compas_cgal.adaptive.neck import NeckSide
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NeckTraversalOrientation
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.adaptive.traversal import CausalNeckTransit
from compas_cgal.adaptive.traversal import DirectedEdgeCursor
from compas_cgal.adaptive.traversal import MatTraversalState
from compas_cgal.adaptive.traversal import NeckSidePosition
from compas_cgal.adaptive.traversal import VisitedEdgeIncidence
from compas_cgal.adaptive.traversal import resolve_neck_side
from compas_cgal.adaptive.traversal import _resolve_next_neck_state
from compas_cgal.adaptive.traversal_graph import TraversalGraph
from compas_cgal.adaptive.traversal_graph import TraversalGraphEdge

L_SHAPE = (
    (0.0, 0.0),
    (6.0, 0.0),
    (6.0, 2.0),
    (2.0, 2.0),
    (2.0, 6.0),
    (0.0, 6.0),
)


def _digest(label: str) -> bytes:
    return hashlib.sha256(label.encode("ascii")).digest()


def _edge(
    *,
    component: str,
    edge: str,
    source: str,
    target: str,
) -> TraversalGraphEdge:
    return TraversalGraphEdge.build(
        component_id=ComponentId(_digest(component)),
        edge_id=_digest(edge),
        branch_id=BranchId(_digest(f"branch:{edge}")),
        source_node_id=_digest(source),
        target_node_id=_digest(target),
        source_cursor_id=_digest(f"cursor:{edge}:source"),
        target_cursor_id=_digest(f"cursor:{edge}:target"),
    )


def _policy() -> TraversalPolicy:
    return TraversalPolicy.build(forward_window=4)


def _axis(
    *,
    tool_radius: float = 0.5,
    station_spacing: float = 0.75,
    max_refinement_depth: int = 32,
) -> MedialAxis:
    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y in L_SHAPE),
    )
    return MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(tool_radius),
        station_spacing=Spacing.build(station_spacing),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=max_refinement_depth,
    )


def _neck_policy() -> NeckPolicy:
    user_cap = EngagementCap.build(math.radians(120.0))
    return NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
        effective_caps={
            (neck_class, passage_state): EngagementCap.build(
                math.radians(90.0 - 10.0 * passage_state.rank),
            )
            for neck_class in range(2)
            for passage_state in ACTIVE_PASSAGE_STATES
        },
    )


def _candidate_policy() -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.5),
        spatial_refinement_levels=2,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=2,
        phase_count=4,
        minimum_guide_radius=GuideRadius.build(0.125),
        minimum_progress=Spacing.build(0.25),
    )


def _state() -> MatTraversalState:
    axis = _axis()
    inventory = NeckInventory.build(
        axis=axis,
        policy=_neck_policy(),
    )
    entry_edge = axis.edges[0]
    return MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=_policy(),
        entry_edge_id=entry_edge.identity,
        entry_node_id=entry_edge.source.identity,
    )


def _terminal_candidate(state: MatTraversalState):
    active = state.active_cursor
    span = MiddleCurveSpan.build(
        axis=state.authority.axis,
        cursor_before=active.cursor,
        cursor_limit=active.terminal_cursor,
    )
    user_cap = state.authority.inventory.policy.user_cap
    full_cap = FullCapDecision.build(
        user_cap=user_cap,
        effective_cap=user_cap,
    )
    return next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=_candidate_policy(),
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=state.neck_scope,
            effective_cap_decision=full_cap,
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.spatial_progress == span.reported_length
    )


def test_exact_route_orders_chain_and_branch_from_shared_node_ids() -> None:
    """Discover each branch from authenticated incidence, never coordinates.

    The root edge reaches one degree-three node. Both remaining edges must
    start at that exact shared node, and their relative order must follow the
    canonical branch identities regardless of input record order.
    """
    root = _edge(
        component="component:a",
        edge="edge:root",
        source="node:entry",
        target="node:branch",
    )
    left = _edge(
        component="component:a",
        edge="edge:left",
        source="node:branch",
        target="node:left",
    )
    right = _edge(
        component="component:a",
        edge="edge:right",
        source="node:branch",
        target="node:right",
    )
    graph = TraversalGraph.build(
        authority_digest=_digest("mat:branch"),
        edges=(right, root, left),
    )

    route = graph.route(
        policy=_policy(),
        entry_edge_id=root.edge_id,
        entry_node_id=root.source_node_id,
    )
    repeated = TraversalGraph.build(
        authority_digest=_digest("mat:branch"),
        edges=(left, right, root),
    ).route(
        policy=_policy(),
        entry_edge_id=root.edge_id,
        entry_node_id=root.source_node_id,
    )

    assert route == repeated
    assert route[0].edge_id == root.edge_id
    assert route[0].entry_node_id == root.source_node_id
    assert route[0].exit_node_id == root.target_node_id
    assert {step.edge_id for step in route} == {root.edge_id, left.edge_id, right.edge_id}
    assert all(step.entry_node_id == root.target_node_id for step in route[1:])
    assert tuple(step.branch_id for step in route[1:]) == _policy().order_branches(
        (left.branch_id, right.branch_id),
    )
    assert len({step.canonical_bytes for step in route}) == 3


def test_exact_route_retains_cycle_incidences_and_entry_component_first() -> None:
    """Visit a cycle edge once and preserve both of its directed incidences.

    A second disconnected component is deliberately assigned a lexically
    earlier identity. Entry authority still places the seeded component first;
    only remaining components use policy order.
    """
    cycle_component = "component:z-cycle"
    first = _edge(
        component=cycle_component,
        edge="edge:cycle-0",
        source="node:cycle-0",
        target="node:cycle-1",
    )
    second = _edge(
        component=cycle_component,
        edge="edge:cycle-1",
        source="node:cycle-1",
        target="node:cycle-2",
    )
    closing = _edge(
        component=cycle_component,
        edge="edge:cycle-2",
        source="node:cycle-2",
        target="node:cycle-0",
    )
    isolated = _edge(
        component="component:a-isolated",
        edge="edge:isolated",
        source="node:isolated-0",
        target="node:isolated-1",
    )
    graph = TraversalGraph.build(
        authority_digest=_digest("mat:cycle-and-component"),
        edges=(isolated, closing, second, first),
    )

    route = graph.route(
        policy=_policy(),
        entry_edge_id=first.edge_id,
        entry_node_id=first.source_node_id,
    )

    assert tuple(step.component_id for step in route[:3]) == (first.component_id,) * 3
    assert route[-1].component_id == isolated.component_id
    assert {step.edge_id for step in route[:3]} == {first.edge_id, second.edge_id, closing.edge_id}
    assert {(step.edge_id, step.entry_node_id) for step in route[:3]} == {
        (first.edge_id, first.source_node_id),
        (second.edge_id, second.source_node_id),
        (closing.edge_id, closing.source_node_id),
    }
    incidences = {(step.edge_id, node_id) for step in route[:3] for node_id in (step.entry_node_id, step.exit_node_id)}
    assert len(incidences) == 6


def test_exact_route_rejects_coordinate_plausible_but_disconnected_node_identity() -> None:
    """Reject a relabelled chain even if reporting points could coincide.

    Traversal records intentionally contain no coordinates. Replacing the
    shared node identity therefore disconnects the declared component and
    must fail instead of being repaired by geometric proximity.
    """
    first = _edge(
        component="component:chain",
        edge="edge:chain-0",
        source="node:chain-0",
        target="node:shared",
    )
    relabelled = _edge(
        component="component:chain",
        edge="edge:chain-1",
        source="node:coordinate-equal-but-foreign",
        target="node:chain-2",
    )

    with pytest.raises(InvalidTraversalGraphError, match="connected"):
        TraversalGraph.build(
            authority_digest=_digest("mat:broken-chain"),
            edges=(first, relabelled),
        )


def test_causal_transit_binds_one_of_three_exact_sides_and_passage_state() -> None:
    """Carry an exact three-way partition identity into oriented neck policy.

    The production plateau neck is not binary. A transit must bind two owned
    side records, derive orientation from their stable identities, and select
    exactly the matching independent passage state.
    """
    axis = _axis()
    inventory = NeckInventory.build(
        axis=axis,
        policy=_neck_policy(),
    )
    neck = inventory.necks[0]
    assert len(neck.sides) == 3
    source_side = neck.sides[0]
    target_side = neck.sides[1]

    transit = CausalNeckTransit.build(
        neck=neck,
        source_side=source_side,
        target_side=target_side,
    )
    expected_orientation = NeckTraversalOrientation.FORWARD if source_side.identity < target_side.identity else NeckTraversalOrientation.REVERSE
    passages = tuple(passage for classified in inventory.necks for passage in (classified.forward, classified.reverse))

    assert transit.orientation is expected_orientation
    assert transit.scope == transit.passage(passages).scope
    assert transit.passage(passages).state.rank == 0
    assert (
        transit.canonical_bytes
        == CausalNeckTransit.build(
            neck=neck,
            source_side=source_side,
            target_side=target_side,
        ).canonical_bytes
    )

    relabelled = NeckSide.build(
        neck_evidence_digest=neck.evidence_digest,
        partition_ordinal=target_side.partition_ordinal,
        edge_ids=source_side.edge_ids,
    )
    with pytest.raises(InvalidCausalNeckTransitError, match="owned side"):
        CausalNeckTransit.build(
            neck=neck,
            source_side=source_side,
            target_side=relabelled,
        )


def test_global_state_seeds_one_cursor_per_edge_with_exact_owner_identity() -> None:
    """Seed immutable graph state without copying physical generation state.

    The MAT certificate, neck inventory, deterministic route, one cursor per
    edge, entry incidence, and initial exact side of every production neck
    must all contribute to the state identity.
    """
    state = _state()
    repeated = MatTraversalState.seed(
        axis=state.authority.axis,
        inventory=state.authority.inventory,
        policy=state.authority.policy,
        entry_edge_id=state.authority.route[0].edge_id,
        entry_node_id=state.authority.route[0].entry_node_id,
    )

    assert state.canonical_bytes == repeated.canonical_bytes
    assert state.digest == repeated.digest
    assert len(state.cursors) == len(state.authority.axis.edges)
    assert tuple(cursor.route_step for cursor in state.cursors) == state.authority.route
    assert state.active_route_index == 0
    assert state.active_cursor.cursor_identity == state.authority.route[0].initial_cursor_id
    assert len(state.visited_incidences) == 1
    assert all(position.side is not None for position in state.neck_positions)
    assert type(state.neck_scope) is NoNeckScope
    with pytest.raises(NonterminalMatTraversalError):
        state.require_terminal()

    foreign_axis = _axis()
    foreign_inventory = NeckInventory.build(
        axis=foreign_axis,
        policy=_neck_policy(),
    )
    with pytest.raises(InvalidMatTraversalStateError, match="exact MAT owner"):
        MatTraversalState.seed(
            axis=state.authority.axis,
            inventory=foreign_inventory,
            policy=_policy(),
            entry_edge_id=state.authority.route[0].edge_id,
            entry_node_id=state.authority.route[0].entry_node_id,
        )


def test_seed_builds_each_cold_authority_and_cursor_index_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Prevent duplicate full builds before measuring candidate throughput.

    Graph projection, causal-frontier construction, and the complete native
    cursor index have graph lifetime. Seeding may invoke each exactly once;
    an immutable child must reuse that authority rather than rematch every
    edge against every native sample.
    """
    axis = _axis()
    inventory = NeckInventory.build(
        axis=axis,
        policy=_neck_policy(),
    )
    entry_edge = axis.edges[0]
    projection_calls = 0
    frontier_calls = 0
    neck_index_calls = 0
    cursor_index_calls = 0
    original_from_axis = TraversalGraph.from_axis
    original_frontiers = traversal_module._build_neck_frontiers
    original_neck_indexes = traversal_module._build_neck_lookups
    original_cursor_build = DirectedEdgeCursor.build_all

    def counted_from_axis(cls, *, axis, sample_index):
        nonlocal projection_calls
        projection_calls += 1
        assert cls is TraversalGraph
        return original_from_axis(
            axis=axis,
            sample_index=sample_index,
        )

    def counted_frontiers(*, axis, inventory, route):
        nonlocal frontier_calls
        frontier_calls += 1
        return original_frontiers(
            axis=axis,
            inventory=inventory,
            route=route,
        )

    def counted_neck_indexes(*, axis, inventory):
        nonlocal neck_index_calls
        neck_index_calls += 1
        return original_neck_indexes(
            axis=axis,
            inventory=inventory,
        )

    def counted_cursor_build(cls, *, sample_index, route):
        nonlocal cursor_index_calls
        cursor_index_calls += 1
        assert cls is DirectedEdgeCursor
        return original_cursor_build(
            sample_index=sample_index,
            route=route,
        )

    monkeypatch.setattr(
        TraversalGraph,
        "from_axis",
        classmethod(counted_from_axis),
    )
    monkeypatch.setattr(
        traversal_module,
        "_build_neck_frontiers",
        counted_frontiers,
    )
    monkeypatch.setattr(
        traversal_module,
        "_build_neck_lookups",
        counted_neck_indexes,
    )
    monkeypatch.setattr(
        DirectedEdgeCursor,
        "build_all",
        classmethod(counted_cursor_build),
    )

    state = MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=_policy(),
        entry_edge_id=entry_edge.identity,
        entry_node_id=entry_edge.source.identity,
    )
    candidate = _terminal_candidate(state)
    child = state.advance(candidate)

    assert (
        projection_calls,
        frontier_calls,
        neck_index_calls,
        cursor_index_calls,
    ) == (1, 1, 1, 1)
    assert child.authority is state.authority


def test_authority_binds_refinement_dependent_cursor_inventory() -> None:
    """Separate refinement-invariant MAT proof from traversal identity.

    Finer proposal sampling must preserve the exact topology certificate while
    changing the complete internal cursor inventory and traversal authority.
    The retained lookup views must also remain immutable.
    """
    coarse_axis = _axis(station_spacing=0.75)
    fine_axis = _axis(station_spacing=0.25)
    same_output_policy_axis = _axis(
        station_spacing=0.75,
        max_refinement_depth=33,
    )
    coarse_entry = coarse_axis.edges[0]
    fine_entry = fine_axis.edges[0]
    coarse = MatTraversalState.seed(
        axis=coarse_axis,
        inventory=NeckInventory.build(
            axis=coarse_axis,
            policy=_neck_policy(),
        ),
        policy=_policy(),
        entry_edge_id=coarse_entry.identity,
        entry_node_id=coarse_entry.source.identity,
    )
    fine = MatTraversalState.seed(
        axis=fine_axis,
        inventory=NeckInventory.build(
            axis=fine_axis,
            policy=_neck_policy(),
        ),
        policy=_policy(),
        entry_edge_id=fine_entry.identity,
        entry_node_id=fine_entry.source.identity,
    )
    same_output_entry = same_output_policy_axis.edges[0]
    same_output_policy = MatTraversalState.seed(
        axis=same_output_policy_axis,
        inventory=NeckInventory.build(
            axis=same_output_policy_axis,
            policy=_neck_policy(),
        ),
        policy=_policy(),
        entry_edge_id=same_output_entry.identity,
        entry_node_id=same_output_entry.source.identity,
    )

    assert coarse_axis.mat_certificate == fine_axis.mat_certificate
    assert coarse_axis.mat_certificate == same_output_policy_axis.mat_certificate
    assert coarse.authority.graph.canonical_bytes == fine.authority.graph.canonical_bytes
    assert len(coarse_axis.samples) < len(fine_axis.samples)
    assert coarse.authority.sample_index.canonical_bytes != fine.authority.sample_index.canonical_bytes
    assert coarse.authority.sample_index.canonical_bytes == same_output_policy.authority.sample_index.canonical_bytes
    assert coarse.authority.digest != fine.authority.digest
    assert coarse.authority.digest != same_output_policy.authority.digest
    with pytest.raises(TypeError):
        operator.setitem(
            coarse.authority.sample_index.sample_by_cursor_id,
            coarse.active_cursor.cursor_identity,
            coarse.active_cursor.cursor,
        )


def test_candidate_advance_mutates_one_cursor_then_activates_next_route_edge() -> None:
    """Commit one terminal candidate without aliasing any other graph cursor.

    A repeated candidate must fail on the terminal authoritative cursor, and
    branch activation must be impossible before the current edge terminates.
    """
    state = _state()
    candidate = _terminal_candidate(state)

    with pytest.raises(NonterminalMatTraversalError, match="active cursor"):
        state.activate_next()

    child = state.advance(candidate)
    changed = tuple(
        index
        for index, (before, after) in enumerate(
            zip(state.cursors, child.cursors, strict=True),
        )
        if before.canonical_bytes != after.canonical_bytes
    )

    assert changed == (0,)
    assert state.active_cursor.terminal is False
    assert child.active_cursor.terminal is True
    assert child.active_cursor.accepted_candidate_count == 1
    assert len(child.visited_incidences) == 2
    with pytest.raises(TerminalTraversalCursorError):
        child.advance(candidate)

    next_state = child.activate_next()
    assert next_state.active_route_index == 1
    assert next_state.active_cursor.route_step == next_state.authority.route[1]
    assert len(next_state.visited_incidences) == 3
    assert child.active_route_index == 0


def test_zero_guide_candidate_advances_the_closed_traversal_union() -> None:
    """A proved terminal link reaches its native endpoint without a circle."""
    axis = _axis(tool_radius=1.0)
    edge = axis.edge_by_id[axis.zero_guide_inventory.runs[0].edge_id]
    inventory = NeckInventory.build(
        axis=axis,
        policy=_neck_policy(),
    )
    state = MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=_policy(),
        entry_edge_id=edge.identity,
        entry_node_id=edge.source.identity,
    )
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=state.active_cursor.cursor,
        cursor_limit=state.active_cursor.terminal_cursor,
    )
    candidate = enumerate_zero_guide_link_candidates(
        span=span,
        policy=_candidate_policy(),
        neck_scope=state.neck_scope,
        effective_cap_decision=FullCapDecision.build(
            user_cap=inventory.policy.user_cap,
            effective_cap=inventory.policy.user_cap,
        ),
        makes_cursor_terminal_at_limit=True,
    )[0]

    child = state.advance(candidate)

    assert type(candidate) is ZeroGuideLinkCandidate
    assert candidate.spatial_progress == span.reported_length
    assert child.active_cursor.terminal
    assert child.active_cursor.cursor == state.active_cursor.terminal_cursor
    assert child.active_cursor.accepted_candidate_count == 1


def test_reverse_route_edge_advances_to_its_exact_entry_endpoint() -> None:
    """Consume the first route edge when graph discovery enters at its target.

    Seeding the real L-pocket parabola from its exact target node reverses the
    native sample order. The global ledger must accept a positive-progress
    terminal candidate and land on ordinal zero without changing route
    causality.
    """
    axis = _axis()
    inventory = NeckInventory.build(
        axis=axis,
        policy=_neck_policy(),
    )
    entry_edge = axis.edges[0]
    state = MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=_policy(),
        entry_edge_id=entry_edge.identity,
        entry_node_id=entry_edge.target.identity,
    )
    active = state.active_cursor
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=active.cursor,
        cursor_limit=active.terminal_cursor,
    )
    candidate = _terminal_candidate(state)
    child = state.advance(candidate)

    assert span.ordinal_step == -1
    assert active.cursor.ordinal_on_edge > active.terminal_cursor.ordinal_on_edge
    assert candidate.spatial_progress > 0
    assert child.active_cursor.terminal
    assert child.active_cursor.cursor == active.terminal_cursor


def test_clearance_leaf_retains_proof_carrying_exhausted_cursor() -> None:
    """Terminalize a tool-radius leaf without forging its native endpoint.

    The adopted L pocket contains a line leaf whose exact MAT endpoint has
    clearance equal to the tool radius. No positive-radius MATHSM circle can
    reach that point, so exhaustive enumeration marks its greatest feasible
    station terminal. Global traversal must retain that derived station while
    keeping the distinct native endpoint as the immutable route bound.
    """
    axis = _axis()
    inventory = NeckInventory.build(
        axis=axis,
        policy=_neck_policy(),
    )
    edge = next(edge for edge in axis.edges if edge.source.point == Point2[WorldXY].build(5.0, 1.0) and edge.target.point == Point2[WorldXY].build(5.5, 1.5))
    state = MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=_policy(),
        entry_edge_id=edge.identity,
        entry_node_id=edge.source.identity,
    )
    active = state.active_cursor
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=active.cursor,
        cursor_limit=active.terminal_cursor,
    )
    candidate = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=_candidate_policy(),
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=state.neck_scope,
            effective_cap_decision=FullCapDecision.build(
                user_cap=inventory.policy.user_cap,
                effective_cap=inventory.policy.user_cap,
            ),
            makes_cursor_terminal_at_limit=True,
        )
        if candidate.traversal_decision.makes_cursor_terminal
    )

    child = state.advance(candidate)

    assert candidate.spatial_progress < span.reported_length
    assert child.active_cursor.terminal
    assert type(child.active_cursor.cursor) is ExhaustedCandidateCursor
    assert child.active_cursor.cursor.candidate == candidate
    assert child.active_cursor.cursor.cursor_identity == candidate.traversal_decision.cursor_after
    assert child.active_cursor.cursor.point == candidate.middle_point
    assert child.active_cursor.terminal_cursor == active.terminal_cursor
    assert child.active_cursor.cursor != child.active_cursor.terminal_cursor
    assert child.activate_next().active_route_index == 1


def test_intermediate_native_limit_remains_a_continuable_sample() -> None:
    """Retain exact native identity when a candidate reaches a window limit.

    An intermediate forward-window endpoint is not an interior derived cursor
    and is not the route terminal. Its accepted candidate must advance onto
    the owned native sample so the next span can continue from that exact
    algebraic parameter without rematching coordinates.
    """
    axis = _axis()
    edge = next(edge for edge in axis.edges if len(tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)) >= 3)
    samples = tuple(sample for sample in axis.samples if sample.edge_id == edge.identity)
    state = MatTraversalState.seed(
        axis=axis,
        inventory=NeckInventory.build(
            axis=axis,
            policy=_neck_policy(),
        ),
        policy=_policy(),
        entry_edge_id=edge.identity,
        entry_node_id=edge.source.identity,
    )
    active = state.active_cursor
    intermediate = samples[1]
    span = MiddleCurveSpan.build(
        axis=axis,
        cursor_before=active.cursor,
        cursor_limit=intermediate,
    )
    candidate = next(
        candidate
        for candidate in enumerate_middle_curve_candidates(
            span=span,
            policy=_candidate_policy(),
            circle_orientation=CircleOrientation.COUNTERCLOCKWISE,
            neck_scope=state.neck_scope,
            effective_cap_decision=FullCapDecision.build(
                user_cap=state.authority.inventory.policy.user_cap,
                effective_cap=state.authority.inventory.policy.user_cap,
            ),
            makes_cursor_terminal_at_limit=False,
        )
        if candidate.spatial_progress == span.reported_length
    )

    child = state.advance(candidate)

    assert not child.active_cursor.terminal
    assert child.active_cursor.accepted_candidate_count == 1
    assert type(child.active_cursor.cursor) is MatSample
    assert child.active_cursor.cursor == intermediate
    assert (
        MiddleCurveSpan.build(
            axis=axis,
            cursor_before=child.active_cursor.cursor,
            cursor_limit=samples[2],
        ).reported_length
        > 0
    )


def test_route_activation_retains_source_side_until_scoped_candidate_commit() -> None:
    """Carry a plateau side through its locus before exposing one transit.

    The deterministic L route reaches the second plateau locus and then its
    one-edge leaf. Branch activation may expose the causal scope, but the
    authoritative side must remain the source until a candidate commits.
    """
    state = _state()
    for _ in range(3):
        assert type(state.neck_scope) is NoNeckScope
        state = state.advance(_terminal_candidate(state)).activate_next()

    transit = state.pending_transit
    assert transit is not None
    assert state.active_route_index == 3
    assert state.neck_scope == transit.scope
    assert state.neck_positions[1].side == transit.source_side
    assert state.active_cursor.route_step.edge_id in transit.target_side.edge_ids
    assert state.active_cursor.accepted_candidate_count == 0


def test_state_rejects_forged_side_frontier_and_accepts_complete_terminal_snapshot() -> None:
    """Validate causal side state against the precomputed exact route frontier.

    Replacing one owned side with another is internally well-typed but
    historically false and must fail. A structurally complete all-edge
    terminal snapshot must satisfy `require_terminal`.
    """
    state = _state()
    neck = state.authority.inventory.necks[0]
    forged_positions = (
        replace(state.neck_positions[0], side=neck.sides[1]),
        *state.neck_positions[1:],
    )
    with pytest.raises(InvalidMatTraversalStateError, match="causal route frontier"):
        replace(
            state,
            neck_positions=forged_positions,
        )

    terminal_cursors = tuple(
        DirectedEdgeCursor(
            cursor.route_step,
            cursor.terminal_cursor,
            cursor.terminal_cursor,
            1,
            True,
        )
        for cursor in state.cursors
    )
    visited = tuple(
        VisitedEdgeIncidence.build(
            edge_id=step.edge_id,
            node_id=node_id,
        )
        for step in state.authority.route
        for node_id in (step.entry_node_id, step.exit_node_id)
    )
    terminal = MatTraversalState(
        state.authority,
        terminal_cursors,
        None,
        visited,
        state.authority.neck_frontiers[-1].positions_after,
        None,
    )

    terminal.require_terminal()
    assert terminal.active_route_index is None
    assert len(terminal.visited_incidences) == 2 * len(terminal.cursors)


def test_side_resolution_rejects_same_union_overlap_and_multiple_active_necks() -> None:
    """Reject topology that a flattened cut union cannot distinguish.

    First, add one central edge to a leaf partition while preserving the exact
    union; side lookup becomes ambiguous. Then present a well-typed synthetic
    side frontier whose next production edge changes both necks; one operation
    cannot silently choose or merge those restrictions.
    """
    state = _state()
    first_neck, second_neck = state.authority.inventory.necks
    duplicated_edge = first_neck.sides[0].edge_ids[0]
    overlapping_side = NeckSide.build(
        neck_evidence_digest=first_neck.evidence_digest,
        partition_ordinal=first_neck.sides[1].partition_ordinal,
        edge_ids=tuple(
            sorted(
                (
                    duplicated_edge,
                    *first_neck.sides[1].edge_ids,
                )
            )
        ),
    )
    same_union_neck = replace(
        first_neck,
        sides=(
            first_neck.sides[0],
            overlapping_side,
            first_neck.sides[2],
        ),
    )
    with pytest.raises(AmbiguousNeckSideError, match="unique"):
        resolve_neck_side(
            neck=same_union_neck,
            edge_id=duplicated_edge,
        )

    target_step = next(step for step in state.authority.route if step.edge_id in first_neck.sides[0].edge_ids and step.edge_id in second_neck.sides[2].edge_ids)
    synthetic_positions = (
        NeckSidePosition.build(
            neck=first_neck,
            side=first_neck.sides[1],
        ),
        NeckSidePosition.build(
            neck=second_neck,
            side=second_neck.sides[0],
        ),
    )
    with pytest.raises(OverlappingNeckTransitError):
        _resolve_next_neck_state(
            axis=state.authority.axis,
            inventory=state.authority.inventory,
            positions=synthetic_positions,
            target_step=target_step,
        )
