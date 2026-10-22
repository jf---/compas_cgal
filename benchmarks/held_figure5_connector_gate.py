"""Connector acceptance inside the refinement: tier 1, tier 1b and the exact convex-corner refusal predicates."""

from __future__ import annotations

from collections.abc import Iterator
from collections.abc import Sequence
from fractions import Fraction

import numpy as np

from benchmarks.held_connector_engagement import ConnectorEngagementLemma
from benchmarks.held_connector_engagement import LemmaCircle
from benchmarks.held_connector_engagement import _bisect
from benchmarks.held_connector_engagement import _CandidateIndex
from benchmarks.held_connector_engagement import _DyadicGrid
from benchmarks.held_connector_engagement import foot_squared_distance
from benchmarks.held_connector_engagement import prior_circle_indices
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_trace import EngagementFootprint
from benchmarks.held_figure5_engagement_trace import PieceBounds
from benchmarks.held_figure5_engagement_trace import floor_reach
from benchmarks.held_figure5_engagement_trace import rational_point
from benchmarks.held_figure5_engagement_trace import tier_1b_kept
from benchmarks.held_figure5_engagement_trace import tier_1b_reach
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


def exterior_angle_exceeds_cap(boundary: tuple[Point2[WorldXY], ...], vertex: int, cos_sq_half_cap: Fraction) -> bool:
    """Exact: ``vertex`` is convex (CCW left turn) and its exterior angle beta exceeds the cap.

    With incoming side e0, outgoing side e1, ``D = e0 . e1`` and
    ``L = |e0|^2 |e1|^2``, ``cos beta = D / sqrt(L)``. The cap enters as its
    exact surrogate ``C = cos^2(cap / 2) = 1 - s / 4`` (``s`` the binary64
    chord surrogate, injected exactly), so ``cos(cap)`` is ``m = 2C - 1``
    and, both angles lying in [0, pi], beta > cap iff ``D < m sqrt(L)``:
    decided by signs and ``D^2`` against ``m^2 L``, all rational.
    """
    before, at, after = (boundary[(vertex + k) % len(boundary)] for k in (-1, 0, 1))
    e0x, e0y = Fraction(float(at.x)) - Fraction(float(before.x)), Fraction(float(at.y)) - Fraction(float(before.y))
    e1x, e1y = Fraction(float(after.x)) - Fraction(float(at.x)), Fraction(float(after.y)) - Fraction(float(at.y))
    if e0x * e1y - e0y * e1x <= 0:
        return False
    dot = e0x * e1x + e0y * e1y
    lengths = (e0x * e0x + e0y * e0y) * (e1x * e1x + e1y * e1y)
    m = 2 * cos_sq_half_cap - 1
    if m > 0:
        return dot <= 0 or dot * dot < m * m * lengths
    if m == 0:
        return dot < 0
    return dot < 0 and dot * dot > m * m * lengths


def _outward_normal(boundary: tuple[Point2[WorldXY], ...], side: int) -> tuple[Fraction, Fraction]:
    """The CCW loop's outward side normal, unnormalised: the edge direction turned clockwise."""
    a, b = rational_point(boundary[side]), rational_point(boundary[(side + 1) % len(boundary)])
    return b[1] - a[1], a[0] - b[0]


def _in_wedge(point: tuple[Fraction, Fraction], apex: tuple[Fraction, Fraction], normals: tuple[tuple[Fraction, Fraction], ...]) -> bool:
    dx, dy = point[0] - apex[0], point[1] - apex[1]
    return all(dx * nx + dy * ny <= 0 for nx, ny in normals)


def _crossed_vertices(sides: int, start: Figure5CounterclockwiseCircle, end: Figure5CounterclockwiseCircle) -> Iterator[int]:
    """The boundary vertices the connector start -> end passes or reaches, in CCW order (as ``_transition`` walks them)."""
    if start.boundary_site == end.boundary_site and start is not end:
        return
    start_side, end_side = int(start.boundary_site.segment_id), int(end.boundary_site.segment_id)
    start_key, end_key = (start_side, start.boundary_site.parameter), (end_side, end.boundary_site.parameter)
    if start is end or end_key < start_key or start_side != end_side:
        side = start_side
        while True:
            side = (side + 1) % sides
            yield side
            if side == end_side:
                return


class _ConnectorGate:
    """Refinement acceptance for one connector: tier 1, else tier 1b, exactly as the post-gate decides it.

    ``ConnectorEngagementCertificate`` decides connector ``i`` from circles
    ``0..i``, the connector's samples, the tool and the cap only; every test is
    exact and homogeneous, so a grid over this connector's inputs decides as
    the post-gate's grid over the whole path does (both keep
    ``MAX_BISECTION_DEPTH`` spare bits, so the bisection points are the same
    rationals, and the candidate filter sees the same binary64 values). An
    accepted connector is therefore certified on the final path, whose circles
    ``0..i`` are emitted before it and never change. Tier 1b puts only the
    prior circles near the connector on the grid (``_nearby``), a superset of
    every circle the candidate filter can return, so it decides as over all.
    """

    def __init__(self, boundary: tuple[Point2[WorldXY], ...], tool_radius: ToolRadius, cap: EngagementCap) -> None:
        self._boundary = boundary
        self._lengths = _side_lengths(boundary)
        self._tool = tool_radius
        self._cos_sq_half_cap = ConnectorEngagementLemma.build(tool_radius, cap).cos_sq_half_cap
        # Per vertex, decided on first query: a splice continuation crosses few vertices.
        self._sharp: dict[int, bool] = {}
        self._centres: list[tuple[float, float]] = []
        self._reach: list[float] = []
        self.footprints: list[EngagementFootprint] = []

    def transition(self, predecessor: Figure5CounterclockwiseCircle, successor: Figure5CounterclockwiseCircle) -> Figure5BoundaryTransition:
        """The connector ``Figure5BoundaryConnectors.build`` emits between these circles."""
        return _transition(self._boundary, self._lengths, predecessor, successor)

    def sharp_vertex(self, emitted: Sequence[Figure5CounterclockwiseCircle], successor: Figure5CounterclockwiseCircle) -> int | None:
        """The convex vertex over the cap whose floor the connector ``emitted[-1] -> successor`` cannot escape, if any.

        Exact, all on rationals or integers: the connector crosses exactly one
        vertex V; V is convex with exterior angle beta over the cap
        (``exterior_angle_exceeds_cap``); its contacts lie on V's incoming and
        outgoing sides; every emitted guide disk within ``rho + 2r`` of V lies
        in V's tangent wedge T, and so does every non-adjacent polygon edge
        within ``2r`` of V (``_floor_hypothesis``).

        Why that suffices (docs/held_connector_engagement.md, "Convex
        corners"). For y in T, ``|V + r e - y|^2 >= r^2 + |y - V|^2`` for every
        e in V's normal cone N, so a cut whose guide points lie in T never
        reaches a rim point ``V + r e``, e in int N, which only points within
        2r of V can. The edge condition gives ``P ∩ B(V, 2r) ⊂ T`` exactly
        (B(V, 2r) minus T is connected, no boundary edge crosses it, and P is
        locally T at the convex vertex V). Circles a continued refinement could
        still insert cannot matter: on V's outgoing side they are cut after the
        arrival; on its incoming side, a guide circle inside the centre domain
        P has its points within 2r of V in T. A continuation the refusal
        forgoes would therefore need a guide circle crossing the centre-domain
        boundary within 2r of V, a gouge. No float is read here, and none of
        those future circles enters the decision. Emitted circles are checked
        as they are: a guide disk leaving T (the tooth fixture's R = 1.8
        circle, which crosses side 2) can pre-cut the corner, so that
        connector is left to the gate.
        """
        predecessor = emitted[-1]
        crossed = tuple(_crossed_vertices(len(self._boundary), predecessor, successor))
        if len(crossed) != 1 or not self._is_sharp(crossed[0]):
            return None
        vertex = crossed[0]
        incoming = (vertex - 1) % len(self._boundary)
        if int(predecessor.boundary_site.segment_id) != incoming or int(successor.boundary_site.segment_id) != vertex:
            return None
        self.footprints.append(EngagementFootprint(len(emitted), vertex=vertex))
        return vertex if self._floor_hypothesis(vertex, emitted) else None

    def _is_sharp(self, vertex: int) -> bool:
        if vertex not in self._sharp:
            self._sharp[vertex] = exterior_angle_exceeds_cap(self._boundary, vertex, self._cos_sq_half_cap)
        return self._sharp[vertex]

    def _floor_hypothesis(self, vertex: int, emitted: Sequence[Figure5CounterclockwiseCircle]) -> bool:
        """Exact: every emitted guide disk within ``rho + 2r`` of V, and every non-adjacent edge within ``2r``, lies in T."""
        sides = len(self._boundary)
        at = rational_point(self._boundary[vertex])
        normals = tuple(_outward_normal(self._boundary, side) for side in ((vertex - 1) % sides, vertex))
        reach = floor_reach(self._tool)
        for side in range(sides):
            if side in ((vertex - 1) % sides, vertex):
                continue
            a, b = rational_point(self._boundary[side]), rational_point(self._boundary[(side + 1) % sides])
            if foot_squared_distance(a, b, at) <= reach * reach and not (_in_wedge(a, at, normals) and _in_wedge(b, at, normals)):
                return False
        for circle in emitted:
            centre, rho = rational_point(circle.center), Fraction(float(circle.radius.value))
            dx, dy = centre[0] - at[0], centre[1] - at[1]
            if dx * dx + dy * dy > (rho + reach) * (rho + reach):
                continue
            for nx, ny in normals:
                signed = dx * nx + dy * ny
                if signed > 0 or signed * signed < rho * rho * (nx * nx + ny * ny):
                    return False
        return True

    def certifies(self, emitted: Sequence[Figure5CounterclockwiseCircle], transition: Figure5BoundaryTransition) -> bool:
        """Exact: the connector leaving ``emitted[-1]`` is certified by tier 1 (``emitted[-1]``) or tier 1b (any of ``emitted``)."""
        predecessor = emitted[-1]
        grid = _DyadicGrid(self._inputs((predecessor,), transition))
        pieces = self._pieces(grid, transition)
        if not pieces:
            return True
        own = LemmaCircle(grid.point(predecessor.center), grid.length(float(predecessor.radius.value)))
        lemma = ConnectorEngagementLemma(grid.length(float(self._tool.value)), self._cos_sq_half_cap)
        if all(lemma.piece_certified_by(a, b, own) for a, b in pieces):
            return True
        circles = self._nearby(emitted, transition)
        grid = _DyadicGrid(self._inputs(circles, transition))
        lemma = ConnectorEngagementLemma(grid.length(float(self._tool.value)), self._cos_sq_half_cap)
        exact = tuple(LemmaCircle(grid.point(circle.center), grid.length(float(circle.radius.value))) for circle in circles)
        index = _CandidateIndex(circles, self._tool, grid)
        prior = prior_circle_indices(len(circles) - 1, len(circles))
        return all(
            lemma.piece_certified_by(a, b, exact[-1]) or _bisect(lemma, exact, index, prior, piece, a, b) is not None
            for piece, (a, b) in enumerate(self._pieces(grid, transition))
        )

    def _nearby(self, emitted: Sequence[Figure5CounterclockwiseCircle], transition: Figure5BoundaryTransition) -> tuple[Figure5CounterclockwiseCircle, ...]:
        """The emitted circles any sub-piece of ``transition`` could ask about, in emission order, ``emitted[-1]`` last.

        ``emitted`` only appends during one refinement, so the binary64 centres
        and reaches are cached across calls. The decision's piece bounds are
        recorded as its footprint (``EngagementTrace``).
        """
        for circle in emitted[len(self._reach) :]:
            self._centres.append((float(circle.center.x), float(circle.center.y)))
            self._reach.append(tier_1b_reach(circle, self._tool))
        count = len(emitted)
        centres = np.array(self._centres[:count], dtype=np.float64).reshape(-1, 2)
        reach = np.array(self._reach[:count], dtype=np.float64)
        pieces = PieceBounds.of([(float(point.x), float(point.y)) for point in transition.samples])
        self.footprints.append(EngagementFootprint(count, pieces=pieces))
        near = tier_1b_kept(centres, reach, pieces).any(axis=1)
        near[count - 1] = True
        return tuple(emitted[int(k)] for k in np.nonzero(near)[0])

    def _inputs(self, circles: Sequence[Figure5CounterclockwiseCircle], transition: Figure5BoundaryTransition) -> Iterator[float]:
        yield float(self._tool.value)
        for circle in circles:
            yield from (float(circle.center.x), float(circle.center.y), float(circle.radius.value))
        for point in transition.samples:
            yield from (float(point.x), float(point.y))

    @staticmethod
    def _pieces(grid: _DyadicGrid, transition: Figure5BoundaryTransition) -> list[tuple[tuple[int, int], tuple[int, int]]]:
        vertices = [grid.point(point) for point in transition.samples]
        return [(a, b) for a, b in zip(vertices, vertices[1:]) if a != b]
