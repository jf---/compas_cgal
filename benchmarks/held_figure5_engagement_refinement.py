"""Bounded engagement refinement proposals on the existing polygon path.

This is an approximate spacing experiment, not machining qualification or MAT
construction. Source circles keep their contacts and family attribution while
their normals and corner radii are repaired. Added circles interpolate source
intervals; original_indices locates repaired source stations. Containment must
be assessed independently.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace

from benchmarks.held_actual_stock_decision import ActualStockDecision
from benchmarks.held_actual_stock_decision import EntryCircle
from benchmarks.held_actual_stock_decision import read_disk
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _validate_boundary
from benchmarks.held_figure5_connector_gate import _ConnectorGate as _ConnectorGate
from benchmarks.held_figure5_connector_gate import exterior_angle_exceeds_cap as exterior_angle_exceeds_cap
from benchmarks.held_figure5_engagement_circles import circle_on_side as circle_on_side
from benchmarks.held_figure5_engagement_circles import concave_corner_bridges as concave_corner_bridges
from benchmarks.held_figure5_engagement_circles import midpoint as midpoint
from benchmarks.held_figure5_engagement_errors import ActualStockEngagementResolutionError as ActualStockEngagementResolutionError
from benchmarks.held_figure5_engagement_errors import ConnectorEngagementResolutionError as ConnectorEngagementResolutionError
from benchmarks.held_figure5_engagement_errors import ConvexCornerEngagementError as ConvexCornerEngagementError
from benchmarks.held_figure5_engagement_errors import Figure5EngagementCosineRangeError as Figure5EngagementCosineRangeError
from benchmarks.held_figure5_engagement_errors import Figure5EngagementFragmentationError as Figure5EngagementFragmentationError
from benchmarks.held_figure5_engagement_errors import Figure5EngagementPrefixError as Figure5EngagementPrefixError
from benchmarks.held_figure5_engagement_errors import Figure5EngagementResolutionError as Figure5EngagementResolutionError
from benchmarks.held_figure5_engagement_errors import Figure5EngagementTraceError as Figure5EngagementTraceError
from benchmarks.held_figure5_engagement_trace import DecisionSite
from benchmarks.held_figure5_engagement_trace import EngagementFootprint
from benchmarks.held_figure5_engagement_trace import EngagementRefinementKey
from benchmarks.held_figure5_engagement_trace import EngagementTrace
from benchmarks.held_figure5_engagement_trace import _tier_1b_neighbours as _tier_1b_neighbours
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_standard_placement import ExactEngagementCosineRangeError
from benchmarks.held_standard_placement import PredecessorEngagementCosine
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

MAX_REFINED_CIRCLES = 20000  # Explicit resource ceiling for the Figure 5 experiment.


@dataclass(frozen=True)
class Figure5EngagementRefinement:
    """Refined motion proposals with original guide stations distinguished."""

    circles: tuple[Figure5CounterclockwiseCircle, ...]
    original_indices: tuple[int, ...]


@dataclass(frozen=True)
class TracedEngagementRefinement(Figure5EngagementRefinement):
    """A refinement with the trace of its decisions: the type every refinement of this module returns.

    ``trace`` (not part of equality) records which emitted circles each
    decision could read (``EngagementTrace``). The coverage repair requires
    this type (a splice stops its continuation early from the trace), so a
    refinement assembled by hand cannot be repaired by accident.
    """

    trace: EngagementTrace = field(compare=False, repr=False)


def _transition_engagement(
    predecessor: Figure5CounterclockwiseCircle,
    target: Figure5CounterclockwiseCircle,
    tool_radius: ToolRadius,
    source_interval: int,
) -> PredecessorEngagementCosine:
    """Exact engagement of one emitted transition; placement-model failures become named refusals."""
    try:
        return exact_predecessor_engagement(_paper_candidate(predecessor), _paper_candidate(target), tool_radius)
    except StandardPlacementFragmentationError as error:
        raise Figure5EngagementFragmentationError(f"Fragmented transition at source interval {source_interval}: {error}") from error
    except ExactEngagementCosineRangeError as error:
        raise Figure5EngagementCosineRangeError(f"Engagement cosine out of range at source interval {source_interval}: {error}") from error


@dataclass
class _Emission:
    """The refinement's growing emission and its budgets; circles only ever append."""

    boundary: tuple[Point2[WorldXY], ...]
    tool_radius: ToolRadius
    cap: EngagementCap
    gate: _ConnectorGate
    stock: ActualStockDecision
    site: DecisionSite
    max_depth: int
    max_circles: int
    corner_approaches: bool
    emitted: list[Figure5CounterclockwiseCircle]

    def stock_refuses(self, target: Figure5CounterclockwiseCircle, *, seam: bool) -> bool:
        """The stock decision on ``target``, its footprint recorded whatever the verdict."""
        position = len(self.emitted)
        read = read_disk(target)
        refused = self.stock.exceeded(read, position)
        self.gate.footprints.append(EngagementFootprint(position, read=read, site=DecisionSite.CLOSE_SEAM if seam else self.site, refused=refused))
        return refused

    def advance(self, targets: tuple[Figure5CounterclockwiseCircle, ...], depth: int, source_interval: int, *, seam: bool = False) -> None:
        """Emit ``targets`` in order, inserting midpoint circles until each transition is accepted.

        A transition is accepted iff its circle pair is within the cap, its
        target is within the cap on its circle view (``ActualStockDecision``)
        and its connector is certified (tier 1 or 1b). Any refusal inserts the
        midpoint circle, depth + 1. A connector that is the first arrival at a
        convex vertex over the cap is refused first, by name
        (``require_no_sharp_corner``); ``seam`` marks seam repair circles.
        """
        pending = [(target, depth) for target in reversed(targets)]
        while pending:
            target, level = pending.pop()
            self.require_no_sharp_corner(target, seam=seam)
            predecessor = self.emitted[-1]
            engagement = _transition_engagement(predecessor, target, self.tool_radius, source_interval)
            exhausted: Figure5EngagementResolutionError | None = None
            if engagement.within(self.cap):
                if self.stock_refuses(target, seam=seam):
                    exhausted = ActualStockEngagementResolutionError(len(self.emitted), seam, f"source interval {source_interval}")
                elif self.gate.certifies(self.emitted, self.gate.transition(predecessor, target)):
                    if len(self.emitted) >= self.max_circles:
                        raise Figure5EngagementResolutionError(f"Circle budget exhausted at source interval {source_interval}.")
                    self.emitted.append(target)
                    self.stock.append(target)
                    continue
                else:
                    exhausted = ConnectorEngagementResolutionError(len(self.emitted) - 1, False, f"source interval {source_interval}")
                if level >= self.max_depth:
                    raise exhausted
            elif level >= self.max_depth:
                raise Figure5EngagementResolutionError(f"Depth exhausted at source interval {source_interval}: {math.degrees(float(engagement.reported_angle())):.6f} degrees.")
            try:
                halfway = midpoint(predecessor, target, self.boundary, corner_approaches=self.corner_approaches)
            except Figure5EngagementResolutionError as error:
                if exhausted is None:
                    raise
                # The repair has no circle left to insert (a vanishing corner-limited circle at a convex vertex).
                if isinstance(exhausted, ActualStockEngagementResolutionError):
                    raise ActualStockEngagementResolutionError(len(self.emitted), seam, f"source interval {source_interval}, repair midpoint refused: {error}") from error
                raise ConnectorEngagementResolutionError(len(self.emitted) - 1, False, f"source interval {source_interval}, repair midpoint refused: {error}") from error
            pending.append((target, level + 1))
            pending.append((halfway, level + 1))

    def require_no_sharp_corner(self, successor: Figure5CounterclockwiseCircle, *, seam: bool) -> None:
        """Refuse, before any midpoint, the first arrival at a convex vertex whose exterior angle exceeds the cap.

        ``_ConnectorGate.sharp_vertex`` states the exact conditions and why
        no in-domain continuation can certify the arrival: halving would only
        drive the midpoint circles toward zero radius at the vertex.
        """
        vertex = self.gate.sharp_vertex(self.emitted, successor)
        if vertex is not None:
            raise ConvexCornerEngagementError(len(self.emitted) - 1, seam, vertex)

    def close_seam(self, source_interval: int) -> None:
        """Certify the seam connector (last circle back to circles[0]), appending midpoint circles as new last circles.

        circles[0] is never re-cut and the closing pair (last, circles[0]) stays
        report-only (user rulings); a ``COINCIDENT`` seam has no connector.
        """
        rounds = 0
        while True:
            last, first = self.emitted[-1], self.emitted[0]
            if last is not first and last.contact_point == first.contact_point:
                return
            self.require_no_sharp_corner(first, seam=True)
            if self.gate.certifies(self.emitted, self.gate.transition(last, first)):
                return
            if rounds >= self.max_depth:
                raise ConnectorEngagementResolutionError(len(self.emitted) - 1, True, f"{rounds} seam rounds")
            rounds += 1
            try:
                halfway = midpoint(last, first, self.boundary, corner_approaches=self.corner_approaches)
            except Figure5EngagementResolutionError as error:
                raise ConnectorEngagementResolutionError(len(self.emitted) - 1, True, f"repair midpoint refused: {error}") from error
            self.advance((halfway,), rounds, source_interval, seam=True)


def refine_figure5_engagement(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    part: PartPolygon,
    max_depth: int = 16,
    max_circles: int = MAX_REFINED_CIRCLES,
    corner_approaches: bool = False,
    close_seam: bool = False,
) -> TracedEngagementRefinement:
    """Repair corner radii and resolve gaps, preserving source order and q.

    ``boundary`` is the placement domain; ``part`` is the material domain P
    the stock decision reads (``benchmarks/held_material_domain.py``). A
    transition is accepted iff its circle pair is within the cap (Held's
    predecessor model), its target is within the cap at all 45 loop stations
    on P minus every earlier circle's swept annulus
    (``benchmarks/held_actual_stock_decision.py``), and its connector, as
    ``Figure5BoundaryConnectors.build`` emits it, is certified by tier 1 or 1b
    of the connector-engagement certificate
    (``benchmarks/held_connector_engagement.py``). Any refusal
    inserts the midpoint circle and retries, depth-bounded. With
    ``close_seam`` the seam connector back to circles[0] is certified too,
    repaired by appending midpoint circles as new last circles; circles[0] is
    never re-cut and the closing pair stays report-only. The post-gate
    ``ConnectorEngagementCertificate`` then passes on the result (with the
    seam only when ``close_seam``).

    Raises:
        Figure5EngagementResolutionError: Empty input, refinement exhaustion,
            or a midpoint whose circle collapses onto its contact (a convex
            corner apex admits no positive circle).
        Figure5EngagementFragmentationError: A transition's swept disks have
            no real intersection in the exact engagement model.
        Figure5EngagementCosineRangeError: A transition's exact engagement
            cosine left [-1, 1].
        ConnectorEngagementResolutionError: A connector (the seam with
            ``close_seam``) is still uncertified when the depth runs out.
        ActualStockEngagementResolutionError: A target is still over the cap
            on its circle view when the repair runs out.
        AmbiguousFigure5BoundaryPathError: A connector has no positive CCW
            progress (as ``Figure5BoundaryConnectors.build`` would refuse it).

    Circle and connector containment are not implied by this model-only step.
    """
    if not circles or len(boundary) < 3 or max_depth < 0 or max_circles < len(circles):
        raise Figure5EngagementResolutionError("Refinement requires circles, a polygon, and sufficient non-negative budgets.")
    _validate_boundary(boundary)
    first = circle_on_side(circles[0], circles[0].boundary_site, circles[0].contact_point, circles[0].radius, boundary)
    key = EngagementRefinementKey(boundary, tool_radius, cap, max_depth, corner_approaches, close_seam, part)
    return continue_figure5_engagement(
        TracedEngagementRefinement((first,), (0,), EngagementTrace(key, ())),
        circles[1:],
        boundary,
        tool_radius,
        cap,
        part=part,
        max_depth=max_depth,
        max_circles=max_circles,
        corner_approaches=corner_approaches,
        close_seam=close_seam,
    )


def continue_figure5_engagement(
    prefix: TracedEngagementRefinement,
    successors: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    part: PartPolygon,
    max_depth: int = 16,
    max_circles: int = MAX_REFINED_CIRCLES,
    corner_approaches: bool = False,
    close_seam: bool = False,
) -> TracedEngagementRefinement:
    """Resume ``refine_figure5_engagement`` after ``prefix``, the refinement of the earlier sources.

    The refinement is a sequential process whose whole state at a source
    boundary is the emitted circle list (tier 1b reads every prior circle), so
    continuing from the exact emission of sources ``0..k`` with the remaining
    sources equals ``refine_figure5_engagement`` on all of them. ``prefix``
    must end at its last source (no seam circles after it); its circles are
    already in emitted form and are never re-cut.

    Raises:
        Figure5EngagementPrefixError: ``prefix`` is empty, its indices do not
            start at 0 and increase, or circles follow its last source.
        Figure5EngagementTraceError: ``prefix`` carries no trace, or one made
            under other settings (``close_seam`` aside: a prefix ends before
            any seam decision).
        Figure5EngagementResolutionError, ConnectorEngagementResolutionError
        and the other errors of ``refine_figure5_engagement``.
    """
    continuation = EngagementContinuation(
        prefix,
        len(successors),
        boundary,
        tool_radius,
        cap,
        part=part,
        max_depth=max_depth,
        max_circles=max_circles,
        corner_approaches=corner_approaches,
        close_seam=close_seam,
    )
    for successor in successors:
        continuation.emit(successor)
    return continuation.finish()


class EngagementContinuation:
    """``continue_figure5_engagement`` one source at a time, for callers that stop it early.

    ``emit`` refines the interval to the next source exactly as the full
    refinement does; ``finish`` closes the seam (with ``close_seam``) and
    returns the refinement with its trace; ``adopt`` instead ends it with a
    reference refinement's remaining circles, for a caller that has proved
    they are what the continuation would emit (the coverage repair's splice).
    """

    def __init__(
        self,
        prefix: TracedEngagementRefinement,
        successors: int,
        boundary: tuple[Point2[WorldXY], ...],
        tool_radius: ToolRadius,
        cap: EngagementCap,
        *,
        part: PartPolygon,
        max_depth: int,
        max_circles: int,
        corner_approaches: bool,
        close_seam: bool,
    ) -> None:
        """Validate ``prefix`` and the budgets for ``successors`` more sources.

        Raises:
            Figure5EngagementPrefixError: see ``continue_figure5_engagement``.
            Figure5EngagementTraceError: see ``continue_figure5_engagement``.
            Figure5EngagementResolutionError: a degenerate polygon or budget.
        """
        indices = prefix.original_indices
        if not indices or indices[0] != 0 or indices[-1] != len(prefix.circles) - 1 or any(a >= b for a, b in zip(indices, indices[1:])):
            raise Figure5EngagementPrefixError(f"prefix indices {indices[:3]}... do not locate {len(prefix.circles)} emitted circles ending at a source")
        if len(boundary) < 3 or max_depth < 0 or max_circles < len(prefix.circles) + successors:
            raise Figure5EngagementResolutionError("Refinement requires circles, a polygon, and sufficient non-negative budgets.")
        _validate_boundary(boundary)
        self.key = EngagementRefinementKey(boundary, tool_radius, cap, max_depth, corner_approaches, close_seam, part)
        # The prefix ends at a source, before any seam decision, so its decisions never read close_seam.
        if not isinstance(prefix, TracedEngagementRefinement) or replace(prefix.trace.key, close_seam=close_seam) != self.key:
            raise Figure5EngagementTraceError("the prefix carries no trace, or one made under other settings")
        self.stock = ActualStockDecision.open(part, tool_radius, cap, EntryCircle(prefix.circles[0]), prefix.circles[1:])
        site = DecisionSite.CONTINUATION if prefix.trace.footprints else DecisionSite.ADVANCE
        gate = _ConnectorGate(boundary, tool_radius, cap)
        self._emission = _Emission(boundary, tool_radius, cap, gate, self.stock, site, max_depth, max_circles, corner_approaches, list(prefix.circles))
        self._prior = tuple(footprint for footprint in prefix.trace.footprints if footprint.position < len(prefix.circles))
        self.original_indices = list(indices)
        self._intervals = len(indices) - 1

    @property
    def emitted(self) -> Sequence[Figure5CounterclockwiseCircle]:
        """The circles emitted so far (read only)."""
        return self._emission.emitted

    def emit(self, source: Figure5CounterclockwiseCircle) -> None:
        """Refine the interval to ``source`` (placed on its side first) and emit it."""
        boundary = self._emission.boundary
        successor = circle_on_side(source, source.boundary_site, source.contact_point, source.radius, boundary)
        bridges = concave_corner_bridges(self._emission.emitted[-1], successor, boundary) if self.key.corner_approaches else ()
        self._emission.advance((*bridges, successor), 0, self._intervals)
        self._intervals += 1
        self.original_indices.append(len(self._emission.emitted) - 1)

    def finish(self) -> TracedEngagementRefinement:
        """Close the seam (with ``close_seam``) and return the refinement with its trace."""
        if self.key.close_seam:
            self._emission.close_seam(self._intervals)
        return self._result((), (), ())

    def last_stock_flip(self, reference: TracedEngagementRefinement, circles_from: int, decisions: tuple[EngagementFootprint, ...]) -> int:
        """The last of ``reference``'s stock ``decisions`` whose verdict differs on this emission followed by ``reference``'s circles from ``circles_from``, or -1.

        ``decisions`` are at positions ``>= circles_from``, in position order. Each is
        re-decided on the view it would have after adoption. The view is cut back
        afterwards, so the continuation can go on.
        """
        size = self.stock.size
        delta = size - circles_from
        last, tail = -1, circles_from
        for decision in decisions:
            while tail < decision.position:
                self.stock.append(reference.circles[tail])
                tail += 1
            if decision.read is not None and self.stock.exceeded(decision.read, decision.position + delta) != decision.refused:
                last = decision.position
        self.stock.truncate(size)
        return last

    def adopt(self, reference: TracedEngagementRefinement, circles_from: int, sources_from: int) -> TracedEngagementRefinement:
        """End with ``reference``'s circles from ``circles_from`` and its sources from ``sources_from`` on.

        The caller's obligation: the continuation from here would emit
        exactly those circles (``reference``'s decisions after
        ``circles_from`` are unchanged by the circles emitted here). Indices
        and footprints shift by the emitted-count difference.

        Raises:
            Figure5EngagementTraceError: ``reference`` has no trace, or one
                made under other settings.
        """
        if not isinstance(reference, TracedEngagementRefinement) or reference.trace.key != self.key:
            raise Figure5EngagementTraceError("the reference carries no trace, or one made under other settings")
        delta = len(self._emission.emitted) - circles_from
        indices = tuple(index + delta for index in reference.original_indices[sources_from:])
        footprints = tuple(footprint.shifted(delta) for footprint in reference.trace.footprints if footprint.position >= circles_from)
        return self._result(reference.circles[circles_from:], indices, footprints)

    def _result(self, tail: tuple[Figure5CounterclockwiseCircle, ...], indices: tuple[int, ...], footprints: tuple[EngagementFootprint, ...]) -> TracedEngagementRefinement:
        trace = EngagementTrace(self.key, (*self._prior, *self._emission.gate.footprints, *footprints), self.stock.counters)
        return TracedEngagementRefinement((*self._emission.emitted, *tail), (*self.original_indices, *indices), trace)
