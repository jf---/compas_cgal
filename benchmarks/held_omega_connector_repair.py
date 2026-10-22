"""Repair of Ω* connectors the exact engagement gate refuses, by refinement along ∂Ω*.

User ruling (polygon path, applied to Ω*): a connector without a tier-1 or
tier-1b certificate is repaired by refinement, not reported. For each refused
connector ``i`` (circles ``i -> i+1``) the refinement's own midpoint circle
(:func:`benchmarks.held_omega_refinement.omega_midpoint_circle`: an exact ∂Ω*
site strictly between the two sites, contained in Ω* exactly) is inserted, and
again between the new neighbours, until every sub-connector is accepted:

- its circle pair is within the cap by the exact standard-model decision
  (:func:`benchmarks.held_omega_refinement.transition_engagement`), re-decided
  for every inserted circle;
- its exact connector is certified by tier 1 or 1b
  (:class:`benchmarks.held_omega_connector_gate.OmegaConnectorGate`, the
  post-gate's decision).

The seam (last circle back to circle 0) is repaired by appending midpoint
circles as new last circles; circle 0 is never re-cut and the closing pair
stays report-only, as on the polygon path.

Gouge-free by construction: every inserted circle is contained in Ω* exactly
and every connector is an exact ∂Ω* chain.

Why one pass suffices. Connectors certified on the input stay certified:
an insertion between ``i`` and ``i+1`` leaves every other connector's two
circles adjacent and only adds priors to connectors after it (more tier-1b
candidates; tier 1 reads only the circle just cut). Their verdicts are carried
over with circle indices renumbered; the certificate they state is unchanged.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from dataclasses import replace

from benchmarks.held_connector_engagement_exact import ExactConnectorAssessment
from benchmarks.held_connector_engagement_exact import ExactConnectorVerdict
from benchmarks.held_connector_engagement_exact import ExactRefusal
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_declared_clearance import CircleClearanceOracle
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connector_engagement import exact_omega_connector
from benchmarks.held_omega_connector_engagement import omega_connector_path
from benchmarks.held_omega_connector_engagement import prior_circle
from benchmarks.held_omega_connector_gate import OmegaConnectorGate
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_connectors import connect
from benchmarks.held_omega_convex_corner import sharp_vertex
from benchmarks.held_omega_refinement import OmegaRefinementError
from benchmarks.held_omega_refinement import omega_midpoint_circle
from benchmarks.held_omega_refinement import transition_engagement
from benchmarks.held_omega_sites import OmegaLoop
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius

# Explicit resource ceilings. Depth: midpoint halvings of one refused
# connector's ∂Ω* chain (2^-16 of a connector is far below every corpus tool
# radius, where the lemma's own circle covers any sweep it can reach). Circles:
# inserted over the whole path, the polygon path's <= 2 % budget of the largest
# corpus path (Monstera, ~4400 circles) with a 2x margin, as a hard ceiling.
MAX_OMEGA_CONNECTOR_REPAIR_DEPTH = 16
MAX_OMEGA_CONNECTOR_REPAIR_CIRCLES = 200


class OmegaConnectorRepairError(OmegaRefinementError):
    """Base failure of Ω* connector repair."""


class OmegaConnectorRepairDepthError(OmegaConnectorRepairError):
    """A refused connector is not repaired within the midpoint depth; the last refusals are named."""

    def __init__(self, connector: int, seam: bool, refusals: tuple[ExactRefusal, ...], detail: str) -> None:
        super().__init__(f"Connector {connector}{' (seam)' if seam else ''} not repaired: {detail}")
        self.connector = connector
        self.seam = seam
        self.refusals = refusals


class OmegaConnectorRepairBudgetError(OmegaConnectorRepairError):
    """Repair would insert more circles than the global budget allows; ``connector`` is the input connector being repaired when it ran out."""

    def __init__(self, connector: int, seam: bool, budget: int) -> None:
        super().__init__(f"Connector repair budget of {budget} circles exhausted at connector {connector}{' (seam)' if seam else ''}.")
        self.connector = connector
        self.seam = seam
        self.budget = budget


class OmegaConnectorRepairContractError(OmegaConnectorRepairError):
    """The repaired verdicts do not pair one-to-one with the connectors of the repaired circles."""


class OmegaConvexCornerEngagementError(OmegaConnectorRepairError):
    """A connector first reaches a convex ∂Ω* vertex whose exterior angle exceeds the cap: a true violation no refinement removes.

    Raised before any midpoint is tried on it (:func:`benchmarks.held_omega_convex_corner.sharp_vertex`).
    ``piece`` is the outgoing piece; the vertex is its source.
    """

    def __init__(self, connector: int, seam: bool, piece: int) -> None:
        super().__init__(f"Connector {connector}{' (seam)' if seam else ''} first reaches the convex vertex at the source of piece {piece}, over the cap.")
        self.connector = connector
        self.seam = seam
        self.piece = piece


@dataclass(frozen=True)
class OmegaConnectorRepairAction:
    """Input connector ``connector`` was refused (its first blocking ``refusals``) and repaired with ``inserted`` circles."""

    connector: int
    seam: bool
    refusals: tuple[ExactRefusal, ...]
    inserted: int


@dataclass(frozen=True)
class OmegaConnectorRepair:
    """Emitted circles whose every connector is certified, their connectors and the verdicts.

    ``verdicts[j]`` certifies ``connectors.connectors[j]`` on ``circles``
    (tier 1, tier 1b or a stationary fan; never refused). ``positions[k]`` is
    the index in ``circles`` of input circle ``k``. ``actions`` lists the input
    connectors that needed repair. ``gate_seconds`` is the input assessment
    (the gate), ``repair_seconds`` the insertion pass.
    """

    circles: tuple[OmegaCircle, ...]
    positions: tuple[int, ...]
    connectors: OmegaConnectors
    verdicts: tuple[ExactConnectorVerdict, ...]
    actions: tuple[OmegaConnectorRepairAction, ...]
    gate_seconds: float
    repair_seconds: float

    @property
    def inserted(self) -> int:
        return sum(action.inserted for action in self.actions)


def _renumbered(verdict: ExactConnectorVerdict, position: int, circle_index: list[int]) -> ExactConnectorVerdict:
    subs = tuple(replace(sub, circle=circle_index[sub.circle]) for sub in verdict.sub_sweeps)
    return replace(verdict, connector=position, sub_sweeps=subs)


def repair_omega_connectors(
    circles: tuple[OmegaCircle, ...],
    loop: OmegaLoop,
    clearance: CircleClearanceOracle,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    max_depth: int = MAX_OMEGA_CONNECTOR_REPAIR_DEPTH,
    max_circles: int = MAX_OMEGA_CONNECTOR_REPAIR_CIRCLES,
) -> OmegaConnectorRepair:
    """Gate every connector of the emitted ``circles`` exactly, then repair each refused one by ∂Ω* midpoint circles.

    ``circles`` are the path's circles in cut order, as emitted (after every
    pass that removes circles: a removal can leave a connector unproved).

    Raises:
        OmegaConnectorRepairDepthError: A refused connector is not repaired within ``max_depth``.
        OmegaConnectorRepairBudgetError: More than ``max_circles`` circles would be inserted.
        OmegaConvexCornerEngagementError: A refused sub-connector first reaches a convex vertex over the cap.
        OmegaEngagementFragmentationError: See the engagement model.
        OmegaEngagementCosineRangeError: See the engagement model.
        DegenerateOmegaMidpointError: See :func:`benchmarks.held_omega_refinement.omega_midpoint_site`.
        NoContainedCircleError: A midpoint site admits no contained circle.
        OmegaCornerFirstArrivalError: The circles do not lap ∂Ω* once in order (a caller contract breach).
        OmegaConnectorRepairContractError: The verdicts do not pair with the repaired connectors.
    """
    gate = OmegaConnectorGate(tool_radius, cap)
    started = time.perf_counter()
    connectors = OmegaConnectors.build(loop.pieces, circles)
    path = omega_connector_path(circles, connectors)
    assessment = ExactConnectorAssessment.build(path, gate.lemma.tool_radius, cap.chord_ratio)
    gate_seconds = time.perf_counter() - started
    started = time.perf_counter()
    if not assessment.refused:
        return OmegaConnectorRepair(circles, tuple(range(len(circles))), connectors, assessment.verdicts, (), gate_seconds, time.perf_counter() - started)
    count = len(circles)
    emitted: list[OmegaCircle] = [circles[0]]
    priors = [path.circles[0]]
    circle_index = [0]
    verdicts: list[ExactConnectorVerdict] = []
    actions: list[OmegaConnectorRepairAction] = []
    budget = max_circles
    for i, verdict in enumerate(assessment.verdicts):
        successor = circles[(i + 1) % count]
        if verdict.certified:
            verdicts.append(_renumbered(verdict, len(emitted) - 1, circle_index))
        else:
            inserted = _repair(i, verdict, successor, emitted, priors, verdicts, loop, clearance, tool_radius, cap, gate, max_depth, budget, max_circles)
            budget -= inserted
            actions.append(OmegaConnectorRepairAction(i, verdict.seam, verdict.refusals, inserted))
        if not verdict.seam:
            if verdict.certified:
                emitted.append(successor)
                priors.append(path.circles[i + 1])
            circle_index.append(len(emitted) - 1)
    repaired = tuple(emitted)
    final = OmegaConnectors.build(loop.pieces, repaired)
    _require_paired(verdicts, final, assessment.verdicts[-1].seam)
    return OmegaConnectorRepair(repaired, tuple(circle_index), final, tuple(verdicts), tuple(actions), gate_seconds, time.perf_counter() - started)


def _require_paired(verdicts: list[ExactConnectorVerdict], connectors: OmegaConnectors, seam: bool) -> None:
    """Fail loud unless verdict ``j`` labels connector ``j`` of the repaired circles, with the input's seam flag on the last one only.

    Raises:
        OmegaConnectorRepairContractError: A count, label or seam-flag mismatch.
    """
    count = len(connectors.connectors)
    if len(verdicts) != count:
        raise OmegaConnectorRepairContractError(f"{len(verdicts)} verdicts for {count} repaired connectors")
    for j, verdict in enumerate(verdicts):
        expected_seam = seam and j == count - 1
        if verdict.connector != j or verdict.seam is not expected_seam:
            raise OmegaConnectorRepairContractError(f"verdict {j} labels connector {verdict.connector} (seam {verdict.seam}), expected {j} (seam {expected_seam})")


def _repair(
    connector: int,
    verdict: ExactConnectorVerdict,
    successor: OmegaCircle,
    emitted: list[OmegaCircle],
    priors: list[PriorCircle],
    verdicts: list[ExactConnectorVerdict],
    loop: OmegaLoop,
    clearance: CircleClearanceOracle,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    gate: OmegaConnectorGate,
    max_depth: int,
    budget: int,
    max_circles: int,
) -> int:
    """Insert midpoint circles before ``successor`` until every sub-connector is accepted; returns how many.

    On a non-seam connector ``successor`` is emitted last; on the seam it is
    circle 0, which closes the loop and is not emitted again. The input
    connector is not decided again: ``verdict`` already refused it.
    ``budget`` is what remains of the global ``max_circles``.
    """
    seam = verdict.seam
    pending: list[tuple[OmegaCircle, int, bool]] = [(successor, 0, True)]
    inserted = 0
    last_refusals = verdict.refusals
    while pending:
        target, depth, refused = pending.pop()
        closing = seam and target is successor
        predecessor = emitted[-1]
        pair_within = closing or transition_engagement(predecessor, target, tool_radius, connector).within(cap)
        if pair_within and not refused:
            exact = exact_omega_connector(predecessor, target, connect(loop.pieces, predecessor, target))
            decided = gate.decide(priors, exact, len(emitted) - 1, seam=closing)
            if decided.certified:
                verdicts.append(decided)
                if not closing:
                    emitted.append(target)
                    priors.append(prior_circle(target))
                continue
            last_refusals = decided.refusals
        corner = sharp_vertex(loop, emitted, target, gate.lemma.cos_sq_half_cap, gate.lemma.tool_radius)
        if corner is not None:
            raise OmegaConvexCornerEngagementError(connector, seam, corner)
        if depth >= max_depth:
            detail = f"depth {max_depth} exhausted ({'connector refused' if pair_within else 'pair over the cap'})"
            raise OmegaConnectorRepairDepthError(connector, seam, last_refusals, detail)
        if inserted >= budget:
            raise OmegaConnectorRepairBudgetError(connector, seam, max_circles)
        midpoint = omega_midpoint_circle(loop, predecessor, target, clearance)
        inserted += 1
        pending.append((target, depth + 1, False))
        pending.append((midpoint, depth + 1, False))
    return inserted
