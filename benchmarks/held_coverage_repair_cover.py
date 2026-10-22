"""Witness-to-candidate matching and the greedy set cover of one repair round.

Membership is decided exactly (``contains_exactly``, ``Fraction`` arithmetic on
the injected binary64 centre and radii); the outward-rounded boxes are a filter
(a superset), never a decision. The greedy cover compares integer counts and
indices only. Design: docs/held_coverage_preserving_selection.md.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Sequence
from fractions import Fraction

import numpy as np
import numpy.typing as npt

from benchmarks.held_coverage_repair_errors import CoverageNotRepairableBySelectionError
from benchmarks.held_coverage_repair_errors import CoverageRepairRefusedError
from benchmarks.held_coverage_repair_splice import emitted_form
from benchmarks.held_coverage_repair_splice import splice_source
from benchmarks.held_coverage_repair_types import BoundaryOrderedEmission
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_coverage_repair_types import RepairInsertion
from benchmarks.held_coverage_repair_types import RepairRefusal
from benchmarks.held_coverage_repair_types import Splice
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementResolutionError
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


def contains_exactly(circle: Figure5CounterclockwiseCircle, tool_radius: ToolRadius, point: tuple[Fraction, Fraction]) -> bool:
    """Exact membership of ``point`` in the swept annulus max(R - r, 0) <= |x - c| <= R + r (closed)."""
    cx, cy = Fraction(float(circle.center.x)), Fraction(float(circle.center.y))
    big_r, r = Fraction(float(circle.radius.value)), Fraction(float(tool_radius.value))
    dx, dy = point[0] - cx, point[1] - cy
    distance2 = dx * dx + dy * dy
    inner = max(big_r - r, Fraction(0))
    return inner * inner <= distance2 <= (big_r + r) * (big_r + r)


def outward_boxes(circles: Sequence[Figure5CounterclockwiseCircle], tool_radius: ToolRadius) -> npt.NDArray[np.float64]:
    """Boxes (xlo, xhi, ylo, yhi) that contain each swept disk of radius R + r.

    Every binary64 operation is correctly rounded (half an ulp), so one
    ``nextafter`` step outward after each rounded operation yields a superset.
    A filter only: membership is decided by ``contains_exactly``.
    """
    centres = np.array([(float(c.center.x), float(c.center.y)) for c in circles], dtype=np.float64).reshape(-1, 2)
    reach = np.array([float(c.radius.value) for c in circles], dtype=np.float64) + float(tool_radius.value)
    reach = np.nextafter(reach, np.inf)
    lo = np.nextafter(centres - reach[:, None], -np.inf)
    hi = np.nextafter(centres + reach[:, None], np.inf)
    return np.stack((lo[:, 0], hi[:, 0], lo[:, 1], hi[:, 1]), axis=1)


def witness_bracket(value: Fraction) -> tuple[float, float]:
    """Doubles ``low <= value <= high``: ``float()`` of a Fraction rounds to nearest, so one step outward brackets it."""
    nearest = float(value)
    return math.nextafter(nearest, -math.inf), math.nextafter(nearest, math.inf)


def candidates_for(
    witness: UncoveredWitness, dropped: Sequence[int], forms: Sequence[Figure5CounterclockwiseCircle], boxes: npt.NDArray[np.float64], tool_radius: ToolRadius
) -> tuple[int, ...]:
    """The dropped sources whose emitted circle contains ``witness`` exactly, in ``dropped`` order.

    Args:
        witness: An exact uncovered point of the certificate.
        dropped: Source indices of the dropped sources.
        forms: ``emitted_form`` of each dropped source, aligned with ``dropped``.
        boxes: ``outward_boxes(forms, tool_radius)``; a filter only.
        tool_radius: Tool radius r.

    Returns:
        Each ``dropped[k]`` whose box overlaps the witness bracket and whose
        swept annulus contains the witness (``contains_exactly``).
    """
    point = witness.uncovered_point
    (xlo, xhi), (ylo, yhi) = witness_bracket(point[0]), witness_bracket(point[1])
    hits = np.nonzero((boxes[:, 0] <= xhi) & (boxes[:, 1] >= xlo) & (boxes[:, 2] <= yhi) & (boxes[:, 3] >= ylo))[0]
    return tuple(dropped[int(k)] for k in hits if contains_exactly(forms[int(k)], tool_radius, (point[0], point[1])))


def _choose(
    sources: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
    emission: BoundaryOrderedEmission,
    members: dict[int, frozenset[int]],
    open_witnesses: frozenset[int],
    round_index: int,
    refusals: list[RepairRefusal],
    refused: dict[tuple[int, tuple[int, int]], tuple[BoundaryOrderedEmission, Figure5EngagementResolutionError]],
) -> tuple[Splice, frozenset[int]] | None:
    """Greedy step: most open witnesses, then fewest added circles, then lowest source index.

    A candidate whose splice raises ``Figure5EngagementResolutionError`` (the
    cap's depth, the circle budget, or a wrapped placement-model failure) is
    refused: recorded once per (candidate, kept interval) and round, and the
    next candidate is tried. The refusal is reused only while the emission
    is the one it was decided on: the splice re-refines from its slot to the
    end, so any other insertion changes its prior circles (tier 1b) or its
    tail, and the candidate is re-spliced.
    """
    best: tuple[tuple[int, int, int], Splice, frozenset[int]] | None = None
    for candidate in sorted(members):
        covered = members[candidate] & open_witnesses
        if not covered:
            continue
        slot = bisect.bisect_left(emission.selected, candidate)
        interval = (emission.selected[slot - 1], emission.selected[slot]) if 0 < slot < len(emission.selected) else (emission.selected[0], emission.selected[-1])
        previous = refused.get((candidate, interval))
        if previous is not None and previous[0] is emission:
            continue
        try:
            splice = splice_source(sources, boundary, settings, emission, candidate)
        except Figure5EngagementResolutionError as error:
            if previous is None:
                refusals.append(RepairRefusal.build(candidate, interval, round_index, f"{type(error).__name__}: {error}"))
            refused[(candidate, interval)] = (emission, error)
            continue
        key = (-len(covered), splice.added, candidate)
        if best is None or key < best[0]:
            best = (key, splice, covered)
    return None if best is None else (best[1], best[2])


def repair_round(
    sources: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
    emission: BoundaryOrderedEmission,
    witnesses: tuple[UncoveredWitness, ...],
    forms: dict[int, Figure5CounterclockwiseCircle],
    round_index: int,
    refusals: list[RepairRefusal],
) -> tuple[BoundaryOrderedEmission, tuple[RepairInsertion, ...]]:
    """One repair round: match every witness to the dropped sources containing it, then greedy-cover them.

    Each greedy step splices the candidate covering the most open witnesses,
    then the fewest added circles, then the lowest source index
    (``splice_source``), until every witness of the round is covered.

    Args:
        sources: Every source circle in boundary order.
        boundary: The CCW centre-domain loop Omega.
        settings: The refinement settings of the emission.
        emission: The current selection and its refinement.
        witnesses: The certificate's witnesses of this round (non-empty).
        forms: Cache of ``emitted_form`` per source index; filled in place.
        round_index: This round's 1-based index, recorded in the ledger.
        refusals: The run's refusal ledger; this round's refusals are appended.

    Returns:
        The emission after every insertion of the round, and the insertions
        in order.

    Raises:
        CoverageNotRepairableBySelectionError: a witness lies in no dropped source.
        CoverageRepairRefusedError: every candidate of an open witness was refused.
        RepairOrderViolationError: a candidate is outside its boundary slot.
    """
    chosen = set(emission.selected)
    dropped = tuple(index for index in range(len(sources)) if index not in chosen)
    for index in dropped:
        if index not in forms:
            forms[index] = emitted_form(sources[index], boundary)
    dropped_forms = tuple(forms[index] for index in dropped)
    boxes = outward_boxes(dropped_forms, settings.tool_radius)
    containing = tuple(candidates_for(witness, dropped, dropped_forms, boxes, settings.tool_radius) for witness in witnesses)
    orphans = tuple(witness for witness, found in zip(witnesses, containing) if not found)
    if orphans:
        raise CoverageNotRepairableBySelectionError(orphans)
    members: dict[int, set[int]] = {}
    for witness_index, found in enumerate(containing):
        for candidate in found:
            members.setdefault(candidate, set()).add(witness_index)
    frozen = {candidate: frozenset(found) for candidate, found in members.items()}
    open_witnesses = frozenset(range(len(witnesses)))
    insertions: list[RepairInsertion] = []
    refused: dict[tuple[int, tuple[int, int]], tuple[BoundaryOrderedEmission, Figure5EngagementResolutionError]] = {}
    while open_witnesses:
        picked = _choose(sources, boundary, settings, emission, frozen, open_witnesses, round_index, refusals, refused)
        if picked is None:
            # Every candidate of every open witness was refused (each open
            # witness has a candidate, or it would be an orphan above).
            named = min(open_witnesses)
            cause = next(error for (candidate, _), (_, error) in refused.items() if candidate in containing[named])
            raise CoverageRepairRefusedError(witnesses[named], refusals) from cause
        splice, covered = picked
        emission = splice.emission
        insertions.append(RepairInsertion.build(splice.candidate, witnesses[min(covered)], round_index))
        frozen.pop(splice.candidate)
        open_witnesses = open_witnesses - covered
    return emission, tuple(insertions)
