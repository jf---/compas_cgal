"""Engagement spacing and refinement of Ω* paths: boundary-order selection, exact containment, ∂Ω* midpoints.

The Ω* counterparts of :func:`benchmarks.held_boundary_order_spacing.select_boundary_ordered_sources`
and :func:`benchmarks.held_figure5_engagement_refinement.refine_figure5_engagement`.
Engagement is the same exact standard-model decision
(:func:`benchmarks.held_standard_placement.exact_predecessor_engagement`) on the
binary64 circles. What changes is the domain:

- every emitted circle is contained in Ω* exactly
  (:func:`benchmarks.held_omega_circles.contain`, replacing the polygon's
  corner-radius limit);
- a midpoint is an exact ∂Ω* site half-way (by reported length) along the
  CCW chain between the two sites, strictly between them by exact comparison;
- there are no concave-corner bridges: ∂Ω* has no reflex vertex (Ω* has an
  exterior r-disk at every boundary point), so the polygon path's
  ``corner_approaches`` has nothing to do here.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.errors import BenchmarkError
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_circles import contain
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_connectors import OmegaConnectorPiece
from benchmarks.held_omega_connectors import _chain
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_omega_sites import SiteKey
from benchmarks.held_omega_sites import SiteOutsidePieceError
from benchmarks.held_omega_sites import _Placement
from benchmarks.held_standard_placement import ExactEngagementCosineRangeError
from benchmarks.held_standard_placement import PredecessorEngagementCosine
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# Explicit resource ceilings, as the polygon refinement's.
MAX_OMEGA_REFINED_CIRCLES = 20000
MAX_OMEGA_REFINEMENT_DEPTH = 16


class OmegaRefinementError(BenchmarkError):
    """Base failure of Ω* engagement refinement."""


class InvalidOmegaSpacingInputError(OmegaRefinementError):
    """Circle order and non-empty source-run attribution must correspond."""


class OmegaEngagementResolutionError(OmegaRefinementError):
    """A transition cannot be refined within the depth or circle budget."""


class OmegaEngagementFragmentationError(OmegaEngagementResolutionError):
    """A transition's swept disks have no real intersection in the exact engagement model."""


class OmegaEngagementCosineRangeError(OmegaEngagementResolutionError):
    """A transition's exact engagement cosine left [-1, 1]."""


class DegenerateOmegaMidpointError(OmegaRefinementError):
    """No representable site lies strictly between two distinct sites (their chain is below binary64 resolution)."""


def _engagement(predecessor: OmegaCircle, target: OmegaCircle, tool_radius: ToolRadius, interval: int) -> PredecessorEngagementCosine:
    try:
        return exact_predecessor_engagement(predecessor.paper_candidate(), target.paper_candidate(), tool_radius)
    except StandardPlacementFragmentationError as error:
        raise OmegaEngagementFragmentationError(f"Fragmented transition at source interval {interval}: {error}") from error
    except ExactEngagementCosineRangeError as error:
        raise OmegaEngagementCosineRangeError(f"Engagement cosine out of range at source interval {interval}: {error}") from error


def select_omega_boundary_sources(
    circles: tuple[OmegaCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
) -> tuple[int, ...]:
    """Linear look-ahead selection in emitted order (the polygon path's rule, on Ω* circles).

    Keeps the endpoints and the last representative of an unseen run; keeps the
    current circle when skipping it would exceed the cap or leave the model's
    overlap domain.

    Raises:
        InvalidOmegaSpacingInputError: If attribution does not match the circles.
    """
    if not circles or len(circles) != len(circle_sources) or any(not runs or len(set(runs)) != len(runs) for runs in circle_sources):
        raise InvalidOmegaSpacingInputError("Each ordered circle needs distinct, non-empty source-run attribution.")
    remaining = Counter(run for runs in circle_sources for run in runs)
    represented = set(circle_sources[0])
    selected = [0]
    remaining.subtract(circle_sources[0])
    for index in range(1, len(circles) - 1):
        remaining.subtract(circle_sources[index])
        required = any(run not in represented and remaining[run] == 0 for run in circle_sources[index])
        if not required:
            try:
                engagement = exact_predecessor_engagement(circles[selected[-1]].paper_candidate(), circles[index + 1].paper_candidate(), tool_radius)
            except StandardPlacementFragmentationError:
                required = True
            else:
                required = not engagement.within(cap)
        if required:
            selected.append(index)
            represented.update(circle_sources[index])
    if len(circles) > 1:
        selected.append(len(circles) - 1)
    return tuple(selected)


def _lower_key(sub: OmegaConnectorPiece, start: OmegaSite) -> SiteKey:
    return start.key if sub.index == start.index and sub.start == start.point() else SiteKey(sub.index, _Placement.SOURCE)


def _strictly_inside(candidate: OmegaSite, lower: SiteKey, upper: SiteKey | None) -> bool:
    """Exact: the candidate lies strictly after the sub-piece's start and strictly before its end.

    ``upper`` is the end site's key when the sub-piece ends there; otherwise
    the sub-piece ends at its piece's target, which a projected site is
    already strictly before (:meth:`OmegaSite.project` decides it).
    """
    if candidate.at_vertex or not lower < candidate.key:
        return False
    return upper is None or candidate.key < upper


def _probe(sub: OmegaConnectorPiece, fraction: float) -> ExactPoint2[WorldXY]:
    """A rational probe whose projection lies at about ``fraction`` of the sub-piece (proposal; checked exactly)."""
    support = sub.piece.support
    if isinstance(support, OffsetLineSupport):
        begin, finish = float(support.parameter(sub.start)), float(support.parameter(sub.end))
        along = Fraction(begin + fraction * (finish - begin))
        foot = support.anchor.offset(support.direction, along)
        return ExactPoint2[WorldXY].build(foot.x, foot.y)
    cx, cy = float(support.centre.x), float(support.centre.y)
    sx, sy = sub.start.as_floats()
    ex, ey = sub.end.as_floats()
    begin = math.atan2(sy - cy, sx - cx)
    sweep = (support.sense * (math.atan2(ey - cy, ex - cx) - begin)) % math.tau
    if sub.full or sweep == 0.0:
        sweep = math.tau if sub.full else sweep
    angle = begin + support.sense * fraction * sweep
    radius = float(support.radius)
    return ExactPoint2[WorldXY].from_binary64(cx + radius * math.cos(angle), cy + radius * math.sin(angle))


def omega_midpoint_site(loop: OmegaLoop, start: OmegaSite, end: OmegaSite) -> OmegaSite:
    """An exact ∂Ω* site strictly between two sites, near half the reported chain length.

    Raises:
        DegenerateOmegaMidpointError: If no tried probe lands strictly inside (a sub-ulp chain).
    """
    if start.same_position(end):
        return start
    chain = _chain(loop.pieces, start, end)
    lengths = [sub.reported_length() for sub in chain]
    total = sum(lengths)
    for share in (0.5, 0.25, 0.75):
        target = share * total
        walked = 0.0
        for sub, length in zip(chain, lengths):
            if walked + length < target and sub is not chain[-1]:
                walked += length
                continue
            fraction = 0.5 if length <= 0.0 else min(max((target - walked) / length, 0.0), 1.0)
            for local in (fraction, 0.5):
                try:
                    candidate = OmegaSite.project(loop.frames[sub.index], _probe(sub, local))
                except SiteOutsidePieceError:
                    continue
                upper = end.key if sub is chain[-1] and not end.at_vertex else None
                if _strictly_inside(candidate, _lower_key(sub, start), upper):
                    return candidate
            break
    raise DegenerateOmegaMidpointError(f"No representable site strictly between pieces {start.index} and {end.index}.")


@dataclass(frozen=True)
class OmegaRefinement:
    """Refined Ω* motion; ``original_indices`` locate the input circles among the emitted ones."""

    circles: tuple[OmegaCircle, ...]
    original_indices: tuple[int, ...]
    limited: int


def refine_omega_engagement(
    circles: tuple[OmegaCircle, ...],
    loop: OmegaLoop,
    clearance: DeclaredClearance,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    max_depth: int = MAX_OMEGA_REFINEMENT_DEPTH,
    max_circles: int = MAX_OMEGA_REFINED_CIRCLES,
) -> OmegaRefinement:
    """Contain every circle in Ω* exactly, then subdivide every over-cap transition along ∂Ω*.

    Raises:
        OmegaEngagementResolutionError: Empty input, depth or circle budget exhausted.
        OmegaEngagementFragmentationError: See the engagement model.
        OmegaEngagementCosineRangeError: See the engagement model.
        NoContainedCircleError: A site admits no contained circle.
        DegenerateOmegaMidpointError: See :func:`omega_midpoint_site`.
    """
    if not circles or max_depth < 0 or max_circles < len(circles):
        raise OmegaEngagementResolutionError("Refinement requires circles and sufficient non-negative budgets.")
    domain = loop.domain
    targets = tuple(contain(circle, clearance, domain) for circle in circles)
    emitted = [targets[0]]
    original = [0]
    for interval, successor in enumerate(targets[1:]):
        pending = [(successor, 0)]
        while pending:
            target, depth = pending.pop()
            predecessor = emitted[-1]
            engagement = _engagement(predecessor, target, tool_radius, interval)
            if engagement.within(cap):
                if len(emitted) >= max_circles:
                    raise OmegaEngagementResolutionError(f"Circle budget exhausted at source interval {interval}.")
                emitted.append(target)
                continue
            if depth >= max_depth:
                raise OmegaEngagementResolutionError(f"Depth exhausted at source interval {interval}: {math.degrees(float(engagement.reported_angle())):.6f} degrees.")
            site = omega_midpoint_site(loop, predecessor.site, target.site)
            radius = (float(predecessor.radius.value) + float(target.radius.value)) / 2
            provenance = (predecessor.run_id, predecessor.station_ordinal, predecessor.source_boundary_site)
            midpoint = contain(place_local(site, radius, provenance), clearance, domain)
            pending.append((target, depth + 1))
            pending.append((midpoint, depth + 1))
        original.append(len(emitted) - 1)
    return OmegaRefinement(tuple(emitted), tuple(original), sum(1 for circle in emitted if circle.limited))


def refine_omega_boundary_ordered(
    circles: tuple[OmegaCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    reached_run_ids: frozenset[GuideRunId],
    loop: OmegaLoop,
    clearance: DeclaredClearance,
    tool_radius: ToolRadius,
    cap: EngagementCap,
) -> tuple[tuple[int, ...], OmegaRefinement]:
    """Select in emitted order, then refine (the polygon path's boundary spacing on Ω*).

    Raises:
        InvalidOmegaSpacingInputError: If the selection loses a reached run.
    """
    selected = select_omega_boundary_sources(circles, circle_sources, tool_radius, cap)
    represented = frozenset(run for index in selected for run in circle_sources[index])
    if represented != reached_run_ids:
        raise InvalidOmegaSpacingInputError("Boundary spacing source attribution differs from the reached-run ledger.")
    return selected, refine_omega_engagement(tuple(circles[index] for index in selected), loop, clearance, tool_radius, cap)
