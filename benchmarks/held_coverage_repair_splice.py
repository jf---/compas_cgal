"""Splice one dropped source into its exact boundary-order slot and re-refine from its slot on.

Cyclic boundary order is compared exactly on ``(segment, Fraction parameter)``.
The continuation stops early, exactly, once the circles it changed can no
longer affect any later decision of the previous refinement
(``continue_with_early_stop``). Design: docs/held_coverage_preserving_selection.md.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass

from benchmarks.held_coverage_repair_errors import InvalidCoverageRepairInputError
from benchmarks.held_coverage_repair_errors import RepairOrderViolationError
from benchmarks.held_coverage_repair_types import BoundaryOrderedEmission
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_coverage_repair_types import SiteKey
from benchmarks.held_coverage_repair_types import Splice
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import EngagementContinuation
from benchmarks.held_figure5_engagement_refinement import Figure5EngagementTraceError
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_figure5_engagement_refinement import circle_on_side
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


def _site_key(site: ProjectionBoundarySite) -> SiteKey:
    return int(site.segment_id), site.parameter


def cyclically_between(first: SiteKey, middle: SiteKey, last: SiteKey) -> bool:
    """Exact cyclic CCW betweenness first <= middle <= last on one loop.

    Equal endpoints denote a fan (no progress), so only that site is between.
    """
    if first < last:
        return first <= middle <= last
    if first > last:
        return middle >= first or middle <= last
    return middle == first


def emitted_form(source: Figure5CounterclockwiseCircle, boundary: tuple[Point2[WorldXY], ...]) -> Figure5CounterclockwiseCircle:
    """The circle exactly as ``refine_figure5_engagement`` emits it for this source."""
    return circle_on_side(source, source.boundary_site, source.contact_point, source.radius, boundary)


@dataclass(frozen=True)
class ContinuedSplice:
    """A spliced refinement and how many of the later sources were re-refined before the early stop adopted the rest."""

    refinement: TracedEngagementRefinement
    rerefined: int


def _same_source(emitted: Figure5CounterclockwiseCircle, source: Figure5CounterclockwiseCircle) -> bool:
    """``emitted`` is the refinement's form of ``source``: ``circle_on_side`` keeps everything but centre and radius."""
    return (emitted.run_id, emitted.station_ordinal, emitted.boundary_site, emitted.contact_point) == (
        source.run_id,
        source.station_ordinal,
        source.boundary_site,
        source.contact_point,
    )


def continue_with_early_stop(
    reference: TracedEngagementRefinement,
    slot: int,
    inserted: Figure5CounterclockwiseCircle,
    later: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
) -> ContinuedSplice:
    """The refinement of ``reference``'s sources with ``inserted`` before source ``slot``, stopping as soon as that is exact.

    ``reference`` refines sources ``S`` (``later`` is ``S[slot:]``) under
    ``settings`` (with another ``close_seam`` its tail differs, so it is
    never adopted). The continuation keeps ``reference``'s circles up to source
    ``slot - 1`` and refines ``inserted`` and then ``later``. After each later
    source whose interval re-emits exactly ``reference``'s circles for it,
    with ``X`` the circles emitted since the slot and ``Y`` ``reference``'s
    circles they replace (the common re-emitted tail removed from both), it
    adopts ``reference``'s remaining circles when (i) no connector decision of
    ``reference`` after this source could read a circle of ``X`` or ``Y``
    (``EngagementTrace.last_influenced``: neither the tier-1b superset
    filter keeps it for any piece of that decision's connector, nor is it
    within ``rho + 2 r`` of a convex-corner floor vertex), (ii) every stock
    decision of ``reference`` after this source whose read set an annulus of
    ``X`` or ``Y`` can reach keeps its verdict when re-decided on the view it
    would have after adoption (``EngagementContinuation.last_stock_flip``), and
    (iii) the circle budget still holds with the count shifted. A dropped ``Y``
    matters for (ii): without its annulus a later circle's stock is larger.

    Why that is exact. The refinement's state at a source boundary is the
    emitted list, and here both lists end in the same circle and continue
    with the same sources. Every later decision reads its predecessor and
    target (equal, inductively), and beyond them only the circles listed in
    its footprint: by (i) those sequences are equal in both runs, so the grid,
    the candidate index and every exact test see the same inputs and decide
    the same. A stock decision's verdict is a function of its target, P and
    the annuli meeting its read set: one no changed circle reaches has the same
    inputs, and one that a changed circle reaches was re-decided by (ii) on
    exactly the view the continuation would have, every earlier verdict being
    equal. So each stock decision takes the same branch. Decisions without a
    footprint (tier 1, the pair cap, the midpoint) read nothing else. By (iii)
    no budget check flips. So the continuation would emit exactly
    ``reference``'s remaining circles.

    Raises:
        InvalidCoverageRepairInputError: ``slot`` is not within ``1..len(sources)``, or
            ``later`` are not ``reference``'s sources from ``slot`` on.
        Figure5EngagementTraceError: ``reference`` has no trace, or one made under other settings.
        Figure5EngagementResolutionError: the cap's depth or the circle budget refuses the continuation.
    """
    indices = reference.original_indices
    if not 1 <= slot <= len(indices):
        raise InvalidCoverageRepairInputError(f"slot {slot} is not within 1..{len(indices)}: the prefix must keep at least the first source")
    if not isinstance(reference, TracedEngagementRefinement):
        raise Figure5EngagementTraceError(f"the reference is a {type(reference).__name__}, which carries no trace")
    if len(indices) != slot + len(later) or not all(_same_source(reference.circles[index], source) for index, source in zip(indices[slot:], later)):
        raise InvalidCoverageRepairInputError(f"the {len(later)} later sources are not the reference's sources from slot {slot} on")
    trace = reference.trace
    start = indices[slot - 1]
    continuation = EngagementContinuation(
        TracedEngagementRefinement(reference.circles[: start + 1], indices[:slot], trace),
        len(later) + 1,
        boundary,
        settings.tool_radius,
        settings.cap,
        part=settings.part,
        max_depth=settings.max_depth,
        max_circles=settings.max_circles,
        corner_approaches=settings.corner_approaches,
        close_seam=settings.close_seam,
    )
    continuation.emit(inserted)
    # A reference refined under another seam setting has another tail: there is nothing to adopt.
    adoptable = trace.key == continuation.key
    common, limit = 0, None
    for offset, source in enumerate(later):
        mark = len(continuation.emitted)
        continuation.emit(source)
        begin, end = indices[slot + offset - 1], indices[slot + offset] + 1
        fresh = continuation.emitted[mark:]
        if not adoptable or tuple(fresh) != reference.circles[begin + 1 : end]:
            common, limit = 0, None
            continue
        common += len(fresh)
        if limit is None:
            changed = (*continuation.emitted[start + 1 : len(continuation.emitted) - common], *reference.circles[start + 1 : end - common])
            limit = max(trace.last_influenced(changed, end), continuation.last_stock_flip(reference, end, trace.stock_decisions_read_by(changed, end)))
        if limit < end and len(reference.circles) + len(continuation.emitted) - end <= settings.max_circles:
            return ContinuedSplice(continuation.adopt(reference, end, slot + offset + 1), offset + 1)
    return ContinuedSplice(continuation.finish(), len(later))


def splice_source(
    sources: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
    emission: BoundaryOrderedEmission,
    candidate: int,
) -> Splice:
    """Insert dropped ``candidate`` at its boundary-order slot and re-refine from (s-, q) to the end.

    Equal to a full ``refine_figure5_engagement`` on the augmented selection:
    the emission up to s- is decided by the sources before the slot alone and
    is kept verbatim; the refinement then continues with q and the later
    sources until the exact early stop (``continue_with_early_stop``) proves
    the rest of the old emission unchanged. A local re-refinement of only
    (s-, q) and (q, s+) is not equal: tier 1b certifies a connector from any
    prior circle, so the window needs the circles before s-, and the
    connectors after s+ gain q's circles as priors and may need fewer
    midpoints; the early stop decides exactly when they no longer do.

    Raises:
        InvalidCoverageRepairInputError: ``candidate`` is already selected or
            the emission does not index its selection.
        RepairOrderViolationError: the candidate's exact site is not between
            its kept neighbours, or it lies outside the selected span.
        Figure5EngagementResolutionError: the cap's depth or the circle budget
            refuses the splice.
        Figure5EngagementTraceError: the emission's refinement carries no
            trace, or one made under other settings.
    """
    selected, refinement = emission.selected, emission.refinement
    if candidate in selected or len(refinement.original_indices) != len(selected):
        raise InvalidCoverageRepairInputError(f"candidate #{candidate} is selected, or the emission does not index its {len(selected)} selected sources")
    slot = bisect.bisect_left(selected, candidate)
    if slot == 0 or slot == len(selected):
        raise RepairOrderViolationError(
            candidate,
            selected[0],
            selected[-1],
            (_site_key(sources[selected[0]].boundary_site), _site_key(sources[candidate].boundary_site), _site_key(sources[selected[-1]].boundary_site)),
        )
    before, after = selected[slot - 1], selected[slot]
    sites = (_site_key(sources[before].boundary_site), _site_key(sources[candidate].boundary_site), _site_key(sources[after].boundary_site))
    if not cyclically_between(*sites):
        raise RepairOrderViolationError(candidate, before, after, sites)
    spliced = continue_with_early_stop(refinement, slot, sources[candidate], tuple(sources[index] for index in selected[slot:]), boundary, settings).refinement
    augmented = (*selected[:slot], candidate, *selected[slot:])
    added = len(spliced.circles) - len(refinement.circles)
    return Splice.build(BoundaryOrderedEmission.build(augmented, spliced), candidate, (before, after), added)
