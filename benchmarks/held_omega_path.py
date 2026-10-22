"""Held reference paths on Ω*: hypothesis contacts on ∂Ω*, exact boundary order, lane placement.

The Ω* counterpart of :func:`benchmarks.held_figure5_path.build_held_reference_path`,
beside it (the polygon-Ω path is unchanged). Every decision that the polygon
path made on the float offset polygon is made here on Ω*
(docs/held_exact_centre_domain.md, "Integration"):

- **contacts**: every hypothesis contact is decided onto ∂Ω* exactly
  (``CentreDomain.decide_contacts``, the admissible gap as the distance bound);
- **ownership**: a station's hypotheses are grouped by cyclically adjacent
  ∂Ω* pieces (not polygon sides), the nearest of each group owns it;
- **order**: circles sort by the exact site key along ∂Ω*
  (:mod:`benchmarks.held_omega_sites`), the seed is the site exactly nearest the
  published start;
- **circles**: binary64 on the inward normal, safe against their own support
  (:func:`benchmarks.held_omega_circles.place_local`); global containment is
  enforced on the emitted circles by the refinement.

A contact decided exactly onto a convex vertex of Ω* admits no tangent circle
(no circle through the apex of a convex wedge lies in it); such hypotheses are
recorded as refused, never silently dropped. The raw guide and its
hypotheses still read the projection ring: they propose, Ω* decides.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.errors import BenchmarkError
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_exact_centre_domain import ContactBeyondGapError
from benchmarks.held_exact_centre_domain import DecidedContact
from benchmarks.held_figure5_path import _canonical_indices
from benchmarks.held_figure5_path import _cyclic_boundary_neighborhoods
from benchmarks.held_figure5_path import _dyadic_compare
from benchmarks.held_figure5_path import _dyadic_minimum
from benchmarks.held_figure5_path import _dyadic_squared_distance
from benchmarks.held_figure5_path import _PlacementLane
from benchmarks.held_figure5_path import continuous_placement_lanes
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_standard_placement import OrderedPlacementSpan
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacement
from benchmarks.held_standard_placement import _coincident_centres
from benchmarks.held_standard_placement import predecessor_overlap_margin
from benchmarks.held_standard_placement import select_next_in_span
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


class OmegaPathIntegrationError(BenchmarkError):
    """Base failure while integrating hypotheses into an Ω* path."""


class MissingOmegaGuideRunError(OmegaPathIntegrationError):
    """A guide run has no admissible contact on Ω*, or placement lost one."""


class AmbiguousOmegaStartError(OmegaPathIntegrationError):
    """The published start has no unique nearest ∂Ω* site, or exceeds its evidence bound."""


@dataclass(frozen=True)
class OmegaContactSites:
    """Hypotheses with an admissible ∂Ω* decision (parallel tuples), and every refusal with its reason.

    Sites are built from the decisions on demand (:meth:`site`): ownership
    reads only the decided piece.
    """

    loop: OmegaLoop
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...]
    decisions: tuple[DecidedContact, ...]
    beyond_gap: tuple[ProjectionAdmissibleBoundaryHypothesis, ...]
    at_convex_vertex: tuple[ProjectionAdmissibleBoundaryHypothesis, ...]

    def site(self, index: int) -> OmegaSite:
        """The exact site of hypothesis ``index``."""
        return OmegaSite.decided(self.loop, self.decisions[index])


def decide_omega_sites(
    loop: OmegaLoop,
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    gap: ExactMillimetre,
    *,
    qualify: bool,
) -> OmegaContactSites:
    """Decide every hypothesis contact onto ∂Ω* in one native call.

    With ``qualify`` a contact beyond ``gap`` is refused and recorded; without
    it, it fails (the polygon path's bounded projection). A contact at a
    convex vertex is refused and recorded either way.

    Raises:
        ContactBeyondGapError: Without ``qualify``, the first contact beyond the gap.
        MissingOmegaGuideRunError: A run keeps no admissible contact.
        AmbiguousContactProjectionError: A contact has two exactly nearest ∂Ω* points.
    """
    points = [(float(hypothesis.contact_point.x), float(hypothesis.contact_point.y)) for hypothesis in hypotheses]
    decisions = loop.domain.decide_contacts(points, gap).contacts
    kept: list[ProjectionAdmissibleBoundaryHypothesis] = []
    sites: list[DecidedContact] = []
    beyond: list[ProjectionAdmissibleBoundaryHypothesis] = []
    vertex: list[ProjectionAdmissibleBoundaryHypothesis] = []
    for position, (hypothesis, decision) in enumerate(zip(hypotheses, decisions, strict=True)):
        if decision is None:
            if not qualify:
                raise ContactBeyondGapError(f"Contact {points[position]} of run {int(hypothesis.run_id)} is beyond the gap.")
            beyond.append(hypothesis)
            continue
        if decision.is_vertex and loop.convex_sources[decision.piece.index]:
            vertex.append(hypothesis)
            continue
        kept.append(hypothesis)
        sites.append(decision)
    missing = {hypothesis.run_id for hypothesis in hypotheses} - {hypothesis.run_id for hypothesis in kept}
    if missing:
        raise MissingOmegaGuideRunError(f"No admissible Ω* contact for guide runs {sorted(map(int, missing))}.")
    return OmegaContactSites(loop, tuple(kept), tuple(sites), tuple(beyond), tuple(vertex))


def _junction_bridged_indices(piece_kinds: tuple[PieceKind, ...]) -> tuple[tuple[int, ...], int]:
    """Ownership index of every ∂Ω* piece, junction arcs bridged, and the bridged count.

    Non-junction pieces are numbered in loop order; a junction arc takes the
    index of the non-junction piece before it, cyclically. Decided by piece
    kind only.

    Raises:
        OmegaPathIntegrationError: ∂Ω* has no non-junction piece.
    """
    first = next((index for index, kind in enumerate(piece_kinds) if kind is not PieceKind.JUNCTION_ARC), None)
    if first is None:
        raise OmegaPathIntegrationError("∂Ω* has no offset piece to own a station.")
    bridged = [0] * len(piece_kinds)
    current = -1
    for step in range(len(piece_kinds)):
        index = (first + step) % len(piece_kinds)
        if piece_kinds[index] is not PieceKind.JUNCTION_ARC:
            current += 1
        bridged[index] = current
    return tuple(bridged), current + 1


def _owner_indices(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    pieces: tuple[int, ...],
    piece_kinds: tuple[PieceKind, ...],
    published_start_center: Point2[WorldXY],
) -> tuple[int, ...]:
    """The polygon path's neighbourhood ownership over ∂Ω* pieces, junction arcs bridged.

    A junction arc offsets a G1 junction of the declared boundary (chord down
    to zero); hypotheses almost never decide onto one, so a station straddling
    it lands on the pieces either side. Bridging keeps that station one
    neighbourhood, as the polygon rule's projection ring (which has no junction
    element) does. Pieces separated by an offset piece stay separate.
    """
    bridged, piece_count = _junction_bridged_indices(piece_kinds)
    pieces = tuple(bridged[piece] for piece in pieces)
    by_station: dict[tuple[GuideRunId, GuideRunStationOrdinal], list[int]] = {}
    for index, hypothesis in enumerate(hypotheses):
        by_station.setdefault((hypothesis.run_id, hypothesis.station_ordinal), []).append(index)
    owners: list[int] = []
    for station_indices in by_station.values():
        for neighborhood in _cyclic_boundary_neighborhoods({pieces[index] for index in station_indices}, piece_count):
            candidates = tuple(index for index in station_indices if pieces[index] in neighborhood)
            distances = tuple(_dyadic_squared_distance(hypotheses[index].middle_point, hypotheses[index].boundary_footpoint) for index in candidates)
            minimum = _dyadic_minimum(distances)
            owners.extend(index for index, distance in zip(candidates, distances, strict=True) if _dyadic_compare(distance, minimum) == 0)
    start_distances = tuple(_dyadic_squared_distance(published_start_center, hypothesis.contact_point) for hypothesis in hypotheses)
    start_minimum = _dyadic_minimum(start_distances)
    for index, distance in enumerate(start_distances):
        if _dyadic_compare(distance, start_minimum) == 0 and hypotheses[index] not in (hypotheses[owner] for owner in owners):
            owners.append(index)
    return tuple(owners)


def _compare_circles(first: OmegaCircle, second: OmegaCircle) -> int:
    """Emission order: exact ∂Ω* site, then run, station and radius."""
    by_site = first.site.key.compare(second.site.key)
    if by_site != 0:
        return by_site
    left = (int(first.run_id), int(first.station_ordinal), float(first.radius.value))
    right = (int(second.run_id), int(second.station_ordinal), float(second.radius.value))
    return (left > right) - (left < right)


circle_order_key = functools.cmp_to_key(_compare_circles)


def _compare_stations(first: OmegaCircle, second: OmegaCircle) -> int:
    left = (int(first.station_ordinal), first.source_boundary_site.parameter)
    right = (int(second.station_ordinal), second.source_boundary_site.parameter)
    if left != right:
        return (left > right) - (left < right)
    by_site = first.site.key.compare(second.site.key)
    if by_site != 0:
        return by_site
    return (float(first.radius.value) > float(second.radius.value)) - (float(first.radius.value) < float(second.radius.value))


_station_key = functools.cmp_to_key(_compare_stations)


def _paper(circle: OmegaCircle) -> PaperCircleCandidate:
    return circle.paper_candidate()


def phase_from_start(circles: tuple[OmegaCircle, ...], start: Point2[WorldXY], start_bound: Millimetre) -> tuple[OmegaCircle, ...]:
    """Sort by the exact ∂Ω* order and rotate to the site exactly nearest the published start.

    Raises:
        AmbiguousOmegaStartError: Nearest site beyond ``start_bound``, or two distinct nearest sites.
    """
    ordered = sorted(circles, key=circle_order_key)
    origin = ExactPoint2[WorldXY].from_binary64(float(start.x), float(start.y))
    origin_box = (Interval.of(origin.x), Interval.of(origin.y))
    boxes = []
    for circle in ordered:
        x, y = circle.site.point_box()
        boxes.append((x - origin_box[0]).square() + (y - origin_box[1]).square())
    # Filter: only sites whose enclosure can reach the least upper end compete exactly.
    ceiling = min(box.hi for box in boxes)
    competing = [index for index, box in enumerate(boxes) if box.lo <= ceiling]
    exact = {index: ordered[index].site.squared_distance_to(origin) for index in competing}
    minimum = exact[competing[0]]
    for index in competing[1:]:
        if exact[index].compare(minimum) < 0:
            minimum = exact[index]
    bound = Fraction(float(start_bound))
    if minimum.compare(bound * bound) > 0:
        raise AmbiguousOmegaStartError("The ∂Ω* site nearest the published start exceeds the start evidence bound.")
    nearest = [index for index in competing if exact[index].compare(minimum) == 0]
    if any(not ordered[index].site.same_position(ordered[nearest[0]].site) for index in nearest):
        raise AmbiguousOmegaStartError("The published start has no unique nearest ∂Ω* site.")
    seed = nearest[0]
    return tuple(ordered[seed:] + ordered[:seed])


def _merge_identical_circles(placed: list[tuple[OmegaCircle, tuple[GuideRunId, ...]]]) -> tuple[tuple[OmegaCircle, ...], dict[int, tuple[GuideRunId, ...]]]:
    """One circle per distinct binary64 circle, carrying every source run that placed it.

    On Ω* the placed circle is the canonical identity: hypotheses that differ
    only along a ∂Ω* normal decide onto one exact site and, at one guide
    radius, place one circle. The representative is the first in emission
    order (exact site, then run, station, radius).
    """
    groups: dict[tuple[float, float, float], list[tuple[OmegaCircle, tuple[GuideRunId, ...]]]] = {}
    for circle, runs in placed:
        groups.setdefault((float(circle.center.x), float(circle.center.y), float(circle.radius.value)), []).append((circle, runs))
    circles: list[OmegaCircle] = []
    runs_of: dict[int, tuple[GuideRunId, ...]] = {}
    for group in groups.values():
        representative = min((circle for circle, _ in group), key=circle_order_key)
        circles.append(representative)
        runs_of[id(representative)] = tuple(sorted({run for _, runs in group for run in runs}, key=int))
    return tuple(circles), runs_of


def _fold_coincident_centres(lanes: tuple[_PlacementLane[OmegaCircle], ...], runs_of: dict[int, tuple[GuideRunId, ...]]) -> tuple[_PlacementLane[OmegaCircle], ...]:
    """Fold a lane circle into its predecessor when both share one centre exactly.

    A zero restart is the one the exact standard placement cannot place
    (``StandardPlacementResolutionError``); identical circles are already one
    circle, so what folds here is a coincident centre at another radius. Distinct
    centres, sub-ulp apart included, stay: every placement decision is exact. The
    predecessor stays and carries the folded circle's runs (``runs_of`` is
    updated in place).
    """
    merged: list[_PlacementLane[OmegaCircle]] = []
    for lane in lanes:
        kept = [lane.circles[0]]
        for circle in lane.circles[1:]:
            if _coincident_centres(kept[-1].paper_candidate(), circle.paper_candidate()):
                runs_of[id(kept[-1])] = tuple(sorted({*runs_of[id(kept[-1])], *runs_of[id(circle)]}, key=int))
                continue
            kept.append(circle)
        merged.append(_PlacementLane(lane.source_segment, lane.direction, tuple(kept)))
    return tuple(merged)


def _split_at_gaps(lanes: tuple[_PlacementLane[OmegaCircle], ...], tool_radius: ToolRadius) -> tuple[_PlacementLane[OmegaCircle], ...]:
    """Split every lane where two consecutive circles leave the paper overlap domain.

    On Ω* a lane can lose circles the polygon lane kept (contacts refused at a
    convex vertex, circles merged into another run's identical circle), and
    standard placement refuses a successor beyond the overlap bound. The split
    uses the rule the lane joiner and the placement use
    (:func:`benchmarks.held_standard_placement.predecessor_overlap_margin` > 0,
    an exact sign), so every resulting lane is one the placement accepts.
    """
    split: list[_PlacementLane[OmegaCircle]] = []
    for lane in lanes:
        start = 0
        for index in range(1, len(lane.circles)):
            if predecessor_overlap_margin(lane.circles[index - 1].paper_candidate(), lane.circles[index].paper_candidate(), tool_radius) <= 0:
                split.append(_PlacementLane(lane.source_segment, lane.direction, lane.circles[start:index]))
                start = index
        split.append(_PlacementLane(lane.source_segment, lane.direction, lane.circles[start:]))
    return tuple(split)


@dataclass(frozen=True)
class OmegaPlacementRecord:
    """One predecessor placement and every guide run it represents."""

    site: OmegaSite
    source_run_ids: tuple[GuideRunId, ...]
    placement: StandardPlacement


@dataclass(frozen=True)
class OmegaHypothesisPath:
    """Placed circles in phased ∂Ω* order, with provenance and every refused contact."""

    circles: tuple[OmegaCircle, ...]
    circle_sources: tuple[tuple[GuideRunId, ...], ...]
    placements: tuple[OmegaPlacementRecord, ...]
    reached_run_ids: frozenset[GuideRunId]
    candidate_count: int
    canonical_candidate_count: int
    placement_lane_count: int
    contacts: OmegaContactSites


def build_omega_hypothesis_path(
    *,
    contacts: OmegaContactSites,
    reached_run_ids: frozenset[GuideRunId],
    published_start_center: Point2[WorldXY],
    published_start_radius: Millimetre,
    start_evidence_bound: Millimetre,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    preserve_source_runs: bool,
) -> OmegaHypothesisPath:
    """Own, canonicalise, place on Ω* and thin by standard placement along lanes (the polygon path's stages).

    Raises:
        OmegaPathIntegrationError: No hypotheses, or a published start radius off the tool evidence.
        MissingOmegaGuideRunError: Placement did not preserve every reached run.
        AmbiguousOmegaStartError: See :func:`phase_from_start`.
    """
    hypotheses = contacts.hypotheses
    if not hypotheses or not reached_run_ids:
        raise OmegaPathIntegrationError("An Ω* path requires hypotheses and a reached-run ledger.")
    radius_gap = abs(Fraction(float(published_start_radius)) - Fraction(float(tool_radius.value)))
    if radius_gap > Fraction(float(start_evidence_bound)):
        raise OmegaPathIntegrationError("The published start radius does not match the tool evidence.")
    owners = _owner_indices(
        hypotheses, tuple(decision.piece.index for decision in contacts.decisions), tuple(piece.kind for piece in contacts.loop.pieces), published_start_center
    )
    owned = tuple(hypotheses[index] for index in owners)
    canonical = _canonical_indices(owned)
    placed: list[tuple[OmegaCircle, tuple[GuideRunId, ...]]] = []
    for index, runs in canonical:
        hypothesis = owned[index]
        provenance = (hypothesis.run_id, hypothesis.station_ordinal, hypothesis.projected_boundary_site)
        placed.append((place_local(contacts.site(owners[index]), float(hypothesis.guide_radius.value), provenance), runs))
    circles, runs_of = _merge_identical_circles(placed)
    ordered = phase_from_start(circles, published_start_center, start_evidence_bound)
    lanes = _split_at_gaps(_fold_coincident_centres(continuous_placement_lanes(ordered, tool_radius, _station_key, _paper), runs_of), tool_radius)
    selected: list[OmegaCircle] = []
    placements: list[OmegaPlacementRecord] = []
    for lane in lanes:
        lane_circles = lane.circles
        candidates = tuple(circle.paper_candidate() for circle in lane_circles)
        selected.append(lane_circles[0])
        stops: list[int] = []
        represented = set(runs_of[id(lane_circles[0])])
        if preserve_source_runs:
            for position, circle in enumerate(lane_circles[1:], start=1):
                if any(run not in represented for run in runs_of[id(circle)]):
                    stops.append(position)
                    represented.update(runs_of[id(circle)])
        stops.append(len(candidates) - 1)
        cursor = 0
        for stop in stops:
            if cursor >= stop:
                continue
            span_start = cursor
            span = OrderedPlacementSpan.build(candidates[span_start : stop + 1], tool_radius, cap)
            while cursor < stop:
                placement = select_next_in_span(span, cursor - span_start)
                cursor += placement.candidate_advance
                chosen = lane_circles[cursor]
                selected.append(chosen)
                placements.append(OmegaPlacementRecord(chosen.site, runs_of[id(chosen)], placement))
    selected_runs = frozenset(run for circle in selected for run in runs_of[id(circle)])
    if selected_runs != reached_run_ids:
        raise MissingOmegaGuideRunError("Engagement placement did not preserve every reached guide run.")
    # A circle reached by two lanes is emitted once (the polygon path keyed the same circle by origin).
    unique = tuple({id(circle): circle for circle in selected}.values())
    phased = phase_from_start(unique, published_start_center, start_evidence_bound)
    return OmegaHypothesisPath(
        phased,
        tuple(runs_of[id(circle)] for circle in phased),
        tuple(placements),
        reached_run_ids,
        len(hypotheses),
        len(canonical),
        len(lanes),
        contacts,
    )


def omega_reference_path(
    loop: OmegaLoop,
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    reached_run_ids: frozenset[GuideRunId],
    *,
    gap: ExactMillimetre,
    published_start_center: Point2[WorldXY],
    published_start_radius: Millimetre,
    start_evidence_bound: Millimetre,
    tool_radius: ToolRadius,
    cap: EngagementCap,
    qualify_contacts: bool,
    preserve_source_runs: bool,
) -> OmegaHypothesisPath:
    """Decide the contacts on ∂Ω*, then place (see :func:`build_omega_hypothesis_path`)."""
    contacts = decide_omega_sites(loop, hypotheses, gap, qualify=qualify_contacts)
    return build_omega_hypothesis_path(
        contacts=contacts,
        reached_run_ids=reached_run_ids,
        published_start_center=published_start_center,
        published_start_radius=published_start_radius,
        start_evidence_bound=start_evidence_bound,
        tool_radius=tool_radius,
        cap=cap,
        preserve_source_runs=preserve_source_runs,
    )
