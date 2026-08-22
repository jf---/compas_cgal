"""Chain ordering: reach the next chain through stock the earlier chains already cleared.

Every generator in this package so far machines each skeleton chain as an island:
plunge, cut, retract, rapid across at clearance height, plunge again. On
`rect_12x8` that is FIVE material entries for five chains, and each one is a
full-immersion cut into virgin stock -- 360 degrees, bounded below by the geometry
in `docs/loop_radius_degeneracy.md` and unfixable by any choice of radius. Four of
those five entries do not have to exist. The rectangle's four corner chains each
run from a corner tip inward to the SPINE, and by the time the spine has been
machined their inner ends sit in cleared material.

This module orders and orients the chains so that every chain after the first
begins where the tool already is, in stock that is already open, and links into it
with a CUT MOVE AT CUTTING DEPTH instead of a retract, a rapid and a plunge.

THE RULE IS THE SAME EXACT PREDICATE, NOT A DISTANCE HEURISTIC
--------------------------------------------------------------
A link is taken at cutting depth only when EVERY evaluated position along it
reports the engagement cap not exceeded, decided by `_stock_2.engagement_at`'s
exact `cap_exceeded` boolean against the depleting stock. "Already cleared" is
therefore a measured property of this pocket at this moment, not an assumption
about rectangles. Where no orientation of any remaining chain has an admissible
link, the tool retracts and plunges exactly as before -- and the chain is named in
`OrderedToolpathResult.isolated_chains` with the engagement that refused it, so
"this pocket needs five entries" is a reported measurement rather than a silent
default.

ORIENTATION IS PART OF THE DECISION, NOT A CONVENTION
------------------------------------------------------
A chain has two ends and they are not interchangeable. A corner chain entered at
its tip starts in virgin material; the same chain entered at its spine end starts
in cleared material and walks outward. Both orientations of every remaining chain
are offered to the link test, so the walk direction is chosen by the same
predicate that chooses the link.

WHAT THIS COSTS, AND THE ONE THING IT ADDS
-------------------------------------------
Walking a corner chain from the spine towards the tip means its BRIDGES advance
into virgin material rather than out of it. The bridge between two machining
circles was never cap-checked by any generator in this package -- the older
modules say so explicitly -- and on a tip-ward walk that omission shows up as
slotting motions. So this module checks the bridge with the same predicate it
checks the link: a candidate advance is admissible only when its machining circle
AND the bridge that reaches it are both under the cap. That shortens some advances
and lengthens the path; the trade is measured rather than assumed.

WHAT THIS GUARANTEES, EXACTLY
-----------------------------
Engagement <= `tea_cap_deg`, decided by an exact predicate, at each EVALUATED tool
position -- now including positions along the bridges and the inter-chain links,
which is strictly more of the path than the other generators evaluate. It remains
NOT a continuous guarantee between evaluated positions, no certificate is produced
and none is returned.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from compas.geometry import Polygon

from compas_cgal import _coverage_2
from compas_cgal import _stock_2
from compas_cgal.engagement_rho_toolpath import DEGENERATE_LOOP_TOOL_RADII
from compas_cgal.engagement_rho_toolpath import DeclinedRegion
from compas_cgal.engagement_rho_toolpath import EmptyRhoGuideError
from compas_cgal.engagement_rho_toolpath import RhoToolpathResult
from compas_cgal.engagement_rho_toolpath import _as_guide_station
from compas_cgal.engagement_rho_toolpath import _declined_regions
from compas_cgal.engagement_rho_toolpath import _first_trochoidal_station
from compas_cgal.engagement_rho_toolpath import _loop_operation
from compas_cgal.engagement_rho_toolpath import _radius_step_admissible
from compas_cgal.engagement_rho_toolpath import _rho_stations
from compas_cgal.engagement_rho_toolpath import _RhoStation
from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import MAX_ADVANCE_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import POLYLINE_SAMPLES_PER_RADIAN
from compas_cgal.engagement_toolpath import UnavoidableEngagementWarning
from compas_cgal.engagement_toolpath import _guide_chains
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.engagement_toolpath import _line_operation
from compas_cgal.engagement_toolpath import _Regulation
from compas_cgal.engagement_toolpath import _station_is_admissible
from compas_cgal.engagement_toolpath import _tessellate
from compas_cgal.engagement_toolpath import _unit_tangent
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import _polygon_to_ccw_vertices

# Spacing of the evaluated positions along a straight move, in TOOL RADII. A
# capsule swept by a tool of radius r is covered by disks whose centres are 2r
# apart; sampling at r therefore puts a full tool diameter of overlap between
# consecutive evaluated positions, so a pocket of material narrower than the tool
# cannot sit unseen between two probes. Halving it again doubles the cost and can
# only find material a half-radius across that the disks already straddle.
LINK_PROBE_SPACING_TOOL_RADII = 1.0

# Evaluated positions on the shortest possible straight move. Two -- the endpoints
# -- is what the spacing rule degenerates to for a move shorter than one tool
# radius, and a move that short cannot hide anything the endpoints miss.
MIN_LINK_PROBES = 2

# Fraction of the engagement cap a STRAIGHT move is held to. A straight cut is a
# transfer through material the loops already cleared: its nominal engagement is
# zero, and whatever it does engage is material the stepover control never
# budgeted for. For a straight cut into a straight wall the engagement angle theta
# and the radial depth a_e are related exactly by a_e = r * (1 - cos(theta / 2)),
# so half the cap allows a link about a quarter of the regulated radial bite -- a
# generous allowance for clipping the corner of a cleared corridor, and a strictly
# tighter test than the one the machining circles face.
#
# WHY A STRAIGHT MOVE NEEDS ITS OWN THRESHOLD AT ALL. Checking links and bridges
# against the FULL cap was measured and is not enough: on `rect_12x8` it left
# eight straight moves engaging between half the cap and the cap, which is what a
# slot is. A circle at that engagement unloads once per turn; a straight move
# never unloads. ENGINEERING JUDGEMENT, matching
# `benchmarks.survey.SLOT_ENGAGEMENT_FRACTION` so the generator is held to the
# criterion it is measured by -- deliberately duplicated rather than imported,
# because `compas_cgal` must not depend on `benchmarks`.
STRAIGHT_MOVE_CAP_FRACTION = 0.5


@dataclass(frozen=True)
class IsolatedChain:
    """A chain the tool had to retract and plunge into, and the measurement that forced it.

    A generator that quietly plunges wherever it likes hides the cost of its own
    ordering. This is the record that stops that: it names a chain the ordering
    could NOT reach through cleared stock, with the engagement that refused the
    best link on offer, so "this pocket needs N entries" is evidence.

    Attributes:
        path_index: Index of the chain, matching its operations' `path_index`.
        entry: Tool-centre point the chain was plunged at, in world XY.
        best_link_length: Length of the shortest link that was tried, in model
            units, or ``0.0`` for the first chain, which has nothing to link from.
        reason: Why the link was refused, in one phrase.
    """

    path_index: int
    entry: Tuple[float, float]
    best_link_length: float
    reason: str


@dataclass
class OrderedToolpathResult(RhoToolpathResult):
    """A `RhoToolpathResult` that also says which chains could not be reached through cut stock.

    Attributes:
        isolated_chains: Every chain the ordering had to enter with a fresh plunge.
            One of these is unavoidable -- the first cut of the operation is always
            into virgin stock -- so a path with exactly one is optimal on this
            criterion, and every further one is a measured claim that no admissible
            link existed.
    """

    isolated_chains: Tuple[IsolatedChain, ...] = ()


@dataclass(frozen=True)
class _OrientedChain:
    """One skeleton chain in one of its two walk directions.

    Attributes:
        stations: The chain's stations in walk order, entries already resolved.
        origin: Index of the chain in the guide's own emission order, so a
            reversed and a forward copy of the same chain are known to be the same
            chain.
        reversed_walk: Whether this is the guide's order or its reverse.
        first: Index of the first station clearing the degeneracy floor.
        last: Index of the last station clearing it.
    """

    stations: List[_RhoStation]
    origin: int
    reversed_walk: bool
    first: int
    last: int

    @property
    def entry_point(self) -> Tuple[float, float]:
        """Tool-centre point this orientation would begin machining at."""
        return self.stations[self.first].entry

    @property
    def exit_point(self) -> Tuple[float, float]:
        """Tool-centre point this orientation would finish at."""
        return self.stations[self.last].entry


def _reversed_guide_chain(stations: Sequence[_GuideStation]) -> List[_GuideStation]:
    """The same guide stations walked the other way.

    The guide tangent is a DIRECTION and reverses with the walk; the radius and the
    turn direction are properties of the station and do not. Reversing before
    `_rho_stations` rather than after is what keeps the entry-radius construction
    consistent, because that construction reads the tangent and the local radius
    slope, both of which change sign together.

    Args:
        stations: The chain's guide stations in emission order.

    Returns:
        The stations in reverse order, each with its tangent reversed.
    """
    return [
        _GuideStation(
            cx=station.cx,
            cy=station.cy,
            radius=station.radius,
            clockwise=station.clockwise,
            tx=-station.tx,
            ty=-station.ty,
        )
        for station in reversed(stations)
    ]


def _orientations(guide_chain: Sequence[_GuideStation], origin: int, tool_radius: float) -> List[_OrientedChain]:
    """Both walk directions of one chain, or none when it carries no trochoid.

    Args:
        guide_chain: The chain's guide stations in emission order.
        origin: The chain's index in the guide's emission order.
        tool_radius: Tool radius, which fixes the degeneracy floor.

    Returns:
        Zero or two orientations.
    """
    floor = DEGENERATE_LOOP_TOOL_RADII * tool_radius
    oriented: List[_OrientedChain] = []
    for reversed_walk, source in ((False, list(guide_chain)), (True, _reversed_guide_chain(guide_chain))):
        stations = _rho_stations(source)
        first = _first_trochoidal_station(stations, tool_radius)
        if first is None:
            return []
        last = len(stations) - 1
        while last > first and stations[last].radius <= floor:
            last -= 1
        oriented.append(_OrientedChain(stations=stations, origin=origin, reversed_walk=reversed_walk, first=first, last=last))
    return oriented


def _straight_probe_positions(start: Tuple[float, float], end: Tuple[float, float], tool_radius: float) -> List[Tuple[float, float]]:
    """Tool-centre positions at which a straight move is decided.

    Spaced `LINK_PROBE_SPACING_TOOL_RADII` tool radii apart, endpoints included, so
    consecutive tool disks overlap by a full diameter. SAMPLED, NOT BOUNDED: the
    answer at each position is exact, and nothing here bounds what happens between
    two of them.

    Args:
        start: Where the move begins.
        end: Where it ends.
        tool_radius: Tool radius, which sets the spacing.

    Returns:
        The evaluated positions, start first and end last.
    """
    length = math.hypot(end[0] - start[0], end[1] - start[1])
    count = max(MIN_LINK_PROBES, int(math.ceil(length / (LINK_PROBE_SPACING_TOOL_RADII * tool_radius))) + 1)
    return [
        (
            start[0] + (end[0] - start[0]) * index / (count - 1),
            start[1] + (end[1] - start[1]) * index / (count - 1),
        )
        for index in range(count)
    ]


def _straight_move_ratio(regulation: _Regulation) -> float:
    """The exact rational surrogate a STRAIGHT move is decided against.

    BOUNDARY (`docs/exactness.md`): the transcendental threshold crosses into
    exact-land once, here, as the squared-chord surrogate `4*sin^2(theta/2)`, by
    the same `_stock_2.cap_chord_ratio` the cap itself uses. It is a stricter
    threshold than `regulation.cap_ratio`, never a looser one.

    Args:
        regulation: The validated parameters, for the caller's cap angle.

    Returns:
        The surrogate for `STRAIGHT_MOVE_CAP_FRACTION` of the cap.
    """
    return _stock_2.cap_chord_ratio(STRAIGHT_MOVE_CAP_FRACTION * regulation.cap_angle)


def _straight_move_is_admissible(
    stock: Stock,
    start: Tuple[float, float],
    end: Tuple[float, float],
    regulation: _Regulation,
    centre_domain: "_coverage_2.ExactRegion2",
) -> bool:
    """Whether a straight move is LEGAL and stays a link rather than a slot.

    TWO CONDITIONS, AND CONTAINMENT COMES FIRST. A move whose every position is
    under the engagement threshold can still be illegal: engagement says how much
    material the cutter meets, not whether the cutter may be there at all. Between
    two chains the straight line can leave the region a tool of this radius may be
    centred in -- measured on `L_shape`, where an engagement-only test passed a
    link that gouged. `ExactRegion2.contains` is exact point location in the exact
    reachable centre domain, so legality is decided by a predicate rather than by
    a distance with a tolerance.

    The same exact `cap_exceeded` predicate the machining circles are decided by,
    applied to the move the other generators never checked, and at the STRICTER
    `STRAIGHT_MOVE_CAP_FRACTION` threshold -- a straight move that engages like a
    machining circle is a slot, because it never unloads. ``gap_close_ratio``
    stays at zero: gap closure is the pessimism a CONTINUOUS certifier needs, and
    this makes no between-position claim to bridge.

    Args:
        stock: The current stock (unmodified by this call).
        start: Where the move begins.
        end: Where it ends.
        regulation: The validated parameters.
        centre_domain: Exact region a cutter of this radius may be centred in.

    Returns:
        ``True`` if every evaluated position is inside the centre domain and none
        reports the straight-move threshold exceeded.
    """
    raw = stock.raw
    ratio = _straight_move_ratio(regulation)
    for probe_x, probe_y in _straight_probe_positions(start, end, regulation.tool_radius):
        if not centre_domain.contains(probe_x, probe_y):
            return False
        _total, _max_run, exceeded = _stock_2.engagement_at(raw, probe_x, probe_y, regulation.tool_radius, ratio, 0.0)
        if exceeded:
            return False
    return True


def _largest_admissible_advance(
    stock: Stock,
    stations: List[_RhoStation],
    origin_index: int,
    window_end: int,
    regulation: _Regulation,
    centre_domain: "_coverage_2.ExactRegion2",
) -> Tuple[int, bool]:
    """SCAN the window from the furthest end for the largest advance whose LOOP AND BRIDGE comply.

    This is `compas_cgal.engagement_rho_toolpath._largest_admissible_advance` with
    one condition added: the bridge that reaches the candidate is decided by the
    same exact predicate as the candidate's machining circle. The older generators
    check only the circle, which is safe while the walk advances OUT of cleared
    material and unsafe the moment it advances INTO virgin material -- which is
    exactly what a chain entered from its wide end does.

    Still a scan and never a bisection: engagement is not monotone in the advance,
    so the admissible advances are not an up-set in the station index, and only a
    walk from the furthest candidate downward returns the furthest admissible one
    under an arbitrary pass/fail pattern.

    Args:
        stock: The current stock (unmodified by this call).
        stations: The chain's ordered stations.
        origin_index: Index of the last accepted station.
        window_end: Highest candidate index the search may consider (inclusive).
        regulation: The validated parameters.

    Returns:
        ``(index, forced)``: the chosen station index, and whether it was taken
        despite no candidate being admissible.
    """
    origin = stations[origin_index]
    for index in range(window_end, origin_index, -1):
        candidate = stations[index]
        if not _radius_step_admissible(origin, candidate, regulation.tool_radius):
            continue
        advance = _unit_tangent(origin.cx, origin.cy, candidate.cx, candidate.cy)
        if not _station_is_admissible(stock, _as_guide_station(candidate), advance, regulation.tool_radius, regulation.cap_ratio):
            continue
        if _straight_move_is_admissible(stock, origin.entry, candidate.entry, regulation, centre_domain):
            return index, False
    minimum = origin_index + 1
    for index in range(minimum, window_end + 1):
        if _radius_step_admissible(origin, stations[index], regulation.tool_radius):
            return index, True
    return minimum, True


def _walk_chain(
    stock: Stock,
    chain: _OrientedChain,
    path_index: int,
    regulation: _Regulation,
    centre_domain: "_coverage_2.ExactRegion2",
    cut_z: float,
    operations: List[ToolpathOperation],
) -> Tuple[Tuple[float, float], int]:
    """Emit one oriented chain's machining circles and bridges, depleting as it goes.

    Emits NO plunge and NO retract: entering and leaving is the caller's decision,
    because whether this chain can be reached through cleared stock is a property
    of the ordering rather than of the chain.

    Args:
        stock: The depleting stock, mutated in place.
        chain: The oriented chain to machine.
        path_index: Path index stamped on this chain's operations.
        regulation: The validated parameters.
        cut_z: Cutting-plane height.
        operations: Output stream, appended to in place.

    Returns:
        ``(exit_point, forced_advances)``.
    """
    stations = chain.stations
    forced_advances = 0
    index = chain.first
    while True:
        station = stations[index]
        entry_x, entry_y = station.entry
        operations.append(_loop_operation(station, cut_z, path_index))
        stock.subtract_arc_sweep_local(station.cx, station.cy, entry_x, entry_y, entry_x, entry_y, station.clockwise, regulation.tool_radius)
        if index == chain.last:
            return (entry_x, entry_y), forced_advances
        next_index, forced = _largest_admissible_advance(
            stock,
            stations,
            index,
            min(index + regulation.window, chain.last),
            regulation,
            centre_domain,
        )
        forced_advances += int(forced)
        next_x, next_y = stations[next_index].entry
        operations.append(_line_operation((entry_x, entry_y), (next_x, next_y), cut_z, cut_z, OperationType.CUT, path_index))
        stock.subtract_capsule_quad(entry_x, entry_y, next_x, next_y, regulation.tool_radius)
        index = next_index


def _choose_next(
    stock: Stock,
    pending: List[List[_OrientedChain]],
    current_exit: Optional[Tuple[float, float]],
    regulation: _Regulation,
    centre_domain: "_coverage_2.ExactRegion2",
) -> Tuple[int, _OrientedChain, bool, float]:
    """Pick the next chain and the end to enter it from.

    Every remaining chain is offered in BOTH orientations. An orientation whose
    link from the current position is admissible under the exact cap predicate is
    preferred over any that is not, and among equals the shortest link wins --
    length breaks ties, it never overrides the predicate.

    Args:
        stock: The current stock (unmodified by this call).
        pending: Remaining chains, each as its list of orientations.
        current_exit: Where the tool is, or ``None`` before the first chain.
        regulation: The validated parameters.

    Returns:
        ``(pending_index, orientation, linked, link_length)``, where *linked* says
        whether the link may be cut at depth rather than flown at clearance.
    """
    best: Optional[Tuple[bool, float, int, _OrientedChain]] = None
    for pending_index, orientations in enumerate(pending):
        for orientation in orientations:
            entry = orientation.entry_point
            if current_exit is None:
                length = 0.0
                linked = False
            else:
                length = math.hypot(entry[0] - current_exit[0], entry[1] - current_exit[1])
                linked = _straight_move_is_admissible(stock, current_exit, entry, regulation, centre_domain)
            key = (not linked, length, pending_index, orientation)
            if best is None or key[:3] < best[:3]:
                best = key
    assert best is not None  # `pending` is non-empty at every call site
    not_linked, length, pending_index, orientation = best
    return pending_index, orientation, not not_linked, length


def chain_ordered_toolpath(
    polygon: Polygon,
    tool_diameter: float,
    tea_cap_deg: float,
    *,
    holes: Optional[List[Polygon]] = None,
    guide_step_tool_diameters: float = GUIDE_STEP_TOOL_DIAMETERS,
    max_advance_tool_diameters: float = MAX_ADVANCE_TOOL_DIAMETERS,
    radial_clearance: Optional[float] = None,
    climb: bool = True,
    cut_z: float = 0.0,
    clearance_z: Optional[float] = None,
    max_passes: int = 1000,
    samples_per_radian: float = POLYLINE_SAMPLES_PER_RADIAN,
) -> OrderedToolpathResult:
    """Rho-regulated pocketing that reaches each chain through stock it has already cut.

    Same guide, same degeneracy floor, same external-tangent linking and same exact
    cap predicate as
    `compas_cgal.engagement_rho_toolpath.rho_regulated_toolpath`. What changes is
    the sequence: chains are ordered and ORIENTED so each one after the first begins
    where the tool already is, and the move into it is cut at depth whenever every
    evaluated position on it is under the cap. Only a chain no admissible link
    reaches costs a retract, a rapid and a plunge, and each of those is named in the
    result.

    The bridges inside a chain are cap-checked too, which the other generators do
    not do. That is not gold-plating: a chain entered at its wide end advances INTO
    virgin material, and an unchecked bridge there is a slot.

    WHAT THIS GUARANTEES, EXACTLY: engagement <= *tea_cap_deg*, decided by an exact
    predicate, at each EVALUATED tool position -- the loop probe rings, the bridges,
    and the inter-chain links. Not a continuous guarantee between evaluated
    positions; no certificate is produced or returned.

    Args:
        polygon: Outer pocket boundary `Polygon` in the world XY plane.
        tool_diameter: Tool diameter; the tool radius is half of this.
        tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``.
        holes: Optional island polygons strictly inside *polygon*.
        guide_step_tool_diameters: Guide station spacing in tool diameters.
        max_advance_tool_diameters: Largest advance considered, in tool diameters.
        radial_clearance: Safety clearance subtracted from each loop's available
            radius. Defaults to ``RADIAL_CLEARANCE_FRACTION * tool_diameter``.
        climb: ``True`` for climb milling (clockwise loops).
        cut_z: Z-height of the cutting plane.
        clearance_z: Z-height for rapid travel between chains.
        max_passes: Maximum number of skeleton chains the guide may emit.
        samples_per_radian: Tessellation density of the returned polyline.

    Returns:
        OrderedToolpathResult: The operation stream, the visualisation polyline,
        `declined_regions` for the material below the degeneracy floor, and
        `isolated_chains` for every chain that still needed its own plunge.

    Raises:
        InvalidEngagementCapDegreesError: If *tea_cap_deg* is not in ``(0, 180]``.
        NonPositiveToolDiameterError: If *tool_diameter* is not strictly positive.
        InvalidGuideResolutionError: If the guide step or advance bound leaves no
            integer window to scan.
        InvalidClearanceHeightError: If *clearance_z* is not above *cut_z*.
        EmptyGuideError: If the pocket admits no gouge-free trochoid at this tool.
        EmptyRhoGuideError: If no chain carries a non-degenerate machining circle.
        ImplausibleClearanceSlopeError: If a chain's radius grows faster along its
            own arc length than a 1-Lipschitz clearance function can.

    Warns:
        UnavoidableEngagementWarning: When advances had to be forced past a
            refusing predicate, chains needed their own plunge, or guide stations
            were declined for being below the degeneracy floor.
    """
    regulation = _Regulation.build(
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        guide_step_tool_diameters=guide_step_tool_diameters,
        max_advance_tool_diameters=max_advance_tool_diameters,
        radial_clearance=radial_clearance,
        cut_z=cut_z,
        clearance_z=clearance_z,
    )
    guide_chains = _guide_chains(
        polygon,
        tool_diameter,
        regulation.guide_step,
        regulation.radial_clearance,
        climb,
        max_passes,
        holes,
    )

    pending: List[List[_OrientedChain]] = []
    declined: List[DeclinedRegion] = []
    for origin, guide_chain in enumerate(guide_chains):
        orientations = _orientations(guide_chain, origin, regulation.tool_radius)
        if orientations:
            pending.append(orientations)
        else:
            declined.extend(_declined_regions(_rho_stations(guide_chain), origin, regulation.tool_radius))
    if not pending:
        raise EmptyRhoGuideError(
            f"The straight-skeleton guide produced {len(guide_chains)} chain(s) for tool_diameter={tool_diameter!r}, but no station on any of them "
            f"has a gouge-free radius above the {DEGENERATE_LOOP_TOOL_RADII * regulation.tool_radius!r} degeneracy floor."
        )

    # The exact region a cutter of this radius may be centred in. Built once: it
    # depends only on the boundary and the tool, never on the depleting stock.
    centre_domain = _coverage_2.ReachableDomain2(
        _polygon_to_ccw_vertices(polygon),
        [_polygon_to_ccw_vertices(hole) for hole in (holes or [])],
        regulation.tool_radius,
    ).center_domain()

    stock = Stock(polygon, holes=holes)
    operations: List[ToolpathOperation] = []
    isolated: List[IsolatedChain] = []
    forced_advances = 0
    current_exit: Optional[Tuple[float, float]] = None
    path_index = 0

    while pending:
        pending_index, chain, linked, link_length = _choose_next(stock, pending, current_exit, regulation, centre_domain)
        entry = chain.entry_point
        if linked and current_exit is not None:
            # Cut into the next chain at depth. Every evaluated position on this
            # move is under the cap, so it is a machining move like any other and
            # it removes what it passes through.
            operations.append(_line_operation(current_exit, entry, cut_z, cut_z, OperationType.CUT, path_index))
            stock.subtract_capsule_quad(current_exit[0], current_exit[1], entry[0], entry[1], regulation.tool_radius)
        else:
            if current_exit is not None:
                operations.append(_line_operation(current_exit, current_exit, cut_z, regulation.clearance_z, OperationType.RETRACT, path_index))
                operations.append(_line_operation(current_exit, entry, regulation.clearance_z, regulation.clearance_z, OperationType.LINK, path_index))
            operations.append(_line_operation(entry, entry, regulation.clearance_z, cut_z, OperationType.PLUNGE, path_index))
            isolated.append(
                IsolatedChain(
                    path_index=path_index,
                    entry=entry,
                    best_link_length=link_length,
                    reason=(
                        "first cut of the operation, which is always into virgin stock"
                        if current_exit is None
                        else "no orientation of any remaining chain had a link that stays inside the centre domain and under the straight-move threshold"
                    ),
                )
            )
        current_exit, chain_forced = _walk_chain(stock, chain, path_index, regulation, centre_domain, cut_z, operations)
        forced_advances += chain_forced
        declined.extend(_declined_regions(chain.stations, path_index, regulation.tool_radius))
        pending.pop(pending_index)
        path_index += 1

    operations.append(_line_operation(current_exit, current_exit, cut_z, regulation.clearance_z, OperationType.RETRACT, path_index - 1))

    if forced_advances or declined or len(isolated) > 1:
        warnings.warn(
            f"{forced_advances} advance(s) were taken past a refusing exact cap predicate at "
            f"tea_cap_deg={tea_cap_deg}. {len(isolated)} chain(s) needed their own plunge into virgin stock; one of "
            "those is the operation's first cut and unavoidable, and any others are named in "
            f"`OrderedToolpathResult.isolated_chains` with the link that was refused. {len(declined)} run(s) of guide "
            "were declined for being below the degeneracy floor and are named in "
            "`OrderedToolpathResult.declined_regions`.",
            UnavoidableEngagementWarning,
            stacklevel=2,
        )

    return OrderedToolpathResult(
        operations=operations,
        polyline=_tessellate(operations, samples_per_radian),
        declined_regions=tuple(declined),
        isolated_chains=tuple(isolated),
    )
