"""Spiral-ramp entry: open the hole before the first full-radius loop, in the cut plane.

STATUS: A PROOF ARTIFACT AND A PARTIAL MITIGATION, NOT THE ENTRY STRATEGY. This
module exists because building it is what established the 240-degree floor below,
and because the concentric-tangent argument it forced -- concentric circles have
no common tangent, so ramp turns can only be linked radially, at a right angle,
through uncut material -- is a geometric finding worth keeping. It is NOT
registered in `benchmarks.gate.GATE_GENERATORS` and should not be: on the corpus
it fails six criteria where `rho_regulated_toolpath` fails four, and by the
structural cap below it can never reach the one criterion it was built to fix.
Use it to measure the entry, not to machine a pocket.

`compas_cgal.engagement_rho_toolpath.rho_regulated_toolpath` plunges at the first
machining circle's ENTRY POINT and then runs that circle at its full gouge-free
radius. Entering solid stock that way is a full-immersion cut -- 360 degrees, for
any generator, at any non-degenerate radius -- and the engagement STEP from it to
the motion after it is the largest number in the whole CUT group. This module is
the same walk with a different entry: plunge at the machining circle's CENTRE, then
climb to the full radius through a ramp of concentric circles, each chosen by the
same exact `cap_exceeded` predicate.

WHY A TRUE HELIX IS NOT WHAT THIS EMITS
---------------------------------------
The obvious answer is a helical ramp: orbit while descending in Z, so each turn
takes a shallow axial bite. THE MEASUREMENT MODEL CANNOT REPRESENT IT.
`benchmarks.depletion._replay_kind` classifies operations against ONE inferred cut
plane and raises `UnreplayableOperationError` on any move with both Z change and
XY travel -- "ramped 3D cutting is outside the cut-plane depletion model". A
helical arc emitted as a `Circle` above the cut plane is classified
`ReplayKind.RAPID`: it would remove nothing in the replay, so the hole it bored
would be invisible to the coverage grid and the first loop would still measure 360
degrees. A helix is therefore not a thing this framework can generate OR score,
and pretending otherwise would produce a path whose reported engagement is a
fiction.

What CAN be done in one plane is the in-plane equivalent of a helical bore: open
the hole with concentric full-width turns instead of shallow ones. That is what
this module does, and it is a genuinely different motion -- more full-width cuts,
no axial ramping -- so it is named for what it is.

THE EXACT BOUND THIS BUYS, AND WHERE IT STOPS
---------------------------------------------
With the plunge at the CENTRE, the void is a disk of radius ``r`` concentric with
every ramp circle, so a tool centred at radius ``rho`` has the void subtend
``2*acos((rho^2 + r^2 - V^2) / (2*rho*r))`` at its centre. At ``V = r`` that
collapses to

    peak engagement = 360 - 2*acos(rho / 2r)

which is INCREASING in rho, so over the non-degenerate range ``rho > r`` its
infimum is at ``rho -> r``:

    360 - 2*acos(1/2) = 360 - 120 = 240 degrees.

**No non-degenerate first loop after a single plunge can engage less than 240
degrees.** Measured against `_stock_2.engagement_at` on `rect_12x8` at ten radii,
the closed form and the exact predicate agree to the last reported digit; the
first admissible rung, 1.048 r, measures 243.20 degrees against the advance-only
generator's 360.00.

The consequence for the engagement STEP is a sharp bound rather than a fix. The
motion after the entry respects the cap, so the step is at least ``240 - cap``,
and at a 120 degree cap that infimum is EXACTLY 120 -- reachable only in the
double limit of a first loop at the degeneracy boundary and a second loop exactly
at the cap, both of which are open. Measured on `rect_12x8` the ramp descends
243.20 -> 118.56 -> 113.91 -> 119.27 -> 116.39 -> 114.23 -> 104.10, so the largest
step is 124.64 degrees against the plunge entry's 329.49. Better by a factor of
2.6, and still not under the criterion, because it cannot be.

THE IN-PLANE APPROACH IS CAPPED AT 180 DEGREES, WHATEVER THE RADIUS
-------------------------------------------------------------------
The 240 figure is the floor over the NON-DEGENERATE range. Taken over ALL
positive radii the same expression gives

    lim (rho -> 0) of 360 - 2*acos(rho / 2r)  =  360 - 2*acos(0)  =  180 degrees,

so no concentric turn after a plunge engages below 180 degrees at ANY radius,
degenerate or not. Deliberately boring a tiny first turn to widen the hole does
not rescue the construction; it is not a limitation of this implementation but of
the geometry, and it means a concentric ramp cannot meet any cap below 180.

WHAT AN ENTRY WOULD HAVE TO OPEN
--------------------------------
Rearranging the same law of cosines, a void of radius V concentric with a first
loop of radius rho meets a cap of theta only while

    V >= sqrt(rho^2 + r^2 - 2*rho*r*cos((360 - theta) / 2))

At a 120 degree cap that is 1.7321*r for rho = r, 2.179*r for rho = 1.5*r, and
3.604*r for the spine's rho = 2.998*r. A plunge opens V = r. It is short by a
factor of 1.73 in the very best case, which is the whole problem in one number:
the hole has to be OPENED, not merely started, and nothing in the cut-plane
vocabulary opens it without a full-immersion motion of some kind.

WHAT IT COSTS
-------------
Concentric circles have no common tangent -- a line at distance ``rho`` from the
shared centre cannot also be at distance ``rho'`` -- so consecutive ramp turns can
only be linked by a RADIAL move, which meets both circles at a right angle. Every
ramp step therefore costs two tangent breaks, and the ramp needs seven turns on
`rect_12x8`. That is the trade this module exists to measure, and
`spiral_entry_toolpath` reports it rather than burying it: the tangent-continuous
alternative would have to advance the ramp ALONG the guide, which grows the radius
at the clearance slope and leaves an uncut wedge behind the ramp that only a
second pass can clear.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Tuple

from compas.geometry import Polygon

from compas_cgal.engagement_rho_toolpath import DEGENERATE_LOOP_TOOL_RADII
from compas_cgal.engagement_rho_toolpath import DeclinedRegion
from compas_cgal.engagement_rho_toolpath import RhoToolpathResult
from compas_cgal.engagement_rho_toolpath import _as_guide_station
from compas_cgal.engagement_rho_toolpath import _declined_regions
from compas_cgal.engagement_rho_toolpath import _first_trochoidal_station
from compas_cgal.engagement_rho_toolpath import _largest_admissible_advance
from compas_cgal.engagement_rho_toolpath import _loop_operation
from compas_cgal.engagement_rho_toolpath import _rho_stations
from compas_cgal.engagement_rho_toolpath import _RhoStation
from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import MAX_ADVANCE_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import POLYLINE_SAMPLES_PER_RADIAN
from compas_cgal.engagement_toolpath import UnavoidableEngagementWarning
from compas_cgal.engagement_toolpath import _guide_chains
from compas_cgal.engagement_toolpath import _line_operation
from compas_cgal.engagement_toolpath import _Regulation
from compas_cgal.engagement_toolpath import _station_is_admissible
from compas_cgal.engagement_toolpath import _tessellate
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

# Hard ceiling on ramp turns at one entry, so a ramp that cannot make progress
# raises instead of spinning. Every turn must strictly enlarge the emitted radius
# on a grid whose rungs are one guide step apart, so the reachable count is
# bounded by the station's radius over the guide step -- forty at the shipped
# defaults on the widest gate pocket. This is twice that: a ramp that reaches it
# is not converging and the run is a defect, not a slow success.
MAX_RAMP_TURNS = 80


class RampCannotProgressError(RuntimeError):
    """A spiral entry ramp reached its turn ceiling without arriving at the station's full radius."""


@dataclass(frozen=True)
class _RampOutcome:
    """What one entry ramp emitted.

    Attributes:
        turns: Machining circles the ramp emitted, the final full-radius one
            included.
        entry: Tool-centre point the ramp plunged at -- the circle's CENTRE, which
            is the whole point of the construction.
        exit: Tool-centre point the ramp finished on, where the chain walk takes
            over.
        forced_turns: Turns emitted at radii the exact cap predicate refuses. The
            first is expected and bounded below by 240 degrees; more than one means
            the ramp is climbing through material the cap cannot reach.
        entry_peak_deg: REPORTING ONLY. Peak engagement the first turn measured,
            so a caller can compare it against the 360 degrees a plunge entry
            gives without replaying the path.
    """

    turns: int
    entry: Tuple[float, float]
    exit: Tuple[float, float]
    forced_turns: int
    entry_peak_deg: float


def _ramp_rungs(station: _RhoStation, regulation: _Regulation) -> List[float]:
    """Candidate ramp radii at one station, LARGEST FIRST.

    The rungs sit one guide step apart below the station's gouge-free radius and
    stop at the degeneracy floor, so no rung the ramp can pick is a bore. Rung 0
    is the station's own radius bit-for-bit, so arriving at it is an identity test
    rather than a proximity one.

    Args:
        station: The entry station.
        regulation: The validated parameters.

    Returns:
        The candidate radii in descending order.
    """
    floor = DEGENERATE_LOOP_TOOL_RADII * regulation.tool_radius
    rungs: List[float] = []
    radius = station.radius
    while radius > floor:
        rungs.append(radius)
        radius -= regulation.guide_step
    return rungs


def _emit_ramp(
    stock: Stock,
    station: _RhoStation,
    path_index: int,
    regulation: _Regulation,
    cut_z: float,
    operations: List[ToolpathOperation],
) -> _RampOutcome:
    """Plunge at the circle's centre and climb to its full radius through concentric turns.

    Each turn takes the LARGEST rung the exact predicate accepts, scanning
    downward -- the same order argument the advance search uses, and for the same
    reason: engagement is not monotone in the radius, so the admissible rungs are
    not an up-set and only a scan from the top returns the largest admissible one.

    The first turn is expected to be refused at every rung: with only a plunge
    hole open, the smallest non-degenerate radius still engages 240 degrees (see
    the module docstring's bound). It is emitted, counted, and reported rather
    than replaced by something that would misreport itself.

    Args:
        stock: The depleting stock, mutated in place.
        station: The entry station, carrying its full gouge-free radius.
        path_index: Path index stamped on this chain's operations.
        regulation: The validated parameters.
        cut_z: Cutting-plane height.
        operations: Output stream, appended to in place.

    Returns:
        What the ramp emitted.

    Raises:
        RampCannotProgressError: If the ramp did not arrive at the station's full
            radius within `MAX_RAMP_TURNS`.
    """
    rungs = _ramp_rungs(station, regulation)
    # Plunge at the CENTRE. This is the one geometric difference from
    # `engagement_rho_toolpath`, and it is what makes the hole concentric with
    # every ramp turn instead of sitting on the rim of the first one.
    operations.append(_line_operation((station.cx, station.cy), (station.cx, station.cy), regulation.clearance_z, cut_z, OperationType.PLUNGE, path_index))
    stock.subtract_capsule_quad(station.cx, station.cy, station.cx, station.cy, regulation.tool_radius)

    emitted = 0.0
    forced_turns = 0
    entry_peak = 0.0
    exit_point = (station.cx, station.cy)
    for turn in range(MAX_RAMP_TURNS):
        chosen: Optional[float] = None
        for radius in rungs:
            if radius <= emitted:
                continue
            candidate = _turn_station(station, radius)
            if _station_is_admissible(stock, _as_guide_station(candidate), candidate.entry_tangent, regulation.tool_radius, regulation.cap_ratio):
                chosen = radius
                break
        forced = chosen is None
        if forced:
            # Nothing complies. Take the smallest non-degenerate rung: it is the
            # gentlest motion that is still a trochoid rather than a bore, and the
            # bound in the module docstring says 240 degrees is the floor here.
            reachable = [radius for radius in rungs if radius > emitted]
            if not reachable:
                break
            chosen = min(reachable)
            forced_turns += 1
        candidate = _turn_station(station, chosen)
        loop_entry = candidate.entry
        if emitted > 0.0:
            # Concentric turns admit no common tangent, so the link is RADIAL and
            # meets both circles at a right angle. Two tangent breaks per step is
            # the measured cost of this entry; the module docstring carries it.
            operations.append(_line_operation(exit_point, loop_entry, cut_z, cut_z, OperationType.CUT, path_index))
            stock.subtract_capsule_quad(exit_point[0], exit_point[1], loop_entry[0], loop_entry[1], regulation.tool_radius)
        if turn == 0:
            entry_peak = _peak_engagement_deg(stock, candidate, regulation)
        operations.append(_loop_operation(candidate, cut_z, path_index))
        stock.subtract_arc_sweep_local(
            candidate.cx,
            candidate.cy,
            loop_entry[0],
            loop_entry[1],
            loop_entry[0],
            loop_entry[1],
            candidate.clockwise,
            regulation.tool_radius,
        )
        emitted = chosen
        exit_point = loop_entry
        if chosen == station.radius:
            return _RampOutcome(
                turns=turn + 1,
                entry=(station.cx, station.cy),
                exit=exit_point,
                forced_turns=forced_turns,
                entry_peak_deg=entry_peak,
            )
    raise RampCannotProgressError(
        f"The spiral entry ramp at ({station.cx!r}, {station.cy!r}) did not reach the station's gouge-free radius "
        f"{station.radius!r} within {MAX_RAMP_TURNS} turns. Every turn must strictly enlarge the emitted radius on a "
        f"grid of {regulation.guide_step!r}-wide rungs, so this means the ramp stopped enlarging it."
    )


def _turn_station(station: _RhoStation, radius: float) -> _RhoStation:
    """The entry station at one ramp radius, keeping its centre and turn direction."""
    return _RhoStation(
        cx=station.cx,
        cy=station.cy,
        radius=radius,
        clockwise=station.clockwise,
        tx=station.tx,
        ty=station.ty,
        wx=station.wx,
        wy=station.wy,
    )


def _peak_engagement_deg(stock: Stock, station: _RhoStation, regulation: _Regulation) -> float:
    """Largest engaged-run angle REPORTED over one circle's probe ring, in degrees.

    REPORTING, NOT DECIDING (`docs/exactness.md`, the deciding/reporting split).
    Nothing is gated on this number: it travels on `_RampOutcome` so a caller can
    see what the entry cost without replaying the path.
    """
    from compas_cgal import _stock_2
    from compas_cgal.engagement_toolpath import _probe_positions

    peak = 0.0
    for probe_x, probe_y in _probe_positions(_as_guide_station(station), station.entry_tangent):
        _total, max_run, _exceeded = _stock_2.engagement_at(stock.raw, probe_x, probe_y, regulation.tool_radius, regulation.cap_ratio, 0.0)
        peak = max(peak, max_run)
    return math.degrees(peak)


def _machine_chain(
    stock: Stock,
    stations: List[_RhoStation],
    path_index: int,
    regulation: _Regulation,
    cut_z: float,
    operations: List[ToolpathOperation],
) -> Optional[Tuple[Tuple[float, float], Tuple[float, float], int, float]]:
    """Ramp in at the first trochoidal station, then walk the chain as the rho generator does.

    Args:
        stock: The depleting stock, mutated in place.
        stations: The chain's ordered stations.
        path_index: Path index stamped on this chain's operations.
        regulation: The validated parameters.
        cut_z: Cutting-plane height.
        operations: Output stream, appended to in place.

    Returns:
        ``(entry, exit, forced, entry_peak_deg)``, or ``None`` when no station on
        the chain clears the degeneracy floor.

    Raises:
        RampCannotProgressError: If the entry ramp did not converge.
    """
    first = _first_trochoidal_station(stations, regulation.tool_radius)
    if first is None:
        return None
    last = len(stations) - 1
    while last > first and stations[last].radius <= DEGENERATE_LOOP_TOOL_RADII * regulation.tool_radius:
        last -= 1

    ramp = _emit_ramp(stock, stations[first], path_index, regulation, cut_z, operations)
    forced = ramp.forced_turns
    exit_x, exit_y = ramp.exit

    index = first
    while index != last:
        next_index, advance_forced = _largest_admissible_advance(
            stock,
            stations,
            index,
            min(index + regulation.window, last),
            regulation,
        )
        forced += int(advance_forced)
        next_x, next_y = stations[next_index].entry
        operations.append(_line_operation((exit_x, exit_y), (next_x, next_y), cut_z, cut_z, OperationType.CUT, path_index))
        stock.subtract_capsule_quad(exit_x, exit_y, next_x, next_y, regulation.tool_radius)
        station = stations[next_index]
        operations.append(_loop_operation(station, cut_z, path_index))
        stock.subtract_arc_sweep_local(station.cx, station.cy, next_x, next_y, next_x, next_y, station.clockwise, regulation.tool_radius)
        exit_x, exit_y = next_x, next_y
        index = next_index

    operations.append(_line_operation((exit_x, exit_y), (exit_x, exit_y), cut_z, regulation.clearance_z, OperationType.RETRACT, path_index))
    return ramp.entry, (exit_x, exit_y), forced, ramp.entry_peak_deg


def spiral_entry_toolpath(
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
) -> RhoToolpathResult:
    """Rho-regulated pocketing whose chains ramp in instead of plunging into a full-radius loop.

    Identical to `compas_cgal.engagement_rho_toolpath.rho_regulated_toolpath` in
    every respect but the entry: the plunge goes at the first machining circle's
    CENTRE rather than on its rim, and the radius climbs to the station's full
    gouge-free value through concentric turns chosen by the exact cap predicate.

    WHAT THIS GUARANTEES, EXACTLY: engagement <= *tea_cap_deg*, decided by an exact
    predicate, at each EVALUATED tool position -- unchanged, and no certificate is
    produced or returned. What changes is the ENTRY: measured on `rect_12x8` with a
    2 mm tool at a 120 degree cap, the first machining circle engages 243.20 degrees
    instead of 360.00, and the largest engagement step falls from 329.49 to 124.64
    degrees.

    IT DOES NOT REACH A 120 DEGREE STEP, AND CANNOT. With a plunge hole of radius
    ``r``, a loop at radius ``rho`` engages ``360 - 2*acos(rho / 2r)``, which is
    increasing in rho and equals 240 degrees at ``rho = r`` -- the degeneracy
    boundary. Every non-degenerate first loop therefore engages at least 240
    degrees, the motion after it respects the cap, and the step is at least
    ``240 - cap``. At a 120 degree cap that infimum is exactly the criterion, and
    it is approached only from above.

    IT COSTS TANGENT BREAKS. Concentric circles have no common tangent, so ramp
    turns are linked radially, at a right angle, two breaks per turn.

    Args:
        polygon: Outer pocket boundary `Polygon` in the world XY plane.
        tool_diameter: Tool diameter; the tool radius is half of this.
        tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``.
        holes: Optional island polygons strictly inside *polygon*.
        guide_step_tool_diameters: Guide station spacing in tool diameters, and
            the spacing of the ramp's radius rungs.
        max_advance_tool_diameters: Largest advance considered, in tool diameters.
        radial_clearance: Safety clearance subtracted from each loop's available
            radius. Defaults to ``RADIAL_CLEARANCE_FRACTION * tool_diameter``.
        climb: ``True`` for climb milling (clockwise loops).
        cut_z: Z-height of the cutting plane.
        clearance_z: Z-height for rapid travel between chains.
        max_passes: Maximum number of skeleton chains the guide may emit.
        samples_per_radian: Tessellation density of the returned polyline.

    Returns:
        RhoToolpathResult: The operation stream, the visualisation polyline, and
        `RhoToolpathResult.declined_regions` naming every stretch of guide left
        unmachined because no circle there would be a trochoid.

    Raises:
        InvalidEngagementCapDegreesError: If *tea_cap_deg* is not in ``(0, 180]``.
        NonPositiveToolDiameterError: If *tool_diameter* is not strictly positive.
        InvalidGuideResolutionError: If the guide step or advance bound leaves no
            integer window to scan.
        InvalidClearanceHeightError: If *clearance_z* is not above *cut_z*.
        EmptyGuideError: If the pocket admits no gouge-free trochoid at this tool.
        RampCannotProgressError: If an entry ramp did not reach its station's full
            radius within `MAX_RAMP_TURNS`.
        ImplausibleClearanceSlopeError: If a chain's radius grows faster along its
            own arc length than a 1-Lipschitz clearance function can.

    Warns:
        UnavoidableEngagementWarning: When ramp turns had to be emitted over the
            cap, or guide stations were declined for being below the degeneracy
            floor (never silently).
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
    chains = _guide_chains(
        polygon,
        tool_diameter,
        regulation.guide_step,
        regulation.radial_clearance,
        climb,
        max_passes,
        holes,
    )

    stock = Stock(polygon, holes=holes)
    operations: List[ToolpathOperation] = []
    forced = 0
    declined: List[DeclinedRegion] = []
    entry_peaks: List[float] = []
    last_exit: Optional[Tuple[float, float]] = None
    path_index = 0

    for guide_chain in chains:
        stations = _rho_stations(guide_chain)
        pending = len(operations)
        outcome = _machine_chain(stock, stations, path_index, regulation, cut_z, operations)
        declined.extend(_declined_regions(stations, path_index, regulation.tool_radius))
        if outcome is None:
            continue
        entry, chain_exit, chain_forced, entry_peak = outcome
        if last_exit is not None:
            operations.insert(pending, _line_operation(last_exit, entry, regulation.clearance_z, regulation.clearance_z, OperationType.LINK, path_index))
        forced += chain_forced
        entry_peaks.append(entry_peak)
        last_exit = chain_exit
        path_index += 1

    if last_exit is None:
        from compas_cgal.engagement_rho_toolpath import EmptyRhoGuideError

        raise EmptyRhoGuideError(
            f"The straight-skeleton guide produced {len(chains)} chain(s) for tool_diameter={tool_diameter!r}, but no station on any of them "
            f"has a gouge-free radius above the {DEGENERATE_LOOP_TOOL_RADII * regulation.tool_radius!r} degeneracy floor."
        )

    if forced or declined:
        worst = max(entry_peaks) if entry_peaks else 0.0
        warnings.warn(
            f"{forced} machining circle(s) were emitted at positions where the exact cap predicate reports "
            f"tea_cap_deg={tea_cap_deg} exceeded. The worst chain entry measured {worst:.2f} degrees; the floor for a "
            "non-degenerate first loop after a single plunge is 240 degrees, so an entry at or near it is the "
            f"construction working rather than failing. {len(declined)} run(s) of guide were declined for being below "
            "the degeneracy floor and are returned as `RhoToolpathResult.declined_regions`.",
            UnavoidableEngagementWarning,
            stacklevel=2,
        )

    return RhoToolpathResult(
        operations=operations,
        polyline=_tessellate(operations, samples_per_radian),
        declined_regions=tuple(declined),
    )
