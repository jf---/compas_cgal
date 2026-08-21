"""Engagement-regulated trochoidal pocketing with a SECOND knob: the loop radius.

`compas_cgal.engagement_toolpath.engagement_controlled_toolpath` regulates the
ADVANCE only. Its loop radius is whatever the straight-skeleton guide derived from
the clearance at that station -- the largest gouge-free circle -- so where no
advance is admissible it emits that maximal circle anyway and reports it. Measured
on a 20x12 pocket with a 2 mm tool, that pins the worst loop engagement away from
the chain entries at 125.8 deg for a 40, a 60 and an 80 deg cap alike.

READ THAT NUMBER PRECISELY: it is the worst over ALL circles away from the chain
entries, and every circle carrying it is one the advance-only generator already
FORCED -- refused by its own predicate and emitted with a warning because no
shorter advance exists there. Over the circles it ACCEPTS, the same walk finds
43.4 / 62.5 / 81.7 deg at those three caps. So the radius knob is not what stops a
generator from claiming a cap it does not hold; that is the probe ring's job
(`LOOP_PROBE_COUNT`). What the radius knob buys is a smaller circle where a
maximal one is genuinely inadmissible, which converts forced circles into cutting
ones: at a 40 deg cap it takes the forced count from 88 to 44 and the accepted
count from 272 to 316.

This module adds the missing knob. At each station the loop radius is chosen by a
LADDER SEARCH over the same exact `cap_exceeded` predicate: the candidate radii
are a descending integer-indexed grid below the station's clearance-derived
maximum, and the search takes the LARGEST ADMISSIBLE one.

WHY A LADDER AND NEVER A BISECTION
----------------------------------
ENGAGEMENT IS NOT MONOTONE IN THE LOOP RADIUS. Measured at one mid-path station of
a 20x12 pocket (centre 14.0, 6.0, tool 2 mm, eight maximal loops swept behind it at
a half-tool-diameter advance), the peak engagement over the loop falls to 69.2 deg
at radius 4.148 and then RISES again as the circle shrinks further -- 72.4 deg at
4.098, 77.0 deg at 4.048, 85.0 deg at 3.998 -- so one guide step of extra retreat
carries the loop back over a 70 deg cap it had just met. The admissible set is
therefore not an up-set in the rung index, and a bisection -- which is only correct
on one -- would converge on a rung that is not the largest admissible one while
still looking like it worked: on that state the scan returns rung 17 (radius 4.148)
and a bisection returns rung 36 (radius 3.198), a circle nearly a millimetre
smaller. `_largest_admissible_radius` consequently SCANS the ladder from the top
and returns the first rung that passes, which is the largest admissible radius
under any pass/fail pattern whatsoever.

This is the second instance of the same structure in this repository: `spacing`
does not order engagement either (`benchmarks.figure6`, "Spacing does not order
engagement"), which is why the constant-spacing baseline is a minimum over a
brute-force sweep rather than a bisection. Both are recorded so the assumption is
not quietly reintroduced.

WHAT THIS GUARANTEES, EXACTLY
-----------------------------
Engagement <= `tea_cap_deg`, decided by an exact predicate, **at each evaluated
tool position**. That is the whole claim, and it is the same claim
`engagement_toolpath` makes -- no more. It is NOT a continuous guarantee between
evaluated positions, no certificate is produced, none is returned, and no
`MotionWitness` / `CapRefutation` object exists on this path. The bridge cuts
between machining circles are not regulated at all, by either generator.

Chain-entry loops remain a full slot by construction: the tool meeting virgin
stock is surrounded by material at every radius, so no rung is admissible there,
the maximal circle is emitted, and it is counted and warned about rather than
hidden.

WHY MORE PASSES, AND WHY THE PATH GETS LONGER
---------------------------------------------
A loop smaller than the station's maximum sweeps a narrower annulus, so it leaves
the outer band of that station's reachable material uncut. The chain is therefore
walked REPEATEDLY: each sweep opens the band further outward, and a station is
finished only once its maximal circle is emitted -- which is what the unregulated
walk does on its first and only pass. Tighter caps force smaller first loops, hence
more sweeps, hence LONGER paths. That is the physics of taking a smaller radial bite
and it is not tuned away here.

RELATION TO THE ADVANCE-ONLY GENERATOR
--------------------------------------
`engagement_controlled_toolpath` is untouched and remains the entry point for
existing callers. This is a SECOND public entry point rather than a flag on the
first, because the two do not differ by a parameter: this one emits several passes
per skeleton chain, plunges and retracts once per pass, and carries its own
diagnostic counts. A boolean that silently multiplied the number of passes in the
returned operation stream would be a worse API than a separate name.

Where every station's maximal circle is admissible on the first pass -- any cap
loose enough that the advance regulation alone already complies -- the ladder
returns rung 0 everywhere, the second sweep finds nothing to do, and the emitted
stream is the advance-only generator's stream.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from dataclasses import replace
from typing import List
from typing import Optional
from typing import Tuple

from compas.geometry import Polygon

from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import MAX_ADVANCE_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import POLYLINE_SAMPLES_PER_RADIAN
from compas_cgal.engagement_toolpath import UnavoidableEngagementWarning
from compas_cgal.engagement_toolpath import _guide_chains
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.engagement_toolpath import _largest_admissible_advance
from compas_cgal.engagement_toolpath import _line_operation
from compas_cgal.engagement_toolpath import _loop_operation
from compas_cgal.engagement_toolpath import _probe_positions
from compas_cgal.engagement_toolpath import _Regulation
from compas_cgal.engagement_toolpath import _station_is_admissible
from compas_cgal.engagement_toolpath import _tessellate
from compas_cgal.engagement_toolpath import _unit_tangent
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

# Spacing of the radius ladder, in tool diameters. Deliberately the SAME grid the
# advance search runs on (`GUIDE_STEP_TOOL_DIAMETERS`), because the two knobs are
# the same physical quantity measured in perpendicular directions: an advance of
# one station is a radial depth of cut of one station TANGENTIALLY along the
# guide, and a rung of the ladder is a radial depth of cut of one station
# OUTWARD from it -- the next sweep's loop reaches exactly one rung further than
# the last one did. Sharing the grid means the two searches quantise the cap with
# the same resolution: by the derivation recorded on GUIDE_STEP_TOOL_DIAMETERS,
# one step is worth at most ~7 deg of engagement in the textbook radial-immersion
# relation, spent in the safe direction because a rung is only ever taken when its
# evaluated positions passed.
RADIUS_LADDER_STEP_TOOL_DIAMETERS = GUIDE_STEP_TOOL_DIAMETERS

# How far below the station's maximal radius the ladder reaches, in tool
# diameters. Same bound and same derivation as `MAX_ADVANCE_TOOL_DIAMETERS`:
# TEA = 2*acos(1 - ae/r) saturates at a full turn when the radial depth of cut ae
# reaches 2r = D, so retreating the radius by more than one tool diameter cannot
# convert a slotting loop into an admissible one -- beyond that depth the loop is
# either already in swept void, in which case the scan stopped at a higher rung,
# or it is surrounded by material, in which case no rung helps and the station is
# forced. Extending the ladder past D would buy nothing but predicate calls.
RADIUS_LADDER_SPAN_TOOL_DIAMETERS = MAX_ADVANCE_TOOL_DIAMETERS

# Rungs on the ladder, span over spacing. Rung 0 is the station's own maximal
# radius -- the ONE the advance-only generator always emits -- so a ladder whose
# rung 0 is admissible reproduces that generator exactly.
RADIUS_LADDER_RUNGS = int(RADIUS_LADDER_SPAN_TOOL_DIAMETERS / RADIUS_LADDER_STEP_TOOL_DIAMETERS)

# Rung index of the station's maximal, clearance-derived radius. Named because it
# is load-bearing in three places: it is the first rung the scan tries, it is the
# rung whose emission FINISHES a station, and it is the rung emitted when the
# scan finds nothing admissible.
FULL_RADIUS_RUNG = 0

# Passes one skeleton chain may take before the walk is declared broken. This is a
# BUDGET, not the termination rule, and it is deliberately not dressed up as a
# theorem. The walk DOES terminate structurally: `_radius_ladder` offers a station
# only radii strictly above the largest it has already emitted, so a station can
# emit at most `RADIUS_LADDER_RUNGS` times, and a pass that emits nothing ends the
# chain -- but that only bounds the passes by stations x rungs, which is far too
# large to catch anything. The value below is an empirical guard: the worst
# requirement measured over 6x4, 10x6, 12x8 and 20x12 pockets with a 2 mm tool at
# caps from 20 to 170 degrees is TWO passes on any chain, so a chain reaching forty
# has a station re-climbing rungs it already emitted. Exhausting the budget RAISES
# rather than silently truncating the chain, because a truncated chain is
# unmachined material presented as a finished toolpath.
MAX_RADIAL_SWEEPS_PER_CHAIN = RADIUS_LADDER_RUNGS


class RadialSweepBudgetExceededError(RuntimeError):
    """A skeleton chain did not finish within `MAX_RADIAL_SWEEPS_PER_CHAIN` sweeps."""


@dataclass(frozen=True)
class _SweepOutcome:
    """What one pass over a skeleton chain emitted and how hard it had to work.

    Attributes:
        operations: The pass's operation stream -- plunge, machining circles and
            bridges, retract -- or empty when the pass found nothing left to cut.
        entry: Tool-centre point the pass plunged at, or ``None`` when empty.
        exit: Tool-centre point the pass retracted from, or ``None`` when empty.
        forced_radii: Machining circles emitted at the station's maximal radius
            because NO rung of the ladder was admissible.
        forced_advances: Advances taken past a refusing predicate, as counted by
            `_largest_admissible_advance`.
        first_loop_forced: Whether the pass's FIRST machining circle was one of
            the forced ones. On a chain's first pass that circle is the entry into
            virgin stock, which is a full slot for any generator at any radius.
    """

    operations: List[ToolpathOperation]
    entry: Optional[Tuple[float, float]]
    exit: Optional[Tuple[float, float]]
    forced_radii: int
    forced_advances: int
    first_loop_forced: bool


def _radius_ladder(full_radius: float, emitted_radius: float, ladder_step: float) -> List[float]:
    """Candidate loop radii at one station, largest first.

    Rung ``k`` is ``full_radius - k * ladder_step``. The ladder stops at the
    largest radius already emitted at this station, because a loop at or below it
    sweeps an annulus this station has already swept and would remove nothing:
    that floor is what makes every emission strictly increase the station's
    cleared radius by at least one rung, which is what makes the sweep loop
    terminate structurally rather than by a budget.

    Rung 0 is always offered when anything is left, even when it lies below one
    ladder step, so a station whose clearance admits only a hair of a circle still
    gets the circle the advance-only generator would have emitted there.

    Args:
        full_radius: The station's clearance-derived maximal radius.
        emitted_radius: Largest radius already emitted at this station; ``0.0``
            before the first pass.
        ladder_step: Spacing between rungs in model units.

    Returns:
        The candidate radii in descending order, empty when nothing is left.
    """
    if not full_radius > emitted_radius:
        return []
    ladder = [full_radius]
    for rung in range(1, RADIUS_LADDER_RUNGS):
        radius = full_radius - rung * ladder_step
        if radius <= emitted_radius:
            break
        ladder.append(radius)
    return ladder


def _loop_reaches_material(stock: Stock, station: _GuideStation, advance: Tuple[float, float], tool_radius: float) -> bool:
    """Whether this machining circle still bites into uncut stock.

    WHY THE QUESTION HAS TO BE ASKED. In the steady trochoidal regime the material
    a loop meets is a crescent at its OUTER rim, so a loop drawn inward is not a
    lighter cut, it is NO cut: it spins inside the annulus its predecessors already
    swept. A cap test alone reads those idle loops as admissible -- they engage
    nothing, so they exceed nothing -- and a scan that takes the largest of them
    never moves the frontier, so the next station admits a still smaller loop, and
    the ladder walks itself to the bottom one rung per station before falling back
    on the maximal circle in virgin stock. Measured on a 10x6 pocket at a 40 deg
    cap, that produced 360 deg loops the advance-only generator never emits. This
    condition is what stops it.

    THE POINT ASKED ABOUT is the loop's deepest reach at each evaluated position:
    the tool centre pushed one tool radius further out along the ray from the
    station centre, i.e. the point at distance ``radius + tool_radius`` from the
    centre. That is where a trochoidal loop takes its bite, so material there is
    exactly what "this loop still has something to cut" means. Membership is
    `Stock.contains`, exact point location in the exact arrangement -- no tolerance
    and no reported double. Like the cap probes themselves it is SAMPLED at the
    same handful of positions, so it decides where the loop is asked about, never
    how the answer is computed.

    Args:
        stock: The current stock (unmodified by this call).
        station: The station under test, carrying the candidate radius.
        advance: Unit direction of travel into this station.
        tool_radius: Tool radius.

    Returns:
        ``True`` if at least one evaluated position reaches uncut stock.
    """
    for px, py in _probe_positions(station, advance):
        ux, uy = _unit_tangent(station.cx, station.cy, px, py)
        if stock.contains(px + tool_radius * ux, py + tool_radius * uy):
            return True
    return False


def _largest_admissible_radius(
    stock: Stock,
    station: _GuideStation,
    emitted_radius: float,
    advance: Tuple[float, float],
    regulation: _Regulation,
) -> Tuple[Optional[int], bool]:
    """SCAN the radius ladder for the largest radius that both complies and cuts.

    THE SCAN IS THE ALGORITHM, and it must not become a bisection. Engagement is
    not monotone in the loop radius (module docstring, with the measurement), so
    the admissible rungs are not an up-set and the pass/fail pattern down the
    ladder can alternate. Walking from rung 0 downward and returning the FIRST
    pass is correct under any pattern: rung indices order the radii strictly
    downward, so the first rung that passes carries the largest admissible radius,
    whatever the rungs below it do.

    A rung is admissible on TWO exact conditions, and both are load-bearing:

    1. no evaluated position exceeds the cap, and
    2. some evaluated position still reaches uncut stock --
       `_loop_reaches_material` carries why. Condition 1 alone is satisfied by
       every loop small enough to spin inside already-swept void, which is how a
       scan that omits condition 2 walks itself to the bottom of the ladder.

    Rung 0, the station's maximal circle, is exempt from condition 2: it is the
    circle the advance-only generator emits unconditionally and the one whose
    emission FINISHES a station, so when it complies it is taken whether or not
    there is anything left there to cut.

    EVERY rung is decided at `LOOP_PROBE_ANGLES_DEG`, the uniform advance-phased
    ring the advance-only generator uses. That the SAME probe set serves both the
    maximal circle and a reduced one is now a property of the set rather than an
    assumption about it: a reduced circle sits inside its station's clearance disk
    where material can lie on any side, and a uniform ring makes no assumption
    about which side that is. The triple this ring replaced did -- it looked only
    at the advance-facing half -- and the constant's comment in
    `compas_cgal.engagement_toolpath` records the measurement that falsified it.

    Density is bounded by cost, not by belief: `LOOP_PROBE_COUNT` positions are
    evaluated per rung, and a rung that fails is abandoned at the first refusing
    probe, so only a rung that PASSES pays for the whole ring.

    The scan is also the cheap order for the common case. A cap loose enough that
    the maximal circle already complies costs one rung -- the same evaluation the
    advance-only generator makes -- and only a station whose maximal circle is
    refused pays for the descent.

    Args:
        stock: The current stock (unmodified by this call).
        station: The station under test, carrying its maximal radius.
        emitted_radius: Largest radius already emitted at this station.
        advance: Unit direction of travel into this station, for probe placement.
        regulation: The validated parameters, for the tool radius and the exact
            rational cap surrogate.

    Returns:
        ``(rung, forced)``. ``rung`` is the index into the ladder, or ``None``
        when the station has nothing left to cut. ``forced`` is ``True`` when no
        rung was admissible and `FULL_RADIUS_RUNG` is returned anyway, which is
        the virgin-stock and neck regime: refusing to cut is not an option there,
        so the maximal circle is emitted and counted.
    """
    ladder = _radius_ladder(station.radius, emitted_radius, regulation.guide_step)
    if not ladder:
        return None, False
    for rung, radius in enumerate(ladder):
        candidate = replace(station, radius=radius)
        if not _station_is_admissible(stock, candidate, advance, regulation.tool_radius, regulation.cap_ratio):
            continue
        if rung == FULL_RADIUS_RUNG or _loop_reaches_material(stock, candidate, advance, regulation.tool_radius):
            return rung, False
    return FULL_RADIUS_RUNG, True


def _radial_sweep(
    stock: Stock,
    stations: List[_GuideStation],
    emitted_radii: List[float],
    finished: List[bool],
    path_index: int,
    regulation: _Regulation,
    cut_z: float,
) -> _SweepOutcome:
    """Walk one skeleton chain once, cutting only what is still left to cut.

    A station is FINISHED once its maximal circle has been emitted: that circle
    sweeps the whole annulus the clearance allows there, which is precisely what
    the advance-only generator emits on its single pass, so nothing further is
    reachable at that station and later sweeps skip it. Stations the advance
    search jumped OVER are finished too, on the advance bound's own coverage
    argument (`MAX_ADVANCE_TOOL_DIAMETERS`): two maximal circles no more than a
    tool diameter apart sweep overlapping annuli, so the ground between them is
    covered and re-walking it on the next sweep would emit circles that remove
    nothing.

    That coverage argument only holds behind a MAXIMAL circle. After a reduced
    one the walk therefore advances a single station and marks nothing finished --
    the same minimum advance the advance-only generator takes when its predicate
    refuses, and the conservative choice in exactly the heavily loaded regime where
    the radius had to be cut back.

    Args:
        stock: The depleting stock, mutated in place.
        stations: The chain's ordered stations.
        emitted_radii: Per station, the largest radius emitted there so far;
            updated in place.
        finished: Per station, whether it has nothing left to cut; updated in
            place.
        path_index: Path index stamped on this pass's operations.
        regulation: The validated parameters.
        cut_z: Cutting-plane height.

    Returns:
        The pass's outcome, whose `operations` are empty when it found nothing.
    """
    last = len(stations) - 1
    operations: List[ToolpathOperation] = []
    entry: Optional[Tuple[float, float]] = None
    previous_entry: Optional[Tuple[float, float]] = None
    previous_index: Optional[int] = None
    forced_radii = 0
    forced_advances = 0
    first_loop_forced = False

    index = 0
    while True:
        station = stations[index]
        rung: Optional[int] = None
        forced = False
        if not finished[index]:
            if previous_index is None:
                advance = (station.tx, station.ty)
            else:
                origin = stations[previous_index]
                advance = _unit_tangent(origin.cx, origin.cy, station.cx, station.cy)
            rung, forced = _largest_admissible_radius(stock, station, emitted_radii[index], advance, regulation)

        if rung is not None:
            radius = station.radius - rung * regulation.guide_step
            loop_station = replace(station, radius=radius)
            loop_entry = loop_station.entry
            if previous_entry is None:
                entry = loop_entry
                operations.append(_line_operation(loop_entry, loop_entry, regulation.clearance_z, cut_z, OperationType.PLUNGE, path_index))
                first_loop_forced = forced
            else:
                operations.append(_line_operation(previous_entry, loop_entry, cut_z, cut_z, OperationType.CUT, path_index))
                stock.subtract_capsule_quad(previous_entry[0], previous_entry[1], loop_entry[0], loop_entry[1], regulation.tool_radius)
            operations.append(_loop_operation(loop_station, cut_z, path_index))
            stock.subtract_arc_sweep_local(
                loop_station.cx,
                loop_station.cy,
                loop_entry[0],
                loop_entry[1],
                loop_entry[0],
                loop_entry[1],
                loop_station.clockwise,
                regulation.tool_radius,
            )
            emitted_radii[index] = radius
            finished[index] = rung == FULL_RADIUS_RUNG
            forced_radii += int(forced)
            previous_entry = loop_entry
            previous_index = index

        if index == last:
            break
        if rung == FULL_RADIUS_RUNG:
            next_index, advance_forced = _largest_admissible_advance(
                stock,
                stations,
                index,
                min(index + regulation.window, last),
                regulation.tool_radius,
                regulation.cap_ratio,
            )
            forced_advances += int(advance_forced)
            if not advance_forced:
                for skipped in range(index + 1, next_index):
                    finished[skipped] = True
            index = next_index
        else:
            index += 1

    if previous_entry is not None:
        operations.append(_line_operation(previous_entry, previous_entry, cut_z, regulation.clearance_z, OperationType.RETRACT, path_index))

    return _SweepOutcome(
        operations=operations,
        entry=entry,
        exit=previous_entry,
        forced_radii=forced_radii,
        forced_advances=forced_advances,
        first_loop_forced=first_loop_forced,
    )


@dataclass(frozen=True)
class _ChainOutcome:
    """Aggregate of every pass over one skeleton chain.

    Attributes:
        sweeps: Passes that emitted something.
        entry_forced: Whether the chain's first machining circle -- the tool
            meeting virgin stock -- had to be emitted over the cap.
        forced_radii: Machining circles other than that entry emitted over the cap
            because no rung of the ladder was admissible.
        forced_advances: Advances taken past a refusing predicate.
        exit: Tool-centre point the chain's last pass retracted from, or ``None``.
    """

    sweeps: int
    entry_forced: bool
    forced_radii: int
    forced_advances: int
    exit: Optional[Tuple[float, float]]


def _machine_chain_radially(
    stock: Stock,
    stations: List[_GuideStation],
    path_index: int,
    regulation: _Regulation,
    cut_z: float,
    last_exit: Optional[Tuple[float, float]],
    operations: List[ToolpathOperation],
) -> _ChainOutcome:
    """Sweep one skeleton chain until every station has had its maximal circle.

    Each pass is linked to the previous one by a clearance-height `LINK`, exactly
    as consecutive chains are, so the returned stream is a continuous motion
    sequence with no unexplained jumps for a post-processor to guess at.

    Termination is structural, not budgeted: `_radius_ladder` offers a station only
    radii strictly above the largest it has already emitted, so each of a station's
    emissions climbs at least one rung of a ladder with `RADIUS_LADDER_RUNGS`
    rungs, and a pass that emits nothing ends the chain. That argument bounds the
    passes by stations times rungs, which is finite but useless as a guard, so
    `MAX_RADIAL_SWEEPS_PER_CHAIN` is an empirical budget on top of it -- twenty
    times the worst requirement measured on any chain -- and it raises.

    Args:
        stock: The depleting stock, mutated in place.
        stations: The chain's ordered stations.
        path_index: Path index stamped on this chain's operations.
        regulation: The validated parameters.
        cut_z: Cutting-plane height.
        last_exit: Tool-centre point the previous emitted pass retracted from, or
            ``None`` when nothing has been emitted yet.
        operations: Output stream, appended to in place.

    Returns:
        The chain's aggregate outcome.

    Raises:
        RadialSweepBudgetExceededError: If the chain did not finish within
            `MAX_RADIAL_SWEEPS_PER_CHAIN` passes.
    """
    emitted_radii = [0.0] * len(stations)
    finished = [False] * len(stations)
    sweeps = 0
    entry_forced = False
    forced_radii = 0
    forced_advances = 0

    for _ in range(MAX_RADIAL_SWEEPS_PER_CHAIN):
        outcome = _radial_sweep(stock, stations, emitted_radii, finished, path_index, regulation, cut_z)
        forced_advances += outcome.forced_advances
        if not outcome.operations:
            break
        if sweeps == 0:
            entry_forced = outcome.first_loop_forced
            forced_radii += outcome.forced_radii - int(outcome.first_loop_forced)
        else:
            forced_radii += outcome.forced_radii
        if last_exit is not None:
            operations.append(_line_operation(last_exit, outcome.entry, regulation.clearance_z, regulation.clearance_z, OperationType.LINK, path_index))
        operations.extend(outcome.operations)
        last_exit = outcome.exit
        sweeps += 1
    else:
        raise RadialSweepBudgetExceededError(
            f"Skeleton chain {path_index} still had material to cut after {MAX_RADIAL_SWEEPS_PER_CHAIN} sweeps. "
            "Each sweep must lift some station at least one rung of a "
            f"{RADIUS_LADDER_RUNGS}-rung ladder, so this means a station is re-climbing rungs it already emitted."
        )

    return _ChainOutcome(sweeps=sweeps, entry_forced=entry_forced, forced_radii=forced_radii, forced_advances=forced_advances, exit=last_exit)


def radius_regulated_toolpath(
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
) -> ToolpathResult:
    """Trochoidal pocketing whose ADVANCE and LOOP RADIUS are both regulated.

    Walks the straight-skeleton guide chain by chain, and each chain repeatedly.
    At every station the loop radius is the largest one on a descending
    integer-indexed ladder below the station's clearance-derived maximum whose
    evaluated tool positions all report the engagement cap NOT exceeded against
    the depleting exact stock; the advance to the next station is chosen as in
    `compas_cgal.engagement_toolpath.engagement_controlled_toolpath`, which this
    function leaves untouched.

    The ladder is SCANNED, never bisected. Engagement is not monotone in the loop
    radius -- measured at one mid-path station of a 20x12 pocket, the peak falls to
    69.2 deg at radius 4.148 and rises again to 72.4 deg at 4.098 and 85.0 deg at
    3.998 -- so the admissible rungs are not an up-set and a bisection would
    silently return a rung that is not the largest admissible one.

    WHAT THIS GUARANTEES, EXACTLY: engagement <= *tea_cap_deg*, decided by an exact
    predicate, at each EVALUATED tool position -- the loop entry point and the
    `LOOP_PROBE_ANGLES_DEG` ring on each accepted machining circle. It is NOT a
    continuous guarantee between evaluated positions -- raising `LOOP_PROBE_COUNT`
    narrows that gap and does not close it -- and it says nothing about the bridge
    cuts, which neither generator regulates. No certificate is produced or
    returned.

    A smaller loop leaves the outer band of its station uncut, so the chain is
    swept again until every station has emitted its maximal circle. Tighter caps
    therefore produce MORE passes and LONGER paths; that is the cost of a smaller
    radial bite, and it is not tuned away.

    Args:
        polygon: Outer pocket boundary `Polygon` in the world XY plane.
        tool_diameter: Tool diameter; the tool radius is half of this.
        tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``. Converted
            once to the exact rational surrogate ``4*sin^2(theta/2)`` the exact
            predicate consumes.
        holes: Optional island polygons strictly inside *polygon*.
        guide_step_tool_diameters: Guide station spacing in tool diameters; also
            the resolution of every advance the search can pick AND the spacing of
            the radius ladder.
        max_advance_tool_diameters: Largest advance the search considers, in tool
            diameters.
        radial_clearance: Safety clearance subtracted from each loop's available
            radius. Defaults to ``RADIAL_CLEARANCE_FRACTION * tool_diameter``.
        climb: ``True`` for climb milling (clockwise loops), ``False`` for
            conventional milling.
        cut_z: Z-height of the cutting plane.
        clearance_z: Z-height for rapid travel between passes. Defaults to
            ``cut_z + CLEARANCE_RISE_TOOL_DIAMETERS * tool_diameter``.
        max_passes: Maximum number of skeleton chains the guide may emit.
        samples_per_radian: Tessellation density of the returned polyline.

    Returns:
        ToolpathResult: The typed operation stream -- `PLUNGE`, `CUT` machining
        circles and bridge lines, `RETRACT`, and clearance-height `LINK` moves
        between passes -- plus the tessellated visualisation polyline. The same
        types the other generators emit, so `audit_toolpath_engagement` and every
        downstream consumer work unchanged. Note that a chain contributes one
        plunge/retract pair PER PASS, not one in total.

    Raises:
        InvalidEngagementCapDegreesError: If *tea_cap_deg* is not in ``(0, 180]``.
        NonPositiveToolDiameterError: If *tool_diameter* is not strictly positive.
        InvalidGuideResolutionError: If the guide step or advance bound leaves no
            integer bracket to bisect.
        InvalidClearanceHeightError: If *clearance_z* is not above *cut_z*.
        EmptyGuideError: If the pocket admits no gouge-free trochoid at this tool.
        InvalidPolygonError: If *polygon* or a hole is non-planar or degenerate.
        RadialSweepBudgetExceededError: If a chain did not finish within
            `MAX_RADIAL_SWEEPS_PER_CHAIN` passes.

    Warns:
        UnavoidableEngagementWarning: When machining circles had to be emitted at
            positions the exact cap predicate refuses at every rung of the ladder
            (never silently).
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
    last_exit: Optional[Tuple[float, float]] = None
    entry_slots = 0
    forced_radii = 0
    forced_advances = 0
    sweeps = 0
    for path_index, stations in enumerate(chains):
        outcome = _machine_chain_radially(stock, stations, path_index, regulation, cut_z, last_exit, operations)
        entry_slots += int(outcome.entry_forced)
        forced_radii += outcome.forced_radii
        forced_advances += outcome.forced_advances
        sweeps += outcome.sweeps
        last_exit = outcome.exit

    if entry_slots or forced_radii or forced_advances:
        warnings.warn(
            f"{entry_slots + forced_radii} machining circle(s) were emitted at positions where the exact cap predicate "
            f"reports tea_cap_deg={tea_cap_deg} exceeded at EVERY rung of the {RADIUS_LADDER_RUNGS}-rung radius ladder: "
            f"{entry_slots} chain-entry loop(s), where the tool first meets virgin stock and a first cut is a full slot "
            f"by construction, and {forced_radii} other station(s) whose surrounding material no loop radius can escape. "
            f"{forced_advances} advance(s) were also taken past a refusing predicate. Cutting less is not an option "
            f"there, so they are emitted and counted rather than dropped. The path took {sweeps} pass(es) over "
            f"{len(chains)} skeleton chain(s).",
            UnavoidableEngagementWarning,
            stacklevel=2,
        )

    return ToolpathResult(operations=operations, polyline=_tessellate(operations, samples_per_radian))
