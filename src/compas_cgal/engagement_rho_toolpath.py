"""Rho-regulated trochoidal pocketing: the loop radius and the advance are chosen together.

`compas_cgal.engagement_toolpath.engagement_controlled_toolpath` regulates the
ADVANCE and takes whatever loop radius the straight-skeleton guide derived from
the clearance at that station.
`compas_cgal.engagement_radial_toolpath.radius_regulated_toolpath` adds a radius
LADDER below that maximum. Neither of them bounds how fast the emitted radius may
CHANGE from one machining circle to the next, and neither notices that a station
whose gouge-free radius has fallen to or below the tool radius cannot carry a
trochoid at all. This module regulates those two things, on the same guide, with
the same exact predicate.

WHAT THIS GUARANTEES, EXACTLY
-----------------------------
Engagement <= `tea_cap_deg`, decided by an exact predicate, **at each evaluated
tool position** -- the same claim the other two generators make, and no more. It
is NOT a continuous guarantee between evaluated positions, no certificate is
produced, none is returned, and no `MotionWitness` / `CapRefutation` object
exists on this path. What is ADDED here is not a stronger engagement claim; it is
two structural properties of the emitted circle sequence, each of which is a
consequence of an emission rule rather than of a measurement:

1. every emitted machining circle has ``rho > r``, so it sweeps an annulus with
   an uncut core rather than a filled disk, and
2. consecutive machining circles differ in radius by at most
   ``MAX_LOOP_RADIUS_STEP_TOOL_RADII * r``, so their swept annuli overlap.

THE THREE RULES
---------------
DEGENERACY FLOOR. A tool of radius ``r`` running a circle of radius ``rho``
sweeps the radii ``[rho - r, rho + r]`` about the centre, so the sweep has an
uncut core exactly when ``rho > r``. At or below that the "trochoid" is a bore
wearing a circle's name (`benchmarks.quality.DEGENERATE_LOOP_RATIO`, and
`docs/loop_radius_degeneracy.md` for where the boundary comes from). Stations
whose gouge-free radius sits at or below it are NOT loop sites here, and the
material at them is left to the neighbouring loops' sweeps. That is a REAL COST,
paid deliberately and reported rather than hidden: on `rect_12x8` it takes the
uncut fraction from 0.0080 to 0.0295, and `UnavoidableEngagementWarning` carries
how many stations were treated that way. See "WHAT THIS DOES NOT FIX" for why no
motion in this vocabulary covers that material without misreporting itself.

RADIUS-STEP BOUND. Two machining circles at radii ``rho`` and ``rho'`` sweep
annuli that overlap only while ``|rho - rho'| < 2r``, so
`MAX_LOOP_RADIUS_STEP_TOOL_RADII` is the overlap limit itself and not a tuned
number. The advance search will not step over a radius change larger than it.

TANGENT-CONTINUOUS LINKING. This is the rule the other two generators are missing
and the reason the radius-step bound is safe to add. With the entry placed one
loop radius along the guide NORMAL, the chord between two entries on a straight
guide is ``advance * t + (rho' - rho) * n``, so it leaves the guide tangent by
``atan(|d rho| / advance) = atan(d rho / d s)`` -- the arctangent of the CLEARANCE
SLOPE, which does not depend on the advance. Shortening the advance shrinks
numerator and denominator together. Measured on `rect_12x8` with a 2 mm tool: all
24 of the advance-only generator's tangent breaks are the four corner chains'
six junctions each, every one of them at the 35.26 degrees that is
``atan(sin 45 deg)`` on a corner bisector, while the constant-clearance spine has
none. **Bounding the radius step without fixing this makes matters worse**: it
forces shorter advances, which emits more bridges, each still tilted by the same
angle. So the entry here is placed at the EXTERNAL COMMON TANGENT of the two
circles instead -- `_entry_normal` -- which makes the chord perpendicular to the
entry radius, hence tangent to both circles, for any radius difference.

WHAT THIS DOES NOT FIX, AND WHY
-------------------------------
A chain entering virgin stock is a full-immersion cut at EVERY non-degenerate
radius: a loop escapes 360 degrees only while its far point stays within ``2r`` of
the plunge hole, which is ``2 * rho < 2r``, which is degeneracy. The crossover is
exact and measured at ``rho = r`` (`docs/loop_radius_degeneracy.md`). No radius
this module may emit changes it, so the engagement STEP from an entry loop to the
motion after it is at least ``360 - cap``. That is a property of entering solid
material with a circular motion, not a defect this or any radius rule can remove.

The material at a convex corner is the same kind of fact. Its reachable centres
run out at a clearance of exactly ``r``, so every motion that touches the corner
material is full immersion, and the vocabulary offers only two ways to spend it:
a circle, which is then a bore, or a straight move, which is then a slot. THREE
VARIANTS WERE MEASURED on `rect_12x8` (2 mm tool, 120 degree cap), and each pays
in a different currency:

| variant                                | uncut  | degenerate | slotting | tangent breaks |
| -------------------------------------- | ------ | ---------- | -------- | -------------- |
| bore the corner (`engagement_toolpath`) | 0.0080 | 4          | 0        | 24             |
| leave it (THIS MODULE)                  | 0.0295 | 0          | 0        | 0              |
| run out to it on a straight move        | 0.0033 | 0          | 12       | 4              |

This module takes the middle row: it never emits a motion that misreports what it
is, and it says how much it left. Clearing a corner narrower than a trochoid is a
second operation with a smaller tool, not a rule this generator can be given.

WHY A LADDER AND NEVER A BISECTION
----------------------------------
Engagement is not monotone in the advance, any more than it is in the loop radius
(`compas_cgal.engagement_radial_toolpath`, with the measurement). The admissible
advances are therefore not an up-set in the station index, and a bisection -- which
is only correct on one -- can converge on a station that is not the furthest
admissible while still looking like it worked. `_largest_admissible_advance`
consequently SCANS the window from the furthest station downward and returns the
first that passes, which is the furthest admissible advance under any pass/fail
pattern whatsoever. The argument is about the ORDER of the scan and nothing else.
`tests/test_engagement_rho_toolpath.py` pins the DIRECTION with a pass/fail
pattern a bisection gets wrong and an inverted scan gets backwards.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Tuple

import numpy as np
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Polygon

from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import MAX_ADVANCE_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import POLYLINE_SAMPLES_PER_RADIAN
from compas_cgal.engagement_toolpath import DegenerateMachiningCircleError
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
from compas_cgal.toolpath import ToolpathResult

# A machining circle must sweep an ANNULUS, not a filled disk. The tool covers
# the radii [rho - r, rho + r] about the loop centre, so the sweep has an uncut
# core exactly when rho > r; at or below it the cutter passes over its own centre
# and never unloads. This is the qualitative boundary at which the annulus loses
# its hole -- physics, not a tuned constant -- and it is the same boundary
# `benchmarks.quality.DEGENERATE_LOOP_RATIO` measures against.
DEGENERATE_LOOP_TOOL_RADII = 1.0

# Largest radius change the advance search will step over, in tool radii. Two
# machining circles at radii rho and rho' sweep the annuli [rho - r, rho + r] and
# [rho' - r, rho' + r] about their centres; those overlap radially only while
# |rho - rho'| < 2r. The bound IS the overlap limit, so a step at or under it
# leaves consecutive circles sharing swept material and a step over it opens a
# band neither of them touches. Nothing is tuned: raising it past 2 would admit
# provably disjoint annuli.
#
# IT DOES NOT BIND AT THE SHIPPED DEFAULTS, and saying so is the difference
# between a guard and a claim. A clearance function is 1-Lipschitz, so over the
# widest advance the search may consider -- MAX_ADVANCE_TOOL_DIAMETERS = 1.0 tool
# diameters, which is 2r -- the radius cannot change by more than 2r, which is
# this bound. The measured improvement in `max_loop_radius_step` on the corpus
# therefore comes from the DEGENERACY FLOOR removing the near-zero corner-tip
# circles, not from this rule ever refusing a candidate. The rule earns its place
# on a coarser guide or a wider advance window, and
# `tests/test_engagement_rho_toolpath.py` exercises it there.
MAX_LOOP_RADIUS_STEP_TOOL_RADII = 2.0

# Largest clearance slope the external-tangent construction accepts, as
# d(rho)/d(s) along the guide. The external common tangent of two circles of
# radii rho and rho' whose centres are D apart exists only while
# |rho - rho'| <= D, i.e. while the slope is at most one. A clearance function is
# 1-Lipschitz -- moving a unit along the guide cannot change the distance to the
# nearest wall by more than a unit -- so a guide that reports more than this is
# not reporting a clearance, and the generator raises rather than clamping a
# number it cannot explain.
MAX_CLEARANCE_SLOPE = 1.0


class NonTrochoidalChainError(RuntimeError):
    """No station on a skeleton chain clears the degeneracy floor, so the chain admits no trochoid."""


class ImplausibleClearanceSlopeError(ValueError):
    """A guide chain's radius grows faster along its own arc length than a clearance function can."""


class EmptyRhoGuideError(RuntimeError):
    """The straight-skeleton guide yielded no chain with a non-degenerate machining circle."""


@dataclass(frozen=True)
class _RhoStation:
    """One station on a guide chain, carrying the entry the external tangent puts on it.

    Distinct from `compas_cgal.engagement_toolpath._GuideStation`, which places the
    entry one loop radius along the guide NORMAL. That placement is correct only
    where the radius is constant; where it changes, the chord between two entries
    tilts off the guide tangent by the arctangent of the clearance slope. This one
    places the entry on the external common tangent instead, so the chord is
    perpendicular to the entry radius and therefore tangent to the circle.

    Attributes:
        cx: X coordinate of the machining-circle centre, on the guide.
        cy: Y coordinate of the machining-circle centre, on the guide.
        radius: Gouge-free machining-circle radius derived by the guide.
        clockwise: Turn direction of the machining circle.
        tx: X component of the unit guide tangent.
        ty: Y component of the unit guide tangent.
        wx: X component of the unit ENTRY RADIUS direction, centre to entry.
        wy: Y component of the unit ENTRY RADIUS direction, centre to entry.
    """

    cx: float
    cy: float
    radius: float
    clockwise: bool
    tx: float
    ty: float
    wx: float
    wy: float

    @property
    def entry(self) -> Tuple[float, float]:
        """Tool-centre point where the bridge meets the machining circle."""
        return self.cx + self.radius * self.wx, self.cy + self.radius * self.wy

    @property
    def entry_tangent(self) -> Tuple[float, float]:
        """Unit travel direction on the circle at `entry`, in the loop's turn direction.

        A quarter turn from the entry radius, the way the circle is traversed:
        clockwise motion turns the outward radius to ``(wy, -wx)``, counterclockwise
        to ``(-wy, wx)``. On a constant-radius guide this reduces to the guide
        tangent, which is what the advance-only generator declares.
        """
        if self.clockwise:
            return self.wy, -self.wx
        return -self.wy, self.wx


def _entry_normal(station: _GuideStation, slope: float) -> Tuple[float, float]:
    """Unit entry-radius direction that puts the bridge on the external common tangent.

    THE CONSTRUCTION. A line at distance ``rho`` from one centre and ``rho'`` from
    the next, with both centres on the same side, has a unit normal ``n``
    satisfying ``n . (c' - c) = rho' - rho``; the tangency points are
    ``c - rho * n`` and ``c' - rho' * n``, and the chord between them is
    perpendicular to ``n`` by construction. Writing ``w = -n`` for the outward
    entry radius and decomposing on the guide frame gives
    ``w = -slope * t + sqrt(1 - slope^2) * m``, where ``m`` is the guide normal on
    the loop's material side and ``slope`` is ``d(rho)/d(s)``.

    At ``slope == 0`` this returns exactly the guide normal, so a constant-radius
    chain gets the advance-only generator's entry placement unchanged and the two
    agree wherever the older model was already right.

    Args:
        station: The station, for its guide tangent and turn direction.
        slope: Rate of change of the gouge-free radius along the guide's arc
            length at this station.

    Returns:
        The unit entry-radius direction, from the centre towards the entry point.

    Raises:
        ImplausibleClearanceSlopeError: If ``|slope|`` exceeds
            `MAX_CLEARANCE_SLOPE`, where the external tangent does not exist and
            the guide is not reporting a 1-Lipschitz clearance.
    """
    if not abs(slope) <= MAX_CLEARANCE_SLOPE:
        raise ImplausibleClearanceSlopeError(
            f"Guide station at ({station.cx!r}, {station.cy!r}) reports a radius slope of {slope!r} along its own arc "
            f"length, above the {MAX_CLEARANCE_SLOPE} a 1-Lipschitz clearance function can have. The external common "
            "tangent between consecutive machining circles does not exist at that slope."
        )
    material_x, material_y = (-station.ty, station.tx) if station.clockwise else (station.ty, -station.tx)
    lateral = math.sqrt(1.0 - slope * slope)
    return (
        -slope * station.tx + lateral * material_x,
        -slope * station.ty + lateral * material_y,
    )


def _rho_stations(stations: List[_GuideStation]) -> List[_RhoStation]:
    """Re-express one guide chain with external-tangent entries.

    The slope at a station is the central difference of the gouge-free radius
    over the guide's own arc length, one-sided at the chain ends -- the same local
    estimate the guide already uses for its tangent, and exact wherever the slope
    is constant, which is every straight chain of a polygonal pocket.

    Args:
        stations: The chain's ordered guide stations.

    Returns:
        The same stations carrying their entry-radius directions.

    Raises:
        ImplausibleClearanceSlopeError: If a station's radius slope exceeds
            `MAX_CLEARANCE_SLOPE`.
    """
    count = len(stations)
    rho_stations: List[_RhoStation] = []
    for index, station in enumerate(stations):
        before = stations[max(0, index - 1)]
        after = stations[min(count - 1, index + 1)]
        travel = math.hypot(after.cx - before.cx, after.cy - before.cy)
        # A single-station chain, or two coincident centres, has no direction in
        # which the radius could be said to change: the entry falls back to the
        # guide normal, which is what `slope == 0` returns.
        slope = 0.0 if travel == 0.0 else (after.radius - before.radius) / travel
        entry_x, entry_y = _entry_normal(station, slope)
        rho_stations.append(
            _RhoStation(
                cx=station.cx,
                cy=station.cy,
                radius=station.radius,
                clockwise=station.clockwise,
                tx=station.tx,
                ty=station.ty,
                wx=entry_x,
                wy=entry_y,
            )
        )
    return rho_stations


def _as_guide_station(station: _RhoStation) -> _GuideStation:
    """Adapt a `_RhoStation` to the probe-ring helpers, which take the older type.

    `_probe_positions` and `_station_is_admissible` read a station's centre,
    radius and ENTRY. The entry is the one field whose definition differs between
    the two types, and the difference is the point of this module, so it is
    carried across rather than recomputed: the returned `_GuideStation` is given
    the tangent whose guide normal reproduces this station's entry exactly.

    Args:
        station: The station to adapt.

    Returns:
        A `_GuideStation` with the same centre, radius, turn direction and entry.
    """
    # `_GuideStation.entry` offsets along (-ty, tx) when clockwise and (ty, -tx)
    # otherwise. Inverting that for a required entry direction w gives the tangent
    # below, so the adapted station's entry is w to the bit.
    tangent = (station.wy, -station.wx) if station.clockwise else (-station.wy, station.wx)
    return _GuideStation(
        cx=station.cx,
        cy=station.cy,
        radius=station.radius,
        clockwise=station.clockwise,
        tx=tangent[0],
        ty=tangent[1],
    )


def _first_trochoidal_station(stations: List[_RhoStation], tool_radius: float) -> Optional[int]:
    """Index of the first station clearing the degeneracy floor, or ``None`` if none does.

    Args:
        stations: The chain's ordered stations.
        tool_radius: Tool radius.

    Returns:
        The index, or ``None`` when every station on the chain is sub-threshold.
    """
    floor = DEGENERATE_LOOP_TOOL_RADII * tool_radius
    for index, station in enumerate(stations):
        if station.radius > floor:
            return index
    return None


@dataclass(frozen=True)
class DeclinedRegion:
    """A stretch of guide the generator refused to machine, and where it is.

    SILENT UNDER-CUT IS THE ONE FAILURE A ROUGHING GENERATOR MUST NOT HAVE. The
    degeneracy floor means some reachable material is left standing, and a count
    buried in a warning string is not something a caller can act on: it cannot be
    handed to a corner-clearing operation, drawn on a plot, or asserted on in a
    test. So the declined stretch is returned as data, with the geometry needed to
    machine it with a smaller tool.

    A region is a MAXIMAL RUN of consecutive sub-threshold stations, so a chain
    that runs out at both ends yields two of these rather than one span hiding a
    machined middle.

    Attributes:
        path_index: Index of the chain this run belongs to, matching the
            `path_index` stamped on that chain's operations.
        first_center: Centre of the run's first station, in world XY.
        last_center: Centre of the run's last station, in world XY.
        station_count: How many consecutive guide stations the run covers.
        largest_gouge_free_radius: The biggest machining-circle radius any station
            in the run admits. It is at or below the tool radius by construction --
            that is why the run was declined -- and it says how much smaller a tool
            would have to be to trochoid here: one whose radius is below it.
    """

    path_index: int
    first_center: Tuple[float, float]
    last_center: Tuple[float, float]
    station_count: int
    largest_gouge_free_radius: float


@dataclass
class RhoToolpathResult(ToolpathResult):
    """A `ToolpathResult` that also says WHICH reachable material it declined.

    Additive: every consumer of a `ToolpathResult` keeps working on one of these
    unchanged, and a caller that does not care never sees the extra field.

    Attributes:
        declined_regions: Every maximal run of guide stations left unmachined
            because no circle there would be a trochoid. Empty when the whole
            guide was machined.
    """

    declined_regions: Tuple[DeclinedRegion, ...] = ()


def _declined_regions(stations: List[_RhoStation], path_index: int, tool_radius: float) -> List[DeclinedRegion]:
    """Every maximal run of sub-threshold stations on one chain.

    Args:
        stations: The chain's ordered stations.
        path_index: Index stamped on this chain's operations.
        tool_radius: Tool radius, which fixes the degeneracy floor.

    Returns:
        One region per maximal run, in station order.
    """
    floor = DEGENERATE_LOOP_TOOL_RADII * tool_radius
    regions: List[DeclinedRegion] = []
    run: List[_RhoStation] = []
    for station in list(stations) + [None]:  # type: ignore[list-item]
        if station is not None and station.radius <= floor:
            run.append(station)
            continue
        if run:
            regions.append(
                DeclinedRegion(
                    path_index=path_index,
                    first_center=(run[0].cx, run[0].cy),
                    last_center=(run[-1].cx, run[-1].cy),
                    station_count=len(run),
                    largest_gouge_free_radius=max(entry.radius for entry in run),
                )
            )
            run = []
    return regions


def _radius_step_admissible(origin: _RhoStation, candidate: _RhoStation, tool_radius: float) -> bool:
    """Whether two consecutive machining circles' swept annuli still overlap radially.

    Args:
        origin: The last emitted station.
        candidate: The station under test.
        tool_radius: Tool radius, which scales the bound.

    Returns:
        ``True`` while the radius change is within
        `MAX_LOOP_RADIUS_STEP_TOOL_RADII` tool radii.
    """
    return abs(candidate.radius - origin.radius) <= MAX_LOOP_RADIUS_STEP_TOOL_RADII * tool_radius


def _largest_admissible_advance(
    stock: Stock,
    stations: List[_RhoStation],
    origin_index: int,
    window_end: int,
    regulation: _Regulation,
) -> Tuple[int, bool]:
    """SCAN the station window from the furthest end for the largest admissible advance.

    THE SCAN IS THE ALGORITHM AND IT MUST NOT BECOME A BISECTION. Engagement is
    not monotone in the advance, so the admissible advances are not an up-set in
    the station index and the pass/fail pattern down the window can alternate.
    Walking from the furthest candidate DOWNWARD and returning the first that
    passes is correct under any pattern, because station indices order the
    advances strictly and the first that passes therefore carries the largest
    admissible one. That argument is about the ORDER of the walk and nothing else.

    A candidate is admissible on TWO conditions, in this order:

    1. the radius step to it is within `MAX_LOOP_RADIUS_STEP_TOOL_RADII`, which is
       arithmetic on two doubles the guide already produced, and
    2. no evaluated position on its machining circle exceeds the cap, which is the
       exact `cap_exceeded` predicate at each of the probe positions.

    The cheap condition is tested first so the exact predicate is only paid for on
    candidates that could be taken. That is an evaluation ORDER, not a relaxation:
    a candidate refused by condition 1 would have to be refused anyway.

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
        if _station_is_admissible(stock, _as_guide_station(candidate), advance, regulation.tool_radius, regulation.cap_ratio):
            return index, False
    # No admissible advance exists -- the tool is entering virgin stock or crossing
    # a neck, where any motion at all exceeds the cap. Refusing to advance means
    # refusing to machine, so the minimum advance is taken and the caller is told.
    # The radius-step bound is NOT relaxed here: it is what keeps the emitted
    # sequence's annuli overlapping, and a forced advance that opened a gap would
    # trade a reported cap exceedance for an unreported patch of uncut material.
    minimum = origin_index + 1
    for index in range(minimum, window_end + 1):
        if _radius_step_admissible(origin, stations[index], regulation.tool_radius):
            return index, True
    return minimum, True


def _loop_operation(station: _RhoStation, cut_z: float, path_index: int) -> ToolpathOperation:
    """Build the `CUT` machining-circle operation for one accepted station.

    The circle's frame x-axis points at the entry tool-centre position, so
    ``Circle.point_at(0)`` is the entry point -- the convention the existing
    generators use and the one the audit's `subtract_arc_sweep` replay relies on.
    The declared tangents are the circle's OWN tangent at the entry
    (`_RhoStation.entry_tangent`), which is what makes the junction with a bridge
    on the external common tangent G1 rather than merely close.

    Args:
        station: The accepted station, carrying its radius and entry.
        cut_z: Cutting-plane height.
        path_index: Path index stamped on the operation.

    Returns:
        The machining-circle operation.

    Raises:
        DegenerateMachiningCircleError: If the station's radius vanished, leaving
            the entry point on the centre and the frame with no x-axis.
    """
    entry_x, entry_y = station.entry
    unit_x, unit_y = _unit_tangent(station.cx, station.cy, entry_x, entry_y)
    if (unit_x, unit_y) == (0.0, 0.0):
        raise DegenerateMachiningCircleError(
            f"Guide station at ({station.cx!r}, {station.cy!r}) has radius {station.radius!r}: its machining circle collapses to its centre and has no start point."
        )
    frame = Frame([station.cx, station.cy, cut_z], [unit_x, unit_y, 0.0], [-unit_y, unit_x, 0.0])
    tangent_x, tangent_y = station.entry_tangent
    tangent = np.array([tangent_x, tangent_y, 0.0], dtype=np.float64)
    return ToolpathOperation(
        geometry=Circle(station.radius, frame=frame),
        operation=OperationType.CUT,
        path_index=path_index,
        clockwise=station.clockwise,
        start_tangent=tangent,
        end_tangent=tangent,
    )


@dataclass(frozen=True)
class _ChainOutcome:
    """What one walked chain emitted and where it left the tool.

    Attributes:
        entry: Tool-centre point the chain plunged at.
        exit: Tool-centre point the chain retracted from.
        entry_over_cap: Whether the chain's first machining circle was emitted at
            positions the exact cap predicate refuses -- the tool meeting virgin
            stock, which is a full-immersion cut for any generator at any
            non-degenerate radius.
        forced_advances: Advances taken past a refusing predicate.
        sub_threshold_stations: Stations skipped because their gouge-free radius
            was at or below the degeneracy floor.
    """

    entry: Tuple[float, float]
    exit: Tuple[float, float]
    entry_over_cap: bool
    forced_advances: int
    sub_threshold_stations: int


def _machine_chain(
    stock: Stock,
    stations: List[_RhoStation],
    path_index: int,
    regulation: _Regulation,
    cut_z: float,
    operations: List[ToolpathOperation],
) -> Optional[_ChainOutcome]:
    """Walk one skeleton chain end to end, emitting and depleting as it goes.

    The walk starts at the first station clearing the degeneracy floor and ends at
    the last one. Sub-threshold stations at either end carry no machining circle:
    a circle there is a bore, and emitting one would report a trochoid the motion
    is not. The material at them is left to the neighbouring loops' sweeps, and
    the count of skipped stations is returned so the caller can say how much of
    the guide was treated that way rather than discovering it from a coverage
    grid.

    Args:
        stock: The depleting stock, mutated in place.
        stations: The chain's ordered stations.
        path_index: Path index stamped on this chain's operations.
        regulation: The validated parameters.
        cut_z: Cutting-plane height.
        operations: Output stream, appended to in place.

    Returns:
        The chain's outcome, or ``None`` when no station clears the floor.
    """
    first = _first_trochoidal_station(stations, regulation.tool_radius)
    if first is None:
        return None
    last = len(stations) - 1
    while last > first and stations[last].radius <= DEGENERATE_LOOP_TOOL_RADII * regulation.tool_radius:
        last -= 1
    skipped = first + (len(stations) - 1 - last)

    entry_station = stations[first]
    entry_x, entry_y = entry_station.entry
    operations.append(_line_operation((entry_x, entry_y), (entry_x, entry_y), regulation.clearance_z, cut_z, OperationType.PLUNGE, path_index))
    entry_over_cap = not _station_is_admissible(
        stock,
        _as_guide_station(entry_station),
        entry_station.entry_tangent,
        regulation.tool_radius,
        regulation.cap_ratio,
    )

    forced_advances = 0
    index = first
    while True:
        station = stations[index]
        loop_entry_x, loop_entry_y = station.entry
        operations.append(_loop_operation(station, cut_z, path_index))
        stock.subtract_arc_sweep_local(
            station.cx,
            station.cy,
            loop_entry_x,
            loop_entry_y,
            loop_entry_x,
            loop_entry_y,
            station.clockwise,
            regulation.tool_radius,
        )
        if index == last:
            break
        next_index, forced = _largest_admissible_advance(
            stock,
            stations,
            index,
            min(index + regulation.window, last),
            regulation,
        )
        forced_advances += int(forced)
        next_x, next_y = stations[next_index].entry
        operations.append(_line_operation((loop_entry_x, loop_entry_y), (next_x, next_y), cut_z, cut_z, OperationType.CUT, path_index))
        stock.subtract_capsule_quad(loop_entry_x, loop_entry_y, next_x, next_y, regulation.tool_radius)
        index = next_index

    exit_x, exit_y = stations[last].entry
    operations.append(_line_operation((exit_x, exit_y), (exit_x, exit_y), cut_z, regulation.clearance_z, OperationType.RETRACT, path_index))
    return _ChainOutcome(
        entry=(entry_x, entry_y),
        exit=(exit_x, exit_y),
        entry_over_cap=entry_over_cap,
        forced_advances=forced_advances,
        sub_threshold_stations=skipped,
    )


def rho_regulated_toolpath(
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
    """Trochoidal pocketing whose loop radius sequence is regulated, not just its advance.

    Walks the straight-skeleton guide chain by chain, as
    `compas_cgal.engagement_toolpath.engagement_controlled_toolpath` does, with
    three differences: a station whose gouge-free radius is at or below the tool
    radius carries no machining circle, the advance may not step over a radius
    change of more than `MAX_LOOP_RADIUS_STEP_TOOL_RADII` tool radii, and the
    bridge between two circles runs along their external common tangent so the
    junctions are G1 at any radius difference.

    WHAT THIS GUARANTEES, EXACTLY: engagement <= *tea_cap_deg*, decided by an
    exact predicate, at each EVALUATED tool position -- the loop entry and the
    probe ring on each accepted machining circle. It is NOT a continuous guarantee
    between evaluated positions, and no certificate is produced or returned. The
    added guarantees are structural rather than measured: every emitted circle has
    ``rho > r``, and consecutive circles' swept annuli overlap.

    Where the predicate refuses -- the tool meeting virgin stock, or a neck no
    advance can cross under the cap -- the circle is emitted anyway, because
    refusing to advance means refusing to machine, and an
    `UnavoidableEngagementWarning` reports how many and of which kind. A chain
    entry is a full-immersion cut at every non-degenerate radius, so it is
    expected there and counted rather than hidden.

    Args:
        polygon: Outer pocket boundary `Polygon` in the world XY plane.
        tool_diameter: Tool diameter; the tool radius is half of this.
        tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``.
        holes: Optional island polygons strictly inside *polygon*.
        guide_step_tool_diameters: Guide station spacing in tool diameters; also
            the resolution of every advance the search can pick.
        max_advance_tool_diameters: Largest advance considered, in tool diameters.
        radial_clearance: Safety clearance subtracted from each loop's available
            radius. Defaults to ``RADIAL_CLEARANCE_FRACTION * tool_diameter``.
        climb: ``True`` for climb milling (clockwise loops).
        cut_z: Z-height of the cutting plane.
        clearance_z: Z-height for rapid travel between chains.
        max_passes: Maximum number of skeleton chains the guide may emit.
        samples_per_radian: Tessellation density of the returned polyline.

    Returns:
        RhoToolpathResult: The typed operation stream and the tessellated
        visualisation polyline, in the same types the other generators emit, plus
        `RhoToolpathResult.declined_regions` naming every stretch of guide left
        unmachined because no circle there would be a trochoid. That field is the
        contract against silent under-cut: the material is reachable, this
        generator will not take it, and a caller can hand the regions to a
        corner-clearing pass with a smaller tool.

    Raises:
        InvalidEngagementCapDegreesError: If *tea_cap_deg* is not in ``(0, 180]``.
        NonPositiveToolDiameterError: If *tool_diameter* is not strictly positive.
        InvalidGuideResolutionError: If the guide step or advance bound leaves no
            integer window to scan.
        InvalidClearanceHeightError: If *clearance_z* is not above *cut_z*.
        EmptyGuideError: If the pocket admits no gouge-free trochoid at this tool.
        EmptyRhoGuideError: If the guide produced chains but no station on any of
            them clears the degeneracy floor, so the pocket admits circles that
            are all bores.
        ImplausibleClearanceSlopeError: If a chain's radius grows faster along its
            own arc length than a 1-Lipschitz clearance function can.

    Warns:
        UnavoidableEngagementWarning: When machining circles had to be emitted at
            positions the exact cap predicate refuses, or when guide stations were
            skipped for being below the degeneracy floor (never silently).
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
    over_cap_entries = 0
    forced_advances = 0
    sub_threshold_stations = 0
    non_trochoidal_chains = 0
    declined: List[DeclinedRegion] = []
    last_exit: Optional[Tuple[float, float]] = None
    path_index = 0

    for guide_chain in chains:
        stations = _rho_stations(guide_chain)
        pending = len(operations)
        outcome = _machine_chain(stock, stations, path_index, regulation, cut_z, operations)
        declined.extend(_declined_regions(stations, path_index, regulation.tool_radius))
        if outcome is None:
            non_trochoidal_chains += 1
            sub_threshold_stations += len(stations)
            continue
        if last_exit is not None:
            # Rebase: the chain wrote its plunge first, and the link into it has to
            # precede that plunge in the stream.
            operations.insert(pending, _line_operation(last_exit, outcome.entry, regulation.clearance_z, regulation.clearance_z, OperationType.LINK, path_index))
        over_cap_entries += int(outcome.entry_over_cap)
        forced_advances += outcome.forced_advances
        sub_threshold_stations += outcome.sub_threshold_stations
        last_exit = outcome.exit
        path_index += 1

    if last_exit is None:
        raise EmptyRhoGuideError(
            f"The straight-skeleton guide produced {len(chains)} chain(s) for tool_diameter={tool_diameter!r}, but no station on any of them "
            f"has a gouge-free radius above the {DEGENERATE_LOOP_TOOL_RADII * regulation.tool_radius!r} degeneracy floor: every machining circle "
            "the pocket admits at this tool size would be a bore rather than a trochoid."
        )

    if over_cap_entries or forced_advances or declined:
        warnings.warn(
            f"{over_cap_entries} chain-entry loop(s) were emitted at positions where the exact cap predicate reports "
            f"tea_cap_deg={tea_cap_deg} exceeded -- entering virgin stock is a full-immersion cut at every "
            f"non-degenerate radius -- alongside {forced_advances} forced minimum advance(s) where not even the "
            f"shortest admissible advance is under the cap. {sub_threshold_stations} guide station(s) across "
            f"{non_trochoidal_chains} wholly sub-threshold chain(s) and the ends of the rest carry no machining "
            "circle, their gouge-free radius being at or below the tool radius; the material at them is left to the "
            f"neighbouring sweeps and is not certified covered here. The {len(declined)} declined run(s) are returned "
            "as `RhoToolpathResult.declined_regions` with the geometry a corner-clearing pass needs.",
            UnavoidableEngagementWarning,
            stacklevel=2,
        )

    return RhoToolpathResult(
        operations=operations,
        polyline=_tessellate(operations, samples_per_radian),
        declined_regions=tuple(declined),
    )
