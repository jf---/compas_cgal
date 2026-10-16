"""Narrow-station rules, ported from the rho-regulated generator into typed code (Ruling 1.1-C7).

Three rules, each a port of `compas_cgal.engagement_rho_toolpath`. The suite
holds them against rho's public outputs (`tests/pipeline/test_against_rho.py`);
bit-equality with rho's privates is a logged scratch probe (Ruling 1.1-B-oracle):

* DEGENERACY FLOOR. A tool of radius ``r`` on a circle of radius ``rho`` sweeps
  ``[rho - r, rho + r]``, an annulus with an uncut core exactly when
  ``rho > r``. The floor is decided first on the exact guide station radius;
  only stations above it are authored as `IncidentLoop`, and the floor is then
  re-decided on the authored radius (Ruling 1.1-M3), because authoring rounds
  the radius onto its chart quantum and may land on ``r``.
* RADIUS-STEP BOUND. Consecutive annuli overlap radially only while
  ``|rho - rho'| <= 2r``; the advance never steps over that, forced or not,
  and where nothing in the window keeps it the chain is split instead
  (Ruling 1.1-B-unbridgeable). Decided exactly over the rationals of the doubles.
* EXTERNAL-TANGENT ENTRY. The entry radius ``w = -slope * t + sqrt(1 - slope^2) * m``
  aims the bridge between consecutive entries at the circles' external common
  tangent, so junctions are approximately tangent by construction. That is not
  certified: no predicate checks the junction tangency, and certifying it is
  Task 1.6 work.

`advance_ladder` is the narrow instance of `compas_cgal.pipeline.ladder`: a
furthest-first scan over the floor-clearing stations of the window with injected
conditions, the radius-step bound always first; when none passes, the minimum
advance keeping the bound is taken and returned as a `ForcedAdvance` record.

Every station index a ledger record carries is in guide emission order
(`OrientedChain.guide_index`); walk order stays inside the walk.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Self

from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import EmptyLadderWindowError
from compas_cgal.pipeline.errors import ImplausibleClearanceSlopeError
from compas_cgal.pipeline.errors import SubFloorStationError
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.ladder import CandidateCondition
from compas_cgal.pipeline.ladder import NoAdmissibleCandidate
from compas_cgal.pipeline.ladder import furthest_first
from compas_cgal.pipeline.result import DeclineCause
from compas_cgal.pipeline.result import DeclinedRun
from compas_cgal.pipeline.result import ForcedAdvance
from compas_cgal.stock import Stock
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# The floor in tool radii: rho > r is the boundary at which the swept annulus
# loses its hole (physics, not a tuned value; benchmarks.quality's
# DEGENERATE_LOOP_RATIO measures against the same boundary).
DEGENERATE_LOOP_TOOL_RADII = 1.0

# Largest radius change an advance may step over, in tool radii: annuli
# [rho - r, rho + r] and [rho' - r, rho' + r] overlap radially only while
# |rho - rho'| <= 2r. The bound is the overlap limit itself.
MAX_LOOP_RADIUS_STEP_TOOL_RADII = 2

# Largest clearance slope d(rho)/d(s) the external tangent accepts: the external
# common tangent of two circles exists only while |rho - rho'| <= D, i.e. slope
# <= 1, and a clearance function is 1-Lipschitz, so a larger slope is not a
# clearance and the guide is refused rather than clamped.
MAX_CLEARANCE_SLOPE = 1.0

# Largest advance considered, in tool diameters: the advance is the radial depth
# of cut in the steady regime, and TEA saturates at a full turn at ae = 2r = D.
MAX_ADVANCE_TOOL_DIAMETERS = 1.0


def advance_window(guide_step_tool_diameters: float) -> int:
    """The advance window in whole guide stations.

    Floor division of the two doubles, as the frozen generators compute it
    (``1.0 // 0.025 == 39.0``, one station short of the nominal 40): kept
    verbatim so the ported walk scans the same window.

    Args:
        guide_step_tool_diameters: Guide spacing in tool diameters.

    Returns:
        The window.
    """
    return int(MAX_ADVANCE_TOOL_DIAMETERS // guide_step_tool_diameters)


@dataclass(frozen=True)
class NarrowStation:
    """A guide station carrying the external-tangent entry direction.

    Attributes:
        centre: Machining-circle centre, on the guide.
        radius: Gouge-free guide radius (the exact station radius).
        clockwise: Turn direction.
        tangent: Unit guide tangent in the walk direction.
        entry_direction: Unit entry-radius direction ``w``, centre to entry.
    """

    centre: Point2[WorldXY]
    radius: GuideRadius
    clockwise: bool
    tangent: Direction2[WorldXY]
    entry_direction: Direction2[WorldXY]

    @property
    def entry(self) -> tuple[float, float]:
        """The double-constructed entry ``c + rho * w`` (placement; the exact entry is the authored loop's)."""
        return self.centre.x + self.radius.value * self.entry_direction.x, self.centre.y + self.radius.value * self.entry_direction.y

    @property
    def entry_tangent(self) -> tuple[float, float]:
        """Unit travel direction at the entry, a quarter turn from ``w`` in the turn direction."""
        wx, wy = self.entry_direction.x, self.entry_direction.y
        if self.clockwise:
            return wy, -wx
        return -wy, wx


def entry_normal(station: GuideStation, slope: float) -> tuple[float, float]:
    """Unit entry-radius direction putting the bridge on the external common tangent.

    ``w = -slope * t + sqrt(1 - slope^2) * m`` with ``m`` the guide normal on the
    loop's material side; at ``slope == 0`` it is the guide normal exactly.

    A CONSTRUCTION only: plausibility of the slope is decided exactly by
    `clearance_slope` before this runs. The rounded slope of a plausible pair can
    still exceed one by an ulp, so the lateral term is guarded against a negative
    radicand (a construction guard, never a decision).

    Args:
        station: The station (tangent and turn direction).
        slope: ``d(rho)/d(s)`` at the station, decided plausible.

    Returns:
        ``(wx, wy)``.
    """
    tx, ty = station.tangent.x, station.tangent.y
    material_x, material_y = (-ty, tx) if station.clockwise else (ty, -tx)
    lateral = math.sqrt(max(0.0, 1.0 - slope * slope))
    return -slope * tx + lateral * material_x, -slope * ty + lateral * material_y


def clearance_slope(before: GuideStation, after: GuideStation) -> float:
    """The radius slope between two stations, refused exactly when no clearance can have it.

    The DECISION ``|rho_after - rho_before| <= MAX_CLEARANCE_SLOPE * |c_after - c_before|``
    is taken over the rationals of the doubles (squared, so no root is formed).
    The returned slope is the double construction rho uses: the radius
    difference over the double ``hypot`` of the centre difference, ``0`` when the
    centres coincide (a single station, or a pair with equal radius).

    Args:
        before: The earlier station.
        after: The later station.

    Returns:
        The slope.

    Raises:
        ImplausibleClearanceSlopeError: The pair's radius changes faster than its
            centres move (or the centres coincide with different radii).
    """
    rise = Fraction(after.radius.value) - Fraction(before.radius.value)
    run_x = Fraction(after.centre.x) - Fraction(before.centre.x)
    run_y = Fraction(after.centre.y) - Fraction(before.centre.y)
    if rise * rise > Fraction(MAX_CLEARANCE_SLOPE) ** 2 * (run_x * run_x + run_y * run_y):
        raise ImplausibleClearanceSlopeError(
            f"Guide stations at ({before.centre.x!r}, {before.centre.y!r}) and ({after.centre.x!r}, {after.centre.y!r}) change radius by "
            f"{float(rise)!r} over a centre distance of {math.hypot(float(run_x), float(run_y))!r}, above the {MAX_CLEARANCE_SLOPE} a 1-Lipschitz "
            "clearance can have; the external common tangent does not exist."
        )
    travel = math.hypot(after.centre.x - before.centre.x, after.centre.y - before.centre.y)
    return 0.0 if travel == 0.0 else (after.radius.value - before.radius.value) / travel


def narrow_stations(chain: GuideChain) -> tuple[NarrowStation, ...]:
    """Re-express a chain with external-tangent entries.

    The slope is the central difference of the radius over the guide's arc
    length, one-sided at the ends; a single station or coincident neighbours
    have no slope direction and read ``0``.

    Args:
        chain: The chain, in walk order.

    Returns:
        One narrow station per guide station.

    Raises:
        ImplausibleClearanceSlopeError: A slope exceeds `MAX_CLEARANCE_SLOPE`.
    """
    stations = chain.stations
    count = len(stations)
    out: list[NarrowStation] = []
    for index, station in enumerate(stations):
        slope = clearance_slope(stations[max(0, index - 1)], stations[min(count - 1, index + 1)])
        wx, wy = entry_normal(station, slope)
        out.append(
            NarrowStation(
                centre=station.centre,
                radius=station.radius,
                clockwise=station.clockwise,
                tangent=station.tangent,
                entry_direction=Direction2[WorldXY].build(wx, wy),
            )
        )
    return tuple(out)


def clears_floor(radius: float, tool: ToolRadius) -> bool:
    """The degeneracy floor ``rho > r`` on an exact double radius.

    Args:
        radius: Station or authored radius, mm.
        tool: Tool radius.

    Returns:
        True iff the circle sweeps an annulus with an uncut core.
    """
    return radius > DEGENERATE_LOOP_TOOL_RADII * tool.value


def station_loop(station: NarrowStation, tool: ToolRadius) -> IncidentLoop:
    """Author the station's loop with an exactly incident entry at the external-tangent phase.

    Args:
        station: A station whose radius clears the floor.
        tool: Tool radius.

    Returns:
        The authored loop (its radius may still fall onto the floor; callers
        re-decide the floor on it).

    Raises:
        SubFloorStationError: The station radius does not clear the floor, so
            authoring would make a bore (and may not be authorable at all).
    """
    if not clears_floor(station.radius.value, tool):
        raise SubFloorStationError(f"station at ({station.centre.x!r}, {station.centre.y!r}) has radius {station.radius.value!r} <= the floor {tool.value!r}.")
    phase = Radian(math.atan2(station.entry_direction.y, station.entry_direction.x))
    return IncidentLoop.build(station.centre, station.radius, phase)


def radius_step_admissible(origin: GuideRadius, candidate: GuideRadius, tool: ToolRadius) -> bool:
    """``|rho' - rho| <= 2r``, decided exactly over the rationals of the three doubles.

    Args:
        origin: Radius advanced from.
        candidate: Radius advanced to.
        tool: Tool radius.

    Returns:
        True iff consecutive annuli still overlap radially.
    """
    return abs(Fraction(candidate.value) - Fraction(origin.value)) <= MAX_LOOP_RADIUS_STEP_TOOL_RADII * Fraction(tool.value)


@dataclass(frozen=True)
class NarrowCandidate:
    """A floor-clearing station offered to the ladder.

    Attributes:
        index: Station index in walk order.
        station: The station.
        loop: Its authored loop (radius above the floor).
    """

    index: int
    station: NarrowStation
    loop: IncidentLoop


@dataclass(frozen=True)
class RadiusStepCondition:
    """Ladder condition: the radius-step bound from origin to candidate.

    Attributes:
        tool: Tool radius.
    """

    tool: ToolRadius

    def __call__(self, stock: Stock, origin: NarrowCandidate, candidate: NarrowCandidate, /) -> bool:
        """Whether the radius step stays within the bound (*stock* unused)."""
        del stock
        return radius_step_admissible(origin.station.radius, candidate.station.radius, self.tool)


@dataclass(frozen=True)
class OrientedChain:
    """One chain in one walk direction, with its floor decided and its loops authored.

    Attributes:
        stations: Narrow stations in walk order.
        loops: The authored loop per station, ``None`` where the floor declines it.
        origin: The chain's emission index.
        reversed_walk: Whether the walk runs against the guide's order.
        first: First floor-clearing station.
        last: Last floor-clearing station.
        declined: Maximal declined runs, in guide emission order.
        tool: Tool radius the floor was decided for.
    """

    stations: tuple[NarrowStation, ...]
    loops: tuple[IncidentLoop | None, ...]
    origin: int
    reversed_walk: bool
    first: int
    last: int
    declined: tuple[DeclinedRun, ...]
    tool: ToolRadius

    @classmethod
    def build(cls, chain: GuideChain, tool: ToolRadius) -> Self | None:
        """Decide the floor on *chain* and author its loops.

        Args:
            chain: The chain in its walk direction.
            tool: Tool radius.

        Returns:
            The oriented chain, or ``None`` if no station clears the floor.

        Raises:
            ImplausibleClearanceSlopeError: A slope exceeds `MAX_CLEARANCE_SLOPE`.
        """
        stations = narrow_stations(chain)
        loops = tuple(_floor_loop(station, tool) for station in stations)
        kept = [index for index, loop in enumerate(loops) if loop is not None]
        if not kept:
            return None
        return cls(
            stations=stations,
            loops=loops,
            origin=chain.origin,
            reversed_walk=chain.reversed_walk,
            first=kept[0],
            last=kept[-1],
            declined=_declined_runs(loops, chain.origin, chain.reversed_walk),
            tool=tool,
        )

    @property
    def entry_station(self) -> NarrowStation:
        """The station this orientation begins machining at."""
        return self.stations[self.first]

    def guide_index(self, index: int) -> int:
        """The guide-emission-order index of the walk-order station *index* (one numbering for the ledger)."""
        return len(self.stations) - 1 - index if self.reversed_walk else index

    def candidate(self, index: int) -> NarrowCandidate:
        """The floor-clearing station at *index* as a ladder candidate.

        Args:
            index: Station index in walk order.

        Returns:
            The candidate.

        Raises:
            SubFloorStationError: The station is declined by the floor.
        """
        loop = self.loops[index]
        if loop is None:
            raise SubFloorStationError(f"station {index} of chain {self.origin} is declined by the floor.")
        return NarrowCandidate(index=index, station=self.stations[index], loop=loop)


def _floor_loop(station: NarrowStation, tool: ToolRadius) -> IncidentLoop | None:
    """The authored loop, or ``None`` if the floor declines the station or its authored radius."""
    if not clears_floor(station.radius.value, tool):
        return None
    loop = station_loop(station, tool)
    if not clears_floor(loop.radius.value, tool):
        return None
    return loop


def _declined_runs(loops: tuple[IncidentLoop | None, ...], chain: int, reversed_walk: bool) -> tuple[DeclinedRun, ...]:
    """Maximal runs of declined stations, indexed and ordered in guide emission order."""
    in_guide_order = tuple(reversed(loops)) if reversed_walk else loops
    runs: list[DeclinedRun] = []
    start: int | None = None
    for index, loop in enumerate(in_guide_order):
        if loop is None and start is None:
            start = index
        elif loop is not None and start is not None:
            runs.append(DeclinedRun(chain=chain, first_station=start, last_station=index - 1))
            start = None
    if start is not None:
        runs.append(DeclinedRun(chain=chain, first_station=start, last_station=len(in_guide_order) - 1))
    return tuple(runs)


def orientations(chain: GuideChain, tool: ToolRadius) -> tuple[OrientedChain, ...]:
    """The chain's walk directions that carry a trochoid.

    Both directions normally; the authored floor is per direction (the entry
    phase differs), so a chain whose only stations sit on the floor's quantum
    may be offered in one.

    Args:
        chain: The chain in the guide's direction.
        tool: Tool radius.

    Returns:
        Forward first, then reversed, each only if it clears the floor.
    """
    built = (OrientedChain.build(chain, tool), OrientedChain.build(chain.reversed(), tool))
    return tuple(oriented for oriented in built if oriented is not None)


@dataclass(frozen=True)
class Advance:
    """The ladder's answer when some station keeps the radius-step bound.

    Attributes:
        index: Station advanced to, in walk order.
        forced: The forced-decision record when no candidate passed every condition.
    """

    index: int
    forced: ForcedAdvance | None


@dataclass(frozen=True)
class ChainSplit:
    """The ladder's answer when no floor-clearing station of the window keeps the radius-step bound.

    The walk never breaks the bound (Ruling 1.1-B-unbridgeable): it stops at the
    origin, the uncut gap is recorded, and the chain is re-entered at `resume`
    as a forced entry.

    Attributes:
        origin: Station the walk stops at, in walk order.
        resume: Nearest floor-clearing station beyond *origin*, in walk order.
        gap: The declined record of the gap, cause `DeclineCause.RADIUS_STEP`.
    """

    origin: int
    resume: int
    gap: DeclinedRun


def advance_ladder(
    stock: Stock,
    chain: OrientedChain,
    origin: int,
    window_end: int,
    conditions: tuple[CandidateCondition[NarrowCandidate], ...],
) -> Advance | ChainSplit:
    """The furthest admissible advance from *origin*, the forced minimum one, or a split.

    The radius-step bound is structural (Ruling 1.1-B-boundinvariant): it is
    always evaluated first, before the injected *conditions*, so no caller can
    drop it. Only floor-clearing stations are candidates. When none passes every
    condition, the nearest floor-clearing station keeping the bound is taken and
    recorded as `ForcedAdvance`; when none keeps the bound either, the chain is
    split there (`ChainSplit`), never bridged across the bound.

    Args:
        stock: The current stock (never mutated here).
        chain: The oriented chain.
        origin: Station advanced from (floor-clearing), in walk order.
        window_end: Highest candidate index (inclusive), in walk order.
        conditions: Further ladder conditions, in evaluation order.

    Returns:
        The advance or the split.

    Raises:
        EmptyLadderWindowError: ``window_end <= origin`` or *origin* is the
            chain's last floor-clearing station (nothing to advance to).
        SubFloorStationError: *origin* is declined by the floor.
    """
    if window_end <= origin or origin >= chain.last:
        raise EmptyLadderWindowError(f"advance window ({origin}, {window_end}] on chain {chain.origin} is empty.")
    start = chain.candidate(origin)
    window = [index for index in range(window_end, origin, -1) if chain.loops[index] is not None]
    bounded = (RadiusStepCondition(chain.tool), *conditions)
    chosen = furthest_first(stock, start, tuple(chain.candidate(index) for index in window), bounded)
    if not isinstance(chosen, NoAdmissibleCandidate):
        return Advance(index=chosen.index, forced=None)
    for index in reversed(window):
        if radius_step_admissible(start.station.radius, chain.stations[index].radius, chain.tool):
            forced = ForcedAdvance(chain=chain.origin, reversed_walk=chain.reversed_walk, origin_station=chain.guide_index(origin), station=chain.guide_index(index))
            return Advance(index=index, forced=forced)
    resume = next(index for index in range(origin + 1, chain.last + 1) if chain.loops[index] is not None)
    first, last = sorted((chain.guide_index(origin), chain.guide_index(resume)))
    return ChainSplit(origin=origin, resume=resume, gap=DeclinedRun(chain=chain.origin, first_station=first, last_station=last, cause=DeclineCause.RADIUS_STEP))
