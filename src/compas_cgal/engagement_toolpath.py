"""Engagement-controlled trochoidal pocketing: advance chosen by an exact per-position predicate.

This generator replaces the geometric stepover proxy of
`compas_cgal.toolpath.trochoidal_mat_toolpath_circular` with a measured
engagement criterion. It walks the same straight-skeleton guide, but the advance
between machining circles is not a dialled-in distance: it is the largest advance
whose evaluated tool positions all pass the exact `_stock_2.engagement_at`
cap predicate against the depleting exact stock.

WHAT THIS GUARANTEES, EXACTLY
-----------------------------
Engagement <= `tea_cap_deg`, decided by an exact predicate, **at each evaluated
tool position**. That is the whole claim. It is NOT a continuous guarantee
between evaluated positions: the tool centre traverses a full circle and a bridge
segment between one evaluated position and the next, and nothing here bounds what
happens in between. No certificate is produced, none is returned, and no
`MotionWitness` / `CapRefutation` object exists on this path.

Two honest comparisons:

- Against Held's trochoidal engagement control, which bisects a floating-point
  engagement expression to a tolerance of 1e-3: **per evaluated position this is
  stronger** -- each verdict is an exact predicate on the exact arrangement, so
  there is no tolerance and no floating-point decision anywhere in the accept /
  reject path.
- Against a continuous partition of the motion (the certified path in
  `compas_cgal.adaptive`, which this module deliberately does not use): **this is
  weaker**. A continuous partition bounds every centre on the motion; this bounds
  only the positions it evaluated.

DESIGN
------
1. The guide comes from the existing generator run at a fine, uniform pitch:
   its emitted `Circle` cut operations are, per `path_index` and in emission
   order, the ordered skeleton-chain stations with their exact clearance-derived
   loop radii. Reusing it keeps the gouge-free radius derivation, the chain
   extraction, and the chain ordering in the one C++ implementation that already
   owns them, and gives this module an advance grid whose resolution is an
   integer count of stations rather than a float.
2. The walk over each chain is greedy: from the accepted station, bisect the
   integer station window for the largest admissible advance, emit the machining
   circle, deplete the stock, and continue to the chain's last station.
3. Candidate acceptance is decided at `LOOP_PROBE_ANGLES_DEG` positions on the
   candidate loop plus the loop's entry point (K = 4 evaluated positions), each
   decided by the exact `cap_exceeded` boolean.

Coverage is not part of the accept/reject rule and is not certified here: the
advance bound `MAX_ADVANCE_TOOL_DIAMETERS` keeps consecutive machining circles
overlapping, and the consequence is measured on a grid rather than proved.

Per the exact-kernel boundary doctrine (`docs/exactness.md`), the transcendental
cap crosses into exact-land exactly once, in `_cap_surrogate`, as the rational
squared-chord surrogate ``4*sin^2(theta/2)``. Probe *placement* is a modelling
choice evaluated in doubles -- choosing where to look is not a geometric decision;
what is found there is decided exactly.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

from compas_cgal import _stock_2  # type: ignore
from compas_cgal.stock import Stock
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

# Guide station spacing as a fraction of the tool diameter -- also the resolution
# of every advance this generator can choose. Derivation: one station is
# D/40 = r/20 of radial depth of cut. In the textbook radial-immersion relation
# TEA = 2*acos(1 - ae/r) the sensitivity at half-radius immersion is
# dTEA/dae = 2/(r*sqrt(1-(1-ae/r)^2)) ~ 2.31/r, so one station is worth at most
# ~0.12 rad ~ 7 deg of engagement. That is the cap headroom the quantisation
# leaves unused; it is spent in the safe direction, because the search only ever
# returns an advance whose evaluated positions passed.
GUIDE_STEP_TOOL_DIAMETERS = 0.025

# Largest advance the search will consider, in tool diameters. It carries two
# independent derivations that land on the same number:
#
# ENGAGEMENT. In the steady trochoidal regime the advance IS the radial depth of
# cut ae, and TEA = 2*acos(1 - ae/r) saturates at a full turn when ae = 2r = D. No
# advance at or beyond one tool diameter can satisfy any cap, so bisecting past it
# is wasted work, not lost capability.
#
# COVERAGE. Nothing in the accept/reject rule mentions coverage -- a candidate
# whose probes all sit in void is accepted whatever lies between it and the last
# accepted circle -- so the advance bound is what stops the walk from buying low
# engagement by skipping material. A machining circle of radius R sweeps the
# annulus [R-r, R+r], of width 2r = D; consecutive circles therefore always overlap
# radially while the advance stays at or below D. Measured consequence on a
# 6x4 and a 10x6 pocket with a 2 mm tool: no residual material survives further
# than 0.29 mm from a wall, against a tool radius of 1.0 mm, and less residue than
# the unregulated generator leaves at a comparable stepover.
MAX_ADVANCE_TOOL_DIAMETERS = 1.0

# Probe angles on a candidate machining circle, in degrees, measured from the
# ADVANCE direction (previous accepted centre -> candidate centre). With the loop
# entry point these are the K = 4 evaluated positions per candidate.
#
# Why these and not a uniform sweep: in the steady regime the material a new loop
# meets is a crescent at its outer rim whose radial thickness varies like
# a*cos(phi), phi measured from the advance direction. Engagement therefore peaks
# at phi = 0 -- that is where the loop bites deepest into uncut material -- and
# the backward half (cos phi < 0) lies inside the union of the preceding loops'
# swept annuli. The 0 deg probe catches the peak. The +/-60 deg flanks sit on the
# half-depth contour (cos 60 deg = 1/2): far enough off-axis to catch the two
# regimes the idealisation misses -- guide curvature rotating the crescent off
# the nominal advance direction, and clearance growing along the guide (corner
# spokes) lifting fresh material onto the flanks -- while staying out of the
# provably-swept backward half. Sampling 16+ positions per candidate would make
# generation cost dominate without testing a materially different regime.
LOOP_PROBE_ANGLES_DEG = (-60.0, 0.0, 60.0)

# Rise of the rapid-travel plane above the cutting plane, in tool diameters, used
# when the caller supplies no explicit clearance height. Scale-free rather than an
# absolute millimetre value; real fixturing overrides it.
CLEARANCE_RISE_TOOL_DIAMETERS = 1.0

# The guide is requested with its own stepover set to this multiple of its pitch,
# so the pitch (not the stepover) governs the station spacing and the guide grid
# comes out uniform. Two is enough: the clearance function along a skeleton chain
# is 1-Lipschitz, so the radius growth per station never exceeds the advance and
# advance + growth <= 2 * pitch.
GUIDE_STEPOVER_ADVANCE_MULTIPLE = 2.0

# The guide's own tessellated polyline is discarded -- only its typed operations
# are read -- so it is requested at the coarsest density the API accepts.
GUIDE_SAMPLES_PER_RADIAN = 1.0

# Tessellation density of the returned visualisation polyline, matching the
# default of `trochoidal_mat_toolpath_circular` so both generators visualise at
# the same fidelity.
POLYLINE_SAMPLES_PER_RADIAN = 10.0

# Resolution floor of the advance bisection, in guide stations. The search bracket
# is a pair of integer station indices, so halving terminates when the bracket
# collapses to adjacent indices -- a floor of exactly one station. There is no
# float tolerance in the loop and no iteration budget to tune: the bracket is at
# most MAX_ADVANCE_TOOL_DIAMETERS / GUIDE_STEP_TOOL_DIAMETERS wide, so bisection
# terminates in at most ceil(log2(window)) steps by construction.
ADVANCE_SEARCH_INDEX_FLOOR = 1

# Guide tangent used for a chain of a single station, where there is nowhere to
# advance to and the tangent is geometrically undefined. Fixed deterministically at
# +X: the tangent only orients the entry point around a circle that is traversed in
# full either way, so the choice moves where the loop starts and never what it
# removes. Pinning it keeps generation reproducible.
SINGLE_STATION_TANGENT = (1.0, 0.0)


class InvalidEngagementCapDegreesError(ValueError):
    """The engagement cap lies outside the admissible half-turn range ``(0, 180]`` degrees."""


class NonPositiveToolDiameterError(ValueError):
    """The tool diameter is not strictly positive (a real cutter has a positive radius)."""


class InvalidGuideResolutionError(ValueError):
    """The guide step or advance bound leaves the advance search no integer bracket to bisect."""


class InvalidClearanceHeightError(ValueError):
    """The rapid-travel plane is not strictly above the cutting plane."""


class EmptyGuideError(RuntimeError):
    """The straight-skeleton guide yielded no machining stations (pocket too small for the tool)."""


class DegenerateMachiningCircleError(ValueError):
    """A guide station's machining-circle radius vanished, so the circle has no reconstructable start point."""


class UnavoidableEngagementWarning(UserWarning):
    """Machining circles were emitted at positions the exact cap predicate reports as exceeding the cap."""


@dataclass(frozen=True)
class _GuideStation:
    """One ordered station on a skeleton-chain guide.

    Attributes:
        cx: X coordinate of the machining-circle centre (on the guide).
        cy: Y coordinate of the machining-circle centre (on the guide).
        radius: Machining-circle radius, derived by the C++ guide from the exact
            clearance at this centre less the tool radius and the radial
            clearance, so a tool centred on this circle stays gouge-free.
        clockwise: Turn direction of the machining circle (climb vs conventional).
        tx: X component of the unit guide tangent at this station.
        ty: Y component of the unit guide tangent at this station.
    """

    cx: float
    cy: float
    radius: float
    clockwise: bool
    tx: float
    ty: float

    @property
    def entry(self) -> tuple[float, float]:
        """Tool-centre point where the bridge meets the machining circle.

        Offset one loop radius along the guide normal, with the normal picked so
        that the circle's tangent there is the guide tangent: the bridge into and
        out of the loop is then tangent-continuous on a straight guide segment
        (and approximately so through a turn, the same G1 model the existing
        generator uses).
        """
        nx, ny = (-self.ty, self.tx) if self.clockwise else (self.ty, -self.tx)
        return self.cx + self.radius * nx, self.cy + self.radius * ny


def _cap_surrogate(tea_cap_deg: float) -> float:
    """Convert the engagement cap in degrees to its exact rational chord surrogate.

    BOUNDARY (`docs/exactness.md`, boundary doctrine): this is the ONE place the
    caller's transcendental cap crosses into exact-land, and it crosses as the
    dimensionless rational ``4*sin^2(theta/2)``. Delegating the conversion to
    `_stock_2.cap_chord_ratio` keeps it byte-identical to the surrogate the C++
    certifier uses. The sub-ulp gap between the surrogate and the angle the caller
    typed is documented API semantics, never an in-core correction constant.

    Args:
        tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``.

    Returns:
        The surrogate ``4 * sin(theta/2)**2`` in ``(0, 4]``.

    Raises:
        InvalidEngagementCapDegreesError: If *tea_cap_deg* is NaN or outside
            ``(0, 180]``. A single engaged run subtends at most a half turn before
            the ``> pi`` case is an exact orientation verdict, so a larger cap is
            meaningless.
    """
    if not (tea_cap_deg > 0.0 and tea_cap_deg <= 180.0):
        raise InvalidEngagementCapDegreesError(f"tea_cap_deg must be in (0, 180]; got {tea_cap_deg!r}.")
    return _stock_2.cap_chord_ratio(math.radians(tea_cap_deg))


def _unit_tangent(ax: float, ay: float, bx: float, by: float) -> tuple[float, float]:
    """Unit vector from ``a`` to ``b``, or the zero vector when the two coincide.

    The zero return is a genuine "no direction here" and every caller treats it as
    such: a vertical plunge or retract gets a zero tangent (the same convention the
    C++ generator writes for a degenerate tangent), a one-station chain falls back
    to `SINGLE_STATION_TANGENT`, and a vanished loop radius raises
    `DegenerateMachiningCircleError` rather than silently building a broken frame.
    """
    dx, dy = bx - ax, by - ay
    length = math.hypot(dx, dy)
    if length == 0.0:
        return 0.0, 0.0
    return dx / length, dy / length


def _guide_chains(
    polygon: Polygon,
    tool_diameter: float,
    guide_step: float,
    radial_clearance: float,
    climb: bool,
    max_passes: int,
    holes: list[Polygon] | None,
) -> list[list[_GuideStation]]:
    """Ordered skeleton-chain stations, read off the existing generator at a fine pitch.

    `trochoidal_mat_toolpath_circular` already extracts the interior straight
    skeleton, splits it into chains, orders the chains, and derives every
    machining-circle radius from an exact clearance query at its own centre. Its
    emitted `Circle` cut operations are therefore exactly the ordered station
    sequence this module needs, grouped by ``path_index`` in emission order.
    Reading them back is reuse of a validated walk; rebuilding chain adjacency
    from the unordered vertex set that `polygon_skeleton_clearance` returns would
    mean inventing a connectivity heuristic the C++ side does not need.

    The guide is requested with linking, leads, and the final retract switched off
    so the operation stream contains nothing but the stations and their bridges.

    Args:
        polygon: Outer pocket boundary in the world XY plane.
        tool_diameter: Tool diameter.
        guide_step: Uniform station spacing along the guide.
        radial_clearance: Safety clearance subtracted from each available radius.
        climb: ``True`` for climb milling (clockwise loops).
        max_passes: Maximum number of skeleton chains the guide may emit.
        holes: Optional island polygons strictly inside *polygon*.

    Returns:
        One list of `_GuideStation` per skeleton chain, each in walk order.

    Raises:
        EmptyGuideError: If the guide emitted no machining circles at all.
    """
    guide = trochoidal_mat_toolpath_circular(
        polygon,
        tool_diameter=tool_diameter,
        stepover=GUIDE_STEPOVER_ADVANCE_MULTIPLE * guide_step,
        pitch=guide_step,
        radial_clearance=radial_clearance,
        max_passes=max_passes,
        climb=climb,
        holes=holes,
        link_paths=False,
        optimize_order=True,
        lead_in=0.0,
        lead_out=0.0,
        retract_at_end=False,
        samples_per_radian=GUIDE_SAMPLES_PER_RADIAN,
    )

    # Group by RUNS of equal path_index rather than by value: emission order is
    # the machining order the guide chose, and a run boundary is a chain boundary.
    runs: list[list[tuple[float, float, float, bool]]] = []
    current_index = None
    for op in guide.operations:
        geometry = op.geometry
        if not isinstance(geometry, Circle):
            continue
        if op.path_index != current_index:
            runs.append([])
            current_index = op.path_index
        center = geometry.frame.point
        runs[-1].append((float(center[0]), float(center[1]), float(geometry.radius), bool(op.clockwise)))

    chains: list[list[_GuideStation]] = []
    for run in runs:
        count = len(run)
        stations: list[_GuideStation] = []
        for k, (cx, cy, radius, clockwise) in enumerate(run):
            # Central difference where it exists, one-sided at the chain ends: the
            # tangent only orients the loop entry point, so a local estimate is
            # what the model needs.
            before = run[max(0, k - 1)]
            after = run[min(count - 1, k + 1)]
            tx, ty = _unit_tangent(before[0], before[1], after[0], after[1])
            if (tx, ty) == (0.0, 0.0):
                tx, ty = SINGLE_STATION_TANGENT
            stations.append(_GuideStation(cx=cx, cy=cy, radius=radius, clockwise=clockwise, tx=tx, ty=ty))
        chains.append(stations)

    if not chains:
        raise EmptyGuideError(
            f"The straight-skeleton guide produced no machining circles for tool_diameter={tool_diameter!r}; the pocket admits no gouge-free trochoid at this tool size."
        )
    return chains


def _probe_positions(station: _GuideStation, advance: tuple[float, float]) -> list[tuple[float, float]]:
    """Tool-centre positions at which a machining circle is decided.

    `LOOP_PROBE_ANGLES_DEG` rotates the *advance* direction around the loop. The
    station's entry point is appended because it is the terminus of the bridge cut
    that immediately precedes the loop, and it is evaluated against the stock
    BEFORE the bridge is removed -- so that probe measures the material the linking
    cut itself runs into.

    Args:
        station: The station whose machining circle is being decided.
        advance: Unit direction of travel into this station.

    Returns:
        The evaluated tool-centre positions, entry point first.
    """
    dx, dy = advance
    positions = [station.entry]
    for degrees in LOOP_PROBE_ANGLES_DEG:
        angle = math.radians(degrees)
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        ux, uy = dx * cos_a - dy * sin_a, dx * sin_a + dy * cos_a
        positions.append((station.cx + station.radius * ux, station.cy + station.radius * uy))
    return positions


def _station_is_admissible(stock: Stock, station: _GuideStation, advance: tuple[float, float], tool_radius: float, cap_ratio: float) -> bool:
    """Whether every evaluated position on this machining circle is under the cap.

    Each position is decided by `_stock_2.engagement_at`'s ``cap_exceeded``
    boolean, which is an exact per-run predicate on the exact arrangement -- not a
    comparison of the reported doubles. ``gap_close_ratio`` stays at zero: gap
    closure is the pessimism a CONTINUOUS certifier needs to bridge its station
    spacing, and this generator makes no between-position claim to bridge.

    The station is measured against the stock as it stands, i.e. WITHOUT first
    removing the bridge that would carry the tool to it. That is deliberate and
    conservative in the safe direction: the un-removed bridge sliver leaves more
    material at the probes, so a station is refused slightly earlier than a
    bridge-aware model would refuse it, never later.

    Args:
        stock: The current stock (unmodified by this call).
        station: The station under test.
        advance: Unit direction of travel into this station.
        tool_radius: Tool radius.
        cap_ratio: The exact rational cap surrogate from `_cap_surrogate`.

    Returns:
        ``True`` if no evaluated position reports the cap exceeded.
    """
    raw = stock.raw
    for px, py in _probe_positions(station, advance):
        _total_tea, _max_run_tea, cap_exceeded = _stock_2.engagement_at(raw, px, py, tool_radius, cap_ratio, 0.0)
        if cap_exceeded:
            return False
    return True


def _largest_admissible_advance(
    stock: Stock,
    stations: list[_GuideStation],
    origin_index: int,
    window_end: int,
    tool_radius: float,
    cap_ratio: float,
) -> tuple[int, bool]:
    """Bisect the station window for the largest advance the cap predicate accepts.

    Engagement is taken to be monotone non-decreasing in advance distance: a
    longer advance exposes a thicker crescent of uncut material at the loop's
    forward rim, and (where clearance grows along the guide) a larger loop radius
    reaching further outward. This monotonicity is a stated MODELLING ASSUMPTION
    of the search, not a theorem -- what is exact is each individual verdict.
    Under it, bisection on the integer station bracket returns the largest
    admissible station index; without it, it returns an admissible one.

    Termination is structural: the bracket is a pair of integers narrowing by
    halving, so it collapses at `ADVANCE_SEARCH_INDEX_FLOOR` -- one guide station.
    No float tolerance participates.

    Args:
        stock: The current stock.
        stations: The chain's ordered stations.
        origin_index: Index of the last accepted station.
        window_end: Highest candidate index the search may consider (inclusive).
        tool_radius: Tool radius.
        cap_ratio: The exact rational cap surrogate.

    Returns:
        ``(index, forced)``: the chosen station index, and whether it was taken
        despite the predicate refusing it because no admissible advance exists.
    """
    origin = stations[origin_index]
    lo = origin_index + ADVANCE_SEARCH_INDEX_FLOOR
    hi = window_end
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        candidate = stations[mid]
        advance = _unit_tangent(origin.cx, origin.cy, candidate.cx, candidate.cy)
        if _station_is_admissible(stock, candidate, advance, tool_radius, cap_ratio):
            best = mid
            lo = mid + ADVANCE_SEARCH_INDEX_FLOOR
        else:
            hi = mid - ADVANCE_SEARCH_INDEX_FLOOR
    if best is None:
        # No admissible advance exists -- the regime where the tool is entering
        # virgin stock or crossing a neck, where any motion at all exceeds the cap.
        # Refusing to advance would mean refusing to machine, so the minimum
        # advance is emitted and the caller is told, never silently.
        return origin_index + ADVANCE_SEARCH_INDEX_FLOOR, True
    return best, False


def _line_operation(
    start: tuple[float, float],
    end: tuple[float, float],
    z_start: float,
    z_end: float,
    operation: OperationType,
    path_index: int,
) -> ToolpathOperation:
    """Build a `ToolpathOperation` for a straight move, tangents filled for XY travel."""
    tx, ty = _unit_tangent(start[0], start[1], end[0], end[1])
    tangent = np.array([tx, ty, 0.0], dtype=np.float64)
    return ToolpathOperation(
        geometry=Line([start[0], start[1], z_start], [end[0], end[1], z_end]),
        operation=operation,
        path_index=path_index,
        clockwise=False,
        start_tangent=tangent,
        end_tangent=tangent,
    )


def _loop_operation(station: _GuideStation, cut_z: float, path_index: int) -> ToolpathOperation:
    """Build the `CUT` machining-circle operation for one accepted station.

    The circle's frame x-axis points at the entry tool-centre position, so
    ``Circle.point_at(0)`` is the entry point -- the same convention the existing
    generator's circle reconstruction uses, which is what makes the audit's
    `subtract_arc_sweep` replay land on the same swept annulus.

    Raises:
        DegenerateMachiningCircleError: If the station's radius vanished, leaving
            the entry point on the centre and the frame with no x-axis.
    """
    ex, ey = station.entry
    ux, uy = _unit_tangent(station.cx, station.cy, ex, ey)
    if (ux, uy) == (0.0, 0.0):
        raise DegenerateMachiningCircleError(
            f"Guide station at ({station.cx!r}, {station.cy!r}) has radius {station.radius!r}: its machining circle collapses to its centre and has no start point."
        )
    frame = Frame([station.cx, station.cy, cut_z], [ux, uy, 0.0], [-uy, ux, 0.0])
    # Traversing the circle in its turn direction, the tangent at the entry point
    # is the guide tangent (rotating the entry normal back by a quarter turn),
    # which is what makes bridge-loop-bridge tangent-continuous on a straight guide.
    tangent = np.array([station.tx, station.ty, 0.0], dtype=np.float64)
    return ToolpathOperation(
        geometry=Circle(station.radius, frame=frame),
        operation=OperationType.CUT,
        path_index=path_index,
        clockwise=station.clockwise,
        start_tangent=tangent,
        end_tangent=tangent,
    )


def _tessellate(operations: list[ToolpathOperation], samples_per_radian: float) -> np.ndarray:
    """Sample the operation stream into one Nx3 visualisation polyline.

    Args:
        operations: The emitted operations, in order.
        samples_per_radian: Angular sampling density for circular motions.

    Returns:
        An ``(N, 3)`` ``float64`` array of points, consecutive duplicates dropped.
    """
    points: list[tuple[float, float, float]] = []
    # Every circular motion here is a full turn, so one sample count serves them all.
    samples = max(2, int(math.ceil(2.0 * math.pi * samples_per_radian)))

    def push(point) -> None:
        candidate = (float(point[0]), float(point[1]), float(point[2]))
        if not points or points[-1] != candidate:
            points.append(candidate)

    for op in operations:
        geometry = op.geometry
        if isinstance(geometry, Circle):
            for i in range(samples + 1):
                push(geometry.point_at(i / samples))
        else:
            push(geometry.start)
            push(geometry.end)

    if not points:
        return np.empty((0, 3), dtype=np.float64)
    return np.asarray(points, dtype=np.float64, order="C")


def _machine_chain(
    stock: Stock,
    stations: list[_GuideStation],
    path_index: int,
    tool_radius: float,
    cap_ratio: float,
    window: int,
    cut_z: float,
    clearance_z: float,
    operations: list[ToolpathOperation],
) -> tuple[int, int]:
    """Walk one skeleton chain end to end, emitting and depleting as it goes.

    Every chain is walked to its LAST station: the advance is at least one guide
    station, so the walk always terminates there rather than stopping early and
    silently leaving the chain's far end unmachined.

    The chain's FIRST machining circle has no advance to search -- it is where the
    tool enters -- but its probe positions are still evaluated, against the guide
    tangent as the direction of travel, so that every emitted machining circle has
    an exact verdict attached and an over-cap entry is counted rather than assumed.

    Args:
        stock: The depleting stock, mutated in place.
        stations: The chain's ordered stations.
        path_index: Path index stamped on this chain's operations.
        tool_radius: Tool radius.
        cap_ratio: The exact rational cap surrogate.
        window: Maximum advance in guide stations.
        cut_z: Cutting-plane height.
        clearance_z: Rapid-travel height.
        operations: Output stream, appended to in place.

    Returns:
        ``(over_cap_entries, forced_advances)``: whether this chain's entry loop
        exceeded the cap (0 or 1), and how many advances had to be forced past a
        refusing predicate.
    """
    last = len(stations) - 1
    entry_x, entry_y = stations[0].entry
    operations.append(_line_operation((entry_x, entry_y), (entry_x, entry_y), clearance_z, cut_z, OperationType.PLUNGE, path_index))

    entry_station = stations[0]
    over_cap_entries = int(not _station_is_admissible(stock, entry_station, (entry_station.tx, entry_station.ty), tool_radius, cap_ratio))

    forced_advances = 0
    index = 0
    while True:
        station = stations[index]
        ex, ey = station.entry
        operations.append(_loop_operation(station, cut_z, path_index))
        stock.subtract_arc_sweep(station.cx, station.cy, ex, ey, ex, ey, station.clockwise, tool_radius)
        if index == last:
            break

        next_index, forced = _largest_admissible_advance(stock, stations, index, min(index + window, last), tool_radius, cap_ratio)
        forced_advances += int(forced)
        nx, ny = stations[next_index].entry
        operations.append(_line_operation((ex, ey), (nx, ny), cut_z, cut_z, OperationType.CUT, path_index))
        stock.subtract_capsule(ex, ey, nx, ny, tool_radius)
        index = next_index

    exit_x, exit_y = stations[last].entry
    operations.append(_line_operation((exit_x, exit_y), (exit_x, exit_y), cut_z, clearance_z, OperationType.RETRACT, path_index))
    return over_cap_entries, forced_advances


def engagement_controlled_toolpath(
    polygon: Polygon,
    tool_diameter: float,
    tea_cap_deg: float,
    *,
    holes: list[Polygon] | None = None,
    guide_step_tool_diameters: float = GUIDE_STEP_TOOL_DIAMETERS,
    max_advance_tool_diameters: float = MAX_ADVANCE_TOOL_DIAMETERS,
    radial_clearance: float | None = None,
    climb: bool = True,
    cut_z: float = 0.0,
    clearance_z: float | None = None,
    max_passes: int = 1000,
    samples_per_radian: float = POLYLINE_SAMPLES_PER_RADIAN,
) -> ToolpathResult:
    """Trochoidal pocketing whose advance is regulated by an exact engagement predicate.

    Walks the straight-skeleton guide chain by chain. At every step the advance to
    the next machining circle is the largest one, on the guide's integer station
    grid, whose evaluated tool positions all report the engagement cap NOT exceeded
    against the depleting exact stock.

    WHAT THIS GUARANTEES, EXACTLY: engagement <= *tea_cap_deg*, decided by an exact
    predicate, at each EVALUATED tool position -- the loop entry point and the
    `LOOP_PROBE_ANGLES_DEG` probes on each accepted machining circle. It is NOT a
    continuous guarantee between evaluated positions: nothing here bounds
    engagement at the tool centres that lie between two probes, along the rest of a
    machining circle, or in the interior of a bridge cut. No certificate is
    produced or returned. Held's trochoidal engagement control bisects a
    floating-point engagement expression to a 1e-3 tolerance; per evaluated
    position this is stronger, because each verdict is exact with no tolerance
    anywhere in the accept/reject path. It is weaker than a continuous partition of
    each motion, which is not used here.

    Every emitted machining circle has its evaluated positions decided, including
    each chain's entry loop, which has no advance to search. Where the predicate
    refuses -- the tool meeting virgin stock, or a neck no advance can cross under
    the cap -- the circle is emitted anyway, because refusing to advance means
    refusing to machine, and an `UnavoidableEngagementWarning` reports how many and
    of which kind. Nothing over the cap is emitted silently.

    Args:
        polygon: Outer pocket boundary `Polygon` in the world XY plane.
        tool_diameter: Tool diameter; the tool radius is half of this.
        tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``. Converted
            once, at `_cap_surrogate`, to the exact rational surrogate
            ``4*sin^2(theta/2)`` that the exact predicate consumes.
        holes: Optional island polygons strictly inside *polygon*.
        guide_step_tool_diameters: Guide station spacing in tool diameters; also
            the resolution of every advance the search can pick.
        max_advance_tool_diameters: Largest advance the search considers, in tool
            diameters.
        radial_clearance: Safety clearance subtracted from each loop's available
            radius. Defaults to ``RADIAL_CLEARANCE_FRACTION * tool_diameter``.
        climb: ``True`` for climb milling (clockwise loops), ``False`` for
            conventional milling.
        cut_z: Z-height of the cutting plane.
        clearance_z: Z-height for rapid travel between chains. Defaults to
            ``cut_z + CLEARANCE_RISE_TOOL_DIAMETERS * tool_diameter``.
        max_passes: Maximum number of skeleton chains the guide may emit.
        samples_per_radian: Tessellation density of the returned polyline.

    Returns:
        ToolpathResult: The typed operation stream -- `PLUNGE`, `CUT` machining
        circles and bridge lines, `RETRACT`, and clearance-height `LINK` moves
        between chains -- plus the tessellated visualisation polyline. The same
        types the existing generator emits, so `audit_toolpath_engagement` and
        every downstream consumer work unchanged.

    Raises:
        InvalidEngagementCapDegreesError: If *tea_cap_deg* is not in ``(0, 180]``.
        NonPositiveToolDiameterError: If *tool_diameter* is not strictly positive.
        InvalidGuideResolutionError: If the guide step or advance bound leaves no
            integer bracket to bisect.
        InvalidClearanceHeightError: If *clearance_z* is not above *cut_z*.
        EmptyGuideError: If the pocket admits no gouge-free trochoid at this tool.
        InvalidPolygonError: If *polygon* or a hole is non-planar or degenerate.

    Warns:
        UnavoidableEngagementWarning: When machining circles had to be emitted at
            positions the exact cap predicate refuses (never silently).
    """
    cap_ratio = _cap_surrogate(tea_cap_deg)

    if not tool_diameter > 0.0:
        raise NonPositiveToolDiameterError(f"tool_diameter must be strictly positive; got {tool_diameter!r}.")
    if not guide_step_tool_diameters > 0.0:
        raise InvalidGuideResolutionError(f"guide_step_tool_diameters must be strictly positive; got {guide_step_tool_diameters!r}.")
    if not max_advance_tool_diameters > 0.0:
        raise InvalidGuideResolutionError(f"max_advance_tool_diameters must be strictly positive; got {max_advance_tool_diameters!r}.")

    window = int(max_advance_tool_diameters // guide_step_tool_diameters)
    if window < 2 * ADVANCE_SEARCH_INDEX_FLOOR:
        raise InvalidGuideResolutionError(
            f"guide_step_tool_diameters={guide_step_tool_diameters!r} against max_advance_tool_diameters="
            f"{max_advance_tool_diameters!r} gives an advance window of {window} station(s); at least "
            f"{2 * ADVANCE_SEARCH_INDEX_FLOOR} are needed for the bisection to have a bracket. Use a finer guide step."
        )

    if radial_clearance is None:
        radial_clearance = RADIAL_CLEARANCE_FRACTION * tool_diameter
    if clearance_z is None:
        clearance_z = cut_z + CLEARANCE_RISE_TOOL_DIAMETERS * tool_diameter
    if not clearance_z > cut_z:
        raise InvalidClearanceHeightError(f"clearance_z must be strictly above cut_z; got clearance_z={clearance_z!r}, cut_z={cut_z!r}.")

    tool_radius = 0.5 * tool_diameter
    chains = _guide_chains(
        polygon,
        tool_diameter,
        guide_step_tool_diameters * tool_diameter,
        radial_clearance,
        climb,
        max_passes,
        holes,
    )

    stock = Stock(polygon, holes=holes)
    operations: list[ToolpathOperation] = []
    over_cap_entries = 0
    forced_advances = 0
    for path_index, stations in enumerate(chains):
        if path_index > 0:
            previous_exit = chains[path_index - 1][-1].entry
            next_entry = stations[0].entry
            operations.append(_line_operation(previous_exit, next_entry, clearance_z, clearance_z, OperationType.LINK, path_index))
        chain_entries, chain_forced = _machine_chain(
            stock,
            stations,
            path_index,
            tool_radius,
            cap_ratio,
            window,
            cut_z,
            clearance_z,
            operations,
        )
        over_cap_entries += chain_entries
        forced_advances += chain_forced

    if over_cap_entries or forced_advances:
        warnings.warn(
            f"{over_cap_entries + forced_advances} machining circle(s) were emitted at positions where the exact cap "
            f"predicate reports tea_cap_deg={tea_cap_deg} exceeded: {over_cap_entries} chain-entry loop(s), where the "
            f"tool first meets virgin stock and a first cut is a full slot by construction, and {forced_advances} "
            "forced minimum advance(s), where not even a one-station advance is admissible. Advancing less is not an "
            "option there, so they are emitted and counted rather than dropped.",
            UnavoidableEngagementWarning,
            stacklevel=2,
        )

    return ToolpathResult(operations=operations, polyline=_tessellate(operations, samples_per_radian))
