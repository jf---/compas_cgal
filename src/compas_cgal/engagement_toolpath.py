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
between evaluated positions: the tool centre traverses the arc between one
evaluated position and the next, and a whole bridge segment between loops, and
nothing here bounds what happens in between. Raising `LOOP_PROBE_COUNT` narrows
that gap; it does not close it, and the residue is measurable rather than
hypothetical -- see the constant's comment. No certificate is produced, none is
returned, and no `MotionWitness` / `CapRefutation` object exists on this path.

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
   candidate loop plus the loop's entry point (K = ``LOOP_PROBE_COUNT + 1``
   evaluated positions), each decided by the exact `cap_exceeded` boolean.

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

# Tool-centre positions evaluated on a candidate machining circle, spaced
# uniformly around the loop from the ADVANCE direction (previous accepted centre
# -> candidate centre). With the loop's entry point these are the
# K = LOOP_PROBE_COUNT + 1 evaluated positions per candidate.
#
# THE DERIVATION THIS REPLACES WAS FALSIFIED, and the falsification is recorded
# here rather than quietly dropped. Until 2026-08-21 the probes were the triple
# (-60, 0, +60) deg, justified by the steady-regime idealisation in which the
# material a new loop meets is a crescent at its forward rim of radial thickness
# ~a*cos(phi): engagement peaks at phi = 0, and "the backward half lies inside
# the union of the preceding loops' swept annuli". The second half of that claim
# is FALSE. Measured on a 20x12 pocket, 2 mm tool, by walking every machining
# circle the generator ACCEPTED at 32 tool-centre positions on the same depleting
# stock: at an 80 deg cap all 120 over-cap positions lie between -30 and -150 deg
# from the advance direction and NONE lies at 0, +30, +60, +90 or +120; at a
# 100 deg cap all 24 do. On the worst accepted circle at an 80 deg cap the peak is
# 101.7 deg at -135 deg from the advance -- and re-measuring the OLD three probes
# against that same stock returns 39.1 deg, so the gap was probe PLACEMENT, not a
# difference between the generation and audit depletion models. The forward peak
# is real (54 of 90 circles peak at phi = 0); it is simply not where the cap is
# broken.
#
# WHY UNIFORM AND NEVER ONE-SIDED. The loaded quadrant is the TRAILING-lateral
# one, on the side fixed by the loop's turn direction, and it mirrors exactly with
# the milling direction: the same pocket and cap that puts 4/16/52/36/12 over-cap
# positions in the (-30,-60,-90,-120,-150) bins under climb milling puts
# 12/36/52/16/4 in the (+150,+120,+90,+60,+30) bins under conventional. Any
# placement biased to one side is therefore tuned to one winding and blind on the
# other, which is a worse failure than the one being fixed. A uniform ring is the
# only winding-agnostic placement, and it needs no empirical tuning to stay
# correct when the guide curves or the clearance grows.
#
# WHY THIS COUNT: MEASURED CONVERGENCE, NOT A DERIVATION. There is a geometric
# floor -- consecutive evaluated positions are a chord 2*R*sin(pi/K) apart, so
# their tool disks only overlap at all while that chord stays under 2*r, which on
# the largest loops in this corpus (R = 4.998, r = 1.0) needs K >= 16 -- but the
# measurement says the floor is not enough, so the shipped value is the measured
# one. Worst engagement an INDEPENDENT 60-position walk (phase-offset by half a
# step, so it shares no grid with the probes) finds on ACCEPTED circles, 20x12,
# 2 mm tool:
#
#   K       3(old)   8     12    16    24    32    40    48
#   cap 40   86.7  71.7  55.3  57.8  49.3  43.4  43.4  43.7
#   cap 80  101.7  99.2  87.3  88.6  82.8  81.7  82.2  82.4
#
# The last count at which either column moves is 32; 40 and 48 buy nothing and
# cost 24% more. Generation on that pocket goes 1.35 s -> 5.54 s at a 40 deg cap
# and 0.19 s -> 0.39 s at a 120 deg cap.
#
# WHAT THIS DOES NOT DO. Raising the density NARROWS the sampling gap; it does
# not close it. Nothing here bounds engagement between two evaluated positions,
# and the residue is measurable: at a 40 deg cap on that pocket the independent
# walk still finds 43.4 deg on a circle every probe accepted.
LOOP_PROBE_COUNT = 32

# The probe angles themselves, in degrees from the advance direction. Derived
# from the count rather than written out so the two can never disagree.
LOOP_PROBE_ANGLES_DEG = tuple(360.0 * index / LOOP_PROBE_COUNT for index in range(LOOP_PROBE_COUNT))

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
class _Regulation:
    """Validated, unit-resolved parameters shared by the engagement-regulated generators.

    Built only through `build`, which owns the whole error model of the parameter
    seam: every raw constructor field is already validated and already in absolute
    units, so a walk that holds one of these never re-checks a caller's numbers and
    never re-derives a default. The two generators in this package
    (`engagement_controlled_toolpath` and
    `compas_cgal.engagement_radial_toolpath.radius_regulated_toolpath`) share it so
    that a parameter can only ever mean the same thing to both.

    Attributes:
        cap_ratio: The exact rational cap surrogate ``4*sin^2(theta/2)`` from
            `_cap_surrogate` -- the only form of the cap that reaches a predicate.
        tool_radius: Tool radius in model units.
        guide_step: Guide station spacing in model units.
        radial_clearance: Safety clearance subtracted from each available radius.
        clearance_z: Rapid-travel plane height.
        window: Largest advance the search may consider, in whole guide stations.
    """

    cap_ratio: float
    tool_radius: float
    guide_step: float
    radial_clearance: float
    clearance_z: float
    window: int

    @classmethod
    def build(
        cls,
        *,
        tool_diameter: float,
        tea_cap_deg: float,
        guide_step_tool_diameters: float,
        max_advance_tool_diameters: float,
        radial_clearance: float | None,
        cut_z: float,
        clearance_z: float | None,
    ) -> "_Regulation":
        """Validate the caller's parameters once and resolve them to absolute units.

        Args:
            tool_diameter: Tool diameter; the tool radius is half of this.
            tea_cap_deg: Engagement-angle cap in degrees, in ``(0, 180]``.
            guide_step_tool_diameters: Guide station spacing in tool diameters.
            max_advance_tool_diameters: Largest advance considered, in tool diameters.
            radial_clearance: Safety clearance, or ``None`` for
                ``RADIAL_CLEARANCE_FRACTION * tool_diameter``.
            cut_z: Z-height of the cutting plane.
            clearance_z: Rapid-travel height, or ``None`` for
                ``cut_z + CLEARANCE_RISE_TOOL_DIAMETERS * tool_diameter``.

        Returns:
            The validated parameters.

        Raises:
            InvalidEngagementCapDegreesError: If *tea_cap_deg* is not in ``(0, 180]``.
            NonPositiveToolDiameterError: If *tool_diameter* is not strictly positive.
            InvalidGuideResolutionError: If the guide step or advance bound leaves
                the advance search no integer bracket to bisect.
            InvalidClearanceHeightError: If *clearance_z* is not above *cut_z*.
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

        return cls(
            cap_ratio=cap_ratio,
            tool_radius=0.5 * tool_diameter,
            guide_step=guide_step_tool_diameters * tool_diameter,
            radial_clearance=radial_clearance,
            clearance_z=clearance_z,
            window=window,
        )


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

    `LOOP_PROBE_ANGLES_DEG` rotates the *advance* direction around the loop, so
    the returned ring is uniform and its phase is the direction of travel. The
    station's entry point is prepended because it is the terminus of the bridge
    cut that immediately precedes the loop, and it is evaluated against the stock
    BEFORE the bridge is removed -- so that probe measures the material the linking
    cut itself runs into.

    SAMPLED, NOT BOUNDED. These positions are where the loop is ASKED about; the
    answer at each of them is exact. The constant's comment carries what the
    density buys and what it does not: a denser ring narrows the gap between two
    evaluated positions, and nothing here closes it.

    Args:
        station: The station whose machining circle is being decided.
        advance: Unit direction of travel into this station.

    Returns:
        The evaluated tool-centre positions, entry point first, then the ring in
        increasing angle from the advance direction.
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
        stock.subtract_arc_sweep_local(station.cx, station.cy, ex, ey, ex, ey, station.clockwise, tool_radius)
        if index == last:
            break

        next_index, forced = _largest_admissible_advance(stock, stations, index, min(index + window, last), tool_radius, cap_ratio)
        forced_advances += int(forced)
        nx, ny = stations[next_index].entry
        operations.append(_line_operation((ex, ey), (nx, ny), cut_z, cut_z, OperationType.CUT, path_index))
        # The bridge's swept capsule, removed as two exact end disks plus the
        # rectangle between them rather than as a chain of hundreds of disks.
        # Same under-covering contract, same slack budget -- but a bounded number
        # of curves per bridge, so neither the bridge's length nor the arrangement
        # it has already accumulated drives the cost.
        stock.subtract_capsule_quad(ex, ey, nx, ny, tool_radius)
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
    `LOOP_PROBE_ANGLES_DEG` ring on each accepted machining circle. It is NOT a
    continuous guarantee between evaluated positions: nothing here bounds
    engagement at the tool centres that lie between two probes, along the rest of a
    machining circle, or in the interior of a bridge cut. Raising
    `LOOP_PROBE_COUNT` NARROWS that gap and does not close it -- on a 20x12 pocket
    with a 2 mm tool at a 40 deg cap, an independent 60-position walk still finds
    43.4 deg on a circle every probe accepted. No certificate is produced or
    returned. Held's trochoidal engagement control bisects a floating-point
    engagement expression to a 1e-3 tolerance; per evaluated position this is
    stronger, because each verdict is exact with no tolerance anywhere in the
    accept/reject path. It is weaker than a continuous partition of each motion,
    which is not used here.

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
    regulation = _Regulation.build(
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        guide_step_tool_diameters=guide_step_tool_diameters,
        max_advance_tool_diameters=max_advance_tool_diameters,
        radial_clearance=radial_clearance,
        cut_z=cut_z,
        clearance_z=clearance_z,
    )
    cap_ratio = regulation.cap_ratio
    window = regulation.window
    clearance_z = regulation.clearance_z
    tool_radius = regulation.tool_radius
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
