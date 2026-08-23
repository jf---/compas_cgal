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
maximum, and the search takes the LARGEST ADMISSIBLE one. Where the coarse grid
finds nothing it is rescanned at `RADIUS_LADDER_SUBDIVISIONS` times the
resolution, because the coarse spacing is the ADVANCE quantisation and can step
straight over a band of radii that comply; and where even that finds nothing, the
circle emitted is the gentlest of the refused candidates BY MEASUREMENT, never
the first or last by index.

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

WHERE NO RADIUS COMPLIES
------------------------
Some stations cannot be cut within the cap at any radius: a chain entry into
virgin stock is surrounded by material at every radius, and so is a corner tip
whose whole clearance disk is narrower than the tool. Refusing to cut there is
not an option -- the guide runs through it -- so a circle the predicate refuses is
emitted, counted, and warned about rather than hidden.

WHICH circle is a measurement, not an index. `_least_bad_rung` ranks the refused
candidates by their REPORTED peak engagement and returns the mildest; ties go to
the maximal circle, so virgin-stock entries still emit a full slot and still
finish in one pass. This is the one place in the module where a reported double
influences a choice, and it is admissible because the exact predicate has already
refused every option being ranked -- the ordering cannot promote one back into the
guarantee.

IT IS NOT FREE, and the cost is structural rather than incidental. A circle
smaller than the station's maximum does not FINISH the station, so the maximal
circle still follows on a later sweep: cutting back where nothing complies buys a
gentler worst circle at the price of an extra circle and an extra sweep. Measured
on a 6x4 pocket at a 40 deg cap -- three tool diameters wide, so the largest loop
the clearance allows is the tool radius itself -- ranking alone takes the worst
machining circle from 86.4 to 54.5 deg while taking the circles the audit finds
over the cap from 34 to 41. `RADIUS_LADDER_REFINEMENT_MARGIN` is where that trade
is balanced, and its comment carries the table.

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
from typing import FrozenSet
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
from compas_cgal.engagement_toolpath import _measured_peak_engagement
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
# rung whose emission FINISHES a station, and it is the rung the ranking falls back
# to when every candidate measures the same.
FULL_RADIUS_RUNG = 0

# Passed as `_radius_ladder`'s minimum radius where the caller wants none. The
# coarse ladder has no floor of its own: it stops on the largest radius already
# emitted, and its last rung is whatever the maximal radius leaves modulo the
# step. Named so that "this ladder is unfloored" reads as a decision rather than
# as a bare zero someone might mistake for a tolerance.
NO_RADIUS_FLOOR = 0.0

# Sub-intervals each coarse ladder interval is split into for the SECOND scan, run
# only where the first one found nothing admissible.
#
# WHY A SECOND SCAN EXISTS AT ALL. The coarse spacing above is the ADVANCE
# quantisation, whose derivation prices one step against the tool radius. That
# pricing is wrong for the radius knob at a station whose maximal radius is itself
# a fraction of the tool radius, and it is wrong by enough to step over the answer.
# Measured on the 20x12 pocket at a 60 deg cap, station (18.482, 10.482), maximal
# radius 0.5156, coarse step 0.05: rung 6 (radius 0.2156) measures 61.3 deg and
# still cuts, rung 7 (radius 0.1656) measures 5.9 deg and cuts NOTHING -- one step
# crosses from over-cap to idle. Sampling that same interval finely shows a band of
# radii from 0.1719 to 0.2123, about four fifths of a coarse step wide, where the
# loop BOTH complies and cuts. The coarse ladder straddles it, reports the station
# unsalvageable, and forces a 107 deg circle where a 59 deg one exists.
#
# WHY THIS COUNT: MEASURED, NOT DERIVED -- the same footing as `LOOP_PROBE_COUNT`,
# and for the same reason: what has to be resolved is the width of a compliant band
# in a non-monotone function, which no closed form in this module bounds. Worst
# machining-circle engagement an independent 32-position walk finds away from the
# chain entries on the 20x12 pocket, against the subdivision count:
#
# Measured on the 20x12 pocket at a 60 deg cap, tool diameter 2.0, over the 244
# machining circles that are NOT chain entries, each probed at 16 positions:
#
#   N                        1      2      4      8     16
#   worst engagement, deg  88.6   88.6   88.6   88.6   88.6
#   circles over the cap     12     12      8      8      8
#
# The peak does not move; the COUNT is what the resolution buys, and it converges
# at N = 4. Eight is one doubling of headroom past that, at no measured cost.
#
# READ THE FIRST COLUMN CAREFULLY: N = 1 is not "refinement off". The rescan
# branch runs whenever `_least_bad_rung` leaves the station within
# `RADIUS_LADDER_REFINEMENT_MARGIN` of the cap, and no constant disables it --
# setting this to 1 only collapses the refined ladder onto the coarse one while
# the RANKING still applies, which is why N = 1 and N = 2 agree. Deleting the
# branch outright is a different experiment and a much worse one: 126.1 deg worst
# engagement, measured separately.
#
# AND DO NOT KEEP THIS WHILE DROPPING `_least_bad_rung`. Measured on 6x4 at a
# 40 deg cap, refinement without the ranking takes the over-cap count from 34 to
# 1 and emits a 213.6 deg circle -- a slotting cut that breaks tools. The refined
# scan rescues a station by DEFERRING it, and the maximal circle that must still
# follow lands in a worse state; the ranking is what caps the tail, the
# refinement is what caps the count, and neither is safe alone. Both together
# UNGATED is also worse than the ranking alone on both pockets, which is why the
# gate exists.
#
# TERMINATION IS THE INTEGER, not a tolerance: the refined ladder is a fixed-length
# list of `(rungs - 1) * N + 1` radii built once, and the scan walks it. No
# convergence test, no float comparison, nothing to tune.
#
# THE REFINED LADDER SPANS EXACTLY THE COARSE ONE -- from the top coarse rung to
# the bottom one, no further. It resolves WITHIN the set of radii the module
# already considers worth emitting; it does not extend that set downward. That
# floor is load-bearing rather than incidental: at station (18.941, 10.941) on the
# same pocket, whose coarse ladder is the two rungs 0.0568 and 0.0068, a fine sweep
# does eventually find a complying radius -- at 0.0011, a circle roughly a
# thousandth of the tool radius. Emitting that is not a lighter cut, it is a
# degenerate motion that also leaves the station unfinished and buys another sweep.
# Lowering the minimum useful circle is a separate change with its own evidence.
RADIUS_LADDER_SUBDIVISIONS = 8

# Smallest radius the refined ladder offers below rung 0, as a multiple of the
# COARSE ladder step. A loop of radius R opens an annulus of width 2R around the
# bore it sits in, so a floor of half a step is exactly "a rung must open at least
# one coarse step of annulus" -- one quantum of radial depth of cut, the same
# quantity `RADIUS_LADDER_STEP_TOOL_DIAMETERS` is measured in. Below that a circle
# removes less than the search can resolve while still leaving the station
# unfinished, which costs a whole extra sweep to come back for the rest.
#
# MEASURED, on the 20x12 pocket at a 60 deg cap with the gate below in place:
# floors of 0.25 and 0.5 are indistinguishable (8 machining circles over the cap),
# and a floor of one whole step loses a rescue (12). So the rationale and the
# measurement agree, and the value sits inside the flat region rather than on its
# edge.
RADIUS_LADDER_FLOOR_STEPS = 0.5

# How far over the cap a station's gentlest available circle may measure and still
# earn the second, finer scan -- as a multiple of the cap the caller asked for.
#
# WHAT IT GATES, AND WHAT IT CANNOT. This decides how hard to LOOK, never what is
# found: every verdict on every candidate is the same exact `cap_exceeded`
# predicate whether the gate opened or not, and a station the gate skips is
# emitted exactly as it would have been with no refinement at all. It is a
# comparison of two doubles -- a REPORTED engagement against the caller's own
# transcendental cap -- and that is admissible precisely because no emission
# depends on it.
#
# WHY A GATE AND NOT ALWAYS. The two halves of this search pull against each
# other, and the gate is where they are balanced. Cutting a station back below its
# maximal circle does not FINISH it (`_RadiusChoice.finishes`), so the maximal
# circle still has to follow on a later sweep: a cut-back at a station where
# nothing complies buys a gentler worst circle and pays for it with an extra
# circle and an extra sweep. Refining unconditionally takes that trade everywhere,
# including at stations whose rescue radius is a small fraction of the coarse rung
# above it, where the extra sweeps cost more engagement than the rescue saves.
#
# MEASURED, on the 6x4 pocket at a 40 deg cap -- the hard case, a pocket three tool
# diameters wide where the largest loop the clearance allows is the tool radius
# itself. Worst machining-circle engagement away from the chain entries, circles
# the dense audit finds over the cap, and total cutting length:
#
#   gate      1.25    1.4     1.5    1.75    2.0
#   worst     54.5   61.1    76.3    76.3  110.0
#   over cap    40     23      18      17     11
#   length     413    593     661     704    835
#
# against 86.4 deg / 34 circles / 295 for this generator with no refinement and no
# ranking at all. 1.4 is the smallest gate at which BOTH quality columns beat that
# baseline; below it the count regresses, above it the worst engagement climbs back
# and the path keeps growing. On the 20x12 pocket at a 60 deg cap every gate from
# 1.1 up gives the same 88.6 deg / 12 circles / 3382, against 126.1 / 12 / 3236 for
# the same baseline -- so that pocket does not constrain the value and this one
# fixes it.
#
# THE PATH GETS LONGER, and that is the trade being made rather than an oversight:
# 6x4 at a 40 deg cap goes from 295 to 593 units of cutting travel. A smaller
# radial bite taken more times is what respecting the cap costs, which is the same
# statement the module docstring makes about tighter caps.
RADIUS_LADDER_REFINEMENT_MARGIN = 1.4

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


@dataclass
class RadialToolpathResult(ToolpathResult):
    """A `ToolpathResult` that also says WHICH of its circles the cap predicate refused.

    Every consumer of a `ToolpathResult` keeps working on one of these unchanged;
    the extra field is additive, and a caller that does not care about it never
    sees it.

    WHY IT IS ON THE RESULT AND NOT ONLY IN THE WARNING. `UnavoidableEngagementWarning`
    reports HOW MANY circles were emitted over the cap, and until now that was the
    only way to find out. A count cannot answer the question that actually matters
    once the search picks a reduced radius where nothing complies: is a given
    over-cap circle one the search HAD NO CHOICE about, or one it chose badly?
    Those two were indistinguishable while every forced emission was the station's
    maximal circle, and they are not any more. Recovering the answer by parsing a
    warning string would be a worse API than a field, and re-deriving it by
    replaying the path would be a second implementation of the search.

    Attributes:
        forced_loops: Indices into `operations` of the machining circles emitted
            at positions the exact cap predicate REFUSES at every candidate radius
            -- chain entries into virgin stock included. These are the circles the
            accompanying `UnavoidableEngagementWarning` counts. An over-cap circle
            outside this set is a search defect; one inside it is the least-bad
            available cut, and `_least_bad_rung` carries how "least bad" is
            decided.
    """

    forced_loops: FrozenSet[int] = frozenset()


@dataclass(frozen=True)
class _SweepOutcome:
    """What one pass over a skeleton chain emitted and how hard it had to work.

    Attributes:
        operations: The pass's operation stream -- plunge, machining circles and
            bridges, retract -- or empty when the pass found nothing left to cut.
        entry: Tool-centre point the pass plunged at, or ``None`` when empty.
        exit: Tool-centre point the pass retracted from, or ``None`` when empty.
        forced_radii: Machining circles emitted over the cap because NO candidate
            radius was admissible at their station.
        forced_loops: WHERE those circles are -- indices into `operations`. The
            count alone cannot distinguish a circle the search had no choice about
            from one it chose badly, and that distinction is the whole content of
            the ranking, so it is carried rather than recomputed.
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
    forced_loops: Tuple[int, ...]
    forced_advances: int
    first_loop_forced: bool


@dataclass(frozen=True)
class _GentlestRung:
    """The mildest of a station's already-refused candidate radii, with its measurement.

    The peak travels with the rung because the two have exactly one consumer each
    and both are reporting quantities: the rung says which circle to emit, the
    peak says how far over the cap that circle sits, which is what
    `_largest_admissible_radius` reads to decide whether looking harder is worth
    it. Recomputing the peak at that call site would mean walking the probe ring a
    second time for a number that was just measured.

    Attributes:
        rung: Index into the ladder the ranking was run on.
        peak: The reported engaged-run angle in radians that ranked it there;
            ``0.0`` for an empty candidate set.
    """

    rung: int
    peak: float


@dataclass(frozen=True)
class _RadiusChoice:
    """The loop radius picked at one station, and how the search arrived at it.

    A radius rather than a rung index, because the refined ladder's radii do not
    sit on integer multiples of the guide step and reconstructing them from an
    index would mean the caller re-deriving the search's own grid.

    Attributes:
        radius: The chosen loop radius in model units.
        finishes: Whether this is the station's maximal circle. Only that circle
            sweeps the whole annulus the clearance allows there, so only its
            emission leaves the station with nothing further to cut.
        forced: Whether the exact cap predicate refused every candidate and this
            one was ranked out of the refused set by `_least_bad_rung`.
    """

    radius: float
    finishes: bool
    forced: bool

    @classmethod
    def of(cls, radius: float, station: _GuideStation, *, forced: bool) -> "_RadiusChoice":
        """Build a choice, deriving `finishes` from the station it was chosen at.

        The maximal-circle test is float IDENTITY, not proximity, and it is sound
        because both ladders carry the station's own ``radius`` value through
        unmodified at rung 0 (`_radius_ladder`, `_refined_radius_ladder`). Deriving
        the flag here rather than at each call site keeps that invariant in one
        place.

        Args:
            radius: The chosen loop radius, taken from a ladder built at *station*.
            station: The station the radius was chosen at.
            forced: Whether the exact cap predicate refused every candidate.

        Returns:
            The choice.
        """
        return cls(radius=radius, finishes=radius == station.radius, forced=forced)


def _radius_ladder(full_radius: float, emitted_radius: float, ladder_step: float, rungs: int, floor_radius: float) -> List[float]:
    """Candidate loop radii at one station, largest first.

    Rung ``k`` is ``full_radius - k * ladder_step``. The ladder stops at the
    largest radius already emitted at this station, because a loop at or below it
    sweeps an annulus this station has already swept and would remove nothing:
    that floor is what makes every emission strictly increase the station's
    cleared radius, which is what makes the sweep loop terminate structurally
    rather than by a budget. It stops as well at *floor_radius*, the smallest
    radius the caller considers worth emitting at all, and after *rungs* entries.

    Rung 0 is always offered when anything is left, even when it lies below one
    ladder step or below *floor_radius*, so a station whose clearance admits only
    a hair of a circle still gets the circle the advance-only generator would have
    emitted there.

    Args:
        full_radius: The station's clearance-derived maximal radius.
        emitted_radius: Largest radius already emitted at this station; ``0.0``
            before the first pass.
        ladder_step: Spacing between rungs in model units.
        rungs: Hard cap on the number of entries, the ladder's only termination
            bound that is not a comparison between two radii.
        floor_radius: Smallest radius offered below rung 0; ``0.0`` for no floor.

    Returns:
        The candidate radii in descending order, empty when nothing is left.
    """
    if not full_radius > emitted_radius:
        return []
    ladder = [full_radius]
    for rung in range(1, rungs):
        radius = full_radius - rung * ladder_step
        if radius <= emitted_radius or radius < floor_radius:
            break
        ladder.append(radius)
    return ladder


def _refined_radius_ladder(full_radius: float, emitted_radius: float, ladder_step: float) -> List[float]:
    """The same ladder at `RADIUS_LADDER_SUBDIVISIONS` times the resolution.

    Same span, same rung 0 -- the station's maximal radius bit-for-bit, so the
    maximal-circle test that finishes a station stays an identity -- and the coarse
    rungs are a subset, sitting at indices that are multiples of the subdivision
    count. The extra rungs are the ones the coarse grid stepped over.

    THE FLOOR IS THE REFINED STEP ITSELF, and it is the whole reason this is a
    ladder builder rather than a subdivision of the coarse list. Two measurements
    on the 20x12 pocket at a 60 deg cap fix it from both sides:

    - station (18.800, 10.800), maximal radius 0.1980. The largest radius that
      both complies and cuts is 0.0429. The coarse ladder's rungs there are
      0.1980, 0.1480, 0.0980, 0.0480 -- it stops one rung ABOVE the answer, on the
      arbitrary remainder of the maximal radius modulo the step, and forces a
      117.8 deg circle. So the refinement has to be allowed BELOW the coarse
      ladder's last rung; subdividing its intervals is not enough.
    - station (18.941, 10.941), maximal radius 0.0568. Here a radius that complies
      and cuts also exists -- at 0.00142, a circle a thousandth of the tool radius.
      Emitting that is not a lighter cut, it is a degenerate motion that removes a
      1.4 micron annulus, leaves the station unfinished, and buys another sweep to
      come back for the rest. So the refinement must NOT run to zero.

    One refined step is the natural stop between those two: it is the smallest
    radius this grid can distinguish from no circle at all, so a rung below it is
    a radius the search could not have resolved in the first place.

    Args:
        full_radius: The station's clearance-derived maximal radius.
        emitted_radius: Largest radius already emitted at this station.
        ladder_step: The COARSE spacing; the refined one is derived from it.

    Returns:
        The refined candidate radii in descending order.
    """
    return _radius_ladder(
        full_radius,
        emitted_radius,
        ladder_step / RADIUS_LADDER_SUBDIVISIONS,
        rungs=RADIUS_LADDER_RUNGS * RADIUS_LADDER_SUBDIVISIONS,
        floor_radius=ladder_step * RADIUS_LADDER_FLOOR_STEPS,
    )


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


def _least_bad_rung(
    stock: Stock,
    ladder: List[float],
    station: _GuideStation,
    advance: Tuple[float, float],
    regulation: _Regulation,
) -> "_GentlestRung":
    """Rank the ALREADY-REFUSED rungs and return the gentlest one.

    WHEN THIS RUNS the exact predicate has refused every rung: there is no
    admissible radius at this station, and refusing to cut is not one of the
    options -- the material is in the way of a guide the walk has to get past. The
    only question left is WHICH inadmissible circle to emit, and it has to be
    answered by measurement, because no index rule answers it.

    WHY NOT RUNG 0 (what this replaces). Returning the station's maximal circle
    was returning the WORST candidate wherever the load falls off with the radius.
    Measured on the 20x12 pocket at a 60 deg cap, station (18.482, 10.482),
    maximal radius 0.5156: the ladder's reported peaks descend 107.0, 99.0, 91.2,
    82.8, 75.2, 68.6, 61.3 deg over rungs 0 to 6, so rung 0 was emitted at 107 deg
    where rung 6 was available at 61.3.

    WHY NOT THE SMALLEST RADIUS EITHER, and why not any index rule. Engagement is
    not monotone in the radius -- the property this whole module is built around --
    so the gentlest rung is generally INTERIOR. On a station whose swept void is
    NARROWER THAN THE TOOL, where no radius gets the tool clear of the material
    (2 mm tool, void radius 0.8, maximal radius 1.5), the reported peak falls from
    300.9 deg at rung 0 to 253.7 deg at rung 18 and RISES back to a full 360 deg by
    rung 26: first, last, smallest and largest all miss it.
    `tests/test_engagement_radial_toolpath.py` pins that state.

    WHY A DOUBLE MAY DECIDE THIS ONE THING. Everything the cap governs is settled
    before this function is called, by the exact `cap_exceeded` predicate, and it
    said NO to every candidate here. This ranking cannot promote a refused
    candidate to an accepted one -- the caller flags whatever comes back as forced
    either way -- so the reported `max_run_tea` doubles are ordering options that
    are already outside the guarantee. That is the deciding/reporting split of
    `docs/exactness.md` used exactly as written, and it is the ONE place in this
    module where a reported number influences a choice.

    CANDIDATES are the rungs that still reach uncut stock, plus rung 0
    unconditionally: a rung that cuts nothing is not a lesser evil, it is a wasted
    motion that also leaves the station unfinished, and rung 0 is the rung whose
    emission FINISHES the station (`_radial_sweep`), so it must always be
    reachable. Ties go to the LOWEST rung index -- the largest radius -- so a
    chain entry into virgin stock, where every rung measures a full turn, still
    emits the maximal circle and still finishes in one pass, exactly as before.

    Args:
        stock: The current stock (unmodified by this call).
        ladder: The station's candidate radii, largest first, as built by
            `_radius_ladder`; must be non-empty.
        station: The station under test, carrying its maximal radius.
        advance: Unit direction of travel into this station, for probe placement.
        regulation: The validated parameters.

    Returns:
        The gentlest candidate: its index into *ladder*, and the reported peak
        engagement that ranked it there.
    """
    best = _GentlestRung(rung=FULL_RADIUS_RUNG, peak=0.0)
    seen = False
    for rung, radius in enumerate(ladder):
        candidate = replace(station, radius=radius)
        if rung != FULL_RADIUS_RUNG and not _loop_reaches_material(stock, candidate, advance, regulation.tool_radius):
            continue
        peak = _measured_peak_engagement(stock, candidate, advance, regulation.tool_radius, regulation.cap_ratio)
        if not seen or peak < best.peak:
            best, seen = _GentlestRung(rung=rung, peak=peak), True
    return best


def _first_admissible_rung(
    stock: Stock,
    ladder: List[float],
    station: _GuideStation,
    advance: Tuple[float, float],
    regulation: _Regulation,
) -> Optional[int]:
    """SCAN one ladder from the top for the largest radius that both complies and cuts.

    THE SCAN IS THE ALGORITHM, and it must not become a bisection. Engagement is
    not monotone in the loop radius (module docstring, with the measurement), so
    the admissible rungs are not an up-set and the pass/fail pattern down the
    ladder can alternate. Walking from rung 0 downward and returning the FIRST
    pass is correct under any pattern: rung indices order the radii strictly
    downward, so the first rung that passes carries the largest admissible radius,
    whatever the rungs below it do. That argument is about the ORDER of the list
    and nothing else, which is why it holds unchanged on a refined ladder.

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

    Args:
        stock: The current stock (unmodified by this call).
        ladder: Candidate radii, largest first; rung 0 must be the station's
            maximal radius.
        station: The station under test, carrying its maximal radius.
        advance: Unit direction of travel into this station, for probe placement.
        regulation: The validated parameters, for the tool radius and the exact
            rational cap surrogate.

    Returns:
        The index into *ladder* of the largest admissible radius, or ``None`` when
        no rung is admissible.
    """
    for rung, radius in enumerate(ladder):
        candidate = replace(station, radius=radius)
        if not _station_is_admissible(stock, candidate, advance, regulation.tool_radius, regulation.cap_ratio):
            continue
        if rung == FULL_RADIUS_RUNG or _loop_reaches_material(stock, candidate, advance, regulation.tool_radius):
            return rung
    return None


def _largest_admissible_radius(
    stock: Stock,
    station: _GuideStation,
    emitted_radius: float,
    advance: Tuple[float, float],
    regulation: _Regulation,
) -> Optional[_RadiusChoice]:
    """Choose this station's loop radius: coarse scan, then refined scan, then ranking.

    THREE STAGES, in the order that keeps the common case cheap.

    1. SCAN THE COARSE LADDER. A cap loose enough that the maximal circle already
       complies costs one rung -- the same single evaluation the advance-only
       generator makes -- and only a station whose maximal circle is refused pays
       for the descent.
    2. RANK THE REFUSED RUNGS. `_least_bad_rung` measures the coarse ladder and
       returns the gentlest circle available there. This runs BEFORE the finer
       scan, not after it, because its measurement is also what decides whether
       the finer scan is worth running: a station whose gentlest circle is a whole
       multiple of the cap over it is not one rung short of complying, it is in a
       regime the cap does not reach.
    3. SCAN THE REFINED LADDER, if stage 2 came back within
       `RADIUS_LADDER_REFINEMENT_MARGIN` of the cap. The coarse spacing is the
       ADVANCE quantisation, and it is too coarse for the radius knob at a station
       whose maximal radius is itself a fraction of the tool radius: one step
       there can cross straight over a band of radii that both comply and cut.
       Where that band exists this stage finds it and the station is not forced at
       all; where it does not, stage 2's answer stands.

    NOTHING HERE SKIPS THE STATION. Leaving a guide station uncut leaves material
    behind, so when no candidate complies a refused circle is emitted, flagged
    forced, and counted.

    WHERE THE DOUBLES ARE. Two of them, both in stage 2 and 3's plumbing and
    neither in a verdict: the ranking orders candidates the exact predicate has
    already refused, and the margin decides how hard to keep looking. Every
    accept / reject on every candidate, coarse or refined, is the same exact
    `cap_exceeded` predicate.

    Args:
        stock: The current stock (unmodified by this call).
        station: The station under test, carrying its maximal radius.
        emitted_radius: Largest radius already emitted at this station.
        advance: Unit direction of travel into this station, for probe placement.
        regulation: The validated parameters, for the tool radius and the exact
            rational cap surrogate.

    Returns:
        The chosen radius and its provenance, or ``None`` when the station has
        nothing left to cut.
    """
    ladder = _radius_ladder(station.radius, emitted_radius, regulation.guide_step, rungs=RADIUS_LADDER_RUNGS, floor_radius=NO_RADIUS_FLOOR)
    if not ladder:
        return None
    rung = _first_admissible_rung(stock, ladder, station, advance, regulation)
    if rung is not None:
        return _RadiusChoice.of(ladder[rung], station, forced=False)
    gentlest = _least_bad_rung(stock, ladder, station, advance, regulation)
    if gentlest.peak <= RADIUS_LADDER_REFINEMENT_MARGIN * regulation.cap_angle:
        refined = _refined_radius_ladder(station.radius, emitted_radius, regulation.guide_step)
        rung = _first_admissible_rung(stock, refined, station, advance, regulation)
        if rung is not None:
            return _RadiusChoice.of(refined[rung], station, forced=False)
    return _RadiusChoice.of(ladder[gentlest.rung], station, forced=True)


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
    forced_loops: List[int] = []
    forced_advances = 0
    first_loop_forced = False

    index = 0
    while True:
        station = stations[index]
        choice: Optional[_RadiusChoice] = None
        if not finished[index]:
            if previous_index is None:
                advance = (station.tx, station.ty)
            else:
                origin = stations[previous_index]
                advance = _unit_tangent(origin.cx, origin.cy, station.cx, station.cy)
            choice = _largest_admissible_radius(stock, station, emitted_radii[index], advance, regulation)

        if choice is not None:
            loop_station = replace(station, radius=choice.radius)
            loop_entry = loop_station.entry
            if previous_entry is None:
                entry = loop_entry
                operations.append(_line_operation(loop_entry, loop_entry, regulation.clearance_z, cut_z, OperationType.PLUNGE, path_index))
                first_loop_forced = choice.forced
            else:
                operations.append(_line_operation(previous_entry, loop_entry, cut_z, cut_z, OperationType.CUT, path_index))
                stock.subtract_capsule_quad(previous_entry[0], previous_entry[1], loop_entry[0], loop_entry[1], regulation.tool_radius)
            if choice.forced:
                forced_loops.append(len(operations))
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
            emitted_radii[index] = choice.radius
            finished[index] = choice.finishes
            forced_radii += int(choice.forced)
            previous_entry = loop_entry
            previous_index = index

        if index == last:
            break
        if choice is not None and choice.finishes:
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
        forced_loops=tuple(forced_loops),
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
            because no candidate radius was admissible.
        forced_loops: Where every one of this chain's forced circles is, the entry
            included, as indices into the WHOLE path's operation stream rather
            than into one pass's.
        forced_advances: Advances taken past a refusing predicate.
        exit: Tool-centre point the chain's last pass retracted from, or ``None``.
    """

    sweeps: int
    entry_forced: bool
    forced_radii: int
    forced_loops: Tuple[int, ...]
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
    forced_loops: List[int] = []
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
        # Rebase the pass's own indices onto the path stream, AFTER the link is in
        # place: the pass numbered its operations from its own plunge.
        base = len(operations)
        forced_loops.extend(base + local for local in outcome.forced_loops)
        operations.extend(outcome.operations)
        last_exit = outcome.exit
        sweeps += 1
    else:
        raise RadialSweepBudgetExceededError(
            f"Skeleton chain {path_index} still had material to cut after {MAX_RADIAL_SWEEPS_PER_CHAIN} sweeps. "
            "Each sweep must lift some station at least one rung of a "
            f"{RADIUS_LADDER_RUNGS}-rung ladder, so this means a station is re-climbing rungs it already emitted."
        )

    return _ChainOutcome(
        sweeps=sweeps,
        entry_forced=entry_forced,
        forced_radii=forced_radii,
        forced_loops=tuple(forced_loops),
        forced_advances=forced_advances,
        exit=last_exit,
    )


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

    Where no radius on that ladder complies, the ladder is rescanned at
    `RADIUS_LADDER_SUBDIVISIONS` times the resolution, and if that finds nothing
    either, the circle emitted is the gentlest of the refused candidates by
    MEASURED engagement -- reported, counted, and warned about, never presented as
    meeting the cap. Cutting back where nothing complies does not finish the
    station, so it costs an extra circle and an extra sweep; that trade is
    balanced by `RADIUS_LADDER_REFINEMENT_MARGIN` and measured in its comment.

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
        RadialToolpathResult: The typed operation stream -- `PLUNGE`, `CUT`
        machining circles and bridge lines, `RETRACT`, and clearance-height `LINK`
        moves between passes -- plus the tessellated visualisation polyline, and
        `forced_loops`, the indices of the circles no candidate radius could bring
        under the cap. A `ToolpathResult` in every other respect, so
        `audit_toolpath_engagement` and every downstream consumer work unchanged.
        Note that a chain contributes one plunge/retract pair PER PASS, not one in
        total.

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
    forced_loops: List[int] = []
    forced_advances = 0
    sweeps = 0
    for path_index, stations in enumerate(chains):
        outcome = _machine_chain_radially(stock, stations, path_index, regulation, cut_z, last_exit, operations)
        entry_slots += int(outcome.entry_forced)
        forced_radii += outcome.forced_radii
        forced_loops.extend(outcome.forced_loops)
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

    return RadialToolpathResult(
        operations=operations,
        polyline=_tessellate(operations, samples_per_radian),
        forced_loops=frozenset(forced_loops),
    )
