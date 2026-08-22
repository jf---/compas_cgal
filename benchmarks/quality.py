"""Machining quality in four groups: is the path VALID, GOOD, FAST, and KIND to the tool.

Nothing in this repository asked whether a generated path was worth running. The
certifier answers "is this motion inside the cap?", and a path that plunges into
a corner and cuts straight back out of it answers that question perfectly well.
This module asks the other question, and it separates four kinds of answer that
should never be averaged together:

* `ElementaryQuality` -- VALIDITY. A defect here makes the path wrong, not merely
  poor: it gouges, it leaves material, it emits motions that remove nothing.
  Binary or exact wherever the exact kernel can decide it.
* `CutQuality` -- CUT MECHANICS. What the cutting edge experiences.
* `SpeedQuality` -- what the path costs a machine that has to execute it.
* `LongevityQuality` -- what it costs the tool.

TWO PROXIES THIS PROJECT RESTS ON, AND WHAT THEY ARE PROXIES FOR.

ENGAGEMENT ANGLE IS A PROXY; CHIP THICKNESS IS THE LOAD. The chip is
`h(phi) = f_z sin(phi)` and the engaged arc reaches `theta`, so the maximum
undeformed chip thickness is `h_ex = f_z sin(min(theta, 90 deg))`. It rises to a
PLATEAU at 90 degrees and never falls again: tightening a cap from 172 to 130
degrees changes `h_ex` NOT AT ALL, while below 90 degrees the same feed produces
a thinner and thinner chip until, under a material-dependent `h_min`, the edge
ploughs and RUBS instead of cutting -- which wears a tool faster than a heavier
cut does. That is why this module reports a LOWER bound on the load
(`low_chip_thickness_ratio`) as well as an upper one; nothing else in the suite
expresses that a cut can be too light.

PATH LENGTH IS A PROXY; FEED-LIMITED TIME IS THE COST. Feed through a curve is
bounded by `v <= sqrt(a_max / kappa)`, so a shorter path with sharper curvature
can take LONGER. `SpeedQuality` reports curvature and its discontinuities from
geometry alone; the time itself needs a `MachineModel` and is not offered
without one.

GEOMETRY-DERIVED VERSUS MODEL-DERIVED. Every field of the four groups is
computable from the toolpath plus the exact stock -- ours, reproducible here.
Everything needing a workpiece material or a machine's dynamics lives in
`material_outcome`, `machine_outcome`, and `tool_life_outcome`, each of which
takes its model as an argument and raises a NAMED error when it is absent. A
guessed coefficient is never substituted, because a number that depends on a
calibration is only as good as the calibration.

EXACT VERSUS SAMPLED, within the geometry-derived half. Effects on the stock are
exact set operations (`Stock.exactly_equals`) and radius comparisons are exact
arithmetic. Everything read ALONG a motion is sampled at
`QUALITY_SAMPLES_PER_MOTION` positions and everything read over an AREA is
sampled on a grid; each individual query is an exact predicate, the aggregation
is not. A sampled field can demonstrate a defect and can never certify its
absence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from benchmarks.coverage import COVERAGE_GRID_SAMPLES
from benchmarks.coverage import measure_coverage
from benchmarks.errors import MissingMachineModelError
from benchmarks.errors import MissingMaterialModelError
from benchmarks.errors import ZeroLengthToolpathError
from benchmarks.models import MachineModel
from benchmarks.models import MaterialModel
from benchmarks.pathmetrics import entry_cut_indices
from benchmarks.spec import PocketSpec
from benchmarks.survey import QUALITY_SAMPLES_PER_MOTION
from benchmarks.survey import MotionKind
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.survey import RapidMotion
from benchmarks.survey import survey_path
from compas_cgal.toolpath import ToolpathResult

# A machining loop must sweep an ANNULUS. A tool of radius r running a circle of
# radius rho covers the radii [rho - r, rho + r] about the loop centre, so the
# swept region has an uncut core exactly when rho > r. At rho <= r the tool
# passes over its own loop centre: the sweep is a filled disk, the cutter never
# leaves the material at the middle of the loop, and the "trochoid" is a bore
# wearing a circle's name. THIS IS PHYSICS, not a tuned constant -- it is the
# qualitative boundary at which the annulus loses its hole, and the test is an
# exact comparison of two doubles.
DEGENERATE_LOOP_RATIO = 1.0

# Above the degeneracy boundary a loop has a hole, but a hole a fraction of a
# millimetre across buys no trochoidal relief worth the arc length. A loop is
# called MARGINAL up to this multiple of the tool radius. ENGINEERING JUDGEMENT
# and nothing more: at 1.5 r the uncut core is half a tool radius across, which
# is the smallest core that still lets the cutter unload measurably between
# passes. Deliberately kept a separate field from `degenerate_loops` so a reader
# can discard the judgement and keep the physics.
MARGINAL_LOOP_RATIO = 1.5

# Chip thickness plateaus once the engaged arc passes a quarter turn: the chip is
# h(phi) = f_z sin(phi) and phi reaches 90 degrees, where sin is 1. Past that,
# more engagement adds engaged arc at a thickness already at maximum.
CHIP_PLATEAU_DEG = 90.0

# A full turn of engaged rim: the cutter is surrounded, which is a plunge and
# the most immersed it can be. Named so `textbook_engagement_deg` can clamp to
# it rather than carrying a bare 360.0.
FULL_TURN_DEG = 360.0

# Engagement bands for the time-at-engagement histogram, in degrees. Six equal
# 30-degree bands span the cap's legal range of (0, 180], and a seventh catches
# the entry cuts, which reach a full turn for any generator entering solid stock
# and would otherwise silently fall off the top of the histogram.
ENGAGEMENT_BANDS_DEG: Tuple[float, ...] = (0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 360.0)

# Percentile reported beside the maximum. The maximum is one position; the 95th
# says whether the load is a spike or the path's normal state.
ENGAGEMENT_PERCENTILE = 95.0

# Percentile that stands for the LIGHT end of the load distribution. The raw
# minimum is worthless here: every path contains positions where the rim grazes
# material tangentially and removes nothing, so the minimum chip thickness reads
# ~1e-8 on all of them. The 5th percentile of the engaged cut length is the
# lightest load the tool actually spends time at.
LOW_CHIP_PERCENTILE = 5.0

# Coincidence tolerance for endpoints and tangent directions, as a fraction of
# the tool radius and as a dot-product slack. Continuity is a question about
# doubles a generator computed, not an exact-kernel predicate: two endpoints the
# generator intended to coincide can differ in the last bits. A thousandth of the
# tool radius is far below any real discontinuity and far above rounding.
CONTINUITY_TOOL_RADIUS_FRACTION = 1e-3

# Tangents are called continuous when their dot product is within this of one --
# about 0.8 degrees. Below that a junction is a corner the machine must
# decelerate through.
TANGENT_CONTINUITY_SLACK = 1e-4

# Half-width of the steady band, as a fraction of the DESIGN radial immersion
# the engagement cap implies. Chatter stability is a function of radial
# immersion a_e/D: the stability lobes for a given tool and spindle shift as
# immersion changes, so a path whose immersion wanders crosses between lobe
# regions and can go unstable even where every instantaneous engagement sits
# under the cap.
#
# The band is deliberately scaled by the CAP and not by the path's own median.
# A median-relative band degenerates: on a path that spends most of its length
# barely touching material the median immersion is near zero, a +/-10% band
# around it is near zero wide, and the metric reports a catastrophic-looking
# number that measures the median's smallness rather than the cut's steadiness.
# Scaling by the commissioned load keeps one fixed absolute width for every path
# measured at the same cap, which is also what makes two generators comparable.
IMMERSION_STEADY_BAND_FRACTION = 0.10

# Percentile of the CUT-BLOCK-LENGTH distribution reported as its short tail.
# A control reads one block at a time, so the short tail is what starves its
# look-ahead; a mean hides it completely.
SHORT_BLOCK_PERCENTILE = 5.0

# The middle of a distribution. Named so the percentile helper's call sites read
# as intent rather than as a literal.
MEDIAN_PERCENTILE = 50.0


@dataclass(frozen=True)
class ElementaryQuality:
    """Validity: defects that make a path WRONG rather than merely poor.

    Attributes:
        uncut_fraction: Share of the tool-reachable material left behind.
            SAMPLED on a grid.
        gouge_free: Whether no sampled cutter centre left the legal centre
            domain. SAMPLED: a violation is certain, an absence is not.
        gouging_motions: Cut motions with at least one sampled centre outside.
        rapid_safety: Whether no rapid travels in XY on the cutting plane.
        unsafe_rapids: Rapids that do.
        continuity_breaks: Junctions where one motion does not end where the
            next begins.
        zero_length_motions: Operations of no length at all.
        degenerate_loops: Loops with ``rho <= r`` -- a disk sweep, not a
            trochoid. EXACT, and physics rather than judgement.
        marginal_loops: Loops with ``r < rho <= MARGINAL_LOOP_RATIO * r`` -- a
            hole too small to be useful. ENGINEERING JUDGEMENT.
        redundant_operations: Cut motions that leave the stock exactly
            unchanged. EXACT.
        recut_fraction: Share of the swept area that met no fresh material --
            re-cut chips and already-cleared space together. Held discusses chip
            re-cutting explicitly and this is its geometric shadow, but READ THE
            NUMBER CAREFULLY: a trochoidal path scores near 0.9 BY CONSTRUCTION,
            because consecutive loops are meant to overlap and most of a loop's
            sweep passes through space the previous loop cleared. It is
            REPORTED, never gated, and a change in it between two paths on the
            same pocket is the only reading of it that means anything.
    """

    uncut_fraction: float
    gouge_free: bool
    gouging_motions: int
    rapid_safety: bool
    unsafe_rapids: int
    continuity_breaks: int
    zero_length_motions: int
    degenerate_loops: int
    marginal_loops: int
    redundant_operations: int
    recut_fraction: float


@dataclass(frozen=True)
class CutQuality:
    """Cut mechanics: what the cutting edge experiences.

    Every engagement statistic is LENGTH-WEIGHTED -- a long loop at 90 degrees
    counts for more than a short bridge at 90 degrees, because the edge spends
    more of its life there. Sample-counting would weight a 0.1 mm link exactly
    as heavily as a 7 mm loop.

    Attributes:
        max_engagement_deg: Largest sampled engaged run. SAMPLED lower bound.
        cap_exceedances: Cut motions where the exact cap predicate fired at some
            sampled position. Includes the entry cut after each plunge, which is
            a full slot for any generator entering solid stock without a helical
            or pre-drilled entry. SAMPLED lower bound on true exceedance.
        engagement_p95_deg: `ENGAGEMENT_PERCENTILE` of the length-weighted
            engagement distribution.
        engagement_variance_deg2: Length-weighted variance of engagement over
            the cutting motions -- steadiness of load.
        max_chip_thickness_ratio: Largest ``h_ex / f_z`` over the path, which is
            ``sin(min(theta, 90 deg))``. Dimensionless and geometry-derived; the
            length needs a `MaterialModel`.
        low_chip_thickness_ratio: The 5th percentile of ``h_ex / f_z`` over the
            ENGAGED cut length. The LOWER bound the cap alone cannot express,
            and the one that decides whether the edge cuts or rubs. A
            percentile rather than the raw minimum on purpose: the minimum is
            always a position where the rim merely grazes material and removes
            nothing, so it reads 1e-8 on every path and says nothing about any
            of them.
        max_engagement_gradient_deg_per_length: Steepest change of engagement
            per unit travel WITHIN a motion. Load shock, which damages an edge
            more than a high steady load does.
        max_engagement_step_deg: Largest change between two motions the tool
            runs back to back without leaving the cut. The same shock across a
            junction, where the gradient is undefined because the travel is zero.
        slotting_motions: Straight cut motions engaging past
            `benchmarks.survey.SLOT_ENGAGEMENT_FRACTION` of the cap.
        immersion_steady_fraction: Share of cut length whose radial immersion
            lies within `IMMERSION_STEADY_BAND_FRACTION` of the DESIGN immersion
            either side of the length-weighted median. One means a single-lobe
            cut; low means the path sweeps across stability lobes.
        immersion_at_design_fraction: Share of cut length cutting within the
            same band of the DESIGN immersion the cap implies. Steadiness at a
            tenth of the intended load is stable and slow; this separates them.
        immersion_excursions: Times the immersion leaves that band along the
            path. Variance measures spread; this measures how OFTEN the cut is
            disturbed, which is what a spindle feels.
        mean_radial_depth: Length-weighted mean radial depth of cut,
            ``a_e = r (1 - cos(theta/2))``. Removal rate per unit travel, which
            in a fixed-depth plane is what material removal rate reduces to.
        radial_depth_variance: Its length-weighted variance.
        wall_scallop_height: Deepest residue against a wall. SAMPLED on a grid.
        loop_radius_cv: Coefficient of variation of the guide radius over the
            emitted loops -- the standard deviation over the mean. A LARGE VALUE
            INDICATES THE CENTRE LOCUS, NOT THE REGULATOR: Held's Figure 5 puts
            the machining-circle centres on a "middle" curve BETWEEN the medial
            axis and the boundary, where clearance varies smoothly, while
            placing them on skeleton stations puts them where clearance is
            maximal and collapses abruptly at branches and corners.
        max_loop_radius_step: Largest ``|rho_i - rho_(i+1)|`` between two loops
            the tool runs BACK TO BACK, in tool radii. The same defect measured
            locally instead of globally, and the sharper of the two: a single
            jump is invisible in a coefficient of variation. Scoped to a
            continuous run, never across a retract -- a change in guide radius
            the cutter performs in the air is not a load step.
    """

    max_engagement_deg: float
    cap_exceedances: int
    engagement_p95_deg: float
    engagement_variance_deg2: float
    max_chip_thickness_ratio: float
    low_chip_thickness_ratio: float
    max_engagement_gradient_deg_per_length: float
    max_engagement_step_deg: float
    slotting_motions: int
    immersion_steady_fraction: float
    immersion_at_design_fraction: float
    immersion_excursions: int
    mean_radial_depth: float
    radial_depth_variance: float
    wall_scallop_height: float
    loop_radius_cv: float
    max_loop_radius_step: float


@dataclass(frozen=True)
class SpeedQuality:
    """What the path costs a machine that has to execute it.

    Attributes:
        cutting_length: Length of the cut-plane motions.
        air_length: Length of the motions that remove nothing.
        air_fraction: Air over total.
        max_curvature: Largest path curvature, in reciprocal length. The feed
            ceiling ``sqrt(a_max / kappa)`` is tightest here.
        tangent_breaks: Junctions where the travel direction jumps -- G1
            discontinuities, at which a machine must come off feed.
        curvature_breaks: Junctions where curvature jumps but the tangent does
            not -- G2 discontinuities, which cost acceleration rather than a
            full stop.
        direction_reversals: Junctions where the tool turns back on itself.
        retract_count: Upward clearance moves.
        reentry_count: Downward plunges, each of which is a re-entry.
    """

    cutting_length: float
    air_length: float
    air_fraction: float
    max_curvature: float
    tangent_breaks: int
    curvature_breaks: int
    direction_reversals: int
    retract_count: int
    reentry_count: int


@dataclass(frozen=True)
class LongevityQuality:
    """What the path costs the tool.

    Attributes:
        material_entries: Cut motions that are the first contact after a plunge.
            Each is an impact, and impacts are the leading cause of edge
            chipping.
        engagement_length_histogram: Cut length spent in each band of
            `ENGAGEMENT_BANDS_DEG`, as ``(low_deg, high_deg, length)``.
            Cumulative damage, which a peak alone does not describe.
        cut_air_alternations: Transitions between cutting and not cutting. Each
            is a thermal cycle for a coated edge.
        alternations_per_length: The same, per unit of total path.
    """

    material_entries: int
    engagement_length_histogram: Tuple[Tuple[float, float, float], ...]
    cut_air_alternations: int
    alternations_per_length: float


@dataclass(frozen=True)
class ProgramQuality:
    """Whether a control can execute the path at the speed it was programmed for.

    Every field is a property of the PATH ALONE -- no material, no machine, no
    calibration -- so it is reproducible from this repository and defensible
    without declaring an assumption. What it deliberately does NOT say is
    whether a PARTICULAR control starves: that verdict needs a feed rate and an
    interpolation period, which are machine properties and live in
    `MachineModel`.

    The failure this group exists to catch is DATA STARVATION. A control reads
    the program one block at a time; when blocks are shorter than it can consume
    at the commanded feed, look-ahead empties and the machine slows below the
    programmed rate no matter what the geometry allows. A trochoidal path is
    exactly the shape that provokes it -- hundreds of short arcs where a contour
    path would have a handful of long ones -- so a trochoidal generator that is
    never measured here is shipping an unmeasured risk.

    Attributes:
        block_count: Cut motions plus rapids -- every block the control reads.
        cut_blocks: Cut motions alone.
        min_block_length: Shortest cut motion.
        median_block_length: The middle one. Reported beside the mean's absence
            deliberately: block-length distributions here are heavily skewed and
            a mean describes neither tail.
        short_block_length: The `SHORT_BLOCK_PERCENTILE` tail -- the length
            below which that share of BLOCKS falls, counted per block rather
            than weighted by length, because starvation is a per-block event.
        block_length_cv: Coefficient of variation of cut-block length. A uniform
            program is one a control can pipeline; a mixed one is not.
        blocks_per_unit_length: Cut blocks per unit of cut length. The direct
            measure of program density, comparable across pocket sizes.
        arc_length_fraction: Share of cut length carried by arcs and loops
            rather than straight blocks. Arcs let look-ahead breathe: one G2
            block spans what a chorded polyline would spend dozens on.
    """

    block_count: int
    cut_blocks: int
    min_block_length: float
    median_block_length: float
    short_block_length: float
    block_length_cv: float
    blocks_per_unit_length: float
    arc_length_fraction: float


@dataclass(frozen=True)
class PathQuality:
    """The four groups, plus what they were measured on.

    Attributes:
        elementary: Validity.
        cut: Cut mechanics.
        speed: Machine cost.
        longevity: Tool cost.
        program: Whether a control can execute it at the programmed feed.
        cut_operations: Cut-plane motions in the path.
        path_length: Analytic length of every operation.
    """

    elementary: ElementaryQuality
    cut: CutQuality
    speed: SpeedQuality
    longevity: LongevityQuality
    program: ProgramQuality
    cut_operations: int
    path_length: float


@dataclass(frozen=True)
class MaterialOutcome:
    """Cut mechanics as lengths, which needs a feed and an edge radius.

    Attributes:
        max_chip_thickness_mm: ``h_ex`` at its worst.
        low_chip_thickness_mm: ``h_ex`` at the 5th percentile of the engaged
            cut length -- the light end, robust to grazing positions.
        rubbing_length: Cut length where ``h_ex`` falls under ``h_min`` while
            the tool is still engaged -- ploughing, not cutting.
        rubbing_fraction: That length over the engaged cutting length.
    """

    max_chip_thickness_mm: float
    low_chip_thickness_mm: float
    rubbing_length: float
    rubbing_fraction: float


@dataclass(frozen=True)
class MachineOutcome:
    """Time, which needs feed rates and an acceleration limit.

    Attributes:
        cutting_seconds: Time on the cut-plane motions, each segment held to
            ``min(feed, sqrt(a_max / kappa))``.
        air_seconds: Time on the rapids, at the traverse rate.
        total_seconds: The two together.
        mean_cutting_feed_mm_per_s: Cutting length over cutting time -- the feed
            the machine actually averages, against the one programmed.
        feed_utilisation: That mean over the programmed feed, in ``[0, 1]``.
    """

    cutting_seconds: float
    air_seconds: float
    total_seconds: float
    mean_cutting_feed_mm_per_s: float
    feed_utilisation: float


@dataclass(frozen=True)
class ToolLifeOutcome:
    """A Taylor life estimate, which needs both models.

    Attributes:
        cutting_speed_m_per_min: Surface speed ``V`` from the spindle and tool.
        taylor_life_minutes: ``T`` from ``V * T^n = C``.
        cutting_minutes: Time the tool spends in the cut on this path.
        life_fraction_consumed: Cutting time over Taylor life -- how much of one
            edge this single pocket costs.
    """

    cutting_speed_m_per_min: float
    taylor_life_minutes: float
    cutting_minutes: float
    life_fraction_consumed: float


def measure_quality(
    spec: PocketSpec,
    result: ToolpathResult,
    *,
    samples_per_motion: int = QUALITY_SAMPLES_PER_MOTION,
    grid: int = COVERAGE_GRID_SAMPLES,
) -> PathQuality:
    """Measure *result* as a way of machining *spec*, in five groups.

    Args:
        spec: The pocket, tool, and engagement cap the path was generated for.
        result: The generated toolpath.
        samples_per_motion: Cutter positions probed along each cut motion.
        grid: Coverage samples along the pocket's longer bounding-box side.

    Returns:
        The four groups.

    Raises:
        InvalidMotionSampleCountError: *samples_per_motion* is below one.
        InvalidGridResolutionError: *grid* is below one.
        CoarseCoverageGridError: *grid* cannot resolve residue at this tool size.
        EmptyReachableRegionError: No grid sample landed in the reachable region.
        ZeroLengthToolpathError: The path has no length, so no fraction exists.
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        UnsampleableMotionError: A cut motion carries no cutter-centre path.
        UnmeasurableOperationLengthError: An operation's length is undefined.
    """
    survey = survey_path(spec, result, samples_per_motion=samples_per_motion)
    if survey.total_length <= 0.0:
        raise ZeroLengthToolpathError(f"{spec.name}: the toolpath sums to zero length, so no length fraction is defined.")
    coverage = measure_coverage(spec, survey.final_stock, grid=grid)
    entries = entry_cut_indices(result)
    return PathQuality(
        elementary=_elementary(spec, survey, coverage.uncut_fraction, coverage.remaining_area),
        cut=_cut(spec, survey, coverage.wall_scallop_height),
        speed=_speed(survey),
        longevity=_longevity(survey, len(entries)),
        program=_program(survey),
        cut_operations=len(survey.motions),
        path_length=survey.total_length,
    )


def _elementary(spec: PocketSpec, survey: PathSurvey, uncut_fraction: float, remaining_area: float) -> ElementaryQuality:
    """Reduce a survey to the validity group.

    `recut_fraction` compares the area the tool SWEPT, which is closed-form
    arithmetic over the primitives, against the area it actually REMOVED, which
    is the pocket's area less what the grid still finds standing. Sweeping more
    than you remove is either re-cutting your own chips or cutting air.

    Args:
        spec: The instance, for the pocket area.
        survey: The replay's findings.
        uncut_fraction: From the coverage grid.
        remaining_area: Estimated area still standing, from the coverage grid.

    Returns:
        The validity group.
    """
    loops = _loop_radii(survey.motions)
    swept = survey.swept_area
    removed = max(0.0, abs(spec.polygon.area) - remaining_area)
    return ElementaryQuality(
        uncut_fraction=uncut_fraction,
        gouge_free=not any(motion.gouges for motion in survey.motions),
        gouging_motions=sum(1 for motion in survey.motions if motion.gouges),
        rapid_safety=not any(rapid.horizontal_at_cut_plane for rapid in survey.rapids),
        unsafe_rapids=sum(1 for rapid in survey.rapids if rapid.horizontal_at_cut_plane),
        continuity_breaks=_continuity_breaks(survey, spec.tool_radius),
        zero_length_motions=sum(1 for motion in survey.motions if motion.length == 0.0) + sum(1 for rapid in survey.rapids if rapid.length == 0.0),
        degenerate_loops=sum(1 for radius in loops if radius <= DEGENERATE_LOOP_RATIO * spec.tool_radius),
        marginal_loops=sum(1 for radius in loops if DEGENERATE_LOOP_RATIO * spec.tool_radius < radius <= MARGINAL_LOOP_RATIO * spec.tool_radius),
        redundant_operations=sum(1 for motion in survey.motions if not motion.removes_material),
        recut_fraction=0.0 if swept <= 0.0 else max(0.0, 1.0 - removed / swept),
    )


def _cut(spec: PocketSpec, survey: PathSurvey, wall_scallop_height: float) -> CutQuality:
    """Reduce a survey to the cut-mechanics group.

    Args:
        spec: The instance, for the tool radius.
        survey: The replay's findings.
        wall_scallop_height: From the coverage grid.

    Returns:
        The cut-mechanics group.
    """
    weighted = _weighted_engagement(survey.motions)
    radii = _loop_radii(survey.motions)
    runs = _loop_runs(survey.motions, survey.rapids)
    engaged = [(value, weight) for value, weight in weighted if value > 0.0]
    immersions = [(radial_immersion(value), weight) for value, weight in weighted]
    depths = [(spec.tool_diameter * value, weight) for value, weight in immersions]
    design = design_immersion(spec.tea_cap_deg)
    band = IMMERSION_STEADY_BAND_FRACTION * design
    return CutQuality(
        max_engagement_deg=max((value for value, _ in weighted), default=0.0),
        cap_exceedances=sum(1 for motion in survey.motions if motion.cap_exceeded),
        engagement_p95_deg=_weighted_percentile(weighted, ENGAGEMENT_PERCENTILE),
        engagement_variance_deg2=_weighted_variance(weighted),
        max_chip_thickness_ratio=max((chip_thickness_ratio_from_rim(value) for value, _ in weighted), default=0.0),
        low_chip_thickness_ratio=_weighted_percentile([(chip_thickness_ratio_from_rim(value), weight) for value, weight in engaged], LOW_CHIP_PERCENTILE),
        max_engagement_gradient_deg_per_length=_max_engagement_gradient(survey.motions),
        max_engagement_step_deg=_max_engagement_step(survey.motions),
        slotting_motions=sum(1 for motion in survey.motions if motion.slot_exceeded),
        immersion_steady_fraction=_immersion_steady_fraction(immersions, band),
        immersion_at_design_fraction=_immersion_at_design_fraction(immersions, design, band),
        immersion_excursions=_immersion_excursions(immersions, band),
        mean_radial_depth=_weighted_mean(depths),
        radial_depth_variance=_weighted_variance(depths),
        wall_scallop_height=wall_scallop_height,
        loop_radius_cv=_loop_radius_cv(radii),
        max_loop_radius_step=_max_loop_radius_step(runs, spec.tool_radius),
    )


def _speed(survey: PathSurvey) -> SpeedQuality:
    """Reduce a survey to the machine-cost group.

    Args:
        survey: The replay's findings.

    Returns:
        The machine-cost group.
    """
    tangent_breaks, curvature_breaks, reversals = _junction_breaks(survey.motions)
    return SpeedQuality(
        cutting_length=survey.cut_length,
        air_length=survey.air_length,
        air_fraction=survey.air_length / survey.total_length,
        max_curvature=max((motion.curvature for motion in survey.motions), default=0.0),
        tangent_breaks=tangent_breaks,
        curvature_breaks=curvature_breaks,
        direction_reversals=reversals,
        retract_count=survey.retracts,
        reentry_count=survey.plunges,
    )


def _longevity(survey: PathSurvey, material_entries: int) -> LongevityQuality:
    """Reduce a survey to the tool-cost group.

    Args:
        survey: The replay's findings.
        material_entries: Cut motions that are the first contact after a plunge.

    Returns:
        The tool-cost group.
    """
    alternations = _cut_air_alternations(survey)
    return LongevityQuality(
        material_entries=material_entries,
        engagement_length_histogram=_engagement_histogram(survey.motions),
        cut_air_alternations=alternations,
        alternations_per_length=alternations / survey.total_length,
    )


def chip_thickness_ratio(engagement_deg: float) -> float:
    """Maximum undeformed chip thickness over the feed per tooth, ``h_ex / f_z``.

    The chip at angular position phi is ``f_z sin(phi)`` and the engaged arc
    reaches *engagement_deg*, so the maximum is ``sin(min(theta, 90 deg))``. It
    rises as a sine to a PLATEAU at a quarter turn and NEVER FALLS: past 90
    degrees more engagement adds engaged arc at a thickness already at its
    maximum. Writing ``min(sin(theta), 1)`` instead would make the curve fall
    again past 90 degrees and tell a reader that heavier engagement means a
    lighter chip.

    Args:
        engagement_deg: Engaged arc of the cutter rim, in degrees.

    Returns:
        The ratio in ``[0, 1]``.
    """
    return math.sin(math.radians(min(max(engagement_deg, 0.0), CHIP_PLATEAU_DEG)))


@dataclass(frozen=True)
class CalibratedOutcome:
    """The three outcomes that only exist once a calibration is declared.

    Kept OUT of `PathQuality` on purpose. Everything in `PathQuality` is derived
    from the toolpath and the exact stock, so it is reproducible from this
    repository alone. Everything here additionally depends on a workpiece
    material or a machine's dynamics, and is therefore only as good as the
    `MaterialModel` and `MachineModel` it was measured against. Merging the two
    would make a calibrated estimate indistinguishable from a computed fact in
    every table that prints them.

    Any figure or report built on these MUST name the models used.

    Attributes:
        material: Chip thickness as a length, and the rubbing share.
        machine: Feed-limited cycle time. A LOWER bound -- see `machine_outcome`.
        tool_life: A Taylor life estimate and the fraction this path consumes.
    """

    material: MaterialOutcome
    machine: MachineOutcome
    tool_life: ToolLifeOutcome


def measure_outcomes(survey: PathSurvey, *, material: MaterialModel, machine: MachineModel) -> CalibratedOutcome:
    """Measure the calibration-dependent outcomes for an already-surveyed path.

    Takes a `PathSurvey` rather than a spec and a result so that a caller
    measuring both `PathQuality` and `CalibratedOutcome` replays the path once.

    Args:
        survey: The replay's findings, from `survey_path`.
        material: The workpiece and cutting-edge coefficients. Required.
        machine: The feed, rapid, acceleration and spindle limits. Required.

    Returns:
        The three calibrated outcomes.

    Raises:
        MissingMaterialModelError: *material* is None.
        MissingMachineModelError: *machine* is None.
    """
    return CalibratedOutcome(
        material=material_outcome(survey, material),
        machine=machine_outcome(survey, machine),
        tool_life=tool_life_outcome(survey, material, machine),
    )


def material_outcome(survey: PathSurvey, material: Optional[MaterialModel]) -> MaterialOutcome:
    """Chip thickness as a length, and how much of the cut is rubbing.

    Args:
        survey: The replay's findings.
        material: The coefficients. Required.

    Returns:
        The material-derived outcome.

    Raises:
        MissingMaterialModelError: *material* is None. No default is substituted:
            `h_min` decides whether the edge cuts or ploughs, and a guessed one
            would make a fabricated verdict look like a measurement.
    """
    if material is None:
        raise MissingMaterialModelError("max_chip_thickness_mm and rubbing_fraction need a MaterialModel; pass one rather than accepting a guessed feed and edge radius.")
    weighted = _weighted_engagement(survey.motions)
    engaged = [(value, weight) for value, weight in weighted if value > 0.0]
    engaged_length = sum(weight for _, weight in engaged)
    floor = material.min_chip_thickness_mm
    rubbing = sum(weight for value, weight in engaged if material.feed_per_tooth_mm * chip_thickness_ratio_from_rim(value) < floor)
    return MaterialOutcome(
        max_chip_thickness_mm=material.feed_per_tooth_mm * max((chip_thickness_ratio_from_rim(value) for value, _ in weighted), default=0.0),
        low_chip_thickness_mm=material.feed_per_tooth_mm
        * _weighted_percentile([(chip_thickness_ratio_from_rim(value), weight) for value, weight in engaged], LOW_CHIP_PERCENTILE),
        rubbing_length=rubbing,
        rubbing_fraction=0.0 if engaged_length <= 0.0 else rubbing / engaged_length,
    )


def machine_outcome(survey: PathSurvey, machine: Optional[MachineModel]) -> MachineOutcome:
    """Feed-limited time, integrating the curvature ceiling along the path.

    Each cut motion is held to ``min(feed, sqrt(a_max / kappa))`` for its whole
    length, which is exact for a circle (constant curvature) and exact for a
    straight move (no bound). It IGNORES the acceleration ramps into and out of
    each motion and the jerk cost of the tangent breaks `SpeedQuality` counts, so
    it is a LOWER BOUND on cycle time -- an optimistic one, and the honest
    direction for a claim about a path being fast.

    Args:
        survey: The replay's findings.
        machine: The limits. Required.

    Returns:
        The machine-derived outcome.

    Raises:
        MissingMachineModelError: *machine* is None.
    """
    if machine is None:
        raise MissingMachineModelError("Cycle time needs a MachineModel; feed, rapid rate and the acceleration limit are properties of a machine, not of the path.")
    cutting = sum(motion.length / machine.curvature_limited_feed_mm_per_s(motion.curvature) for motion in survey.motions if motion.length > 0.0)
    air = survey.air_length / machine.rapid_mm_per_s
    mean_feed = 0.0 if cutting <= 0.0 else survey.cut_length / cutting
    return MachineOutcome(
        cutting_seconds=cutting,
        air_seconds=air,
        total_seconds=cutting + air,
        mean_cutting_feed_mm_per_s=mean_feed,
        feed_utilisation=0.0 if machine.feed_mm_per_s <= 0.0 else mean_feed / machine.feed_mm_per_s,
    )


def tool_life_outcome(survey: PathSurvey, material: Optional[MaterialModel], machine: Optional[MachineModel]) -> ToolLifeOutcome:
    """A Taylor tool-life estimate, ``V * T^n = C``.

    Args:
        survey: The replay's findings.
        material: The Taylor coefficients. Required.
        machine: The spindle and feed limits. Required.

    Returns:
        The life estimate.

    Raises:
        MissingMaterialModelError: *material* is None.
        MissingMachineModelError: *machine* is None.
    """
    if material is None:
        raise MissingMaterialModelError("A Taylor life estimate needs a MaterialModel for the exponent n and the constant C.")
    if machine is None:
        raise MissingMachineModelError("A Taylor life estimate needs a MachineModel for the spindle speed that sets the cutting speed V.")
    speed = math.pi * survey.spec.tool_diameter * machine.spindle_rpm / 1000.0
    life = (material.taylor_constant_m_per_min / speed) ** (1.0 / material.taylor_exponent)
    minutes = machine_outcome(survey, machine).cutting_seconds / 60.0
    return ToolLifeOutcome(
        cutting_speed_m_per_min=speed,
        taylor_life_minutes=life,
        cutting_minutes=minutes,
        life_fraction_consumed=0.0 if life <= 0.0 else minutes / life,
    )


def _program(survey: PathSurvey) -> ProgramQuality:
    """Reduce a survey to the program-feasibility group.

    Block lengths are counted PER BLOCK, not weighted by length. Starvation is a
    per-block event -- the control stalls once per short block regardless of how
    little distance that block covers -- so length-weighting would hide exactly
    the population that causes it.

    Args:
        survey: The replay's findings.

    Returns:
        The program-feasibility group.
    """
    lengths = [motion.length for motion in survey.motions if motion.length > 0.0]
    counted = [(length, 1.0) for length in lengths]
    mean = _weighted_mean(counted)
    arc_length = sum(motion.length for motion in survey.motions if motion.kind in (MotionKind.LOOP, MotionKind.ARC))
    return ProgramQuality(
        block_count=len(survey.motions) + len(survey.rapids),
        cut_blocks=len(survey.motions),
        min_block_length=min(lengths, default=0.0),
        median_block_length=_weighted_percentile(counted, MEDIAN_PERCENTILE),
        short_block_length=_weighted_percentile(counted, SHORT_BLOCK_PERCENTILE),
        block_length_cv=0.0 if mean <= 0.0 else math.sqrt(_weighted_variance(counted)) / mean,
        blocks_per_unit_length=0.0 if survey.cut_length <= 0.0 else len(survey.motions) / survey.cut_length,
        arc_length_fraction=0.0 if survey.cut_length <= 0.0 else arc_length / survey.cut_length,
    )


def textbook_engagement_deg(rim_arc_deg: float) -> float:
    """The textbook engagement angle for a rim arc this kernel reported.

    THE TWO ANGLES ARE NOT THE SAME AND DIFFER BY A FACTOR OF TWO. The exact
    kernel's `max_run_tea` is the ENGAGED ARC OF THE CUTTER RIM, which reaches a
    full turn when the tool is surrounded -- which is why a plunge measures 360
    degrees here. The textbook milling angle is measured from entry to exit and
    tops out at 180 degrees for a full slot. Writing one where the other belongs
    is the single easiest way to get every chip-load number on this page wrong.

    DERIVED, THEN MEASURED. For a cut into a straight wall the rim arc satisfies
    ``a_e = r (1 - cos(theta_rim / 2))``, while the textbook relation is
    ``a_e / D = (1 - cos theta_tb) / 2``; equating them gives
    ``theta_tb = theta_rim / 2``. Confirmed against `_stock_2.engagement_at` on a
    half-plane at tool radius 1 (2026-08-22): radial depths of 0, 0.5, 1.0, 1.5
    and 2.0 measured rim arcs of 0, 120, 180, 240 and 360 degrees, and
    ``r (1 - cos(theta_rim / 2))`` reproduces every one of them exactly.
    `tests/benchmarks/test_quality.py` pins that table against the kernel, so the
    conversion cannot be "corrected" back on algebra alone.

    THE LIMIT OF IT. The equality is exact for a straight wall and approximate
    for the curved, partly-cleared material a trochoid actually meets. It is the
    standard one-sided assumption every mechanistic chip model makes, and it is
    an assumption, not a measurement.

    Args:
        rim_arc_deg: The engaged rim arc in degrees, as the kernel reports it.

    Returns:
        The textbook engagement angle in degrees, in ``[0, 180]``.
    """
    return 0.5 * max(0.0, min(rim_arc_deg, FULL_TURN_DEG))


def radial_immersion(rim_arc_deg: float) -> float:
    """Radial immersion ``a_e / D`` for an engaged rim arc of *rim_arc_deg*.

    ``a_e = r (1 - cos(theta_rim / 2))``, so against the diameter the tool size
    cancels and the ratio is ``(1 - cos(theta_rim / 2)) / 2``. Equivalently it is
    the textbook ``(1 - cos theta_tb) / 2`` evaluated at
    `textbook_engagement_deg`.

    It anchors where the KERNEL puts it, not where the textbook angle would: a
    rim arc of 180 degrees is HALF immersion (``a_e = r``), and a full turn is a
    slot (``a_e = D``). Reading 180 degrees as a slot would double every load
    number on the page, and reading a full turn through ``(1 - cos theta) / 2``
    would report a PLUNGE as zero immersion.

    Args:
        rim_arc_deg: The engaged rim arc in degrees.

    Returns:
        The immersion ratio, in ``[0, 1]``.
    """
    return 0.5 * (1.0 - math.cos(math.radians(textbook_engagement_deg(rim_arc_deg))))


def chip_thickness_ratio_from_rim(rim_arc_deg: float) -> float:
    """``h_ex / f_z`` for an engaged rim arc, converting the angle first.

    The metrics measure rim arcs; `chip_thickness_ratio` is stated in the
    textbook angle. This is the only place the two meet, so a caller cannot
    accidentally feed one to the other.

    Args:
        rim_arc_deg: The engaged rim arc in degrees.

    Returns:
        The ratio in ``[0, 1]``.
    """
    return chip_thickness_ratio(textbook_engagement_deg(rim_arc_deg))


def design_immersion(tea_cap_deg: float) -> float:
    """Radial immersion a cut running exactly at the engagement cap would carry.

    The commissioned load, and a property of the cap alone.

    Args:
        tea_cap_deg: The engagement cap in degrees.

    Returns:
        The design immersion ratio, in ``[0, 1]``.
    """
    return radial_immersion(tea_cap_deg)


def _immersion_steady_fraction(immersions: Sequence[Tuple[float, float]], band: float) -> float:
    """Share of cut length whose immersion sits within *band* of the median.

    A DISPERSION measure: it asks whether the cut is steady, wherever it sits,
    which is the property that keeps a spindle inside one stability lobe.

    Args:
        immersions: Radial immersion paired with the cut length it stands for.
        band: Half-width of the steady band, from `design_immersion`.

    Returns:
        The share, or zero when there is no cut length or no band.
    """
    total = sum(weight for _, weight in immersions)
    if total <= 0.0 or band <= 0.0:
        return 0.0
    median = _weighted_percentile(immersions, MEDIAN_PERCENTILE)
    return sum(weight for value, weight in immersions if abs(value - median) <= band) / total


def _immersion_at_design_fraction(immersions: Sequence[Tuple[float, float]], design: float, band: float) -> float:
    """Share of cut length actually cutting at the commissioned load.

    Distinct from `_immersion_steady_fraction`: a path can be perfectly steady
    at a tenth of the intended load, which is stable and slow, and this is the
    field that separates the two. It is the direct measure of whether an
    engagement-controlled path delivers the engagement it was asked for.

    Args:
        immersions: Radial immersion paired with the cut length it stands for.
        design: The design immersion, from `design_immersion`.
        band: Half-width of the band around it.

    Returns:
        The share, or zero when there is no cut length or no band.
    """
    total = sum(weight for _, weight in immersions)
    if total <= 0.0 or band <= 0.0:
        return 0.0
    return sum(weight for value, weight in immersions if abs(value - design) <= band) / total


def _immersion_excursions(immersions: Sequence[Tuple[float, float]], band: float) -> int:
    """Times the immersion leaves the steady band along the path.

    Counts inside-to-outside transitions only, and takes its initial state from
    the first sample, so a path that STARTS outside the band is not charged an
    excursion for merely beginning there.

    Args:
        immersions: Radial immersion paired with the cut length it stands for,
            in path order.
        band: Half-width of the steady band.

    Returns:
        The number of departures from the band.
    """
    if band <= 0.0:
        return 0
    median = _weighted_percentile(immersions, MEDIAN_PERCENTILE)
    inside_previous: Optional[bool] = None
    excursions = 0
    for value, _ in immersions:
        inside = abs(value - median) <= band
        if inside_previous and not inside:
            excursions += 1
        inside_previous = inside
    return excursions


def _weighted_engagement(motions: Sequence[MotionQuality]) -> List[Tuple[float, float]]:
    """Every sampled engagement paired with the path length it stands for.

    A motion of length L probed at k positions gives each position ``L / k``, so
    a 7 mm loop outweighs a 0.1 mm link instead of counting equally with it.

    Args:
        motions: The cut motions.

    Returns:
        ``(engagement_deg, length)`` pairs.
    """
    weighted: List[Tuple[float, float]] = []
    for motion in motions:
        if not motion.samples:
            continue
        share = motion.length / len(motion.samples)
        weighted.extend((sample.engagement_deg, share) for sample in motion.samples)
    return weighted


def _weighted_mean(values: Sequence[Tuple[float, float]]) -> float:
    """Length-weighted mean, or ``0.0`` for no weight."""
    total = sum(weight for _, weight in values)
    if total <= 0.0:
        return 0.0
    return sum(value * weight for value, weight in values) / total


def _weighted_variance(values: Sequence[Tuple[float, float]]) -> float:
    """Length-weighted variance, or ``0.0`` for no weight."""
    total = sum(weight for _, weight in values)
    if total <= 0.0:
        return 0.0
    mean = _weighted_mean(values)
    return sum(weight * (value - mean) ** 2 for value, weight in values) / total


def _weighted_percentile(values: Sequence[Tuple[float, float]], percentile: float) -> float:
    """The length-weighted *percentile* of *values*.

    The value below which that share of the CUT LENGTH lies, taken at the first
    sample whose cumulative weight reaches the target so the number reported is
    one that was measured rather than interpolated between two that were.

    Args:
        values: ``(value, weight)`` pairs, in any order.
        percentile: The percentile in ``(0, 100]``.

    Returns:
        The value, or ``0.0`` for an empty population.
    """
    if not values:
        return 0.0
    total = sum(weight for _, weight in values)
    if total <= 0.0:
        return max(value for value, _ in values)
    target = percentile / 100.0 * total
    cumulative = 0.0
    for value, weight in sorted(values):
        cumulative += weight
        if cumulative >= target:
            return value
    return max(value for value, _ in values)


def _max_engagement_gradient(motions: Sequence[MotionQuality]) -> float:
    """Steepest change of engagement per unit travel within any motion.

    Args:
        motions: The cut motions.

    Returns:
        Degrees per unit length; ``0.0`` when no motion carries two samples.
    """
    steepest = 0.0
    for motion in motions:
        for previous, current in zip(motion.samples, motion.samples[1:]):
            travel = current.distance - previous.distance
            if travel > 0.0:
                steepest = max(steepest, abs(current.engagement_deg - previous.engagement_deg) / travel)
    return steepest


def _max_engagement_step(motions: Sequence[MotionQuality]) -> float:
    """Largest peak-engagement change between back-to-back cut motions.

    Consecutive means ADJACENT IN THE TOOLPATH. A plunge, a retract, or a link
    between two cut motions shows up as a gap in the operation indices, and it is
    also the tool leaving the material -- so the load picked up afterwards is an
    entry, not a step.

    Args:
        motions: The cut motions, in toolpath order.

    Returns:
        The largest step in degrees; ``0.0`` when no two motions are adjacent.
    """
    step = 0.0
    for previous, current in zip(motions, motions[1:]):
        if current.index == previous.index + 1:
            step = max(step, abs(current.peak_engagement_deg - previous.peak_engagement_deg))
    return step


def _engagement_histogram(motions: Sequence[MotionQuality]) -> Tuple[Tuple[float, float, float], ...]:
    """Cut length spent in each band of `ENGAGEMENT_BANDS_DEG`.

    Args:
        motions: The cut motions.

    Returns:
        One ``(low_deg, high_deg, length)`` per band, in ascending order.
    """
    weighted = _weighted_engagement(motions)
    bands: List[Tuple[float, float, float]] = []
    for low, high in zip(ENGAGEMENT_BANDS_DEG, ENGAGEMENT_BANDS_DEG[1:]):
        bands.append((low, high, sum(weight for value, weight in weighted if low <= value < high)))
    top = ENGAGEMENT_BANDS_DEG[-1]
    if bands:
        low, high, length = bands[-1]
        bands[-1] = (low, high, length + sum(weight for value, weight in weighted if value >= top))
    return tuple(bands)


def _cut_air_alternations(survey: PathSurvey) -> int:
    """Transitions between cutting and not cutting, over the whole operation stream.

    A cut motion that engages nothing counts as air: the edge is not loaded, so
    it is not a thermal cycle whatever the operation is labelled.

    Args:
        survey: The replay's findings.

    Returns:
        The number of transitions.
    """
    states = [(motion.index, motion.is_engaged) for motion in survey.motions] + [(rapid.index, False) for rapid in survey.rapids]
    ordered = [engaged for _, engaged in sorted(states)]
    return sum(1 for previous, current in zip(ordered, ordered[1:]) if previous != current)


def _continuity_breaks(survey: PathSurvey, tool_radius: float) -> int:
    """Junctions between cut motions where one does not end where the next begins.

    Only the CUT motions are walked. A rapid deliberately jumps in XY -- that is
    what a rapid is for -- so counting its endpoints as breaks would report the
    path's design as a defect.

    Args:
        survey: The replay's findings.
        tool_radius: Tool radius, which scales the coincidence tolerance.

    Returns:
        The number of breaks.
    """
    tolerance = CONTINUITY_TOOL_RADIUS_FRACTION * tool_radius
    breaks = 0
    for previous, current in zip(survey.motions, survey.motions[1:]):
        if current.index != previous.index + 1:
            continue
        if math.hypot(current.start[0] - previous.end[0], current.start[1] - previous.end[1]) > tolerance:
            breaks += 1
    return breaks


def _junction_breaks(motions: Sequence[MotionQuality]) -> Tuple[int, int, int]:
    """Tangent breaks, curvature breaks, and reversals between adjacent cut motions.

    A junction is counted ONCE, in the most severe class it belongs to: a
    reversal is also a tangent break, and a tangent break is also a curvature
    break, so reporting all three for one corner would treble-count it.

    Args:
        motions: The cut motions, in toolpath order.

    Returns:
        ``(tangent_breaks, curvature_breaks, direction_reversals)``.

    """
    tangent = 0
    curvature = 0
    reversals = 0
    for previous, current in zip(motions, motions[1:]):
        if current.index != previous.index + 1:
            continue
        dot = previous.end_tangent[0] * current.start_tangent[0] + previous.end_tangent[1] * current.start_tangent[1]
        if dot < 0.0:
            reversals += 1
            tangent += 1
        elif dot < 1.0 - TANGENT_CONTINUITY_SLACK:
            tangent += 1
        elif previous.curvature != current.curvature:
            curvature += 1
    return tangent, curvature, reversals


def _loop_radii(motions: Sequence[MotionQuality]) -> List[float]:
    """Guide radii of the closed machining circles, in toolpath order."""
    return [motion.loop_radius for motion in motions if motion.kind is MotionKind.LOOP and motion.loop_radius is not None]


def _loop_runs(motions: Sequence[MotionQuality], rapids: Sequence[RapidMotion]) -> List[List[float]]:
    """Guide radii grouped into runs the tool performs without leaving material.

    A radius STEP is a load change the cutter actually passes through, so two
    loops may only be compared when the tool stayed down between them. A retract
    and its following plunge end a run: the cutter climbs to clearance height,
    rapids across, and re-enters somewhere else, and the difference in guide
    radius across that gap is not a step the machine ever performs.

    Splitting on rapid indices is what `_max_engagement_step` achieves with its
    ``index + 1`` test. That test cannot be reused here because loops within one
    chain are separated by the bridge motions between them and so are never
    index-adjacent; the run boundary, not adjacency, is the right notion.

    Args:
        motions: The cut motions, in toolpath order.
        rapids: The non-cutting motions, whose indices mark the breaks.

    Returns:
        One list of radii per continuous run; runs of fewer than two loops are
        kept so a caller can count them.
    """
    breaks = sorted(rapid.index for rapid in rapids)
    runs: List[List[float]] = [[]]
    position = 0
    for motion in motions:
        while position < len(breaks) and breaks[position] < motion.index:
            position += 1
            if runs[-1]:
                runs.append([])
        if motion.kind is MotionKind.LOOP and motion.loop_radius is not None:
            runs[-1].append(motion.loop_radius)
    return [run for run in runs if run]


def _loop_radius_cv(radii: Sequence[float]) -> float:
    """Coefficient of variation of the guide radius over the emitted loops.

    Args:
        radii: Loop radii, in toolpath order.

    Returns:
        Standard deviation over mean; ``0.0`` for fewer than two loops or a
        vanishing mean.
    """
    if len(radii) < 2:
        return 0.0
    mean = sum(radii) / len(radii)
    if mean <= 0.0:
        return 0.0
    variance = sum((radius - mean) ** 2 for radius in radii) / len(radii)
    return math.sqrt(variance) / mean


def _max_loop_radius_step(runs: Sequence[Sequence[float]], tool_radius: float) -> float:
    """Largest jump in guide radius between loops the tool runs back to back.

    Compared WITHIN a run and never across one, so a retract to clearance height
    followed by a plunge into a different chain cannot manufacture a step. An
    earlier form here compared consecutive entries of the flat radius list and
    did exactly that: on a 12x8 rectangle it reported 2.982 tool radii for the
    jump from the spine's last loop to a corner chain's first, which the cutter
    performs in the air.

    Args:
        runs: Loop radii grouped by `_loop_runs`.
        tool_radius: Tool radius, which normalises the step.

    Returns:
        The largest step; ``0.0`` when no run holds two loops.
    """
    if tool_radius <= 0.0:
        return 0.0
    steps = [abs(later - earlier) for run in runs for earlier, later in zip(run, run[1:])]
    if not steps:
        return 0.0
    return max(steps) / tool_radius
