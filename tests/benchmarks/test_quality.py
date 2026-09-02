"""The machining-quality gate: what a toolpath must satisfy to be worth running.

THE ASSERTIONS BELOW ARE WRITTEN AT THE STANDARD A CAM ENGINEER DEMANDS, NOT AT
WHAT TODAY'S GENERATORS ACHIEVE, AND MOST OF THEM FAIL. That is the point of the
file. A gate calibrated to pass on today's paths would certify a path that
plunges into a corner, cuts straight back out of it, and leaves four loops of a
sixtieth of the tool radius behind -- all of which passed every test that existed
before `benchmarks/quality.py`, because none of them asked whether the path was
any good, only whether it was legal.

Each gate test evaluates EVERY criterion and reports all of them together,
measured against required, instead of stopping at the first. A quality gate whose
failure message is one line is a quality gate nobody can act on.

The machinery tests above the gate hold the measurement itself: they run on
synthetic operation streams, cost milliseconds, and would catch a metric that
silently stopped measuring long before the slow gate noticed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

import benchmarks.quality as quality_module
from benchmarks.coverage import CoarseCoverageGridError
from benchmarks.coverage import CoverageEstimate
from benchmarks.coverage import measure_coverage
from benchmarks.errors import EmptyReachableRegionError
from benchmarks.errors import InvalidGridResolutionError
from benchmarks.errors import InvalidMachineModelError
from benchmarks.errors import InvalidMaterialModelError
from benchmarks.errors import InvalidMotionSampleCountError
from benchmarks.errors import MissingMachineModelError
from benchmarks.errors import MissingMaterialModelError
from benchmarks.errors import ZeroLengthToolpathError
from benchmarks.families.analytic import rectangle
from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_GENERATOR_NAMES
from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.models import MachineModel
from benchmarks.models import MaterialModel
from benchmarks.quality import CHIP_PLATEAU_DEG
from benchmarks.quality import IMMERSION_STEADY_BAND_FRACTION
from benchmarks.quality import FULL_TURN_DEG
from benchmarks.quality import CalibratedOutcome
from benchmarks.quality import design_immersion
from benchmarks.quality import measure_outcomes
from benchmarks.quality import chip_thickness_ratio_from_rim
from benchmarks.quality import radial_immersion
from benchmarks.quality import textbook_engagement_deg
from benchmarks.quality import DEGENERATE_LOOP_RATIO
from benchmarks.quality import MARGINAL_LOOP_RATIO
from benchmarks.quality import PathQuality
from benchmarks.quality import chip_thickness_ratio
from benchmarks.quality import machine_outcome
from benchmarks.quality import material_outcome
from benchmarks.quality import measure_quality
from benchmarks.quality import tool_life_outcome
from benchmarks.quality_observations import OperationPair
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import assess_path_quality
from benchmarks.spec import PocketSpec
from benchmarks.survey import _swept_area
from benchmarks.survey import MotionKind
from benchmarks.survey import PathSurvey
from benchmarks.survey import survey_path
from benchmarks.units import OperationIndex
from compas_cgal import _stock_2
from compas_cgal.engagement import _cap_chord_ratio
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

# ---------------------------------------------------------------------------
# What a good path must satisfy. Every threshold is derived from machining or
# from the geometry of the sweep, never read off a measurement.
# ---------------------------------------------------------------------------

# ELEMENTARY -- validity. Each of these is a defect that makes the path WRONG.
REQUIRED_UNCUT_FRACTION = 0.0
REQUIRED_GOUGING_MOTIONS = 0
REQUIRED_UNSAFE_RAPIDS = 0
REQUIRED_CONTINUITY_BREAKS = 0
REQUIRED_ZERO_LENGTH_MOTIONS = 0

# A loop that sweeps a disk instead of an annulus is a bore: the cutter never
# disengages, which is the one behaviour a trochoidal loop exists to prevent. It
# is not a small loop, it is the wrong motion. PHYSICS, not judgement.
REQUIRED_DEGENERATE_LOOPS = 0

# An operation leaving the stock exactly unchanged is cutting time, tool wear and
# program length bought for nothing. Decided exactly, so there is no quantity of
# it that a measurement uncertainty could excuse.
REQUIRED_REDUNDANT_OPERATIONS = 0

# CUT -- mechanics. The cap is the whole contract: a cut motion over it is the
# tool taking more load than the program allows. This count INCLUDES the entry
# cut after each plunge, which is a full slot for any generator entering solid
# stock without a helical or pre-drilled entry. That is a real defect and not an
# exemption -- a machinist handed a program that full-slots on entry sends it
# back -- so the gate counts it and the failure message says how many.
REQUIRED_CAP_EXCEEDANCES = 0

# A straight transfer between machining circles travels through material the
# loops already cleared. One engaging past half the cap is cutting on a motion
# nothing regulates.
REQUIRED_SLOTTING_MOTIONS = 0

# The load may not change by more than the ENTIRE budgeted load between two
# motions the tool runs back to back without leaving the cut. A swing wider than
# the cap is the shock the cap exists to bound, arriving as a step instead of a
# level. Expressed as a multiple of the instance's own cap so it scales with the
# program rather than with this pocket.
REQUIRED_ENGAGEMENT_STEP_CAP_MULTIPLE = 1.0

# A radius jump larger than the largest advance the generator can take is
# STRUCTURAL PROOF that the machining-circle centre locus jumped rather than
# tracked a curve, and it is derived rather than chosen. The clearance function
# along a skeleton chain is 1-Lipschitz (stated at
# `compas_cgal.engagement_toolpath.GUIDE_STEPOVER_ADVANCE_MULTIPLE`), so the
# radius change between consecutive stations of one chain cannot exceed the
# advance, and the advance cannot exceed `MAX_ADVANCE_TOOL_DIAMETERS` = 1.0 tool
# diameters = 2 tool radii. A larger step did not come from advancing along a
# chain. Held's Figure 5 avoids these jumps by putting the centres on a "middle"
# curve between the medial axis and the boundary, where clearance varies
# smoothly, instead of on the skeleton where it collapses at every branch.
REQUIRED_MAX_LOOP_RADIUS_STEP_TOOL_RADII = 2.0

# SPEED -- a tangent discontinuity is a corner the machine must come off feed
# for, and a trochoidal generator that stitches its loops to its bridges properly
# has none. This repository already treats tangent continuity as a design goal
# (`docs/toolpath_tangent_continuity.md`), so the gate holds it to that goal
# rather than to a softer number.
REQUIRED_TANGENT_BREAKS = 0

# NOT GATED, and each for a stated reason:
#   air_fraction          -- depends on the machine's rapid rate against its feed
#                            rate and on how many chains the skeleton carries,
#                            neither a property of the path.
#   recut_fraction        -- a trochoidal path scores ~0.9 BY CONSTRUCTION,
#                            because consecutive loops are meant to overlap.
#   marginal_loops        -- engineering judgement, deliberately not a gate.
#   material_entries      -- one per chain is unavoidable without a helical
#                            entry, and the cap gate already counts that defect.
#   curvature_breaks      -- a G2 break costs acceleration, not a stop, and no
#                            defensible count exists for a trochoidal path.
#   loop_radius_cv        -- REPORTED, and deliberately not gated, though it was
#                            asked for. On a rectangle the medial axis carries a
#                            spine of constant clearance and four branches on
#                            which clearance falls to zero at the corners, so a
#                            coefficient of variation near 0.5 is INTRINSIC to
#                            placing centres on the skeleton at all. The number
#                            that would make a threshold defensible is the one a
#                            middle-curve implementation produces, and it does
#                            not exist yet; picking one now would be reading the
#                            threshold off current behaviour, which is the single
#                            thing this file refuses to do. `max_loop_radius_step`
#                            gates the same defect from a DERIVED bound instead.

# ---------------------------------------------------------------------------
# Synthetic fixtures for the machinery tests.
# ---------------------------------------------------------------------------

SYNTHETIC_TOOL_DIAMETER = 2.0
SYNTHETIC_CAP_DEG = 120.0

# `MAX_CELL_TOOL_RADIUS_FRACTION` requires a cell at or below a tenth of the 1.0
# tool radius, so a 12-wide pocket needs at least 120 samples and a 6-wide one 60.
SYNTHETIC_GRID = 120
SMALL_POCKET_GRID = 80

# Few probes: the machinery tests assert classifications, not engagement values.
FAST_SAMPLES = 4

_SYNTHETIC_POLYGON = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
SYNTHETIC = PocketSpec.build(
    name="synthetic_12x8",
    family="analytic",
    polygon=_SYNTHETIC_POLYGON,
    tool_diameter=SYNTHETIC_TOOL_DIAMETER,
    tea_cap_deg=SYNTHETIC_CAP_DEG,
)


def _op(geometry, operation: OperationType, path_index: int = 0) -> ToolpathOperation:
    """One operation, on chain *path_index* unless a test cares otherwise."""
    return ToolpathOperation(geometry=geometry, operation=operation, path_index=path_index)


def _plunge(x: float, y: float, z_top: float = 4.0) -> ToolpathOperation:
    """A downward bore from the clearance height to the cutting plane."""
    return _op(Line([x, y, z_top], [x, y, 0.0]), OperationType.PLUNGE)


def _retract(x: float, y: float, z_top: float = 4.0) -> ToolpathOperation:
    """An upward rapid to the clearance height."""
    return _op(Line([x, y, 0.0], [x, y, z_top]), OperationType.RETRACT)


def _circle(x: float, y: float, radius: float, path_index: int = 0) -> ToolpathOperation:
    """A closed machining circle at the cutting plane."""
    return _op(Circle(radius, frame=Frame([x, y, 0.0])), OperationType.CUT, path_index)


def _cut_line(x0: float, y0: float, x1: float, y1: float, path_index: int = 0) -> ToolpathOperation:
    """A straight cut at the cutting plane."""
    return _op(Line([x0, y0, 0.0], [x1, y1, 0.0]), OperationType.CUT, path_index)


def _link_line(x0: float, y0: float, x1: float, y1: float) -> ToolpathOperation:
    """A straight LINK at the cutting plane, which the cut-plane model treats as a cut."""
    return _op(Line([x0, y0, 0.0], [x1, y1, 0.0]), OperationType.LINK)


def _result(operations: list) -> ToolpathResult:
    """A synthetic toolpath; the polyline is unused by every metric under test."""
    return ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))


def _survey(operations: list, samples: int = FAST_SAMPLES):
    """Survey a synthetic path on the synthetic pocket."""
    return survey_path(SYNTHETIC, _result(operations), samples_per_motion=samples)


# ---------------------------------------------------------------------------
# The chip-thickness physics. This is the page's centrepiece and the one formula
# a plausible-looking mistake survives review in.
# ---------------------------------------------------------------------------


def test_chip_thickness_plateaus_past_a_quarter_turn_and_never_falls() -> None:
    """``h_ex / f_z = sin(min(theta, 90 deg))`` -- a PLATEAU, not a sine over the whole range.

    Writing ``min(sin(theta), 1)`` instead makes the curve fall again past 90
    degrees, which would tell a reader that heavier engagement produces a lighter
    chip. It does not: past a quarter turn the chip is already at full thickness
    and more engagement adds engaged arc at that same thickness. This is the
    regression test for that error.
    """
    for degrees in (90.0, 91.0, 120.0, 150.0, 179.0, 180.0):
        assert chip_thickness_ratio(degrees) == pytest.approx(1.0), f"h_ex/f_z must stay at 1.0 at {degrees} deg"


def test_the_kernel_anchors_the_rim_arc_to_radial_depth_relation() -> None:
    """The measurement that decides the angle convention, pinned against the kernel.

    `_stock_2.engagement_at` reports the ENGAGED ARC OF THE RIM, which reaches a
    full turn when the cutter is surrounded. Feeding a half-plane at known radial
    depths and reading the arc back is the only way to know that, and it settles
    a question algebra alone got wrong once already: the textbook milling angle
    tops out at 180 degrees for a slot, this one tops out at 360.

    Without this test `radial_immersion` can be "corrected" to the textbook
    ``(1 - cos theta) / 2``, which reports a PLUNGE -- a full turn of engaged rim
    -- as ZERO immersion.
    """
    tool_radius = 1.0
    material = Polygon([[-20.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 20.0, 0.0], [-20.0, 20.0, 0.0]])
    stock = Stock(material, [])
    ratio = _cap_chord_ratio(math.radians(120.0))
    # (radial depth into the material, the rim arc the kernel must report)
    for depth, expected_rim_deg in ((0.0, 0.0), (0.5, 120.0), (1.0, 180.0), (1.5, 240.0), (2.0, 360.0)):
        centre_y = depth - tool_radius
        _total, run, _exceeded = _stock_2.engagement_at(stock.raw, 0.0, centre_y, tool_radius, ratio, 0.0)
        measured = math.degrees(run)
        assert measured == pytest.approx(expected_rim_deg, abs=1e-6), f"depth {depth}: kernel reports {measured} deg"
        assert 2.0 * tool_radius * radial_immersion(measured) == pytest.approx(depth, abs=1e-9)


def test_the_rim_arc_is_twice_the_textbook_engagement_angle() -> None:
    """A full turn of rim is a slot; half a turn is half immersion."""
    assert textbook_engagement_deg(360.0) == pytest.approx(180.0)
    assert textbook_engagement_deg(180.0) == pytest.approx(90.0)
    assert textbook_engagement_deg(0.0) == pytest.approx(0.0)
    assert radial_immersion(360.0) == pytest.approx(1.0)
    assert radial_immersion(180.0) == pytest.approx(0.5)
    assert chip_thickness_ratio_from_rim(180.0) == pytest.approx(1.0)
    assert chip_thickness_ratio_from_rim(FULL_TURN_DEG) == pytest.approx(1.0)


def test_chip_thickness_rises_as_a_sine_below_the_plateau() -> None:
    assert chip_thickness_ratio(0.0) == pytest.approx(0.0)
    assert chip_thickness_ratio(30.0) == pytest.approx(0.5)
    assert chip_thickness_ratio(60.0) == pytest.approx(math.sqrt(3.0) / 2.0)
    assert chip_thickness_ratio(CHIP_PLATEAU_DEG) == pytest.approx(1.0)


def test_chip_thickness_is_monotone_non_decreasing_across_the_whole_range() -> None:
    """The property the plateau bug broke, stated directly rather than at samples."""
    values = [chip_thickness_ratio(float(degrees)) for degrees in range(0, 361)]
    assert all(later >= earlier - 1e-12 for earlier, later in zip(values, values[1:]))


# ---------------------------------------------------------------------------
# Elementary: the exact tests.
# ---------------------------------------------------------------------------


def test_a_loop_no_wider_than_the_tool_is_a_bore_and_not_a_trochoid() -> None:
    """Physics, not judgement: at rho <= r the swept annulus has no hole at all."""
    radius = DEGENERATE_LOOP_RATIO * SYNTHETIC.tool_radius
    below = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 0.9 * radius)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    above = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 1.1 * radius)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert below.elementary.degenerate_loops == 1
    assert above.elementary.degenerate_loops == 0


def test_a_loop_with_a_hole_too_small_to_use_is_marginal_and_not_degenerate() -> None:
    """The judgement is a SEPARATE count, so a reader can discard it and keep the physics."""
    radius = 0.5 * (DEGENERATE_LOOP_RATIO + MARGINAL_LOOP_RATIO) * SYNTHETIC.tool_radius
    quality = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, radius)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.elementary.degenerate_loops == 0
    assert quality.elementary.marginal_loops == 1


def test_engagement_cannot_decide_redundancy_in_either_direction() -> None:
    """Why the exact test exists: engagement measures the RIM, redundancy the REMOVAL.

    A loop of radius 1.5 leaves an uncut core of radius 0.5 inside its swept
    annulus. A short cut across that core removes it while the cutter's rim runs
    entirely through cleared space -- zero engagement, material removed. Cutting
    it again is zero engagement, nothing removed. Two motions, the same measured
    engagement, opposite verdicts: no engagement threshold separates them.
    """
    operations = [_plunge(3.0, 0.0), _circle(0.0, 0.0, 1.5), _cut_line(-0.1, 0.0, 0.1, 0.0), _cut_line(-0.1, 0.0, 0.1, 0.0)]
    survey = _survey(operations)
    clears_the_core, repeats_it = survey.motions[1], survey.motions[2]
    assert clears_the_core.peak_engagement_deg == pytest.approx(0.0)
    assert repeats_it.peak_engagement_deg == pytest.approx(0.0)
    assert clears_the_core.removes_material is True
    assert repeats_it.removes_material is False


def test_a_cut_that_leaves_the_legal_centre_domain_is_a_gouge() -> None:
    """The cutter centre must stay a tool radius clear of every wall."""
    inside = _survey([_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 4.0, 0.0)])
    outside = _survey([_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 5.9, 0.0)])
    assert inside.motions[0].gouges is False
    assert outside.motions[0].gouges is True


def test_a_break_between_consecutive_cuts_is_counted() -> None:
    joined = measure_quality(
        SYNTHETIC, _result([_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 1.0, 0.0), _cut_line(1.0, 0.0, 2.0, 0.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID
    )
    broken = measure_quality(
        SYNTHETIC, _result([_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 1.0, 0.0), _cut_line(2.0, 0.0, 3.0, 0.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID
    )
    assert joined.elementary.continuity_breaks == 0
    assert broken.elementary.continuity_breaks == 1


def test_swept_area_is_the_closed_form_for_each_primitive() -> None:
    """A capsule, an annulus, and -- below the degeneracy boundary -- a disk."""
    tool = 1.0
    assert _swept_area(Line([0.0, 0.0, 0.0], [5.0, 0.0, 0.0]), tool) == pytest.approx(2.0 * tool * 5.0 + math.pi * tool**2)
    assert _swept_area(Circle(3.0, frame=Frame([0.0, 0.0, 0.0])), tool) == pytest.approx(4.0 * math.pi * 3.0 * tool)
    assert _swept_area(Circle(0.5, frame=Frame([0.0, 0.0, 0.0])), tool) == pytest.approx(math.pi * (0.5 + tool) ** 2)


def test_the_uncut_denominator_is_what_the_tool_can_reach_and_not_the_pocket() -> None:
    """The load-bearing definition: a corner a round tool cannot enter is not residue.

    The pocket is cleared by five capsules whose union covers the reachable
    region -- the centre domain's four edges plus its midline, each swept at the
    tool radius. Every reachable sample is then gone while the pocket's sharp
    corners still hold material, which is the only configuration that tells the
    two possible denominators apart.
    """
    spec = rectangle(width=6.0, height=4.0, tool_diameter=SYNTHETIC_TOOL_DIAMETER, tea_cap_deg=SYNTHETIC_CAP_DEG)
    stock = Stock(spec.polygon, [])
    assert measure_coverage(spec, stock, grid=SMALL_POCKET_GRID).uncut_fraction == pytest.approx(1.0)

    for x0, y0, x1, y1 in ((-2.0, -1.0, 2.0, -1.0), (-2.0, 0.0, 2.0, 0.0), (-2.0, 1.0, 2.0, 1.0), (-2.0, -1.0, -2.0, 1.0), (2.0, -1.0, 2.0, 1.0)):
        stock.subtract_capsule(x0, y0, x1, y1, spec.tool_radius)

    assert measure_coverage(spec, stock, grid=SMALL_POCKET_GRID).uncut_fraction == pytest.approx(0.0)
    assert stock.contains(-2.99, -1.99) is True


# ---------------------------------------------------------------------------
# Cut, speed and longevity reductions.
# ---------------------------------------------------------------------------


def test_the_engagement_histogram_accounts_for_the_whole_cut_length() -> None:
    """Every probed position is banded exactly once, so the bands sum to the cut."""
    survey = _survey([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _cut_line(0.0, 0.0, 3.0, 0.0)])
    quality = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _cut_line(0.0, 0.0, 3.0, 0.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    banded = sum(length for _low, _high, length in quality.longevity.engagement_length_histogram)
    assert banded == pytest.approx(survey.cut_length)


def test_a_straight_cut_has_no_curvature_and_a_loop_has_its_reciprocal_radius() -> None:
    survey = _survey([_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 3.0, 0.0), _circle(0.0, 0.0, 2.0)])
    assert survey.motions[0].curvature == pytest.approx(0.0)
    assert survey.motions[1].curvature == pytest.approx(0.5)


def test_a_reversal_is_counted_once_and_not_also_as_a_tangent_break() -> None:
    """A corner belongs to the most severe class it is in, or one corner counts three times."""
    quality = measure_quality(
        SYNTHETIC, _result([_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 3.0, 0.0), _cut_line(3.0, 0.0, 0.0, 0.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID
    )
    assert quality.speed.direction_reversals == 1
    assert quality.speed.tangent_breaks == 1
    assert quality.speed.curvature_breaks == 0


def test_air_counts_the_retract_and_never_the_plunge() -> None:
    """A plunge bores its own disk, so it is cutting even though it makes no XY progress."""
    operations = [_plunge(0.0, 0.0, z_top=4.0), _circle(0.0, 0.0, 2.0), _retract(0.0, 0.0, z_top=4.0)]
    quality = measure_quality(SYNTHETIC, _result(operations), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    plunge_length, retract_length, circle_length = 4.0, 4.0, 4.0 * math.pi
    assert quality.path_length == pytest.approx(plunge_length + retract_length + circle_length)
    assert quality.speed.air_fraction == pytest.approx(retract_length / (plunge_length + retract_length + circle_length))


def test_attributed_gate_observations_match_old_reducers_on_a_synthetic_path() -> None:
    result = _result(
        [
            _plunge(0.0, 0.0),
            _circle(0.0, 0.0, 1.0),
            _link_line(1.0, 0.0, 2.0, 0.0),
            _circle(-1.0, 0.0, 3.0),
        ]
    )
    quality = measure_quality(SYNTHETIC, result, samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    survey = survey_path(SYNTHETIC, result, samples_per_motion=FAST_SAMPLES)
    coverage = measure_coverage(SYNTHETIC, survey.final_stock, grid=SYNTHETIC_GRID)
    snapshot = snapshot_toolpath(result)
    assessment = assess_path_quality(SYNTHETIC, snapshot, survey, coverage)

    _assert_quality_parity(SYNTHETIC, quality, assessment, snapshot, survey)
    engagement_maximum = assessment.attribution.max_engagement_step
    loop_maximum = assessment.attribution.max_loop_radius_step
    assert engagement_maximum is not None
    assert engagement_maximum.pair == OperationPair.build(previous=OperationIndex(2), current=OperationIndex(3), operation_count=4)
    assert loop_maximum is not None
    assert loop_maximum.pair == OperationPair.build(previous=OperationIndex(1), current=OperationIndex(3), operation_count=4)


def test_quality_parity_rejects_a_wrong_in_bounds_maximum_pair() -> None:
    result = _result(
        [
            _plunge(0.0, 0.0),
            _circle(0.0, 0.0, 1.0),
            _link_line(1.0, 0.0, 2.0, 0.0),
            _circle(-1.0, 0.0, 3.0),
        ]
    )
    quality = measure_quality(SYNTHETIC, result, samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    survey = survey_path(SYNTHETIC, result, samples_per_motion=FAST_SAMPLES)
    coverage = measure_coverage(SYNTHETIC, survey.final_stock, grid=SYNTHETIC_GRID)
    snapshot = snapshot_toolpath(result)
    assessment = assess_path_quality(SYNTHETIC, snapshot, survey, coverage)
    maximum = assessment.attribution.max_engagement_step
    assert maximum is not None
    object.__setattr__(
        maximum,
        "pair",
        OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=4),
    )

    with pytest.raises(AssertionError):
        _assert_quality_parity(SYNTHETIC, quality, assessment, snapshot, survey)


# ---------------------------------------------------------------------------
# The model boundary: a coefficient is never guessed.
# ---------------------------------------------------------------------------


def test_a_model_derived_metric_refuses_to_guess_its_coefficients() -> None:
    """The whole point of the split: no material model, no chip thickness in millimetres."""
    survey = _survey([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0)])
    with pytest.raises(MissingMaterialModelError):
        material_outcome(survey, None)
    with pytest.raises(MissingMachineModelError):
        machine_outcome(survey, None)
    with pytest.raises(MissingMaterialModelError):
        tool_life_outcome(survey, None, MachineModel.build())
    with pytest.raises(MissingMachineModelError):
        tool_life_outcome(survey, MaterialModel.build(), None)


def test_an_out_of_range_coefficient_fails_at_construction() -> None:
    with pytest.raises(InvalidMaterialModelError):
        MaterialModel.build(feed_per_tooth_mm=0.0)
    with pytest.raises(InvalidMaterialModelError):
        MaterialModel.build(taylor_exponent=1.0)
    with pytest.raises(InvalidMaterialModelError):
        MaterialModel.build(min_chip_fraction_of_edge_radius=1.5)
    with pytest.raises(InvalidMachineModelError):
        MachineModel.build(max_acceleration_mm_per_s2=float("nan"))


def test_the_feed_ceiling_is_the_normal_acceleration_bound_clamped_to_the_feed() -> None:
    """``v <= sqrt(a_max / kappa)``, and a straight move is not bound by it."""
    machine = MachineModel.build(feed_rate_mm_per_min=6000.0, max_acceleration_mm_per_s2=2000.0)
    assert machine.feed_mm_per_s == pytest.approx(100.0)
    assert machine.curvature_limited_feed_mm_per_s(0.0) == pytest.approx(100.0)
    assert machine.curvature_limited_feed_mm_per_s(0.02) == pytest.approx(100.0)
    assert machine.curvature_limited_feed_mm_per_s(2.0) == pytest.approx(math.sqrt(1000.0))


def test_the_rubbing_floor_comes_from_the_edge_radius() -> None:
    material = MaterialModel.build(edge_radius_mm=0.005, min_chip_fraction_of_edge_radius=0.2)
    assert material.min_chip_thickness_mm == pytest.approx(0.001)


# ---------------------------------------------------------------------------
# Radial immersion. The second formula on this page a plausible-looking mistake
# survives review in, and one that DID survive: nothing pinned it until here.
# ---------------------------------------------------------------------------


def test_radial_immersion_anchors_where_the_kernel_puts_it() -> None:
    """``a_e / D = (1 - cos(theta_rim / 2)) / 2`` -- pinned at the arcs the kernel reports.

    NOT the textbook ``(1 - cos theta) / 2``. That relation is stated in the
    entry-to-exit angle, which tops out at 180 degrees for a slot; this kernel
    reports the ENGAGED RIM ARC, which reaches a full turn when the cutter is
    surrounded. The two differ by a factor of two, and
    `test_the_kernel_anchors_the_rim_arc_to_radial_depth_relation` measures which
    one this code is in rather than arguing it.

    So: a rim arc of 180 degrees is HALF immersion, and a full turn is the slot.
    """
    assert radial_immersion(0.0) == pytest.approx(0.0)
    assert radial_immersion(120.0) == pytest.approx(0.25)
    assert radial_immersion(180.0) == pytest.approx(0.5)
    assert radial_immersion(FULL_TURN_DEG) == pytest.approx(1.0)


def test_radial_immersion_never_falls_back_as_the_cutter_engages_further() -> None:
    """Monotone across the whole 0-360 range the engagement metric can produce.

    The textbook cosine turns back down past 180 degrees and would report a
    surrounded cutter -- a plunge -- as a light cut, or at 360 degrees as no cut
    at all. Halving the angle first keeps the relation inside the cosine's
    monotone quarter.
    """
    previous = -1.0
    for degrees in range(0, 361, 5):
        value = radial_immersion(float(degrees))
        assert value >= previous
        previous = value
    assert radial_immersion(FULL_TURN_DEG) == pytest.approx(1.0)


def test_the_design_immersion_is_a_property_of_the_cap_alone() -> None:
    """The load the program commissions, before any path is generated."""
    assert design_immersion(120.0) == pytest.approx(radial_immersion(120.0))
    assert design_immersion(120.0) == pytest.approx(0.25)
    assert design_immersion(FULL_TURN_DEG) == pytest.approx(1.0)


def test_the_steady_band_does_not_collapse_on_a_path_that_barely_cuts() -> None:
    """Regression: the band is scaled by the CAP, never by the path's own median.

    Three of these four identical circles retrace the first one's annulus and so
    cut nothing, which puts the median immersion at zero. A band taken as a
    fraction of that median is zero wide, and the metric then reports a
    catastrophic-looking zero that measures the median's smallness rather than
    the cut's steadiness. Scaled by the cap the band keeps its width and the
    path reads as what it is: entirely steady, at no load at all.
    """
    path = [_plunge(0.0, 0.0), _circle(0.0, 0.0, 1.5), _circle(0.0, 0.0, 1.5), _circle(0.0, 0.0, 1.5), _circle(0.0, 0.0, 1.5)]
    quality = measure_quality(SYNTHETIC, _result(path), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.cut.immersion_steady_fraction > 0.0
    assert IMMERSION_STEADY_BAND_FRACTION > 0.0


def test_a_steady_cut_and_a_correct_one_are_different_measurements() -> None:
    """Steadiness at a tenth of the commissioned load is stable, and useless.

    The path above is as steady as a path can be -- most of it cuts nothing at
    all, at a perfectly constant nothing -- while carrying none of the load the
    cap commissioned. One field must say so while the other does not, or the
    group cannot tell a good path from an idle one.
    """
    path = [_plunge(0.0, 0.0), _circle(0.0, 0.0, 1.5), _circle(0.0, 0.0, 1.5), _circle(0.0, 0.0, 1.5), _circle(0.0, 0.0, 1.5)]
    quality = measure_quality(SYNTHETIC, _result(path), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.cut.immersion_steady_fraction > quality.cut.immersion_at_design_fraction


def test_a_radius_jump_the_cutter_performs_in_the_air_is_not_a_load_step() -> None:
    """Regression: guide-radius steps are scoped to a continuous run.

    Two chains, each a pair of near-identical loops, with a retract and a plunge
    between them. Inside a chain the radius barely moves; across the gap it
    triples. The cutter is at clearance height for that change, so it is not a
    load step, and an earlier form that zipped the flat radius list reported it
    as the worst one on the path -- on every pocket in the gate corpus.

    `_max_engagement_step` already draws this distinction for engagement; this
    pins the same rule for radius.
    """
    path = [
        _plunge(-3.0, 0.0),
        _circle(-3.0, 0.0, 1.4),
        _cut_line(-3.0, 0.0, -2.6, 0.0),
        _circle(-2.6, 0.0, 1.5),
        _retract(-2.6, 0.0),
        _plunge(3.0, 0.0),
        _circle(3.0, 0.0, 4.4),
        _cut_line(3.0, 0.0, 3.4, 0.0),
        _circle(3.4, 0.0, 4.5),
    ]
    quality = measure_quality(SYNTHETIC, _result(path), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    within_run = 0.1 / SYNTHETIC.tool_radius
    across_the_gap = 2.9 / SYNTHETIC.tool_radius
    assert quality.cut.max_loop_radius_step == pytest.approx(within_run, abs=1e-6)
    assert quality.cut.max_loop_radius_step < across_the_gap


def test_two_chains_linked_at_cutting_depth_are_still_two_guides() -> None:
    """Regression: a run ends at a CHANGE OF CHAIN, not only at a retract.

    A generator that orders its chains so the tool never lifts links them at
    cutting depth. Splitting runs on rapids alone then puts the whole path in one
    run, and the radius change from one chain's last loop to the next chain's
    first is charged as a load step -- 3.988 tool radii on `rect_20x12`, against
    1.024 once chains are honoured.

    There is no rapid anywhere in this path. The two loops on chain 0 sit at
    1.40 and 1.45; chain 1 opens at 4.40. Only the within-chain 0.05 may count.
    """
    path = [
        _plunge(-3.0, 0.0),
        _circle(-3.0, 0.0, 1.40, path_index=0),
        _cut_line(-3.0, 0.0, -2.6, 0.0, path_index=0),
        _circle(-2.6, 0.0, 1.45, path_index=0),
        _cut_line(-2.6, 0.0, 3.0, 0.0, path_index=1),
        _circle(3.0, 0.0, 4.40, path_index=1),
        _cut_line(3.0, 0.0, 3.4, 0.0, path_index=1),
        _circle(3.4, 0.0, 4.45, path_index=1),
    ]
    quality = measure_quality(SYNTHETIC, _result(path), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.speed.retract_count == 0
    assert quality.cut.max_loop_radius_step == pytest.approx(0.05 / SYNTHETIC.tool_radius, abs=1e-6)


# ---------------------------------------------------------------------------
# Program feasibility: whether a control can execute the path at speed.
# ---------------------------------------------------------------------------


def test_a_path_of_arcs_and_a_path_of_lines_are_told_apart() -> None:
    """Arc fraction is what lets a control's look-ahead breathe, so it is measured."""
    arcs = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    lines = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _cut_line(-2.0, 0.0, 2.0, 0.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert arcs.program.arc_length_fraction == pytest.approx(1.0)
    assert lines.program.arc_length_fraction == pytest.approx(0.0)


def test_every_block_the_control_reads_is_counted_including_the_rapids() -> None:
    """A program is what the machine reads, not only what it cuts."""
    quality = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _retract(0.0, 0.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.program.cut_blocks == 1
    assert quality.program.block_count > quality.program.cut_blocks


def test_the_short_block_tail_is_reported_beside_the_median_it_hides_behind() -> None:
    """Starvation is a per-block event, so a mean would bury the population that causes it."""
    path = [_plunge(-5.0, 0.0), _cut_line(-5.0, 0.0, 5.0, 0.0)] + [_cut_line(x / 10.0, 1.0, (x + 1) / 10.0, 1.0) for x in range(-20, 0)]
    quality = measure_quality(SYNTHETIC, _result(path), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.program.cut_blocks == 21
    assert quality.program.short_block_length <= quality.program.median_block_length
    assert quality.program.min_block_length <= quality.program.short_block_length


# ---------------------------------------------------------------------------
# The calibration boundary.
# ---------------------------------------------------------------------------


def test_the_calibrated_outcomes_are_kept_out_of_the_computed_quality() -> None:
    """A calibrated estimate must never print in the same table as a computed fact.

    `PathQuality` is derived from the toolpath and the exact stock and is
    reproducible from this repository alone; `CalibratedOutcome` additionally
    depends on a material and a machine. Merging them would make the two
    indistinguishable to every reader downstream.
    """
    survey = _survey([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0)])
    quality = measure_quality(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    outcome = measure_outcomes(survey, material=MaterialModel.build(), machine=MachineModel.build())
    assert isinstance(outcome, CalibratedOutcome)
    assert not any("outcome" in name for name in vars(quality))
    assert outcome.machine.total_seconds > 0.0
    assert outcome.tool_life.life_fraction_consumed >= 0.0


# ---------------------------------------------------------------------------
# Failing loudly rather than measuring the wrong thing.
# ---------------------------------------------------------------------------


def test_a_grid_too_coarse_to_resolve_residue_fails_loudly() -> None:
    """Under-reporting residue silently is the one thing a coverage metric must not do."""
    with pytest.raises(CoarseCoverageGridError):
        measure_coverage(SYNTHETIC, Stock(SYNTHETIC.polygon, []), grid=10)


def test_a_grid_below_one_sample_is_refused() -> None:
    with pytest.raises(InvalidGridResolutionError):
        measure_coverage(SYNTHETIC, Stock(SYNTHETIC.polygon, []), grid=0)


def test_a_motion_probed_at_no_positions_is_refused() -> None:
    with pytest.raises(InvalidMotionSampleCountError):
        survey_path(SYNTHETIC, _result([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0)]), samples_per_motion=0)


def test_a_toolpath_with_no_length_is_refused_rather_than_divided_by() -> None:
    with pytest.raises(ZeroLengthToolpathError):
        measure_quality(SYNTHETIC, _result([]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)


def test_the_named_empty_reachable_error_exists_for_a_pocket_a_grid_cannot_see() -> None:
    """Pinned so the failure mode keeps its name even while no corpus instance hits it."""
    assert issubclass(EmptyReachableRegionError, Exception)


# ---------------------------------------------------------------------------
# The gate.
# ---------------------------------------------------------------------------


def _expected_maximum_pairs(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> tuple[OperationPair | None, OperationPair | None]:
    """Select the first maximum pairs directly from source observations."""
    operation_count = len(snapshot)
    engagement_pair: OperationPair | None = None
    engagement_value = 0.0
    for previous, current in zip(survey.motions, survey.motions[1:]):
        if current.index != previous.index + 1:
            continue
        value = abs(current.peak_engagement_deg - previous.peak_engagement_deg)
        if value > engagement_value:
            engagement_value = value
            engagement_pair = OperationPair.build(
                previous=OperationIndex(previous.index),
                current=OperationIndex(current.index),
                operation_count=operation_count,
            )

    chain_by_index = {int(operation.ordinal): operation.path_index for operation in snapshot}
    rapid_indices = sorted(rapid.index for rapid in survey.rapids)
    rapid_position = 0
    current_chain: int | None = None
    previous_loop: tuple[int, float] | None = None
    loop_pair: OperationPair | None = None
    loop_value = 0.0
    for motion in survey.motions:
        boundary = False
        while rapid_position < len(rapid_indices) and rapid_indices[rapid_position] < motion.index:
            rapid_position += 1
            boundary = True
        motion_chain = chain_by_index.get(motion.index, current_chain)
        if current_chain is not None and motion_chain != current_chain:
            boundary = True
        current_chain = motion_chain
        if boundary:
            previous_loop = None
        if motion.kind is not MotionKind.LOOP or motion.loop_radius is None:
            continue
        if previous_loop is not None:
            value = abs(motion.loop_radius - previous_loop[1]) / spec.tool_radius
            if value > loop_value:
                loop_value = value
                loop_pair = OperationPair.build(
                    previous=OperationIndex(previous_loop[0]),
                    current=OperationIndex(motion.index),
                    operation_count=operation_count,
                )
        previous_loop = (motion.index, motion.loop_radius)
    return engagement_pair, loop_pair


def _assert_quality_parity(
    spec: PocketSpec,
    quality: PathQuality,
    assessment: PathQualityAssessment,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> None:
    """Prove the additive observations reproduce every existing gate input."""
    old_and_new = (
        (quality.elementary.uncut_fraction, assessment.uncut_fraction.measured),
        (quality.elementary.gouging_motions, assessment.gouging_motions.measured),
        (quality.elementary.unsafe_rapids, assessment.unsafe_rapids.measured),
        (quality.elementary.continuity_breaks, assessment.continuity_breaks.measured),
        (quality.elementary.zero_length_motions, assessment.zero_length_motions.measured),
        (quality.elementary.degenerate_loops, assessment.degenerate_loops.measured),
        (quality.elementary.redundant_operations, assessment.redundant_operations.measured),
        (quality.cut.cap_exceedances, assessment.cap_exceedances.measured),
        (quality.cut.slotting_motions, assessment.slotting_motions.measured),
        (quality.cut.max_engagement_step_deg, assessment.max_engagement_step.measured),
        (quality.cut.max_loop_radius_step, assessment.max_loop_radius_step.measured),
        (quality.speed.tangent_breaks, assessment.tangent_breaks.measured),
    )
    assert all(old == new for old, new in old_and_new)

    required = (
        assessment.uncut_fraction.required,
        assessment.gouging_motions.required,
        assessment.unsafe_rapids.required,
        assessment.continuity_breaks.required,
        assessment.zero_length_motions.required,
        assessment.degenerate_loops.required,
        assessment.redundant_operations.required,
        assessment.cap_exceedances.required,
        assessment.slotting_motions.required,
        assessment.max_engagement_step.required,
        assessment.max_loop_radius_step.required,
        assessment.tangent_breaks.required,
    )
    assert required == (
        REQUIRED_UNCUT_FRACTION,
        REQUIRED_GOUGING_MOTIONS,
        REQUIRED_UNSAFE_RAPIDS,
        REQUIRED_CONTINUITY_BREAKS,
        REQUIRED_ZERO_LENGTH_MOTIONS,
        REQUIRED_DEGENERATE_LOOPS,
        REQUIRED_REDUNDANT_OPERATIONS,
        REQUIRED_CAP_EXCEEDANCES,
        REQUIRED_SLOTTING_MOTIONS,
        REQUIRED_ENGAGEMENT_STEP_CAP_MULTIPLE * spec.tea_cap_deg,
        REQUIRED_MAX_LOOP_RADIUS_STEP_TOOL_RADII,
        REQUIRED_TANGENT_BREAKS,
    )

    attribution = assessment.attribution
    expected_engagement_pair, expected_loop_pair = _expected_maximum_pairs(spec, snapshot, survey)
    actual_engagement_pair = None if attribution.max_engagement_step is None else attribution.max_engagement_step.pair
    actual_loop_pair = None if attribution.max_loop_radius_step is None else attribution.max_loop_radius_step.pair
    assert actual_engagement_pair == expected_engagement_pair
    assert actual_loop_pair == expected_loop_pair
    count_and_sources = (
        (assessment.gouging_motions.measured, attribution.gouging_operations),
        (assessment.unsafe_rapids.measured, attribution.unsafe_rapid_operations),
        (assessment.continuity_breaks.measured, attribution.continuity_break_pairs),
        (assessment.zero_length_motions.measured, attribution.zero_length_operations),
        (assessment.degenerate_loops.measured, attribution.degenerate_loop_operations),
        (assessment.redundant_operations.measured, attribution.redundant_operations),
        (assessment.cap_exceedances.measured, attribution.cap_exceeded_operations),
        (assessment.slotting_motions.measured, attribution.slotting_operations),
        (assessment.tangent_breaks.measured, attribution.tangent_break_pairs),
    )
    assert all(count == len(sources) for count, sources in count_and_sources)
    assert attribution.uncut_operations == ()

    for old_maximum, observation in (
        (quality.cut.max_engagement_step_deg, attribution.max_engagement_step),
        (quality.cut.max_loop_radius_step, attribution.max_loop_radius_step),
    ):
        if old_maximum == 0.0:
            assert observation is None
        else:
            assert observation is not None
            assert observation.value == old_maximum


def _violations(spec: PocketSpec, quality: PathQuality) -> list:
    """Every gate criterion *quality* fails, as ``measured against required`` lines.

    Args:
        spec: The instance, whose cap scales the engagement-step criterion.
        quality: The measured quality.

    Returns:
        One line per violated criterion, empty when the path is machinable.
    """
    step_limit = REQUIRED_ENGAGEMENT_STEP_CAP_MULTIPLE * spec.tea_cap_deg
    checks = (
        ("elementary", "uncut fraction", quality.elementary.uncut_fraction, REQUIRED_UNCUT_FRACTION, "{:.6f}"),
        ("elementary", "gouging motions", quality.elementary.gouging_motions, REQUIRED_GOUGING_MOTIONS, "{:g}"),
        ("elementary", "unsafe rapids", quality.elementary.unsafe_rapids, REQUIRED_UNSAFE_RAPIDS, "{:g}"),
        ("elementary", "continuity breaks", quality.elementary.continuity_breaks, REQUIRED_CONTINUITY_BREAKS, "{:g}"),
        ("elementary", "zero-length motions", quality.elementary.zero_length_motions, REQUIRED_ZERO_LENGTH_MOTIONS, "{:g}"),
        ("elementary", "degenerate loops", quality.elementary.degenerate_loops, REQUIRED_DEGENERATE_LOOPS, "{:g}"),
        ("elementary", "redundant operations", quality.elementary.redundant_operations, REQUIRED_REDUNDANT_OPERATIONS, "{:g}"),
        ("cut", "cap exceedances", quality.cut.cap_exceedances, REQUIRED_CAP_EXCEEDANCES, "{:g}"),
        ("cut", "slotting motions", quality.cut.slotting_motions, REQUIRED_SLOTTING_MOTIONS, "{:g}"),
        ("cut", "max engagement step (deg)", quality.cut.max_engagement_step_deg, step_limit, "{:.2f}"),
        ("cut", "max loop radius step (tool radii)", quality.cut.max_loop_radius_step, REQUIRED_MAX_LOOP_RADIUS_STEP_TOOL_RADII, "{:.3f}"),
        ("speed", "tangent breaks", quality.speed.tangent_breaks, REQUIRED_TANGENT_BREAKS, "{:g}"),
    )
    return [
        f"[{group}] {name}: measured {fmt.format(measured)}, required <= {fmt.format(required)}" for group, name, measured, required, fmt in checks if not measured <= required
    ]


def _report(pocket: str, generator: str, quality: PathQuality, violations: list) -> str:
    """A failure message a machinist could act on: the whole measurement, then the verdict.

    Args:
        pocket: Instance name.
        generator: Generator name.
        quality: The measured quality.
        violations: Lines from `_violations`.

    Returns:
        The formatted report.
    """
    groups = (
        (
            "ELEMENTARY",
            (
                ("uncut fraction", f"{quality.elementary.uncut_fraction:.6f}"),
                ("gouge free", str(quality.elementary.gouge_free)),
                ("rapid safety", str(quality.elementary.rapid_safety)),
                ("continuity breaks", str(quality.elementary.continuity_breaks)),
                ("zero-length motions", str(quality.elementary.zero_length_motions)),
                ("degenerate loops (rho <= r)", str(quality.elementary.degenerate_loops)),
                ("marginal loops (judgement)", str(quality.elementary.marginal_loops)),
                ("redundant operations", str(quality.elementary.redundant_operations)),
                ("recut fraction (reported)", f"{quality.elementary.recut_fraction:.4f}"),
            ),
        ),
        (
            "CUT",
            (
                ("max engagement (deg)", f"{quality.cut.max_engagement_deg:.2f}"),
                ("cap exceedances", str(quality.cut.cap_exceedances)),
                ("engagement p95 (deg)", f"{quality.cut.engagement_p95_deg:.2f}"),
                ("engagement variance (deg^2)", f"{quality.cut.engagement_variance_deg2:.1f}"),
                ("max chip thickness h_ex/f_z", f"{quality.cut.max_chip_thickness_ratio:.4f}"),
                ("low chip thickness h_ex/f_z", f"{quality.cut.low_chip_thickness_ratio:.4f}"),
                ("max engagement gradient (deg/mm)", f"{quality.cut.max_engagement_gradient_deg_per_length:.1f}"),
                ("max engagement step (deg)", f"{quality.cut.max_engagement_step_deg:.2f}"),
                ("slotting motions", str(quality.cut.slotting_motions)),
                ("immersion steady fraction", f"{quality.cut.immersion_steady_fraction:.4f}"),
                ("immersion at design load", f"{quality.cut.immersion_at_design_fraction:.4f}"),
                ("immersion excursions", str(quality.cut.immersion_excursions)),
                ("mean radial depth", f"{quality.cut.mean_radial_depth:.4f}"),
                ("radial depth variance", f"{quality.cut.radial_depth_variance:.4f}"),
                ("wall scallop height", f"{quality.cut.wall_scallop_height:.4f}"),
                ("loop radius CV (reported)", f"{quality.cut.loop_radius_cv:.4f}"),
                ("max loop radius step (r)", f"{quality.cut.max_loop_radius_step:.3f}"),
            ),
        ),
        (
            "SPEED",
            (
                ("cutting length", f"{quality.speed.cutting_length:.3f}"),
                ("air length", f"{quality.speed.air_length:.3f}"),
                ("air fraction (reported)", f"{quality.speed.air_fraction:.4f}"),
                ("max curvature (1/mm)", f"{quality.speed.max_curvature:.3f}"),
                ("tangent breaks", str(quality.speed.tangent_breaks)),
                ("curvature breaks (reported)", str(quality.speed.curvature_breaks)),
                ("direction reversals", str(quality.speed.direction_reversals)),
                ("retracts", str(quality.speed.retract_count)),
                ("re-entries", str(quality.speed.reentry_count)),
            ),
        ),
        (
            "PROGRAM",
            (
                ("blocks (cut + rapid)", str(quality.program.block_count)),
                ("cut blocks", str(quality.program.cut_blocks)),
                ("min block length", f"{quality.program.min_block_length:.5f}"),
                ("median block length", f"{quality.program.median_block_length:.5f}"),
                ("p05 block length", f"{quality.program.short_block_length:.5f}"),
                ("block length CV", f"{quality.program.block_length_cv:.4f}"),
                ("blocks per unit length", f"{quality.program.blocks_per_unit_length:.4f}"),
                ("arc length fraction", f"{quality.program.arc_length_fraction:.4f}"),
            ),
        ),
        (
            "LONGEVITY",
            (
                ("material entries", str(quality.longevity.material_entries)),
                ("cut/air alternations", str(quality.longevity.cut_air_alternations)),
                ("alternations per unit length", f"{quality.longevity.alternations_per_length:.5f}"),
            ),
        ),
    )
    lines = [f"\n{generator} on {pocket}:"]
    for title, rows in groups:
        lines.append(f"  {title}")
        lines.extend(f"    {name:<34} {value}" for name, value in rows)
    lines.append("  ENGAGEMENT HISTOGRAM (cut length per band)")
    lines.extend(f"    {low:>5.0f}-{high:<5.0f} deg{'':<16} {length:.3f}" for low, high, length in quality.longevity.engagement_length_histogram)
    lines.append(f"  NOT MACHINABLE -- {len(violations)} criteria failed:")
    lines.extend(f"    {line}" for line in violations)
    return "\n".join(lines) + "\n"


QUALITY_GATE_CASES = tuple(
    pytest.param(
        tea_cap_deg,
        generator_name,
        pocket_name,
        id=f"{cap_id}-{generator_name}-{pocket_name}",
    )
    for cap_id, tea_cap_deg in GATE_CAP_CASES
    for generator_name in GATE_GENERATOR_NAMES
    for pocket_name in GATE_POCKET_NAMES
)


@pytest.mark.parametrize(
    ("tea_cap_deg", "generator_name", "pocket_name"),
    QUALITY_GATE_CASES,
)
def test_the_generated_path_is_worth_running(
    tea_cap_deg: GateCapDegrees,
    generator_name: str,
    pocket_name: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every machining-quality criterion, on one pocket and one generator.

    All criteria are evaluated and reported together, so the failure message
    carries the full four-group measurement and the next person reads a diagnosis
    rather than a boolean.
    """
    spec = gate_pocket(pocket_name, tea_cap_deg=tea_cap_deg)
    result = GATE_GENERATORS[generator_name](spec)
    snapshot = snapshot_toolpath(result)
    surveys: list[PathSurvey] = []
    coverages: list[CoverageEstimate] = []
    original_survey_path = quality_module.survey_path
    original_measure_coverage = quality_module.measure_coverage

    def capture_survey(spec_arg: PocketSpec, result_arg: ToolpathResult, *, samples_per_motion: int) -> PathSurvey:
        survey = original_survey_path(spec_arg, result_arg, samples_per_motion=samples_per_motion)
        surveys.append(survey)
        return survey

    def capture_coverage(spec_arg: PocketSpec, stock: Stock, *, grid: int) -> CoverageEstimate:
        coverage = original_measure_coverage(spec_arg, stock, grid=grid)
        coverages.append(coverage)
        return coverage

    monkeypatch.setattr(quality_module, "survey_path", capture_survey)
    monkeypatch.setattr(quality_module, "measure_coverage", capture_coverage)
    quality = measure_quality(spec, result)
    assert len(surveys) == 1
    assert len(coverages) == 1
    assessment = assess_path_quality(spec, snapshot, surveys[0], coverages[0])
    _assert_quality_parity(spec, quality, assessment, snapshot, surveys[0])
    violations = _violations(spec, quality)
    assert not violations, _report(pocket_name, generator_name, quality, violations)
