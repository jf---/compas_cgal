"""The gate rules (a)-(f) of Task 0.3, each pinned on a synthetic path.

Every test here was written red against the pre-0.3 gate and states the rule it
pins in its docstring. The fixtures are the synthetic 12x8 pocket and helpers of
`tests.benchmarks.test_quality`, so a path here means what it means there.
"""

from __future__ import annotations

import dataclasses
import math

import pytest
from compas.geometry import Circle
from compas.geometry import Frame

from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.quality import measure_quality
from benchmarks.errors import MissingLoopSeamStationError
from benchmarks.quality_observations import CRITERION_NAMES
from benchmarks.quality_observations import _junction_step
from benchmarks.quality_observations import reduce_quality_evidence
from benchmarks.survey import MotionKind
from benchmarks.survey import survey_path
from benchmarks.units import OperationIndex
from compas_cgal import _stock_2
from compas_cgal.engagement import _cap_chord_ratio
from compas_cgal.stock import Stock
from tests.benchmarks.test_quality import FAST_SAMPLES
from tests.benchmarks.test_quality import SYNTHETIC
from tests.benchmarks.test_quality import SYNTHETIC_CAP_DEG
from tests.benchmarks.test_quality import SYNTHETIC_GRID
from tests.benchmarks.test_quality import _circle
from tests.benchmarks.test_quality import _cut_line
from tests.benchmarks.test_quality import _link_line
from tests.benchmarks.test_quality import _plunge
from tests.benchmarks.test_quality import _result
from tests.benchmarks.test_quality import _retract
from tests.benchmarks.test_quality import _survey
from tests.benchmarks.test_quality import _violations

# Probes along the wall-following cut of the rule (c) fixture: enough interior
# stations that the precondition (every interior reading strictly between half
# the cap and the cap) is a statement about the cut and not about one sample.
WALL_FOLLOW_SAMPLES = 8


def _evidence(operations: list):
    result = _result(operations)
    survey = survey_path(SYNTHETIC, result, samples_per_motion=FAST_SAMPLES)
    return survey, reduce_quality_evidence(SYNTHETIC, snapshot_toolpath(result), survey, grid=SYNTHETIC_GRID)


# --- (a) redundant = a CUT that removes nothing; a LINK at depth is a transit --


def test_a_link_at_depth_through_cleared_stock_is_a_transit_and_not_redundant() -> None:
    """R1c: staying down through cleared material is legal; it is reported as air at depth."""
    operations = [_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _link_line(2.0, 0.0, 1.5, 0.0)]
    quality = measure_quality(SYNTHETIC, _result(operations), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.elementary.redundant_operations == 0
    assert quality.speed.air_at_depth_length == pytest.approx(0.5)


def test_a_cut_that_removes_nothing_is_still_redundant() -> None:
    """The label decides (Ruling 0.3-C11): a CUT owns removal, so an empty CUT is waste."""
    operations = [_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _cut_line(2.0, 0.0, 1.5, 0.0)]
    quality = measure_quality(SYNTHETIC, _result(operations), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.elementary.redundant_operations == 1
    assert quality.speed.air_at_depth_length == 0.0


# --- (b) plunges into material <= centre-domain components --------------------


def test_a_second_plunge_into_virgin_stock_is_counted_and_gated() -> None:
    """One connected centre domain admits one plunge into material; the second is a defect."""
    operations = [_plunge(-3.0, 0.0), _circle(-3.0, 0.0, 1.5), _retract(-1.5, 0.0), _plunge(3.0, 0.0), _circle(3.0, 0.0, 1.5)]
    survey, evidence = _evidence(operations)
    assert evidence.coverage.centre_domain_components == 1
    assert survey.material_plunge_indices == (OperationIndex(0), OperationIndex(3))
    assert evidence.path_quality.elementary.material_plunges == 2
    assert evidence.assessment.material_plunges.required == 1
    assert any(line.startswith("[elementary] material plunges: measured 2, required <= 1") for line in _violations(evidence.assessment))


def test_a_single_plunge_into_one_centre_domain_is_within_the_rule() -> None:
    _survey_record, evidence = _evidence([_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0)])
    assert evidence.path_quality.elementary.material_plunges == 1
    assert not any("material plunges" in line for line in _violations(evidence.assessment))


def test_a_plunge_into_cleared_stock_is_not_a_material_plunge() -> None:
    """The plunge disk at (2, 0) lies inside the disk of radius 3 the first plunge and loop cleared."""
    operations = [_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _retract(2.0, 0.0), _plunge(2.0, 0.0), _cut_line(2.0, 0.0, 2.5, 0.0)]
    survey = _survey(operations)
    assert survey.plunge_indices == (OperationIndex(0), OperationIndex(3))
    assert survey.material_plunge_indices == (OperationIndex(0),)


# --- (c) slotting = a straight motion over the cap -------------------------------


def test_a_straight_cut_under_the_cap_is_not_slotting() -> None:
    """v-T2: a straight motion is judged by the same cap as a loop, not by half of it."""
    operations = [_plunge(-4.0, 0.0), _cut_line(-4.0, 0.0, 4.0, 0.0), _cut_line(4.0, 0.0, 4.0, 1.0), _cut_line(4.0, 1.0, -2.0, 1.0)]
    survey = _survey(operations, samples=WALL_FOLLOW_SAMPLES)
    wall_follow = survey.motions[2]
    # It stops at x = -2: run on to x = -4 and its last station meets the slot's
    # start corner, which reads exactly the cap (120 deg, an exceedance).
    interior = [sample.engagement_deg for sample in wall_follow.samples[1:-1]]
    assert all(SYNTHETIC_CAP_DEG / 2 < engagement < SYNTHETIC_CAP_DEG for engagement in interior), interior  # fixture precondition
    assert not wall_follow.cap_exceeded
    assert wall_follow.slot_exceeded is False


def test_slotting_is_the_straight_subset_of_cap_exceedances() -> None:
    """A slot over the cap is slotting; a loop over the cap is a cap exceedance only."""
    survey = _survey([_plunge(-4.0, 0.0), _cut_line(-4.0, 0.0, 4.0, 0.0), _circle(0.0, 2.5, 1.5)], samples=WALL_FOLLOW_SAMPLES)
    slot = survey.motions[0]
    assert slot.cap_exceeded and slot.slot_exceeded
    for motion in survey.motions:
        assert motion.slot_exceeded == (motion.kind is MotionKind.LINE and motion.cap_exceeded)


# --- (d) degenerate loops are reported, not gated ------------------------------


def test_a_degenerate_loop_under_the_cap_is_reported_but_not_a_violation() -> None:
    """g6-T4: a rho <= r circle in cleared stock is legitimate cleanup; the cap decides."""
    operations = [_plunge(0.0, 0.0), _circle(0.0, 0.0, 2.0), _cut_line(2.0, 0.0, 0.5, 0.0), _circle(0.0, 0.0, 0.5)]
    _survey_record, evidence = _evidence(operations)
    assert evidence.path_quality.elementary.degenerate_loops == 1
    assert "degenerate loops" not in CRITERION_NAMES
    assert not any("degenerate" in line for line in _violations(evidence.assessment))


# --- (e) a plunge's XY must equal the next motion's start ----------------------


def test_a_plunge_that_lands_away_from_the_next_cut_is_misplaced() -> None:
    """v-T9: a plunge -> cut teleport is otherwise unchecked."""
    entry = Circle(2.0, frame=Frame([0.0, 0.0, 0.0])).point_at(0.0)
    assert (float(entry[0]), float(entry[1])) == (2.0, 0.0)  # fixture precondition: bit-exact entry
    good = measure_quality(SYNTHETIC, _result([_plunge(2.0, 0.0), _circle(0.0, 0.0, 2.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    bad = measure_quality(SYNTHETIC, _result([_plunge(1.0, 0.0), _circle(0.0, 0.0, 2.0)]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert good.elementary.misplaced_plunges == 0
    assert bad.elementary.misplaced_plunges == 1


def test_a_plunge_not_followed_by_a_cut_motion_is_misplaced() -> None:
    """A plunge whose next operation is not a cut-plane motion has no start to land on."""
    _survey_record, evidence = _evidence([_plunge(2.0, 0.0), _circle(0.0, 0.0, 2.0), _retract(2.0, 0.0), _plunge(2.0, 0.0), _retract(2.0, 0.0)])
    assert evidence.assessment.attribution.misplaced_plunge_operations == (OperationIndex(3),)
    assert any(line.startswith("[elementary] misplaced plunges: measured 1") for line in _violations(evidence.assessment))


# --- (f) the engagement step is read at the junction ---------------------------


def test_the_engagement_step_is_read_at_the_junction_and_not_from_peaks() -> None:
    """v-T9: a step is what the edge feels crossing a junction, not two motions' peaks.

    The junction value is the next motion's first station with a non-empty
    prefix (index 1; station 0 reads 0 deg by construction) against the
    previous motion's last station (Ruling 0.3-f-station).
    """
    operations = [_plunge(-4.0, 0.0), _cut_line(-4.0, 0.0, 0.0, 0.0), _cut_line(0.0, 0.0, 0.0, 3.0)]
    survey, evidence = _evidence(operations)
    first, second = survey.motions
    expected = abs(second.samples[1].engagement_deg - first.samples[-1].engagement_deg)
    assert expected != pytest.approx(abs(second.peak_engagement_deg - first.peak_engagement_deg))  # fixture precondition
    assert float(evidence.assessment.max_engagement_step.measured) == pytest.approx(expected)


def _rim_deg(stock: Stock, x: float, y: float) -> float:
    """Static engaged rim arc at (x, y) on *stock*, read by the kernel directly."""
    return math.degrees(_stock_2.engagement_at(stock.raw, x, y, SYNTHETIC.tool_radius, _cap_chord_ratio(SYNTHETIC.tea_cap_rad), 0.0)[1])


def test_the_junction_after_a_loop_is_read_at_the_seam_with_the_full_turn_cleared() -> None:
    """Ruling 0.3-f-loop: a loop ends where it began, after its whole annulus is swept.

    The loop side's oracle replays the stock independently of the survey -- the
    plunge disk, then the loop's full annulus -- and reads the kernel statically
    at the seam. The next motion's side is its station 1 (Ruling
    0.3-f-station, unchanged). The loop's last per-motion station sits a sample
    spacing before the seam and reads otherwise.
    """
    operations = [_plunge(2.0, 0.0), _circle(0.0, 0.0, 2.0), _cut_line(2.0, 0.0, 3.5, 0.0)]
    survey, evidence = _evidence(operations)
    loop, line = survey.motions

    stock = Stock(SYNTHETIC.polygon, [])
    stock.subtract_disk(2.0, 0.0, SYNTHETIC.tool_radius)
    stock.subtract_circle_sweep(0.0, 0.0, 2.0, SYNTHETIC.tool_radius)
    loop_end = _rim_deg(stock, 2.0, 0.0)
    expected = abs(line.samples[1].engagement_deg - loop_end)

    assert expected != pytest.approx(abs(line.samples[1].engagement_deg - loop.samples[-1].engagement_deg))  # fixture precondition
    assert float(evidence.assessment.max_engagement_step.measured) == pytest.approx(expected)


def test_a_loop_record_without_its_seam_station_is_refused_by_name() -> None:
    """Fail loud (R2-6): a loop's junction is its seam; reading its last per-motion station instead would reopen I1."""
    operations = [_plunge(2.0, 0.0), _circle(0.0, 0.0, 2.0), _cut_line(2.0, 0.0, 3.5, 0.0)]
    survey, _evidence_record = _evidence(operations)
    loop, line = survey.motions
    assert loop.kind is MotionKind.LOOP and loop.exit_sample is not None  # the survey always records the seam
    with pytest.raises(MissingLoopSeamStationError, match=f"motion {loop.index}"):
        _junction_step(dataclasses.replace(loop, exit_sample=None), line)


def test_a_line_junction_reads_its_last_station_by_definition() -> None:
    """A line's hand-over station is its last one; no seam record is needed or read."""
    operations = [_plunge(0.0, 0.0), _cut_line(0.0, 0.0, 2.0, 0.0), _cut_line(2.0, 0.0, 2.0, 1.5)]
    survey, _evidence_record = _evidence(operations)
    first, second = survey.motions
    assert first.kind is MotionKind.LINE
    expected = abs(second.samples[1].engagement_deg - first.samples[-1].engagement_deg)
    assert _junction_step(dataclasses.replace(first, exit_sample=None), second) == expected
