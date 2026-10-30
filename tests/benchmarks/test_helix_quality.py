"""The quality reduction scores a HELIX entry (Task 1.2, set C deferred items Q-4..Q-11, set B M7 and C4/C5).

Rulings in force (W/progress.md, USER 2026-10-16):

* Q3(a): `material_plunges` counts PLUNGE + HELIX into material; the key is kept.
* Q3(b): `ramp_angle_violations` is a new GATED criterion, decided exactly by
  `Helix.ramp_within` against the machine's ramp surrogate.
* E8 / N7: the helix is not cap-sampled; rule (e) covers the helix exit.
* D4: `max_loop_radius_step` is a <= 2 r verdict row, no longer ratcheted.
* Q6 / M7: the helix is cutting; its length is in the cutting length.

HELIX operations come from the real emitter (`tests.benchmarks.helix_fixtures`).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import InvalidMachineModelError
from benchmarks.errors import MissingMachineModelError
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.models import GATE_MACHINE
from benchmarks.models import MachineModel
from benchmarks.quality import helix_space_curvature
from benchmarks.quality import machine_outcome
from benchmarks.quality import ramp_angle_violations
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import assess_path_quality
from benchmarks.quality_observations import reduce_quality_evidence
from benchmarks.stage_acceptance import STAGE_ACCEPTANCE
from benchmarks.stage_acceptance import VectorRow
from benchmarks.survey import PathSurvey
from benchmarks.survey import survey_path
from compas_cgal.helix import Helix
from compas_cgal.helix_units import DEFAULT_MAX_RAMP_ANGLE_DEG
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import line_operation
from compas_cgal.pipeline.emit import loop_operation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY
from tests.benchmarks.helix_fixtures import CUT_Z
from tests.benchmarks.helix_fixtures import MEMBER_RADIUS
from tests.benchmarks.helix_fixtures import POCKET
from tests.benchmarks.helix_fixtures import RAMP_DEG
from tests.benchmarks.helix_fixtures import TOP_Z
from tests.benchmarks.helix_fixtures import emitted
from tests.benchmarks.helix_fixtures import first_member
from tests.benchmarks.helix_fixtures import make_helix
from tests.benchmarks.helix_fixtures import result
from tools.gate_measure import CRITERIA
from tools.gate_measure import MeasuredCriterion
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import load_gate_vectors

ROOT = Path(__file__).resolve().parents[2]
SAMPLES = 8
# 12 mm / 120 = 0.1 mm cells: the coarsest grid that still resolves residue at r = 1 (CoarseCoverageGridError).
COARSE_GRID = 120

# A ramp cap just below the fixture's 2 deg: the fixture helix descends 4 mm in
# ceil(4 / (2 pi tan 2deg)) = 19 turns on rho_h = 1, i.e. pitch 4/19 = 0.2105,
# while 2 pi tan(1.9 deg) = 0.2084. One violation, decided exactly.
STEEPER_THAN_CAP_DEG = 1.9

# The re-baseline record of the "1.2" rules bump (Ruling 1.2-D2), as the
# committed `benchmarks/results/gate_vectors.json` must carry it.
REBASELINE_1_2 = {
    "from": "0.3",
    "approval": "Ruling 1.2-D2: USER-APPROVED 2026-10-16",
    "criteria": ["material_plunges", "ramp_angle_violations", "max_engagement_step_deg", "cap_exceedances", "misplaced_plunges"],
    "retired": ["max_loop_radius_step"],
}


def _evidence(operations: list[ToolpathOperation]) -> tuple[PathSurvey, QualityEvidence]:
    path = result(operations)
    survey = survey_path(POCKET, path, samples_per_motion=SAMPLES)
    return survey, reduce_quality_evidence(POCKET, snapshot_toolpath(path), survey, grid=COARSE_GRID)


def _loop(centre: tuple[float, float], radius: float, *, clockwise: bool = True) -> ToolpathOperation:
    loop = IncidentLoop.build(Point2[WorldXY].build(*centre), GuideRadius.build(radius), Radian(0.0))
    return loop_operation(loop, CutZ.build(CUT_Z), 0, clockwise=clockwise)


def _retract(helix: Helix) -> ToolpathOperation:
    entry = helix.entry_xy
    return line_operation(entry, entry, CUT_Z, TOP_Z, OperationType.RETRACT, 0)


def _criterion(name: str) -> MeasuredCriterion:
    (found,) = [criterion for criterion in CRITERIA if criterion.spec.name == name]
    return found


# --- Q-4 / Q-5 / Q-6: the helix partition and material_plunges (Q3(a)) -------


def test_a_helix_into_material_is_a_material_plunge() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), first_member(helix)])
    assert evidence.assessment.attribution.material_plunge_operations == (0,)
    assert evidence.assessment.material_plunges.measured == 1
    assert evidence.path_quality.elementary.material_plunges == 1


def test_a_helix_into_its_own_bore_is_not_a_material_plunge() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), first_member(helix), _retract(helix), emitted(helix), first_member(helix)])
    assert evidence.assessment.attribution.material_plunge_operations == (0,)


def test_the_survey_binding_refuses_a_survey_that_drops_its_helix() -> None:
    helix = make_helix()
    path = result([emitted(helix), first_member(helix)])
    survey = survey_path(POCKET, path, samples_per_motion=SAMPLES)
    with pytest.raises(InvalidHeldPathEvidenceError, match="helix"):
        reduce_quality_evidence(POCKET, snapshot_toolpath(path), dataclasses.replace(survey, helices=()), grid=COARSE_GRID)


def test_the_survey_binding_refuses_a_helix_observation_off_its_snapshot_entry() -> None:
    helix = make_helix()
    path = result([emitted(helix), first_member(helix)])
    survey = survey_path(POCKET, path, samples_per_motion=SAMPLES)
    (observed,) = survey.helices
    moved = dataclasses.replace(survey, helices=(dataclasses.replace(observed, start=(observed.start[0] + 1.0, observed.start[1])),))
    with pytest.raises(InvalidHeldPathEvidenceError, match="helix"):
        reduce_quality_evidence(POCKET, snapshot_toolpath(path), moved, grid=COARSE_GRID)


# --- Q-7: a helix is no loop-radius-step operand, and it splits a run ---------


def test_the_helix_is_not_a_loop_radius_step_operand() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), first_member(helix)])
    assert evidence.assessment.attribution.max_loop_radius_step is None


def test_a_helix_splits_a_loop_run() -> None:
    helix = make_helix()
    before = _loop((3.0, 0.0), 2.0)
    _survey, evidence = _evidence([before, emitted(helix), first_member(helix)])
    assert evidence.assessment.attribution.max_loop_radius_step is None


# --- Q-8: continuity, tangent and the junction step at (helix, next) ---------


def test_the_helix_to_first_member_junction_is_a_measured_engagement_step() -> None:
    helix = make_helix()
    survey, evidence = _evidence([emitted(helix), first_member(helix)])
    (observed,) = survey.helices
    (member,) = survey.motions
    step = evidence.assessment.attribution.max_engagement_step
    assert step is not None
    assert (step.pair.previous, step.pair.current) == (0, 1)
    assert step.value == abs(member.samples[1].engagement_deg - observed.exit_sample.engagement_deg)


def test_the_helix_exit_is_g1_with_the_first_member() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), first_member(helix)])
    attribution = evidence.assessment.attribution
    assert attribution.continuity_break_pairs == ()
    assert attribution.tangent_break_pairs == ()
    # G2 only: the footprint (rho_h = r) and the member (1.5 r) differ in curvature.
    assert [(pair.previous, pair.current) for pair in attribution.curvature_break_pairs] == [(0, 1)]


def test_a_member_that_starts_off_the_helix_exit_is_a_continuity_break() -> None:
    helix = make_helix()
    entry = helix.entry_xy
    away = _loop((entry.x - MEMBER_RADIUS, entry.y + 0.5), MEMBER_RADIUS)
    _survey, evidence = _evidence([emitted(helix), away])
    assert [(pair.previous, pair.current) for pair in evidence.assessment.attribution.continuity_break_pairs] == [(0, 1)]


def test_a_reversed_member_is_a_tangent_break_at_the_helix_exit() -> None:
    helix = make_helix(clockwise=True)
    reversed_member = loop_operation(
        IncidentLoop.build(Point2[WorldXY].build(helix.entry_xy.x - MEMBER_RADIUS, helix.entry_xy.y), GuideRadius.build(MEMBER_RADIUS), Radian(0.0)),
        CutZ.build(CUT_Z),
        0,
        clockwise=False,
    )
    _survey, evidence = _evidence([emitted(helix), reversed_member])
    attribution = evidence.assessment.attribution
    assert [(pair.previous, pair.current) for pair in attribution.reversal_pairs] == [(0, 1)]
    assert [(pair.previous, pair.current) for pair in attribution.tangent_break_pairs] == [(0, 1)]


# --- Q-9 / N7: rule (e) covers the helix exit ---------------------------------


def test_a_helix_whose_exit_is_the_next_start_is_not_misplaced() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), first_member(helix)])
    assert evidence.assessment.attribution.misplaced_plunge_operations == ()


def test_a_helix_whose_exit_is_not_the_next_start_is_misplaced() -> None:
    helix = make_helix()
    entry = helix.entry_xy
    away = _loop((entry.x - MEMBER_RADIUS, entry.y + 0.5), MEMBER_RADIUS)
    _survey, evidence = _evidence([emitted(helix), away])
    assert evidence.assessment.attribution.misplaced_plunge_operations == (0,)


def test_a_helix_followed_by_no_cut_motion_is_misplaced() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), _retract(helix)])
    assert evidence.assessment.attribution.misplaced_plunge_operations == (0,)


# --- Q-10 / M7: the helix is cutting -------------------------------------------


def test_the_helix_length_is_cutting_length_and_the_lengths_are_consistent() -> None:
    helix = make_helix()
    survey, evidence = _evidence([emitted(helix), first_member(helix)])
    (member,) = survey.motions
    assert survey.helix_length == float(helix.length())
    assert survey.cut_length == member.length
    assert survey.total_length == survey.cut_length + survey.helix_length + survey.air_length
    assert evidence.path_quality.speed.cutting_length == member.length + float(helix.length())
    assert evidence.path_quality.speed.air_fraction == survey.air_length / survey.total_length


def test_reported_cells_account_for_the_helix() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), first_member(helix)])
    assert evidence.path_quality.speed.reentry_count == 1
    assert evidence.path_quality.longevity.material_entries == 1


def test_a_helix_into_material_is_a_cutting_state_for_alternations() -> None:
    helix = make_helix()
    _survey, evidence = _evidence([emitted(helix), _retract(helix)])
    # helix (cutting) -> retract (air): one transition.
    assert evidence.path_quality.longevity.cut_air_alternations == 1


# --- C4 / C5 / Q3(b): the machine ramp limit and ramp_angle_violations --------


def test_the_machine_ramp_limit_defaults_to_the_placeholder_and_builds_its_surrogate() -> None:
    machine = MachineModel.build()
    assert machine.max_ramp_angle_deg == DEFAULT_MAX_RAMP_ANGLE_DEG
    assert machine.helix_pitch_per_radius == HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(DEFAULT_MAX_RAMP_ANGLE_DEG))
    assert GATE_MACHINE == MachineModel.build()


@pytest.mark.parametrize("angle", [0.0, -1.0, 90.0, float("nan"), float("inf")])
def test_the_machine_refuses_a_ramp_angle_outside_the_open_quarter_turn(angle: float) -> None:
    with pytest.raises(InvalidMachineModelError):
        MachineModel.build(max_ramp_angle_deg=angle)


def test_the_helix_ramp_verdict_is_the_exact_surrogate_check() -> None:
    helix = make_helix()
    survey, _evidence_record = _evidence([emitted(helix), first_member(helix)])
    (observed,) = survey.helices
    at_cap = MachineModel.build(max_ramp_angle_deg=RAMP_DEG)
    steeper = MachineModel.build(max_ramp_angle_deg=STEEPER_THAN_CAP_DEG)
    assert observed.ramp_within_limit(at_cap) is helix.ramp_within(at_cap.helix_pitch_per_radius) is True
    assert observed.ramp_within_limit(steeper) is helix.ramp_within(steeper.helix_pitch_per_radius) is False


def test_ramp_angle_violations_counts_the_helices_steeper_than_the_machine_allows() -> None:
    helix = make_helix()
    survey, _evidence_record = _evidence([emitted(helix), first_member(helix)])
    assert ramp_angle_violations(survey, MachineModel.build(max_ramp_angle_deg=RAMP_DEG)) == ()
    assert ramp_angle_violations(survey, MachineModel.build(max_ramp_angle_deg=STEEPER_THAN_CAP_DEG)) == (0,)


def test_ramp_angle_violations_needs_a_machine() -> None:
    helix = make_helix()
    survey, _evidence_record = _evidence([emitted(helix), first_member(helix)])
    with pytest.raises(MissingMachineModelError):
        ramp_angle_violations(survey, None)


def test_a_path_without_a_helix_has_no_ramp_violation() -> None:
    survey, _evidence_record = _evidence([_loop((0.0, 0.0), 1.5)])
    assert ramp_angle_violations(survey, GATE_MACHINE) == ()


# --- the gate criteria, rows and the 1.2 rebaseline record --------------------


def test_ramp_angle_violations_is_a_gated_criterion_read_against_the_gate_machine() -> None:
    criterion = _criterion("ramp_angle_violations")
    assert (criterion.spec.direction, criterion.spec.decimals, criterion.spec.required) == ("lower", 0, "0")
    helix = make_helix()
    survey, evidence = _evidence([emitted(helix), first_member(helix)])
    assert criterion.read is not None
    assert criterion.read(evidence, survey) == len(ramp_angle_violations(survey, GATE_MACHINE)) == 0


def test_the_loop_step_is_a_coverage_verdict_not_a_ratchet() -> None:
    criterion = _criterion("max_loop_radius_step")
    assert (criterion.spec.direction, criterion.spec.required) == ("reported", None)
    rows = [row for row in STAGE_ACCEPTANCE["stage-1"] if isinstance(row, VectorRow) and row.criterion == "max_loop_radius_step"]
    assert [(row.case, row.relation, row.bound) for row in rows] == [("*/pipeline/*", "<=", "2")]


def test_the_1_2_stage_rows_are_in_force() -> None:
    rows = {(row.case, row.criterion, row.relation, row.bound) for row in STAGE_ACCEPTANCE["stage-1"] if isinstance(row, VectorRow)}
    assert ("cap-120-default/pipeline/*", "cap_exceedances", "==", "0") in rows
    assert ("cap-120-default/pipeline/*", "max_engagement_step_deg", "<=", "120") in rows
    assert ("cap-120-default/pipeline/*", "ramp_angle_violations", "==", "0") in rows
    assert ("cap-120-default/pipeline/*", "material_plunges", "==", "1") in rows


def test_the_1_2_rebaseline_record_validates() -> None:
    gated = tuple(criterion.spec.name for criterion in CRITERIA if criterion.spec.direction == "lower")
    record = load_gate_vectors(ROOT / GATE_VECTORS_PATH).rebaseline
    assert record is not None
    assert record.to_json() == REBASELINE_1_2
    assert set(record.criteria) <= set(gated)
    assert record.retired == ("max_loop_radius_step",)


def test_assess_path_quality_and_the_reduction_agree_on_a_helix_path() -> None:
    helix = make_helix()
    path = result([emitted(helix), first_member(helix)])
    survey = survey_path(POCKET, path, samples_per_motion=SAMPLES)
    evidence = reduce_quality_evidence(POCKET, snapshot_toolpath(path), survey, grid=COARSE_GRID)
    assert assess_path_quality(POCKET, snapshot_toolpath(path), survey, evidence.coverage) == evidence.assessment


# --- fix round 1 (W/task-1.2-g2-review.md I1, I2, M1) ---------------------------


def test_the_survey_binding_refuses_a_helix_observation_whose_helix_is_not_its_snapshot() -> None:
    """I1: the gated ramp verdict reads `HelixObservation.helix`, so the binding checks the whole helix."""
    helix = make_helix()
    path = result([emitted(helix), first_member(helix)])
    survey = survey_path(POCKET, path, samples_per_motion=SAMPLES)
    (observed,) = survey.helices
    # Same exact entry P = (1, 0), other footprint: centre (0.1, 0), radius 0.9.
    other = make_helix(centre=(0.1, 0.0), radius=0.9)
    assert (other.entry_xy.x, other.entry_xy.y) == (helix.entry_xy.x, helix.entry_xy.y)
    swapped = dataclasses.replace(survey, helices=(dataclasses.replace(observed, helix=other),))
    with pytest.raises(InvalidHeldPathEvidenceError, match="helix"):
        reduce_quality_evidence(POCKET, snapshot_toolpath(path), swapped, grid=COARSE_GRID)


def test_the_machine_outcome_times_the_helix_at_its_curvature_limited_feed() -> None:
    """I2: cutting seconds cover the helix like `cutting_length` does, at ``min(feed, sqrt(a / kappa))`` on its space curvature."""
    helix = make_helix()
    survey, _evidence_record = _evidence([emitted(helix), first_member(helix)])
    (observed,) = survey.helices
    (member,) = survey.motions
    outcome = machine_outcome(survey, GATE_MACHINE)
    member_seconds = member.length / GATE_MACHINE.curvature_limited_feed_mm_per_s(member.curvature)
    helix_seconds = observed.length / GATE_MACHINE.curvature_limited_feed_mm_per_s(helix_space_curvature(observed))
    expected = member_seconds + helix_seconds
    assert 0.0 < helix_space_curvature(observed) < observed.curvature
    assert outcome.cutting_seconds == expected
    assert outcome.mean_cutting_feed_mm_per_s == (survey.cut_length + survey.helix_length) / expected


def test_a_helix_into_its_own_bore_is_cut_length_and_air_for_alternations_like_a_loop_that_removes_nothing() -> None:
    """M1: one convention with 1.1's non-removing loop: a feed move (cutting length), unloaded (air for alternations)."""
    helix = make_helix()
    survey, evidence = _evidence([emitted(helix), first_member(helix), _retract(helix), emitted(helix), first_member(helix)])
    first, second = survey.helices
    assert (first.into_material, second.into_material) == (True, False)
    assert evidence.path_quality.speed.cutting_length == survey.cut_length + first.length + second.length
    motion_states = [(motion.index, motion.is_engaged) for motion in survey.motions]
    rapid_states = [(rapid.index, False) for rapid in survey.rapids]
    helix_states = [(helix.index, helix.into_material) for helix in survey.helices]
    states = sorted(motion_states + rapid_states + helix_states)
    engaged = [state for _index, state in states]
    assert evidence.path_quality.longevity.cut_air_alternations == sum(1 for previous, current in zip(engaged, engaged[1:]) if previous != current)
    assert engaged[3] is False
