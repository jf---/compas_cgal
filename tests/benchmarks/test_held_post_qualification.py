from __future__ import annotations

from dataclasses import replace
import re
from typing import cast

import pytest

import benchmarks.held_post_qualification as post_qualification_module
from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import HeldPathNotEligibleForPostQualificationError
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_post_qualification import HeldPostQualificationCandidate
from benchmarks.held_post_qualification import post_qualification_failures
from benchmarks.held_post_qualification import require_post_qualification_candidate
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.quality import PathQuality
from benchmarks.quality import _cut
from benchmarks.quality import _elementary
from benchmarks.quality import _longevity
from benchmarks.quality import _program
from benchmarks.quality import _speed
from benchmarks.quality_observations import CRITERION_NAMES
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import _build_record
from benchmarks.quality_observations import assess_path_quality
from benchmarks.survey import EngagementSample
from benchmarks.survey import MotionKind
from benchmarks.survey import PathSurvey
from benchmarks.survey import RapidMotion
from benchmarks.units import OperationIndex
from benchmarks.units import Seconds
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement import EngagementReport
from compas_cgal.engagement import OperationEngagement
from compas_cgal.toolpath import OperationType
from tests.benchmarks.test_quality_observations import _motion
from tests.benchmarks.test_quality_observations import _snapshot
from tests.benchmarks.test_quality_observations import _survey


def _sample(*, engagement: float = 20.0, cap: bool = False, inside: bool = True) -> EngagementSample:
    return EngagementSample(0.0, Point2[WorldXY].build(0.0, 0.0), engagement, cap, inside)


def _quality(
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    coverage: CoverageEstimate,
) -> QualityEvidence:
    assessment = assess_path_quality(survey.spec, snapshot, survey, coverage)
    path_quality = PathQuality(
        elementary=_elementary(survey.spec, survey, coverage.remaining_area, assessment),
        cut=_cut(survey.spec, survey, coverage.wall_scallop_height, assessment),
        speed=_speed(survey, assessment),
        longevity=_longevity(survey, 1),
        program=_program(survey),
        cut_operations=len(survey.motions),
        path_length=survey.total_length,
    )
    return _build_record(
        QualityEvidence,
        {"path_quality": path_quality, "assessment": assessment, "coverage": coverage},
    )


def _characterization(
    *,
    criterion: str | None = None,
    unresolved: bool = False,
    all_open: bool = False,
) -> HeldFigure5Characterization:
    case = load_held_reference_case("figure5")
    if criterion == "degenerate loops":
        circle_indices = {0}
    elif criterion == "max loop radius step (tool radii)" or all_open:
        circle_indices = {0, 1}
    else:
        circle_indices = set()
    snapshot = _snapshot(
        3,
        operation_roles={2: OperationType.RETRACT},
        circle_indices=circle_indices,
        horizontal_role_indices={2},
    )
    first = _motion(0, samples=(_sample(),))
    second = _motion(1, start=(1.0, 0.0), end=(2.0, 0.0), samples=(_sample(),))
    rapid = RapidMotion(2, OperationType.RETRACT, MotionKind.LINE, 1.0, False)
    coverage = CoverageEstimate(10, 10, 1.0, 10, 0, 0, 0.0)

    if criterion == "uncut fraction" or all_open:
        coverage = CoverageEstimate(10, 10, 1.0, 10, 1, 1, 0.0)
    if criterion == "gouging motions" or all_open:
        first = replace(first, samples=(_sample(inside=False),))
    if criterion == "unsafe rapids" or all_open:
        rapid = replace(rapid, horizontal_at_cut_plane=True)
    if criterion == "continuity breaks" or all_open:
        second = replace(second, start=(4.0, 0.0))
    if criterion == "zero-length motions" or all_open:
        first = replace(first, length=0.0)
    if criterion == "degenerate loops" or all_open:
        first = replace(first, kind=MotionKind.LOOP, loop_radius=0.5, curvature=2.0)
    if criterion == "redundant operations" or all_open:
        first = replace(first, removes_material=False)
    if criterion == "cap exceedances" or all_open:
        first = replace(first, samples=(_sample(cap=True, inside=not all_open),), cap_exceeded=True)
    if criterion == "slotting motions" or all_open:
        first = replace(first, slot_exceeded=True)
    if criterion == "max engagement step (deg)" or all_open:
        first = replace(first, samples=(_sample(engagement=0.0, cap=all_open, inside=not all_open),), cap_exceeded=all_open)
        second = replace(second, samples=(_sample(engagement=81.0),))
    if criterion == "max loop radius step (tool radii)" or all_open:
        first_radius = 0.5 if all_open else 1.1
        second_radius = 4.0 if all_open else 4.2
        first = replace(first, kind=MotionKind.LOOP, loop_radius=first_radius, curvature=1.0 / first_radius)
        second = replace(second, kind=MotionKind.LOOP, loop_radius=second_radius, curvature=1.0 / second_radius)
    if criterion == "tangent breaks" or all_open:
        second = replace(second, start_tangent=(-1.0, 0.0))

    survey = _survey(
        (first, second),
        (rapid,),
        spec=case.pocket_spec(),
        retract_indices=(OperationIndex(2),),
        source_snapshot=snapshot,
    )
    demonstrated = any(sample.cap_exceeded for motion in survey.motions for sample in motion.samples)
    unresolved_index = 1 if unresolved or all_open else None
    operations = [
        OperationEngagement(
            index,
            operation.operation,
            0.0,
            index not in ({0} if demonstrated else set()) | ({unresolved_index} if unresolved_index is not None else set()),
            0 if index == 2 else 3,
        )
        for index, operation in enumerate(snapshot)
    ]
    cap_violations = int(demonstrated) + int(unresolved_index is not None)
    audit = EngagementReport(2.0, 1.3962634015954636, operations, 0.0, cap_violations, 0)
    quality = _quality(snapshot, survey, coverage)
    return HeldFigure5Characterization.build(
        case=case,
        snapshot=snapshot,
        audit=audit,
        survey=survey,
        quality=quality,
        generation_seconds=Seconds(1.0),
        audit_seconds=Seconds(2.0),
        survey_seconds=Seconds(3.0),
        reduction_seconds=Seconds(4.0),
    )


@pytest.mark.parametrize("criterion", CRITERION_NAMES)
def test_each_open_quality_criterion_refuses_candidate(criterion: str) -> None:
    characterization = _characterization(criterion=criterion)

    assert any(failure.startswith(f"{criterion}: measured=") for failure in post_qualification_failures(characterization))
    with pytest.raises(HeldPathNotEligibleForPostQualificationError, match=re.escape(criterion)):
        require_post_qualification_candidate(characterization)


def test_witnessed_exceedance_refuses_both_reachable_open_conditions() -> None:
    characterization = _characterization(criterion="cap exceedances")

    assert post_qualification_failures(characterization) == (
        "demonstrated engagement exceedances: 1",
        "cap exceedances: measured=1, required=0, outcome=failure_observed",
    )
    with pytest.raises(HeldPathNotEligibleForPostQualificationError):
        require_post_qualification_candidate(characterization)


def test_unresolved_motion_refuses_candidate() -> None:
    characterization = _characterization(unresolved=True)

    assert post_qualification_failures(characterization) == ("unresolved TEA-audited operations: 1",)
    with pytest.raises(HeldPathNotEligibleForPostQualificationError):
        require_post_qualification_candidate(characterization)


def test_every_open_condition_has_stable_order_and_complete_error_context() -> None:
    characterization = _characterization(all_open=True)
    expected = (
        "demonstrated engagement exceedances: 1",
        "unresolved TEA-audited operations: 1",
        "uncut fraction: measured=0.1, required=0.0, outcome=failure_observed",
        "gouging motions: measured=1, required=0, outcome=failure_observed",
        "unsafe rapids: measured=1, required=0, outcome=outside_declared_tolerance",
        "continuity breaks: measured=1, required=0, outcome=outside_declared_tolerance",
        "zero-length motions: measured=1, required=0, outcome=criterion_violated",
        "degenerate loops: measured=1, required=0, outcome=criterion_violated",
        "redundant operations: measured=1, required=0, outcome=criterion_violated",
        "cap exceedances: measured=1, required=0, outcome=failure_observed",
        "slotting motions: measured=1, required=0, outcome=failure_observed",
        "max engagement step (deg): measured=81.0, required=80.0, outcome=failure_observed",
        "max loop radius step (tool radii): measured=3.5, required=2.0, outcome=criterion_violated",
        "tangent breaks: measured=1, required=0, outcome=outside_declared_tolerance",
    )
    failures = post_qualification_failures(characterization)

    assert failures == expected
    with pytest.raises(HeldPathNotEligibleForPostQualificationError) as exc_info:
        HeldPostQualificationCandidate.build(characterization)
    assert (
        str(exc_info.value)
        == "Held Figure 5 is not eligible for post qualification: " + "; ".join(failures) + ". Engagement counts: certified=0, demonstrated_exceeded=1, unresolved=1."
    )


def test_closed_characterization_builds_snapshot_bound_candidate() -> None:
    characterization = _characterization()

    assert post_qualification_failures(characterization) == ()
    candidate = HeldPostQualificationCandidate.build(characterization)

    assert candidate.characterization is characterization
    assert candidate.snapshot is characterization.snapshot
    assert candidate.characterization.case_name == "figure5"
    assert candidate.characterization.tool_diameter == 2.0
    assert candidate.characterization.tea_cap == 80.0


def test_direct_candidate_construction_is_disabled() -> None:
    with pytest.raises(TypeError):
        HeldPostQualificationCandidate()


def test_functional_helper_delegates_to_candidate_factory(monkeypatch: pytest.MonkeyPatch) -> None:
    characterization = _characterization()
    sentinel = cast(HeldPostQualificationCandidate, object())
    calls: list[HeldFigure5Characterization] = []

    def build(value: HeldFigure5Characterization) -> HeldPostQualificationCandidate:
        calls.append(value)
        return sentinel

    monkeypatch.setattr(post_qualification_module.HeldPostQualificationCandidate, "build", build)

    assert require_post_qualification_candidate(characterization) is sentinel
    assert calls == [characterization]


def test_candidate_factory_consumes_failure_collector_once(monkeypatch: pytest.MonkeyPatch) -> None:
    characterization = _characterization()
    calls: list[HeldFigure5Characterization] = []

    def failures(value: HeldFigure5Characterization) -> tuple[str, ...]:
        calls.append(value)
        return ()

    monkeypatch.setattr(post_qualification_module, "post_qualification_failures", failures)

    candidate = HeldPostQualificationCandidate.build(characterization)
    assert candidate.characterization is characterization
    assert calls == [characterization]
