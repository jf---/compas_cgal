from __future__ import annotations

from dataclasses import replace
from typing import Any
from typing import cast

import numpy as np
import pytest
from compas.geometry import Line

from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import ContradictoryEngagementEvidenceError
from benchmarks.errors import ContradictoryPathQualityEvidenceError
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import UnexpectedHeldPathCaseError
from benchmarks.held_path_evidence import EngagementDispositionCounts
from benchmarks.held_path_evidence import EngagementExceedanceWitness
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.quality import PathQuality
from benchmarks.quality import _cut
from benchmarks.quality import _elementary
from benchmarks.quality import _longevity
from benchmarks.quality import _program
from benchmarks.quality import _speed
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import MeasuredStep
from benchmarks.quality_observations import OperationPair
from benchmarks.quality_observations import _build_record
from benchmarks.quality_observations import assess_path_quality
from benchmarks.survey import EngagementSample
from benchmarks.survey import MotionKind
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.survey import RapidMotion
from benchmarks.units import OperationIndex
from benchmarks.units import MotionCount
from benchmarks.units import Seconds
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.engagement import EngagementReport
from compas_cgal.engagement import OperationEngagement
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult


def _snapshot() -> tuple[HeldOperationSnapshot, ...]:
    roles = (OperationType.PLUNGE, OperationType.CUT, OperationType.CUT, OperationType.RETRACT)
    positive_x = Direction3[WorldXYZ].build(1.0, 0.0, 0.0)
    positive_z = Direction3[WorldXYZ].build(0.0, 0.0, 1.0)
    negative_z = Direction3[WorldXYZ].build(0.0, 0.0, -1.0)
    records: list[HeldOperationSnapshot] = []
    for index, role in enumerate(roles):
        if role is OperationType.PLUNGE:
            start = Point3[WorldXYZ].build(0.0, 0.0, 1.0)
            end = Point3[WorldXYZ].build(0.0, 0.0, 0.0)
            tangent = negative_z
        elif role is OperationType.RETRACT:
            start = Point3[WorldXYZ].build(2.0, 0.0, 0.0)
            end = Point3[WorldXYZ].build(2.0, 0.0, 1.0)
            tangent = positive_z
        else:
            start = Point3[WorldXYZ].build(float(index - 1), 0.0, 0.0)
            end = Point3[WorldXYZ].build(float(index), 0.0, 0.0)
            tangent = positive_x
        records.append(
            HeldLineSnapshot.build(
                ordinal=OperationIndex(index),
                operation=role,
                path_index=0,
                clockwise=False,
                start=start,
                end=end,
                start_tangent=tangent,
                end_tangent=tangent,
            )
        )
    return tuple(records)


def _sample(*, cap: bool = False, engagement: float = 0.0, x: float = 1.0, y: float = 0.0) -> EngagementSample:
    return EngagementSample(0.0, Point2[WorldXY].build(x, y), engagement, cap, True)


def _motion(index: int, samples: tuple[EngagementSample, ...], *, removes: bool = True) -> MotionQuality:
    return MotionQuality(
        index=index,
        operation=OperationType.CUT,
        kind=MotionKind.LINE,
        length=1.0,
        loop_radius=None,
        curvature=0.0,
        swept_area=1.0,
        start=(float(index - 1), 0.0),
        end=(float(index), 0.0),
        start_tangent=(1.0, 0.0),
        end_tangent=(1.0, 0.0),
        samples=samples,
        cap_exceeded=any(sample.cap_exceeded for sample in samples),
        slot_exceeded=False,
        removes_material=removes,
    )


def _survey(snapshot: tuple[HeldOperationSnapshot, ...], second_samples: tuple[EngagementSample, ...] = ()) -> PathSurvey:
    motions = (_motion(1, (_sample(engagement=20.0),)), _motion(2, second_samples, removes=False))
    rapid = RapidMotion(3, OperationType.RETRACT, MotionKind.LINE, 1.0, False)
    return PathSurvey(
        spec=load_held_reference_case("figure5").pocket_spec(),
        source_snapshot=snapshot,
        motions=motions,
        rapids=(rapid,),
        plunges=1,
        retracts=1,
        plunge_indices=(OperationIndex(0),),
        retract_indices=(OperationIndex(3),),
        final_stock=cast(Stock, object()),
        total_length=3.0,
        cut_length=2.0,
        air_length=1.0,
        plunge_swept_area=1.0,
    )


def _quality(snapshot: tuple[HeldOperationSnapshot, ...], survey: PathSurvey) -> QualityEvidence:
    coverage = CoverageEstimate(10, 10, 1.0, 10, 0, 0, 0.0)
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
    return _build_record(QualityEvidence, {"path_quality": path_quality, "assessment": assessment, "coverage": coverage})


def _audit(snapshot: tuple[HeldOperationSnapshot, ...], *, second_certified: bool = True, cap_violations: int = 0) -> EngagementReport:
    operations = [OperationEngagement(index, operation.operation, 0.0, index != 2 or second_certified, 0 if index in {0, 3} else 3) for index, operation in enumerate(snapshot)]
    return EngagementReport(2.0, 1.3962634015954636, operations, 0.0, cap_violations, 0)


def _inputs(*, samples: tuple[EngagementSample, ...] = (), second_certified: bool = True, cap_violations: int = 0) -> dict[str, object]:
    snapshot = _snapshot()
    survey = _survey(snapshot, samples)
    return {
        "case": load_held_reference_case("figure5"),
        "snapshot": snapshot,
        "audit": _audit(snapshot, second_certified=second_certified, cap_violations=cap_violations),
        "survey": survey,
        "quality": _quality(snapshot, survey),
        "generation_seconds": Seconds(1.0),
        "audit_seconds": Seconds(2.0),
        "survey_seconds": Seconds(3.0),
        "reduction_seconds": Seconds(4.0),
    }


def _build(values: dict[str, object]) -> HeldFigure5Characterization:
    return HeldFigure5Characterization.build(**values)  # type: ignore[arg-type]


def test_unmeasured_certified_row_is_excluded_not_certified() -> None:
    result = _build(_inputs())
    assert result.reference_primitive_count == 31
    assert result.projection_vertex_count == 65
    assert result.engagement.tea_audited == (OperationIndex(1), OperationIndex(2))
    assert result.engagement.certified == (OperationIndex(1), OperationIndex(2))
    assert result.engagement.excluded == (OperationIndex(0), OperationIndex(3))


def test_five_witnesses_reduce_to_one_demonstrated_operation() -> None:
    samples = tuple(_sample(cap=True, engagement=10.0, x=float(index)) for index in range(5))
    result = _build(_inputs(samples=samples, second_certified=False, cap_violations=1))
    assert len(result.witnesses) == 5
    assert result.engagement.demonstrated_exceeded == (OperationIndex(2),)


def test_certified_witness_overlap_is_contradictory() -> None:
    with pytest.raises(ContradictoryEngagementEvidenceError):
        _build(_inputs(samples=(_sample(cap=True),), cap_violations=0))


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "reordered", "wrong-kind", "out-of-bounds"])
def test_audit_rows_must_match_snapshot_exactly(mutation: str) -> None:
    values = _inputs()
    audit = cast(EngagementReport, values["audit"])
    if mutation == "missing":
        audit.operations.pop()
    elif mutation == "duplicate":
        audit.operations[2] = audit.operations[1]
    elif mutation == "reordered":
        audit.operations[1], audit.operations[2] = audit.operations[2], audit.operations[1]
    elif mutation == "wrong-kind":
        audit.operations[1] = replace(audit.operations[1], operation=OperationType.LINK)
    else:
        audit.operations[1] = replace(audit.operations[1], op_index=9)
    with pytest.raises(ContradictoryEngagementEvidenceError):
        _build(values)


def test_audited_indices_must_equal_survey_motion_indices() -> None:
    values = _inputs()
    audit = cast(EngagementReport, values["audit"])
    audit.operations[1] = replace(audit.operations[1], stations=0)
    with pytest.raises(ContradictoryEngagementEvidenceError):
        _build(values)


def test_plunge_and_repeated_cleared_cut_have_distinct_semantics() -> None:
    result = _build(_inputs())
    assert OperationIndex(0) in result.engagement.excluded
    assert OperationIndex(2) in result.engagement.tea_audited
    assert result.sampled_material_contact_operations == 1


def test_reporting_values_do_not_decide_disposition() -> None:
    baseline = _inputs()
    changed = _inputs()
    audit = cast(EngagementReport, changed["audit"])
    audit.operations[1] = replace(audit.operations[1], max_tea=999.0)
    survey = cast(PathSurvey, changed["survey"])
    object.__setattr__(survey.motions[0].samples[0], "engagement_deg", -999.0)
    changed["quality"] = _quality(cast(tuple[HeldOperationSnapshot, ...], changed["snapshot"]), survey)
    assert _build(baseline).engagement == _build(changed).engagement


def test_requires_exact_canonical_figure5_case() -> None:
    values = _inputs()
    values["case"] = replace(load_held_reference_case("figure5"), title="not canonical")
    with pytest.raises(UnexpectedHeldPathCaseError):
        _build(values)


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf"), True, 1])
def test_requires_non_negative_finite_typed_timings(value: object) -> None:
    values = _inputs()
    values["audit_seconds"] = value
    with pytest.raises(InvalidHeldPathEvidenceError):
        _build(values)


def test_witness_coordinates_may_be_outside_domain_but_must_be_finite() -> None:
    outside = _inputs(samples=(_sample(cap=True, x=-1_000_000.0, y=1_000_000.0),), second_certified=False, cap_violations=1)
    assert _build(outside).witnesses[0].position == Point2[WorldXY].build(-1_000_000.0, 1_000_000.0)

    invalid = _inputs(samples=(_sample(cap=True),), second_certified=False, cap_violations=1)
    survey = cast(PathSurvey, invalid["survey"])
    object.__setattr__(survey.motions[1].samples[0].position, "x", float("nan"))
    with pytest.raises(InvalidHeldPathEvidenceError):
        _build(invalid)


def test_path_quality_must_equal_all_assessment_values() -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    changed_elementary = replace(quality.path_quality.elementary, redundant_operations=99)
    changed_quality = replace(quality.path_quality, elementary=changed_elementary)
    values["quality"] = _build_record(QualityEvidence, {"path_quality": changed_quality, "assessment": quality.assessment, "coverage": quality.coverage})
    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _build(values)


@pytest.mark.parametrize(
    ("group", "field"),
    [
        ("elementary", "uncut_fraction"),
        ("elementary", "gouging_motions"),
        ("elementary", "unsafe_rapids"),
        ("elementary", "continuity_breaks"),
        ("elementary", "zero_length_motions"),
        ("elementary", "degenerate_loops"),
        ("elementary", "redundant_operations"),
        ("cut", "cap_exceedances"),
        ("cut", "slotting_motions"),
        ("cut", "max_engagement_step_deg"),
        ("cut", "max_loop_radius_step"),
        ("speed", "tangent_breaks"),
    ],
)
def test_each_path_quality_projection_is_bound_to_the_assessment(group: str, field: str) -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    group_value = getattr(quality.path_quality, group)
    changed_group = replace(group_value, **{field: getattr(group_value, field) + 1})
    changed_quality = replace(quality.path_quality, **{group: changed_group})
    values["quality"] = _build_record(
        QualityEvidence,
        {"path_quality": changed_quality, "assessment": quality.assessment, "coverage": quality.coverage},
    )
    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _build(values)


@pytest.mark.parametrize("mutation", ["duplicate-name", "wrong-evidence"])
def test_criterion_vocabulary_is_closed(mutation: str) -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    if mutation == "duplicate-name":
        object.__setattr__(quality.assessment.gouging_motions, "name", "uncut fraction")
    else:
        object.__setattr__(quality.assessment.gouging_motions, "evidence", "derived_geometry")
    with pytest.raises(InvalidHeldPathEvidenceError):
        _build(values)


@pytest.mark.parametrize("field", ["required", "outcome"])
def test_assessment_factory_rejects_forged_criterion_contract(field: str) -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    forged: object = MotionCount(99) if field == "required" else "failure_observed"
    object.__setattr__(quality.assessment.gouging_motions, field, forged)
    with pytest.raises(InvalidHeldPathEvidenceError):
        _build(values)


def test_assessment_factory_rejects_forged_pair_attribution() -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    pair = OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=4)
    object.__setattr__(quality.assessment.attribution, "continuity_break_pairs", (pair,))
    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _build(values)


def test_assessment_factory_rejects_forged_extremum_attribution() -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    pair = OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=4)
    forged = MeasuredStep.build(value=99.0, pair=pair, unit="degrees")
    object.__setattr__(quality.assessment.attribution, "max_engagement_step", forged)
    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _build(values)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("nx", 10.0),
        ("reachable_samples", True),
        ("cell_area", "bad"),
        ("wall_scallop_height", None),
    ],
)
def test_malformed_coverage_fields_raise_named_error(field: str, value: object) -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    malformed = replace(quality.coverage, **{field: cast(Any, value)})
    values["quality"] = _build_record(
        QualityEvidence,
        {"path_quality": quality.path_quality, "assessment": quality.assessment, "coverage": malformed},
    )
    with pytest.raises(InvalidHeldPathEvidenceError):
        _build(values)


def test_coverage_wall_scallop_must_match_path_quality() -> None:
    values = _inputs()
    quality = cast(QualityEvidence, values["quality"])
    forged_coverage = replace(quality.coverage, wall_scallop_height=1.0)
    values["quality"] = _build_record(
        QualityEvidence,
        {"path_quality": quality.path_quality, "assessment": quality.assessment, "coverage": forged_coverage},
    )
    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _build(values)


@pytest.mark.parametrize("mutation", ["missing", "overlap", "out-of-bounds", "reordered", "wrong-role", "wrong-source"])
def test_survey_requires_complete_ordered_source_partition(mutation: str) -> None:
    values = _inputs()
    survey = cast(PathSurvey, values["survey"])
    if mutation == "missing":
        object.__setattr__(survey, "motions", survey.motions[:1])
    elif mutation == "overlap":
        object.__setattr__(survey, "plunge_indices", (OperationIndex(0), OperationIndex(1)))
        object.__setattr__(survey, "plunges", 2)
    elif mutation == "out-of-bounds":
        object.__setattr__(survey, "plunge_indices", (OperationIndex(9),))
    elif mutation == "reordered":
        object.__setattr__(survey, "motions", tuple(reversed(survey.motions)))
    elif mutation == "wrong-role":
        object.__setattr__(survey, "motions", (replace(survey.motions[0], operation=OperationType.LINK), survey.motions[1]))
    else:
        object.__setattr__(survey, "source_snapshot", survey.source_snapshot[:-1])
    with pytest.raises(InvalidHeldPathEvidenceError):
        _build(values)


def test_witness_factory_validates_bounds_and_copies_typed_position() -> None:
    point = Point2[WorldXY].build(-100.0, 100.0)
    first = EngagementExceedanceWitness.build(operation_index_=OperationIndex(0), position=point, operation_count=4)
    last = EngagementExceedanceWitness.build(operation_index_=OperationIndex(3), position=point, operation_count=4)
    assert (first.operation_index, last.operation_index) == (OperationIndex(0), OperationIndex(3))
    assert first.position == point and first.position is not point
    with pytest.raises(InvalidHeldPathEvidenceError):
        EngagementExceedanceWitness.build(operation_index_=OperationIndex(4), position=point, operation_count=4)
    object.__setattr__(point, "x", float("nan"))
    with pytest.raises(InvalidHeldPathEvidenceError):
        EngagementExceedanceWitness.build(operation_index_=OperationIndex(0), position=point, operation_count=4)


@pytest.mark.parametrize("mutation", ["overlap", "non-exhaustive"])
def test_disposition_factory_rejects_invalid_partitions(mutation: str) -> None:
    kwargs = {
        "operation_count": 3,
        "tea_audited": frozenset({OperationIndex(1), OperationIndex(2)}),
        "excluded": frozenset({OperationIndex(0)}),
        "certified": frozenset({OperationIndex(1)}),
        "demonstrated_exceeded": frozenset(),
        "unresolved": frozenset({OperationIndex(2)}),
    }
    if mutation == "overlap":
        kwargs["demonstrated_exceeded"] = frozenset({OperationIndex(1)})
    else:
        kwargs["unresolved"] = frozenset()
    with pytest.raises(InvalidHeldPathEvidenceError):
        EngagementDispositionCounts.build(**kwargs)  # type: ignore[arg-type]


def test_cap_violation_aggregate_must_match_noncertified_dispositions() -> None:
    with pytest.raises(ContradictoryEngagementEvidenceError):
        _build(_inputs(second_certified=False, cap_violations=0))


def test_inputs_are_not_retained_or_aliased() -> None:
    values = _inputs()
    result = _build(values)
    before = repr(result)
    audit = cast(EngagementReport, values["audit"])
    audit.operations.clear()
    survey = cast(PathSurvey, values["survey"])
    object.__setattr__(survey, "final_stock", object())
    assert repr(result) == before
    assert not hasattr(result, "audit")
    assert not hasattr(result, "survey")
    assert not hasattr(result, "case")


def test_source_operation_arrays_can_change_without_mutating_characterization() -> None:
    tangent = np.array([1.0, 0.0, 0.0])
    operations = [
        ToolpathOperation(Line([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]), OperationType.PLUNGE, 0),
        ToolpathOperation(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), OperationType.CUT, 0, start_tangent=tangent, end_tangent=tangent),
        ToolpathOperation(Line([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]), OperationType.CUT, 0, start_tangent=tangent, end_tangent=tangent),
        ToolpathOperation(Line([2.0, 0.0, 0.0], [2.0, 0.0, 1.0]), OperationType.RETRACT, 0),
    ]
    source_polyline = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    source = ToolpathResult(operations=operations, polyline=source_polyline)
    snapshot = snapshot_toolpath(source)
    survey = _survey(snapshot)
    values = _inputs()
    values["snapshot"] = snapshot
    values["audit"] = _audit(snapshot)
    values["survey"] = survey
    values["quality"] = _quality(snapshot, survey)
    characterization = _build(values)
    before = characterization.snapshot

    tangent[0] = -1.0
    source_polyline[0, 0] = 99.0
    cast(Line, source.operations[1].geometry).start.x = 99.0

    assert characterization.snapshot == before
    assert characterization.snapshot[1].start == Point3[WorldXYZ].build(0.0, 0.0, 0.0)


def test_direct_construction_is_disabled() -> None:
    with pytest.raises(TypeError):
        EngagementExceedanceWitness()
    with pytest.raises(TypeError):
        EngagementDispositionCounts()
    with pytest.raises(TypeError):
        HeldFigure5Characterization()
