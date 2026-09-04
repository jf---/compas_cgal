"""Validated immutable evidence boundary for Held--Pfeiffer Figure 5."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Type
from typing import TypeVar
from typing import cast

from typing_extensions import Self

from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import ContradictoryEngagementEvidenceError
from benchmarks.errors import ContradictoryPathQualityEvidenceError
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import UnexpectedHeldPathCaseError
from benchmarks.held_path_snapshot import HeldArcSnapshot
from benchmarks.held_path_snapshot import HeldCircleSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.quality import PathQuality
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import QualityEvidence
from benchmarks.survey import EngagementSample
from benchmarks.survey import MotionKind
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.survey import RapidMotion
from benchmarks.units import Degrees
from benchmarks.units import OperationIndex
from benchmarks.units import Seconds
from benchmarks.units import SquareMillimetre
from benchmarks.units import UnitFraction
from benchmarks.units import closed_unit_fraction
from benchmarks.units import operation_index
from benchmarks.units import seconds_value
from benchmarks.units import square_millimetres_value
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement import EngagementReport
from compas_cgal.engagement import OperationEngagement

RecordT = TypeVar("RecordT")
ZERO_COUNT = 0


def _build_record(record_type: Type[RecordT], values: dict[str, object]) -> RecordT:
    record = object.__new__(record_type)
    for name, value in values.items():
        object.__setattr__(record, name, value)
    return record


def _ordered_indices(values: object, operation_count: int, *, name: str) -> tuple[OperationIndex, ...]:
    if type(values) is not tuple:
        raise InvalidHeldPathEvidenceError(f"{name} must be one immutable tuple.")
    checked = tuple(operation_index(value, operation_count=operation_count) for value in cast(tuple[int, ...], values))
    if checked != tuple(sorted(checked)) or len(set(checked)) != len(checked):
        raise InvalidHeldPathEvidenceError(f"{name} must contain unique indices in source order.")
    return checked


def _snapshot_kind(operation: HeldOperationSnapshot) -> MotionKind:
    if type(operation) is HeldLineSnapshot:
        return MotionKind.LINE
    if type(operation) is HeldArcSnapshot:
        return MotionKind.ARC
    if type(operation) is HeldCircleSnapshot:
        return MotionKind.LOOP
    raise InvalidHeldPathEvidenceError("snapshot contains an unsupported operation record.")


@dataclass(frozen=True, init=False)
class EngagementExceedanceWitness:
    """One retained exact-predicate exceedance at a world-XY position."""

    operation_index: OperationIndex
    position: Point2[WorldXY]

    def __init__(self) -> None:
        raise TypeError("EngagementExceedanceWitness must be created with EngagementExceedanceWitness.build().")

    @classmethod
    def build(cls, *, operation_index_: OperationIndex, position: Point2[WorldXY], operation_count: int) -> Self:
        checked_index = operation_index(operation_index_, operation_count=operation_count)
        if type(position) is not Point2:
            raise InvalidHeldPathEvidenceError("witness position must be one typed world-XY point.")
        try:
            checked_position = Point2[WorldXY].build(position.x, position.y)
        except (TypeError, ValueError, OverflowError) as error:
            raise InvalidHeldPathEvidenceError("witness position must contain two finite world-XY coordinates.") from error
        return _build_record(cls, {"operation_index": checked_index, "position": checked_position})


@dataclass(frozen=True, init=False)
class EngagementDispositionCounts:
    """Exhaustive guarded-audit disposition indices."""

    tea_audited: tuple[OperationIndex, ...]
    excluded: tuple[OperationIndex, ...]
    certified: tuple[OperationIndex, ...]
    demonstrated_exceeded: tuple[OperationIndex, ...]
    unresolved: tuple[OperationIndex, ...]
    tea_audited_count: int
    excluded_count: int
    certified_count: int
    demonstrated_exceeded_count: int
    unresolved_count: int

    def __init__(self) -> None:
        raise TypeError("EngagementDispositionCounts must be created with EngagementDispositionCounts.build().")

    @classmethod
    def build(
        cls,
        *,
        operation_count: int,
        tea_audited: frozenset[OperationIndex],
        excluded: frozenset[OperationIndex],
        certified: frozenset[OperationIndex],
        demonstrated_exceeded: frozenset[OperationIndex],
        unresolved: frozenset[OperationIndex],
    ) -> Self:
        if type(operation_count) is not int or operation_count <= ZERO_COUNT:
            raise InvalidHeldPathEvidenceError("engagement dispositions require a positive operation count.")
        named = {
            "tea_audited": tea_audited,
            "excluded": excluded,
            "certified": certified,
            "demonstrated_exceeded": demonstrated_exceeded,
            "unresolved": unresolved,
        }
        checked: dict[str, frozenset[OperationIndex]] = {}
        for name, values in named.items():
            if type(values) is not frozenset:
                raise InvalidHeldPathEvidenceError(f"{name} must be an immutable operation-index set.")
            checked[name] = frozenset(operation_index(value, operation_count=operation_count) for value in values)
        all_indices = frozenset(OperationIndex(index) for index in range(operation_count))
        audited = checked["tea_audited"]
        excluded_set = checked["excluded"]
        certified_set = checked["certified"]
        demonstrated = checked["demonstrated_exceeded"]
        unresolved_set = checked["unresolved"]
        if audited & excluded_set or audited | excluded_set != all_indices:
            raise InvalidHeldPathEvidenceError("TEA-audited and excluded indices must exactly partition every source operation.")
        if any(left & right for left, right in ((certified_set, demonstrated), (certified_set, unresolved_set), (demonstrated, unresolved_set))):
            raise InvalidHeldPathEvidenceError("certified, demonstrated, and unresolved dispositions must be disjoint.")
        if certified_set | demonstrated | unresolved_set != audited:
            raise InvalidHeldPathEvidenceError("engagement dispositions must exactly partition TEA-audited operations.")
        record_values: dict[str, object] = {name: tuple(sorted(indices)) for name, indices in checked.items()}
        record_values.update({f"{name}_count": len(indices) for name, indices in checked.items()})
        return _build_record(cls, record_values)


@dataclass(frozen=True, init=False)
class HeldFigure5Characterization:
    """Complete immutable Figure 5 evidence consumed by reporting and gates."""

    case_name: str
    tool_diameter: Millimetre
    tea_cap: Degrees
    snapshot: tuple[HeldOperationSnapshot, ...]
    engagement: EngagementDispositionCounts
    witnesses: tuple[EngagementExceedanceWitness, ...]
    source_operation_count: int
    tea_audited_operation_count: int
    excluded_operation_count: int
    sampled_material_contact_operations: int
    generation_seconds: Seconds
    audit_seconds: Seconds
    survey_seconds: Seconds
    reduction_seconds: Seconds
    path_quality: PathQuality
    assessment: PathQualityAssessment
    coverage_nx: int
    coverage_ny: int
    coverage_cell_area: SquareMillimetre
    coverage_reachable_samples: int
    coverage_uncut_reachable_samples: int
    coverage_remaining_samples: int
    coverage_wall_scallop_height: Millimetre
    coverage_uncut_fraction: UnitFraction
    coverage_remaining_area: SquareMillimetre

    def __init__(self) -> None:
        raise TypeError("HeldFigure5Characterization must be created with HeldFigure5Characterization.build().")

    @classmethod
    def build(
        cls,
        *,
        case: HeldReferenceCase,
        snapshot: tuple[HeldOperationSnapshot, ...],
        audit: EngagementReport,
        survey: PathSurvey,
        quality: QualityEvidence,
        generation_seconds: Seconds,
        audit_seconds: Seconds,
        survey_seconds: Seconds,
        reduction_seconds: Seconds,
    ) -> Self:
        canonical_case = load_held_reference_case("figure5")
        if type(case) is not HeldReferenceCase or case != canonical_case:
            raise UnexpectedHeldPathCaseError("Held path characterization requires the exact committed Figure 5 case.")
        operation_count = _validate_snapshot(snapshot)
        timings = {
            "generation_seconds": _timing(generation_seconds, name="generation duration"),
            "audit_seconds": _timing(audit_seconds, name="guarded-audit duration"),
            "survey_seconds": _timing(survey_seconds, name="survey duration"),
            "reduction_seconds": _timing(reduction_seconds, name="quality-reduction duration"),
        }
        _validate_audit(case, snapshot, audit)
        motion_indices, excluded_indices = _validate_survey(case, snapshot, survey)
        tea_audited = frozenset(OperationIndex(row.op_index) for row in audit.operations if row.stations > ZERO_COUNT)
        if tea_audited != motion_indices:
            raise ContradictoryEngagementEvidenceError("TEA-audited operation indices disagree with survey motion indices.")
        witnesses = _witnesses(survey, operation_count)
        certified = frozenset(OperationIndex(row.op_index) for row in audit.operations if row.stations > ZERO_COUNT and row.cap_certified)
        demonstrated = frozenset(witness.operation_index for witness in witnesses)
        if certified & demonstrated:
            raise ContradictoryEngagementEvidenceError("A certified operation has a sampled exact-predicate exceedance witness.")
        unresolved = tea_audited - certified - demonstrated
        excluded = frozenset(OperationIndex(index) for index in range(operation_count)) - tea_audited
        if excluded != excluded_indices:
            raise InvalidHeldPathEvidenceError("TEA-audit-excluded indices disagree with the complete survey partition.")
        if len(excluded) != len(survey.rapids) + survey.plunges:
            raise InvalidHeldPathEvidenceError("TEA-audit-excluded count must equal survey rapid and plunge counts.")
        if type(audit.cap_violations) is not int or audit.cap_violations != len(demonstrated | unresolved):
            raise ContradictoryEngagementEvidenceError("audit cap-violation count disagrees with non-certified dispositions.")
        engagement = EngagementDispositionCounts.build(
            operation_count=operation_count,
            tea_audited=tea_audited,
            excluded=excluded,
            certified=certified,
            demonstrated_exceeded=demonstrated,
            unresolved=unresolved,
        )
        path_quality, assessment, coverage = _validate_quality(quality, operation_count, survey)
        sampled_contact = sum(motion.is_engaged for motion in survey.motions)
        values: dict[str, object] = {
            "case_name": case.name,
            "tool_diameter": Millimetre(2.0 * float(case.tool_radius.value)),
            "tea_cap": Degrees(float(case.tea_cap)),
            "snapshot": tuple(snapshot),
            "engagement": engagement,
            "witnesses": witnesses,
            "source_operation_count": operation_count,
            "tea_audited_operation_count": len(engagement.tea_audited),
            "excluded_operation_count": len(engagement.excluded),
            "sampled_material_contact_operations": sampled_contact,
            **timings,
            "path_quality": path_quality,
            "assessment": assessment,
            "coverage_nx": coverage.nx,
            "coverage_ny": coverage.ny,
            "coverage_cell_area": square_millimetres_value(coverage.cell_area, name="coverage cell area"),
            "coverage_reachable_samples": coverage.reachable_samples,
            "coverage_uncut_reachable_samples": coverage.uncut_reachable_samples,
            "coverage_remaining_samples": coverage.remaining_samples,
            "coverage_wall_scallop_height": Millimetre(coverage.wall_scallop_height),
            "coverage_uncut_fraction": closed_unit_fraction(coverage.uncut_fraction, name="coverage uncut fraction"),
            "coverage_remaining_area": square_millimetres_value(coverage.remaining_area, name="coverage remaining area"),
        }
        return _build_record(cls, values)


def _timing(value: Seconds, *, name: str) -> Seconds:
    if type(value) is not float:
        raise InvalidHeldPathEvidenceError(f"{name} must be one typed Seconds value.")
    return seconds_value(value, name=name)


def _validate_snapshot(snapshot: object) -> int:
    if type(snapshot) is not tuple or not snapshot:
        raise InvalidHeldPathEvidenceError("Figure 5 characterization requires a non-empty immutable snapshot.")
    typed = cast(tuple[HeldOperationSnapshot, ...], snapshot)
    for expected, operation in enumerate(typed):
        _snapshot_kind(operation)
        if operation.ordinal != expected:
            raise InvalidHeldPathEvidenceError("snapshot ordinals must cover source operations in order.")
    return len(typed)


def _validate_audit(case: HeldReferenceCase, snapshot: tuple[HeldOperationSnapshot, ...], audit: object) -> None:
    if type(audit) is not EngagementReport or type(audit.operations) is not list:
        raise InvalidHeldPathEvidenceError("guarded audit must be one exact EngagementReport with an operation list.")
    if audit.tool_diameter != 2.0 * float(case.tool_radius.value) or audit.tea_cap != math.radians(float(case.tea_cap)):
        raise ContradictoryEngagementEvidenceError("guarded audit tool or cap disagrees with Figure 5.")
    if len(audit.operations) != len(snapshot):
        raise ContradictoryEngagementEvidenceError("guarded audit must retain one ordered row per source operation.")
    for expected, (row, operation) in enumerate(zip(audit.operations, snapshot)):
        if type(row) is not OperationEngagement or row.op_index != expected or row.operation is not operation.operation:
            raise ContradictoryEngagementEvidenceError("guarded audit row sequence disagrees with the source snapshot.")
        if type(row.stations) is not int or row.stations < ZERO_COUNT or type(row.cap_certified) is not bool:
            raise InvalidHeldPathEvidenceError("guarded audit rows require non-negative station counts and exact certification flags.")


def _validate_survey(
    case: HeldReferenceCase,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: object,
) -> tuple[frozenset[OperationIndex], frozenset[OperationIndex]]:
    if type(survey) is not PathSurvey:
        raise InvalidHeldPathEvidenceError("survey must be one exact PathSurvey.")
    canonical_spec = case.pocket_spec()
    if survey.spec != canonical_spec or type(survey.source_snapshot) is not tuple or survey.source_snapshot != snapshot:
        raise InvalidHeldPathEvidenceError("survey case or source snapshot disagrees with Figure 5 characterization.")
    if type(survey.motions) is not tuple or type(survey.rapids) is not tuple:
        raise InvalidHeldPathEvidenceError("survey motion observations must be immutable tuples.")
    operation_count = len(snapshot)
    motion_indices = _ordered_indices(tuple(motion.index for motion in survey.motions), operation_count, name="survey motions")
    rapid_indices = _ordered_indices(tuple(rapid.index for rapid in survey.rapids), operation_count, name="survey rapids")
    plunge_indices = _ordered_indices(survey.plunge_indices, operation_count, name="survey plunges")
    retract_indices = _ordered_indices(survey.retract_indices, operation_count, name="survey retracts")
    if type(survey.plunges) is not int or survey.plunges != len(plunge_indices):
        raise InvalidHeldPathEvidenceError("survey plunge count disagrees with retained indices.")
    if type(survey.retracts) is not int or survey.retracts != len(retract_indices):
        raise InvalidHeldPathEvidenceError("survey retract count disagrees with retained indices.")
    motion_set = frozenset(motion_indices)
    rapid_set = frozenset(rapid_indices)
    plunge_set = frozenset(plunge_indices)
    retract_set = frozenset(retract_indices)
    if not retract_set <= rapid_set:
        raise InvalidHeldPathEvidenceError("every retract must retain one rapid observation.")
    pure_rapid = rapid_set - retract_set
    categories = (motion_set, pure_rapid, plunge_set, retract_set)
    for position, category in enumerate(categories):
        if any(category & later for later in categories[position + 1 :]):
            raise InvalidHeldPathEvidenceError("survey source-operation categories must be disjoint.")
    all_indices = frozenset(OperationIndex(index) for index in range(operation_count))
    if frozenset().union(*categories) != all_indices:
        raise InvalidHeldPathEvidenceError("survey categories must exactly partition the source snapshot.")
    for motion in survey.motions:
        if type(motion) is not MotionQuality:
            raise InvalidHeldPathEvidenceError("survey motions must contain exact MotionQuality records.")
        source = snapshot[motion.index]
        if motion.operation is not source.operation or motion.kind is not _snapshot_kind(source):
            raise InvalidHeldPathEvidenceError("survey motion role or primitive kind disagrees with the source snapshot.")
        if type(motion.samples) is not tuple:
            raise InvalidHeldPathEvidenceError("survey samples must be retained in immutable tuples.")
    for rapid in survey.rapids:
        if type(rapid) is not RapidMotion:
            raise InvalidHeldPathEvidenceError("survey rapids must contain exact RapidMotion records.")
        source = snapshot[rapid.index]
        if rapid.operation is not source.operation or rapid.kind is not _snapshot_kind(source):
            raise InvalidHeldPathEvidenceError("survey rapid role or primitive kind disagrees with the source snapshot.")
    return motion_set, all_indices - motion_set


def _witnesses(survey: PathSurvey, operation_count: int) -> tuple[EngagementExceedanceWitness, ...]:
    witnesses: list[EngagementExceedanceWitness] = []
    for motion in survey.motions:
        checked_index = operation_index(motion.index, operation_count=operation_count)
        for sample in motion.samples:
            if type(sample) is not EngagementSample or type(sample.cap_exceeded) is not bool:
                raise InvalidHeldPathEvidenceError("survey samples must be exact EngagementSample records with exact verdicts.")
            if sample.cap_exceeded:
                witnesses.append(
                    EngagementExceedanceWitness.build(
                        operation_index_=checked_index,
                        position=sample.position,
                        operation_count=operation_count,
                    )
                )
    return tuple(witnesses)


def _validate_quality(
    quality: object,
    operation_count: int,
    survey: PathSurvey,
) -> tuple[PathQuality, PathQualityAssessment, CoverageEstimate]:
    if type(quality) is not QualityEvidence:
        raise InvalidHeldPathEvidenceError("quality must be one validated QualityEvidence record.")
    if type(quality.path_quality) is not PathQuality or type(quality.assessment) is not PathQualityAssessment:
        raise InvalidHeldPathEvidenceError("quality evidence requires exact aggregate and assessment records.")
    if type(quality.coverage) is not CoverageEstimate:
        raise InvalidHeldPathEvidenceError("quality evidence requires one exact coverage record.")
    assessment = quality.assessment
    canonical_assessment = PathQualityAssessment.build(
        spec=survey.spec,
        snapshot=survey.source_snapshot,
        survey=survey,
        uncut_fraction=assessment.uncut_fraction,
        gouging_motions=assessment.gouging_motions,
        unsafe_rapids=assessment.unsafe_rapids,
        continuity_breaks=assessment.continuity_breaks,
        zero_length_motions=assessment.zero_length_motions,
        degenerate_loops=assessment.degenerate_loops,
        redundant_operations=assessment.redundant_operations,
        cap_exceedances=assessment.cap_exceedances,
        slotting_motions=assessment.slotting_motions,
        max_engagement_step=assessment.max_engagement_step,
        max_loop_radius_step=assessment.max_loop_radius_step,
        tangent_breaks=assessment.tangent_breaks,
        attribution=assessment.attribution,
    )
    if canonical_assessment != assessment:
        raise ContradictoryPathQualityEvidenceError("quality assessment differs from its canonical validated reconstruction.")
    if assessment.attribution.operation_count != operation_count:
        raise InvalidHeldPathEvidenceError("quality attribution operation count disagrees with the source snapshot.")
    path_quality = quality.path_quality
    if path_quality.cut_operations != len(survey.motions) or path_quality.path_length != survey.total_length:
        raise ContradictoryPathQualityEvidenceError("PathQuality source counts or path length disagree with the supplied survey.")
    projections = (
        path_quality.elementary.uncut_fraction,
        path_quality.elementary.gouging_motions,
        path_quality.elementary.unsafe_rapids,
        path_quality.elementary.continuity_breaks,
        path_quality.elementary.zero_length_motions,
        path_quality.elementary.degenerate_loops,
        path_quality.elementary.redundant_operations,
        path_quality.cut.cap_exceedances,
        path_quality.cut.slotting_motions,
        path_quality.cut.max_engagement_step_deg,
        path_quality.cut.max_loop_radius_step,
        path_quality.speed.tangent_breaks,
    )
    measured = tuple(getattr(canonical_assessment, name).measured for name in _CRITERION_FIELDS)
    if projections != measured:
        raise ContradictoryPathQualityEvidenceError("PathQuality values disagree with the twelve assessed criterion values.")
    coverage = quality.coverage
    integer_scalars = (coverage.nx, coverage.ny, coverage.reachable_samples, coverage.uncut_reachable_samples, coverage.remaining_samples)
    if any(type(value) is not int or value < ZERO_COUNT for value in integer_scalars):
        raise InvalidHeldPathEvidenceError("coverage counts must be exact non-negative integers.")
    if coverage.nx == ZERO_COUNT or coverage.ny == ZERO_COUNT or coverage.reachable_samples == ZERO_COUNT:
        raise InvalidHeldPathEvidenceError("coverage dimensions and reachable sample count must be positive.")
    if coverage.uncut_reachable_samples > coverage.reachable_samples:
        raise InvalidHeldPathEvidenceError("uncut reachable samples cannot exceed reachable samples.")
    scalar_values = (
        coverage.cell_area,
        coverage.wall_scallop_height,
    )
    if any(type(value) is not float for value in scalar_values):
        raise InvalidHeldPathEvidenceError("coverage measurements must be exact floating-point values.")
    if any(not math.isfinite(value) for value in scalar_values):
        raise InvalidHeldPathEvidenceError("coverage scalars must be finite.")
    if coverage.cell_area <= 0.0 or coverage.wall_scallop_height < 0.0:
        raise InvalidHeldPathEvidenceError("coverage cell area must be positive and wall scallop height non-negative.")
    if float(assessment.uncut_fraction.measured) != coverage.uncut_fraction:
        raise ContradictoryPathQualityEvidenceError("coverage uncut fraction disagrees with the assessed value.")
    if path_quality.cut.wall_scallop_height != coverage.wall_scallop_height:
        raise ContradictoryPathQualityEvidenceError("coverage wall scallop height disagrees with PathQuality.")
    return path_quality, canonical_assessment, coverage


_CRITERION_FIELDS = (
    "uncut_fraction",
    "gouging_motions",
    "unsafe_rapids",
    "continuity_breaks",
    "zero_length_motions",
    "degenerate_loops",
    "redundant_operations",
    "cap_exceedances",
    "slotting_motions",
    "max_engagement_step",
    "max_loop_radius_step",
    "tangent_breaks",
)
