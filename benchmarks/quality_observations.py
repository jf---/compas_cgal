"""Typed, operation-attributed observations for the existing path-quality gate.

This module is additive.  It reproduces the twelve decisions currently made in
``benchmarks.quality`` but does not redirect that judge; parity and convergence
remain separate review gates.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Generic
from typing import Literal
from typing import Optional
from typing import Sequence
from typing import Type
from typing import TypeVar
from typing import Union
from typing import cast
from typing import overload

from typing_extensions import Self
from typing_extensions import TypeAlias

from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import ContradictoryPathQualityEvidenceError
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import UnreplayableOperationError
from benchmarks.held_path_snapshot import HeldArcSnapshot
from benchmarks.held_path_snapshot import HeldCircleSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality import CONTINUITY_TOOL_RADIUS_FRACTION
from benchmarks.quality import DEGENERATE_LOOP_RATIO
from benchmarks.quality import TANGENT_CONTINUITY_SLACK
from benchmarks.spec import PocketSpec
from benchmarks.survey import MotionKind
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.units import Degrees
from benchmarks.units import MotionCount
from benchmarks.units import OperationIndex
from benchmarks.units import ToolRadiusMultiple
from benchmarks.units import UnitFraction
from benchmarks.units import closed_unit_fraction
from benchmarks.units import degrees_value
from benchmarks.units import motion_count
from benchmarks.units import operation_index
from benchmarks.units import tool_radius_multiple
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.replay_classification import CutPlaneRampError
from compas_cgal.replay_classification import OffPlaneReplayCurveError
from compas_cgal.replay_classification import ReplayCategory
from compas_cgal.replay_classification import classify_line_replay
from compas_cgal.replay_classification import classify_planar_replay
from compas_cgal.replay_classification import line_cut_height_anchor
from compas_cgal.replay_classification import minimum_cut_height
from compas_cgal.toolpath import OperationType

CriterionName: TypeAlias = Literal[
    "uncut fraction",
    "gouging motions",
    "unsafe rapids",
    "continuity breaks",
    "zero-length motions",
    "degenerate loops",
    "redundant operations",
    "cap exceedances",
    "slotting motions",
    "max engagement step (deg)",
    "max loop radius step (tool radii)",
    "tangent breaks",
]
EvidenceKind: TypeAlias = Literal[
    "guarded_replay_certificate",
    "sampled_exact_predicate",
    "exact_depletion_replay",
    "sampled_diagnostic",
    "derived_geometry",
    "tolerance_diagnostic",
]
ObservationOutcome: TypeAlias = Literal[
    "no_failure_observed",
    "failure_observed",
    "criterion_satisfied",
    "criterion_violated",
    "within_declared_tolerance",
    "outside_declared_tolerance",
]

CRITERION_NAMES: tuple[CriterionName, ...] = (
    "uncut fraction",
    "gouging motions",
    "unsafe rapids",
    "continuity breaks",
    "zero-length motions",
    "degenerate loops",
    "redundant operations",
    "cap exceedances",
    "slotting motions",
    "max engagement step (deg)",
    "max loop radius step (tool radii)",
    "tangent breaks",
)

EVIDENCE_BY_CRITERION: dict[CriterionName, EvidenceKind] = {
    "uncut fraction": "sampled_diagnostic",
    "gouging motions": "sampled_diagnostic",
    "unsafe rapids": "tolerance_diagnostic",
    "continuity breaks": "tolerance_diagnostic",
    "zero-length motions": "derived_geometry",
    "degenerate loops": "derived_geometry",
    "redundant operations": "exact_depletion_replay",
    "cap exceedances": "sampled_exact_predicate",
    "slotting motions": "sampled_exact_predicate",
    "max engagement step (deg)": "sampled_diagnostic",
    "max loop radius step (tool radii)": "derived_geometry",
    "tangent breaks": "tolerance_diagnostic",
}

CountCriterionName: TypeAlias = Literal[
    "gouging motions",
    "unsafe rapids",
    "continuity breaks",
    "zero-length motions",
    "degenerate loops",
    "redundant operations",
    "cap exceedances",
    "slotting motions",
    "tangent breaks",
]

ZERO_FLOAT = 0.0
ZERO_COUNT = 0
UNIT_DOT = 1.0
NEXT_OPERATION_OFFSET = 1
REQUIRED_FRACTION = UnitFraction(ZERO_FLOAT)
REQUIRED_COUNT = MotionCount(ZERO_COUNT)
REQUIRED_LOOP_STEP = ToolRadiusMultiple(2.0)

FRACTION_NAMES = frozenset({"uncut fraction"})
COUNT_NAMES = frozenset(
    {
        "gouging motions",
        "unsafe rapids",
        "continuity breaks",
        "zero-length motions",
        "degenerate loops",
        "redundant operations",
        "cap exceedances",
        "slotting motions",
        "tangent breaks",
    }
)
DEGREES_NAMES = frozenset({"max engagement step (deg)"})
TOOL_RADIUS_MULTIPLE_NAMES = frozenset({"max loop radius step (tool radii)"})
SAMPLED_EVIDENCE = frozenset({"sampled_diagnostic", "sampled_exact_predicate"})
EXACT_OR_DERIVED_EVIDENCE = frozenset({"exact_depletion_replay", "derived_geometry"})

CriterionT = TypeVar("CriterionT")
StepUnitT = TypeVar("StepUnitT", Degrees, ToolRadiusMultiple)


def _build_record(record_type: Type[CriterionT], values: dict[str, object]) -> CriterionT:
    record = object.__new__(record_type)
    for name, value in values.items():
        object.__setattr__(record, name, value)
    return record


def _criterion_name(value: object, allowed: frozenset[str]) -> CriterionName:
    if not isinstance(value, str) or value not in allowed:
        raise InvalidHeldPathEvidenceError(f"criterion name {value!r} does not belong to this unit domain.")
    return cast(CriterionName, value)


def _evidence(name: CriterionName, value: object) -> EvidenceKind:
    expected = EVIDENCE_BY_CRITERION[name]
    if value != expected:
        raise InvalidHeldPathEvidenceError(f"{name} requires evidence kind {expected!r}, got {value!r}.")
    return expected


def _outcome(evidence: EvidenceKind, satisfied: bool) -> ObservationOutcome:
    if evidence in SAMPLED_EVIDENCE:
        return "no_failure_observed" if satisfied else "failure_observed"
    if evidence in EXACT_OR_DERIVED_EVIDENCE:
        return "criterion_satisfied" if satisfied else "criterion_violated"
    if evidence == "tolerance_diagnostic":
        return "within_declared_tolerance" if satisfied else "outside_declared_tolerance"
    raise InvalidHeldPathEvidenceError("guarded replay evidence is not valid for a path-quality criterion.")


def _validated_outcome(evidence: EvidenceKind, measured: float, required: float, value: object) -> ObservationOutcome:
    expected = _outcome(evidence, measured <= required)
    if value != expected:
        raise InvalidHeldPathEvidenceError(f"outcome {value!r} contradicts the measured and required values; expected {expected!r}.")
    return expected


@dataclass(frozen=True, init=False)
class FractionCriterion:
    name: CriterionName
    measured: UnitFraction
    required: UnitFraction
    evidence: EvidenceKind
    outcome: ObservationOutcome

    def __init__(self) -> None:
        raise TypeError("FractionCriterion must be created with FractionCriterion.build().")

    @classmethod
    def build(cls, *, name: CriterionName, measured: UnitFraction, required: UnitFraction, evidence: EvidenceKind, outcome: ObservationOutcome) -> Self:
        checked_name = _criterion_name(name, FRACTION_NAMES)
        checked_measured = closed_unit_fraction(measured, name=f"{checked_name} measured")
        checked_required = closed_unit_fraction(required, name=f"{checked_name} required")
        if checked_required != REQUIRED_FRACTION:
            raise InvalidHeldPathEvidenceError("uncut fraction requires an exact zero threshold.")
        checked_evidence = _evidence(checked_name, evidence)
        checked_outcome = _validated_outcome(checked_evidence, checked_measured, checked_required, outcome)
        return _build_record(cls, {"name": checked_name, "measured": checked_measured, "required": checked_required, "evidence": checked_evidence, "outcome": checked_outcome})


@dataclass(frozen=True, init=False)
class CountCriterion:
    name: CriterionName
    measured: MotionCount
    required: MotionCount
    evidence: EvidenceKind
    outcome: ObservationOutcome

    def __init__(self) -> None:
        raise TypeError("CountCriterion must be created with CountCriterion.build().")

    @classmethod
    def build(
        cls,
        *,
        name: CountCriterionName,
        measured: MotionCount,
        required: MotionCount,
        evidence: EvidenceKind,
        outcome: ObservationOutcome,
        attribution_count: int,
    ) -> Self:
        checked_name = _criterion_name(name, COUNT_NAMES)
        checked_measured = motion_count(measured, name=f"{checked_name} measured")
        checked_required = motion_count(required, name=f"{checked_name} required")
        checked_attribution = motion_count(attribution_count, name=f"{checked_name} attribution count")
        if checked_required != REQUIRED_COUNT:
            raise InvalidHeldPathEvidenceError(f"{checked_name} requires an exact zero threshold.")
        if checked_measured != checked_attribution:
            raise InvalidHeldPathEvidenceError(f"{checked_name} measured count disagrees with its attribution cardinality.")
        checked_evidence = _evidence(checked_name, evidence)
        checked_outcome = _validated_outcome(checked_evidence, checked_measured, checked_required, outcome)
        return _build_record(cls, {"name": checked_name, "measured": checked_measured, "required": checked_required, "evidence": checked_evidence, "outcome": checked_outcome})


@dataclass(frozen=True, init=False)
class DegreesCriterion:
    name: CriterionName
    measured: Degrees
    required: Degrees
    evidence: EvidenceKind
    outcome: ObservationOutcome

    def __init__(self) -> None:
        raise TypeError("DegreesCriterion must be created with DegreesCriterion.build().")

    @classmethod
    def build(cls, *, name: CriterionName, measured: Degrees, required: Degrees, evidence: EvidenceKind, outcome: ObservationOutcome) -> Self:
        checked_name = _criterion_name(name, DEGREES_NAMES)
        checked_measured = degrees_value(measured, name=f"{checked_name} measured")
        checked_required = degrees_value(required, name=f"{checked_name} required")
        if checked_measured < ZERO_FLOAT or checked_required <= ZERO_FLOAT:
            raise InvalidHeldPathEvidenceError("engagement steps must be non-negative and their case-cap threshold positive.")
        checked_evidence = _evidence(checked_name, evidence)
        checked_outcome = _validated_outcome(checked_evidence, checked_measured, checked_required, outcome)
        return _build_record(cls, {"name": checked_name, "measured": checked_measured, "required": checked_required, "evidence": checked_evidence, "outcome": checked_outcome})


@dataclass(frozen=True, init=False)
class ToolRadiusMultipleCriterion:
    name: CriterionName
    measured: ToolRadiusMultiple
    required: ToolRadiusMultiple
    evidence: EvidenceKind
    outcome: ObservationOutcome

    def __init__(self) -> None:
        raise TypeError("ToolRadiusMultipleCriterion must be created with ToolRadiusMultipleCriterion.build().")

    @classmethod
    def build(
        cls,
        *,
        name: CriterionName,
        measured: ToolRadiusMultiple,
        required: ToolRadiusMultiple,
        evidence: EvidenceKind,
        outcome: ObservationOutcome,
    ) -> Self:
        checked_name = _criterion_name(name, TOOL_RADIUS_MULTIPLE_NAMES)
        checked_measured = tool_radius_multiple(measured, name=f"{checked_name} measured")
        checked_required = tool_radius_multiple(required, name=f"{checked_name} required")
        if checked_required != REQUIRED_LOOP_STEP:
            raise InvalidHeldPathEvidenceError("loop-radius step requires an exact two-tool-radii threshold.")
        checked_evidence = _evidence(checked_name, evidence)
        checked_outcome = _validated_outcome(checked_evidence, checked_measured, checked_required, outcome)
        return _build_record(cls, {"name": checked_name, "measured": checked_measured, "required": checked_required, "evidence": checked_evidence, "outcome": checked_outcome})


@dataclass(frozen=True, init=False)
class OperationPair:
    previous: OperationIndex
    current: OperationIndex

    def __init__(self) -> None:
        raise TypeError("OperationPair must be created with OperationPair.build().")

    @classmethod
    def build(cls, *, previous: OperationIndex, current: OperationIndex, operation_count: int) -> Self:
        checked_previous = operation_index(previous, operation_count=operation_count)
        checked_current = operation_index(current, operation_count=operation_count)
        if checked_previous >= checked_current:
            raise InvalidHeldPathEvidenceError("an operation pair must follow forward source order.")
        return _build_record(cls, {"previous": checked_previous, "current": checked_current})


@dataclass(frozen=True, init=False)
class MeasuredStep(Generic[StepUnitT]):
    value: StepUnitT
    pair: OperationPair
    unit: Literal["degrees", "tool_radius_multiple"]

    def __init__(self) -> None:
        raise TypeError("MeasuredStep must be created with MeasuredStep.build().")

    @classmethod
    @overload
    def build(cls, *, value: Degrees, pair: OperationPair, unit: Literal["degrees"]) -> MeasuredStep[Degrees]: ...

    @classmethod
    @overload
    def build(
        cls,
        *,
        value: ToolRadiusMultiple,
        pair: OperationPair,
        unit: Literal["tool_radius_multiple"],
    ) -> MeasuredStep[ToolRadiusMultiple]: ...

    @classmethod
    def build(
        cls,
        *,
        value: Union[Degrees, ToolRadiusMultiple],
        pair: OperationPair,
        unit: Literal["degrees", "tool_radius_multiple"],
    ) -> Union[MeasuredStep[Degrees], MeasuredStep[ToolRadiusMultiple]]:
        if type(pair) is not OperationPair:
            raise InvalidHeldPathEvidenceError("a measured step requires one validated operation pair.")
        if unit == "degrees":
            checked_degrees = degrees_value(value, name="measured engagement step")
            if checked_degrees < ZERO_FLOAT:
                raise InvalidHeldPathEvidenceError("measured engagement step must be non-negative.")
            return _build_record(cls, {"value": checked_degrees, "pair": pair, "unit": unit})
        if unit == "tool_radius_multiple":
            checked_multiple = tool_radius_multiple(value, name="measured loop-radius step")
            return _build_record(cls, {"value": checked_multiple, "pair": pair, "unit": unit})
        raise InvalidHeldPathEvidenceError(f"measured step unit {unit!r} is not supported.")


StepFinding: TypeAlias = Union[MeasuredStep[Degrees], MeasuredStep[ToolRadiusMultiple]]


@dataclass(frozen=True, init=False)
class PathQualityAttribution:
    operation_count: int
    uncut_operations: tuple[OperationIndex, ...]
    gouging_operations: tuple[OperationIndex, ...]
    unsafe_rapid_operations: tuple[OperationIndex, ...]
    continuity_break_pairs: tuple[OperationPair, ...]
    zero_length_operations: tuple[OperationIndex, ...]
    degenerate_loop_operations: tuple[OperationIndex, ...]
    redundant_operations: tuple[OperationIndex, ...]
    cap_exceeded_operations: tuple[OperationIndex, ...]
    slotting_operations: tuple[OperationIndex, ...]
    max_engagement_step: Optional[MeasuredStep[Degrees]]
    engagement_step_failure_pairs: tuple[OperationPair, ...]
    max_loop_radius_step: Optional[MeasuredStep[ToolRadiusMultiple]]
    loop_radius_step_failure_pairs: tuple[OperationPair, ...]
    tangent_break_pairs: tuple[OperationPair, ...]
    curvature_break_pairs: tuple[OperationPair, ...]
    reversal_pairs: tuple[OperationPair, ...]

    def __init__(self) -> None:
        raise TypeError("PathQualityAttribution must be created with PathQualityAttribution.build().")

    @classmethod
    def build(
        cls,
        *,
        operation_count: int,
        uncut_operations: tuple[OperationIndex, ...],
        gouging_operations: tuple[OperationIndex, ...],
        unsafe_rapid_operations: tuple[OperationIndex, ...],
        continuity_break_pairs: tuple[OperationPair, ...],
        zero_length_operations: tuple[OperationIndex, ...],
        degenerate_loop_operations: tuple[OperationIndex, ...],
        redundant_operations: tuple[OperationIndex, ...],
        cap_exceeded_operations: tuple[OperationIndex, ...],
        slotting_operations: tuple[OperationIndex, ...],
        max_engagement_step: Optional[MeasuredStep[Degrees]],
        engagement_step_failure_pairs: tuple[OperationPair, ...],
        max_loop_radius_step: Optional[MeasuredStep[ToolRadiusMultiple]],
        loop_radius_step_failure_pairs: tuple[OperationPair, ...],
        tangent_break_pairs: tuple[OperationPair, ...],
        curvature_break_pairs: tuple[OperationPair, ...],
        reversal_pairs: tuple[OperationPair, ...],
    ) -> Self:
        if type(operation_count) is not int or operation_count < ZERO_COUNT:
            raise InvalidHeldPathEvidenceError("operation count must be an exact non-negative integer.")
        component_values = (
            uncut_operations,
            gouging_operations,
            unsafe_rapid_operations,
            continuity_break_pairs,
            zero_length_operations,
            degenerate_loop_operations,
            redundant_operations,
            cap_exceeded_operations,
            slotting_operations,
            max_engagement_step,
            engagement_step_failure_pairs,
            max_loop_radius_step,
            loop_radius_step_failure_pairs,
            tangent_break_pairs,
            curvature_break_pairs,
            reversal_pairs,
        )
        if operation_count == ZERO_COUNT and any(component_values):
            raise InvalidHeldPathEvidenceError("zero operations require a wholly empty standalone attribution.")
        if uncut_operations:
            raise InvalidHeldPathEvidenceError("uncut stock is spatial aggregate evidence and cannot name operations.")
        index_fields = {
            "uncut_operations": uncut_operations,
            "gouging_operations": gouging_operations,
            "unsafe_rapid_operations": unsafe_rapid_operations,
            "zero_length_operations": zero_length_operations,
            "degenerate_loop_operations": degenerate_loop_operations,
            "redundant_operations": redundant_operations,
            "cap_exceeded_operations": cap_exceeded_operations,
            "slotting_operations": slotting_operations,
        }
        pair_fields = {
            "continuity_break_pairs": continuity_break_pairs,
            "engagement_step_failure_pairs": engagement_step_failure_pairs,
            "loop_radius_step_failure_pairs": loop_radius_step_failure_pairs,
            "tangent_break_pairs": tangent_break_pairs,
            "curvature_break_pairs": curvature_break_pairs,
            "reversal_pairs": reversal_pairs,
        }
        checked: dict[str, object] = {"operation_count": operation_count}
        for name, indices in index_fields.items():
            checked[name] = _indices(indices, operation_count) if indices else ()
        for name, pairs in pair_fields.items():
            checked[name] = _validated_pairs(pairs, operation_count)
        for name, step in (
            ("max_engagement_step", max_engagement_step),
            ("max_loop_radius_step", max_loop_radius_step),
        ):
            if step is not None and type(step) is not MeasuredStep:
                raise InvalidHeldPathEvidenceError(f"{name} must be one validated measured step or None.")
            if step is not None:
                expected_unit = "degrees" if name == "max_engagement_step" else "tool_radius_multiple"
                if getattr(step, "unit", None) != expected_unit:
                    raise InvalidHeldPathEvidenceError(f"{name} requires a measured step in {expected_unit!r}.")
                checked_pair = OperationPair.build(
                    previous=step.pair.previous,
                    current=step.pair.current,
                    operation_count=operation_count,
                )
                if name == "max_engagement_step":
                    checked[name] = MeasuredStep.build(value=cast(Degrees, step.value), pair=checked_pair, unit="degrees")
                else:
                    checked[name] = MeasuredStep.build(
                        value=cast(ToolRadiusMultiple, step.value),
                        pair=checked_pair,
                        unit="tool_radius_multiple",
                    )
            else:
                checked[name] = None
        tangent_pairs = set(tangent_break_pairs)
        reversal_pair_set = set(reversal_pairs)
        curvature_pairs = set(curvature_break_pairs)
        if not reversal_pair_set <= tangent_pairs:
            raise InvalidHeldPathEvidenceError("every reversal pair must also be attributed as a tangent break.")
        if not tangent_pairs.isdisjoint(curvature_pairs):
            raise InvalidHeldPathEvidenceError("tangent-break and curvature-break pairs must be disjoint.")
        return _build_record(cls, checked)


@dataclass(frozen=True, init=False)
class PathQualityAssessment:
    uncut_fraction: FractionCriterion
    gouging_motions: CountCriterion
    unsafe_rapids: CountCriterion
    continuity_breaks: CountCriterion
    zero_length_motions: CountCriterion
    degenerate_loops: CountCriterion
    redundant_operations: CountCriterion
    cap_exceedances: CountCriterion
    slotting_motions: CountCriterion
    max_engagement_step: DegreesCriterion
    max_loop_radius_step: ToolRadiusMultipleCriterion
    tangent_breaks: CountCriterion
    attribution: PathQualityAttribution

    def __init__(self) -> None:
        raise TypeError("PathQualityAssessment must be created with PathQualityAssessment.build().")

    @classmethod
    def build(
        cls,
        *,
        spec: PocketSpec,
        snapshot: tuple[HeldOperationSnapshot, ...],
        survey: PathSurvey,
        uncut_fraction: FractionCriterion,
        gouging_motions: CountCriterion,
        unsafe_rapids: CountCriterion,
        continuity_breaks: CountCriterion,
        zero_length_motions: CountCriterion,
        degenerate_loops: CountCriterion,
        redundant_operations: CountCriterion,
        cap_exceedances: CountCriterion,
        slotting_motions: CountCriterion,
        max_engagement_step: DegreesCriterion,
        max_loop_radius_step: ToolRadiusMultipleCriterion,
        tangent_breaks: CountCriterion,
        attribution: PathQualityAttribution,
    ) -> Self:
        if type(attribution) is not PathQualityAttribution:
            raise InvalidHeldPathEvidenceError("an assessment requires one validated attribution record.")
        _validate_survey_binding(spec, snapshot, survey)
        cardinalities = {
            "gouging_motions": len(attribution.gouging_operations),
            "unsafe_rapids": len(attribution.unsafe_rapid_operations),
            "continuity_breaks": len(attribution.continuity_break_pairs),
            "zero_length_motions": len(attribution.zero_length_operations),
            "degenerate_loops": len(attribution.degenerate_loop_operations),
            "redundant_operations": len(attribution.redundant_operations),
            "cap_exceedances": len(attribution.cap_exceeded_operations),
            "slotting_motions": len(attribution.slotting_operations),
            "tangent_breaks": len(attribution.tangent_break_pairs),
        }
        supplied_counts = {
            "gouging_motions": gouging_motions,
            "unsafe_rapids": unsafe_rapids,
            "continuity_breaks": continuity_breaks,
            "zero_length_motions": zero_length_motions,
            "degenerate_loops": degenerate_loops,
            "redundant_operations": redundant_operations,
            "cap_exceedances": cap_exceedances,
            "slotting_motions": slotting_motions,
            "tangent_breaks": tangent_breaks,
        }
        count_names: dict[str, CountCriterionName] = {
            "gouging_motions": "gouging motions",
            "unsafe_rapids": "unsafe rapids",
            "continuity_breaks": "continuity breaks",
            "zero_length_motions": "zero-length motions",
            "degenerate_loops": "degenerate loops",
            "redundant_operations": "redundant operations",
            "cap_exceedances": "cap exceedances",
            "slotting_motions": "slotting motions",
            "tangent_breaks": "tangent breaks",
        }
        checked_counts: dict[str, CountCriterion] = {}
        for field, criterion in supplied_counts.items():
            if type(criterion) is not CountCriterion:
                raise InvalidHeldPathEvidenceError(f"{field} must be one validated count criterion record.")
            checked_counts[field] = CountCriterion.build(
                name=cast(CountCriterionName, criterion.name),
                measured=criterion.measured,
                required=criterion.required,
                evidence=criterion.evidence,
                outcome=criterion.outcome,
                attribution_count=criterion.measured,
            )
        canonical_counts = {field: _count(count_names[field], cardinalities[field]) for field in supplied_counts}
        if checked_counts != canonical_counts:
            raise ContradictoryPathQualityEvidenceError("count criteria disagree with their canonical thresholds, outcomes, or source attribution.")

        if type(uncut_fraction) is not FractionCriterion:
            raise InvalidHeldPathEvidenceError("uncut_fraction must be one validated fraction criterion record.")
        checked_uncut = FractionCriterion.build(
            name=uncut_fraction.name,
            measured=uncut_fraction.measured,
            required=uncut_fraction.required,
            evidence=uncut_fraction.evidence,
            outcome=uncut_fraction.outcome,
        )
        uncut_evidence = EVIDENCE_BY_CRITERION["uncut fraction"]
        canonical_uncut = FractionCriterion.build(
            name="uncut fraction",
            measured=checked_uncut.measured,
            required=REQUIRED_FRACTION,
            evidence=uncut_evidence,
            outcome=_outcome(uncut_evidence, checked_uncut.measured <= REQUIRED_FRACTION),
        )

        engagement_value = ZERO_FLOAT if attribution.max_engagement_step is None else attribution.max_engagement_step.value
        if type(max_engagement_step) is not DegreesCriterion:
            raise InvalidHeldPathEvidenceError("max_engagement_step must be one validated degrees criterion record.")
        checked_engagement = DegreesCriterion.build(
            name=max_engagement_step.name,
            measured=max_engagement_step.measured,
            required=max_engagement_step.required,
            evidence=max_engagement_step.evidence,
            outcome=max_engagement_step.outcome,
        )
        engagement_evidence = EVIDENCE_BY_CRITERION["max engagement step (deg)"]
        canonical_engagement = DegreesCriterion.build(
            name="max engagement step (deg)",
            measured=Degrees(engagement_value),
            required=Degrees(spec.tea_cap_deg),
            evidence=engagement_evidence,
            outcome=_outcome(engagement_evidence, engagement_value <= spec.tea_cap_deg),
        )

        loop_value = ZERO_FLOAT if attribution.max_loop_radius_step is None else attribution.max_loop_radius_step.value
        if type(max_loop_radius_step) is not ToolRadiusMultipleCriterion:
            raise InvalidHeldPathEvidenceError("max_loop_radius_step must be one validated tool-radius-multiple criterion record.")
        checked_loop = ToolRadiusMultipleCriterion.build(
            name=max_loop_radius_step.name,
            measured=max_loop_radius_step.measured,
            required=max_loop_radius_step.required,
            evidence=max_loop_radius_step.evidence,
            outcome=max_loop_radius_step.outcome,
        )
        loop_evidence = EVIDENCE_BY_CRITERION["max loop radius step (tool radii)"]
        canonical_loop = ToolRadiusMultipleCriterion.build(
            name="max loop radius step (tool radii)",
            measured=ToolRadiusMultiple(loop_value),
            required=REQUIRED_LOOP_STEP,
            evidence=loop_evidence,
            outcome=_outcome(loop_evidence, loop_value <= REQUIRED_LOOP_STEP),
        )
        if (checked_uncut, checked_engagement, checked_loop) != (canonical_uncut, canonical_engagement, canonical_loop):
            raise ContradictoryPathQualityEvidenceError("fraction or maximum criterion disagrees with its canonical threshold, outcome, or source attribution.")

        expected_attribution = _attribution_from_survey(snapshot, survey)
        if attribution != expected_attribution:
            raise ContradictoryPathQualityEvidenceError("path-quality attribution disagrees with its exact source observations.")
        values: dict[str, object] = {
            "uncut_fraction": canonical_uncut,
            **canonical_counts,
            "max_engagement_step": canonical_engagement,
            "max_loop_radius_step": canonical_loop,
            "attribution": attribution,
        }
        return _build_record(cls, values)


def _indices(values: Sequence[int], operation_count: int) -> tuple[OperationIndex, ...]:
    typed = tuple(operation_index(value, operation_count=operation_count) for value in values)
    if len(set(typed)) != len(typed):
        raise InvalidHeldPathEvidenceError("operation attribution must contain unique indices.")
    return typed


def _validated_pairs(values: Sequence[OperationPair], operation_count: int) -> tuple[OperationPair, ...]:
    checked: list[OperationPair] = []
    for value in values:
        if type(value) is not OperationPair:
            raise InvalidHeldPathEvidenceError("pair attribution requires validated operation pairs.")
        checked.append(OperationPair.build(previous=value.previous, current=value.current, operation_count=operation_count))
    if len(set(checked)) != len(checked):
        raise InvalidHeldPathEvidenceError("pair attribution must contain unique operation pairs.")
    return tuple(checked)


def _pair(previous: MotionQuality, current: MotionQuality, operation_count: int) -> OperationPair:
    return OperationPair.build(previous=OperationIndex(previous.index), current=OperationIndex(current.index), operation_count=operation_count)


def _adjacent_motion_pairs(motions: Sequence[MotionQuality]) -> Sequence[tuple[MotionQuality, MotionQuality]]:
    return [(previous, current) for previous, current in zip(motions, motions[1:]) if current.index == previous.index + NEXT_OPERATION_OFFSET]


def _snapshot_kind(operation: HeldOperationSnapshot) -> MotionKind:
    if type(operation) is HeldLineSnapshot:
        return MotionKind.LINE
    if type(operation) is HeldArcSnapshot:
        return MotionKind.ARC
    if type(operation) is HeldCircleSnapshot:
        return MotionKind.LOOP
    raise InvalidHeldPathEvidenceError("source snapshot contains an unsupported operation record.")


def _snapshot_cut_height(snapshot: tuple[HeldOperationSnapshot, ...]) -> float:
    line_anchors: list[Millimetre] = []
    curve_heights: list[Millimetre] = []
    for operation in snapshot:
        if isinstance(operation, HeldLineSnapshot):
            anchor = line_cut_height_anchor(
                operation.operation,
                start_z=Millimetre(float(operation.start.z)),
                end_z=Millimetre(float(operation.end.z)),
                xy_travel=Millimetre(
                    math.hypot(
                        float(operation.end.x) - float(operation.start.x),
                        float(operation.end.y) - float(operation.start.y),
                    )
                ),
            )
            if anchor is not None:
                line_anchors.append(anchor)
        elif operation.operation in {OperationType.CUT, OperationType.LEAD_IN, OperationType.LEAD_OUT}:
            curve_heights.append(Millimetre(float(operation.centre.z)))
    return float(minimum_cut_height(line_anchors if line_anchors else curve_heights))


def _snapshot_replay_category(
    operation: HeldOperationSnapshot,
    cut_height: float,
) -> ReplayCategory:
    try:
        if isinstance(operation, HeldLineSnapshot):
            z_start = float(operation.start.z)
            z_end = float(operation.end.z)
            return classify_line_replay(
                operation.operation,
                start_z=Millimetre(z_start),
                end_z=Millimetre(z_end),
                xy_travel=Millimetre(
                    math.hypot(
                        float(operation.end.x) - float(operation.start.x),
                        float(operation.end.y) - float(operation.start.y),
                    )
                ),
                cut_z=Millimetre(cut_height),
            )
        return classify_planar_replay(
            operation.operation,
            motion_z=Millimetre(float(operation.centre.z)),
            cut_z=Millimetre(cut_height),
        )
    except (CutPlaneRampError, OffPlaneReplayCurveError) as error:
        replay_error = UnreplayableOperationError(f"Operation {operation.ordinal} ({operation.operation.value}) {error}.")
        replay_error.__cause__ = error
        raise InvalidHeldPathEvidenceError(str(replay_error)) from replay_error


def _validate_ordered_indices(values: Sequence[int], operation_count: int, *, name: str) -> tuple[OperationIndex, ...]:
    checked = _indices(values, operation_count)
    if tuple(checked) != tuple(sorted(checked)):
        raise InvalidHeldPathEvidenceError(f"{name} must retain source operation order.")
    return checked


def _validate_survey_binding(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> None:
    if survey.spec is not spec:
        raise InvalidHeldPathEvidenceError("path survey must retain the exact PocketSpec supplied for assessment.")
    if type(snapshot) is not tuple or type(survey.source_snapshot) is not tuple or survey.source_snapshot != snapshot:
        raise InvalidHeldPathEvidenceError("path survey source snapshot disagrees with the supplied operation snapshot.")
    operation_count = len(snapshot)
    if operation_count == ZERO_COUNT:
        raise InvalidHeldPathEvidenceError("a Figure 5 quality assessment requires at least one operation.")
    for expected, operation in enumerate(snapshot):
        if operation.ordinal != expected:
            raise InvalidHeldPathEvidenceError("snapshot ordinals must cover source operations in order.")

    motion_indices = _validate_ordered_indices([motion.index for motion in survey.motions], operation_count, name="survey motions")
    rapid_indices = _validate_ordered_indices([rapid.index for rapid in survey.rapids], operation_count, name="survey rapids")
    plunge_indices = _validate_ordered_indices(survey.plunge_indices, operation_count, name="survey plunges")
    retract_indices = _validate_ordered_indices(survey.retract_indices, operation_count, name="survey retracts")
    if type(survey.plunges) is not int or survey.plunges != len(plunge_indices):
        raise InvalidHeldPathEvidenceError("plunge count disagrees with retained plunge operation indices.")
    if type(survey.retracts) is not int or survey.retracts != len(retract_indices):
        raise InvalidHeldPathEvidenceError("retract count disagrees with retained retract operation indices.")

    motion_set = set(motion_indices)
    rapid_set = set(rapid_indices)
    plunge_set = set(plunge_indices)
    retract_set = set(retract_indices)
    if not retract_set <= rapid_set:
        raise InvalidHeldPathEvidenceError("every retained retract must have one non-cutting motion observation.")
    pure_rapid_set = rapid_set - retract_set
    categories = (motion_set, pure_rapid_set, plunge_set, retract_set)
    for position, observed_category in enumerate(categories):
        if any(observed_category & later for later in categories[position + NEXT_OPERATION_OFFSET :]):
            raise InvalidHeldPathEvidenceError("motion, rapid, plunge, and retract observations must be disjoint.")
    if set().union(*categories) != set(range(operation_count)):
        raise InvalidHeldPathEvidenceError("motion, rapid, plunge, and retract observations must exactly partition the operation snapshot.")

    cut_height = _snapshot_cut_height(snapshot)
    expected_categories: dict[str, set[OperationIndex]] = {"motion": set(), "rapid": set(), "plunge": set(), "retract": set()}
    for index, operation in enumerate(snapshot):
        expected_category = _snapshot_replay_category(operation, cut_height)
        expected_categories[expected_category].add(OperationIndex(index))
    observed_categories: dict[str, set[OperationIndex]] = {
        "motion": motion_set,
        "rapid": pure_rapid_set,
        "plunge": plunge_set,
        "retract": retract_set,
    }
    if observed_categories != expected_categories:
        raise InvalidHeldPathEvidenceError("survey categories disagree with the source snapshot's cut-plane replay classification.")

    for motion in survey.motions:
        source = snapshot[motion.index]
        if motion.operation is not source.operation or motion.kind is not _snapshot_kind(source):
            raise InvalidHeldPathEvidenceError("cut-motion role or primitive kind disagrees with the source snapshot.")
    for rapid in survey.rapids:
        source = snapshot[rapid.index]
        if rapid.operation is not source.operation or rapid.kind is not _snapshot_kind(source):
            raise InvalidHeldPathEvidenceError("rapid-motion role or primitive kind disagrees with the source snapshot.")


def _maximum_step(steps: Sequence[tuple[float, OperationPair]], unit: Literal["degrees", "tool_radius_multiple"]) -> Optional[StepFinding]:
    winner: Optional[tuple[float, OperationPair]] = None
    for step in steps:
        if winner is None or step[0] > winner[0]:
            winner = step
    if winner is None or winner[0] <= ZERO_FLOAT:
        return None
    if unit == "degrees":
        return MeasuredStep[Degrees].build(value=Degrees(winner[0]), pair=winner[1], unit="degrees")
    return MeasuredStep[ToolRadiusMultiple].build(value=ToolRadiusMultiple(winner[0]), pair=winner[1], unit="tool_radius_multiple")


def _loop_steps(
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    tool_radius: float,
    operation_count: int,
) -> list[tuple[float, OperationPair]]:
    rapid_indices = sorted(rapid.index for rapid in survey.rapids)
    chain_of = {int(operation.ordinal): operation.path_index for operation in snapshot}
    runs: list[list[tuple[int, float]]] = [[]]
    position = ZERO_COUNT
    chain: Optional[int] = None
    for motion in survey.motions:
        boundary = False
        while position < len(rapid_indices) and rapid_indices[position] < motion.index:
            position += NEXT_OPERATION_OFFSET
            boundary = True
        current_chain = chain_of.get(motion.index, chain)
        if chain is not None and current_chain != chain:
            boundary = True
        chain = current_chain
        if boundary and runs[-1]:
            runs.append([])
        if motion.kind is MotionKind.LOOP and motion.loop_radius is not None:
            runs[-1].append((motion.index, motion.loop_radius))
    steps: list[tuple[float, OperationPair]] = []
    for run in runs:
        for earlier, later in zip(run, run[1:]):
            value = abs(later[1] - earlier[1]) / tool_radius
            pair = OperationPair.build(previous=OperationIndex(earlier[0]), current=OperationIndex(later[0]), operation_count=operation_count)
            steps.append((value, pair))
    return steps


def _attribution_from_survey(
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> PathQualityAttribution:
    operation_count = len(snapshot)
    spec = survey.spec
    gouging = _indices(
        [motion.index for motion in survey.motions if any(not sample.inside_centre_domain for sample in motion.samples)],
        operation_count,
    )
    unsafe = _indices([rapid.index for rapid in survey.rapids if rapid.horizontal_at_cut_plane], operation_count)
    zero_length = _indices(
        sorted([motion.index for motion in survey.motions if motion.length == ZERO_FLOAT] + [rapid.index for rapid in survey.rapids if rapid.length == ZERO_FLOAT]),
        operation_count,
    )
    degenerate = _indices(
        [
            motion.index
            for motion in survey.motions
            if motion.kind is MotionKind.LOOP and motion.loop_radius is not None and motion.loop_radius <= DEGENERATE_LOOP_RATIO * spec.tool_radius
        ],
        operation_count,
    )
    redundant = _indices([motion.index for motion in survey.motions if not motion.removes_material], operation_count)
    cap_exceeded = _indices([motion.index for motion in survey.motions if any(sample.cap_exceeded for sample in motion.samples)], operation_count)
    slotting = _indices([motion.index for motion in survey.motions if motion.slot_exceeded], operation_count)

    continuity: list[OperationPair] = []
    tangent: list[OperationPair] = []
    curvature: list[OperationPair] = []
    reversals: list[OperationPair] = []
    engagement_steps: list[tuple[float, OperationPair]] = []
    continuity_tolerance = CONTINUITY_TOOL_RADIUS_FRACTION * spec.tool_radius
    for previous, current in _adjacent_motion_pairs(survey.motions):
        pair = _pair(previous, current, operation_count)
        if math.hypot(current.start[0] - previous.end[0], current.start[1] - previous.end[1]) > continuity_tolerance:
            continuity.append(pair)
        dot = previous.end_tangent[0] * current.start_tangent[0] + previous.end_tangent[1] * current.start_tangent[1]
        if dot < ZERO_FLOAT:
            reversals.append(pair)
            tangent.append(pair)
        elif dot < UNIT_DOT - TANGENT_CONTINUITY_SLACK:
            tangent.append(pair)
        elif previous.curvature != current.curvature:
            curvature.append(pair)
        previous_peak = max((sample.engagement_deg for sample in previous.samples), default=ZERO_FLOAT)
        current_peak = max((sample.engagement_deg for sample in current.samples), default=ZERO_FLOAT)
        engagement_steps.append((abs(current_peak - previous_peak), pair))

    loop_steps = _loop_steps(snapshot, survey, spec.tool_radius, operation_count)
    maximum_engagement = cast(Optional[MeasuredStep[Degrees]], _maximum_step(engagement_steps, "degrees"))
    maximum_loop = cast(Optional[MeasuredStep[ToolRadiusMultiple]], _maximum_step(loop_steps, "tool_radius_multiple"))
    engagement_failures = tuple(pair for value, pair in engagement_steps if value > spec.tea_cap_deg)
    loop_failures = tuple(pair for value, pair in loop_steps if value > REQUIRED_LOOP_STEP)
    return PathQualityAttribution.build(
        operation_count=operation_count,
        uncut_operations=(),
        gouging_operations=gouging,
        unsafe_rapid_operations=unsafe,
        continuity_break_pairs=tuple(continuity),
        zero_length_operations=zero_length,
        degenerate_loop_operations=degenerate,
        redundant_operations=redundant,
        cap_exceeded_operations=cap_exceeded,
        slotting_operations=slotting,
        max_engagement_step=maximum_engagement,
        engagement_step_failure_pairs=engagement_failures,
        max_loop_radius_step=maximum_loop,
        loop_radius_step_failure_pairs=loop_failures,
        tangent_break_pairs=tuple(tangent),
        curvature_break_pairs=tuple(curvature),
        reversal_pairs=tuple(reversals),
    )


def _count(
    name: CountCriterionName,
    cardinality: int,
) -> CountCriterion:
    evidence = EVIDENCE_BY_CRITERION[name]
    return CountCriterion.build(
        name=name,
        measured=MotionCount(cardinality),
        required=REQUIRED_COUNT,
        evidence=evidence,
        outcome=_outcome(evidence, cardinality <= REQUIRED_COUNT),
        attribution_count=cardinality,
    )


def assess_path_quality(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    coverage: CoverageEstimate,
) -> PathQualityAssessment:
    """Reproduce the current twelve gate decisions with source attribution."""
    _validate_survey_binding(spec, snapshot, survey)
    attribution = _attribution_from_survey(snapshot, survey)
    max_engagement_value = Degrees(ZERO_FLOAT if attribution.max_engagement_step is None else attribution.max_engagement_step.value)
    max_loop_value = ToolRadiusMultiple(ZERO_FLOAT if attribution.max_loop_radius_step is None else attribution.max_loop_radius_step.value)
    uncut = closed_unit_fraction(coverage.uncut_fraction, name="uncut fraction measured")
    uncut_evidence = EVIDENCE_BY_CRITERION["uncut fraction"]
    engagement_evidence = EVIDENCE_BY_CRITERION["max engagement step (deg)"]
    loop_evidence = EVIDENCE_BY_CRITERION["max loop radius step (tool radii)"]
    return PathQualityAssessment.build(
        spec=spec,
        uncut_fraction=FractionCriterion.build(
            name="uncut fraction",
            measured=uncut,
            required=REQUIRED_FRACTION,
            evidence=uncut_evidence,
            outcome=_outcome(uncut_evidence, uncut <= REQUIRED_FRACTION),
        ),
        gouging_motions=_count("gouging motions", len(attribution.gouging_operations)),
        unsafe_rapids=_count("unsafe rapids", len(attribution.unsafe_rapid_operations)),
        continuity_breaks=_count("continuity breaks", len(attribution.continuity_break_pairs)),
        zero_length_motions=_count("zero-length motions", len(attribution.zero_length_operations)),
        degenerate_loops=_count("degenerate loops", len(attribution.degenerate_loop_operations)),
        redundant_operations=_count("redundant operations", len(attribution.redundant_operations)),
        cap_exceedances=_count("cap exceedances", len(attribution.cap_exceeded_operations)),
        slotting_motions=_count("slotting motions", len(attribution.slotting_operations)),
        max_engagement_step=DegreesCriterion.build(
            name="max engagement step (deg)",
            measured=max_engagement_value,
            required=Degrees(spec.tea_cap_deg),
            evidence=engagement_evidence,
            outcome=_outcome(engagement_evidence, max_engagement_value <= spec.tea_cap_deg),
        ),
        max_loop_radius_step=ToolRadiusMultipleCriterion.build(
            name="max loop radius step (tool radii)",
            measured=max_loop_value,
            required=REQUIRED_LOOP_STEP,
            evidence=loop_evidence,
            outcome=_outcome(loop_evidence, max_loop_value <= REQUIRED_LOOP_STEP),
        ),
        tangent_breaks=_count("tangent breaks", len(attribution.tangent_break_pairs)),
        attribution=attribution,
        snapshot=snapshot,
        survey=survey,
    )
