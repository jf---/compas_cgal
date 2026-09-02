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

from typing_extensions import Self
from typing_extensions import TypeAlias

from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import InvalidHeldPathEvidenceError
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

    def __init__(self) -> None:
        raise TypeError("MeasuredStep must be created with MeasuredStep.build().")

    @classmethod
    def build(cls, *, value: StepUnitT, pair: OperationPair, unit: Literal["degrees", "tool_radius_multiple"]) -> Self:
        if type(pair) is not OperationPair:
            raise InvalidHeldPathEvidenceError("a measured step requires one validated operation pair.")
        if unit == "degrees":
            checked_degrees = degrees_value(value, name="measured engagement step")
            if checked_degrees < ZERO_FLOAT:
                raise InvalidHeldPathEvidenceError("measured engagement step must be non-negative.")
            return _build_record(cls, {"value": checked_degrees, "pair": pair})
        checked_multiple = tool_radius_multiple(value, name="measured loop-radius step")
        return _build_record(cls, {"value": checked_multiple, "pair": pair})


StepFinding: TypeAlias = Union[MeasuredStep[Degrees], MeasuredStep[ToolRadiusMultiple]]


@dataclass(frozen=True, init=False)
class PathQualityAttribution:
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
        checked: dict[str, object] = {}
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
            checked[name] = step
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
        expected_names = {
            "uncut_fraction": "uncut fraction",
            "gouging_motions": "gouging motions",
            "unsafe_rapids": "unsafe rapids",
            "continuity_breaks": "continuity breaks",
            "zero_length_motions": "zero-length motions",
            "degenerate_loops": "degenerate loops",
            "redundant_operations": "redundant operations",
            "cap_exceedances": "cap exceedances",
            "slotting_motions": "slotting motions",
            "max_engagement_step": "max engagement step (deg)",
            "max_loop_radius_step": "max loop radius step (tool radii)",
            "tangent_breaks": "tangent breaks",
        }
        values: dict[str, object] = {
            "uncut_fraction": uncut_fraction,
            "gouging_motions": gouging_motions,
            "unsafe_rapids": unsafe_rapids,
            "continuity_breaks": continuity_breaks,
            "zero_length_motions": zero_length_motions,
            "degenerate_loops": degenerate_loops,
            "redundant_operations": redundant_operations,
            "cap_exceedances": cap_exceedances,
            "slotting_motions": slotting_motions,
            "max_engagement_step": max_engagement_step,
            "max_loop_radius_step": max_loop_radius_step,
            "tangent_breaks": tangent_breaks,
        }
        for field, expected_name in expected_names.items():
            criterion = values[field]
            if not isinstance(criterion, (FractionCriterion, CountCriterion, DegreesCriterion, ToolRadiusMultipleCriterion)):
                raise InvalidHeldPathEvidenceError(f"{field} must be one validated criterion record.")
            if criterion.name != expected_name:
                raise InvalidHeldPathEvidenceError(f"{field} carries criterion {criterion.name!r}, expected {expected_name!r}.")
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
        for field, cardinality in cardinalities.items():
            criterion = cast(CountCriterion, values[field])
            if criterion.measured != cardinality:
                raise InvalidHeldPathEvidenceError(f"{field} disagrees with the assessment attribution cardinality.")
        engagement_value = ZERO_FLOAT if attribution.max_engagement_step is None else attribution.max_engagement_step.value
        if max_engagement_step.measured != engagement_value:
            raise InvalidHeldPathEvidenceError("maximum engagement step disagrees with its attributed extremum.")
        loop_value = ZERO_FLOAT if attribution.max_loop_radius_step is None else attribution.max_loop_radius_step.value
        if max_loop_radius_step.measured != loop_value:
            raise InvalidHeldPathEvidenceError("maximum loop-radius step disagrees with its attributed extremum.")
        values["attribution"] = attribution
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
    operation_count = len(snapshot)
    for expected, operation in enumerate(snapshot):
        if operation.ordinal != expected:
            raise InvalidHeldPathEvidenceError("snapshot ordinals must cover source operations in order.")
    all_observed = [motion.index for motion in survey.motions] + [rapid.index for rapid in survey.rapids]
    if all_observed and operation_count == ZERO_COUNT:
        raise InvalidHeldPathEvidenceError("survey operations require a non-empty operation snapshot.")
    _indices(all_observed, operation_count) if all_observed else ()

    gouging = _indices([motion.index for motion in survey.motions if motion.gouges], operation_count)
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
        engagement_steps.append((abs(current.peak_engagement_deg - previous.peak_engagement_deg), pair))

    loop_steps = _loop_steps(snapshot, survey, spec.tool_radius, operation_count)
    maximum_engagement = cast(Optional[MeasuredStep[Degrees]], _maximum_step(engagement_steps, "degrees"))
    maximum_loop = cast(Optional[MeasuredStep[ToolRadiusMultiple]], _maximum_step(loop_steps, "tool_radius_multiple"))
    max_engagement_value = Degrees(ZERO_FLOAT if maximum_engagement is None else maximum_engagement.value)
    max_loop_value = ToolRadiusMultiple(ZERO_FLOAT if maximum_loop is None else maximum_loop.value)
    engagement_failures = tuple(pair for value, pair in engagement_steps if value > spec.tea_cap_deg)
    loop_failures = tuple(pair for value, pair in loop_steps if value > REQUIRED_LOOP_STEP)

    attribution = PathQualityAttribution.build(
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
    uncut = closed_unit_fraction(coverage.uncut_fraction, name="uncut fraction measured")
    uncut_evidence = EVIDENCE_BY_CRITERION["uncut fraction"]
    engagement_evidence = EVIDENCE_BY_CRITERION["max engagement step (deg)"]
    loop_evidence = EVIDENCE_BY_CRITERION["max loop radius step (tool radii)"]
    return PathQualityAssessment.build(
        uncut_fraction=FractionCriterion.build(
            name="uncut fraction",
            measured=uncut,
            required=REQUIRED_FRACTION,
            evidence=uncut_evidence,
            outcome=_outcome(uncut_evidence, uncut <= REQUIRED_FRACTION),
        ),
        gouging_motions=_count("gouging motions", len(gouging)),
        unsafe_rapids=_count("unsafe rapids", len(unsafe)),
        continuity_breaks=_count("continuity breaks", len(continuity)),
        zero_length_motions=_count("zero-length motions", len(zero_length)),
        degenerate_loops=_count("degenerate loops", len(degenerate)),
        redundant_operations=_count("redundant operations", len(redundant)),
        cap_exceedances=_count("cap exceedances", len(cap_exceeded)),
        slotting_motions=_count("slotting motions", len(slotting)),
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
        tangent_breaks=_count("tangent breaks", len(tangent)),
        attribution=attribution,
    )
