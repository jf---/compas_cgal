"""Criterion and record types of the path-quality gate.

The thirteen criterion names, their evidence kinds, and the validated records
(criteria, operation pairs, measured steps, attribution) that the assessment
and reductions in `benchmarks.quality_observations` build. Split out of that
module.
"""

from __future__ import annotations

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

from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.quality_binding import ZERO_COUNT
from benchmarks.quality_binding import _indices
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
    "redundant operations",
    "material plunges",
    "misplaced plunges",
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
    "redundant operations",
    "material plunges",
    "misplaced plunges",
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
    "redundant operations": "exact_depletion_replay",
    "material plunges": "exact_depletion_replay",
    "misplaced plunges": "derived_geometry",
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
    "redundant operations",
    "material plunges",
    "misplaced plunges",
    "cap exceedances",
    "slotting motions",
    "tangent breaks",
]


ZERO_FLOAT = 0.0

REQUIRED_FRACTION = UnitFraction(ZERO_FLOAT)
REQUIRED_COUNT = MotionCount(ZERO_COUNT)
REQUIRED_LOOP_STEP = ToolRadiusMultiple(2.0)

# The one count criterion whose threshold is not zero: a connected centre
# domain admits exactly one plunge into material (gate rule b), so the
# threshold is the number of components, which is at least one for any
# machinable pocket.
COMPONENT_BOUNDED_COUNT_NAME = "material plunges"
MINIMUM_CENTRE_DOMAIN_COMPONENTS = 1

FRACTION_NAMES = frozenset({"uncut fraction"})
COUNT_NAMES = frozenset(
    {
        "gouging motions",
        "unsafe rapids",
        "continuity breaks",
        "zero-length motions",
        "redundant operations",
        "material plunges",
        "misplaced plunges",
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
        if checked_name == COMPONENT_BOUNDED_COUNT_NAME:
            if checked_required < MINIMUM_CENTRE_DOMAIN_COMPONENTS:
                raise InvalidHeldPathEvidenceError(f"{checked_name} requires a threshold of at least one centre-domain component.")
        elif checked_required != REQUIRED_COUNT:
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
    air_at_depth_operations: tuple[OperationIndex, ...]
    material_plunge_operations: tuple[OperationIndex, ...]
    misplaced_plunge_operations: tuple[OperationIndex, ...]
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
        air_at_depth_operations: tuple[OperationIndex, ...],
        material_plunge_operations: tuple[OperationIndex, ...],
        misplaced_plunge_operations: tuple[OperationIndex, ...],
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
            air_at_depth_operations,
            material_plunge_operations,
            misplaced_plunge_operations,
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
            "air_at_depth_operations": air_at_depth_operations,
            "material_plunge_operations": material_plunge_operations,
            "misplaced_plunge_operations": misplaced_plunge_operations,
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
        if not set(redundant_operations).isdisjoint(air_at_depth_operations):
            raise InvalidHeldPathEvidenceError("a motion that removes nothing is either redundant (a CUT) or air at depth, never both.")
        return _build_record(cls, checked)


def _validated_pairs(values: Sequence[OperationPair], operation_count: int) -> tuple[OperationPair, ...]:
    checked: list[OperationPair] = []
    for value in values:
        if type(value) is not OperationPair:
            raise InvalidHeldPathEvidenceError("pair attribution requires validated operation pairs.")
        checked.append(OperationPair.build(previous=value.previous, current=value.current, operation_count=operation_count))
    if len(set(checked)) != len(checked):
        raise InvalidHeldPathEvidenceError("pair attribution must contain unique operation pairs.")
    return tuple(checked)
