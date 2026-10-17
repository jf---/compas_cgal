"""Typed, operation-attributed observations for the existing path-quality gate.

It makes the thirteen gated decisions of the machining-quality gate (rules
a-f of Task 0.3: redundant counts only CUT motions, plunges into material are
bounded by the centre-domain components, a plunge must land on the next
motion's start, slotting is the straight subset of cap exceedances, degenerate
loops are reported and not gated, and the engagement step is read at the
junction) and attributes each to source operations.

A helical entry (Task 1.2, rules "1.2") is an entry event, not a motion: it is
never cap-sampled (E8), it counts as a plunge into material when its exact
bore met the stock (Q3(a)), rule (e) covers its exit (N7), and its exit pairs
with the next motion for continuity, tangency and the junction step. It is no
loop-radius-step operand and it splits a loop run. Its ramp verdict needs a
machine and is `benchmarks.quality.ramp_angle_violations`, outside the
thirteen.
The criterion and record types live in `benchmarks.quality_criteria`; the
snapshot dispatch and survey binding in `benchmarks.quality_binding`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal
from typing import Optional
from typing import Sequence
from typing import cast

from typing_extensions import Self

from benchmarks.coverage import COVERAGE_GRID_SAMPLES
from benchmarks.coverage import CoverageEstimate
from benchmarks.coverage import measure_coverage
from benchmarks.errors import ContradictoryPathQualityEvidenceError
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import MissingLoopSeamStationError
from benchmarks.errors import ZeroLengthToolpathError
from benchmarks.held_path_snapshot import HeldHelixSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality import CONTINUITY_TOOL_RADIUS_FRACTION
from benchmarks.quality import DEGENERATE_LOOP_RATIO
from benchmarks.quality import TANGENT_CONTINUITY_SLACK
from benchmarks.quality import PathQuality
from benchmarks.quality import _cut
from benchmarks.quality import _elementary
from benchmarks.quality import _longevity
from benchmarks.quality import _program
from benchmarks.quality import _speed
from benchmarks.quality_binding import NEXT_OPERATION_OFFSET
from benchmarks.quality_binding import ZERO_COUNT
from benchmarks.quality_binding import HelixHasNoMotionKindError as HelixHasNoMotionKindError
from benchmarks.quality_binding import _indices
from benchmarks.quality_binding import _snapshot_cut_height as _snapshot_cut_height
from benchmarks.quality_binding import _snapshot_replay_category as _snapshot_replay_category
from benchmarks.quality_binding import _validate_survey_binding
from benchmarks.quality_binding import snapshot_motion_kind as snapshot_motion_kind
from benchmarks.quality_criteria import COMPONENT_BOUNDED_COUNT_NAME
from benchmarks.quality_criteria import CRITERION_NAMES as CRITERION_NAMES
from benchmarks.quality_criteria import EVIDENCE_BY_CRITERION as EVIDENCE_BY_CRITERION
from benchmarks.quality_criteria import REQUIRED_COUNT
from benchmarks.quality_criteria import REQUIRED_FRACTION
from benchmarks.quality_criteria import REQUIRED_LOOP_STEP
from benchmarks.quality_criteria import ZERO_FLOAT
from benchmarks.quality_criteria import CountCriterion as CountCriterion
from benchmarks.quality_criteria import CountCriterionName
from benchmarks.quality_criteria import DegreesCriterion as DegreesCriterion
from benchmarks.quality_criteria import FractionCriterion as FractionCriterion
from benchmarks.quality_criteria import MeasuredStep as MeasuredStep
from benchmarks.quality_criteria import OperationPair as OperationPair
from benchmarks.quality_criteria import PathQualityAttribution as PathQualityAttribution
from benchmarks.quality_criteria import StepFinding
from benchmarks.quality_criteria import ToolRadiusMultipleCriterion as ToolRadiusMultipleCriterion
from benchmarks.quality_criteria import _build_record as _build_record
from benchmarks.quality_criteria import _outcome
from benchmarks.spec import PocketSpec
from benchmarks.survey import MotionKind
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.units import Degrees
from benchmarks.units import MotionCount
from benchmarks.units import OperationIndex
from benchmarks.units import ToolRadiusMultiple
from benchmarks.units import closed_unit_fraction
from compas_cgal.toolpath import OperationType

UNIT_DOT = 1.0


# Index offset from a motion to the station that holds its junction value: the
# first station with a non-empty swept prefix. Station 0 reads the stock minus
# the tool's own start disk, which is 0 deg by construction (Ruling
# 0.3-f-station), so it says nothing about the junction.
JUNCTION_STATION = 1


@dataclass(frozen=True, init=False)
class PathQualityAssessment:
    uncut_fraction: FractionCriterion
    gouging_motions: CountCriterion
    unsafe_rapids: CountCriterion
    continuity_breaks: CountCriterion
    zero_length_motions: CountCriterion
    redundant_operations: CountCriterion
    material_plunges: CountCriterion
    misplaced_plunges: CountCriterion
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
        redundant_operations: CountCriterion,
        material_plunges: CountCriterion,
        misplaced_plunges: CountCriterion,
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
            "redundant_operations": len(attribution.redundant_operations),
            "material_plunges": len(attribution.material_plunge_operations),
            "misplaced_plunges": len(attribution.misplaced_plunge_operations),
            "cap_exceedances": len(attribution.cap_exceeded_operations),
            "slotting_motions": len(attribution.slotting_operations),
            "tangent_breaks": len(attribution.tangent_break_pairs),
        }
        supplied_counts = {
            "gouging_motions": gouging_motions,
            "unsafe_rapids": unsafe_rapids,
            "continuity_breaks": continuity_breaks,
            "zero_length_motions": zero_length_motions,
            "redundant_operations": redundant_operations,
            "material_plunges": material_plunges,
            "misplaced_plunges": misplaced_plunges,
            "cap_exceedances": cap_exceedances,
            "slotting_motions": slotting_motions,
            "tangent_breaks": tangent_breaks,
        }
        count_names: dict[str, CountCriterionName] = {
            "gouging_motions": "gouging motions",
            "unsafe_rapids": "unsafe rapids",
            "continuity_breaks": "continuity breaks",
            "zero_length_motions": "zero-length motions",
            "redundant_operations": "redundant operations",
            "material_plunges": "material plunges",
            "misplaced_plunges": "misplaced plunges",
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
        # The material-plunge threshold is the centre domain's component count,
        # which the coverage pass measures (one predicate per reduction); it is
        # carried by the supplied criterion, validated there as a count of at
        # least one, and bound to the coverage record by the evidence boundary
        # (`benchmarks.held_path_evidence`). Every other count threshold is zero.
        required_counts = {
            field: checked_counts["material_plunges"].required if count_names[field] == COMPONENT_BOUNDED_COUNT_NAME else REQUIRED_COUNT for field in supplied_counts
        }
        canonical_counts = {field: _count(count_names[field], cardinalities[field], required_counts[field]) for field in supplied_counts}
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
    # A rapid splits a run, and so does a helix: it is an entry, never a loop
    # step operand (it is not among the motions at all), and the loops either
    # side of it are not run back to back.
    rapid_indices = sorted([rapid.index for rapid in survey.rapids] + [helix.index for helix in survey.helices])
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
    # Rule (a): the label decides (Ruling 0.3-C11). A CUT owns removal, so a CUT
    # that removes nothing is waste; a LINK or lead at the cut plane that removes
    # nothing is a stay-down transit, reported as air at depth and still judged
    # by the cap like every motion.
    redundant = _indices([motion.index for motion in survey.motions if not motion.removes_material and motion.operation is OperationType.CUT], operation_count)
    air_at_depth = _indices([motion.index for motion in survey.motions if not motion.removes_material and motion.operation is not OperationType.CUT], operation_count)
    # Q3(a): a plunge into material is a PLUNGE or a HELIX whose exact bore met
    # the stock's interior; the key "material plunges" is kept.
    material_plunges = _indices(sorted([*survey.material_plunge_indices, *(helix.index for helix in survey.helices if helix.into_material)]), operation_count)
    misplaced_plunges = _misplaced_plunges(snapshot, survey, operation_count)
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
        junction = _junction_step(previous, current)
        if junction is not None:
            engagement_steps.append((junction, pair))
    motion_at = {motion.index: motion for motion in survey.motions}
    for helix in survey.helices:
        following = motion_at.get(helix.index + NEXT_OPERATION_OFFSET)
        if following is None:
            continue
        pair = OperationPair.build(previous=OperationIndex(helix.index), current=OperationIndex(following.index), operation_count=operation_count)
        if math.hypot(following.start[0] - helix.end[0], following.start[1] - helix.end[1]) > continuity_tolerance:
            continuity.append(pair)
        dot = helix.end_tangent[0] * following.start_tangent[0] + helix.end_tangent[1] * following.start_tangent[1]
        if dot < ZERO_FLOAT:
            reversals.append(pair)
            tangent.append(pair)
        elif dot < UNIT_DOT - TANGENT_CONTINUITY_SLACK:
            tangent.append(pair)
        elif helix.curvature != following.curvature:
            curvature.append(pair)
        if len(following.samples) > JUNCTION_STATION:
            engagement_steps.append((abs(following.samples[JUNCTION_STATION].engagement_deg - helix.exit_sample.engagement_deg), pair))

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
        air_at_depth_operations=air_at_depth,
        material_plunge_operations=material_plunges,
        misplaced_plunge_operations=misplaced_plunges,
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


def _junction_step(previous: MotionQuality, current: MotionQuality) -> Optional[float]:
    """The engagement change across one junction, or None when it is not read.

    Gate rule (f): the value is the next motion's first station with a
    non-empty prefix (`JUNCTION_STATION`) against the station where the
    previous motion ends (`MotionQuality.exit_sample`: a loop's seam after the
    full turn, Ruling 0.3-f-loop; a line's or arc's last station). A motion
    with no such station (a zero-length move, which is one station with an
    empty prefix) has no junction reading, and neither does a previous motion
    with no stations. REPORTING degrees (Ruling 0.3-C6): the step is a sampled
    diagnostic; the cap verdicts at those stations stay exact.

    The leaving station is chosen by the previous motion's kind, never by
    whether a field happens to be set: a line or an arc hands over at its last
    station by definition; a loop hands over at its seam, which only
    `MotionQuality.exit_sample` carries. A loop record without it is refused
    rather than read at its last per-motion station (that reading is the I1
    defect Ruling 0.3-f-loop fixed).

    Args:
        previous: The motion the tool leaves.
        current: The motion it runs next without leaving the cut.

    Returns:
        ``|current[1] - exit(previous)|`` in degrees, or None.

    Raises:
        MissingLoopSeamStationError: *previous* is a loop with no seam station.
    """
    if not previous.samples or len(current.samples) <= JUNCTION_STATION:
        return None
    if previous.kind is MotionKind.LOOP:
        if previous.exit_sample is None:
            raise MissingLoopSeamStationError(f"motion {previous.index} is a loop without its seam station (exit_sample); its junction is read at the seam (Ruling 0.3-f-loop).")
        leaving = previous.exit_sample
    else:
        leaving = previous.samples[-1]
    return abs(current.samples[JUNCTION_STATION].engagement_deg - leaving.engagement_deg)


def _misplaced_plunges(
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    operation_count: int,
) -> tuple[OperationIndex, ...]:
    """Plunges whose foot, and helices whose exit, is not the next motion's start (gate rule e, N7).

    Exact: the replay injects both binary64 pairs as themselves, so equal
    doubles are the same exact point and unequal ones are not; there is no
    tolerance. A plunge or helix whose next operation is not a cut-plane motion
    has no start to land on and is misplaced. A helix exits at its footprint's
    exact entry ``P`` (the snapshot record's ``entry``); the next motion's start
    is read from its ``entry_xy`` (Ruling 1.1-C14-A).

    Args:
        snapshot: The source operations.
        survey: The validated survey.
        operation_count: Operations in the path.

    Returns:
        The misplaced plunges, in source order.

    Raises:
        InvalidHeldPathEvidenceError: A surveyed plunge is not a line operation.
    """
    starts = {motion.index: motion.start for motion in survey.motions}
    misplaced: list[int] = []
    for index in survey.plunge_indices:
        plunge = snapshot[index]
        if type(plunge) is not HeldLineSnapshot:
            raise InvalidHeldPathEvidenceError(f"plunge {index} is not a line operation.")
        foot = (float(plunge.end.x), float(plunge.end.y))
        if starts.get(int(index) + NEXT_OPERATION_OFFSET) != foot:
            misplaced.append(int(index))
    for helix in survey.helices:
        record = snapshot[helix.index]
        if not isinstance(record, HeldHelixSnapshot):
            raise InvalidHeldPathEvidenceError(f"helix {helix.index} is not a helix record.")
        if starts.get(helix.index + NEXT_OPERATION_OFFSET) != (float(record.entry.x), float(record.entry.y)):
            misplaced.append(helix.index)
    return _indices(sorted(misplaced), operation_count)


def _count(
    name: CountCriterionName,
    cardinality: int,
    required: MotionCount = REQUIRED_COUNT,
) -> CountCriterion:
    evidence = EVIDENCE_BY_CRITERION[name]
    return CountCriterion.build(
        name=name,
        measured=MotionCount(cardinality),
        required=required,
        evidence=evidence,
        outcome=_outcome(evidence, cardinality <= required),
        attribution_count=cardinality,
    )


def assess_path_quality(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    coverage: CoverageEstimate,
) -> PathQualityAssessment:
    """Make the thirteen gate decisions with source attribution."""
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
        redundant_operations=_count("redundant operations", len(attribution.redundant_operations)),
        material_plunges=_count("material plunges", len(attribution.material_plunge_operations), MotionCount(coverage.centre_domain_components)),
        misplaced_plunges=_count("misplaced plunges", len(attribution.misplaced_plunge_operations)),
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


@dataclass(frozen=True, init=False)
class QualityEvidence:
    """One canonical quality reduction and the evidence it was built from."""

    path_quality: PathQuality
    assessment: PathQualityAssessment
    coverage: CoverageEstimate

    def __init__(self) -> None:
        raise TypeError("QualityEvidence is returned only by reduce_quality_evidence().")


def _material_entry_count(survey: PathSurvey) -> int:
    """Count first cut contacts after plunges and helical entries from the validated survey."""
    plunge_indices = set(survey.plunge_indices) | {helix.index for helix in survey.helices}
    motion_indices = {motion.index for motion in survey.motions}
    armed = False
    entries = 0
    for index in range(len(survey.source_snapshot)):
        if index in plunge_indices:
            armed = True
        elif index in motion_indices and armed:
            entries += NEXT_OPERATION_OFFSET
            armed = False
    return entries


def reduce_quality_evidence(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    *,
    grid: int = COVERAGE_GRID_SAMPLES,
) -> QualityEvidence:
    """Reduce one validated survey through the canonical quality assessment."""
    if survey.total_length <= ZERO_FLOAT:
        raise ZeroLengthToolpathError(f"{spec.name}: the toolpath sums to zero length, so no length fraction is defined.")
    _validate_survey_binding(spec, snapshot, survey)
    coverage = measure_coverage(spec, survey.final_stock, grid=grid)
    assessment = assess_path_quality(spec, snapshot, survey, coverage)
    quality = PathQuality(
        elementary=_elementary(spec, survey, coverage.remaining_area, assessment),
        cut=_cut(spec, survey, coverage.wall_scallop_height, assessment),
        speed=_speed(survey, assessment),
        longevity=_longevity(survey, _material_entry_count(survey)),
        program=_program(survey),
        cut_operations=len(survey.motions),
        path_length=survey.total_length,
    )
    return _build_record(
        QualityEvidence,
        {
            "path_quality": quality,
            "assessment": assessment,
            "coverage": coverage,
        },
    )
