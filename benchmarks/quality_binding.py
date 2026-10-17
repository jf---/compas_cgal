"""Snapshot dispatch and survey binding for the path-quality gate.

Binds a `PathSurvey` to the exact operation snapshot it was taken from: the
survey kind of each snapshot record, the cut-plane replay category, and the
helix round trip. Split out of `benchmarks.quality_observations`.
"""

from __future__ import annotations

import math
from typing import NoReturn
from typing import Sequence
from typing import get_args

from typing_extensions import Never

from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import UnreplayableOperationError
from benchmarks.held_path_snapshot import HeldArcSnapshot
from benchmarks.held_path_snapshot import HeldCircleSnapshot
from benchmarks.held_path_snapshot import HeldHelixSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.helix_observation import HelixObservation
from benchmarks.spec import PocketSpec
from benchmarks.survey import MotionKind
from benchmarks.survey import PathSurvey
from benchmarks.units import OperationIndex
from benchmarks.units import operation_index
from compas_cgal.replay_classification import CutPlaneRampError
from compas_cgal.replay_classification import OffPlaneHelixError
from compas_cgal.replay_classification import OffPlaneReplayCurveError
from compas_cgal.replay_classification import ReplayCategory
from compas_cgal.replay_classification import classify_helix_replay
from compas_cgal.replay_classification import classify_line_replay
from compas_cgal.replay_classification import classify_planar_replay
from compas_cgal.replay_classification import line_cut_height_anchor
from compas_cgal.replay_classification import minimum_cut_height
from compas_cgal.toolpath import OperationType
from compas_cgal.units import Millimetre

ZERO_COUNT = 0

NEXT_OPERATION_OFFSET = 1

# Every cut-plane replay category, one partition cell each (5-way since Task
# 1.2), read off the Literal so a new category cannot be left out of the binding.
REPLAY_CATEGORIES: tuple[ReplayCategory, ...] = get_args(ReplayCategory)


def _indices(values: Sequence[int], operation_count: int) -> tuple[OperationIndex, ...]:
    typed = tuple(operation_index(value, operation_count=operation_count) for value in values)
    if len(set(typed)) != len(typed):
        raise InvalidHeldPathEvidenceError("operation attribution must contain unique indices.")
    return typed


class HelixHasNoMotionKindError(InvalidHeldPathEvidenceError):
    """A helix source was asked for a `MotionKind`; a survey routes a helix as an entry event, never a motion or rapid."""


SNAPSHOT_RECORD_TYPES = (HeldLineSnapshot, HeldArcSnapshot, HeldCircleSnapshot, HeldHelixSnapshot)


def _unsupported_snapshot(operation: Never) -> NoReturn:
    raise InvalidHeldPathEvidenceError(f"source snapshot contains an unsupported operation record {type(operation).__name__!r}.")


def snapshot_motion_kind(operation: HeldOperationSnapshot) -> MotionKind:
    """The survey `MotionKind` of one snapshot record; the one kind dispatch for snapshot consumers.

    Raises:
        HelixHasNoMotionKindError: *operation* is a helix (a survey entry event).
        InvalidHeldPathEvidenceError: *operation* is no exact snapshot record.
    """
    if type(operation) not in SNAPSHOT_RECORD_TYPES:
        raise InvalidHeldPathEvidenceError("source snapshot contains an unsupported operation record.")
    if isinstance(operation, HeldLineSnapshot):
        return MotionKind.LINE
    if isinstance(operation, HeldArcSnapshot):
        return MotionKind.ARC
    if isinstance(operation, HeldCircleSnapshot):
        return MotionKind.LOOP
    if isinstance(operation, HeldHelixSnapshot):
        raise HelixHasNoMotionKindError(f"source operation {operation.ordinal} is a helical entry; it has no survey motion kind.")
    _unsupported_snapshot(operation)


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
        elif isinstance(operation, HeldHelixSnapshot):
            # Independent cut-plane evidence, like a plunge foot (parity with
            # `replay_classification.infer_operation_cut_height`).
            line_anchors.append(Millimetre(float(operation.bottom_z)))
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
        if isinstance(operation, HeldHelixSnapshot):
            return classify_helix_replay(
                operation.operation,
                bottom_z=Millimetre(float(operation.bottom_z)),
                cut_z=Millimetre(cut_height),
            )
        if isinstance(operation, (HeldArcSnapshot, HeldCircleSnapshot)):
            return classify_planar_replay(
                operation.operation,
                motion_z=Millimetre(float(operation.centre.z)),
                cut_z=Millimetre(cut_height),
            )
        _unsupported_snapshot(operation)
    except (CutPlaneRampError, OffPlaneReplayCurveError, OffPlaneHelixError) as error:
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
    material_plunge_indices = _validate_ordered_indices(survey.material_plunge_indices, operation_count, name="survey material plunges")
    if not set(material_plunge_indices) <= set(plunge_indices):
        raise InvalidHeldPathEvidenceError("every plunge into material must be one of the survey's plunges.")
    helix_indices = _validate_ordered_indices([helix.index for helix in survey.helices], operation_count, name="survey helices")

    motion_set = set(motion_indices)
    rapid_set = set(rapid_indices)
    plunge_set = set(plunge_indices)
    retract_set = set(retract_indices)
    helix_set = set(helix_indices)
    if not retract_set <= rapid_set:
        raise InvalidHeldPathEvidenceError("every retained retract must have one non-cutting motion observation.")
    pure_rapid_set = rapid_set - retract_set
    categories = (motion_set, pure_rapid_set, plunge_set, retract_set, helix_set)
    for position, observed_category in enumerate(categories):
        if any(observed_category & later for later in categories[position + NEXT_OPERATION_OFFSET :]):
            raise InvalidHeldPathEvidenceError("motion, rapid, plunge, retract, and helix observations must be disjoint.")
    if set().union(*categories) != set(range(operation_count)):
        raise InvalidHeldPathEvidenceError("motion, rapid, plunge, retract, and helix observations must exactly partition the operation snapshot.")

    cut_height = _snapshot_cut_height(snapshot)
    expected_categories: dict[ReplayCategory, set[OperationIndex]] = {category: set() for category in REPLAY_CATEGORIES}
    for index, operation in enumerate(snapshot):
        expected_category = _snapshot_replay_category(operation, cut_height)
        expected_categories[expected_category].add(OperationIndex(index))
    observed_categories: dict[ReplayCategory, set[OperationIndex]] = {
        "motion": motion_set,
        "rapid": pure_rapid_set,
        "plunge": plunge_set,
        "retract": retract_set,
        "helix": helix_set,
    }
    if observed_categories != expected_categories:
        raise InvalidHeldPathEvidenceError("survey categories disagree with the source snapshot's cut-plane replay classification.")

    for motion in survey.motions:
        source = snapshot[motion.index]
        if motion.operation is not source.operation or motion.kind is not snapshot_motion_kind(source):
            raise InvalidHeldPathEvidenceError("cut-motion role or primitive kind disagrees with the source snapshot.")
    for rapid in survey.rapids:
        source = snapshot[rapid.index]
        if rapid.operation is not source.operation or rapid.kind is not snapshot_motion_kind(source):
            raise InvalidHeldPathEvidenceError("rapid-motion role or primitive kind disagrees with the source snapshot.")
    for helix in survey.helices:
        _validate_helix_binding(helix, snapshot[helix.index])


def _validate_helix_binding(helix: HelixObservation, source: HeldOperationSnapshot) -> None:
    """Bind one helix observation to its snapshot record: the same exact helix, entry, turns and pitch.

    Exact: the observed `Helix` must equal the record's exact round trip
    (`HeldHelixSnapshot.helix`: centre, radius, chart, entry, top and bottom
    Z, turns, turn sense), because the gated ramp verdict reads it; the entry
    doubles are the snapshot's exact ``P``; the pitch is a `Fraction` on both
    sides.

    Raises:
        InvalidHeldPathEvidenceError: The record is not a helix, or the
            observation's helix, ends, turns or pitch disagree with it.
    """
    if not isinstance(source, HeldHelixSnapshot):
        raise InvalidHeldPathEvidenceError(f"survey helix {helix.index} has no helix record in the source snapshot.")
    entry = (float(source.entry.x), float(source.entry.y))
    if helix.helix != source.helix():
        raise InvalidHeldPathEvidenceError(f"survey helix {helix.index} is not its source helix record's exact helix.")
    if helix.start != entry or helix.end != entry or helix.turns != source.turns or helix.pitch != source.pitch:
        raise InvalidHeldPathEvidenceError(f"survey helix {helix.index} disagrees with its source helix record (entry, turns or pitch).")
