from __future__ import annotations

import math
from dataclasses import FrozenInstanceError
from dataclasses import replace
from typing import Any
from typing import cast

import pytest
from compas.geometry import Line
from compas.geometry import Polygon
from compas.tolerance import TOL

import benchmarks.depletion as depletion
import benchmarks.quality_observations as quality_observations
from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import ContradictoryPathQualityEvidenceError
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import UnreplayableOperationError
from benchmarks.held_path_snapshot import HeldCircleSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality import CONTINUITY_TOOL_RADIUS_FRACTION
from benchmarks.quality import TANGENT_CONTINUITY_SLACK
from benchmarks.quality_observations import CRITERION_NAMES
from benchmarks.quality_observations import EVIDENCE_BY_CRITERION
from benchmarks.quality_observations import CountCriterion
from benchmarks.quality_observations import DegreesCriterion
from benchmarks.quality_observations import FractionCriterion
from benchmarks.quality_observations import MeasuredStep
from benchmarks.quality_observations import OperationPair
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import PathQualityAttribution
from benchmarks.quality_observations import ToolRadiusMultipleCriterion
from benchmarks.quality_observations import assess_path_quality
from benchmarks.spec import PocketSpec
from benchmarks.survey import EngagementSample
from benchmarks.survey import MotionKind
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.survey import RapidMotion
from benchmarks.units import Degrees
from benchmarks.units import MotionCount
from benchmarks.units import OperationIndex
from benchmarks.units import ToolRadiusMultiple
from benchmarks.units import UnitFraction
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

EXPECTED_NAMES = (
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

EXPECTED_EVIDENCE = {
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

SPEC = PocketSpec.build(
    name="quality-observations",
    family="analytic",
    polygon=Polygon([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 20.0, 0.0], [0.0, 20.0, 0.0]]),
    tool_diameter=2.0,
    tea_cap_deg=120.0,
)


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


def _sample(engagement: float, *, cap: bool = False, inside: bool = True) -> EngagementSample:
    return EngagementSample(distance=0.0, position=_point(0.0, 0.0), engagement_deg=engagement, cap_exceeded=cap, inside_centre_domain=inside)


def _motion(
    index: int,
    *,
    operation: OperationType = OperationType.CUT,
    kind: MotionKind = MotionKind.LINE,
    length: float = 1.0,
    radius: float | None = None,
    start: tuple[float, float] = (0.0, 0.0),
    end: tuple[float, float] = (1.0, 0.0),
    start_tangent: tuple[float, float] = (1.0, 0.0),
    end_tangent: tuple[float, float] = (1.0, 0.0),
    samples: tuple[EngagementSample, ...] = (),
    slot: bool = False,
    removes: bool = True,
) -> MotionQuality:
    return MotionQuality(
        index=index,
        operation=operation,
        kind=kind,
        length=length,
        loop_radius=radius,
        curvature=0.0 if radius is None else 1.0 / radius,
        swept_area=1.0,
        start=start,
        end=end,
        start_tangent=start_tangent,
        end_tangent=end_tangent,
        samples=samples,
        cap_exceeded=any(sample.cap_exceeded for sample in samples),
        slot_exceeded=slot,
        removes_material=removes,
    )


def _rapid(
    index: int,
    *,
    operation: OperationType,
    length: float,
    horizontal_at_cut_plane: bool,
    kind: MotionKind = MotionKind.LINE,
) -> RapidMotion:
    rapid = RapidMotion(
        index=index,
        operation=operation,
        kind=kind,
        length=length,
        horizontal_at_cut_plane=horizontal_at_cut_plane,
    )
    return rapid


def _snapshot(
    count: int,
    path_indices: dict[int, int] | None = None,
    operation_roles: dict[int, OperationType] | None = None,
    circle_indices: set[int] | None = None,
    heights: dict[int, float] | None = None,
    horizontal_role_indices: set[int] | None = None,
) -> tuple[HeldOperationSnapshot, ...]:
    paths = path_indices or {}
    roles = operation_roles or {}
    direction = Direction3[WorldXYZ].build(1.0, 0.0, 0.0)
    circles = circle_indices or set()
    operation_heights = heights or {}
    horizontal_roles = horizontal_role_indices or set()
    y_direction = Direction3[WorldXYZ].build(0.0, 1.0, 0.0)
    positive_z = Direction3[WorldXYZ].build(0.0, 0.0, 1.0)
    negative_z = Direction3[WorldXYZ].build(0.0, 0.0, -1.0)
    snapshots: list[HeldOperationSnapshot] = []
    for index in range(count):
        if index in circles:
            snapshots.append(
                HeldCircleSnapshot.build(
                    ordinal=OperationIndex(index),
                    operation=roles.get(index, OperationType.CUT),
                    path_index=paths.get(index, 0),
                    clockwise=False,
                    centre=Point3[WorldXYZ].build(float(index), 0.0, operation_heights.get(index, 0.0)),
                    xaxis=direction,
                    yaxis=y_direction,
                    radius=Millimetre(1.0),
                    start_tangent=y_direction,
                    end_tangent=y_direction,
                )
            )
        else:
            role = roles.get(index, OperationType.CUT)
            height = operation_heights.get(index, 0.0)
            if role is OperationType.PLUNGE and index not in horizontal_roles:
                start = Point3[WorldXYZ].build(float(index), 0.0, height + 1.0)
                end = Point3[WorldXYZ].build(float(index), 0.0, height)
                tangent = negative_z
            elif role is OperationType.RETRACT and index not in horizontal_roles:
                start = Point3[WorldXYZ].build(float(index), 0.0, height)
                end = Point3[WorldXYZ].build(float(index), 0.0, height + 1.0)
                tangent = positive_z
            else:
                start = Point3[WorldXYZ].build(float(index), 0.0, height)
                end = Point3[WorldXYZ].build(float(index + 1), 0.0, height)
                tangent = direction
            snapshots.append(
                HeldLineSnapshot.build(
                    ordinal=OperationIndex(index),
                    operation=role,
                    path_index=paths.get(index, 0),
                    clockwise=False,
                    start=start,
                    end=end,
                    start_tangent=tangent,
                    end_tangent=tangent,
                )
            )
    return tuple(snapshots)


def _survey(
    motions: tuple[MotionQuality, ...],
    rapids: tuple[RapidMotion, ...] = (),
    *,
    spec: PocketSpec = SPEC,
    plunge_indices: tuple[OperationIndex, ...] = (),
    retract_indices: tuple[OperationIndex, ...] = (),
    source_snapshot: tuple[HeldOperationSnapshot, ...],
) -> PathSurvey:
    survey = PathSurvey(
        spec=spec,
        source_snapshot=source_snapshot,
        motions=motions,
        rapids=rapids,
        plunges=len(plunge_indices),
        retracts=len(retract_indices),
        plunge_indices=plunge_indices,
        retract_indices=retract_indices,
        final_stock=cast(Stock, object()),
        total_length=sum(motion.length for motion in motions) + sum(rapid.length for rapid in rapids),
        cut_length=sum(motion.length for motion in motions),
        air_length=sum(rapid.length for rapid in rapids),
        plunge_swept_area=0.0,
    )
    return survey


def _coverage(uncut: int = 0) -> CoverageEstimate:
    return CoverageEstimate(nx=10, ny=10, cell_area=1.0, reachable_samples=10, uncut_reachable_samples=uncut, remaining_samples=uncut, wall_scallop_height=0.0)


def _expected_count_sources(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> tuple[tuple[object, ...], ...]:
    operation_count = len(snapshot)
    continuity: list[OperationPair] = []
    tangent: list[OperationPair] = []
    tolerance = CONTINUITY_TOOL_RADIUS_FRACTION * spec.tool_radius
    for previous, current in zip(survey.motions, survey.motions[1:]):
        if current.index != previous.index + 1:
            continue
        pair = OperationPair.build(
            previous=OperationIndex(previous.index),
            current=OperationIndex(current.index),
            operation_count=operation_count,
        )
        if math.hypot(current.start[0] - previous.end[0], current.start[1] - previous.end[1]) > tolerance:
            continuity.append(pair)
        dot = previous.end_tangent[0] * current.start_tangent[0] + previous.end_tangent[1] * current.start_tangent[1]
        if dot < 1.0 - TANGENT_CONTINUITY_SLACK:
            tangent.append(pair)
    return (
        tuple(OperationIndex(motion.index) for motion in survey.motions if any(not sample.inside_centre_domain for sample in motion.samples)),
        tuple(OperationIndex(rapid.index) for rapid in survey.rapids if rapid.horizontal_at_cut_plane),
        tuple(continuity),
        tuple(
            OperationIndex(index)
            for index in sorted([motion.index for motion in survey.motions if motion.length == 0.0] + [rapid.index for rapid in survey.rapids if rapid.length == 0.0])
        ),
        tuple(
            OperationIndex(motion.index)
            for motion in survey.motions
            if motion.kind is MotionKind.LOOP and motion.loop_radius is not None and motion.loop_radius <= spec.tool_radius
        ),
        tuple(OperationIndex(motion.index) for motion in survey.motions if not motion.removes_material),
        tuple(OperationIndex(motion.index) for motion in survey.motions if any(sample.cap_exceeded for sample in motion.samples)),
        tuple(OperationIndex(motion.index) for motion in survey.motions if motion.slot_exceeded),
        tuple(tangent),
    )


def _line_snapshot(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    operation: OperationType,
    *,
    ordinal: int = 0,
) -> HeldLineSnapshot:
    return HeldLineSnapshot.build(
        ordinal=OperationIndex(ordinal),
        operation=operation,
        path_index=0,
        clockwise=False,
        start=Point3[WorldXYZ].build(*start),
        end=Point3[WorldXYZ].build(*end),
        start_tangent=None,
        end_tangent=None,
    )


def _expected_maximum_pairs(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> tuple[OperationPair | None, OperationPair | None]:
    operation_count = len(snapshot)
    engagement_pair: OperationPair | None = None
    engagement_value = 0.0
    for previous, current in zip(survey.motions, survey.motions[1:]):
        if current.index != previous.index + 1:
            continue
        previous_peak = max((sample.engagement_deg for sample in previous.samples), default=0.0)
        current_peak = max((sample.engagement_deg for sample in current.samples), default=0.0)
        value = abs(current_peak - previous_peak)
        if value > engagement_value:
            engagement_value = value
            engagement_pair = OperationPair.build(
                previous=OperationIndex(previous.index),
                current=OperationIndex(current.index),
                operation_count=operation_count,
            )

    chain_by_index = {int(operation.ordinal): operation.path_index for operation in snapshot}
    rapid_indices = sorted(rapid.index for rapid in survey.rapids)
    rapid_position = 0
    current_chain: int | None = None
    previous_loop: tuple[int, float] | None = None
    loop_pair: OperationPair | None = None
    loop_value = 0.0
    for motion in survey.motions:
        boundary = False
        while rapid_position < len(rapid_indices) and rapid_indices[rapid_position] < motion.index:
            rapid_position += 1
            boundary = True
        motion_chain = chain_by_index[motion.index]
        if current_chain is not None and motion_chain != current_chain:
            boundary = True
        current_chain = motion_chain
        if boundary:
            previous_loop = None
        if motion.kind is not MotionKind.LOOP or motion.loop_radius is None:
            continue
        if previous_loop is not None:
            value = abs(motion.loop_radius - previous_loop[1]) / spec.tool_radius
            if value > loop_value:
                loop_value = value
                loop_pair = OperationPair.build(
                    previous=OperationIndex(previous_loop[0]),
                    current=OperationIndex(motion.index),
                    operation_count=operation_count,
                )
        previous_loop = (motion.index, motion.loop_radius)
    return engagement_pair, loop_pair


def _assert_exact_count_sources(spec: PocketSpec, snapshot: tuple[HeldOperationSnapshot, ...], survey: PathSurvey, assessment: PathQualityAssessment) -> None:
    attribution = assessment.attribution
    actual = (
        attribution.gouging_operations,
        attribution.unsafe_rapid_operations,
        attribution.continuity_break_pairs,
        attribution.zero_length_operations,
        attribution.degenerate_loop_operations,
        attribution.redundant_operations,
        attribution.cap_exceeded_operations,
        attribution.slotting_operations,
        attribution.tangent_break_pairs,
    )
    assert actual == _expected_count_sources(spec, snapshot, survey)
    expected_engagement, expected_loop = _expected_maximum_pairs(spec, snapshot, survey)
    actual_engagement = None if attribution.max_engagement_step is None else attribution.max_engagement_step.pair
    actual_loop = None if attribution.max_loop_radius_step is None else attribution.max_loop_radius_step.pair
    assert (actual_engagement, actual_loop) == (expected_engagement, expected_loop)


def _assess(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    coverage: CoverageEstimate,
) -> PathQualityAssessment:
    assessment = assess_path_quality(spec, snapshot, survey, coverage)
    _assert_exact_count_sources(spec, snapshot, survey, assessment)
    return assessment


@pytest.mark.parametrize(
    ("start", "end", "operation", "expected_category", "expected_replay"),
    [
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), OperationType.CUT, "motion", depletion.ReplayKind.CUT),
        ((0.0, 0.0, 2.0), (1.0, 0.0, 2.0), OperationType.LINK, "rapid", depletion.ReplayKind.RAPID),
        ((0.0, 0.0, 2.0), (0.0, 0.0, 0.0), OperationType.PLUNGE, "plunge", depletion.ReplayKind.PLUNGE),
        ((0.0, 0.0, 0.0), (0.0, 0.0, 2.0), OperationType.RETRACT, "retract", depletion.ReplayKind.RAPID),
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), OperationType.PLUNGE, "rapid", depletion.ReplayKind.RAPID),
        ((0.0, 0.0, 0.0), (0.0, 0.0, TOL.absolute), OperationType.CUT, "motion", depletion.ReplayKind.CUT),
        ((0.0, 0.0, 0.0), (0.0, 0.0, 2.0 * TOL.absolute), OperationType.CUT, "rapid", depletion.ReplayKind.RAPID),
        (
            (0.0, 0.0, 2.0 * TOL.absolute),
            (TOL.absolute, 0.0, 0.0),
            OperationType.CUT,
            "plunge",
            depletion.ReplayKind.PLUNGE,
        ),
        ((0.0, 0.0, TOL.absolute), (1.0, 0.0, TOL.absolute), OperationType.CUT, "motion", depletion.ReplayKind.CUT),
        (
            (0.0, 0.0, 2.0 * TOL.absolute),
            (1.0, 0.0, 2.0 * TOL.absolute),
            OperationType.CUT,
            "rapid",
            depletion.ReplayKind.RAPID,
        ),
    ],
)
def test_toolpath_and_snapshot_use_one_replay_classification_decision(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    operation: OperationType,
    expected_category: depletion.ReplayCategory,
    expected_replay: depletion.ReplayKind,
) -> None:
    toolpath_operation = ToolpathOperation(geometry=Line(start, end), operation=operation, path_index=0)
    snapshot = _line_snapshot(start, end, operation)

    assert depletion._replay_kind(0, toolpath_operation, 0.0) is expected_replay
    assert quality_observations._snapshot_replay_category(snapshot, 0.0) == expected_category


def test_toolpath_and_snapshot_share_invalid_ramp_classification() -> None:
    ordinal = 7
    start = (0.0, 0.0, 2.0 * TOL.absolute)
    end = (2.0 * TOL.absolute, 0.0, 0.0)
    toolpath_operation = ToolpathOperation(geometry=Line(start, end), operation=OperationType.CUT, path_index=0)
    snapshot = _line_snapshot(start, end, OperationType.CUT, ordinal=ordinal)

    with pytest.raises(UnreplayableOperationError, match=r"Operation 7 \(cut\)"):
        depletion._replay_kind(ordinal, toolpath_operation, 0.0)
    with pytest.raises(InvalidHeldPathEvidenceError) as exc_info:
        quality_observations._snapshot_replay_category(snapshot, 0.0)
    assert isinstance(exc_info.value.__cause__, UnreplayableOperationError)
    assert "Operation 7 (cut)" in str(exc_info.value)


def test_snapshot_curve_cannot_authorize_itself_below_line_cut_plane() -> None:
    snapshot = _snapshot(2, circle_indices={1}, heights={1: -1.0})
    cut_height = quality_observations._snapshot_cut_height(snapshot)

    assert cut_height == 0.0
    with pytest.raises(InvalidHeldPathEvidenceError) as exc_info:
        quality_observations._snapshot_replay_category(snapshot[1], cut_height)
    assert isinstance(exc_info.value.__cause__, UnreplayableOperationError)


def test_quality_assessment_rejects_foreign_survey_spec_binding() -> None:
    foreign = PocketSpec.build(
        name="foreign-quality-observations",
        family="analytic",
        polygon=SPEC.polygon,
        tool_diameter=SPEC.tool_diameter,
        tea_cap_deg=SPEC.tea_cap_deg,
    )
    snapshot = _snapshot(1)
    survey = _survey((_motion(0),), spec=foreign, source_snapshot=snapshot)

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, snapshot, survey, _coverage())


@pytest.mark.parametrize("case", ["omitted-operation", "reordered-operations", "wrong-primitive-kind", "wrong-operation-role"])
def test_quality_assessment_rejects_snapshot_survey_binding(case: str) -> None:
    snapshot = _snapshot(2 if case in {"omitted-operation", "reordered-operations"} else 1)
    if case == "omitted-operation":
        motions = (_motion(0),)
    elif case == "reordered-operations":
        motions = (_motion(1), _motion(0))
    elif case == "wrong-primitive-kind":
        motions = (_motion(0, kind=MotionKind.LOOP, radius=2.0),)
    else:
        motions = (_motion(0, operation=OperationType.LINK),)
    survey = _survey(motions, source_snapshot=snapshot)

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, snapshot, survey, _coverage())


def test_quality_assessment_rejects_structurally_foreign_source_snapshot_binding() -> None:
    retained = _snapshot(1)
    supplied = _snapshot(1, path_indices={0: 1})
    survey = _survey((_motion(0),), source_snapshot=retained)

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, supplied, survey, _coverage())


def test_quality_assessment_rejects_wrong_rapid_primitive_kind_binding() -> None:
    snapshot = _snapshot(2, operation_roles={1: OperationType.LINK}, heights={1: 2.0})
    rapid = _rapid(index=1, operation=OperationType.LINK, kind=MotionKind.LOOP, length=1.0, horizontal_at_cut_plane=False)
    survey = _survey((_motion(0),), (rapid,), source_snapshot=snapshot)

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, snapshot, survey, _coverage())


def test_quality_assessment_rejects_cut_plane_cut_classified_as_rapid_binding() -> None:
    snapshot = _snapshot(1)
    rapid = _rapid(index=0, operation=OperationType.CUT, length=1.0, horizontal_at_cut_plane=True)
    survey = _survey((), (rapid,), source_snapshot=snapshot)

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, snapshot, survey, _coverage())


def test_quality_assessment_rejects_clearance_link_classified_as_motion_binding() -> None:
    snapshot = _snapshot(2, operation_roles={1: OperationType.LINK}, heights={1: 2.0})
    survey = _survey(
        (_motion(0), _motion(1, operation=OperationType.LINK)),
        source_snapshot=snapshot,
    )

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, snapshot, survey, _coverage())


def test_quality_assessment_requires_complete_motion_rapid_plunge_retract_partition() -> None:
    snapshot = _snapshot(
        4,
        operation_roles={
            0: OperationType.PLUNGE,
            2: OperationType.RETRACT,
            3: OperationType.LINK,
        },
        heights={3: 2.0},
    )
    rapids = (
        _rapid(index=2, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),
        _rapid(index=3, operation=OperationType.LINK, length=1.0, horizontal_at_cut_plane=False),
    )
    incomplete = _survey(
        (_motion(1),),
        rapids,
        source_snapshot=snapshot,
        retract_indices=(OperationIndex(2),),
    )

    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, snapshot, incomplete, _coverage())


def test_quality_assessment_accepts_exact_motion_rapid_plunge_retract_partition() -> None:
    snapshot = _snapshot(
        4,
        operation_roles={
            0: OperationType.PLUNGE,
            2: OperationType.RETRACT,
            3: OperationType.LINK,
        },
        heights={3: 2.0},
    )
    rapids = (
        _rapid(index=2, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),
        _rapid(index=3, operation=OperationType.LINK, length=1.0, horizontal_at_cut_plane=False),
    )
    survey = _survey(
        (_motion(1),),
        rapids,
        source_snapshot=snapshot,
        plunge_indices=(OperationIndex(0),),
        retract_indices=(OperationIndex(2),),
    )

    assessment = _assess(SPEC, snapshot, survey, _coverage())

    assert survey.plunge_indices == (OperationIndex(0),)
    assert survey.retract_indices == (OperationIndex(2),)
    assert assessment.attribution.operation_count == 4


def test_criterion_vocabulary_thresholds_evidence_and_outcomes_are_closed() -> None:
    assert CRITERION_NAMES == EXPECTED_NAMES
    assert EVIDENCE_BY_CRITERION == EXPECTED_EVIDENCE

    fraction = FractionCriterion.build(name="uncut fraction", measured=UnitFraction(0.0), required=UnitFraction(0.0), evidence="sampled_diagnostic", outcome="no_failure_observed")
    count = CountCriterion.build(
        name="gouging motions", measured=MotionCount(0), required=MotionCount(0), evidence="sampled_diagnostic", outcome="no_failure_observed", attribution_count=0
    )
    degrees = DegreesCriterion.build(
        name="max engagement step (deg)", measured=Degrees(120.0), required=Degrees(120.0), evidence="sampled_diagnostic", outcome="no_failure_observed"
    )
    radii = ToolRadiusMultipleCriterion.build(
        name="max loop radius step (tool radii)", measured=ToolRadiusMultiple(2.0), required=ToolRadiusMultiple(2.0), evidence="derived_geometry", outcome="criterion_satisfied"
    )

    assert (fraction.required, count.required, degrees.required, radii.required) == (0.0, 0, 120.0, 2.0)
    with pytest.raises(FrozenInstanceError):
        fraction.measured = UnitFraction(1.0)  # type: ignore[misc]


@pytest.mark.parametrize(
    ("factory", "kwargs"),
    [
        (
            FractionCriterion.build,
            {"name": "uncut fraction", "measured": UnitFraction(math.nan), "required": UnitFraction(0.0), "evidence": "sampled_diagnostic", "outcome": "failure_observed"},
        ),
        (
            FractionCriterion.build,
            {"name": "uncut fraction", "measured": UnitFraction(0.0), "required": UnitFraction(0.0), "evidence": "derived_geometry", "outcome": "criterion_satisfied"},
        ),
        (
            FractionCriterion.build,
            {
                "name": "max engagement step (deg)",
                "measured": UnitFraction(0.0),
                "required": UnitFraction(0.0),
                "evidence": "sampled_diagnostic",
                "outcome": "no_failure_observed",
            },
        ),
        (
            CountCriterion.build,
            {
                "name": "gouging motions",
                "measured": MotionCount(-1),
                "required": MotionCount(0),
                "evidence": "sampled_diagnostic",
                "outcome": "failure_observed",
                "attribution_count": 0,
            },
        ),
        (
            CountCriterion.build,
            {
                "name": "gouging motions",
                "measured": MotionCount(1),
                "required": MotionCount(0),
                "evidence": "sampled_diagnostic",
                "outcome": "failure_observed",
                "attribution_count": 0,
            },
        ),
        (
            DegreesCriterion.build,
            {"name": "max engagement step (deg)", "measured": Degrees(121.0), "required": Degrees(120.0), "evidence": "sampled_diagnostic", "outcome": "no_failure_observed"},
        ),
        (
            ToolRadiusMultipleCriterion.build,
            {
                "name": "max loop radius step (tool radii)",
                "measured": ToolRadiusMultiple(-1.0),
                "required": ToolRadiusMultiple(2.0),
                "evidence": "derived_geometry",
                "outcome": "criterion_satisfied",
            },
        ),
    ],
)
def test_criterion_factories_reject_invalid_evidence(factory: object, kwargs: dict[str, object]) -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        cast(object, factory)(**kwargs)  # type: ignore[operator]


def test_operation_pairs_reject_invalid_indices() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        OperationPair.build(previous=OperationIndex(2), current=OperationIndex(2), operation_count=3)
    with pytest.raises(InvalidHeldPathEvidenceError):
        OperationPair.build(previous=OperationIndex(1), current=OperationIndex(3), operation_count=3)


def test_measured_step_factory_rejects_open_unit_discriminator() -> None:
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)

    with pytest.raises(InvalidHeldPathEvidenceError):
        MeasuredStep.build(value=ToolRadiusMultiple(1.0), pair=pair, unit=cast(Any, "seconds"))


def test_measured_step_factory_retains_closed_unit_discriminator() -> None:
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)
    degrees = MeasuredStep.build(value=Degrees(1.0), pair=pair, unit="degrees")
    radius = MeasuredStep.build(value=ToolRadiusMultiple(1.0), pair=pair, unit="tool_radius_multiple")

    assert degrees.unit == "degrees"
    assert radius.unit == "tool_radius_multiple"


def test_attribution_factory_rejects_duplicate_indices() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        PathQualityAttribution.build(
            operation_count=2,
            uncut_operations=(),
            gouging_operations=(OperationIndex(0), OperationIndex(0)),
            unsafe_rapid_operations=(),
            continuity_break_pairs=(),
            zero_length_operations=(),
            degenerate_loop_operations=(),
            redundant_operations=(),
            cap_exceeded_operations=(),
            slotting_operations=(),
            max_engagement_step=None,
            engagement_step_failure_pairs=(),
            max_loop_radius_step=None,
            loop_radius_step_failure_pairs=(),
            tangent_break_pairs=(),
            curvature_break_pairs=(),
            reversal_pairs=(),
        )


def _attribution(
    *,
    operation_count: int = 2,
    gouging_operations: tuple[OperationIndex, ...] = (),
    max_engagement_step: MeasuredStep[Degrees] | None = None,
    engagement_step_failure_pairs: tuple[OperationPair, ...] = (),
    max_loop_radius_step: MeasuredStep[ToolRadiusMultiple] | None = None,
    loop_radius_step_failure_pairs: tuple[OperationPair, ...] = (),
    tangent_break_pairs: tuple[OperationPair, ...] = (),
    curvature_break_pairs: tuple[OperationPair, ...] = (),
    reversal_pairs: tuple[OperationPair, ...] = (),
) -> PathQualityAttribution:
    return PathQualityAttribution.build(
        operation_count=operation_count,
        uncut_operations=(),
        gouging_operations=gouging_operations,
        unsafe_rapid_operations=(),
        continuity_break_pairs=(),
        zero_length_operations=(),
        degenerate_loop_operations=(),
        redundant_operations=(),
        cap_exceeded_operations=(),
        slotting_operations=(),
        max_engagement_step=max_engagement_step,
        engagement_step_failure_pairs=engagement_step_failure_pairs,
        max_loop_radius_step=max_loop_radius_step,
        loop_radius_step_failure_pairs=loop_radius_step_failure_pairs,
        tangent_break_pairs=tangent_break_pairs,
        curvature_break_pairs=curvature_break_pairs,
        reversal_pairs=reversal_pairs,
    )


def _rebuild_attribution(
    assessment: PathQualityAssessment,
    *,
    operation_count: int,
    **changes: object,
) -> PathQualityAttribution:
    values = {name: value for name, value in vars(assessment.attribution).items() if name != "operation_count"}
    values.update(changes)
    return PathQualityAttribution.build(operation_count=operation_count, **values)  # type: ignore[arg-type]


def _rebuild_assessment(
    assessment: PathQualityAssessment,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
    *,
    spec: PocketSpec,
    attribution: PathQualityAttribution,
    **changes: object,
) -> PathQualityAssessment:
    values = {
        name: getattr(assessment, name)
        for name in (
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
    }
    values.update(changes)
    return PathQualityAssessment.build(
        spec=spec,
        snapshot=snapshot,
        survey=survey,
        attribution=attribution,
        **values,  # type: ignore[arg-type]
    )


def _forge_record(record: object, **changes: object) -> object:
    forged = object.__new__(type(record))
    for name, value in vars(record).items():
        object.__setattr__(forged, name, value)
    for name, value in changes.items():
        object.__setattr__(forged, name, value)
    return forged


def test_attribution_operation_count_is_an_exact_non_boolean_integer() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(operation_count=cast(Any, True))
    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(operation_count=-1)


def test_zero_operation_attribution_is_valid_only_when_wholly_empty() -> None:
    empty = _attribution(operation_count=0)
    assert empty.operation_count == 0

    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(operation_count=0, gouging_operations=(OperationIndex(0),))


def test_zero_operation_stream_is_ineligible_for_quality_assessment() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        assess_path_quality(SPEC, (), _survey((), source_snapshot=()), _coverage())


@pytest.mark.parametrize(
    ("slot", "unit"),
    [("engagement", "tool_radius_multiple"), ("loop", "degrees")],
)
def test_attribution_factory_revalidates_measured_step_unit_on_slot_installation(slot: str, unit: str) -> None:
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)
    forged = object.__new__(MeasuredStep)
    object.__setattr__(forged, "value", 1.0)
    object.__setattr__(forged, "pair", pair)
    object.__setattr__(forged, "unit", unit)

    with pytest.raises(InvalidHeldPathEvidenceError):
        if slot == "engagement":
            _attribution(max_engagement_step=cast(Any, forged))
        else:
            _attribution(max_loop_radius_step=cast(Any, forged))


def test_attribution_factory_revalidates_forged_same_unit_nan_step() -> None:
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)
    forged = MeasuredStep.build(value=Degrees(1.0), pair=pair, unit="degrees")
    object.__setattr__(forged, "value", math.nan)

    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(max_engagement_step=forged)


@pytest.mark.parametrize(
    ("field", "changes", "error"),
    [
        ("uncut_fraction", {"required": UnitFraction(1.0)}, InvalidHeldPathEvidenceError),
        ("gouging_motions", {"required": MotionCount(1)}, InvalidHeldPathEvidenceError),
        ("max_engagement_step", {"required": Degrees(119.0)}, ContradictoryPathQualityEvidenceError),
        ("max_loop_radius_step", {"required": ToolRadiusMultiple(1.0)}, InvalidHeldPathEvidenceError),
        ("uncut_fraction", {"evidence": "derived_geometry"}, InvalidHeldPathEvidenceError),
        ("gouging_motions", {"evidence": "derived_geometry"}, InvalidHeldPathEvidenceError),
        ("max_engagement_step", {"evidence": "derived_geometry"}, InvalidHeldPathEvidenceError),
        ("max_loop_radius_step", {"evidence": "sampled_diagnostic"}, InvalidHeldPathEvidenceError),
        ("uncut_fraction", {"outcome": "failure_observed"}, InvalidHeldPathEvidenceError),
        ("gouging_motions", {"outcome": "failure_observed"}, InvalidHeldPathEvidenceError),
        ("max_engagement_step", {"outcome": "failure_observed"}, InvalidHeldPathEvidenceError),
        ("max_loop_radius_step", {"outcome": "criterion_violated"}, InvalidHeldPathEvidenceError),
    ],
)
def test_assessment_rejects_forged_installed_criterion_contradiction(
    field: str,
    changes: dict[str, object],
    error: type[Exception],
) -> None:
    snapshot = _snapshot(2, circle_indices={0, 1})
    survey = _survey(
        (
            _motion(0, kind=MotionKind.LOOP, radius=1.0, samples=(_sample(0.0),)),
            _motion(1, kind=MotionKind.LOOP, radius=2.0, samples=(_sample(10.0),)),
        ),
        source_snapshot=snapshot,
    )
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    forged = _forge_record(getattr(assessment, field), **changes)

    with pytest.raises(error):
        _rebuild_assessment(
            assessment,
            snapshot,
            survey,
            spec=SPEC,
            attribution=assessment.attribution,
            **{field: forged},
        )


def test_assessment_rejects_installed_count_criterion_field_identity_contradiction() -> None:
    snapshot = _snapshot(1)
    survey = _survey((_motion(0),), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(
            assessment,
            snapshot,
            survey,
            spec=SPEC,
            attribution=assessment.attribution,
            gouging_motions=assessment.unsafe_rapids,
        )


def test_direct_assessment_builder_rejects_foreign_spec_binding() -> None:
    snapshot = _snapshot(1)
    survey = _survey((_motion(0),), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    foreign = PocketSpec.build(
        name="foreign-direct-builder",
        family=SPEC.family,
        polygon=SPEC.polygon,
        tool_diameter=SPEC.tool_diameter,
        tea_cap_deg=SPEC.tea_cap_deg,
    )
    values = {
        name: getattr(assessment, name)
        for name in (
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
    }

    with pytest.raises(InvalidHeldPathEvidenceError):
        PathQualityAssessment.build(
            spec=foreign,
            snapshot=snapshot,
            survey=survey,
            attribution=assessment.attribution,
            **values,  # type: ignore[arg-type,call-arg]
        )


def test_assessment_rejects_contradiction_of_violated_engagement_maximum_without_failure_pair() -> None:
    snapshot = _snapshot(2)
    survey = _survey((_motion(0, samples=(_sample(0.0),)), _motion(1, samples=(_sample(130.0),))), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    attribution = _rebuild_attribution(assessment, operation_count=2, engagement_step_failure_pairs=())

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(assessment, snapshot, survey, spec=SPEC, attribution=attribution)


def test_assessment_rejects_contradiction_of_satisfied_engagement_maximum_with_failure_pair() -> None:
    snapshot = _snapshot(2)
    survey = _survey((_motion(0, samples=(_sample(0.0),)), _motion(1, samples=(_sample(10.0),))), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    assert assessment.attribution.max_engagement_step is not None
    attribution = _rebuild_attribution(
        assessment,
        operation_count=2,
        engagement_step_failure_pairs=(assessment.attribution.max_engagement_step.pair,),
    )

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(assessment, snapshot, survey, spec=SPEC, attribution=attribution)


def test_assessment_rejects_contradiction_of_violated_maximum_absent_from_failure_pairs() -> None:
    snapshot = _snapshot(3)
    survey = _survey(
        (
            _motion(0, samples=(_sample(0.0),)),
            _motion(1, samples=(_sample(130.0),)),
            _motion(2, samples=(_sample(0.0),)),
        ),
        source_snapshot=snapshot,
    )
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    assert assessment.attribution.engagement_step_failure_pairs
    attribution = _rebuild_attribution(
        assessment,
        operation_count=3,
        engagement_step_failure_pairs=(assessment.attribution.engagement_step_failure_pairs[-1],),
    )

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(assessment, snapshot, survey, spec=SPEC, attribution=attribution)


def test_assessment_rejects_contradiction_of_in_bounds_nonadjacent_engagement_pair() -> None:
    snapshot = _snapshot(3)
    survey = _survey((_motion(0), _motion(1), _motion(2)), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(2), operation_count=3)
    maximum = MeasuredStep.build(value=Degrees(1.0), pair=pair, unit="degrees")
    attribution = _rebuild_attribution(assessment, operation_count=3, max_engagement_step=maximum)
    criterion = DegreesCriterion.build(
        name="max engagement step (deg)",
        measured=Degrees(1.0),
        required=Degrees(SPEC.tea_cap_deg),
        evidence="sampled_diagnostic",
        outcome="no_failure_observed",
    )

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(
            assessment,
            snapshot,
            survey,
            spec=SPEC,
            attribution=attribution,
            max_engagement_step=criterion,
        )


@pytest.mark.parametrize("boundary", ["rapid", "path-chain"])
def test_assessment_rejects_contradiction_of_loop_pair_crossing_rapid_or_path_chain_boundary(boundary: str) -> None:
    paths = {0: 0, 2: 0} if boundary == "rapid" else {0: 0, 1: 1}
    roles = {1: OperationType.RETRACT} if boundary == "rapid" else {}
    last_index = 2 if boundary == "rapid" else 1
    snapshot = _snapshot(3 if boundary == "rapid" else 2, paths, roles, circle_indices={0, last_index})
    motions = (
        _motion(0, kind=MotionKind.LOOP, radius=1.0),
        _motion(2 if boundary == "rapid" else 1, kind=MotionKind.LOOP, radius=2.0),
    )
    rapids = (_rapid(index=1, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),) if boundary == "rapid" else ()
    retracts = (OperationIndex(1),) if boundary == "rapid" else ()
    survey = _survey(motions, rapids, source_snapshot=snapshot, retract_indices=retracts)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    pair = OperationPair.build(
        previous=OperationIndex(0),
        current=OperationIndex(2 if boundary == "rapid" else 1),
        operation_count=len(snapshot),
    )
    maximum = MeasuredStep.build(value=ToolRadiusMultiple(1.0), pair=pair, unit="tool_radius_multiple")
    attribution = _rebuild_attribution(assessment, operation_count=len(snapshot), max_loop_radius_step=maximum)
    criterion = ToolRadiusMultipleCriterion.build(
        name="max loop radius step (tool radii)",
        measured=ToolRadiusMultiple(1.0),
        required=ToolRadiusMultiple(2.0),
        evidence="derived_geometry",
        outcome="criterion_satisfied",
    )

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(
            assessment,
            snapshot,
            survey,
            spec=SPEC,
            attribution=attribution,
            max_loop_radius_step=criterion,
        )


@pytest.mark.parametrize(
    ("field", "wrong_sources"),
    [
        ("gouging_operations", (OperationIndex(1),)),
        ("unsafe_rapid_operations", (OperationIndex(2),)),
        (
            "continuity_break_pairs",
            (OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=4),),
        ),
        ("zero_length_operations", (OperationIndex(1),)),
        ("degenerate_loop_operations", (OperationIndex(1),)),
        ("redundant_operations", (OperationIndex(1),)),
        ("cap_exceeded_operations", (OperationIndex(1),)),
        ("slotting_operations", (OperationIndex(1),)),
        (
            "tangent_break_pairs",
            (OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=4),),
        ),
    ],
)
def test_assessment_rejects_same_cardinality_wrong_exact_sources(
    field: str,
    wrong_sources: tuple[object, ...],
) -> None:
    snapshot = _snapshot(
        4,
        operation_roles={3: OperationType.PLUNGE},
        circle_indices={0},
        horizontal_role_indices={3},
    )
    motions = (
        _motion(
            0,
            kind=MotionKind.LOOP,
            length=0.0,
            radius=SPEC.tool_radius,
            end=(0.0, 0.0),
            samples=(_sample(130.0, cap=True, inside=False),),
            slot=True,
            removes=False,
        ),
        _motion(1, start=(2.0, 0.0), end=(1.0, 0.0), start_tangent=(0.0, 1.0)),
        _motion(2, start=(1.0, 0.0)),
    )
    rapid = _rapid(index=3, operation=OperationType.PLUNGE, length=1.0, horizontal_at_cut_plane=True)
    survey = _survey(motions, (rapid,), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    attribution = _rebuild_attribution(
        assessment,
        operation_count=4,
        **{field: wrong_sources},
    )

    with pytest.raises(ContradictoryPathQualityEvidenceError):
        _rebuild_assessment(assessment, snapshot, survey, spec=SPEC, attribution=attribution)


def test_attribution_factory_revalidates_engagement_maximum_pair_bounds() -> None:
    larger_stream_pair = OperationPair.build(previous=OperationIndex(2), current=OperationIndex(3), operation_count=4)
    maximum = MeasuredStep[Degrees].build(value=Degrees(1.0), pair=larger_stream_pair, unit="degrees")

    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(max_engagement_step=maximum)


def test_attribution_factory_revalidates_loop_maximum_pair_bounds() -> None:
    larger_stream_pair = OperationPair.build(previous=OperationIndex(2), current=OperationIndex(3), operation_count=4)
    maximum = MeasuredStep[ToolRadiusMultiple].build(value=ToolRadiusMultiple(1.0), pair=larger_stream_pair, unit="tool_radius_multiple")

    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(max_loop_radius_step=maximum)


def test_attribution_factory_rejects_reversal_without_tangent_break() -> None:
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)

    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(reversal_pairs=(pair,))


def test_attribution_factory_rejects_tangent_and_curvature_overlap() -> None:
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)

    with pytest.raises(InvalidHeldPathEvidenceError):
        _attribution(tangent_break_pairs=(pair,), curvature_break_pairs=(pair,))


def test_simple_findings_retain_exact_unique_operation_indices() -> None:
    motions = (
        _motion(0, length=0.0, samples=(_sample(300.0, cap=False, inside=False),), removes=False),
        _motion(2, kind=MotionKind.LOOP, radius=SPEC.tool_radius, samples=(_sample(10.0, cap=True),), slot=True),
        _motion(4, length=0.0),
    )
    rapids = (
        _rapid(index=1, operation=OperationType.PLUNGE, length=0.0, horizontal_at_cut_plane=True),
        _rapid(index=3, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),
    )
    snapshot = _snapshot(
        5,
        operation_roles={1: OperationType.PLUNGE, 3: OperationType.RETRACT},
        circle_indices={2},
        horizontal_role_indices={1},
    )
    survey = _survey(motions, rapids, source_snapshot=snapshot, retract_indices=(OperationIndex(3),))
    assessment = _assess(SPEC, snapshot, survey, _coverage(uncut=2))

    assert assessment.uncut_fraction.measured == UnitFraction(0.2)
    assert assessment.attribution.uncut_operations == ()
    assert assessment.attribution.gouging_operations == (OperationIndex(0),)
    assert assessment.attribution.unsafe_rapid_operations == (OperationIndex(1),)
    assert assessment.attribution.zero_length_operations == (OperationIndex(0), OperationIndex(1), OperationIndex(4))
    assert assessment.attribution.degenerate_loop_operations == (OperationIndex(2),)
    assert assessment.attribution.redundant_operations == (OperationIndex(0),)
    assert assessment.attribution.cap_exceeded_operations == (OperationIndex(2),)
    assert assessment.attribution.slotting_operations == (OperationIndex(2),)
    assert assessment.cap_exceedances.measured == MotionCount(1)
    assert tuple(
        criterion.required
        for criterion in (
            assessment.uncut_fraction,
            assessment.gouging_motions,
            assessment.unsafe_rapids,
            assessment.continuity_breaks,
            assessment.zero_length_motions,
            assessment.degenerate_loops,
            assessment.redundant_operations,
            assessment.cap_exceedances,
            assessment.slotting_motions,
            assessment.max_engagement_step,
            assessment.max_loop_radius_step,
            assessment.tangent_breaks,
        )
    ) == (0.0, 0, 0, 0, 0, 0, 0, 0, 0, SPEC.tea_cap_deg, 2.0, 0)


def _source_attributions(assessment: PathQualityAssessment) -> dict[str, tuple[object, ...]]:
    attribution = assessment.attribution
    return {
        "gouging_operations": attribution.gouging_operations,
        "unsafe_rapid_operations": attribution.unsafe_rapid_operations,
        "continuity_break_pairs": attribution.continuity_break_pairs,
        "zero_length_operations": attribution.zero_length_operations,
        "degenerate_loop_operations": attribution.degenerate_loop_operations,
        "redundant_operations": attribution.redundant_operations,
        "cap_exceeded_operations": attribution.cap_exceeded_operations,
        "slotting_operations": attribution.slotting_operations,
        "tangent_break_pairs": attribution.tangent_break_pairs,
    }


@pytest.mark.parametrize(
    ("mutation", "expected_field"),
    [
        ("inside-centre-domain", "gouging_operations"),
        ("horizontal-at-cut-plane", "unsafe_rapid_operations"),
        ("endpoint-displacement", "continuity_break_pairs"),
        ("zero-length", "zero_length_operations"),
        ("loop-radius", "degenerate_loop_operations"),
        ("removes-material", "redundant_operations"),
        ("sample-cap", "cap_exceeded_operations"),
        ("slot-exceeded", "slotting_operations"),
        ("tangent-dot", "tangent_break_pairs"),
    ],
)
def test_raw_observation_sensitivity_changes_exactly_one_count_source(mutation: str, expected_field: str) -> None:
    snapshot = _snapshot(
        4,
        operation_roles={3: OperationType.PLUNGE},
        circle_indices={0},
        horizontal_role_indices={3},
    )
    motions = [
        _motion(0, kind=MotionKind.LOOP, radius=2.0, end=(0.0, 0.0), samples=(_sample(10.0),)),
        _motion(1, start=(0.0, 0.0), end=(1.0, 0.0), samples=(_sample(10.0),)),
        _motion(2, start=(1.0, 0.0), samples=(_sample(10.0),)),
    ]
    rapids = [_rapid(index=3, operation=OperationType.PLUNGE, length=1.0, horizontal_at_cut_plane=False)]
    baseline = _assess(SPEC, snapshot, _survey(tuple(motions), tuple(rapids), source_snapshot=snapshot), _coverage())

    if mutation == "inside-centre-domain":
        motions[2] = replace(motions[2], samples=(replace(motions[2].samples[0], inside_centre_domain=False),))
    elif mutation == "horizontal-at-cut-plane":
        rapids[0] = replace(rapids[0], horizontal_at_cut_plane=True)
    elif mutation == "endpoint-displacement":
        motions[1] = replace(motions[1], start=(2.0, 0.0))
    elif mutation == "zero-length":
        motions[2] = replace(motions[2], length=0.0)
    elif mutation == "loop-radius":
        motions[0] = replace(motions[0], loop_radius=SPEC.tool_radius)
    elif mutation == "removes-material":
        motions[2] = replace(motions[2], removes_material=False)
    elif mutation == "sample-cap":
        motions[2] = replace(motions[2], samples=(replace(motions[2].samples[0], cap_exceeded=True),))
    elif mutation == "slot-exceeded":
        motions[2] = replace(motions[2], slot_exceeded=True)
    else:
        motions[1] = replace(motions[1], start_tangent=(-1.0, 0.0))

    changed = _assess(SPEC, snapshot, _survey(tuple(motions), tuple(rapids), source_snapshot=snapshot), _coverage())
    baseline_sources = _source_attributions(baseline)
    changed_sources = _source_attributions(changed)
    changed_fields = {name for name in baseline_sources if baseline_sources[name] != changed_sources[name]}
    assert changed_fields == {expected_field}
    assert len(set(baseline_sources[expected_field]) ^ set(changed_sources[expected_field])) == 1


def test_nonleading_engagement_sample_changes_the_causal_maximum_pair_and_value() -> None:
    snapshot = _snapshot(3)
    baseline_motions = (
        _motion(0, samples=(_sample(0.0),)),
        _motion(1, samples=(_sample(20.0),)),
        _motion(2, samples=(_sample(15.0), _sample(15.0))),
    )
    changed_motions = (*baseline_motions[:2], replace(baseline_motions[2], samples=(_sample(15.0), _sample(50.0))))

    baseline = _assess(SPEC, snapshot, _survey(baseline_motions, source_snapshot=snapshot), _coverage())
    changed = _assess(SPEC, snapshot, _survey(changed_motions, source_snapshot=snapshot), _coverage())

    assert baseline.attribution.max_engagement_step == MeasuredStep.build(
        value=Degrees(20.0),
        pair=OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=3),
        unit="degrees",
    )
    assert changed.attribution.max_engagement_step == MeasuredStep.build(
        value=Degrees(30.0),
        pair=OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=3),
        unit="degrees",
    )


@pytest.mark.parametrize("mutation", ["radius", "path-boundary"])
def test_loop_radius_or_path_boundary_changes_the_causal_maximum_pair(mutation: str) -> None:
    baseline_snapshot = _snapshot(3, circle_indices={0, 1, 2})
    baseline_motions = (
        _motion(0, kind=MotionKind.LOOP, radius=2.0),
        _motion(1, kind=MotionKind.LOOP, radius=5.0),
        _motion(2, kind=MotionKind.LOOP, radius=6.0),
    )
    baseline = _assess(SPEC, baseline_snapshot, _survey(baseline_motions, source_snapshot=baseline_snapshot), _coverage())

    if mutation == "radius":
        changed_snapshot = baseline_snapshot
        changed_motions = (baseline_motions[0], replace(baseline_motions[1], loop_radius=3.0), baseline_motions[2])
    else:
        changed_snapshot = _snapshot(3, path_indices={1: 1}, circle_indices={0, 1, 2})
        changed_motions = baseline_motions
    changed = _assess(SPEC, changed_snapshot, _survey(changed_motions, source_snapshot=changed_snapshot), _coverage())

    assert baseline.attribution.max_loop_radius_step == MeasuredStep.build(
        value=ToolRadiusMultiple(3.0),
        pair=OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=3),
        unit="tool_radius_multiple",
    )
    if mutation == "radius":
        assert changed.attribution.max_loop_radius_step == MeasuredStep.build(
            value=ToolRadiusMultiple(3.0),
            pair=OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=3),
            unit="tool_radius_multiple",
        )
    else:
        assert changed.attribution.max_loop_radius_step is None


def test_junction_findings_preserve_adjacency_and_classification() -> None:
    tolerance = CONTINUITY_TOOL_RADIUS_FRACTION * SPEC.tool_radius
    motions = (
        _motion(0, end=(0.0, 0.0), end_tangent=(1.0, 0.0)),
        _motion(1, start=(tolerance, 0.0), end_tangent=(-1.0, 0.0), samples=(_sample(10.0),)),
        _motion(2, start=(2.0 + tolerance, 0.0), start_tangent=(1.0, 0.0), end_tangent=(1.0, 0.0), samples=(_sample(40.0),)),
        _motion(4, start=(99.0, 0.0), start_tangent=(0.0, 1.0), samples=(_sample(80.0),)),
    )
    snapshot = _snapshot(5, operation_roles={3: OperationType.LINK}, heights={3: 2.0})
    survey = _survey(motions, (_rapid(index=3, operation=OperationType.LINK, length=1.0, horizontal_at_cut_plane=False),), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())

    assert assessment.attribution.continuity_break_pairs == (OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=5),)
    assert assessment.attribution.reversal_pairs == (OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=5),)
    assert assessment.attribution.tangent_break_pairs == assessment.attribution.reversal_pairs
    assert assessment.attribution.curvature_break_pairs == ()


def test_engagement_steps_keep_first_maximum_tie_and_every_failure_pair() -> None:
    motions = (
        _motion(0, samples=(_sample(0.0),)),
        _motion(1, samples=(_sample(130.0),)),
        _motion(2, samples=(_sample(0.0),)),
        _motion(4, samples=(_sample(400.0),)),
    )
    snapshot = _snapshot(5, operation_roles={3: OperationType.LINK}, heights={3: 2.0})
    survey = _survey(motions, (_rapid(index=3, operation=OperationType.LINK, length=1.0, horizontal_at_cut_plane=False),), source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    pair_01 = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=5)
    pair_12 = OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=5)

    assert assessment.max_engagement_step.measured == Degrees(130.0)
    assert assessment.attribution.max_engagement_step is not None
    assert assessment.attribution.max_engagement_step.pair == pair_01
    assert assessment.attribution.engagement_step_failure_pairs == (pair_01, pair_12)


def test_zero_steps_have_no_invented_maximum_pair() -> None:
    motions = (_motion(0, samples=(_sample(20.0),)), _motion(1, samples=(_sample(20.0),)))
    snapshot = _snapshot(2)
    survey = _survey(motions, source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())

    assert assessment.max_engagement_step.measured == Degrees(0.0)
    assert assessment.attribution.max_engagement_step is None


def test_non_reversal_tangent_break_excludes_curvature_break() -> None:
    motions = (
        _motion(0, radius=2.0, end_tangent=(1.0, 0.0)),
        _motion(1, radius=3.0, start_tangent=(0.0, 1.0)),
    )
    snapshot = _snapshot(2)
    survey = _survey(motions, source_snapshot=snapshot)
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    expected = (OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2),)

    assert assessment.attribution.tangent_break_pairs == expected
    assert assessment.attribution.reversal_pairs == ()
    assert assessment.attribution.curvature_break_pairs == ()


def test_loop_steps_split_on_rapids_and_path_indices() -> None:
    motions = (
        _motion(0, kind=MotionKind.LOOP, radius=1.0),
        _motion(1, kind=MotionKind.LINE),
        _motion(2, kind=MotionKind.LOOP, radius=4.0),
        _motion(4, kind=MotionKind.LOOP, radius=9.0),
        _motion(5, kind=MotionKind.LOOP, radius=1.0),
        _motion(6, kind=MotionKind.LOOP, radius=4.0),
        _motion(7, kind=MotionKind.LOOP, radius=1.0),
    )
    rapids = (_rapid(index=3, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),)
    snapshot = _snapshot(
        8,
        {0: 0, 1: 0, 2: 0, 4: 0, 5: 1, 6: 1, 7: 1},
        {3: OperationType.RETRACT},
        circle_indices={0, 2, 4, 5, 6, 7},
    )
    survey = _survey(motions, rapids, source_snapshot=snapshot, retract_indices=(OperationIndex(3),))
    assessment = _assess(SPEC, snapshot, survey, _coverage())
    pair_02 = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(2), operation_count=8)
    pair_56 = OperationPair.build(previous=OperationIndex(5), current=OperationIndex(6), operation_count=8)
    pair_67 = OperationPair.build(previous=OperationIndex(6), current=OperationIndex(7), operation_count=8)

    assert assessment.max_loop_radius_step.measured == ToolRadiusMultiple(3.0)
    assert assessment.attribution.max_loop_radius_step is not None
    assert assessment.attribution.max_loop_radius_step.pair == pair_02
    assert assessment.attribution.loop_radius_step_failure_pairs == (pair_02, pair_56, pair_67)
