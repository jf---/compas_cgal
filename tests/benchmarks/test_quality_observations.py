from __future__ import annotations

import math
from dataclasses import FrozenInstanceError
from typing import cast

import pytest
from compas.geometry import Polygon

from benchmarks.coverage import CoverageEstimate
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.quality import CONTINUITY_TOOL_RADIUS_FRACTION
from benchmarks.quality_observations import CRITERION_NAMES
from benchmarks.quality_observations import EVIDENCE_BY_CRITERION
from benchmarks.quality_observations import CountCriterion
from benchmarks.quality_observations import DegreesCriterion
from benchmarks.quality_observations import FractionCriterion
from benchmarks.quality_observations import MeasuredStep
from benchmarks.quality_observations import OperationPair
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
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType

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
        operation=OperationType.CUT,
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


def _snapshot(count: int, path_indices: dict[int, int] | None = None) -> tuple[HeldLineSnapshot, ...]:
    paths = path_indices or {}
    direction = Direction3[WorldXYZ].build(1.0, 0.0, 0.0)
    return tuple(
        HeldLineSnapshot.build(
            ordinal=OperationIndex(index),
            operation=OperationType.CUT,
            path_index=paths.get(index, 0),
            clockwise=False,
            start=Point3[WorldXYZ].build(float(index), 0.0, 0.0),
            end=Point3[WorldXYZ].build(float(index + 1), 0.0, 0.0),
            start_tangent=direction,
            end_tangent=direction,
        )
        for index in range(count)
    )


def _survey(motions: tuple[MotionQuality, ...], rapids: tuple[RapidMotion, ...] = ()) -> PathSurvey:
    return PathSurvey(
        spec=SPEC,
        motions=motions,
        rapids=rapids,
        plunges=0,
        retracts=0,
        final_stock=cast(Stock, object()),
        total_length=sum(motion.length for motion in motions) + sum(rapid.length for rapid in rapids),
        cut_length=sum(motion.length for motion in motions),
        air_length=sum(rapid.length for rapid in rapids),
        plunge_swept_area=0.0,
    )


def _coverage(uncut: int = 0) -> CoverageEstimate:
    return CoverageEstimate(nx=10, ny=10, cell_area=1.0, reachable_samples=10, uncut_reachable_samples=uncut, remaining_samples=uncut, wall_scallop_height=0.0)


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
    max_engagement_step: MeasuredStep[Degrees] | None = None,
    max_loop_radius_step: MeasuredStep[ToolRadiusMultiple] | None = None,
    tangent_break_pairs: tuple[OperationPair, ...] = (),
    curvature_break_pairs: tuple[OperationPair, ...] = (),
    reversal_pairs: tuple[OperationPair, ...] = (),
) -> PathQualityAttribution:
    return PathQualityAttribution.build(
        operation_count=operation_count,
        uncut_operations=(),
        gouging_operations=(),
        unsafe_rapid_operations=(),
        continuity_break_pairs=(),
        zero_length_operations=(),
        degenerate_loop_operations=(),
        redundant_operations=(),
        cap_exceeded_operations=(),
        slotting_operations=(),
        max_engagement_step=max_engagement_step,
        engagement_step_failure_pairs=(),
        max_loop_radius_step=max_loop_radius_step,
        loop_radius_step_failure_pairs=(),
        tangent_break_pairs=tangent_break_pairs,
        curvature_break_pairs=curvature_break_pairs,
        reversal_pairs=reversal_pairs,
    )


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
        RapidMotion(index=1, operation=OperationType.LINK, length=0.0, horizontal_at_cut_plane=True),
        RapidMotion(index=3, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),
    )
    assessment = assess_path_quality(SPEC, _snapshot(5), _survey(motions, rapids), _coverage(uncut=2))

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


def test_junction_findings_preserve_adjacency_and_classification() -> None:
    tolerance = CONTINUITY_TOOL_RADIUS_FRACTION * SPEC.tool_radius
    motions = (
        _motion(0, end=(0.0, 0.0), end_tangent=(1.0, 0.0)),
        _motion(1, start=(tolerance, 0.0), end_tangent=(-1.0, 0.0), samples=(_sample(10.0),)),
        _motion(2, start=(2.0 + tolerance, 0.0), start_tangent=(1.0, 0.0), end_tangent=(1.0, 0.0), samples=(_sample(40.0),)),
        _motion(4, start=(99.0, 0.0), start_tangent=(0.0, 1.0), samples=(_sample(80.0),)),
    )
    assessment = assess_path_quality(SPEC, _snapshot(5), _survey(motions), _coverage())

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
    assessment = assess_path_quality(SPEC, _snapshot(5), _survey(motions), _coverage())
    pair_01 = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=5)
    pair_12 = OperationPair.build(previous=OperationIndex(1), current=OperationIndex(2), operation_count=5)

    assert assessment.max_engagement_step.measured == Degrees(130.0)
    assert assessment.attribution.max_engagement_step is not None
    assert assessment.attribution.max_engagement_step.pair == pair_01
    assert assessment.attribution.engagement_step_failure_pairs == (pair_01, pair_12)


def test_zero_steps_have_no_invented_maximum_pair() -> None:
    motions = (_motion(0, samples=(_sample(20.0),)), _motion(1, samples=(_sample(20.0),)))
    assessment = assess_path_quality(SPEC, _snapshot(2), _survey(motions), _coverage())

    assert assessment.max_engagement_step.measured == Degrees(0.0)
    assert assessment.attribution.max_engagement_step is None


def test_non_reversal_tangent_break_excludes_curvature_break() -> None:
    motions = (
        _motion(0, radius=2.0, end_tangent=(1.0, 0.0)),
        _motion(1, radius=3.0, start_tangent=(0.0, 1.0)),
    )
    assessment = assess_path_quality(SPEC, _snapshot(2), _survey(motions), _coverage())
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
    rapids = (RapidMotion(index=3, operation=OperationType.RETRACT, length=1.0, horizontal_at_cut_plane=False),)
    assessment = assess_path_quality(SPEC, _snapshot(8, {0: 0, 1: 0, 2: 0, 4: 0, 5: 1, 6: 1, 7: 1}), _survey(motions, rapids), _coverage())
    pair_02 = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(2), operation_count=8)
    pair_56 = OperationPair.build(previous=OperationIndex(5), current=OperationIndex(6), operation_count=8)
    pair_67 = OperationPair.build(previous=OperationIndex(6), current=OperationIndex(7), operation_count=8)

    assert assessment.max_loop_radius_step.measured == ToolRadiusMultiple(3.0)
    assert assessment.attribution.max_loop_radius_step is not None
    assert assessment.attribution.max_loop_radius_step.pair == pair_02
    assert assessment.attribution.loop_radius_step_failure_pairs == (pair_02, pair_56, pair_67)
