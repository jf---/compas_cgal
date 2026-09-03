from typing import Literal
from typing import Optional

from typing_extensions import assert_type

from benchmarks.coverage import CoverageEstimate
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality_observations import CountCriterion
from benchmarks.quality_observations import DegreesCriterion
from benchmarks.quality_observations import FractionCriterion
from benchmarks.quality_observations import MeasuredStep
from benchmarks.quality_observations import OperationPair
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import ToolRadiusMultipleCriterion
from benchmarks.quality_observations import assess_path_quality
from benchmarks.spec import PocketSpec
from benchmarks.survey import PathSurvey
from benchmarks.units import Degrees
from benchmarks.units import OperationIndex
from benchmarks.units import ToolRadiusMultiple


def _contract(spec: PocketSpec, snapshot: tuple[HeldOperationSnapshot, ...], survey: PathSurvey, coverage: CoverageEstimate) -> None:
    assessment = assert_type(assess_path_quality(spec, snapshot, survey, coverage), PathQualityAssessment)
    assert_type(assessment.uncut_fraction, FractionCriterion)
    assert_type(assessment.gouging_motions, CountCriterion)
    assert_type(assessment.max_engagement_step, DegreesCriterion)
    assert_type(assessment.max_loop_radius_step, ToolRadiusMultipleCriterion)
    assert_type(assessment.attribution.max_engagement_step, Optional[MeasuredStep[Degrees]])
    assert_type(assessment.attribution.max_loop_radius_step, Optional[MeasuredStep[ToolRadiusMultiple]])
    pair = OperationPair.build(previous=OperationIndex(0), current=OperationIndex(1), operation_count=2)
    degrees = assert_type(MeasuredStep.build(value=Degrees(1.0), pair=pair, unit="degrees"), MeasuredStep[Degrees])
    radii = assert_type(
        MeasuredStep.build(value=ToolRadiusMultiple(1.0), pair=pair, unit="tool_radius_multiple"),
        MeasuredStep[ToolRadiusMultiple],
    )
    assert_type(degrees.unit, Literal["degrees", "tool_radius_multiple"])
    assert_type(radii.unit, Literal["degrees", "tool_radius_multiple"])

    MeasuredStep.build(value=Degrees(1.0), pair=pair, unit="tool_radius_multiple")  # type: ignore[call-overload]
    MeasuredStep.build(value=ToolRadiusMultiple(1.0), pair=pair, unit="degrees")  # type: ignore[call-overload]
