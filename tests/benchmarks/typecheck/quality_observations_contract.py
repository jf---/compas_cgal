from typing import Optional

from typing_extensions import assert_type

from benchmarks.coverage import CoverageEstimate
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality_observations import CountCriterion
from benchmarks.quality_observations import DegreesCriterion
from benchmarks.quality_observations import FractionCriterion
from benchmarks.quality_observations import MeasuredStep
from benchmarks.quality_observations import PathQualityAssessment
from benchmarks.quality_observations import ToolRadiusMultipleCriterion
from benchmarks.quality_observations import assess_path_quality
from benchmarks.spec import PocketSpec
from benchmarks.survey import PathSurvey
from benchmarks.units import Degrees
from benchmarks.units import ToolRadiusMultiple


def _contract(spec: PocketSpec, snapshot: tuple[HeldOperationSnapshot, ...], survey: PathSurvey, coverage: CoverageEstimate) -> None:
    assessment = assert_type(assess_path_quality(spec, snapshot, survey, coverage), PathQualityAssessment)
    assert_type(assessment.uncut_fraction, FractionCriterion)
    assert_type(assessment.gouging_motions, CountCriterion)
    assert_type(assessment.max_engagement_step, DegreesCriterion)
    assert_type(assessment.max_loop_radius_step, ToolRadiusMultipleCriterion)
    assert_type(assessment.attribution.max_engagement_step, Optional[MeasuredStep[Degrees]])
    assert_type(assessment.attribution.max_loop_radius_step, Optional[MeasuredStep[ToolRadiusMultiple]])
