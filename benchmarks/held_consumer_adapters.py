"""Value-only adapters from Figure 5 evidence ports to existing consumers."""

from benchmarks.coverage import COVERAGE_GRID_SAMPLES
from benchmarks.coverage import minimum_coverage_grid
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import reduce_quality_evidence
from benchmarks.spec import PocketSpec
from benchmarks.survey import PathSurvey
from compas_cgal.engagement import EngagementReport
from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.toolpath import ToolpathResult


def audit_figure5_engagement(
    spec: PocketSpec,
    result: ToolpathResult,
) -> EngagementReport:
    """Run the existing guarded audit with the fixed Figure 5 spec values."""
    return audit_toolpath_engagement(
        spec.polygon,
        result,
        spec.tool_diameter,
        spec.tea_cap_rad,
        list(spec.holes),
    )


def reduce_figure5_quality(
    spec: PocketSpec,
    snapshot: tuple[HeldOperationSnapshot, ...],
    survey: PathSurvey,
) -> QualityEvidence:
    """Reduce Figure 5 at the smallest protected coverage grid."""
    grid = max(COVERAGE_GRID_SAMPLES, minimum_coverage_grid(spec))
    return reduce_quality_evidence(spec, snapshot, survey, grid=grid)
