"""Value-only adapters from Figure 5 evidence ports to existing consumers."""

from benchmarks.spec import PocketSpec
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
