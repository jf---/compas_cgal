"""Re-measure the g6 census claims: publisher Figure 5 vs our figure5_refined (orchestrator check, not committed)."""
import math
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks import held_toolpath_benchmark as hb
from compas_cgal.motion import EngagementCap

pub = load_figure5_publisher_path()
r_pub = [float(t.radius) for t in pub.turns]
print(f"publisher: turns={len(r_pub)} min_radius={min(r_pub):.4f} turn_length={sum(abs(float(t.signed_sweep)) * float(t.radius) for t in pub.turns):.1f} mm")

case = load_held_reference_case("figure5")
r = float(case.tool_radius.value)
guide = build_figure5_raw_guide(case)
components = figure7_inward_offset(case).components
initial = hb.build_figure5_approximate_path(guide, components)
cap = EngagementCap.build(math.radians(float(case.tea_cap)))
refined = refine_figure5_engagement(initial.path.circles, components[0], case.tool_radius, cap)
radii = [float(c.radius.value) for c in refined.circles]
print(f"ours: tool r={r} initial_circles={len(initial.path.circles)} refined_circles={len(radii)} "
      f"rho<=r={sum(x <= r for x in radii)} rho<=0.1r={sum(x <= 0.1 * r for x in radii)} min_rho={min(radii):.5f} "
      f"circle_length={sum(2 * math.pi * x for x in radii):.1f} mm")
