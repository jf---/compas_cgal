"""Measure the off-live analytic Figure 5 path against Held's published path (orchestrator check, not committed)."""
import math
from benchmarks.held_figure5_analytic import AnalyticFigure5Input, build_analytic_figure5_path
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.adaptive.motion import EngagementCap

case = load_held_reference_case("figure5")
r = float(case.tool_radius.value)
path = build_analytic_figure5_path(AnalyticFigure5Input.build(boundary=case.boundary, tool_radius=case.tool_radius, cap=EngagementCap.build(math.radians(float(case.tea_cap)))))
radii = [float(c.guide_radius.value) for c in path.circles]
pub = load_figure5_publisher_path()
print(f"cap={float(case.tea_cap)} deg, tool r={r}")
print(f"analytic: circles={len(radii)} rho<=r={sum(x <= r for x in radii)} min_rho={min(radii):.4f} max_rho={max(radii):.4f} "
      f"circle_length={sum(2 * math.pi * x for x in radii):.1f} mm  termination={path.termination_reason} remaining_progress={float(path.remaining_progress):.2f}")
print(f"publisher: turns={len(pub.turns)} min_radius={min(float(t.radius) for t in pub.turns):.4f}")
