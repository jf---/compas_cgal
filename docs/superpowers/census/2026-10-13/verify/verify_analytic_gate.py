"""Run the off-live analytic Held lap on gate-shaped pockets (orchestrator check, not committed)."""
import math
import traceback
from benchmarks.held_figure5_analytic import AnalyticFigure5Input, build_analytic_figure5_path
from benchmarks.held_reference_geometry import ReferenceBoundary, ReferenceLine
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Millimetre, Point2, ToolRadius, WorldXY

POCKETS = {
    "rect_12x8": [(0, 0), (12, 0), (12, 8), (0, 8)],
    "L_shape": [(0, 0), (12, 0), (12, 6), (6, 6), (6, 10), (0, 10)],
}
tool = ToolRadius.build(Millimetre(1.0))
for name, ring in POCKETS.items():
    pts = [Point2[WorldXY].build(float(x), float(y)) for x, y in ring]
    lines = [ReferenceLine.build(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    boundary = ReferenceBoundary.build(lines, tool, Millimetre(0.01))
    for cap_deg in (80.0, 120.0):
        try:
            path = build_analytic_figure5_path(AnalyticFigure5Input.build(boundary=boundary, tool_radius=tool, cap=EngagementCap.build(math.radians(cap_deg))))
            radii = [float(c.guide_radius.value) for c in path.circles]
            print(f"{name} cap {cap_deg:.0f}: circles={len(radii)} rho<=r={sum(x <= 1.0 for x in radii)} min_rho={min(radii):.4f} max_rho={max(radii):.4f} "
                  f"circle_len={sum(2 * math.pi * x for x in radii):.1f} remaining={float(path.remaining_progress):.2f} term={path.termination_reason}")
        except Exception as error:  # probe: report the failure mode, do not hide it
            print(f"{name} cap {cap_deg:.0f}: FAILED {type(error).__name__}: {error}")
