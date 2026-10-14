"""Orientation and corner-rounding controls for the analytic Held lap failure (orchestrator check)."""
import math
from benchmarks.held_figure5_analytic import AnalyticFigure5Input, build_analytic_figure5_path
from benchmarks.held_reference_geometry import ReferenceArc, ReferenceBoundary, ReferenceLine
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Millimetre, Point2, Radian, ToolRadius, WorldXY

tool = ToolRadius.build(Millimetre(1.0))
P = lambda x, y: Point2[WorldXY].build(float(x), float(y))
cap = EngagementCap.build(math.radians(80.0))

def run(name, prims):
    try:
        path = build_analytic_figure5_path(AnalyticFigure5Input.build(boundary=ReferenceBoundary.build(prims, tool, Millimetre(0.01)), tool_radius=tool, cap=cap))
        radii = [float(c.guide_radius.value) for c in path.circles]
        print(f"{name}: circles={len(radii)} rho<=r={sum(x <= 1.0 for x in radii)} min_rho={min(radii):.4f} circle_len={sum(2 * math.pi * x for x in radii):.1f}")
    except Exception as error:
        print(f"{name}: FAILED {type(error).__name__}: {error}")

fig5 = load_held_reference_case("figure5").boundary
area = sum((float(p.start.x) * float(p.end.y) - float(p.end.x) * float(p.start.y)) for p in fig5.primitives) / 2
print(f"figure5 boundary chord-area sign: {'CCW' if area > 0 else 'CW'}")
cw = [(0, 0), (0, 8), (12, 8), (12, 0)]
run("rect_12x8 CW", [ReferenceLine.build(P(*cw[i]), P(*cw[(i + 1) % 4])) for i in range(4)])
# Rounded rectangle 12x8, corner radius 2 (> r): CCW lines + quarter arcs.
R = 2.0
prims = [
    ReferenceLine.build(P(R, 0), P(12 - R, 0)), ReferenceArc.build(P(12 - R, 0), P(12, R), P(12 - R, R), Radian(math.pi / 2)),
    ReferenceLine.build(P(12, R), P(12, 8 - R)), ReferenceArc.build(P(12, 8 - R), P(12 - R, 8), P(12 - R, 8 - R), Radian(math.pi / 2)),
    ReferenceLine.build(P(12 - R, 8), P(R, 8)), ReferenceArc.build(P(R, 8), P(0, 8 - R), P(R, 8 - R), Radian(math.pi / 2)),
    ReferenceLine.build(P(0, 8 - R), P(0, R)), ReferenceArc.build(P(0, R), P(R, 0), P(R, R), Radian(math.pi / 2)),
]
run("rounded rect 12x8 R2 CCW", prims)
