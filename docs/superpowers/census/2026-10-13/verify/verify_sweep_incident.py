"""verify_sweep.py with loops AUTHORED through IncidentLoop: does a plunge at a loop's own entry, or the same loop from another start, still change the stock? Orchestrator check."""
import math
import numpy as np
from compas.geometry import Circle, Frame, Polygon
from compas_cgal.engagement import _subtract_operation
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.stock import Stock
from compas_cgal.units import GuideRadius, Point2, Radian, WorldXY
from compas_cgal.toolpath import OperationType, ToolpathOperation

square = Polygon([[0, 0, 0], [40, 0, 0], [40, 40, 0], [0, 40, 0]])
r, rho, c = 1.0, 2.5, (20.0, 20.0)

def loop(authored, angle):
    ux, uy = math.cos(angle), math.sin(angle)
    frame = Frame([authored.centre.x, authored.centre.y, 0.0], [ux, uy, 0.0], [-uy, ux, 0.0])
    return ToolpathOperation(Circle(authored.radius.value, frame=frame), OperationType.CUT, 0, True, np.zeros(3), np.zeros(3))

plunge_changes = resweep_changes = 0
trials = 12
for k in range(trials):
    a = 2 * math.pi * k / trials
    authored = IncidentLoop.build(Point2[WorldXY].build(*c), GuideRadius.build(rho), Radian(a))
    stock = Stock(square)
    _subtract_operation(stock, loop(authored, authored.phase), r)
    probe = stock.clone(); probe.subtract_disk(authored.entry.x, authored.entry.y, r)
    plunge_changes += int(not probe.exactly_equals(stock))
    probe = stock.clone(); _subtract_operation(probe, loop(authored, authored.phase + 1.0), r)
    resweep_changes += int(not probe.exactly_equals(stock))
print(f"plunge at own entry changes stock: {plunge_changes}/{trials}; same loop from another start changes stock: {resweep_changes}/{trials}")
