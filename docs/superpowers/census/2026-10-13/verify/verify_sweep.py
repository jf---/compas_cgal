"""After replaying a full loop, does a plunge disk at its own entry, or the same loop again, still change the stock? Orchestrator check."""
import math
import numpy as np
from compas.geometry import Circle, Frame, Polygon
from compas_cgal.engagement import _subtract_operation
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType, ToolpathOperation

square = Polygon([[0, 0, 0], [40, 0, 0], [40, 40, 0], [0, 40, 0]])
r, rho, c = 1.0, 2.5, (20.0, 20.0)

def loop(angle):
    ux, uy = math.cos(angle), math.sin(angle)
    frame = Frame([c[0], c[1], 0.0], [ux, uy, 0.0], [-uy, ux, 0.0])
    return ToolpathOperation(Circle(rho, frame=frame), OperationType.CUT, 0, True, np.zeros(3), np.zeros(3))

plunge_changes = resweep_changes = 0
trials = 12
for k in range(trials):
    a = 2 * math.pi * k / trials
    stock = Stock(square)
    _subtract_operation(stock, loop(a), r)
    entry = (c[0] + rho * math.cos(a), c[1] + rho * math.sin(a))
    probe = stock.clone(); probe.subtract_disk(entry[0], entry[1], r)
    plunge_changes += int(not probe.exactly_equals(stock))
    probe = stock.clone(); _subtract_operation(probe, loop(a + 1.0), r)
    resweep_changes += int(not probe.exactly_equals(stock))
print(f"plunge at own entry changes stock: {plunge_changes}/{trials}; same loop from another start changes stock: {resweep_changes}/{trials}")
