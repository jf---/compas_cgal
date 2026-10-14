"""Does the gate score a straight slot into solid stock at 180 deg (physical) or 360 deg (direction-blind)? Orchestrator check."""
import math
import numpy as np
from compas.geometry import Line, Polygon
from benchmarks.spec import PocketSpec
from benchmarks.survey import survey_path
from compas_cgal.toolpath import OperationType, ToolpathOperation, ToolpathResult

square = Polygon([[0, 0, 0], [40, 0, 0], [40, 40, 0], [0, 40, 0]])
spec = PocketSpec.build(name="probe", family="analytic", polygon=square, tool_diameter=2.0, tea_cap_deg=120.0, params={})
t = np.array([1.0, 0.0, 0.0])
ops = [
    ToolpathOperation(Line([10, 20, 5], [10, 20, 0]), OperationType.PLUNGE, 0, False, np.zeros(3), np.zeros(3)),
    ToolpathOperation(Line([10, 20, 0], [20, 20, 0]), OperationType.CUT, 0, False, t, t),
    ToolpathOperation(Line([20, 20, 0], [20, 20, 5]), OperationType.RETRACT, 0, False, np.zeros(3), np.zeros(3)),
]
survey = survey_path(spec, ToolpathResult(operations=ops, polyline=np.array([[10, 20, 0], [20, 20, 0]], dtype=float)))
for motion in survey.motions:
    degs = [round(s.engagement_deg, 1) for s in motion.samples]
    print("straight slot into solid stock after a plunge, sampled engagement deg:", degs[:3], "...", degs[-3:])
