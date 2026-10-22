"""The engagement generator's visualisation tessellation names the geometry it cannot sample.

Its stream holds only machining circles and lines. Any other operation geometry
(a helix, an arc) has no ``start``/``end`` pair that describes it, so it must
be refused loudly rather than drawn as a chord or crash on a missing attribute.
"""

from __future__ import annotations

import pytest
from compas.geometry import Line

from compas_cgal.engagement_toolpath import UntessellatedGeometryError
from compas_cgal.engagement_toolpath import _tessellate
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

# Angular sampling density passed to the tessellation, samples per radian.
SAMPLES_PER_RADIAN = 4.0


def _helix_operation() -> ToolpathOperation:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(1.0), Millimetre(2.0)), GuideRadius(Millimetre(1.0)), Radian(0.0))
    plane = CutPlane.build(CutZ.build(0.0), ClearanceZ.build(4.0))
    helix = Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(2.0)), clockwise=False)
    return ToolpathOperation(geometry=helix, operation=OperationType.HELIX, path_index=0)


def test_a_line_tessellates_to_its_endpoints() -> None:
    line = ToolpathOperation(geometry=Line([0.0, 0.0, 1.0], [2.0, 0.0, 1.0]), operation=OperationType.CUT, path_index=0)
    assert _tessellate([line], SAMPLES_PER_RADIAN).tolist() == [[0.0, 0.0, 1.0], [2.0, 0.0, 1.0]]


def test_a_helix_is_refused_by_name() -> None:
    assert isinstance(_helix_operation().geometry, Helix)
    with pytest.raises(UntessellatedGeometryError):
        _tessellate([_helix_operation()], SAMPLES_PER_RADIAN)
