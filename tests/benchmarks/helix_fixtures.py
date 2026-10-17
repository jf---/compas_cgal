"""Shared HELIX fixtures for the Task 1.2 consumer tests (set B): helices built by the real emitter.

A plain helper module, not a conftest, so no existing test's fixtures change.
"""

from __future__ import annotations

import numpy as np
from compas.geometry import Polygon

from benchmarks.spec import PocketSpec
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import helix_operation
from compas_cgal.pipeline.emit import loop_operation
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

TOOL_DIAMETER = 2.0
TOOL_RADIUS = TOOL_DIAMETER / 2.0
CAP_DEG = 120.0
TOP_Z = 4.0
CUT_Z = 0.0
RAMP_DEG = 2.0
# Footprint radius = tool radius: the bore is the disk of radius 2r (C10 target).
FOOTPRINT_RADIUS = TOOL_RADIUS
# First pencil member: tangent to the footprint at P, larger by 0.5 r.
MEMBER_RADIUS = 1.5

POCKET = PocketSpec.build(
    name="helix_12x8",
    family="analytic",
    polygon=Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]]),
    tool_diameter=TOOL_DIAMETER,
    tea_cap_deg=CAP_DEG,
)


def make_helix(clockwise: bool = True, *, centre: tuple[float, float] = (0.0, 0.0), radius: float = FOOTPRINT_RADIUS, phase: float = 0.0) -> Helix:
    """A helix over the fixture cut plane at the fixture ramp cap."""
    footprint = IncidentLoop.build(Point2[WorldXY].build(*centre), GuideRadius.build(radius), Radian(phase))
    plane = CutPlane.build(CutZ.build(CUT_Z), ClearanceZ.build(TOP_Z))
    return Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(RAMP_DEG)), clockwise=clockwise)


def emitted(helix: Helix) -> ToolpathOperation:
    """The HELIX operation the pipeline emitter produces."""
    return helix_operation(helix, 0)


def first_member(helix: Helix) -> ToolpathOperation:
    """The first pencil member: centre on the ray from P through c_h, entry P (phase-0 helices only)."""
    entry = helix.entry_xy
    loop = IncidentLoop.build(Point2[WorldXY].build(entry.x - MEMBER_RADIUS, entry.y), GuideRadius.build(MEMBER_RADIUS), Radian(0.0))
    return loop_operation(loop, CutZ.build(CUT_Z), 0, clockwise=helix.clockwise)


def result(operations: list[ToolpathOperation]) -> ToolpathResult:
    """A toolpath result; the polyline is unused by every consumer under test."""
    return ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))
