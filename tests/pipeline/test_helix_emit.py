"""Emitting a helix (Task 1.2): the HELIX operation and its visualisation polyline.

`emit.helix_operation` builds the operation directly (amendment M5: it never
passes through `PathState.emit_motion`, which relabels by removal), and
`emit.tessellate` samples it for reporting along its exact descent.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Line

from compas_cgal.helix import HELIX_CLOSING_TURNS
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline import emit
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

BOTTOM_Z = 0.0
TOP_Z = 4.0
CENTRE = (1.0, 2.0)
FOOTPRINT_RADIUS = 1.0
RAMP_DEG = 2.0
# Sampled points are doubles of cos/sin rotations: a few ulps off the exact
# footprint circle at radius ~1 mm. 1e-9 mm is far above that and far below any
# visible error.
RADIUS_SAMPLE_TOLERANCE_MM = 1e-9


def _helix(clockwise: bool) -> Helix:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(CENTRE[0]), Millimetre(CENTRE[1])), GuideRadius(Millimetre(FOOTPRINT_RADIUS)), Radian(0.3))
    plane = CutPlane.build(CutZ.build(BOTTOM_Z), ClearanceZ.build(TOP_Z))
    return Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(RAMP_DEG)), clockwise=clockwise)


@pytest.mark.parametrize("clockwise", [False, True])
def test_helix_operation_is_a_helix_op_with_the_exact_entry(clockwise: bool) -> None:
    helix = _helix(clockwise)
    op = emit.helix_operation(helix, path_index=3)
    assert op.operation is OperationType.HELIX
    assert op.geometry is helix
    assert op.path_index == 3
    assert op.clockwise is clockwise
    assert op.entry_xy == (helix.entry_xy.x, helix.entry_xy.y)


@pytest.mark.parametrize("clockwise", [False, True])
def test_helix_operation_tangents_are_the_unit_xy_travel_direction_at_the_entry(clockwise: bool) -> None:
    """Repo convention (`_tangent_array`, `line_operation`, toolpath.cpp): unit XY travel, z = 0; the slope lives in `Helix.pitch`."""
    helix = _helix(clockwise)
    op = emit.helix_operation(helix, path_index=0)
    assert op.start_tangent is not None and op.end_tangent is not None
    np.testing.assert_array_equal(op.start_tangent, op.end_tangent)
    assert float(op.start_tangent[2]) == 0.0
    assert math.isclose(float(np.linalg.norm(op.start_tangent)), 1.0)
    radial = np.array([helix.entry_xy.x - CENTRE[0], helix.entry_xy.y - CENTRE[1], 0.0])
    assert abs(float(np.dot(radial, op.start_tangent))) < RADIUS_SAMPLE_TOLERANCE_MM
    turn = float(np.cross(radial, op.start_tangent)[2])
    assert bool(turn < 0.0) is clockwise


@pytest.mark.parametrize("clockwise", [False, True])
def test_tessellate_samples_the_helix_descent_from_top_to_bottom_at_the_entry(clockwise: bool) -> None:
    helix = _helix(clockwise)
    polyline = emit.tessellate([emit.helix_operation(helix, path_index=0)])
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    assert tuple(polyline[0]) == (entry[0], entry[1], TOP_Z)
    assert tuple(polyline[-1][2:]) == (BOTTOM_Z,)
    assert tuple(polyline[-1]) == (helix.exit_xy.x, helix.exit_xy.y, BOTTOM_Z)
    ramp_segments = math.ceil(helix.turns * math.tau * emit.POLYLINE_SAMPLES_PER_RADIAN)
    closing_segments = math.ceil(HELIX_CLOSING_TURNS * math.tau * emit.POLYLINE_SAMPLES_PER_RADIAN)
    # Ramp then the flat closing turn; the ramp's exact end is the closing turn's start (one shared point).
    assert len(polyline) == ramp_segments + 1 + closing_segments
    assert np.all(np.diff(polyline[: ramp_segments + 1, 2]) < 0.0)
    assert np.all(polyline[ramp_segments:, 2] == BOTTOM_Z)
    radii = np.hypot(polyline[:, 0] - CENTRE[0], polyline[:, 1] - CENTRE[1])
    assert np.all(np.abs(radii - float(helix.footprint.radius.value)) < RADIUS_SAMPLE_TOLERANCE_MM)
    second = polyline[1]
    turn = (entry[0] - CENTRE[0]) * (second[1] - CENTRE[1]) - (entry[1] - CENTRE[1]) * (second[0] - CENTRE[0])
    assert bool(turn < 0.0) is clockwise


def test_tessellate_joins_a_helix_to_the_motion_after_it() -> None:
    helix = _helix(False)
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    after = ToolpathOperation(geometry=Line([entry[0], entry[1], BOTTOM_Z], [entry[0] + 1.0, entry[1], BOTTOM_Z]), operation=OperationType.CUT, path_index=0)
    polyline = emit.tessellate([emit.helix_operation(helix, path_index=0), after])
    assert tuple(polyline[-1]) == (entry[0] + 1.0, entry[1], BOTTOM_Z)


def test_tessellate_refuses_an_unknown_geometry_by_name() -> None:
    op = ToolpathOperation(geometry=object(), operation=OperationType.CUT, path_index=0)
    with pytest.raises(emit.UnhandledOperationGeometryError):
        emit.tessellate([op])


def test_arc_geometry_is_the_arc_operation_geometry() -> None:
    centre = Point2[WorldXY].build(0.0, 0.0)
    start = Point2[WorldXY].build(1.0, 0.0)
    end = Point2[WorldXY].build(0.0, 1.0)
    geometry = emit.arc_geometry(centre, start, end, CutZ.build(BOTTOM_Z))
    op = emit.arc_operation(centre, start, end, True, CutZ.build(BOTTOM_Z), 0, OperationType.LINK)
    assert isinstance(op.geometry, Arc)
    assert geometry.__data__ == op.geometry.__data__
