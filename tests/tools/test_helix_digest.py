"""`pipeline_timing.ops_digest` covers HELIX content deterministically (Task 1.2).

A helix carries no compas ``__data__``; its digest content is its exact
defining fields, read off the op `emit.helix_operation` emits. Every other geometry digests exactly as before, so the
committed digests of output without a helix are unchanged.
"""

from __future__ import annotations

import hashlib
from dataclasses import replace

import numpy as np
from compas.geometry import Circle
from compas.geometry import Line

from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import helix_operation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY
from tools.pipeline_timing import ops_digest
from tools.pipeline_timing_results import DIGEST_HEX_CHARS

RAMP_DEG = 2.0


def _helix(centre_x: float = 1.0, top_z: float = 4.0, clockwise: bool = False, radius: float = 1.0, phase: float = 0.0) -> Helix:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(centre_x), Millimetre(2.0)), GuideRadius(Millimetre(radius)), Radian(phase))
    plane = CutPlane.build(CutZ.build(0.0), ClearanceZ.build(top_z))
    return Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(RAMP_DEG)), clockwise=clockwise)


def _result(*operations: ToolpathOperation) -> ToolpathResult:
    return ToolpathResult(operations=list(operations), polyline=np.empty((0, 3)))


def _committed_formula(result: ToolpathResult) -> str:
    digest = hashlib.sha256()
    for operation in result.operations:
        geometry = operation.geometry
        assert not isinstance(geometry, Helix), "the committed formula predates the helix"
        digest.update(repr((operation.operation, type(geometry).__name__, geometry.__data__)).encode())
    return digest.hexdigest()[:DIGEST_HEX_CHARS]


def _planar() -> ToolpathResult:
    return _result(
        ToolpathOperation(geometry=Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), operation=OperationType.CUT, path_index=0),
        ToolpathOperation(geometry=Circle(2.0), operation=OperationType.CUT, path_index=0, entry_xy=(2.0, 0.0)),
    )


def test_a_helix_digests_deterministically() -> None:
    assert ops_digest(_result(helix_operation(_helix(), 0))) == ops_digest(_result(helix_operation(_helix(), 0)))


def test_the_helix_content_moves_the_digest() -> None:
    base = ops_digest(_result(helix_operation(_helix(), 0)))
    assert ops_digest(_result(helix_operation(_helix(centre_x=1.5), 0))) != base
    assert ops_digest(_result(helix_operation(_helix(top_z=8.0), 0))) != base
    assert ops_digest(_result(helix_operation(_helix(clockwise=True), 0))) != base
    assert ops_digest(_result(helix_operation(_helix(radius=0.75), 0))) != base
    assert ops_digest(_result(helix_operation(_helix(phase=1.0), 0))) != base
    assert ops_digest(_result(helix_operation(replace(_helix(), turns=_helix().turns + 1), 0))) != base


def test_geometry_without_a_helix_digests_exactly_as_committed() -> None:
    assert ops_digest(_planar()) == _committed_formula(_planar())
