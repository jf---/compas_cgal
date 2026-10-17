"""The helix carries one flat closing turn at the cut plane in the same op (Task 1.2, ruling Q6 shape).

The HELIX operation is the ramp turns followed by ONE flat full turn at
``bottom_z`` on the footprint circle. The model carries it (REPORTING length
includes it) and `emit.tessellate` samples ramp then closing turn with both
ends at the exact entry ``P``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from compas_cgal.helix import HELIX_CLOSING_TURNS
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline import emit
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


def _helix(clockwise: bool) -> Helix:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(CENTRE[0]), Millimetre(CENTRE[1])), GuideRadius(Millimetre(FOOTPRINT_RADIUS)), Radian(0.3))
    plane = CutPlane.build(CutZ.build(BOTTOM_Z), ClearanceZ.build(TOP_Z))
    return Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(RAMP_DEG)), clockwise=clockwise)


def test_the_reported_length_is_the_ramp_plus_one_flat_footprint_turn() -> None:
    helix = _helix(False)
    circumference = 2.0 * math.pi * helix.footprint.radius.value
    assert helix.closing_turn_length() == circumference
    assert helix.ramp_length() == helix.turns * math.hypot(circumference, float(helix.pitch))
    assert helix.length() == helix.ramp_length() + helix.closing_turn_length()


@pytest.mark.parametrize("clockwise", [False, True])
def test_tessellate_samples_the_ramp_then_one_flat_closing_turn_with_exact_ends(clockwise: bool) -> None:
    helix = _helix(clockwise)
    polyline = emit.tessellate([emit.helix_operation(helix, path_index=0)])
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    ramp_segments = math.ceil(helix.turns * math.tau * emit.POLYLINE_SAMPLES_PER_RADIAN)
    closing_segments = math.ceil(HELIX_CLOSING_TURNS * math.tau * emit.POLYLINE_SAMPLES_PER_RADIAN)
    # The ramp's exact end is the closing turn's exact start: one shared point.
    assert len(polyline) == ramp_segments + 1 + closing_segments
    assert tuple(polyline[0]) == (entry[0], entry[1], TOP_Z)
    assert tuple(polyline[ramp_segments]) == (helix.exit_xy.x, helix.exit_xy.y, BOTTOM_Z)
    assert tuple(polyline[-1]) == (helix.exit_xy.x, helix.exit_xy.y, BOTTOM_Z)
    assert np.all(np.diff(polyline[: ramp_segments + 1, 2]) < 0.0)
    assert np.all(polyline[ramp_segments:, 2] == BOTTOM_Z)
