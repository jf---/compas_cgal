"""`OperationType.HELIX` and its replay rule (Task 1.2; amendment A0 requires an explicit rule).

A helix is never AIR: replayed as a rapid, depletion would skip the bore. It
has its own rule and category, is not cap-sampled (like a plunge, its next
motion is), and its bottom Z is independent cut-plane evidence.
"""

from __future__ import annotations

from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

import pytest

import compas_cgal.replay_classification as rc
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.toolpath import UnknownNativeOperationCodeError
from compas_cgal.toolpath import decode_native_operation
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

HEAD_CODES = (
    OperationType.CUT,
    OperationType.LEAD_IN,
    OperationType.LEAD_OUT,
    OperationType.LINK,
    OperationType.RETRACT,
    OperationType.PLUNGE,
)
BOTTOM_Z = 0.0
TOP_Z = 4.0
ABOVE_LINE_Z = 5.0


def _helix() -> Helix:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(1.0), Millimetre(2.0)), GuideRadius(Millimetre(1.0)), Radian(0.0))
    plane = CutPlane.build(CutZ.build(BOTTOM_Z), ClearanceZ.build(TOP_Z))
    return Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(2.0)), clockwise=False)


def _helix_operation(helix: Helix, operation: OperationType = OperationType.HELIX) -> ToolpathOperation:
    entry = helix.entry_xy
    return ToolpathOperation(geometry=helix, operation=operation, path_index=0, clockwise=helix.clockwise, entry_xy=(entry.x, entry.y))


def test_native_codes_decode_to_the_six_native_kinds_only() -> None:
    assert OperationType.HELIX.value == "helix"
    assert tuple(decode_native_operation(code) for code in range(len(HEAD_CODES))) == HEAD_CODES


@pytest.mark.parametrize("code", [len(HEAD_CODES), -1, 99])
def test_a_stray_native_code_raises_a_named_error(code: int) -> None:
    with pytest.raises(UnknownNativeOperationCodeError):
        decode_native_operation(code)


def test_helix_has_its_own_replay_rule_and_it_is_not_air() -> None:
    assert rc.replay_rule(OperationType.HELIX) is rc.ReplayRule.HELIX
    assert rc.ReplayRule.HELIX is not rc.ReplayRule.AIR


def test_a_helix_is_excluded_from_engaged_sampling_like_a_plunge() -> None:
    assert not rc.rule_is_engaged(rc.ReplayRule.HELIX)
    assert OperationType.HELIX not in rc.AUDIT_ENGAGED


def test_a_helix_replays_as_the_helix_category() -> None:
    assert rc.classify_operation_replay(_helix_operation(_helix()), Millimetre(BOTTOM_Z)) == "helix"


def test_a_helix_off_the_cut_plane_raises_a_named_error() -> None:
    with pytest.raises(rc.OffPlaneHelixError):
        rc.classify_operation_replay(_helix_operation(_helix()), Millimetre(BOTTOM_Z - 1.0))


def test_the_helix_bottom_anchors_the_inferred_cut_plane() -> None:
    above = ToolpathOperation(geometry=Line([0.0, 0.0, ABOVE_LINE_Z], [1.0, 0.0, ABOVE_LINE_Z]), operation=OperationType.CUT, path_index=0)
    assert rc.infer_operation_cut_height([_helix_operation(_helix())]) == BOTTOM_Z
    assert rc.infer_operation_cut_height([above, _helix_operation(_helix())]) == BOTTOM_Z


def test_a_helix_geometry_under_another_operation_is_refused() -> None:
    with pytest.raises(rc.HelixGeometryError):
        rc.classify_operation_replay(_helix_operation(_helix(), OperationType.CUT), Millimetre(BOTTOM_Z))
    with pytest.raises(rc.HelixGeometryError):
        rc.infer_operation_cut_height([_helix_operation(_helix(), OperationType.CUT)])


def test_a_helix_operation_without_helix_geometry_is_refused() -> None:
    line = ToolpathOperation(geometry=Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), operation=OperationType.HELIX, path_index=0)
    circle = ToolpathOperation(geometry=Circle(1.0, frame=Frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])), operation=OperationType.HELIX, path_index=0)
    for operation in (line, circle):
        with pytest.raises(rc.HelixGeometryError):
            rc.classify_operation_replay(operation, Millimetre(BOTTOM_Z))
