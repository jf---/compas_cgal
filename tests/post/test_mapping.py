"""Operation-to-motion mapping on small hand-made operation streams."""

from enum import Enum
from fractions import Fraction

import numpy as np
import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

from compas_cgal.pipeline.emit import arc_operation
from compas_cgal.post.errors import DiscontinuousPathError
from compas_cgal.post.errors import MissingCircleEntryError
from compas_cgal.post.errors import MissingEntryPlungeError
from compas_cgal.post.errors import MissingFinalRetractError
from compas_cgal.post.errors import MissingStartTangentError
from compas_cgal.post.errors import NonVerticalMoveError
from compas_cgal.post.errors import ToolDiameterMismatchError
from compas_cgal.post.errors import UnmappedHeightError
from compas_cgal.post.errors import UnsupportedOperationError
from compas_cgal.post.errors import WindingDisagreementError
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.mapping import SourcePath
from compas_cgal.post.mapping import map_operations
from compas_cgal.post.program import CentreOffset
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import ProgramPoint
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.post.support import mm

SOURCE = SourcePath(cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0), tool=ToolRadius.build(1.0))


class FutureOperationType(str, Enum):
    """An operation the post has no mapping for yet (Task 1.2's helix)."""

    HELIX = "helix"


def _plunge(x: float = 0.0, y: float = 0.0, top: float = 4.0) -> ToolpathOperation:
    return ToolpathOperation(geometry=Line([x, y, top], [x, y, 0.0]), operation=OperationType.PLUNGE, path_index=0)


def _retract(x: float = 0.0, y: float = 0.0) -> ToolpathOperation:
    return ToolpathOperation(geometry=Line([x, y, 0.0], [x, y, 4.0]), operation=OperationType.RETRACT, path_index=0)


def _loop(clockwise: bool, entry: tuple[float, float] | None = (0.0, 0.0), centre_x: float = 1.0, travels_clockwise: bool | None = None) -> ToolpathOperation:
    """A unit circle entered on its west point; the tangent follows *travels_clockwise* (default: the flag)."""
    frame = Frame([centre_x, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    sense = clockwise if travels_clockwise is None else travels_clockwise
    tangent = np.array([0.0, 1.0 if sense else -1.0, 0.0])
    return ToolpathOperation(
        geometry=Circle(1.0, frame=frame), operation=OperationType.CUT, path_index=0, clockwise=clockwise, start_tangent=tangent, end_tangent=tangent, entry_xy=entry
    )


def test_loop_becomes_two_half_arcs(setup: PostSetup) -> None:
    posted = map_operations([_plunge(), _loop(clockwise=True), _retract()], setup, SOURCE)
    motions = posted.program.motions
    assert posted.program.approach_z == 12
    assert [m.kind for m in motions] == [MotionKind.RAPID, MotionKind.LINEAR, MotionKind.CLOCKWISE_ARC, MotionKind.CLOCKWISE_ARC, MotionKind.RAPID]
    assert motions[0].target == ProgramPoint(mm(18), mm(12), mm(12))
    assert motions[1].target == ProgramPoint(mm(18), mm(12), mm(Fraction(-3, 2))) and motions[1].feed == 100
    assert motions[2].target == ProgramPoint(mm(24), mm(12), mm(Fraction(-3, 2))) and motions[2].centre == CentreOffset(mm(3), mm(0))
    assert motions[3].target == ProgramPoint(mm(18), mm(12), mm(Fraction(-3, 2))) and motions[3].centre == CentreOffset(mm(-3), mm(0))
    assert [(s.operation_index, s.first_motion, s.motion_count) for s in posted.spans] == [(0, 1, 1), (1, 2, 2), (2, 4, 1)]


def test_counter_clockwise_loop_is_g3(setup: PostSetup) -> None:
    posted = map_operations([_plunge(), _loop(clockwise=False), _retract()], setup, SOURCE)
    assert [m.kind for m in posted.program.motions[2:4]] == [MotionKind.COUNTERCLOCKWISE_ARC] * 2


def test_clockwise_arc_travels_the_short_way(setup: PostSetup) -> None:
    # emit.py stores a CW quarter turn as its CCW complement (3pi/2); direction must come from `clockwise`.
    origin = Point2[WorldXY].build(0.0, 0.0)
    arc = arc_operation(Point2[WorldXY].build(1.0, 0.0), origin, Point2[WorldXY].build(1.0, 1.0), True, CutZ.build(0.0), 0, OperationType.LINK)
    posted = map_operations([_plunge(), arc, _retract(1.0, 1.0)], setup, SOURCE)
    motion = posted.program.motions[2]
    assert motion.kind is MotionKind.CLOCKWISE_ARC
    assert motion.target == ProgramPoint(mm(21), mm(15), mm(Fraction(-3, 2)))
    assert motion.centre == CentreOffset(mm(3), mm(0))


def test_unknown_operation_raises(setup: PostSetup) -> None:
    helix = ToolpathOperation(geometry=Line([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]), operation=FutureOperationType.HELIX, path_index=0)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedOperationError, match="helix"):
        map_operations([_plunge(), helix, _retract()], setup, SOURCE)


def test_helix_on_a_non_helix_geometry_raises(setup: PostSetup) -> None:
    helix = ToolpathOperation(geometry=Line([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]), operation=OperationType.HELIX, path_index=0)
    with pytest.raises(UnsupportedOperationError, match="HELIX"):
        map_operations([_plunge(), helix, _retract()], setup, SOURCE)


def test_circle_plunge_raises(setup: PostSetup) -> None:
    circular = _loop(clockwise=True)
    circular.operation = OperationType.PLUNGE
    with pytest.raises(UnsupportedOperationError):
        map_operations([_plunge(), circular, _retract()], setup, SOURCE)


def test_circle_without_exact_entry_raises(setup: PostSetup) -> None:
    with pytest.raises(MissingCircleEntryError):
        map_operations([_plunge(), _loop(clockwise=True, entry=None), _retract()], setup, SOURCE)


def test_height_off_both_planes_raises(setup: PostSetup) -> None:
    with pytest.raises(UnmappedHeightError):
        map_operations([_plunge(top=3.0), _loop(clockwise=True), _retract()], setup, SOURCE)


def test_gap_between_operations_raises(setup: PostSetup) -> None:
    with pytest.raises(DiscontinuousPathError):
        map_operations([_plunge(), _loop(clockwise=True, entry=(0.5, 0.0), centre_x=1.5), _retract(0.5, 0.0)], setup, SOURCE)


def test_sloped_plunge_raises(setup: PostSetup) -> None:
    sloped = ToolpathOperation(geometry=Line([0.0, 0.0, 4.0], [0.5, 0.0, 0.0]), operation=OperationType.PLUNGE, path_index=0)
    with pytest.raises(NonVerticalMoveError):
        map_operations([sloped, _retract(0.5, 0.0)], setup, SOURCE)


def test_program_must_start_with_a_plunge(setup: PostSetup) -> None:
    with pytest.raises(MissingEntryPlungeError):
        map_operations([_loop(clockwise=True), _retract()], setup, SOURCE)


def test_program_must_end_with_a_retract(setup: PostSetup) -> None:
    with pytest.raises(MissingFinalRetractError):
        map_operations([_plunge(), _loop(clockwise=True)], setup, SOURCE)


def test_cutter_must_match_the_scaled_source_tool(setup: PostSetup) -> None:
    other = SourcePath(cut=SOURCE.cut, clearance=SOURCE.clearance, tool=ToolRadius.build(2.0))
    with pytest.raises(ToolDiameterMismatchError):
        map_operations([_plunge(), _loop(clockwise=True), _retract()], setup, other)


@pytest.mark.parametrize("clockwise", [True, False])
def test_flag_disagreeing_with_the_tangent_raises(setup: PostSetup, clockwise: bool) -> None:
    flipped = _loop(clockwise=clockwise, travels_clockwise=not clockwise)
    with pytest.raises(WindingDisagreementError):
        map_operations([_plunge(), flipped, _retract()], setup, SOURCE)


def test_defaulted_flag_on_a_clockwise_arc_raises(setup: PostSetup) -> None:
    # An emitter that forgets `clockwise` gets the dataclass default False; the tangent still says CW.
    origin = Point2[WorldXY].build(0.0, 0.0)
    arc = arc_operation(Point2[WorldXY].build(1.0, 0.0), origin, Point2[WorldXY].build(1.0, 1.0), True, CutZ.build(0.0), 0, OperationType.LINK)
    arc.clockwise = False
    with pytest.raises(WindingDisagreementError):
        map_operations([_plunge(), arc, _retract(1.0, 1.0)], setup, SOURCE)


def test_circular_operation_without_tangent_raises(setup: PostSetup) -> None:
    loop = _loop(clockwise=True)
    loop.start_tangent = None
    with pytest.raises(MissingStartTangentError):
        map_operations([_plunge(), loop, _retract()], setup, SOURCE)
