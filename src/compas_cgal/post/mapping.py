"""Map generated `ToolpathOperation`s onto program motions: the post's one boundary from doubles to exact.

Each double is injected exactly (`Fraction(float)`), scaled and shifted in
`Fraction`, then quantised once onto the output grid. Heights are not scaled:
the source has exactly two planes (cut, clearance) and each maps onto the
setup's program plane.

Direction comes from `ToolpathOperation.clockwise`, cross-checked exactly
against the sign of ``cross(centre -> start, start_tangent)`` because the flag
has a silent default (False). The geometry encodes
counter-clockwise: a `Circle` has no sense (frame z is +Z), and `emit.py`
stores a clockwise `Arc` as its counter-clockwise complement, so reading the
sense or sweep from geometry inverts every clockwise motion. Arc ends are read
as points (`point_at(0)`, `point_at(1)`), never as angles.

Mapping (`testcut-emission-audit.md` section 4):

* PLUNGE -> G1 in Z at the plunge feed; RETRACT -> G0 in Z.
* LINK line at clearance -> G0; LINK line on the cut plane -> G1 at the cut feed.
* CUT, LEAD_IN, LEAD_OUT, and LINK arcs -> G1 (line), G2/G3 (arc), or two
  G2/G3 half-arcs from `entry_xy` (full circle; never a 360-degree block).
* HELIX -> a G1 at the plunge feed from clearance straight down to the ramp
  start ``PostSetup.ramp_start_z`` (in air), then LinuxCNC helical
  interpolation at the ramp feed from there to the cut plane. The ramp turn
  count is decided in PROGRAM units (`program_ramp_turns`): the source count
  is for source heights, and the post does not scale heights. Each ramp turn
  is two G2/G3 half-turns (entry -> exact antipode -> entry, like a full
  circle) whose Z is interpolated exactly (`Fraction`) and quantised once per
  block; then `HELIX_CLOSING_TURNS` flat closing turn (two half-arcs) at the
  cut plane, same operation (Q6 ruling). Footprint centre, radius, entry and
  sense are the pipeline helix's; sense from ``Helix.clockwise`` only,
  cross-checked against the operation's flag and start tangent.
* The mapper synthesises the first motion, a G0 to clearance above the entry
  (the entry PLUNGE's top, or the entry HELIX's ``entry_xy``), after the
  declared Z-only approach.

Approach rule for a HELIX entry (orchestrator ruling, G3 round 1): the G0
reaches clearance above the entry (the helix top, ``Helix.top_z``, is the
source clearance and maps onto it), then the G1 descends to the ramp start.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any
from typing import Never
from typing import NoReturn

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from compas_cgal.helix import HELIX_CLOSING_TURNS
from compas_cgal.helix import Helix
from compas_cgal.helix import ramp_turns
from compas_cgal.post.errors import DiscontinuousPathError
from compas_cgal.post.errors import HelixOperationMismatchError
from compas_cgal.post.errors import InexactProgramHelixError
from compas_cgal.post.errors import MissingCircleEntryError
from compas_cgal.post.errors import MissingEntryPlungeError
from compas_cgal.post.errors import MissingFinalRetractError
from compas_cgal.post.errors import MissingStartTangentError
from compas_cgal.post.errors import NonCentreCuttingHelixError
from compas_cgal.post.errors import NonVerticalMoveError
from compas_cgal.post.errors import ToolDiameterMismatchError
from compas_cgal.post.errors import UnmappedHeightError
from compas_cgal.post.errors import UnsupportedOperationError
from compas_cgal.post.errors import WindingDisagreementError
from compas_cgal.post.grid import quantise
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.program import CentreOffset
from compas_cgal.post.program import Motion
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import OperationSpan
from compas_cgal.post.program import PostedProgram
from compas_cgal.post.program import Program
from compas_cgal.post.program import ProgramPoint
from compas_cgal.post.quantities import MillimetrePerMinute
from compas_cgal.post.quantities import ProgramMillimetre
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import ToolRadius

# Source points are XYZ triples; XY points are pairs.
XY = tuple[Fraction, Fraction]


@dataclass(frozen=True)
class SourcePath:
    """The generated path's own planes and tool, in source units.

    Attributes:
        cut: The source cut plane.
        clearance: The source clearance plane.
        tool: The source tool radius.
    """

    cut: CutZ
    clearance: ClearanceZ
    tool: ToolRadius


def _unknown_operation(operation: Never) -> NoReturn:
    raise UnsupportedOperationError(f"no mapping for operation {operation!r}.")


def _exact_xy(point: Any) -> XY:
    return Fraction(float(point[0])), Fraction(float(point[1]))


def program_ramp_turns(helix: Helix, setup: PostSetup, top: ProgramMillimetre, bottom: ProgramMillimetre) -> int:
    """The ramp turn count of *helix* in program units, by the core rule `compas_cgal.helix.ramp_turns`.

    The count is decided on the exact program radius ``scale * rho_h`` (before
    quantisation); `compas_cgal.post.checks` certifies the ramp on the
    emitted, quantised radius. Near a tie (the exact pitch within a quantum
    of ``K * r``) the two can disagree: the post then refuses the program
    (`HelixRampAngleError`) instead of adding a turn. That is fail-closed.

    Args:
        helix: The pipeline helix (source units).
        setup: Supplies the scale and the ramp cap ``K``.
        top: Ramp start, mm.
        bottom: Cut plane, mm.

    Returns:
        The program-unit ramp turn count.

    Raises:
        InexactProgramHelixError: The program radius is not exactly a double (`GuideRadius` carries a double).
    """
    radius = setup.scale * Fraction(helix.footprint.radius.value)
    if Fraction(float(radius)) != radius:
        raise InexactProgramHelixError(f"program helix radius {radius} mm is not exactly a double.")
    return ramp_turns(ExactMillimetre(top - bottom), GuideRadius(Millimetre(float(radius))), setup.ramp_cap)


def _check_winding(index: int, op: ToolpathOperation, centre: XY, start: XY) -> None:
    """Refuse a circular op whose flag disagrees with its own start tangent (exact, on the op's doubles).

    Raises:
        MissingStartTangentError: No start tangent.
        WindingDisagreementError: The tangent is radial, or turns the other way than `clockwise` says.
    """
    tangent = op.start_tangent
    if tangent is None:
        raise MissingStartTangentError(f"operation {index}: a circular motion without a start tangent.")
    tangent_x, tangent_y = Fraction(float(tangent[0])), Fraction(float(tangent[1]))
    turn = (start[0] - centre[0]) * tangent_y - (start[1] - centre[1]) * tangent_x
    if turn == 0 or (turn < 0) != op.clockwise:
        sense = "radial" if turn == 0 else ("clockwise" if turn < 0 else "counter-clockwise")
        raise WindingDisagreementError(f"operation {index}: clockwise={op.clockwise} but the start tangent is {sense}.")


def _entry_xy(first: ToolpathOperation) -> XY:
    """Where the entry starts in XY: the PLUNGE line's top, or the HELIX's exact entry.

    Raises:
        UnsupportedOperationError: The entry PLUNGE is not a line, or the entry HELIX carries no `Helix`.
    """
    geometry = first.geometry
    if first.operation is OperationType.HELIX:
        if not isinstance(geometry, Helix):
            raise UnsupportedOperationError(f"the entry HELIX is a {type(geometry).__name__}, not a helix.")
        return Fraction(geometry.entry_xy.x), Fraction(geometry.entry_xy.y)
    if not isinstance(geometry, Line):
        raise UnsupportedOperationError(f"the entry PLUNGE is a {type(geometry).__name__}, not a line.")
    return _exact_xy(geometry.start)


class _Mapper:
    """Accumulates motions while tracking the emitted (quantised) tool position."""

    def __init__(self, setup: PostSetup, source: SourcePath) -> None:
        self.setup = setup
        self.source = source
        self.motions: list[Motion] = []

    @property
    def position(self) -> ProgramPoint:
        return self.motions[-1].target

    def xy(self, xy: XY) -> tuple[ProgramMillimetre, ProgramMillimetre]:
        shift_x, shift_y = self.setup.xy_shift
        digits = self.setup.decimals
        return ProgramMillimetre(quantise(self.setup.scale * xy[0] + shift_x, digits)), ProgramMillimetre(quantise(self.setup.scale * xy[1] + shift_y, digits))

    def height(self, z: float, index: int) -> ProgramMillimetre:
        if z == self.source.cut.value:
            return self.setup.cut_z
        if z == self.source.clearance.value:
            return self.setup.clearance_z
        raise UnmappedHeightError(f"operation {index}: Z {z!r} is neither the cut plane {self.source.cut.value!r} nor clearance {self.source.clearance.value!r}.")

    def point(self, point: Any, index: int) -> ProgramPoint:
        x, y = self.xy(_exact_xy(point))
        return ProgramPoint(x, y, self.height(float(point[2]), index))

    def expect_at(self, point: ProgramPoint, index: int) -> None:
        if point != self.position:
            raise DiscontinuousPathError(f"operation {index} starts at {point}, the program left the tool at {self.position}.")

    def arc_to(self, kind: MotionKind, end: XY, centre: XY, z: ProgramMillimetre, feed: MillimetrePerMinute) -> None:
        here = self.position
        centre_x, centre_y = self.xy(centre)
        end_x, end_y = self.xy(end)
        offset = CentreOffset(ProgramMillimetre(centre_x - here.x), ProgramMillimetre(centre_y - here.y))
        self.motions.append(Motion(kind, ProgramPoint(end_x, end_y, z), offset, feed))

    def vertical(self, index: int, op: ToolpathOperation, kind: MotionKind, feed: MillimetrePerMinute | None) -> None:
        line = op.geometry
        if not isinstance(line, Line):
            raise UnsupportedOperationError(f"operation {index}: {op.operation!r} on a {type(line).__name__}.")
        if _exact_xy(line.start) != _exact_xy(line.end):
            raise NonVerticalMoveError(f"operation {index}: {op.operation!r} moves in XY.")
        self.expect_at(self.point(line.start, index), index)
        self.motions.append(Motion(kind, self.point(line.end, index), None, feed))

    def link(self, index: int, op: ToolpathOperation) -> None:
        line = op.geometry
        if not isinstance(line, Line):
            self.feed(index, op)
            return
        start, end = self.point(line.start, index), self.point(line.end, index)
        if start.z != end.z:
            raise UnsupportedOperationError(f"operation {index}: a LINK between planes.")
        self.expect_at(start, index)
        if end.z == self.setup.clearance_z:
            self.motions.append(Motion(MotionKind.RAPID, end, None, None))
        else:
            self.motions.append(Motion(MotionKind.LINEAR, end, None, self.setup.cutting.cut_feed))

    def feed(self, index: int, op: ToolpathOperation) -> None:
        geometry = op.geometry
        kind = MotionKind.CLOCKWISE_ARC if op.clockwise else MotionKind.COUNTERCLOCKWISE_ARC
        if isinstance(geometry, Line):
            self.expect_at(self.point(geometry.start, index), index)
            self.motions.append(Motion(MotionKind.LINEAR, self.point(geometry.end, index), None, self.setup.cutting.cut_feed))
        elif isinstance(geometry, Circle):
            if op.entry_xy is None:
                raise MissingCircleEntryError(f"operation {index}: a full circle without its exact entry.")
            z = self.height(float(geometry.frame.point[2]), index)
            centre, entry = _exact_xy(geometry.frame.point), (Fraction(op.entry_xy[0]), Fraction(op.entry_xy[1]))
            _check_winding(index, op, centre, entry)
            x, y = self.xy(entry)
            self.expect_at(ProgramPoint(x, y, z), index)
            feed = self.setup.cutting.cut_feed
            self.arc_to(kind, (2 * centre[0] - entry[0], 2 * centre[1] - entry[1]), centre, z, feed)
            self.arc_to(kind, entry, centre, z, feed)
        elif isinstance(geometry, Arc):
            z = self.height(float(geometry.frame.point[2]), index)
            start, centre = _exact_xy(geometry.point_at(0.0)), _exact_xy(geometry.frame.point)
            _check_winding(index, op, centre, start)
            x, y = self.xy(start)
            self.expect_at(ProgramPoint(x, y, z), index)
            self.arc_to(kind, _exact_xy(geometry.point_at(1.0)), centre, z, self.setup.cutting.cut_feed)
        else:
            raise UnsupportedOperationError(f"operation {index}: geometry {type(geometry).__name__}.")

    def helix(self, index: int, op: ToolpathOperation) -> None:
        """G1 in air to the ramp start, program-unit ramp turns as helical half-turns, then the flat closing turn, at the ramp feed.

        Raises:
            UnsupportedOperationError: The geometry is not a `Helix`.
            MissingCircleEntryError: No `entry_xy`.
            MissingStartTangentError: No start tangent.
            WindingDisagreementError: The flag disagrees with the start tangent.
            HelixOperationMismatchError: The operation's entry or sense is not its helix's.
            InexactProgramHelixError: The program-unit turn count cannot be decided exactly.
        """
        helix = op.geometry
        if not isinstance(helix, Helix):
            raise UnsupportedOperationError(f"operation {index}: HELIX on a {type(helix).__name__}.")
        if op.entry_xy is None:
            raise MissingCircleEntryError(f"operation {index}: a helix without its exact entry.")
        entry = (Fraction(op.entry_xy[0]), Fraction(op.entry_xy[1]))
        centre = (Fraction(helix.footprint.centre.x), Fraction(helix.footprint.centre.y))
        _check_winding(index, op, centre, entry)
        if entry != (Fraction(helix.entry_xy.x), Fraction(helix.entry_xy.y)) or op.clockwise != helix.clockwise:
            raise HelixOperationMismatchError(f"operation {index}: the operation's entry/sense is not its helix's.")
        x, y = self.xy(entry)
        self.expect_at(ProgramPoint(x, y, self.height(helix.top_z.value, index)), index)
        top, bottom = self.setup.ramp_start_z, self.height(helix.bottom_z.value, index)
        self.motions.append(Motion(MotionKind.LINEAR, ProgramPoint(x, y, top), None, self.setup.cutting.plunge_feed))
        kind = MotionKind.CLOCKWISE_ARC if helix.clockwise else MotionKind.COUNTERCLOCKWISE_ARC
        antipode = (2 * centre[0] - entry[0], 2 * centre[1] - entry[1])
        feed = self.setup.ramp_feed
        halves = 2 * program_ramp_turns(helix, self.setup, top, bottom)
        for half in range(1, halves + 1):
            z = ProgramMillimetre(quantise(top + (bottom - top) * Fraction(half, halves), self.setup.decimals))
            self.arc_to(kind, antipode if half % 2 else entry, centre, z, feed)
        for _ in range(HELIX_CLOSING_TURNS):
            self.arc_to(kind, antipode, centre, bottom, feed)
            self.arc_to(kind, entry, centre, bottom, feed)

    def operation(self, index: int, op: ToolpathOperation) -> None:
        # Exhaustive by identity: mypy narrows the final branch to Never, so a new
        # `OperationType` member fails type checking here and raises at run time.
        operation = op.operation
        if operation is OperationType.PLUNGE:
            self.vertical(index, op, MotionKind.LINEAR, self.setup.cutting.plunge_feed)
        elif operation is OperationType.RETRACT:
            self.vertical(index, op, MotionKind.RAPID, None)
        elif operation is OperationType.LINK:
            self.link(index, op)
        elif operation is OperationType.CUT or operation is OperationType.LEAD_IN or operation is OperationType.LEAD_OUT:
            self.feed(index, op)
        elif operation is OperationType.HELIX:
            self.helix(index, op)
        else:
            _unknown_operation(operation)


def map_operations(operations: Sequence[ToolpathOperation], setup: PostSetup, source: SourcePath) -> PostedProgram:
    """Map a generated path onto program motions, with a span per operation.

    Args:
        operations: The path, starting with its PLUNGE or HELIX entry and ending with its RETRACT.
        setup: The setup.
        source: The path's planes and tool.

    Returns:
        The program and its trace back to *operations*.

    Raises:
        ToolDiameterMismatchError: ``setup.cutter.diameter != 2 * scale * source tool radius``.
        NonCentreCuttingHelixError: A HELIX with a cutter that is not centre-cutting.
        MissingEntryPlungeError: The path starts with neither a PLUNGE nor a HELIX.
        MissingFinalRetractError: The path does not end with a RETRACT.
        UnsupportedOperationError: An operation type, or operation/geometry pairing, with no mapping.
        UnmappedHeightError: A height off both source planes.
        NonVerticalMoveError: A PLUNGE or RETRACT moving in XY.
        MissingCircleEntryError: A circle without `entry_xy`.
        DiscontinuousPathError: An operation not starting at the emitted tool position.
        MissingStartTangentError: A circular operation without a start tangent.
        WindingDisagreementError: A circular operation whose flag disagrees with its tangent.
        HelixOperationMismatchError: A HELIX operation whose entry or sense is not its helix's.
        InexactProgramHelixError: A HELIX whose program-unit turn count cannot be decided exactly.
    """
    scaled_tool = 2 * setup.scale * Fraction(source.tool.value)
    if setup.cutter.diameter != scaled_tool:
        raise ToolDiameterMismatchError(f"cutter {setup.cutter.diameter} mm is not the source tool scaled: {scaled_tool} mm.")
    if not setup.cutter.centre_cutting and any(op.operation is OperationType.HELIX for op in operations):
        raise NonCentreCuttingHelixError("a HELIX entry ramps on the end face: the cutter must be centre-cutting (pending the user's answer).")
    if not operations or not (operations[0].operation is OperationType.PLUNGE or operations[0].operation is OperationType.HELIX):
        raise MissingEntryPlungeError("the path must start with its PLUNGE or HELIX entry.")
    if operations[-1].operation is not OperationType.RETRACT:
        raise MissingFinalRetractError("the path must end with its RETRACT.")
    mapper = _Mapper(setup, source)
    entry_x, entry_y = mapper.xy(_entry_xy(operations[0]))
    mapper.motions.append(Motion(MotionKind.RAPID, ProgramPoint(entry_x, entry_y, setup.clearance_z), None, None))
    spans: list[OperationSpan] = []
    for index, op in enumerate(operations):
        first = len(mapper.motions)
        mapper.operation(index, op)
        spans.append(OperationSpan(index, first, len(mapper.motions) - first))
    return PostedProgram(Program(setup.cutting.spindle_rpm, setup.spindle_dwell, setup.clearance_z, tuple(mapper.motions)), tuple(spans))
