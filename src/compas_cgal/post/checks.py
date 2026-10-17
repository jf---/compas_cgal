"""Exact checks over a parsed program: each failure mode raises its own named error.

Run on what was written (parse the text back first), never on what the mapper
meant to write. Every comparison is on `Fraction`s read from the emitted
decimals.

There is no stock model here, so "a rapid in material" is not decidable. The
rapid rule is stated without it: every G0 ends on the clearance plane, and a G0
that moves in XY also starts on it. That leaves exactly two named exemptions
for a G0 that starts below clearance:

* the declared approach (`Program.approach_z`, a Z-only rapid from the unknown
  start position, which must itself end on the clearance plane), and
* the vertical retract: a Z-only G0 up to the clearance plane.

A Z-only G1 is either a plunge (descending onto the cut plane, counted
against `PostSetup.accepted_plunges`) or an air descent (descending to a
height at or above the stock top `STOCK_TOP_Z`, as before a HELIX ramp; not
a plunge, and a feed move, so no rapid exemption).

The exemption list is unchanged by the HELIX entry: a helix is feed motion.

A helical block is a G2/G3 whose Z changes. It must descend between the two
planes; each maximal run of helical blocks (the ramp of one HELIX entry) is
checked as a whole on the emitted decimals:

* structure (`HelixStructureError`): one sense, one exact centre, an even
  number of blocks whose every second block ends exactly at the run's start
  XY (each pair is then exactly one full turn: two same-sense sweeps in
  (0, 2 pi) that return to their start sum to 2 pi), ending on the cut plane;
* descent (`HelixDescentMismatchError`): with ``m`` blocks from ``Z_0`` to
  ``Z_m`` the uniform half-turn descent is ``delta = (Z_m - Z_0) / m``; every
  block descends (``d_k < 0``) and ``|d_k - delta| <= q``, ``q`` the output
  quantum. Derivation: the mapper writes ``Z_k = quantise(Z_0 + k delta)``
  (``Z_0`` and ``Z_m`` are grid values), so each emitted Z is within ``q/2``
  of the uniform helix and a difference of two is within ``q``;
* ramp (`HelixRampAngleError`): each full turn's emitted descent ``D`` on the
  smaller of its emitted squared radii ``r^2`` satisfies ``D <= exact(K) * r``,
  decided as ``D^2 <= K^2 r^2`` (both sides non-negative; no square root, no
  ``atan``). ``K`` is `PostSetup.ramp_cap`; the smaller radius is the steeper,
  conservative reading. Certified per full turn, which is exact because the
  turn sweeps exactly 2 pi. A single half-turn block is not certified on its
  own: its descent is within ``q`` of half the turn's (the descent bound), so
  over its length ``pi r`` its slope exceeds the turn's by at most
  ``q / (pi r)``, about 1e-5 at r = 3 mm against tan(2 deg) ~ 0.035.

Per block, a helical arc gets the same radius, envelope and tool-grown stock
checks as any arc (in XY; its Z is monotone between two checked end points).
"""

from fractions import Fraction

from compas_cgal.post.arc_geometry import Vector
from compas_cgal.post.arc_geometry import radius_within_tolerance
from compas_cgal.post.arc_geometry import sweep_contains
from compas_cgal.post.errors import ArcRadiusMismatchError
from compas_cgal.post.errors import DegenerateArcError
from compas_cgal.post.errors import HelixDescentMismatchError
from compas_cgal.post.errors import HelixRampAngleError
from compas_cgal.post.errors import HelixStructureError
from compas_cgal.post.errors import MalformedProgramError
from compas_cgal.post.errors import OffPlaneFeedError
from compas_cgal.post.errors import OutOfEnvelopeError
from compas_cgal.post.errors import PocketOutsideStockError
from compas_cgal.post.errors import UnacceptedPlungeError
from compas_cgal.post.errors import UnsafeRapidError
from compas_cgal.post.grid import quantum
from compas_cgal.post.machine import STOCK_TOP_Z
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post.program import Motion
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import Program
from compas_cgal.post.program import ProgramPoint
from compas_cgal.post.quantities import ProgramMillimetre

# Axis directions an arc can reach its X/Y extremes along, with the axis they bound
# and whether that is the axis's max side.
_EXTREMES = (
    (Vector(Fraction(1), Fraction(0)), "x", True),
    (Vector(Fraction(-1), Fraction(0)), "x", False),
    (Vector(Fraction(0), Fraction(1)), "y", True),
    (Vector(Fraction(0), Fraction(-1)), "y", False),
)


def _inside(point: ProgramPoint, envelope: TravelEnvelope) -> bool:
    return all(low <= value <= high for value, (low, high) in ((point.x, envelope.x), (point.y, envelope.y), (point.z, envelope.z)))


def _check_envelope_point(point: ProgramPoint, envelope: TravelEnvelope, index: int) -> None:
    if not _inside(point, envelope):
        raise OutOfEnvelopeError(f"motion {index}: ({point.x}, {point.y}, {point.z}) is outside the travel envelope.")


Range = tuple[ProgramMillimetre, ProgramMillimetre]


def _arc_extent_within(centre: Vector, to_start: Vector, to_end: Vector, clockwise: bool, radius_squared: Fraction, inflation: Fraction, box: tuple[Range, Range]) -> bool:
    """Whether the arc, grown by *inflation* (a tool radius, or 0), stays inside the XY *box* at each swept axis extreme.

    ``centre + (r + inflation) <= high`` is ``room = high - centre - inflation >= 0``
    and ``r**2 <= room**2``: exact, no square root.
    """
    for direction, axis, is_max in _EXTREMES:
        if not sweep_contains(to_start, to_end, direction, clockwise):
            continue
        offset = centre.x if axis == "x" else centre.y
        low, high = box[0] if axis == "x" else box[1]
        room = (high - offset if is_max else offset - low) - inflation
        if room < 0 or radius_squared > room * room:
            return False
    return True


def _check_stock_disc(point: ProgramPoint, setup: PostSetup, index: int) -> None:
    """The tool disc at a cut-plane point lies inside the blank in XY (closed sides)."""
    radius = setup.cutter.diameter / 2
    (x_low, x_high), (y_low, y_high) = setup.stock.x_range, setup.stock.y_range
    if not (x_low <= point.x - radius and point.x + radius <= x_high and y_low <= point.y - radius and point.y + radius <= y_high):
        raise PocketOutsideStockError(f"motion {index}: the tool disc at ({point.x}, {point.y}) leaves the {setup.stock.material} blank.")


def _check_arc(motion: Motion, start: ProgramPoint, setup: PostSetup, index: int) -> None:
    """Radius agreement, non-degeneracy, and the swept extent against the envelope and, tool-grown, the blank."""
    if motion.centre is None:
        raise MalformedProgramError(f"motion {index}: an arc without I/J.")
    centre_x, centre_y = start.x + motion.centre.i, start.y + motion.centre.j
    to_start = Vector(start.x - centre_x, start.y - centre_y)
    to_end = Vector(motion.target.x - centre_x, motion.target.y - centre_y)
    start_squared = to_start.x**2 + to_start.y**2
    end_squared = to_end.x**2 + to_end.y**2
    if start_squared == 0 or end_squared == 0 or (start.x, start.y) == (motion.target.x, motion.target.y):
        raise DegenerateArcError(f"motion {index}: zero radius, or an arc ending where it starts.")
    if not radius_within_tolerance(start_squared, end_squared, setup.arc_tolerance):
        raise ArcRadiusMismatchError(f"motion {index}: r_start^2 = {start_squared}, r_end^2 = {end_squared} differ by more than {setup.arc_tolerance} mm in radius.")
    clockwise = motion.kind is MotionKind.CLOCKWISE_ARC
    radius_squared = max(start_squared, end_squared)
    centre = Vector(centre_x, centre_y)
    if not _arc_extent_within(centre, to_start, to_end, clockwise, radius_squared, Fraction(0), (setup.envelope.x, setup.envelope.y)):
        raise OutOfEnvelopeError(f"motion {index}: the arc's swept extent leaves the travel envelope.")
    if not _arc_extent_within(centre, to_start, to_end, clockwise, radius_squared, setup.cutter.diameter / 2, (setup.stock.x_range, setup.stock.y_range)):
        raise PocketOutsideStockError(f"motion {index}: the arc's swept tool extent leaves the {setup.stock.material} blank.")


def _check_rapid(motion: Motion, start: ProgramPoint, setup: PostSetup, index: int) -> None:
    clearance = setup.clearance_z
    if motion.target.z != clearance:
        raise UnsafeRapidError(f"motion {index}: G0 ends at Z {motion.target.z}, not on the clearance plane {clearance}.")
    moves_xy = (start.x, start.y) != (motion.target.x, motion.target.y)
    if moves_xy and start.z != clearance:
        raise UnsafeRapidError(f"motion {index}: G0 moves in XY from Z {start.z}, below the clearance plane {clearance}.")


def _is_helical(motion: Motion, start: ProgramPoint) -> bool:
    """Whether a motion is a helical block: an arc whose Z changes."""
    return motion.kind is not MotionKind.LINEAR and motion.kind is not MotionKind.RAPID and start.z != motion.target.z


def _is_plunge(motion: Motion, start: ProgramPoint, setup: PostSetup, index: int) -> bool:
    """Whether a feed motion is a plunge; raise unless it is a plunge, a cut on the cut plane, or a descending helical block."""
    cut = setup.cut_z
    moves_xy = (start.x, start.y) != (motion.target.x, motion.target.y)
    if _is_helical(motion, start):
        if not cut <= motion.target.z < start.z <= setup.clearance_z:
            raise OffPlaneFeedError(f"motion {index}: a helical {motion.kind.value} must descend between the planes (Z {start.z} -> {motion.target.z}).")
        return False
    if motion.kind is not MotionKind.LINEAR or moves_xy:
        if start.z != cut or motion.target.z != cut:
            raise OffPlaneFeedError(f"motion {index}: {motion.kind.value} moves in XY off the cut plane {cut} (Z {start.z} -> {motion.target.z}).")
        return False
    if not start.z > motion.target.z:
        raise OffPlaneFeedError(f"motion {index}: a Z-only G1 must descend (Z {start.z} -> {motion.target.z}).")
    if motion.target.z == cut:
        return True
    if motion.target.z < STOCK_TOP_Z:
        raise OffPlaneFeedError(f"motion {index}: a Z-only G1 must end on the cut plane {cut} or in air at or above the stock top (Z {start.z} -> {motion.target.z}).")
    return False


def _arc_centre(motion: Motion, start: ProgramPoint) -> tuple[Fraction, Fraction]:
    if motion.centre is None:
        raise MalformedProgramError("an arc without I/J.")
    return start.x + motion.centre.i, start.y + motion.centre.j


def _squared_radius(point: ProgramPoint, centre: tuple[Fraction, Fraction]) -> Fraction:
    return (point.x - centre[0]) ** 2 + (point.y - centre[1]) ** 2


def _check_helix_run(blocks: list[tuple[ProgramPoint, Motion]], first: int, setup: PostSetup) -> None:
    """Structure, uniform descent and ramp of one run of helical blocks (module docstring).

    Args:
        blocks: The run's (start, motion) pairs, in order.
        first: Motion index of the run's first block.
        setup: The setup.

    Raises:
        HelixStructureError: Mixed senses or centres, a part turn, or a run ending off the cut plane.
        HelixDescentMismatchError: A block not descending, or off the uniform pitch by more than a quantum.
        HelixRampAngleError: A turn steeper than ``setup.ramp_cap``.
    """
    origin = blocks[0][0]
    centre = _arc_centre(blocks[0][1], origin)
    if len(blocks) % 2:
        raise HelixStructureError(f"motions {first}+: {len(blocks)} helical half-turns are not whole turns.")
    for offset, (start, motion) in enumerate(blocks):
        if motion.kind is not blocks[0][1].kind or _arc_centre(motion, start) != centre:
            raise HelixStructureError(f"motion {first + offset}: the helix changes sense or centre.")
        if offset % 2 and (motion.target.x, motion.target.y) != (origin.x, origin.y):
            raise HelixStructureError(f"motion {first + offset}: a helical turn does not return to its start.")
    bottom = blocks[-1][1].target.z
    if bottom != setup.cut_z:
        raise HelixStructureError(f"motions {first}+: the helix ends at Z {bottom}, not on the cut plane {setup.cut_z}.")
    delta = (bottom - origin.z) / len(blocks)
    step = quantum(setup.decimals)
    for offset, (start, motion) in enumerate(blocks):
        descent = motion.target.z - start.z
        if not descent < 0 or abs(descent - delta) > step:
            raise HelixDescentMismatchError(f"motion {first + offset}: Z step {descent} off the uniform {delta} by more than {step} mm.")
    cap = setup.ramp_cap.exact
    for turn in range(0, len(blocks), 2):
        start, end = blocks[turn][0], blocks[turn + 1][1].target
        drop = start.z - end.z
        radius_squared = min(_squared_radius(start, centre), _squared_radius(blocks[turn][1].target, centre))
        if drop * drop > cap * cap * radius_squared:
            raise HelixRampAngleError(f"motion {first + turn}: a turn drops {drop} mm on r^2 = {radius_squared} mm^2, beyond the ramp cap K = {cap}.")


def check_program(program: Program, setup: PostSetup) -> None:
    """Refuse *program* unless every check holds.

    Args:
        program: A parsed program.
        setup: The setup it was posted for.

    Raises:
        UnsafeRapidError: A rapid outside the rapid rule and its two exemptions.
        OffPlaneFeedError: A feed motion that is neither a cut on the cut plane nor a plunge onto it.
        DegenerateArcError: A zero-radius or closed arc.
        ArcRadiusMismatchError: An arc's radii differ by more than the setup's arc tolerance.
        OutOfEnvelopeError: A point or an arc's extent leaves the travel envelope.
        PocketOutsideStockError: A feed motion's tool disc, or an arc's tool-grown extent, leaves the blank in XY.
        UnacceptedPlungeError: More plunges than `PostSetup.accepted_plunges`.
        HelixStructureError: A helical run that is not whole same-sense turns about one centre onto the cut plane.
        HelixDescentMismatchError: A helical block off the run's uniform descent by more than one quantum.
        HelixRampAngleError: A helical turn steeper than `PostSetup.ramp_cap`.
    """
    clearance = setup.clearance_z
    if program.approach_z != clearance:
        raise UnsafeRapidError(f"the approach rapid ends at Z {program.approach_z}, not on the clearance plane {clearance}.")
    if not program.motions or program.motions[0].kind is not MotionKind.RAPID:
        raise UnsafeRapidError("the first motion after the approach must be the rapid to clearance above the entry.")
    low, high = setup.envelope.z
    if not low <= program.approach_z <= high:
        raise OutOfEnvelopeError(f"the approach Z {program.approach_z} is outside the travel envelope.")
    first = program.motions[0]
    _check_rapid(first, ProgramPoint(first.target.x, first.target.y, program.approach_z), setup, 0)
    _check_envelope_point(first.target, setup.envelope, 0)
    plunges = 0
    runs: list[tuple[int, list[tuple[ProgramPoint, Motion]]]] = []
    for index in range(1, len(program.motions)):
        start, motion = program.motions[index - 1].target, program.motions[index]
        _check_envelope_point(motion.target, setup.envelope, index)
        if motion.kind is MotionKind.RAPID:
            _check_rapid(motion, start, setup, index)
            continue
        if _is_plunge(motion, start, setup, index):
            plunges += 1
        _check_stock_disc(motion.target, setup, index)
        if motion.kind is not MotionKind.LINEAR:
            _check_arc(motion, start, setup, index)
        if _is_helical(motion, start):
            if runs and runs[-1][0] + len(runs[-1][1]) == index:
                runs[-1][1].append((start, motion))
            else:
                runs.append((index, [(start, motion)]))
    for run_start, blocks in runs:
        _check_helix_run(blocks, run_start, setup)
    if plunges > setup.accepted_plunges:
        raise UnacceptedPlungeError(f"{plunges} plunges, {setup.accepted_plunges} accepted.")
