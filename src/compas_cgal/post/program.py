"""The emitted program as data: the vocabulary the mapper builds, the writer prints and the parser reads back.

The program is stateless per block: every motion carries X, Y and Z, every feed
motion carries F, every arc carries I and J. The one exception is the declared
approach, a Z-only rapid from the unknown start position, held apart as
`Program.approach_z`.
"""

from dataclasses import dataclass
from enum import Enum

from compas_cgal.post.quantities import MillimetrePerMinute
from compas_cgal.post.quantities import ProgramMillimetre
from compas_cgal.post.quantities import RevolutionsPerMinute
from compas_cgal.post.quantities import Seconds

# LinuxCNC modal safety words, in emission order: mm, absolute distance, feed
# per minute, cutter compensation off, canned cycle cancelled. Every modal group
# a motion depends on is set, so nothing rests on what an earlier program left.
LINUXCNC_MODAL_WORDS: tuple[str, ...] = ("G21", "G90", "G94", "G40", "G80")

# Arc plane and centre addressing: XY plane, incremental I/J.
XY_INCREMENTAL_ARC_WORDS: tuple[str, ...] = ("G17", "G91.1")

# The work offset word (Q9).
G54_WORD = "G54"

# The one header the parser accepts: the LinuxCNC / G17 G91.1 / G54 setup.
HEADER_WORDS: tuple[str, ...] = (*LINUXCNC_MODAL_WORDS, *XY_INCREMENTAL_ARC_WORDS, G54_WORD)

# Dwell word (G4) and its seconds letter (P).
DWELL = "G4"
DWELL_LETTER = "P"

# Spindle on clockwise; spindle off; program end.
SPINDLE_ON = "M3"
SPINDLE_OFF = "M5"
PROGRAM_END = "M30"


class MotionKind(Enum):
    """The four motion words."""

    RAPID = "G0"
    LINEAR = "G1"
    CLOCKWISE_ARC = "G2"
    COUNTERCLOCKWISE_ARC = "G3"


@dataclass(frozen=True)
class ProgramPoint:
    """A tool-centre position in the work frame, mm."""

    x: ProgramMillimetre
    y: ProgramMillimetre
    z: ProgramMillimetre


@dataclass(frozen=True)
class CentreOffset:
    """An arc centre relative to the arc's start (G91.1 I/J), mm."""

    i: ProgramMillimetre
    j: ProgramMillimetre


@dataclass(frozen=True)
class Motion:
    """One motion block.

    Attributes:
        kind: The motion word.
        target: Where it ends.
        centre: I/J from the block's start; arcs only.
        feed: Feed rate; feed motions (G1, G2, G3) only.
    """

    kind: MotionKind
    target: ProgramPoint
    centre: CentreOffset | None
    feed: MillimetrePerMinute | None


@dataclass(frozen=True)
class Program:
    """A whole program between the header and the program end.

    Attributes:
        spindle_rpm: The S word of the spindle-on block.
        dwell: The spindle-up dwell before any motion.
        approach_z: The declared approach: a Z-only rapid from the unknown start.
        motions: The motion blocks, in order.
    """

    spindle_rpm: RevolutionsPerMinute
    dwell: Seconds
    approach_z: ProgramMillimetre
    motions: tuple[Motion, ...]


@dataclass(frozen=True)
class OperationSpan:
    """The motions one source operation became.

    Attributes:
        operation_index: Index into the source operations.
        first_motion: Index of its first motion in `Program.motions`.
        motion_count: How many motions it became.
    """

    operation_index: int
    first_motion: int
    motion_count: int


@dataclass(frozen=True)
class PostedProgram:
    """A mapped program and its trace back to the source operations."""

    program: Program
    spans: tuple[OperationSpan, ...]
