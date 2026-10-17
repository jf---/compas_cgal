"""Print a `Program` as G-code text in the emitted dialect, exactly (no rounding happens here).

The writer prints faithfully; the grammar (which words a block may carry) is
enforced where text is read back, in `compas_cgal.post.parser`.
"""

from typing import Never
from typing import NoReturn

from compas_cgal.post.errors import UnsupportedArcModeError
from compas_cgal.post.errors import UnsupportedDialectError
from compas_cgal.post.grid import format_decimal
from compas_cgal.post.machine import ArcCentreMode
from compas_cgal.post.machine import Dialect
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.program import DWELL
from compas_cgal.post.program import DWELL_LETTER
from compas_cgal.post.program import LINUXCNC_MODAL_WORDS
from compas_cgal.post.program import PROGRAM_END
from compas_cgal.post.program import SPINDLE_OFF
from compas_cgal.post.program import SPINDLE_ON
from compas_cgal.post.program import XY_INCREMENTAL_ARC_WORDS
from compas_cgal.post.program import Motion
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import Program


def _unsupported_dialect(dialect: Never) -> NoReturn:
    raise UnsupportedDialectError(f"no writer for dialect {dialect!r}.")


def _unsupported_arc_mode(mode: Never) -> NoReturn:
    raise UnsupportedArcModeError(f"no writer for arc mode {mode!r}.")


def _modal_words(dialect: Dialect) -> tuple[str, ...]:
    # Exhaustive by identity: mypy narrows the fall-through to Never.
    if dialect is Dialect.LINUXCNC:
        return LINUXCNC_MODAL_WORDS
    _unsupported_dialect(dialect)


def _arc_words(mode: ArcCentreMode) -> tuple[str, ...]:
    if mode is ArcCentreMode.XY_INCREMENTAL_IJ:
        return XY_INCREMENTAL_ARC_WORDS
    _unsupported_arc_mode(mode)


def _motion_line(motion: Motion, setup: PostSetup) -> str:
    digits = setup.decimals
    target = motion.target
    words = [motion.kind.value, f"X{format_decimal(target.x, digits)}", f"Y{format_decimal(target.y, digits)}", f"Z{format_decimal(target.z, digits)}"]
    if motion.centre is not None:
        words += [f"I{format_decimal(motion.centre.i, digits)}", f"J{format_decimal(motion.centre.j, digits)}"]
    if motion.feed is not None:
        words.append(f"F{motion.feed}")
    return " ".join(words)


def write_program(program: Program, setup: PostSetup) -> str:
    """The program text, one block per line, newline-terminated.

    Args:
        program: The program; every coordinate already on the setup's grid.
        setup: Supplies the dialect, arc mode, work offset and grid digits.

    Returns:
        The text.

    Raises:
        OffGridValueError: A coordinate is not on the grid.
        UnsupportedDialectError: A dialect other than LinuxCNC.
        UnsupportedArcModeError: An arc mode other than G17 with incremental I/J.
    """
    header = " ".join((*_modal_words(setup.dialect), *_arc_words(setup.arc_mode), setup.work_offset.value))
    lines = [
        header,
        f"S{program.spindle_rpm} {SPINDLE_ON}",
        f"{DWELL} {DWELL_LETTER}{format_decimal(program.dwell, setup.decimals)}",
        f"{MotionKind.RAPID.value} Z{format_decimal(program.approach_z, setup.decimals)}",
    ]
    lines += [_motion_line(motion, setup) for motion in program.motions]
    lines += [SPINDLE_OFF, PROGRAM_END]
    return "\n".join(lines) + "\n"
