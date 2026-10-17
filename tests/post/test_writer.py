"""The writer dispatches exhaustively on dialect and arc mode, and emits the full modal safety block."""

from enum import Enum
from fractions import Fraction

import pytest

from compas_cgal.post.errors import UnsupportedArcModeError
from compas_cgal.post.errors import UnsupportedDialectError
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.program import Program
from compas_cgal.post.quantities import RevolutionsPerMinute
from compas_cgal.post.quantities import Seconds
from compas_cgal.post.writer import write_program
from tests.post.support import mm
from tests.post.support import rebuild

EMPTY = Program(spindle_rpm=RevolutionsPerMinute(10000), dwell=Seconds(Fraction(2)), approach_z=mm(12), motions=())


class FanucDialect(Enum):
    """A controller the writer has no dialect for."""

    FANUC = "fanuc"


class AbsoluteCentres(Enum):
    """An arc mode the writer has no mapping for."""

    XY_ABSOLUTE_IJ = "G17 G90.1"


def test_header_sets_every_modal_group(setup: PostSetup) -> None:
    header = write_program(EMPTY, setup).splitlines()[0].split()
    assert header == ["G21", "G90", "G94", "G40", "G80", "G17", "G91.1", "G54"]


def test_spindle_up_dwell_precedes_any_motion(setup: PostSetup) -> None:
    assert write_program(EMPTY, setup).splitlines()[1:4] == ["S10000 M3", "G4 P2.0000", "G0 Z12.0000"]


def test_unknown_dialect_raises(setup: PostSetup) -> None:
    with pytest.raises(UnsupportedDialectError):
        write_program(EMPTY, rebuild(setup, dialect=FanucDialect.FANUC))


def test_unknown_arc_mode_raises(setup: PostSetup) -> None:
    with pytest.raises(UnsupportedArcModeError):
        write_program(EMPTY, rebuild(setup, arc_mode=AbsoluteCentres.XY_ABSOLUTE_IJ))
