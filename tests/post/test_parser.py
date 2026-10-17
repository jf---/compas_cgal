"""The parser reads only the emitted dialect and raises on anything else."""

from fractions import Fraction

import pytest

from compas_cgal.post.errors import MalformedProgramError
from compas_cgal.post.errors import UnsupportedWordError
from compas_cgal.post.parser import parse_program
from compas_cgal.post.program import MotionKind

GOOD = """G21 G90 G94 G40 G80 G17 G91.1 G54
S10000 M3
G4 P2.0000
G0 Z12.0000
G0 X1.0000 Y2.0000 Z12.0000
G1 X1.0000 Y2.0000 Z-1.5000 F100
G2 X3.0000 Y2.0000 Z-1.5000 I1.0000 J0.0000 F600
G0 X3.0000 Y2.0000 Z12.0000
M5
M30
"""


def test_parses_the_emitted_dialect() -> None:
    program = parse_program(GOOD)
    assert program.spindle_rpm == 10000
    assert program.approach_z == 12
    assert program.dwell == 2
    assert [m.kind for m in program.motions] == [MotionKind.RAPID, MotionKind.LINEAR, MotionKind.CLOCKWISE_ARC, MotionKind.RAPID]
    arc = program.motions[2]
    assert arc.centre is not None and (arc.centre.i, arc.centre.j) == (Fraction(1), Fraction(0)) and arc.feed == 600


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("G21 G90", "G20 G90"),  # inch mode
        ("G94", "G93"),  # inverse-time feed
        ("G40", "G41"),  # cutter compensation on
        ("G91.1", "G90.1"),  # absolute arc centres
        ("S10000 M3", "S10000 M4"),  # reverse spindle
        ("G1 X1.0000", "N10 G1 X1.0000"),  # block numbers
        ("M5", "M8"),  # coolant
        ("G0 Z12.0000\n", "G0 Z12.0000\nT1\n"),  # tool change
        ("G2 X3.0000", "G2 X3.0000 (c)"),  # comment
    ],
)
def test_any_other_word_raises(old: str, new: str) -> None:
    with pytest.raises(UnsupportedWordError):
        parse_program(GOOD.replace(old, new, 1))


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("G21 G90 G94 G40 G80 G17 G91.1 G54", "G21 G90 G17 G91.1 G54"),  # no feed mode, compensation or canned-cycle cancel
        ("G4 P2.0000\n", ""),  # no spindle-up dwell
        ("G2 X3.0000 Y2.0000 Z-1.5000 I1.0000 J0.0000 F600", "G2 X3.0000 Y2.0000 Z-1.5000 F600"),  # arc without centre
        ("G1 X1.0000 Y2.0000 Z-1.5000 F100", "G1 X1.0000 Y2.0000 Z-1.5000"),  # feed block without F
        ("G0 X3.0000 Y2.0000 Z12.0000", "G0 X3.0000 Y2.0000 Z12.0000 F600"),  # rapid with F
        ("G1 X1.0000 Y2.0000", "G1 Y2.0000 X1.0000"),  # word order
        ("M30\n", ""),  # no program end
        ("X1.0000 Y2.0000 Z12.0000", "X1.0000 Y2.0000 Z1e1"),  # exponent
    ],
)
def test_malformed_structure_raises(old: str, new: str) -> None:
    with pytest.raises(MalformedProgramError):
        parse_program(GOOD.replace(old, new, 1))
