"""Read back exactly the dialect the writer emits, into exact `Fraction` coordinates.

Two failure modes, checked in this order over the whole text: a word outside
the dialect (`UnsupportedWordError`), then dialect words in a structure the
writer never produces (`MalformedProgramError`). Nothing is skipped, defaulted
or repaired.
"""

import re

from compas_cgal.post.errors import MalformedProgramError
from compas_cgal.post.errors import UnsupportedWordError
from compas_cgal.post.grid import parse_decimal
from compas_cgal.post.grid import parse_integer
from compas_cgal.post.program import DWELL
from compas_cgal.post.program import DWELL_LETTER
from compas_cgal.post.program import HEADER_WORDS
from compas_cgal.post.program import PROGRAM_END
from compas_cgal.post.program import SPINDLE_OFF
from compas_cgal.post.program import SPINDLE_ON
from compas_cgal.post.program import CentreOffset
from compas_cgal.post.program import Motion
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import Program
from compas_cgal.post.program import ProgramPoint
from compas_cgal.post.quantities import MillimetrePerMinute
from compas_cgal.post.quantities import ProgramMillimetre
from compas_cgal.post.quantities import RevolutionsPerMinute
from compas_cgal.post.quantities import Seconds

_WORD = re.compile(r"([A-Z])(\S+)")

# Code words the dialect knows, by letter.
_CODES = {
    "G": frozenset(word[1:] for word in (*HEADER_WORDS, DWELL)) | frozenset(kind.value[1:] for kind in MotionKind),
    "M": frozenset(code[1:] for code in (SPINDLE_ON, SPINDLE_OFF, PROGRAM_END)),
}

# Value words the dialect knows.
_VALUE_LETTERS = frozenset("SFXYZIJ" + DWELL_LETTER)

# Letter sequence of each motion block after its G word.
_MOTION_LETTERS = {
    MotionKind.RAPID: "XYZ",
    MotionKind.LINEAR: "XYZF",
    MotionKind.CLOCKWISE_ARC: "XYZIJF",
    MotionKind.COUNTERCLOCKWISE_ARC: "XYZIJF",
}

# Blocks around the motions: header, spindle on, dwell, approach before; spindle off, end after.
_LEADING_BLOCKS = 4
_TRAILING_BLOCKS = 2

Word = tuple[str, str]


def _words(line: str, number: int) -> list[Word]:
    words: list[Word] = []
    for token in line.split():
        match = _WORD.fullmatch(token)
        if match is None:
            raise UnsupportedWordError(f"line {number}: {token!r} is not a word of the dialect.")
        letter, value = match.group(1), match.group(2)
        if (letter in _CODES and value not in _CODES[letter]) or (letter not in _CODES and letter not in _VALUE_LETTERS):
            raise UnsupportedWordError(f"line {number}: {token!r} is outside the emitted dialect.")
        words.append((letter, value))
    if not words:
        raise MalformedProgramError(f"line {number} is empty.")
    return words


def _codes(words: list[Word]) -> tuple[str, ...]:
    return tuple(letter + value for letter, value in words)


def _motion(words: list[Word], number: int) -> Motion:
    head = words[0][0] + words[0][1]
    kinds = {kind.value: kind for kind in MotionKind}
    if head not in kinds:
        raise MalformedProgramError(f"line {number} does not start with a motion word.")
    kind = kinds[head]
    letters = "".join(letter for letter, _ in words[1:])
    if letters != _MOTION_LETTERS[kind]:
        raise MalformedProgramError(f"line {number}: {head} carries {letters!r}, the writer emits {_MOTION_LETTERS[kind]!r}.")
    values = {letter: value for letter, value in words[1:]}
    target = ProgramPoint(*(ProgramMillimetre(parse_decimal(values[axis])) for axis in "XYZ"))
    centre = CentreOffset(ProgramMillimetre(parse_decimal(values["I"])), ProgramMillimetre(parse_decimal(values["J"]))) if "I" in values else None
    feed = MillimetrePerMinute(parse_integer(values["F"])) if "F" in values else None
    return Motion(kind, target, centre, feed)


def parse_program(text: str) -> Program:
    """Parse a program the writer could have produced.

    Args:
        text: The program text.

    Returns:
        The program.

    Raises:
        UnsupportedWordError: A word outside the dialect.
        MalformedProgramError: The blocks are not header, spindle on, dwell, approach,
            motions, spindle off, program end, each in the writer's form.
    """
    blocks = [_words(line, number) for number, line in enumerate(text.splitlines(), start=1)]
    if len(blocks) < _LEADING_BLOCKS + _TRAILING_BLOCKS:
        raise MalformedProgramError(f"{len(blocks)} blocks cannot hold header, spindle, approach and program end.")
    header, spindle, dwell, approach = blocks[0], blocks[1], blocks[2], blocks[3]
    if _codes(header) != HEADER_WORDS:
        raise MalformedProgramError(f"header {_codes(header)} is not {HEADER_WORDS}.")
    if len(spindle) != 2 or spindle[0][0] != "S" or _codes(spindle[1:]) != (SPINDLE_ON,):
        raise MalformedProgramError(f"block 2 {_codes(spindle)} is not 'S<rpm> {SPINDLE_ON}'.")
    if len(dwell) != 2 or _codes(dwell[:1]) != (DWELL,) or dwell[1][0] != DWELL_LETTER:
        raise MalformedProgramError(f"block 3 {_codes(dwell)} is not the spindle-up dwell '{DWELL} {DWELL_LETTER}<s>'.")
    if len(approach) != 2 or _codes(approach[:1]) != (MotionKind.RAPID.value,) or approach[1][0] != "Z":
        raise MalformedProgramError(f"block 4 {_codes(approach)} is not the Z-only approach rapid.")
    if [_codes(block) for block in blocks[-_TRAILING_BLOCKS:]] != [(SPINDLE_OFF,), (PROGRAM_END,)]:
        raise MalformedProgramError(f"the program does not end with {SPINDLE_OFF} then {PROGRAM_END}.")
    motions = tuple(_motion(block, number) for number, block in enumerate(blocks[_LEADING_BLOCKS:-_TRAILING_BLOCKS], start=_LEADING_BLOCKS + 1))
    return Program(
        RevolutionsPerMinute(parse_integer(spindle[0][1])),
        Seconds(parse_decimal(dwell[1][1])),
        ProgramMillimetre(parse_decimal(approach[1][1])),
        motions,
    )
