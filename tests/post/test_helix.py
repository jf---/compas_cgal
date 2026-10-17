"""The HELIX entry: helical G2/G3 ramp turns plus one flat closing turn, round-tripped writer -> parser -> checks.

The helix comes from `compas_cgal.pipeline.emit.helix_operation`, the emitter
the generator uses. Source planes: cut 0, clearance 4, tool radius 1 (x3 onto
the test-cut setup: clearance 12, cut -1.5, 6 mm cutter). Footprint radius 1 ->
3 mm in the program. The post feeds down in air to the ramp start (1 mm above
the top) and ramps 2.5 mm from there; at the 2-degree placeholder cap
``K = 2*pi*tan(2 deg)`` that is ``ceil(2.5 / (3 K)) = ceil(3.79..) = 4`` turns,
decided in program units whatever the source helix's count.
"""

import dataclasses
import math
from fractions import Fraction
from typing import Any

import pytest
from compas.geometry import Line

from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import helix_operation
from compas_cgal.post.checks import check_program
from compas_cgal.post.errors import HelixDescentMismatchError
from compas_cgal.post.errors import HelixOperationMismatchError
from compas_cgal.post.errors import HelixRampAngleError
from compas_cgal.post.errors import HelixStructureError
from compas_cgal.post.errors import RampStartAboveClearanceError
from compas_cgal.post.errors import MalformedProgramError
from compas_cgal.post.errors import NonPositiveQuantityError
from compas_cgal.post.errors import OutOfEnvelopeError
from compas_cgal.post.errors import PocketOutsideStockError
from compas_cgal.post.errors import UnsupportedOperationError
from compas_cgal.post.errors import UnsupportedWordError
from compas_cgal.post.errors import WindingDisagreementError
from compas_cgal.post.fail_closed import post_operations
from compas_cgal.post.grid import quantise
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import StockBlank
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post import mapping
from compas_cgal.post.mapping import SourcePath
from compas_cgal.post.mapping import map_operations
from compas_cgal.post.parser import parse_program
from compas_cgal.post.program import Motion
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import Program
from compas_cgal.post.program import ProgramPoint
from compas_cgal.post.quantities import ProgramMillimetre
from compas_cgal.post.testcut_defaults import TESTCUT_RAMP_FEED
from compas_cgal.post.writer import write_program
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.post.support import mm
from tests.post.support import rebuild

SOURCE = SourcePath(cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0), tool=ToolRadius.build(1.0))

# Program-unit ramp turns at the defaults (module docstring).
PROGRAM_TURNS = 4

# `Helix.build`'s count for the same footprint in SOURCE units: ceil(4 / K) = 19.
SOURCE_UNIT_TURNS = 19

# A ramp start 11.5 mm above the top: a 13 mm program ramp needs ceil(13 / (3 K)) = 20 turns, so the
# source count (19) is too steep there. Heights are not scaled by the post while radii scale x3.
HIGH_RAMP_START = Fraction(23, 2)

# Entry directions: on an axis, and generic (an irrational-looking double phase).
PHASES = (0.0, 2.0)

# A Z perturbation far above the output quantum (0.0001 mm) and far below a half-turn's descent.
Z_DEFECT = Fraction(1, 100)


def _helix(clockwise: bool, phase: float, turns: int = SOURCE_UNIT_TURNS) -> Helix:
    footprint = IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(1.0), Radian(phase))
    return Helix(footprint=footprint, top_z=SOURCE.clearance, bottom_z=SOURCE.cut, turns=turns, clockwise=clockwise)


def _stream(helix: Helix) -> list[ToolpathOperation]:
    entry = helix.entry_xy
    retract = ToolpathOperation(geometry=Line([entry.x, entry.y, 0.0], [entry.x, entry.y, 4.0]), operation=OperationType.RETRACT, path_index=0)
    return [helix_operation(helix, 0), retract]


def _parsed(helix: Helix, setup: PostSetup) -> Program:
    return parse_program(post_operations(_stream(helix), setup, SOURCE))


def _with_motion(program: Program, index: int, motion: Motion) -> Program:
    motions = list(program.motions)
    motions[index] = motion
    return dataclasses.replace(program, motions=tuple(motions))


@pytest.mark.parametrize("phase", PHASES)
@pytest.mark.parametrize("clockwise", [True, False])
def test_helix_round_trips_writer_parser_checks(setup: PostSetup, clockwise: bool, phase: float) -> None:
    helix = _helix(clockwise, phase)
    posted = map_operations(_stream(helix), setup, SOURCE)
    text = post_operations(_stream(helix), setup, SOURCE)
    assert parse_program(text) == posted.program
    motions = posted.program.motions
    halves = 2 * PROGRAM_TURNS
    shift_x, shift_y = setup.xy_shift
    entry = (mm(quantise(3 * Fraction(helix.entry_xy.x) + shift_x, setup.decimals)), mm(quantise(3 * Fraction(helix.entry_xy.y) + shift_y, setup.decimals)))
    antipode = (mm(quantise(-3 * Fraction(helix.entry_xy.x) + shift_x, setup.decimals)), mm(quantise(-3 * Fraction(helix.entry_xy.y) + shift_y, setup.decimals)))
    assert motions[0] == Motion(MotionKind.RAPID, ProgramPoint(entry[0], entry[1], setup.clearance_z), None, None)
    assert setup.ramp_start_z == 1
    assert motions[1] == Motion(MotionKind.LINEAR, ProgramPoint(entry[0], entry[1], setup.ramp_start_z), None, setup.cutting.plunge_feed)
    kind = MotionKind.CLOCKWISE_ARC if clockwise else MotionKind.COUNTERCLOCKWISE_ARC
    ramp, closing = motions[2 : 2 + halves], motions[2 + halves : 4 + halves]
    assert all(m.kind is kind and m.feed == setup.ramp_feed for m in (*ramp, *closing))
    descent = setup.cut_z - setup.ramp_start_z
    assert [m.target.z for m in ramp] == [quantise(setup.ramp_start_z + descent * Fraction(k, halves), setup.decimals) for k in range(1, halves + 1)]
    assert [(m.target.x, m.target.y) for m in ramp] == [antipode, entry] * PROGRAM_TURNS
    assert [m.target for m in closing] == [ProgramPoint(antipode[0], antipode[1], setup.cut_z), ProgramPoint(entry[0], entry[1], setup.cut_z)]
    assert motions[4 + halves] == Motion(MotionKind.RAPID, ProgramPoint(entry[0], entry[1], setup.clearance_z), None, None)
    assert [(s.operation_index, s.first_motion, s.motion_count) for s in posted.spans] == [(0, 1, halves + 3), (1, halves + 4, 1)]


def test_program_turns_follow_the_cap_not_the_source_count(setup: PostSetup) -> None:
    for source_turns in (1, SOURCE_UNIT_TURNS, 40):
        assert mapping.program_ramp_turns(_helix(True, 0.0, source_turns), setup, setup.ramp_start_z, setup.cut_z) == PROGRAM_TURNS


def test_a_helix_entry_needs_no_plunge(setup: PostSetup) -> None:
    # The one Z-only G1 is the air descent to the ramp start, not a plunge.
    program = _parsed(_helix(True, 0.0), rebuild(setup, accepted_plunges=0))
    assert [m.target.z for m in program.motions if m.kind is MotionKind.LINEAR] == [setup.ramp_start_z]


@pytest.mark.parametrize("start", [Fraction(0), Fraction(-1)])
def test_ramp_start_below_the_top_is_refused(setup: PostSetup, start: Fraction) -> None:
    with pytest.raises(NonPositiveQuantityError):
        rebuild(setup, ramp_start_above_top=start)


@pytest.mark.parametrize("start", [Fraction(12), Fraction(13)])
def test_ramp_start_at_or_above_clearance_is_refused(setup: PostSetup, start: Fraction) -> None:
    with pytest.raises(RampStartAboveClearanceError):
        rebuild(setup, ramp_start_above_top=start)


def test_ramp_feed_is_a_placeholder_setup_field(setup: PostSetup) -> None:
    assert setup.ramp_feed == TESTCUT_RAMP_FEED
    with pytest.raises(NonPositiveQuantityError):
        rebuild(setup, ramp_feed=0)


@pytest.mark.parametrize("block", [2, 5, 2 * PROGRAM_TURNS])  # ramp blocks are motions 2..9; 9 ends on the cut plane
def test_mutated_helix_z_is_refused(setup: PostSetup, block: int) -> None:
    program = _parsed(_helix(True, 2.0), setup)
    motion = program.motions[block]
    lifted = dataclasses.replace(motion, target=dataclasses.replace(motion.target, z=ProgramMillimetre(motion.target.z + Z_DEFECT)))
    with pytest.raises(HelixDescentMismatchError):
        check_program(parse_program(write_program(_with_motion(program, block, lifted), setup)), setup)


def test_flipped_helical_block_is_refused(setup: PostSetup) -> None:
    program = _parsed(_helix(True, 0.0), setup)
    flipped = dataclasses.replace(program.motions[3], kind=MotionKind.COUNTERCLOCKWISE_ARC)
    with pytest.raises(HelixStructureError):
        check_program(_with_motion(program, 3, flipped), setup)


@pytest.mark.parametrize("clockwise", [True, False])
def test_flipped_operation_sense_is_refused(setup: PostSetup, clockwise: bool) -> None:
    stream = _stream(_helix(clockwise, 2.0))
    stream[0].clockwise = not clockwise
    with pytest.raises(WindingDisagreementError):
        map_operations(stream, setup, SOURCE)


def test_operation_disagreeing_with_its_helix_is_refused(setup: PostSetup) -> None:
    stream = _stream(_helix(True, 2.0))
    stream[0].geometry = _helix(False, 2.0)
    with pytest.raises(HelixOperationMismatchError):
        map_operations(stream, setup, SOURCE)


def test_operation_entry_off_its_helix_is_refused(setup: PostSetup) -> None:
    stream = _stream(_helix(True, 2.0))
    assert stream[0].entry_xy is not None
    stream[0].entry_xy = (stream[0].entry_xy[0], math.nextafter(stream[0].entry_xy[1], math.inf))
    with pytest.raises(HelixOperationMismatchError):
        map_operations(stream, setup, SOURCE)


def test_source_turn_count_ramps_too_steep_in_the_program(setup: PostSetup, monkeypatch: pytest.MonkeyPatch) -> None:
    # Why the recount exists: the source helix's own `Helix.build` count, forced onto a 13 mm program ramp.
    source = Helix.build(_helix(True, 0.0).footprint, CutPlane.build(SOURCE.cut, SOURCE.clearance), setup.ramp_cap, True)
    assert source.turns == SOURCE_UNIT_TURNS
    high = rebuild(setup, ramp_start_above_top=HIGH_RAMP_START)
    post_operations(_stream(source), high, SOURCE)

    def source_count(helix: Helix, *_: Any) -> int:
        return helix.turns

    monkeypatch.setattr(mapping, "program_ramp_turns", source_count)
    with pytest.raises(HelixRampAngleError):
        post_operations(_stream(source), high, SOURCE)


def test_ramp_is_decided_on_the_setup_surrogate(setup: PostSetup) -> None:
    # K = 0.2: ceil(2.5 / 0.6) = 5 turns, each 0.5 mm on r = 3 mm.
    program = _parsed(_helix(True, 0.0), rebuild(setup, ramp_cap=HelixPitchPerRadius.build(0.2)))
    assert sum(m.kind is MotionKind.CLOCKWISE_ARC for m in program.motions) == 2 * 5 + 2  # ramp + closing half-turns


def test_helix_bulge_outside_the_envelope_is_refused(setup: PostSetup) -> None:
    # Phase 0: every emitted point has y = 12; the arcs bulge to y = 9 and 15.
    box = TravelEnvelope.build(x=(Fraction(0), Fraction(40)), y=(Fraction(5), Fraction(14)), z=(Fraction(-5), Fraction(20)))
    with pytest.raises(OutOfEnvelopeError):
        post_operations(_stream(_helix(True, 0.0)), rebuild(setup, envelope=box), SOURCE)


def test_helix_tool_extent_outside_the_blank_is_refused(setup: PostSetup) -> None:
    # Phase 0: tool-grown arcs reach y = 15 + 3 = 18; a blank ending at y = 17 cannot hold them.
    tight = StockBlank.build(material="wax", min_corner=(Fraction(-12), Fraction(-13)), size_x=Fraction(60), size_y=Fraction(30), size_z=Fraction(15))
    with pytest.raises(PocketOutsideStockError):
        post_operations(_stream(_helix(True, 0.0)), rebuild(setup, stock=tight), SOURCE)


def test_entry_plunge_must_still_be_a_line(setup: PostSetup) -> None:
    stream = _stream(_helix(True, 0.0))
    stream[0].operation = OperationType.PLUNGE
    with pytest.raises(UnsupportedOperationError):
        map_operations(stream, setup, SOURCE)


HELICAL = """G21 G90 G94 G40 G80 G17 G91.1 G54
S10000 M3
G4 P2.0000
G0 Z12.0000
G0 X21.0000 Y12.0000 Z12.0000
G2 X15.0000 Y12.0000 Z11.6786 I-3.0000 J0.0000 F300
M5
M30
"""


def test_parser_reads_a_helical_block() -> None:
    motion = parse_program(HELICAL).motions[1]
    assert motion.kind is MotionKind.CLOCKWISE_ARC and motion.target.z == Fraction("11.6786")


@pytest.mark.parametrize("new", ["I-3.0000 J0.0000 K0.5000 F300", "I-3.0000 J0.0000 R3.0000 F300"])
def test_parser_refuses_other_helical_words(new: str) -> None:
    with pytest.raises(UnsupportedWordError):
        parse_program(HELICAL.replace("I-3.0000 J0.0000 F300", new))


@pytest.mark.parametrize("new", ["I-3.0000 J0.0000 P2 F300", "Z11.6786 J0.0000 I-3.0000 F300"])
def test_parser_refuses_other_helical_structures(new: str) -> None:
    with pytest.raises(MalformedProgramError):
        parse_program(HELICAL.replace("Z11.6786 I-3.0000 J0.0000 F300" if new.startswith("Z") else "I-3.0000 J0.0000 F300", new))
