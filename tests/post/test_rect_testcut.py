"""Round trip writer -> parser -> checks on the real rect_12x8 120-degree pipeline output, scaled x3.

The program is generated here, never read from a committed fixture: a change
to the generator shows up as a change in what the post must accept. The entry
is whatever the generator emits (1.1: a PLUNGE; 1.2: a HELIX); every test
reads it from the stream.
"""

import math
from fractions import Fraction

import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from compas_cgal.helix import Helix
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.post.checks import check_program
from compas_cgal.post.errors import OutOfEnvelopeError
from compas_cgal.post.fail_closed import post_operations
from compas_cgal.post.grid import quantise
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post.mapping import SourcePath
from compas_cgal.post.mapping import map_operations
from compas_cgal.post.parser import parse_program
from compas_cgal.post.program import CentreOffset
from compas_cgal.post.program import Motion
from compas_cgal.post.program import MotionKind
from compas_cgal.post.program import PostedProgram
from compas_cgal.post.program import ProgramPoint
from compas_cgal.post.quantities import ProgramMillimetre
from compas_cgal.post.writer import write_program
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from tests.post.support import mm
from tests.post.support import rebuild

# Arc/Circle operations on this output, every one clockwise (climb): 31 on the 1.2 generator (HELIX
# entry), pinned 2026-10-17 on the landed generate_narrow. 1.1 (PLUNGE entry) had 26
# (testcut-emission-audit.md section 1); that entry is no longer emitted here.
RECT_CIRCULAR_OPERATIONS = 31


@pytest.fixture(scope="module")
def posted(rect_result: PipelineResult, setup: PostSetup, rect_source: SourcePath) -> PostedProgram:
    return map_operations(rect_result.operations, setup, rect_source)


def _xy(setup: PostSetup, x: Fraction, y: Fraction) -> tuple[ProgramMillimetre, ProgramMillimetre]:
    shift_x, shift_y = setup.xy_shift
    return mm(quantise(setup.scale * x + shift_x, setup.decimals)), mm(quantise(setup.scale * y + shift_y, setup.decimals))


def _expected_helix(op: ToolpathOperation, helix: Helix, start: ProgramPoint, setup: PostSetup) -> list[Motion]:
    """G1 to the ramp start, program-unit ramp turns as half-turns with exactly interpolated, once-quantised Z, then one flat closing turn."""
    kind = MotionKind.CLOCKWISE_ARC if op.clockwise else MotionKind.COUNTERCLOCKWISE_ARC
    cx, cy = Fraction(helix.footprint.centre.x), Fraction(helix.footprint.centre.y)
    ex, ey = Fraction(helix.entry_xy.x), Fraction(helix.entry_xy.y)
    qcx, qcy = _xy(setup, cx, cy)
    entry, antipode = _xy(setup, ex, ey), _xy(setup, 2 * cx - ex, 2 * cy - ey)
    top = setup.ramp_start_z
    descent = top - setup.cut_z
    # Oracle of `Helix.build`'s rule in program units: the fewest turns with descent / n <= K * r.
    halves = 2 * math.ceil(descent / (setup.ramp_cap.exact * setup.scale * Fraction(helix.footprint.radius.value)))
    ends = [(antipode if k % 2 else entry, mm(quantise(top - descent * Fraction(k, halves), setup.decimals))) for k in range(1, halves + 1)]
    ends += [(antipode, setup.cut_z), (entry, setup.cut_z)]
    start = ProgramPoint(entry[0], entry[1], top)
    motions = [Motion(MotionKind.LINEAR, start, None, setup.cutting.plunge_feed)]
    for (x, y), z in ends:
        motions.append(Motion(kind, ProgramPoint(x, y, z), CentreOffset(mm(qcx - start.x), mm(qcy - start.y)), setup.ramp_feed))
        start = ProgramPoint(x, y, z)
    return motions


def _expected(op: ToolpathOperation, start: ProgramPoint, setup: PostSetup) -> list[Motion]:
    """The operation quantised once, independently of the mapper: the round trip's reference."""
    geometry = op.geometry
    z = setup.cut_z
    if isinstance(geometry, Helix):
        return _expected_helix(op, geometry, start, setup)
    if isinstance(geometry, Line):
        x, y = _xy(setup, Fraction(float(geometry.end[0])), Fraction(float(geometry.end[1])))
        if op.operation is OperationType.RETRACT:
            return [Motion(MotionKind.RAPID, ProgramPoint(x, y, setup.clearance_z), None, None)]
        feed = setup.cutting.plunge_feed if op.operation is OperationType.PLUNGE else setup.cutting.cut_feed
        return [Motion(MotionKind.LINEAR, ProgramPoint(x, y, z), None, feed)]
    kind = MotionKind.CLOCKWISE_ARC if op.clockwise else MotionKind.COUNTERCLOCKWISE_ARC
    assert isinstance(geometry, (Arc, Circle))
    cx, cy = Fraction(float(geometry.frame.point[0])), Fraction(float(geometry.frame.point[1]))
    qcx, qcy = _xy(setup, cx, cy)
    if isinstance(geometry, Circle):
        assert op.entry_xy is not None
        ex, ey = Fraction(op.entry_xy[0]), Fraction(op.entry_xy[1])
        ends = [_xy(setup, 2 * cx - ex, 2 * cy - ey), _xy(setup, ex, ey)]
    else:
        assert isinstance(geometry, Arc)
        end = geometry.point_at(1.0)
        ends = [_xy(setup, Fraction(float(end[0])), Fraction(float(end[1])))]
    motions = []
    for x, y in ends:
        motions.append(Motion(kind, ProgramPoint(x, y, z), CentreOffset(mm(qcx - start.x), mm(qcy - start.y)), setup.cutting.cut_feed))
        start = ProgramPoint(x, y, z)
    return motions


def test_text_round_trips_exactly(posted: PostedProgram, setup: PostSetup) -> None:
    assert parse_program(write_program(posted.program, setup)) == posted.program


def test_every_operation_matches_its_quantisation(posted: PostedProgram, rect_result: PipelineResult, setup: PostSetup) -> None:
    motions = posted.program.motions
    assert len(posted.spans) == len(rect_result.operations)
    for span, op in zip(posted.spans, rect_result.operations):
        start = motions[span.first_motion - 1].target
        assert list(motions[span.first_motion : span.first_motion + span.motion_count]) == _expected(op, start, setup), span


def test_checks_accept_the_parsed_program(posted: PostedProgram, setup: PostSetup) -> None:
    check_program(parse_program(write_program(posted.program, setup)), setup)


def test_every_clockwise_circular_operation_is_g2(posted: PostedProgram, rect_result: PipelineResult) -> None:
    assert len([op for op in rect_result.operations if isinstance(op.geometry, (Arc, Circle))]) == RECT_CIRCULAR_OPERATIONS
    circular = [index for index, op in enumerate(rect_result.operations) if isinstance(op.geometry, (Arc, Circle, Helix))]
    assert all(rect_result.operations[index].clockwise for index in circular)
    for index in circular:
        span = posted.spans[index]
        kinds = {motion.kind for motion in posted.program.motions[span.first_motion : span.first_motion + span.motion_count]}
        # A HELIX span opens with its Z-only G1 air descent to the ramp start.
        assert kinds - ({MotionKind.LINEAR} if isinstance(rect_result.operations[index].geometry, Helix) else set()) == {MotionKind.CLOCKWISE_ARC}, index
    assert not any(motion.kind is MotionKind.COUNTERCLOCKWISE_ARC for motion in posted.program.motions)


def test_first_rapid_reaches_clearance_above_the_entry(posted: PostedProgram, rect_result: PipelineResult, setup: PostSetup) -> None:
    entry = rect_result.operations[0]
    assert entry.operation is OperationType.HELIX and isinstance(entry.geometry, Helix)
    x, y = _xy(setup, Fraction(entry.geometry.entry_xy.x), Fraction(entry.geometry.entry_xy.y))
    assert posted.program.approach_z == setup.clearance_z
    first = posted.program.motions[0]
    assert first == Motion(MotionKind.RAPID, ProgramPoint(x, y, setup.clearance_z), None, None)
    assert posted.spans[0].first_motion == 1
    assert posted.program.motions[1].kind is MotionKind.LINEAR


def test_fail_closed_post_returns_checked_text(rect_result: PipelineResult, setup: PostSetup, rect_source: SourcePath) -> None:
    text = post_operations(rect_result.operations, setup, rect_source)
    assert text.startswith("G21 G90 G94 G40 G80 G17 G91.1 G54\nS10000 M3\nG4 P2.0000\nG0 Z12.0000\n")
    assert text.endswith("M5\nM30\n")


def test_envelope_violation_raises(rect_result: PipelineResult, setup: PostSetup, rect_source: SourcePath) -> None:
    # The swept pocket reaches x = 35.994 mm in the work frame; a 30 mm travel cannot hold it.
    short = TravelEnvelope.build(x=(Fraction(-10), Fraction(30)), y=(Fraction(-10), Fraction(100)), z=(Fraction(-10), Fraction(50)))
    with pytest.raises(OutOfEnvelopeError):
        post_operations(rect_result.operations, rebuild(setup, envelope=short), rect_source)


def test_helix_entry_rect_posts_with_no_accepted_plunge(rect_result: PipelineResult, setup: PostSetup, rect_source: SourcePath) -> None:
    # The HELIX entry is feed motion, not a plunge, and the rect has no other entry: nothing to accept.
    # The over-allowance refusal itself is pinned in test_checks.test_second_plunge_raises.
    assert rect_result.operations[0].operation is OperationType.HELIX
    assert not any(op.operation is OperationType.PLUNGE for op in rect_result.operations)
    post_operations(rect_result.operations, rebuild(setup, accepted_plunges=0), rect_source)
