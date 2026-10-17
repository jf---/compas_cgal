"""Inject one defect per check into a small emitted program and assert the named raise."""

from fractions import Fraction

import pytest

from compas_cgal.post.checks import check_program
from compas_cgal.post.errors import ArcRadiusMismatchError
from compas_cgal.post.errors import DegenerateArcError
from compas_cgal.post.errors import OffPlaneFeedError
from compas_cgal.post.errors import OutOfEnvelopeError
from compas_cgal.post.errors import PocketOutsideStockError
from compas_cgal.post.errors import UnacceptedPlungeError
from compas_cgal.post.errors import UnsafeRapidError
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import StockBlank
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post.parser import parse_program
from tests.post.support import rebuild

# Clearance 12, cut -1.5 (the test-cut defaults). A 4 mm circle about (10, 10) in two CW halves.
GOOD = """G21 G90 G94 G40 G80 G17 G91.1 G54
S10000 M3
G4 P2.0000
G0 Z12.0000
G0 X14.0000 Y10.0000 Z12.0000
G1 X14.0000 Y10.0000 Z-1.5000 F100
G2 X6.0000 Y10.0000 Z-1.5000 I-4.0000 J0.0000 F600
G2 X14.0000 Y10.0000 Z-1.5000 I4.0000 J0.0000 F600
G1 X16.0000 Y10.0000 Z-1.5000 F600
G0 X16.0000 Y10.0000 Z12.0000
M5
M30
"""


def test_good_program_passes(setup: PostSetup) -> None:
    check_program(parse_program(GOOD), setup)


def test_retract_exemption_is_vertical_only(setup: PostSetup) -> None:
    # GOOD's final G0 starts on the cut plane and passes as a vertical retract; the same G0 with an XY step does not.
    text = GOOD.replace("G0 X16.0000 Y10.0000 Z12.0000", "G0 X17.0000 Y10.0000 Z12.0000")
    with pytest.raises(UnsafeRapidError):
        check_program(parse_program(text), setup)


def test_rapid_across_below_clearance_raises(setup: PostSetup) -> None:
    text = GOOD.replace("G1 X16.0000 Y10.0000 Z-1.5000 F600", "G0 X16.0000 Y10.0000 Z-1.5000")
    with pytest.raises(UnsafeRapidError):
        check_program(parse_program(text), setup)


def test_approach_below_clearance_raises(setup: PostSetup) -> None:
    with pytest.raises(UnsafeRapidError):
        check_program(parse_program(GOOD.replace("G0 Z12.0000", "G0 Z5.0000")), setup)


def test_feed_above_cut_plane_raises(setup: PostSetup) -> None:
    text = GOOD.replace("G1 X16.0000 Y10.0000 Z-1.5000 F600", "G1 X16.0000 Y10.0000 Z-1.0000 F600")
    with pytest.raises(OffPlaneFeedError):
        check_program(parse_program(text), setup)


def test_perturbed_arc_raises(setup: PostSetup) -> None:
    text = GOOD.replace("G2 X6.0000 Y10.0000", "G2 X5.9900 Y10.0000")
    with pytest.raises(ArcRadiusMismatchError):
        check_program(parse_program(text), setup)


def test_rounding_level_arc_passes(setup: PostSetup) -> None:
    text = GOOD.replace("G2 X6.0000 Y10.0000", "G2 X5.9999 Y10.0001")
    check_program(parse_program(text), setup)


def test_closed_arc_raises(setup: PostSetup) -> None:
    text = GOOD.replace("G2 X6.0000 Y10.0000 Z-1.5000 I-4.0000", "G2 X14.0000 Y10.0000 Z-1.5000 I-4.0000")
    with pytest.raises(DegenerateArcError):
        check_program(parse_program(text), setup)


def test_second_plunge_raises(setup: PostSetup) -> None:
    text = GOOD.replace(
        "G0 X16.0000 Y10.0000 Z12.0000\n",
        "G0 X16.0000 Y10.0000 Z12.0000\nG0 X20.0000 Y10.0000 Z12.0000\nG1 X20.0000 Y10.0000 Z-1.5000 F100\nG0 X20.0000 Y10.0000 Z12.0000\n",
    )
    check_program(parse_program(text), rebuild(setup, accepted_plunges=2))
    with pytest.raises(UnacceptedPlungeError):
        check_program(parse_program(text), setup)


def test_target_outside_envelope_raises(setup: PostSetup) -> None:
    box = TravelEnvelope.build(x=(Fraction(0), Fraction(15)), y=(Fraction(0), Fraction(20)), z=(Fraction(-5), Fraction(20)))
    with pytest.raises(OutOfEnvelopeError):
        check_program(parse_program(GOOD), rebuild(setup, envelope=box))


def test_arc_bulge_outside_envelope_raises(setup: PostSetup) -> None:
    # Every emitted point has y = 10; the second CW half, (6,10) to (14,10), bulges to y = 14.
    box = TravelEnvelope.build(x=(Fraction(0), Fraction(20)), y=(Fraction(5), Fraction(13)), z=(Fraction(-5), Fraction(20)))
    with pytest.raises(OutOfEnvelopeError):
        check_program(parse_program(GOOD), rebuild(setup, envelope=box))


def test_cut_outside_the_blank_raises(setup: PostSetup) -> None:
    # The line to x = 16 with a 3 mm tool radius reaches exactly x = 19, the narrow blank's closed max side.
    narrow = StockBlank.build(material="wax", min_corner=(Fraction(-12), Fraction(-13)), size_x=Fraction(31), size_y=Fraction(50), size_z=Fraction(15))
    check_program(parse_program(GOOD), rebuild(setup, stock=narrow))
    # The second half-circle bulges to y = 14, so the tool reaches y = 17 > 13.
    tight = StockBlank.build(material="wax", min_corner=(Fraction(-12), Fraction(-13)), size_x=Fraction(31), size_y=Fraction(26), size_z=Fraction(15))
    with pytest.raises(PocketOutsideStockError):
        check_program(parse_program(GOOD), rebuild(setup, stock=tight))
