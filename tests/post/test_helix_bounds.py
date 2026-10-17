"""Boundaries of the HELIX checks on hand-written emitted programs: air descent, descent bound, ramp cap, setup gates.

One ramp turn about (18, 12) on r = 3 mm: entry P = (21, 12), antipode
A = (15, 12), two helical half-turns from the ramp start to the cut plane
(-1.5), then the flat closing turn. The ramp cap is set to exact doubles
(``K = 1/2``, ``K = 1``) so the boundaries fall on the 0.0001 mm grid.
"""

from fractions import Fraction

import pytest

from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.post.checks import check_program
from compas_cgal.post.errors import HelixDescentMismatchError
from compas_cgal.post.errors import HelixRampAngleError
from compas_cgal.post.errors import NonCentreCuttingHelixError
from compas_cgal.post.errors import OffGridValueError
from compas_cgal.post.errors import OffPlaneFeedError
from compas_cgal.post.errors import UnacceptedPlungeError
from compas_cgal.post.machine import CutterSpec
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.mapping import map_operations
from compas_cgal.post.parser import parse_program
from tests.post.support import rebuild
from tests.post.test_checks import GOOD
from tests.post.test_helix import SOURCE
from tests.post.test_helix import _helix
from tests.post.test_helix import _stream

# K = 1/2 on r = 3 mm admits exactly 1.5 mm of descent per turn.
HALF = HelixPitchPerRadius.build(0.5)
# K = 1 on r = 3 mm admits 3 mm per turn: wide enough that only the descent bound decides.
ONE = HelixPitchPerRadius.build(1.0)

PLUNGE = "G1 X14.0000 Y10.0000 Z-1.5000 F100\n"


def _turn(top: str, mid: str, antipode_x: str = "15.0000", back_i: str = "3.0000") -> str:
    """One ramp turn from *top* through *mid* to the cut plane, plus the flat closing turn."""
    return f"""G21 G90 G94 G40 G80 G17 G91.1 G54
S10000 M3
G4 P2.0000
G0 Z12.0000
G0 X21.0000 Y12.0000 Z12.0000
G1 X21.0000 Y12.0000 Z{top} F100
G2 X{antipode_x} Y12.0000 Z{mid} I-3.0000 J0.0000 F300
G2 X21.0000 Y12.0000 Z-1.5000 I{back_i} J0.0000 F300
G2 X15.0000 Y12.0000 Z-1.5000 I-3.0000 J0.0000 F300
G2 X21.0000 Y12.0000 Z-1.5000 I3.0000 J0.0000 F300
G0 X21.0000 Y12.0000 Z12.0000
M5
M30
"""


# --- the air-descent rule: a Z-only G1 ends on the cut plane (a plunge) or at/above the stock top Z0 ---


@pytest.mark.parametrize("depth", ["-1.0000", "-0.0001"])
def test_air_descent_into_the_stock_is_refused(setup: PostSetup, depth: str) -> None:
    text = GOOD.replace(PLUNGE, f"G1 X14.0000 Y10.0000 Z{depth} F100\n{PLUNGE}")
    with pytest.raises(OffPlaneFeedError):
        check_program(parse_program(text), setup)


def test_air_descent_to_the_stock_top_is_not_a_plunge(setup: PostSetup) -> None:
    # One accepted plunge: if the descent to Z0 counted, two plunges would raise.
    text = GOOD.replace(PLUNGE, f"G1 X14.0000 Y10.0000 Z0.0000 F100\n{PLUNGE}")
    check_program(parse_program(text), setup)
    with pytest.raises(UnacceptedPlungeError):
        check_program(parse_program(text), rebuild(setup, accepted_plunges=0))


# --- the descent bound |d_k - delta| <= q ---


def test_descent_exactly_one_quantum_off_uniform_passes(setup: PostSetup) -> None:
    # delta = -1.25; the uniform midpoint is -0.25, one quantum above it gives |d - delta| = q.
    check_program(parse_program(_turn("1.0000", "-0.2499")), rebuild(setup, ramp_cap=ONE))


def test_descent_two_quanta_off_uniform_is_refused(setup: PostSetup) -> None:
    with pytest.raises(HelixDescentMismatchError):
        check_program(parse_program(_turn("1.0000", "-0.2498")), rebuild(setup, ramp_cap=ONE))


# --- the ramp cap D^2 <= K^2 r^2 ---


def test_ramp_exactly_at_the_cap_passes(setup: PostSetup) -> None:
    check_program(parse_program(_turn("0.0000", "-0.7500")), rebuild(setup, ramp_cap=HALF))


def test_ramp_one_quantum_over_the_cap_is_refused(setup: PostSetup) -> None:
    with pytest.raises(HelixRampAngleError):
        check_program(parse_program(_turn("0.0001", "-0.7500")), rebuild(setup, ramp_cap=HALF))


def test_ramp_is_certified_on_the_smaller_radius(setup: PostSetup) -> None:
    # Antipode at r = 3.001: D = 1.5003 exceeds K * 3 = 1.5 but not K * 3.001 = 1.5005.
    text = _turn("0.0003", "-0.7498", antipode_x="14.9990", back_i="3.0010")
    with pytest.raises(HelixRampAngleError):
        check_program(parse_program(text), rebuild(setup, ramp_cap=HALF))


# --- setup gates ---


def test_off_grid_ramp_start_is_refused(setup: PostSetup) -> None:
    with pytest.raises(OffGridValueError):
        rebuild(setup, ramp_start_above_top=Fraction(1, 10**5))


def test_helix_entry_needs_a_centre_cutting_cutter(setup: PostSetup) -> None:
    side_only = rebuild(setup, accepted_plunges=0, cutter=CutterSpec.build(Fraction(6), 3, False))
    with pytest.raises(NonCentreCuttingHelixError):
        map_operations(_stream(_helix(True, 0.0)), side_only, SOURCE)
