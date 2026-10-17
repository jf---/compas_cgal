"""The first test cut's setup: every bracketed scope default, each a placeholder pending the user's answer to Qn.

Held apart from `compas_cgal.post.machine` because these values change with
every user answer, while the setup types and their invariants do not.
"""

from fractions import Fraction

from compas_cgal.helix_units import DEFAULT_MAX_RAMP_ANGLE_DEG
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.post.machine import ArcCentreMode
from compas_cgal.post.machine import CutterSpec
from compas_cgal.post.machine import CuttingParameters
from compas_cgal.post.machine import Dialect
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import StockBlank
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post.machine import WorkOffset

# Q1 placeholder: LinuxCNC reads coordinates as doubles, so it does not round a 0.0001 mm grid; the
# machine's real input resolution is pending the user's answer to Q1.
TESTCUT_CONTROLLER_RESOLUTION = Fraction(1, 10000)
# Q3 placeholder: LinuxCNC's default arc-radius tolerance, mm.
TESTCUT_ARC_TOLERANCE = Fraction(2, 1000)
# Q4 placeholder: 6 mm, 3-flute, centre-cutting carbide (the gate's Ø 2 tool at x3).
TESTCUT_CUTTER_DIAMETER = Fraction(6)
TESTCUT_CUTTER_FLUTES = 3
# Q5 placeholder: unitless gate pocket x3 -> 36 x 24 mm rect_12x8.
TESTCUT_SCALE = Fraction(3)
# Q6 placeholder: 6061 aluminium blank 60 x 50 x 15 mm (wax or MDF dry run first).
TESTCUT_STOCK_MATERIAL = "6061 aluminium"
TESTCUT_STOCK_SIZE = (Fraction(60), Fraction(50), Fraction(15))
# Q6/Q9 placeholder: the 36 x 24 mm pocket centred on the 60 x 50 blank, with the work origin at the
# pocket's min corner, puts the blank's min corner at (-(60-36)/2, -(50-24)/2) = (-12, -13) mm.
TESTCUT_STOCK_MIN_CORNER = (Fraction(-12), Fraction(-13))
# Q7 placeholder: one Z level at 0.25 D = 1.5 mm.
TESTCUT_AXIAL_DEPTH = Fraction(3, 2)
# Q8 placeholders: 10000 rpm, 600 mm/min cut, 100 mm/min plunge, clearance 2 D = 12 mm above the stock top.
TESTCUT_SPINDLE_RPM = 10000
TESTCUT_CUT_FEED = 600
TESTCUT_PLUNGE_FEED = 100
# Q8 placeholder: helix entry feed, half the cut feed, mm/min.
TESTCUT_RAMP_FEED = TESTCUT_CUT_FEED // 2
# Ramp-cap placeholder: the plan's D1 default angle (`DEFAULT_MAX_RAMP_ANGLE_DEG`, 2 deg) as its exact
# surrogate K = 2*pi*tan(alpha), pending the cutter datasheet (Q4).
# Q8 placeholder: the helix ramp starts 1 mm above the stock top (the tool feeds down to it in air), mm.
TESTCUT_RAMP_START_ABOVE_TOP = Fraction(1)
TESTCUT_RAMP_CAP = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(DEFAULT_MAX_RAMP_ANGLE_DEG))
# Q8 placeholder: 2 s for the spindle to reach speed before the first motion, s.
TESTCUT_SPINDLE_DWELL = Fraction(2)
TESTCUT_CLEARANCE_ABOVE_TOP = 2 * TESTCUT_CUTTER_DIAMETER
# Output quantum (scope "e.g. 0.001 mm", not a user question): 4 digits = 0.0001 mm. The scope's 0.001 mm
# fails ARC_RADIUS_GAP_QUANTA_SQUARED against Q3's 0.002 mm (8e-6 > 4e-6 mm^2); 0.0001 mm gives 8e-8.
TESTCUT_DECIMALS = 4


def build_testcut_defaults(*, envelope: TravelEnvelope, pocket_min_corner: tuple[Fraction, Fraction], accepted_plunges: int) -> PostSetup:
    """The first test cut's setup: every bracketed scope default, each a placeholder pending the user's answer.

    Q1 LinuxCNC, 0.0001 mm resolution; Q2 G17 + G91.1; Q3 0.002 mm; Q4 6 mm,
    3 flutes, centre-cutting; Q5 x3; Q6 6061, 60 x 50 x 15 mm, pocket centred;
    Q7 1.5 mm; Q8 10000 rpm, 2 s dwell, 600/100/300 mm/min (cut/plunge/ramp),
    clearance 12 mm above the top, helix ramp from 1 mm above the top; ramp cap 2 deg (D1); Q9 G54 with XY at the pocket's min corner. Q10 (*envelope*) has
    no default.

    Args:
        envelope: Q10, the machine travel box in the work frame. Required.
        pocket_min_corner: The pocket's min corner in source units (a fact of
            the pocket, not a default): Q9 places it at the work origin.
        accepted_plunges: The plunges the user accepts (not a default).

    Returns:
        The validated setup.
    """
    return PostSetup.build(
        dialect=Dialect.LINUXCNC,
        controller_resolution=TESTCUT_CONTROLLER_RESOLUTION,
        arc_mode=ArcCentreMode.XY_INCREMENTAL_IJ,
        arc_tolerance=TESTCUT_ARC_TOLERANCE,
        cutter=CutterSpec.build(TESTCUT_CUTTER_DIAMETER, TESTCUT_CUTTER_FLUTES, True),
        scale=TESTCUT_SCALE,
        stock=StockBlank.build(TESTCUT_STOCK_MATERIAL, TESTCUT_STOCK_MIN_CORNER, *TESTCUT_STOCK_SIZE),
        axial_depth=TESTCUT_AXIAL_DEPTH,
        cutting=CuttingParameters.build(TESTCUT_SPINDLE_RPM, TESTCUT_CUT_FEED, TESTCUT_PLUNGE_FEED),
        spindle_dwell=TESTCUT_SPINDLE_DWELL,
        clearance_above_top=TESTCUT_CLEARANCE_ABOVE_TOP,
        work_offset=WorkOffset.G54,
        xy_shift=(-TESTCUT_SCALE * pocket_min_corner[0], -TESTCUT_SCALE * pocket_min_corner[1]),
        envelope=envelope,
        decimals=TESTCUT_DECIMALS,
        accepted_plunges=accepted_plunges,
        ramp_feed=TESTCUT_RAMP_FEED,
        ramp_cap=TESTCUT_RAMP_CAP,
        ramp_start_above_top=TESTCUT_RAMP_START_ABOVE_TOP,
    )
