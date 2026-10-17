"""The helix turn-count rule as a public pure function (Task 1.2; the post maps program units through it).

``ramp_turns`` is the least ``n >= 1`` with ``descent / n <= exact(K) * rho_h``,
over the rationals; `Helix.build` uses it, so the two cannot drift.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from compas_cgal.helix import Helix
from compas_cgal.helix import ramp_turns
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

DEPTH_MM = Fraction(4)
# Integer-ratio boundary (as pinned for Helix.build): K * rho_h divides the
# descent exactly, so the least count meets the bound with equality.
INTEGER_RATIO_CASES = ((1.0, 4), (0.5, 8))
# Probe P3 at rho_h = 1, depth 4 mm.
P3_TURNS_BY_ANGLE_DEG = {1.0: 37, 2.0: 19, 3.0: 13, 5.0: 8}


@pytest.mark.parametrize(("pitch_per_radius", "turns"), INTEGER_RATIO_CASES)
def test_the_least_count_meets_an_integer_ratio_with_equality(pitch_per_radius: float, turns: int) -> None:
    limit = HelixPitchPerRadius.build(pitch_per_radius)
    radius = GuideRadius.build(1.0)
    assert ramp_turns(ExactMillimetre(DEPTH_MM), radius, limit) == turns
    assert DEPTH_MM / turns == limit.exact * Fraction(radius.value)
    assert DEPTH_MM / (turns - 1) > limit.exact * Fraction(radius.value)


@pytest.mark.parametrize(("angle_deg", "turns"), sorted(P3_TURNS_BY_ANGLE_DEG.items()))
def test_the_p3_table(angle_deg: float, turns: int) -> None:
    limit = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(angle_deg))
    assert ramp_turns(ExactMillimetre(DEPTH_MM), GuideRadius.build(1.0), limit) == turns


def test_a_shallow_descent_is_one_turn() -> None:
    assert ramp_turns(ExactMillimetre(Fraction(1, 1000)), GuideRadius.build(1.0), HelixPitchPerRadius.build(1.0)) == 1


def test_helix_build_uses_the_rule() -> None:
    limit = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(2.0))
    footprint = IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(1.25), Radian(0.4))
    helix = Helix.build(footprint, CutPlane.build(CutZ.build(-1.5), ClearanceZ.build(2.5)), limit, clockwise=True)
    assert helix.turns == ramp_turns(helix.descent, footprint.radius, limit)
