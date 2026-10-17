"""The helical-entry units are frozen dataclasses with `build` (Task 1.2 amendment M1)."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest

from compas_cgal.helix_units import DEFAULT_MAX_RAMP_ANGLE_DEG
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import InvalidHelixPitchError
from compas_cgal.helix_units import InvalidRampAngleError
from compas_cgal.helix_units import RampAngleDegrees


@pytest.mark.parametrize("value", [0.0, -0.0, -1.0, math.nan, math.inf, -math.inf, True])
def test_helix_pitch_per_radius_refuses_nonpositive_and_nonfinite(value: float) -> None:
    with pytest.raises(InvalidHelixPitchError):
        HelixPitchPerRadius.build(value)
    with pytest.raises(InvalidHelixPitchError):
        HelixPitchPerRadius(value)


def test_the_pitch_surrogate_is_the_double_injected_exactly() -> None:
    angle = RampAngleDegrees.build(2.0)
    surrogate = HelixPitchPerRadius.from_ramp_angle(angle)
    double = 2.0 * math.pi * math.tan(math.radians(2.0))
    assert surrogate.value == double
    assert surrogate.exact == Fraction(double)


@pytest.mark.parametrize("value", [0.0, -2.0, 90.0, 91.0, math.nan, math.inf, False])
def test_a_ramp_angle_outside_the_open_quarter_turn_is_refused(value: float) -> None:
    with pytest.raises(InvalidRampAngleError):
        RampAngleDegrees.build(value)


def test_the_default_ramp_angle_is_the_d1_placeholder() -> None:
    assert DEFAULT_MAX_RAMP_ANGLE_DEG == 2.0
    assert RampAngleDegrees.build(DEFAULT_MAX_RAMP_ANGLE_DEG).value == 2.0
