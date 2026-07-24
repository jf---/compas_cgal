import math
import struct
from fractions import Fraction
from typing import assert_type

import numpy as np
import pytest

from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import DegenerateCircleMotionError
from compas_cgal.adaptive.errors import DegenerateSegmentMotionError
from compas_cgal.adaptive.errors import InvalidEngagementCapError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


def test_exact_segment_rejects_zero_progress() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)

    with pytest.raises(DegenerateSegmentMotionError):
        ExactSegmentMotion.build(point, point)


def test_exact_segment_retains_binary64_endpoints() -> None:
    start = Point2[WorldXY].build(0.1, -0.2)
    end = Point2[WorldXY].build(2.5, 3.75)

    motion = ExactSegmentMotion.build(start, end)

    assert motion.start is start
    assert motion.end is end


def test_exact_circle_derives_exact_squared_radius_from_phase() -> None:
    center = Point2[WorldXY].build(2.0, 3.0)
    phase = Vector2[WorldXY].build(3.0, 4.0)

    motion = ExactCircleMotion.build(center, phase, clockwise=True)

    assert motion.squared_radius == Fraction(25, 1)
    assert_type(motion.squared_radius, SquaredMillimetre)
    assert motion.clockwise is True


def test_exact_circle_uses_binary64_rationals_not_rounded_decimal_math() -> None:
    phase = Vector2[WorldXY].build(0.1, 0.2)

    motion = ExactCircleMotion.build(Point2[WorldXY].build(0.0, 0.0), phase, clockwise=False)

    expected = Fraction.from_float(0.1) ** 2 + Fraction.from_float(0.2) ** 2
    assert motion.squared_radius == expected


def test_exact_circle_rejects_zero_phase() -> None:
    with pytest.raises(DegenerateCircleMotionError):
        ExactCircleMotion.build(
            Point2[WorldXY].build(0.0, 0.0),
            Vector2[WorldXY].build(0.0, 0.0),
            clockwise=False,
        )


def test_raw_motion_constructors_cannot_bypass_invariants() -> None:
    point = Point2[WorldXY].build(0.0, 0.0)
    zero = Vector2[WorldXY].build(0.0, 0.0)

    with pytest.raises(DegenerateSegmentMotionError):
        ExactSegmentMotion(point, point)
    with pytest.raises(DegenerateCircleMotionError):
        ExactCircleMotion(point, zero, False)
    with pytest.raises(DegenerateCircleMotionError):
        ExactCircleMotion.build(point, Vector2[WorldXY].build(1.0, 0.0), clockwise=1)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "theta",
    [
        0.0,
        -0.0,
        -1.0,
        math.inf,
        math.nan,
        math.nextafter(math.pi, math.inf),
        float.fromhex("0x0.0000000000001p-1022"),
        True,
    ],
)
def test_engagement_cap_rejects_angles_outside_native_domain(theta: float) -> None:
    with pytest.raises(InvalidEngagementCapError):
        EngagementCap.build(theta)


def test_engagement_cap_owns_native_surrogate_bytes() -> None:
    theta = math.radians(87.0)

    cap = EngagementCap.build(theta)
    native_ratio = _stock_2.cap_chord_ratio(theta)

    assert cap.chord_ratio == native_ratio
    assert cap.chord_ratio_bytes == struct.pack(">d", native_ratio)


def test_pi_cap_has_exact_full_chord_surrogate() -> None:
    assert _stock_2.cap_chord_ratio(cap_radians=math.pi) == 4.0


def test_native_cap_comparison_orders_adjacent_binary64_surrogates() -> None:
    upper_ratio = _stock_2.cap_chord_ratio(math.pi / 2.0)
    lower_ratio = math.nextafter(upper_ratio, 0.0)

    assert _stock_2.cap_chord_ratio_le(lhs=lower_ratio, rhs=upper_ratio)
    assert _stock_2.cap_chord_ratio_le(lhs=upper_ratio, rhs=upper_ratio)
    assert not _stock_2.cap_chord_ratio_le(lhs=upper_ratio, rhs=lower_ratio)


@pytest.mark.parametrize(
    "theta",
    [
        0.0,
        -0.0,
        -1.0,
        math.inf,
        math.nan,
        math.nextafter(math.pi, math.inf),
        float.fromhex("0x0.0000000000001p-1022"),
    ],
)
def test_native_cap_conversion_rejects_unrepresentable_domain(theta: float) -> None:
    with pytest.raises(ValueError):
        _stock_2.cap_chord_ratio(cap_radians=theta)


def test_segment_certifier_reuses_native_cap_conversion() -> None:
    square = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    stock = _stock_2.Stock2(square, [])
    underflowing_cap = float.fromhex("0x0.0000000000001p-1022")

    with pytest.raises(ValueError, match="representable chord ratio"):
        _stock_2.certify_segment_tea(stock, 0.5, 0.5, 1.5, 0.5, 0.1, underflowing_cap)


@pytest.mark.parametrize("ratio", [0.0, -0.0, -1.0, math.inf, math.nan, math.nextafter(4.0, math.inf)])
def test_native_cap_comparison_rejects_values_outside_surrogate_domain(ratio: float) -> None:
    with pytest.raises(ValueError):
        _stock_2.cap_chord_ratio_le(ratio, 1.0)
    with pytest.raises(ValueError):
        _stock_2.cap_chord_ratio_le(1.0, ratio)
