"""The side probe's step h (side_probe_2.h): (w, w + h·dir] must meet no support, in both filter modes.

A support crossed at exactly τ = 1 has a certified enclosure whose lower end is
1.0: the LAZY pre-pass must still halve, or the probe lands on that support.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from compas_cgal import _coverage_certificate_2 as native


def _wire(value: Fraction | int) -> tuple[str, str]:
    value = Fraction(value)
    return str(value.numerator), str(value.denominator)


def _point(x: Fraction | int, y: Fraction | int) -> tuple[tuple[str, str], tuple[str, str]]:
    return _wire(x), _wire(y)


def _step(circles: list[tuple[int, int, int]], exact_only: bool) -> Fraction:
    numerator, denominator = native.contract_probe_probe_step(_point(0, 0), _point(1, 0), [tuple(_wire(v) for v in c) for c in circles], exact_only)
    return Fraction(int(numerator), int(denominator))


@pytest.mark.parametrize("exact_only", [False, True])
def test_a_support_crossed_at_exactly_one_halves_the_step(exact_only: bool) -> None:
    # The circle about (2, 0) of radius 1 meets the ray w + τ (1, 0) at τ = 1 and τ = 3.
    assert _step([(2, 0, 1)], exact_only) == Fraction(1, 2)


@pytest.mark.parametrize("exact_only", [False, True])
def test_a_support_beyond_one_keeps_the_full_step(exact_only: bool) -> None:
    # About (3, 0), radius 1: τ = 2 and 4, so h = 1 already clears it.
    assert _step([(3, 0, 1)], exact_only) == Fraction(1)


@pytest.mark.parametrize("exact_only", [False, True])
def test_a_support_just_below_one_takes_a_power_of_two_below_it(exact_only: bool) -> None:
    # About (7/4, 0), radius 1: τ = 3/4. The exact path halves to 1/2; the LAZY
    # pre-pass may stop lower (any power of two below the certified lower end).
    numerator, denominator = native.contract_probe_probe_step(_point(0, 0), _point(1, 0), [(_wire(Fraction(7, 4)), _wire(0), _wire(1))], exact_only)
    step = Fraction(int(numerator), int(denominator))
    assert 0 < step < Fraction(3, 4)
    assert step.numerator == 1 and step.denominator & (step.denominator - 1) == 0


@pytest.mark.parametrize("exact_only", [False, True])
def test_a_crossing_just_beyond_w_is_decided_exactly(exact_only: bool) -> None:
    # About (1 + 2^-60, 0), radius 1: τ = 2^-60 and 2 + 2^-60. No binary64
    # enclosure separates the small root from 0, so LAZY must take it exactly;
    # the step halves to 2^-61.
    centre = 1 + Fraction(1, 2**60)
    numerator, denominator = native.contract_probe_probe_step(_point(0, 0), _point(1, 0), [(_wire(centre), _wire(0), _wire(1))], exact_only)
    assert Fraction(int(numerator), int(denominator)) == Fraction(1, 2**61)
