"""The public atlas helpers of `compas_cgal.incident_loop` that consumers share."""

from __future__ import annotations

import dataclasses
import math
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import strategies as st

from compas_cgal.incident_loop import INCIDENCE_CHART_COUNT
from compas_cgal.incident_loop import INCIDENCE_CHART_DENOMINATOR
from compas_cgal.incident_loop import ChartParameterRangeError
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.incident_loop import LoopEntryIncidenceError
from compas_cgal.incident_loop import chart_phase
from compas_cgal.incident_loop import chart_triple
from compas_cgal.incident_loop import exact_double
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

CHARTS = st.integers(min_value=0, max_value=INCIDENCE_CHART_COUNT - 1)
NUMERATORS = st.integers(min_value=0, max_value=INCIDENCE_CHART_DENOMINATOR)
PHASES = st.floats(min_value=-20.0, max_value=20.0, allow_nan=False, allow_infinity=False)


def test_the_chart_count_is_the_four_quarters() -> None:
    assert INCIDENCE_CHART_COUNT == 4


@given(CHARTS, NUMERATORS)
def test_chart_triple_is_a_reduced_pythagorean_triple(chart: int, numerator: int) -> None:
    a, b, m = chart_triple(chart, numerator)
    assert a * a + b * b == m * m
    assert math.gcd(a, b) == 1
    assert m > 0


@given(CHARTS, NUMERATORS)
def test_chart_phase_is_atan2_of_the_triple(chart: int, numerator: int) -> None:
    a, b, _ = chart_triple(chart, numerator)
    assert chart_phase(chart, numerator) == math.atan2(b, a)


@given(PHASES)
def test_build_stores_the_phase_derived_from_its_chart_parameter(phase: float) -> None:
    loop = IncidentLoop.build(Point2[WorldXY](Millimetre(1.25), Millimetre(-3.5)), GuideRadius(Millimetre(2.0)), Radian(phase))
    assert loop.phase == chart_phase(loop.chart, loop.chart_numerator)
    assert loop.certifies_chart()


@given(PHASES)
def test_build_entry_direction_is_the_chart_triple(phase: float) -> None:
    loop = IncidentLoop.build(Point2[WorldXY](Millimetre(0.0), Millimetre(0.0)), GuideRadius(Millimetre(3.0)), Radian(phase))
    a, b, m = chart_triple(loop.chart, loop.chart_numerator)
    rho = Fraction(loop.radius.value)
    assert (Fraction(loop.entry.x) - Fraction(loop.centre.x)) * m == a * rho
    assert (Fraction(loop.entry.y) - Fraction(loop.centre.y)) * m == b * rho


def test_certifies_chart_refuses_a_contradicting_chart_field() -> None:
    loop = IncidentLoop.build(Point2[WorldXY](Millimetre(0.0), Millimetre(0.0)), GuideRadius(Millimetre(2.0)), Radian(0.7))
    rotated = dataclasses.replace(loop, chart=(loop.chart + 1) % INCIDENCE_CHART_COUNT)
    assert not rotated.certifies_chart()


def test_certifies_chart_refuses_a_phase_not_derived_from_the_chart() -> None:
    loop = IncidentLoop.build(Point2[WorldXY](Millimetre(0.0), Millimetre(0.0)), GuideRadius(Millimetre(2.0)), Radian(0.7))
    nudged = dataclasses.replace(loop, phase=Radian(math.nextafter(float(loop.phase), math.inf)))
    assert not nudged.certifies_chart()


@pytest.mark.parametrize(("chart", "numerator"), [(-1, 0), (INCIDENCE_CHART_COUNT, 0), (0, -1), (0, INCIDENCE_CHART_DENOMINATOR + 1)])
def test_chart_parameter_out_of_range_is_refused_by_name(chart: int, numerator: int) -> None:
    with pytest.raises(ChartParameterRangeError):
        chart_triple(chart, numerator)
    with pytest.raises(ChartParameterRangeError):
        chart_phase(chart, numerator)


def test_exact_double_returns_the_equal_double() -> None:
    assert exact_double(Fraction(3, 2**40)) == 3 * 2.0**-40


def test_exact_double_refuses_a_rational_that_is_not_a_double() -> None:
    with pytest.raises(LoopEntryIncidenceError):
        exact_double(Fraction(1, 3))
