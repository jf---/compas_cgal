"""The exact decimal grid: half-even quantisation, exact formatting and parsing."""

from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import strategies as st

from compas_cgal.post.errors import MalformedProgramError
from compas_cgal.post.errors import OffGridValueError
from compas_cgal.post.grid import format_decimal
from compas_cgal.post.grid import parse_decimal
from compas_cgal.post.grid import quantise
from compas_cgal.post.quantities import DecimalPlaces

FOUR = DecimalPlaces(4)


def test_quantise_rounds_half_even() -> None:
    assert quantise(Fraction(5, 100000), FOUR) == 0
    assert quantise(Fraction(15, 100000), FOUR) == Fraction(2, 10000)
    assert quantise(Fraction(-15, 100000), FOUR) == Fraction(-2, 10000)


def test_format_is_fixed_point() -> None:
    assert format_decimal(Fraction(-18), FOUR) == "-18.0000"
    assert format_decimal(Fraction(-3, 10000), FOUR) == "-0.0003"
    assert format_decimal(Fraction(0), FOUR) == "0.0000"


def test_off_grid_value_raises() -> None:
    with pytest.raises(OffGridValueError):
        format_decimal(Fraction(1, 3), FOUR)


@pytest.mark.parametrize("text", ["1e3", "1/2", ".5", "nan", "1.", "+-1"])
def test_parse_rejects_what_we_never_emit(text: str) -> None:
    with pytest.raises(MalformedProgramError):
        parse_decimal(text)


@given(st.fractions())
def test_format_parse_round_trip(value: Fraction) -> None:
    on_grid = quantise(value, FOUR)
    assert parse_decimal(format_decimal(on_grid, FOUR)) == on_grid
