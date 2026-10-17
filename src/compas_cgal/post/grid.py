"""The exact decimal grid of the emitted program: quantise, format, parse.

Quantisation is the post's one declared boundary injection: an exact rational
is rounded half-even onto the 10**-decimals mm grid once, and everything
downstream (text, parse, checks) is exact on that grid.
"""

import re
from fractions import Fraction

from compas_cgal.post.errors import MalformedProgramError
from compas_cgal.post.errors import OffGridValueError
from compas_cgal.post.quantities import DecimalPlaces

# Exactly the decimal forms the writer produces: optional minus, digits, and a
# fractional part of digits. No exponent, no leading '.', no '+'.
_DECIMAL = re.compile(r"-?\d+(\.\d+)?")

# Exactly the integer forms the writer produces for S and F words.
_INTEGER = re.compile(r"\d+")


def quantum(decimals: DecimalPlaces) -> Fraction:
    """The grid step, mm."""
    return Fraction(1, 10**decimals)


def quantise(value: Fraction, decimals: DecimalPlaces) -> Fraction:
    """Round *value* half-even onto the grid.

    Args:
        value: Exact value, mm.
        decimals: Grid digits.

    Returns:
        The nearest grid value (ties to even).
    """
    step = quantum(decimals)
    return round(value / step) * step


def format_decimal(value: Fraction, decimals: DecimalPlaces) -> str:
    """Write a grid value as fixed-point text, exactly.

    Args:
        value: A value on the grid.
        decimals: Grid digits.

    Returns:
        The text, e.g. ``-18.0000``.

    Raises:
        OffGridValueError: *value* is not on the grid.
    """
    scaled = value * 10**decimals
    if scaled.denominator != 1:
        raise OffGridValueError(f"{value} is not a multiple of 10**-{decimals}.")
    digits = str(abs(scaled.numerator)).rjust(decimals + 1, "0")
    sign = "-" if scaled.numerator < 0 else ""
    return f"{sign}{digits[:-decimals]}.{digits[-decimals:]}"


def parse_decimal(text: str) -> Fraction:
    """Read a decimal the writer could have produced, exactly.

    Raises:
        MalformedProgramError: Any other numeric form.
    """
    if _DECIMAL.fullmatch(text) is None:
        raise MalformedProgramError(f"{text!r} is not a plain decimal.")
    return Fraction(text)


def parse_integer(text: str) -> int:
    """Read a non-negative integer the writer could have produced.

    Raises:
        MalformedProgramError: Any other numeric form.
    """
    if _INTEGER.fullmatch(text) is None:
        raise MalformedProgramError(f"{text!r} is not a plain integer.")
    return int(text)
