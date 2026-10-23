"""Exact squared distances of double point pairs over one power-of-two denominator."""

from __future__ import annotations

from collections.abc import Iterable
from typing import NamedTuple

from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


class CommonScaleSquaredDistances(NamedTuple):
    """Exact squared distances ``numerators[i] / 4**exponent`` in mm², one exponent for all.

    With one shared denominator, the exact order and equality of the squared
    distances are the order and equality of the integer numerators, so an
    argmin is the builtin ``min`` and a tie is ``==``.
    """

    numerators: tuple[int, ...]
    exponent: int

    def exceeds_squared(self, numerator: int, bound: float) -> bool:
        """Exact verdict ``numerator / 4**exponent > bound**2``.

        ``bound`` is a finite double ``n / d`` (``d`` a power of two), so the
        verdict is the integer comparison ``numerator * d**2 > n**2 * 4**exponent``.

        Args:
            numerator: One of ``numerators`` (or any integer on this scale).
            bound: The distance bound in mm; enters exactly.

        Returns:
            True when the squared distance is strictly larger than ``bound**2``.
        """
        bound_numerator, bound_denominator = bound.as_integer_ratio()
        return numerator * bound_denominator * bound_denominator > (bound_numerator * bound_numerator) << (2 * self.exponent)


def common_scale_squared_distances(
    pairs: Iterable[tuple[Point2[WorldXY], Point2[WorldXY]]],
) -> CommonScaleSquaredDistances:
    """Exact squared distance of every double point pair over one power-of-two denominator.

    Every finite double is ``n / 2**k``; the largest ``k`` over all coordinates
    is a common denominator of every coordinate, so every difference is an
    integer over ``2**k`` and every squared distance an integer over ``4**k``.
    """
    ratios = [
        (
            float(first.x).as_integer_ratio(),
            float(first.y).as_integer_ratio(),
            float(second.x).as_integer_ratio(),
            float(second.y).as_integer_ratio(),
        )
        for first, second in pairs
    ]
    scale = max((denominator for row in ratios for _, denominator in row), default=1)
    numerators = []
    for (first_x, first_x_scale), (first_y, first_y_scale), (second_x, second_x_scale), (second_y, second_y_scale) in ratios:
        delta_x = second_x * (scale // second_x_scale) - first_x * (scale // first_x_scale)
        delta_y = second_y * (scale // second_y_scale) - first_y * (scale // first_y_scale)
        numerators.append(delta_x * delta_x + delta_y * delta_y)
    return CommonScaleSquaredDistances(tuple(numerators), scale.bit_length() - 1)
