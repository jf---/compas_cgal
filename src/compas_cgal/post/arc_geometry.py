"""Exact predicates on emitted arcs: radius agreement and swept directions.

Every quantity is a `Fraction` read from the emitted decimals; nothing is
rounded, and no square root is taken: radii are compared squared, in the exact
form below.
"""

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Vector:
    """An exact XY vector, mm."""

    x: Fraction
    y: Fraction


def cross(a: Vector, b: Vector) -> Fraction:
    """z of ``a x b``: positive iff *b* is counter-clockwise of *a*."""
    return a.x * b.y - a.y * b.x


def dot(a: Vector, b: Vector) -> Fraction:
    """The dot product."""
    return a.x * b.x + a.y * b.y


def radius_within_tolerance(start_squared: Fraction, end_squared: Fraction, tolerance: Fraction) -> bool:
    """Whether ``|r_s - r_e| <= tol``, decided exactly from the squared radii.

    ``(r_s - r_e)**2 <= tol**2`` is ``r_s**2 + r_e**2 - tol**2 <= 2 r_s r_e``.
    The right side is non-negative, so a non-positive left side passes; else
    both sides are non-negative and squaring is exact:
    ``(r_s**2 + r_e**2 - tol**2)**2 <= 4 r_s**2 r_e**2``.

    Args:
        start_squared: ``r_s**2``, mm^2.
        end_squared: ``r_e**2``, mm^2.
        tolerance: ``tol``, mm.

    Returns:
        True iff the radii agree within *tolerance*.
    """
    left = start_squared + end_squared - tolerance * tolerance
    if left <= 0:
        return True
    return left * left <= 4 * start_squared * end_squared


def sweep_contains(start: Vector, end: Vector, direction: Vector, clockwise: bool) -> bool:
    """Whether the closed arc swept from direction *start* to *end* contains *direction*.

    All three are nonzero offsets from the arc centre. A clockwise sweep from
    *start* to *end* is the counter-clockwise sweep from *end* to *start*. Two
    collinear same-sense directions are the near-full or near-zero turn a
    controller resolves its own way: they count as containing every direction
    (the conservative reading for an envelope).

    Args:
        start: Centre-to-start offset.
        end: Centre-to-end offset.
        direction: The direction asked about.
        clockwise: Travel sense.

    Returns:
        True iff *direction* is in the swept closed sector.
    """
    first, last = (end, start) if clockwise else (start, end)
    turn = cross(first, last)
    if turn > 0:
        return cross(first, direction) >= 0 and cross(direction, last) >= 0
    if turn < 0:
        return not (cross(last, direction) > 0 and cross(direction, first) > 0)
    if dot(first, last) > 0:
        return True
    return cross(first, direction) >= 0
