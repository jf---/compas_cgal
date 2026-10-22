"""Pins on the tier-1b prefilter's stated contract: it keeps every circle within ``R + H + E`` of the piece midpoint.

``_tier_1b_neighbours`` documents its guarantee as: a circle is kept unless
``|c - M| > R + H + E`` provably, with ``E = u (6 (R + H) + 8 P + 2 C)``. The
error term ``E`` is what the rounding argument needs to cover the candidate
filter's binary64 sub-piece midpoints; the per-op outward nudges alone cover
only rounding at the scale of ``R + H``. These tests check the contract
directly, in exact ``Fraction`` arithmetic: circles are placed in the band
``R + H < |c - M| <= R + H + E_low`` with
``E_low = u (6 R + 8 P + 2 C) <= E`` (dropping the ``6 H`` term), where the
coordinates are large enough that ``E_low`` exceeds the outward nudges by
orders of magnitude. A prefilter that drops ``E`` drops these circles.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_engagement_refinement import _tier_1b_neighbours

# Unit roundoff of binary64 (round to nearest), the u of the docstring's E.
UNIT_ROUNDOFF = Fraction(1, 2**53)
# Tool radius of the fixtures (mm); reach = rho + 2 r as the candidate filter computes it.
TOOL_RADIUS = 1.0
Point = tuple[float, float]


def _reach(rho: float) -> float:
    return rho + 2 * TOOL_RADIUS


def _error_floor(reach: float, a: Point, b: Point, centre: Point) -> Fraction:
    """E_low = u (6 R + 8 P + 2 C), a lower bound on the docstring's E (H >= 0)."""
    extent = max(abs(Fraction(v)) for v in (*a, *b))
    spread = max(abs(Fraction(v)) for v in centre)
    return UNIT_ROUNDOFF * (6 * Fraction(reach) + 8 * extent + 2 * spread)


def _within_contract(centre: Point, reach: float, a: Point, b: Point) -> bool:
    """Exact: |c - M| <= R + H + E_low, H = |B - A| / 2 (irrational, compared by squaring)."""
    mx, my = (Fraction(a[0]) + Fraction(b[0])) / 2, (Fraction(a[1]) + Fraction(b[1])) / 2
    d_sq = (Fraction(centre[0]) - mx) ** 2 + (Fraction(centre[1]) - my) ** 2
    q = (Fraction(b[0]) - Fraction(a[0])) ** 2 + (Fraction(b[1]) - Fraction(a[1])) ** 2
    s = Fraction(reach) + _error_floor(reach, a, b, centre)
    # d <= s + sqrt(q)/2  <=>  d² - s² - q/4 <= s sqrt(q)
    t = d_sq - s * s - q / 4
    return t <= 0 or t * t <= s * s * q


def _beyond_r_plus_h(centre: Point, reach: float, a: Point, b: Point) -> bool:
    """Exact: |c - M| > R + H (the band the error term alone must cover)."""
    mx, my = (Fraction(a[0]) + Fraction(b[0])) / 2, (Fraction(a[1]) + Fraction(b[1])) / 2
    d_sq = (Fraction(centre[0]) - mx) ** 2 + (Fraction(centre[1]) - my) ** 2
    q = (Fraction(b[0]) - Fraction(a[0])) ** 2 + (Fraction(b[1]) - Fraction(a[1])) ** 2
    s = Fraction(reach)
    t = d_sq - s * s - q / 4
    return t > 0 and t * t > s * s * q


def _kept(centre: Point, reach: float, a: Point, b: Point) -> bool:
    centres = np.array([centre], dtype=np.float64)
    return bool(_tier_1b_neighbours(centres, np.array([reach], dtype=np.float64), [a, b])[0])


def _place(a: Point, b: Point, reach: float, angle: float, fraction: float) -> Point:
    """A centre at about R + H + fraction * E_low from M in direction ``angle`` (binary64; checked exactly by the caller)."""
    mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
    half = math.dist(a, b) / 2
    guess = (mx + (reach + half) * math.cos(angle), my + (reach + half) * math.sin(angle))
    band = float(_error_floor(reach, a, b, guess))
    distance = reach + half + fraction * band
    return mx + distance * math.cos(angle), my + distance * math.sin(angle)


def test_a_circle_in_the_error_band_on_the_extension_is_kept() -> None:
    """Piece at x ≈ 4096 along the x axis; the circle on its extension at R + H + 3/4 E_low."""
    a, b = (4096.0, 0.0), (4100.0, 0.0)
    reach = _reach(1.5)
    centre = _place(a, b, reach, 0.0, 0.75)
    assert _beyond_r_plus_h(centre, reach, a, b)
    assert _within_contract(centre, reach, a, b)
    assert _kept(centre, reach, a, b)


_MAGNITUDE = st.sampled_from([64.0, 1024.0, 65536.0, 1048576.0])
_UNIT = st.floats(min_value=-1.0, max_value=1.0, allow_nan=False)


@settings(max_examples=400, deadline=None)
@given(
    _MAGNITUDE,
    _UNIT,
    _UNIT,
    st.floats(min_value=1e-3, max_value=8.0),
    st.floats(min_value=0.0, max_value=2 * math.pi),
    st.floats(min_value=0.25, max_value=6.0),
    st.floats(min_value=0.0, max_value=2 * math.pi),
    st.floats(min_value=0.1, max_value=0.9),
)
def test_every_circle_within_r_plus_h_plus_e_is_kept(magnitude: float, ux: float, uy: float, length: float, heading: float, rho: float, angle: float, fraction: float) -> None:
    """Contract: |c - M| <= R + H + E (exact) implies kept; circles drawn in the band beyond R + H."""
    a = (magnitude * ux, magnitude * uy)
    b = (a[0] + length * math.cos(heading), a[1] + length * math.sin(heading))
    assume(a != b)
    reach = _reach(rho)
    centre = _place(a, b, reach, angle, fraction)
    assume(_within_contract(centre, reach, a, b))
    assert _kept(centre, reach, a, b)
