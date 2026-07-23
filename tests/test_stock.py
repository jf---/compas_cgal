import math
import random
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal import _stock_2

getcontext().prec = 80

SQUARE = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=np.float64)
ISLAND = np.array([[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]], dtype=np.float64)


def cap_ratio(cap):
    """Caller-side exact squared-chord surrogate for an angular cap: 4*sin^2(cap/2)."""
    return 4.0 * math.sin(cap / 2.0) ** 2


def test_stock_init_square():
    stock = _stock_2.Stock2(SQUARE, [])
    assert stock.contains(5.0, 5.0)
    assert not stock.contains(-1.0, 5.0)
    assert not stock.contains(11.0, 5.0)
    assert not stock.is_empty()


def test_stock_init_with_island():
    stock = _stock_2.Stock2(SQUARE, [ISLAND])
    assert stock.contains(2.0, 2.0)
    assert not stock.contains(5.0, 5.0)  # inside island = void


def test_stock_rejects_non_simple():
    bowtie = np.array([[0, 0, 0], [6, 0, 0], [0, 4, 0], [6, 4, 0]], dtype=np.float64)
    with pytest.raises(ValueError, match="simple"):
        _stock_2.Stock2(bowtie, [])


def test_subtract_disk_clears_point():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not stock.contains(5.0, 5.0)
    assert not stock.contains(5.9, 5.0)  # inside disk
    assert stock.contains(6.1, 5.0)  # outside disk
    assert stock.contains(2.0, 2.0)


def test_subtract_capsule_geometry():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(2.0, 5.0, 8.0, 5.0, 0.5)
    assert not stock.contains(5.0, 5.0)  # on the axis
    assert not stock.contains(5.0, 5.4)  # inside half-width
    assert stock.contains(5.0, 5.6)  # outside half-width
    assert not stock.contains(8.4, 5.0)  # inside end cap
    assert stock.contains(8.6, 5.0)  # outside end cap


def test_overlapping_subtractions_are_regularized():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(2.0, 5.0, 8.0, 5.0, 0.5)
    stock.subtract_capsule(2.0, 5.2, 8.0, 5.2, 0.5)  # heavy overlap
    stock.subtract_disk(5.0, 5.0, 0.5)  # fully inside cleared
    assert not stock.contains(5.0, 5.3)
    assert stock.contains(5.0, 6.0)


def test_subtract_degenerate_capsule_is_disk():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(5.0, 5.0, 5.0, 5.0, 1.0)
    assert not stock.contains(5.5, 5.0)
    assert stock.contains(6.5, 5.0)


def test_subtract_arc_sweep_full_circle():
    stock = _stock_2.Stock2(SQUARE, [])
    # tool r=0.4 swept around full circle of radius 2 centered at (5,5)
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, 0.4)
    assert not stock.contains(7.0, 5.0)  # on guide circle
    assert not stock.contains(5.0, 7.0)  # opposite side of circle
    assert not stock.contains(7.3, 5.0)  # within tool band (outer)
    assert not stock.contains(6.7, 5.0)  # within tool band (inner)
    assert stock.contains(5.0, 5.0)  # annulus center survives
    assert stock.contains(7.6, 5.0)  # outside band + slack


def test_subtract_arc_sweep_quarter_arc():
    stock = _stock_2.Stock2(SQUARE, [])
    # CCW quarter arc from (7,5) to (5,7) about (5,5), tool r=0.4
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 5.0, 7.0, False, 0.4)
    mid = (5.0 + 2.0 * math.cos(math.pi / 4), 5.0 + 2.0 * math.sin(math.pi / 4))
    assert not stock.contains(*mid)
    assert stock.contains(5.0, 3.0)  # untouched opposite quadrant
    assert stock.contains(3.0, 5.0)


def test_arc_sweep_under_covers_never_over():
    """Conservatism: nothing farther than tool_radius from the guide arc is removed."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, 0.4)
    rng = np.random.default_rng(3)
    for _ in range(400):
        x, y = rng.uniform(0.2, 9.8, 2)
        d_guide = abs(math.hypot(x - 5.0, y - 5.0) - 2.0)
        if d_guide > 0.4:  # strictly outside the true sweep
            assert stock.contains(x, y), (x, y, d_guide)


# --------------------------------------------------------------------------- #
# engagement_at: exact station TEA                                            #
# --------------------------------------------------------------------------- #


def test_engagement_full_material():
    stock = _stock_2.Stock2(SQUARE, [])
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(2.0 * math.pi, abs=1e-9)
    assert exceeded  # 2*pi run > pi cap, decided exactly


def test_engagement_zero_in_cleared():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 2.0)
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(0.0, abs=1e-9)
    assert not exceeded


def test_engagement_half_plane():
    """Cutter centered on the boundary of a cleared half: TEA ~= pi."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(0.0, 7.0, 10.0, 7.0, 2.0)  # clears the band y in [5, 9]
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(math.pi, abs=0.02)  # chain slack << tolerance
    assert max_run == pytest.approx(math.pi, abs=0.02)


def test_engagement_cap_exceeded_over_half():
    """Run just over pi -> exact ORIENTATION branch (theta > pi): any cap <= pi exceeded.

    The cutter sits on a cleared-half boundary; chain under-coverage leaves the
    full lower semicircle plus slivers engaged, so the run measures ~pi + 4e-4,
    strictly above pi. No cap in the contractual (0, pi] range can fail to be
    exceeded here -- the >pi verdict is topological, from sign of orientation,
    never from the reported double.
    """
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(0.0, 7.0, 10.0, 7.0, 2.0)  # clears the band y in [5, 9]
    for cap_deg in (100.0, 179.0, 180.0):
        _, _, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.radians(cap_deg)))
        assert exceeded, cap_deg


def test_engagement_cap_chord_predicate_both_ways():
    """Sub-pi run -> exact CHORD branch: the same test flips on the threshold.

    Clearing the band y in [4.7, 8.7] engages only the cutter's bottom ~106 deg
    cap (theta < pi), so the squared-chord predicate decides. A cap below the
    run is exceeded, one above it is not -- exercising |pq|^2 > T in both signs.
    """
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(0.0, 6.7, 10.0, 6.7, 2.0)  # clears the band y in [4.7, 8.7]
    total, max_run, _ = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(math.radians(106.26), abs=0.02)
    assert max_run == pytest.approx(math.radians(106.26), abs=0.02)
    _, _, exceeded_tight = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.radians(90.0)))
    assert exceeded_tight
    _, _, not_exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.radians(150.0)))
    assert not not_exceeded


def test_engagement_seam_crossing_run():
    """A run straddling the +x seam is split at (cx+r, cy) and re-joined by the
    exact wrap-around merge; its two arc spans sum across the seam to ~pi."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(3.0, 0.0, 3.0, 10.0, 2.0)  # clears x in [1, 5]; right half engaged
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(math.pi, abs=0.02)
    assert max_run == pytest.approx(math.pi, abs=0.02)  # one merged run, not two
    assert exceeded  # ~pi + slack > pi cap, orientation branch


def test_engagement_exact_pi_run():
    """Exact-pi run -> the sign primitive's ZERO branch (collinear center/start/end).

    Stock is the lower half-plane y <= 5, its top edge passing exactly through the
    cutter centre, so the engaged run is an exact semicircle with antipodal
    endpoints and orientation(center, p, q) == COLLINEAR. At cap == pi (ratio == 4)
    run == cap and it does NOT exceed; any tighter cap does. Exact equality, no
    tolerance.
    """
    lower_half = np.array([[0, 0, 0], [10, 0, 0], [10, 5, 0], [0, 5, 0]], dtype=np.float64)
    stock = _stock_2.Stock2(lower_half, [])
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(math.pi, abs=1e-9)
    assert max_run == pytest.approx(math.pi, abs=1e-9)
    assert not exceeded  # run == pi == cap: ZERO branch, cap_chord_ratio < 4 is False
    _, _, exceeded_tight = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.radians(179.0)))
    assert exceeded_tight  # run pi > 179 deg cap


def test_engagement_tangent_touch():
    """Tangent contact is a zero-measure event and must not alter the 2*pi run.

    A cleared disk internally tangent to the cutter rim at a single point removes
    only that measure-zero point; an externally tangent one removes nothing. Both
    leave the rim fully engaged.
    """
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.25, 0.25)  # inside cutter, tangent to rim at (5, 5.5)
    total, _, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total == pytest.approx(2.0 * math.pi, abs=1e-9)
    assert exceeded

    stock2 = _stock_2.Stock2(SQUARE, [])
    stock2.subtract_disk(6.0, 5.0, 0.5)  # externally tangent to rim at (5.5, 5)
    total2, _, _ = _stock_2.engagement_at(stock2, 5.0, 5.0, 0.5, cap_ratio(math.pi))
    assert total2 == pytest.approx(2.0 * math.pi, abs=1e-9)


def test_engagement_rejects_bad_cap_ratio():
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(ValueError):
        _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, 0.0)
    with pytest.raises(ValueError):
        _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, -1.0)
    with pytest.raises(ValueError):
        _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, 4.0 + 1e-9)


# --------------------------------------------------------------------------- #
# _sign_mixed_radical: exact sign of A + B*sqrt(a) + C*sqrt(b) + D*sqrt(a*b)   #
# --------------------------------------------------------------------------- #


def _ref_sign_two_term(A, B, alpha):
    """Exact sign of A + B*sqrt(alpha) over rationals (Fraction)."""
    sA = (A > 0) - (A < 0)
    sB = (B > 0) - (B < 0)
    if sB == 0:
        return sA
    if alpha == 0:
        return sA
    if sA == 0:
        return sB
    if sA == sB:
        return sA
    disc = A * A - B * B * alpha
    s = (disc > 0) - (disc < 0)
    return s if sA > 0 else -s


def _ref_sign_mixed(A, B, C, D, alpha, beta):
    """Exact sign of A + B*sqrt(a) + C*sqrt(b) + D*sqrt(a*b) over rationals."""
    if beta == 0:
        return _ref_sign_two_term(A, B, alpha)
    su = _ref_sign_two_term(A, B, alpha)
    sw = _ref_sign_two_term(C, D, alpha)
    if su == sw:
        return su
    if su == 0:
        return sw
    if sw == 0:
        return su
    Ap = A * A + B * B * alpha - beta * (C * C + D * D * alpha)
    Bp = 2 * A * B - 2 * beta * C * D
    s = _ref_sign_two_term(Ap, Bp, alpha)
    return s if su > 0 else -s


def _decimal_value(a, b, c, d, alpha, beta):
    """Independent 80-digit numeric evaluation (sign valid when bounded off 0)."""
    return Decimal(a) + Decimal(b) * Decimal(alpha).sqrt() + Decimal(c) * Decimal(beta).sqrt() + Decimal(d) * Decimal(alpha).sqrt() * Decimal(beta).sqrt()


# (a, b, c, d, alpha, beta, expected sign) — hand-verified, including exact zeros.
CURATED = [
    (1, 0, 0, 0, 0, 0, 1),
    (-3, 0, 0, 0, 0, 0, -1),
    (0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 0, 2, 0, 1),  # sqrt(2)
    (0, -1, 0, 0, 2, 0, -1),
    (-1, 1, 0, 0, 2, 0, 1),  # -1 + sqrt(2) > 0
    (-2, 1, 0, 0, 2, 0, -1),  # -2 + sqrt(2) < 0
    (0, 0, 1, 0, 0, 3, 1),  # sqrt(3)
    (1, 1, -1, 0, 2, 3, 1),  # 1 + sqrt(2) - sqrt(3) > 0
    (0, 1, -1, 0, 2, 3, -1),  # sqrt(2) - sqrt(3) < 0
    (1, -1, -1, 0, 2, 3, -1),  # 1 - sqrt(2) - sqrt(3) < 0
    (0, 0, 0, 1, 2, 3, 1),  # sqrt(6)
    (0, 0, 0, -1, 2, 3, -1),
    (0, 1, 0, -1, 2, 3, -1),  # sqrt(2) - sqrt(6) < 0
    (5, 1, 1, 1, 2, 3, 1),
    (3, -1, 0, 0, 9, 0, 0),  # 3 - sqrt(9) = 0
    (0, 1, -1, 0, 2, 2, 0),  # sqrt(2) - sqrt(2) = 0
    (2, 0, 0, -1, 1, 4, 0),  # 2 - sqrt(4) = 0            (general branch, compare EQUAL)
    (-2, 0, 0, 1, 4, 1, 0),  # -2 + sqrt(4) = 0           (general branch, compare EQUAL)
    (0, 1, -0.5, 0, 2, 8, 0),  # sqrt(2) - 0.5*sqrt(8) = 0 (u extended, opposite-sign EQUAL)
]


@pytest.mark.parametrize("a,b,c,d,alpha,beta,expected", CURATED)
def test_sign_mixed_radical_curated(a, b, c, d, alpha, beta, expected):
    assert _stock_2._sign_mixed_radical(a, b, c, d, alpha, beta) == expected


def test_sign_mixed_radical_randomized():
    """C++ exact primitive vs a Fraction-exact oracle and an 80-digit numeric check."""
    rng = random.Random(20260723)
    roots = [0, 1, 2, 3, 5, 6, 7, 8, 10, 11]
    for _ in range(3000):
        a, b, c, d = (rng.randint(-6, 6) for _ in range(4))
        alpha = rng.choice(roots)
        beta = rng.choice(roots)
        got = _stock_2._sign_mixed_radical(a, b, c, d, alpha, beta)
        exact = _ref_sign_mixed(*(Fraction(v) for v in (a, b, c, d, alpha, beta)))
        assert got == exact, (a, b, c, d, alpha, beta, got, exact)
        val = _decimal_value(a, b, c, d, alpha, beta)
        if abs(val) > Decimal("1e-12"):  # bounded off zero -> numeric sign is trustworthy
            assert got == (1 if val > 0 else -1), (a, b, c, d, alpha, beta, got, val)


def test_sign_mixed_radical_dyadic():
    """Non-integer (exactly-representable) coefficients: FT(double) stays exact."""
    rng = random.Random(4242)
    roots = [0, 1, 2, 3, 5, 6, 7]
    for _ in range(2000):
        a, b, c, d = (rng.randint(-8, 8) / 4.0 for _ in range(4))
        alpha = rng.choice(roots)
        beta = rng.choice(roots)
        got = _stock_2._sign_mixed_radical(a, b, c, d, alpha, beta)
        exact = _ref_sign_mixed(*(Fraction(v) for v in (a, b, c, d, alpha, beta)))
        assert got == exact, (a, b, c, d, alpha, beta, got, exact)
