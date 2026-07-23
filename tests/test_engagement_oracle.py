"""Monte-Carlo raster differential oracle for the exact station-TEA kernel.

This is Task 7 of SP1: a *verification net*, not a unit test of one function. It
builds a second, deliberately independent, model of the stock -- a pure-numpy
boolean occupancy grid -- and drives it through the *same* sequence of material
subtractions as the exact CGAL kernel (`compas_cgal._stock_2`). The two models
share no geometry code: the exact kernel is an ``Epeck`` circle-segment
arrangement with one-root point coordinates, the oracle is a bitmap with per-cell
Euclidean distance tests. When two independent implementations of the same
geometry agree, that is strong evidence neither harbours a systematic error the
other would also make -- the differential signal a single-implementation test
cannot give.

Two properties are cross-checked:

* **Tool-engagement angle (TEA).** For a probe circle of radius ``r`` centred at
  ``(cx, cy)``, the exact kernel reports the total angular measure of the rim
  still lying in material. The oracle estimates the same quantity by Monte-Carlo:
  the fraction of ``N_ANGLE`` equally spaced rim samples whose grid cell is
  material, times ``2*pi``. `test_tea_agrees_with_raster_oracle` asserts the two
  agree within a tolerance derived entirely from the oracle's discretisation (see
  that test's docstring for the budget).
* **Growth-bound soundness.** The certifier bounds how fast TEA can change as the
  probe centre moves by a distance ``d`` (the Task 5 Lipschitz-style lemma).
  `test_growth_bound_sound` samples the exact kernel at pairs of nearby centres
  and asserts the observed change never exceeds the raw lemma bound.

Chord-ratio API (Task 4 amendment): `_stock_2.engagement_at` consumes the
caller's exact rational surrogate ``4*sin^2(cap/2)`` for the transcendental
angular cap, not the cap in radians. Both tests use the full-engagement cap
``cap = pi``, whose surrogate is exactly ``4.0`` -- the whole engaged rim counts,
so the reported ``total_tea`` is the uncapped engaged angle the oracle measures.
"""

import math

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

# --- Oracle discretisation -------------------------------------------------
# Grid cell size on a 10-unit-wide pocket: 500x500 cells. H sets both the
# spatial resolution of the occupancy grid and the boundary-cell ambiguity of
# the TEA estimate (derived in test_tea_agrees_with_raster_oracle).
H = 0.02
# Rim samples per probe circle. 1440 == 0.25 deg spacing; the Monte-Carlo TEA
# estimate is granular to one step, a hard floor of 2*pi/N_ANGLE radians.
N_ANGLE = 1440

# Full-engagement angular cap and its exact rational chord surrogate. cap in
# (0, pi] maps to 4*sin^2(cap/2) in (0, 4]; cap = pi => surrogate = 4.0 (no cap:
# the whole engaged rim counts, so total_tea is the raw uncapped engaged angle).
CAP = math.pi
CAP_CHORD_RATIO = 4.0 * math.sin(CAP / 2.0) ** 2

# Empirical TEA slack absorbing (a) the exact kernel's certified chain
# under-coverage -- capsule/arc subtraction chains remove <= CHAIN_SLACK_FRACTION
# * r = 1e-4 * r LESS than the ideal swept region (a one-sided conservative bias
# so the removal certificate never over-claims), leaving the exact stock with
# marginally more material and hence marginally higher TEA than the oracle's
# ideal-capsule raster; and (b) residual boundary quantisation when several
# engaged arcs are active beyond the single-run 4*H/r term.
CHAIN_SLACK_TEA = 0.01

# Fraction of `tol` the worst observed exact-vs-oracle gap is pinned under, so an
# erosion of the agreement headroom (observed ~0.80) surfaces as a failure rather
# than silently consuming the raw tolerance. Not a geometric tolerance -- a
# regression sentinel on the margin between observation and the derived bound.
TEA_AGREEMENT_HEADROOM = 0.90

# Pure numerical floors (NOT geometric tolerances): guard a degenerate
# zero-length segment against divide-by-zero, and give the abs() float
# comparison of the growth bound a few ulp of slack.
_SEGMENT_L2_FLOOR = 1e-30
_FLOAT_SLACK = 1e-9


class RasterStock:
    """Independent occupancy-grid stock model -- same subtraction semantics, no CGAL.

    A boolean bitmap of the pocket: ``material[j, i]`` is True while the cell
    whose centre is ``((i + 0.5) * h, (j + 0.5) * h)`` still holds stock. Each
    ``subtract_*`` clears every cell whose centre lies inside the removed region
    via an exact per-cell distance test -- the *ideal* swept area with no
    under-coverage, so the oracle is a faithful ground truth against which the
    exact kernel's conservatively under-covering chains are compared.
    """

    def __init__(self, size: float = 10.0, h: float = H) -> None:
        n = int(round(size / h))
        self.h = h
        xs = (np.arange(n) + 0.5) * h
        self.gx, self.gy = np.meshgrid(xs, xs)
        self.material = np.ones_like(self.gx, dtype=bool)

    def subtract_disk(self, cx: float, cy: float, r: float) -> None:
        """Clear every cell whose centre lies within ``r`` of ``(cx, cy)``."""
        self.material &= (self.gx - cx) ** 2 + (self.gy - cy) ** 2 > r * r

    def subtract_capsule(self, x0: float, y0: float, x1: float, y1: float, r: float) -> None:
        """Clear every cell within ``r`` of the segment ``(x0, y0)-(x1, y1)``."""
        dx, dy = x1 - x0, y1 - y0
        l2 = max(dx * dx + dy * dy, _SEGMENT_L2_FLOOR)
        t = np.clip(((self.gx - x0) * dx + (self.gy - y0) * dy) / l2, 0.0, 1.0)
        px, py = x0 + t * dx, y0 + t * dy
        self.material &= (self.gx - px) ** 2 + (self.gy - py) ** 2 > r * r

    def tea_at(self, cx: float, cy: float, r: float) -> float:
        """Monte-Carlo tool-engagement angle: material rim fraction times ``2*pi``."""
        ang = np.linspace(0.0, 2.0 * math.pi, N_ANGLE, endpoint=False)
        x = cx + r * np.cos(ang)
        y = cy + r * np.sin(ang)
        i = np.clip((x / self.h).astype(int), 0, self.material.shape[1] - 1)
        j = np.clip((y / self.h).astype(int), 0, self.material.shape[0] - 1)
        return float(self.material[j, i].mean()) * 2.0 * math.pi


SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


@pytest.mark.parametrize("seed", range(8))
def test_tea_agrees_with_raster_oracle(seed: int) -> None:
    """Exact station TEA agrees with the independent raster oracle.

    Drive both the exact kernel and the numpy occupancy grid through an identical
    random chain of six capsule subtractions, then compare their TEA at 40 random
    probe centres. Two implementations sharing no geometry code agreeing to
    ``tol`` is the differential evidence this net exists to provide.

    Tolerance budget (``tol``), every term in radians of engagement:

    * Rim quantisation ``2*pi/N_ANGLE``: the oracle samples the rim at that
      spacing (~4.4 mrad at N_ANGLE=1440), so its estimate is granular to one
      step regardless of geometry.
    * Boundary-cell ambiguity ``4*H/r`` (a TRANSVERSAL / generic-position budget):
      at each end of an engaged arc that crosses a rasterised material boundary at
      a nonzero angle, whether a rim sample reads material is decided by an
      ``H``-sized cell. Both the sample's cell-centre position and the rasterised
      subtraction boundary quantise by up to ``H``, leaving the true arc endpoint
      uncertain by up to ``2*H/r`` rad; an engaged arc has two endpoints, so up to
      ``4*H/r`` rad (0.16 rad at H=0.02, r=0.5) -- the dominant term. This is a
      transversal-crossing bound. GRAZING INCIDENCE -- a material boundary running
      nearly tangent to the rim -- is the KNOWN gap it does not cover: there a
      whole ``O(sqrt(H/r))`` run of rim samples flips together, an error that can
      exceed ``4*H/r``. This suite samples random probe centres, so grazing
      alignments are measure-zero and do not arise; ``CHAIN_SLACK_TEA`` below
      carries what residual multi-run quantisation the random draws do produce.
    * Chain under-coverage + multi-run residual ``CHAIN_SLACK_TEA``: the exact
      kernel's subtraction chains under-cover the ideal sweep by
      <= CHAIN_SLACK_FRACTION * r = 1e-4 * r (leaving slightly more material, so
      slightly higher exact TEA), and several arcs may be engaged at once beyond
      the single-run boundary budget; ``CHAIN_SLACK_TEA = 0.01`` rad absorbs both.

    => ``tol = 2*pi/N_ANGLE + 4*H/r + CHAIN_SLACK_TEA``.

    HEADROOM PIN. The largest exact-vs-oracle disagreement observed over this
    suite is ``~0.80 * tol`` (measured 2026-07-23 on the post-merge-repair build,
    all 8 seeds). ``TEA_AGREEMENT_HEADROOM`` asserts every probe stays under
    ``0.90 * tol`` so a regression that erodes the ~20% headroom -- e.g. a grazing
    alignment creeping in, or the transversal budget being consumed -- surfaces
    here instead of silently riding the raw ``tol`` bound.
    """
    rng = np.random.default_rng(seed)
    exact = Stock(SQUARE)
    raster = RasterStock()
    for _ in range(6):
        x0, y0, x1, y1 = rng.uniform(1.5, 8.5, 4)
        r = float(rng.uniform(0.3, 0.9))
        exact.subtract_capsule(x0, y0, x1, y1, r)
        raster.subtract_capsule(x0, y0, x1, y1, r)
    r_tool = 0.5
    tol = 2 * math.pi / N_ANGLE + 4 * H / r_tool + CHAIN_SLACK_TEA
    max_gap = 0.0
    for _ in range(40):
        cx, cy = rng.uniform(1.0, 9.0, 2)
        total, _, _ = _stock_2.engagement_at(exact.raw, cx, cy, r_tool, CAP_CHORD_RATIO)
        oracle = raster.tea_at(cx, cy, r_tool)
        assert total == pytest.approx(oracle, abs=tol)
        max_gap = max(max_gap, abs(total - oracle))
    assert max_gap <= TEA_AGREEMENT_HEADROOM * tol  # ~20% headroom pinned (observed ~0.80)


def test_growth_bound_sound() -> None:
    """The Task 5 growth-bound lemma is empirically sound at safety factor 1.

    As the probe centre moves by ``d``, TEA cannot change by more than the raw
    lemma bound ``4*asin(min(1, d/2r)) + 2*acos(max(-1, 1 - d/r))``. This asserts
    the *lemma itself* -- the certifier's own factor-2 safety margin is
    deliberately NOT applied here (Task 5 amendment): the net verifies the bound
    is sound before any margin, so that a regression eroding the margin still
    leaves a live, un-padded bound for this test to catch. 200 random
    (centre, displacement) pairs on a cut-up stock; ``d`` in (1e-4, 0.2) keeps
    both inverse-trig arguments inside their domains without relying on the clamp.
    """
    rng = np.random.default_rng(11)
    exact = Stock(SQUARE)
    for _ in range(5):
        x0, y0, x1, y1 = rng.uniform(1.5, 8.5, 4)
        exact.subtract_capsule(x0, y0, x1, y1, float(rng.uniform(0.4, 1.0)))
    r = 0.5
    for _ in range(200):
        cx, cy = rng.uniform(1.5, 8.5, 2)
        d = float(rng.uniform(1e-4, 0.2))
        theta = float(rng.uniform(0, 2 * math.pi))
        t0, _, _ = _stock_2.engagement_at(exact.raw, cx, cy, r, CAP_CHORD_RATIO)
        t1, _, _ = _stock_2.engagement_at(exact.raw, cx + d * math.cos(theta), cy + d * math.sin(theta), r, CAP_CHORD_RATIO)
        bound = 4 * math.asin(min(1.0, d / (2 * r))) + 2 * math.acos(max(-1.0, 1.0 - d / r))
        assert abs(t1 - t0) <= bound + _FLOAT_SLACK


# --------------------------------------------------------------------------- #
# Run-merge discontinuity and its gap-closure repair (Task-5 certificate hole) #
# --------------------------------------------------------------------------- #

# Merge witness stock: lower half-plane (y <= 5) with a thin central void slot
# biting the rim bottom. A radius-0.5 cutter at cy = 5 engages the lower
# semicircle; the slot splits it into two runs, which MERGE into one ~pi run as
# the centre slides in x and the slot leaves the rim near cx ~ 4.71. There the
# growth lemma (bounded above by test_growth_bound_sound for a SINGLE run) is
# violated by an O(1) jump in max_run -- the hole the certificate had to close.
_MERGE_LOWER_HALF = np.array([[0, 0, 0], [10, 0, 0], [10, 5, 0], [0, 5, 0]], dtype=np.float64)
_MERGE_LOCUS_X = 4.7125  # cx where the two runs fuse (max_run jumps ~pi/2 -> pi)
_MERGE_R = 0.5


def _merge_stock():
    stock = _stock_2.Stock2(_MERGE_LOWER_HALF, [])
    stock.subtract_capsule(5.0, 4.30, 5.0, 4.53, 0.05)
    return stock


def _cap_ratio(cap: float) -> float:
    """Caller-side exact squared-chord surrogate 4*sin^2(cap/2) for an angular cap/gamma."""
    return 4.0 * math.sin(cap / 2.0) ** 2


def _growth_bound(d: float, r: float) -> float:
    """Raw factor-1 TEA-growth lemma (mirrors tea_growth_bound in engagement_2.cpp)."""
    return 4.0 * math.asin(min(1.0, d / (2.0 * r))) + 2.0 * math.acos(max(-1.0, 1.0 - d / r))


def _pess_max_run(stock, cx: float, r: float, gamma_ratio: float, iters: int = 34) -> float:
    """Largest PESSIMISTIC run angle at a station (void gaps <= gamma absorbed).

    The reported ``max_run`` stays the TRUE (unclosed) measure by design, so the
    pessimistic span is recovered through the exact DECISION: bisect the cap for the
    ``cap_exceeded`` flip at ``gap_close_ratio = gamma_ratio``. Saturates at ``pi`` --
    a run ``>= pi`` exceeds every cap ``<= pi`` (the exact orientation branch), so a
    merged (``~pi``) pessimistic run reads as ``pi`` here.
    """
    if _stock_2.engagement_at(stock, cx, 5.0, r, _cap_ratio(math.pi), gamma_ratio)[2]:
        return math.pi
    lo, hi = 1e-4, math.pi
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if _stock_2.engagement_at(stock, cx, 5.0, r, _cap_ratio(mid), gamma_ratio)[2]:
            lo = mid  # exceeded -> pessimistic max_run > mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def test_merge_jump_exceeds_growth_lemma() -> None:
    """The run-merge discontinuity is real: raw max_run jumps by more than the growth
    lemma allows.

    An INVERTED assertion -- it documents the HOLE, not the fix. Sweeping the cutter
    centre across the merge locus with raw engagement (``gap_close_ratio = 0``), the
    largest step-to-step change in ``max_run`` is an O(1) jump (~pi/2 rad: two ~pi/2
    runs fusing into one ~pi run) even though the step ``d`` is tiny. The raw lemma
    bound ``4*asin(min(1,d/2r)) + 2*acos(max(-1,1-d/r))`` is O(sqrt d) and cannot cover
    it, so a certificate bridging stations by that lemma alone (the pre-repair Task-5
    design) is unsound across merges. The repair removes the event at the station
    (test_pessimistic_max_run_growth_is_sound;
    tests/test_stock.py::test_certificate_refuses_merge_over_cap).

    The geometry is fixed, not tuned to the outcome: if this jump ever fell within the
    bound the merge hole would be illusory and the whole repair suspect -- a BLOCKED
    signal, not something to soften.
    """
    stock = _merge_stock()
    r = _MERGE_R
    step = 0.005
    xs = [4.700 + i * step for i in range(6)]  # 4.700 .. 4.725, straddling the merge ~4.7125
    runs = [_stock_2.engagement_at(stock, cx, 5.0, r, CAP_CHORD_RATIO)[1] for cx in xs]
    max_jump = max(abs(runs[i + 1] - runs[i]) for i in range(len(runs) - 1))
    bound = _growth_bound(step, r)
    assert max_jump > bound  # ~1.04 rad jump dwarfs the ~0.30 rad O(sqrt d) lemma bound
    assert max_jump > math.radians(45.0)  # a genuine O(1) merge, not boundary quantisation


def test_pessimistic_max_run_growth_is_sound() -> None:
    """The quantity the repaired certificate NOW relies on -- the PESSIMISTIC max_run
    at a station -- grows within the raw lemma across the merge, where the true max_run
    does not.

    For probe pairs straddling the merge locus at spacing ``d``, the certifier measures
    each station with void gaps ``<= gamma = 2*GROWTH(d/2)`` pre-absorbed. Both probes
    then see the closing gap as already shut, report the merged ~pi run, and the
    pessimistic max_run changes by ~0 across the pair -- inside the raw lemma bound. The
    TRUE max_run over the same pair jumps by ~pi/2 (the discontinuity of
    test_merge_jump_exceeds_growth_lemma), so the lemma bounds this quantity only AFTER
    the pessimism is applied. Sound at safety factor 1 (raw lemma, no guard margin),
    mirroring test_growth_bound_sound. If the pessimistic delta ever exceeded the bound,
    the certificate's new load-bearing quantity would itself be unbounded -- report
    BLOCKED, do not widen the bound.
    """
    stock = _merge_stock()
    r, m = _MERGE_R, _MERGE_LOCUS_X
    for d in (0.01, 0.02, 0.04):
        gamma = min(2.0 * _growth_bound(0.5 * d, r), math.pi)  # certifier's guard for hs = d/2
        gr = _cap_ratio(gamma)
        a, b = m - 0.5 * d, m + 0.5 * d
        pa = _pess_max_run(stock, a, r, gr)
        pb = _pess_max_run(stock, b, r, gr)
        bound = _growth_bound(d, r)
        assert abs(pa - pb) <= bound + _FLOAT_SLACK  # pessimism keeps the growth lemma-bounded
        # nonvacuous: the TRUE max_run across the SAME pair exceeds the lemma (the jump absorbed)
        ta = _stock_2.engagement_at(stock, a, 5.0, r, CAP_CHORD_RATIO)[1]
        tb = _stock_2.engagement_at(stock, b, 5.0, r, CAP_CHORD_RATIO)[1]
        assert abs(ta - tb) > bound
