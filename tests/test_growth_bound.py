"""Falsification harness for the analytic TEA-growth bound.

``tea_growth_bound(d, r)`` is the ONLY thing standing between the guarded-station
method and an unsound certificate: ``certify_segment_tea`` subtracts
``2 * tea_growth_bound(half_spacing, r)`` from the cap and concludes that two
passing stations certify everything between them. If the bound under-estimates
true growth anywhere, that conclusion is false.

These tests measure true growth with the EXACT kernel and compare it against the
SHIPPED bound -- ``_stock_2.tea_growth_bound``, never a local copy of the formula
-- across the boundary shapes the stock model actually produces. They are
deliberately adversarial.

What a verdict means here, precisely:

* A **red** is unconditional evidence: a measured cutter displacement of exactly
  ``d`` grew the largest engaged run by more than the bound permits, so the bound
  is not an upper bound and any certificate resting on it is unsound.
* A **green** is bounded evidence: no violation was found among the probed
  centres and the sampled displacement directions. Sampling can only ever be a
  LOWER bound on the worst growth, so each green test additionally pins the
  configuration it claims to measure (the probe really straddles the feature, the
  measured growth really is non-zero) -- a green that came from an empty arc list
  or an uncut stock would be worthless.

Every measurement is the RAW factor-1 lemma against the RAW true growth: the
certifier's own factor-2 safety margin is deliberately not applied, so a
regression eroding that margin still leaves a live, un-padded bound here, and
gap-closure pessimism is switched off (``GAP_CLOSE_NONE``) because the run-merge
discontinuity is a separate, already-repaired hole (see
``tests/test_engagement_oracle.py::test_merge_jump_exceeds_growth_lemma``) and
must not be confused with a growth-lemma failure.

Scale invariance: ``tea_growth_bound`` depends on its arguments only through
``d / r``, and the true growth of every feature modelled here depends only on
``d / r`` and ``rho / r``. One tool radius with parametrised ratios therefore
covers the whole family; a second radius would add no information.

One test here is NOT about the bound. Measuring growth surfaced a separate,
independent defect in the kernel's engagement harvest -- stations at which
``engagement_at`` reports an engaged run larger than a full turn, which is
geometrically impossible -- so `test_reported_engagement_never_exceeds_a_full_turn`
pins it and `_max_run_tea` refuses such a reading at the point of measurement.
Without that split, a growth number computed from an impossible engagement would
be read as a growth-bound failure it is not.
"""

from __future__ import annotations

import math

import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

# Tool radius r, in model units, for every probe in this module (see the module
# docstring on scale invariance).
TOOL_RADIUS = 0.5

# Engagement cap handed to `engagement_at` as its exact squared-chord surrogate
# 4*sin^2(cap/2): cap = pi/2 (90 deg) -> ratio 2.0, inside the required (0, 4].
# The cap selects only the `cap_exceeded` verdict, which this module never reads:
# `max_run_tea` is the TRUE largest engaged run at every admissible cap
# (engagement_2.cpp, step 4a REPORTING), so the choice cannot bias a measurement.
CAP_RATIO_90 = 4.0 * math.sin(math.pi / 4.0) ** 2

# Gap-closure surrogate 4*sin^2(gamma/2) for gamma = 0: absorb no void gap, i.e.
# measure the TRUE runs. The growth lemma is a statement about the true runs; the
# pessimistic (gap-closed) runs are the run-merge repair and are measured
# elsewhere. Mixing the two would hide a growth failure behind the repair.
GAP_CLOSE_NONE = 0.0

# Directions sampled on the ring of radius d around a probe centre. 360 == 1 deg
# spacing AND divisible by 4, so the four cardinal directions -- where the
# extremal growth of every feature in this module lies (radially out of a void,
# straight into a wall) -- are sampled EXACTLY rather than merely approached.
PROBE_DIRECTIONS = 360

# A full turn of engaged rim (radians): the reading when the whole cutter rim is
# in material.
FULL_TURN = 2.0 * math.pi

# Slack for comparing a measured TEA against its closed-form value. Reported spans
# are doubles assembled by atan2 from exact one-root arc endpoints; the kernel's
# own zone-vs-overlay study bounds that representation artefact at <= 1e-15 rad
# (docs/superpowers/state/engagement-zone-divergence.md). 1e-12 rad == 6e-11 deg
# keeps three decades of headroom over the artefact while still catching any real
# geometric error, which is orders of magnitude larger. Used ONLY on the
# configuration pins, never on the bound comparison -- that one is exact-vs-exact
# and takes no tolerance.
TEA_REPORTING_SLACK = 1e-12

# Ambient block half-width, in model units: 12x the tool radius, so the block's
# own straight walls stay >= 4.5 units clear of every probe below and never
# contribute an engaged run. Each test then measures exactly ONE boundary feature.
BLOCK_HALF_WIDTH = 12.0 * TOOL_RADIUS

_BLOCK = Polygon(
    [
        (-BLOCK_HALF_WIDTH, -BLOCK_HALF_WIDTH, 0),
        (BLOCK_HALF_WIDTH, -BLOCK_HALF_WIDTH, 0),
        (BLOCK_HALF_WIDTH, BLOCK_HALF_WIDTH, 0),
        (-BLOCK_HALF_WIDTH, BLOCK_HALF_WIDTH, 0),
    ]
)

# Half-plane stock: the same block cut off at y = 0, leaving one straight wall.
_HALF_PLANE = Polygon(
    [
        (-BLOCK_HALF_WIDTH, -BLOCK_HALF_WIDTH, 0),
        (BLOCK_HALF_WIDTH, -BLOCK_HALF_WIDTH, 0),
        (BLOCK_HALF_WIDTH, 0, 0),
        (-BLOCK_HALF_WIDTH, 0, 0),
    ]
)

# Probe centre heights above the wall at y = 0, in model units, covering the whole
# engagement range of the cutter against a half-plane: 0.0 centred ON the wall
# (TEA = pi), 0.1 and 0.25 partially engaged, 0.5 == r exactly tangent from the
# clear side (TEA = 0 -- the configuration where clause (b) of the lemma is
# asymptotically TIGHT, growth 2*acos(1 - d/r)), 0.75 == 1.5r fully clear.
WALL_PROBE_HEIGHTS = (0.0, 0.1, 0.25, 0.5, 0.75)

# Clearance past a tangency, in model units (0.2r): far enough that the feature is
# unambiguously out of contact, near enough to stay a local probe.
CLEAR_MARGIN = 0.2 * TOOL_RADIUS

# Depleted-stock sweep (the shipped-audit regime): parallel slots cut by the SAME
# radius the query uses. Six passes of a full tool diameter step-over 0.8 == 1.6r
# < 2r, so consecutive capsules OVERLAP and the union is one connected pocket with
# the scalloped, part-circular boundary a real roughing pass leaves behind.
SWEEP_PASSES = 6
SWEEP_STEPOVER = 1.6 * TOOL_RADIUS
SWEEP_FIRST_X = -2.0
SWEEP_Y0 = -1.0
SWEEP_Y1 = 1.0
# Probe heights along each pass, in model units: both capsule ends, the midpoint,
# and the quarter points.
SWEEP_PROBE_HEIGHTS = (-1.0, -0.5, 0.0, 0.5, 1.0)

# Void radius for the impossible-engagement witnesses below, in model units
# (0.2r -- small enough that the rim only ever takes a bite out of an otherwise
# fully engaged turn, which is what makes the surplus turn unmistakable).
FULL_TURN_WITNESS_VOID_RADIUS = 0.2 * TOOL_RADIUS

# Cutter-centre stations at which `engagement_at` reports an engaged run LARGER
# than a full turn against `_void_stock(FULL_TURN_WITNESS_VOID_RADIUS)`. Verbatim
# doubles, not rounded decimals: the defect is representation-sensitive, and a
# neighbouring station a milliradian away reports correctly.
#
# All three are configurations where the void's crossing of the rim coincides with
# the rim's x-extreme ``(cx - r, cy)`` -- the point `make_x_monotone_2` splits the
# cutter circle at. The first two were found on the ``d == rho`` probe ring of
# test_bound_holds_against_a_void_smaller_than_the_tool, where a base offset of
# ``s = r`` puts the rim x-extreme exactly ``rho`` from the void centre in EVERY
# direction; the third was constructed independently from that condition. At each
# the report exceeds the closed-form engagement by EXACTLY one full turn
# (2*pi to 1e-14), which is what the harvest's ``if (span <= 0) span += 2*pi``
# normalisation produces when applied to the sub-ulp sub-arc such a coincidence
# leaves behind (engagement_2.cpp, engaged_arcs_zone).
FULL_TURN_WITNESS_STATIONS = (
    (0.5981627183447664, -0.019080899537654468),
    (0.5961261695938319, -0.02756373558169998),
    (0.5714142842854285, 0.07),
)


def _max_run_tea(stock: Stock, x: float, y: float, r: float = TOOL_RADIUS) -> float:
    """Exact largest contiguous engaged run at a station, in radians (reporting).

    Args:
        stock: The stock region to measure against (frozen: nothing is cut here).
        x: X coordinate of the cutter centre.
        y: Y coordinate of the cutter centre.
        r: Tool radius.

    Returns:
        The largest contiguous engaged-run angle in radians; ``0.0`` when the rim
        lies nowhere in material.
    """
    _total, max_run, _exceeded = _stock_2.engagement_at(stock.raw, x, y, r, CAP_RATIO_90, GAP_CLOSE_NONE)
    # Contract check on the kernel's OWN answer -- not a tolerance, and nothing is
    # clamped or filtered. A contiguous run on a closed rim cannot exceed a full
    # turn, so a larger reading is the harvest defect
    # `test_reported_engagement_never_exceeds_a_full_turn` pins, and every growth
    # computed from it is meaningless. Failing HERE, naming the station, stops that
    # reading from surfacing downstream as a growth-bound violation it is not.
    assert max_run <= FULL_TURN + TEA_REPORTING_SLACK, (
        f"engagement_at reported max_run_tea = {max_run:.9f} rad > 2*pi at station ({x!r}, {y!r}), r={r!r}: "
        f"a contiguous run cannot exceed a full turn -- kernel harvest defect, not a growth-bound failure"
    )
    return float(max_run)


def _worst_growth_over(stock: Stock, x0: float, y0: float, d: float, samples: int = PROBE_DIRECTIONS) -> float:
    """Largest TEA increase from ``(x0, y0)`` to any probed point exactly *d* away.

    The quantity ``tea_growth_bound(d, r)`` claims to bound. Sampling the ring of
    displacement directions gives a LOWER bound on the true worst growth, so a
    returned value above the bound falsifies it outright, while a value below it
    only means no violation was found at this angular resolution.

    Args:
        stock: The stock region to measure against.
        x0: X coordinate of the base station.
        y0: Y coordinate of the base station.
        d: Centre-travel distance (the displacement magnitude).
        samples: Number of equally spaced displacement directions.

    Returns:
        The largest observed increase in the largest engaged run, in radians;
        ``0.0`` when no probed direction increases it.
    """
    base = _max_run_tea(stock, x0, y0)
    worst = 0.0
    for i in range(samples):
        a = 2.0 * math.pi * i / samples
        worst = max(worst, _max_run_tea(stock, x0 + d * math.cos(a), y0 + d * math.sin(a)) - base)
    return worst


def _void_stock(rho: float) -> Stock:
    """Block with one circular void of radius *rho* centred at the origin.

    Args:
        rho: Void radius in model units.

    Returns:
        A `Stock` whose only boundary feature within reach of the probes is the
        void's circular wall.
    """
    stock = Stock(_BLOCK)
    stock.subtract_disk(0.0, 0.0, rho)
    return stock


def _sweep_stock() -> Stock:
    """Block depleted by parallel same-radius passes (the shipped-audit regime).

    Returns:
        A `Stock` whose pocket boundary is built exclusively from arcs of radius
        ``TOOL_RADIUS`` and the straight flanks between them -- the only shapes
        ``Stock2::subtract_capsule`` can produce at that radius.
    """
    stock = Stock(_BLOCK)
    for k in range(SWEEP_PASSES):
        x = SWEEP_FIRST_X + k * SWEEP_STEPOVER
        stock.subtract_capsule(x, SWEEP_Y0, x, SWEEP_Y1, TOOL_RADIUS)
    return stock


def _small_void_probe_offsets(rho: float) -> tuple[float, ...]:
    """Probe centre offsets along +x for a void strictly smaller than the tool.

    Covers the complete contact life of the rim against such a void: concentric,
    centre on the void rim, rim entering the void, rim at its DEEPEST inside it,
    rim externally tangent, and fully clear. The two interior offsets
    (``r - rho`` and ``r``) are what make the case non-vacuous: at every other
    offset the void is either swallowed whole by the cutter disk or out of
    contact, both of which read a full turn and cannot grow.

    Args:
        rho: Void radius, strictly less than ``TOOL_RADIUS``.

    Returns:
        Offsets from the void centre, in model units.
    """
    return (
        0.0,
        rho,
        TOOL_RADIUS - rho,
        TOOL_RADIUS,
        rho + TOOL_RADIUS,
        rho + TOOL_RADIUS + CLEAR_MARGIN,
    )


@pytest.mark.parametrize("station", FULL_TURN_WITNESS_STATIONS)
def test_reported_engagement_never_exceeds_a_full_turn(station):
    """A reported engaged run larger than a full turn is impossible -- and the kernel reports one.

    NOT a growth-bound test. `engagement_at`'s ``total_tea`` and ``max_run_tea``
    are the REPORTING half of the deciding/reporting split (engagement_2.cpp, step
    4a): doubles summed from ``atan2`` spans of exact one-root arc endpoints. Their
    ceiling is structural -- the rim is a closed curve of angular measure ``2*pi``,
    so no contiguous run on it can measure more, whatever the stock looks like.

    At the stations below the kernel reports the true engagement plus exactly one
    full turn. Surfaced while measuring growth against a void smaller than the tool
    (a growth of 6.62 rad is impossible on its face: growth is a difference of two
    quantities in ``[0, 2*pi]``), so it is pinned separately here rather than left
    to masquerade as a growth-bound violation. Consumers inherit it: the audit's
    reported ``max_tea`` is this number.

    The cap DECISION is not implicated by this evidence -- it runs on exact
    predicates over the run ENDPOINTS (step 4b), never on these doubles -- but a
    harvest that mis-assembles a run's span is not obviously assembling its
    endpoints correctly either, and nothing here establishes that it does.
    """
    cx, cy = station
    stock = _void_stock(FULL_TURN_WITNESS_VOID_RADIUS)
    total, max_run, _exceeded = _stock_2.engagement_at(stock.raw, cx, cy, TOOL_RADIUS, CAP_RATIO_90, GAP_CLOSE_NONE)

    # Configuration pin: the station really is in the regime it claims -- the cutter
    # straddles the void wall, so SOME rim is engaged and some is not.
    assert total > 0.0

    assert max_run <= FULL_TURN + TEA_REPORTING_SLACK, (
        f"station ({cx!r}, {cy!r}): max_run_tea {max_run:.9f} rad exceeds a full turn by {max_run - FULL_TURN:.9f} rad (total_tea {total:.9f})"
    )


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
def test_bound_holds_against_a_straight_wall(d):
    """The half-plane case: clause (b) is asymptotically tight here, so this is the floor.

    A straight wall is the shape clause (b) was derived from -- the deepest first
    bite over travel ``d`` reaches radial depth ``d`` and cuts a chord spanning
    ``2*acos(1 - d/r)``. Growth here therefore approaches the bound from below as
    ``d`` shrinks (the clause-(a) drift term is the whole remaining slack) without
    ever crossing it. If this case ever went red the lemma would be wrong at its
    own defining geometry, not merely incomplete on curved features.
    """
    stock = Stock(_HALF_PLANE)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)

    # Configuration pin: the probe set really straddles the wall. The deepest
    # probe sits ON it (half the rim in material) and the clear probe is out of
    # contact -- so a growth of zero below could only be a real measurement.
    engagements = [_max_run_tea(stock, 0.0, y) for y in WALL_PROBE_HEIGHTS]
    assert max(engagements) == pytest.approx(math.pi, abs=TEA_REPORTING_SLACK)
    assert min(engagements) == 0.0

    worst = max(_worst_growth_over(stock, 0.0, y, d) for y in WALL_PROBE_HEIGHTS)
    assert worst > 0.0, f"straight wall: no probe grew at all at d={d} -- the harness measured nothing"
    assert worst <= bound, f"straight wall: growth {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
@pytest.mark.parametrize("rho", [0.1, 0.25])
def test_bound_holds_against_a_void_smaller_than_the_tool(d, rho):
    """rho < r: the rim can never sit inside the void, so the feature behaves convexly.

    A void strictly smaller than the tool cannot contain the cutter, so there is no
    tangency for the rim to emerge from: the void merely takes a bite out of an
    otherwise fully engaged rim, and that bite shrinks and grows continuously with
    the centre. Growth is bounded by clause (a) drift of the two bite endpoints,
    well inside the lemma.

    ``rho == d`` is the one case that does not reach its growth comparison: the
    ``s = r`` probe's displacement ring then puts the rim's x-extreme exactly
    ``rho`` from the void centre in every direction, and some of those stations
    trip the impossible-engagement defect
    (`test_reported_engagement_never_exceeds_a_full_turn`). `_max_run_tea` refuses
    the reading, so that parametrisation reports the harvest defect by name
    instead of a growth number derived from it. The probe stays: it is the offset
    at which the rim cuts deepest into the void, and dropping it would make this
    test vacuous -- every other offset reads a full turn at both ends and cannot
    grow at all.
    """
    assert rho < TOOL_RADIUS  # the regime this test claims to cover
    stock = _void_stock(rho)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)

    # Configuration pin: the probe set really straddles the void wall. At s = r the
    # rim passes through the void centre and is at its deepest inside it, so the
    # engaged run must fall short of a full turn; past external tangency the void
    # is out of contact and the run must be a full turn exactly.
    assert _max_run_tea(stock, TOOL_RADIUS, 0.0) < FULL_TURN - TEA_REPORTING_SLACK
    assert _max_run_tea(stock, rho + TOOL_RADIUS + CLEAR_MARGIN, 0.0) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)

    worst = max(_worst_growth_over(stock, s, 0.0, d) for s in _small_void_probe_offsets(rho))
    assert worst > 0.0, f"rho={rho}: no probe grew at all at d={d} -- the harness measured nothing"
    assert worst <= bound, f"rho={rho}: growth {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
def test_bound_holds_against_a_void_the_size_of_the_tool(d):
    """rho == r: the tool sitting in the hole it just cut -- the regime's singular limit.

    This is a plunge followed by a departure, the most ordinary motion pair a
    toolpath contains: ``subtract_disk(c, r)`` then a cut starting at ``c``. The
    cutter is internally tangent to its own hole, so the base station reads TEA = 0,
    and the emerged half-angle at travel ``delta`` satisfies

        cos(psi) = 1 - (2*delta*rho + delta**2) / (2*r*(rho - r + delta))

    which at ``rho == r`` collapses to ``cos(psi) = -delta/(2r)``: the run appears at
    just over ``pi`` for ANY positive travel, however small. It is the ``rho -> r+``
    limit of the concave family below, where ``true/bound = sqrt(rho/(rho - r))``
    diverges -- so this is the WORST case of the family, not a safe boundary of it.

    That contradicts the premise that the shipped audit is protected by every
    ``Stock2::subtract_*`` removing disks of the query's own radius: ``rho == r``
    does not make the failing regime unreachable, it makes it maximal. What keeps
    ``test_bound_holds_on_stock_the_shipped_audit_actually_produces`` green is
    where its probes sit, not the equal radii.
    """
    rho = TOOL_RADIUS
    stock = _void_stock(rho)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)

    # Configuration pin: the cutter really is sitting in its own hole (rim exactly
    # on the void wall, nothing engaged) and there really is material to emerge
    # into one displacement away.
    assert _max_run_tea(stock, 0.0, 0.0) == 0.0
    assert _max_run_tea(stock, d, 0.0) > 0.0

    worst = _worst_growth_over(stock, 0.0, 0.0, d)
    assert worst <= bound, f"void rho=r: growth {worst:.6f} exceeds bound {bound:.6f} at d={d} (ratio {worst / bound:.2f}x)"


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
@pytest.mark.parametrize("rho_over_r", [1.01, 1.1, 1.3, 2.0, 4.0])
def test_bound_holds_against_a_concave_void_larger_than_the_tool(d, rho_over_r):
    """rho > r: the rim can sit INSIDE the void and emerge.

    This is the regime the audit falsified. At internal tangency s0 = rho - r the
    emerged run half-angle satisfies

        psi**2 = (2*d*rho + d**2) / (r * (rho - r + d))

    so growth ~ 2*sqrt(2*d*rho / (r*(rho - r))), which exceeds the current
    4*asin(d/2r) + 2*acos(1 - d/r) ~ 2*sqrt(2*d/r) by a factor
    sqrt(rho/(rho - r)) -- without limit as rho -> r+, and above 1 for EVERY
    rho > r, so no void larger than the tool is safe, only less unsafe.
    """
    rho = rho_over_r * TOOL_RADIUS
    stock = _void_stock(rho)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)

    # Configuration pin: the base station really is the internal tangency (cutter
    # wholly inside the void, nothing engaged) and material really is reachable
    # one displacement away -- otherwise the growth below would measure nothing.
    assert _max_run_tea(stock, rho - TOOL_RADIUS, 0.0) == 0.0
    assert _max_run_tea(stock, rho - TOOL_RADIUS + d, 0.0) > 0.0

    worst = _worst_growth_over(stock, rho - TOOL_RADIUS, 0.0, d)
    assert worst <= bound, f"concave void rho/r={rho_over_r}: growth {worst:.6f} exceeds bound {bound:.6f} at d={d} (ratio {worst / bound:.2f}x)"


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
def test_bound_holds_on_stock_the_shipped_audit_actually_produces(d):
    """Voids cut by the SAME radius as the query, probed on the pocket's flank.

    Every ``Stock2::subtract_*`` removes a union of disks of the subtraction
    radius, so within ``audit_toolpath_engagement`` the void arcs always have
    rho == r. This test pins the behaviour on the flank of such a pocket, where
    the probes stand on the boundary rather than inside the void, so a future
    change to the audit's wiring that moves this geometry fails here rather than
    in a certificate.

    It does NOT certify the regime: equal radii are the singular limit of the
    concave family, and ``test_bound_holds_against_a_void_the_size_of_the_tool``
    shows the same rho == r stock falsifying the bound by an order of magnitude
    once the probe sits INSIDE the void instead of on its wall. The green here is
    about where these probes are, not about rho == r being safe.
    """
    stock = _sweep_stock()
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)

    # Configuration pin: the stock really was cut (an uncut block engages the whole
    # rim everywhere and could never grow) and the probes really stand on the
    # pocket flank, partially engaged rather than buried or clear.
    flank_x = SWEEP_FIRST_X + TOOL_RADIUS
    assert 0.0 < _max_run_tea(stock, flank_x, SWEEP_Y1) < FULL_TURN

    worst = 0.0
    for k in range(SWEEP_PASSES):
        for y in SWEEP_PROBE_HEIGHTS:
            worst = max(worst, _worst_growth_over(stock, SWEEP_FIRST_X + k * SWEEP_STEPOVER + TOOL_RADIUS, y, d))
    assert worst > 0.0, f"depleted stock: no probe grew at all at d={d} -- the harness measured nothing"
    assert worst <= bound, f"depleted stock: growth {worst:.6f} exceeds bound {bound:.6f} at d={d}"
