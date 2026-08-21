"""Falsification harness for the certificate's interior bound.

``swept_run_bound(stock, cx, cy, r, d)`` is what now stands between
``certify_segment_tea`` and an unsound certificate. It claims to bound the
largest contiguous engaged run that ANY cutter centre within ``d`` of
``(cx, cy)`` can see; ``certify_recursive`` evaluates it at both stations of a
span at half-station-spacing ``d`` and, since every centre on the span lies
within ``d`` of the nearer station, concludes that a bound under the cap
certifies the whole span. If the bound under-estimates the true reachable run
anywhere, that conclusion is false.

**This file used to measure ``tea_growth_bound`` and it FALSIFIED it** -- 18
certificate-critical violations by ratios up to 24x, on concave features
(``rho > r``) and on the tool sitting in the hole it just cut (``rho == r``). The
false certificates that followed are committed as
`tests/test_false_certificate.py`. The measurements are kept below, verbatim
where they still apply, but they are re-pointed at the bound that now carries the
obligation. ``tea_growth_bound`` still exists and is still exported -- deleting it
would delete the record of why the certifier changed shape -- but it certifies
nothing, so it has no soundness obligation and this file no longer asserts one
against it.

THE OBLIGATION CHANGED SHAPE WITH THE BOUND, and the difference matters. The old
lemma bounded GROWTH, a DIFFERENCE of two engagements, so the harness measured
differences. The new bound bounds the RUN ITSELF over a whole disk of centres, so
the harness measures the ABSOLUTE largest run reachable within ``d`` -- base
station included, since ``|S - S| = 0 <= d``. That is a strictly stronger thing to
ask of a stock configuration, and it is the thing the certificate actually needs.

What a verdict means here, precisely:

* A **red** is unconditional evidence: some cutter centre within exactly ``d`` of
  the probe sees a larger engaged run than the bound permits, so the bound is not
  an upper bound and the certificate's PROOF does not close.
* A **green** is bounded evidence: no violation was found among the probed
  centres and the sampled displacement directions. Sampling can only ever be a
  LOWER bound on the worst reachable run, so each green test additionally pins the
  configuration it claims to measure (the probe really straddles the feature, the
  measured run really is non-zero) -- a green that came from an empty arc list or
  an uncut stock would be worthless.

A BOUND THAT ALWAYS SATURATED TO A FULL TURN would pass every soundness
assertion in this file and be worthless, so `test_bound_tracks_the_half_plane_closed_form`
pins it from ABOVE against a closed form. Several probes below legitimately DO
saturate -- a void smaller than the tool leaves material wrapping the station in
every direction, and a full turn is then the true answer, not a capitulation --
which is exactly why the anti-capitulation pin is stated separately on geometry
whose answer is small.

Scale invariance: ``swept_run_bound`` depends on its arguments only through
``d / r`` and the stock's shape relative to ``r``, and the true reachable run of
every feature modelled here depends only on ``d / r`` and ``rho / r``. One tool
radius with parametrised ratios therefore covers the whole family; a second radius
would add no information.

One test here is NOT about the bound. Measuring engagement surfaced a separate,
independent defect in the kernel's engagement harvest -- stations at which
``engagement_at`` reports an engaged run larger than a full turn, which is
geometrically impossible -- so `test_reported_engagement_never_exceeds_a_full_turn`
pins it and `_max_run_tea` refuses such a reading at the point of measurement.
Without that split, a bound comparison computed from an impossible engagement
would be read as a bound failure it is not. That defect is untouched by the
interior-bound repair and those reds are still red.
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
# measure the TRUE runs. The interior bound is a statement about the true runs;
# the pessimistic (gap-closed) runs are the station test's merge repair and are
# measured elsewhere. Mixing the two would hide a bound failure behind the repair.
GAP_CLOSE_NONE = 0.0

# Directions sampled on the ring of radius d around a probe centre. 360 == 1 deg
# spacing AND divisible by 4, so the four cardinal directions -- where the
# extremal engagement of every feature in this module lies (radially out of a
# void, straight into a wall) -- are sampled EXACTLY rather than merely approached.
PROBE_DIRECTIONS = 360

# A full turn of engaged rim (radians): the reading when the whole cutter rim is
# in material, and the value `swept_run_bound` saturates to when it cannot bound
# the swept region below one.
FULL_TURN = 2.0 * math.pi

# Slack for comparing a measured TEA against its closed-form value. Reported spans
# are doubles assembled by atan2 from exact one-root arc endpoints; the kernel's
# own zone-vs-overlay study bounds that representation artefact at <= 1e-15 rad
# (docs/superpowers/state/engagement-zone-divergence.md). 1e-12 rad == 6e-11 deg
# keeps three decades of headroom over the artefact while still catching any real
# geometric error, which is orders of magnitude larger. Used ONLY on the
# configuration pins, never on the bound comparison -- that one is measurement
# against the shipped bound and takes no tolerance.
TEA_REPORTING_SLACK = 1e-12

# Ceiling on how far `swept_run_bound` may exceed the closed-form worst reachable
# run, in radians PER UNIT of probe travel d. It is dominated by the bound's own
# angular-transfer term 2*asin(d / (r - d)) ~ 2*d/r = 4*d at r = 0.5; the
# station-frame corner read-off of a component's angular extent adds well under
# d again. Measured maximum over every half-plane row below: 3.90. 8.0 is 2x that
# -- loose enough not to be a knife edge, tight enough that a bound which merely
# doubled, or one that saturated to a full turn, fails here.
HALF_PLANE_SLACK_PER_TRAVEL = 8.0

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
# clear side (TEA = 0), 0.75 == 1.5r fully clear.
WALL_PROBE_HEIGHTS = (0.0, 0.1, 0.25, 0.5, 0.75)

# Wall heights at which the bound is pinned from ABOVE against the closed form.
# 0.0 is excluded on purpose and NOT as an inconvenience: a centre exactly on the
# wall has material on the station's own row, so the swept material's bounding box
# straddles the station and the bound saturates by construction (engagement_2.cpp,
# component_angular_bound). That configuration is already pi-engaged, which no cap
# in (0, pi] can certify anyway, so the saturation costs nothing -- but it is not
# a tightness measurement and must not be read as one.
WALL_TIGHTNESS_HEIGHTS = (0.1, 0.25, 0.5, 0.75)

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

# The witness that retired the growth lemma, in its shortest form: an annular RIB
# whose centreline circle has the TOOL's own radius, thickness 0.008 r. Stated
# here as well as in `tests/test_false_certificate.py` -- there it drives the
# certificate-level verdict, here it pins the single geometric fact the repair
# turns on, namely that a rib wrapping the rim saturates the bound. The two must
# be able to fail independently.
RIB_THICKNESS = 0.004
RIB_FACETS = 256

# Travels at which the rib's fully immersed centre IS reachable from the station
# at ``(d, 0)`` -- the certifier's own root geometry, since a span centred on the
# rib centre puts both stations exactly one half-spacing from it. Floored at
# 0.003: below ``RIB_THICKNESS / 2 = 0.002`` the station itself is buried in the
# rib and reads a full turn, which measures the rib's thickness rather than the
# bound. 0.1 down to 0.0125 is the certifier's whole refinement ladder for the
# witness motions.
RIB_REACHING_TRAVELS = (1e-1, 5e-2, 2.5e-2, 1.25e-2, 3e-3)

# Station offset for the rib's anti-capitulation pin, in model units: the same
# half-spacing the un-refined witness motion uses.
RIB_STATION_OFFSET = 0.0125

# Half-width of the rib's cap-violating window, in model units: the exact oracle
# reports the 90 deg cap exceeded across ``|x| <= 0.002823``
# (`tests/test_false_certificate.py`). A travel that cannot bring the centre
# inside it cannot reach a violation.
RIB_IMMERSION_HALF_WIDTH = 0.002823

# Travels too short to reach that window from `RIB_STATION_OFFSET`, and the cap
# the bound must stay under there: 90 deg, the canonical roughing cap the
# certificate-level witnesses use. Measured bounds 0.658 / 0.490 / 0.339 rad
# against a 1.571 rad cap -- a 2.4x margin at the worst of them.
RIB_UNREACHING_TRAVELS = (2e-3, 1e-3, 1e-4)
RIB_UNREACHED_CAP = 0.5 * math.pi

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
    # `test_reported_engagement_never_exceeds_a_full_turn` pins, and every bound
    # comparison computed from it is meaningless. Failing HERE, naming the station,
    # stops that reading from surfacing downstream as a bound failure it is not.
    assert max_run <= FULL_TURN + TEA_REPORTING_SLACK, (
        f"engagement_at reported max_run_tea = {max_run:.9f} rad > 2*pi at station ({x!r}, {y!r}), r={r!r}: "
        f"a contiguous run cannot exceed a full turn -- kernel harvest defect, not an interior-bound failure"
    )
    return float(max_run)


def _swept_bound(stock: Stock, x: float, y: float, d: float, r: float = TOOL_RADIUS) -> float:
    """The SHIPPED interior bound at one station, never a local copy of the formula.

    Args:
        stock: The stock region to bound against.
        x: X coordinate of the station.
        y: Y coordinate of the station.
        d: Centre-travel radius the bound must cover.
        r: Tool radius.

    Returns:
        Upper bound (radians) on the largest engaged run at any centre within
        ``d`` of ``(x, y)``.
    """
    return float(_stock_2.swept_run_bound(stock.raw, x, y, r, d))


def _worst_run_within(stock: Stock, x0: float, y0: float, d: float, samples: int = PROBE_DIRECTIONS) -> float:
    """Largest engaged run measured at any probed centre within *d* of ``(x0, y0)``.

    The quantity ``swept_run_bound(stock, x0, y0, r, d)`` claims to bound. The base
    station is included (``|S - S| = 0 <= d``, so the bound covers it too) and the
    ring of displacement directions is sampled, which gives a LOWER bound on the
    true worst reachable run -- so a returned value above the bound falsifies it
    outright, while a value below it only means no violation was found at this
    angular resolution.

    Args:
        stock: The stock region to measure against.
        x0: X coordinate of the base station.
        y0: Y coordinate of the base station.
        d: Centre-travel distance (the displacement magnitude).
        samples: Number of equally spaced displacement directions.

    Returns:
        The largest observed engaged run, in radians.
    """
    worst = _max_run_tea(stock, x0, y0)
    for i in range(samples):
        a = 2.0 * math.pi * i / samples
        worst = max(worst, _max_run_tea(stock, x0 + d * math.cos(a), y0 + d * math.sin(a)))
    return worst


def _half_plane_worst_run(h: float, d: float) -> float:
    """Closed-form largest engaged run reachable from centre height *h* within *d*.

    Against material ``y <= 0`` the engaged run at centre height ``s`` is
    ``2*acos(s / r)``, monotone decreasing in ``s``, so over a disk of centres of
    radius ``d`` the largest run is attained at the lowest reachable centre
    ``h - d`` -- straight at the wall.

    Args:
        h: Height of the base station above the wall, in model units.
        d: Centre-travel radius.

    Returns:
        The exact largest reachable engaged run, in radians.
    """
    deepest = h - d
    if deepest >= TOOL_RADIUS:
        return 0.0
    if deepest <= -TOOL_RADIUS:
        return FULL_TURN
    return 2.0 * math.acos(deepest / TOOL_RADIUS)


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


def _ngon(radius: float, sides: int) -> Polygon:
    """Regular polygon of *sides* vertices inscribed in a circle about the origin.

    Args:
        radius: Circumradius in model units.
        sides: Number of vertices.

    Returns:
        A closed `Polygon` in the world XY plane, wound counterclockwise.
    """
    return Polygon([(radius * math.cos(2.0 * math.pi * i / sides), radius * math.sin(2.0 * math.pi * i / sides), 0.0) for i in range(sides)])


def _rib_stock() -> Stock:
    """Annular rib of `RIB_THICKNESS` whose centreline circle has radius `TOOL_RADIUS`.

    Returns:
        A `Stock` whose material is the annulus between radii
        ``TOOL_RADIUS -/+ RIB_THICKNESS / 2`` -- what two concentric passes leave
        when the step-over overshoots.
    """
    return Stock(
        _ngon(TOOL_RADIUS + 0.5 * RIB_THICKNESS, RIB_FACETS),
        holes=[_ngon(TOOL_RADIUS - 0.5 * RIB_THICKNESS, RIB_FACETS)],
    )


def _small_void_probe_offsets(rho: float) -> tuple[float, ...]:
    """Probe centre offsets along +x for a void strictly smaller than the tool.

    Covers the complete contact life of the rim against such a void: concentric,
    centre on the void rim, rim entering the void, rim at its DEEPEST inside it,
    rim externally tangent, and fully clear. The two interior offsets
    (``r - rho`` and ``r``) are what make the case non-vacuous: at every other
    offset the void is either swallowed whole by the cutter disk or out of
    contact, both of which read a full turn.

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


# --------------------------------------------------------------------------- #
# The bound in isolation, against geometry whose answer is a closed form        #
# --------------------------------------------------------------------------- #


def test_bound_is_exactly_zero_when_no_material_lies_in_the_swept_annulus():
    """No material within reach of any displaced rim -> the bound is EXACTLY zero.

    Not "small": zero. Every rim point of every reachable centre lies in the
    annulus ``A(S, r - d, r + d)``, so an annulus empty of material means no
    reachable centre touches anything, and the bound reports that outright rather
    than falling back to its angular-transfer term (which measures the spread of a
    run that does not exist).
    """
    stock = Stock(_HALF_PLANE)
    d = 1e-3
    clear = TOOL_RADIUS + 1.5 * d  # rim of every reachable centre stays above the wall

    # Configuration pin: the wall really is out of reach here, and really is in
    # reach a little lower -- so the zero below is a measurement, not an empty stock.
    assert _worst_run_within(stock, 0.0, clear, d) == 0.0
    assert _max_run_tea(stock, 0.0, TOOL_RADIUS - d) > 0.0

    assert _swept_bound(stock, 0.0, clear, d) == 0.0


def test_bound_saturates_when_the_station_is_buried_in_material():
    """Cutter buried in virgin stock -> a full turn, which is also the truth there.

    The rim is wholly in material, so the largest run IS ``2*pi`` and the bound
    cannot be smaller. It is reached through the wrap branch: the swept material
    surrounds the station, so its bounding box holds the station and the bound
    saturates (engagement_2.cpp, component_angular_bound).
    """
    stock = Stock(_BLOCK)
    d = 1e-2

    # Configuration pin: the reachable runs really are a full turn, so a saturated
    # bound here is tight rather than a capitulation.
    assert _worst_run_within(stock, 0.0, 0.0, d) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)

    assert _swept_bound(stock, 0.0, 0.0, d) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
@pytest.mark.parametrize("h", WALL_TIGHTNESS_HEIGHTS)
def test_bound_tracks_the_half_plane_closed_form(h, d):
    """ANTI-CAPITULATION PIN: the bound is squeezed between the closed form and it plus O(d).

    Every other assertion in this file is one-sided -- ``worst <= bound`` -- and a
    bound that returned ``2*pi`` unconditionally would satisfy all of them while
    certifying nothing. This test is the other side, stated against geometry whose
    answer is exact: over a half-plane the largest run reachable from height ``h``
    within ``d`` is ``2*acos((h - d) / r)``, attained straight at the wall.

    The upper pin is ``closed form + HALF_PLANE_SLACK_PER_TRAVEL * d``, dominated
    by the bound's own angular-transfer term ``2*asin(d / (r - d))``. It therefore
    tightens as the certifier refines, which is the property that makes the
    certifier converge instead of refusing everything.
    """
    stock = Stock(_HALF_PLANE)
    truth = _half_plane_worst_run(h, d)
    bound = _swept_bound(stock, 0.0, h, d)

    # Configuration pin: the closed form really is what the kernel measures, so a
    # comparison against it measures the BOUND and not a modelling mistake.
    assert _worst_run_within(stock, 0.0, h, d) == pytest.approx(truth, abs=TEA_REPORTING_SLACK)

    assert bound >= truth, f"h={h}, d={d}: bound {bound:.6f} is BELOW the reachable run {truth:.6f} -- not an upper bound"
    assert bound <= truth + HALF_PLANE_SLACK_PER_TRAVEL * d, (
        f"h={h}, d={d}: bound {bound:.6f} exceeds the closed form {truth:.6f} by {bound - truth:.6f} rad, "
        f"more than {HALF_PLANE_SLACK_PER_TRAVEL} rad per unit travel -- the bound has gone slack"
    )


@pytest.mark.parametrize("d", RIB_REACHING_TRAVELS)
def test_bound_saturates_on_the_rib_that_retired_the_growth_lemma(d):
    """The single geometric fact the repair turns on: a rib wrapping the rim saturates.

    An annular rib of the tool's own radius surrounds the station, so the material
    inside the swept annulus reaches every direction from it and no angular extent
    below a full turn can be claimed. That is what makes `certify_segment_tea`
    refuse the motions of `tests/test_false_certificate.py` -- and it is what the
    retired growth lemma could not see, because the rib creates engagement that
    does not grow from anything either station measures: the station reads 0.32 rad
    at ``d = 0.0125`` while the fully immersed centre sits one half-spacing away.

    The station is placed at ``(d, 0)``, which is the certifier's OWN root geometry
    for the witness motion -- a span centred on the rib centre puts both stations
    exactly one half-spacing from it. Pinned separately from the certificate-level
    tests so the bound's behaviour and the certifier's verdict can fail
    independently.
    """
    stock = _rib_stock()

    # Configuration pin: the rib really is built and really is straddled -- the
    # station is partially engaged, and the rib centre one travel away is a full
    # immersion. Both are what make a saturated bound the TRUE answer here.
    assert 0.0 < _max_run_tea(stock, d, 0.0) < FULL_TURN - TEA_REPORTING_SLACK
    assert _max_run_tea(stock, 0.0, 0.0) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)

    assert _swept_bound(stock, d, 0.0, d) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)


@pytest.mark.parametrize("d", RIB_UNREACHING_TRAVELS)
def test_bound_does_not_refuse_a_rib_the_travel_cannot_reach(d):
    """ANTI-CAPITULATION PIN on the rib: saturating is a MEASUREMENT, not a reflex.

    The same rib, the same station, but a travel too short for the immersed centre
    to be reachable: the cap-violating window is ``|x| <= 0.002823`` and the closest
    centre a travel of at most `RIB_STATION_OFFSET` - d can reach stays outside it.
    No centre on such a span violates, and the bound says so -- it drops to well
    under the 90 deg cap the witnesses are certified against instead of refusing
    on sight of a ring.

    Without this, `test_bound_saturates_on_the_rib_that_retired_the_growth_lemma`
    would be satisfied by a bound that returned a full turn for any stock
    containing a hole, which certifies nothing and repairs nothing. It is the
    bound-level statement of the certificate-level green control
    ``test_a_motion_clear_of_the_rib_centre_is_soundly_certified``.
    """
    stock = _rib_stock()

    # Configuration pin: the immersion really is out of reach, so a bound below the
    # cap is correct here and not a missed violation.
    assert RIB_STATION_OFFSET - d > RIB_IMMERSION_HALF_WIDTH
    assert _worst_run_within(stock, RIB_STATION_OFFSET, 0.0, d) < RIB_UNREACHED_CAP

    bound = _swept_bound(stock, RIB_STATION_OFFSET, 0.0, d)
    assert bound < RIB_UNREACHED_CAP, f"bound {bound:.6f} refuses a rib no reachable centre engages past {RIB_UNREACHED_CAP:.6f} rad, at d={d}"


@pytest.mark.parametrize("d", [0.0, 0.5 * TOOL_RADIUS, TOOL_RADIUS, 2.0 * TOOL_RADIUS])
def test_bound_saturates_when_the_annulus_degenerates(d):
    """Travels the construction cannot bound report a full turn, not a comfortable number.

    ``d == 0`` leaves the annulus with empty interior -- it cannot see the rim at
    all, so ``0.0`` would be a lie -- and ``d >= r/2`` breaks the angular-transfer
    step, which needs ``d < r - d``. Both saturate, which forces the certifier to
    refine or refuse. A silent small answer at either would be a false certificate
    generator.

    The certifier never reaches the ``d == 0`` case as an interior claim: a span of
    zero length IS its single station, which the exact station predicate decides
    outright (engagement_2.cpp, interior_run_within_cap).
    """
    stock = Stock(_HALF_PLANE)

    # Configuration pin: this station is genuinely engaged, so a full turn is not
    # simply what an untouched stock would return.
    assert 0.0 < _max_run_tea(stock, 0.0, 0.25) < FULL_TURN

    assert _swept_bound(stock, 0.0, 0.25, d) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)


@pytest.mark.parametrize(
    "cx, cy, r, d",
    [
        (float("nan"), 0.0, TOOL_RADIUS, 1e-3),
        (0.0, float("inf"), TOOL_RADIUS, 1e-3),
        (0.0, 0.0, 0.0, 1e-3),
        (0.0, 0.0, -TOOL_RADIUS, 1e-3),
        (0.0, 0.0, float("inf"), 1e-3),
        (0.0, 0.0, TOOL_RADIUS, float("nan")),
        (0.0, 0.0, TOOL_RADIUS, -1e-9),
    ],
)
def test_bound_rejects_non_physical_arguments(cx, cy, r, d):
    """The seam contract: refuse before any exact injection, never guess.

    A non-finite coordinate is not a rational and has no exact image at the seam, a
    non-positive radius is not a cutter, and a negative travel is not a shorter
    motion but a nonsensical one. Each is refused by name rather than folded into a
    confident number -- the defect the sibling ``_sign_mixed_radical`` binding
    records, where an unguarded caller got a confident sign for ``sqrt(-1)``.
    """
    stock = Stock(_BLOCK)
    with pytest.raises(ValueError):
        _stock_2.swept_run_bound(stock.raw, cx, cy, r, d)


# --------------------------------------------------------------------------- #
# A separate, independent kernel defect -- NOT about the interior bound         #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("station", FULL_TURN_WITNESS_STATIONS)
def test_reported_engagement_never_exceeds_a_full_turn(station):
    """A reported engaged run larger than a full turn is impossible -- and the kernel reports one.

    NOT a bound test. `engagement_at`'s ``total_tea`` and ``max_run_tea`` are the
    REPORTING half of the deciding/reporting split (engagement_2.cpp, step 4a):
    doubles summed from ``atan2`` spans of exact one-root arc endpoints. Their
    ceiling is structural -- the rim is a closed curve of angular measure ``2*pi``,
    so no contiguous run on it can measure more, whatever the stock looks like.

    At the stations below the kernel reports the true engagement plus exactly one
    full turn. Surfaced while measuring engagement against a void smaller than the
    tool (a growth of 6.62 rad is impossible on its face), so it is pinned
    separately here rather than left to masquerade as a bound failure. Consumers
    inherit it: the audit's reported ``max_tea`` is this number.

    Untouched by the interior-bound repair, which changes what the certifier
    CONCLUDES and not how the harvest assembles a span. The cap DECISION is not
    implicated by this evidence either -- it runs on exact predicates over the run
    ENDPOINTS (step 4b), never on these doubles -- but a harvest that mis-assembles
    a run's span is not obviously assembling its endpoints correctly, and nothing
    here establishes that it does.
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


# --------------------------------------------------------------------------- #
# The bound against the boundary shapes the stock model actually produces       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
def test_bound_holds_against_a_straight_wall(d):
    """The half-plane case: the simplest feature, and the one every pocket flank reduces to.

    A straight wall is where the interior bound must be TIGHT rather than merely
    true -- it is the shape ordinary cutting presents to the tool, so a bound that
    went slack here would refuse the motions the certifier exists to pass.
    `test_bound_tracks_the_half_plane_closed_form` pins that tightness; this test
    pins soundness over the whole probe set with the ring actually swept.
    """
    stock = Stock(_HALF_PLANE)

    # Configuration pin: the probe set really straddles the wall. The deepest
    # probe sits ON it (half the rim in material) and the clear probe is out of
    # contact -- so a zero below could only be a real measurement.
    engagements = [_max_run_tea(stock, 0.0, y) for y in WALL_PROBE_HEIGHTS]
    assert max(engagements) == pytest.approx(math.pi, abs=TEA_REPORTING_SLACK)
    assert min(engagements) == 0.0

    for y in WALL_PROBE_HEIGHTS:
        worst = _worst_run_within(stock, 0.0, y, d)
        bound = _swept_bound(stock, 0.0, y, d)
        assert worst <= bound, f"straight wall y={y}: reachable run {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
@pytest.mark.parametrize("rho", [0.1, 0.25])
def test_bound_holds_against_a_void_smaller_than_the_tool(d, rho):
    """rho < r: the rim can never sit inside the void, so material wraps the station.

    A void strictly smaller than the tool cannot contain the cutter: it merely
    takes a bite out of an otherwise fully engaged rim. Material therefore reaches
    the station from every direction and the bound saturates -- correctly, because
    the reachable run genuinely IS a full turn at five of the six offsets probed.
    A saturated bound here is a tight one, not a capitulation.

    ``rho == d`` is the one case that does not reach its bound comparison: the
    ``s = r`` probe's displacement ring then puts the rim's x-extreme exactly
    ``rho`` from the void centre in every direction, and some of those stations
    trip the impossible-engagement defect
    (`test_reported_engagement_never_exceeds_a_full_turn`). `_max_run_tea` refuses
    the reading, so that parametrisation reports the harvest defect by name instead
    of a bound comparison derived from it. THAT RED IS NOT AN INTERIOR-BOUND
    FAILURE and the repair does not touch it. The probe stays: it is the offset at
    which the rim cuts deepest into the void, and dropping it would make this test
    vacuous.
    """
    assert rho < TOOL_RADIUS  # the regime this test claims to cover
    stock = _void_stock(rho)

    # Configuration pin: the probe set really straddles the void wall. At s = r the
    # rim passes through the void centre and is at its deepest inside it, so the
    # engaged run must fall short of a full turn; past external tangency the void
    # is out of contact and the run must be a full turn exactly.
    assert _max_run_tea(stock, TOOL_RADIUS, 0.0) < FULL_TURN - TEA_REPORTING_SLACK
    assert _max_run_tea(stock, rho + TOOL_RADIUS + CLEAR_MARGIN, 0.0) == pytest.approx(FULL_TURN, abs=TEA_REPORTING_SLACK)

    for s in _small_void_probe_offsets(rho):
        worst = _worst_run_within(stock, s, 0.0, d)
        bound = _swept_bound(stock, s, 0.0, d)
        assert worst <= bound, f"rho={rho}, s={s}: reachable run {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
def test_bound_holds_against_a_void_the_size_of_the_tool(d):
    """rho == r: the tool sitting in the hole it just cut -- where the growth lemma failed hardest.

    This is a plunge followed by a departure, the most ordinary motion pair a
    toolpath contains: ``subtract_disk(c, r)`` then a cut starting at ``c``. The
    cutter is internally tangent to its own hole, so the base station reads
    TEA = 0, and the emerged half-angle at travel ``delta`` satisfies
    ``cos(psi) = -delta / (2r)``: the run appears at just over ``pi`` for ANY
    positive travel, however small. Growth is therefore discontinuous at zero and
    NO growth lemma can bound it -- ``tea_growth_bound`` was falsified here by an
    order of magnitude.

    The interior bound is not troubled by it, because it never asks how a run
    grows: material fills the annulus all the way round the station, so the bound
    saturates at a full turn and the certifier refines instead of concluding.
    """
    rho = TOOL_RADIUS
    stock = _void_stock(rho)

    # Configuration pin: the cutter really is sitting in its own hole (rim exactly
    # on the void wall, nothing engaged) and there really is material to emerge
    # into one displacement away -- so the discontinuity is live.
    assert _max_run_tea(stock, 0.0, 0.0) == 0.0
    assert _max_run_tea(stock, d, 0.0) > math.pi

    worst = _worst_run_within(stock, 0.0, 0.0, d)
    bound = _swept_bound(stock, 0.0, 0.0, d)
    assert worst <= bound, f"void rho=r: reachable run {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
@pytest.mark.parametrize("rho_over_r", [1.01, 1.1, 1.3, 2.0, 4.0])
def test_bound_holds_against_a_concave_void_larger_than_the_tool(d, rho_over_r):
    """rho > r: the regime that falsified the growth lemma, at every ratio it falsified.

    At internal tangency ``s0 = rho - r`` the emerged run half-angle satisfies
    ``psi**2 = (2*d*rho + d**2) / (r * (rho - r + d))``, so the run grows like
    ``2*sqrt(2*d*rho / (r*(rho - r)))`` -- exceeding
    ``4*asin(d/2r) + 2*acos(1 - d/r)`` by ``sqrt(rho/(rho - r))``, without limit as
    ``rho -> r+`` and above 1 for EVERY ``rho > r``. No void larger than the tool
    was safe under the old lemma, only less unsafe.

    The interior bound covers the same configurations by construction: the emerged
    run lies inside the swept annulus whether or not it grew from anything, so
    bounding what the annulus can hold bounds it. Measured slack at
    ``rho/r = 4, d = 1e-4``: 0.9%, where the old lemma under-estimated by 14%.
    """
    rho = rho_over_r * TOOL_RADIUS
    stock = _void_stock(rho)

    # Configuration pin: the base station really is the internal tangency (cutter
    # wholly inside the void, nothing engaged) and material really is reachable
    # one displacement away -- otherwise the comparison below would measure nothing.
    assert _max_run_tea(stock, rho - TOOL_RADIUS, 0.0) == 0.0
    assert _max_run_tea(stock, rho - TOOL_RADIUS + d, 0.0) > 0.0

    worst = _worst_run_within(stock, rho - TOOL_RADIUS, 0.0, d)
    bound = _swept_bound(stock, rho - TOOL_RADIUS, 0.0, d)
    assert worst <= bound, f"concave void rho/r={rho_over_r}: reachable run {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
def test_bound_holds_on_stock_the_shipped_audit_actually_produces(d):
    """Voids cut by the SAME radius as the query, probed on the pocket's flank.

    Every ``Stock2::subtract_*`` removes a union of disks of the subtraction
    radius, so within ``audit_toolpath_engagement`` the void arcs always have
    ``rho == r`` and the boundary is the scalloped mixture of same-radius arcs and
    straight flanks a real roughing pass leaves. This is the regime the certifier
    ships into, so a future change to the audit's wiring that moves this geometry
    fails here rather than in a certificate.
    """
    stock = _sweep_stock()

    # Configuration pin: the stock really was cut (an uncut block engages the whole
    # rim everywhere) and the probes really stand on the pocket flank, partially
    # engaged rather than buried or clear.
    flank_x = SWEEP_FIRST_X + TOOL_RADIUS
    assert 0.0 < _max_run_tea(stock, flank_x, SWEEP_Y1) < FULL_TURN

    for k in range(SWEEP_PASSES):
        for y in SWEEP_PROBE_HEIGHTS:
            x = SWEEP_FIRST_X + k * SWEEP_STEPOVER + TOOL_RADIUS
            worst = _worst_run_within(stock, x, y, d)
            bound = _swept_bound(stock, x, y, d)
            assert worst <= bound, f"depleted stock ({x}, {y}): reachable run {worst:.6f} exceeds bound {bound:.6f} at d={d}"
