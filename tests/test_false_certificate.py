"""Certificate-level search for a FALSE TEA certificate -- and the witness it found.

`certify_segment_tea` states an unconditional contract (`src/engagement_2.h`):
"Certify TEA(P) <= cap_radians for EVERY cutter center P on the segment", with
``cap_certified`` "true iff no cutter center on the motion can exceed the cap".
`tests/test_growth_bound.py` showed the analytic guard that contract rests on is
not an upper bound, so the certificate's PROOF does not close. That left the only
question that matters open: is any VERDICT actually wrong?

It is. The witness below is a motion `certify_segment_tea` CERTIFIES at
cap = 90 deg while the cutter is, part-way along it, immersed a FULL TURN
(``max_run_tea`` = 2*pi, `cap_exceeded` exactly true at 91 of 401 probed centres).
The certifier never sees it: it reports ``max_tea`` = 0.322 rad, twenty times
under the truth.

THE CONSTRUCTION -- an annular rib whose mean radius is the tool radius.

Two boundary features, one inside the other: a circular void of radius
``r - tau/2`` and, outside it, material removed again beyond ``r + tau/2``,
leaving a thin annular RIB of thickness ``tau`` whose centreline is a circle of
the TOOL's own radius. This is what two concentric passes leave behind when the
step-over overshoots by ``tau`` -- `test_machined_stock_certified_motion_has_no_cap_violating_centre`
builds exactly that, with `subtract_disk` + `subtract_arc_sweep` and nothing else.
The rib is built three independent ways below -- through `Stock`'s boundary/hole
constructor, by material removal alone, and as one plain simple polygon -- so no
verdict here can be blamed on a constructor.

With the cutter concentric with the rib its whole rim lies in material, so
TEA = 2*pi. Move the centre off by ``a`` and the rim cuts across the rib instead,
leaving two SHORT arcs of roughly ``tau / a`` radians each -- vanishing as the rib
thins. So TEA along a motion through the rib centre is small, then a full turn,
then small again: the interior maximum is NOT at an endpoint, and the certifier
measures only endpoints.

WHY THE ENDPOINT-ATTAINMENT ARGUMENT DOES NOT COVER THIS. The single-void argument
(TEA monotone in the distance ``s`` to the void centre, ``s`` convex along a line,
hence the maximum is at an endpoint) needs TEA to be a function of ONE distance.
Here it is a function of two -- distance to the inner wall and to the outer wall --
and the engaged arc is the region BETWEEN them, which is maximal in the middle.
Two features are enough to break it, and the break is not marginal: the interior
value is the largest TEA that exists.

THE ORACLE IS EXACT, so a red here cannot be a sampling artifact. Each probed
centre is judged by `engagement_at`'s ``cap_exceeded``, the same exact
orientation + squared-chord predicate the certifier's own stations use
(engagement_2.cpp, step 4b) -- never a comparison of reported doubles. The
certifier's claim is UNIVERSAL over the segment, so ONE exactly-violating centre
refutes it outright. Sampling density therefore governs only whether a violation
is FOUND, never whether a found one is real: refining the sample grid can turn a
green into a red, never a red into a false alarm. ``max_run_tea`` appears below
only to say HOW BADLY, never as evidence -- it is a REPORTING double subject to
the known one-full-turn harvest defect
(`test_growth_bound.py::test_reported_engagement_never_exceeds_a_full_turn`).

HOW FAR IT GOES. Open the rib into a 135 deg sector and the flanking stations stop
touching it at all: the certifier then returns ``max_tea = 0.0`` -- "the cutter
never contacted material anywhere on this motion" -- for a motion on which it is,
at one point, 135 deg engaged. A guard added to a growth bound cannot repair a
verdict drawn from two measurements that are both identically zero, which is why
`test_certified_motion_whose_stations_report_no_contact_at_all_has_no_cap_violating_centre`
constrains the repair more tightly than the rest.

TWO GREEN CONTROLS keep the reds from being vacuous: the certifier is SOUND on a
motion of this same stock that never violates, and it correctly REFUSES a motion
whose refinement happens to land a station on the violation. The defect is a blind
spot between stations, not a certifier that says yes to everything.
"""

from __future__ import annotations

import math

from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

# Tool radius r, in model units, shared with `tests/test_growth_bound.py`. The
# whole construction is scale-free -- rib thickness, motion length and the guard
# all depend on their arguments only through a ratio to r -- and the witness was
# measured identical (certified, interior fully immersed, one station pair) at
# r = 0.05, 0.1, 0.5, 2.0 and 10.0, so one radius covers the family.
TOOL_RADIUS = 0.5

# Engagement cap for every certificate below: 90 deg, the canonical roughing cap
# and the one `tests/test_growth_bound.py` measures against. Not a special value:
# the witness was reproduced for caps from 0.05 rad to pi by scaling the rib
# thickness with the station spacing the cap forces (see the module report).
CAP_RADIANS = math.pi / 2.0

# Exact squared-chord surrogate 4*sin^2(cap/2) for CAP_RADIANS, the form
# `engagement_at` takes (engagement_2.h, boundary doctrine). This is the FULL cap
# -- the oracle asks whether the TRUE cap is exceeded, never a guarded one.
CAP_CHORD_RATIO = 4.0 * math.sin(0.5 * CAP_RADIANS) ** 2

# Gap-closure surrogate for gamma = 0: the oracle absorbs no void gap, so it reads
# the TRUE runs. Gap-closure pessimism is the certifier's internal merge repair and
# inflating the oracle with it would understate the violation, not overstate it.
GAP_CLOSE_NONE = 0.0

# Rib thickness tau, in model units (0.008 r). Ceiling: a station at half-spacing
# `a` reads two engaged arcs of about tau / a radians, and it passes only while
# that stays under the guarded cap, so tau <~ a * (cap - tea_guard(a, r)) =
# 0.0125 * 0.5745 = 0.00718 for the short motion below. 0.004 is 56% of that
# ceiling -- an ~1.8x margin, so the witness is a region of the parameter space
# and not a knife edge. Measured: reds over tau/r in [0.001, 0.008] against
# motion lengths from 0.0125 r to 0.1 r.
RIB_THICKNESS = 0.004

# Sides of the regular polygons approximating the rib walls. The faceting error
# r * (1 - cos(pi / n)) = 3.76e-5 is 0.94% of RIB_THICKNESS, so the rib is uniform
# to within 1% and every margin claimed here clears the faceting by 50x or more.
RIB_FACETS = 256

# Angular extent of the OPEN rib -- an annular sector rather than a closed ring,
# because the rib need not close for the certificate to fail. 135 deg is 1.5x the
# cap, so the sector's own extent IS the engaged run at its centre and the
# violation is read straight off the geometry. Measured: reds down to 93.6 deg,
# i.e. the moment the sector outruns the cap; 72 deg is soundly certified.
SECTOR_EXTENT = 0.75 * math.pi

# Vertices per sector wall, set so the angular step matches the closed ring's
# 2*pi / RIB_FACETS and therefore its 0.94% faceting error exactly.
SECTOR_FACETS = round(RIB_FACETS * SECTOR_EXTENT / (2.0 * math.pi))

# Short motion, in model units (0.05 r), centred on the rib centre. Two bounds fix
# it: it must be SHORT enough that the certifier can certify without refining at
# all (the guarded cap is positive only below 0.0561 = 0.112 r at this cap), and
# LONG enough that the cap-violating window -- measured at |x| <= 0.002823, i.e.
# 22.6% of this motion -- falls strictly between the two stations. 0.025 is 45% of
# the certifiable ceiling, comfortably inside both.
SHORT_MOTION = 0.025

# A long, obliquely-oriented motion through the rib centre, 0.75 units = 1.5 tool
# diameters, with a direction (0.8, 0.6) chosen so it is not axis-aligned. Its job
# is to show the defect survives ADAPTIVE REFINEMENT: the certifier bisects it to
# 35 stations, settling at a leaf spacing of 0.75/16 = 0.046875 -- the coarsest
# spacing at which any guarded cap exists -- which is 8.3x the 0.005646-wide
# violating window, so the window still fits between two stations.
LONG_MOTION_START = (-0.4, -0.3)
LONG_MOTION_END = (0.2, 0.15)

# Alignment sweep: the same 0.6-long crossing motion slid along its own direction
# through one whole leaf spacing and a bit more, to show the red is not an artifact
# of where the dyadic stations happen to fall. 33 offsets over +/- 0.15 samples the
# 0.0375 leaf spacing four times over.
SWEEP_MOTION = 0.6
SWEEP_OFFSETS = 33
SWEEP_OFFSET_SPAN = 0.15

# Probes along a motion for the exact oracle. The violating window spans 22.6% of
# SHORT_MOTION, so any grid of a handful of points already lands in it; 400 is
# chosen for a legible count (91 hits) rather than for detection, and the oracle's
# exactness means a denser grid could only ever raise that count.
SCAN_SAMPLES = 400

# Ambient block half-width for the machined build, in model units (12 r): the
# block's straight walls stay >= 5 units clear of the rib, so they never reach the
# cutter rim and the measurement sees the rib alone.
BLOCK_HALF_WIDTH = 12.0 * TOOL_RADIUS

_BLOCK = Polygon(
    [
        (-BLOCK_HALF_WIDTH, -BLOCK_HALF_WIDTH, 0),
        (BLOCK_HALF_WIDTH, -BLOCK_HALF_WIDTH, 0),
        (BLOCK_HALF_WIDTH, BLOCK_HALF_WIDTH, 0),
        (-BLOCK_HALF_WIDTH, BLOCK_HALF_WIDTH, 0),
    ]
)

# A full turn of engaged rim (radians): the reading when the whole cutter rim is in
# material, which is exactly what the cutter reads at the rib centre.
FULL_TURN = 2.0 * math.pi

# Slack for comparing a REPORTING span against its closed-form value, carried over
# verbatim from `tests/test_growth_bound.py`: reported spans are doubles assembled
# by atan2 from exact one-root arc endpoints and the kernel's own zone-vs-overlay
# study bounds that representation artefact at <= 1e-15 rad. Used ONLY on the
# configuration pins; no verdict below rests on a double comparison.
TEA_REPORTING_SLACK = 1e-12


def _ngon(radius: float, sides: int) -> Polygon:
    """Regular polygon of *sides* vertices inscribed in a circle about the origin.

    Args:
        radius: Circumradius in model units.
        sides: Number of vertices.

    Returns:
        A closed `Polygon` in the world XY plane, wound counterclockwise.
    """
    return Polygon([(radius * math.cos(2.0 * math.pi * i / sides), radius * math.sin(2.0 * math.pi * i / sides), 0.0) for i in range(sides)])


def _rib_stock(thickness: float = RIB_THICKNESS) -> Stock:
    """Annular rib of the given thickness whose centreline circle has radius `TOOL_RADIUS`.

    Built through the `Stock` boundary/hole constructor -- the shortest statement
    of the shape. `_machined_rib_stock` builds the same rib by removal only.

    Args:
        thickness: Radial thickness of the rib in model units.

    Returns:
        A `Stock` whose material is the annulus between radii
        ``TOOL_RADIUS -/+ thickness / 2``.
    """
    return Stock(
        _ngon(TOOL_RADIUS + 0.5 * thickness, RIB_FACETS),
        holes=[_ngon(TOOL_RADIUS - 0.5 * thickness, RIB_FACETS)],
    )


def _machined_rib_stock(thickness: float = RIB_THICKNESS) -> Stock:
    """The same rib left behind in a solid block by two removals and nothing else.

    A bore of radius ``TOOL_RADIUS - thickness / 2`` opens the inside; one full
    circular contour pass with a cutter of radius `TOOL_RADIUS`, guided at radius
    ``TOOL_RADIUS + thickness / 2 + TOOL_RADIUS``, clears the outside. What remains
    between them is the rib -- the ordinary consequence of a step-over that
    overshoots by ``thickness``. `subtract_arc_sweep` lays its guide as a chain of
    disks whose under-coverage is bounded at ``1e-4 * tool_radius`` = 5e-5, i.e.
    1.3% of the rib, so the machined rib matches the ideal one well inside every
    margin claimed here.

    Args:
        thickness: Radial thickness of the rib in model units.

    Returns:
        A `Stock` carrying the same rib, reachable by material removal alone.
    """
    stock = Stock(_BLOCK)
    stock.subtract_disk(0.0, 0.0, TOOL_RADIUS - 0.5 * thickness)
    guide = TOOL_RADIUS + 0.5 * thickness + TOOL_RADIUS
    stock.subtract_arc_sweep(0.0, 0.0, guide, 0.0, guide, 0.0, False, TOOL_RADIUS)
    return stock


def _sector_rib_stock(extent: float = SECTOR_EXTENT, thickness: float = RIB_THICKNESS) -> Stock:
    """An OPEN rib: one annular sector, stated as a single simple polygon.

    The third independent construction route -- no holes, no material removal,
    just a boundary. It is also the weakest shape that still breaks the
    certificate: an arc of leftover material of the tool's own radius, spanning a
    little more than the cap.

    Args:
        extent: Angular extent of the sector in radians, centred on the +x axis.
        thickness: Radial thickness of the rib in model units.

    Returns:
        A `Stock` whose material is the annular sector of mean radius
        `TOOL_RADIUS`, traced outer wall then inner wall.
    """
    outer, inner = TOOL_RADIUS + 0.5 * thickness, TOOL_RADIUS - 0.5 * thickness
    steps = SECTOR_FACETS
    points = [(outer * math.cos(-0.5 * extent + extent * i / steps), outer * math.sin(-0.5 * extent + extent * i / steps), 0.0) for i in range(steps + 1)]
    points += [(inner * math.cos(-0.5 * extent + extent * (steps - i) / steps), inner * math.sin(-0.5 * extent + extent * (steps - i) / steps), 0.0) for i in range(steps + 1)]
    return Stock(Polygon(points))


def _sample(stock: Stock, x: float, y: float) -> tuple[float, bool]:
    """Reported largest engaged run and the EXACT cap verdict at one cutter centre.

    Args:
        stock: The stock region to measure against (frozen: nothing is cut here).
        x: X coordinate of the cutter centre.
        y: Y coordinate of the cutter centre.

    Returns:
        ``(max_run_tea, cap_exceeded)``. The first is a REPORTING double, quoted in
        failure messages only; the second is the exact predicate that decides
        everything here.
    """
    _total, max_run, exceeded = _stock_2.engagement_at(stock.raw, x, y, TOOL_RADIUS, CAP_CHORD_RATIO, GAP_CLOSE_NONE)
    return float(max_run), bool(exceeded)


def _violating_centres(stock: Stock, start: tuple[float, float], end: tuple[float, float], samples: int = SCAN_SAMPLES) -> list[tuple[float, float]]:
    """Probed centres on a motion whose EXACT cap verdict is "exceeded".

    Every entry is a counterexample to `certify_segment_tea`'s universal claim over
    the same segment, decided by the same exact predicate the certifier's stations
    use. An empty list is bounded evidence (nothing found at this resolution); a
    non-empty one is unconditional.

    Args:
        stock: The stock region to measure against.
        start: ``(x, y)`` of the motion start.
        end: ``(x, y)`` of the motion end.
        samples: Number of sub-intervals of the motion to probe.

    Returns:
        The violating cutter centres, in order along the motion.
    """
    x0, y0 = start
    x1, y1 = end
    hits = []
    for i in range(samples + 1):
        t = i / samples
        x, y = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
        if _sample(stock, x, y)[1]:
            hits.append((x, y))
    return hits


def _certify(stock: Stock, start: tuple[float, float], end: tuple[float, float]) -> tuple[float, bool, int]:
    """Run the shipped certificate over a motion.

    Args:
        stock: The stock region to certify against.
        start: ``(x, y)`` of the motion start.
        end: ``(x, y)`` of the motion end.

    Returns:
        ``(max_tea, cap_certified, stations)`` exactly as `certify_segment_tea`
        returns it.
    """
    return _stock_2.certify_segment_tea(stock.raw, start[0], start[1], end[0], end[1], TOOL_RADIUS, CAP_RADIANS)


def _pin_rib(stock: Stock, half_spacing: float) -> None:
    """Assert the rib really is in the regime the tests below claim to measure.

    Without this a broken construction could produce a green (or a red) that means
    nothing: an empty stock reads zero everywhere, and a rib the cutter cannot
    straddle reads a full turn everywhere.

    Args:
        stock: The rib stock under test.
        half_spacing: Distance from the rib centre at which the flanking stations
            of the motion under test sit.

    Raises:
        AssertionError: If the cutter is not fully immersed at the rib centre, or
            if the stations at ``half_spacing`` are not partially engaged.
    """
    centre_run, centre_exceeded = _sample(stock, 0.0, 0.0)
    assert centre_run == FULL_TURN, f"rib centre should immerse the whole rim: max_run_tea {centre_run!r} != 2*pi"
    assert centre_exceeded, "rib centre should exceed a 90 deg cap -- the exact oracle disagrees, so the rib is not built"
    for sign in (-1.0, 1.0):
        run, exceeded = _sample(stock, sign * half_spacing, 0.0)
        assert 0.0 < run < FULL_TURN - TEA_REPORTING_SLACK, f"station at x={sign * half_spacing!r} should straddle the rib, not miss it or be buried: max_run_tea {run!r}"
        assert not exceeded, f"station at x={sign * half_spacing!r} should be under the 90 deg cap: max_run_tea {run!r}"


def test_certified_short_motion_has_no_cap_violating_centre():
    """A motion the certifier certifies WITHOUT refining at all, whose middle is fully immersed.

    The minimal witness: the segment is already short enough for a positive guarded
    cap, so `certify_segment_tea` visits exactly ONE station pair (``stations == 1``)
    and returns immediately. No refinement, no bisection, nothing adaptive -- the
    verdict is the two endpoint measurements and the analytic guard, which is the
    certificate in its purest form.

    Measured: stations at x = -/+0.0125 read max_run_tea 0.32190 against a guarded
    cap of 0.57449, so both pass; the centre reads a FULL TURN and `cap_exceeded`
    is exactly true across |x| <= 0.002823, 22.6% of the motion. The certifier
    returns ``cap_certified = True`` with ``max_tea = 0.32190``.
    """
    stock = _rib_stock()
    half = 0.5 * SHORT_MOTION
    start, end = (-half, 0.0), (half, 0.0)
    _pin_rib(stock, half)

    max_tea, certified, stations = _certify(stock, start, end)

    # Configuration pin: this is the un-refined regime the docstring describes. If
    # the certifier ever starts bisecting here, the test is measuring a different
    # thing and must be re-derived rather than believed.
    assert stations == 1, f"expected a single un-refined station pair, got {stations}"

    violations = _violating_centres(stock, start, end)
    assert not (certified and violations), (
        f"certified a motion with {len(violations)} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{start} -> {end}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}, reported max_tea={max_tea:.6f}; "
        f"first violating centre {violations[0] if violations else None}, worst reading {_sample(stock, 0.0, 0.0)[0]:.6f} rad at the rib centre"
    )


def test_certified_long_motion_has_no_cap_violating_centre():
    """The same rib, but a long oblique motion the certifier must adaptively refine to reach.

    1.5 tool diameters long and not axis-aligned, so the verdict is the product of
    real bisection rather than one lucky station pair. Measured: 35 stations, i.e.
    the certifier drove the spacing down to 0.75/16 = 0.046875 -- the coarsest
    spacing at which a positive guarded cap exists at this cap -- and stopped,
    because at that spacing both flanking stations pass. The violating window is
    0.005646 wide, 8.3x smaller, so it sits entirely between two stations.

    This is the shape of the defect: refinement does not converge onto the
    violation, it converges onto the spacing at which the guard first admits a
    verdict, and a feature narrower than that spacing is invisible.
    """
    stock = _rib_stock()
    _pin_rib(stock, 0.5 * SHORT_MOTION)

    max_tea, certified, stations = _certify(stock, LONG_MOTION_START, LONG_MOTION_END)

    # Configuration pin: the certifier really did refine, so this is not the
    # single-station-pair case re-tested under another name.
    assert stations > 1, f"expected adaptive refinement on a long motion, got {stations} station(s)"

    violations = _violating_centres(stock, LONG_MOTION_START, LONG_MOTION_END)
    assert not (certified and violations), (
        f"certified a REFINED motion with {len(violations)} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{LONG_MOTION_START} -> {LONG_MOTION_END}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}, "
        f"reported max_tea={max_tea:.6f} against a true {_sample(stock, 0.0, 0.0)[0]:.6f} rad at the rib centre"
    )


def test_no_alignment_of_a_crossing_motion_is_falsely_certified():
    """Sliding the motion along its own direction does not rescue the verdict.

    The stations of an adaptive bisection sit at dyadic fractions of the segment,
    so WHERE a motion starts decides whether one lands on the violation. If the red
    only appeared at a hand-picked offset it would be an alignment artifact worth
    discounting. It is not: sliding a 0.6-long crossing motion through four whole
    leaf spacings, a large majority of offsets are certified while the rib centre
    they pass over is fully immersed. Measured on the +/-0.15 sweep: 24 of 33, each
    after 31 to 39 stations of genuine refinement.
    """
    stock = _rib_stock()
    _pin_rib(stock, 0.5 * SHORT_MOTION)

    falsely_certified = []
    for i in range(SWEEP_OFFSETS):
        offset = -SWEEP_OFFSET_SPAN + 2.0 * SWEEP_OFFSET_SPAN * i / (SWEEP_OFFSETS - 1)
        start = (-0.5 * SWEEP_MOTION + offset, 0.0)
        end = (0.5 * SWEEP_MOTION + offset, 0.0)
        # Configuration pin: the motion really does pass over the violating window,
        # so a "certified" verdict here is a claim about geometry it crosses.
        assert start[0] < 0.0 < end[0], f"offset {offset!r} moves the motion off the rib centre"
        _max_tea, certified, stations = _certify(stock, start, end)
        if certified:
            falsely_certified.append((offset, stations))

    assert not falsely_certified, (
        f"{len(falsely_certified)} of {SWEEP_OFFSETS} alignments of a {SWEEP_MOTION}-long motion over a fully immersed rib centre "
        f"were certified at cap={CAP_RADIANS:.6f}; first {falsely_certified[0]} as (offset, stations)"
    )


def test_machined_stock_certified_motion_has_no_cap_violating_centre():
    """The rib reached by material removal alone -- so the red is not a constructor artifact.

    `_rib_stock` states the rib through `Stock`'s boundary/hole constructor, which
    invites the objection that no toolpath could ever produce it. This test removes
    the objection: the same rib is left in a solid block by a bore and ONE full
    circular contour pass with a cutter of the query's own radius -- the ordinary
    consequence of a step-over that overshoots. Nothing else is touched.

    Measured: stations at x = -/+0.0125 read 0.32149 (against 0.32190 for the ideal
    rib -- the chain-of-disks under-coverage in `subtract_arc_sweep`, well inside
    every margin), the centre reads a full turn, and the certifier again returns
    ``cap_certified = True`` after a single station pair.
    """
    stock = _machined_rib_stock()
    half = 0.5 * SHORT_MOTION
    start, end = (-half, 0.0), (half, 0.0)
    _pin_rib(stock, half)

    max_tea, certified, stations = _certify(stock, start, end)

    violations = _violating_centres(stock, start, end)
    assert not (certified and violations), (
        f"certified a motion over a rib built by removal only, with {len(violations)} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{start} -> {end}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}, reported max_tea={max_tea:.6f}"
    )


def test_certified_motion_whose_stations_report_no_contact_at_all_has_no_cap_violating_centre():
    """The sharpest form: every station reports ZERO engagement, and the certificate is still false.

    Open the rib into a 135 deg sector and the flanking stations stop touching it
    altogether -- the rim of a cutter half a spacing off centre crosses the rib's
    radial band only near +/-90 deg from the sector, which the sector does not
    reach. Measured: both stations read ``total_tea`` = 0, `certify_segment_tea`
    returns ``max_tea = 0.0`` with ``cap_certified = True``, and mid-motion the
    exact oracle reads 2.35619 rad = 3*pi/4, the sector's own extent, exceeding the
    90 deg cap over 64 of 401 probed centres.

    So the certifier does not merely under-estimate the engagement here; it
    certifies a motion it believes never touches material, over stock the cutter is
    at one point 135 deg engaged in. No guard on a growth bound can repair a
    verdict drawn from two measurements that are both identically zero -- which is
    why this witness constrains the repair more tightly than the others.
    """
    stock = _sector_rib_stock()
    half = 0.5 * SHORT_MOTION
    start, end = (-half, 0.0), (half, 0.0)

    # Configuration pins: the stations really are out of contact, and the sector's
    # centre really is engaged over its whole extent (a closed-form value, so a
    # mis-built sector cannot masquerade as a violation).
    for sign in (-1.0, 1.0):
        run, exceeded = _sample(stock, sign * half, 0.0)
        assert run == 0.0 and not exceeded, f"station at x={sign * half!r} should be clear of the sector, got max_run_tea {run!r}"
    centre_run, centre_exceeded = _sample(stock, 0.0, 0.0)
    assert abs(centre_run - SECTOR_EXTENT) <= TEA_REPORTING_SLACK, f"sector centre should engage exactly its own extent {SECTOR_EXTENT!r}, got {centre_run!r}"
    assert centre_exceeded, "sector centre should exceed the 90 deg cap -- the exact oracle disagrees, so the sector is not built"

    max_tea, certified, stations = _certify(stock, start, end)

    violations = _violating_centres(stock, start, end)
    assert not (certified and violations), (
        f"certified a motion reporting NO contact (max_tea={max_tea:.6f}) with {len(violations)} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{start} -> {end}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}; "
        f"true engagement at the sector centre {centre_run:.6f} rad"
    )


def test_a_motion_clear_of_the_rib_centre_is_soundly_certified():
    """Green control: on this same stock the certifier is right when it says yes.

    A short motion far from the rib centre never immerses the cutter -- the rim
    crosses the rib twice and reads about 0.02 rad -- and the certifier certifies
    it. Without this control the reds above would be consistent with a harness that
    calls every certificate false.
    """
    stock = _rib_stock()
    start, end = (0.2, 0.0), (0.25, 0.0)

    max_tea, certified, stations = _certify(stock, start, end)
    violations = _violating_centres(stock, start, end)

    # Configuration pin: the cutter really is in contact along this motion, so the
    # sound verdict is about a measured engagement and not about an empty rim.
    assert 0.0 < max_tea < CAP_RADIANS, f"expected a small but non-zero engagement along {start} -> {end}, got max_tea={max_tea!r}"
    assert certified, f"expected a sound certificate for {start} -> {end} (stations={stations}), got cap_certified=False"
    assert not violations, f"the oracle found {len(violations)} violating centres on a motion that should never violate: first {violations[0]}"


def test_a_motion_whose_station_lands_on_the_violation_is_refused():
    """Green control: when refinement DOES land a station on the rib centre, the certifier refuses.

    The symmetric 0.1-long motion bisects onto the rib centre, so a station reads
    the full turn and the guarded test fails all the way to the floor. The defect
    is therefore a blind spot BETWEEN stations, not a certifier that cannot see a
    violation it measures -- which is exactly why the repair belongs in the
    guard/refinement, not in the station predicate.
    """
    stock = _rib_stock()
    start, end = (-0.05, 0.0), (0.05, 0.0)

    max_tea, certified, stations = _certify(stock, start, end)
    violations = _violating_centres(stock, start, end)

    # Configuration pin: there really is something to refuse.
    assert violations, "control is vacuous: the oracle found no violating centre on this motion"
    assert not certified, f"expected refusal for {start} -> {end} (stations={stations}, max_tea={max_tea:.6f}), got cap_certified=True"
