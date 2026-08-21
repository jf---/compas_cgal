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
The rib is built below in three different TOPOLOGIES over two API paths -- a
boundary with a hole, a boolean difference, and one plain simple polygon; the
first and third both enter through `Stock.__init__` -- so no verdict here rests on
a single construction route. The removal-only route is the one that answers
"could a toolpath produce this", and it does.

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
to say HOW BADLY and as a CONFIGURATION PIN alongside the exact predicate, never
as a verdict -- it is a REPORTING double subject to the known one-full-turn
harvest defect
(`test_growth_bound.py::test_reported_engagement_never_exceeds_a_full_turn`),
which does fire elsewhere in this very geometry (9.23669 rad at a tangency
configuration one ulp from a clean reading). Every pin that reads it is paired
with the exact `cap_exceeded` and sits far from that tangency.

A THIRD ORACLE agrees. Dense `Stock.contains` sampling of the cutter rim shares no
code with the engagement harvest -- it is a point-in-region query on the
arrangement, not an arc harvest -- and it reproduces every reading quoted here to
within its own sampling resolution: at the spiral witness 1.899093 rad against
`engagement_at`'s 1.899052 (20,000 rim samples, 3.1e-4 rad resolution), and 0 of
20,011 rim points outside material at the annular rib's centre, confirming the 2*pi
is genuine full immersion rather than the harvest defect.

HOW FAR IT GOES. Open the rib into a 135 deg sector and the flanking stations stop
touching it at all: the certifier then returns ``max_tea = 0.0`` -- "the cutter
never contacted material anywhere on this motion" -- for a motion on which it is,
at one point, 135 deg engaged. A guard added to a growth bound cannot repair a
verdict drawn from two measurements that are both identically zero, which is why
`test_certified_motion_whose_stations_report_no_contact_at_all_has_no_cap_violating_centre`
constrains the repair more tightly than the rest.

IT IS NOT A COINCIDENCE, AND THIS IS THE POINT THAT SIZES THE REPAIR. A rib of
CONSTANT radius has to match the tool radius to within about 2.5 * tau (2% of r
here) for any of this to happen, which invites the reading that the failure is a
codimension-1 accident. It is not. Put the rib on a SPIRAL -- centreline
``rho(theta) = r + k*theta``, so its radius SWEEPS THROUGH the tool radius instead
of matching it -- and no coincidence is engineered anywhere: a spiral necessarily
contains a point where its local curvature radius equals the cutter's, whatever
that cutter's radius is. The result is a REGION of falsely-certified motions of
positive area (measured: 12 of 33 probe centres on a 3x11 grid) and every motion
direction through the witness centre is certified at ``stations = 1`` (6 of 6
tried). So a thin leftover sliver along ANY spiral or ramped contour pass carries
this failure. The repair must be sized against "whenever a thin leftover curves
with the tool", not against a 1% radius match.

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
#
# TOLERANCE ON THE RIB'S RADIUS. A rib of CONSTANT radius must match the tool
# radius for any of this to bite, but the window is 2.5 * tau wide (2% of r), not
# tau/2: measured red at every mismatch from 0 to 2.5 tau with the motion held at
# the origin, sound from 3 tau. Its SHAPE depends on where the motion sits -- with
# the motion recentred on the mismatch there is a REFUSAL GAP at 1.25-1.5 tau,
# with red resuming at 2.0-2.5 tau. That gap is placement luck, not detection: the
# bisection happens to land a station on the violation there. `_spiral_rib_stock`
# removes the question entirely by letting the radius sweep.
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
# all, and LONG enough that the cap-violating window -- measured at
# |x| <= 0.002823, i.e. 22.6% of this motion -- falls strictly between the two
# stations. The upper bound is where the guarded cap `cap - tea_guard(L/2, r)`
# reaches zero, obtained by BISECTING that condition: 0.056111 = 0.112 r at
# cap = pi/2, and 0.178464 at cap = pi. Do NOT substitute the small-`d`
# asymptotic `r * (cap/4)^2` for it -- that drops the endpoint-drift term
# `4*asin(d/2r)`, which is not negligible at these spacings, and over-estimates
# the certifiable range by 37% at cap = pi/2 and 73% at cap = pi. 0.025 is 45% of
# the bisected ceiling, comfortably inside both bounds.
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

# Probes along a motion for the exact oracle, used for the COUNT in failure
# messages. Detection does not depend on it: liveness is established by probing the
# segment's closest point to the rib centre, which is derived rather than searched
# for (`_closest_point_on`) and so cannot be stepped over by a sample grid.
SCAN_SAMPLES = 400

# Centre of curvature the closed and open ribs are built around. Every TEA peak on
# a motion past those ribs sits at the motion's closest approach to this point.
RIB_CENTRE = (0.0, 0.0)

# Fraction of the closed-form engaged run a probe must reach before the green
# control will believe the cutter is in contact. A cutter centred `s` from the rib
# centre (s >> tau) crosses the rib in two arcs of tau / s radians each -- measured
# within 2% to 4.5% of that closed form over s in [0.15, 0.3]. Half of it is the
# floor: a 2x margin over the closed form, and eight decades above the 1.6e-07 a
# hairline rib returns, which is what a bare `run > 0` would have accepted.
CONTACT_FLOOR_FRACTION = 0.5

# --- The spiral rib: the same failure with NO radius coincidence -------------
#
# Centreline `rho(theta) = TOOL_RADIUS + SPIRAL_GROWTH * theta`, walls at
# +/- RIB_THICKNESS/2. SPIRAL_GROWTH = 0.005 makes the radius sweep 0.488 -> 0.512
# over the arc below, i.e. +/-2.4% of r: the tool radius is crossed on the way
# past rather than matched. Larger growth rates were measured too -- 0.01 (+/-4.8%)
# and 0.02 (+/-9.6%) are both still falsely certified at `stations = 1`; 0.04
# (+/-19%) is where the rib stops tracking the rim long enough and the certifier
# starts refusing. So the failure needs the sliver to curve WITH the tool, not to
# match it.
SPIRAL_GROWTH = 0.005

# Half the spiral's angular extent (radians). 2.4 rad each way spans 4.8 rad of
# arc -- long enough to carry the engaged run of 1.9 rad found below with room to
# spare, short enough that the two turns never approach each other.
SPIRAL_HALF_TURN = 2.4

# Vertices per spiral wall, set so the angular step matches the closed ring's
# 2*pi / RIB_FACETS: the measured sagitta is 3.75e-5, 0.94% of RIB_THICKNESS,
# identical to the ring's.
SPIRAL_STEPS = round(RIB_FACETS * 2.0 * SPIRAL_HALF_TURN / (2.0 * math.pi))

# Centre of the spiral witness motion, in model units. It is where the cutter rim
# best osculates the spiral rib, LOCATED BY SEARCH over the centre plane and
# stated here so the test is deterministic. Nothing about it is fine-tuned:
# `test_no_spiral_probe_centre_is_falsely_certified` shows the falsely-certified
# set has positive area, and this centre is picked from it for its MARGINS rather
# than its peak -- both stations sit 19% under the guarded cap while the interior
# runs 21% over the cap.
#
# NOT the closed-form osculating centre of the centreline. Bisecting
# `Rc(theta) = (rho^2 + k^2)^(3/2) / (rho^2 + 2*k^2)` to `Rc = r` gives
# theta* = 0.0049984 and the centre (2.49976e-05, 4.99969e-03), where the same
# centred 0.025 motion is REFUSED at 10 stations. That refusal is a property of
# THAT PLACEMENT, not of the curvature-matching point, and it is NOT detection: its
# two endpoint stations read 0.708505 and 0.712094, both ABOVE the guarded cap of
# 0.574493, so no certificate is available at that spacing and refinement is forced
# until a station happens to land on the violation. A repair may NOT infer from it
# that the curvature-matching region is handled -- move the motion 0.0016 in y and
# the same region certifies.
SPIRAL_PROBE_CENTRE = (0.0, 0.0034)

# Probe-centre region sweep for the spiral: a 3 x 11 grid deliberately spanning
# sound, refused AND falsely-certified centres (y from 0 to 0.010 crosses the whole
# osculating band), so the count it reports is an honest fraction and not a window
# chosen to be red. Measured: 12 red, 12 soundly certified, 9 refused.
SPIRAL_SWEEP_X = (-0.004, 0.0, 0.004)
SPIRAL_SWEEP_ROWS = 11
SPIRAL_SWEEP_Y_STEP = 0.001

# Directions probed through SPIRAL_PROBE_CENTRE, in radians. Chosen to include the
# axis-aligned pair and three obliques, none of them aligned with the spiral's
# tangent at the witness. Measured: all six certified at `stations = 1`.
SPIRAL_DIRECTIONS = (0.0, math.pi / 6.0, math.pi / 4.0, math.pi / 3.0, math.pi / 2.0, 2.0 * math.pi / 3.0)

# Liveness floor for the spiral probe-centre grid: how many of its motions must
# carry a real cap violation before the sweep's verdict means anything. The grid
# deliberately spans dead centres as well as live ones, so an exact count would be
# brittle; measured 21 of 33 live, and 15 leaves 29% of slack for faceting or
# kernel drift while still catching total degeneration -- a rib thinned to a
# hairline and a spiral grown to 0.04 both measure 0 of 33.
SPIRAL_SWEEP_LIVE_FLOOR = 15

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


def _spiral_rib_stock(growth: float = SPIRAL_GROWTH, thickness: float = RIB_THICKNESS) -> Stock:
    """A rib whose radius SWEEPS THROUGH the tool radius instead of matching it.

    The walls are two Archimedean spirals ``rho(theta) = TOOL_RADIUS + growth *
    theta +/- thickness / 2``. Because the centreline radius varies monotonically,
    the rib necessarily contains a point where its local curvature radius equals
    the cutter's, whatever the cutter's radius is -- so this construction removes
    the "the rib has to match the tool" objection to `_rib_stock` entirely, and
    with it the reading that the failure is a codimension-1 accident.

    Args:
        growth: Radius gained per radian of the centreline spiral.
        thickness: Radial thickness of the rib in model units.

    Returns:
        A `Stock` whose material is the spiral sliver, traced outer wall then
        inner wall as one simple polygon.
    """
    angles = [-SPIRAL_HALF_TURN + 2.0 * SPIRAL_HALF_TURN * i / SPIRAL_STEPS for i in range(SPIRAL_STEPS + 1)]

    def wall(theta: float, offset: float) -> tuple[float, float, float]:
        rho = TOOL_RADIUS + growth * theta + offset
        return (rho * math.cos(theta), rho * math.sin(theta), 0.0)

    points = [wall(theta, 0.5 * thickness) for theta in angles]
    points += [wall(theta, -0.5 * thickness) for theta in reversed(angles)]
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


def _closest_point_on(start: tuple[float, float], end: tuple[float, float], target: tuple[float, float]) -> tuple[float, float]:
    """The point of a motion nearest *target*, computed by projection and clamped to the segment.

    For every rib here TEA decreases with the distance to the rib's centre of
    curvature, so this is where the motion's TEA peaks -- the interior minimum of a
    convex distance, which is exactly the point endpoint attainment does not reach.
    Deriving it beats scanning for it: it needs one oracle call instead of a grid,
    and it cannot be missed by a sample grid stepping over a narrow window.

    Args:
        start: ``(x, y)`` of the motion start.
        end: ``(x, y)`` of the motion end.
        target: ``(x, y)`` the rib is centred on.

    Returns:
        The ``(x, y)`` on the segment closest to *target*.
    """
    dx, dy = end[0] - start[0], end[1] - start[1]
    length_sq = dx * dx + dy * dy
    t = ((target[0] - start[0]) * dx + (target[1] - start[1]) * dy) / length_sq
    t = min(1.0, max(0.0, t))
    return (start[0] + t * dx, start[1] + t * dy)


def _pin_cap_violation_at(stock: Stock, point: tuple[float, float], what: str) -> float:
    """LIVENESS: assert the EXACT oracle finds the cap exceeded at *point*.

    Decided by `engagement_at`'s ``cap_exceeded`` alone -- the certifier is not
    consulted, so this cannot be satisfied by the verdict the test is about to
    challenge. Its job is to fail LOUDLY, naming the construction, when a witness
    stops being a witness: a rib thinned to nothing or a spiral growing too fast
    both leave the certifier with nothing to be wrong about, and a red that goes
    green that way would read as a repair.

    Args:
        stock: The stock region to measure against.
        point: ``(x, y)`` cutter centre that must violate.
        what: Name of the construction, for the failure message.

    Returns:
        The reported largest engaged run at *point* (for failure messages only).

    Raises:
        AssertionError: If the exact oracle does not report the cap exceeded there.
    """
    run, exceeded = _sample(stock, *point)
    assert exceeded, f"{what}: the exact oracle finds NO cap violation at {point} (max_run_tea {run!r}) -- the construction is dead, so nothing here tests the certifier"
    return run


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
        half_spacing: Distance from the rib centre at which the rib is probed for
            the straddling regime. For the un-refined witnesses this is also the
            motion's own half-spacing; for the refined ones the certifier's leaf
            half-spacing differs (0.0234 and 0.0187), and the pin's job there is
            only to establish that the rib is built and that the origin -- which
            lies on every red motion -- genuinely violates.

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

    # 1. LIVENESS -- exact oracle only, the certifier is not consulted. The motion's
    #    closest approach to the rib centre genuinely exceeds the cap, so there IS
    #    something for the certificate to be wrong about.
    peak = _pin_cap_violation_at(stock, _closest_point_on(start, end, RIB_CENTRE), "annular rib, short motion")

    # 2. VERDICT -- what the repair must change.
    max_tea, certified, stations = _certify(stock, start, end)

    # Configuration pin: this is the un-refined regime the docstring describes. If
    # the certifier ever starts bisecting here, the test is measuring a different
    # thing and must be re-derived rather than believed.
    assert stations == 1, f"expected a single un-refined station pair, got {stations}"

    assert not certified, (
        f"certified a motion whose interior reaches {peak:.6f} rad, with {len(_violating_centres(stock, start, end))} of {SCAN_SAMPLES + 1} "
        f"probed centres exactly over the cap: {start} -> {end}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, "
        f"stations={stations}, reported max_tea={max_tea:.6f}"
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

    # 1. LIVENESS -- exact oracle only. Derived rather than scanned for, which
    #    matters here: the violating window is 0.75% of this motion's length, so a
    #    sample grid is a poor instrument for establishing that it exists at all.
    peak = _pin_cap_violation_at(stock, _closest_point_on(LONG_MOTION_START, LONG_MOTION_END, RIB_CENTRE), "annular rib, long oblique motion")

    # 2. VERDICT.
    max_tea, certified, stations = _certify(stock, LONG_MOTION_START, LONG_MOTION_END)

    # Configuration pin: the certifier really did refine, so this is not the
    # single-station-pair case re-tested under another name.
    assert stations > 1, f"expected adaptive refinement on a long motion, got {stations} station(s)"

    assert not certified, (
        f"certified a REFINED motion whose interior reaches {peak:.6f} rad, with "
        f"{len(_violating_centres(stock, LONG_MOTION_START, LONG_MOTION_END))} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{LONG_MOTION_START} -> {LONG_MOTION_END}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}, reported max_tea={max_tea:.6f}"
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

        # 1. LIVENESS, per motion, exact oracle only. Both endpoints lie at y = 0
        #    and straddle x = 0, so the rib centre is EXACTLY on this segment -- no
        #    projection, no rounding -- and the oracle says it violates. One cheap
        #    call per motion establishes that every "certified" below is a claim
        #    about geometry that genuinely breaks the cap.
        assert start[0] < 0.0 < end[0] and start[1] == end[1] == RIB_CENTRE[1], f"offset {offset!r} moves the motion off the rib centre"
        _pin_cap_violation_at(stock, RIB_CENTRE, f"annular rib, alignment offset {offset!r}")

        # 2. VERDICT.
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

    # 1. LIVENESS -- exact oracle only. Independent of `_rib_stock`: if the two
    #    removals ever stop leaving a rib, this fails naming the machined build
    #    rather than passing as a repair.
    peak = _pin_cap_violation_at(stock, _closest_point_on(start, end, RIB_CENTRE), "machined rib (bore + contour pass)")

    # 2. VERDICT.
    max_tea, certified, stations = _certify(stock, start, end)

    assert not certified, (
        f"certified a motion over a rib built by removal only, whose interior reaches {peak:.6f} rad, with "
        f"{len(_violating_centres(stock, start, end))} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
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

    # Configuration pin: the stations really are out of contact, which is the whole
    # point of this witness and is what makes it un-repairable by a guard.
    for sign in (-1.0, 1.0):
        run, exceeded = _sample(stock, sign * half, 0.0)
        assert run == 0.0 and not exceeded, f"station at x={sign * half!r} should be clear of the sector, got max_run_tea {run!r}"

    # 1. LIVENESS -- exact oracle only, plus a closed-form shape pin: the sector's
    #    engaged run at its centre must be exactly its own angular extent, so a
    #    sector that degenerated to a hairline (or closed into a ring) fails here
    #    naming the construction instead of passing as a repair.
    peak = _pin_cap_violation_at(stock, _closest_point_on(start, end, RIB_CENTRE), f"{math.degrees(SECTOR_EXTENT):.1f} deg sector rib")
    assert abs(peak - SECTOR_EXTENT) <= TEA_REPORTING_SLACK, f"sector centre should engage exactly its own extent {SECTOR_EXTENT!r}, got {peak!r}"

    # 2. VERDICT.
    max_tea, certified, stations = _certify(stock, start, end)

    assert not certified, (
        f"certified a motion reporting NO contact (max_tea={max_tea:.6f}) whose interior reaches {peak:.6f} rad, with "
        f"{len(_violating_centres(stock, start, end))} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{start} -> {end}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}"
    )


def test_certified_spiral_rib_motion_has_no_cap_violating_centre():
    """No radius coincidence at all: the rib's radius sweeps past the tool's, and the certificate is still false.

    `_rib_stock` needs its constant radius to match the tool's within ~2.5 tau,
    which invites the reading that this whole failure is a codimension-1 accident.
    A spiral rib settles that: its centreline radius runs 0.488 -> 0.512 across the
    sliver, so it CROSSES the tool radius on the way past. Every spiral does --
    which is why a thin leftover along any spiral or ramped contour pass carries
    this failure, whatever the cutter's radius.

    Measured, on the 0.025-long motion through `SPIRAL_PROBE_CENTRE`:

    | quantity | value |
    | --- | --- |
    | station at x = -0.0125 | ``max_run_tea`` 0.447584, `cap_exceeded` False |
    | station at x = +0.0125 | ``max_run_tea`` 0.462523, `cap_exceeded` False |
    | guarded cap | 0.574493 -- both stations 19% under it |
    | motion peak | 1.910381 rad at x = 6.25e-05, i.e. 1.22x the cap |
    | `certify_segment_tea` | ``max_tea`` 0.462523, `cap_certified` **True**, ``stations`` 1 |
    | exact oracle | 44 of 401 probed centres over the cap |

    Cross-checked against the independent `Stock.contains` rim oracle, which shares
    no code with the engagement harvest: 1.899093 rad at the centre against
    `engagement_at`'s 1.899052, and 0.447677 / 0.462757 at the two stations against
    0.447584 / 0.462523 -- agreement to within the 3.1e-4 rad resolution of 20,000
    rim samples.

    Note this is a HARDER case for the repair than the annular rib, not an easier
    one: the stations here read 0.45-0.46 rather than 0.32, so simply tightening
    the guard would not catch it either.
    """
    stock = _spiral_rib_stock()
    px, py = SPIRAL_PROBE_CENTRE
    half = 0.5 * SHORT_MOTION
    start, end = (px - half, py), (px + half, py)

    # Configuration pin: both stations really are in contact and really do pass the
    # guarded cap, so this is a genuine blind-spot verdict, not an out-of-contact
    # one (that form is pinned separately by the sector test).
    for probe in (start, end):
        run, exceeded = _sample(stock, *probe)
        assert 0.0 < run < CAP_RADIANS, f"station {probe} should be engaged but under the cap, got max_run_tea {run!r}"
        assert not exceeded, f"station {probe} should be under the cap, got max_run_tea {run!r}"

    # 1. LIVENESS -- exact oracle only. A rib thinned to a hairline or a spiral
    #    grown too fast leaves the certifier nothing to be wrong about; both fail
    #    HERE, naming the spiral, rather than turning this red green.
    peak = _pin_cap_violation_at(stock, _closest_point_on(start, end, SPIRAL_PROBE_CENTRE), f"spiral rib (growth {SPIRAL_GROWTH!r})")

    # 2. VERDICT.
    max_tea, certified, stations = _certify(stock, start, end)

    # Configuration pin: still the un-refined regime, so the verdict is two exact
    # measurements plus the analytic guard and nothing adaptive.
    assert stations == 1, f"expected a single un-refined station pair, got {stations}"

    assert not certified, (
        f"certified a motion over a SPIRAL rib -- no radius coincidence -- whose interior reaches {peak:.6f} rad, with "
        f"{len(_violating_centres(stock, start, end))} of {SCAN_SAMPLES + 1} probed centres exactly over the cap: "
        f"{start} -> {end}, r={TOOL_RADIUS}, cap={CAP_RADIANS:.6f}, stations={stations}, reported max_tea={max_tea:.6f}"
    )


def test_no_spiral_probe_centre_is_falsely_certified():
    """The falsely-certified set on a spiral rib has positive AREA, and every direction through it fails.

    `SPIRAL_PROBE_CENTRE` was located by search, which invites the objection that a
    searched point proves nothing about ordinary geometry. This test answers it in
    two directions at once.

    ACROSS THE PLANE: a 3x11 grid of probe centres spanning y = 0 to 0.010 -- the
    whole osculating band, deliberately including centres that are soundly
    certified and centres that are refused -- yields 12 false certificates, 12
    sound certificates and 9 refusals. A measure-zero coincidence cannot occupy a
    third of a grid laid across it.

    ACROSS DIRECTION: all six directions probed through the witness centre (0, 30,
    45, 60, 90, 120 deg) are certified at ``stations = 1``. The failure does not
    need the motion to be aligned with anything.

    LIVENESS IS SEPARATE FROM THE VERDICT here, and deliberately so. This sweep's
    only claim would otherwise be "nothing was certified", which a DEAD
    construction satisfies perfectly: measured, a rib thinned to 1e-9 and a spiral
    grown at 0.04 each yield 0 of 33 live motions and 0 false certificates, so the
    test would pass while the certifier remained exactly as unsound. The live count
    below is taken from the exact oracle alone and floored, so degeneration fails
    the test instead of silencing it.
    """
    stock = _spiral_rib_stock()
    half = 0.5 * SHORT_MOTION
    motions = [((px - half, py), (px + half, py)) for px in SPIRAL_SWEEP_X for py in (row * SPIRAL_SWEEP_Y_STEP for row in range(SPIRAL_SWEEP_ROWS))]

    # 1. LIVENESS -- exact oracle only. The grid deliberately spans dead centres as
    #    well as live ones, so this is a floor rather than a count.
    live_grid = [motion for motion in motions if _violating_centres(stock, *motion, samples=100)]
    assert len(live_grid) >= SPIRAL_SWEEP_LIVE_FLOOR, (
        f"only {len(live_grid)} of {len(motions)} grid motions carry a real cap violation (floor {SPIRAL_SWEEP_LIVE_FLOOR}) -- "
        f"the spiral rib is degenerate, so a clean sweep below would mean nothing"
    )
    _pin_cap_violation_at(stock, SPIRAL_PROBE_CENTRE, f"spiral rib (growth {SPIRAL_GROWTH!r})")

    # 2. VERDICT, over the same grid plus a fan of directions through the witness.
    falsely_certified = []
    for start, end in live_grid:
        _max_tea, certified, stations = _certify(stock, start, end)
        if certified:
            falsely_certified.append((round(start[0] + half, 4), round(start[1], 4), stations))

    px, py = SPIRAL_PROBE_CENTRE
    for angle in SPIRAL_DIRECTIONS:
        dx, dy = half * math.cos(angle), half * math.sin(angle)
        motion = ((px - dx, py - dy), (px + dx, py + dy))
        # Every direction passes through the witness centre, which liveness above
        # already showed violates -- so each is a live motion by construction.
        _max_tea, certified, stations = _certify(stock, *motion)
        if certified:
            falsely_certified.append((f"dir {math.degrees(angle):.0f}deg", stations))

    assert not falsely_certified, (
        f"{len(falsely_certified)} false certificates over a spiral rib: {len(live_grid)} live motions of a "
        f"{len(SPIRAL_SWEEP_X)}x{SPIRAL_SWEEP_ROWS} probe-centre grid plus {len(SPIRAL_DIRECTIONS)} directions "
        f"through {SPIRAL_PROBE_CENTRE}; first {falsely_certified[0]}"
    )


def test_a_motion_clear_of_the_rib_centre_is_soundly_certified():
    """Green control: on this same stock the certifier is right when it says yes.

    A short motion far from the rib centre never immerses the cutter -- the rim
    crosses the rib twice and reads about 0.02 rad -- and the certifier certifies
    it. Without this control the reds above would be consistent with a harness that
    calls every certificate false.

    Its liveness check is the mirror image of theirs, and asked of the exact oracle
    rather than of ``max_tea``: the cutter must be genuinely IN CONTACT along the
    motion. A control that went green because its stock had vanished would prove
    nothing, and reading contact off the certifier's own report would make the
    check circular.
    """
    stock = _rib_stock()
    start, end = (0.2, 0.0), (0.25, 0.0)

    # 1. LIVENESS -- exact oracle only, mirrored: the cutter is in contact all along
    #    this motion, and NO centre on it violates. "In contact" is measured against
    #    the closed form, not against zero: a hairline rib still returns a non-zero
    #    reading (1.6e-07 at thickness 1e-9), so `run > 0` would let this control
    #    go green on a stock that had effectively vanished.
    for probe in (start, _closest_point_on(start, end, RIB_CENTRE), end):
        run, exceeded = _sample(stock, *probe)
        floor = CONTACT_FLOOR_FRACTION * RIB_THICKNESS / math.hypot(probe[0] - RIB_CENTRE[0], probe[1] - RIB_CENTRE[1])
        assert run >= floor, f"control is vacuous: cutter barely in contact at {probe} (max_run_tea {run!r} < {floor!r}) -- the rib is not built"
        assert not exceeded, f"control is mis-stated: {probe} exceeds the cap (max_run_tea {run!r}), so this motion is not a sound-certificate case"
    violations = _violating_centres(stock, start, end)
    assert not violations, f"the oracle found {len(violations)} violating centres on a motion that should never violate: first {violations[0]}"

    # 2. VERDICT -- and here the certifier is RIGHT.
    max_tea, certified, stations = _certify(stock, start, end)
    assert certified, f"expected a sound certificate for {start} -> {end} (stations={stations}, max_tea={max_tea:.6f}), got cap_certified=False"


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

    # 1. LIVENESS -- exact oracle only: there really is something to refuse. Without
    #    it, a stock that had degenerated to nothing would still be "refused" for
    #    entirely the wrong reason and this control would read as evidence.
    peak = _pin_cap_violation_at(stock, _closest_point_on(start, end, RIB_CENTRE), "annular rib, station-on-violation control")

    # 2. VERDICT -- and here the certifier is RIGHT.
    max_tea, certified, stations = _certify(stock, start, end)
    assert not certified, f"expected refusal for {start} -> {end} over an interior reaching {peak:.6f} rad (stations={stations}, max_tea={max_tea:.6f}), got cap_certified=True"
