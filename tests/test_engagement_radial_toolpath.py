"""Radius-regulated trochoidal generator tests.

`compas_cgal.engagement_radial_toolpath` adds a second regulated knob -- the loop
radius -- to the advance regulation of `compas_cgal.engagement_toolpath`. These
tests hold the three things that makes or breaks:

1. the radius search is a SCAN of a ladder and must never become a bisection,
   because engagement is not monotone in the loop radius;
2. a rung is admissible only if it also still CUTS, because a loop drawn inward
   into already-swept void complies with any cap while removing nothing; and
3. a cap loose enough for the advance regulation alone must come out byte-identical
   to the advance-only generator, so the new knob is additive rather than a rewrite.

Pockets are sized so a generate-plus-audit round trip stays tractable: the audit
replays every motion against a depleting exact stock, which dominates the cost.
"""

import math

import pytest
from compas.geometry import Circle, Line, Polygon

from compas_cgal.engagement import _subtract_operation, audit_toolpath_engagement
from compas_cgal.engagement_radial_toolpath import (
    FULL_RADIUS_RUNG,
    RADIUS_LADDER_RUNGS,
    _largest_admissible_radius,
    _loop_reaches_material,
    _radius_ladder,
    radius_regulated_toolpath,
)
from compas_cgal.engagement_toolpath import (
    GUIDE_STEP_TOOL_DIAMETERS,
    InvalidEngagementCapDegreesError,
    InvalidGuideResolutionError,
    NonPositiveToolDiameterError,
    _cap_surrogate,
    _guide_chains,
    _GuideStation,
    _Regulation,
    _station_is_admissible,
    engagement_controlled_toolpath,
)
from compas_cgal.stock import Stock
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION, OperationType, ToolpathResult

# Same 6x4 pocket the advance-only generator's tests use, for the same reason: it
# exercises the full topology (a central skeleton chain plus four corner spokes)
# at a cost that keeps a full audit affordable.
POCKET = Polygon([[0, 0, 0], [6, 0, 0], [6, 4, 0], [0, 4, 0]])
TOOL_DIAMETER = 2.0
TOOL_RADIUS = 0.5 * TOOL_DIAMETER

# A cap the advance regulation alone already meets on this pocket, so every
# station's maximal circle is admissible and the ladder never leaves rung 0.
LOOSE_CAP_DEG = 120.0

# A cap it does not, so the ladder actually descends.
TIGHT_CAP_DEG = 40.0

# The pocket and station the non-monotonicity pin is built on. Twenty by twelve
# with a 2 mm tool puts the central skeleton chain at y = 6 with a maximal loop
# radius of 4.998, which is the regime the finding was measured in.
PIN_POCKET = Polygon([[0, 0, 0], [20, 0, 0], [20, 12, 0], [0, 12, 0]])
PIN_STATION_X = 14.0
PIN_STATION_Y = 6.0
PIN_FULL_RADIUS = 4.998
# Three maximal loops immediately behind the station, at the guide's own step, is
# the smallest depletion that reproduces the measured pattern: enough swept
# annulus for the middle of the ladder to fall into void, not so much that the
# whole ladder does.
PIN_PRECEDING_CENTRES = (12.5, 13.0, 13.5)
# The cap at which the pattern bites. Chosen because the ladder's engagement dips
# to ~59 deg around rung 8 and rises again below it, so only a narrow band of
# rungs complies and the admissible set is provably not an up-set.
PIN_CAP_DEG = 60.0


def _regulation(cap_deg):
    """The validated parameter object the search helpers take."""
    return _Regulation.build(
        tool_diameter=TOOL_DIAMETER,
        tea_cap_deg=cap_deg,
        guide_step_tool_diameters=GUIDE_STEP_TOOL_DIAMETERS,
        max_advance_tool_diameters=1.0,
        radial_clearance=None,
        cut_z=0.0,
        clearance_z=None,
    )


def _pin_state():
    """The stock, station and travel direction the non-monotonicity pin is read on."""
    stock = Stock(PIN_POCKET)
    for cx in PIN_PRECEDING_CENTRES:
        stock.subtract_annulus(cx, PIN_STATION_Y, PIN_FULL_RADIUS - TOOL_RADIUS, PIN_FULL_RADIUS + TOOL_RADIUS)
    station = _GuideStation(cx=PIN_STATION_X, cy=PIN_STATION_Y, radius=PIN_FULL_RADIUS, clockwise=True, tx=1.0, ty=0.0)
    return stock, station, (1.0, 0.0)


def _rung_flags(stock, station, advance, regulation):
    """Per rung, whether the cap predicate accepts that radius. Descending radius."""
    flags = []
    for radius in _radius_ladder(station.radius, 0.0, regulation.guide_step):
        candidate = _GuideStation(cx=station.cx, cy=station.cy, radius=radius, clockwise=station.clockwise, tx=station.tx, ty=station.ty)
        flags.append(_station_is_admissible(stock, candidate, advance, regulation.tool_radius, regulation.cap_ratio))
    return flags


def _bisect_rung(flags):
    """Smallest accepting rung a BISECTION would report, assuming an up-set.

    Written out here rather than imported because no such function exists in the
    package: this is the search the module must never contain, kept in the test so
    the difference between it and the scan is asserted rather than asserted about.
    """
    lo, hi, best = 0, len(flags) - 1, None
    while lo <= hi:
        mid = (lo + hi) // 2
        if flags[mid]:
            best = mid
            hi = mid - 1
        else:
            lo = mid + 1
    return best


def _op_signature(result: ToolpathResult):
    """Hashable, exact-float signature of the emitted operation sequence."""
    signature = []
    for op in result.operations:
        geometry = op.geometry
        if isinstance(geometry, Circle):
            center = geometry.frame.point
            key = ("circle", float(geometry.radius), tuple(float(v) for v in center), tuple(float(v) for v in geometry.point_at(0.0)))
        else:
            key = ("line", tuple(float(v) for v in geometry.start), tuple(float(v) for v in geometry.end))
        signature.append((op.operation.value, op.path_index, op.clockwise, key))
    return tuple(signature)


def _cut_length(result: ToolpathResult) -> float:
    """Total XY length of the cutting motions (circles plus bridges)."""
    total = 0.0
    for op in result.operations:
        if op.operation is not OperationType.CUT:
            continue
        geometry = op.geometry
        if isinstance(geometry, Circle):
            total += 2.0 * math.pi * float(geometry.radius)
        elif isinstance(geometry, Line):
            total += math.hypot(float(geometry.end[0]) - float(geometry.start[0]), float(geometry.end[1]) - float(geometry.start[1]))
    return total


def _chain_entry_indices(result: ToolpathResult):
    """Index of the machining circle each chain is ENTERED on, once per chain.

    The radial generator plunges once per PASS, and only a chain's first pass
    meets virgin stock, so the later plunges are re-entries into material this
    path has already opened and their circles stay in every statistic below.
    """
    entries = set()
    seen = set()
    for index, op in enumerate(result.operations):
        if op.operation is OperationType.PLUNGE and op.path_index not in seen:
            seen.add(op.path_index)
            entries.add(index + 1)
    return entries


def _replay(result: ToolpathResult, polygon: Polygon) -> Stock:
    """Depleting replay of a whole path, for the coverage check."""
    stock = Stock(polygon)
    for op in result.operations:
        geometry = op.geometry
        if isinstance(geometry, Line):
            z_start, z_end = float(geometry.start[2]), float(geometry.end[2])
            if z_start != z_end:
                if z_end < z_start:  # plunge bores a full disk at the cutting plane
                    stock.subtract_disk(float(geometry.end[0]), float(geometry.end[1]), TOOL_RADIUS)
                continue
            if z_start != 0.0:  # clearance-height rapid, cuts nothing
                continue
        _subtract_operation(stock, op, TOOL_RADIUS)
    return stock


def test_radius_ladder_descends_from_the_maximal_radius():
    ladder = _radius_ladder(PIN_FULL_RADIUS, 0.0, GUIDE_STEP_TOOL_DIAMETERS * TOOL_DIAMETER)

    assert ladder[FULL_RADIUS_RUNG] == PIN_FULL_RADIUS, "rung 0 must be the station's own maximal circle"
    assert len(ladder) == RADIUS_LADDER_RUNGS
    assert all(ladder[k] > ladder[k + 1] for k in range(len(ladder) - 1)), "rung index must order the radii strictly downward"


def test_radius_ladder_stops_at_what_the_station_already_cut():
    step = GUIDE_STEP_TOOL_DIAMETERS * TOOL_DIAMETER
    already = PIN_FULL_RADIUS - 3.5 * step

    ladder = _radius_ladder(PIN_FULL_RADIUS, already, step)

    # Only radii that reach beyond what this station already swept are offered:
    # anything at or below removes nothing, and offering it would let a pass
    # "progress" without cutting, which is what makes the sweep loop terminate.
    assert ladder == [PIN_FULL_RADIUS, PIN_FULL_RADIUS - step, PIN_FULL_RADIUS - 2.0 * step, PIN_FULL_RADIUS - 3.0 * step]
    assert _radius_ladder(PIN_FULL_RADIUS, PIN_FULL_RADIUS, step) == [], "a finished station must offer no rung at all"


def test_engagement_is_not_monotone_in_the_loop_radius():
    """THE PIN: the admissible rungs are not an up-set, so a bisection is unsound.

    Shrinking a machining circle does not monotonically lighten it. On the state
    built here the cap predicate accepts a band of rungs in the middle of the
    ladder and refuses rungs on BOTH sides of it, so "admissible" is not upward
    closed in the rung index and the bisection that would be correct on an up-set
    reports a different answer -- here, no answer at all, which would force the
    maximal circle the scan avoids.

    If this test ever fails because the pattern became monotone, that is a finding
    about the geometry, not licence to bisect: the scan is correct either way and
    the bisection is only correct on one of the two.
    """
    regulation = _regulation(PIN_CAP_DEG)
    stock, station, advance = _pin_state()

    flags = _rung_flags(stock, station, advance, regulation)
    accepting = [rung for rung, ok in enumerate(flags) if ok]
    assert accepting, "the pin state must have at least one admissible rung, else it pins nothing"

    largest = accepting[0]
    refused_below = [rung for rung in range(largest + 1, len(flags)) if not flags[rung]]
    assert refused_below, f"admissible rungs {accepting} form an up-set here -- this state no longer pins non-monotonicity"

    rung, forced = _largest_admissible_radius(stock, station, 0.0, advance, regulation)
    assert not forced
    assert rung == largest, "the scan must return the LARGEST admissible radius, i.e. the smallest accepting rung"
    assert rung > FULL_RADIUS_RUNG, "the pin is vacuous unless the maximal circle itself is refused"

    assert _bisect_rung(flags) != rung, "a bisection agreed with the scan here, so this state no longer distinguishes them"


def test_a_rung_is_only_admissible_if_it_still_cuts():
    """A loop drawn inside already-swept void complies with every cap and cuts nothing.

    Without the cutting condition the scan would take such a loop -- it engages
    nothing, so it exceeds nothing -- the frontier would not move, and the ladder
    would walk itself down one rung per station. This checks the condition
    directly on a loop deep inside the swept annulus.
    """
    regulation = _regulation(PIN_CAP_DEG)
    advance = (1.0, 0.0)

    # A disk of swept void around the station, and a loop whose whole swept
    # annulus -- centre path plus the tool radius on either side -- fits inside it.
    void_radius = 4.0
    stock = Stock(PIN_POCKET)
    stock.subtract_annulus(PIN_STATION_X, PIN_STATION_Y, 0.0, void_radius)
    idle = _GuideStation(cx=PIN_STATION_X, cy=PIN_STATION_Y, radius=void_radius - 2.0 * TOOL_RADIUS, clockwise=True, tx=1.0, ty=0.0)

    assert _station_is_admissible(stock, idle, advance, regulation.tool_radius, regulation.cap_ratio), "an idle loop trivially satisfies the cap"
    assert not _loop_reaches_material(stock, idle, advance, regulation.tool_radius), "an idle loop must not count as admissible"

    # ... while the maximal circle at the same station, which reaches past the
    # void, does still cut -- so the condition rejects idleness, not smallness.
    biting = _GuideStation(cx=PIN_STATION_X, cy=PIN_STATION_Y, radius=PIN_FULL_RADIUS, clockwise=True, tx=1.0, ty=0.0)
    assert _loop_reaches_material(stock, biting, advance, regulation.tool_radius)


def test_loose_cap_reproduces_the_advance_only_generator():
    """A cap the advance regulation alone meets must come out byte-identical.

    The new knob is additive: where every station's maximal circle is admissible
    the ladder never leaves rung 0, the second pass finds nothing to do, and the
    emitted stream is the one `engagement_controlled_toolpath` produces.
    """
    radial = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, LOOSE_CAP_DEG)
    advance_only = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, LOOSE_CAP_DEG)

    assert _op_signature(radial) == _op_signature(advance_only)


def _guide_maximal_radii(polygon: Polygon):
    """Every clearance-derived maximal loop radius the guide offers on *polygon*.

    Read from the guide the generator itself walks, so the set is complete rather
    than whichever subset one generator's advance happened to land on.
    """
    chains = _guide_chains(polygon, TOOL_DIAMETER, GUIDE_STEP_TOOL_DIAMETERS * TOOL_DIAMETER, RADIAL_CLEARANCE_FRACTION * TOOL_DIAMETER, True, 1000, None)
    return {station.radius for chain in chains for station in chain}


def _over_cap_circles(polygon: Polygon, result: ToolpathResult, cap_deg: float):
    """Radii of the machining circles the audit measures over *cap_deg*, entries aside.

    The audit walks twenty stations per circle against the generator's four, so
    this is a genuinely independent read of what the generator accepted.
    """
    cap_rad = math.radians(cap_deg)
    report = audit_toolpath_engagement(polygon, result, TOOL_DIAMETER, cap_rad)
    entries = _chain_entry_indices(result)
    return [
        float(result.operations[e.op_index].geometry.radius)
        for e in report.operations
        if e.op_index not in entries and e.max_tea > cap_rad and isinstance(result.operations[e.op_index].geometry, Circle)
    ]


def test_tight_cap_regulates_the_radius_and_leaves_fewer_circles_over_the_cap():
    """At a cap the advance alone cannot meet, the ladder must fire and must help.

    Three claims, because none of them alone rules out a degenerate path: radii
    BELOW the guide's maximal ones are emitted (the knob is used); the audit finds
    strictly fewer machining circles over the cap than the advance-only generator
    leaves there; and every circle it still finds over the cap is at a station's
    MAXIMAL radius -- that is, one the ladder could not improve on because no rung
    was admissible there, never a reduced circle the ladder chose badly.
    """
    radial = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)
    advance_only = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)

    maximal = _guide_maximal_radii(POCKET)
    emitted = {float(op.geometry.radius) for op in radial.operations if isinstance(op.geometry, Circle)}
    assert emitted - maximal, "no reduced radius was ever emitted, so the ladder never left rung 0"

    radial_over = _over_cap_circles(POCKET, radial, TIGHT_CAP_DEG)
    advance_over = _over_cap_circles(POCKET, advance_only, TIGHT_CAP_DEG)
    assert advance_over, "baseline is unexpectedly already under the cap, so there is nothing to improve"
    assert len(radial_over) < len(advance_over), f"radial {len(radial_over)} vs advance-only {len(advance_over)} machining circles measured over the cap"

    reduced_over = sorted(radius for radius in radial_over if radius not in maximal)
    assert not reduced_over, f"reduced circles measured over the cap: {reduced_over}"


def test_no_machining_circle_away_from_a_chain_entry_is_a_full_slot():
    """Regression: the ladder must not walk itself down and then jump back to maximal.

    A loop that spins in swept void moves no frontier, so if such loops were
    admissible the ladder would drop one rung per station until it bottomed out and
    the station was forced back to its maximal circle -- by then far from anything
    cleared, i.e. a full slot. Measured on 10x6 at a 40 deg cap that produced
    360 deg circles the advance-only generator never emits. Only a chain's virgin
    stock entry may be a full slot.
    """
    cap_rad = math.radians(TIGHT_CAP_DEG)
    result = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)
    report = audit_toolpath_engagement(POCKET, result, TOOL_DIAMETER, cap_rad)
    entries = _chain_entry_indices(result)

    slots = {
        e.op_index: math.degrees(e.max_tea)
        for e in report.operations
        if e.op_index not in entries and isinstance(result.operations[e.op_index].geometry, Circle) and e.max_tea >= math.radians(180.0)
    }
    assert not slots, f"machining circles at or beyond a half turn away from a chain entry: {slots}"


def test_tight_cap_takes_more_than_one_pass_over_some_chain():
    """A reduced loop leaves an outer band, so the chain must be swept again.

    One plunge per chain means the radius was never reduced anywhere; more plunges
    than chains is the multi-pass mechanism working.
    """
    result = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)

    plunges = sum(1 for op in result.operations if op.operation is OperationType.PLUNGE)
    chains = len({op.path_index for op in result.operations if op.operation is OperationType.PLUNGE})
    assert plunges > chains, f"{plunges} plunge(s) over {chains} chain(s): no chain needed a second pass"

    retracts = sum(1 for op in result.operations if op.operation is OperationType.RETRACT)
    assert plunges == retracts, "every pass must be left by a retract"


def test_radial_path_clears_everything_the_tool_can_reach():
    """Regulating the radius must not buy low engagement by leaving material.

    A smaller loop sweeps a narrower annulus, so the outer band of its station
    survives the pass. What makes that safe is that a station is finished only when
    its MAXIMAL circle has been emitted, so the sweeps keep going until the band is
    open. This checks the consequence directly: replay the path and confirm nothing
    survives that a round tool of this radius could have reached. Residue nearer the
    wall than the tool radius is the round tool's inherent corner residue.
    """
    result = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)
    stock = _replay(result, POCKET)

    width, height = 6.0, 4.0
    columns, rows = 60, 40
    survivors = []
    for row in range(rows):
        y = (row + 0.5) * height / rows
        for column in range(columns):
            x = (column + 0.5) * width / columns
            wall_distance = min(x, width - x, y, height - y)
            if wall_distance > TOOL_RADIUS and stock.contains(x, y):
                survivors.append((round(x, 3), round(y, 3), round(wall_distance, 3)))
    assert not survivors, f"material left further than the tool radius from the wall: {survivors[:10]}"


def test_tighter_cap_is_not_shorter():
    tight = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, 90.0)
    loose = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, 170.0)
    # A tighter cap forces shorter advances and smaller loops, hence more machining
    # circles and more passes over the same guide, hence at least as much cutting
    # travel. Smaller loops clearing less per pass is the physics, not a defect.
    assert _cut_length(tight) >= _cut_length(loose)


def test_generation_is_deterministic():
    first = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)
    second = radius_regulated_toolpath(POCKET, TOOL_DIAMETER, TIGHT_CAP_DEG)
    assert _op_signature(first) == _op_signature(second)


@pytest.mark.parametrize("bad_cap_deg", [0.0, -1.0, 180.001, 360.0, float("nan")])
def test_rejects_out_of_range_cap(bad_cap_deg):
    with pytest.raises(InvalidEngagementCapDegreesError):
        radius_regulated_toolpath(POCKET, TOOL_DIAMETER, bad_cap_deg)


@pytest.mark.parametrize("bad_diameter", [0.0, -2.0])
def test_rejects_non_positive_tool_diameter(bad_diameter):
    with pytest.raises(NonPositiveToolDiameterError):
        radius_regulated_toolpath(POCKET, bad_diameter, LOOSE_CAP_DEG)


def test_rejects_guide_resolution_coarser_than_the_advance_bound():
    with pytest.raises(InvalidGuideResolutionError):
        radius_regulated_toolpath(POCKET, TOOL_DIAMETER, LOOSE_CAP_DEG, guide_step_tool_diameters=1.0)


def test_cap_surrogate_is_the_only_transcendental_crossing():
    # The cap reaches the predicate as the rational chord surrogate 4*sin^2(t/2)
    # and nothing else does; a half turn is the largest legal request and maps to 4.
    assert _cap_surrogate(180.0) == pytest.approx(4.0)
    assert 0.0 < _cap_surrogate(TIGHT_CAP_DEG) < _cap_surrogate(LOOSE_CAP_DEG) < 4.0
