"""Generic-double witness for the station lane, with a wall-clock ceiling.

The station lane's existing coverage centres the tool at (5.0, 5.0) with radius
0.5 on a 10x10 integer square. On fixtures like that every square root is a
perfect square, CORE's approximation error is exactly zero, and the refinement
path is never entered. This file supplies the missing case: coordinates with no
structure, on geometry that actually enters the lane's one-root algebra, with a
bound on the clock so a regression is attributable to this lane rather than
surfacing as a suite that quietly takes an hour.

Vehicle: ``segment_station_cap_exceeded_exact``. It calls ``classify_station_cell``
unconditionally (``segment_oracle.cpp:647``) and only afterwards raises on an
``UNRESOLVED`` decision, so the station work always happens. The full-circle
entry point is not usable here: ``construct_full_circle_cell_authority`` sits
behind the ``cap_exceeded`` early return at ``circle_oracle.cpp:782``, so a
cap-tripping fixture never reaches the station lane at all.

What makes this fixture stressing, measured on this tree (Apple M1 Max, macOS
15.3.2, CPython 3.12.13, 2026-09-08; warm-up run discarded, minimum of N):

| axis held / varied                       | cheap    | expensive | ratio  |
|------------------------------------------|----------|-----------|--------|
| tool circle crosses the void-disk arc     | 0.298 ms | 10.413 ms | 34.9x  |
| same, on an all-integer fixture           | 0.226 ms |  6.452 ms | 28.5x  |
| generic vs integer, both arc-crossing     | 6.452 ms | 10.413 ms |  1.61x |
| generic vs integer, both clear of the arc | 0.226 ms |  0.298 ms |  1.32x |

So genericity alone is worth only ~1.6x on this lane, not the 10^4-10^7x that
``docs/number_types.md`` records elsewhere. The 30x lever is whether the tool
circle at the station **properly crosses the curved stock boundary**, which is
what puts the ``Sqrt_extension`` one-root algebra on the critical path. Both
properties are therefore asserted below: a fixture that loses either one passes
in a third of a millisecond and certifies nothing.

Sweeping the tool radius on a fixed segment reproduces the crossing band to its
endpoints: 0.550 -> 0.336 ms and 0.600 -> 10.477 ms across a predicted lower
tangency at 0.5785; 2.150 -> 9.915 ms and 2.200 -> 0.312 ms across a predicted
upper tangency at 2.1809. The tool radius / half-segment-length ratio, by
contrast, is **not** causal here: holding arc-crossing fixed, the cost is flat
from a ratio of 0.04 to 1.04.
"""

from __future__ import annotations

import math
import time

import numpy as np

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2


# No coordinate here is dyadic with a small denominator, and no two are related
# by a rational scale factor. Digits of pi, e, phi and sqrt(2) are used purely
# because they are structureless, not for any mathematical property.
GENERIC_POCKET = np.array(
    [
        [0.3141, 0.2718, 0.0],
        [9.7183, 0.4142, 0.0],
        [9.4949, 9.6180, 0.0],
        [0.5772, 9.3010, 0.0],
    ],
    dtype=np.float64,
)
GENERIC_DISK_CENTER_X = 4.7312
GENERIC_DISK_CENTER_Y = 5.2891
GENERIC_DISK_RADIUS = 1.3797

GENERIC_START_X = 2.7183
GENERIC_START_Y = 3.1416
GENERIC_END_X = 7.3891
GENERIC_END_Y = 6.2832
GENERIC_STATION_NUMERATOR = 3
GENERIC_STATION_DENOMINATOR = 7
GENERIC_CAP_CHORD_RATIO = 3.1416

# Arc-crossing: the tool circle at the station meets the void disk's boundary in
# two points, so the decision needs one-root algebra. Measured 10.413 ms.
CROSSING_TOOL_RADIUS = 0.9137
# Clear of the arc: the same station and stock, a radius below the lower
# tangency, so the curved boundary never enters the decision. Measured 0.298 ms.
# Present only as the denominator of the separation guard below.
CLEAR_TOOL_RADIUS = 0.3137

# Wall-clock ceiling, seconds. The plan's rule is max(2.0, 4 x slowest observed);
# the slowest of 60 timed runs on this tree was 11.5 ms, so the 2.0 s floor
# dominates. The floor is deliberate: several agents run `pytest -n auto` in this
# worktree concurrently and a tight ceiling would flake on scheduling alone. At
# 10.4 ms measured this leaves ~190x of headroom, which the separation guard
# below covers from the other side.
STATION_SECONDS_CEILING = 2.0

# Proper-crossing margin, as a fraction of the distance to each tangency. The
# measured transition is sharp -- the cost flips between radii 0.055 apart at the
# lower tangency -- so the guard demands the fixture sit 20% clear of both
# tangencies rather than merely inside the band. The fixture's actual margins are
# 1.58x the lower tangency and 0.42x the upper.
CROSSING_TANGENCY_MARGIN = 0.2

# Repeats behind each timing used in the separation guard. The reported figure is
# the MINIMUM: scheduling noise on a loaded machine only ever adds time, so a
# minimum over a handful of runs is robust where a single run is not.
SEPARATION_TIMING_REPEATS = 5

# The arc-crossing fixture must cost at least this many times the arc-clear one.
# Measured 28.95x-34.76x over eight min-of-5 trials, so a floor of 5x carries
# 5.8x headroom against the worst trial while still failing loudly if the lane
# stops doing the exact work -- a short-circuit collapses the ratio towards 1x,
# which the ceiling alone can never detect because a wrong answer arrives fast.
MINIMUM_CROSSING_SEPARATION = 5.0

# A binary64 is exactly a dyadic rational. Multiplying by 2**DYADIC_GUARD_BITS
# and asking for a non-integer rejects any coordinate expressible with that few
# fractional bits, which is what an "integer-ish" fixture edit would produce.
DYADIC_GUARD_BITS = 8


def _generic_stock() -> _stock_2.Stock2:
    stock = _stock_2.Stock2(GENERIC_POCKET, [])
    stock.subtract_disk(
        GENERIC_DISK_CENTER_X,
        GENERIC_DISK_CENTER_Y,
        GENERIC_DISK_RADIUS,
    )
    return stock


def _station_gap() -> float:
    """Distance from the interpolated station to the void disk's centre."""
    parameter = GENERIC_STATION_NUMERATOR / GENERIC_STATION_DENOMINATOR
    station_x = GENERIC_START_X + parameter * (GENERIC_END_X - GENERIC_START_X)
    station_y = GENERIC_START_Y + parameter * (GENERIC_END_Y - GENERIC_START_Y)
    return math.hypot(
        station_x - GENERIC_DISK_CENTER_X,
        station_y - GENERIC_DISK_CENTER_Y,
    )


def _tangency_radii() -> tuple[float, float]:
    """The two tool radii at which the tool circle is tangent to the void disk.

    Returns:
        The internal and external tangency radii. Between them the two circles
        cross in two points; outside them one contains the other or they are
        disjoint, and the curved boundary drops out of the decision.
    """
    gap = _station_gap()
    return abs(GENERIC_DISK_RADIUS - gap), GENERIC_DISK_RADIUS + gap


def _decide(tool_radius: float) -> bool:
    return _continuous_tea_2.segment_station_cap_exceeded_exact(
        _generic_stock(),
        GENERIC_START_X,
        GENERIC_START_Y,
        GENERIC_END_X,
        GENERIC_END_Y,
        GENERIC_STATION_NUMERATOR,
        GENERIC_STATION_DENOMINATOR,
        tool_radius,
        GENERIC_CAP_CHORD_RATIO,
    )


def _fastest_decision_seconds(tool_radius: float) -> float:
    """Time the station decision, discarding a warm-up run.

    The stock is rebuilt outside every timed region, so only the station lane is
    measured.
    """
    _decide(tool_radius)
    timings = []
    for _ in range(SEPARATION_TIMING_REPEATS):
        stock = _generic_stock()
        start = time.perf_counter()
        _continuous_tea_2.segment_station_cap_exceeded_exact(
            stock,
            GENERIC_START_X,
            GENERIC_START_Y,
            GENERIC_END_X,
            GENERIC_END_Y,
            GENERIC_STATION_NUMERATOR,
            GENERIC_STATION_DENOMINATOR,
            tool_radius,
            GENERIC_CAP_CHORD_RATIO,
        )
        timings.append(time.perf_counter() - start)
    return min(timings)


def test_segment_station_decides_on_generic_doubles_within_budget() -> None:
    stock = _generic_stock()

    start = time.perf_counter()
    exceeded = _continuous_tea_2.segment_station_cap_exceeded_exact(
        stock,
        GENERIC_START_X,
        GENERIC_START_Y,
        GENERIC_END_X,
        GENERIC_END_Y,
        GENERIC_STATION_NUMERATOR,
        GENERIC_STATION_DENOMINATOR,
        CROSSING_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    )
    elapsed = time.perf_counter() - start

    # The lane reached a definite disposition. An UNRESOLVED station raises
    # IncompleteSegmentOracleError instead of returning, so a bool IS the proof
    # that classify_station_cell decided rather than gave up.
    assert exceeded is True
    assert elapsed < STATION_SECONDS_CEILING, f"segment station decision took {elapsed:.3f}s on generic doubles, ceiling {STATION_SECONDS_CEILING:.1f}s"


def test_the_witness_fixture_properly_crosses_the_curved_boundary() -> None:
    """Guard the guard: the geometry that makes this fixture stressing.

    Measured on this tree, the arc-crossing/arc-clear split is 34.9x while
    genericity alone is 1.61x. A later edit that slid the tool radius past either
    tangency would leave every other assertion in this file green while the
    one-root algebra stopped running -- exactly the failure this witness exists
    to prevent, one level up.
    """
    internal_tangency, external_tangency = _tangency_radii()

    assert CROSSING_TOOL_RADIUS > (1.0 + CROSSING_TANGENCY_MARGIN) * internal_tangency, (
        f"tool radius {CROSSING_TOOL_RADIUS} is within "
        f"{CROSSING_TANGENCY_MARGIN:.0%} of internal tangency {internal_tangency:.4f}; "
        "the tool circle no longer properly crosses the void-disk arc"
    )
    assert CROSSING_TOOL_RADIUS < (1.0 - CROSSING_TANGENCY_MARGIN) * external_tangency, (
        f"tool radius {CROSSING_TOOL_RADIUS} is within "
        f"{CROSSING_TANGENCY_MARGIN:.0%} of external tangency {external_tangency:.4f}; "
        "the tool circle no longer properly crosses the void-disk arc"
    )
    assert CLEAR_TOOL_RADIUS < internal_tangency, (
        f"the separation guard's control radius {CLEAR_TOOL_RADIUS} crosses the arc (internal tangency {internal_tangency:.4f}); it is no longer a control"
    )


def test_the_witness_fixture_is_not_an_integer_fixture() -> None:
    """Guard the guard: a later edit must not quietly round the fixture off.

    An integer or small-dyadic fixture makes every square root a perfect square,
    so CORE's approximation error is exactly zero and the refinement path is
    never entered.
    """
    coordinates = [
        *GENERIC_POCKET[:, 0].tolist(),
        *GENERIC_POCKET[:, 1].tolist(),
        GENERIC_DISK_CENTER_X,
        GENERIC_DISK_CENTER_Y,
        GENERIC_DISK_RADIUS,
        GENERIC_START_X,
        GENERIC_START_Y,
        GENERIC_END_X,
        GENERIC_END_Y,
        CROSSING_TOOL_RADIUS,
        CLEAR_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    ]
    scale = float(2**DYADIC_GUARD_BITS)
    for value in coordinates:
        assert value != round(value), f"{value} is an integer fixture"
        assert (value * scale) != round(value * scale), f"{value} is a dyadic fixture with at most {DYADIC_GUARD_BITS} fractional bits"


def test_the_crossing_fixture_costs_what_the_exact_algebra_costs() -> None:
    """The lane did the work, not merely finished in time.

    A ceiling can only fail upwards. If a change short-circuits the station
    decision -- returning fast and wrong -- every other assertion here still
    passes. This one fails instead: the arc-crossing fixture must remain
    materially more expensive than the same station with the curved boundary out
    of reach.
    """
    crossing_seconds = _fastest_decision_seconds(CROSSING_TOOL_RADIUS)
    clear_seconds = _fastest_decision_seconds(CLEAR_TOOL_RADIUS)

    separation = crossing_seconds / clear_seconds
    assert separation > MINIMUM_CROSSING_SEPARATION, (
        f"the arc-crossing station cost {crossing_seconds * 1e3:.3f} ms against "
        f"{clear_seconds * 1e3:.3f} ms clear of the arc, a separation of "
        f"{separation:.1f}x below the {MINIMUM_CROSSING_SEPARATION:.0f}x floor; "
        "the one-root algebra is no longer on the critical path"
    )
