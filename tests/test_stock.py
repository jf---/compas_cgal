import math
import multiprocessing
import random
from decimal import Decimal, getcontext
from fractions import Fraction
from multiprocessing.connection import Connection

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
# subtract_annulus: the exact swept region of a full circular guide           #
# --------------------------------------------------------------------------- #

# Full-circle guide used by every annulus test: tool radius 0.5 swept about the
# circle of radius 2 centred at (5, 5). Both radii and both of their sums are
# exact in binary floating point (2 - 0.5 = 1.5, 2 + 0.5 = 2.5), so the double
# API and `subtract_arc_sweep` -- which forms the same two bounds in exact
# rational arithmetic -- describe the identical region, with no rounding to
# reason about. Radii whose sum is NOT exactly representable would differ by an
# ulp, which is a property of the caller's arithmetic, not of the annulus.
ANNULUS_CENTER = (5.0, 5.0)
ANNULUS_GUIDE_RADIUS = 2.0
ANNULUS_TOOL_RADIUS = 0.5
ANNULUS_INNER = ANNULUS_GUIDE_RADIUS - ANNULUS_TOOL_RADIUS
ANNULUS_OUTER = ANNULUS_GUIDE_RADIUS + ANNULUS_TOOL_RADIUS

# Under-coverage of the disk chain, from its own construction in src/stock_2.cpp:
# spacing s = 2*r*sqrt(CHAIN_SLACK_FRACTION) gives a sagitta shortfall of
# s^2/(4r) = r*CHAIN_SLACK_FRACTION. At r = 0.5 that is 5.0e-5 mm.
CHAIN_SLACK_FRACTION = 1e-4
CHAIN_SAGITTA = ANNULUS_TOOL_RADIUS * CHAIN_SLACK_FRACTION
# Ring probes per turn: comfortably denser than the 1260 disks the chain places
# over the full turn, so probes land between disks rather than on them.
CHAIN_PROBE_COUNT = 4001

# Exact-centred reference chain: chord bound and center budget for
# `subtract_exact_full_circle`, which parametrises the guide circle rationally.
EXACT_CHAIN_MAX_CHORD = 0.1
EXACT_CHAIN_CENTER_LIMIT = 4096


def _annulus_stock():
    """Virgin square with the exact swept annulus of the reference motion removed."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_annulus(*ANNULUS_CENTER, ANNULUS_INNER, ANNULUS_OUTER)
    return stock


def _chain_stock():
    """Virgin square with the same full turn removed by the under-covering disk chain.

    Driven as two clockwise half arcs rather than one full circle: the full-circle
    branch now takes the exact annulus, so this is how the still-live chain path is
    exercised over the identical region.
    """
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 3.0, 5.0, True, ANNULUS_TOOL_RADIUS)  # lower half
    stock.subtract_arc_sweep(5.0, 5.0, 3.0, 5.0, 7.0, 5.0, True, ANNULUS_TOOL_RADIUS)  # upper half
    return stock


def _full_circle_arc_sweep_stock():
    """Virgin square depleted through the public full-circle `subtract_arc_sweep`."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, ANNULUS_TOOL_RADIUS)
    return stock


def test_subtract_annulus_removes_the_band_and_nothing_else():
    stock = _annulus_stock()
    assert not stock.contains(7.0, 5.0)  # mid-band, on the guide circle
    assert not stock.contains(5.0, 7.0)  # mid-band, quarter turn away
    assert not stock.contains(6.55, 5.0)  # just inside the inner bound
    assert not stock.contains(7.45, 5.0)  # just inside the outer bound
    assert stock.contains(5.0, 5.0)  # island inside the hole survives
    assert stock.contains(6.45, 5.0)  # strictly inside inner_radius
    assert stock.contains(7.55, 5.0)  # strictly outside outer_radius


def test_annulus_removes_the_slivers_the_disk_chain_retains():
    """The chain UNDER-covers the swept region; the annulus IS the swept region.

    The chain's shortfall is the sagitta between consecutive tool disks,
    `CHAIN_SAGITTA` = 5.0e-5 mm here -- material the model kept although the tool
    had in fact removed it, which reads high on every later engagement query.
    Probing a ring at half that depth inside the outer bound, the chain leaves
    material at the mid-angles between its disks and the annulus leaves none.

    Note this is deliberately NOT phrased as exact set inclusion of the annulus
    result in the chain result. That inclusion is FALSE, and not because of the
    annulus: the chain's centres are `cx + guide_r*cos(a)` in double arithmetic
    and land up to 7.0e-16 mm OFF the guide circle (measured over these 1260
    disks), so the chain also removes a sub-femtometre crescent outside the true
    swept region. Neither set contains the other; what is true, and what matters
    at 11 orders of magnitude larger, is the shortfall asserted here.
    `test_annulus_is_the_certified_full_circle_sweep_region` states the exact
    inclusion, against a chain whose centres ARE exactly on the guide circle.
    """
    annulus = _annulus_stock()
    chain = _chain_stock()
    probe_radius = ANNULUS_OUTER - 0.5 * CHAIN_SAGITTA
    angles = [2.0 * math.pi * i / CHAIN_PROBE_COUNT for i in range(CHAIN_PROBE_COUNT)]
    ring = [(5.0 + probe_radius * math.cos(a), 5.0 + probe_radius * math.sin(a)) for a in angles]
    assert not any(annulus.contains(x, y) for x, y in ring)
    assert any(chain.contains(x, y) for x, y in ring)


def test_annulus_is_the_certified_full_circle_sweep_region():
    """The wired-in full-turn path removes exactly what the depletion oracle certifies.

    `subtract_exact_full_circle` places its tool disks at EXACT rational points on
    the guide circle and is proved to under-cover the swept annulus
    (`exact_full_circle_undercover_holds`). Against that exactly-centred chain the
    inclusion does hold exactly, in the direction the certificate states: the
    annulus removes everything the chain removes, and strictly more.
    """
    swept = _full_circle_arc_sweep_stock()
    exact_chain = _stock_2.Stock2(SQUARE, [])
    exact_chain.subtract_exact_full_circle(
        *ANNULUS_CENTER,
        ANNULUS_GUIDE_RADIUS,
        0.0,
        False,
        ANNULUS_TOOL_RADIUS,
        EXACT_CHAIN_MAX_CHORD,
        EXACT_CHAIN_CENTER_LIMIT,
    )
    assert swept.is_subset_of(exact_chain)
    assert not exact_chain.is_subset_of(swept)


def test_zero_inner_radius_is_exactly_a_disk():
    annulus = _stock_2.Stock2(SQUARE, [])
    annulus.subtract_annulus(*ANNULUS_CENTER, 0.0, ANNULUS_OUTER)
    disk = _stock_2.Stock2(SQUARE, [])
    disk.subtract_disk(*ANNULUS_CENTER, ANNULUS_OUTER)
    assert annulus.exactly_equals(disk)


# Floor for the arrangement-size win, well under the measured ratio so it
# ratchets without tracking CGAL's exact arrangement bookkeeping.
ANNULUS_VERTEX_REDUCTION_FLOOR = 100


def test_annulus_arrangement_is_orders_smaller_than_the_chain():
    """Two boundary circles instead of a 1260-disk union, in arrangement features.

    Measured: 8 vertices after the annulus against 2536 after the chain, a 317x
    reduction. Arrangement size is what makes every later `engagement_at` query
    expensive, so this is the number the whole change exists to move.
    """
    annulus_vertices = _annulus_stock().arrangement_stats()[0]
    chain_vertices = _chain_stock().arrangement_stats()[0]
    assert annulus_vertices * ANNULUS_VERTEX_REDUCTION_FLOOR <= chain_vertices, (annulus_vertices, chain_vertices)


def test_full_circle_arc_sweep_uses_the_exact_annulus():
    """The wired-in path: a full-circle guide removes the annulus, not a disk chain.

    Exact equality is available here only because ``2 +/- 0.5`` is exact in binary
    floating point, so the double radii this test passes and the exact rational
    bounds `subtract_arc_sweep` forms are the same two numbers.
    """
    assert _full_circle_arc_sweep_stock().exactly_equals(_annulus_stock())


def test_partial_arc_sweep_still_uses_the_disk_chain():
    """A quarter arc is untouched by this change: it sweeps one quadrant, not a ring."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 5.0, 7.0, False, ANNULUS_TOOL_RADIUS)
    assert not stock.contains(7.0, 5.0)  # arc start is swept
    # Mid-band in the opposite quadrant: material iff the annulus was NOT taken.
    opposite = 5.0 - ANNULUS_GUIDE_RADIUS / math.sqrt(2.0)
    assert stock.contains(opposite, opposite)
    assert stock.contains(5.0, 3.0)  # mid-band, quarter turn back from the start
    assert stock.contains(3.0, 5.0)  # mid-band, quarter turn past the end


@pytest.mark.parametrize(
    "inner, outer",
    [
        (2.5, 2.5),  # empty annulus
        (2.5, 1.5),  # inverted
        (-0.1, 2.5),  # negative inner
        (1.5, 0.0),  # non-positive outer
    ],
)
def test_annulus_rejects_malformed_radii(inner, outer):
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(_stock_2.InvalidAnnulusRadiiError):
        stock.subtract_annulus(*ANNULUS_CENTER, inner, outer)


@pytest.mark.parametrize("bad", [math.inf, -math.inf, math.nan])
def test_annulus_rejects_non_finite_input(bad):
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(_stock_2.NonFiniteAnnulusInputError):
        stock.subtract_annulus(5.0, 5.0, ANNULUS_INNER, bad)
    with pytest.raises(_stock_2.NonFiniteAnnulusInputError):
        stock.subtract_annulus(bad, 5.0, ANNULUS_INNER, ANNULUS_OUTER)


def test_stock_wrapper_exposes_subtract_annulus():
    from compas.geometry import Polygon

    from compas_cgal.stock import Stock

    stock = Stock(Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]))
    stock.subtract_annulus(*ANNULUS_CENTER, ANNULUS_INNER, ANNULUS_OUTER)
    assert not stock.contains(7.0, 5.0)
    assert stock.contains(5.0, 5.0)


# --------------------------------------------------------------------------- #
# engagement_at: exact station TEA                                            #
# --------------------------------------------------------------------------- #

# Fixed 18-vertex arrangement recovered from the scale-invariance property's
# minimal hanging example. A healthy exact-station query on this tiny input is
# sub-second; ten seconds leaves ample cold-process headroom while converting a
# native GIL-holding regression into a finite, attributable test failure.
EXACT_STATION_QUERY_BUDGET_SECONDS = 10.0
PATHOLOGICAL_STOCK_BOUNDARY = np.array(
    [[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]],
    dtype=np.float64,
)
PATHOLOGICAL_STATION_Y = -0.375
PATHOLOGICAL_TOOL_RADIUS = 1.0
PATHOLOGICAL_CIRCLES = (
    (-3.375, 1.0),
    (-2.5, 2.5),
    (-1.625, 2.25),
)
PATHOLOGICAL_MATERIAL_X = 2.5
PATHOLOGICAL_VOID_POINT = (0.0, 0.0)
PATHOLOGICAL_CAP_DEGREES = 120.0
PATHOLOGICAL_LOCAL_DISK_RADIUS = 0.125


def _pathological_stock():
    """Build the multi-outer-CCB stock recovered from the native stall."""
    stock = _stock_2.Stock2(PATHOLOGICAL_STOCK_BOUNDARY, [])
    first_x = PATHOLOGICAL_CIRCLES[0][0]
    stock.subtract_disk(
        first_x,
        PATHOLOGICAL_STATION_Y,
        PATHOLOGICAL_TOOL_RADIUS,
    )
    for index, (center_x, guide_radius) in enumerate(PATHOLOGICAL_CIRCLES):
        rim_x = center_x + guide_radius
        stock.subtract_arc_sweep(
            center_x,
            PATHOLOGICAL_STATION_Y,
            rim_x,
            PATHOLOGICAL_STATION_Y,
            rim_x,
            PATHOLOGICAL_STATION_Y,
            False,
            PATHOLOGICAL_TOOL_RADIUS,
        )
        if index + 1 < len(PATHOLOGICAL_CIRCLES):
            next_x = PATHOLOGICAL_CIRCLES[index + 1][0]
            stock.subtract_capsule(
                center_x,
                PATHOLOGICAL_STATION_Y,
                next_x,
                PATHOLOGICAL_STATION_Y,
                PATHOLOGICAL_TOOL_RADIUS,
            )
    return stock


def _query_pathological_exact_station(sender: Connection) -> None:
    """Query the original zero-engagement station in a killable process."""
    stock = _pathological_stock()

    query_x = PATHOLOGICAL_CIRCLES[-1][0]
    result = _stock_2.engagement_at(
        stock,
        query_x,
        PATHOLOGICAL_STATION_Y,
        PATHOLOGICAL_TOOL_RADIUS,
        cap_ratio(math.radians(PATHOLOGICAL_CAP_DEGREES)),
    )
    sender.send(result)
    sender.close()


def _query_pathological_overlap_station(sender: Connection) -> None:
    """Query a cutter whose rim is mostly in material on the pathological stock."""
    stock = _pathological_stock()
    result = _stock_2.engagement_at(
        stock,
        PATHOLOGICAL_MATERIAL_X,
        PATHOLOGICAL_STATION_Y,
        PATHOLOGICAL_TOOL_RADIUS,
        cap_ratio(math.radians(120.0)),
    )
    sender.send(result)
    sender.close()


def _query_pathological_contains(sender: Connection) -> None:
    """Query a point strictly inside the removed annulus."""
    stock = _pathological_stock()
    sender.send(stock.contains(*PATHOLOGICAL_VOID_POINT))
    sender.close()


def _query_after_stock_mutation(sender: Connection) -> None:
    """Build a locator, mutate its arrangement, then query the new epoch."""
    stock = _stock_2.Stock2(SQUARE, [])
    before = stock.contains(5.0, 5.0)
    stock.subtract_disk(5.0, 5.0, 1.0)
    sender.send((before, stock.contains(5.0, 5.0)))
    sender.close()


def _query_pathological_local_depletion(sender: Connection) -> None:
    """Compare local and global removal inside the pathological void face."""
    reference = _pathological_stock()
    local = _pathological_stock()
    reference.subtract_disk(*PATHOLOGICAL_VOID_POINT, PATHOLOGICAL_LOCAL_DISK_RADIUS)
    local.subtract_disk_local(*PATHOLOGICAL_VOID_POINT, PATHOLOGICAL_LOCAL_DISK_RADIUS)
    sender.send(
        (
            local.exactly_equals(reference),
            reference.exactly_equals(local),
            local.representation_is_valid(),
        )
    )
    sender.close()


def _bounded_pathological_query(target):
    """Return one native result or fail if the child stalls."""
    context = multiprocessing.get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(
        target=target,
        args=(sender,),
    )
    process.start()
    sender.close()
    try:
        process.join(EXACT_STATION_QUERY_BUDGET_SECONDS)
        if process.is_alive():
            process.terminate()
            process.join()
            pytest.fail(f"exact stock query did not terminate within {EXACT_STATION_QUERY_BUDGET_SECONDS:.0f} seconds")
        assert process.exitcode == 0
        assert receiver.poll(), "exact stock child returned no result"
        return receiver.recv()
    finally:
        if process.is_alive():
            process.terminate()
            process.join()
        receiver.close()


def test_exact_station_query_terminates_on_depleted_dyadic_stock():
    """The original zero-engagement station must terminate and remain clear."""
    assert _bounded_pathological_query(_query_pathological_exact_station) == (
        0.0,
        0.0,
        False,
    )


def test_exact_station_reports_material_on_multi_outer_ccb_stock():
    """A locator must not certify zero when the cutter rim crosses material."""
    total, max_run, exceeded = _bounded_pathological_query(_query_pathological_overlap_station)
    assert total > math.pi
    assert max_run > math.radians(PATHOLOGICAL_CAP_DEGREES)
    assert exceeded


def test_contains_returns_void_on_multi_outer_ccb_stock():
    """Point containment must terminate and classify a strict annulus point."""
    assert _bounded_pathological_query(_query_pathological_contains) is False


def test_point_location_rebuilds_after_stock_mutation():
    """A query accelerator from an earlier arrangement must be detached."""
    assert _bounded_pathological_query(_query_after_stock_mutation) == (True, False)


def test_local_depletion_terminates_on_multi_outer_ccb_stock():
    """Local removal must preserve the reference set on a multi-outer-CCB face."""
    assert _bounded_pathological_query(_query_pathological_local_depletion) == (
        True,
        True,
        True,
    )


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


@pytest.mark.parametrize("invalid", [math.nan, math.inf, -math.inf])
def test_engagement_native_boundary_raises_named_errors(invalid):
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(_stock_2.NonFiniteEngagementInputError):
        _stock_2.engagement_at(stock, invalid, 5.0, 0.5, 2.0)
    with pytest.raises(_stock_2.NonPositiveEngagementToolRadiusError):
        _stock_2.engagement_at(stock, 5.0, 5.0, 0.0, 2.0)
    with pytest.raises(_stock_2.InvalidEngagementCapRatioError):
        _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, 0.0)
    with pytest.raises(_stock_2.InvalidEngagementGapRatioError):
        _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, 2.0, 5.0)


@pytest.mark.parametrize("invalid", [math.nan, math.inf, -math.inf])
def test_segment_and_radical_boundaries_raise_named_errors(invalid):
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(_stock_2.NonFiniteSegmentCertificationInputError):
        _stock_2.certify_segment_tea(stock, invalid, 0.0, 1.0, 0.0, 0.5, 1.0)
    with pytest.raises(_stock_2.NonPositiveSegmentCertificationToolRadiusError):
        _stock_2.certify_segment_tea(stock, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0)
    with pytest.raises(_stock_2.InvalidSegmentCertificationCapError):
        _stock_2.certify_segment_tea(stock, 0.0, 0.0, 1.0, 0.0, 0.5, 0.0)
    with pytest.raises(_stock_2.NonFiniteMixedRadicalInputError):
        _stock_2._sign_mixed_radical(invalid, 0.0, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(_stock_2.InvalidMixedRadicalRootError):
        _stock_2._sign_mixed_radical(0.0, 0.0, 0.0, 0.0, -1.0, 0.0)


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


# --------------------------------------------------------------------------- #
# certify_segment_tea: exact-station TEA cap certificate along a linear motion #
# --------------------------------------------------------------------------- #


def test_certify_segment_in_open_field():
    """Rim never touches material along the motion -> TEA == 0, cap certified."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(1.0, 5.0, 9.0, 5.0, 1.5)  # wide cleared corridor
    max_tea, ok, n = _stock_2.certify_segment_tea(stock, 2.0, 5.0, 8.0, 5.0, 0.5, math.radians(120))
    assert ok
    assert max_tea == pytest.approx(0.0, abs=1e-9)


def test_certify_segment_flags_slotting():
    """Virgin slotting move: TEA ~ a full turn >> 120 deg cap -> not certified."""
    stock = _stock_2.Stock2(SQUARE, [])
    max_tea, ok, n = _stock_2.certify_segment_tea(stock, 2.0, 5.0, 8.0, 5.0, 0.5, math.radians(120))
    assert not ok
    assert max_tea > math.radians(170)


def test_certify_stations_refine_near_feature():
    """Motion leaving a cleared pocket into material forces refinement + fails."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 2.0)  # cleared pocket in the middle
    max_tea, ok, n = _stock_2.certify_segment_tea(stock, 5.0, 5.0, 9.0, 5.0, 0.5, math.radians(120))
    assert not ok  # exits the cleared disk into full slotting
    assert n > 2  # refinement actually happened


def test_certify_flags_unclosable_guard_margin():
    """Cap set barely above a near-uniform corridor's station TEA: no station
    exceeds the cap, yet the analytic guard margin (>= ~0.13 rad down at the
    spacing floor) can never be closed, so the certificate conservatively reports
    False -- the guarded-threshold discipline refusing a false pass."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(0.0, 6.7, 10.0, 6.7, 2.0)  # clears band y in [4.7, 8.7]
    r, y = 0.5, 5.0
    # self-calibrate: the largest station TEA the motion's dyadic stations hit
    tau = max(_stock_2.engagement_at(stock, x, y, r, cap_ratio(math.pi))[1] for x in (3.0, 4.0, 5.0, 6.0, 7.0))
    cap = tau + 0.02  # barely above every station's TEA
    assert cap <= math.pi
    _, ok, _ = _stock_2.certify_segment_tea(stock, 3.0, y, 7.0, y, r, cap)
    assert not ok


# --------------------------------------------------------------------------- #
# gap-closure pessimism: the certificate's run-merge repair                    #
# --------------------------------------------------------------------------- #

# Stock = lower half-plane (y <= 5) with a thin vertical void slot at x = 5 that
# bites the cutter rim's bottom. A cutter of radius 0.5 centred at cy = 5 engages
# the lower semicircle; the slot splits that engaged region into two maximal runs
# separated by a thin void gap. Sliding the centre in x slides the slot off the
# rim -- the two runs MERGE into one ~pi run when it no longer bites -- so max_run
# jumps discontinuously by O(1) within an arbitrarily small step. That jump is the
# configuration the O(sqrt d) growth lemma cannot bridge (demonstrated in
# tests/test_engagement_oracle.py::test_merge_jump_exceeds_growth_lemma) and that
# gap-closure pessimism repairs.
LOWER_HALF = np.array([[0, 0, 0], [10, 0, 0], [10, 5, 0], [0, 5, 0]], dtype=np.float64)


def _slotted_lower_half():
    """Lower half-plane with a thin central void slot reaching the rim bottom."""
    stock = _stock_2.Stock2(LOWER_HALF, [])
    stock.subtract_capsule(5.0, 4.30, 5.0, 4.53, 0.05)
    return stock


def _growth_bound(d, r):
    """Raw factor-1 TEA-growth lemma (mirrors tea_growth_bound in engagement_2.cpp)."""
    return 4.0 * math.asin(min(1.0, d / (2.0 * r))) + 2.0 * math.acos(max(-1.0, 1.0 - d / r))


def test_gap_closure_flips_cap_decision():
    """Pessimism unit test: absorbing a thin void gap between two engaged runs flips
    the exact cap decision.

    At cx = 5 the central slot splits the lower semicircle into two symmetric ~84
    deg runs separated by a ~11 deg void gap. With gap_close_ratio = 0 the cap
    decision sees only the larger single run (~84 deg); closing the gap (gamma = 30
    deg, comfortably above the ~11 deg gap) absorbs it into one ~180 deg pessimistic
    run. For a cap BETWEEN the single-run and merged spans the verdict flips
    False -> True -- the exact mechanism that lets a station account for a merge that
    could complete before the next station. Reported max_run stays the TRUE (unclosed)
    measure regardless of gamma (the deciding/reporting split).
    """
    stock = _slotted_lower_half()
    r, cx, cy = 0.5, 5.0, 5.0
    total, max_run, _ = _stock_2.engagement_at(stock, cx, cy, r, cap_ratio(math.pi))
    assert total == pytest.approx(2.0 * max_run, rel=0.05)  # two symmetric runs
    assert max_run < math.radians(90.0)

    cap = math.radians(135.0)  # between the single run (~84) and the merged span (~180)
    assert max_run < cap
    _, max_open, exceeded_open = _stock_2.engagement_at(stock, cx, cy, r, cap_ratio(cap), 0.0)
    assert not exceeded_open  # gap open: largest run ~84 deg < 135 deg cap

    gamma = cap_ratio(math.radians(30.0))  # 30 deg > the ~11 deg gap -> absorbed
    _, max_closed, exceeded_closed = _stock_2.engagement_at(stock, cx, cy, r, cap_ratio(cap), gamma)
    assert exceeded_closed  # gap closed: merged ~180 deg pessimistic run > 135 deg cap
    # deciding/reporting split: gamma changes only the DECISION, never reported max_run
    assert max_closed == pytest.approx(max_open, abs=1e-9)


def test_certificate_refuses_merge_over_cap():
    """Load-bearing soundness: a linear motion crossing the run-merge locus with a cap
    below the merged span is NOT certified.

    The cutter centre sweeps in x across the slot's disengagement (cx ~ 4.71), where
    the two engaged runs fuse into one ~180 deg (== pi) run. A station in the merged
    stretch genuinely exceeds a 150 deg cap while a split station further along is
    under it, so the guarded stations straddle the O(1) merge jump. Gap-closure
    pessimism pre-absorbs the closing gap at the stations, so the merge is seen and
    certify_segment_tea reports the motion uncertified -- the certificate the pre-repair
    growth lemma could unsoundly pass (see the jump vs bound in
    test_engagement_oracle.py::test_merge_jump_exceeds_growth_lemma).
    """
    stock = _slotted_lower_half()
    r, cap = 0.5, math.radians(150.0)
    # the merged stretch (cx <= ~4.71) holds a full pi run -- a genuine > 150 deg violation
    _, merged_run, _ = _stock_2.engagement_at(stock, 4.68, 5.0, r, cap_ratio(math.pi))
    assert merged_run > cap
    # a split station further along is under the cap (two runs, largest < 150 deg)
    _, split_run, _ = _stock_2.engagement_at(stock, 4.85, 5.0, r, cap_ratio(math.pi))
    assert split_run < cap

    _, certified, stations = _stock_2.certify_segment_tea(stock, 4.65, 5.0, 4.85, 5.0, r, cap)
    assert not certified
    assert stations > 2  # refinement engaged, not a single-shot verdict


def test_gap_close_ratio_zero_preserves_true_verdict():
    """gap_close_ratio = 0 (the default) closes no gap: the cap decision is the true
    per-run verdict, bit-for-bit the pre-pessimism result. Explicit 0.0 and the default
    must agree, and both must see only the true ~84 deg runs (not exceeded at 135 deg)."""
    stock = _slotted_lower_half()
    r, cap = 0.5, cap_ratio(math.radians(135.0))
    default = _stock_2.engagement_at(stock, 5.0, 5.0, r, cap)
    explicit_zero = _stock_2.engagement_at(stock, 5.0, 5.0, r, cap, 0.0)
    assert default == explicit_zero
    assert not default[2]  # true largest run ~84 deg < 135 deg cap


def test_gap_close_ratio_rejects_out_of_range():
    """gap_close_ratio is the chord surrogate 4*sin^2(gamma/2) for 0 <= gamma <= pi, so
    it must lie in [0, 4]; values outside raise (validated before exact injection)."""
    stock = _slotted_lower_half()
    r, cap = 0.5, cap_ratio(math.pi)
    with pytest.raises(ValueError):
        _stock_2.engagement_at(stock, 5.0, 5.0, r, cap, -1e-9)
    with pytest.raises(ValueError):
        _stock_2.engagement_at(stock, 5.0, 5.0, r, cap, 4.0 + 1e-9)


def _thin_slotted_lower_half():
    """Lower half-plane with a THIN central void slot: the two engaged runs are separated
    by a ~4.6 deg void gap.

    The gap is deliberately below the certifier's floor-level guard
    gamma_floor = 2*GROWTH(hs_floor) ~ 7.4 deg (hs_floor = STATION_FLOOR_FRACTION * r / 2 =
    2.5e-4 at r = 0.5), so gap-closure pessimism keeps the two runs merged at EVERY
    refinement level down to the floor -- the certifier can never resolve the gap as
    genuinely separate, and the guarded station test refuses the motion.
    """
    stock = _stock_2.Stock2(LOWER_HALF, [])
    stock.subtract_capsule(5.0, 4.44, 5.0, 4.53, 0.02)
    return stock


def test_certifier_gap_wiring_regression_guard():
    """Regression guard: the certifier's gap-closure gamma-WIRING must not be silently
    removed.

    Reverting ``gap_close_ratio -> 0`` in ``certify_segment_tea``'s deciding branch leaves
    every other test green (test_certificate_refuses_merge_over_cap refuses at gamma = 0
    anyway, via the guard margin), so nothing else catches a wiring removal. This test
    closes that gap by STATION-ANCHORING the merge witness: the motion is centred so its
    MIDPOINT cutter-centre (5, 5) -- ALWAYS a station of the dyadic refinement's first
    bisection -- sits on the split witness (two ~88 deg runs + a ~4.6 deg void gap).

    With the wiring present, every station's PESSIMISTIC run absorbs the thin gap (below
    gamma_floor, so closed at every level to the floor) and reads ~pi > the guarded cap, so
    the motion is refused. Remove the wiring and every station reads only its true ~88 deg
    run, both endpoints clear the root guarded cap, the root early-returns, and the motion
    is FALSELY certified -- flipping this assertion.

    BRANCH: this witness is the accepted SUB-FLOOR / OVER-CONSERVATIVE-refusal branch. The
    ~4.6 deg gap is below the refinement floor's guard resolution, so the certifier cannot
    rule out a merge and conservatively refuses (the safe direction). The motion is in fact
    merge-free over [4.99, 5.01] -- all TRUE max_run ~88 deg < 170 deg cap -- so the guard is
    on the WIRING, not on a genuine engagement violation. Ablation confirmed 2026-07-23:
    reverting the wiring flips certify to True (this test the sole failure); wiring restored,
    suite green.
    """
    stock = _thin_slotted_lower_half()
    r, cap = 0.5, math.radians(170.0)
    x0, x1, y = 4.99, 5.01, 5.0

    # midpoint (5,5) is the anchored split witness: two symmetric runs, largest TRUE run < cap
    total, max_run, _ = _stock_2.engagement_at(stock, 5.0, y, r, cap_ratio(math.pi))
    assert total == pytest.approx(2.0 * max_run, rel=0.05)
    assert max_run < cap  # no TRUE violation anywhere -- the motion is genuinely merge-free

    # the certifier's guarded cap lands strictly between the single-run and merged spans, so
    # gamma = 0 clears the station while gamma > 0 (gap absorbed) does not
    hs = 0.5 * (x1 - x0)
    cap_guarded = cap - 2.0 * _growth_bound(hs, r)
    assert max_run < cap_guarded < math.pi

    # discrimination probe (no production edit): at the anchored midpoint the exact station
    # verdict flips False -> True between gap_close_ratio 0 and the certifier's guard surrogate
    guard_ratio = cap_ratio(min(2.0 * _growth_bound(hs, r), math.pi))
    _, _, exceeded_unwired = _stock_2.engagement_at(stock, 5.0, y, r, cap_ratio(cap_guarded), 0.0)
    _, _, exceeded_wired = _stock_2.engagement_at(stock, 5.0, y, r, cap_ratio(cap_guarded), guard_ratio)
    assert not exceeded_unwired and exceeded_wired  # only the pessimism makes the station exceed

    # the guard: with the gamma-wiring the anchored merge is seen at every level -> refused
    _, certified, _ = _stock_2.certify_segment_tea(stock, x0, y, x1, y, r, cap)
    assert certified is False
