"""In-process engagement: stations inside a motion against stock minus its own swept prefix."""

from __future__ import annotations

import math

import numpy as np
import pytest

from compas_cgal import _stock_2

FULL_TURN_CAP_RATIO = 4.0
NO_GAP_CLOSURE_RATIO = 0.0
# radians; accumulated atan2 reporting only, never a decision.
REPORTING_ANGLE_ABS_TOLERANCE = 2e-12
# Exact cap surrogate (boundary doctrine): 4*sin^2(cap/2) for a 179 degree cap, injected once.
CAP_179_DEG_RATIO = 4.0 * math.sin(math.radians(179.0) / 2.0) ** 2
# The slow arc sweep under-covers by at most 1e-4 tool radii per its construction
# contract; at this bounded 2R guide case 0.1 degree isolates that residual from a
# wrong leading/trailing half.
SLOW_ARC_ORACLE_REPORTING_TOLERANCE = math.radians(0.1)
SLOW_ARC_ORACLE_APPROACH_ANGLE = 1e-6  # radians; resolves the tangent direction in binary64 (oracle only).


def _rectangle_stock() -> _stock_2.Stock2:
    return _stock_2.Stock2(np.array([[-8.0, -8.0], [8.0, -8.0], [8.0, 8.0], [-8.0, 8.0]], dtype=float), [])


def test_line_prefix_zero_is_the_static_station_and_interior_is_the_leading_half() -> None:
    stock = _rectangle_stock()
    before = stock.clone()
    zero = _stock_2.engagement_at(stock, -2.0, 0.0, 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)

    samples = _stock_2.engagement_line_prefixes(stock, -2.0, 0.0, 2.0, 0.0, [0.0, 0.25, 0.5, 1.0], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)

    assert samples[0] == zero
    for total, run, _exceeded in samples[1:]:
        assert math.isclose(total, math.pi, abs_tol=REPORTING_ANGLE_ABS_TOLERANCE)
        assert math.isclose(run, math.pi, abs_tol=REPORTING_ANGLE_ABS_TOLERANCE)
    assert stock.exactly_equals(before)  # the query never depletes; the replay owns depletion


def test_line_prefix_slot_cap_verdict_is_exact_at_the_half_turn() -> None:
    stock = _rectangle_stock()

    within = _stock_2.engagement_line_prefixes(stock, -2.0, 0.0, 2.0, 0.0, [0.5], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)
    over = _stock_2.engagement_line_prefixes(stock, -2.0, 0.0, 2.0, 0.0, [0.5], 1.0, CAP_179_DEG_RATIO, NO_GAP_CLOSURE_RATIO)

    assert within[0][2] is False
    assert over[0][2] is True


@pytest.mark.parametrize("clockwise", [False, True])
@pytest.mark.parametrize("guide_radius", [0.5, 1.0, 2.0])
def test_full_turn_prefix_starts_full_and_ends_empty(clockwise: bool, guide_radius: float) -> None:
    stock = _rectangle_stock()
    before = stock.clone()

    samples = _stock_2.engagement_arc_prefixes(
        stock,
        0.0,
        0.0,
        guide_radius,
        0.0,
        guide_radius,
        0.0,
        clockwise,
        [0.0, 0.25, 0.5, 0.75, 0.999, 1.0],
        1.0,
        FULL_TURN_CAP_RATIO,
        NO_GAP_CLOSURE_RATIO,
    )

    assert math.isclose(samples[0][1], math.tau, abs_tol=REPORTING_ANGLE_ABS_TOLERANCE)
    assert all(0.0 <= run <= math.pi + REPORTING_ANGLE_ABS_TOLERANCE for _total, run, _exceeded in samples[1:5])
    assert math.isclose(samples[-1][1], 0.0, abs_tol=REPORTING_ANGLE_ABS_TOLERANCE)
    assert stock.exactly_equals(before)


def test_arc_prefix_beyond_a_half_turn_uses_the_union_sector() -> None:
    # p = 0.75 of a CCW full turn sweeps 3*pi/2, so the sector test is the union.
    # The station at 270 degrees on a 2R guide in solid stock engages exactly its
    # leading half: the start disk is disjoint from it.
    stock = _rectangle_stock()

    [(_total, run, _exceeded)] = _stock_2.engagement_arc_prefixes(stock, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0, False, [0.75], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)

    assert math.isclose(run, math.pi, abs_tol=REPORTING_ANGLE_ABS_TOLERANCE)


def test_partial_arc_prefix_matches_the_slow_swept_stock_oracle() -> None:
    stock = _rectangle_stock()
    slow_prefix = stock.clone()
    guide_radius = 2.0
    half_parameter = 0.5

    samples = _stock_2.engagement_arc_prefixes(
        stock, 0.0, 0.0, guide_radius, 0.0, 0.0, guide_radius, False, [0.0, half_parameter, 1.0], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO
    )

    prefix_angle = half_parameter * math.pi / 2.0
    prior_angle = prefix_angle - SLOW_ARC_ORACLE_APPROACH_ANGLE
    slow_prefix.subtract_arc_sweep(0.0, 0.0, guide_radius, 0.0, guide_radius * math.cos(prior_angle), guide_radius * math.sin(prior_angle), False, 1.0)
    slow = _stock_2.engagement_at(slow_prefix, guide_radius * math.cos(prefix_angle), guide_radius * math.sin(prefix_angle), 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)
    assert len(samples) == 3
    assert math.isclose(samples[1][1], slow[1], abs_tol=SLOW_ARC_ORACLE_REPORTING_TOLERANCE)


def test_prefix_parameters_must_be_ordered_and_in_unit_interval() -> None:
    stock = _rectangle_stock()

    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_line_prefixes(stock, 0.0, 0.0, 1.0, 0.0, [0.5, 0.25], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_line_prefixes(stock, 0.0, 0.0, 1.0, 0.0, [1.5], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_line_prefixes(stock, 0.0, 0.0, 0.0, 0.0, [0.5], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_arc_prefixes(stock, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, [0.5], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)


NEAR_SEAM_PARAMETER = 1.0 - 2.0**-40
# radians, reporting: at a nearly full turn only a sliver wedge of rim is unswept.
NEAR_SEAM_SLIVER = 1e-3


@pytest.mark.parametrize("clockwise", [False, True])
def test_full_turn_station_just_before_the_seam_is_beyond_a_half_turn(clockwise: bool) -> None:
    # The prefix extent (at most / beyond a half turn) is decided by an exact
    # orientation on the exact chart station, never by a double swept angle, so a
    # station a hair before the seam reads like a nearly complete turn.
    stock = _rectangle_stock()

    [(_total, run, _exceeded)] = _stock_2.engagement_arc_prefixes(
        stock, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0, clockwise, [NEAR_SEAM_PARAMETER], 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO
    )

    assert run <= NEAR_SEAM_SLIVER


PREDECESSOR_GUIDE_RADIUS = 2.0
CANDIDATE_GUIDE_RADIUS = 2.0
TOOL_RADIUS = 1.0
HELD_STOCK_HALF_WIDTH = 20.0
CROSS_CHECK_STATIONS = 64
# radians; reporting comparison of two independently rounded acos/atan2 reports.
BOUND_REPORTING_TOLERANCE = 1e-9


def _held_stock() -> _stock_2.Stock2:
    w = HELD_STOCK_HALF_WIDTH
    stock = _stock_2.Stock2(np.array([[-w, -w], [w, -w], [w, w], [-w, w]], dtype=float), [])
    stock.subtract_disk(0.0, 0.0, PREDECESSOR_GUIDE_RADIUS + TOOL_RADIUS)
    return stock


@pytest.mark.parametrize("clockwise", [False, True])
@pytest.mark.parametrize(("candidate_x", "cap_deg", "bound_exceeds"), [(1.0, 110.0, False), (1.5, 140.0, False), (1.5, 120.0, True)])
def test_in_process_loop_engagement_respects_the_held_whole_orbit_bound(clockwise: bool, candidate_x: float, cap_deg: float, bound_exceeds: bool) -> None:
    contour = _stock_2.HeldDiskContour2((0.0, 0.0), PREDECESSOR_GUIDE_RADIUS, TOOL_RADIUS)
    ratio = 4.0 * math.sin(math.radians(cap_deg) / 2.0) ** 2
    bound, exceeded = contour.engagement_bound((candidate_x, 0.0), CANDIDATE_GUIDE_RADIUS, ratio)
    assert exceeded is bound_exceeds

    start = (candidate_x - CANDIDATE_GUIDE_RADIUS, 0.0)  # inside the cleared disk
    parameters = [i / CROSS_CHECK_STATIONS for i in range(1, CROSS_CHECK_STATIONS)]
    samples = _stock_2.engagement_arc_prefixes(_held_stock(), candidate_x, 0.0, *start, *start, clockwise, parameters, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO)

    assert max(run for _total, run, _exceeded in samples) <= bound + BOUND_REPORTING_TOLERANCE
    if bound_exceeds:
        assert any(e for _total, _run, e in samples)  # the bound is tight: the query sees the excess too
    else:
        assert not any(e for _total, _run, e in samples)  # exact: bound within cap => every station within cap


ASYMMETRIC_PREFIX_STATIONS = 16


def _asymmetric_stock() -> _stock_2.Stock2:
    """Solid stock with two cleared disks off every symmetry line of the guide circles below.

    A clockwise and a counter-clockwise prefix through this stock meet different
    material, so a reversed sweep cannot pass for the right one.
    """
    stock = _rectangle_stock()
    stock.subtract_disk(-1.0, 2.6, 1.1)
    stock.subtract_disk(2.5, -1.5, 0.7)
    return stock


def _slow_oracle_run(guide_radius: float, signed_angle: float, clockwise: bool) -> float:
    """Sweep the stock to just before the station at *signed_angle*, then read it statically."""
    direction = -1.0 if clockwise else 1.0
    prior = signed_angle - direction * SLOW_ARC_ORACLE_APPROACH_ANGLE
    slow = _asymmetric_stock()
    slow.subtract_arc_sweep(0.0, 0.0, guide_radius, 0.0, guide_radius * math.cos(prior), guide_radius * math.sin(prior), clockwise, 1.0)
    return _stock_2.engagement_at(slow, guide_radius * math.cos(signed_angle), guide_radius * math.sin(signed_angle), 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO)[1]


@pytest.mark.parametrize("clockwise", [False, True])
@pytest.mark.parametrize("guide_radius", [0.5, 1.0, 2.0])
def test_full_turn_prefix_follows_the_travel_direction_on_asymmetric_stock(clockwise: bool, guide_radius: float) -> None:
    parameters = [i / ASYMMETRIC_PREFIX_STATIONS for i in range(1, ASYMMETRIC_PREFIX_STATIONS)]
    samples = _stock_2.engagement_arc_prefixes(
        _asymmetric_stock(), 0.0, 0.0, guide_radius, 0.0, guide_radius, 0.0, clockwise, parameters, 1.0, FULL_TURN_CAP_RATIO, NO_GAP_CLOSURE_RATIO
    )

    direction = -1.0 if clockwise else 1.0
    for parameter, (_total, run, _exceeded) in zip(parameters, samples, strict=True):
        expected = _slow_oracle_run(guide_radius, direction * math.tau * parameter, clockwise)
        assert math.isclose(run, expected, abs_tol=SLOW_ARC_ORACLE_REPORTING_TOLERANCE), f"p={parameter}: {math.degrees(run)} vs oracle {math.degrees(expected)}"


@pytest.mark.parametrize("clockwise", [False, True])
def test_open_arc_prefix_follows_the_travel_direction_on_asymmetric_stock(clockwise: bool) -> None:
    # Three quarters of a turn in the travel direction from (2, 0); the same end
    # point is a quarter turn the other way, so a reversed sweep is detectable.
    guide_radius = 2.0
    direction = -1.0 if clockwise else 1.0
    sweep = 1.5 * math.pi
    end = direction * sweep
    parameters = [i / ASYMMETRIC_PREFIX_STATIONS for i in range(1, ASYMMETRIC_PREFIX_STATIONS)]
    samples = _stock_2.engagement_arc_prefixes(
        _asymmetric_stock(),
        0.0,
        0.0,
        guide_radius,
        0.0,
        guide_radius * math.cos(end),
        guide_radius * math.sin(end),
        clockwise,
        parameters,
        1.0,
        FULL_TURN_CAP_RATIO,
        NO_GAP_CLOSURE_RATIO,
    )

    for parameter, (_total, run, _exceeded) in zip(parameters, samples, strict=True):
        expected = _slow_oracle_run(guide_radius, direction * sweep * parameter, clockwise)
        assert math.isclose(run, expected, abs_tol=SLOW_ARC_ORACLE_REPORTING_TOLERANCE), f"p={parameter}: {math.degrees(run)} vs oracle {math.degrees(expected)}"
