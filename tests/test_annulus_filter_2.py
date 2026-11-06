"""The interval filter of the rim engine against the exact engine it must never contradict (unit U3 of the actual-stock refinement design).

The filter's contract is one-sided. A certified station is within the cap on the exact
engine, on the view that certified it, and the filter's chord bound is at or above the
exact engine's largest run there. A station it does not certify goes to the exact path,
so the decision always equals the exact decision. Integer fixtures reach the measure-zero
contacts; generic doubles are the witness that the fallback does not stall.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _stock_2

CAP_RATIO = _stock_2.cap_chord_ratio(math.radians(80.0))
LOOP_PARAMETERS = [k / 45.0 for k in range(45)]
STAR_DIRECTIONS = ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))
UNCERTIFIED = 3
# Reported spans are to_double of exact endpoints (chord squared units on the tool radius).
REPORTED_CHORD_AGREEMENT = 1e-9
# A generic-double decision that stalls is the failure this bound exists for (seconds, wall, under any load).
GENERIC_DECISION_WALL_S = 30.0
# How far below a chord enclosure's upper bound the cap walk starts (ulps of the cap); the enclosure is narrower.
CAP_WALK_ULPS = 256
Annulus = tuple[float, float, float]


def rows(points: list[tuple[float, float]]) -> np.ndarray:
    return np.array([(x, y, 0.0) for x, y in points], dtype=np.float64)


def square(half: float) -> np.ndarray:
    return rows([(-half, -half), (half, -half), (half, half), (-half, half)])


def open_view(polygon: np.ndarray, annuli: list[Annulus], tool: float) -> _stock_2.AnnulusView2:
    view = _stock_2.AnnulusView2.open(polygon, tool)
    for row in annuli:
        view.append(*row)
    return view


def tier_range(tier: int, end: int) -> tuple[int, int]:
    return {0: (max(0, end - 1), end), 1: (max(0, end - 16), end)}.get(tier, (0, end))


def assert_filter_sound(
    polygon: np.ndarray, annuli: list[Annulus], tool: float, circle: tuple[float, float, float, float, float], cap: float = CAP_RATIO
) -> list[tuple[int, float]]:
    """Every certified station is within on its view with its chord bound above the exact run; the decision equals the exact one."""
    view = open_view(polygon, annuli, tool)
    end = len(annuli)
    report = view.circle_filter_report(0, end, *circle, LOOP_PARAMETERS, cap)
    assert len(report) == len(LOOP_PARAMETERS)
    args = (*circle, False, LOOP_PARAMETERS, cap, 0.0)
    for tier in {row[0] for row in report}:
        first, last = tier_range(tier, end)
        exact = view.circle_prefixes(first, last, *args)
        for (station_tier, bound), (_total, span, exceeded) in zip(report, exact):
            if station_tier != tier or tier == UNCERTIFIED:
                continue
            assert not exceeded, "a certified station is within on the view that certified it"
            chord = 4.0 * math.sin(span / 2.0) ** 2 * tool * tool
            assert chord <= bound + REPORTED_CHORD_AGREEMENT, "M+ contains every exact run"
    assert view.circle_exceeded_filtered(0, end, *args) == view.circle_exceeded(0, end, *args)
    return report


# ---- the two legs of the cap test and the own prefix ---------------------------------


def test_a_run_past_a_half_turn_is_never_certified_by_its_small_chord() -> None:
    # A disk 1.9 away with radius 1.05 covers about 36 deg of the rim. The run left is about 324 deg: its chord squared is about 0.4, below the cap ratio 1.65, and it is exceeded.
    ex, ey = 1.0, 0.0
    annuli = [(ex + 1.9, ey, 0.05)]
    view = open_view(square(20.0), annuli, 1.0)
    args = (0.0, 0.0, 1.0, ex, ey)
    report = assert_filter_sound(square(20.0), annuli, 1.0, args)
    exact = view.circle_prefixes(0, 1, *args, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    assert exact[0][2], "station 0 keeps more than half the rim"
    assert 4.0 * math.sin(exact[0][1] / 2.0) ** 2 < CAP_RATIO, "the chord alone would pass"
    assert report[0][0] == UNCERTIFIED


def test_a_rim_with_no_cover_at_all_is_a_closed_loop_and_never_certified() -> None:
    report = assert_filter_sound(square(20.0), [], 1.0, (0.0, 0.0, 2.0, 2.0, 0.0))
    assert report[0][0] == UNCERTIFIED


def test_a_rim_inside_an_annulus_certifies_with_no_material() -> None:
    report = assert_filter_sound(square(20.0), [(0.0, 0.0, 4.0)], 1.0, (0.0, 0.0, 4.0, 4.0, 0.0))
    assert report[0] == (0, 0.0)


def test_the_own_prefix_certifies_stations_the_static_rim_exceeds() -> None:
    # Stations k > 0 of a concentric successor on a wide-open part: the previous annulus leaves the outer half of each rim.
    annuli = [(0.0, 0.0, 3.0)]
    circle = (0.0, 0.0, 3.5, 3.5, 0.0)
    report = assert_filter_sound(square(30.0), annuli, 1.0, circle)
    view = open_view(square(30.0), annuli, 1.0)
    static = view.station(0, 1, 3.5, 0.0, CAP_RATIO, 0.0)
    assert static[2], "the static rim at the entry exceeds"
    exact = view.circle_prefixes(0, 1, *circle, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    certified = [k for k, row in enumerate(report) if row[0] != UNCERTIFIED and k > 0]
    assert certified, "the clip has to do the certifying here"
    assert not any(exact[k][2] for k in certified)


def test_a_guide_smaller_than_the_tool_keeps_the_prefix_sound() -> None:
    assert_filter_sound(square(20.0), [(0.25, 0.0, 0.5)], 1.0, (0.0, 0.0, 0.5, 0.5, 0.0))


# ---- station 0 at a wall (ruling M3: a sub-ulp crossing is real stock) -----------------


@pytest.mark.parametrize(
    ("entry_x", "exact_exceeded", "tier"),
    [
        pytest.param(math.nextafter(9.0, 10.0), False, 0, id="crosses_by_one_ulp_which_this_wall_certifies"),
        pytest.param(9.0, True, UNCERTIFIED, id="tangent_keeps_one_run"),
        pytest.param(math.nextafter(9.0, 0.0), True, UNCERTIFIED, id="clear_of_the_wall_by_one_ulp_keeps_one_run"),
    ],
)
def test_a_sub_ulp_wall_on_exact_coordinates_keeps_the_exact_verdict(entry_x: float, exact_exceeded: bool, tier: int) -> None:
    # The merged run is 100 deg, over the 80 deg cap, and 50 deg per side when the wall splits it. On an axis-aligned
    # wall the enclosure is tight enough to see even a one-ulp crossing.
    annuli = [(entry_x - 1.0, 0.0, 0.8125)]
    view = open_view(square(10.0), annuli, 1.0)
    args = (entry_x - 2.0, 0.0, 2.0, entry_x, 0.0, False, [0.0], CAP_RATIO, 0.0)
    assert view.circle_filter_report(0, 1, *args[:5], [0.0], CAP_RATIO)[0][0] == tier
    assert view.circle_exceeded_filtered(0, 1, *args) is exact_exceeded


def test_a_crossing_inside_the_enclosure_falls_through_to_the_exact_verdict() -> None:
    # A sloped wall on large coordinates: the crossing is real and the enclosure straddles it, so the merged 100 deg run
    # stands and the filter cannot certify. The exact path sees the split (ruling M3: no tolerance, no near-tangent class).
    wall_a = (1004.0441298821165, -2000.7879140748912)
    wall_b = (1035.857917041054, 2000.786933132295)
    y0 = 1.8837966023380506
    entry_x = 1018.9659725372238
    polygon = rows([(-1000.1, -2000.3), wall_a, wall_b, (-1000.2, 2000.9)])
    view = open_view(polygon, [(entry_x - 1.0, y0, 0.8125)], 1.0)
    args = (entry_x - 2.0, y0, 2.0, entry_x, y0, False, [0.0], CAP_RATIO, 0.0)
    assert view.circle_filter_report(0, 1, *args[:5], [0.0], CAP_RATIO)[0][0] == UNCERTIFIED
    assert view.circle_exceeded_filtered(0, 1, *args) is False
    assert view.circle_exceeded(0, 1, *args) is False
    assert view.filter_counters.exact_stations == 1
    assert view.filter_counters.exact_station_zero == 1


def test_wall_tangency_whose_merged_run_is_within_needs_no_exact_arithmetic() -> None:
    # The same wall with a 60 deg merged run: the filter decides, whichever side of the wall the entry is on.
    entry_x = math.nextafter(9.0, 10.0)
    annuli = [(entry_x - 1.0, 0.0, 1.5)]
    view = open_view(square(10.0), annuli, 1.0)
    args = (entry_x - 2.0, 0.0, 2.0, entry_x, 0.0, False, [0.0], CAP_RATIO, 0.0)
    assert view.circle_exceeded_filtered(0, 1, *args) is False
    assert view.filter_counters.exact_stations == 0
    assert view.filter_counters.certified_last_one == 1


def test_wall_cover_certifies_a_rim_mostly_outside_the_part() -> None:
    # A rim centred on a straight wall keeps 180 deg, and with an annulus over most of that it is within only thanks to the wall.
    polygon = square(10.0)
    annuli = [(9.0, 4.0, 3.5)]
    report = assert_filter_sound(polygon, annuli, 1.0, (8.0, 0.0, 2.0, 10.0, 0.0))
    assert report[0][0] != UNCERTIFIED


def test_clockwise_polygon_gives_the_same_certificates() -> None:
    annuli = [(9.0, 4.0, 3.5)]
    circle = (8.0, 0.0, 2.0, 10.0, 0.0)
    ccw = open_view(square(10.0), annuli, 1.0).circle_filter_report(0, 1, *circle, LOOP_PARAMETERS, CAP_RATIO)
    cw = open_view(square(10.0)[::-1], annuli, 1.0).circle_filter_report(0, 1, *circle, LOOP_PARAMETERS, CAP_RATIO)
    assert [row[0] for row in ccw] == [row[0] for row in cw]


def test_rim_through_a_vertex_is_sound() -> None:
    polygon = rows([(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (0.0, 0.0), (-5.0, 5.0)])
    assert_filter_sound(polygon, [(0.0, -3.0, 1.0)], 1.0, (0.0, -3.0, 2.0, 0.0, -1.0))


def test_rim_outside_the_part_decides_through_the_exact_path() -> None:
    # Nothing is material, and the filter never proves a rim outside P without a wall cover.
    view = open_view(square(5.0), [], 1.0)
    args = (20.0, 0.0, 2.0, 22.0, 0.0, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    assert view.circle_exceeded_filtered(0, 0, *args) is False
    assert view.filter_counters.exact_stations > 0


# ---- nested views and fall-through ----------------------------------------------------


def test_a_station_escalates_only_as_far_as_it_needs() -> None:
    # The last annulus covers nothing near the rim; an earlier one covers it.
    annuli = [(0.0, 0.0, 4.0), (30.0, 30.0, 1.0)]
    view = open_view(square(40.0), annuli, 1.0)
    report = view.circle_filter_report(0, 2, 0.0, 0.0, 4.0, 4.0, 0.0, [0.0], CAP_RATIO)
    assert report[0][0] == 1
    deep = [(0.0, 0.0, 4.0)] + [(30.0, 30.0, 1.0)] * 16
    view = open_view(square(40.0), deep, 1.0)
    report = view.circle_filter_report(0, len(deep), 0.0, 0.0, 4.0, 4.0, 0.0, [0.0], CAP_RATIO)
    assert report[0][0] == 2


def test_forced_fall_through_decides_every_station_exactly_and_equally() -> None:
    annuli = [(0.0, 0.0, 3.0)]
    args = (0.0, 0.0, 3.5, 3.5, 0.0, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    view = open_view(square(30.0), annuli, 1.0)
    free = view.circle_exceeded_filtered(0, 1, *args)
    view.reset_filter_counters()
    view.force_filter_fallthrough(True)
    forced = view.circle_exceeded_filtered(0, 1, *args)
    assert forced == free == view.circle_exceeded(0, 1, *args)
    assert view.filter_counters.certified_last_one == 0
    assert view.filter_counters.exact_stations > 0


@pytest.mark.parametrize(
    "overrides",
    [
        pytest.param({"clockwise": True}, id="clockwise"),
        pytest.param({"gap": 0.5}, id="gap_closure"),
        pytest.param({"parameters": [0.0, 0.5, 1.0]}, id="closed_end_of_the_loop"),
    ],
)
def test_inputs_the_filter_does_not_take_go_wholly_to_the_exact_path(overrides: dict[str, object]) -> None:
    annuli = [(0.0, 0.0, 3.0)]
    view = open_view(square(30.0), annuli, 1.0)
    clockwise = overrides.get("clockwise", False)
    parameters = overrides.get("parameters", LOOP_PARAMETERS)
    gap = overrides.get("gap", 0.0)
    args = (0.0, 0.0, 3.5, 3.5, 0.0, clockwise, parameters, CAP_RATIO, gap)
    assert view.circle_exceeded_filtered(0, 1, *args) == view.circle_exceeded(0, 1, *args)  # type: ignore[arg-type]
    assert view.filter_counters.exact_decisions == 1
    assert view.filter_counters.stations == 0


@pytest.mark.parametrize(
    ("args", "error"),
    [
        pytest.param((0.0, 0.0, -1.0, 1.0, 0.0), ValueError, id="negative_guide_radius"),
        pytest.param((0.0, 0.0, 1.0, 0.0, 0.0), ValueError, id="entry_on_the_guide_centre"),
        pytest.param((math.inf, 0.0, 1.0, 1.0, 0.0), ValueError, id="non_finite_guide"),
    ],
)
def test_invalid_circles_raise_what_the_exact_path_raises(args: tuple[float, float, float, float, float], error: type[Exception]) -> None:
    view = open_view(square(5.0), [], 1.0)
    with pytest.raises(error) as filtered:
        view.circle_exceeded_filtered(0, 0, *args, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    with pytest.raises(error) as exact:
        view.circle_exceeded(0, 0, *args, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    assert type(filtered.value) is type(exact.value)


@pytest.mark.parametrize("parameters", [[0.5, 0.25], [-0.1], [float("nan")], [1.5]])
def test_invalid_parameters_raise_what_the_exact_path_raises(parameters: list[float]) -> None:
    view = open_view(square(5.0), [], 1.0)
    args = (0.0, 0.0, 1.0, 1.0, 0.0, False, parameters, CAP_RATIO, 0.0)
    with pytest.raises(ValueError) as filtered:
        view.circle_exceeded_filtered(0, 0, *args)
    with pytest.raises(ValueError) as exact:
        view.circle_exceeded(0, 0, *args)
    assert type(filtered.value) is type(exact.value)


def test_invalid_cap_ratio_raises_on_a_cached_parameter_list() -> None:
    view = open_view(square(5.0), [], 1.0)
    args = (0.0, 0.0, 1.0, 1.0, 0.0, False, LOOP_PARAMETERS)
    view.circle_exceeded_filtered(0, 0, *args, CAP_RATIO, 0.0)
    with pytest.raises(ValueError):
        view.circle_exceeded_filtered(0, 0, *args, 9.0, 0.0)


def test_truncation_reaches_the_filter() -> None:
    annuli = [(0.0, 0.0, 4.0), (0.0, 0.0, 4.0)]
    view = open_view(square(30.0), annuli, 1.0)
    view.truncate(0)
    args = (0.0, 0.0, 4.0, 4.0, 0.0, False, [0.0], CAP_RATIO, 0.0)
    assert view.circle_exceeded_filtered(0, 0, *args) is True


def test_generic_double_wall_station_decides_quickly_and_equals_exact() -> None:
    entry_x = 8.999999999999998
    annuli = [(entry_x - 1.0000000000000002, 1e-3, 0.8125000000000001), (3.3333333333333335, -2.718281828459045, 1.4142135623730951)]
    view = open_view(square(10.0), annuli, 1.0)
    args = (entry_x - 2.0000000000000004, 1e-3, 2.0000000000000004, entry_x, 0.0, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    start = time.perf_counter()
    assert view.circle_exceeded_filtered(0, 2, *args) == view.circle_exceeded(0, 2, *args)
    assert time.perf_counter() - start < GENERIC_DECISION_WALL_S


# ---- filtered == exact on random stock ----------------------------------------------

integer_annuli = st.lists(st.tuples(st.integers(-9, 9), st.integers(-9, 9), st.integers(0, 6)), max_size=8)
integer_radii = st.lists(st.integers(1, 12), min_size=8, max_size=8)
generic = st.floats(min_value=-9.0, max_value=9.0, allow_nan=False, allow_infinity=False, allow_subnormal=False)
generic_annuli = st.lists(st.tuples(generic, generic, st.floats(min_value=0.0, max_value=6.0, allow_subnormal=False)), max_size=24)
generic_radii = st.lists(st.floats(min_value=0.5, max_value=12.0, allow_subnormal=False), min_size=8, max_size=8)
phase_angle = st.floats(min_value=-3.0, max_value=3.0, allow_subnormal=False)
guide_radius = st.floats(min_value=0.25, max_value=6.0, allow_subnormal=False)
tool_radius = st.floats(0.25, 3.0, allow_subnormal=False)
stream_guides = st.lists(st.floats(min_value=0.1, max_value=6.0, allow_subnormal=False), min_size=1, max_size=20)
spacing_ratio = st.floats(min_value=0.1, max_value=0.9, allow_subnormal=False)
axis_direction = st.sampled_from([(1, 0), (0, 1), (-1, 0), (0, -1)])


def star(radii: list[float]) -> np.ndarray:
    return rows([(dx * radius, dy * radius) for (dx, dy), radius in zip(STAR_DIRECTIONS, radii)])


@settings(max_examples=120, deadline=None)
@given(integer_radii, integer_annuli, st.integers(-9, 9), st.integers(-9, 9), st.integers(1, 6), axis_direction)
def test_h1_integer_fixtures_filter_is_sound(radii: list[int], annuli: list[tuple[int, int, int]], x: int, y: int, rho: int, direction: tuple[int, int]) -> None:
    rows_ = [(float(a), float(b), float(c)) for a, b, c in annuli]
    circle = (float(x), float(y), float(rho), float(x + rho * direction[0]), float(y + rho * direction[1]))
    assert_filter_sound(star([float(r) for r in radii]), rows_, 1.0, circle)


@settings(max_examples=80, deadline=None)
@given(
    generic_radii,
    generic_annuli,
    generic,
    generic,
    st.floats(min_value=0.25, max_value=6.0, allow_subnormal=False),
    st.floats(min_value=-3.0, max_value=3.0, allow_subnormal=False),
    st.floats(0.25, 3.0, allow_subnormal=False),
)
def test_h1_generic_doubles_filter_is_sound(radii: list[float], annuli: list[Annulus], x: float, y: float, rho: float, phase: float, tool: float) -> None:
    circle = (x, y, rho, x + rho * math.cos(phase), y + rho * math.sin(phase))
    assert_filter_sound(star(radii), annuli, tool, circle)


@settings(max_examples=60, deadline=None)
@given(
    generic_radii,
    st.lists(st.floats(min_value=0.1, max_value=6.0, allow_subnormal=False), min_size=1, max_size=20),
    st.floats(min_value=-3.0, max_value=3.0, allow_subnormal=False),
    st.floats(min_value=0.1, max_value=0.9, allow_subnormal=False),
)
def test_h1_a_held_like_stream_of_touching_circles_is_sound(radii: list[float], guides: list[float], phase: float, spacing: float) -> None:
    """Successive circles sharing a contact direction: the nearly-coincident boundaries the workloads produce."""
    annuli: list[Annulus] = []
    position = 0.0
    for guide in guides:
        annuli.append((position * math.cos(phase), position * math.sin(phase), guide))
        position += spacing * guide
    guide = guides[-1] + spacing
    cx, cy = position * math.cos(phase), position * math.sin(phase)
    assert_filter_sound(star(radii), annuli, 1.0, (cx, cy, guide, cx + guide * math.cos(phase), cy + guide * math.sin(phase)))


# ---- thresholds inside the enclosures: left to the exact path, never decided on a midpoint ---------


def _station_is_sound(view: _stock_2.AnnulusView2, entry: tuple[float, float], guide: tuple[float, float, float], end: int, cap: float) -> None:
    circle = (*guide, *entry)
    report = view.circle_filter_report(0, end, *circle, [0.0], cap)
    exact = view.circle_prefixes(0, end, *circle, False, [0.0], cap, 0.0)
    if report[0][0] != UNCERTIFIED:
        assert not exact[0][2], "certified inside the enclosure of the deciding quantity"
    assert view.circle_exceeded_filtered(0, end, *circle, False, [0.0], cap, 0.0) == exact[0][2]


@pytest.mark.parametrize("angle", [0.1 + 0.15 * k for k in range(40)])
def test_caps_within_the_chord_enclosure_are_left_to_the_exact_path(angle: float) -> None:
    """The cap walks ulp by ulp through the enclosure of one run's chord squared."""
    entry = (3.3 + math.cos(angle), -2.7 + math.sin(angle))
    centre = (entry[0] + 1.0000000000000002 * math.cos(angle), entry[1] + 1.0000000000000002 * math.sin(angle))
    view = open_view(square(40.0), [(*centre, 0.8125000000000001)], 1.0)
    guide = (entry[0] - 2.0 * math.cos(angle), entry[1] - 2.0 * math.sin(angle), 2.0)
    bound = view.circle_filter_report(0, 1, *guide, *entry, [0.0], 3.9)[0][1]
    assert bound > 0.0
    cap = bound
    for _ in range(CAP_WALK_ULPS):
        cap = math.nextafter(cap, 0.0)
    for _ in range(CAP_WALK_ULPS + CAP_WALK_ULPS // 8):
        _station_is_sound(view, entry, guide, 1, cap)
        cap = math.nextafter(cap, 4.0)


@pytest.mark.parametrize("cap", [4.0, CAP_RATIO])
@pytest.mark.parametrize("angle", [0.2, 0.9, 1.6, 2.3, 3.0, 3.7, 4.4, 5.1, 5.8])
def test_a_run_of_about_half_a_turn_is_left_to_the_exact_path(angle: float, cap: float) -> None:
    """Orientation inside its enclosure: on the cap of a half turn a run a few ulps over is exceeded with a chord below the cap."""
    entry = (3.3 + math.cos(angle), -2.7 + math.sin(angle))
    distance = 1.5
    centre = (entry[0] + distance * math.cos(angle), entry[1] + distance * math.sin(angle))
    guide_radius = math.sqrt(distance * distance + 1.0) - 1.0
    for _ in range(24):
        guide_radius = math.nextafter(guide_radius, 0.0)
    guide = (entry[0] - 2.0 * math.cos(angle), entry[1] - 2.0 * math.sin(angle), 2.0)
    for _ in range(48):
        view = open_view(square(40.0), [(*centre, guide_radius)], 1.0)
        _station_is_sound(view, entry, guide, 1, cap)
        guide_radius = math.nextafter(guide_radius, 4.0)


def test_a_tiny_cover_whose_enclosed_ends_cross_is_not_read_as_nearly_the_whole_rim() -> None:
    """A rim crossing a disk by ulps just off an x extreme: the arc's two ends share one key enclosure, and the rest of the rim is a closed loop."""
    for angle in (1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, math.pi - 1e-8, math.pi - 1e-7, 1e-8 - math.pi):
        for centre_distance in (3.0, 7.3):
            entry = (3.3, -2.7)
            centre = (entry[0] + centre_distance * math.cos(angle), entry[1] + centre_distance * math.sin(angle))
            outer_needed = math.hypot(centre[0] - entry[0], centre[1] - entry[1]) - 1.0
            guide_radius = outer_needed - 1.0
            for _ in range(8):
                guide_radius = math.nextafter(guide_radius, 0.0)
            guide = (entry[0] - 2.0, entry[1], 2.0)
            for _ in range(16):
                view = open_view(square(40.0), [(*centre, guide_radius)], 1.0)
                _station_is_sound(view, entry, guide, 1, CAP_RATIO)
                guide_radius = math.nextafter(guide_radius, 40.0)


def test_forced_fall_through_applies_the_own_prefix_clip() -> None:
    """Stations whose static rim exceeds and whose clipped rim is within: an exact fallback that skips the clip would refuse."""
    annuli = [(-0.4, 0.0, 4.0)]
    circle = (0.0, 0.0, 4.2, -4.2, 0.0)
    view = open_view(square(30.0), annuli, 1.0)
    static = [view.station(0, 1, 4.2 * math.cos(2.0 * math.pi * k / 45.0), 4.2 * math.sin(2.0 * math.pi * k / 45.0), CAP_RATIO, 0.0)[2] for k in range(1, 45)]
    assert any(static)
    view.force_filter_fallthrough(True)
    assert view.circle_exceeded_filtered(0, 1, *circle, False, LOOP_PARAMETERS, CAP_RATIO, 0.0) is False
    assert view.circle_exceeded(0, 1, *circle, False, LOOP_PARAMETERS, CAP_RATIO, 0.0) is False


def test_a_new_parameter_list_gets_its_own_rotations() -> None:
    annuli = [(0.0, 0.0, 3.0)]
    circle = (0.0, 0.0, 3.5, 3.5, 0.0)
    view = open_view(square(30.0), annuli, 1.0)
    args = (*circle, False)
    view.circle_exceeded_filtered(0, 1, *args, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    assert len(view.circle_filter_report(0, 1, *circle, [0.0, 0.5], CAP_RATIO)) == 2
    fresh = open_view(square(30.0), annuli, 1.0)
    assert view.circle_filter_report(0, 1, *circle, [0.0, 0.5], CAP_RATIO) == fresh.circle_filter_report(0, 1, *circle, [0.0, 0.5], CAP_RATIO)
