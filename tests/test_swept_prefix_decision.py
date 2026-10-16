"""Decision-only swept-prefix queries: the verdict is exactly ``any`` over the prefix call's ``cap_exceeded`` flags.

`engagement_{line,arc,circle}_exceeded` classify the same stations with the
same exact predicates as `engagement_{line,arc,circle}_prefixes`; only the
probing order (last station first) and the early exit differ, and "some
station exceeds" does not depend on order.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _stock_2

TOOL_RADIUS = 1.0
NO_GAP_CLOSURE_RATIO = 0.0
# The gate's 45 samples per motion: k / 45 for k = 0..45 (open) and k < 45 (loops).
SAMPLES = 45
MOTION_PARAMETERS = [k / SAMPLES for k in range(SAMPLES + 1)]
LOOP_PARAMETERS = [k / SAMPLES for k in range(SAMPLES)]
CAPS_DEG = (40.0, 90.0, 120.0, 179.0)


def _ratio(cap_deg: float) -> float:
    return _stock_2.cap_chord_ratio(math.radians(cap_deg))


def _rectangle() -> _stock_2.Stock2:
    return _stock_2.Stock2(np.array([[-6.0, -4.0], [6.0, -4.0], [6.0, 4.0], [-6.0, 4.0]], dtype=float), [])


def _partly_cleared() -> _stock_2.Stock2:
    """A rectangle with a cleared disk and a cleared slot: engagement varies along every fixture motion."""
    stock = _rectangle()
    stock.subtract_disk(-2.0, 0.0, 1.0)
    stock.subtract_capsule(-2.0, 0.0, 1.5, 1.2, 1.0)
    return stock


def _generic() -> _stock_2.Stock2:
    """Generic doubles (no dyadic structure) in the stock, the cut and the motions: the witness against identities."""
    stock = _stock_2.Stock2(np.array([[-5.913, -3.871], [6.127, -4.093], [5.771, 4.219], [-6.233, 3.917]], dtype=float), [])
    stock.subtract_disk(-1.7310000000000001, 0.3170000000000002, 1.0730000000000002)
    return stock


STOCKS = {"virgin": _rectangle, "partly_cleared": _partly_cleared, "generic": _generic}

LINES = [(-2.0, 0.0, 2.0, 0.0), (-2.0, 0.0, -2.0, 2.5), (-1.7310000000000001, 0.3170000000000002, 2.0930000000000004, -1.1130000000000002), (1.5, 1.2, -3.3, 1.9)]
ARCS = [
    # (guide_cx, guide_cy, start_x, start_y, end_x, end_y, clockwise)
    (-2.0, 0.0, -1.0, 0.0, -2.0, 1.0, False),
    (-2.0, 0.0, -1.0, 0.0, -2.0, 1.0, True),
    (-1.7310000000000001, 0.3170000000000002, -0.6410000000000002, 0.5230000000000001, -2.9130000000000003, -0.4110000000000001, True),
    (0.5, 0.5, 2.0, 0.5, 2.0, 0.5, False),
]
CIRCLES = [
    # (guide_cx, guide_cy, guide_radius, entry_x, entry_y, clockwise)
    (-2.0, 0.0, 0.6, -1.4, 0.0, True),
    (-2.0, 0.0, 0.6, -1.4, 0.0, False),
    (0.3170000000000002, -0.7310000000000001, 0.9130000000000001, 1.2300000000000002, -0.7310000000000001, True),
    (1.5, 1.2, 1.1, 1.5, 2.3, False),
]


def _any_exceeded(samples: list[tuple[float, float, bool]]) -> bool:
    return any(exceeded for _total, _run, exceeded in samples)


@pytest.mark.parametrize("cap_deg", CAPS_DEG)
@pytest.mark.parametrize("stock_name", STOCKS)
@pytest.mark.parametrize("line", LINES)
def test_line_decision_is_any_of_the_prefix_flags(stock_name: str, line: tuple[float, float, float, float], cap_deg: float) -> None:
    stock = STOCKS[stock_name]()
    before = stock.clone()
    ratio = _ratio(cap_deg)
    expected = _any_exceeded(_stock_2.engagement_line_prefixes(stock, *line, MOTION_PARAMETERS, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO))
    assert _stock_2.engagement_line_exceeded(stock, *line, MOTION_PARAMETERS, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO) is expected
    assert stock.exactly_equals(before)


@pytest.mark.parametrize("cap_deg", CAPS_DEG)
@pytest.mark.parametrize("stock_name", STOCKS)
@pytest.mark.parametrize("arc", ARCS)
def test_arc_decision_is_any_of_the_prefix_flags(stock_name: str, arc: tuple[float, float, float, float, float, float, bool], cap_deg: float) -> None:
    stock = STOCKS[stock_name]()
    before = stock.clone()
    ratio = _ratio(cap_deg)
    expected = _any_exceeded(_stock_2.engagement_arc_prefixes(stock, *arc, MOTION_PARAMETERS, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO))
    assert _stock_2.engagement_arc_exceeded(stock, *arc, MOTION_PARAMETERS, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO) is expected
    assert stock.exactly_equals(before)


@pytest.mark.parametrize("cap_deg", CAPS_DEG)
@pytest.mark.parametrize("stock_name", STOCKS)
@pytest.mark.parametrize("circle", CIRCLES)
def test_circle_decision_is_any_of_the_prefix_flags(stock_name: str, circle: tuple[float, float, float, float, float, bool], cap_deg: float) -> None:
    stock = STOCKS[stock_name]()
    before = stock.clone()
    ratio = _ratio(cap_deg)
    expected = _any_exceeded(_stock_2.engagement_circle_prefixes(stock, *circle, LOOP_PARAMETERS, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO))
    assert _stock_2.engagement_circle_exceeded(stock, *circle, LOOP_PARAMETERS, TOOL_RADIUS, ratio, NO_GAP_CLOSURE_RATIO) is expected
    assert stock.exactly_equals(before)


def test_the_fixtures_reach_both_verdicts() -> None:
    """The agreement tests are not vacuous: every motion kind is refused somewhere and admitted somewhere."""
    seen: dict[str, set[bool]] = {"line": set(), "arc": set(), "circle": set()}
    for make in STOCKS.values():
        stock = make()
        for cap_deg in CAPS_DEG:
            ratio = _ratio(cap_deg)
            for line in LINES:
                seen["line"].add(_stock_2.engagement_line_exceeded(stock, *line, MOTION_PARAMETERS, TOOL_RADIUS, ratio))
            for arc in ARCS:
                seen["arc"].add(_stock_2.engagement_arc_exceeded(stock, *arc, MOTION_PARAMETERS, TOOL_RADIUS, ratio))
            for circle in CIRCLES:
                seen["circle"].add(_stock_2.engagement_circle_exceeded(stock, *circle, LOOP_PARAMETERS, TOOL_RADIUS, ratio))
    assert seen == {"line": {False, True}, "arc": {False, True}, "circle": {False, True}}


def test_interior_exceedances_are_found_when_both_ends_are_under_the_cap() -> None:
    """Both ends sit in cleared slots and only interior stations cross the web between them: the scan must reach past the end-first probe."""
    stock = _rectangle()
    stock.subtract_capsule(-4.0, 0.0, -2.0, 0.0, 1.25)
    stock.subtract_capsule(2.0, 0.0, 4.0, 0.0, 1.25)  # a solid web left across the slot, |x| < 0.75
    ratio = _ratio(40.0)
    parameters = [0.0, 0.25, 0.5, 0.75, 1.0]
    flags = [exceeded for _t, _r, exceeded in _stock_2.engagement_line_prefixes(stock, -3.0, 0.0, 3.0, 0.0, parameters, TOOL_RADIUS, ratio)]
    assert flags[0] is False and flags[-1] is False and any(flags)
    assert _stock_2.engagement_line_exceeded(stock, -3.0, 0.0, 3.0, 0.0, parameters, TOOL_RADIUS, ratio) is True


def test_no_station_is_no_exceedance() -> None:
    stock = _rectangle()
    assert _stock_2.engagement_line_exceeded(stock, -2.0, 0.0, 2.0, 0.0, [], TOOL_RADIUS, _ratio(40.0)) is False


@pytest.mark.parametrize("parameters", [[0.5, 0.25], [-0.1], [1.5], [math.nan]])
def test_the_decision_refuses_the_inputs_the_prefix_call_refuses(parameters: list[float]) -> None:
    stock = _rectangle()
    ratio = _ratio(120.0)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_line_prefixes(stock, -2.0, 0.0, 2.0, 0.0, parameters, TOOL_RADIUS, ratio)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_line_exceeded(stock, -2.0, 0.0, 2.0, 0.0, parameters, TOOL_RADIUS, ratio)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_arc_exceeded(stock, -2.0, 0.0, -1.0, 0.0, -2.0, 1.0, False, parameters, TOOL_RADIUS, ratio)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_circle_exceeded(stock, -2.0, 0.0, 0.6, -1.4, 0.0, True, parameters, TOOL_RADIUS, ratio)


def test_a_degenerate_motion_is_refused_like_the_prefix_call() -> None:
    stock = _rectangle()
    ratio = _ratio(120.0)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_line_exceeded(stock, 1.0, 1.0, 1.0, 1.0, MOTION_PARAMETERS, TOOL_RADIUS, ratio)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_arc_exceeded(stock, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, False, MOTION_PARAMETERS, TOOL_RADIUS, ratio)
    with pytest.raises(_stock_2.InvalidSweptPrefixInputError):
        _stock_2.engagement_circle_exceeded(stock, 1.0, 1.0, 0.5, 1.0, 1.0, True, LOOP_PARAMETERS, TOOL_RADIUS, ratio)


# ---- per-station agreement (Hypothesis): the decision-only station verdict skips the
# clip when the static verdict is WITHIN_CAP (clip monotonicity, swept_prefix_engagement_2.h);
# every single station must still read exactly what the full classification reads.

COORDINATE = st.floats(min_value=-4.5, max_value=4.5, allow_nan=False, allow_infinity=False)
CUT = st.tuples(COORDINATE, COORDINATE, COORDINATE, COORDINATE, st.floats(min_value=0.3, max_value=1.6))
PARAMETER = st.one_of(st.sampled_from([0.0, 1.0, 0.5]), st.floats(min_value=0.0, max_value=1.0, allow_nan=False))
CAP_RATIO = st.sampled_from([_ratio(cap) for cap in (20.0, 40.0, 90.0, 120.0, 179.0, 180.0)])
GAP_RATIO = st.sampled_from([0.0, _ratio(5.0), _ratio(30.0)])


def _cut_stock(cuts: list[tuple[float, float, float, float, float]]) -> _stock_2.Stock2:
    stock = _rectangle()
    for x0, y0, x1, y1, radius in cuts:
        if (x0, y0) == (x1, y1):
            stock.subtract_disk(x0, y0, radius)
        else:
            stock.subtract_capsule(x0, y0, x1, y1, radius)
    return stock


@settings(max_examples=150, deadline=None)
@given(st.lists(CUT, max_size=3), COORDINATE, COORDINATE, COORDINATE, COORDINATE, PARAMETER, CAP_RATIO, GAP_RATIO)
def test_every_line_station_reads_the_full_classification(cuts, x0, y0, x1, y1, parameter, ratio, gap) -> None:  # noqa: ANN001
    assume((x0, y0) != (x1, y1))
    stock = _cut_stock(cuts)
    expected = _stock_2.engagement_line_prefixes(stock, x0, y0, x1, y1, [parameter], TOOL_RADIUS, ratio, gap)[0][2]
    assert _stock_2.engagement_line_exceeded(stock, x0, y0, x1, y1, [parameter], TOOL_RADIUS, ratio, gap) is expected


@settings(max_examples=150, deadline=None)
@given(
    st.lists(CUT, max_size=3),
    COORDINATE,
    COORDINATE,
    st.floats(min_value=0.2, max_value=2.5),
    st.floats(min_value=0.0, max_value=6.2),
    st.floats(min_value=0.05, max_value=6.2),
    st.booleans(),
    PARAMETER,
    CAP_RATIO,
    GAP_RATIO,
)
def test_every_arc_station_reads_the_full_classification(cuts, cx, cy, radius, phase, sweep, clockwise, parameter, ratio, gap) -> None:  # noqa: ANN001
    stock = _cut_stock(cuts)
    sx, sy = cx + radius * math.cos(phase), cy + radius * math.sin(phase)
    turn = -sweep if clockwise else sweep
    ex, ey = cx + radius * math.cos(phase + turn), cy + radius * math.sin(phase + turn)
    assume((sx, sy) != (cx, cy) and (ex, ey) != (cx, cy))
    args = (stock, cx, cy, sx, sy, ex, ey, clockwise, [parameter], TOOL_RADIUS, ratio, gap)
    try:
        expected = _stock_2.engagement_arc_prefixes(*args)[0][2]
    except _stock_2.ArcPrefixExtentAmbiguousError:
        with pytest.raises(_stock_2.ArcPrefixExtentAmbiguousError):
            _stock_2.engagement_arc_exceeded(*args)
        return
    assert _stock_2.engagement_arc_exceeded(*args) is expected


@settings(max_examples=150, deadline=None)
@given(
    st.lists(CUT, max_size=3),
    COORDINATE,
    COORDINATE,
    st.floats(min_value=0.2, max_value=2.5),
    st.floats(min_value=0.0, max_value=6.2),
    st.booleans(),
    PARAMETER,
    CAP_RATIO,
    GAP_RATIO,
)
def test_every_circle_station_reads_the_full_classification(cuts, cx, cy, radius, phase, clockwise, parameter, ratio, gap) -> None:  # noqa: ANN001
    stock = _cut_stock(cuts)  # parameter 1.0 (the closing station) is drawn: PARAMETER samples it
    ex, ey = cx + radius * math.cos(phase), cy + radius * math.sin(phase)
    assume((ex, ey) != (cx, cy))
    args = (stock, cx, cy, radius, ex, ey, clockwise, [parameter], TOOL_RADIUS, ratio, gap)
    assert _stock_2.engagement_circle_exceeded(*args) is _stock_2.engagement_circle_prefixes(*args)[0][2]


# ---- the clip flip: static CAP_EXCEEDED, in-process WITHIN_CAP. The decision call must clip
# there; a shortcut that answers "exceeded" from the static verdict alone fails these.
# Generic doubles, deep in virgin stock: every station is fully immersed at rest (more than a
# tool radius from the boundary), so its static verdict is a full turn, exceeded at any cap.
FLIP_STOCK_OUTLINE = [[-5.913, -3.871], [6.127, -4.093], [5.771, 4.219], [-6.233, 3.917]]
FLIP_CENTRE = (0.3170000000000002, -0.7310000000000001)
FLIP_GUIDE_RADIUS = 0.2130000000000001
# A line's in-process rim in virgin stock is exactly the leading half turn: only the
# half-turn cap (chord ratio 4) admits it, so that is where a line can flip.
HALF_TURN_RATIO = 4.0


def _flip_stock() -> _stock_2.Stock2:
    return _stock_2.Stock2(np.array(FLIP_STOCK_OUTLINE, dtype=float), [])


def _fully_immersed_at_rest(stock: _stock_2.Stock2, x: float, y: float) -> bool:
    total, run, exceeded = _stock_2.engagement_at(stock, x, y, TOOL_RADIUS, _ratio(179.0))
    return exceeded and run == pytest.approx(math.tau)


def _flip_cases() -> list[tuple[str, tuple[object, ...], list[float], list[tuple[float, float]]]]:
    cx, cy = FLIP_CENTRE
    rho = FLIP_GUIDE_RADIUS
    x0, y0, x1, y1 = -1.7310000000000001, 0.3170000000000002, 1.3130000000000002, -0.4110000000000001
    line_parameters = [0.5, 1.0]
    circle_parameters = [0.5, 0.9, 1.0]  # 1.0: the closing station, the whole rim behind the full turn
    arc_parameters = [0.3, 0.7, 1.0]
    on_guide = [(cx + rho * math.cos(a), cy + rho * math.sin(a)) for a in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0)]
    return [
        ("line", (x0, y0, x1, y1), line_parameters, [(x0 + p * (x1 - x0), y0 + p * (y1 - y0)) for p in line_parameters]),
        ("circle", (cx, cy, rho, cx + rho, cy, True), circle_parameters, on_guide),
        ("arc", (cx, cy, cx + rho, cy, cx, cy + rho, False), arc_parameters, on_guide),
    ]


FLIP_RATIOS = {"line": HALF_TURN_RATIO, "circle": _ratio(120.0), "arc": _ratio(170.0)}


@pytest.mark.parametrize(("kind", "motion", "parameters", "stations"), _flip_cases(), ids=lambda case: case if isinstance(case, str) else "")
def test_the_decision_clips_where_the_static_verdict_flips(kind: str, motion: tuple[object, ...], parameters: list[float], stations: list[tuple[float, float]]) -> None:
    stock = _flip_stock()
    ratio = FLIP_RATIOS[kind]
    assert all(_fully_immersed_at_rest(stock, x, y) for x, y in stations)  # static: exceeded at every station
    full = getattr(_stock_2, f"engagement_{kind}_prefixes")(stock, *motion, parameters, TOOL_RADIUS, ratio)
    assert not _any_exceeded(full)  # in process: within at every station -- the flip
    for parameter in parameters:  # station by station, and as one motion
        assert getattr(_stock_2, f"engagement_{kind}_exceeded")(stock, *motion, [parameter], TOOL_RADIUS, ratio) is False
    assert getattr(_stock_2, f"engagement_{kind}_exceeded")(stock, *motion, parameters, TOOL_RADIUS, ratio) is _any_exceeded(full)
