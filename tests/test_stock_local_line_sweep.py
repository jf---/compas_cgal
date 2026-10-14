"""Local line and circle sweeps equal their global twins, exactly and canonically.

The local path edits the stock's own arrangement around the removed region and
floods it from a seed, stopping only at edges that lie on a SUB-CURVE of the
region's boundary (docs/local_depletion.md). The gate is the double one the page
requires: exact point-set equality with the global difference AND the
`General_polygon_set_2` representation invariant, after EVERY removal, plus
equal arrangement sizes.
"""

import time

from compas.geometry import Polygon
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.stock import Stock

TOOL_RADII = (0.5, 1.0)
# Dyadic lattice: quarter-millimetre coordinates are exact doubles, so collinear
# sub-segments built on it are EXACTLY collinear -- the shared-support case.
LATTICE = st.integers(-4, 52).map(lambda q: q / 4.0)
# Pythagorean directions (scaled by small dyadic lengths) and a generic oblique one.
PYTHAGOREAN_STEPS = ((3, 4), (4, -3), (-5, 12), (12, 5), (1, 0), (0, 1), (8, 15))
EXAMPLES = 60
# Wall-clock ceiling, seconds, for two generic-double oblique sweeps on a 12x8
# stock: the exact path takes milliseconds; a filter collapse into bignum
# evaluation of an identically-zero decision takes seconds to minutes.
GENERIC_DOUBLE_CEILING_S = 1.0


def _rect() -> Stock:
    return Stock(Polygon([[0, 0, 0], [12, 0, 0], [12, 8, 0], [0, 8, 0]]))


@st.composite
def _segments(draw: st.DrawFn) -> tuple[float, float, float, float]:
    kind = draw(st.sampled_from(("axis", "pythagorean", "oblique", "collinear", "point")))
    x0, y0 = draw(LATTICE), draw(LATTICE)
    if kind == "axis":
        if draw(st.booleans()):
            return (x0, y0, draw(LATTICE), y0)
        return (x0, y0, x0, draw(LATTICE))
    if kind == "pythagorean":
        dx, dy = draw(st.sampled_from(PYTHAGOREAN_STEPS))
        scale = draw(st.sampled_from((0.25, 0.5, 1.0, -0.75)))
        return (x0, y0, x0 + scale * dx, y0 + scale * dy)
    if kind == "oblique":
        fx = draw(st.floats(-1.0, 13.0, allow_nan=False))
        fy = draw(st.floats(-1.0, 9.0, allow_nan=False))
        return (x0, y0, fx, fy)
    if kind == "collinear":
        dx, dy = draw(st.integers(-6, 6)), draw(st.integers(-6, 6))
        t0, t1 = draw(st.sampled_from((0.0, 0.25, 0.5))), draw(st.sampled_from((0.75, 1.0)))
        return (x0 + t0 * dx, y0 + t0 * dy, x0 + t1 * dx, y0 + t1 * dy)
    return (x0, y0, x0, y0)


@st.composite
def _circles(draw: st.DrawFn) -> tuple[float, float, float]:
    return (draw(LATTICE), draw(LATTICE), draw(st.sampled_from((0.25, 0.5, 1.0, 1.5, 2.5))))


def _assert_twins(local: Stock, global_: Stock) -> None:
    assert local.exactly_equals(global_)
    assert local.representation_is_valid()
    assert local.arrangement_stats() == global_.arrangement_stats()


@settings(max_examples=EXAMPLES, deadline=None)
@given(segments=st.lists(_segments(), min_size=1, max_size=6), repeat=st.booleans(), radius=st.sampled_from(TOOL_RADII))
def test_local_line_sweep_equals_global(segments: list[tuple[float, float, float, float]], repeat: bool, radius: float) -> None:
    local, global_ = _rect(), _rect()
    for segment in segments + (segments[:1] if repeat else []):
        local.subtract_line_sweep_local(*segment, radius)
        global_.subtract_line_sweep(*segment, radius)
        _assert_twins(local, global_)


@settings(max_examples=EXAMPLES, deadline=None)
@given(circles=st.lists(_circles(), min_size=1, max_size=5), repeat=st.booleans(), radius=st.sampled_from(TOOL_RADII))
def test_local_circle_sweep_equals_global(circles: list[tuple[float, float, float]], repeat: bool, radius: float) -> None:
    local, global_ = _rect(), _rect()
    for cx, cy, rho in circles + (circles[:1] if repeat else []):
        local.subtract_circle_sweep_local(cx, cy, rho, radius)
        global_.subtract_circle_sweep(cx, cy, rho, radius)
        _assert_twins(local, global_)


@settings(max_examples=EXAMPLES, deadline=None)
@given(
    operations=st.lists(st.one_of(_segments().map(lambda s: ("line", s)), _circles().map(lambda c: ("circle", c))), min_size=2, max_size=7),
    radius=st.sampled_from(TOOL_RADII),
)
def test_mixed_local_sweeps_equal_global(operations: list[tuple[str, tuple[float, ...]]], radius: float) -> None:
    local, global_ = _rect(), _rect()
    for kind, values in operations:
        if kind == "line":
            local.subtract_line_sweep_local(*values, radius)
            global_.subtract_line_sweep(*values, radius)
        else:
            local.subtract_circle_sweep_local(*values, radius)
            global_.subtract_circle_sweep(*values, radius)
        _assert_twins(local, global_)


def test_four_wall_runs_locally_equal_the_global_runs() -> None:
    local, global_ = Stock(Polygon([[0, 0, 0], [20, 0, 0], [20, 12, 0], [0, 12, 0]])), Stock(Polygon([[0, 0, 0], [20, 0, 0], [20, 12, 0], [0, 12, 0]]))
    for y in (1.5, 3.5, 8.5, 10.5):
        local.subtract_line_sweep_local(2.0, y, 18.0, y, 1.0)
        global_.subtract_line_sweep(2.0, y, 18.0, y, 1.0)
        _assert_twins(local, global_)


def test_a_generic_double_oblique_sweep_stays_fast() -> None:
    """Witness against identically-zero exact decisions stalling the filter."""
    local, global_ = _rect(), _rect()
    start = time.perf_counter()
    for segment in ((0.3141592653589793, 1.2718281828459045, 11.1, 6.7071067811865475), (11.1, 6.7071067811865475, 1.4142135623730951, 7.3)):
        local.subtract_line_sweep_local(*segment, 1.0)
        global_.subtract_line_sweep(*segment, 1.0)
    elapsed = time.perf_counter() - start
    _assert_twins(local, global_)
    assert elapsed < GENERIC_DOUBLE_CEILING_S
