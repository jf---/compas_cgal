"""The annulus store behind ``AnnulusView2``: append, truncate, candidate selection and the tiered views (unit U2 of the actual-stock refinement design).

The selection contract: every annulus that meets a query, by its exact outward box, is
returned once, whatever was appended and truncated before. The engine contract: a view
opened empty and fed by ``append`` decides what ``AnnulusView2.build`` decides on the
same rows, and what ``Stock2`` decides on the same replay.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _stock_2

CAP_RATIO = _stock_2.cap_chord_ratio(math.radians(80.0))
LOOP_PARAMETERS = [k / 45.0 for k in range(45)]
SQUARE = np.array([(-20.0, -20.0, 0.0), (20.0, -20.0, 0.0), (20.0, 20.0, 0.0), (-20.0, 20.0, 0.0)], dtype=np.float64)
# Largest coordinate of the query box and of the annulus centres in the Hypothesis streams (mm), past the part so the grid's border cells are exercised.
STREAM_REACH = 40.0
TIER_LAST = 16


def outward_box(cx: float, cy: float, rho: float, tool: float) -> tuple[float, float, float, float]:
    """The annulus's box, one ulp wider than any rounding of ``c +- (rho + r)``: the brute-force side of the selection test."""
    outer = math.nextafter(rho + tool, math.inf)
    return (
        math.nextafter(cx - outer, -math.inf),
        math.nextafter(cy - outer, -math.inf),
        math.nextafter(cx + outer, math.inf),
        math.nextafter(cy + outer, math.inf),
    )


def boxes_meet(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])


def open_view(annuli: list[tuple[float, float, float]], tool: float = 1.0) -> _stock_2.AnnulusView2:
    view = _stock_2.AnnulusView2.open(SQUARE, tool)
    for row in annuli:
        view.append(*row)
    return view


coordinate = st.floats(min_value=-STREAM_REACH, max_value=STREAM_REACH, allow_subnormal=False)
guide_radius = st.floats(min_value=0.0, max_value=12.0, allow_subnormal=False)
annulus = st.tuples(coordinate, coordinate, guide_radius)
near = st.floats(-8.0, 8.0, allow_subnormal=False)


@settings(max_examples=120, deadline=None)
@given(
    st.lists(annulus, max_size=60),
    st.integers(0, 60),
    st.lists(annulus, max_size=30),
    st.tuples(coordinate, coordinate),
    st.floats(min_value=0.0, max_value=20.0, allow_subnormal=False),
    st.integers(0, 90),
    st.integers(0, 90),
)
def test_h2_selection_finds_every_box_meeting_annulus_after_truncation(
    first_batch: list[tuple[float, float, float]],
    keep: int,
    second_batch: list[tuple[float, float, float]],
    query_centre: tuple[float, float],
    query_half: float,
    first: int,
    end: int,
) -> None:
    view = open_view(first_batch)
    kept = first_batch[: min(keep, len(first_batch))]
    view.truncate(len(kept))
    for row in second_batch:
        view.append(*row)
    rows = kept + second_batch
    assert view.size == len(rows)
    low, high = sorted((min(first, len(rows)), min(end, len(rows))))
    box = (query_centre[0] - query_half, query_centre[1] - query_half, query_centre[0] + query_half, query_centre[1] + query_half)
    found = view.candidates(low, high, *box)
    assert len(found) == len(set(found)), "each candidate once"
    assert all(low <= index < high for index in found), "the range filters at query time"
    must = {index for index in range(low, high) if boxes_meet(outward_box(*rows[index], 1.0), box)}
    assert must <= set(found), "a box that meets the query is never dropped"


def test_truncate_then_append_reuses_the_index_without_a_stale_hit() -> None:
    view = open_view([(10.0, 10.0, 1.0), (-10.0, -10.0, 1.0)])
    view.truncate(1)
    assert view.candidates(0, 1, -11.0, -11.0, -9.0, -9.0) == []
    assert view.stale_references > 0 or view.store_rebuilds > 0
    view.append(5.0, 5.0, 1.0)
    # Index 1 now names the new annulus: its old cell still references the index, its box does not meet.
    assert view.candidates(0, 2, -11.0, -11.0, -9.0, -9.0) == []
    assert view.candidates(0, 2, 4.0, 4.0, 6.0, 6.0) == [1]


def test_stale_references_compact_lazily() -> None:
    view = open_view([(float(k % 7), float(k % 5), 1.0) for k in range(200)])
    before = view.store_rebuilds
    view.truncate(150)
    assert view.store_rebuilds == before, "a small truncation leaves its references to the range filter"
    assert view.stale_references > 0
    view.truncate(10)
    assert view.store_rebuilds == before + 1, "stale references outnumbering live ones rebuild the grid"
    assert view.stale_references == 0
    assert sorted(view.candidates(0, 10, -30.0, -30.0, 30.0, 30.0)) == list(range(10)), "a rebuilt grid keeps every live annulus"


def test_an_index_reused_after_truncation_has_no_stale_exact_annulus() -> None:
    """The exact path caches an annulus's exact form by index: a truncation has to forget it."""
    view = open_view([(0.0, 0.0, 4.0)])
    args = (0.0, 0.0, 4.0, 4.0, 0.0, False, [0.0], CAP_RATIO, 0.0)
    assert view.circle_exceeded(0, 1, *args) is False
    view.truncate(0)
    view.append(100.0, 100.0, 1.0)
    assert view.circle_exceeded(0, 1, *args) is True, "the new index 0 is far away: the rim is a closed loop of material"
    rows = [(0.0, 0.0, 4.0), (2.5, 0.0, 1.5)]
    view = open_view(rows[:1])
    tangent = (6.0, 0.0, 1.0, 5.0, 0.0, False, [0.0], CAP_RATIO, 0.0)
    view.circle_prefixes(0, 1, *tangent)
    view.truncate(0)
    view.append(*rows[1])
    built = _stock_2.AnnulusView2.build(SQUARE, np.array(rows[1:], dtype=np.float64), 1.0)
    assert view.circle_prefixes(0, 1, *tangent) == built.circle_prefixes(0, 1, *tangent)


def test_the_tiers_are_the_last_one_and_the_last_sixteen() -> None:
    covering = [(0.0, 0.0, 4.0)]
    near = (0.0, 0.0, 4.0, 4.0, 0.0, [0.0], CAP_RATIO)
    far = (30.0, 30.0, 1.0)
    for trailing, tier in ((0, 0), (1, 1), (15, 1), (16, 2)):
        report = open_view(covering + [far] * trailing, tool=1.0).circle_filter_report(0, 1 + trailing, *near)
        assert report[0][0] == tier, f"{trailing} later emissions"
    last = open_view([far, *covering], tool=1.0).circle_filter_report(0, 2, *near)
    assert last[0][0] == 0, "the last emission alone is the first tier"
    sixteenth = open_view([far, *covering, *([far] * 15)], tool=1.0).circle_filter_report(0, 17, *near)
    assert sixteenth[0][0] == 1, "the sixteenth emission from the end is in the second tier"


def test_truncation_beyond_the_store_raises() -> None:
    view = open_view([(0.0, 0.0, 1.0)])
    with pytest.raises(_stock_2.InvalidAnnulusViewRangeError):
        view.truncate(2)


def test_append_on_a_built_view_raises_its_own_error() -> None:
    view = _stock_2.AnnulusView2.build(SQUARE, np.array([(0.0, 0.0, 1.0)], dtype=np.float64), 1.0)
    with pytest.raises(_stock_2.AnnulusViewNotAppendableError):
        view.append(1.0, 1.0, 1.0)
    with pytest.raises(_stock_2.AnnulusViewNotAppendableError):
        view.truncate(0)


@pytest.mark.parametrize(("row", "error"), [((float("nan"), 0.0, 1.0), "NonFiniteAnnulusViewInputError"), ((0.0, 0.0, -1.0), "InvalidAnnulusViewRadiiError")])
def test_invalid_append_raises(row: tuple[float, float, float], error: str) -> None:
    with pytest.raises(getattr(_stock_2, error)):
        open_view([]).append(*row)


def test_coincident_centres_stay_separate_annuli() -> None:
    """Concentric circles share a centre: a centre-keyed index would lose one."""
    view = open_view([(0.0, 0.0, 2.0), (0.0, 0.0, 6.0)])
    assert sorted(view.candidates(0, 2, -1.0, -1.0, 1.0, 1.0)) == [0, 1]


@settings(max_examples=60, deadline=None)
@given(st.lists(st.tuples(st.integers(-9, 9), st.integers(-9, 9), st.integers(0, 6)), max_size=12), st.integers(-9, 9), st.integers(-9, 9), st.integers(1, 6))
def test_opened_view_decides_what_a_built_view_and_the_arrangement_decide(annuli: list[tuple[int, int, int]], x: int, y: int, rho: int) -> None:
    rows = [(float(a), float(b), float(c)) for a, b, c in annuli]
    opened = open_view(rows)
    built = _stock_2.AnnulusView2.build(SQUARE, np.array(rows, dtype=np.float64).reshape(-1, 3), 1.0)
    stock = _stock_2.Stock2(SQUARE, [])
    for row in rows:
        stock.subtract_circle_sweep(*row, 1.0)
    args = (float(x), float(y), float(rho), float(x + rho), float(y), False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    verdict = opened.circle_exceeded(0, len(rows), *args)
    assert verdict == built.circle_exceeded(0, len(rows), *args)
    assert all(a and b for a, b in opened.circle_matches_stock(stock, 0, len(rows), *args))


@settings(max_examples=40, deadline=None)
@given(
    st.lists(
        st.tuples(st.floats(-8.0, 8.0, allow_subnormal=False), st.floats(-8.0, 8.0, allow_subnormal=False), st.floats(0.0, 6.0, allow_subnormal=False)), min_size=1, max_size=24
    ),
    near,
    near,
    st.floats(0.25, 6.0, allow_subnormal=False),
)
def test_h2_nested_views_only_shrink_the_max_run_and_the_tiered_decision_is_the_full_view(annuli: list[tuple[float, float, float]], x: float, y: float, rho: float) -> None:
    view = open_view(annuli)
    end = len(annuli)
    args = (x, y, rho, x + rho, y, False, LOOP_PARAMETERS, CAP_RATIO, 0.0)
    last_one = view.circle_prefixes(end - 1, end, *args)
    last_sixteen = view.circle_prefixes(max(0, end - TIER_LAST), end, *args)
    every = view.circle_prefixes(0, end, *args)
    for wide, middle, narrow in zip(last_one, last_sixteen, every):
        assert not narrow[2] or (middle[2] and wide[2]), "a pass on a larger view is a pass on the smaller"
        assert wide[1] + 1e-9 >= middle[1] and middle[1] + 1e-9 >= narrow[1]
    assert view.circle_exceeded_tiered(0, end, *args) == view.circle_exceeded(0, end, *args)
