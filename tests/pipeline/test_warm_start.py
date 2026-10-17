"""Warm-started evaluation order of the at-rest ring (Task 1.2 lever L6): cost only, never a verdict."""

from __future__ import annotations

import math
from collections import Counter

import pytest
from compas.geometry import Polygon
from hypothesis import event
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _stock_2
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline import admissibility as admissibility_module
from compas_cgal.pipeline.admissibility import GAP_CLOSE_RATIO
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import ring_positions
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.ladder import WarmStart
from compas_cgal.stock import Stock
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
ADVANCE = Direction2[WorldXY].build(1.0, 0.0)


class _CountingStock2:
    """`_stock_2` as `admissibility` sees it, counting the at-rest ring's native probes."""

    def __init__(self) -> None:
        self.calls: Counter[str] = Counter()

    def __getattr__(self, name: str) -> object:
        return getattr(_stock_2, name)

    def engagement_at(self, *args: object) -> tuple[float, float, bool]:
        self.calls["engagement_at"] += 1
        return _stock_2.engagement_at(*args)  # type: ignore[arg-type]


def _admissibility(cap_deg: float) -> Admissibility:
    return Admissibility.build(TOOL, cap_ratio(Radian(math.radians(cap_deg))), centre_domain(RECT_12X8, [], TOOL))


def _entry_cleared_stock() -> Stock:
    """The loop's entry sits in a cleared disk; the ring first bites five positions round (fixture probed: index 5 at 120 deg)."""
    stock = Stock(RECT_12X8)
    stock.subtract_disk(2.0, 0.0, 2.0)
    return stock


def _loop() -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(2.0), Radian(0.0))


def test_a_repeated_ring_refusal_is_asked_at_its_refusing_position_first(monkeypatch: pytest.MonkeyPatch) -> None:
    """The second identical refusal costs one native probe: the ring starts where it last refused."""
    counting = _CountingStock2()
    monkeypatch.setattr(admissibility_module, "_stock_2", counting)
    adm = _admissibility(120.0)
    stock, loop = _entry_cleared_stock(), _loop()
    assert not adm.loop_within_cap_at_rest(stock, loop, ADVANCE)
    cold = counting.calls["engagement_at"]
    assert cold > 1  # non-vacuous: the natural order does not refuse at its first position
    counting.calls.clear()
    assert not adm.loop_within_cap_at_rest(stock, loop, ADVANCE)
    assert counting.calls["engagement_at"] == 1


@given(verdicts=st.lists(st.booleans(), max_size=40), hint=st.integers(min_value=-3, max_value=43))
def test_warm_start_asks_each_probe_once_hint_first_and_finds_a_refusal_iff_one_exists(verdicts: list[bool], hint: int) -> None:
    order = WarmStart(hint=hint)
    asked: list[int] = []

    def refuses(index: int) -> bool:
        asked.append(index)
        return verdicts[index]

    found = order.first_refusal(len(verdicts), refuses)
    assert len(asked) == len(set(asked))
    if 0 <= hint < len(verdicts):
        assert asked[0] == hint
    if found is None:
        assert not any(verdicts)
        assert sorted(asked) == list(range(len(verdicts)))
        assert order.hint == hint
    else:
        assert verdicts[found] and asked[-1] == found and order.hint == found
    warm_order = ([hint] if 0 <= hint < len(verdicts) else []) + [index for index in range(len(verdicts)) if index != hint]
    assert asked == warm_order[: len(asked)]


@given(
    verdicts=st.lists(st.booleans(), min_size=1, max_size=40),
    hint=st.integers(min_value=-3, max_value=43),
    rewrites=st.lists(st.integers(min_value=-3, max_value=43), min_size=40, max_size=40),
)
def test_warm_start_survives_its_hint_being_rewritten_mid_evaluation(verdicts: list[bool], hint: int, rewrites: list[int]) -> None:
    """An interleaved evaluation rewriting the shared hint mid-scan neither skips nor repeats a probe."""
    order = WarmStart(hint=hint)
    asked: list[int] = []

    def refuses(index: int) -> bool:
        asked.append(index)
        order.hint = rewrites[len(asked) - 1]  # another evaluation sharing the hint moved it
        return verdicts[index]

    found = order.first_refusal(len(verdicts), refuses)
    assert len(asked) == len(set(asked))
    if found is None:
        assert not any(verdicts) and sorted(asked) == list(range(len(verdicts)))
    else:
        assert verdicts[found] and asked[-1] == found


def _natural_order_ring(adm: Admissibility, stock: Stock, loop: IncidentLoop, advance: Direction2[WorldXY]) -> bool:
    """HEAD's ring verdict: every position in natural order, entry first."""
    raw, ratio = stock.raw, float(adm.cap_ratio)
    return not any(_stock_2.engagement_at(raw, x, y, TOOL.value, ratio, GAP_CLOSE_RATIO)[2] for x, y in ring_positions(loop, advance))


SHARED_120 = _admissibility(120.0)
SHARED_40 = _admissibility(40.0)

_coordinate = st.floats(min_value=-1.5, max_value=1.5, allow_nan=False, allow_infinity=False)
_candidate = st.tuples(
    st.floats(min_value=-3.0, max_value=3.0, allow_nan=False, allow_infinity=False),
    _coordinate,
    st.floats(min_value=1.05, max_value=2.4, allow_nan=False, allow_infinity=False),
    st.floats(min_value=0.0, max_value=6.28, allow_nan=False, allow_infinity=False),
    st.floats(min_value=0.0, max_value=6.28, allow_nan=False, allow_infinity=False),
)


# Outcomes the equivalence property saw, so its non-vacuity is asserted, not hoped for.
RING_OUTCOMES: Counter[bool] = Counter()


@settings(max_examples=100, deadline=None)
@given(
    cleared=st.lists(st.tuples(st.floats(-5.0, 5.0), st.floats(-3.0, 3.0), st.floats(0.5, 3.0)), max_size=3),
    candidates=st.lists(_candidate, min_size=1, max_size=6),
    hint=st.integers(min_value=-1, max_value=40),
    cap_40=st.booleans(),
    clear_own=st.booleans(),
)
def test_the_warm_started_ring_decides_what_the_natural_order_ring_decides(
    cleared: list[tuple[float, float, float]], candidates: list[tuple[float, float, float, float, float]], hint: int, cap_40: bool, clear_own: bool
) -> None:
    """Equivalence over many candidates, with the hint carried across them and seeded anywhere (out of range included)."""
    stock = Stock(RECT_12X8)
    for cx, cy, radius in cleared:
        stock.subtract_disk(cx, cy, radius)
    if clear_own:  # each candidate's own annulus cut first: the admitting side, often
        for cx, cy, rho, _phase, _heading in candidates:
            stock.subtract_circle_sweep_local(cx, cy, rho, TOOL.value)
    adm = SHARED_40 if cap_40 else SHARED_120
    adm.ring_order.hint = hint
    for cx, cy, rho, phase, heading in candidates:
        loop = IncidentLoop.build(Point2[WorldXY].build(cx, cy), GuideRadius.build(rho), Radian(phase))
        advance = Direction2[WorldXY].build(math.cos(heading), math.sin(heading))
        before = adm.ring_order.hint
        warm = adm.loop_within_cap_at_rest(stock, loop, advance)
        assert warm == _natural_order_ring(adm, stock, loop, advance)
        event("admit" if warm else "refuse")
        RING_OUTCOMES[warm] += 1
        if warm:
            assert adm.ring_order.hint == before
        else:
            x, y = ring_positions(loop, advance)[adm.ring_order.hint]
            assert _stock_2.engagement_at(stock.raw, x, y, TOOL.value, float(adm.cap_ratio), GAP_CLOSE_RATIO)[2]


def test_the_equivalence_property_exercised_both_admit_and_refuse() -> None:
    """M3: the property above is run here and must have compared both outcomes (Hypothesis also reports them as events)."""
    RING_OUTCOMES.clear()
    test_the_warm_started_ring_decides_what_the_natural_order_ring_decides()
    assert RING_OUTCOMES[True] > 0 and RING_OUTCOMES[False] > 0, dict(RING_OUTCOMES)


def test_the_warm_start_is_not_part_of_the_admissibility_value() -> None:
    """Two admissibilities differing only in their hint compare equal: the hint is cost, not a verdict input."""
    a = _admissibility(120.0)
    b = Admissibility(tool=a.tool, cap_ratio=a.cap_ratio, domain=a.domain)
    b.ring_order.hint = 7
    assert a == b


def test_an_admitted_ring_asks_every_position_and_keeps_the_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    """The admitting side of the equivalence, non-vacuously: the loop's own annulus is cut, so no position exceeds."""
    counting = _CountingStock2()
    monkeypatch.setattr(admissibility_module, "_stock_2", counting)
    adm, loop = _admissibility(120.0), _loop()
    stock = Stock(RECT_12X8)
    stock.subtract_circle_sweep_local(loop.centre.x, loop.centre.y, loop.radius.value, TOOL.value)
    adm.ring_order.hint = 5
    assert adm.loop_within_cap_at_rest(stock, loop, ADVANCE)
    assert counting.calls["engagement_at"] == len(ring_positions(loop, ADVANCE))
    assert adm.ring_order.hint == 5
    assert _natural_order_ring(adm, stock, loop, ADVANCE)
