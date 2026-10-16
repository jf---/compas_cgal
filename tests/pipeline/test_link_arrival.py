"""A link-reached chain's root loop is part of the link's admissibility (A-I1).

`walk_chain(..., emit_start_loop=True)` cuts the arriving root loop with no
ladder in front of it, so `_link_admissible` is that loop's only cap decision:
a link whose pieces pass but whose arriving loop exceeds the cap is refused,
and the chain falls through to another offer or a forced entry.
"""

from __future__ import annotations

import math
from types import SimpleNamespace

import pytest
from compas.geometry import Polygon

from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline import generate
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.linking import StraightLink
from compas_cgal.stock import Stock
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

TOOL = ToolRadius.build(1.0)
RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
CAP_DEG = 40.0
# A cleared disk 2.5 tool radii wide: the link chord and a loop centred on it sweep inside it.
CLEARED_CENTRE = (-3.0, 0.0)
CLEARED_RADIUS = 2.5
# The arriving loop: guide radius 0.5, phase pi puts its entry on the cleared (-x) side.
LOOP_RADIUS = 0.5
LOOP_PHASE = Radian(math.pi)
# Centred at x = -1.5, the loop's entry (-2, 0) is cleared but its sweep reaches
# x = 0, past the cleared disk's edge at x = -0.5.
OVER_CAP_CENTRE = (-1.5, 0.0)
DEPARTURE_OFFSET = 0.5


def _path() -> generate.PathState:
    stock = Stock(RECT_12X8)
    stock.subtract_disk(*CLEARED_CENTRE, CLEARED_RADIUS)
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(math.radians(CAP_DEG))), centre_domain(RECT_12X8, [], TOOL))
    return generate.PathState(stock=stock, admissibility=admissibility, cut=CutZ.build(0.0), clearance=ClearanceZ.build(5.0))


def _arrival(cx: float, cy: float) -> tuple[StraightLink, SimpleNamespace]:
    loop = IncidentLoop.build(Point2[WorldXY].build(cx, cy), GuideRadius.build(LOOP_RADIUS), LOOP_PHASE)
    chain = SimpleNamespace(first=0, loops=[loop], entry_station=SimpleNamespace(clockwise=True))
    departure = Point2[WorldXY].build(CLEARED_CENTRE[0] + DEPARTURE_OFFSET, CLEARED_CENTRE[1])
    return StraightLink(next_loop=loop, departure=departure), chain


def test_a_link_into_an_over_cap_root_loop_is_refused() -> None:
    path = _path()
    link, chain = _arrival(*OVER_CAP_CENTRE)
    # The premise: every piece passes; only the arriving loop exceeds the cap.
    assert all(generate._piece_admissible(path.admissibility, path.stock, piece, True) for piece in link.pieces)
    assert not path.admissibility.loop_within_cap(path.stock, chain.loops[0], True)
    assert not generate._link_admissible(path, link, chain)  # type: ignore[arg-type]


def test_a_link_into_a_cleared_root_loop_is_admitted() -> None:
    path = _path()
    link, chain = _arrival(*CLEARED_CENTRE)
    assert generate._link_admissible(path, link, chain)  # type: ignore[arg-type]


def test_a_refused_arrival_falls_through_to_a_forced_entry(monkeypatch: pytest.MonkeyPatch) -> None:
    path = _path()
    link, chain = _arrival(*OVER_CAP_CENTRE)
    monkeypatch.setattr(generate, "link_between", lambda *_args, **_kwargs: link)
    current = generate.WalkEnd(loop=chain.loops[0], clockwise=True, resume=None)
    choice = generate.choose_next(path, current, [(chain,)])  # type: ignore[list-item]
    assert choice.link is None
