"""A chain's turn direction is uniform: a mixed-turn chain raises by name, never a bare assert.

`link_between` returns `NoLink` for loops turning opposite ways before it looks
at ``tangent_entries``, so a consecutive-station link exists only on a
uniform-turn chain. Both sites that build one -- the ladder's bridge condition
and the walk's emitted link -- raise `MixedTurnChainError` instead.
"""

from __future__ import annotations

import dataclasses

import pytest
from compas.geometry import Polygon

from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import MixedTurnChainError
from compas_cgal.pipeline.generate import LoopConditions
from compas_cgal.pipeline.generate import PathState
from compas_cgal.pipeline.generate import enter_first_cut
from compas_cgal.pipeline.generate import walk_chain
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.stock import Stock
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)
# Stations 1.2 apart along x, guide radius 2: three floor-clearing loops inside RECT_12X8.
STATION_PITCH = 1.2
STATION_X0 = -4.2
GUIDE_RADIUS = 2.0
# A permissive cap (radians, below pi) so only the turn rule is under test.
PERMISSIVE_CAP = Radian(3.0)
WALK_WINDOW = 3


def _chain(turns: list[bool]) -> OrientedChain:
    stations = tuple(
        GuideStation(
            centre=Point2[WorldXY].build(STATION_X0 + STATION_PITCH * index, 0.0),
            radius=GuideRadius.build(GUIDE_RADIUS),
            tangent=Direction2[WorldXY].build(1.0, 0.0),
            clockwise=clockwise,
        )
        for index, clockwise in enumerate(turns)
    )
    chain = OrientedChain.build(GuideChain(stations=stations, origin=0), TOOL)
    assert chain is not None
    return chain


def _path() -> tuple[PathState, Admissibility]:
    admissibility = Admissibility.build(TOOL, cap_ratio(PERMISSIVE_CAP), centre_domain(RECT_12X8, [], TOOL))
    return PathState(stock=Stock(RECT_12X8), admissibility=admissibility, cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0)), admissibility


def test_the_bridge_condition_raises_on_a_mixed_turn_pair() -> None:
    chain = _chain([True, True, True])
    _state, admissibility = _path()
    origin, candidate = chain.candidate(chain.first), chain.candidate(chain.first + 1)
    flipped = dataclasses.replace(candidate, station=dataclasses.replace(candidate.station, clockwise=not candidate.station.clockwise))
    with pytest.raises(MixedTurnChainError):
        LoopConditions(admissibility).bridge_admissible(Stock(RECT_12X8), origin, flipped)


def test_the_walk_raises_on_a_mixed_turn_chain(monkeypatch: pytest.MonkeyPatch) -> None:
    """The walk's own link: the ring refuses, so the advance is forced without asking the bridge."""
    chain = _chain([True, False, True])
    state, admissibility = _path()
    monkeypatch.setattr(LoopConditions, "loop_admissible_at_rest", lambda *_args: False)
    enter_first_cut(state, chain, chain.first)
    with pytest.raises(MixedTurnChainError):
        walk_chain(state, chain, chain.first, WALK_WINDOW, LoopConditions(admissibility), emit_start_loop=False)
