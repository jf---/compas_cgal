"""Indexed placement lanes equal the pre-lever ``_continuous_placement_lanes`` (verbatim oracle).

The overlap margin is replaced, on both sides, by a recorded deterministic stub
over lattice inputs (signs -, 0, +), so the property sees ties, rejections and
merges often, and pins that the margin is evaluated on the same pairs in the
same order; the full-path digests cover the real margin.
"""

from __future__ import annotations

from fractions import Fraction
from types import ModuleType

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import tests.benchmarks.held_figure5_path_pre_lever_oracle as oracle
from benchmarks import held_figure5_path
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

SEGMENTS = 3


@st.composite
def _circles(draw: st.DrawFn) -> tuple[Figure5CounterclockwiseCircle, ...]:
    circles = []
    for _ in range(draw(st.integers(min_value=0, max_value=40))):
        segment = draw(st.integers(min_value=0, max_value=SEGMENTS - 1))
        site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(segment), Fraction(draw(st.integers(min_value=0, max_value=7)), 8), SEGMENTS)
        circles.append(
            Figure5CounterclockwiseCircle(
                GuideRunId(draw(st.integers(min_value=0, max_value=5))),
                GuideRunStationOrdinal(draw(st.integers(min_value=0, max_value=5))),
                site,
                site,
                Point2[WorldXY].build(draw(st.integers(min_value=0, max_value=3)) / 2, draw(st.integers(min_value=0, max_value=3)) / 2),
                Point2[WorldXY].build(0.0, 0.0),
                GuideRadius.build(draw(st.sampled_from((0.5, 1.0)))),
            )
        )
    return tuple(circles)


def _stubbed(module: ModuleType, monkeypatch: pytest.MonkeyPatch, calls: list[tuple[object, object]]) -> None:
    def margin(tail: Figure5CounterclockwiseCircle, head: Figure5CounterclockwiseCircle, tool_radius: ToolRadius) -> float:
        calls.append((tail, head))
        return float((int(tail.station_ordinal) * 7 + int(head.station_ordinal) * 3 + int(tail.run_id)) % 3 - 1)

    monkeypatch.setattr(module, "predecessor_overlap_margin", margin)
    monkeypatch.setattr(module, "_paper_candidate", lambda circle: circle)


@settings(max_examples=1500, deadline=None)
@given(_circles())
def test_indexed_lanes_equal_the_pre_lever_function(circles: tuple[Figure5CounterclockwiseCircle, ...]) -> None:
    with pytest.MonkeyPatch.context() as monkeypatch:
        old_calls: list[tuple[object, object]] = []
        new_calls: list[tuple[object, object]] = []
        _stubbed(oracle, monkeypatch, old_calls)
        _stubbed(held_figure5_path, monkeypatch, new_calls)
        tool_radius = ToolRadius.build(1.0)
        old = oracle._continuous_placement_lanes(circles, tool_radius)
        new = held_figure5_path._continuous_placement_lanes(circles, tool_radius)
    # The oracle's _PlacementLane is its own class: compare the full repr.
    assert repr(new) == repr(old)
    assert new_calls == old_calls
