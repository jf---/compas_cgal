"""The U2 gate harness on a small stream, plus its report contract."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from benchmarks.held_annulus_store_gate import StoreCandidateMissError
from benchmarks.held_annulus_store_gate import StoreGateReport
from benchmarks.held_annulus_store_gate import ViewBoundBreakError
from benchmarks.held_annulus_store_gate import nested_breaks
from benchmarks.held_annulus_store_gate import read_box
from benchmarks.held_annulus_store_gate import selection_misses
from benchmarks.held_annulus_store_gate import store_gate
from compas_cgal import _stock_2
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import ReplayWorkload

CAP_RATIO = _stock_2.cap_chord_ratio(math.radians(80.0))


@pytest.fixture
def spiral() -> ReplayWorkload:
    circles = tuple(
        CircleSweep(0.25 * math.cos(k / 3.0), 0.25 * math.sin(k / 3.0), 0.5 + 0.2 * k, 0.25 * math.cos(k / 3.0) + 0.5 + 0.2 * k, 0.25 * math.sin(k / 3.0)) for k in range(24)
    )
    polygon = np.array([(-12.0, -12.0, 0.0), (12.0, -12.0, 0.0), (12.0, 12.0, 0.0), (-12.0, 12.0, 0.0)], dtype=np.float64)
    return ReplayWorkload(name="spiral", polygon=polygon, tool_radius=1.0, cap_chord_ratio=CAP_RATIO, circles=circles, segments=(), generation_s=0.0)


def test_read_box_holds_every_rim_including_the_entry_rim() -> None:
    circle = CircleSweep(1.0, 2.0, 3.0, 9.0, 2.0)
    box = read_box(circle, 1.0)
    assert box[2] >= 10.0
    assert box[0] <= 1.0 - 9.0


def test_selection_has_no_miss_on_a_stream(spiral: ReplayWorkload) -> None:
    misses, candidates = selection_misses(spiral)
    assert misses == 0
    assert candidates > 0


def test_nested_views_have_no_break_on_a_stream(spiral: ReplayWorkload) -> None:
    stations, breaks, first = nested_breaks(spiral, 1)
    assert stations == (len(spiral.circles) - 1) * 45
    assert breaks == 0
    assert first == []


def test_gate_report_round_trips_and_carries_its_load(spiral: ReplayWorkload) -> None:
    report = store_gate(spiral, every=4, repeats=1)
    report.require()
    record = json.loads(json.dumps(report.to_json()))
    assert record["schema"] == "held-annulus-store-gate/1"
    assert record["candidate_misses"] == 0
    assert record["stale_references"] > 0
    assert len(record["load_before"]) == 3


def _report(misses: int, breaks: int) -> StoreGateReport:
    return StoreGateReport("x", 1, 1, misses, 0.0, 0.0, 0.0, 0.0, 0, 0, breaks, (), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))


def test_require_raises_on_a_miss_and_on_a_break() -> None:
    with pytest.raises(StoreCandidateMissError):
        _report(1, 0).require()
    with pytest.raises(ViewBoundBreakError):
        _report(0, 1).require()
