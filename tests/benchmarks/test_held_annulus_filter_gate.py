"""The U3 gate harness on a small stream, plus its report contract."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from benchmarks.held_annulus_filter_gate import FilteredDecisionMismatchError
from benchmarks.held_annulus_filter_gate import FilterGateReport
from benchmarks.held_annulus_filter_gate import FilterUnsoundError
from benchmarks.held_annulus_filter_gate import filter_gate
from compas_cgal import _stock_2
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import ReplayWorkload

CAP_RATIO = _stock_2.cap_chord_ratio(math.radians(80.0))


@pytest.fixture
def stream() -> ReplayWorkload:
    circles = [CircleSweep(0.0, 0.0, 0.5, 0.5, 0.0)]
    for k in range(1, 30):
        rho = 0.5 + 0.1 * k
        circles.append(CircleSweep(0.0, 0.0, rho, rho, 0.0))
    polygon = np.array([(-12.0, -12.0, 0.0), (12.0, -12.0, 0.0), (12.0, 12.0, 0.0), (-12.0, 12.0, 0.0)], dtype=np.float64)
    return ReplayWorkload(name="stream", polygon=polygon, tool_radius=1.0, cap_chord_ratio=CAP_RATIO, circles=tuple(circles), segments=(), generation_s=0.0)


def test_gate_on_a_stream_is_sound_and_equal_to_the_exact_decision(stream: ReplayWorkload) -> None:
    report = filter_gate(stream, passes=1, ratio_every=7)
    report.require()
    assert report.decisions == len(stream.circles) - 1
    assert report.stations == report.decisions * 45
    assert report.soundness_compared > 0
    assert report.certified_last_one + report.certified_last_sixteen + report.certified_range + report.exact_stations == report.stations
    assert report.arrangement_over_filtered > 0.0


def test_gate_report_round_trips_and_carries_its_load(stream: ReplayWorkload) -> None:
    record = json.loads(json.dumps(filter_gate(stream, passes=1, ratio_every=7, check_soundness=False).to_json()))
    assert record["schema"] == "held-annulus-filter-gate/1"
    assert record["soundness_compared"] == 0
    assert len(record["load_after"]) == 3


def _report(unsound: int, breaks: int, mismatches: int) -> FilterGateReport:
    return FilterGateReport("x", 1, 1, 1, 0, 0, 0, 0, 0, 1, unsound, breaks, mismatches, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))


@pytest.mark.parametrize(("unsound", "breaks", "mismatches", "error"), [(1, 0, 0, FilterUnsoundError), (0, 1, 0, FilterUnsoundError), (0, 0, 1, FilteredDecisionMismatchError)])
def test_require_raises_on_each_failure_mode(unsound: int, breaks: int, mismatches: int, error: type[Exception]) -> None:
    with pytest.raises(error):
        _report(unsound, breaks, mismatches).require()
