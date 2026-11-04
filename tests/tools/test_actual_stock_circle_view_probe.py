"""Circle-view design probe (ported from build/actual-stock-design/probe_circles_only.py): rows, gap statistics, record round trip."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from compas_cgal import _stock_2
from tools.actual_stock_circle_view_probe import SCHEMA
from tools.actual_stock_circle_view_probe import MalformedProbeRecordError
from tools.actual_stock_circle_view_probe import ProbeRecord
from tools.actual_stock_circle_view_probe import probe
from tools.actual_stock_circle_view_probe import render_markdown
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import ReplayWorkload
from tools.actual_stock_profile import SegmentSweep

TOOL = 1.0
CAP_DEG = 80.0


def _workload() -> ReplayWorkload:
    """Circle 0 big; circle 1 inside a band its two connectors cleared; circle 2 repeats circle 0 (reads only cut material)."""
    band = (SegmentSweep(0, 4.0, 0.5, 12.0, 0.5), SegmentSweep(0, 4.0, -0.5, 12.0, -0.5))
    circles = (CircleSweep(0.0, 0.0, 3.0, 3.0, 0.0), CircleSweep(8.0, 0.0, 0.25, 8.25, 0.0), CircleSweep(0.0, 0.0, 3.0, 3.0, 0.0))
    polygon = np.array([(-10.0, -10.0, 0.0), (50.0, -10.0, 0.0), (50.0, 10.0, 0.0), (-10.0, 10.0, 0.0)], dtype=np.float64)
    return ReplayWorkload("synthetic", polygon, TOOL, _stock_2.cap_chord_ratio(math.radians(CAP_DEG)), circles, band, 0.0)


def test_probe_rows_compare_the_full_window_on_both_views() -> None:
    record = probe(_workload(), stride=1)
    assert [r.circle for r in record.rows] == [0, 1, 2]
    plunge, cleared, repeat = record.rows
    assert plunge.actual_exceeded and plunge.circle_view_exceeded
    assert not cleared.actual_exceeded and cleared.circle_view_exceeded
    assert cleared.circle_view_deg > cleared.actual_deg
    # Circle 0's annulus is proven exactly to miss circle 1's window; only its two connectors reach it.
    assert cleared.circles == 0 and cleared.sweeps == 2
    # Circle 1 (outer radius 1.25 at x = 8) misses circle 2's window (radius 4.5 at the origin); circle 0 and both connectors reach it.
    assert repeat.circles == 1 and repeat.sweeps == 3
    assert repeat.circle_view_deg >= repeat.actual_deg
    assert record.circle_view_only_refusals == (1,)
    assert record.gap_max_deg == pytest.approx(cleared.circle_view_deg - cleared.actual_deg)


def test_probe_stride_samples_every_kth_circle() -> None:
    assert [r.circle for r in probe(_workload(), stride=2).rows] == [0, 2]


def test_record_round_trips_through_json() -> None:
    record = probe(_workload(), stride=1)
    loaded = ProbeRecord.from_json(json.loads(json.dumps(record.to_json())))
    assert loaded == record
    assert "circle-view-only refusals: [1]" in render_markdown(loaded)


def test_record_with_a_foreign_schema_is_refused() -> None:
    raw = probe(_workload(), stride=1).to_json()
    with pytest.raises(MalformedProbeRecordError):
        ProbeRecord.from_json({**raw, "schema": SCHEMA + "-other"})
    with pytest.raises(MalformedProbeRecordError):
        ProbeRecord.from_json({k: v for k, v in raw.items() if k != "rows"})
