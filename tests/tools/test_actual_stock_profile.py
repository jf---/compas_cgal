"""Contracts of tools.actual_stock_profile: the replay measurements on a small exact workload, the record schema and the generated markdown."""

from __future__ import annotations

import dataclasses
import json
import math

import numpy as np
import pytest

from benchmarks.performance_record import BuildIdentity
from benchmarks.performance_record import MachineIdentity
from compas_cgal import _stock_2
from tools.actual_stock_profile import ENGINE_BUDGETS
from tools.actual_stock_profile import GENERATED_BY
from tools.actual_stock_profile import JSON_PATH
from tools.actual_stock_profile import MD_PATH
from tools.actual_stock_profile import MONSTERA
from tools.actual_stock_profile import SCHEMA
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import CpuMs
from tools.actual_stock_profile import CpuS
from tools.actual_stock_profile import EngineBudget
from tools.actual_stock_profile import MalformedProfileRecordError
from tools.actual_stock_profile import ProfileRecord
from tools.actual_stock_profile import ReplayWorkload
from tools.actual_stock_profile import SegmentSweep
from tools.actual_stock_profile import render_markdown
from tools.actual_stock_profile import run_budget
from tools.actual_stock_profile import run_sequential
from tools.actual_stock_profile import run_windows

# A 24 mm square pocket, 1 mm tool, guide radius 2 mm, centres 1.5 mm apart along x: every circle overlaps its predecessors.
SIDE = 24.0
TOOL = 1.0
GUIDE = 2.0
STEP = 1.5
CIRCLE_COUNT = 8
CAP_DEG = 80.0
GENEROUS_BUDGET = EngineBudget(total_s=CpuS(3600.0), per_sweep_ms=None, per_query_ms=None)
ZERO_BUDGET = EngineBudget(total_s=CpuS(0.0), per_sweep_ms=None, per_query_ms=None)


def _workload() -> ReplayWorkload:
    polygon = np.array([(0.0, 0.0, 0.0), (SIDE, 0.0, 0.0), (SIDE, SIDE, 0.0), (0.0, SIDE, 0.0)], dtype=np.float64)
    centre_y = SIDE / 2
    circles = tuple(CircleSweep(4.0 + STEP * i, centre_y, GUIDE, 4.0 + STEP * i, centre_y - GUIDE) for i in range(CIRCLE_COUNT))
    segments = tuple(SegmentSweep(i, a.entry_x, a.entry_y, b.entry_x, b.entry_y) for i, (a, b) in enumerate(zip(circles, circles[1:])))
    return ReplayWorkload("square", polygon, TOOL, _stock_2.cap_chord_ratio(math.radians(CAP_DEG)), circles, segments, 0.0)


@pytest.fixture(scope="module")
def record() -> ProfileRecord:
    w = _workload()
    return ProfileRecord(
        generated_by=GENERATED_BY,
        schema=SCHEMA,
        machine=MachineIdentity("arm64", "Darwin", "test cpu", "3.12"),
        build=BuildIdentity("0" * 40, False, "0", "2026-10-21T00:00:00+00:00"),
        workload=w.name,
        circles=len(w.circles),
        segments=len(w.segments),
        stations=45,
        generation_s=0.0,
        load_start=(1.0, 2.0, 3.0),
        load_end=(1.0, 2.0, 3.0),
        sequential=run_sequential(w, CIRCLE_COUNT),
        windows=run_windows(w, 1),
        budget_runs=(dataclasses.replace(run_budget(w, GENEROUS_BUDGET), workload=MONSTERA, budget=ENGINE_BUDGETS[MONSTERA]),),
    )


def test_sequential_replay_queries_each_circle_before_its_own_sweep(record: ProfileRecord) -> None:
    rows = record.sequential
    assert [r.circle for r in rows] == list(range(CIRCLE_COUNT))
    assert [r.segments for r in rows] == [1] * (CIRCLE_COUNT - 1) + [0]
    # The first circle meets uncut stock: its decision runs on the bare polygon.
    assert rows[0].halfedges == 8
    assert rows[-1].halfedges > rows[0].halfedges


def test_window_views_nest_and_never_break_the_bound(record: ProfileRecord) -> None:
    assert len(record.windows) == CIRCLE_COUNT - 1
    for row in record.windows:
        assert row.last1.sweeps <= row.last16.sweeps <= row.exact.sweeps
        assert not row.exact.query.exceeded or (row.last16.query.exceeded and row.last1.query.exceeded)


def test_budget_run_completes_under_a_generous_budget() -> None:
    w = _workload()
    run = run_budget(w, GENEROUS_BUDGET)
    assert run.completed
    assert (run.sweeps_done, run.queries_done) == (len(w.circles) + len(w.segments), len(w.circles))
    assert run.failures() == []


def test_budget_run_stops_at_the_budget_and_is_red() -> None:
    run = run_budget(_workload(), ZERO_BUDGET)
    assert not run.completed
    assert (run.sweeps_done, run.queries_done) == (0, 1)
    assert run.failures() and "stopped at the budget" in run.failures()[0]


def test_per_item_budget_lines_are_checked(record: ProfileRecord) -> None:
    run = dataclasses.replace(record.budget_runs[0], sweep_cpu_s=CpuS(1.0), query_cpu_s=CpuS(1.0), completed=True)
    failures = run.failures()
    assert any("per sweep" in f for f in failures) and any("per query" in f for f in failures)


def test_record_round_trips_through_json(record: ProfileRecord) -> None:
    text = record.to_json()
    assert ProfileRecord.from_json(text).to_json() == text


def test_wrong_schema_is_refused(record: ProfileRecord) -> None:
    raw = json.loads(record.to_json())
    raw["schema"] = "actual-stock-profile/0"
    with pytest.raises(MalformedProfileRecordError, match="schema"):
        ProfileRecord.from_json(json.dumps(raw))


def test_missing_field_is_refused(record: ProfileRecord) -> None:
    raw = json.loads(record.to_json())
    del raw["sequential"][0]["query"]["exceeded"]
    with pytest.raises(MalformedProfileRecordError, match="fields"):
        ProfileRecord.from_json(json.dumps(raw))


def test_markdown_is_marked_generated_and_has_every_table(record: ProfileRecord) -> None:
    md = render_markdown(record)
    assert md.startswith(f"<!-- {GENERATED_BY} -->")
    for heading in ("**Engine budget, measured**", "**Sequential replay**", "**Windowed local stocks**", "**Last-16 bound vs exact**", "**Labelled-engine projection"):
        assert heading in md


def test_committed_markdown_is_rendered_from_the_committed_record() -> None:
    record = ProfileRecord.from_json(JSON_PATH.read_text(encoding="utf-8"))
    assert MD_PATH.read_text(encoding="utf-8") == render_markdown(record)


def test_engine_budgets_are_the_approved_ones() -> None:
    assert ENGINE_BUDGETS[MONSTERA] == EngineBudget(total_s=CpuS(5.0), per_sweep_ms=CpuMs(0.25), per_query_ms=CpuMs(0.6))
    assert ENGINE_BUDGETS["figure5_refined"].total_s == 2.5
