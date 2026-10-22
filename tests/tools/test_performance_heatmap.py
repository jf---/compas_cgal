"""The performance figures' arithmetic: buckets partition the profile, lever cells read named code, renames fail loud."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import pytest

import benchmarks.held_figure5_path as figure5_path
import benchmarks.held_figure5_raw_guide as raw_guide
from tools.performance_callpaths import CALLPATH_BATCH
from tools.performance_callpaths import CALLPATH_HYPOTHESES
from tools.performance_callpaths import CLEARANCE_DECODE
from tools.performance_callpaths import CLEARANCE_MODULE
from tools.performance_callpaths import CLEARANCE_QUERY
from tools.performance_callpaths import FRACTION_MODULE
from tools.performance_callpaths import FRACTION_OPERATOR
from tools.performance_callpaths import GCD
from tools.performance_callpaths import GCD_MODULE
from tools.performance_callpaths import callpath_rows
from tools.performance_heatmap import HELD_CASE_LOAD
from tools.performance_heatmap import HELD_CLEARANCE
from tools.performance_heatmap import HELD_HYPOTHESES
from tools.performance_heatmap import HELD_PLACEMENT
from tools.performance_heatmap import Bucket
from tools.performance_heatmap import NegativeBucketError
from tools.performance_heatmap import held_buckets
from tools.performance_heatmap import stock_buckets
from tools.performance_levers import Lever
from tools.performance_levers import LeverDoubleCountError
from tools.performance_levers import UnknownLeverFunctionError
from tools.performance_levers import lever_share
from tools.performance_levers import load_levers
from tools.performance_stats import FunctionKey
from tools.performance_stats import FunctionRow
from tools.performance_stats import LoadedStats
from tools.performance_stats import MissingProfileFunctionError

GUIDE = "/repo/benchmarks/held_figure5_raw_guide.py"
CASES = "/repo/benchmarks/held_reference_cases.py"
PATH = "/repo/benchmarks/held_figure5_path.py"
STOCK = "/repo/src/compas_cgal/stock.py"
ADMISSIBILITY = "/repo/src/compas_cgal/pipeline/admissibility.py"
TOTAL_S = 10.0

PLACEMENT: FunctionKey = (PATH, 493, HELD_PLACEMENT)
HYPOTHESES: FunctionKey = (GUIDE, 393, HELD_HYPOTHESES)
CLEARANCE: FunctionKey = (GUIDE, 373, HELD_CLEARANCE)
LOAD: FunctionKey = (CASES, 887, HELD_CASE_LOAD)


def _row(self_s: float, cumulative_s: float, callers: Mapping[FunctionKey, object] | None = None) -> FunctionRow:
    return (1, 1, self_s, cumulative_s, dict(callers or {}))


def _held(clearance_s: float = 5.0) -> LoadedStats:
    return LoadedStats(
        stats={
            LOAD: _row(0.0, 2.0),
            PLACEMENT: _row(0.5, 7.0),
            HYPOTHESES: _row(0.0, 6.0, {PLACEMENT: None}),
            CLEARANCE: _row(0.0, clearance_s, {HYPOTHESES: None}),
        },
        total_tt=TOTAL_S,
    )


def _stock(*, with_depletion: bool = True) -> LoadedStats:
    stats: dict[FunctionKey, FunctionRow] = {(ADMISSIBILITY, 10, "exceeds"): _row(6.0, 6.0)}
    if with_depletion:
        stats[(STOCK, 241, "subtract_line_sweep_local")] = _row(3.0, 3.0)
    return LoadedStats(stats=stats, total_tt=TOTAL_S)


def test_held_buckets_partition_the_profile() -> None:
    shares = held_buckets(_held()).shares
    assert shares == pytest.approx({Bucket.CLEARANCE: 50.0, Bucket.CASE_LOAD: 20.0, Bucket.HYPOTHESES: 20.0, Bucket.ORCHESTRATION: 10.0})
    assert sum(shares.values()) == pytest.approx(100.0)


def test_stock_buckets_partition_the_profile() -> None:
    shares = stock_buckets(_stock()).shares
    assert shares == pytest.approx({Bucket.PREDICATES: 60.0, Bucket.DEPLETION: 30.0, Bucket.ORCHESTRATION: 10.0})
    assert sum(shares.values()) == pytest.approx(100.0)


def test_held_bucket_function_renamed_raises() -> None:
    stats = _held()
    renamed = {(key if key != CLEARANCE else (GUIDE, 373, f"{HELD_CLEARANCE}_renamed")): row for key, row in stats.stats.items()}
    with pytest.raises(MissingProfileFunctionError, match=HELD_CLEARANCE):
        held_buckets(LoadedStats(stats=renamed, total_tt=TOTAL_S))


def test_held_bucket_functions_are_on_the_generation_path() -> None:
    # Static call graph of the generator: placement calls the guide-wide
    # hypotheses builder, which calls the native batch. A reroute that leaves a
    # bucket function off the path fails here, before any profile is read.
    assert HELD_PLACEMENT in figure5_path.build_figure5_approximate_path.__code__.co_names
    assert HELD_HYPOTHESES in getattr(figure5_path, HELD_PLACEMENT).__code__.co_names
    assert HELD_CLEARANCE in getattr(raw_guide, HELD_HYPOTHESES).__code__.co_names


def test_callpath_functions_are_on_the_generation_path() -> None:
    assert CALLPATH_HYPOTHESES in figure5_path.build_held_reference_path.__code__.co_names
    assert CALLPATH_BATCH in getattr(raw_guide, CALLPATH_HYPOTHESES).__code__.co_names


def test_callpath_rows_read_the_named_functions() -> None:
    # The table and the buckets read the same two functions; `_held` carries them.
    assert (CALLPATH_HYPOTHESES, CALLPATH_BATCH) == (HELD_HYPOTHESES, HELD_CLEARANCE)
    held = _held()
    stats = dict(held.stats)
    stats[(PATH, 600, "build_figure5_approximate_path")] = _row(0.0, 7.0)
    stats[(f"/repo/src/compas_cgal/{CLEARANCE_MODULE}", 1, CLEARANCE_QUERY)] = _row(4.0, 4.0)
    stats[(f"/repo/src/compas_cgal/{CLEARANCE_MODULE}", 2, CLEARANCE_DECODE)] = _row(0.5, 0.5)
    stats[(f"/lib/{FRACTION_MODULE}", 3, FRACTION_OPERATOR)] = _row(0.1, 0.1)
    stats[(GCD_MODULE, 0, GCD)] = _row(0.1, 0.1)
    rows = callpath_rows(LoadedStats(stats=stats, total_tt=TOTAL_S))
    assert [row.cumulative for row in rows[:3]] == ["7.0 s", "6.0 s", "5.0 s"]
    assert CALLPATH_HYPOTHESES in rows[1].path and CALLPATH_BATCH in rows[2].path


def test_stock_bucket_without_any_function_raises() -> None:
    with pytest.raises(MissingProfileFunctionError, match="stock.py"):
        stock_buckets(_stock(with_depletion=False))


def test_non_nested_subtrees_raise_instead_of_a_negative_bucket() -> None:
    with pytest.raises(NegativeBucketError):
        held_buckets(_held(clearance_s=9.0))


def _lever(functions: list[dict[str, object]]) -> Lever:
    return Lever.build({"key": "probe", "label": "probe", "status": "proposed", "functions": functions, "estimated": {}})


def test_lever_share_sums_named_functions() -> None:
    lever = _lever(
        [
            {"module": "benchmarks/held_reference_cases.py", "function": "load_held_reference_case", "time": "cumulative", "workloads": ["held_figure5"]},
            {"module": "benchmarks/held_figure5_path.py", "function": "build_held_reference_path", "time": "self", "workloads": ["held_figure5"]},
        ]
    )
    assert lever_share(lever, "held_figure5", _held()) == pytest.approx(25.0)
    assert lever_share(lever, "held_crossed_skis", _held()) is None


def test_lever_share_rejects_double_counting() -> None:
    lever = _lever(
        [
            {"module": "benchmarks/held_figure5_raw_guide.py", "function": HELD_HYPOTHESES, "time": "cumulative", "workloads": ["held_figure5"]},
            {"module": "benchmarks/held_figure5_raw_guide.py", "function": HELD_CLEARANCE, "time": "cumulative", "workloads": ["held_figure5"]},
        ]
    )
    with pytest.raises(LeverDoubleCountError):
        lever_share(lever, "held_figure5", _held())


def test_lever_function_missing_from_profile_raises() -> None:
    lever = _lever([{"module": "benchmarks/held_figure5_raw_guide.py", "function": "_boundary_arclength", "time": "cumulative", "workloads": ["held_figure5"]}])
    with pytest.raises(MissingProfileFunctionError, match="_boundary_arclength"):
        lever_share(lever, "held_figure5", _held())


def test_lever_function_missing_from_source_raises() -> None:
    with pytest.raises(UnknownLeverFunctionError, match="no_such_function"):
        _lever([{"module": "benchmarks/held_figure5_raw_guide.py", "function": "no_such_function", "time": "cumulative", "workloads": []}])


def test_committed_lever_mapping_names_existing_code(tmp_path: Path) -> None:
    levers = load_levers()
    assert levers, "the committed mapping is empty"
    assert all(lever.status.value in {"wired", "built", "proposed", "design"} for lever in levers)
    broken = json.loads((Path(__file__).resolve().parents[2] / "benchmarks/results/performance_levers.json").read_text())
    broken["levers"][0]["functions"][0]["function"] = "renamed_away"
    path = tmp_path / "levers.json"
    path.write_text(json.dumps(broken))
    with pytest.raises(UnknownLeverFunctionError, match="renamed_away"):
        load_levers(path)
