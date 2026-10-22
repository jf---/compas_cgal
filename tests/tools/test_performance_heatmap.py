"""The performance figures' arithmetic: buckets partition the profile, lever cells read named code, renames fail loud."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import pytest

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

PLACEMENT: FunctionKey = (PATH, 493, "build_held_reference_path")
HYPOTHESES: FunctionKey = (GUIDE, 393, "build_distance_admissible_hypotheses")
CLEARANCE: FunctionKey = (GUIDE, 373, "_native_station_boundary")
LOAD: FunctionKey = (CASES, 887, "load_held_reference_case")


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
    renamed = {(key if key != CLEARANCE else (GUIDE, 373, "_native_station_boundary_lookup")): row for key, row in stats.stats.items()}
    with pytest.raises(MissingProfileFunctionError, match="_native_station_boundary"):
        held_buckets(LoadedStats(stats=renamed, total_tt=TOTAL_S))


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
            {"module": "benchmarks/held_figure5_raw_guide.py", "function": "build_distance_admissible_hypotheses", "time": "cumulative", "workloads": ["held_figure5"]},
            {"module": "benchmarks/held_figure5_raw_guide.py", "function": "_native_station_boundary", "time": "cumulative", "workloads": ["held_figure5"]},
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
