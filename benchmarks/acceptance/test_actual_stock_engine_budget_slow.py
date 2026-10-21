"""Actual-stock engine budget: full replay of every sweep plus one 45-station decision per circle, one core.

Not collected by pytest (benchmarks/acceptance); run with `pixi run
actual-stock-engine-budget`, about two minutes while over budget. The
budgets are user-approved (docs/held_actual_stock_engine.md section 1,
2026-10-21) and live in ``tools.actual_stock_profile.ENGINE_BUDGETS``:
Monstera ≤ 5 s (≤ 0.25 ms per sweep, ≤ 0.6 ms per query), Figure 5 ≤ 2.5 s.

Each run is a fresh interpreter (``tools.actual_stock_profile
--budget-worker``), serial, timed in process CPU seconds and stopped as soon
as its CPU passes the budget: past that point the verdict cannot change, so
an over-budget engine fails in budget time instead of the hour today's full
replay takes. The median run over ``RUNS_PER_WORKLOAD`` decides; every run
and its load average is printed.
"""

from __future__ import annotations

import statistics

import pytest

from tools.actual_stock_profile import ENGINE_BUDGETS
from tools.actual_stock_profile import BudgetRun
from tools.actual_stock_profile import measure_budget_in_fresh_process

# Serial fresh-process runs per workload; odd, so the median is a measured run.
RUNS_PER_WORKLOAD = 3


@pytest.mark.parametrize("workload", list(ENGINE_BUDGETS))
def test_actual_stock_engine_is_within_budget(workload: str) -> None:
    runs: list[BudgetRun] = [measure_budget_in_fresh_process(workload) for _ in range(RUNS_PER_WORKLOAD)]
    detail = "\n".join(f"  run {i + 1}: {run.describe()}" for i, run in enumerate(runs))
    print(f"\n{workload}\n{detail}")
    median = sorted(runs, key=lambda run: run.total_cpu_s)[RUNS_PER_WORKLOAD // 2]
    assert median.total_cpu_s == statistics.median(run.total_cpu_s for run in runs)
    failures = median.failures()
    assert not failures, f"{workload}: median run over the engine budget: " + "; ".join(failures) + f"\n{detail}"
