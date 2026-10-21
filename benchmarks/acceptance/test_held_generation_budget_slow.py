"""Held generation-time budgets: median of three serial fresh-process generations per workload.

Not collected by pytest (benchmarks/acceptance); run with `pixi run
held-generation-budget`, about five minutes while Monstera is over budget.
Generation is ``GeneratedToolpath.generation_seconds``: guide, inward offset,
initial placement, engagement refinement and connectors; case load excluded.
Each run is a fresh interpreter (``tools.held_path_digest --worker``), run
one after the other, never in parallel, so no run shares a warm cache or a
core with another. The three runs must also agree on every digest.
"""

from __future__ import annotations

import os
import statistics

import pytest

from tools.held_path_digest import WorkloadMeasurement
from tools.held_path_digest import measure_in_fresh_process

# Serial fresh-process generations per workload; odd, so the median is a measured run.
RUNS_PER_WORKLOAD = 3
# Generation budgets, seconds, quiet single core. Source: the "lever budget
# (initial / generation)" column of the hypotheses-and-selection lever plan,
# .superpowers/sdd/2026-10-13-product-pipeline/lever-hypotheses-selection-plan.md,
# "Time budget for this lever" (2026-10-21), built from the measured components
# of the contact-projection lever after Python steps P1-P7 and C1.
FIGURE5_GENERATION_BUDGET_S = 1.7
CROSSED_SKIS_GENERATION_BUDGET_S = 2.2
MONSTERA_GENERATION_BUDGET_S = 9.0
GENERATION_BUDGETS_S: dict[str, float] = {
    "figure5_refined": FIGURE5_GENERATION_BUDGET_S,
    "figure8_crossed_skis_boundary": CROSSED_SKIS_GENERATION_BUDGET_S,
    "figure8_monstera_boundary_qualified_corners": MONSTERA_GENERATION_BUDGET_S,
}


def _describe(runs: list[WorkloadMeasurement], loads: list[tuple[float, float, float]]) -> str:
    lines = []
    for index, (run, load) in enumerate(zip(runs, loads, strict=True)):
        stages = " ".join(f"{stage} {seconds:.2f}" for stage, seconds in run.stage_seconds.items())
        lines.append(f"  run {index + 1}: {run.generation_seconds:.2f} s [{stages}]  load avg before {load[0]:.1f} / {load[1]:.1f} / {load[2]:.1f}")
    return "\n".join(lines)


@pytest.mark.parametrize("workload", list(GENERATION_BUDGETS_S))
def test_held_generation_is_within_budget(workload: str) -> None:
    budget = GENERATION_BUDGETS_S[workload]
    runs: list[WorkloadMeasurement] = []
    loads: list[tuple[float, float, float]] = []
    for _ in range(RUNS_PER_WORKLOAD):
        loads.append(os.getloadavg())
        runs.append(measure_in_fresh_process(workload))
    detail = _describe(runs, loads)
    print(f"\n{workload}\n{detail}")
    assert len({run.digest for run in runs}) == 1, f"{workload}: fresh-process runs disagree on the generated path:\n" + "\n".join(repr(run.digest) for run in runs)
    median = statistics.median(run.generation_seconds for run in runs)
    assert median <= budget, f"{workload}: median generation {median:.2f} s > budget {budget} s over {RUNS_PER_WORKLOAD} serial fresh-process runs:\n{detail}"
