"""Held Ω* generation-time budgets: minimum process CPU time of five serial fresh-process generations per workload.

The exact product path (``benchmarks.held_omega_toolpath.generate_omega``) is
held to the polygon path's budgets and gate rules
(``test_held_generation_budget_slow.py``). It is expected RED until the levers
of the Ω* profile land (.superpowers/sdd/2026-10-13-product-pipeline/
op-omega-profile-report.md: 1.75-2.9x over on 2026-10-22); the test documents
the target. Generation time is a gate: slow is a red finding, never a baseline.

Not collected by pytest (benchmarks/acceptance); run with ``pytest
benchmarks/acceptance/test_held_omega_generation_budget_slow.py -q -s`` through
the editable build, about five minutes while over budget. The gate is process
CPU time (``time.process_time``) of the generation, user ruling 2026-10-22: the
stages of ``OmegaGeneratedToolpath.stage_seconds`` (guide, Ω* domain, initial
placement, refinement, NOTHING-circle removal, connector repair) plus the
collector's exit collection; case load excluded. What is counted and why:
``tools/held_generation_clock.py``. The gated value is the minimum CPU over the
runs (second user ruling 2026-10-22): interference from other sessions only
adds time. Median and spread of CPU, and wall time with the load average, are
printed beside it and not gated.

Budget verdicts count only from quiet-machine runs (user ruling D2,
2026-10-22): before every run the 1-minute load average is read, and above
``QUIET_LOAD_AVERAGE_MAX`` (2) the test raises ``MachineNotQuietError`` instead
of measuring. Not a skip, and no verdict: a loaded machine at the start is
refused in fixture setup (pytest reports ERROR); a load that rises between runs
is refused in the test body (pytest reports FAILED, by the same named error,
never by the budget assertion). Under load, measure indicatively with
``python -m tools.held_generation_clock --indicative --domain omega
<workload>...``, which prints the same numbers and never compares them with a
budget.

Each run is a fresh interpreter (``tools.held_generation_clock --domain
omega``), run one after the other, never in parallel, so no run shares a warm
cache with another. The runs must also agree on every digest.
"""

from __future__ import annotations

import statistics

import pytest

from benchmarks.held_omega_toolpath import PathDomain
from tools.held_generation_clock import ClockedGeneration
from tools.held_generation_clock import CpuSummary
from tools.held_generation_clock import MeasuredRun
from tools.held_generation_clock import measure_clocked_in_fresh_process
from tools.held_generation_clock import measure_gated_runs
from tools.held_generation_clock import require_quiet_machine

# Serial fresh-process generations per workload; the gate reads their minimum
# CPU (user ruling 2026-10-22); odd, so the reported median is a measured run.
# The polygon gate's count.
RUNS_PER_WORKLOAD = 5
# Generation budgets, CPU seconds, one core: the polygon path's budgets
# (test_held_generation_budget_slow.py), whose source is the "lever budget
# (initial / generation)" column of the hypotheses-and-selection lever plan,
# .superpowers/sdd/2026-10-13-product-pipeline/lever-hypotheses-selection-plan.md
# (2026-10-21). Ω* is the exact replacement of that path and is held to the same
# target (op-omega-profile-report.md, "Ω* budget test: recommendation").
FIGURE5_OMEGA_GENERATION_BUDGET_S = 1.7
CROSSED_SKIS_OMEGA_GENERATION_BUDGET_S = 2.2
MONSTERA_OMEGA_GENERATION_BUDGET_S = 9.0
OMEGA_GENERATION_BUDGETS_S: dict[str, float] = {
    "figure5_refined": FIGURE5_OMEGA_GENERATION_BUDGET_S,
    "figure8_crossed_skis_boundary": CROSSED_SKIS_OMEGA_GENERATION_BUDGET_S,
    "figure8_monstera_boundary_qualified_corners": MONSTERA_OMEGA_GENERATION_BUDGET_S,
}


def _measure_omega(workload: str) -> ClockedGeneration:
    return measure_clocked_in_fresh_process(workload, PathDomain.OMEGA)


def _describe(measured: list[MeasuredRun]) -> str:
    lines = []
    for index, entry in enumerate(measured):
        run, load = entry.generation, entry.load_before
        clock = run.clock
        stages = " ".join(f"{stage} {seconds:.2f}" for stage, seconds in clock.stage_cpu_seconds.items())
        lines.append(
            f"  run {index + 1}: CPU {clock.cpu_seconds:.2f} s [{stages} exit-collect {clock.exit_cpu_seconds:.2f}]"
            f"  wall {clock.wall_seconds:.2f} s  worker CPU {run.worker_cpu_seconds:.2f} s"
            f"  load avg before {load[0]:.1f} / {load[1]:.1f} / {load[2]:.1f}"
        )
    return "\n".join(lines)


@pytest.fixture
def quiet_machine() -> None:
    """Refuse a loaded machine before any measurement (user ruling D2, 2026-10-22)."""
    require_quiet_machine()


@pytest.mark.parametrize("workload", list(OMEGA_GENERATION_BUDGETS_S))
@pytest.mark.usefixtures("quiet_machine")
def test_held_omega_generation_is_within_budget(workload: str) -> None:
    budget = OMEGA_GENERATION_BUDGETS_S[workload]
    measured = measure_gated_runs(workload, RUNS_PER_WORKLOAD, measure=_measure_omega)
    runs = [entry.generation for entry in measured]
    assert all(run.domain is PathDomain.OMEGA for run in runs), f"{workload}: a run was not generated on Ω*."
    detail = _describe(measured)
    median_wall = statistics.median(run.clock.wall_seconds for run in runs)
    cpu = CpuSummary.build([run.clock.cpu_seconds for run in runs])
    path = runs[0].digest
    print(
        f"\n{workload} (omega): min CPU {cpu.minimum:.2f} s (gated, budget {budget} s); median {cpu.median:.2f}, max {cpu.maximum:.2f}, "
        f"spread {cpu.spread:.2f} s; median wall {median_wall:.2f} s (not gated); "
        f"circles {path.circle_count} {path.circles[:12]}, connectors {path.connector_count} {path.connectors[:12]}\n{detail}"
    )
    assert len({run.digest for run in runs}) == 1, f"{workload}: fresh-process Ω* runs disagree on the generated path:\n" + "\n".join(repr(run.digest) for run in runs)
    assert cpu.gated <= budget, f"{workload}: minimum Ω* generation CPU {cpu.gated:.2f} s > budget {budget} s over {RUNS_PER_WORKLOAD} serial fresh-process runs:\n{detail}"
