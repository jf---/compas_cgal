"""The certificate budget gate: process CPU, minimum of 5 fresh runs, quiet machine only, covered production path.

User rulings (2026-10-22): the certificate budget gates process CPU, not wall
time; the gated value is the minimum over 5 fresh-process runs; a verdict
counts only from a quiet machine (``MachineNotQuietError`` above the quiet
load), otherwise the numbers are labelled indicative and carry no verdict.
The gate measures the production workload only, and only a covered,
self-consistent certificate is timed against the budget.
"""

from __future__ import annotations

import dataclasses
import json
import subprocess
import sys
from collections.abc import Iterator
from typing import Any

import pytest

from benchmarks.held_omega_toolpath import OMEGA_WORKLOADS
from tools.coverage_certificate_budget_gate import BUDGET_RUNS
from tools.coverage_certificate_budget_gate import BUDGET_S
from tools.coverage_certificate_budget_gate import BUDGET_WORKER_MARKER
from tools.coverage_certificate_budget_gate import PRODUCTION_WORKLOAD
from tools.coverage_certificate_budget_gate import BudgetCoverageError
from tools.coverage_certificate_budget_gate import BudgetLimitError
from tools.coverage_certificate_budget_gate import BudgetRun
from tools.coverage_certificate_budget_gate import BudgetRunCountError
from tools.coverage_certificate_budget_gate import BudgetRunMismatchError
from tools.coverage_certificate_budget_gate import BudgetVerdict
from tools.coverage_certificate_budget_gate import BudgetWorkerError
from tools.coverage_certificate_budget_gate import BudgetWorkloadError
from tools.coverage_certificate_budget_gate import CoverageVerdict
from tools.coverage_certificate_budget_gate import ProductionBudget
from tools.coverage_certificate_budget_gate import measure_budget_run_in_fresh_process
from tools.coverage_certificate_budget_gate import measure_budget_runs
from tools.coverage_certificate_budget_gate import render_budget
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import MachineNotQuietError

QUIET = 1.0
LOADED = 450.0
# A wall time far above the budget: a wall-time gate would say OVER.
SLOW_WALL_S = 9.0
EDGES = 1_210_313


def _run(cpu_s: float, wall_s: float = SLOW_WALL_S, load: float = QUIET, certificate_load: float = QUIET) -> BudgetRun:
    return BudgetRun(
        cpu_s=cpu_s,
        wall_s=wall_s,
        load_1min=load,
        certificate_load_1min=(certificate_load, certificate_load),
        verdict=CoverageVerdict.COVERED,
        witnesses=0,
        edges=EDGES,
    )


def _quiet_budget(cpus: tuple[float, ...]) -> ProductionBudget:
    return ProductionBudget.build(tuple(_run(c) for c in cpus), indicative=False)


def test_the_production_workload_is_the_omega_star_monstera_path() -> None:
    assert PRODUCTION_WORKLOAD == "figure8_monstera_boundary_qualified_corners"
    assert PRODUCTION_WORKLOAD in OMEGA_WORKLOADS


def test_the_verdict_gates_the_minimum_process_cpu_not_wall_time() -> None:
    within = _quiet_budget((5.4, 5.3, BUDGET_S - 0.1, 5.2, 5.6))
    assert within.verdict is BudgetVerdict.WITHIN
    over = ProductionBudget.build(tuple(_run(BUDGET_S + 0.1, wall_s=1.0) for _ in range(BUDGET_RUNS)), indicative=False)
    assert over.verdict is BudgetVerdict.OVER


def test_the_gate_takes_exactly_the_ruled_number_of_runs() -> None:
    assert BUDGET_RUNS == 5
    with pytest.raises(BudgetRunCountError):
        _quiet_budget((4.0,) * (BUDGET_RUNS - 1))


def test_a_gated_record_with_a_loaded_run_is_refused() -> None:
    runs = tuple(_run(4.0, load=QUIET) for _ in range(BUDGET_RUNS - 1)) + (_run(4.0, load=QUIET_LOAD_AVERAGE_MAX + 1),)
    with pytest.raises(MachineNotQuietError):
        ProductionBudget.build(runs, indicative=False)


def test_a_gated_record_with_a_loaded_certificate_is_refused() -> None:
    runs = tuple(_run(4.0) for _ in range(BUDGET_RUNS - 1)) + (_run(4.0, certificate_load=QUIET_LOAD_AVERAGE_MAX + 1),)
    with pytest.raises(MachineNotQuietError):
        ProductionBudget.build(runs, indicative=False)


@pytest.mark.parametrize("indicative", [False, True])
def test_an_uncovered_run_is_refused(indicative: bool) -> None:
    uncovered = dataclasses.replace(_run(4.0), verdict=CoverageVerdict.UNCOVERED, witnesses=3, edges=None)
    with pytest.raises(BudgetCoverageError):
        ProductionBudget.build(tuple(_run(4.0) for _ in range(BUDGET_RUNS - 1)) + (uncovered,), indicative=indicative)


def test_runs_that_disagree_on_the_certificate_are_refused() -> None:
    other = dataclasses.replace(_run(4.0), edges=EDGES + 1)
    with pytest.raises(BudgetRunMismatchError):
        ProductionBudget.build(tuple(_run(4.0) for _ in range(BUDGET_RUNS - 1)) + (other,), indicative=False)


def test_an_indicative_measurement_carries_no_verdict() -> None:
    budget = ProductionBudget.build(tuple(_run(4.0, load=LOADED, certificate_load=LOADED) for _ in range(BUDGET_RUNS)), indicative=True)
    assert budget.verdict is BudgetVerdict.INDICATIVE
    text = "\n".join(render_budget(budget))
    assert "INDICATIVE" in text
    assert "within the" not in text
    assert "OVER the" not in text


def test_the_verdict_line_names_the_production_workload_and_process_cpu() -> None:
    text = "\n".join(render_budget(_quiet_budget((4.9, 5.2, 5.3, 5.1, 5.0))))
    verdict_line = next(line for line in text.splitlines() if f"within the {BUDGET_S:g} s budget" in line)
    assert PRODUCTION_WORKLOAD in verdict_line
    assert "process CPU" in verdict_line


def _loads(values: list[float]) -> Iterator[tuple[float, float, float]]:
    for value in values:
        yield (value, value, value)


def test_gated_measurement_refuses_a_loaded_machine_before_measuring() -> None:
    loads = _loads([QUIET, QUIET, LOADED, QUIET, QUIET])
    measured: list[float] = []

    def measure(load: float) -> BudgetRun:
        measured.append(load)
        return _run(4.0, load=load)

    with pytest.raises(MachineNotQuietError):
        measure_budget_runs(indicative=False, load_reader=lambda: next(loads), measure=measure)
    assert measured == [QUIET, QUIET]


def test_gated_measurement_stops_at_a_run_whose_certificate_ran_under_load() -> None:
    loads = _loads([QUIET] * BUDGET_RUNS)
    measured: list[float] = []

    def measure(load: float) -> BudgetRun:
        measured.append(load)
        return _run(4.0, load=load, certificate_load=LOADED if len(measured) == 2 else QUIET)

    with pytest.raises(MachineNotQuietError):
        measure_budget_runs(indicative=False, load_reader=lambda: next(loads), measure=measure)
    assert len(measured) == 2


def test_indicative_measurement_runs_at_any_load_and_records_it() -> None:
    loads = _loads([LOADED] * BUDGET_RUNS)
    budget = measure_budget_runs(indicative=True, load_reader=lambda: next(loads), measure=lambda load: _run(4.0, load=load, certificate_load=LOADED))
    assert budget.indicative
    assert [run.load_1min for run in budget.runs] == [LOADED] * BUDGET_RUNS


def test_from_raw_round_trips_and_re_validates() -> None:
    budget = _quiet_budget((4.9, 5.2, 5.3, 5.1, 5.0))
    raw = json.loads(json.dumps(dataclasses.asdict(budget)))
    assert ProductionBudget.from_raw(raw) == budget
    raw["runs"][0]["load_1min"] = LOADED
    with pytest.raises(MachineNotQuietError):
        ProductionBudget.from_raw(raw)


def test_from_raw_refuses_another_workload_or_budget() -> None:
    raw = json.loads(json.dumps(dataclasses.asdict(_quiet_budget((4.9, 5.2, 5.3, 5.1, 5.0)))))
    with pytest.raises(BudgetWorkloadError):
        ProductionBudget.from_raw({**raw, "workload": "figure8_upper_boundary"})
    with pytest.raises(BudgetLimitError):
        ProductionBudget.from_raw({**raw, "budget_s": 2 * BUDGET_S})


# A worker row with distinct CPU and wall times and certificate-time loads.
WORKER_CPU_S = 4.25
WORKER_WALL_S = 7.75
WORKER_LOADS = [0.5, 0.75]


def _worker_row() -> dict[str, Any]:
    return {
        "certificate_cpu_seconds": WORKER_CPU_S,
        "certificate_wall_seconds": WORKER_WALL_S,
        "load_1min_before_after": WORKER_LOADS,
        "verdict": "COVERED",
        "witnesses": 0,
        "edges": EDGES,
    }


def _fake_run(monkeypatch: pytest.MonkeyPatch, returncode: int, stdout: str) -> list[list[str]]:
    calls: list[list[str]] = []

    def fake(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, returncode, stdout=stdout, stderr="worker stderr")

    monkeypatch.setattr(subprocess, "run", fake)  # the gate calls subprocess.run through the module
    return calls


def test_a_fresh_worker_run_reads_the_childs_cpu_not_its_wall(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _fake_run(monkeypatch, 0, "noise\n" + BUDGET_WORKER_MARKER + json.dumps(_worker_row()) + "\n")
    run = measure_budget_run_in_fresh_process(QUIET)
    assert run.cpu_s == WORKER_CPU_S
    assert run.wall_s == WORKER_WALL_S
    assert run.load_1min == QUIET
    assert run.certificate_load_1min == tuple(WORKER_LOADS)
    assert run.verdict is CoverageVerdict.COVERED
    assert calls == [[sys.executable, "-m", "tools.coverage_certificate_budget_gate", "--worker"]]


@pytest.mark.parametrize(
    ("returncode", "markers"),
    [(1, 1), (0, 0), (0, 2)],
    ids=["non-zero exit", "no marker", "two markers"],
)
def test_a_bad_worker_raises(monkeypatch: pytest.MonkeyPatch, returncode: int, markers: int) -> None:
    _fake_run(monkeypatch, returncode, "".join(BUDGET_WORKER_MARKER + json.dumps(_worker_row()) + "\n" for _ in range(markers)))
    with pytest.raises(BudgetWorkerError):
        measure_budget_run_in_fresh_process(QUIET)
