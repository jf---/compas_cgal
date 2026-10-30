"""The T(Omega*) certificate's CPU budget gate on the production Monstera path.

The production workload ``PRODUCTION_WORKLOAD`` (the Omega*-built Monstera path
with exact dOmega* connectors) is certified against T(Omega*) by
``tools.held_omega_coverage_results.measure_workload`` in ``BUDGET_RUNS`` fresh
interpreters (``python -m tools.coverage_certificate_budget_gate --worker``).
The verdict compares the minimum certificate process CPU with ``BUDGET_S``
(user rulings 2026-10-22: process CPU, not wall time; minimum of 5 fresh runs,
since interference only adds time). Gated measurement requires a quiet machine
before every run (``tools.held_generation_clock.require_quiet_machine``,
``MachineNotQuietError`` otherwise); indicative measurement runs at any load
and carries no verdict.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.performance_record import REPOSITORY
from tools.held_generation_clock import INDICATIVE_RUNS_PER_WORKLOAD
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import CpuSummary
from tools.held_generation_clock import LoadAverage
from tools.held_generation_clock import MachineNotQuietError
from tools.held_generation_clock import require_quiet_machine
from tools.held_omega_coverage_results import measure_workload as measure_omega_coverage

# Certificate process-CPU budget on Monstera, s (docs/held_coverage_certificate.md, "Cost model").
BUDGET_S = 5.0
# The workload the budget gates: the Omega*-built Monstera production path.
PRODUCTION_WORKLOAD = "figure8_monstera_boundary_qualified_corners"
# Fresh-process runs per budget measurement; the gate reads their minimum CPU (user ruling 2026-10-22).
BUDGET_RUNS = INDICATIVE_RUNS_PER_WORKLOAD
BUDGET_WORKER_MARKER = "CERTIFICATE_BUDGET_RUN "


class BudgetRunCountError(BenchmarkError):
    """A budget measurement does not hold exactly ``BUDGET_RUNS`` fresh-process runs."""


class BudgetWorkerError(BenchmarkError):
    """A fresh-process budget run exited non-zero or reported no single result line."""


class BudgetCoverageError(BenchmarkError):
    """A budget run's certificate was not COVERED: an uncovered path's time is not the production cost."""


class BudgetRunMismatchError(BenchmarkError):
    """The budget runs disagree on the certificate (witnesses or edges): they did not measure the same work."""


class BudgetWorkloadError(BenchmarkError):
    """A budget record names another workload than ``PRODUCTION_WORKLOAD``."""


class BudgetLimitError(BenchmarkError):
    """A budget record carries another budget than ``BUDGET_S``."""


class CoverageVerdict(str, Enum):
    """The certificate's coverage verdict of one run (``held_omega_coverage_results`` row ``verdict``)."""

    COVERED = "COVERED"
    UNCOVERED = "UNCOVERED"


class BudgetVerdict(Enum):
    """The certificate budget verdict of one production measurement."""

    WITHIN = "within"
    OVER = "OVER"
    INDICATIVE = "INDICATIVE"


@dataclass(frozen=True)
class BudgetRun:
    """One fresh-process certificate run of the production workload.

    Attributes:
        cpu_s: Certificate process CPU, s (the gated clock).
        wall_s: Certificate wall time, s (reported, never gated).
        load_1min: 1-minute load average read by the parent just before starting the worker.
        certificate_load_1min: 1-minute load read by the worker just before and just after the certificate.
        verdict: Coverage verdict of the certificate.
        witnesses: Uncovered witnesses (0 when covered).
        edges: Arrangement edges (``None`` when uncovered).
    """

    cpu_s: float
    wall_s: float
    load_1min: float
    certificate_load_1min: tuple[float, float]
    verdict: CoverageVerdict
    witnesses: int
    edges: int | None

    @property
    def peak_load_1min(self) -> float:
        """The highest of the run's three 1-minute load readings."""
        return max(self.load_1min, *self.certificate_load_1min)

    @classmethod
    def from_raw(cls, raw: dict[str, Any]) -> BudgetRun:
        """Rebuild from ``dataclasses.asdict`` output after a JSON round trip."""
        before, after = raw["certificate_load_1min"]
        return cls(
            cpu_s=raw["cpu_s"],
            wall_s=raw["wall_s"],
            load_1min=raw["load_1min"],
            certificate_load_1min=(before, after),
            verdict=CoverageVerdict(raw["verdict"]),
            witnesses=raw["witnesses"],
            edges=raw["edges"],
        )


def require_quiet_run(run: BudgetRun) -> None:
    """Refuse a gated run any of whose load readings is above the quiet limit.

    Raises:
        MachineNotQuietError: a reading (parent before the worker, worker before or after the certificate) is above ``QUIET_LOAD_AVERAGE_MAX``.
    """
    if run.peak_load_1min > QUIET_LOAD_AVERAGE_MAX:
        raise MachineNotQuietError(
            f"{PRODUCTION_WORKLOAD}: gated run measured at 1-minute load {run.load_1min:.2f} (worker start), "
            f"{run.certificate_load_1min[0]:.2f} / {run.certificate_load_1min[1]:.2f} (certificate before / after) > {QUIET_LOAD_AVERAGE_MAX}; "
            "budget verdicts count only from quiet-machine runs (user ruling D2, 2026-10-22); measure with --indicative"
        )


@dataclass(frozen=True)
class ProductionBudget:
    """The certificate budget measurement of ``PRODUCTION_WORKLOAD``.

    Attributes:
        workload: The production workload measured.
        budget_s: Process-CPU budget, s.
        indicative: Measured without the quiet-machine check; carries no verdict.
        runs: The fresh-process runs, in order.
    """

    workload: str
    budget_s: float
    indicative: bool
    runs: tuple[BudgetRun, ...]

    @classmethod
    def build(cls, runs: Sequence[BudgetRun], *, indicative: bool) -> ProductionBudget:
        """The measurement of ``PRODUCTION_WORKLOAD`` against ``BUDGET_S``, validated.

        Raises:
            BudgetRunCountError: not exactly ``BUDGET_RUNS`` runs.
            BudgetCoverageError: a run's certificate was not COVERED.
            BudgetRunMismatchError: the runs disagree on witnesses or edges.
            MachineNotQuietError: gated, and a run's load reading was above ``QUIET_LOAD_AVERAGE_MAX``.
        """
        if len(runs) != BUDGET_RUNS:
            raise BudgetRunCountError(f"{PRODUCTION_WORKLOAD}: {len(runs)} runs, the budget gate takes exactly {BUDGET_RUNS} fresh-process runs")
        uncovered = [i + 1 for i, run in enumerate(runs) if run.verdict is not CoverageVerdict.COVERED]
        if uncovered:
            raise BudgetCoverageError(f"{PRODUCTION_WORKLOAD}: runs {uncovered} are not COVERED; the budget times the covered production certificate only")
        shapes = {(run.witnesses, run.edges) for run in runs}
        if len(shapes) != 1:
            raise BudgetRunMismatchError(f"{PRODUCTION_WORKLOAD}: runs disagree on (witnesses, edges): {sorted(shapes, key=str)}")
        if not indicative:
            for run in runs:
                require_quiet_run(run)
        return cls(PRODUCTION_WORKLOAD, BUDGET_S, indicative, tuple(runs))

    @classmethod
    def from_raw(cls, raw: dict[str, Any]) -> ProductionBudget:
        """Rebuild from ``dataclasses.asdict`` output, re-validating through :meth:`build`.

        Raises:
            BudgetWorkloadError: the record names another workload.
            BudgetLimitError: the record carries another budget.
            BudgetRunCountError, BudgetCoverageError, BudgetRunMismatchError, MachineNotQuietError: as :meth:`build`.
        """
        if raw["workload"] != PRODUCTION_WORKLOAD:
            raise BudgetWorkloadError(f"budget record of {raw['workload']!r}, the gate measures {PRODUCTION_WORKLOAD!r}; re-measure")
        if raw["budget_s"] != BUDGET_S:
            raise BudgetLimitError(f"budget record against {raw['budget_s']!r} s, the declared budget is {BUDGET_S!r} s; re-measure")
        return cls.build(tuple(BudgetRun.from_raw(run) for run in raw["runs"]), indicative=raw["indicative"])

    @property
    def cpu(self) -> CpuSummary:
        """Process CPU over the runs; the gate reads the minimum."""
        return CpuSummary.build([run.cpu_s for run in self.runs])

    @property
    def verdict(self) -> BudgetVerdict:
        """``INDICATIVE`` when measured under load, else the minimum CPU against the budget."""
        if self.indicative:
            return BudgetVerdict.INDICATIVE
        return BudgetVerdict.WITHIN if self.cpu.gated <= self.budget_s else BudgetVerdict.OVER


def budget_worker_command() -> list[str]:
    """The fresh-interpreter command that certifies the production workload once and reports its row."""
    return [sys.executable, "-m", "tools.coverage_certificate_budget_gate", "--worker"]


def measure_budget_run_in_fresh_process(load_1min: float) -> BudgetRun:
    """Certify ``PRODUCTION_WORKLOAD`` once in a fresh interpreter (same executable and environment).

    Args:
        load_1min: The 1-minute load read just before the run, recorded with it.

    Raises:
        BudgetWorkerError: the worker exited non-zero or reported other than one result line.
    """
    completed = subprocess.run(budget_worker_command(), cwd=REPOSITORY, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise BudgetWorkerError(f"budget worker exited {completed.returncode}; stderr tail:\n{completed.stderr[-4000:]}")
    lines = [line for line in completed.stdout.splitlines() if line.startswith(BUDGET_WORKER_MARKER)]
    if len(lines) != 1:
        raise BudgetWorkerError(f"expected one {BUDGET_WORKER_MARKER.strip()} line, found {len(lines)}; stdout tail:\n{completed.stdout[-2000:]}")
    row = json.loads(lines[0][len(BUDGET_WORKER_MARKER) :])
    before, after = row["load_1min_before_after"]
    return BudgetRun(
        cpu_s=row["certificate_cpu_seconds"],
        wall_s=row["certificate_wall_seconds"],
        load_1min=load_1min,
        certificate_load_1min=(before, after),
        verdict=CoverageVerdict(row["verdict"]),
        witnesses=row["witnesses"],
        edges=row["edges"],
    )


def measure_budget_runs(
    *,
    indicative: bool,
    load_reader: Callable[[], LoadAverage] = os.getloadavg,
    measure: Callable[[float], BudgetRun] = measure_budget_run_in_fresh_process,
) -> ProductionBudget:
    """Measure ``BUDGET_RUNS`` serial fresh-process certificate runs of the production workload.

    Gated (``indicative=False``): each run first requires a quiet machine, and the
    worker's own readings around the certificate must be quiet too. Indicative:
    each run reads the load and runs at any load; the result carries no verdict.

    Raises:
        MachineNotQuietError: gated, and the load was above the quiet limit before a run or
            around its certificate; no further run is measured.
    """
    runs = []
    for _ in range(BUDGET_RUNS):
        load = load_reader() if indicative else require_quiet_machine(load_reader)
        run = measure(load[0])
        if not indicative:
            require_quiet_run(run)
        runs.append(run)
    return ProductionBudget.build(runs, indicative=indicative)


def render_budget(budget: ProductionBudget) -> list[str]:
    """The budget gate's verdict line and run table (production workload only)."""
    cpu = budget.cpu
    loads = [run.peak_load_1min for run in budget.runs]
    summary = f"process CPU min {cpu.minimum:.2f} / median {cpu.median:.2f} / max {cpu.maximum:.2f} s over {len(budget.runs)} fresh-process runs"
    if budget.verdict is BudgetVerdict.INDICATIVE:
        verdict = f"INDICATIVE (measured without the quiet-machine check, 1-minute load {min(loads):.2f}-{max(loads):.2f}; no budget verdict): {summary}."
    else:
        verdict = f"min {cpu.gated:.2f} s, {budget.verdict.value} the {budget.budget_s:g} s budget ({summary}, quiet machine)."
    lines = [
        f"**Certificate budget gate**, production workload `{budget.workload}` (Omega*-built path, exact dOmega* connectors) against T(Omega*), gated on process CPU: {verdict}",
        "",
        "| run | coverage | witnesses | edges | CPU s | wall s (not gated) | 1-min load: worker start, certificate before / after |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    lines += [
        f"| {i + 1} | {run.verdict.value} | {run.witnesses} | {'—' if run.edges is None else run.edges} | {run.cpu_s:.2f} | {run.wall_s:.2f} "
        f"| {run.load_1min:.2f}, {run.certificate_load_1min[0]:.2f} / {run.certificate_load_1min[1]:.2f} |"
        for i, run in enumerate(budget.runs)
    ]
    return lines


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--worker", action="store_true", required=True, help="certify the production workload once in this process and print its row")
    parser.parse_args(argv)
    print(BUDGET_WORKER_MARKER + json.dumps(measure_omega_coverage(PRODUCTION_WORKLOAD)), flush=True)


if __name__ == "__main__":
    main()
