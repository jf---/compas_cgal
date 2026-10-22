"""Slow tier of the exact last-16 certificate (stage C1): the differential against the full exact replay, and the budgets.

Not collected by pytest (benchmarks/acceptance). Two tests:

- **Differential** (Figure 5 and Monstera, any load; correctness, no budget):
  every circle is certified over 8 processes, every ``AUDIT_EVERY``-th circle
  also has both window stocks checked ``exactly_equals`` against their
  local-surgery references and the bound (full window inside last-16,
  verdict monotone), and the full sequential exact replay decides every
  circle. The final verdicts must equal the replay's per circle (refused iff
  exceeded on the replay), on every circle; ``summarize`` raises first on a
  bound break. The replay is single-threaded: about 9 CPU minutes on Figure 5
  and 39 on Monstera (C1, load 10-50).
- **Budgets** (Figure 5 and Monstera, quiet machine only): budget verdicts
  count only from quiet-machine runs (user ruling D2, 2026-10-22), so a
  1-minute load average above ``QUIET_LOAD_AVERAGE_MAX`` refuses in fixture
  setup with ``MachineNotQuietError`` (pytest reports ERROR: never a skip,
  never a budget verdict). Then the per-circle CPU sum (AS2 plan section 3:
  Figure 5 <= 100 s, Monstera <= 200 s; the sum of the per-circle process
  CPU of the certificates) and, on Monstera, the wall time on 8 processes
  (<= 30 s) are asserted on one run. A single run's CPU and wall are each
  >= their minimum over N runs, so a PASS settles the min-of-5 verdict; a
  FAIL is inconclusive (single run; a min-of-5 quiet re-run decides).
"""

from __future__ import annotations

import os
import time

import pytest

from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_actual_stock_last16_batch import certify_batch
from benchmarks.held_actual_stock_last16_batch import render_markdown
from benchmarks.held_actual_stock_last16_batch import summarize
from benchmarks.held_actual_stock_last16_reference import replay_exceeded
from tools.actual_stock_profile import FIGURE5
from tools.actual_stock_profile import MONSTERA
from tools.actual_stock_profile import ReplayWorkload
from tools.held_generation_clock import require_quiet_machine

PROCESSES = 8
# Stage C1 budgets, seconds (AS2 plan section 3).
CPU_BUDGET_S = {FIGURE5: 100.0, MONSTERA: 200.0}
MONSTERA_WALL_BUDGET_S = 30.0
# Appended to every budget failure: one run bounds the min-of-5 from above only.
SINGLE_RUN_NOTE = "single run; a min-of-5 quiet re-run decides"
# Every this many circles also gets the exactly_equals stock audit and the
# bound check: 41 circles on Figure 5, 89 on Monstera.
AUDIT_EVERY = 50


@pytest.fixture
def quiet_machine() -> None:
    """Refuse a loaded machine before a budget measurement (user ruling D2, 2026-10-22)."""
    require_quiet_machine()


@pytest.mark.parametrize("name", [FIGURE5, MONSTERA])
def test_last16_verdicts_equal_full_replay(name: str) -> None:
    workload = ReplayWorkload.build(name)
    load = os.getloadavg()
    wall = time.perf_counter()
    certificates = certify_batch(workload, range(len(workload.circles)), PROCESSES, audit_every=AUDIT_EVERY)
    wall_s = time.perf_counter() - wall
    reference = replay_exceeded(workload, progress=True)
    summary = summarize(workload, certificates, PROCESSES, float("nan"), wall_s, load, reference)
    print("\n" + render_markdown(summary))
    refused = {c.circle: c.verdict is Verdict.REFUSED for c in certificates}
    assert refused.keys() == reference.keys(), f"{name}: certified {len(refused)} circles, replayed {len(reference)}"
    mismatched = sorted(i for i in reference if refused[i] != reference[i])
    assert not mismatched, f"{name}: verdicts differ from the full replay on circles {mismatched[:20]}"


@pytest.mark.parametrize("name", [FIGURE5, MONSTERA])
@pytest.mark.usefixtures("quiet_machine")
def test_last16_budget(name: str) -> None:
    workload = ReplayWorkload.build(name)
    load = os.getloadavg()
    wall = time.perf_counter()
    certificates = certify_batch(workload, range(len(workload.circles)), PROCESSES)
    wall_s = time.perf_counter() - wall
    summary = summarize(workload, certificates, PROCESSES, float("nan"), wall_s, load, None)
    print("\n" + render_markdown(summary))
    assert summary.circle_cpu_s <= CPU_BUDGET_S[name], f"{name}: {summary.circle_cpu_s:.1f} s CPU > {CPU_BUDGET_S[name]} s ({SINGLE_RUN_NOTE})"
    if name == MONSTERA:
        assert wall_s <= MONSTERA_WALL_BUDGET_S, f"{name}: {wall_s:.1f} s wall > {MONSTERA_WALL_BUDGET_S} s on {PROCESSES} processes, load {load} ({SINGLE_RUN_NOTE})"
