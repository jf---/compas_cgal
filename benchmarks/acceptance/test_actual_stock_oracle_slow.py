"""Slow tier of the actual-stock differential oracle (unit U0 of the actual-stock refinement design).

Not collected by pytest (benchmarks/acceptance); run by path. Correctness
only, no budget, any load:

- **D4** (U0 exit gate): per workload, circles accepted on the actual stock
  S_t but refused on the circle view C_t (part polygon minus every prior
  circle's swept annulus, connectors left out) must be 0. Figure 5 and skis
  every circle, Monstera every ``MONSTERA_EVERY``-th. Every
  ``AUDIT_EVERY``-th circle also gets the actual verdict, so the inclusion
  C_t ⊇ S_t is checked as a bound (``BoundBreakError``) beyond the refusals.
- **D2** (the property the U4 exit gate will hold production to; a harness
  only today, since no generation stage decides it): every emitted circle
  t >= 1 is certified on the actual stock and circle 0 is refused. Skis
  every circle.

Design: docs/superpowers/specs/2026-11-04-actual-stock-refinement-design.md,
section 3 (D2, D4) and section 5 (U0).
"""

from __future__ import annotations

import pytest

from benchmarks.held_actual_stock_oracle import SKIS
from benchmarks.held_actual_stock_oracle import circle_view_vs_actual
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_actual_stock_oracle import production_implies_certificate
from benchmarks.held_actual_stock_oracle import render_d2
from benchmarks.held_actual_stock_oracle import render_d4
from tools.actual_stock_profile import FIGURE5
from tools.actual_stock_profile import MONSTERA

# Two worker processes: a correctness tier run on a shared, loaded machine.
PROCESSES = 2
# Monstera's D4 sample stride (design section 3: every 20th circle).
MONSTERA_EVERY = 20
# Every this many circles also decides the actual stock: the C_t ⊇ S_t bound check.
AUDIT_EVERY = 50
EVERY = {FIGURE5: 1, SKIS: 1, MONSTERA: MONSTERA_EVERY}


@pytest.mark.parametrize("name", [SKIS, FIGURE5, MONSTERA])
def test_d4_circle_view_refuses_nothing_the_actual_stock_accepts(name: str) -> None:
    workload = differential_workload(name)
    report = circle_view_vs_actual(workload, range(0, len(workload.circles), EVERY[name]), PROCESSES, audit_every=AUDIT_EVERY)
    print("\n" + render_d4(report))
    assert report.compared == len(range(0, len(workload.circles), EVERY[name]))
    report.require()


@pytest.mark.parametrize("name", [SKIS])
def test_d2_every_circle_after_the_entry_is_certified(name: str) -> None:
    workload = differential_workload(name)
    report = production_implies_certificate(workload, range(len(workload.circles)), PROCESSES)
    print("\n" + render_d2(report))
    report.require()
