"""Slow tier of D1 (engine == Stock2 on the circle view), unit U1 of the actual-stock refinement design.

Not collected by pytest (benchmarks/acceptance); run by path. Correctness only,
no budget, any load. U1's exit gate is the first 200 Figure 5 circles, every
station, on C_t. The full D1 table of design section 3 (Figure 5 and skis every
circle, Monstera every 20th) is reached through the same harness.
"""

from __future__ import annotations

import pytest

from benchmarks.held_actual_stock_oracle import SKIS
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_annulus_view_d1 import d1_report
from tools.actual_stock_profile import FIGURE5
from tools.actual_stock_profile import MONSTERA

U1_GATE_CIRCLES = 200
STATIONS = 45
# Monstera's D1 sample stride (design section 3: every 20th circle).
MONSTERA_EVERY = 20


def test_u1_gate_first_figure5_circles_every_station_equal() -> None:
    report = d1_report(differential_workload(FIGURE5), range(U1_GATE_CIRCLES), last=None)
    assert report.stations == U1_GATE_CIRCLES * STATIONS
    report.require()


@pytest.mark.parametrize(("name", "every"), [(FIGURE5, 1), (SKIS, 1), (MONSTERA, MONSTERA_EVERY)])
def test_d1_full_table(name: str, every: int) -> None:
    workload = differential_workload(name)
    report = d1_report(workload, range(0, len(workload.circles), every), last=None)
    assert report.stations == len(range(0, len(workload.circles), every)) * STATIONS
    report.require()
