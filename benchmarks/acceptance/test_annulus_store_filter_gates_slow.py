"""Slow tier of units U2 and U3 of the actual-stock refinement design: the store and the filter on real paths.

Not collected by pytest (benchmarks/acceptance); run by path. Correctness at any load. The time
ceilings are asserted only on a quiet machine (ruling D2); a loaded run records its ratio and load.
"""

from __future__ import annotations

import os

import pytest

from benchmarks.held_actual_stock_oracle import SKIS
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.held_annulus_filter_gate import filter_gate
from benchmarks.held_annulus_store_gate import QUERY_CEILING_US
from benchmarks.held_annulus_store_gate import store_gate
from tools.actual_stock_profile import FIGURE5
from tools.actual_stock_profile import MONSTERA

# The one-minute load below which a budget verdict counts (ruling D2, same bound as the C1 slow tier's quiet check).
QUIET_LOAD = 8.0
FILTER_PASSES = 7
STORE_REPEATS = 5
# Nesting stride: every circle where the suite can afford it, every 20th on Monstera (design section 3).
MONSTERA_EVERY = 20


def _quiet() -> bool:
    return os.getloadavg()[0] < QUIET_LOAD


@pytest.mark.parametrize(("name", "every"), [(FIGURE5, 1), (MONSTERA, MONSTERA_EVERY)])
def test_u2_store_gate(name: str, every: int) -> None:
    report = store_gate(differential_workload(name), every, STORE_REPEATS)
    report.require()
    if _quiet():
        assert report.query_us_full_store <= QUERY_CEILING_US


@pytest.mark.parametrize("name", [FIGURE5, SKIS, MONSTERA])
def test_u3_filter_gate(name: str) -> None:
    report = filter_gate(differential_workload(name), FILTER_PASSES, ratio_every=20)
    report.require()
    assert report.exact_station_zero == 0, "U-1's station-0 fallback rate is 0"
    if _quiet():
        assert report.added_s <= report.ceiling_s
