"""Slow tier of unit U4 of the actual-stock refinement design: D2 and the call count on real paths.

Not collected by pytest (benchmarks/acceptance); run by path. Correctness at any load: D2
is exact, so the machine's load only changes how long it takes.
"""

from __future__ import annotations

import pytest

from benchmarks.held_actual_stock_oracle import SKIS
from benchmarks.held_refinement_decision_gate import decision_gate
from tools.actual_stock_profile import FIGURE5
from tools.actual_stock_profile import MONSTERA

# Certificate processes per workload; U4's gate is every circle on fig5 and skis, every 20th on Monstera (D2-all is U5's).
PROCESSES = 4
MONSTERA_EVERY = 20


@pytest.mark.parametrize(("name", "every"), [(FIGURE5, 1), (SKIS, 1), (MONSTERA, MONSTERA_EVERY)])
def test_u4_production_implies_certificate(name: str, every: int) -> None:
    decision_gate(name, every, PROCESSES).require()
