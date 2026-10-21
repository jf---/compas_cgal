"""Slow-tier acceptance: every published Held workload generated on Ω* has zero exact gouges.

Not collected by pytest (benchmarks/acceptance is outside ``testpaths``). Run
with the editable-skip invocation and this path. The audit is exact against
the declared boundary (:mod:`benchmarks.held_omega_gouge_audit`): circles are
binary64 and may touch (dist == r); connectors run on ∂Ω* pieces whose exact
clearance is r.
"""

from __future__ import annotations

import pytest

from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_omega_toolpath import OMEGA_WORKLOADS
from benchmarks.held_omega_toolpath import generate_omega


@pytest.mark.parametrize("workload", list(OMEGA_WORKLOADS))
def test_omega_path_has_no_exact_gouge(workload: str) -> None:
    path = generate_omega(workload)
    report = path.audit()
    assert not report.gouging_circles, report.gouging_circles[:5]
    assert not report.gouging_connectors, report.gouging_connectors[:5]
    assert path.connectors.seam is Figure5SeamClosure.CONNECTOR or len(path.connectors.connectors) == len(path.circles) - 1
    assert path.initial.reached_run_ids == frozenset(run for runs in path.initial.circle_sources for run in runs)
