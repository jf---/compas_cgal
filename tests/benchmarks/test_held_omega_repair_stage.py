"""Ω* generation repairs every connector the exact gate refuses, as its own timed stage after the NOTHING pass."""

from __future__ import annotations

from functools import cache

import pytest

from benchmarks.held_omega_connector_engagement import assess_omega_connectors
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import generate_omega

pytestmark = pytest.mark.xdist_group("held_omega_figure5")

# Figure 5 Ω*: one connector refused after the NOTHING pass, repaired by one
# ∂Ω* midpoint circle (or-omega-connector-repair-report.md, acceptance table).
FIGURE5_REPAIR_INSERTED = 1


@cache
def _figure5() -> OmegaGeneratedToolpath:
    return generate_omega("figure5_refined")


def test_generation_times_connector_repair_after_the_nothing_stage() -> None:
    stages = list(_figure5().stage_seconds)
    assert stages.index("refine") < stages.index("nothing") < stages.index("connector_repair")
    assert stages[-1] == "connector_repair"
    assert "connectors" not in stages


def test_generation_emits_the_repaired_circles_and_their_connectors() -> None:
    path = _figure5()
    assert path.circles == path.repair.circles
    assert path.connectors is path.repair.connectors
    assert path.repair.inserted == FIGURE5_REPAIR_INSERTED
    assert len(path.circles) == len(path.nothing.circles) + FIGURE5_REPAIR_INSERTED


def test_every_emitted_connector_is_certified_and_the_repair_verdicts_are_the_post_gate() -> None:
    path = _figure5()
    fresh = assess_omega_connectors(path)
    assert fresh.refused == ()
    assert [verdict.tier for verdict in path.repair.verdicts] == [verdict.tier for verdict in fresh.verdicts]
    report = path.audit()
    assert not report.gouging_circles
    assert not report.gouging_connectors
