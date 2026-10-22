"""The exact NOTHING-circle pass runs inside Ω* generation as its own timed stage, before the connectors."""

from __future__ import annotations

import math
from functools import cache

import pytest

import benchmarks.held_omega_toolpath as held_omega_toolpath
from benchmarks.held_nothing_circles import NothingCirclePass
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.held_omega_toolpath import remove_omega_nothing_circles
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.motion import EngagementCap

pytestmark = pytest.mark.xdist_group("held_omega_figure5")

# Figure 5 Ω* circles the pass proves NOTHING (indices into the refined
# circles; NR's acceptance table, nr-nothing-circles-report.md).
FIGURE5_OMEGA_NOTHING = (851, 867, 884, 892, 900)


@cache
def _figure5() -> OmegaGeneratedToolpath:
    return generate_omega("figure5_refined")


def _cap(path: OmegaGeneratedToolpath) -> EngagementCap:
    return EngagementCap.build(math.radians(float(path.case.tea_cap)))


def test_an_exact_duplicate_omega_circle_is_removed() -> None:
    path = _figure5()
    first, second, third = path.refinement.circles[:3]
    # Premise of the closed pass: removing the duplicate makes third -> first the seam, which must stay within the cap.
    assert exact_predecessor_engagement(third.paper_candidate(), first.paper_candidate(), path.case.tool_radius).within(_cap(path))
    result = remove_omega_nothing_circles((first, second, third, third), path.case.tool_radius, _cap(path))
    assert result.removed == (3,)
    assert result.circles == (first, second, third)


def test_generation_times_the_nothing_stage_between_refine_and_connector_repair() -> None:
    stages = list(_figure5().stage_seconds)
    assert stages.index("refine") < stages.index("nothing") < stages.index("connector_repair")


def test_generation_emits_and_connects_the_surviving_circles() -> None:
    path = _figure5()
    assert path.nothing.removed == FIGURE5_OMEGA_NOTHING
    refined = path.refinement.circles
    assert path.nothing.circles == tuple(circle for k, circle in enumerate(refined) if k not in FIGURE5_OMEGA_NOTHING)
    assert len(path.connectors.connectors) == len(path.circles)
    report = path.audit()
    assert not report.gouging_circles
    assert not report.gouging_connectors


def test_generation_runs_the_nothing_pass_on_a_closed_path(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[bool] = []
    original = held_omega_toolpath.remove_nothing_circles

    def recording(*args: object, closed: bool) -> NothingCirclePass[object]:
        calls.append(closed)
        return original(*args, closed=closed)  # type: ignore[arg-type]

    monkeypatch.setattr(held_omega_toolpath, "remove_nothing_circles", recording)
    path = generate_omega("figure5_refined")
    assert calls == [True]
    assert path.nothing.removed == FIGURE5_OMEGA_NOTHING
