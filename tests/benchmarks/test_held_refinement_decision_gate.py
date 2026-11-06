"""The U4 gate harness: its counts on a real generation and each way it refuses."""

from __future__ import annotations

from dataclasses import replace

import pytest

from benchmarks.held_figure5_engagement_trace import DecisionSite
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_refinement_decision_gate import CallCount
from benchmarks.held_refinement_decision_gate import DecisionCallCountError
from benchmarks.held_refinement_decision_gate import call_count
from benchmarks.held_refinement_decision_gate import generate_with_refinement

GOOD = CallCount(10, 9, 0, 9, 9, 0, 0, (), 9, 0, 405, 400, 5, 0, 0, 0, 0)


@pytest.mark.parametrize("workload", ["figure5_refined", "figure8_crossed_skis_boundary"])
def test_a_real_generation_decides_every_refined_circle_on_the_part(workload: str) -> None:
    path, refinement, cpu = generate_with_refinement(workload)
    assert cpu > 0.0
    assert refinement.trace.key.part == PartPolygon.from_case(path.case), "the generator passes P, not the centre domain"
    assert refinement.trace.key.part.points != refinement.trace.key.boundary
    counts = call_count(refinement, len(path.circles))
    counts.require()
    assert counts.continuation == 0
    assert counts.decisions == counts.advance + counts.close_seam == refinement.trace.filter_counters.decisions
    assert refinement.trace.decisions_at(DecisionSite.ADVANCE) == counts.advance


@pytest.mark.parametrize(
    "broken",
    [
        replace(GOOD, decisions=7, engine_decisions=7, advance=7),
        replace(GOOD, engine_decisions=8),
        replace(GOOD, undecided_circles=(3,)),
        replace(GOOD, advance=0),
        replace(GOOD, seam_circles=1),
    ],
)
def test_the_call_count_refuses_each_breach(broken: CallCount) -> None:
    GOOD.require()
    with pytest.raises(DecisionCallCountError):
        broken.require()
