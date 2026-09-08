"""Wall-clock witnesses for the exact full-motion coverage replay on the Figure 5 draft.

Measured 2026-09-08 before the change: 400 circles with their connectors
replayed in 5.4 s and the full 2,013-circle draft had not finished after
eighteen minutes, because `remaining_material` differenced every sweep one at
a time against a residual that grows with every circle.
"""

from __future__ import annotations

import math
import time

import pytest

from benchmarks.held_exact_motion_coverage import replay_exact_motion_coverage
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal.adaptive.motion import EngagementCap

# One divide-and-conquer union of all sweeps and a single difference should
# replay 400 circles in well under a second on the recording machine; 2 s
# keeps room for CI while staying below the 5.4 s measured before.
PREFIX_400_BUDGET_S = 2.0
# The whole draft: generous against CI, three orders of magnitude below the
# unfinished eighteen minutes.
FULL_DRAFT_BUDGET_S = 60.0


@pytest.fixture(scope="module")
def figure5_draft():
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)
    components = figure7_inward_offset(case).components
    result = build_figure5_approximate_path(guide, components)
    refined = refine_figure5_engagement(result.path.circles, components[0], case.tool_radius, EngagementCap.build(math.radians(float(case.tea_cap))))
    lengths = _side_lengths(components[0])
    transitions = tuple(_transition(components[0], lengths, a, b) for a, b in zip(refined.circles, refined.circles[1:]))
    return case, components[0], refined.circles, transitions


def test_replay_of_400_circles_stays_within_budget(figure5_draft) -> None:
    case, boundary, circles, transitions = figure5_draft
    started = time.perf_counter()
    residual = replay_exact_motion_coverage(boundary, circles[:400], transitions[:399], case.tool_radius)
    seconds = time.perf_counter() - started
    assert residual.component_count() >= 1  # a prefix never clears the pocket
    assert seconds < PREFIX_400_BUDGET_S


def test_replay_of_the_full_draft_stays_within_budget(figure5_draft) -> None:
    case, boundary, circles, transitions = figure5_draft
    started = time.perf_counter()
    replay_exact_motion_coverage(boundary, circles, transitions, case.tool_radius)
    assert time.perf_counter() - started < FULL_DRAFT_BUDGET_S
