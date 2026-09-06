"""Prepared Figure 8 inputs must reach the real circle/connector generator."""

import math

from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from compas_cgal.adaptive.motion import EngagementCap


def test_figure8_upper_emits_connected_capped_standard_draft() -> None:
    case = load_held_reference_case("figure8_upper")
    guide = build_held_reference_raw_guide(case)
    components = figure7_inward_offset(case).components
    source = build_held_reference_path(guide, components).path.circles
    cap = EngagementCap.build(math.radians(80))
    refined = refine_figure5_engagement(source, components[0], case.tool_radius, cap)
    assert len(refined.original_indices) == len(source) > 1
    lengths = _side_lengths(components[0])
    for previous, current in zip(refined.circles, refined.circles[1:]):
        transition = _transition(components[0], lengths, previous, current)
        assert transition.samples[0] == previous.contact_point
        assert transition.samples[-1] == current.contact_point
        assert float(maximum_predecessor_engagement(_paper_candidate(previous), _paper_candidate(current), case.tool_radius)) <= float(cap.theta)
