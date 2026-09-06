"""Expensive four-pocket coverage acceptance; every nonempty residual fails.

Run explicitly with ``pixi run held-coverage-corpus``. These machining workloads
are separate from construction diagrams and intermediate stock illustrations.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Literal
from typing import cast

import pytest

from benchmarks.held_boundary_order_spacing import refine_boundary_ordered_path
from benchmarks.held_contour_bound_spacing import select_contour_bound_path
from benchmarks.held_coverage_path import build_coverage_preserving_path
from benchmarks.held_exact_motion_coverage import require_full_motion_coverage
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

PreparedPaths = tuple[
    HeldReferenceCase,
    tuple[Point2[WorldXY], ...],
    tuple[Figure5CounterclockwiseCircle, ...],
    tuple[Figure5CounterclockwiseCircle, ...],
    Figure5BoundaryPath,
]


@pytest.fixture(scope="module", params=[pytest.param(name, marks=pytest.mark.xdist_group(name)) for name in CANONICAL_CASE_NAMES])
def prepared_paths(request: pytest.FixtureRequest) -> PreparedPaths:
    """Generate each source and its three placements once per pocket worker."""
    case = load_held_reference_case(cast(str, request.param))
    guide = build_held_reference_raw_guide(case)
    components = figure7_inward_offset(case).components
    boundary = components[0]
    source = (
        build_figure5_approximate_path(guide, components)
        if case.name == "figure5"
        else build_held_reference_path(guide, components, qualify_contacts=case.name == "figure8_monstera")
    )
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    if case.name == "figure5":
        baseline = refine_figure5_engagement(source.path.circles, boundary, case.tool_radius, cap)
    else:
        _, baseline = refine_boundary_ordered_path(source, boundary, case.tool_radius, cap, corner_approaches=True)
    _, draft, _ = select_contour_bound_path(source.path.circles, source.circle_sources, boundary, case.tool_radius, cap, corner_approaches=True)
    _, corrected, _, _ = build_coverage_preserving_path(
        source.path.circles, source.circle_sources, boundary, case.projection.points, case.tool_radius, cap, corner_approaches=True
    )
    return case, boundary, baseline.circles, draft.circles, corrected


@pytest.mark.parametrize("algorithm", ("standard", "contour_bound", "coverage_preserving"))
def test_full_declared_pocket_is_covered(prepared_paths: PreparedPaths, algorithm: Literal["standard", "contour_bound", "coverage_preserving"]) -> None:
    """Require continuous cutter coverage, including every emitted connector."""
    case, boundary, baseline, draft, corrected = prepared_paths
    output = Path("build/held-coverage-corpus") / f"{case.name}_{algorithm}.json"
    if algorithm == "coverage_preserving":
        # Preserve the emitted connector chain; rebuilding direct connectors
        # here would test different motion from the deletion algorithm.
        require_full_motion_coverage(case, corrected.circles, corrected.transitions, output)
        return
    circles = baseline if algorithm == "standard" else draft
    lengths = _side_lengths(boundary)
    transitions = tuple(_transition(boundary, lengths, first, second) for first, second in zip(circles, circles[1:]))
    require_full_motion_coverage(case, circles, transitions, output)
