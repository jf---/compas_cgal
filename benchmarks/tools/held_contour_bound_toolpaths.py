"""Generate conservative contour-bound Figure 5/8 drafts from prepared inputs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from benchmarks.held_boundary_order_spacing import refine_boundary_ordered_path
from benchmarks.held_contour_bound_spacing import select_contour_bound_path
from benchmarks.held_coverage_path import build_coverage_preserving_path
from benchmarks.held_exact_motion_coverage import IncompleteMotionCoverageError
from benchmarks.held_exact_motion_coverage import report_full_motion_coverage
from benchmarks.held_figure5_boundary_path import Figure5BoundaryPath
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_motion_coverage import compare_sampled_motion_coverage
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.tools.held_contour_bound_plot import render_contour_bound_comparison
from benchmarks.tools.held_coverage_preserving_plot import render_coverage_preserving_comparison
from compas_cgal.adaptive.motion import EngagementCap


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CANONICAL_CASE_NAMES, required=True)
    parser.add_argument("--check-coverage", action="store_true", help="Also compare sampled coverage; exact full-motion coverage is always checked.")
    parser.add_argument("--qualify-contacts", action="store_true", help="Retain bounded contact representatives and report rejected hypotheses.")
    parser.add_argument("--corner-approaches", action="store_true", help="Partition concave corner approach, rotation and departure.")
    parser.add_argument("--coverage-preserving", action="store_true", help="Reproduce the rejected dense-source preservation experiment; not the selected local-repair approach.")
    args = parser.parse_args()
    if args.case == "figure5" and args.qualify_contacts:
        parser.error("Figure 5 uses its established unprojected source contacts.")
    case = load_held_reference_case(args.case)
    guide = build_held_reference_raw_guide(case)
    components = figure7_inward_offset(case).components
    source = build_figure5_approximate_path(guide, components) if case.name == "figure5" else build_held_reference_path(guide, components, qualify_contacts=args.qualify_contacts)
    print(f"{case.name}: {len(source.path.circles)} source circles; {len(source.rejected_contact_hypotheses)} rejected contact hypotheses", flush=True)
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    if case.name == "figure5":
        baseline = refine_figure5_engagement(source.path.circles, components[0], case.tool_radius, cap)
    else:
        _, baseline = refine_boundary_ordered_path(source, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches)
    selected, refined, bounds = select_contour_bound_path(
        source.path.circles,
        source.circle_sources,
        components[0],
        case.tool_radius,
        cap,
        corner_approaches=args.corner_approaches,
    )
    reached = {run for index in selected for run in source.circle_sources[index]}
    print(f"{case.name}: {len(refined.circles)} circles, {len(reached)} source runs, max bound {max(map(math.degrees, bounds)):.6f} degrees", flush=True)
    lengths = _side_lengths(components[0])
    transitions = tuple(_transition(components[0], lengths, a, b) for a, b in zip(refined.circles, refined.circles[1:]))
    baseline_transitions = tuple(_transition(components[0], lengths, a, b) for a, b in zip(baseline.circles, baseline.circles[1:]))
    if args.coverage_preserving:
        before = Figure5BoundaryPath(refined.circles, transitions, source.reached_run_ids)
        dense, corrected, bounds, stock = build_coverage_preserving_path(
            source.path.circles,
            source.circle_sources,
            components[0],
            case.projection.points,
            case.tool_radius,
            cap,
            corner_approaches=args.corner_approaches,
        )
        output = Path("docs/assets/images") / f"held_{case.name}_coverage_preserving.png"
        render_coverage_preserving_comparison(case, before, corrected, bounds, output)
        print(f"{case.name}: {len(dense.circles)} dense -> {len(corrected.circles)} retained; conservative residual empty={stock.is_empty()}; plot {output}", flush=True)
        # Publish the actual path before the more expensive absolute gate.
        residual = report_full_motion_coverage(case, corrected.circles, corrected.transitions, output.with_suffix(".exact-coverage.json"))
        if not residual.is_empty():
            raise IncompleteMotionCoverageError(f"{case.name}: coverage-preserving path still leaves full-design residual material; see exact-coverage report.")
        return
    output = Path("docs/assets/images") / f"held_{case.name}_contour_bound.png"
    render_contour_bound_comparison(case, baseline.circles, baseline_transitions, refined.circles, transitions, bounds, output)
    if args.check_coverage:
        lost, gained = compare_sampled_motion_coverage(case, baseline.circles, baseline_transitions, refined.circles, transitions, output.with_suffix(".coverage.json"))
        print(f"{case.name}: sampled coverage lost={lost}, gained={gained}; not continuous proof", flush=True)
    baseline_residual = report_full_motion_coverage(case, baseline.circles, baseline_transitions, output.with_suffix(".baseline.exact-coverage.json"))
    draft_residual = report_full_motion_coverage(case, refined.circles, transitions, output.with_suffix(".exact-coverage.json"))
    preserves_coverage = draft_residual.is_subset_of(baseline_residual)
    print(f"{case.name}: exact residual inclusion preserves baseline coverage={preserves_coverage}", flush=True)
    incomplete = []
    if not baseline_residual.is_empty():
        incomplete.append("standard baseline")
    if not draft_residual.is_empty():
        incomplete.append("contour-bound draft")
    if incomplete:
        raise IncompleteMotionCoverageError(f"{case.name}: incomplete full-motion coverage for {', '.join(incomplete)}; see exact-coverage reports.")
    print(output, flush=True)


if __name__ == "__main__":
    main()
