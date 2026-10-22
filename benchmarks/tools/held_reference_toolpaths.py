"""Render standard-model drafts on the prepared Figure 5/8 pocket corpus.

A refined draft whose options name a benchmark workload is that workload's
``held_toolpath_benchmark.generate`` path (refinement with the closed seam,
NOTHING pass, connector post-gate), so a published panel is the generator's
path. ``--certified-repair`` draws the certificate-gated coverage repair
instead, and records the sources it re-inserted.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from benchmarks.held_boundary_order_spacing import refine_boundary_ordered_path
from benchmarks.held_coverage_repair import refine_boundary_ordered_path_certified
from benchmarks.held_exact_motion_coverage import require_full_motion_coverage
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_toolpath_benchmark import certify_connectors
from benchmarks.held_toolpath_benchmark import generate
from benchmarks.tools.held_circle_distribution import render_distribution
from benchmarks.tools.held_contour_contact_plot import render_contour_contacts
from benchmarks.tools.held_sequence_coverage import certify_path_coverage
from benchmarks.tools.held_toolpath_progress import render_path
from compas_cgal.motion import EngagementCap

# The published figure directory; only generator paths (or the certified repair) are written there.
PUBLISHED_IMAGES = Path("docs/assets/images")
# Refined drafts that are benchmark workloads, by (case, spacing, qualify contacts, corner approaches).
GENERATOR_WORKLOADS = {
    ("figure8_upper", "lane", False, False): "figure8_upper_lane",
    ("figure8_upper", "boundary", False, False): "figure8_upper_boundary",
    ("figure8_crossed_skis", "boundary", False, False): "figure8_crossed_skis_boundary",
    ("figure8_monstera", "boundary", True, True): "figure8_monstera_boundary_qualified_corners",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CANONICAL_CASE_NAMES, required=True)
    parser.add_argument("--contour-prefix", type=int, help="Inspect exposed-contour contact corrections on this many emitted circles.")
    parser.add_argument("--stage", choices=("initial", "refined"), default="refined")
    parser.add_argument("--qualify-contacts", action="store_true", help="Record out-of-bound contact alternatives; require a bounded representative for every source run.")
    parser.add_argument("--corner-approaches", action="store_true", help="Partition concave corner approach, rotation and departure before refinement.")
    parser.add_argument("--check-coverage", action="store_true", help="Also run the exact full-motion coverage replay; unworkable at this size, see docs/held_motion_coverage.md.")
    parser.add_argument("--spacing", choices=("lane", "boundary"), default="lane", help="Boundary is an experimental reselection in emitted order, with source-run retention.")
    parser.add_argument("--certified-repair", action="store_true", help="Boundary spacing only: the certificate-gated coverage repair (refine_boundary_ordered_path_certified).")
    parser.add_argument("--output-dir", type=Path, default=PUBLISHED_IMAGES, help="Directory every figure and report is written to.")
    args = parser.parse_args()
    if args.contour_prefix is not None and args.contour_prefix < 2:
        parser.error("Contour prefix requires at least two circles.")
    if args.spacing == "boundary" and args.stage != "refined":
        parser.error("Boundary spacing requires the refined stage so no unresolved gaps are emitted.")
    if args.corner_approaches and args.stage != "refined":
        parser.error("Corner approaches require the refined stage.")
    if args.certified_repair and args.spacing != "boundary":
        parser.error("The certified repair acts on the boundary-spacing selection.")
    workload = GENERATOR_WORKLOADS.get((args.case, args.spacing, args.qualify_contacts, args.corner_approaches)) if args.stage == "refined" else None
    if args.stage == "refined" and workload is None and not args.certified_repair and args.output_dir.resolve() == PUBLISHED_IMAGES.resolve():
        parser.error("These options name no generator workload, so the path is not the generator's; pass --output-dir outside docs/assets/images.")
    case = load_held_reference_case(args.case)
    components = figure7_inward_offset(case).components
    # The generator builds its own guide and source path; build them here only for what reads them.
    needs_source = args.stage != "refined" or workload is None or args.certified_repair or args.qualify_contacts or args.spacing == "lane"
    if needs_source:
        guide = build_held_reference_raw_guide(case)
        print(f"{case.name}: {len(guide.runs)} guide runs, {guide.station_count} stations", flush=True)
        result = build_held_reference_path(guide, components, qualify_contacts=args.qualify_contacts)
        print(f"{case.name}: {len(result.rejected_contact_hypotheses)} rejected contact alternatives retained as evidence", flush=True)
        circles = result.path.circles
        print(f"{case.name}: {len(circles)} initial circles", flush=True)
    inserted_sources: tuple[int, ...] | None = None
    if args.stage == "refined":
        cap = EngagementCap.build(math.radians(float(case.tea_cap)))
        if args.certified_repair:
            certified = refine_boundary_ordered_path_certified(result, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches)
            circles = certified.circles
            inserted_sources = tuple(certified.ledger.inserted_sources)
            print(f"{case.name}: certified repair re-inserted sources {inserted_sources}, {len(circles)} emitted circles", flush=True)
        elif workload is not None:
            # The published path is the benchmark generator's own, every stage included, not a re-derivation.
            circles = generate(workload).circles
            print(f"{case.name}: generator workload {workload}, {len(circles)} emitted circles", flush=True)
        else:
            print(f"{case.name}: no generator workload for these options; refining here without the generator's later stages", flush=True)
            if args.spacing == "boundary":
                _, refined = refine_boundary_ordered_path(result, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches, close_seam=True)
            else:
                refined = refine_figure5_engagement(circles, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches, close_seam=True)
            circles = refined.circles
            certify_connectors(circles, Figure5BoundaryConnectors.build(components[0], circles), case.tool_radius, cap)
        if args.spacing == "lane":
            # Diagnostic of the refinement stage itself (before the generator's NOTHING pass).
            refined = refine_figure5_engagement(result.path.circles, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches, close_seam=True)
            render_distribution(
                case, result, refined, Figure5BoundaryConnectors.build(components[0], refined.circles), args.output_dir / f"held_{case.name}_circle_distribution.png"
            )
    transitions = Figure5BoundaryConnectors.build(components[0], circles).transitions
    suffix = f"{args.stage}_boundary" if args.spacing == "boundary" else args.stage
    if args.qualify_contacts:
        suffix += "_qualified_contacts"
    if args.corner_approaches:
        suffix += "_corner_approaches"
    if args.certified_repair:
        suffix += "_certified_repair"
    output = args.output_dir / f"held_{case.name}_standard_{suffix}.png"
    if args.contour_prefix is not None:
        if args.contour_prefix > len(circles):
            parser.error("Contour prefix exceeds emitted circle count.")
        render_contour_contacts(case, circles[: args.contour_prefix], args.output_dir / f"held_{case.name}_contour_contacts_{args.contour_prefix}.png")
    if args.qualify_contacts:
        output.parent.mkdir(parents=True, exist_ok=True)
        contact_report = {
            "case": case.name,
            "evidence_bound_mm": float(guide.site_budget.admissible_distance_gap),
            "retained_source_runs": len(result.reached_run_ids),
            "rejection_reason": "CGAL offset-contact distance exceeds evidence bound",
            "rejected_alternatives": [
                {
                    "run": int(hypothesis.run_id),
                    "station": int(hypothesis.station_ordinal),
                    "source_side": int(hypothesis.projected_boundary_site.segment_id),
                    "contact_xy_mm": [float(hypothesis.contact_point.x), float(hypothesis.contact_point.y)],
                }
                for hypothesis in result.rejected_contact_hypotheses
            ],
        }
        output.with_suffix(".contacts.json").write_text(json.dumps(contact_report, indent=2, allow_nan=False) + "\n")
    render_path(
        case,
        circles,
        transitions,
        output,
        title=(
            f"{case.name.replace('_', ' ')} · {args.spacing} spacing · certified coverage repair"
            if args.certified_repair
            else f"{case.name.replace('_', ' ')} · standard-model {args.stage} draft · {args.spacing} spacing"
        ),
        scope="Publisher-start polygon guide · standard predecessor model · contour-aware reproduction and containment unqualified.",
        coverage=certify_path_coverage(components[0], circles, case.tool_radius),
        inserted_sources=inserted_sources,
    )
    if args.check_coverage:
        require_full_motion_coverage(case, circles, transitions, output.with_suffix(".exact-coverage.json"))
    else:
        print("exact full-motion coverage NOT checked (set aside 2026-09-08; pass --check-coverage to run the replay)", flush=True)
    print(output, flush=True)


if __name__ == "__main__":
    main()
