"""Render standard-model drafts on the prepared Figure 5/8 pocket corpus."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from benchmarks.held_boundary_order_spacing import refine_boundary_ordered_path
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.tools.held_circle_distribution import render_distribution
from benchmarks.tools.held_toolpath_progress import render_path
from compas_cgal.adaptive.motion import EngagementCap


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CANONICAL_CASE_NAMES, required=True)
    parser.add_argument("--stage", choices=("initial", "refined"), default="refined")
    parser.add_argument("--qualify-contacts", action="store_true", help="Record out-of-bound contact alternatives; require a bounded representative for every source run.")
    parser.add_argument("--corner-approaches", action="store_true", help="Partition concave corner approach, rotation and departure before refinement.")
    parser.add_argument("--spacing", choices=("lane", "boundary"), default="lane", help="Boundary is an experimental reselection in emitted order, with source-run retention.")
    args = parser.parse_args()
    if args.spacing == "boundary" and args.stage != "refined":
        parser.error("Boundary spacing requires the refined stage so no unresolved gaps are emitted.")
    if args.corner_approaches and args.stage != "refined":
        parser.error("Corner approaches require the refined stage.")
    case = load_held_reference_case(args.case)
    guide = build_held_reference_raw_guide(case)
    print(f"{case.name}: {len(guide.runs)} guide runs, {guide.station_count} stations", flush=True)
    components = figure7_inward_offset(case).components
    result = build_held_reference_path(guide, components, qualify_contacts=args.qualify_contacts)
    print(f"{case.name}: {len(result.rejected_contact_hypotheses)} rejected contact alternatives retained as evidence", flush=True)
    circles = result.path.circles
    print(f"{case.name}: {len(circles)} initial circles", flush=True)
    if args.stage == "refined":
        cap = EngagementCap.build(math.radians(float(case.tea_cap)))
        if args.spacing == "boundary":
            selected, refined = refine_boundary_ordered_path(result, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches)
            print(f"{case.name}: {len(selected)} boundary-order sources, {len(refined.circles)} refined circles", flush=True)
        else:
            refined = refine_figure5_engagement(circles, components[0], case.tool_radius, cap, corner_approaches=args.corner_approaches)
            render_distribution(case, result, refined, Path("docs/assets/images") / f"held_{case.name}_circle_distribution.png")
        circles = refined.circles
    lengths = _side_lengths(components[0])
    transitions = tuple(_transition(components[0], lengths, a, b) for a, b in zip(circles, circles[1:]))
    suffix = f"{args.stage}_boundary" if args.spacing == "boundary" else args.stage
    if args.qualify_contacts:
        suffix += "_qualified_contacts"
    if args.corner_approaches:
        suffix += "_corner_approaches"
    output = Path("docs/assets/images") / f"held_{case.name}_standard_{suffix}.png"
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
        title=f"{case.name.replace('_', ' ')} · standard-model {args.stage} draft · {args.spacing} spacing",
        scope="Publisher-start polygon guide · standard predecessor model · contour-aware reproduction and containment unqualified.",
    )
    print(output, flush=True)


if __name__ == "__main__":
    main()
