"""Render standard-model drafts on the prepared Figure 5/8 pocket corpus."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.tools.held_toolpath_progress import render_path
from compas_cgal.adaptive.motion import EngagementCap


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CANONICAL_CASE_NAMES, required=True)
    parser.add_argument("--stage", choices=("initial", "refined"), default="refined")
    args = parser.parse_args()
    case = load_held_reference_case(args.case)
    guide = build_held_reference_raw_guide(case)
    print(f"{case.name}: {len(guide.runs)} guide runs, {guide.station_count} stations", flush=True)
    components = figure7_inward_offset(case).components
    result = build_held_reference_path(guide, components)
    circles = result.path.circles
    print(f"{case.name}: {len(circles)} initial circles", flush=True)
    if args.stage == "refined":
        circles = refine_figure5_engagement(circles, components[0], case.tool_radius, EngagementCap.build(math.radians(float(case.tea_cap)))).circles
    lengths = _side_lengths(components[0])
    transitions = tuple(_transition(components[0], lengths, a, b) for a, b in zip(circles, circles[1:]))
    output = Path("docs/assets/images") / f"held_{case.name}_standard_{args.stage}.png"
    render_path(
        case,
        circles,
        transitions,
        output,
        title=f"{case.name.replace('_', ' ')} · standard-model {args.stage} draft",
        scope="Publisher-start polygon guide · standard predecessor model · contour-aware reproduction and containment unqualified.",
    )
    print(output, flush=True)


if __name__ == "__main__":
    main()
