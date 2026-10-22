"""Generate, cache and render a Held workload as a sequence plot (swept bands in emission order + certificate verdict).

Generation (minutes for Monstera) is cached as ``SequencePlotData`` JSON so a
figure re-renders in seconds. The repaired Monstera workload runs the
certificate-gated repair (``refine_boundary_ordered_path_certified``).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from benchmarks.held_coverage_repair import refine_boundary_ordered_path_certified
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_toolpath_benchmark import GeneratedToolpath
from benchmarks.held_toolpath_benchmark import generate
from benchmarks.tools.held_sequence_coverage import certify_path_coverage
from benchmarks.tools.held_sequence_plot import ColourMode
from benchmarks.tools.held_sequence_plot import SequencePlotData
from benchmarks.tools.held_sequence_plot import emitted_engagement
from benchmarks.tools.held_sequence_plot import plot_data
from benchmarks.tools.held_sequence_plot import render_sequence
from benchmarks.tools.held_sequence_plot import render_side_by_side
from compas_cgal.motion import EngagementCap

REPAIRED_MONSTERA = "figure8_monstera_boundary_qualified_corners_repaired"
TITLES = {
    "figure5_refined": "Figure 5 · emission order of the refined path",
    "figure8_upper_lane": "Upper pocket · emission order of the lane-spaced path",
    "figure8_upper_boundary": "Upper pocket · emission order of the boundary-ordered path",
    "figure8_crossed_skis_boundary": "Crossed skis · emission order of the boundary-ordered path",
    "figure8_monstera_boundary_qualified_corners": "Monstera · emission order of the boundary-ordered path",
    REPAIRED_MONSTERA: "Monstera · emission order of the certificate-repaired path",
}


def _repaired_monstera() -> GeneratedToolpath:
    case = load_held_reference_case("figure8_monstera")
    components = figure7_inward_offset(case).components
    source = build_held_reference_path(build_held_reference_raw_guide(case), components, qualify_contacts=True)
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    certified = refine_boundary_ordered_path_certified(source, components[0], case.tool_radius, cap, corner_approaches=True)
    circles = certified.circles
    return GeneratedToolpath(REPAIRED_MONSTERA, case, components[0], circles, Figure5BoundaryConnectors.build(components[0], circles), {})


def sequence_data(path: GeneratedToolpath, title: str) -> SequencePlotData:
    """Report a generated path as plot doubles: engagement recomputed in emitted order, certificate run on it."""
    engagement, _ = emitted_engagement(path.case, path.circles, path.transitions)
    coverage = certify_path_coverage(path.boundary, path.circles, path.case.tool_radius)
    return plot_data(workload=path.workload, title=title, case=path.case, circles=path.circles, transitions=path.transitions, engagement_deg=engagement, coverage=coverage)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workload", required=True, choices=sorted(TITLES))
    parser.add_argument("--cache", type=Path, default=Path("build/seq-plot/data"), help="Directory of generated sequence data (JSON).")
    parser.add_argument("--regenerate", action="store_true", help="Run the generator even if cached data exists.")
    parser.add_argument("--output", type=Path, default=Path("build/seq-plot/new"))
    parser.add_argument("--colour", choices=[mode.value for mode in ColourMode], default=ColourMode.ORDER.value, help="Ramp span: whole path or each guide run.")
    parser.add_argument("--compare", type=Path, help="Existing PNG to set beside the new render in build/seq-plot/cmp/.")
    args = parser.parse_args()
    cached = args.cache / f"{args.workload}.json"
    if args.regenerate or not cached.exists():
        path = _repaired_monstera() if args.workload == REPAIRED_MONSTERA else generate(args.workload)
        data = sequence_data(path, TITLES[args.workload])
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_text(data.to_json())
    else:
        data = SequencePlotData.from_json(cached.read_text())
    mode = ColourMode(args.colour)
    output = args.output / f"held_{args.workload}_sequence_{mode.value}_swept.png"
    render_sequence(data, output, mode)
    print(output)
    if args.compare is not None:
        comparison = Path("build/seq-plot/cmp") / f"held_{args.workload}_{mode.value}_swept_old_vs_new.png"
        render_side_by_side(args.compare, output, comparison)
        print(comparison)


if __name__ == "__main__":
    main()
