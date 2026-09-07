"""Bounded circle construction directly on native imported curved boundaries."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from benchmarks.errors import BenchmarkError
from benchmarks.held_native_boundary_draft import CHECKPOINT_SECONDS
from benchmarks.held_native_boundary_draft import MINIMUM_DISPLAY_SPAN_MM
from benchmarks.held_native_boundary_draft import _draw_primitive
from benchmarks.held_native_curve_import import import_held_boundary
from benchmarks.held_native_curve_import_plot import _draw_source
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import _circle_geometry_2
from compas_cgal import _coverage_2 as native


class InvalidCurvedCircleDiagnosticLimitError(BenchmarkError):
    """The requested diagnostic prefix must contain at least one native piece."""


def _render(
    case: HeldReferenceCase,
    owner: native.NativeBoundary2 | None,
    events: list[tuple[int, _circle_geometry_2.BoundaryNormalCircleProposal2, float]],
    output: Path,
    *,
    elapsed: float,
    finished: bool,
    failure: RuntimeError | None = None,
    failure_site: tuple[float, float] | None = None,
) -> None:
    """Render retained reporting views; never reconstruct deciding geometry."""
    figure = Figure(figsize=(13, 7), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    overview, zoom = figure.subplots(1, 2)
    for axes in (overview, zoom):
        if owner is None:
            _draw_source(axes, case)
        else:
            for primitive in owner.cycle.primitives:
                _draw_primitive(axes, primitive, color="#34403f", linewidth=1)
        for _, proposal, _ in events:
            if proposal.is_stationary:
                axes.scatter(*proposal.center_mm, color="#c0392b", marker="x", s=22)
            else:
                axes.add_patch(Circle(proposal.center_mm, proposal.guide_radius_mm, fill=False, edgecolor="#087e8b", linewidth=0.8))
            axes.scatter(*proposal.m_mm, color="#a87432", s=13)
        if failure_site is not None:
            axes.scatter(*failure_site, color="#c0392b", marker="X", s=70)
        axes.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    extents = [(p.center_mm[0] + sign * p.guide_radius_mm, p.center_mm[1] + sign * p.guide_radius_mm) for _, p, _ in events for sign in (-1, 1)]
    if failure_site is not None:
        extents.append(failure_site)
    if extents:
        coordinates = np.array(extents)
        low, high = coordinates.min(axis=0), coordinates.max(axis=0)
        margin = max(float(np.max(high - low)), MINIMUM_DISPLAY_SPAN_MM) / 4
        zoom.set(xlim=(low[0] - margin, high[0] + margin), ylim=(low[1] - margin, high[1] + margin))
    overview.set_title("Imported native boundary and sampled circles")
    zoom.set_title("Completed circles · gold dots: native medial contacts")
    total = 0 if owner is None else len(owner.cycle.primitives)
    status = "FAILED CONSTRUCTION" if failure else "CONSTRUCTION DIAGNOSTIC — NOT A PATH"
    figure.suptitle(f"{case.name.replace('_', ' ')} · {status}", x=0.065, ha="left", fontsize=16)
    figure.text(0.065, 0.9, f"{len(events)}/{total} native pieces queried · {sum(p.is_stationary for _, p, _ in events)} stationary · {elapsed:.2f} s", fontsize=11)
    if failure is not None:
        figure.text(0.065, 0.07, f"{type(failure).__name__}: {str(failure)[:140]}", color="#c0392b", fontsize=9)
    figure.text(0.065, 0.025, "Actual tool-radius offset, connectors, engagement, coverage and entry remain unqualified. Reporting coordinates only.", fontsize=9)
    figure.subplots_adjust(left=0.065, right=0.98, top=0.83, bottom=0.15, wspace=0.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    report = {
        "case": case.name,
        "algorithm": "native_curved_boundary_circle_diagnostic",
        "finished_requested_queries": finished and failure is None,
        "all_native_pieces_queried": finished and failure is None and len(events) == total,
        "native_pieces": total,
        "completed_queries": len(events),
        "elapsed_seconds": elapsed,
        "failure": None if failure is None else {"type": type(failure).__name__, "message": str(failure), "source_point_mm": failure_site},
        "events": [
            {
                "piece_index": index,
                "query_seconds": seconds,
                "stationary": proposal.is_stationary,
                "p_mm": proposal.p_mm,
                "q_mm": proposal.q_mm,
                "m_mm": proposal.m_mm,
                "center_mm": proposal.center_mm,
                "guide_radius_mm": proposal.guide_radius_mm,
                "clearance_mm": proposal.clearance_mm,
                "competing_arc_indices": proposal.competing_arc_indices,
            }
            for index, proposal, seconds in events
        ],
        "scope": "Native circles on imported curved boundary. No actual tool-radius offset, connector, engagement, coverage or entry qualification.",
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"PLOT {output}: {len(events)}/{total} queries; failed={failure is not None}", flush=True)


def render_native_curved_circles(case_name: str, *, max_pieces: int | None = None, output: Path | None = None) -> None:
    """Query midpoint circles directly on native pieces, keeping failures visible."""
    if max_pieces is not None and max_pieces < 1:
        raise InvalidCurvedCircleDiagnosticLimitError("max_pieces must be positive when supplied.")
    case = load_held_reference_case(case_name)
    output = output or Path(f"docs/assets/images/held_{case.name}_native_curved_circles.png")
    events: list[tuple[int, _circle_geometry_2.BoundaryNormalCircleProposal2, float]] = []
    owner: native.NativeBoundary2 | None = None
    site: tuple[float, float] | None = None
    started = time.perf_counter()
    checkpoint = started
    try:
        owner = import_held_boundary(case)
        for index, primitive in enumerate(owner.cycle.primitives[:max_pieces]):
            # Reporting sample identifies failed queries; it is never reinjected.
            site = primitive.sample(0.5).reporting_xy_mm
            query_started = time.perf_counter()
            print(f"{case.name}: native piece {index + 1} midpoint query", flush=True)
            proposal = owner.circle_on_piece(index, 0.5, float(case.tool_radius.value))
            now = time.perf_counter()
            events.append((index, proposal, now - query_started))
            if len(events) <= 3 or now - checkpoint >= CHECKPOINT_SECONDS:
                _render(case, owner, events, output, elapsed=now - started, finished=False)
                checkpoint = now
    except (
        native.InvalidNativeBoundaryCurveError,
        native.InvalidNativeBoundaryChainError,
        native.ReachableDomainConstructionError,
        native.InvalidBoundaryCircleContactError,
        native.BoundaryContactConstructionError,
        native.CoverageTransitionError,
        native.InvalidCoverageGeometryError,
        native.InvalidNativeBoundaryMedialInputError,
        native.NoPositiveNativeBoundaryCircleError,
        native.NativeBoundaryMedialConstructionError,
    ) as error:
        _render(case, owner, events, output, elapsed=time.perf_counter() - started, finished=False, failure=error, failure_site=site)
        raise
    _render(case, owner, events, output, elapsed=time.perf_counter() - started, finished=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("all", *CANONICAL_CASE_NAMES), default="figure5")
    parser.add_argument("--max-pieces", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.case == "all" and args.output is not None:
        raise InvalidCurvedCircleDiagnosticLimitError("An explicit output path requires one case, so cases cannot overwrite one another.")
    for name in CANONICAL_CASE_NAMES if args.case == "all" else (args.case,):
        render_native_curved_circles(name, max_pieces=args.max_pieces, output=args.output)


if __name__ == "__main__":
    main()
