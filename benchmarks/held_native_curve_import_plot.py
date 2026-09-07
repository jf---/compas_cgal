"""Compare prepared Held curves with endpoint-preserving native imports."""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from benchmarks.held_native_boundary_draft import _draw_primitive
from benchmarks.held_native_curve_import import import_held_boundary
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceLine
from compas_cgal import _coverage_2 as native

SOURCE_DISPLAY_INTERVALS = 32  # Plot tessellation only; never native input.


def _draw_source(axes: Axes, case: HeldReferenceCase) -> None:
    for primitive in case.boundary.primitives:
        start = float(primitive.start.x), float(primitive.start.y)
        end = float(primitive.end.x), float(primitive.end.y)
        if isinstance(primitive, ReferenceLine):
            axes.plot((start[0], end[0]), (start[1], end[1]), color="#c08a3e", linewidth=2.5)
        else:
            center = float(primitive.centre.x), float(primitive.centre.y)
            radius = math.dist(start, center)
            angle = math.atan2(start[1] - center[1], start[0] - center[0])
            angles = np.linspace(angle, angle + float(primitive.sweep), SOURCE_DISPLAY_INTERVALS + 1)
            axes.plot(center[0] + radius * np.cos(angles), center[1] + radius * np.sin(angles), color="#c08a3e", linewidth=2.5)


def render_native_curve_import(case_name: str, output: Path | None = None) -> None:
    """Render successful import or named failure, then propagate any failure."""
    case = load_held_reference_case(case_name)
    output = output or Path(f"docs/assets/images/held_{case.name}_native_curve_import.png")
    owner: native.NativeBoundary2 | None = None
    failure: RuntimeError | None = None
    started = time.perf_counter()
    try:
        owner = import_held_boundary(case)
    except (native.InvalidNativeBoundaryCurveError, native.InvalidNativeBoundaryChainError, native.ReachableDomainConstructionError) as error:
        failure = error
    elapsed = time.perf_counter() - started
    figure = Figure(figsize=(13, 6), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    overlay, adjustments = figure.subplots(1, 2)
    _draw_source(overlay, case)
    overlay.plot([], [], color="#c08a3e", linewidth=2.5, label="Prepared source curves")
    arc_indices: list[int] = []
    shifts: list[float] = []
    if owner is not None:
        for primitive in owner.cycle.primitives:
            _draw_primitive(overlay, primitive, color="#087e8b", linewidth=0.9)
        overlay.plot([], [], color="#087e8b", linewidth=0.9, label="Imported native cycle")
        for index, curve in enumerate(owner.curves):
            if curve.is_arc:
                arc_indices.append(index)
                shifts.append(curve.center_adjustment_mm)
        adjustments.plot(arc_indices, shifts, color="#087e8b", marker=".", linewidth=0.7)
        adjustments.set_title(f"Maximum centre adjustment {max(shifts, default=0):.4g} mm")
    else:
        overlay.set_title("FAILED IMPORT · prepared source shown only")
        adjustments.text(0.5, 0.5, f"{type(failure).__name__}\n{failure}", ha="center", va="center", wrap=True, transform=adjustments.transAxes, color="#c0392b")
    overlay.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    overlay.legend(fontsize=9)
    adjustments.set(xlabel="Arc index in prepared input chain", ylabel="Centre adjustment / mm")
    adjustments.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    adjustments.grid(alpha=0.2)
    figure.suptitle(f"{case.name.replace('_', ' ')} · native line/arc import", x=0.075, ha="left", fontsize=18)
    imported = 0 if owner is None else len(owner.curves)
    figure.text(0.075, 0.90, f"{len(case.boundary.primitives)} prepared curves · {imported} imported · {len(shifts)} arcs · {elapsed:.3f} s import", fontsize=11)
    figure.text(0.075, 0.03, "Import diagnostic only: no G1 continuity, source-fit, machinability or toolpath acceptance claim.", fontsize=10)
    figure.subplots_adjust(left=0.075, right=0.98, top=0.81, bottom=0.15, wspace=0.3)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    report = {
        "case": case.name,
        "complete_import": failure is None,
        "failure": None if failure is None else {"type": type(failure).__name__, "message": str(failure)},
        "prepared_curves": len(case.boundary.primitives),
        "imported_curves": imported,
        "native_cycle_primitives": 0 if owner is None else len(owner.cycle.primitives),
        "import_seconds": elapsed,
        "arc_center_adjustments": [{"input_index": index, "center_adjustment_mm": shift} for index, shift in zip(arc_indices, shifts)],
        "maximum_center_adjustment_mm": max(shifts, default=0) if owner is not None else None,
        "scope": (
            "Endpoint-preserving native import. Center displacement is reporting only, not a trimmed-source error bound. "
            "No G1, source-fit, machinability or toolpath acceptance claim."
        ),
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"PLOT {output}: {imported} imported curves in {elapsed:.3f} s; complete={failure is None}", flush=True)
    if failure is not None:
        raise failure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("all", *CANONICAL_CASE_NAMES), default="all")
    args = parser.parse_args()
    failures: list[RuntimeError] = []
    for case_name in CANONICAL_CASE_NAMES if args.case == "all" else (args.case,):
        try:
            render_native_curve_import(case_name)
        except (native.InvalidNativeBoundaryCurveError, native.InvalidNativeBoundaryChainError, native.ReachableDomainConstructionError) as error:
            failures.append(error)
    if failures:
        raise failures[0]


if __name__ == "__main__":
    main()
