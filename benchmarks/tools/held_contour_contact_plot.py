"""Inspect native contour contact corrections on an existing circle sequence."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_reference_cases import HeldReferenceCase
from compas_cgal import _stock_2

CONTACT_ARC_DISPLAY_SAMPLES = 64  # Annotation polyline only; not geometric decisions.


class InvalidContourPlotPrefixError(ValueError):
    """A contact plot needs at least one predecessor/successor pair."""


def _xy(circle: Figure5CounterclockwiseCircle) -> tuple[float, float]:
    return float(circle.center.x), float(circle.center.y)


def render_contour_contacts(case: HeldReferenceCase, circles: tuple[Figure5CounterclockwiseCircle, ...], output: Path) -> None:
    """Plot exact CW correction decisions; displayed coordinates are approximate.

    Supplied circles remain unchanged. Fully covered or concentric predecessors
    are counted explicitly; neither implies zero engagement for the successor.
    """
    if len(circles) < 2:
        raise InvalidContourPlotPrefixError("Contour contact plot requires at least two circles.")
    tool = float(case.tool_radius.value)
    contour = _stock_2.HeldDiskContour2(_xy(circles[0]), float(circles[0].radius.value), tool)
    corrections: list[tuple[int, tuple[float, float], tuple[float, float]]] = []
    covered = 0
    concentric = 0
    unchanged = 0
    for index, (previous, candidate) in enumerate(zip(circles, circles[1:])):
        try:
            point, moved = contour.contact_toward(_xy(candidate))
        except _stock_2.NoExposedPredecessorArcError:
            covered += 1
        except _stock_2.UndefinedPredecessorDirectionError:
            concentric += 1
        else:
            if moved:
                center = np.array(_xy(previous))
                direction = np.array(_xy(candidate)) - center
                # Reporting only: the original b is shown next to CGAL's exact
                # correction. These display coordinates never drive placement.
                standard = center + direction / np.linalg.norm(direction) * (float(previous.radius.value) + tool)
                corrections.append((index, (float(standard[0]), float(standard[1])), point))
            else:
                unchanged += 1
        contour.append(_xy(candidate), float(candidate.radius.value))

    boundary = np.array([(float(point.x), float(point.y)) for point in case.projection.points])
    boundary = np.vstack((boundary, boundary[0]))
    figure = Figure(figsize=(12, 6), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    overview, detail = figure.subplots(1, 2)
    for axis in (overview, detail):
        axis.plot(*boundary.T, color="#34403f", linewidth=1)
        axis.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    for circle in circles:
        overview.add_patch(Circle(_xy(circle), float(circle.radius.value), fill=False, edgecolor="#087e8b", linewidth=0.5))
    for _, original, corrected in corrections:
        overview.annotate("", xy=corrected, xytext=original, arrowprops={"arrowstyle": "->", "color": "#c34c32", "lw": 0.7})
    overview.set_title(f"{len(circles):,} existing circles · {len(corrections):,} contact corrections")

    if corrections:
        index, original, corrected = max(corrections, key=lambda row: math.dist(row[1], row[2]))
        previous, candidate = circles[index : index + 2]
        for earlier in circles[:index]:
            detail.add_patch(Circle(_xy(earlier), float(earlier.radius.value) + tool, facecolor="#91d0ce", edgecolor="none", alpha=0.15))
        outer = float(previous.radius.value) + tool
        detail.add_patch(Circle(_xy(previous), outer, fill=False, edgecolor="#087e8b", linewidth=2, label="Predecessor outer circle"))
        detail.add_patch(Circle(_xy(candidate), float(candidate.radius.value), fill=False, edgecolor="#34403f", linestyle="--", label="Candidate machining circle"))
        detail.scatter(*original, color="#c34c32", marker="x", s=65, label="Standard b (covered)")
        detail.scatter(*corrected, color="#6d3b88", s=35, label="Corrected b (exposed)")
        cx, cy = _xy(previous)
        start = math.atan2(original[1] - cy, original[0] - cx)
        end = math.atan2(corrected[1] - cy, corrected[0] - cx)
        angles = np.linspace(start, start - (start - end) % math.tau, CONTACT_ARC_DISPLAY_SAMPLES)
        arc_x, arc_y = cx + outer * np.cos(angles), cy + outer * np.sin(angles)
        detail.plot(arc_x, arc_y, color="#6d3b88", linewidth=2, linestyle=":")
        detail.annotate("", xy=corrected, xytext=(arc_x[-2], arc_y[-2]), arrowprops={"arrowstyle": "->", "color": "#6d3b88"})
        extent = outer + tool
        detail.set_xlim(cx - extent, cx + extent)
        detail.set_ylim(cy - extent, cy + extent)
        detail.set_title(f"Largest displayed correction · predecessor {index + 1}")
        figure.legend(*detail.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.065, 0.93), ncol=4, fontsize=8)
    else:
        detail.text(0.5, 0.5, "No contact corrections in this prefix", transform=detail.transAxes, ha="center")
        detail.set_title("No corrected contact witness")
    figure.suptitle(f"{case.name.replace('_', ' ')} · Held Section 3.1 contact correction", x=0.07, ha="left", fontsize=16)
    figure.text(0.07, 0.065, f"{len(corrections)} corrected · {unchanged} unchanged · {covered} fully covered predecessors · {concentric} concentric directions", fontsize=10)
    figure.text(0.07, 0.03, "Exact CGAL filled-disk contour query on existing standard circles. Placement and engagement are not yet contour-aware.", fontsize=9)
    figure.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.17, wspace=0.25)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    output.with_suffix(".json").write_text(
        json.dumps(
            {
                "case": case.name,
                "circle_count": len(circles),
                "corrected": len(corrections),
                "unchanged": unchanged,
                "fully_covered": covered,
                "concentric": concentric,
                "corrections": corrections,
                "coordinate_units": "mm",
                "coordinate_semantics": "approximate reports of exact native contact decisions",
                "scope": "contact query only; existing standard circles unchanged",
            },
            indent=2,
        )
        + "\n"
    )
