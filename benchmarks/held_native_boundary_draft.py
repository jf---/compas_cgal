"""Bounded native offset-boundary traversal diagnostic; no machining qualification."""

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
from matplotlib.patches import Circle

from benchmarks.errors import BenchmarkError
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import _circle_geometry_2
from compas_cgal import _coverage_2

MAX_SAMPLES_PER_PRIMITIVE = 16  # Keep the construction diagnostic bounded.
ARC_DISPLAY_INTERVALS = 32  # Display tessellation only; native transitions stay opaque.
CHECKPOINT_SECONDS = 30.0  # User-facing visual feedback cadence, not geometry.
MINIMUM_DISPLAY_SPAN_MM = 1.0  # Keep stationary-only checkpoint zooms visible.


class InvalidBoundaryDraftSamplingError(BenchmarkError):
    """Diagnostic sampling must use a bounded power-of-two count."""


class DisconnectedNativeBoundaryDraftError(BenchmarkError):
    """Native transition endpoints do not match their original contacts."""


def _draw_primitive(axes: Axes, primitive: _coverage_2.ReachableBoundaryPrimitive2, *, color: str, linewidth: float) -> None:
    """Project a native primitive to plotting coordinates, never back to geometry."""
    first, last = primitive.start_mm, primitive.end_mm
    if primitive.kind == "line":
        axes.plot((first[0], last[0]), (first[1], last[1]), color=color, linewidth=linewidth)
        return
    center, radius = primitive.arc_center_mm, primitive.arc_radius_mm
    start = math.atan2(first[1] - center[1], first[0] - center[0])
    end = math.atan2(last[1] - center[1], last[0] - center[0])
    sweep = (end - start) % math.tau if primitive.arc_counterclockwise else -((start - end) % math.tau)
    angles = np.linspace(start, start + sweep, ARC_DISPLAY_INTERVALS + 1)
    axes.plot(center[0] + radius * np.cos(angles), center[1] + radius * np.sin(angles), color=color, linewidth=linewidth)


def _write_checkpoint(
    case_name: str,
    boundary: list[tuple[float, float]],
    primitives: list[_coverage_2.ReachableBoundaryPrimitive2],
    proposals: list[_circle_geometry_2.BoundaryNormalCircleProposal2],
    completed: int,
    elapsed: float,
    output: Path,
) -> None:
    """Publish actual completed construction while later native calls still run."""
    figure = Figure(figsize=(13, 7), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    overview, zoom = figure.subplots(1, 2)
    ring = np.array(boundary + boundary[:1])
    for axes in (overview, zoom):
        axes.plot(ring[:, 0], ring[:, 1], color="#34403f", linewidth=1)
        for primitive in primitives:
            _draw_primitive(axes, primitive, color="#bd873b", linewidth=1)
        for proposal in proposals:
            if proposal.is_stationary:
                axes.scatter(*proposal.center_mm, color="#c0392b", marker="x", s=25)
            else:
                axes.add_patch(Circle(proposal.center_mm, proposal.guide_radius_mm, fill=False, edgecolor="#087e8b", linewidth=0.8))
        axes.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    if proposals:
        extents = np.array([(p.center_mm[0] + sign * p.guide_radius_mm, p.center_mm[1] + sign * p.guide_radius_mm) for p in proposals for sign in (-1, 1)])
        low, high = extents.min(axis=0), extents.max(axis=0)
        margin = max(float(np.max(high - low)), MINIMUM_DISPLAY_SPAN_MM) / 4
        zoom.set(xlim=(low[0] - margin, high[0] + margin), ylim=(low[1] - margin, high[1] + margin))
    overview.set_title("Whole pocket")
    zoom.set_title("Completed native circles · red crosses are stationary")
    figure.suptitle(f"{case_name} · PARTIAL CONSTRUCTION — NOT A PATH", fontsize=14)
    figure.text(0.1, 0.9, f"{completed}/{len(primitives)} primitives completed · {len(proposals)} native events · {elapsed:.1f} s", fontsize=11)
    figure.text(0.1, 0.025, "Native circles shown as completed. Transitions, engagement and coverage are unqualified.", fontsize=9)
    figure.subplots_adjust(top=0.84, bottom=0.13)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    report = {
        "case": case_name,
        "complete_construction": False,
        "status": "construction checkpoint; not a path",
        "sampled_primitives": completed,
        "offset_primitives": len(primitives),
        "sampled_events": len(proposals),
        "stationary_events": sum(proposal.is_stationary for proposal in proposals),
        "elapsed_seconds": elapsed,
        "events": [
            {"stationary": p.is_stationary, "center_mm": p.center_mm, "guide_radius_mm": p.guide_radius_mm, "p_mm": p.p_mm, "m_mm": p.m_mm, "q_mm": p.q_mm} for p in proposals
        ],
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"PLOT {output}: {completed}/{len(primitives)} primitives completed", flush=True)


def render_native_boundary_draft(case_name: str, samples_per_primitive: int = 1, output: Path | None = None, *, max_primitives: int | None = None) -> None:
    """Construct native circles and closed CCW transitions at dyadic contacts.

    Sampling controls this diagnostic only. Every contact remains native through
    construction and transition lookup. Reported coordinates are used only in
    PNG/JSON output. Native construction errors propagate; sites are not repaired
    or silently omitted. Stationary events retain their native classification.
    """
    if samples_per_primitive < 1 or samples_per_primitive > MAX_SAMPLES_PER_PRIMITIVE or samples_per_primitive & (samples_per_primitive - 1):
        raise InvalidBoundaryDraftSamplingError(f"samples_per_primitive must be a power of two from 1 to {MAX_SAMPLES_PER_PRIMITIVE}.")
    if max_primitives is not None and max_primitives < 1:
        raise InvalidBoundaryDraftSamplingError("max_primitives must be positive when supplied.")
    case = load_held_reference_case(case_name)
    output = output or Path(f"docs/assets/images/held_{case.name}_native_boundary_draft.png")
    boundary = [(float(point.x), float(point.y)) for point in case.projection.points]
    radius = float(case.tool_radius.value)
    started = time.perf_counter()
    primitives: list[_coverage_2.ReachableBoundaryPrimitive2] = []
    contacts: list[_coverage_2.WorldXYBoundaryPointMm] = []
    proposals: list[_circle_geometry_2.BoundaryNormalCircleProposal2] = []
    primitive_indices: list[int] = []
    transition_pieces: list[_coverage_2.ReachableBoundaryPrimitive2] = []
    failure: Exception | None = None
    sampled_primitives = 0
    bounded_prefix = False
    last_checkpoint = started
    try:
        cycle = _coverage_2.build_center_boundary_cycle(np.array([(x, y, 0.0) for x, y in boundary], dtype=np.float64), [], radius)
        primitives = list(cycle.primitives)
        owner = _circle_geometry_2.BoundaryNormalCircle2(boundary)
        # Include junctions: convex stationary events are absent from interiors.
        # Dyadic midpoint parameters are exactly representable at the native seam.
        bounded_prefix = max_primitives is not None and max_primitives < len(primitives)
        for index, primitive in enumerate(primitives[:max_primitives]):
            sampled = [primitive.start] + [primitive.sample((2 * ordinal + 1) / (2 * samples_per_primitive)) for ordinal in range(samples_per_primitive)]
            for contact in sampled:
                if contacts and (contact == contacts[-1] or contact == contacts[0]):
                    continue
                event_started = time.perf_counter()
                proposal = _coverage_2.boundary_circle_at_contact(owner, contact, radius)
                print(
                    f"{case.name}: event {len(proposals) + 1}, primitive {index + 1}, "
                    f"stationary={proposal.is_stationary}, construction {time.perf_counter() - event_started:.3f} s",
                    flush=True,
                )
                contacts.append(contact)
                proposals.append(proposal)
                primitive_indices.append(index)
            sampled_primitives += 1
            print(f"{case.name}: primitive {index + 1}/{len(primitives)} complete; {len(proposals)} events", flush=True)
            now = time.perf_counter()
            if sampled_primitives <= 3 or now - last_checkpoint >= CHECKPOINT_SECONDS:
                _write_checkpoint(case.name, boundary, primitives, proposals, sampled_primitives, now - started, output)
                last_checkpoint = now
        if not contacts:
            raise DisconnectedNativeBoundaryDraftError("Native offset cycle supplied no sampling contacts.")
        for first, last in zip(contacts, contacts[1:] if bounded_prefix else contacts[1:] + contacts[:1]):
            pieces = cycle.ccw_transition(first, last)
            if not pieces or pieces[0].start != first or pieces[-1].end != last:
                raise DisconnectedNativeBoundaryDraftError("Native CCW transition does not join its unchanged contacts.")
            transition_pieces.extend(pieces)
    except (
        _coverage_2.ReachableDomainConstructionError,
        _coverage_2.ReachableArrangementTopologyError,
        _coverage_2.CoverageTransitionError,
        _coverage_2.InvalidCoverageGeometryError,
        _coverage_2.InvalidBoundaryCircleContactError,
        _coverage_2.BoundaryContactConstructionError,
        _circle_geometry_2.InvalidBoundaryPolygonError,
        _circle_geometry_2.InvalidBoundaryNormalInputError,
        _circle_geometry_2.NoPositiveBoundaryCircleError,
        _circle_geometry_2.BoundaryNormalConstructionError,
        DisconnectedNativeBoundaryDraftError,
    ) as error:
        # Preserve completed construction evidence, then propagate the same error.
        failure = error
    elapsed = time.perf_counter() - started
    stationary = sum(proposal.is_stationary for proposal in proposals)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure = Figure(figsize=(13, 7), facecolor="#faf9f5")
    FigureCanvasAgg(figure)
    motion, medial = figure.subplots(1, 2)
    ring = np.array(boundary + boundary[:1])
    for axes in (motion, medial):
        axes.plot(ring[:, 0], ring[:, 1], color="#34403f", linewidth=1)
        for primitive in primitives:
            _draw_primitive(axes, primitive, color="#bd873b", linewidth=1)
        axes.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    for primitive in transition_pieces:
        _draw_primitive(motion, primitive, color="#bd873b", linewidth=1.2)
    for proposal in proposals:
        if proposal.is_stationary:
            motion.scatter(*proposal.center_mm, color="#c0392b", marker="x", s=25)
            medial.scatter(*proposal.m_mm, color="#c0392b", marker="x", s=25)
        else:
            motion.add_patch(Circle(proposal.center_mm, proposal.guide_radius_mm, fill=False, edgecolor="#087e8b", linewidth=0.6, alpha=0.8))
            medial.scatter(*proposal.m_mm, color="#087e8b", s=9)
    motion.plot([], [], color="#087e8b", label="Native machining circles")
    motion.plot([], [], color="#bd873b", label="Exact offset / CCW connectors")
    motion.scatter([], [], color="#c0392b", marker="x", label="Stationary events")
    motion.legend(fontsize=8)
    motion.set_title("PARTIAL CONSTRUCTION — NOT A PATH" if failure or bounded_prefix else "Sampled circles and unchanged boundary transitions")
    medial.set_title("Native medial contacts; red crosses are stationary events")
    figure.suptitle(f"{case.name.replace('_', ' ')} · native boundary traversal", x=0.065, ha="left", fontsize=18)
    figure.text(
        0.065,
        0.90,
        f"{sampled_primitives}/{len(primitives)} sampled primitives · {len(proposals)} sampled events · {stationary} stationary · {elapsed:.2f} s native construction",
        fontsize=11,
    )
    figure.text(0.065, 0.035, "Construction diagnostic only: no engagement-cap, complete-coverage, entry or emitted-motion qualification.", fontsize=10)
    if failure is not None:
        figure.text(0.065, 0.075, f"FAILED: {type(failure).__name__}: {str(failure)[:140]}", fontsize=9, color="#c0392b")
    figure.subplots_adjust(left=0.065, right=0.975, top=0.83, bottom=0.14, wspace=0.2)
    figure.savefig(output, dpi=180)
    report = {
        "case": case.name,
        "algorithm": "native_boundary_traversal_diagnostic",
        "complete_construction": failure is None and not bounded_prefix,
        "sampled_primitives": sampled_primitives,
        "bounded_prefix": bounded_prefix,
        "failure": None if failure is None else {"type": type(failure).__name__, "message": str(failure)},
        "samples_per_primitive": samples_per_primitive,
        "offset_primitives": len(primitives),
        "sampled_events": len(proposals),
        "stationary_events": stationary,
        "transition_pieces": len(transition_pieces),
        "native_construction_seconds": elapsed,
        "scope": "Opaque native contacts and CCW transitions; all coordinates reporting only. No engagement-cap, complete-coverage, entry or emitted-motion qualification.",
        "events": [
            {
                "primitive_index": index,
                "stationary": proposal.is_stationary,
                "p_mm": proposal.p_mm,
                "m_mm": proposal.m_mm,
                "q_mm": proposal.q_mm,
                "center_mm": proposal.center_mm,
                "guide_radius_mm": proposal.guide_radius_mm,
                "clearance_mm": proposal.clearance_mm,
            }
            for index, proposal in zip(primitive_indices, proposals)
        ],
    }
    output.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if failure is not None:
        raise failure
    print(f"{case.name}: {len(proposals)} events, {stationary} stationary, {elapsed:.2f} s; {output}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CANONICAL_CASE_NAMES, default="figure5")
    parser.add_argument("--samples-per-primitive", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-primitives", type=int)
    args = parser.parse_args()
    render_native_boundary_draft(args.case, args.samples_per_primitive, args.output, max_primitives=args.max_primitives)


if __name__ == "__main__":
    main()
