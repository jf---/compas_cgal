"""Splice A/B on Monstera's first coverage-repair step, for docs/held_coverage_preserving_selection.md.

The coverage repair splices each candidate source into its boundary-order
slot and re-refines from the slot on (``splice_source``). This tool
measures, in one process, what the exact early stop buys:

- A: continue from the slot to the end (``continue_figure5_engagement``
  from the reference's prefix), the splice without the early stop;
- B: ``continue_with_early_stop``, the splice the repair runs.

The emission, the certificate's witnesses and the candidate set are those of
the repair's first round on the published Monstera selection
(``refine_boundary_ordered_path_certified``'s selection and refinement, the
NOTHING pass, ``path_sweeps``, ``CoverageCertificate.build``, then
``repair_round``'s exact witness-to-candidate matching). For every candidate
both splices run ``TIMING_REPEATS`` times, interleaved; the minimum is
reported. Each B is asserted equal to its A (circles, indices and the
decision footprints of the trace). The record also times one full
refinement of the selection and the repair round itself
(``repair_round``), which is the repair's ``repair_s`` for that round.

Modes: ``--measure`` (about a minute) writes
``benchmarks/results/held_splice_ab.json`` and renders; the default
re-renders the Markdown from that JSON. Both files are generated; do not
edit them.
"""

from __future__ import annotations

import argparse
import bisect
import gc
import json
import math
import os
import pathlib
import time
from collections.abc import Sequence
from dataclasses import asdict
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.held_boundary_order_spacing import select_boundary_ordered_sources
from benchmarks.held_coverage_repair import path_sweeps
from benchmarks.held_coverage_repair import remove_emission_nothing_circles
from benchmarks.held_coverage_repair_cover import candidates_for
from benchmarks.held_coverage_repair_cover import outward_boxes
from benchmarks.held_coverage_repair_cover import repair_round
from benchmarks.held_coverage_repair_splice import continue_with_early_stop
from benchmarks.held_coverage_repair_splice import emitted_form
from benchmarks.held_coverage_repair_types import BoundaryOrderedEmission
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_figure5_engagement_refinement import continue_figure5_engagement
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.performance_record import REPOSITORY
from benchmarks.performance_record import build_identity
from benchmarks.performance_record import machine_identity
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.motion import EngagementCap
from tools.held_native_identity import capture_current_native
from tools.held_native_identity import native_identity_line
from tools.held_native_identity import native_record_fields
from tools.held_source_closure import capture
from tools.held_source_closure import digest_contents
from tools.held_source_closure import source_closure

RESULTS = REPOSITORY / "benchmarks" / "results"
JSON_PATH = RESULTS / "held_splice_ab.json"
MARKDOWN_PATH = RESULTS / "held_splice_ab.md"
CASE = "figure8_monstera"
WORKLOAD = "figure8_monstera_boundary_qualified_corners"
# Repeats per measured call; the minimum is reported, the standard estimator
# of a deterministic computation's cost on a shared machine.
TIMING_REPEATS = 3


class SpliceABContractError(BenchmarkError):
    """The early stop disagrees with the full continuation, the selection is already covered, or the code changed mid-run."""


def _footprints(refinement: TracedEngagementRefinement) -> list[tuple[Any, ...]]:
    return [(footprint.position, footprint.vertex) for footprint in refinement.trace.footprints]


def _min_seconds(seconds: list[float]) -> float:
    return min(seconds)


def measure() -> dict[str, Any]:
    """The first repair round on Monstera: A/B per candidate, a full refinement, the round itself."""
    identity, loaded = build_identity(), capture(source_closure())
    native_run = capture_current_native()
    gc.collect()
    load = os.getloadavg()[0]
    case = load_held_reference_case(CASE)
    components = figure7_inward_offset(case).components
    boundary = components[0]
    source = build_held_reference_path(build_held_reference_raw_guide(case), components, qualify_contacts=True)
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    settings = EngagementRefinementSettings.build(case.tool_radius, cap, corner_approaches=True)
    sources = tuple(emitted_form(circle, boundary) for circle in source.path.circles)
    selected = select_boundary_ordered_sources(sources, source.circle_sources, case.tool_radius, cap)
    picked = tuple(sources[index] for index in selected)
    options = {
        "max_depth": settings.max_depth,
        "max_circles": settings.max_circles,
        "corner_approaches": settings.corner_approaches,
        "close_seam": settings.close_seam,
    }
    full_seconds = []
    for _ in range(TIMING_REPEATS):
        started = time.perf_counter()
        reference = refine_figure5_engagement(picked, boundary, case.tool_radius, cap, **options)
        full_seconds.append(time.perf_counter() - started)
    emission = BoundaryOrderedEmission.build(selected, reference)
    nothing = remove_emission_nothing_circles(emission, settings)
    domain = CentreDomain.build(list(boundary))
    cuts, sweeps = path_sweeps(domain, boundary, nothing.circles, close_seam=settings.close_seam)
    try:
        CoverageCertificate.build(domain, case.tool_radius, cuts, sweeps)
    except UncoveredPocketError as error:
        witnesses = error.witnesses
    else:
        raise SpliceABContractError(f"{WORKLOAD}: the selection is covered; there is no repair step to measure")
    chosen = set(selected)
    dropped = tuple(index for index in range(len(sources)) if index not in chosen)
    forms = {index: emitted_form(sources[index], boundary) for index in dropped}
    dropped_forms = tuple(forms[index] for index in dropped)
    boxes = outward_boxes(dropped_forms, case.tool_radius)
    candidates = sorted({candidate for witness in witnesses for candidate in candidates_for(witness, dropped, dropped_forms, boxes, case.tool_radius)})
    rows = []
    for candidate in candidates:
        slot = bisect.bisect_left(selected, candidate)
        later = tuple(sources[index] for index in selected[slot:])
        start = reference.original_indices[slot - 1]
        prefix = TracedEngagementRefinement(reference.circles[: start + 1], reference.original_indices[:slot], reference.trace)
        seconds_a, seconds_b = [], []
        for _ in range(TIMING_REPEATS):
            started = time.perf_counter()
            full = continue_figure5_engagement(prefix, (sources[candidate], *later), boundary, case.tool_radius, cap, **options)
            seconds_a.append(time.perf_counter() - started)
            started = time.perf_counter()
            stopped = continue_with_early_stop(reference, slot, sources[candidate], later, boundary, settings)
            seconds_b.append(time.perf_counter() - started)
        if stopped.refinement != full or _footprints(stopped.refinement) != _footprints(full):
            raise SpliceABContractError(f"{WORKLOAD}: candidate {candidate}: the early stop differs from the full continuation")
        rows.append(
            {
                "candidate": candidate,
                "slot": slot,
                "later_sources": len(later),
                "rerefined": stopped.rerefined,
                "added": len(full.circles) - len(reference.circles),
                "a_seconds": _min_seconds(seconds_a),
                "b_seconds": _min_seconds(seconds_b),
            }
        )
    started = time.perf_counter()
    repaired, inserted = repair_round(sources, boundary, settings, emission, witnesses, dict(forms), 1, [])
    round_seconds = time.perf_counter() - started
    if capture(loaded) != loaded:
        raise SpliceABContractError("the measured sources changed during the run; re-measure once edits settle")
    closure = source_closure()
    return {
        "build": asdict(identity),
        "source_sha256": digest_contents(capture(closure)),
        "sources": list(closure),
        **native_record_fields(native_run),
        "machine": asdict(machine_identity()),
        "workload": WORKLOAD,
        "sources_count": len(sources),
        "selected": len(selected),
        "emitted": len(reference.circles),
        "emitted_after_nothing": len(nothing.circles),
        "witnesses": len(witnesses),
        "full_refinement_seconds": _min_seconds(full_seconds),
        "candidates": rows,
        "repair_round_seconds": round_seconds,
        "repair_round_inserted": [insertion.source for insertion in inserted],
        "repair_round_emitted_after": len(repaired.refinement.circles),
        "load_1min_before_after": [load, os.getloadavg()[0]],
    }


def render(record: dict[str, Any]) -> str:
    """The Markdown fragment the docs page includes."""
    build, rows = record["build"], record["candidates"]
    sum_a, sum_b = sum(row["a_seconds"] for row in rows), sum(row["b_seconds"] for row in rows)
    later = [row["later_sources"] for row in rows]
    rerefined = [row["rerefined"] for row in rows]
    a = [row["a_seconds"] for row in rows]
    b = [row["b_seconds"] for row in rows]
    load = record["load_1min_before_after"]
    lines = [
        "<!-- Generated by `pixi run held-splice-ab` (tools/held_splice_ab_results.py --measure); do not edit. -->",
        "",
        f"Recorded {build['recorded_at']} at commit `{build['commit'][:10]}`{' (dirty tree)' if build['dirty'] else ''}, "
        f"measured sources sha256 `{record['source_sha256'][:12]}`, {native_identity_line(record)}, on {record['machine']['cpu_brand']}, load {load[0]:.0f}–{load[1]:.0f}.",
        "",
        f"`{record['workload']}`: {record['sources_count']:,} sources, {record['selected']:,} selected, {record['emitted']:,} emitted "
        f"({record['emitted_after_nothing']:,} after the NOTHING pass), {record['witnesses']} witnesses, {len(rows)} candidates in the first greedy step. "
        f"A: continue from the slot to the end. B: exact early stop. Each B equals its A (asserted). Minimum of {TIMING_REPEATS} runs.",
        "",
        "| | A: continue to the end | B: exact early stop |",
        "|---|---|---|",
        f"| later sources re-refined | {min(later)}–{max(later)} | {min(rerefined)}–{max(rerefined)} |",
        f"| per candidate | {min(a):.3f}–{max(a):.3f} s | {1000 * min(b):.1f}–{1000 * max(b):.1f} ms |",
        f"| sum over {len(rows)} candidates | {sum_a:.2f} s | {sum_b:.3f} s ({sum_a / sum_b:.0f}×) |",
        "",
        f"Full refinement of the selection: {record['full_refinement_seconds']:.3f} s. Repair round 1 (`repair_round`, the round's `repair_s`): "
        f"{record['repair_round_seconds']:.3f} s, inserted sources {record['repair_round_inserted']}, "
        f"{record['emitted']:,} → {record['repair_round_emitted_after']:,} emitted circles before the NOTHING pass.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--measure", action="store_true", help="re-measure (about a minute) before rendering")
    parser.add_argument("--json", type=pathlib.Path, default=JSON_PATH)
    parser.add_argument("--markdown", type=pathlib.Path, default=MARKDOWN_PATH)
    args = parser.parse_args(argv)
    if args.measure:
        args.json.write_text(json.dumps(measure(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    record = json.loads(args.json.read_text(encoding="utf-8"))
    args.markdown.write_text(render(record), encoding="utf-8")
    print(args.markdown.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
