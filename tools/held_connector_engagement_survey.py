"""Survey the connector-engagement certificate (tiers 1 and 1b) on the five Held workloads.

Generates every published path live (``benchmarks.held_toolpath_benchmark``),
runs ``ConnectorEngagementAssessment`` on it, and writes the per-workload table
to ``benchmarks/results/held_connector_engagement.{json,md}``; the docs page
includes the markdown. The tiers' own clock is reported against
``TIERS_BUDGET_S`` as the best of ``TIMING_REPEATS`` runs (generation is
reported beside it, not gated here). Every count is an exact decision; every
degree value is reporting only.
"""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Sequence
from dataclasses import asdict
from dataclasses import dataclass

from benchmarks.errors import BenchmarkError
from benchmarks.held_connector_engagement import TIERS_BUDGET_S
from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import ConnectorEngagementLemma
from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_connector_engagement import ConnectorVerdict
from benchmarks.held_connector_engagement import LemmaCircle
from benchmarks.held_connector_engagement import foot_squared_distance
from benchmarks.held_connector_engagement import moving_pieces
from benchmarks.held_toolpath_benchmark import GENERATORS
from benchmarks.held_toolpath_benchmark import GeneratedToolpath
from benchmarks.held_toolpath_benchmark import generate
from benchmarks.performance_record import REPOSITORY
from benchmarks.performance_record import build_identity
from benchmarks.performance_record import machine_identity
from compas_cgal.motion import EngagementCap

OUTPUT_MD = REPOSITORY / "benchmarks/results/held_connector_engagement.md"
OUTPUT_JSON = REPOSITORY / "benchmarks/results/held_connector_engagement.json"
# Repeats of the tiers' clock per workload; the minimum is reported, the
# standard estimator of a deterministic computation's cost on a shared machine.
TIMING_REPEATS = 3


class NondeterministicCertificateError(BenchmarkError):
    """Repeated assessments of one generated path disagree."""


@dataclass(frozen=True)
class SeamSummary:
    tier: str
    tier1_bound_deg: float
    circles_used: tuple[int, ...]


@dataclass(frozen=True)
class SurveyRow:
    """One workload's certificate counts; degrees are reporting only."""

    workload: str
    cap_deg: float
    circles: int
    connectors: int
    moving: int
    pieces: int
    tier1_certified: int
    tier1_refused: int
    tier1b_certified: int
    left: tuple[int, ...]
    hole_arc_connectors: int
    tier1_refused_bound_deg: tuple[float, ...]
    seam: SeamSummary | None
    tiers_seconds: float
    generation_seconds: float


def _tier1_bound_deg(path: GeneratedToolpath, verdict: ConnectorVerdict, lemma: ConnectorEngagementLemma) -> float:
    circle = LemmaCircle.from_circle(path.circles[verdict.connector])
    return max(lemma.bound_deg(a, b, circle) for a, b in moving_pieces(path.transitions[verdict.connector]))


def _has_hole_arc(path: GeneratedToolpath, verdict: ConnectorVerdict, lemma: ConnectorEngagementLemma) -> bool:
    """REPORTING: some piece passes closer than rho_i to circle i's centre, with rho_i > r (exact)."""
    circle = LemmaCircle.from_circle(path.circles[verdict.connector])
    if circle.guide_radius <= lemma.tool_radius:
        return False
    rho_sq = circle.guide_radius * circle.guide_radius
    return any(foot_squared_distance(a, b, circle.centre) < rho_sq for a, b in moving_pieces(path.transitions[verdict.connector]))


def survey(workload: str) -> SurveyRow:
    path = generate(workload)
    cap = EngagementCap.build(math.radians(float(path.case.tea_cap)))
    tool = path.case.tool_radius
    runs = [ConnectorEngagementAssessment.build(path.circles, path.connectors, tool, cap) for _ in range(TIMING_REPEATS)]
    assessment = runs[0]
    if any(run.verdicts != assessment.verdicts for run in runs[1:]):
        raise NondeterministicCertificateError(f"{workload}: repeated assessments disagree")
    lemma = ConnectorEngagementLemma.build(tool, cap)
    moving = [verdict for verdict in assessment.verdicts if verdict.tier is not ConnectorTier.STATIONARY]
    tier1_refused = [verdict for verdict in moving if verdict.tier is not ConnectorTier.TIER_1]
    seam_verdict = next((verdict for verdict in assessment.verdicts if verdict.seam), None)
    seam = None
    if seam_verdict is not None and seam_verdict.tier is not ConnectorTier.STATIONARY:
        seam = SeamSummary(seam_verdict.tier.value, _tier1_bound_deg(path, seam_verdict, lemma), tuple(sorted({sub.circle for sub in seam_verdict.sub_pieces})))
    return SurveyRow(
        workload=workload,
        cap_deg=float(path.case.tea_cap),
        circles=len(path.circles),
        connectors=len(path.transitions),
        moving=len(moving),
        pieces=sum(verdict.pieces for verdict in moving),
        tier1_certified=assessment.count(ConnectorTier.TIER_1),
        tier1_refused=len(tier1_refused),
        tier1b_certified=assessment.count(ConnectorTier.TIER_1B),
        left=assessment.refused,
        hole_arc_connectors=sum(1 for verdict in moving if _has_hole_arc(path, verdict, lemma)),
        tier1_refused_bound_deg=tuple(sorted(_tier1_bound_deg(path, verdict, lemma) for verdict in tier1_refused)),
        seam=seam,
        tiers_seconds=min(run.seconds for run in runs),
        generation_seconds=path.generation_seconds,
    )


def _seam_cell(seam: SeamSummary | None) -> str:
    if seam is None:
        return "none (coincident)"
    if seam.tier == ConnectorTier.TIER_1.value:
        return f"tier 1 ({seam.tier1_bound_deg:.2f}°)"
    if seam.tier == ConnectorTier.TIER_1B.value:
        return f"tier 1 refuses (bound {seam.tier1_bound_deg:.2f}°); tier 1b certifies with circles {', '.join(map(str, seam.circles_used))}"
    return f"refused (tier-1 bound {seam.tier1_bound_deg:.2f}°)"


def _bounds_cell(bounds: tuple[float, ...]) -> str:
    if not bounds:
        return "—"
    if len(bounds) == 1:
        return f"{bounds[0]:.2f}"
    return f"{bounds[0]:.2f}–{bounds[-1]:.2f} (median {bounds[len(bounds) // 2]:.2f})"


def markdown(rows: Sequence[SurveyRow]) -> str:
    """The generated results fragment the docs page includes."""
    build, machine = build_identity(), machine_identity()
    over = [row.workload for row in rows if row.tiers_seconds > TIERS_BUDGET_S]
    lines = [
        "<!-- Generated by `pixi run held-connector-engagement` (tools/held_connector_engagement_survey.py); do not edit. -->",
        "",
        f"Commit `{build.commit[:10]}`{' (dirty)' if build.dirty else ''} · {machine.cpu_brand} · recorded {build.recorded_at}.",
        "",
        "| workload | circles | connectors | moving | pieces | tier 1 certified | tier 1b certifies | left | seam "
        "| hole-arc connectors | tier-1 bound at refusals (°) | tiers 1 + 1b s | generation s |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| `{row.workload}` | {row.circles:,} | {row.connectors:,} | {row.moving:,} | {row.pieces:,} | {row.tier1_certified:,} "
            f"| {row.tier1b_certified} of {row.tier1_refused} | **{len(row.left)}** | {_seam_cell(row.seam)} | {row.hole_arc_connectors:,} "
            f"| {_bounds_cell(row.tier1_refused_bound_deg)} | {row.tiers_seconds:.3f} | {row.generation_seconds:.1f} |"
        )
    lines += [
        "",
        f"Cap {', '.join(sorted({f'{row.cap_deg:g}°' for row in rows}))}. Tiers' clock: best of {TIMING_REPEATS}, decisions only (path already generated); "
        f"budget {TIERS_BUDGET_S} s per workload: {'within budget on every workload' if not over else 'OVER BUDGET on ' + ', '.join(over)}.",
        f"Left without a certificate: {sum(len(row.left) for row in rows)} connectors in all.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("workloads", nargs="*", default=list(GENERATORS), help="default: all five")
    args = parser.parse_args(argv)
    rows = [survey(workload) for workload in args.workloads]
    OUTPUT_JSON.write_text(json.dumps([asdict(row) for row in rows], indent=2) + "\n")
    OUTPUT_MD.write_text(markdown(rows))
    print(OUTPUT_MD.read_text())


if __name__ == "__main__":
    main()
