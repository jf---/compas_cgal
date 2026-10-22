"""Slow-tier acceptance: select, certify, repair closes both Monstera slivers.

Not collected by pytest (benchmarks/acceptance is outside ``testpaths``). Run with
the editable-skip invocation and this path. Path generation dominates: the module
fixture took 66-73 s (2026-10-21); the repair itself about 0.15 s (exact early
stop, 2026-10-22).

Pinned (2026-10-21, monstera-acceptance-report.md; counts re-pinned 2026-10-22,
user go): the greedy set cover decides on integers only, so the outcome is
deterministic: one round re-inserts sources 744 (island B) then 1273 (island A),
no refusal, two certificate runs (the failing one and the fresh full post-repair
one), 4438 -> 4440 emitted circles and 4160 connectors.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from benchmarks.held_coverage_repair import path_sweeps
from benchmarks.held_coverage_repair import refine_boundary_ordered_path_certified
from benchmarks.held_coverage_repair_types import CertifiedBoundaryOrderedPath
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# Dropped source indices whose swept annulus touches each island
# (build/held-selection-design/repair_cover_diag.py, reporting only).
ISLAND_A_CANDIDATES = frozenset({1273, 1274, 1275, 1276, 1277, 4035, 4036, 4037, 4038, 4046, 6127, 6128})
ISLAND_B_CANDIDATES = frozenset({738, 744, 746, 6159, 6161})
INSERTED_SOURCES = (744, 1273)  # measured 2026-10-21; integer-deterministic tie-break
# Re-pinned 2026-10-22 after the connector repair (+72 circles): 4438 = 2695 selected + 1743 refinement
# circles; 4440 = 4438 + 2 inserted sources (one circle each); 4160 = 4440 transitions (4439 consecutive
# + seam) - 280 fans. Derivation: docs/superpowers/state/coverage-repair-budget-gate-analysis.md.
EMITTED_BEFORE, EMITTED_AFTER, CONNECTORS_AFTER = 4438, 4440, 4160
REPAIR_BUDGET_S = 0.5  # proposed repair budget on Monstera, s (design doc, Cost)
REPORT = Path(__file__).resolve().parents[2] / "build" / "coverage-repair" / "monstera_repair.json"

Certified = tuple[CertifiedBoundaryOrderedPath, EngagementCap, tuple[Point2[WorldXY], ...], ToolRadius]


@pytest.fixture(scope="module")
def certified() -> Certified:
    case = load_held_reference_case("figure8_monstera")
    guide = build_held_reference_raw_guide(case)
    components = figure7_inward_offset(case).components
    assert len(components) == 1, "the Monstera centre domain is one loop without holes"
    source = build_held_reference_path(guide, components, qualify_contacts=True)
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    return refine_boundary_ordered_path_certified(source, components[0], case.tool_radius, cap, corner_approaches=True), cap, components[0], case.tool_radius


def test_monstera_repaired_path_passes_an_independent_fresh_certificate(certified: Certified) -> None:
    """A certificate run here, outside the repair, on the final emission's own sweeps (seam closed)."""
    path, _, boundary, tool_radius = certified
    cuts, connectors = path_sweeps(CentreDomain.build(list(boundary)), boundary, path.emission.refinement.circles)
    assert (len(cuts), len(connectors)) == (EMITTED_AFTER, CONNECTORS_AFTER)
    fresh = CoverageCertificate.build(CentreDomain.build(list(boundary)), tool_radius, cuts, connectors)
    assert fresh.exposed_edge_count == 0
    assert path.ledger.certificate_runs == path.ledger.rounds + 1 == 2


def test_monstera_repair_reinserts_the_pinned_island_candidates(certified: Certified) -> None:
    ledger = certified[0].ledger
    assert ledger.inserted_sources == INSERTED_SOURCES
    assert INSERTED_SOURCES[0] in ISLAND_B_CANDIDATES and INSERTED_SOURCES[1] in ISLAND_A_CANDIDATES
    assert (ledger.rounds, ledger.refusals) == (1, ())
    assert (ledger.emitted_before, ledger.emitted_after) == (EMITTED_BEFORE, EMITTED_AFTER)


def test_monstera_repair_respects_the_engagement_cap_including_the_seam(certified: Certified) -> None:
    """Every emitted transition, the closing seam (last -> first) included."""
    path, cap, _, tool_radius = certified
    circles = path.emission.refinement.circles
    transitions = (*zip(circles, circles[1:]), (circles[-1], circles[0]))
    over = [k for k, (a, b) in enumerate(transitions) if not exact_predecessor_engagement(_paper_candidate(a), _paper_candidate(b), tool_radius).within(cap)]
    assert not over, f"{len(over)} over-cap transitions of {len(transitions)}, first at {over[:5]}"


def test_monstera_repair_reports_circle_count_and_time(certified: Certified) -> None:
    ledger = certified[0].ledger
    record = {
        "emitted_before": ledger.emitted_before,
        "emitted_after": ledger.emitted_after,
        "emitted_change": ledger.emitted_after - ledger.emitted_before,
        "inserted_sources": list(ledger.inserted_sources),
        "refusals": [[r.candidate, list(r.interval), r.round, r.reason] for r in ledger.refusals],
        "rounds": ledger.rounds,
        "repair_s": ledger.repair_s,
        "certificate_s": ledger.certificate_s,
        "certificate_runs": ledger.certificate_runs,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record))
    assert ledger.repair_s <= REPAIR_BUDGET_S, f"repair took {ledger.repair_s:.3f} s"
