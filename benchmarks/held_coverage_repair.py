"""Select, certify, repair: coverage-preserving boundary-order selection.

Boundary-order selection (``select_boundary_ordered_sources``) and engagement
refinement (``refine_figure5_engagement``) run unchanged. Before every
certificate run the exact NOTHING pass (``benchmarks.held_nothing_circles``,
closed over the seam when ``settings.close_seam``) prunes the refined circles,
so the certificate, the repair's witnesses and the connector post-gate all
decide the circle set that is emitted. A full coverage
certificate then runs on the emitted path, closing seam included (when
``settings.close_seam``). When it
raises ``UncoveredPocketError``, this stage re-inserts dropped source circles
that contain the certificate's exact witnesses, each into its own
boundary-order slot, re-refines from that slot on under the
engagement cap, and runs a FRESH FULL certificate again (user ruling: the
post-repair gate is a full certificate, not windowed re-certification).
Design: docs/held_coverage_preserving_selection.md.

Decisions made here, all exact:

- witness membership: (max(R - r, 0))^2 <= |w - c|^2 <= (R + r)^2 on the
  exactly injected binary64 centre and radii, with ``Fraction`` arithmetic,
  on the circle exactly as the refinement emits it;
- cyclic boundary order of a candidate between its kept neighbours, compared
  on ``(segment, Fraction parameter)``;
- the greedy set cover compares integer counts and indices only.

The outward-rounded candidate boxes are a filter (a superset), never a
decision. The certificate is an injected dependency (``CoverageOracle``) so
contract tests run the real one.
"""

from __future__ import annotations

import time

from benchmarks.held_boundary_order_spacing import InvalidBoundarySpacingInputError
from benchmarks.held_boundary_order_spacing import select_boundary_ordered_sources
from benchmarks.held_connector_engagement import ConnectorEngagementAssessment
from benchmarks.held_connector_engagement import ConnectorEngagementCertificate
from benchmarks.held_connector_engagement import UncertifiedConnectorEngagementError
from benchmarks.held_coverage_repair_cover import repair_round
from benchmarks.held_coverage_repair_errors import UncoveredPocketWithoutWitnessesError
from benchmarks.held_coverage_repair_splice import emitted_form
from benchmarks.held_coverage_repair_types import BoundaryOrderedEmission
from benchmarks.held_coverage_repair_types import CertifiedBoundaryOrderedPath
from benchmarks.held_coverage_repair_types import CoverageOracle
from benchmarks.held_coverage_repair_types import CoverageRepairLedger
from benchmarks.held_coverage_repair_types import EngagementRefinementSettings
from benchmarks.held_coverage_repair_types import RepairInsertion
from benchmarks.held_coverage_repair_types import RepairRefusal
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import HypothesisFigure5Path
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_nothing_circles import NothingCirclePass
from benchmarks.held_nothing_circles import remove_nothing_circles
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


def path_sweeps(
    domain: CentreDomain,
    boundary: tuple[Point2[WorldXY], ...],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    *,
    close_seam: bool = True,
) -> tuple[tuple[FullCircleSweep, ...], tuple[BoundaryConnectorSweep, ...]]:
    """The certificate's swept regions for emitted circles and their CCW connectors.

    The connectors are exactly ``Figure5BoundaryConnectors``' transitions:
    consecutive contacts, plus the seam connector from the last contact back
    to the first unless the generator declares the seam coincident. Fan
    transitions (same site) are stationary and sweep only disk(q, r), inside
    the cut at q, so they contribute no region.

    Args:
        close_seam: False drops the seam connector; it reproduces the open
            path emitted before the seam was closed, whose failure stays a
            tested contract.
    """
    connectors = Figure5BoundaryConnectors.build(boundary, circles)
    pairs = list(zip(circles, circles[1:]))
    if close_seam and connectors.seam is Figure5SeamClosure.CONNECTOR:
        pairs.append((circles[-1], circles[0]))
    cuts = tuple(FullCircleSweep.build(circle.center, circle.radius) for circle in circles)
    sweeps = tuple(
        BoundaryConnectorSweep.build(
            domain,
            CentreBoundarySite.build(0, int(a.boundary_site.segment_id), a.boundary_site.parameter),
            CentreBoundarySite.build(0, int(b.boundary_site.segment_id), b.boundary_site.parameter),
        )
        for a, b in pairs
        if a.boundary_site != b.boundary_site
    )
    return cuts, sweeps


def repair_coverage(
    sources: tuple[Figure5CounterclockwiseCircle, ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
    emission: BoundaryOrderedEmission,
    oracle: CoverageOracle,
) -> CertifiedBoundaryOrderedPath:
    """Certify the emission; on ``UncoveredPocketError`` re-insert covering sources and re-certify in full.

    Every certificate run, and so every repair round's witnesses, decides the
    emission's circles after the exact NOTHING pass
    (``remove_emission_nothing_circles``); the emitted path is that pass's
    result (``CertifiedBoundaryOrderedPath.circles``). A removed circle's
    annulus lies in the union of the kept circles before it, so no witness
    can lie in it.

    Once coverage holds, the connector-engagement post-gate runs on the final
    path (``ConnectorEngagementAssessment``, the decisions of
    ``ConnectorEngagementCertificate``): every cut connector must be certified
    by tier 1 or 1b, the seam included when ``settings.close_seam``.

    Args:
        sources: Every source circle in boundary order, as passed to selection.
        boundary: The CCW centre-domain loop Omega the contacts lie on.
        settings: The refinement settings the emission was refined with;
            ``settings.close_seam`` also selects the swept seam (``path_sweeps``).
        emission: The selection and its refinement.
        oracle: The full coverage certificate.

    Raises:
        CoverageNotRepairableBySelectionError: a witness lies in no dropped source.
        CoverageRepairRefusedError: every candidate for a witness was refused.
        RepairOrderViolationError: a candidate is outside its boundary slot.
        UncoveredPocketWithoutWitnessesError: the oracle reported a pocket
            without a witness, so no round can make progress.
        UncertifiedConnectorEngagementError: a cut connector of the final
            path has no engagement certificate (an emission not refined with
            these settings).
    """
    domain = CentreDomain.build(list(boundary))
    forms: dict[int, Figure5CounterclockwiseCircle] = {}
    refusals: list[RepairRefusal] = []
    insertions: list[RepairInsertion] = []
    repair_s = certificate_s = 0.0
    rounds = runs = 0
    nothing = remove_emission_nothing_circles(emission, settings)
    emitted_before = len(nothing.circles)
    while True:
        cuts, connectors = path_sweeps(domain, boundary, nothing.circles, close_seam=settings.close_seam)
        started = time.perf_counter()
        runs += 1
        try:
            certificate = oracle(domain, settings.tool_radius, cuts, connectors)
        except UncoveredPocketError as error:
            certificate_s += time.perf_counter() - started
            if not error.witnesses:
                raise UncoveredPocketWithoutWitnessesError(error.undecided_edge_count) from error
            witnesses = error.witnesses
        else:
            certificate_s += time.perf_counter() - started
            break
        rounds += 1
        started = time.perf_counter()
        emission, inserted = repair_round(sources, boundary, settings, emission, witnesses, forms, rounds, refusals)
        repair_s += time.perf_counter() - started
        insertions.extend(inserted)
        nothing = remove_emission_nothing_circles(emission, settings)
    _require_connector_certificate(boundary, settings, nothing.circles)
    ledger = CoverageRepairLedger.build(
        insertions,
        refusals,
        rounds=rounds,
        emitted_before=emitted_before,
        emitted_after=len(nothing.circles),
        repair_s=repair_s,
        certificate_s=certificate_s,
        certificate_runs=runs,
    )
    return CertifiedBoundaryOrderedPath.build(sources, emission, nothing, certificate, ledger)


def remove_emission_nothing_circles(emission: BoundaryOrderedEmission, settings: EngagementRefinementSettings) -> NothingCirclePass[Figure5CounterclockwiseCircle]:
    """The exact NOTHING pass on an emission's refined circles; the seam is re-decided when ``settings.close_seam``.

    The emission itself is not pruned: its ``original_indices`` and decision
    trace index the refined circles, and the repair splices into them.
    """
    return remove_nothing_circles(emission.refinement.circles, _paper_candidate, settings.tool_radius, settings.cap, closed=settings.close_seam)


def _require_connector_certificate(boundary: tuple[Point2[WorldXY], ...], settings: EngagementRefinementSettings, circles: tuple[Figure5CounterclockwiseCircle, ...]) -> None:
    """Post-gate: every cut connector certified; the seam is cut only on a closed path.

    Closed: ``ConnectorEngagementCertificate.build`` (seam included). Open: the
    same decisions (``ConnectorEngagementAssessment``) with the uncut seam's
    verdict ignored.
    """
    connectors = Figure5BoundaryConnectors.build(boundary, circles)
    if settings.close_seam:
        ConnectorEngagementCertificate.build(circles, connectors, settings.tool_radius, settings.cap)
        return
    assessment = ConnectorEngagementAssessment.build(circles, connectors, settings.tool_radius, settings.cap)
    refused = tuple(verdict.connector for verdict in assessment.verdicts if not verdict.certified and not verdict.seam)
    if refused:
        raise UncertifiedConnectorEngagementError(refused)


def _refined_selection(
    sources: tuple[Figure5CounterclockwiseCircle, ...],
    selected: tuple[int, ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
) -> BoundaryOrderedEmission:
    refinement = refine_figure5_engagement(
        tuple(sources[index] for index in selected),
        boundary,
        settings.tool_radius,
        settings.cap,
        max_depth=settings.max_depth,
        max_circles=settings.max_circles,
        corner_approaches=settings.corner_approaches,
        close_seam=settings.close_seam,
    )
    return BoundaryOrderedEmission.build(selected, refinement)


def certify_boundary_ordered_selection(
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    circle_sources: tuple[tuple[GuideRunId, ...], ...],
    boundary: tuple[Point2[WorldXY], ...],
    settings: EngagementRefinementSettings,
    oracle: CoverageOracle,
) -> CertifiedBoundaryOrderedPath:
    """Select in boundary order, refine, certify, repair: ``refine_boundary_ordered_path`` plus the coverage gate.

    The selection and the initial refinement are exactly those of
    ``refine_boundary_ordered_path`` (sources are first placed on their side,
    as there).
    """
    sources = tuple(emitted_form(circle, boundary) for circle in circles)
    selected = select_boundary_ordered_sources(sources, circle_sources, settings.tool_radius, settings.cap)
    return repair_coverage(sources, boundary, settings, _refined_selection(sources, selected, boundary, settings), oracle)


def refine_boundary_ordered_path_certified(
    source: HypothesisFigure5Path,
    boundary: tuple[Point2[WorldXY], ...],
    tool_radius: ToolRadius,
    cap: EngagementCap,
    *,
    oracle: CoverageOracle = CoverageCertificate.build,
    corner_approaches: bool = False,
) -> CertifiedBoundaryOrderedPath:
    """``refine_boundary_ordered_path`` followed by the always-on certificate and the repair.

    The reached-run ledger is checked on the selection, before refinement, as
    in ``refine_boundary_ordered_path``: a mismatch fails before any
    certificate runs.

    Raises:
        InvalidBoundarySpacingInputError: the selection loses a reached run.
        CoverageNotRepairableBySelectionError, CoverageRepairRefusedError,
        RepairOrderViolationError, UncoveredPocketWithoutWitnessesError: see
        ``repair_coverage``.
    """
    settings = EngagementRefinementSettings.build(tool_radius, cap, corner_approaches=corner_approaches)
    sources = tuple(emitted_form(circle, boundary) for circle in source.path.circles)
    selected = select_boundary_ordered_sources(sources, source.circle_sources, tool_radius, cap)
    represented = frozenset(run for index in selected for run in source.circle_sources[index])
    if represented != source.reached_run_ids:
        raise InvalidBoundarySpacingInputError("Boundary spacing source attribution differs from the reached-run ledger.")
    return repair_coverage(sources, boundary, settings, _refined_selection(sources, selected, boundary, settings), oracle)
