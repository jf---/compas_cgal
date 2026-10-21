"""The five Held workloads generated on Ω* (exact ∂Ω* contacts, order and connectors).

Each generator mirrors its polygon-Ω counterpart in
:mod:`benchmarks.held_toolpath_benchmark` stage by stage (guide, domain,
initial placement, engagement refinement, connectors) with the same contact
qualification, run preservation and spacing, and replaces the float offset
polygon by Ω*. The polygon's ``corner_approaches`` has no Ω* counterpart (∂Ω*
has no reflex vertex). Select it with ``PathDomain.OMEGA``.
"""

from __future__ import annotations

import enum
import functools
import math
import time
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction

from benchmarks.errors import BenchmarkError
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_path import OmegaHypothesisPath
from benchmarks.held_omega_path import omega_reference_path
from benchmarks.held_omega_refinement import OmegaRefinement
from benchmarks.held_omega_refinement import refine_omega_boundary_ordered
from benchmarks.held_omega_refinement import refine_omega_engagement
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Millimetre


class PathDomain(enum.Enum):
    """Which centre domain the Held generators decide on."""

    POLYGON = "polygon"
    OMEGA = "omega"


class UnknownOmegaWorkloadError(BenchmarkError):
    """The requested workload is not one of the five Held generators."""


class Spacing(enum.Enum):
    """The polygon path's two spacing regimes."""

    LANE = "lane"
    BOUNDARY = "boundary"


@dataclass(frozen=True)
class OmegaWorkload:
    """One workload's polygon-path settings, carried over."""

    case_name: str
    spacing: Spacing
    qualify_contacts: bool
    preserve_source_runs: bool


OMEGA_WORKLOADS: dict[str, OmegaWorkload] = {
    "figure5_refined": OmegaWorkload("figure5", Spacing.LANE, qualify_contacts=False, preserve_source_runs=False),
    "figure8_upper_lane": OmegaWorkload("figure8_upper", Spacing.LANE, qualify_contacts=False, preserve_source_runs=True),
    "figure8_upper_boundary": OmegaWorkload("figure8_upper", Spacing.BOUNDARY, qualify_contacts=False, preserve_source_runs=True),
    "figure8_crossed_skis_boundary": OmegaWorkload("figure8_crossed_skis", Spacing.BOUNDARY, qualify_contacts=False, preserve_source_runs=True),
    "figure8_monstera_boundary_qualified_corners": OmegaWorkload("figure8_monstera", Spacing.BOUNDARY, qualify_contacts=True, preserve_source_runs=True),
}


@dataclass(frozen=True)
class OmegaGeneratedToolpath:
    """One workload's Ω* path: binary64 circles and exact ∂Ω* connectors, with its domain and oracle."""

    workload: str
    case: HeldReferenceCase
    declared: DeclaredBoundary
    loop: OmegaLoop
    clearance: DeclaredClearance
    initial: OmegaHypothesisPath
    refinement: OmegaRefinement
    connectors: OmegaConnectors
    stage_seconds: dict[str, float]

    @property
    def circles(self) -> tuple[OmegaCircle, ...]:
        return self.refinement.circles

    @property
    def generation_seconds(self) -> float:
        return sum(self.stage_seconds.values())

    def audit(self) -> HeldGougeAuditReport:
        """Exact gouge audit of every circle and connector sub-piece against the declared boundary."""
        return audit_omega_path(self.clearance, self.circles, self.connectors.connectors)


@dataclass(frozen=True)
class OmegaPathMetrics:
    """Counts and reported doubles of one Ω* path (the polygon benchmark's metrics, where they apply)."""

    circles: int
    connectors: int
    path_length_mm: float
    max_engagement_deg: float
    unresolved_engagements: int
    seam: str
    limited_circles: int
    arc_connector_pieces: int
    line_connector_pieces: int
    beyond_gap_contacts: int
    convex_vertex_contacts: int
    max_junction_gap_mm: float


def _guide(case: HeldReferenceCase) -> Figure5RawGuide:
    return build_figure5_raw_guide(case) if case.name == "figure5" else build_held_reference_raw_guide(case)


def generate_omega(workload: str) -> OmegaGeneratedToolpath:
    """Run one workload's generator on Ω* stage by stage.

    Raises:
        UnknownOmegaWorkloadError: ``workload`` is not in ``OMEGA_WORKLOADS``.
    """
    if workload not in OMEGA_WORKLOADS:
        raise UnknownOmegaWorkloadError(f"unknown workload {workload!r}; known: {list(OMEGA_WORKLOADS)}")
    spec = OMEGA_WORKLOADS[workload]
    case = load_held_reference_case(spec.case_name)
    seconds: dict[str, float] = {}
    started = time.perf_counter()
    guide = _guide(case)
    seconds["guide"] = time.perf_counter() - started
    started = time.perf_counter()
    tool_radius = Fraction(float(case.tool_radius.value))
    declared = declared_boundary_from_reference(case.name, case.boundary)
    domain = CentreDomain.build(declared, tool_radius)
    loop = OmegaLoop.build(domain)
    clearance = DeclaredClearance.build(declared, tool_radius)
    seconds["omega"] = time.perf_counter() - started
    start, start_radius = case.start_marker, case.start_marker_radius
    if start is None or start_radius is None:
        raise UnknownOmegaWorkloadError(f"{case.name} has no published start circle.")
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    started = time.perf_counter()
    hypotheses = tuple(hypothesis for run in guide.runs for station in run.stations for hypothesis in build_distance_admissible_hypotheses(guide, station))
    initial = omega_reference_path(
        loop,
        hypotheses,
        frozenset(run.run_id for run in guide.runs),
        gap=ExactMillimetre(Fraction(float(guide.site_budget.admissible_distance_gap))),
        published_start_center=start,
        published_start_radius=start_radius,
        start_evidence_bound=Millimetre(float(guide.site_budget.reconstruction_bound) + float(guide.site_budget.projection_bound)),
        tool_radius=case.tool_radius,
        cap=cap,
        qualify_contacts=spec.qualify_contacts,
        preserve_source_runs=spec.preserve_source_runs,
    )
    seconds["initial"] = time.perf_counter() - started
    started = time.perf_counter()
    if spec.spacing is Spacing.BOUNDARY:
        _, refinement = refine_omega_boundary_ordered(initial.circles, initial.circle_sources, initial.reached_run_ids, loop, clearance, case.tool_radius, cap)
    else:
        refinement = refine_omega_engagement(initial.circles, loop, clearance, case.tool_radius, cap)
    seconds["refine"] = time.perf_counter() - started
    started = time.perf_counter()
    connectors = OmegaConnectors.build(loop.pieces, refinement.circles)
    seconds["connectors"] = time.perf_counter() - started
    return OmegaGeneratedToolpath(workload, case, declared, loop, clearance, initial, refinement, connectors, seconds)


def _engagement_deg(path: OmegaGeneratedToolpath, first: OmegaCircle, second: OmegaCircle) -> float | None:
    try:
        return math.degrees(float(maximum_predecessor_engagement(first.paper_candidate(), second.paper_candidate(), path.case.tool_radius)))
    except StandardPlacementFragmentationError:
        return None


def omega_metrics(path: OmegaGeneratedToolpath) -> OmegaPathMetrics:
    """Reported counts and lengths of an Ω* path."""
    circles, connectors = path.circles, path.connectors.connectors
    length = sum(2 * math.pi * float(circle.radius.value) for circle in circles) + sum(float(connector.boundary_progress) for connector in connectors)
    angles = [_engagement_deg(path, first, second) for first, second in zip(circles[:-1], circles[1:], strict=True)]
    resolved = [angle for angle in angles if angle is not None]
    kinds = [piece.piece.kind.value for connector in connectors for piece in connector.pieces]
    return OmegaPathMetrics(
        circles=len(circles),
        connectors=len(connectors),
        path_length_mm=length,
        max_engagement_deg=max(resolved, default=0.0),
        unresolved_engagements=len(angles) - len(resolved),
        seam=path.connectors.seam.value,
        limited_circles=path.refinement.limited,
        arc_connector_pieces=sum(1 for kind in kinds if kind != "offset_line"),
        line_connector_pieces=sum(1 for kind in kinds if kind == "offset_line"),
        beyond_gap_contacts=len(path.initial.contacts.beyond_gap),
        convex_vertex_contacts=len(path.initial.contacts.at_convex_vertex),
        max_junction_gap_mm=max((_junction_gap(circle) for circle in circles if not circle.limited), default=0.0),
    )


def _junction_gap(circle: OmegaCircle) -> float:
    """Reported ``| |c - q| - R |``: how far the binary64 circle passes from its exact site (the circle-connector junction)."""
    x, y = circle.site.point().as_floats()
    return abs(math.hypot(x - float(circle.center.x), y - float(circle.center.y)) - float(circle.radius.value))


OMEGA_GENERATORS: dict[str, Callable[[], OmegaGeneratedToolpath]] = {name: functools.partial(generate_omega, name) for name in OMEGA_WORKLOADS}
