"""Reference-guided Figure 5(a) circle correspondence."""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field
from typing import Literal

import numpy as np
from scipy.optimize import linear_sum_assignment  # type: ignore[import-untyped]

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import AmbiguousFigure5BoundaryPathError
from benchmarks.held_figure5_boundary_path import _canonical_offset_point
from benchmarks.held_figure5_boundary_path import _resolved_site
from benchmarks.held_figure5_boundary_path import _validate_boundary
from benchmarks.held_figure5_path import _CanonicalHypothesis
from benchmarks.held_figure5_path import _canonicalize
from benchmarks.held_figure5_path import _circle_geometry_key
from benchmarks.held_figure5_publisher import Figure5PublisherPathEvidence
from benchmarks.held_figure5_publisher import PublisherTurnOrdinal
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import WorldXY

FIGURE5_REFERENCE_GUIDED_PROVENANCE: Literal["publisher-ordered unchanged repository hypotheses"] = "publisher-ordered unchanged repository hypotheses"
FIGURE5_REFERENCE_GUIDED_CLAIM_BOUNDARY: Literal[
    "reference-guided circle correspondence; no independent path, standard-placement, topology, exact-MAT, or coordinate-parity claim"
] = "reference-guided circle correspondence; no independent path, standard-placement, topology, exact-MAT, or coordinate-parity claim"
FIGURE5_ASSOCIATION_REPORTING_BOUND = Millimetre(0.1502430)
FIGURE5_DISCRETE_ENGAGEMENT_TOLERANCE = Radian(0.001)


class InvalidReferenceGuidedFigure5InputError(BenchmarkError):
    """Inputs do not define the canonical Figure 5 comparison."""


class AmbiguousReferenceGuidedFigure5MarkerError(BenchmarkError):
    """The publisher marker is equidistant from both stream endpoints."""


@dataclass(frozen=True)
class ReferenceGuidedFigure5Circle:
    publisher_turn_ordinal: PublisherTurnOrdinal
    canonical_candidate_ordinal: int
    hypothesis: ProjectionAdmissibleBoundaryHypothesis
    source_run_ids: tuple[GuideRunId, ...]
    assignment_residual: Millimetre
    center_residual: Millimetre
    radius_residual: Millimetre
    circle_locus_residual: Millimetre
    clockwise: Literal[False] = field(default=False, init=False)


@dataclass(frozen=True)
class ReferenceGuidedTransitionBreak:
    after_turn_ordinal: PublisherTurnOrdinal
    reason: str
    continuous_transition_residual: Millimetre


@dataclass(frozen=True)
class ReferenceGuidedCorrespondenceDiagnostics:
    reporting_bound: Millimetre
    center_residuals: tuple[Millimetre, ...]
    radius_residuals: tuple[Millimetre, ...]
    circle_locus_residuals: tuple[Millimetre, ...]
    outlier_ordinals: tuple[PublisherTurnOrdinal, ...]

    @property
    def outlier_count(self) -> int:
        return len(self.outlier_ordinals)


@dataclass(frozen=True)
class ReferenceGuidedEngagementDiagnostics:
    reporting_limit: Radian
    within_family_angles: tuple[Radian, ...]
    within_family_over_tolerance_ordinals: tuple[PublisherTurnOrdinal, ...]
    connector_crossing_angles: tuple[Radian, ...]
    connector_crossing_over_tolerance_ordinals: tuple[PublisherTurnOrdinal, ...]

    @property
    def within_family_over_tolerance_count(self) -> int:
        return len(self.within_family_over_tolerance_ordinals)

    @property
    def connector_crossing_over_tolerance_count(self) -> int:
        return len(self.connector_crossing_over_tolerance_ordinals)


@dataclass(frozen=True)
class PublisherMarkerEndpointEvidence:
    publisher_initial_distance: Millimetre
    publisher_terminal_distance: Millimetre
    interpretation: Literal["initial evidence", "terminal evidence"]


@dataclass(frozen=True)
class ReferenceGuidedMarkerEvidence:
    publisher_initial_distance: Millimetre
    publisher_terminal_distance: Millimetre
    repository_initial_locus_distance: Millimetre
    repository_terminal_locus_distance: Millimetre
    interpretation: Literal["terminal evidence"] = field(default="terminal evidence", init=False)


@dataclass(frozen=True)
class ReferenceGuidedFigure5Correspondence:
    """Publisher order associated with unchanged repository circles."""

    provenance: Literal["publisher-ordered unchanged repository hypotheses"]
    claim_boundary: Literal["reference-guided circle correspondence; no independent path, standard-placement, topology, exact-MAT, or coordinate-parity claim"]
    circles: tuple[ReferenceGuidedFigure5Circle, ...]
    transition_breaks: tuple[ReferenceGuidedTransitionBreak, ...]
    correspondence: ReferenceGuidedCorrespondenceDiagnostics
    engagement: ReferenceGuidedEngagementDiagnostics
    marker: ReferenceGuidedMarkerEvidence
    canonical_candidate_count: int


def _candidate_vector(candidate: _CanonicalHypothesis) -> tuple[float, float, float]:
    hypothesis = candidate.hypothesis
    return float(hypothesis.center.x), float(hypothesis.center.y), float(hypothesis.guide_radius.value)


def _all_raw_hypotheses(guide: Figure5RawGuide) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    return tuple(item for run in guide.runs for station in run.stations for item in build_distance_admissible_hypotheses(guide, station))


def _associate(evidence: Figure5PublisherPathEvidence, canonical: tuple[_CanonicalHypothesis, ...]) -> tuple[ReferenceGuidedFigure5Circle, ...]:
    if not canonical:
        raise InvalidReferenceGuidedFigure5InputError("Reference-guided association requires repository candidates.")
    ordered = tuple(sorted(canonical, key=lambda item: _circle_geometry_key(item.hypothesis)))
    observations = np.asarray([(float(turn.center.x), float(turn.center.y), float(turn.radius)) for turn in evidence.turns])
    candidates = np.asarray([_candidate_vector(candidate) for candidate in ordered])
    delta = observations[:, np.newaxis, :] - candidates[np.newaxis, :, :]
    costs = np.sqrt(np.sum(delta * delta, axis=2))
    observation_ordinals, candidate_ordinals = linear_sum_assignment(costs)
    if tuple(map(int, observation_ordinals)) != tuple(range(len(evidence.turns))):
        raise InvalidReferenceGuidedFigure5InputError("One-to-one association did not cover every publisher turn.")
    circles = []
    for observation_ordinal, candidate_ordinal in zip(observation_ordinals, candidate_ordinals, strict=True):
        row = costs[int(observation_ordinal)]
        selected = float(row[int(candidate_ordinal)])
        dx, dy, dr = delta[int(observation_ordinal), int(candidate_ordinal)]
        center = math.hypot(float(dx), float(dy))
        radius = abs(float(dr))
        circles.append(
            ReferenceGuidedFigure5Circle(
                evidence.turns[int(observation_ordinal)].ordinal,
                int(candidate_ordinal),
                ordered[int(candidate_ordinal)].hypothesis,
                ordered[int(candidate_ordinal)].source_run_ids,
                Millimetre(selected),
                Millimetre(center),
                Millimetre(radius),
                Millimetre(center + radius),
            )
        )
    return tuple(circles)


def _paper_candidate(circle: ReferenceGuidedFigure5Circle) -> PaperCircleCandidate:
    item = circle.hypothesis
    return PaperCircleCandidate.build(center=item.center, guide_radius=item.guide_radius, contact_point=item.contact_point)


def _marker_locus_distance(marker: Point2[WorldXY], circle: ReferenceGuidedFigure5Circle) -> Millimetre:
    item = circle.hypothesis
    return Millimetre(abs(math.dist((float(marker.x), float(marker.y)), (float(item.center.x), float(item.center.y))) - float(item.guide_radius.value)))


def _classify_publisher_marker(marker: Point2[WorldXY], initial: Point2[WorldXY], terminal: Point2[WorldXY], uncertainty: Millimetre) -> PublisherMarkerEndpointEvidence:
    initial_distance = Millimetre(math.dist((float(marker.x), float(marker.y)), (float(initial.x), float(initial.y))))
    terminal_distance = Millimetre(math.dist((float(marker.x), float(marker.y)), (float(terminal.x), float(terminal.y))))
    if abs(initial_distance - terminal_distance) <= 2 * uncertainty:
        raise AmbiguousReferenceGuidedFigure5MarkerError("Publisher marker endpoint separation is inside the two-sided evidence uncertainty.")
    interpretation: Literal["initial evidence", "terminal evidence"] = "terminal evidence" if terminal_distance < initial_distance else "initial evidence"
    return PublisherMarkerEndpointEvidence(initial_distance, terminal_distance, interpretation)


def _transition_residual(component: tuple[Point2[WorldXY], ...], circle: ReferenceGuidedFigure5Circle, bound: Millimetre) -> tuple[Millimetre, str]:
    try:
        site = _resolved_site(component, circle.hypothesis, bound)
        canonical = _canonical_offset_point(component, site)
    except AmbiguousFigure5BoundaryPathError as error:
        return Millimetre(math.inf), str(error)
    q = circle.hypothesis.contact_point
    residual = Millimetre(math.dist((float(q.x), float(q.y)), (float(canonical.x), float(canonical.y))))
    reason = (
        "Repository circle q is not exactly on the radius-one boundary component." if residual > 0 else "Continuous transition lineage is unavailable from the approximate guide."
    )
    return residual, reason


def build_reference_guided_figure5_correspondence(
    guide: Figure5RawGuide,
    inward_components: tuple[tuple[Point2[WorldXY], ...], ...],
) -> ReferenceGuidedFigure5Correspondence:
    if type(guide) is not Figure5RawGuide or guide.case.name != "figure5" or len(inward_components) != 1:
        raise InvalidReferenceGuidedFigure5InputError("Reference-guided correspondence requires canonical Figure 5 and one inward component.")
    marker = guide.case.start_marker
    if marker is None:
        raise InvalidReferenceGuidedFigure5InputError("Reference-guided Figure 5 requires its publisher marker evidence.")
    component = inward_components[0]
    _validate_boundary(component)
    evidence = load_figure5_publisher_path()
    canonical = _canonicalize(_all_raw_hypotheses(guide))
    circles = _associate(evidence, canonical)
    marker_uncertainty = Millimetre(float(guide.site_budget.reconstruction_bound) + float(guide.site_budget.projection_bound))
    endpoint = _classify_publisher_marker(marker, evidence.stream_start, evidence.stream_end, marker_uncertainty)
    if endpoint.interpretation != "terminal evidence":
        raise InvalidReferenceGuidedFigure5InputError("Publisher marker does not identify the terminal stream endpoint.")

    breaks = []
    for index, (start, end) in enumerate(zip(circles, circles[1:], strict=False)):
        start_residual, start_reason = _transition_residual(component, start, guide.site_budget.admissible_distance_gap)
        end_residual, end_reason = _transition_residual(component, end, guide.site_budget.admissible_distance_gap)
        breaks.append(
            ReferenceGuidedTransitionBreak(
                PublisherTurnOrdinal(index), start_reason if start_residual >= end_residual else end_reason, Millimetre(max(start_residual, end_residual))
            )
        )

    candidates = tuple(_paper_candidate(circle) for circle in circles)
    angles = tuple(maximum_predecessor_engagement(start, end, guide.case.tool_radius) for start, end in zip(candidates, candidates[1:], strict=False))
    limit = Radian(math.radians(float(guide.case.tea_cap)) + float(FIGURE5_DISCRETE_ENGAGEMENT_TOLERANCE))
    # A circular turn consumes four cubics; the eleven connector-delimited
    # family transfers each contain more than one turn's cubic budget.
    complex_crossings = {index for index, turn in enumerate(evidence.turns[:-1]) if turn.connector_after.cubic_count > 4}
    within = tuple(angle for index, angle in enumerate(angles) if index not in complex_crossings)
    crossings = tuple(angle for index, angle in enumerate(angles) if index in complex_crossings)
    within_over = tuple(circles[index + 1].publisher_turn_ordinal for index, angle in enumerate(angles) if index not in complex_crossings and angle > limit)
    crossing_over = tuple(circles[index + 1].publisher_turn_ordinal for index, angle in enumerate(angles) if index in complex_crossings and angle > limit)
    locus = tuple(circle.circle_locus_residual for circle in circles)
    outliers = tuple(circle.publisher_turn_ordinal for circle in circles if circle.circle_locus_residual > FIGURE5_ASSOCIATION_REPORTING_BOUND)
    return ReferenceGuidedFigure5Correspondence(
        FIGURE5_REFERENCE_GUIDED_PROVENANCE,
        FIGURE5_REFERENCE_GUIDED_CLAIM_BOUNDARY,
        circles,
        tuple(breaks),
        ReferenceGuidedCorrespondenceDiagnostics(
            FIGURE5_ASSOCIATION_REPORTING_BOUND, tuple(circle.center_residual for circle in circles), tuple(circle.radius_residual for circle in circles), locus, outliers
        ),
        ReferenceGuidedEngagementDiagnostics(limit, within, within_over, crossings, crossing_over),
        ReferenceGuidedMarkerEvidence(
            endpoint.publisher_initial_distance, endpoint.publisher_terminal_distance, _marker_locus_distance(marker, circles[0]), _marker_locus_distance(marker, circles[-1])
        ),
        len(canonical),
    )
