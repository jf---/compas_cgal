import math
from dataclasses import replace

import pytest

from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_path import _canonicalize
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_figure5_reference_guided import FIGURE5_REFERENCE_GUIDED_PROVENANCE
from benchmarks.held_figure5_reference_guided import ReferenceGuidedFigure5Circle
from benchmarks.held_figure5_reference_guided import AmbiguousReferenceGuidedFigure5MarkerError
from benchmarks.held_figure5_reference_guided import ReferenceGuidedFigure5Correspondence
from benchmarks.held_figure5_reference_guided import _associate
from benchmarks.held_figure5_reference_guided import _classify_publisher_marker
from benchmarks.held_figure5_reference_guided import build_reference_guided_figure5_correspondence
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal.adaptive.units import Millimetre


def _geometry_key(candidate: ReferenceGuidedFigure5Circle) -> tuple[float, ...]:
    hypothesis = candidate.hypothesis
    return (
        float(hypothesis.center.x),
        float(hypothesis.center.y),
        float(hypothesis.guide_radius.value),
        float(hypothesis.contact_point.x),
        float(hypothesis.contact_point.y),
    )


def test_reference_guided_correspondence_selects_unchanged_raw_geometry_in_publisher_order() -> None:
    case = load_held_reference_case("figure5")
    guide = build_figure5_raw_guide(case)
    raw = tuple(hypothesis for run in guide.runs for station in run.stations for hypothesis in build_distance_admissible_hypotheses(guide, station))

    result = build_reference_guided_figure5_correspondence(
        guide,
        figure7_inward_offset(case).components,
    )

    raw_geometry = {
        (
            float(hypothesis.center.x),
            float(hypothesis.center.y),
            float(hypothesis.guide_radius.value),
            float(hypothesis.contact_point.x),
            float(hypothesis.contact_point.y),
        )
        for hypothesis in raw
    }
    assert result.provenance == FIGURE5_REFERENCE_GUIDED_PROVENANCE
    assert type(result) is ReferenceGuidedFigure5Correspondence
    assert len(result.circles) == 265
    assert all(not circle.clockwise for circle in result.circles)
    assert all(_geometry_key(circle) in raw_geometry for circle in result.circles)
    assert len({circle.canonical_candidate_ordinal for circle in result.circles}) == 265
    assert tuple(circle.publisher_turn_ordinal for circle in result.circles) == tuple(range(265))
    assert len(result.transition_breaks) == 264
    assert result.canonical_candidate_count == 30684
    assert result.correspondence.outlier_count == 4
    assert max(result.correspondence.circle_locus_residuals) == pytest.approx(0.173275, abs=1e-6)
    assert len(result.correspondence.center_residuals) == 265
    assert len(result.correspondence.radius_residuals) == 265
    assert result.engagement.within_family_over_tolerance_count + result.engagement.connector_crossing_over_tolerance_count == 137
    assert len(result.engagement.connector_crossing_angles) == 11
    assert len(result.transition_breaks) == 264
    assert {transition.reason for transition in result.transition_breaks} == {"Repository circle q is not exactly on the radius-one boundary component."}
    assert result.correspondence.outlier_count == sum(residual > result.correspondence.reporting_bound for residual in result.correspondence.circle_locus_residuals)
    assert result.engagement.within_family_over_tolerance_count == sum(angle > math.radians(80.0) + 0.001 for angle in result.engagement.within_family_angles)
    assert result.marker.publisher_terminal_distance < result.marker.publisher_initial_distance
    assert result.marker.repository_terminal_locus_distance < result.marker.repository_initial_locus_distance
    assert result.claim_boundary == ("reference-guided circle correspondence; no independent path, standard-placement, topology, exact-MAT, or coordinate-parity claim")
    permuted = _associate(load_figure5_publisher_path(), tuple(reversed(_canonicalize(raw))))
    assert tuple(_geometry_key(circle) for circle in permuted) == tuple(_geometry_key(circle) for circle in result.circles)


def test_publisher_marker_classification_is_endpoint_owned_and_fails_a_tie() -> None:
    marker = load_held_reference_case("figure5").start_marker
    assert marker is not None
    initial = type(marker).build(float(marker.x) + 1.0, float(marker.y))
    terminal = type(marker).build(float(marker.x) + 0.25, float(marker.y))

    evidence = _classify_publisher_marker(marker, initial, terminal, Millimetre(0.01))
    reversed_evidence = _classify_publisher_marker(marker, terminal, initial, Millimetre(0.01))

    assert evidence.interpretation == "terminal evidence"
    assert reversed_evidence.interpretation == "initial evidence"
    with pytest.raises(AmbiguousReferenceGuidedFigure5MarkerError):
        _classify_publisher_marker(marker, initial, replace(initial), Millimetre(0.01))
    with pytest.raises(AmbiguousReferenceGuidedFigure5MarkerError):
        _classify_publisher_marker(marker, initial, type(marker).build(float(initial.x) + 0.01, float(initial.y)), Millimetre(0.01))
