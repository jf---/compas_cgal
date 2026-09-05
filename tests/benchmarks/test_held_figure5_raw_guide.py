import math
from fractions import Fraction

import pytest
from compas.tolerance import TOL

import benchmarks.held_figure5_raw_guide as raw_guide_module
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import InadmissibleFigure5BoundarySideError
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundaryVertexId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_raw_guide import build_projection_admissible_hypothesis
from benchmarks.held_figure5_raw_guide import distance_admissible_boundary_sites
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


def _distance(first: Point2[WorldXY], second: Point2[WorldXY]) -> float:
    return math.hypot(
        float(second.x) - float(first.x),
        float(second.y) - float(first.y),
    )


def test_raw_guide_preserves_every_station_and_exposes_explicit_sides() -> None:
    case = load_held_reference_case("figure5")

    guide = build_figure5_raw_guide(case)

    assert guide.station_count == 7075
    assert tuple(run.run_id for run in guide.runs) == tuple(GuideRunId(index) for index in range(len(guide.runs)))
    assert all(tuple(station.ordinal_on_run for station in run.stations) == tuple(GuideRunStationOrdinal(index) for index in range(len(run.stations))) for run in guide.runs)

    station = max(
        (station for run in guide.runs for station in run.stations),
        key=lambda item: len(distance_admissible_boundary_sites(guide, item)),
    )
    sites = distance_admissible_boundary_sites(guide, station)
    candidates = build_distance_admissible_hypotheses(guide, station)

    assert len(sites) >= 2
    assert tuple(candidate.projected_boundary_site for candidate in candidates) == sites
    assert len(set(sites)) == len(sites)
    assert all(0.0 <= float(candidate.boundary_arclength) < float(guide.boundary_length) for candidate in candidates)

    tool_radius = float(case.tool_radius.value)
    for candidate in candidates:
        assert candidate.middle_point == station.middle_point
        assert candidate.distance_admissible_sites == sites
        assert isinstance(candidate.projected_boundary_site.parameter, Fraction)
        assert type(candidate.guide_radius) is GuideRadius
        assert TOL.is_between(
            _distance(candidate.boundary_footpoint, candidate.contact_point),
            tool_radius,
            tool_radius,
        )
        assert TOL.is_between(
            _distance(candidate.contact_point, candidate.center),
            float(candidate.guide_radius.value),
            float(candidate.guide_radius.value),
        )
        expected_diameter = (
            _distance(
                candidate.boundary_footpoint,
                candidate.middle_point,
            )
            - tool_radius
        )
        assert TOL.is_between(
            2.0 * float(candidate.guide_radius.value),
            expected_diameter,
            expected_diameter,
        )


def test_projected_boundary_site_canonicalizes_closed_ring_endpoint() -> None:
    at_last_segment_end = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(64),
        Fraction(1),
        segment_count=65,
    )
    at_first_segment_start = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(0),
        Fraction(0),
        segment_count=65,
    )

    assert at_last_segment_end == at_first_segment_start
    assert at_last_segment_end.segment_id == ProjectionBoundarySegmentId(0)
    assert at_last_segment_end.parameter == 0
    assert at_last_segment_end.vertex_id == ProjectionBoundaryVertexId(0)


def test_hypothesis_batch_projects_boundary_once_per_station(monkeypatch: pytest.MonkeyPatch) -> None:
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    station = guide.runs[41].stations[3]
    original = raw_guide_module._station_records
    calls = 0

    def counted_station_records(*args: object, **kwargs: object) -> object:
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(raw_guide_module, "_station_records", counted_station_records)

    hypotheses = build_distance_admissible_hypotheses(guide, station)

    assert len(hypotheses) >= 2
    assert calls == 1


def test_explicit_hypothesis_factory_rejects_site_outside_distance_budget() -> None:
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    station = guide.runs[0].stations[0]
    admissible = set(distance_admissible_boundary_sites(guide, station))
    rejected = next(
        ProjectionBoundarySite.build(
            ProjectionBoundarySegmentId(ordinal),
            Fraction(1, 2),
            guide.boundary_segment_count,
        )
        for ordinal in range(guide.boundary_segment_count)
        if ProjectionBoundarySite.build(
            ProjectionBoundarySegmentId(ordinal),
            Fraction(1, 2),
            guide.boundary_segment_count,
        )
        not in admissible
    )

    with pytest.raises(InadmissibleFigure5BoundarySideError):
        build_projection_admissible_hypothesis(guide, station, rejected)
