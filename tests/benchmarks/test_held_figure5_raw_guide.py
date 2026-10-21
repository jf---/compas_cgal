import dataclasses
import functools
import math
from collections.abc import Sequence
from fractions import Fraction

import pytest
from compas.tolerance import TOL

import benchmarks.held_figure5_raw_guide as raw_guide_module
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
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
from benchmarks.held_figure5_raw_guide_native import NativeStationBoundaryAnalysis
from benchmarks.held_figure5_raw_guide_native import native_station_boundary_analyses
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.held_boundary_clearance import HeldBoundaryClearanceIndex
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


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
    # Wired native path: one exact batch per guide analyses every station
    # exactly once; the brute-force Fraction oracle is never on the path.
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    original_batch = HeldBoundaryClearanceIndex.distance_admissible_sites
    original_records = raw_guide_module._station_records
    batches: list[int] = []
    oracle_calls = 0

    def counted_batch(self: HeldBoundaryClearanceIndex, stations: Sequence[Point2[WorldXY]], gap: Millimetre) -> object:
        batches.append(len(stations))
        return original_batch(self, stations, gap)

    def counted_station_records(*args: object, **kwargs: object) -> object:
        nonlocal oracle_calls
        oracle_calls += 1
        return original_records(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(HeldBoundaryClearanceIndex, "distance_admissible_sites", counted_batch)
    monkeypatch.setattr(raw_guide_module, "_station_records", counted_station_records)

    station = guide.runs[41].stations[3]
    hypotheses = build_distance_admissible_hypotheses(guide, station)
    assert build_distance_admissible_hypotheses(guide, station) == hypotheses
    assert build_distance_admissible_hypotheses(guide, guide.runs[0].stations[0])

    assert len(hypotheses) >= 2
    assert batches == [guide.station_count]
    assert oracle_calls == 0


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


def test_guide_derived_values_never_write_the_frozen_guide() -> None:
    # Review I1.1: memoised values live in the guide's own mutable memo, never
    # in the frozen instance's __dict__ (no cached_property on frozen types).
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    assert not any(isinstance(value, functools.cached_property) for value in vars(Figure5RawGuide).values())
    first = (guide.rational_boundary_segments, guide.boundary_side_lengths, guide.boundary_side_offsets, guide.station_boundary_analyses)
    assert set(vars(guide)) == {field.name for field in dataclasses.fields(guide)}
    second = (guide.rational_boundary_segments, guide.boundary_side_lengths, guide.boundary_side_offsets, guide.station_boundary_analyses)
    assert all(a is b for a, b in zip(first, second, strict=True))
    twin = dataclasses.replace(guide)
    assert twin == guide and twin.derived is not guide.derived and twin.derived.station_boundary_analyses is None


def test_native_module_runs_the_wired_conversion(monkeypatch: pytest.MonkeyPatch) -> None:
    # Review I6.2: one native analysis path; this module re-shapes the wired batch.
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    calls: list[Figure5RawGuide] = []
    wired = raw_guide_module._native_station_boundary_analyses

    def counted(argument: Figure5RawGuide) -> tuple[tuple[raw_guide_module._StationBoundaryAnalysis, ...], ...]:
        calls.append(argument)
        return wired(argument)

    monkeypatch.setattr(raw_guide_module, "_native_station_boundary_analyses", counted)
    native = native_station_boundary_analyses(guide)
    assert calls == [guide]
    assert native == tuple(
        tuple(NativeStationBoundaryAnalysis(analysis.distance_admissible_sites, analysis.records) for analysis in run) for run in guide.station_boundary_analyses
    )
