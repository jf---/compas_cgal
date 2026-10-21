"""The generator's station analysis is the native batch, and it builds what the Fraction oracle builds."""

from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import benchmarks.held_figure5_raw_guide as raw_guide_module
from benchmarks.held_figure5_raw_guide import Figure5BoundarySiteBudget
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import Figure5RawGuideRun
from benchmarks.held_figure5_raw_guide import Figure5RawGuideStation
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import InvalidFigure5RawGuideError
from benchmarks.held_figure5_raw_guide import _analyze_station_boundary
from benchmarks.held_figure5_raw_guide import _build_projection_admissible_hypothesis
from benchmarks.held_figure5_raw_guide import _native_station_boundary
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_raw_guide import distance_admissible_boundary_sites
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


def _assert_wired_matches_oracle(guide: Figure5RawGuide, *, build_hypotheses: bool) -> int:
    """Compare decisions, admissible records and (optionally) hypotheses per station; return the site count."""
    compared = 0
    for run in guide.runs:
        for station in run.stations:
            oracle = _analyze_station_boundary(guide, station)
            wired = _native_station_boundary(guide, station)
            assert wired.distance_admissible_sites == oracle.distance_admissible_sites, (station.run_id, station.ordinal_on_run)
            oracle_records = {site: (footpoint, distance) for site, footpoint, distance in oracle.records}
            assert tuple(record[0] for record in wired.records) == wired.distance_admissible_sites
            for site, footpoint, distance in wired.records:
                assert oracle_records[site] == (footpoint, distance), (station.run_id, station.ordinal_on_run, site)
            if build_hypotheses:
                expected = tuple(_build_projection_admissible_hypothesis(guide, station, site, oracle) for site in oracle.distance_admissible_sites)
                assert build_distance_admissible_hypotheses(guide, station) == expected, (station.run_id, station.ordinal_on_run)
            compared += len(wired.distance_admissible_sites)
    return compared


def test_wired_hypotheses_equal_fraction_oracle_on_every_figure5_station() -> None:
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))

    compared = _assert_wired_matches_oracle(guide, build_hypotheses=True)

    assert guide.station_count == 7075
    assert compared == 31742


def test_native_batch_runs_once_per_guide(monkeypatch: pytest.MonkeyPatch) -> None:
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    original = raw_guide_module._native_station_boundary_analyses
    calls = 0

    def counted(*args: object, **kwargs: object) -> object:
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(raw_guide_module, "_native_station_boundary_analyses", counted)

    for station in (guide.runs[41].stations[3], guide.runs[0].stations[0], guide.runs[-1].stations[-1]):
        assert build_distance_admissible_hypotheses(guide, station)
        assert distance_admissible_boundary_sites(guide, station)

    assert calls == 1


def test_wired_analysis_rejects_a_station_the_guide_does_not_own() -> None:
    guide = build_figure5_raw_guide(load_held_reference_case("figure5"))
    station = guide.runs[0].stations[0]
    foreign = replace(station, ordinal_on_run=GuideRunStationOrdinal(1))

    with pytest.raises(InvalidFigure5RawGuideError):
        build_distance_admissible_hypotheses(guide, foreign)


# Quarter-unit lattice coordinates make exact ties, vertex projections and
# zero-gap boundaries common; generic floats cover the rest.
_LATTICE = st.integers(min_value=-16, max_value=16).map(lambda value: value / 4)
_GENERIC = st.floats(min_value=-4.0, max_value=4.0, allow_nan=False, allow_infinity=False)
_COORDINATE = st.one_of(_LATTICE, _GENERIC)
_GAP = st.one_of(st.just(0.0), st.sampled_from((0.25, 0.5, 1.0)), st.floats(min_value=0.0, max_value=2.0, allow_nan=False))


def _ring_without_zero_length_segments(points: list[tuple[float, float]]) -> bool:
    return all(start != end for start, end in zip(points, (*points[1:], points[0])))


def _synthetic_guide(ring: list[tuple[float, float]], stations: list[tuple[float, float]], gap: float) -> Figure5RawGuide:
    """Bypass the Held-case factory: both analyses read only these fields."""
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ring)
    raw_stations = tuple(Figure5RawGuideStation(0, ordinal, Point2[WorldXY].build(x, y)) for ordinal, (x, y) in enumerate(stations))  # type: ignore[arg-type]
    case = SimpleNamespace(projection=SimpleNamespace(points=boundary))
    budget = Figure5BoundarySiteBudget(Millimetre(0.0), Millimetre(0.0), Millimetre(gap))
    return Figure5RawGuide(case, budget, Millimetre(1.0), (Figure5RawGuideRun(0, raw_stations),))  # type: ignore[arg-type]


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.filter_too_much])
@given(
    ring=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=3, max_size=9).filter(_ring_without_zero_length_segments),
    stations=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=1, max_size=6),
    gap=_GAP,
)
def test_wired_analysis_matches_fraction_oracle_on_random_rings(
    ring: list[tuple[float, float]],
    stations: list[tuple[float, float]],
    gap: float,
) -> None:
    _assert_wired_matches_oracle(_synthetic_guide(ring, stations, gap), build_hypotheses=False)
