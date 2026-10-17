"""Native distance-admissible boundary sites decide exactly what the Fraction path decides."""

from __future__ import annotations

from fractions import Fraction
from types import SimpleNamespace

import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_raw_guide import Figure5BoundarySiteBudget
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import Figure5RawGuideRun
from benchmarks.held_figure5_raw_guide import Figure5RawGuideStation
from benchmarks.held_figure5_raw_guide import _analyze_station_boundary
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_figure5_raw_guide_native import NativeStationBoundaryAnalysis
from benchmarks.held_figure5_raw_guide_native import native_station_boundary_analyses
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import _stock_2
from compas_cgal.held_boundary_clearance import HeldBoundaryClearanceIndex
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


def _assert_equivalent(guide: Figure5RawGuide, native: tuple[tuple[NativeStationBoundaryAnalysis, ...], ...]) -> int:
    """Compare every station's decisions and admissible exact records; return the site count."""
    assert len(native) == len(guide.runs)
    compared = 0
    for run, native_run in zip(guide.runs, native):
        assert len(native_run) == len(run.stations)
        for station, analysis in zip(run.stations, native_run):
            old = _analyze_station_boundary(guide, station)
            assert analysis.distance_admissible_sites == old.distance_admissible_sites, (station.run_id, station.ordinal_on_run)
            old_records = {site: (footpoint, distance) for site, footpoint, distance in old.records}
            for site, footpoint, distance in analysis.admissible_records:
                assert old_records[site] == (footpoint, distance), (station.run_id, station.ordinal_on_run, site)
            assert tuple(record[0] for record in analysis.admissible_records) == analysis.distance_admissible_sites
            compared += len(analysis.distance_admissible_sites)
    return compared


def test_native_matches_fraction_decisions_on_every_figure5_station() -> None:
    guide = build_held_reference_raw_guide(load_held_reference_case("figure5"))

    compared = _assert_equivalent(guide, native_station_boundary_analyses(guide))

    assert guide.station_count == 7075
    assert compared > guide.station_count


def test_degenerate_boundary_segment_is_rejected_at_build() -> None:
    ring = (Point2[WorldXY].build(0.0, 0.0), Point2[WorldXY].build(0.0, 0.0), Point2[WorldXY].build(1.0, 1.0))
    with pytest.raises(_stock_2.DegenerateHeldBoundarySegmentError):
        HeldBoundaryClearanceIndex.build(ring)


def test_negative_gap_is_rejected() -> None:
    ring = (Point2[WorldXY].build(0.0, 0.0), Point2[WorldXY].build(1.0, 0.0), Point2[WorldXY].build(0.0, 1.0))
    index = HeldBoundaryClearanceIndex.build(ring)
    with pytest.raises(_stock_2.InvalidHeldBoundaryClearanceInputError):
        index.distance_admissible_sites((Point2[WorldXY].build(0.2, 0.2),), Millimetre(-1.0))


# Quarter-unit integer lattice coordinates make exact distance ties, vertex
# projections and zero-gap boundaries common; generic floats cover the rest.
_LATTICE = st.integers(min_value=-16, max_value=16).map(lambda value: value / 4)
_GENERIC = st.floats(min_value=-4.0, max_value=4.0, allow_nan=False, allow_infinity=False)
_COORDINATE = st.one_of(_LATTICE, _GENERIC)
_GAP = st.one_of(st.just(0.0), st.sampled_from((0.25, 0.5, 1.0)), st.floats(min_value=0.0, max_value=2.0, allow_nan=False))


def _ring_without_zero_length_segments(points: list[tuple[float, float]]) -> bool:
    return all(start != end for start, end in zip(points, (*points[1:], points[0])))


def _synthetic_guide(ring: list[tuple[float, float]], stations: list[tuple[float, float]], gap: float) -> Figure5RawGuide:
    """Bypass the Held-case factory: the old analysis reads only these fields."""
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ring)
    raw_stations = tuple(Figure5RawGuideStation(0, ordinal, Point2[WorldXY].build(x, y)) for ordinal, (x, y) in enumerate(stations))
    case = SimpleNamespace(projection=SimpleNamespace(points=boundary))
    budget = Figure5BoundarySiteBudget(Millimetre(0.0), Millimetre(0.0), Millimetre(gap))
    return Figure5RawGuide(case, budget, Millimetre(1.0), (Figure5RawGuideRun(0, raw_stations),))  # type: ignore[arg-type]


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.filter_too_much])
@given(
    ring=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=3, max_size=9).filter(_ring_without_zero_length_segments),
    stations=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=1, max_size=6),
    gap=_GAP,
)
def test_native_matches_fraction_decisions_on_random_rings(
    ring: list[tuple[float, float]],
    stations: list[tuple[float, float]],
    gap: float,
) -> None:
    guide = _synthetic_guide(ring, stations, gap)

    _assert_equivalent(guide, native_station_boundary_analyses(guide))


def test_exact_tie_reports_both_sites() -> None:
    ring = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
    guide = _synthetic_guide(ring, [(2.0, 1.0), (1.0, 1.0), (2.0, 2.0)], 0.0)

    native = native_station_boundary_analyses(guide)

    _assert_equivalent(guide, native)
    assert [len(analysis.distance_admissible_sites) for analysis in native[0]] == [1, 2, 4]
    assert native[0][1].admissible_records[0][2] == Fraction(1)
