"""The trimmed bulk hypothesis path builds exactly what the per-site path built.

The oracle is the pre-change ``_build_projection_admissible_hypothesis`` and
``build_distance_admissible_hypotheses``, kept verbatim. The bulk path reads
each native record directly instead of re-finding it by dataclass equality, and
injects the tool radius once per station; outputs must be equal field for field
(floats with ``==``), and the exact ``clearance**2 <= r**2`` rejection must fire
on the same station with the same error.
"""

from __future__ import annotations

import math
from fractions import Fraction
from types import SimpleNamespace

from hypothesis import HealthCheck
from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_raw_guide import Figure5BoundarySiteBudget
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import Figure5RawGuideRun
from benchmarks.held_figure5_raw_guide import Figure5RawGuideStation
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import InadmissibleFigure5BoundarySideError
from benchmarks.held_figure5_raw_guide import InvalidProjectionAdmissibleHypothesisError
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import _analyze_station_boundary
from benchmarks.held_figure5_raw_guide import _boundary_arclength
from benchmarks.held_figure5_raw_guide import _build_projection_admissible_hypothesis
from benchmarks.held_figure5_raw_guide import _exact_square
from benchmarks.held_figure5_raw_guide import _native_station_boundary
from benchmarks.held_figure5_raw_guide import _StationBoundaryAnalysis
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from compas_cgal.errors import InvalidUnitValueError
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

Outcome = tuple[ProjectionAdmissibleBoundaryHypothesis, ...] | tuple[type[Exception], str]


def _old_build_projection_admissible_hypothesis(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
    site: ProjectionBoundarySite,
    analysis: _StationBoundaryAnalysis,
) -> ProjectionAdmissibleBoundaryHypothesis:
    """Verbatim pre-change ``_build_projection_admissible_hypothesis``."""
    admissible = analysis.distance_admissible_sites
    if site not in admissible:
        raise InadmissibleFigure5BoundarySideError(f"Projected boundary site {int(site.segment_id)}:{site.parameter} is outside the Figure 5 evidence budget.")
    footpoint_exact, squared_clearance = next((footpoint, distance) for record_site, footpoint, distance in analysis.records if record_site == site)

    tool_radius = float(guide.case.tool_radius.value)
    exact_tool_radius = Fraction.from_float(tool_radius)
    if squared_clearance <= exact_tool_radius * exact_tool_radius:
        raise InvalidProjectionAdmissibleHypothesisError("Figure 5 projected-boundary hypothesis has no positive paper-derived radius.")
    middle = station.middle_point
    footpoint = Point2[WorldXY].build(
        float(footpoint_exact[0]),
        float(footpoint_exact[1]),
    )
    clearance = math.sqrt(float(squared_clearance))
    direction_x = (float(middle.x) - float(footpoint.x)) / clearance
    direction_y = (float(middle.y) - float(footpoint.y)) / clearance
    radius = (clearance - tool_radius) / 2.0
    contact = Point2[WorldXY].build(
        float(footpoint.x) + tool_radius * direction_x,
        float(footpoint.y) + tool_radius * direction_y,
    )
    center = Point2[WorldXY].build(
        float(contact.x) + radius * direction_x,
        float(contact.y) + radius * direction_y,
    )
    return ProjectionAdmissibleBoundaryHypothesis(
        station.run_id,
        station.ordinal_on_run,
        site,
        admissible,
        _boundary_arclength(guide, site),
        middle,
        footpoint,
        contact,
        center,
        GuideRadius.build(radius),
    )


def _old_build_distance_admissible_hypotheses(
    guide: Figure5RawGuide,
    station: Figure5RawGuideStation,
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    """Verbatim pre-change ``build_distance_admissible_hypotheses``."""
    analysis = _native_station_boundary(guide, station)
    return tuple(
        _old_build_projection_admissible_hypothesis(
            guide,
            station,
            site,
            analysis,
        )
        for site in analysis.distance_admissible_sites
    )


# Quarter-unit lattice coordinates make exact ties and r**2 == clearance**2 common.
_LATTICE = st.integers(min_value=-16, max_value=16).map(lambda value: value / 4)
_GENERIC = st.floats(min_value=-4.0, max_value=4.0, allow_nan=False, allow_infinity=False)
_COORDINATE = st.one_of(_LATTICE, _GENERIC)
_GAP = st.one_of(st.just(0.0), st.sampled_from((0.25, 0.5, 1.0)), st.floats(min_value=0.0, max_value=2.0, allow_nan=False))
_TOOL_RADIUS = st.one_of(st.sampled_from((0.25, 0.5, 1.0)), st.floats(min_value=1e-3, max_value=2.0, allow_nan=False))


def _ring_without_zero_length_segments(points: list[tuple[float, float]]) -> bool:
    return all(start != end for start, end in zip(points, (*points[1:], points[0])))


def _synthetic_guide(ring: list[tuple[float, float]], stations: list[tuple[float, float]], gap: float, tool_radius: float) -> Figure5RawGuide:
    """Bypass the Held-case factory: hypothesis construction reads only these fields."""
    boundary = tuple(Point2[WorldXY].build(x, y) for x, y in ring)
    raw_stations = tuple(Figure5RawGuideStation(GuideRunId(0), GuideRunStationOrdinal(ordinal), Point2[WorldXY].build(x, y)) for ordinal, (x, y) in enumerate(stations))
    case = SimpleNamespace(projection=SimpleNamespace(points=boundary), tool_radius=SimpleNamespace(value=Millimetre(tool_radius)))
    budget = Figure5BoundarySiteBudget(Millimetre(0.0), Millimetre(0.0), Millimetre(gap))
    length = sum(math.hypot(end[0] - start[0], end[1] - start[1]) for start, end in zip(ring, (*ring[1:], ring[0])))
    return Figure5RawGuide(case, budget, Millimetre(length), (Figure5RawGuideRun(GuideRunId(0), raw_stations),))  # type: ignore[arg-type]


def _outcome(build: object, guide: Figure5RawGuide, station: Figure5RawGuideStation) -> Outcome:
    try:
        return build(guide, station)  # type: ignore[operator, no-any-return]
    # A clearance a hair above r rounds to a non-positive double radius: both
    # paths reject it in ``GuideRadius.build``, and must do so identically.
    except (InvalidProjectionAdmissibleHypothesisError, InadmissibleFigure5BoundarySideError, InvalidUnitValueError) as error:
        return type(error), str(error)


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.filter_too_much])
@given(
    ring=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=3, max_size=9).filter(_ring_without_zero_length_segments),
    stations=st.lists(st.tuples(_COORDINATE, _COORDINATE), min_size=1, max_size=6),
    gap=_GAP,
    tool_radius=_TOOL_RADIUS,
)
@example(ring=[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)], stations=[(2.0, 1.0), (1.0, 1.0), (2.0, 2.0)], gap=0.0, tool_radius=1.0)
@example(ring=[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)], stations=[(2.0, 2.0), (3.0, 2.5)], gap=0.5, tool_radius=0.25)
def test_bulk_hypotheses_equal_the_pre_change_per_site_path(
    ring: list[tuple[float, float]],
    stations: list[tuple[float, float]],
    gap: float,
    tool_radius: float,
) -> None:
    guide = _synthetic_guide(ring, stations, gap, tool_radius)
    for station in guide.runs[0].stations:
        new = _outcome(build_distance_admissible_hypotheses, guide, station)
        old = _outcome(_old_build_distance_admissible_hypotheses, guide, station)
        assert new == old
        if isinstance(new, tuple) and new and isinstance(new[0], ProjectionAdmissibleBoundaryHypothesis):
            assert repr(new) == repr(old)
        oracle = _analyze_station_boundary(guide, station)
        for site in oracle.records[0][0], *oracle.distance_admissible_sites:
            single_new = _outcome(lambda g, s, site=site: _build_projection_admissible_hypothesis(g, s, site, oracle), guide, station)
            single_old = _outcome(lambda g, s, site=site: _old_build_projection_admissible_hypothesis(g, s, site, oracle), guide, station)
            assert single_new == single_old


@settings(max_examples=2000, deadline=None)
@given(st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
@example(0.0)
@example(math.ulp(0.0))
@example(1.0 + math.ulp(1.0))
def test_exact_square_is_the_fraction_square(value: float) -> None:
    exact = Fraction.from_float(value)
    assert _exact_square(value) == exact * exact
