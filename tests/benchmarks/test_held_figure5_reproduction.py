import math
from fractions import Fraction

from compas.geometry import Circle
from compas.tolerance import TOL

from benchmarks.held_figure5_reproduction import APPROXIMATION_PROVENANCE
from benchmarks.held_figure5_reproduction import _boundary_footpoint
from benchmarks.held_figure5_reproduction import build_figure5_station_approximation
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


def _distance(first: Point2[WorldXY], second: Point2[WorldXY]) -> float:
    return math.hypot(
        float(second.x) - float(first.x),
        float(second.y) - float(first.y),
    )


def test_equal_distance_boundary_tie_follows_current_circle_phase() -> None:
    left = (
        (Fraction(-1), Fraction(-1)),
        (Fraction(-1), Fraction(1)),
    )
    right = (
        (Fraction(1), Fraction(-1)),
        (Fraction(1), Fraction(1)),
    )

    footpoint, squared_distance = _boundary_footpoint(
        middle_point=(Fraction(0), Fraction(0)),
        phase_point=(Fraction(1), Fraction(0)),
        boundary_segments=(left, right),
    )

    assert footpoint == (Fraction(1), Fraction(0))
    assert squared_distance == Fraction(1)


def test_figure5_adapter_labels_stations_and_preserves_paper_geometry() -> None:
    case = load_held_reference_case("figure5")
    approximation = build_figure5_station_approximation(case)

    assert approximation.provenance == APPROXIMATION_PROVENANCE
    assert approximation.provenance == "straight-skeleton station approximation"
    assert approximation.stations
    source_stations = tuple(operation for operation in approximation.source_path.operations if isinstance(operation.geometry, Circle))
    assert len(approximation.stations) == len(source_stations)
    assert tuple(station.path_index for station in approximation.stations) == tuple(operation.path_index for operation in source_stations)

    tool_radius = float(case.tool_radius.value)
    for station in approximation.stations:
        assert TOL.is_between(
            _distance(station.boundary_footpoint, station.contact_point),
            tool_radius,
            tool_radius,
        )
        assert TOL.is_between(
            _distance(station.contact_point, station.center),
            float(station.guide_radius),
            float(station.guide_radius),
        )
        expected_diameter = _distance(station.boundary_footpoint, station.middle_point) - tool_radius
        assert TOL.is_between(
            2.0 * float(station.guide_radius),
            expected_diameter,
            expected_diameter,
        )
