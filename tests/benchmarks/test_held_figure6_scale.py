from __future__ import annotations

import pytest

from benchmarks.held_figure5_publisher import Figure5PublisherCubicEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherLineEvidence
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_figure6_publisher import load_figure6_publisher_evidence
from benchmarks.held_figure6_scale import FIGURE5_STANDARD_CAP_DEG
from benchmarks.held_figure6_scale import calibrate_publisher_figure6_length
from benchmarks.held_figure6_scale import log_interpolated_graphical_length
from benchmarks.held_figure6_scale import publisher_figure5_planar_length
from benchmarks.held_figure6_scale import publisher_primitive_length
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


def test_publisher_primitive_length_measures_lines_analytically() -> None:
    primitive = Figure5PublisherLineEvidence.build(0, _point(0.0, 0.0), _point(3.0, 4.0))

    assert float(publisher_primitive_length(primitive)) == pytest.approx(5.0)


def test_publisher_primitive_length_integrates_a_straight_cubic() -> None:
    primitive = Figure5PublisherCubicEvidence.build(
        0,
        _point(0.0, 0.0),
        _point(1.0, 0.0),
        _point(2.0, 0.0),
        _point(3.0, 0.0),
    )

    assert float(publisher_primitive_length(primitive)) == pytest.approx(3.0)


def test_tracked_figure5_stream_has_dimensionless_planar_length() -> None:
    case = load_held_reference_case("figure5")
    path = load_figure5_publisher_path()

    length = publisher_figure5_planar_length(path, case.tool_radius)

    assert float(length) == pytest.approx(9412.051197363458, rel=1e-12)


def test_standard_figure6_series_is_log_interpolated_at_figure5_cap() -> None:
    standard = load_figure6_publisher_evidence().series[0]

    length = log_interpolated_graphical_length(standard, FIGURE5_STANDARD_CAP_DEG)

    assert float(length) == pytest.approx(221.9140179195166, rel=1e-12)


def test_publisher_length_calibration_anchors_figure6_to_figure5() -> None:
    case = load_held_reference_case("figure5")
    path = load_figure5_publisher_path()
    figure6 = load_figure6_publisher_evidence()

    calibration = calibrate_publisher_figure6_length(case, path, figure6)

    assert float(calibration.figure5_path_length) == pytest.approx(9412.051197363458, rel=1e-12)
    assert float(calibration.figure6_graphical_length) == pytest.approx(221.9140179195166, rel=1e-12)
    assert calibration.tool_radius_multiples_per_graphical_unit == pytest.approx(42.41305387376207, rel=1e-12)
    assert calibration.to_tool_radius_multiples(calibration.figure6_graphical_length) == pytest.approx(calibration.figure5_path_length)
