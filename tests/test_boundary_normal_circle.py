"""Native boundary-normal proposals use the true first medial contact."""

import math

import pytest

from compas_cgal import _circle_geometry_2 as native


RECTANGLE = [(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)]
REPORT_ULP_BUDGET = 4  # CORE reports approximate binary64; geometry is checked exactly before reporting.


def _assert_report(actual: float, expected: float) -> None:
    assert actual == pytest.approx(expected, rel=0.0, abs=REPORT_ULP_BUDGET * math.ulp(expected))


def _assert_xy_report(actual: tuple[float, float], expected: tuple[float, float]) -> None:
    for value, reference in zip(actual, expected):
        _assert_report(value, reference)


def test_rectangle_query_constructs_medial_contact_and_q_m_diameter() -> None:
    owner = native.BoundaryNormalCircle2(RECTANGLE)
    circle = owner.query(0, 0.5, 1.0)
    assert circle.p_mm == pytest.approx((5.0, 0.0))
    assert circle.m_mm == pytest.approx((5.0, 3.0))
    assert circle.q_mm == pytest.approx((5.0, 1.0))
    assert circle.center_mm == pytest.approx((5.0, 2.0))
    _assert_report(circle.guide_radius_mm, 1.0)
    _assert_report(circle.clearance_mm, 3.0)
    assert circle.competing_segment_indices == [2]


def test_first_medial_contact_changes_near_rectangle_corner() -> None:
    circle = native.BoundaryNormalCircle2(RECTANGLE).query(0, 0.125, 0.25)
    assert circle.m_mm == pytest.approx((1.25, 1.25))
    _assert_report(circle.clearance_mm, 1.25)
    assert circle.competing_segment_indices == [3]


def test_sloping_segment_competitor_uses_exact_line_distance() -> None:
    triangle = [(0.0, 0.0), (8.0, 0.0), (0.0, 6.0)]
    circle = native.BoundaryNormalCircle2(triangle).query(0, 0.375, 0.25)
    assert circle.m_mm == pytest.approx((3.0, 5.0 / 3.0))
    assert circle.clearance_mm == pytest.approx(5.0 / 3.0)
    assert circle.competing_segment_indices == [1]


def test_reflex_point_wins_when_support_line_foot_is_outside_segment() -> None:
    polygon = [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (2.0, 2.0), (2.0, 6.0), (0.0, 6.0)]
    circle = native.BoundaryNormalCircle2(polygon).query(0, 0.25, 0.25)
    assert circle.p_mm == (1.5, 0.0)
    _assert_xy_report(circle.m_mm, (1.5, 1.0625))
    _assert_report(circle.clearance_mm, 1.0625)
    assert circle.competing_vertex_indices == [3]
    assert circle.competing_segment_indices == []
    _assert_xy_report(circle.center_mm, (1.5, 0.65625))
    _assert_report(circle.guide_radius_mm, 0.40625)


@pytest.mark.parametrize("scale", [1.0 / 1024, 1.0, 1024.0])
def test_binary_scale_preserves_proposal_geometry(scale: float) -> None:
    owner = native.BoundaryNormalCircle2([(x * scale, y * scale) for x, y in RECTANGLE])
    circle = owner.query(0, 0.5, scale)
    assert circle.m_mm == pytest.approx((5 * scale, 3 * scale))
    _assert_report(circle.guide_radius_mm, scale)


def test_clockwise_input_preserves_original_segment_index() -> None:
    owner = native.BoundaryNormalCircle2([(0.0, 0.0), (0.0, 6.0), (10.0, 6.0), (10.0, 0.0)])
    _assert_xy_report(owner.query(3, 0.5, 1.0).m_mm, (5.0, 3.0))


@pytest.mark.parametrize("parameter", [0.0, 1.0])
def test_vertex_query_rejected_explicitly(parameter: float) -> None:
    with pytest.raises(native.BoundaryVertexQueryUnsupportedError):
        native.BoundaryNormalCircle2(RECTANGLE).query(0, parameter, 1.0)


@pytest.mark.parametrize("index,parameter,radius", [(-1, 0.5, 1.0), (4, 0.5, 1.0), (0, -0.1, 1.0), (0, float("nan"), 1.0), (0, 0.5, 0.0)])
def test_invalid_query_has_named_error(index: int, parameter: float, radius: float) -> None:
    with pytest.raises(native.InvalidBoundaryNormalInputError):
        native.BoundaryNormalCircle2(RECTANGLE).query(index, parameter, radius)


@pytest.mark.parametrize("radius", [3.0, 4.0])
def test_no_positive_machining_radius_fails(radius: float) -> None:
    with pytest.raises(native.NoPositiveBoundaryCircleError):
        native.BoundaryNormalCircle2(RECTANGLE).query(0, 0.5, radius)


@pytest.mark.parametrize("polygon", [[(0.0, 0.0), (1.0, 0.0)], [(0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)], [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)]])
def test_invalid_polygon_has_named_error(polygon: list[tuple[float, float]]) -> None:
    with pytest.raises(native.InvalidBoundaryPolygonError):
        native.BoundaryNormalCircle2(polygon)


def test_reflex_vertex_sector_queries_true_medial_branch() -> None:
    polygon = [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (2.0, 2.0), (2.0, 6.0), (0.0, 6.0)]
    owner = native.BoundaryNormalCircle2(polygon)
    circle = owner.query_vertex(3, (-1.0, -1.0), 0.25)
    coordinate = 2.0 * 2.0**0.5 / (1.0 + 2.0**0.5)
    assert circle.m_mm == pytest.approx((coordinate, coordinate))
    assert circle.clearance_mm == pytest.approx(coordinate)
    assert circle.competing_segment_indices == [0, 5]
    assert circle.competing_vertex_indices == []
    assert circle.p_mm == (2.0, 2.0)


@pytest.mark.parametrize("direction,expected", [((-1.0, 0.0), (1.0, 2.0)), ((0.0, -1.0), (2.0, 1.0))])
def test_reflex_sector_endpoints_ignore_same_contact_ties(direction: tuple[float, float], expected: tuple[float, float]) -> None:
    polygon = [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (2.0, 2.0), (2.0, 6.0), (0.0, 6.0)]
    circle = native.BoundaryNormalCircle2(polygon).query_vertex(3, direction, 0.25)
    assert circle.m_mm == expected
    assert circle.clearance_mm == 1.0


@pytest.mark.parametrize("vertex,direction", [(3, (1.0, 0.0)), (3, (-1.0, 1.0)), (3, (0.0, 0.0)), (0, (-1.0, -1.0))])
def test_invalid_reflex_normal_cone_has_named_error(vertex: int, direction: tuple[float, float]) -> None:
    polygon = [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (2.0, 2.0), (2.0, 6.0), (0.0, 6.0)]
    with pytest.raises(native.InvalidBoundaryVertexDirectionError):
        native.BoundaryNormalCircle2(polygon).query_vertex(vertex, direction, 0.25)
