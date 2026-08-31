from __future__ import annotations

import math
import sys

import pytest

from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import InvalidReferenceProjectionError
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import NORMALIZED_ROUNDOFF
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PolygonProjection
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import project_boundary
from benchmarks.held_reference_geometry import reconstruct_cubic
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _world(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


def _identity_transform() -> SourceToWorld:
    return SourceToWorld.build(
        source_origin=PdfPoint2.build(0.0, 0.0),
        world_origin=_world(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0),
    )


def _quarter_circle_source() -> SourceCubic:
    control = 4.0 * (math.sqrt(2.0) - 1.0) / 3.0
    return SourceCubic.build(
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(1.0, control),
        PdfPoint2.build(control, 1.0),
        PdfPoint2.build(0.0, 1.0),
    )


def _non_circular_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(0.0, 0.0),
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(1.0, 2.0),
        PdfPoint2.build(2.0, 2.0),
    )


def _asymmetric_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(0.0, 0.0),
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(2.0, 1.0),
        PdfPoint2.build(2.0, 3.0),
    )


def _semicircle_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(-1.0, 0.0),
        PdfPoint2.build(-1.0, 4.0 / 3.0),
        PdfPoint2.build(1.0, 4.0 / 3.0),
        PdfPoint2.build(1.0, 0.0),
    )


def _angular_reversal_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(1.0, 5.0),
        PdfPoint2.build(5.0, 1.0),
        PdfPoint2.build(0.0, 1.0),
    )


def _unit_tangent(arc: ReferenceArc, *, at_end: bool) -> tuple[float, float]:
    point = arc.end if at_end else arc.start
    radius_x = float(point.x) - float(arc.centre.x)
    radius_y = float(point.y) - float(arc.centre.y)
    direction = 1.0 if float(arc.sweep) > 0.0 else -1.0
    length = math.hypot(radius_x, radius_y)
    return direction * -radius_y / length, direction * radius_x / length


def _vector_residual(left: tuple[float, float], right: tuple[float, float]) -> float:
    return math.hypot(left[0] - right[0], left[1] - right[1])


def _rounded_rectangle_boundary() -> ReferenceBoundary:
    primitives = (
        ReferenceLine.build(_world(1.0, 0.0), _world(3.0, 0.0)),
        ReferenceArc.build(_world(3.0, 0.0), _world(4.0, 1.0), _world(3.0, 1.0), Radian(math.pi / 2.0)),
        ReferenceLine.build(_world(4.0, 1.0), _world(4.0, 3.0)),
        ReferenceArc.build(_world(4.0, 3.0), _world(3.0, 4.0), _world(3.0, 3.0), Radian(math.pi / 2.0)),
        ReferenceLine.build(_world(3.0, 4.0), _world(1.0, 4.0)),
        ReferenceArc.build(_world(1.0, 4.0), _world(0.0, 3.0), _world(1.0, 3.0), Radian(math.pi / 2.0)),
        ReferenceLine.build(_world(0.0, 3.0), _world(0.0, 1.0)),
        ReferenceArc.build(_world(0.0, 1.0), _world(1.0, 0.0), _world(1.0, 1.0), Radian(math.pi / 2.0)),
    )
    return ReferenceBoundary.build(primitives, ToolRadius.build(1.0), Millimetre(0.1))


def _unit_circle_boundary(offset: float = 0.0) -> ReferenceBoundary:
    centre = _world(offset, offset)
    cardinal_points = (
        _world(offset + 1.0, offset),
        _world(offset, offset + 1.0),
        _world(offset - 1.0, offset),
        _world(offset, offset - 1.0),
    )
    arcs = tuple(ReferenceArc.build(start, end, centre, Radian(math.pi / 2.0)) for start, end in zip(cardinal_points, (*cardinal_points[1:], cardinal_points[0])))
    return ReferenceBoundary.build(arcs, ToolRadius.build(1.0), Millimetre(0.1))


def _point_segment_distance(
    point: tuple[float, float],
    start: tuple[float, float],
    end: tuple[float, float],
) -> float:
    segment = end[0] - start[0], end[1] - start[1]
    parameter = ((point[0] - start[0]) * segment[0] + (point[1] - start[1]) * segment[1]) / (segment[0] ** 2 + segment[1] ** 2)
    parameter = min(1.0, max(0.0, parameter))
    closest = start[0] + parameter * segment[0], start[1] + parameter * segment[1]
    return math.hypot(point[0] - closest[0], point[1] - closest[1])


def test_pdf_point_rejects_non_finite_coordinate() -> None:
    with pytest.raises(InvalidPublishedPrimitiveError):
        PdfPoint2.build(float("nan"), 0.0)


def test_source_line_rejects_identical_endpoints() -> None:
    point = PdfPoint2.build(2.0, 3.0)
    with pytest.raises(InvalidPublishedPrimitiveError):
        SourceLine.build(point, point)


def test_reference_boundary_rejects_disconnected_primitives() -> None:
    first = ReferenceLine.build(_world(0.0, 0.0), _world(1.0, 0.0))
    second = ReferenceLine.build(_world(2.0, 0.0), _world(0.0, 0.0))

    with pytest.raises(DisconnectedPublishedBoundaryError):
        ReferenceBoundary.build(
            (first, second),
            ToolRadius.build(1.0),
            Millimetre(0.1),
        )


@pytest.mark.parametrize("offset", (0.0, 2.0**30))
def test_reference_arc_rejects_local_radius_defect_independent_of_translation(offset: float) -> None:
    radius = 2.0**-30
    relative_defect = math.sqrt(sys.float_info.epsilon)

    with pytest.raises(InvalidPublishedPrimitiveError):
        ReferenceArc.build(
            _world(offset, radius),
            _world(offset, -radius * (1.0 + relative_defect)),
            _world(offset, 0.0),
            Radian(math.pi),
        )


def test_quarter_circle_cubic_recovers_one_arc() -> None:
    arcs = reconstruct_cubic(
        _quarter_circle_source(),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert len(arcs) == 1
    assert float(arcs[0].centre.x) == pytest.approx(0.0)
    assert float(arcs[0].centre.y) == pytest.approx(0.0)
    assert float(arcs[0].sweep) == pytest.approx(math.pi / 2.0)


def test_non_circular_cubic_splits_into_g1_arcs() -> None:
    arcs = reconstruct_cubic(
        _non_circular_source(),
        _identity_transform(),
        Millimetre(0.01),
    )

    assert len(arcs) >= 2
    for left, right in zip(arcs, arcs[1:]):
        assert left.end == right.start
        assert _unit_tangent(left, at_end=True) == pytest.approx(_unit_tangent(right, at_end=False))


def test_asymmetric_cubic_reconstructs_without_unequal_radius_failure() -> None:
    arcs = reconstruct_cubic(
        _asymmetric_source(),
        _identity_transform(),
        Millimetre(0.01),
    )

    assert arcs[0].start == _world(0.0, 0.0)
    assert arcs[-1].end == _world(2.0, 3.0)
    for left, right in zip(arcs, arcs[1:]):
        assert _unit_tangent(left, at_end=True) == pytest.approx(_unit_tangent(right, at_end=False))


def test_antiparallel_endpoint_tangents_recover_one_semicircle() -> None:
    arcs = reconstruct_cubic(
        _semicircle_source(),
        _identity_transform(),
        Millimetre(0.02),
    )

    assert len(arcs) == 1
    assert float(arcs[0].centre.x) == pytest.approx(0.0)
    assert float(arcs[0].centre.y) == pytest.approx(0.0)
    assert float(arcs[0].sweep) == pytest.approx(-math.pi)


def test_near_parallel_biarc_preserves_both_authored_endpoint_tangents() -> None:
    tangent_x = math.sqrt(32.0 * sys.float_info.epsilon)
    start_tangent = (tangent_x, math.sqrt(1.0 - tangent_x**2))
    end_tangent = (0.0, 1.0)
    source = SourceCubic.build(
        PdfPoint2.build(-1.0, 0.0),
        PdfPoint2.build(-1.0 + start_tangent[0], start_tangent[1]),
        PdfPoint2.build(1.0 - end_tangent[0], -end_tangent[1]),
        PdfPoint2.build(1.0, 0.0),
    )

    arcs = reconstruct_cubic(source, _identity_transform(), Millimetre(2.0))

    assert _vector_residual(_unit_tangent(arcs[0], at_end=False), start_tangent) <= NORMALIZED_ROUNDOFF
    assert _vector_residual(_unit_tangent(arcs[-1], at_end=True), end_tangent) <= NORMALIZED_ROUNDOFF
    for left, right in zip(arcs, arcs[1:]):
        assert _vector_residual(_unit_tangent(left, at_end=True), _unit_tangent(right, at_end=False)) <= NORMALIZED_ROUNDOFF


def test_circle_candidate_rejects_angular_reversal() -> None:
    arcs = reconstruct_cubic(
        _angular_reversal_source(),
        _identity_transform(),
        Millimetre(10.0),
    )

    assert len(arcs) >= 2


def test_projection_closes_and_meets_chord_bound() -> None:
    projection = project_boundary(_rounded_rectangle_boundary(), Millimetre(0.002))

    assert isinstance(projection, PolygonProjection)
    assert projection.points[0] != projection.points[-1]
    assert projection.observed_deviation <= projection.deviation_limit
    assert len(projection.points) >= 8


def test_projection_reports_independently_measured_emitted_chord_deviation() -> None:
    projection = project_boundary(_unit_circle_boundary(), Millimetre(0.002))
    segments_per_quarter = len(projection.points) // 4
    angle_step = math.pi / (2.0 * segments_per_quarter)
    oracle = 0.0
    for index, start in enumerate(projection.points):
        end = projection.points[(index + 1) % len(projection.points)]
        midpoint_angle = (index + 0.5) * angle_step
        oracle = max(
            oracle,
            _point_segment_distance(
                (math.cos(midpoint_angle), math.sin(midpoint_angle)),
                (float(start.x), float(start.y)),
                (float(end.x), float(end.y)),
            ),
        )

    assert float(projection.observed_deviation) == pytest.approx(oracle, rel=0.0, abs=sys.float_info.epsilon)


def test_projection_rejects_translated_circle_when_emitted_coordinates_exceed_bound() -> None:
    offset = 2.0**40

    with pytest.raises(InvalidReferenceProjectionError):
        project_boundary(_unit_circle_boundary(offset), Millimetre(math.ulp(offset)))


def test_projection_rejects_non_positive_bound() -> None:
    with pytest.raises(InvalidReferenceProjectionError):
        project_boundary(_rounded_rectangle_boundary(), Millimetre(0.0))
