from __future__ import annotations

import math

import pytest

from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import InvalidReferenceProjectionError
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
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


def test_projection_rejects_non_positive_bound() -> None:
    with pytest.raises(InvalidReferenceProjectionError):
        project_boundary(_rounded_rectangle_boundary(), Millimetre(0.0))
