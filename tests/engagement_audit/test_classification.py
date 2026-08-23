from __future__ import annotations

import math

import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

from compas_cgal import _stock_2
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.engagement_audit.classification import classify_operation
from compas_cgal.engagement_audit.errors import ContradictoryOperationRoleError
from compas_cgal.engagement_audit.errors import ContradictoryOperationOrientationError
from compas_cgal.engagement_audit.errors import UnsupportedAuditGeometryError
from compas_cgal.engagement_audit.records import AuthenticatedPlungeOperation
from compas_cgal.engagement_audit.records import AuthenticatedLateralOperation
from compas_cgal.engagement_audit.records import AuthenticatedNonEngagingOperation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

CUT_Z = 0.0
CLEARANCE_Z = 5.0
OPERATION_INDEX = 7


def _cut_plane() -> CutPlane:
    return CutPlane.build(CutZ.build(CUT_Z), ClearanceZ.build(CLEARANCE_Z))


def _line(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    operation: OperationType,
) -> ToolpathOperation:
    return ToolpathOperation(
        geometry=Line(start, end),
        operation=operation,
        path_index=3,
    )


@pytest.mark.parametrize(
    ("start", "end", "operation", "native_type"),
    [
        ((1.0, 2.0, CUT_Z), (1.0, 2.0, CLEARANCE_Z), OperationType.RETRACT, _stock_2.AuditVerticalRetract2),
        ((1.0, 2.0, CLEARANCE_Z), (3.0, 2.0, CLEARANCE_Z), OperationType.LINK, _stock_2.AuditClearanceTransport2),
    ],
)
def test_geometry_proves_non_engaging_motion(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    operation: OperationType,
    native_type: type[object],
) -> None:
    result = classify_operation(
        _line(start, end, operation),
        _cut_plane(),
        operation_index=OPERATION_INDEX,
    )

    assert isinstance(result, AuthenticatedNonEngagingOperation)
    assert isinstance(result.motion, native_type)
    assert result.operation_index == OPERATION_INDEX


def test_native_classified_plunge_retains_cut_plane_endpoint() -> None:
    result = classify_operation(
        _line(
            (1.0, 2.0, CLEARANCE_Z),
            (1.0, 2.0, CUT_Z),
            OperationType.PLUNGE,
        ),
        _cut_plane(),
        operation_index=OPERATION_INDEX,
    )

    assert isinstance(result, AuthenticatedPlungeOperation)
    assert result.endpoint == Point2[WorldXY].build(1.0, 2.0)
    assert result.operation_index == OPERATION_INDEX


def test_cut_height_line_is_supported_lateral_motion() -> None:
    result = classify_operation(
        _line((1.0, 2.0, CUT_Z), (3.0, 2.0, CUT_Z), OperationType.CUT),
        _cut_plane(),
        operation_index=OPERATION_INDEX,
    )

    assert isinstance(result, AuthenticatedLateralOperation)
    assert isinstance(result.motion, _stock_2.AuditSegmentMotion2)


def test_cut_height_arc_is_one_opaque_native_motion() -> None:
    geometry = Arc(
        radius=2.0,
        start_angle=0.0,
        end_angle=1.5 * math.pi,
        frame=Frame([1.0, 2.0, CUT_Z]),
    )
    operation = ToolpathOperation(
        geometry=geometry,
        operation=OperationType.CUT,
        path_index=4,
        clockwise=False,
    )

    result = classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)

    assert isinstance(result, AuthenticatedLateralOperation)
    assert isinstance(result.motion, _stock_2.AuditArcMotion2)


def test_cut_height_circle_is_supported_lateral_motion() -> None:
    operation = ToolpathOperation(
        geometry=Circle(2.0, frame=Frame([1.0, 2.0, CUT_Z])),
        operation=OperationType.CUT,
        path_index=5,
        clockwise=True,
    )

    result = classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)

    assert isinstance(result, AuthenticatedLateralOperation)
    assert isinstance(result.motion, _stock_2.AuditCircleMotion2)


def test_compas_arc_direction_must_map_unambiguously_to_native_orientation() -> None:
    operation = ToolpathOperation(
        geometry=Arc(
            radius=2.0,
            start_angle=0.0,
            end_angle=math.pi / 2.0,
            frame=Frame([1.0, 2.0, CUT_Z]),
        ),
        operation=OperationType.CUT,
        path_index=5,
        clockwise=True,
    )

    with pytest.raises(ContradictoryOperationOrientationError, match="exact arc sweep sign"):
        classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)


def test_retract_label_cannot_hide_cut_height_lateral_motion() -> None:
    operation = _line(
        (1.0, 2.0, CUT_Z),
        (3.0, 2.0, CUT_Z),
        OperationType.RETRACT,
    )

    with pytest.raises(ContradictoryOperationRoleError, match="exact lateral motion"):
        classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)


def test_clearance_transport_requires_both_endpoints_on_clearance_plane() -> None:
    operation = _line(
        (1.0, 2.0, CLEARANCE_Z),
        (3.0, 2.0, 4.0),
        OperationType.LINK,
    )

    with pytest.raises(UnsupportedAuditGeometryError, match="ramp"):
        classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)


def test_off_plane_lateral_motion_is_not_non_engaging() -> None:
    operation = _line((1.0, 2.0, 2.0), (3.0, 2.0, 2.0), OperationType.LINK)

    with pytest.raises(UnsupportedAuditGeometryError, match="declared plane"):
        classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)


def test_tilted_circle_is_unsupported_three_dimensional_geometry() -> None:
    operation = ToolpathOperation(
        geometry=Circle(
            2.0,
            frame=Frame(
                [1.0, 2.0, CUT_Z],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ),
        ),
        operation=OperationType.CUT,
        path_index=6,
    )

    with pytest.raises(UnsupportedAuditGeometryError, match="world"):
        classify_operation(operation, _cut_plane(), operation_index=OPERATION_INDEX)
