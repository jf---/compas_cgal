from __future__ import annotations

import math
from dataclasses import FrozenInstanceError
from typing import Callable
from typing import cast

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

from benchmarks.errors import InvalidHeldOperationSnapshotError
from benchmarks.errors import MutatedHeldToolpathError
from benchmarks.held_path_snapshot import HeldArcSnapshot
from benchmarks.held_path_snapshot import HeldCircleSnapshot
from benchmarks.held_path_snapshot import HeldLineSnapshot
from benchmarks.held_path_snapshot import assert_toolpath_matches_snapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.units import OperationIndex
from compas_cgal.adaptive.units import Direction3
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point3
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import WorldXYZ
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

ResultMutation = Callable[[ToolpathResult], None]


def _operation(
    geometry: Line | Arc | Circle,
    *,
    operation: OperationType = OperationType.CUT,
    path_index: int = 0,
    clockwise: bool = False,
    start_tangent: np.ndarray | None = None,
    end_tangent: np.ndarray | None = None,
) -> ToolpathOperation:
    return ToolpathOperation(
        geometry=geometry,
        operation=operation,
        path_index=path_index,
        clockwise=clockwise,
        start_tangent=start_tangent,
        end_tangent=end_tangent,
    )


def _result(operations: list[ToolpathOperation]) -> ToolpathResult:
    return ToolpathResult(operations=operations, polyline=np.empty((0, 3), dtype=np.float64))


def _unchecked_direction(components: tuple[float, float, float]) -> Direction3[WorldXYZ]:
    direction = object.__new__(Direction3)
    object.__setattr__(direction, "x", components[0])
    object.__setattr__(direction, "y", components[1])
    object.__setattr__(direction, "z", components[2])
    return cast(Direction3[WorldXYZ], direction)


def _representative_result() -> ToolpathResult:
    return _result(
        [
            _operation(
                Line([0.0, 0.0, 1.0], [1.0, 0.0, 1.0]),
                path_index=2,
                start_tangent=np.array([1.0, 0.0, 0.0]),
                end_tangent=np.array([1.0, 0.0, 0.0]),
            ),
            _operation(
                Arc(
                    2.0,
                    0.25,
                    1.25,
                    frame=Frame([2.0, 3.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                ),
                operation=OperationType.LEAD_IN,
                path_index=3,
                clockwise=True,
                start_tangent=np.array([0.0, 1.0, 0.0]),
                end_tangent=np.array([-1.0, 0.0, 0.0]),
            ),
            _operation(
                Circle(
                    3.0,
                    frame=Frame([4.0, 5.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]),
                ),
                operation=OperationType.LINK,
                path_index=4,
                clockwise=True,
            ),
        ]
    )


def test_line_snapshot_retains_complete_operation_behavior() -> None:
    snapshot = snapshot_toolpath(_representative_result())[0]

    assert isinstance(snapshot, HeldLineSnapshot)
    assert snapshot.ordinal == OperationIndex(0)
    assert snapshot.operation is OperationType.CUT
    assert snapshot.path_index == 2
    assert snapshot.clockwise is False
    assert snapshot.start == Point3[WorldXYZ].build(0.0, 0.0, 1.0)
    assert snapshot.end == Point3[WorldXYZ].build(1.0, 0.0, 1.0)
    assert snapshot.start_tangent == Direction3[WorldXYZ].build(1.0, 0.0, 0.0)
    assert snapshot.end_tangent == Direction3[WorldXYZ].build(1.0, 0.0, 0.0)


def test_arc_snapshot_retains_complete_curve_behavior() -> None:
    snapshot = snapshot_toolpath(_representative_result())[1]

    assert isinstance(snapshot, HeldArcSnapshot)
    assert snapshot.ordinal == OperationIndex(1)
    assert snapshot.operation is OperationType.LEAD_IN
    assert snapshot.path_index == 3
    assert snapshot.clockwise is True
    assert snapshot.centre == Point3[WorldXYZ].build(2.0, 3.0, 1.0)
    assert snapshot.xaxis == Direction3[WorldXYZ].build(1.0, 0.0, 0.0)
    assert snapshot.yaxis == Direction3[WorldXYZ].build(0.0, 1.0, 0.0)
    assert snapshot.radius == Millimetre(2.0)
    assert snapshot.start_angle == Radian(0.25)
    assert snapshot.end_angle == Radian(1.25)
    assert snapshot.start_tangent == Direction3[WorldXYZ].build(0.0, 1.0, 0.0)
    assert snapshot.end_tangent == Direction3[WorldXYZ].build(-1.0, 0.0, 0.0)


def test_circle_snapshot_retains_frame_phase() -> None:
    circle = Circle(2.0, frame=Frame([3.0, 4.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]))
    snapshot = snapshot_toolpath(_result([_operation(circle, clockwise=True)]))[0]

    assert isinstance(snapshot, HeldCircleSnapshot)
    assert snapshot.ordinal == OperationIndex(0)
    assert snapshot.operation is OperationType.CUT
    assert snapshot.path_index == 0
    assert snapshot.clockwise is True
    assert snapshot.centre == Point3[WorldXYZ].build(3.0, 4.0, 0.0)
    assert snapshot.xaxis == Direction3[WorldXYZ].build(0.0, 1.0, 0.0)
    assert snapshot.yaxis == Direction3[WorldXYZ].build(-1.0, 0.0, 0.0)
    assert snapshot.radius == Millimetre(2.0)
    assert snapshot.start_tangent is None
    assert snapshot.end_tangent is None


def test_snapshot_retains_no_mutable_ingress_objects() -> None:
    result = _representative_result()
    source_line = cast(Line, result.operations[0].geometry)
    source_tangent = cast(np.ndarray, result.operations[0].start_tangent)
    snapshots = snapshot_toolpath(result)

    source_line.start.x = 99.0
    source_tangent[0] = -1.0

    line = cast(HeldLineSnapshot, snapshots[0])
    assert line.start == Point3[WorldXYZ].build(0.0, 0.0, 1.0)
    assert line.start_tangent == Direction3[WorldXYZ].build(1.0, 0.0, 0.0)
    assert all(not isinstance(value, (Line, Arc, Circle, np.ndarray)) for snapshot in snapshots for value in snapshot.__dict__.values())


def test_structural_comparison_accepts_unchanged_toolpath() -> None:
    result = _representative_result()

    assert_toolpath_matches_snapshot(result, snapshot_toolpath(result))


def _reverse_operations(result: ToolpathResult) -> None:
    result.operations.reverse()


def _remove_operation(result: ToolpathResult) -> None:
    result.operations.pop()


def _change_role(result: ToolpathResult) -> None:
    result.operations[0].operation = OperationType.RETRACT


def _change_path_index(result: ToolpathResult) -> None:
    result.operations[0].path_index = 20


def _change_clockwise(result: ToolpathResult) -> None:
    result.operations[0].clockwise = True


def _change_line_start(result: ToolpathResult) -> None:
    result.operations[0].geometry = Line([-1.0, 0.0, 1.0], [1.0, 0.0, 1.0])


def _change_line_end(result: ToolpathResult) -> None:
    result.operations[0].geometry = Line([0.0, 0.0, 1.0], [2.0, 0.0, 1.0])


def _change_line_z(result: ToolpathResult) -> None:
    result.operations[0].geometry = Line([0.0, 0.0, 2.0], [1.0, 0.0, 1.0])


def _change_arc_centre(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(2.0, 0.25, 1.25, frame=Frame([3.0, 3.0, 1.0]))


def _change_arc_centre_z(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(2.0, 0.25, 1.25, frame=Frame([2.0, 3.0, 2.0]))


def _change_arc_xaxis(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(
        2.0,
        0.25,
        1.25,
        frame=Frame([2.0, 3.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]),
    )


def _change_arc_yaxis(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(
        2.0,
        0.25,
        1.25,
        frame=Frame([2.0, 3.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    )


def _change_arc_radius(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(2.5, 0.25, 1.25, frame=Frame([2.0, 3.0, 1.0]))


def _change_arc_start_angle(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(2.0, 0.5, 1.25, frame=Frame([2.0, 3.0, 1.0]))


def _change_arc_end_angle(result: ToolpathResult) -> None:
    result.operations[1].geometry = Arc(2.0, 0.25, 1.5, frame=Frame([2.0, 3.0, 1.0]))


def _change_start_tangent(result: ToolpathResult) -> None:
    result.operations[0].start_tangent = np.array([0.0, 1.0, 0.0])


def _change_end_tangent(result: ToolpathResult) -> None:
    result.operations[0].end_tangent = None


@pytest.mark.parametrize(
    "mutate",
    [
        _reverse_operations,
        _remove_operation,
        _change_role,
        _change_path_index,
        _change_clockwise,
        _change_line_start,
        _change_line_end,
        _change_line_z,
        _change_arc_centre,
        _change_arc_centre_z,
        _change_arc_xaxis,
        _change_arc_yaxis,
        _change_arc_radius,
        _change_arc_start_angle,
        _change_arc_end_angle,
        _change_start_tangent,
        _change_end_tangent,
    ],
)
def test_structural_comparison_rejects_each_behavior_change(mutate: ResultMutation) -> None:
    baseline = _representative_result()
    expected = snapshot_toolpath(baseline)
    changed = _representative_result()
    mutate(changed)

    with pytest.raises(MutatedHeldToolpathError, match="differs from the characterized operation snapshot"):
        assert_toolpath_matches_snapshot(changed, expected)


@pytest.mark.parametrize("coordinate", [math.nan, math.inf, -math.inf])
def test_snapshot_rejects_non_finite_line_geometry(coordinate: float) -> None:
    result = _result([_operation(Line([0.0, 0.0, 0.0], [coordinate, 1.0, 0.0]))])

    with pytest.raises(InvalidHeldOperationSnapshotError, match="finite"):
        snapshot_toolpath(result)


@pytest.mark.parametrize("radius", [0.0, -1.0, math.nan, math.inf])
def test_circle_factory_rejects_invalid_radius(radius: float) -> None:
    with pytest.raises(InvalidHeldOperationSnapshotError, match="radius"):
        HeldCircleSnapshot.build(
            ordinal=OperationIndex(0),
            operation=OperationType.CUT,
            path_index=0,
            clockwise=False,
            centre=Point3[WorldXYZ].build(0.0, 0.0, 0.0),
            xaxis=Direction3[WorldXYZ].build(1.0, 0.0, 0.0),
            yaxis=Direction3[WorldXYZ].build(0.0, 1.0, 0.0),
            radius=Millimetre(radius),
            start_tangent=None,
            end_tangent=None,
        )


@pytest.mark.parametrize(
    ("xaxis", "yaxis"),
    [
        ((2.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        ((1.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
        ((math.nan, 0.0, 0.0), (0.0, 1.0, 0.0)),
    ],
)
def test_curve_factory_rejects_malformed_frame_axes(
    xaxis: tuple[float, float, float],
    yaxis: tuple[float, float, float],
) -> None:
    with pytest.raises(InvalidHeldOperationSnapshotError, match="axes"):
        HeldCircleSnapshot.build(
            ordinal=OperationIndex(0),
            operation=OperationType.CUT,
            path_index=0,
            clockwise=False,
            centre=Point3[WorldXYZ].build(0.0, 0.0, 0.0),
            xaxis=_unchecked_direction(xaxis),
            yaxis=_unchecked_direction(yaxis),
            radius=Millimetre(1.0),
            start_tangent=None,
            end_tangent=None,
        )


@pytest.mark.parametrize("angle", [math.nan, math.inf, -math.inf])
def test_arc_factory_rejects_invalid_angles(angle: float) -> None:
    with pytest.raises(InvalidHeldOperationSnapshotError, match="angle"):
        HeldArcSnapshot.build(
            ordinal=OperationIndex(0),
            operation=OperationType.CUT,
            path_index=0,
            clockwise=False,
            centre=Point3[WorldXYZ].build(0.0, 0.0, 0.0),
            xaxis=Direction3[WorldXYZ].build(1.0, 0.0, 0.0),
            yaxis=Direction3[WorldXYZ].build(0.0, 1.0, 0.0),
            radius=Millimetre(1.0),
            start_angle=Radian(angle),
            end_angle=Radian(1.0),
            start_tangent=None,
            end_tangent=None,
        )


@pytest.mark.parametrize(
    "tangent",
    [
        np.array([1.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
        np.array([math.nan, 0.0, 0.0]),
        [1.0, 0.0, 0.0],
    ],
)
def test_snapshot_rejects_malformed_tangent(tangent: object) -> None:
    result = _result(
        [
            _operation(
                Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
                start_tangent=cast(np.ndarray, tangent),
            )
        ]
    )

    with pytest.raises(InvalidHeldOperationSnapshotError, match="tangent"):
        snapshot_toolpath(result)


def test_snapshot_rejects_unsupported_geometry() -> None:
    operation = _operation(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]))
    operation.geometry = cast(Line | Arc | Circle, object())

    with pytest.raises(InvalidHeldOperationSnapshotError, match="supported exact primitive"):
        snapshot_toolpath(_result([operation]))


@pytest.mark.parametrize("record", [HeldLineSnapshot, HeldArcSnapshot, HeldCircleSnapshot])
def test_direct_snapshot_construction_is_disabled(record: type[object]) -> None:
    with pytest.raises(TypeError):
        record()


def test_snapshots_are_frozen() -> None:
    snapshot = snapshot_toolpath(_representative_result())[0]

    with pytest.raises(FrozenInstanceError):
        snapshot.clockwise = True  # type: ignore[misc]
