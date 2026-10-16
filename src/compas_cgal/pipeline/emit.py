"""Typed geometry to `ToolpathOperation`s, and the visualisation polyline.

One responsibility: turning decided geometry into operations. Nothing here
decides anything.

* `loop_operation` emits an `IncidentLoop` as a full `Circle` at the cut plane
  and carries the loop's EXACT entry in `ToolpathOperation.entry_xy` (Ruling
  1.1-C14-A): ``Circle.point_at(0)`` is only a double rounding of it.
* `arc_operation` follows the 1.6 convention
  ``Arc(rho, 0, ccw_angle(s -> e) mod 2pi, frame x = s - c)`` with the travel
  direction in ``clockwise``. The frame encodes the start; the angle places the
  end, so ``Arc.point_at(1)`` lands a few ulps off the exact ``end``, and the
  replay, the survey and the generator's piece admissibility all read that
  angle-placed end. An exact arc end is Task 1.6 work.
"""

from __future__ import annotations

import math

import numpy as np
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.errors import ZeroLengthLinkError
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import CutZ
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Visualisation density of circular motions, samples per radian (the frozen
# generators' polyline density, so figures compare like for like).
POLYLINE_SAMPLES_PER_RADIAN = 10.0

# Fewest samples on any circular motion of the polyline: both ends.
MIN_POLYLINE_SAMPLES = 2


def _travel_tangent(dx: float, dy: float, clockwise: bool) -> tuple[float, float]:
    """Unit travel direction at offset ``(dx, dy)`` from a centre, given the turn sense."""
    length = math.hypot(dx, dy)
    ux, uy = dx / length, dy / length
    return (uy, -ux) if clockwise else (-uy, ux)


def _tangent_array(direction: tuple[float, float]) -> np.ndarray:
    return np.array([direction[0], direction[1], 0.0], dtype=np.float64)


def loop_operation(loop: IncidentLoop, cut: CutZ, path_index: int, clockwise: bool) -> ToolpathOperation:
    """One full machining loop at the cut plane, starting at its exact entry.

    Args:
        loop: The authored loop.
        cut: Cut-plane height.
        path_index: Path index stamped on the operation.
        clockwise: Turn direction.

    Returns:
        A `CUT` circle whose frame x-axis points along the realised phase and
        whose `entry_xy` is the loop's exact entry.
    """
    ux, uy = math.cos(loop.phase), math.sin(loop.phase)
    frame = Frame([loop.centre.x, loop.centre.y, cut.value], [ux, uy, 0.0], [-uy, ux, 0.0])
    tangent = _tangent_array(_travel_tangent(loop.entry.x - loop.centre.x, loop.entry.y - loop.centre.y, clockwise))
    return ToolpathOperation(
        geometry=Circle(loop.radius.value, frame=frame),
        operation=OperationType.CUT,
        path_index=path_index,
        clockwise=clockwise,
        start_tangent=tangent,
        end_tangent=tangent,
        entry_xy=(loop.entry.x, loop.entry.y),
    )


def arc_operation(
    centre: Point2[WorldXY],
    start: Point2[WorldXY],
    end: Point2[WorldXY],
    clockwise: bool,
    cut: CutZ,
    path_index: int,
    operation: OperationType,
) -> ToolpathOperation:
    """A circular-arc motion from *start* to *end* about *centre*, travelled in *clockwise*.

    The radius is ``|start - centre|`` (the arc's endpoints are what the replay
    reads); the angle is the counter-clockwise angle from *start* to *end*,
    reduced into ``[0, 2pi)``, and reports the geometry only.

    Args:
        centre: Arc centre.
        start: Where the travel begins.
        end: Where it ends.
        clockwise: Travel direction.
        cut: Cut-plane height.
        path_index: Path index stamped on the operation.
        operation: The operation's label.

    Returns:
        The arc operation.

    Raises:
        ZeroLengthLinkError: *start* equals *end* exactly (no motion to emit).
    """
    if start == end:
        raise ZeroLengthLinkError(f"an arc from {start!r} to itself is no motion.")
    sx, sy = start.x - centre.x, start.y - centre.y
    radius = math.hypot(sx, sy)
    angle = (math.atan2(end.y - centre.y, end.x - centre.x) - math.atan2(sy, sx)) % math.tau
    ux, uy = sx / radius, sy / radius
    frame = Frame([centre.x, centre.y, cut.value], [ux, uy, 0.0], [-uy, ux, 0.0])
    return ToolpathOperation(
        geometry=Arc(radius, 0.0, angle, frame=frame),
        operation=operation,
        path_index=path_index,
        clockwise=clockwise,
        start_tangent=_tangent_array(_travel_tangent(sx, sy, clockwise)),
        end_tangent=_tangent_array(_travel_tangent(end.x - centre.x, end.y - centre.y, clockwise)),
    )


def line_operation(start: Point2[WorldXY], end: Point2[WorldXY], z_start: float, z_end: float, operation: OperationType, path_index: int) -> ToolpathOperation:
    """A straight motion, tangents filled for its XY travel (zero for a vertical move).

    Args:
        start: XY start.
        end: XY end.
        z_start: Start height.
        z_end: End height.
        operation: The operation's label.
        path_index: Path index stamped on the operation.

    Returns:
        The line operation.
    """
    dx, dy = end.x - start.x, end.y - start.y
    length = math.hypot(dx, dy)
    direction = (0.0, 0.0) if length == 0.0 else (dx / length, dy / length)
    tangent = _tangent_array(direction)
    return ToolpathOperation(
        geometry=Line([start.x, start.y, z_start], [end.x, end.y, z_end]),
        operation=operation,
        path_index=path_index,
        clockwise=False,
        start_tangent=tangent,
        end_tangent=tangent.copy(),
    )


def tessellate(operations: list[ToolpathOperation], samples_per_radian: float = POLYLINE_SAMPLES_PER_RADIAN) -> np.ndarray:
    """Sample the operation stream into one Nx3 visualisation polyline.

    Circles and arcs are sampled along their travel direction from their start;
    lines contribute their endpoints. Consecutive duplicates are dropped.

    Args:
        operations: The operations, in order.
        samples_per_radian: Angular density of circular motions.

    Returns:
        An ``(N, 3)`` float64 array.
    """
    points: list[tuple[float, float, float]] = []

    def push(x: float, y: float, z: float) -> None:
        candidate = (float(x), float(y), float(z))
        if not points or points[-1] != candidate:
            points.append(candidate)

    for op in operations:
        geometry = op.geometry
        if isinstance(geometry, Line):
            push(*geometry.start)
            push(*geometry.end)
            continue
        centre = geometry.frame.point
        if isinstance(geometry, Circle):
            span = math.tau
            sx, sy = op.entry_xy if op.entry_xy is not None else (float(geometry.point_at(0.0)[0]), float(geometry.point_at(0.0)[1]))
        else:
            start, end = geometry.point_at(0.0), geometry.point_at(1.0)
            sx, sy = float(start[0]), float(start[1])
            ccw = (math.atan2(float(end[1]) - centre[1], float(end[0]) - centre[0]) - math.atan2(sy - centre[1], sx - centre[0])) % math.tau
            span = (math.tau - ccw) % math.tau if op.clockwise else ccw
        count = max(MIN_POLYLINE_SAMPLES, math.ceil(span * samples_per_radian))
        turn = -span if op.clockwise else span
        for index in range(count + 1):
            angle = turn * index / count
            cosine, sine = math.cos(angle), math.sin(angle)
            push(centre[0] + cosine * (sx - centre[0]) - sine * (sy - centre[1]), centre[1] + sine * (sx - centre[0]) + cosine * (sy - centre[1]), centre[2])
    if not points:
        return np.empty((0, 3), dtype=np.float64)
    return np.asarray(points, dtype=np.float64, order="C")
