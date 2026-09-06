"""Exact design-pocket residual from emitted planar cutter operations."""

from __future__ import annotations

import numpy as np
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from benchmarks.errors import BenchmarkError
from benchmarks.held_exact_motion_coverage import IncompleteMotionCoverageError
from benchmarks.spec import PocketSpec
from compas_cgal import _coverage_2
from compas_cgal.replay_classification import classify_operation_replay
from compas_cgal.replay_classification import infer_operation_cut_height
from compas_cgal.toolpath import ToolpathResult


class UnsupportedCoverageMotionError(BenchmarkError):
    """An emitted motion has no exact planar coverage representation."""


def replay_toolpath_coverage(spec: PocketSpec, result: ToolpathResult) -> _coverage_2.ExactRegion2:
    """Return exact residual of the design pocket, initially entirely uncut.

    Uses the public replay classification. Cutting primitives must additionally
    be exactly planar: its reporting tolerance cannot authorize projected cuts.
    Full circles use their declared radius, independent of frame phase. Partial
    arcs are unsupported; the visualization polyline is never replayed.
    """
    target = _coverage_2.ExactRegion2.from_polygon(
        np.asarray(spec.polygon.points, dtype=np.float64),
        [np.asarray(hole.points, dtype=np.float64) for hole in spec.holes],
    )
    circles: list[tuple[float, float, float]] = []
    segments: list[tuple[float, float, float, float]] = []
    disks: list[tuple[float, float]] = []
    cut_z = infer_operation_cut_height(result.operations)
    for index, operation in enumerate(result.operations):
        geometry = operation.geometry
        category = classify_operation_replay(operation, cut_z)
        if category in {"rapid", "retract"}:
            continue
        if isinstance(geometry, Line):
            first, second = geometry.start, geometry.end
            if category == "plunge":
                if first.x != second.x or first.y != second.y or second.z != cut_z:
                    raise UnsupportedCoverageMotionError(f"Operation {index}: plunge must be vertical and end on the cut plane.")
                disks.append((second.x, second.y))
            else:
                if first.z != cut_z or second.z != cut_z:
                    raise UnsupportedCoverageMotionError(f"Operation {index}: cutting line must lie exactly on the cut plane.")
                if first.x == second.x and first.y == second.y:
                    disks.append((first.x, first.y))
                else:
                    segments.append((first.x, first.y, second.x, second.y))
        elif isinstance(geometry, Arc):
            raise UnsupportedCoverageMotionError(f"Operation {index}: exact partial arc cutter sweep is not implemented.")
        elif isinstance(geometry, Circle):
            frame = geometry.frame
            if frame.xaxis.z != 0.0 or frame.yaxis.z != 0.0 or frame.point.z != cut_z:
                raise UnsupportedCoverageMotionError(f"Operation {index}: circle must lie exactly parallel to world XY on the cut plane.")
            circles.append((frame.point.x, frame.point.y, geometry.radius))
        else:
            raise UnsupportedCoverageMotionError(f"Operation {index}: unsupported cutter geometry {type(geometry).__name__}.")
    return _coverage_2.remaining_material(
        target,
        np.array(circles, dtype=np.float64).reshape((-1, 3)),
        np.array(segments, dtype=np.float64).reshape((-1, 4)),
        np.array(disks, dtype=np.float64).reshape((-1, 2)),
        spec.tool_radius,
    )


def require_toolpath_coverage(spec: PocketSpec, result: ToolpathResult) -> None:
    """Reject any nonempty exact design residual before accepting path metrics."""
    residual = replay_toolpath_coverage(spec, result)
    if not residual.is_empty():
        raise IncompleteMotionCoverageError(f"{spec.name}: {residual.component_count()} exact residual components; path length is unqualified.")
