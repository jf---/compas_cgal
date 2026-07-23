"""Trochoidal toolpath generation along the polygon straight skeleton (CGAL backend).

Guide curves are interior straight-skeleton edges; every trochoid radius
derives from an exact clearance query at its own station center and every
bridge, lead, and flat link is certified against the boundary at tool radius,
so emitted motions are gouge-free by construction. The straight skeleton
coincides with the medial axis only for convex polygons; for non-convex
polygons the guide loci differ near reflex vertices while radii remain exact.

Deliberate scope choices (documented, not accidental): path ordering is a
greedy nearest-neighbour heuristic (optimal ordering is a TSP); junctions are
tangent-continuous (G1), not curvature-continuous; a single cutting depth per
call (compose calls for multi-depth roughing); flat links assume the corridor
between paths has been cleared — use ``clearance_z`` for stock-safe linking.
Clearance queries are exact and O(edges) per query without an acceleration
structure, which is ample for pocket boundaries of hundreds of edges.
"""

import math
import warnings
from dataclasses import dataclass
from enum import Enum

import numpy as np
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon
from compas.tolerance import TOL

from compas_cgal import _toolpath  # type: ignore

# Side-effect import: registers the bound std::vector types (nb::bind_vector)
# that _toolpath signatures resolve at call time.  _toolpath must not declare
# eager vector-typed default arguments, or its import would require these
# registrations to exist first (std::bad_cast at module init).
from compas_cgal import _types_std  # noqa: F401  # type: ignore
from compas_cgal.types import PolylinesNumpy

__all__ = [
    "DegeneratePrimitiveError",
    "InvalidPolygonError",
    "OperationType",
    "ToolpathOperation",
    "ToolpathResult",
    "polygon_skeleton_clearance",
    "trochoidal_mat_toolpath",
    "trochoidal_mat_toolpath_circular",
]

# Default certification margin as a fraction of the tool diameter: keeps the
# safety clearance proportional to the tool (scale-free) instead of a fixed
# absolute value that silently changes meaning between unit systems.
RADIAL_CLEARANCE_FRACTION = 1e-3


class InvalidPolygonError(ValueError):
    """The input polygon violates a precondition (planarity, simplicity, orientation)."""


class DegeneratePrimitiveError(ValueError):
    """The backend emitted a primitive whose geometry cannot be reconstructed."""


class OperationType(str, Enum):
    """Toolpath operation semantics; str-valued for G-code post compatibility."""

    CUT = "cut"
    LEAD_IN = "lead_in"
    LEAD_OUT = "lead_out"
    LINK = "link"
    RETRACT = "retract"
    PLUNGE = "plunge"


_OPERATION_BY_CODE = (
    OperationType.CUT,
    OperationType.LEAD_IN,
    OperationType.LEAD_OUT,
    OperationType.LINK,
    OperationType.RETRACT,
    OperationType.PLUNGE,
)


@dataclass
class ToolpathOperation:
    """A single toolpath motion primitive with operation metadata."""

    geometry: Line | Arc | Circle
    operation: OperationType
    path_index: int
    clockwise: bool = False
    start_tangent: np.ndarray | None = None
    end_tangent: np.ndarray | None = None


@dataclass
class ToolpathResult:
    """Toolpath output: typed operations for G-code + tessellated polyline for visualization."""

    operations: list[ToolpathOperation]
    polyline: np.ndarray


def _polygon_to_ccw_vertices(polygon: Polygon) -> np.ndarray:
    pts = np.asarray(polygon.points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[0] < 3:
        raise InvalidPolygonError("Expected a polygon with at least three points.")

    if pts.shape[1] == 2:
        pts = np.column_stack((pts, np.zeros(len(pts), dtype=np.float64)))
    elif pts.shape[1] != 3:
        raise InvalidPolygonError("Expected points with shape Nx2 or Nx3.")

    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]

    max_abs_z = float(np.abs(pts[:, 2]).max())
    if max_abs_z > TOL.absolute:
        raise InvalidPolygonError(f"Polygon should lie in the world XY plane (max |z| = {max_abs_z:.3e}).")

    # Shoelace signed area: robust orientation without a polygon-normal
    # heuristic, and a direct diagnosis for degenerate/self-cancelling input.
    x, y = pts[:, 0], pts[:, 1]
    signed_area = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
    scale = float(np.abs(pts[:, :2]).max()) or 1.0
    if abs(signed_area) <= TOL.absolute * scale * scale:
        raise InvalidPolygonError("Polygon is degenerate or self-intersecting (vanishing signed area).")
    if signed_area < 0.0:
        pts = pts[::-1]

    return np.asarray(pts, dtype=np.float64, order="C")


def _holes_to_vertex_arrays(holes: list[Polygon] | None) -> list[np.ndarray]:
    if not holes:
        return []
    return [_polygon_to_ccw_vertices(hole) for hole in holes]


def _warn_if_truncated(skipped: int, max_passes: int) -> None:
    if skipped > 0:
        warnings.warn(
            f"max_passes={max_passes} reached: {skipped} skeleton edge chain(s) left unmachined; " "raise max_passes for full coverage.",
            UserWarning,
            stacklevel=3,
        )


def _unit_xy_3(v: np.ndarray) -> np.ndarray | None:
    """Unit vector in XY from *v*, with z=0. None when degenerate."""
    norm = math.hypot(float(v[0]), float(v[1]))
    if norm <= TOL.absolute:
        return None
    inv = 1.0 / norm
    return np.array([float(v[0]) * inv, float(v[1]) * inv, 0.0], dtype=np.float64)


def _ccw_angle_xy(u: np.ndarray, v: np.ndarray) -> float:
    """CCW angle from *u* to *v* in [0, 2π)."""
    angle = math.atan2(
        float(u[0]) * float(v[1]) - float(u[1]) * float(v[0]),
        float(u[0]) * float(v[0]) + float(u[1]) * float(v[1]),
    )
    return angle + 2.0 * math.pi if angle < 0.0 else angle


def _row_to_line(start: np.ndarray, end: np.ndarray) -> Line:
    return Line(start.tolist(), end.tolist())


def _row_to_arc(start: np.ndarray, end: np.ndarray, center: np.ndarray, radius: float, clockwise: bool) -> Arc:
    start_vector = _unit_xy_3(start - center)
    end_vector = _unit_xy_3(end - center)
    if start_vector is None or end_vector is None:
        raise DegeneratePrimitiveError(f"Arc with coincident start/end and center cannot be reconstructed (center={center.tolist()}).")

    if clockwise:
        xaxis = end_vector
        yaxis = np.array([-xaxis[1], xaxis[0], 0.0], dtype=np.float64)
        start_angle = _ccw_angle_xy(xaxis, start_vector)
        end_angle = 0.0
    else:
        xaxis = start_vector
        yaxis = np.array([-xaxis[1], xaxis[0], 0.0], dtype=np.float64)
        start_angle = 0.0
        end_angle = _ccw_angle_xy(xaxis, end_vector)

    frame = Frame(center.tolist(), xaxis.tolist(), yaxis.tolist())
    return Arc(radius, start_angle, end_angle, frame=frame)


def _row_to_circle(center: np.ndarray, start: np.ndarray, radius: float) -> Circle:
    xaxis = _unit_xy_3(start - center)
    if xaxis is None:
        raise DegeneratePrimitiveError(f"Circle with start point on its center cannot be reconstructed (center={center.tolist()}).")
    yaxis = np.array([-xaxis[1], xaxis[0], 0.0], dtype=np.float64)
    frame = Frame(center.tolist(), xaxis.tolist(), yaxis.tolist())
    return Circle(radius, frame=frame)


def _matrices_to_operations(
    meta: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    centers: np.ndarray,
    radii: np.ndarray,
    start_tangents: np.ndarray | None = None,
    end_tangents: np.ndarray | None = None,
) -> list[ToolpathOperation]:
    ops: list[ToolpathOperation] = []
    for i in range(meta.shape[0]):
        path_index = int(meta[i, 0])
        ptype = int(meta[i, 1])
        clockwise = bool(meta[i, 2] > 0.5)
        operation = _OPERATION_BY_CODE[int(meta[i, 3])]

        if ptype == 0:
            geom: Line | Arc | Circle = _row_to_line(starts[i], ends[i])
        elif ptype == 1:
            geom = _row_to_arc(starts[i], ends[i], centers[i], float(radii[i]), clockwise)
        else:
            geom = _row_to_circle(centers[i], starts[i], float(radii[i]))

        st = start_tangents[i] if start_tangents is not None else None
        et = end_tangents[i] if end_tangents is not None else None
        ops.append(
            ToolpathOperation(
                geometry=geom,
                operation=operation,
                path_index=path_index,
                clockwise=clockwise,
                start_tangent=st,
                end_tangent=et,
            )
        )
    return ops


def polygon_skeleton_clearance(polygon: Polygon, holes: list[Polygon] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Interior straight-skeleton vertices with exact boundary clearance.

    The straight skeleton coincides with the medial axis for convex polygons
    only; for non-convex polygons the loci differ near reflex vertices, but
    each returned radius is the exact clearance (largest gouge-free disk) at
    its locus.

    Parameters
    ----------
    polygon
        A `Polygon` lying in the XY plane.
    holes
        Optional hole polygons strictly inside *polygon*.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Interior skeleton points (Mx3) and their exact clearance radii (M,).

    Raises
    ------
    InvalidPolygonError
        If any polygon is non-planar, degenerate, or self-intersecting.
    """
    V = _polygon_to_ccw_vertices(polygon)
    H = _holes_to_vertex_arrays(holes)
    points, radii = _toolpath.polygon_skeleton_clearance(V, H)
    points = np.asarray(points, dtype=np.float64, order="C")
    radii = np.asarray(radii, dtype=np.float64, order="C").reshape(-1)
    return points, radii


def trochoidal_mat_toolpath(
    polygon: Polygon,
    tool_diameter: float,
    stepover: float | None = None,
    pitch: float | None = None,
    min_trochoid_radius: float = 0.0,
    max_trochoid_radius: float | None = None,
    mat_scale: float = 1.0,
    radial_clearance: float | None = None,
    samples_per_cycle: int = 20,
    max_passes: int = 1000,
    climb: bool = True,
    holes: list[Polygon] | None = None,
) -> PolylinesNumpy:
    """Certified gouge-free trochoidal pocketing paths along the straight skeleton.

    Every trochoid radius derives from an exact clearance query at its own
    station center and every bridge move is certified against the boundary at
    tool radius, so emitted tool-center polylines never approach a wall closer
    than the tool radius.

    Parameters
    ----------
    polygon
        A `Polygon` lying in the XY plane.
    tool_diameter
        Tool diameter.
    stepover
        Maximum radial engagement per trochoid cycle (station advance plus
        radius growth). Defaults to ``0.4 * tool_diameter``.
    pitch
        Desired trochoid advance per cycle; the effective advance is capped by
        *stepover*. Defaults to ``0.75 * tool_diameter``.
    min_trochoid_radius
        Minimum desired trochoid radius.
    max_trochoid_radius
        Maximum trochoid radius. Defaults to no cap.
    mat_scale
        Scale factor applied to the clearance-derived available radius.
    radial_clearance
        Safety clearance subtracted from the available radius. Defaults to
        ``RADIAL_CLEARANCE_FRACTION * tool_diameter`` (scale-free).
    samples_per_cycle
        Polyline samples per arc primitive.
    max_passes
        Maximum number of emitted toolpaths; exceeding it warns.
    climb
        ``True`` for climb milling (clockwise arcs), ``False`` for
        conventional milling (counterclockwise arcs).
    holes
        Optional island polygons strictly inside *polygon*; paths keep tool
        clearance from island boundaries as well.

    Returns
    -------
    `compas_cgal.types.PolylinesNumpy`
        List of Nx3 toolpath polylines.

    Raises
    ------
    InvalidPolygonError
        If any polygon is non-planar, degenerate, or self-intersecting.

    Warns
    -----
    UserWarning
        When *max_passes* truncates coverage (never silently).
    """
    if stepover is None:
        stepover = 0.4 * tool_diameter
    if pitch is None:
        pitch = 0.75 * tool_diameter
    if max_trochoid_radius is None:
        max_trochoid_radius = float("inf")
    if radial_clearance is None:
        radial_clearance = RADIAL_CLEARANCE_FRACTION * tool_diameter

    V = _polygon_to_ccw_vertices(polygon)
    H = _holes_to_vertex_arrays(holes)
    paths, skipped = _toolpath.trochoidal_mat_toolpath(
        V,
        H,
        float(tool_diameter),
        float(stepover),
        float(pitch),
        float(min_trochoid_radius),
        float(max_trochoid_radius),
        float(mat_scale),
        float(radial_clearance),
        int(samples_per_cycle),
        int(max_passes),
        bool(climb),
    )
    _warn_if_truncated(int(skipped), max_passes)
    return [np.asarray(path, dtype=np.float64, order="C") for path in paths]


def trochoidal_mat_toolpath_circular(
    polygon: Polygon,
    tool_diameter: float,
    stepover: float | None = None,
    pitch: float | None = None,
    min_trochoid_radius: float = 0.0,
    max_trochoid_radius: float | None = None,
    mat_scale: float = 1.0,
    radial_clearance: float | None = None,
    samples_per_cycle: int = 20,
    max_passes: int = 1000,
    climb: bool = True,
    holes: list[Polygon] | None = None,
    lead_in: float = 0.0,
    lead_out: float = 0.0,
    link_paths: bool = True,
    optimize_order: bool = True,
    cut_z: float = 0.0,
    clearance_z: float | None = None,
    retract_at_end: bool = True,
    samples_per_radian: float = 10.0,
) -> ToolpathResult:
    """Generate linked circular arc toolpath operations.

    Leads are certified against the boundary and shrink (or drop) to stay
    gouge-free. Flat links (no *clearance_z*) that would cross a wall raise
    instead of being emitted; they also assume the corridor between paths has
    already been cleared of stock — provide *clearance_z* for stock-safe
    linking. Path ordering is a greedy nearest-neighbour heuristic.

    Parameters
    ----------
    polygon
        A `Polygon` lying in the XY plane.
    tool_diameter
        Tool diameter.
    stepover
        Maximum radial engagement per trochoid cycle. Defaults to
        ``0.4 * tool_diameter``.
    pitch
        Desired trochoid advance per cycle; capped by *stepover*. Defaults to
        ``0.75 * tool_diameter``.
    min_trochoid_radius
        Minimum desired trochoid radius.
    max_trochoid_radius
        Maximum trochoid radius. Defaults to no cap.
    mat_scale
        Scale factor applied to the clearance-derived available radius.
    radial_clearance
        Safety clearance subtracted from the available radius. Defaults to
        ``RADIAL_CLEARANCE_FRACTION * tool_diameter`` (scale-free).
    samples_per_cycle
        Polyline samples per arc primitive.
    max_passes
        Maximum number of emitted toolpaths; exceeding it warns.
    climb
        ``True`` for climb milling (clockwise arcs), ``False`` for
        conventional milling (counterclockwise arcs).
    holes
        Optional island polygons strictly inside *polygon*.
    lead_in
        Entry line length along the start tangent for each path.
    lead_out
        Exit line length along the end tangent for each path.
    link_paths
        If ``True``, connect consecutive paths with linear link moves.
    optimize_order
        If ``True``, greedily reorder paths by nearest entry.
    cut_z
        Z-height of cutting motions.
    clearance_z
        Safe Z-height for inter-path travel.
    retract_at_end
        If ``True``, add a final retract after the last path.
    samples_per_radian
        Tessellation density for the visualization polyline.

    Returns
    -------
    ToolpathResult
        Operations stream with typed compas geometry and tessellated polyline.

    Raises
    ------
    InvalidPolygonError
        If any polygon is non-planar, degenerate, or self-intersecting.
    ValueError
        If a flat link would gouge the boundary (provide *clearance_z*).

    Warns
    -----
    UserWarning
        When *max_passes* truncates coverage (never silently).
    """
    if stepover is None:
        stepover = 0.4 * tool_diameter
    if pitch is None:
        pitch = 0.75 * tool_diameter
    if max_trochoid_radius is None:
        max_trochoid_radius = float("inf")
    if radial_clearance is None:
        radial_clearance = RADIAL_CLEARANCE_FRACTION * tool_diameter

    has_clearance_z = clearance_z is not None
    effective_clearance_z = float(clearance_z) if clearance_z is not None else 0.0

    V = _polygon_to_ccw_vertices(polygon)
    H = _holes_to_vertex_arrays(holes)
    meta, starts, ends, centers, radii, polyline, start_tangents, end_tangents, skipped = _toolpath.trochoidal_mat_toolpath_circular(
        V,
        H,
        float(tool_diameter),
        float(stepover),
        float(pitch),
        float(min_trochoid_radius),
        float(max_trochoid_radius),
        float(mat_scale),
        float(radial_clearance),
        int(samples_per_cycle),
        int(max_passes),
        bool(climb),
        float(lead_in),
        float(lead_out),
        bool(link_paths),
        bool(optimize_order),
        float(cut_z),
        effective_clearance_z,
        has_clearance_z,
        bool(retract_at_end),
        float(samples_per_radian),
    )
    _warn_if_truncated(int(skipped), max_passes)

    meta = np.asarray(meta, dtype=np.float64)
    if meta.size == 0:
        return ToolpathResult(operations=[], polyline=np.empty((0, 3), dtype=np.float64))

    starts = np.asarray(starts, dtype=np.float64, order="C")
    ends = np.asarray(ends, dtype=np.float64, order="C")
    centers = np.asarray(centers, dtype=np.float64, order="C")
    radii = np.asarray(radii, dtype=np.float64, order="C").reshape(-1)
    polyline = np.asarray(polyline, dtype=np.float64, order="C")
    start_tangents = np.asarray(start_tangents, dtype=np.float64, order="C")
    end_tangents = np.asarray(end_tangents, dtype=np.float64, order="C")

    return ToolpathResult(
        operations=_matrices_to_operations(meta, starts, ends, centers, radii, start_tangents, end_tangents),
        polyline=polyline,
    )
