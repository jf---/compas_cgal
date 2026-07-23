"""Typed `Stock` wrapper over the exact-kernel CGAL boolean backend (`_stock_2`).

`Stock` is the ergonomic Python face of `compas_cgal._stock_2.Stock2`: it accepts
compas `Polygon` boundary and holes, validates planarity and orientation once at
the seam, and forwards the material-removal and membership queries 1:1 to the
exact backend. Per the exact-kernel boundary doctrine, doubles cross into
exact-land only here, at construction, by exact injection -- there is no snapping
or tolerance downstream of the vertex arrays handed to `Stock2`.
"""

import numpy as np
from compas.geometry import Polygon
from compas.tolerance import TOL

from compas_cgal import _stock_2  # type: ignore
from compas_cgal.toolpath import InvalidPolygonError


def _polygon_to_ccw_vertices(polygon: Polygon) -> np.ndarray:
    """Convert a planar XY `Polygon` to a C-ordered counterclockwise vertex array.

    Reimplements the conversion `compas_cgal.toolpath` performs for its own
    inputs (local copy by design -- the toolpath helper is private): rejects
    non-planar or degenerate input and canonicalizes winding to counterclockwise
    via the shoelace signed area, so `Stock2` always receives the CCW outer
    boundary it expects.

    Args:
        polygon: A `Polygon` whose vertices lie in the world XY plane.

    Returns:
        An ``(N, 3)`` C-ordered ``float64`` array of CCW vertices with any
        duplicated closing vertex removed.

    Raises:
        InvalidPolygonError: If the polygon has fewer than three points, is not
            ``Nx2``/``Nx3``, does not lie in the XY plane, or has vanishing
            signed area (degenerate or self-intersecting).
    """
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

    # Shoelace signed area: orientation without a polygon-normal heuristic, and a
    # direct diagnosis of degenerate / self-cancelling input (vanishing area).
    x, y = pts[:, 0], pts[:, 1]
    signed_area = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
    scale = float(np.abs(pts[:, :2]).max()) or 1.0
    if abs(signed_area) <= TOL.absolute * scale * scale:
        raise InvalidPolygonError("Polygon is degenerate or self-intersecting (vanishing signed area).")
    if signed_area < 0.0:
        pts = pts[::-1]

    return np.asarray(pts, dtype=np.float64, order="C")


class Stock:
    """Mutable 2D stock region backed by the exact CGAL boolean engine.

    Wraps `compas_cgal._stock_2.Stock2`. The stock is the still-uncut material in
    the XY plane; `subtract_*` methods remove the swept area of a tool motion and
    `contains` reports whether a point is still material. All queries and removals
    are exact (no tolerance): the wrapper only converts compas `Polygon` inputs to
    the CCW vertex arrays the backend expects and injects scalars as doubles.
    """

    def __init__(self, polygon: Polygon, holes: list[Polygon] | None = None) -> None:
        """Build a stock region from an outer boundary and optional island holes.

        Args:
            polygon: Outer boundary `Polygon` in the world XY plane.
            holes: Optional island polygons strictly inside *polygon*; each is
                converted to CCW vertices independently.

        Raises:
            InvalidPolygonError: If the boundary or any hole is non-planar,
                degenerate, or self-intersecting.
        """
        boundary = _polygon_to_ccw_vertices(polygon)
        hole_arrays = [_polygon_to_ccw_vertices(hole) for hole in holes] if holes else []
        self._raw = _stock_2.Stock2(boundary, hole_arrays)

    @property
    def raw(self) -> "_stock_2.Stock2":
        """The underlying exact-kernel `_stock_2.Stock2` backend instance."""
        return self._raw

    def contains(self, x: float, y: float) -> bool:
        """Return whether the point ``(x, y)`` is still material (inside the stock).

        Args:
            x: X coordinate of the query point.
            y: Y coordinate of the query point.

        Returns:
            ``True`` if the point lies in the current stock region.
        """
        return self._raw.contains(float(x), float(y))

    def is_empty(self) -> bool:
        """Return whether the stock region has been fully depleted.

        Returns:
            ``True`` if no material remains.
        """
        return self._raw.is_empty()

    def subtract_capsule(self, x0: float, y0: float, x1: float, y1: float, radius: float) -> None:
        """Remove the capsule swept by a tool of *radius* moving from start to end.

        Args:
            x0: X coordinate of the segment start (tool center).
            y0: Y coordinate of the segment start (tool center).
            x1: X coordinate of the segment end (tool center).
            y1: Y coordinate of the segment end (tool center).
            radius: Tool radius (capsule half-width).
        """
        self._raw.subtract_capsule(float(x0), float(y0), float(x1), float(y1), float(radius))

    def subtract_disk(self, cx: float, cy: float, radius: float) -> None:
        """Remove the disk swept by a tool of *radius* plunging at ``(cx, cy)``.

        Args:
            cx: X coordinate of the disk center (tool center).
            cy: Y coordinate of the disk center (tool center).
            radius: Tool radius (disk radius).
        """
        self._raw.subtract_disk(float(cx), float(cy), float(radius))

    def subtract_arc_sweep(
        self,
        cx: float,
        cy: float,
        sx: float,
        sy: float,
        ex: float,
        ey: float,
        cw: bool,
        tool_radius: float,
    ) -> None:
        """Remove the area swept by a tool of *tool_radius* along a circular arc.

        Args:
            cx: X coordinate of the arc center.
            cy: Y coordinate of the arc center.
            sx: X coordinate of the arc start (tool center).
            sy: Y coordinate of the arc start (tool center).
            ex: X coordinate of the arc end (tool center).
            ey: Y coordinate of the arc end (tool center).
            cw: ``True`` if the arc runs clockwise from start to end.
            tool_radius: Tool radius (half-width of the swept annulus).
        """
        self._raw.subtract_arc_sweep(
            float(cx),
            float(cy),
            float(sx),
            float(sy),
            float(ex),
            float(ey),
            bool(cw),
            float(tool_radius),
        )
