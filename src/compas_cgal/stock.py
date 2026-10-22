"""Typed `Stock` wrapper over the exact-kernel CGAL boolean backend (`_stock_2`).

`Stock` is the ergonomic Python face of `compas_cgal._stock_2.Stock2`: it accepts
compas `Polygon` boundary and holes, validates planarity and orientation once at
the seam, and forwards the material-removal and membership queries 1:1 to the
exact backend. Per the exact-kernel boundary doctrine, doubles cross into
exact-land only here, at construction, by exact injection -- there is no snapping
or tolerance downstream of the vertex arrays handed to `Stock2`.
"""

from __future__ import annotations

import enum

import numpy as np
from compas.geometry import Polygon
from compas.tolerance import TOL

from compas_cgal import _stock_2
from compas_cgal.toolpath import InvalidPolygonError


class SweepRemoval(enum.Enum):
    """What one motion removes, bracketed exactly between two certified regions.

    The replay subtracts an UNDER region (a subset of the true swept area); a
    certified OVER region contains the true swept area. Both tests are exact
    regularized interior intersections with the stock.

    Attributes:
        REMOVES: The replay's own region meets the stock's interior, so the
            replay changes the stock.
        WITHIN_BAND: The replay's region misses the stock but the over-region
            meets it: any material the tool still cuts lies inside the certified
            band, so neither other answer can be claimed.
        NOTHING: Even the over-region misses the stock: the tool removes nothing.
    """

    REMOVES = "removes"
    WITHIN_BAND = "within_band"
    NOTHING = "nothing"


def _sweep_removal(verdict: "_stock_2.SweepRemovalVerdict") -> SweepRemoval:
    """Map the native verdict onto `SweepRemoval` by name (one vocabulary)."""
    return SweepRemoval[verdict.name]


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

    @classmethod
    def _from_raw(cls, raw: "_stock_2.Stock2") -> "Stock":
        stock = cls.__new__(cls)
        stock._raw = raw
        return stock

    def clone(self) -> "Stock":
        """Return an exact deep copy with independent native ownership."""
        return self._from_raw(self._raw.clone())

    def is_subset_of(self, other: "Stock") -> bool:
        """Return exact regularized-set inclusion in *other*."""
        return self._raw.is_subset_of(other._raw)

    def exactly_equals(self, other: "Stock") -> bool:
        """Return exact regularized-set equality with *other*."""
        return self._raw.exactly_equals(other._raw)

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

    def subtract_capsule_quad(self, x0: float, y0: float, x1: float, y1: float, radius: float) -> None:
        """Remove the same capsule as `subtract_capsule`, as two disks and a rectangle.

        The capsule's true side lines stand off the segment by an irrational
        distance that the circle-segment traits cannot hold, so the region is
        under-covered -- exactly as the disk chain is. The difference is the
        SHAPE of the under-approximation: the two end disks (exact) unioned with
        the rectangle along the segment at a marginally reduced half-width, six
        curves instead of the chain's hundreds of arcs.

        Both halves of the contract are checked as exact rational comparisons in
        the backend, not estimated:

        - nothing outside the true swept capsule is ever removed, and
        - everything within ``radius * (1 - 1e-4)`` of the segment IS removed,
          which is the same slack budget the disk chain documents.

        Args:
            x0: X coordinate of the segment start (tool center).
            y0: Y coordinate of the segment start (tool center).
            x1: X coordinate of the segment end (tool center).
            y1: Y coordinate of the segment end (tool center).
            radius: Tool radius (capsule half-width).

        Raises:
            NonFiniteCapsuleInputError: If any endpoint coordinate or the radius
                is not finite.
            CapsuleQuadCertificateError: If the constructed rectangle's exact
                half-width falls outside its certified band -- a broken internal
                invariant, never a malformed input.
        """
        self._raw.subtract_capsule_quad(float(x0), float(y0), float(x1), float(y1), float(radius))

    def subtract_line_sweep(self, x0: float, y0: float, x1: float, y1: float, radius: float) -> None:
        """Remove the replay's model of the sweep along a segment: exact, or certified under.

        The model is decided EXACTLY on the injected endpoints: when
        ``dx**2 + dy**2`` is the square of a rational the capsule's side lines
        are rational and the region removed IS the swept capsule. Otherwise it
        is ``disk(A, r) | rect | disk(B, r)``: the end disks are exact (the
        tool's endpoint positions are), and the rectangle has half-width ``h``
        with ``(1 - 1e-4) r <= h <= r`` checked exactly, fixed by the
        direction class. A plunge or the next motion's start at an endpoint
        therefore meets nothing the model left; only a thin side band may
        survive.

        Args:
            x0: X coordinate of the segment start (tool center).
            y0: Y coordinate of the segment start (tool center).
            x1: X coordinate of the segment end (tool center).
            y1: Y coordinate of the segment end (tool center).
            radius: Tool radius.

        Raises:
            NonFiniteCapsuleInputError: A coordinate or the radius is not finite.
            InvalidLineSweepRadiusError: The radius is not positive.
            LineSweepBandCertificateError: The rectangle's exact half-width left
                its certified band -- a broken internal invariant.
        """
        self._raw.subtract_line_sweep(float(x0), float(y0), float(x1), float(y1), float(radius))

    def subtract_line_sweep_local(self, x0: float, y0: float, x1: float, y1: float, radius: float) -> None:
        """Remove the same region as `subtract_line_sweep`, by local arrangement surgery.

        Args:
            x0: X coordinate of the segment start (tool center).
            y0: Y coordinate of the segment start (tool center).
            x1: X coordinate of the segment end (tool center).
            y1: Y coordinate of the segment end (tool center).
            radius: Tool radius.

        Raises:
            NonFiniteCapsuleInputError: A coordinate or the radius is not finite.
            InvalidLineSweepRadiusError: The radius is not positive.
            LineSweepBandCertificateError: The rectangle's exact half-width left
                its certified band.
            LocalDepletionEscapedError: The flood escaped the region -- a broken
                internal invariant.
        """
        self._raw.subtract_line_sweep_local(float(x0), float(y0), float(x1), float(y1), float(radius))

    def subtract_circle_sweep(self, cx: float, cy: float, guide_radius: float, tool_radius: float) -> None:
        """Remove the exact annulus a tool sweeps over one full turn of a guide circle.

        The AUTHORED guide radius is injected exactly; the annulus is
        ``[max(rho - r, 0), rho + r]`` formed in the exact kernel, a full disk
        when ``rho <= r``.

        Args:
            cx: X coordinate of the guide centre.
            cy: Y coordinate of the guide centre.
            guide_radius: Authored guide radius ``rho``.
            tool_radius: Tool radius ``r``.

        Raises:
            NonFiniteAnnulusInputError: A coordinate or radius is not finite.
            InvalidAnnulusRadiiError: ``rho < 0`` or ``r <= 0``.
        """
        self._raw.subtract_circle_sweep(float(cx), float(cy), float(guide_radius), float(tool_radius))

    def subtract_circle_sweep_local(self, cx: float, cy: float, guide_radius: float, tool_radius: float) -> None:
        """Remove the same annulus as `subtract_circle_sweep`, by local arrangement surgery.

        Args:
            cx: X coordinate of the guide centre.
            cy: Y coordinate of the guide centre.
            guide_radius: Authored guide radius ``rho``.
            tool_radius: Tool radius ``r``.

        Raises:
            NonFiniteAnnulusInputError: A coordinate or radius is not finite.
            InvalidAnnulusRadiiError: ``rho < 0`` or ``r <= 0``.
        """
        self._raw.subtract_circle_sweep_local(float(cx), float(cy), float(guide_radius), float(tool_radius))

    def line_sweep_removal(self, x0: float, y0: float, x1: float, y1: float, radius: float) -> SweepRemoval:
        """What the line sweep would remove from this stock, without removing it.

        Args:
            x0: X coordinate of the segment start (tool center).
            y0: Y coordinate of the segment start (tool center).
            x1: X coordinate of the segment end (tool center).
            y1: Y coordinate of the segment end (tool center).
            radius: Tool radius.

        Returns:
            The exact three-valued verdict; never WITHIN_BAND for an exact
            (Pythagorean or zero-length) motion.

        Raises:
            NonFiniteCapsuleInputError: A coordinate or the radius is not finite.
            InvalidLineSweepRadiusError: The radius is not positive.
        """
        return _sweep_removal(self._raw.line_sweep_removal(float(x0), float(y0), float(x1), float(y1), float(radius)))

    def circle_sweep_removal(self, cx: float, cy: float, guide_radius: float, tool_radius: float) -> SweepRemoval:
        """What one full machining circle would remove from this stock, without removing it.

        Args:
            cx: X coordinate of the guide centre.
            cy: Y coordinate of the guide centre.
            guide_radius: Authored guide radius ``rho``.
            tool_radius: Tool radius ``r``.

        Returns:
            REMOVES or NOTHING: the annulus is the swept region exactly.

        Raises:
            NonFiniteAnnulusInputError: A coordinate or radius is not finite.
            InvalidAnnulusRadiiError: ``rho < 0`` or ``r <= 0``.
        """
        return _sweep_removal(self._raw.circle_sweep_removal(float(cx), float(cy), float(guide_radius), float(tool_radius)))

    def can_remove_circle(
        self,
        previous: tuple[float, float, float],
        current: tuple[float, float, float],
        following: tuple[float, float, float],
        connector_samples: np.ndarray,
        tool_radius: float,
    ) -> bool:
        """Whether deleting the *current* loop provably loses no material.

        Sufficient, not necessary: True means the material the current loop's
        annulus meets is covered by the two neighbouring loops and the
        preserved connector polyline (swept with `subtract_line_sweep`, which
        never over-covers); False retains the loop without claiming a gap.

        Args:
            previous: ``(cx, cy, rho)`` of the preceding loop, world XY mm.
            current: ``(cx, cy, rho)`` of the loop proposed for deletion.
            following: ``(cx, cy, rho)`` of the succeeding loop.
            connector_samples: ``(N, 2)`` world-XY polyline, ``N >= 2``, of the
                route kept in place of the deleted loop.
            tool_radius: Tool radius.

        Returns:
            True when the deletion is proved to preserve material.

        Raises:
            InvalidCircleRemovalInputError: Malformed samples, non-finite input,
                a negative guide radius, or a non-positive tool radius.
        """
        samples = np.asarray(connector_samples, dtype=np.float64, order="C")
        return bool(
            self._raw.can_remove_circle(
                (float(previous[0]), float(previous[1]), float(previous[2])),
                (float(current[0]), float(current[1]), float(current[2])),
                (float(following[0]), float(following[1]), float(following[2])),
                samples,
                float(tool_radius),
            )
        )

    def subtract_disk(self, cx: float, cy: float, radius: float) -> None:
        """Remove the disk swept by a tool of *radius* plunging at ``(cx, cy)``.

        Args:
            cx: X coordinate of the disk center (tool center).
            cy: Y coordinate of the disk center (tool center).
            radius: Tool radius (disk radius).
        """
        self._raw.subtract_disk(float(cx), float(cy), float(radius))

    def subtract_annulus(self, cx: float, cy: float, inner_radius: float, outer_radius: float) -> None:
        """Remove the exact annulus between *inner_radius* and *outer_radius* about ``(cx, cy)``.

        This is the swept region of a tool of radius ``r`` carried once around a
        circular guide of radius ``rho``, with ``inner_radius = rho - r`` and
        ``outer_radius = rho + r``. The region is removed exactly -- two boundary
        circles, no discretisation -- so nothing is left behind to compensate for.

        Args:
            cx: X coordinate of the annulus center (the guide center).
            cy: Y coordinate of the annulus center (the guide center).
            inner_radius: Inner bound of the band; ``0.0`` removes a full disk.
            outer_radius: Outer bound of the band; must exceed *inner_radius*.

        Raises:
            InvalidAnnulusRadiiError: If the radii do not satisfy
                ``outer_radius > inner_radius >= 0``.
            NonFiniteAnnulusInputError: If any coordinate or radius is not finite.
        """
        self._raw.subtract_annulus(float(cx), float(cy), float(inner_radius), float(outer_radius))

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

        A FULL turn -- ``(sx, sy) == (ex, ey)`` -- removes the swept annulus
        exactly, via `subtract_annulus`. A partial arc is still removed as a
        certified under-covering chain of tool disks, which retains slivers of
        material the tool has in fact removed and so reads high on engagement.

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

    def subtract_disk_local(self, cx: float, cy: float, radius: float) -> None:
        """Remove the same disk as `subtract_disk`, by local arrangement surgery.

        Identical result -- exact point set and canonical representation alike --
        reached by editing the arrangement around the disk instead of overlaying
        the whole stock. The cost then grows with the removed neighbourhood
        rather than with the whole stock, which pays off once the stock has been
        cut up; on a barely-cut stock the two are within noise of each other. See
        `docs/local_depletion.md` for the measured crossover.

        Args:
            cx: X coordinate of the disk center (tool center).
            cy: Y coordinate of the disk center (tool center).
            radius: Tool radius (disk radius).
        """
        self._raw.subtract_disk_local(float(cx), float(cy), float(radius))

    def subtract_annulus_local(self, cx: float, cy: float, inner_radius: float, outer_radius: float) -> None:
        """Remove the same annulus as `subtract_annulus`, by local arrangement surgery.

        Args:
            cx: X coordinate of the annulus center (the guide center).
            cy: Y coordinate of the annulus center (the guide center).
            inner_radius: Inner bound of the band; ``0.0`` removes a full disk.
            outer_radius: Outer bound of the band; must exceed *inner_radius*.

        Raises:
            InvalidAnnulusRadiiError: If the radii do not satisfy
                ``outer_radius > inner_radius >= 0``.
            NonFiniteAnnulusInputError: If any coordinate or radius is not finite.
        """
        self._raw.subtract_annulus_local(float(cx), float(cy), float(inner_radius), float(outer_radius))

    def subtract_arc_sweep_local(
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
        """Remove the same arc sweep as `subtract_arc_sweep`, locally where possible.

        A FULL turn takes the local annulus path. A PARTIAL arc is a disk-chain
        union, whose boundary the local update cannot identify exactly, so it
        falls through to `subtract_arc_sweep` unchanged -- same result, same cost.

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
        self._raw.subtract_arc_sweep_local(
            float(cx),
            float(cy),
            float(sx),
            float(sy),
            float(ex),
            float(ey),
            bool(cw),
            float(tool_radius),
        )

    def representation_is_valid(self) -> bool:
        """Does the backing boolean set still satisfy its representation invariant?

        Every edge must separate faces of DIFFERENT containment and must carry the
        contained side on its left. This is orthogonal to point-set equality: a
        set can be exactly the right region and still be stored non-canonically,
        which would silently degrade every later query.

        Returns:
            ``True`` if the underlying `General_polygon_set_2` is canonical.
        """
        return self._raw.representation_is_valid()

    def arrangement_stats(self) -> tuple[int, int, int]:
        """Feature counts of the underlying exact arrangement.

        Cheap -- pure counters with no exact evaluation -- so this is safe to read
        inside a timed measurement.

        Returns:
            ``(vertices, halfedges, faces)``.
        """
        return self._raw.arrangement_stats()

    def coordinate_digits(self) -> tuple[int, float, int]:
        """Printed decimal length of the arrangement's exact vertex coordinates.

        A proxy for the bit growth of the exact representation: only how these
        numbers grow across successive removals is meaningful, never their
        absolute magnitude.

        Warning:
            This forces exact evaluation of every sampled coordinate, collapsing
            the lazy filter and changing subsequent timings. Never read it inside
            a timed measurement; use a separate diagnostic pass.

        Returns:
            ``(max_digits, mean_digits, sampled)``.
        """
        return self._raw.coordinate_digits()
