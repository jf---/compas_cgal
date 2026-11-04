"""The two polygons of a Held polygon-lane path, typed apart: the material domain and the placement domain.

- ``PartPolygon`` P is the **material domain**: the part's projected boundary
  (``case.projection.points``). Every actual-stock decision (the last-16
  certificate, the circle view C_t, the future refinement engine) decides on
  stock that starts as P.
- ``CentreDomainPolygon`` is the **placement domain**: the one-tool-radius
  inward offset of P (``figure7_inward_offset(case).components[0]``), where
  tool centres, midpoints and connectors live.

Deciding on the offset silently drops the r-wide wall band that station 0
reads on almost every circle (design 2026-11-04, review H1), so the two are
distinct nominal types: ``mypy --strict`` refuses a swap, and
``require_part_material`` checks, exactly, that a stock's rows are P's.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from benchmarks.held_reference_cases import HeldReferenceCase
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# A simple polygon needs this many vertices.
MIN_POLYGON_VERTICES = 3
# Stock2 takes (x, y, z) rows; the material domain is planar at z = 0.
STOCK_ROW_Z = 0.0


class MaterialDomainError(Exception):
    """Base class of the material-domain failures."""


class DegeneratePolygonError(MaterialDomainError):
    """A domain polygon has fewer than three vertices."""


class MaterialDomainMismatchError(MaterialDomainError):
    """A stock's polygon is not the part polygon P (for example the tool-centre offset was passed)."""


def _validated(points: Sequence[Point2[WorldXY]], kind: str) -> tuple[Point2[WorldXY], ...]:
    ordered = tuple(points)
    if len(ordered) < MIN_POLYGON_VERTICES:
        raise DegeneratePolygonError(f"{kind}: {len(ordered)} vertices, a polygon needs {MIN_POLYGON_VERTICES}")
    return ordered


def _rows(points: tuple[Point2[WorldXY], ...]) -> NDArray[np.float64]:
    return np.array([(float(p.x), float(p.y), STOCK_ROW_Z) for p in points], dtype=np.float64).reshape(-1, 3)


@dataclass(frozen=True)
class PartPolygon:
    """The material domain P: the part's projected boundary, in world XY."""

    points: tuple[Point2[WorldXY], ...]

    @classmethod
    def build(cls, points: Sequence[Point2[WorldXY]]) -> PartPolygon:
        """Validate and freeze P.

        Raises:
            DegeneratePolygonError: fewer than three vertices.
        """
        return cls(_validated(points, cls.__name__))

    @classmethod
    def from_case(cls, case: HeldReferenceCase) -> PartPolygon:
        """P of a reference case: ``case.projection.points``, the certificate's material source."""
        return cls.build(case.projection.points)

    @classmethod
    def from_stock_rows(cls, rows: NDArray[np.float64]) -> PartPolygon:
        """P read back from ``Stock2`` input rows ``(x, y, z)``."""
        return cls.build([Point2[WorldXY].build(float(x), float(y)) for x, y, *_z in rows.tolist()])

    def stock_rows(self) -> NDArray[np.float64]:
        """``Stock2`` input rows ``(x, y, 0)``."""
        return _rows(self.points)


@dataclass(frozen=True)
class CentreDomainPolygon:
    """The placement domain: the tool-centre (inward offset) polygon, in world XY. Never a material domain."""

    points: tuple[Point2[WorldXY], ...]

    @classmethod
    def build(cls, points: Sequence[Point2[WorldXY]]) -> CentreDomainPolygon:
        """Validate and freeze the centre domain.

        Raises:
            DegeneratePolygonError: fewer than three vertices.
        """
        return cls(_validated(points, cls.__name__))


def require_part_material(rows: NDArray[np.float64], part: PartPolygon, where: str) -> None:
    """Check exactly (binary64 equality, no tolerance) that stock rows ``rows`` are P's vertices in order.

    Raises:
        MaterialDomainMismatchError: the rows differ from P.
    """
    expected = part.stock_rows()
    if rows.shape != expected.shape or not np.array_equal(rows, expected):
        raise MaterialDomainMismatchError(f"{where}: stock polygon ({rows.shape[0]} vertices) is not the part polygon P ({expected.shape[0]} vertices)")
