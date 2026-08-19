"""Probes for the exact kernel's two hidden cost drivers.

Wall time in an exact-constructions kernel is driven by how many arrangement
features exist AND by how many bits the exact rationals carry, because chained
boolean constructions grow both. Only the first is visible from outside; this
module exposes the second so a measured slowdown can be attributed instead of
guessed at.
"""

from __future__ import annotations

from dataclasses import dataclass

from compas_cgal.stock import Stock


@dataclass(frozen=True)
class ArrangementSize:
    """Feature counts of an exact stock arrangement.

    Attributes:
        vertices: Number of arrangement vertices.
        halfedges: Number of arrangement halfedges.
        faces: Number of arrangement faces.
    """

    vertices: int
    halfedges: int
    faces: int


@dataclass(frozen=True)
class CoordinateDigits:
    """Decimal-length statistics of exact coordinates.

    Attributes:
        max_digits: Longest exact rational encountered.
        mean_digits: Mean length across all sampled coordinate parts.
        sampled: Number of coordinate parts inspected.
    """

    max_digits: int
    mean_digits: float
    sampled: int


def probe_size(stock: Stock) -> ArrangementSize:
    """Read arrangement feature counts. Safe inside a timed run.

    Args:
        stock: The stock to inspect.

    Returns:
        The feature counts.
    """
    vertices, halfedges, faces = stock.arrangement_stats()
    return ArrangementSize(vertices=vertices, halfedges=halfedges, faces=faces)


def probe_digits(stock: Stock) -> CoordinateDigits:
    """Read exact-coordinate digit statistics. NEVER call inside a timed run.

    Calling this forces exact evaluation of each sampled coordinate, which
    collapses the lazy-exact filter and inflates every later operation.

    Args:
        stock: The stock to inspect.

    Returns:
        The digit statistics.
    """
    max_digits, mean_digits, sampled = stock.coordinate_digits()
    return CoordinateDigits(max_digits=max_digits, mean_digits=mean_digits, sampled=sampled)
