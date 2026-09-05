"""Bounded repository path-length measurements for the Figure 6 comparison."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from benchmarks.errors import BenchmarkError
from benchmarks.figure6 import controlled_path
from benchmarks.figure6 import reference_pocket
from benchmarks.pathmetrics import path_length
from benchmarks.units import Degrees
from compas_cgal.adaptive.units import Millimetre

REPOSITORY_FIGURE6_CAPS: tuple[Degrees, ...] = (Degrees(80.0), Degrees(120.0), Degrees(160.0))
COMPLIANCE_UNAUDITED: Literal["compliance unaudited"] = "compliance unaudited"


class InvalidFigure6RepositoryMeasurementError(BenchmarkError):
    """A bounded repository Figure 6 measurement is invalid."""


@dataclass(frozen=True)
class RepositoryFigure6Point:
    requested_cap_deg: Degrees
    path_length_mm: Millimetre
    compliance: Literal["compliance unaudited"]

    @classmethod
    def build(cls, requested_cap_deg: Degrees, path_length_mm: Millimetre) -> RepositoryFigure6Point:
        cap = float(requested_cap_deg)
        length = float(path_length_mm)
        if not math.isfinite(cap) or not math.isfinite(length) or cap <= 0.0 or length <= 0.0:
            raise InvalidFigure6RepositoryMeasurementError("Repository Figure 6 cap and path length must be finite and positive.")
        return cls(Degrees(cap), Millimetre(length), COMPLIANCE_UNAUDITED)


@dataclass(frozen=True)
class RepositoryFigure6Comparison:
    points: tuple[RepositoryFigure6Point, ...]
    constant_spacing_available: Literal[False]

    @classmethod
    def build(cls, points: tuple[RepositoryFigure6Point, ...]) -> RepositoryFigure6Comparison:
        if tuple(point.requested_cap_deg for point in points) != REPOSITORY_FIGURE6_CAPS:
            raise InvalidFigure6RepositoryMeasurementError("Repository Figure 6 requires requested caps 80, 120, and 160 degrees in order.")
        return cls(points, False)


def measure_repository_figure6() -> RepositoryFigure6Comparison:
    """Generate three real paths and measure analytic length without any audit."""
    spec = reference_pocket()
    points = tuple(RepositoryFigure6Point.build(cap, Millimetre(path_length(controlled_path(spec, float(cap))))) for cap in REPOSITORY_FIGURE6_CAPS)
    return RepositoryFigure6Comparison.build(points)
