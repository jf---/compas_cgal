"""The generator's native station-boundary batch, exposed beside the Fraction path.

The conversion is the wired one (``held_figure5_raw_guide._native_station_boundary_analyses``),
not a copy: this module only re-shapes it so the equivalence test and the timing
benchmark can run the native batch and the Fraction oracle over the same guide.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from benchmarks import held_figure5_raw_guide as raw_guide
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import RationalPoint


@dataclass(frozen=True)
class NativeStationBoundaryAnalysis:
    """Admissible sites and their exact records; non-admissible records are never built."""

    distance_admissible_sites: tuple[ProjectionBoundarySite, ...]
    admissible_records: tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]


def native_station_boundary_analyses(
    guide: Figure5RawGuide,
) -> tuple[tuple[NativeStationBoundaryAnalysis, ...], ...]:
    """Run the generator's native batch afresh (no guide memo), indexed [run][ordinal]."""
    return tuple(
        tuple(NativeStationBoundaryAnalysis(analysis.distance_admissible_sites, analysis.records) for analysis in run)
        for run in raw_guide._native_station_boundary_analyses(guide)
    )
