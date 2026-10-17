"""Native exact station-boundary analysis beside the Fraction path in held_figure5_raw_guide.

Not wired into the generator: this module exists so the equivalence test and the
timing benchmark can run both paths over the same guide.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import RationalPoint
from compas_cgal.held_boundary_clearance import HeldBoundaryClearanceIndex
from compas_cgal.units import Millimetre


@dataclass(frozen=True)
class NativeStationBoundaryAnalysis:
    """Admissible sites and their exact records; non-admissible records are never built."""

    distance_admissible_sites: tuple[ProjectionBoundarySite, ...]
    admissible_records: tuple[tuple[ProjectionBoundarySite, RationalPoint, Fraction], ...]


def native_station_boundary_analyses(
    guide: Figure5RawGuide,
) -> tuple[tuple[NativeStationBoundaryAnalysis, ...], ...]:
    """Analyse every station of every run in one native batch, indexed [run][ordinal]."""
    index = HeldBoundaryClearanceIndex.build(guide.case.projection.points)
    stations = tuple(station.middle_point for run in guide.runs for station in run.stations)
    # The Fraction path injects this same double; the kernel injects it exactly.
    gap = Millimetre(float(guide.site_budget.admissible_distance_gap))
    rows = iter(index.distance_admissible_sites(stations, gap))
    segment_count = guide.boundary_segment_count
    analyses: list[tuple[NativeStationBoundaryAnalysis, ...]] = []
    for run in guide.runs:
        run_analyses: list[NativeStationBoundaryAnalysis] = []
        for _ in run.stations:
            records = tuple(
                (
                    ProjectionBoundarySite.build(ProjectionBoundarySegmentId(site.segment), site.parameter, segment_count),
                    (site.footpoint[0], site.footpoint[1]),
                    site.squared_distance,
                )
                for site in next(rows)
            )
            run_analyses.append(NativeStationBoundaryAnalysis(tuple(record[0] for record in records), records))
        analyses.append(tuple(run_analyses))
    return tuple(analyses)
