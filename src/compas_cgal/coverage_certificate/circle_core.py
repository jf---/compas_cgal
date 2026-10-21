"""Typed wrapper over the native circle-circle core (phase 1 of the certificate).

For every boundary circle of every full-circle cut the native core finds the
crossings (certified interval front, Circular_kernel_2 exact back end), sweeps
the cover counters, and on every edge no other sweep covers builds an exact
on-edge witness (P3, confirmed by P4) and an uncovered point beside it (P7),
whose membership in T it decides exactly.
"""

from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction

from compas_cgal import _coverage_certificate_2
from compas_cgal.coverage_certificate.errors import native_contract_error
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.target import VertexDisk
from compas_cgal.coverage_certificate.target import loop_coordinates
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import ToolRadius

ExactPoint = tuple[ExactMillimetre, ExactMillimetre]

# Native PieceKind enum values (coverage_certificate_2/regions_2.h).
_PIECE_KINDS = {0: ExposedPieceKind.ANNULUS_OUTER, 1: ExposedPieceKind.ANNULUS_INNER, 2: ExposedPieceKind.DISK}


class FilterMode(Enum):
    """How the native core takes each filtered decision.

    LAZY is the production mode: certified interval filter, exact predicate on
    overlap. EXACT_ONLY sends every decision straight to the exact predicate;
    it exists so contract tests can hold LAZY to byte-identical output.
    """

    LAZY = "lazy"
    EXACT_ONLY = "exact_only"


@dataclass(frozen=True)
class ExposedCircleEdge:
    """One edge of a full-circle boundary that no other sweep covers, with exact witnesses."""

    circle_index: int
    region_index: int
    piece_kind: ExposedPieceKind
    edge_index: int
    whole_circle: bool
    start_vertex: tuple[int, ...]
    end_vertex: tuple[int, ...]
    parameter: Fraction
    on_edge: ExactPoint
    uncovered: ExactPoint
    uncovered_in_target: bool


@dataclass(frozen=True)
class CircleCoreReport:
    """Result of the circle-circle core; ``counters`` and timings are reporting only."""

    exposed: tuple[ExposedCircleEdge, ...]
    counters: Mapping[str, int]
    region_count: int
    piece_count: int
    elapsed_s: float
    phase_s: tuple[float, float, float, float]

    def canonical_text(self) -> str:
        """Deterministic text of every decision (exposed edges and exact witnesses), excluding counters and timings."""
        lines = [f"regions={self.region_count} pieces={self.piece_count} sub_arcs={self.counters['sub_arcs']}"]
        for e in self.exposed:
            lines.append(
                f"{e.circle_index} {e.region_index} {e.piece_kind.name} {e.edge_index} {int(e.whole_circle)} "
                f"{','.join(map(str, e.start_vertex))} {','.join(map(str, e.end_vertex))} {e.parameter} "
                f"{e.on_edge[0]} {e.on_edge[1]} {e.uncovered[0]} {e.uncovered[1]} {int(e.uncovered_in_target)}"
            )
        return "\n".join(lines) + "\n"


def _fraction(text: tuple[str, str]) -> Fraction:
    return Fraction(int(text[0]), int(text[1]))


def _point(text: tuple[tuple[str, str], tuple[str, str]]) -> ExactPoint:
    return ExactMillimetre(_fraction(text[0])), ExactMillimetre(_fraction(text[1]))


def run_circle_core(
    domain: CentreDomain,
    tool_radius: ToolRadius,
    full_circles: Sequence[FullCircleSweep],
    *,
    vertex_disks: Sequence[VertexDisk] = (),
    filter_mode: FilterMode = FilterMode.LAZY,
    audit: bool = False,
) -> CircleCoreReport:
    """Run the native circle-circle core.

    Args:
        domain: Centre domain Omega; with the tool radius it defines T.
        tool_radius: Tool radius r.
        full_circles: Full-circle cuts; region i is ``full_circles[i]``.
        vertex_disks: Extra one-root target members.
        filter_mode: LAZY (production) or EXACT_ONLY (contract tests).
        audit: Also confirm every edge the counters call covered by pure exact
            membership at its own witness, an oracle independent of the
            counters (slow).

    Raises:
        CoverageCertificateNativeError: an internal contract breach.
    """
    try:
        out = _coverage_certificate_2.run_circle_core(
            [(float(cut.centre.x), float(cut.centre.y)) for cut in full_circles],
            [float(cut.path_radius.value) for cut in full_circles],
            float(tool_radius.value),
            loop_coordinates(domain.outer),
            [loop_coordinates(hole) for hole in domain.holes],
            [disk.wire() for disk in vertex_disks],
            filter_mode is FilterMode.EXACT_ONLY,
            audit,
        )
    except (_coverage_certificate_2.CoverageCertificateNativeError, _coverage_certificate_2.InvalidCertificateInputError) as error:
        raise native_contract_error(error) from error
    exposed = tuple(
        ExposedCircleEdge(
            circle_index=circle,
            region_index=region,
            piece_kind=_PIECE_KINDS[kind],
            edge_index=edge,
            whole_circle=whole,
            start_vertex=tuple(start),
            end_vertex=tuple(end),
            parameter=_fraction(parameter),
            on_edge=_point(on_edge),
            uncovered=_point(uncovered),
            uncovered_in_target=in_target,
        )
        for circle, region, kind, edge, whole, start, end, parameter, on_edge, uncovered, in_target in out["exposed"]
    )
    phase = out["phase_s"]
    return CircleCoreReport(
        exposed=exposed,
        counters=dict(out["counters"]),
        region_count=out["region_count"],
        piece_count=out["piece_count"],
        elapsed_s=out["elapsed_s"],
        phase_s=(phase[0], phase[1], phase[2], phase[3]),
    )
