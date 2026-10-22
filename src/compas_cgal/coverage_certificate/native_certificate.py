"""Typed wrapper over the native full certificate (circles, connector capsules, target pieces).

The native run cuts every swept boundary (cut circles, capsule end disks,
capsule flanks) and every target piece (edge-rectangle flanks, vertex circles)
into the arrangement, sweeps the cover counters along each, and returns one
exact uncovered witness per exposed edge (docs/held_coverage_certificate.md).
``CoverageCertificate.build`` is the public entry point; this module also
exposes the filter mode and the audit oracle for the contract tests.
"""

from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

from compas_cgal import _coverage_certificate_2
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.errors import native_contract_error
from compas_cgal.coverage_certificate.exact_connectors import ExactConnectorSweep
from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
from compas_cgal.coverage_certificate.exact_target import InvalidExactTargetError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.target import VertexDisk
from compas_cgal.coverage_certificate.target import loop_coordinates
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import ToolRadius

ExactPoint = tuple[ExactMillimetre, ExactMillimetre]

# Native CertificatePieceKind values (coverage_certificate_2/certificate_2.h).
_PIECE_KINDS = {
    0: ExposedPieceKind.ANNULUS_OUTER,
    1: ExposedPieceKind.ANNULUS_INNER,
    2: ExposedPieceKind.DISK,
    3: ExposedPieceKind.CAPSULE_CAP,
    4: ExposedPieceKind.CAPSULE_FLANK,
    5: ExposedPieceKind.TARGET_BOUNDARY,
    6: ExposedPieceKind.TARGET_INTERIOR,
    7: ExposedPieceKind.CONNECTOR_ARC,
}


@dataclass(frozen=True)
class ExposedPiece:
    """One exposed edge: an exact proof that T is not covered.

    Attributes:
        piece_kind: Kind of the exposed piece.
        region_index: Full circles, then connectors, then target members
            (1 + domain edge, loop-major) of the region the piece bounds.
        piece_index: Native piece id (boundary circle, then flank piece).
        edge_index: Edge index on the piece.
        parameter: Chart parameter of the exact on-edge witness.
        uncovered: Exact point in T outside every swept region.
    """

    piece_kind: ExposedPieceKind
    region_index: int
    piece_index: int
    edge_index: int
    parameter: Fraction
    uncovered: ExactPoint


@dataclass(frozen=True)
class CertificateReport:
    """Result of the native certificate; ``counters`` and timings are reporting only."""

    exposed: tuple[ExposedPiece, ...]
    undecided: int
    counters: Mapping[str, int]
    region_count: int
    piece_count: int
    edge_count: int
    elapsed_s: float
    phase_s: tuple[float, float, float, float]

    def canonical_text(self) -> str:
        """Deterministic text of every decision (exposed edges, exact witnesses), excluding counters and timings."""
        lines = [f"regions={self.region_count} pieces={self.piece_count} edges={self.edge_count} undecided={self.undecided}"]
        for e in self.exposed:
            lines.append(f"{e.piece_index} {e.region_index} {e.piece_kind.name} {e.edge_index} {e.parameter} {e.uncovered[0]} {e.uncovered[1]}")
        return "\n".join(lines) + "\n"


def _fraction(text: tuple[str, str]) -> Fraction:
    return Fraction(int(text[0]), int(text[1]))


def _wire(value: Fraction) -> tuple[str, str]:
    return str(value.numerator), str(value.denominator)


def run_certificate(
    domain: CentreDomain,
    tool_radius: ToolRadius,
    full_circles: Sequence[FullCircleSweep],
    connectors: Sequence[BoundaryConnectorSweep] = (),
    *,
    vertex_disks: Sequence[VertexDisk] = (),
    exact_target: ExactCentreDomainTarget | None = None,
    exact_connectors: Sequence[ExactConnectorSweep] = (),
    filter_mode: FilterMode = FilterMode.LAZY,
    audit: bool = False,
) -> CertificateReport:
    """Run the native certificate.

    Args:
        domain: Centre domain Omega; with the tool radius it defines T.
        tool_radius: Tool radius r.
        full_circles: Full-circle cuts; region i is ``full_circles[i]``.
        connectors: Boundary connectors; region ``len(full_circles) + j``.
        vertex_disks: Extra one-root target members of the polygonal target
            (their boundary is not cut into the arrangement: a run with them
            reports undecided edges). Not allowed with ``exact_target``.
        exact_target: The exact target T = Omega* (+) D_r. When given, ``domain``
            only carries the connectors' sites and T is the exact target,
            decided exactly both ways (exact_certificate_2.h).
        exact_connectors: Exact dOmega* connectors (``exact_target`` only); region
            ``len(full_circles) + len(connectors) + k``. Their sub-pieces must
            cover dOmega*.
        filter_mode: LAZY (production) or EXACT_ONLY (contract tests).
        audit: Also confirm every edge the counters call covered by pure exact
            membership at its own witness (slow).

    Raises:
        InvalidExactTargetError: vertex disks given beside an exact target.
        InvalidExactTargetError: vertex disks given beside an exact target, or
            exact connectors without one.
        ConnectorBoundaryCoverageError: the exact connectors do not cover dOmega*.
        BoundaryChainError: the exact target's dOmega* pieces do not form closed
            chains with non-reflex junctions (exact connectors only).
        CoincidentTargetSupportError: an irrational connector flank lies on an
            exact connector support.
        CoverageCertificateNativeError: an internal contract breach.
    """
    rows = [(c.start.loop, c.start.segment, _wire(c.start.parameter), c.end.segment, _wire(c.end.parameter)) for c in connectors]
    if exact_target is not None and vertex_disks:
        raise InvalidExactTargetError("vertex disks belong to the exact target; pass them in ExactCentreDomainTarget")
    if exact_target is None and exact_connectors:
        raise InvalidExactTargetError("exact connectors run along dOmega* and need the exact target")
    try:
        if exact_target is not None:
            boundary_rows, piece_rows, disk_rows = exact_target.wire()
            out = _coverage_certificate_2.run_exact_target_certificate(
                [(float(cut.centre.x), float(cut.centre.y)) for cut in full_circles],
                [float(cut.path_radius.value) for cut in full_circles],
                float(tool_radius.value),
                loop_coordinates(domain.outer),
                [loop_coordinates(hole) for hole in domain.holes],
                rows,
                boundary_rows,
                piece_rows,
                disk_rows,
                filter_mode is FilterMode.EXACT_ONLY,
                audit,
                [connector.wire() for connector in exact_connectors],
            )
        else:
            out = _coverage_certificate_2.run_certificate(
                [(float(cut.centre.x), float(cut.centre.y)) for cut in full_circles],
                [float(cut.path_radius.value) for cut in full_circles],
                float(tool_radius.value),
                loop_coordinates(domain.outer),
                [loop_coordinates(hole) for hole in domain.holes],
                rows,
                [disk.wire() for disk in vertex_disks],
                filter_mode is FilterMode.EXACT_ONLY,
                audit,
            )
    except (_coverage_certificate_2.CoverageCertificateNativeError, _coverage_certificate_2.InvalidCertificateInputError) as error:
        raise native_contract_error(error) from error
    exposed = tuple(
        ExposedPiece(
            piece_kind=_PIECE_KINDS[kind],
            region_index=region,
            piece_index=piece,
            edge_index=edge,
            parameter=_fraction(parameter),
            uncovered=(ExactMillimetre(_fraction(point[0])), ExactMillimetre(_fraction(point[1]))),
        )
        for kind, region, piece, edge, parameter, point in out["exposed"]
    )
    phase = out["phase_s"]
    return CertificateReport(
        exposed=exposed,
        undecided=out["undecided"],
        counters=dict(out["counters"]),
        region_count=out["region_count"],
        piece_count=out["piece_count"],
        edge_count=out["edge_count"],
        elapsed_s=out["elapsed_s"],
        phase_s=(phase[0], phase[1], phase[2], phase[3]),
    )
