"""``CoverageCertificate.build``: decide exactly whether the swept regions cover the machinable target T.

Phase 2 (docs/held_coverage_certificate.md, "Implementation status"): cut
circles, connector capsules and the target pieces are all cut into the
arrangement, so on a polygonal centre domain the certificate decides both ways:
it returns a ``CoverageCertificate`` when no edge is exposed, and raises
``UncoveredPocketError`` with one exact witness per exposed edge otherwise.

The exact target T = Omega* (+) D_r of the exact centre domain
(``exact_target``, docs/held_coverage_certificate.md "Exact target Omega*") is
cut in completely: its declared boundary supports and its one-root vertex
circles cut every swept piece, so the certificate decides both ways and never
raises ``UnmodelledCertificatePieceError``. The polygon ``domain`` then only
carries the polygon connectors' sites; exact connectors along dOmega*
(``exact_connectors``) are swept exactly, arcs included
(docs/held_coverage_certificate.md "Exact connectors along dOmega*").

Loose VERTEX_DISK members on the polygonal target (``vertex_disks``) enter T's
membership but are not cut into the arrangement: with them, a proven uncovered
point still raises ``UncoveredPocketError``, and the absence of one raises
``UnmodelledCertificatePieceError`` instead of a pass.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.errors import UnmodelledCertificatePieceError
from compas_cgal.coverage_certificate.exact_connectors import ExactConnectorSweep
from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.target import VertexDisk
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.units import ToolRadius


@dataclass(frozen=True)
class CoverageCertificate:
    """A passed certificate: T is covered by the swept regions. Reporting fields only."""

    region_count: int
    piece_count: int
    edge_count: int
    exposed_edge_count: int
    exact_fallback_count: int
    elapsed_s: float

    @classmethod
    def build(
        cls,
        domain: CentreDomain,
        tool_radius: ToolRadius,
        full_circles: Sequence[FullCircleSweep],
        connectors: Sequence[BoundaryConnectorSweep] = (),
        *,
        vertex_disks: Sequence[VertexDisk] = (),
        exact_target: ExactCentreDomainTarget | None = None,
        exact_connectors: Sequence[ExactConnectorSweep] = (),
    ) -> CoverageCertificate:
        """Decide T subset-of U exactly; never returns partial success.

        Args:
            domain: Centre domain Omega.
            tool_radius: Tool radius r.
            full_circles: Full-circle cuts, region indices 0..n-1.
            connectors: Boundary connectors, region indices n..n+m-1.
            vertex_disks: Extra target members D_r(w) of the polygonal target
                (not cut into the arrangement; see the module docstring).
            exact_target: The exact target T = Omega* (+) D_r; ``domain`` then
                only carries the connectors' sites.
            exact_connectors: Exact dOmega* connectors along the exact target
                (region indices after ``connectors``); their sub-pieces must
                cover dOmega*.

        Raises:
            UncoveredPocketError: an exact point of T lies outside every sweep.
            UnmodelledCertificatePieceError: loose vertex disks were given and
                no uncovered point was found (never with ``exact_target``).
            InvalidExactTargetError: vertex disks given beside an exact target.
            ConnectorBoundaryCoverageError: the exact connectors leave part of
                dOmega* unswept.
            BoundaryChainError: the exact target's dOmega* pieces do not form
                closed chains with non-reflex junctions (exact connectors only).
            CoincidentTargetSupportError: an irrational connector flank lies on
                an exact connector support.
            CoverageCertificateNativeError: an internal contract breach.
        """
        report = run_certificate(domain, tool_radius, full_circles, connectors, vertex_disks=vertex_disks, exact_target=exact_target, exact_connectors=exact_connectors)
        if report.exposed:
            witnesses = [UncoveredWitness(edge.piece_kind, edge.region_index, edge.uncovered) for edge in report.exposed]
            raise UncoveredPocketError(witnesses, report.undecided)
        if vertex_disks:
            raise UnmodelledCertificatePieceError(
                (ExposedPieceKind.TARGET_BOUNDARY,),
                report.undecided,
                f"{len(vertex_disks)} VERTEX_DISK member(s) given; their one-root boundary circles are not cut into the arrangement yet, so coverage of T cannot be certified",
            )
        return cls(
            region_count=report.region_count,
            piece_count=report.piece_count,
            edge_count=report.edge_count,
            exposed_edge_count=len(report.exposed),
            exact_fallback_count=report.counters["exact_fallbacks"],
            elapsed_s=report.elapsed_s,
        )
