"""Exact dOmega* connectors: tool-centre motion along sub-pieces of the exact target's own pieces.

A connector on the exact path runs CCW along dOmega* (docs/held_exact_centre_domain.md):
the rest of one piece, whole pieces, and the head of another, each a sub-piece of
an offset line or an offset circle with one-root end points. The sub-pieces of
all connectors must cover dOmega*; the certificate verifies this exactly
(``ConnectorBoundaryCoverageError`` otherwise), together with the declared
pieces forming closed chains without reflex junctions (``BoundaryChainError``).
Their union is then exactly dOmega* (+) D_r, so the certificate sweeps one
composite region per dOmega* PIECE -- a line piece as the strip between its two
flank lines plus its end disks, an arc piece as the annular sector between radii
R - r and R + r plus its end disks -- and checks their flank segments and
support-circle arcs (``coverage_certificate_2/exact_connectors_2.h``). The
sub-pieces are kept for validation and for naming the connector of a witness.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
from compas_cgal.coverage_certificate.exact_target import InvalidExactTargetError
from compas_cgal.coverage_certificate.exact_target import OneRootPoint
from compas_cgal.coverage_certificate.exact_target import RationalText

ExactPieceRow = tuple[int, bool, list[RationalText], list[RationalText]]


@dataclass(frozen=True)
class ExactConnectorPiece:
    """Part of dOmega* piece ``piece`` (an index into the exact target's pieces) from ``start`` to ``end`` in travel order.

    ``full`` marks a whole closed circle piece travelled from its source back to it.
    """

    piece: int
    start: OneRootPoint
    end: OneRootPoint
    full: bool

    @classmethod
    def build(cls, target: ExactCentreDomainTarget, piece: int, start: OneRootPoint, end: OneRootPoint, *, full: bool = False) -> ExactConnectorPiece:
        """Raises:
        InvalidExactTargetError: the piece index is not one of the target's pieces.
        """
        if isinstance(piece, bool) or not 0 <= piece < len(target.pieces):
            raise InvalidExactTargetError(f"connector sub-piece names dOmega* piece {piece!r} of {len(target.pieces)}")
        return cls(piece, start, end, full)

    def wire(self) -> ExactPieceRow:
        return self.piece, self.full, self.start.wire(), self.end.wire()


@dataclass(frozen=True)
class ExactConnectorSweep:
    """One exact connector: its sub-pieces in travel order (a fan has none and sweeps nothing beyond its circles)."""

    pieces: tuple[ExactConnectorPiece, ...]

    @classmethod
    def build(cls, pieces: Sequence[ExactConnectorPiece]) -> ExactConnectorSweep:
        """Build a connector from its sub-pieces (validated against the target by ``ExactConnectorPiece.build``)."""
        return cls(tuple(pieces))

    def wire(self) -> list[ExactPieceRow]:
        return [piece.wire() for piece in self.pieces]
