"""Exposed-piece kinds and the exact uncovered witness a failed certificate carries."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Millimetre


class ExposedPieceKind(Enum):
    """The boundary piece an exposed edge lies on."""

    ANNULUS_OUTER = "annulus_outer"
    ANNULUS_INNER = "annulus_inner"
    DISK = "disk"
    CAPSULE_CAP = "capsule_cap"
    CAPSULE_FLANK = "capsule_flank"
    TARGET_BOUNDARY = "target_boundary"
    TARGET_INTERIOR = "target_interior"
    CONNECTOR_ARC = "connector_arc"


@dataclass(frozen=True)
class UncoveredWitness:
    """An exact point of T outside every swept region, beside one exposed edge.

    Attributes:
        piece_kind: Kind of the exposed piece.
        region_index: Index into ``full_circles``, then ``connectors``, then
            ``exact_connectors``, then the target members, of the region the
            exposed piece bounds.
        uncovered_point: Exact world-XY point, decided in T and outside every
            swept region by exact predicates (P7).
    """

    piece_kind: ExposedPieceKind
    region_index: int
    uncovered_point: tuple[ExactMillimetre, ExactMillimetre]

    @property
    def uncovered_point_mm(self) -> tuple[Millimetre, Millimetre]:
        """The uncovered point rounded to binary64, for reporting only."""
        return Millimetre(float(self.uncovered_point[0])), Millimetre(float(self.uncovered_point[1]))
