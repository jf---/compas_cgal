"""One named exception per failure mode of the coverage certificate."""

from __future__ import annotations

from collections.abc import Sequence

from compas_cgal import _coverage_certificate_2 as native
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.coverage_certificate.witness import UncoveredWitness


class InvalidCentreDomainError(ValueError):
    """The centre domain violates a polygon invariant (simple loops, outer CCW, holes CW and strictly inside, no zero-length edge)."""


class InvalidCentreBoundarySiteError(ValueError):
    """A boundary site has a parameter outside [0, 1), a negative index, or an index the domain does not have."""


class DegenerateConnectorError(ValueError):
    """A boundary connector makes no CCW progress (its start and end sites coincide)."""


class InvalidVertexDiskError(ValueError):
    """A vertex-disk centre is not a one-root point (non-rational parts or a negative radicand)."""


class CoverageCertificateNativeError(RuntimeError):
    """An internal contract breach in the native certificate; raised, never swallowed."""


class StartWitnessSearchError(CoverageCertificateNativeError):
    """The circle-sweep start-witness search exceeded its theorem bound: a support coincides with a circle unrecognised."""


class SignRuleContractError(CoverageCertificateNativeError):
    """An exact sign rule contradicted its own certified interval enclosure (the rule is wrong)."""


class SeparationBoundExceededError(CoverageCertificateNativeError):
    """An exact search loop ran past its root-separation bound: two values called strictly ordered are equal."""


class CoincidentTargetSupportError(CoverageCertificateNativeError):
    """A swept piece lies exactly on a cut support of the exact target (declared boundary or vertex circle): not modelled."""


def native_contract_error(error: Exception) -> CoverageCertificateNativeError:
    """The named Python error for a native contract breach (most specific class first)."""
    for native_type, python_type in (
        (native.StartWitnessSearchError, StartWitnessSearchError),
        (native.SignRuleContractError, SignRuleContractError),
        (native.SeparationBoundExceededError, SeparationBoundExceededError),
        (native.CoincidentTargetSupportError, CoincidentTargetSupportError),
    ):
        if isinstance(error, native_type):
            return python_type(str(error))
    return CoverageCertificateNativeError(str(error))


class UncoveredPocketError(Exception):
    """The swept regions do not cover the machinable target.

    Every witness is an exact point decided in T and outside every swept region,
    so the failure is a proof, not an estimate.

    Attributes:
        witnesses: One exact uncovered point per exposed edge that the
            certificate proved exposed.
        undecided_edge_count: Exposed sweep edges whose side could not be
            decided against T (only with VERTEX_DISK members; zero otherwise).
    """

    def __init__(self, witnesses: Sequence[UncoveredWitness], undecided_edge_count: int) -> None:
        self.witnesses = tuple(witnesses)
        self.undecided_edge_count = undecided_edge_count
        nearest = ", ".join(f"({x:.4f}, {y:.4f})" for x, y in (w.uncovered_point_mm for w in self.witnesses[:4]))
        super().__init__(f"{len(self.witnesses)} exposed edge(s) leave target material uncut, e.g. at {nearest}; {undecided_edge_count} more exposed edge(s) undecided")


class UnmodelledCertificatePieceError(Exception):
    """The certificate needs a piece type that is not implemented yet, so it cannot decide.

    Raised when VERTEX_DISK members (one-root centres of the exact centre domain
    Omega*) are given and no uncovered point was found: their boundary circles
    are not cut into the arrangement yet, so coverage cannot be certified.
    Never raised for full circles and segment connectors on a polygonal domain.

    Attributes:
        piece_kinds: The piece kinds the decision needs.
        undecided_edge_count: Exposed sweep edges left undecided.
    """

    def __init__(self, piece_kinds: Sequence[ExposedPieceKind], undecided_edge_count: int, reason: str) -> None:
        self.piece_kinds = tuple(piece_kinds)
        self.undecided_edge_count = undecided_edge_count
        names = ", ".join(kind.name for kind in self.piece_kinds)
        super().__init__(f"{reason} (needs {names}; {undecided_edge_count} exposed edge(s) undecided)")
