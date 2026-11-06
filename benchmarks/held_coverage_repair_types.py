"""Value types of the coverage repair: settings, emission, splice, ledger, certified path.

Each type's ``build`` factory validates its invariants and raises one named
``CoverageRepairInvariantError`` per violated invariant
(``held_coverage_repair_invariant_errors``); the raw constructors stay
minimal and unvalidated. Design: docs/held_coverage_preserving_selection.md.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Protocol

from benchmarks.held_coverage_repair_invariant_errors import CertificateRunCountError
from benchmarks.held_coverage_repair_invariant_errors import CertifiedLedgerMismatchError
from benchmarks.held_coverage_repair_invariant_errors import CertifiedNothingMismatchError
from benchmarks.held_coverage_repair_invariant_errors import CertifiedSelectionRangeError
from benchmarks.held_coverage_repair_invariant_errors import DuplicateInsertionError
from benchmarks.held_coverage_repair_invariant_errors import EmissionIndexCountError
from benchmarks.held_coverage_repair_invariant_errors import EmissionIndexOrderError
from benchmarks.held_coverage_repair_invariant_errors import EmittedCountError
from benchmarks.held_coverage_repair_invariant_errors import InvalidRepairTimingError
from benchmarks.held_coverage_repair_invariant_errors import LedgerRoundError
from benchmarks.held_coverage_repair_invariant_errors import NegativeRefinementDepthError
from benchmarks.held_coverage_repair_invariant_errors import NegativeSourceIndexError
from benchmarks.held_coverage_repair_invariant_errors import NonPositiveCircleBudgetError
from benchmarks.held_coverage_repair_invariant_errors import NonPositiveRepairRoundError
from benchmarks.held_coverage_repair_invariant_errors import SpliceSlotError
from benchmarks.held_coverage_repair_invariant_errors import UncertifiedPathError
from benchmarks.held_coverage_repair_invariant_errors import UnorderedSelectionError
from benchmarks.held_coverage_repair_invariant_errors import UntracedRefinementError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_refinement import MAX_REFINED_CIRCLES
from benchmarks.held_figure5_engagement_refinement import TracedEngagementRefinement
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_nothing_circles import NothingCirclePass
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius

DEFAULT_MAX_DEPTH = 16  # refine_figure5_engagement's own default bisection depth.
SiteKey = tuple[int, Fraction]


class CoverageOracle(Protocol):
    """A full coverage certificate: returns on T subset-of U, raises ``UncoveredPocketError`` otherwise."""

    def __call__(
        self,
        domain: CentreDomain,
        tool_radius: ToolRadius,
        full_circles: Sequence[FullCircleSweep],
        connectors: Sequence[BoundaryConnectorSweep],
        /,
    ) -> CoverageCertificate: ...


def _require_source_index(index: int, role: str) -> None:
    if index < 0:
        raise NegativeSourceIndexError(f"{role} source index {index} is negative")


def _require_round(round_index: int, role: str) -> None:
    if round_index < 1:
        raise NonPositiveRepairRoundError(f"{role} names round {round_index}; rounds are 1-based")


def _strictly_increasing(values: Sequence[int]) -> bool:
    return all(a < b for a, b in zip(values, values[1:]))


@dataclass(frozen=True)
class RepairRefusal:
    """One candidate the engagement cap or the circle budget refused; reporting fields.

    Attributes:
        candidate: Source index of the refused candidate.
        interval: Kept source indices (s-, s+) the candidate would split.
        round: Repair round (1-based).
        reason: The refinement's ``Figure5EngagementResolutionError`` message.
    """

    candidate: int
    interval: tuple[int, int]
    round: int
    reason: str

    @classmethod
    def build(cls, candidate: int, interval: tuple[int, int], round: int, reason: str) -> RepairRefusal:
        """Validated refusal.

        Raises:
            NegativeSourceIndexError: ``candidate`` or an interval source is negative.
            NonPositiveRepairRoundError: ``round`` is below one.
        """
        _require_source_index(candidate, "refused candidate")
        for index in interval:
            _require_source_index(index, "refusal interval")
        _require_round(round, "refusal")
        return cls(candidate, interval, round, reason)


@dataclass(frozen=True)
class EngagementRefinementSettings:
    """The refinement arguments every interval of one path is refined with.

    ``close_seam`` declares the path closed: its seam connector (last contact
    back to the first) is cut, so it is refined, swept by the coverage
    certificate and post-gated like every other connector. False is the open
    path (no seam connector anywhere). ``part`` is the material domain P the
    stock decision reads.
    """

    tool_radius: ToolRadius
    cap: EngagementCap
    part: PartPolygon
    max_depth: int = DEFAULT_MAX_DEPTH
    max_circles: int = MAX_REFINED_CIRCLES
    corner_approaches: bool = False
    close_seam: bool = True

    @classmethod
    def build(
        cls,
        tool_radius: ToolRadius,
        cap: EngagementCap,
        *,
        part: PartPolygon,
        max_depth: int = DEFAULT_MAX_DEPTH,
        max_circles: int = MAX_REFINED_CIRCLES,
        corner_approaches: bool = False,
        close_seam: bool = True,
    ) -> EngagementRefinementSettings:
        """Validated settings; ``tool_radius`` and ``cap`` are validated by their own factories.

        Raises:
            NegativeRefinementDepthError: ``max_depth`` is negative.
            NonPositiveCircleBudgetError: ``max_circles`` is below one.
        """
        if max_depth < 0:
            raise NegativeRefinementDepthError(f"max_depth {max_depth} is negative")
        if max_circles < 1:
            raise NonPositiveCircleBudgetError(f"max_circles {max_circles} is below one")
        return cls(tool_radius, cap, part, max_depth, max_circles, corner_approaches, close_seam)


@dataclass(frozen=True)
class BoundaryOrderedEmission:
    """A selection and its refined motion.

    ``refinement.original_indices[j]`` is the emitted index of source
    ``selected[j]``.
    """

    selected: tuple[int, ...]
    refinement: TracedEngagementRefinement

    @classmethod
    def build(cls, selected: Sequence[int], refinement: TracedEngagementRefinement) -> BoundaryOrderedEmission:
        """Validated emission.

        Raises:
            NegativeSourceIndexError: a selected source index is negative.
            UnorderedSelectionError: ``selected`` is not strictly increasing.
            EmissionIndexCountError: ``original_indices`` and ``selected``
                differ in length.
            EmissionIndexOrderError: ``original_indices`` is not strictly
                increasing within ``[0, len(circles))``.
            UntracedRefinementError: ``refinement`` carries no decision trace.
        """
        chosen = tuple(selected)
        if chosen:
            _require_source_index(min(chosen), "selected")
        if not _strictly_increasing(chosen):
            raise UnorderedSelectionError(f"selection {chosen[:8]} (first 8) is not strictly increasing")
        indices = refinement.original_indices
        if len(indices) != len(chosen):
            raise EmissionIndexCountError(f"{len(indices)} original indices for {len(chosen)} selected sources")
        if indices and (indices[0] < 0 or indices[-1] >= len(refinement.circles) or not _strictly_increasing(indices)):
            raise EmissionIndexOrderError(f"original indices are not strictly increasing within {len(refinement.circles)} emitted circles")
        if not isinstance(refinement, TracedEngagementRefinement):
            raise UntracedRefinementError(f"a {type(refinement).__name__} carries no decision trace")
        return cls(chosen, refinement)


@dataclass(frozen=True)
class Splice:
    """One source re-inserted into an emission; ``added`` is the net emitted-circle change."""

    emission: BoundaryOrderedEmission
    candidate: int
    interval: tuple[int, int]
    added: int

    @classmethod
    def build(cls, emission: BoundaryOrderedEmission, candidate: int, interval: tuple[int, int], added: int) -> Splice:
        """Validated splice; ``emission`` is the spliced (augmented) emission.

        Raises:
            SpliceSlotError: ``interval[0], candidate, interval[1]`` are not
                consecutive in ``emission.selected``.
        """
        selected = emission.selected
        slot = selected.index(candidate) if candidate in selected else -1
        if not 0 < slot < len(selected) - 1 or (selected[slot - 1], selected[slot + 1]) != interval:
            raise SpliceSlotError(f"candidate #{candidate} is not between #{interval[0]} and #{interval[1]} in the spliced selection")
        return cls(emission, candidate, interval, added)


@dataclass(frozen=True)
class RepairInsertion:
    """A re-inserted source, the first witness it was chosen for, and its round."""

    source: int
    witness: UncoveredWitness
    round: int

    @classmethod
    def build(cls, source: int, witness: UncoveredWitness, round: int) -> RepairInsertion:
        """Validated insertion.

        Raises:
            NegativeSourceIndexError: ``source`` is negative.
            NonPositiveRepairRoundError: ``round`` is below one.
        """
        _require_source_index(source, "inserted")
        _require_round(round, "insertion")
        return cls(source, witness, round)


@dataclass(frozen=True)
class CoverageRepairLedger:
    """Reporting fields of one select-certify-repair run.

    ``repair_s`` is the repair's own work (candidate search, set cover,
    splices); ``certificate_s`` is the wall time of every certificate run,
    the final fresh full certificate included.
    """

    insertions: tuple[RepairInsertion, ...]
    refusals: tuple[RepairRefusal, ...]
    rounds: int
    emitted_before: int
    emitted_after: int
    repair_s: float
    certificate_s: float
    certificate_runs: int

    @classmethod
    def build(
        cls,
        insertions: Sequence[RepairInsertion],
        refusals: Sequence[RepairRefusal],
        *,
        rounds: int,
        emitted_before: int,
        emitted_after: int,
        repair_s: float,
        certificate_s: float,
        certificate_runs: int,
    ) -> CoverageRepairLedger:
        """Validated ledger.

        Raises:
            CertificateRunCountError: ``rounds`` is negative or
                ``certificate_runs != rounds + 1``.
            LedgerRoundError: insertion rounds are out of order, exceed
                ``rounds``, or leave a round without an insertion; or refusal
                rounds are out of order or outside ``1..rounds``.
            DuplicateInsertionError: a source is inserted twice.
            EmittedCountError: an emitted count is negative, or the count
                changed with zero rounds.
            InvalidRepairTimingError: a timing is negative or not finite.
        """
        if rounds < 0 or certificate_runs != rounds + 1:
            raise CertificateRunCountError(f"{certificate_runs} certificate runs for {rounds} rounds; expected rounds + 1")
        inserted_rounds = [insertion.round for insertion in insertions]
        refused_rounds = [refusal.round for refusal in refusals]
        if sorted(inserted_rounds) != inserted_rounds or sorted(set(inserted_rounds)) != list(range(1, rounds + 1)):
            raise LedgerRoundError(f"insertion rounds {inserted_rounds[:8]} (first 8) do not cover rounds 1..{rounds} in order")
        if sorted(refused_rounds) != refused_rounds or any(not 1 <= r <= rounds for r in refused_rounds):
            raise LedgerRoundError(f"refusal rounds {refused_rounds[:8]} (first 8) are out of order or outside 1..{rounds}")
        sources = [insertion.source for insertion in insertions]
        if len(set(sources)) != len(sources):
            raise DuplicateInsertionError(f"inserted sources {sources[:8]} (first 8) repeat a source")
        if emitted_before < 0 or emitted_after < 0 or (rounds == 0 and emitted_after != emitted_before):
            raise EmittedCountError(f"emitted {emitted_before} -> {emitted_after} over {rounds} rounds")
        if not (math.isfinite(repair_s) and math.isfinite(certificate_s) and repair_s >= 0 and certificate_s >= 0):
            raise InvalidRepairTimingError(f"repair_s {repair_s}, certificate_s {certificate_s}")
        return cls(tuple(insertions), tuple(refusals), rounds, emitted_before, emitted_after, repair_s, certificate_s, certificate_runs)

    @property
    def inserted_sources(self) -> tuple[int, ...]:
        """Re-inserted source indices in insertion order."""
        return tuple(insertion.source for insertion in self.insertions)


@dataclass(frozen=True)
class CertifiedBoundaryOrderedPath:
    """A boundary-ordered path whose final full certificate holds.

    ``emission`` is the repair's working state: the selection and every circle
    its refinement emitted. ``nothing`` is the exact NOTHING pass on those
    circles; its kept circles (``circles``) are the emitted path, the one the
    certificate and the connector post-gate decided.
    """

    sources: tuple[Figure5CounterclockwiseCircle, ...]
    emission: BoundaryOrderedEmission
    nothing: NothingCirclePass[Figure5CounterclockwiseCircle]
    certificate: CoverageCertificate
    ledger: CoverageRepairLedger

    @classmethod
    def build(
        cls,
        sources: Sequence[Figure5CounterclockwiseCircle],
        emission: BoundaryOrderedEmission,
        nothing: NothingCirclePass[Figure5CounterclockwiseCircle],
        certificate: CoverageCertificate,
        ledger: CoverageRepairLedger,
    ) -> CertifiedBoundaryOrderedPath:
        """Validated certified path; ``emission`` and ``ledger`` are validated by their own factories.

        Raises:
            CertifiedSelectionRangeError: a selected index is not a source.
            CertifiedNothingMismatchError: ``nothing.removed`` is not strictly
                increasing within ``1..n-1`` (circle 0 is never removed),
                ``nothing.kept_for_successor_cap`` names a circle that is not
                kept, or ``nothing.circles`` is not the emission's circles minus
                ``nothing.removed``.
            CertifiedLedgerMismatchError: ``ledger.emitted_after`` is not the
                emitted count, or an inserted source is not selected.
            UncertifiedPathError: the certificate has exposed edges.
        """
        if emission.selected and max(emission.selected) >= len(sources):
            raise CertifiedSelectionRangeError(f"selected source #{max(emission.selected)} is not one of {len(sources)} sources")
        refined = emission.refinement.circles
        removed = nothing.removed
        if not _strictly_increasing(removed) or (removed and (removed[0] < 1 or removed[-1] >= len(refined))):
            raise CertifiedNothingMismatchError(f"removed indices {removed[:8]} (first 8) are not strictly increasing within 1..{len(refined) - 1} (circle 0 is never removed)")
        dropped = frozenset(removed)
        if not set(nothing.kept_for_successor_cap) <= set(range(len(refined))) - dropped:
            raise CertifiedNothingMismatchError(f"kept-for-cap indices {nothing.kept_for_successor_cap[:8]} (first 8) are not kept refined circles")
        if nothing.circles != tuple(circle for k, circle in enumerate(refined) if k not in dropped):
            raise CertifiedNothingMismatchError(f"{len(nothing.circles)} kept circles are not the {len(refined)} refined circles minus {len(removed)} removed")
        emitted = len(nothing.circles)
        if ledger.emitted_after != emitted or not set(ledger.inserted_sources) <= set(emission.selected):
            raise CertifiedLedgerMismatchError(f"ledger emitted_after {ledger.emitted_after} vs {emitted} emitted; inserted {ledger.inserted_sources[:8]} (first 8)")
        if certificate.exposed_edge_count != 0:
            raise UncertifiedPathError(f"certificate has {certificate.exposed_edge_count} exposed edge(s)")
        return cls(tuple(sources), emission, nothing, certificate, ledger)

    @property
    def circles(self) -> tuple[Figure5CounterclockwiseCircle, ...]:
        """The emitted path: the refined circles the NOTHING pass kept."""
        return self.nothing.circles
