"""Failure modes of the coverage repair, one named exception each.

The base ``CoverageRepairError`` and the value-type invariant violations live
in ``held_coverage_repair_invariant_errors``.

Design: docs/held_coverage_preserving_selection.md.
"""

from __future__ import annotations

from collections.abc import Sequence

from benchmarks.held_coverage_repair_invariant_errors import CoverageRepairError
from benchmarks.held_coverage_repair_types import RepairRefusal
from benchmarks.held_coverage_repair_types import SiteKey
from compas_cgal.coverage_certificate.witness import UncoveredWitness


class InvalidCoverageRepairInputError(CoverageRepairError):
    """The repair received an emission inconsistent with its sources (a contract breach by the caller)."""


class CoverageNotRepairableBySelectionError(CoverageRepairError):
    """A witness lies in no dropped source: placement or the connector chain is at fault, not selection.

    Attributes:
        witnesses: Every witness of the round that no dropped source contains.
    """

    def __init__(self, witnesses: Sequence[UncoveredWitness]) -> None:
        self.witnesses = tuple(witnesses)
        shown = ", ".join(f"{w.piece_kind.name} ({w.uncovered_point_mm[0]:.4f}, {w.uncovered_point_mm[1]:.4f})" for w in self.witnesses[:4])
        super().__init__(f"{len(self.witnesses)} witness(es) lie in no dropped source, e.g. {shown}")


class CoverageRepairRefusedError(CoverageRepairError):
    """Every candidate that contains a witness was refused by the engagement cap or the circle budget.

    Attributes:
        witness: The witness left without an admissible candidate.
        refusals: Every refusal recorded by the repair so far.
    """

    def __init__(self, witness: UncoveredWitness, refusals: Sequence[RepairRefusal]) -> None:
        self.witness = witness
        self.refusals = tuple(refusals)
        names = ", ".join(f"#{r.candidate}" for r in self.refusals)
        x, y = witness.uncovered_point_mm
        super().__init__(f"every candidate covering ({x:.4f}, {y:.4f}) was refused: {names}")


class RepairOrderViolationError(CoverageRepairError):
    """A candidate's exact boundary site is not between its kept neighbours in cyclic CCW order."""

    def __init__(self, candidate: int, predecessor: int, successor: int, sites: tuple[SiteKey, SiteKey, SiteKey]) -> None:
        self.candidate = candidate
        self.predecessor = predecessor
        self.successor = successor
        self.sites = sites
        super().__init__(f"candidate #{candidate} at site {sites[1]} is not between #{predecessor} {sites[0]} and #{successor} {sites[2]}")


class UncoveredPocketWithoutWitnessesError(CoverageRepairError):
    """The oracle raised ``UncoveredPocketError`` with no witness, so no repair can make progress.

    The real certificate never does this (it raises only on exposed edges,
    each with a witness); the ``CoverageOracle`` protocol permits it.

    Attributes:
        undecided_edge_count: The oracle's undecided exposed-edge count.
    """

    def __init__(self, undecided_edge_count: int) -> None:
        self.undecided_edge_count = undecided_edge_count
        super().__init__(f"the coverage oracle reported an uncovered pocket without witnesses ({undecided_edge_count} undecided edge(s)); repair cannot progress")
