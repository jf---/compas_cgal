"""Invariant violations of the coverage-repair value types, one named exception each.

Raised by the ``build`` factories in ``held_coverage_repair_types``; the raw
dataclass constructors stay unvalidated (bypass-safe). ``CoverageRepairError``,
the base of every coverage-repair failure, lives here so the value types can
raise without importing ``held_coverage_repair_errors`` (which imports them).
Design: docs/held_coverage_preserving_selection.md.
"""

from __future__ import annotations

from benchmarks.errors import BenchmarkError


class CoverageRepairError(BenchmarkError):
    """Base class of the coverage-repair failure modes."""


class CoverageRepairInvariantError(CoverageRepairError):
    """A coverage-repair value type was built with a violated invariant (a contract breach by the caller)."""


class NegativeRefinementDepthError(CoverageRepairInvariantError):
    """``EngagementRefinementSettings.max_depth`` is negative; the bisection depth is a count."""


class NonPositiveCircleBudgetError(CoverageRepairInvariantError):
    """``EngagementRefinementSettings.max_circles`` is below one; a path emits at least one circle."""


class NegativeSourceIndexError(CoverageRepairInvariantError):
    """A source index (selection, insertion, refusal) is negative."""


class UnorderedSelectionError(CoverageRepairInvariantError):
    """``BoundaryOrderedEmission.selected`` is not strictly increasing, so it is not a boundary-order selection."""


class EmissionIndexCountError(CoverageRepairInvariantError):
    """``refinement.original_indices`` does not index every selected source exactly once."""


class EmissionIndexOrderError(CoverageRepairInvariantError):
    """``refinement.original_indices`` is not strictly increasing within the emitted circles."""


class UntracedRefinementError(CoverageRepairInvariantError):
    """The emission's refinement carries no decision trace (not a ``TracedEngagementRefinement``), so it cannot be spliced."""


class SpliceSlotError(CoverageRepairInvariantError):
    """A splice's candidate does not sit between its interval's sources, consecutively, in the spliced selection."""


class NonPositiveRepairRoundError(CoverageRepairInvariantError):
    """An insertion or refusal names a repair round below one (rounds are 1-based)."""


class CertificateRunCountError(CoverageRepairInvariantError):
    """A ledger's ``rounds`` is negative or ``certificate_runs != rounds + 1``.

    Each round follows one failed certificate, and exactly one final
    certificate passes.
    """


class LedgerRoundError(CoverageRepairInvariantError):
    """A ledger's round attribution is inconsistent.

    An insertion or refusal round exceeds ``rounds``, rounds are out of order,
    or a round inserted nothing (every round inserts at least one source or
    raises).
    """


class DuplicateInsertionError(CoverageRepairInvariantError):
    """A ledger re-inserts the same source twice; an inserted source is selected and never dropped again."""


class EmittedCountError(CoverageRepairInvariantError):
    """A ledger's emitted circle count is negative, or changes although no round ran."""


class InvalidRepairTimingError(CoverageRepairInvariantError):
    """A ledger's ``repair_s`` or ``certificate_s`` is negative, infinite or NaN."""


class CertifiedSelectionRangeError(CoverageRepairInvariantError):
    """A certified path selects a source index outside its ``sources``."""


class CertifiedLedgerMismatchError(CoverageRepairInvariantError):
    """A certified path's ledger disagrees with its emission.

    ``emitted_after`` is not the emitted circle count, or an inserted source is
    not selected.
    """


class UncertifiedPathError(CoverageRepairInvariantError):
    """A certified path carries a certificate with exposed edges."""
