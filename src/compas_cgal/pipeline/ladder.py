"""One finite furthest-first scan over candidates with injected, ordered conditions (Ruling 1.1-R13).

WHY A SCAN AND NEVER A BISECTION. Engagement is not monotone in the advance (nor
in the loop radius), so the admissible candidates are not an up-set and the
pass/fail pattern down a window can alternate. Walking the candidates from the
furthest and returning the first that passes every condition is correct under
any pattern; a bisection is correct only on an up-set.

Per candidate the conditions are evaluated in their given order and the first
refusal rejects that candidate; the scan then moves to the next candidate and
never stops early (Ruling 1.1-R13-scan). The order is an evaluation order, not a
relaxation: cheap structural conditions first, exact predicates after.

When nothing passes, the scan returns `NoAdmissibleCandidate` and takes no
fallback (amendment D2); a caller that must still advance (the narrow forced
advance) decides that itself and records it.

`WarmStart` is the same principle one level down: the evaluation order of a
conjunction of indexed exact probes, warm-started at the probe that refused
last. It changes which probe is asked first, never which probes exist or how
each is decided, so the conjunction's verdict is order-free (Task 1.2 lever L6).
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol
from typing import TypeVar

from compas_cgal.stock import Stock

CandidateT = TypeVar("CandidateT")
CandidateT_contra = TypeVar("CandidateT_contra", contravariant=True)


class CandidateCondition(Protocol[CandidateT_contra]):
    """One admissibility condition: does *candidate* pass, advancing from *origin* on *stock*?

    Implementations never mutate *stock*.
    """

    def __call__(self, stock: Stock, origin: CandidateT_contra, candidate: CandidateT_contra, /) -> bool: ...


@dataclass(frozen=True)
class NoAdmissibleCandidate:
    """The scan found no candidate passing every condition.

    Attributes:
        scanned: How many candidates were evaluated.
    """

    scanned: int


def furthest_first(
    stock: Stock,
    origin: CandidateT,
    candidates: Sequence[CandidateT],
    conditions: tuple[CandidateCondition[CandidateT], ...],
) -> CandidateT | NoAdmissibleCandidate:
    """The first candidate, in the given furthest-first order, that passes every condition.

    Args:
        stock: The current stock (never mutated here).
        origin: Where the advance starts.
        candidates: Candidates ordered furthest first.
        conditions: Conditions in evaluation order.

    Returns:
        The furthest admissible candidate, or `NoAdmissibleCandidate`.
    """
    for candidate in candidates:
        if all(condition(stock, origin, candidate) for condition in conditions):
            return candidate
    return NoAdmissibleCandidate(scanned=len(candidates))


@dataclass
class WarmStart:
    """Evaluation-order hint for a conjunction of ``count`` indexed probes: ask the last refuser first.

    A COST device, never a decision. "Some probe refuses" is a disjunction of
    the same per-probe verdicts, so any order of asking them gives the same
    answer; the hint only moves the index that refused last to the front and
    keeps the rest in their natural order. Every probe is asked at most once
    per evaluation, and an evaluation that finds no refusal has asked all of
    them. Refusals repeat at the same index across consecutive candidates of a
    scan (the ladder's candidates are near neighbours), so a refusal usually
    costs one probe.

    The hint is mutable state shared by the evaluations that use this
    instance; reading or updating it never changes a verdict.

    Attributes:
        hint: The index asked first while it lies in ``range(count)``.
    """

    hint: int = 0

    def first_refusal(self, count: int, refuses: Callable[[int], bool]) -> int | None:
        """The first refusing index in warm-start order, or None when no probe refuses.

        Args:
            count: Number of probes, indexed ``0 .. count - 1``.
            refuses: The exact per-probe verdict: True iff probe *index* refuses.

        Returns:
            A refusing index (and the hint is moved to it), or None iff no
            index in ``range(count)`` refuses.
        """
        hint = self.hint  # read once: an interleaved evaluation may move it while this one runs
        first: list[int] = [hint] if 0 <= hint < count else []
        order: list[int] = [*first, *(index for index in range(count) if index != hint)]
        for index in order:
            if refuses(index):
                self.hint = index
                return index
        return None
