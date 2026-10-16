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
"""

from __future__ import annotations

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
