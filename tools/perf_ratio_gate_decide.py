"""The ratio gate's one decision: a workload's anchored ratio ``A_w * rho_w`` against ``1 +- m_w``.

The gate (`tools.perf_ratio_gate_verdict.gate`) decides each session workload with
`decide`; the anchor chain (`tools.perf_ratio_gate_chain.link_ratios`)
re-decides every promoted link with the same function, so a link replays only if
the gate would have promoted it. Both read ``rho_w`` through `exact_ratio`. This
module sits below both: it depends on the record only, and imports neither.

Decisions are exact: ``A_w`` is a `Fraction`, and ``rho_w`` and ``m_w`` enter as the
rationals their recorded doubles are. No float product or float bound decides.
"""

from __future__ import annotations

from enum import Enum
from fractions import Fraction
from typing import NewType

from tools.perf_ratio_record import WorkloadRatio

# Cap on a workload's pairs: 4x the calibrated N, three extensions; bounds a red workload's extra wall
# at 15 pairs (omega monstera, the costliest, ~27 s/pair on M1 Max -> ~7 min).
MAX_PAIRS = 20

# A dimensionless allowance on a workload's ratio: the gate passes ``A_w * rho_w <= 1 + m_w``.
Margin = NewType("Margin", float)
# A workload's anchor ratio A_w, exact: the reference's CPU time over the anchor commit's, a product of recorded doubles.
Anchor = NewType("Anchor", Fraction)
# A workload's ``rho_w`` as the exact rational its recorded double is (`exact_ratio`).
ExactRatio = NewType("ExactRatio", Fraction)
# The anchor of a workload at its anchor commit itself: no drift carried.
UNIT_ANCHOR = Anchor(Fraction(1))


class Verdict(str, Enum):
    """One workload's gate verdict."""

    PASS = "PASS"
    IMPROVED = "IMPROVED"
    REGRESSION = "REGRESSION"


def exact_ratio(run: WorkloadRatio) -> ExactRatio:
    """*run*'s ``rho_w``, exactly: the rational the record's ratio double is.

    The decided quotient is the record's ``ratio`` (``min head / min ref`` rounded once to a
    double, as `tools.perf_ratio_record` stores and re-derives it), read exactly; it is not
    ``Fraction(min head) / Fraction(min ref)``. Every decision and every anchor reads ``rho_w``
    here, so the gate and the chain replay decide on one quantity.
    """
    return ExactRatio(Fraction(run.ratio))


def anchored(run: WorkloadRatio, anchor: Anchor) -> Fraction:
    """The exact ``A_w * rho_w`` of *run* (`exact_ratio`)."""
    return anchor * exact_ratio(run)


def decide(run: WorkloadRatio, margin: Margin, anchor: Anchor) -> Verdict:
    """The verdict on *run*'s anchored ratio ``A_w * rho_w`` against *margin*, exactly; both bounds pass.

    Args:
        run: The workload's pairs.
        margin: The workload's calibrated ``m_w``; the bounds are the exact rationals ``1 +- m_w``.
        anchor: The workload's exact anchor ``A_w`` (`UNIT_ANCHOR` on an empty chain).

    Returns:
        REGRESSION when ``A_w * rho_w > 1 + m_w``, IMPROVED when ``A_w * rho_w < 1 - m_w``, else PASS.
    """
    product = anchored(run, anchor)
    exact_margin = Fraction(margin)
    if product > 1 + exact_margin:
        return Verdict.REGRESSION
    if product < 1 - exact_margin:
        return Verdict.IMPROVED
    return Verdict.PASS
