"""Wall-clock witnesses for the native curved-circle query on the prepared Held pockets.

The slowest pieces measured on 2026-09-08 (Monstera 159/85/160 at 2.6-3.1 s,
upper 75 at 0.4 s) were not near-degenerate: their segment-interior tests
failed CORE's floating filter on generic values because the ray parameter
was formed by a division, and each failure cost one root-bound evaluation.
"""

import time

import pytest

from benchmarks.held_native_curve_import import import_held_boundary
from benchmarks.held_reference_cases import load_held_reference_case

# One order of magnitude above the repaired regime (tens of milliseconds for
# 300 competitors) and one below the measured failure (2.6-3.1 s).
SLOW_PIECE_BUDGET_S = 0.25


@pytest.mark.parametrize(
    "case_name,piece",
    [("figure8_monstera", 159), ("figure8_monstera", 85), ("figure8_monstera", 160), ("figure8_upper", 75)],
)
def test_slowest_measured_pieces_stay_within_budget(case_name: str, piece: int) -> None:
    case = load_held_reference_case(case_name)
    owner = import_held_boundary(case)
    started = time.perf_counter()
    proposal = owner.circle_on_piece(piece, 0.5, float(case.tool_radius.value))
    seconds = time.perf_counter() - started
    assert proposal.clearance_mm >= float(case.tool_radius.value)
    assert seconds < SLOW_PIECE_BUDGET_S
