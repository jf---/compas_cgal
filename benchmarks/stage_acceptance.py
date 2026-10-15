"""The acceptance rows each product-pipeline stage must hold before it exits.

A stage exits when `pixi run red-manifest` passes, `pixi run gate-ratchet`
passes, and its exit-gate test (`tests/benchmarks/test_stage_exit.py`) passes.
That test reads the rows below: a `VectorRow` is a predicate on the pinned gate
vectors (`benchmarks/results/gate_vectors.json`) read by criterion name and
bound whatever their rules version, a `TestRow` names an
acceptance test that must exist and must not be a declared red, so the
red-manifest run is what proves it green.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal
from typing import Union

# Matches every case key of the gate vectors.
ALL_CASES = "*"

Relation = Literal["<=", "==", "<"]


@dataclass(frozen=True)
class VectorRow:
    """A predicate on the pinned gate vectors.

    Attributes:
        case: A case key, or `ALL_CASES`.
        criterion: The criterion's name in the vectors.
        relation: How the cell compares with *bound*.
        bound: A decimal string.
        why: The acceptance clause the row stands for.
    """

    case: str
    criterion: str
    relation: Relation
    bound: str
    why: str


@dataclass(frozen=True)
class TestRow:
    """An acceptance test that must be collected and green.

    Attributes:
        nodeid: ``path::function`` of the test.
        why: The acceptance clause the row stands for.
    """

    __test__ = False  # a data record, not a pytest class

    nodeid: str
    why: str


AcceptanceRow = Union[VectorRow, TestRow]


class UnknownAcceptanceCriterionError(Exception):
    """A `VectorRow` names a criterion the current gate vectors do not carry."""


# Rows bind to the current vectors by criterion name and bound, with no rules
# version pin (Ruling 0.3-exit-unpinned): a rules bump that changes a row's
# criterion must edit the row in the same commit.
STAGE_ACCEPTANCE: dict[str, tuple[AcceptanceRow, ...]] = {
    "stage-0": (
        TestRow("tests/benchmarks/test_survey.py::test_a_straight_slot_into_solid_stock_reads_a_half_turn", "0.1: a slot reads a half turn in process"),
        TestRow("tests/benchmarks/test_survey.py::test_the_slot_half_turn_is_an_exact_cap_verdict", "0.1: the half turn is decided exactly"),
        TestRow("tests/test_exact_sweep_idempotence.py::test_a_loop_re_swept_from_any_start_changes_nothing", "0.2: exact idempotent loop re-sweeps"),
        TestRow("tests/test_exact_sweep_idempotence.py::test_a_plunge_at_a_loops_own_entry_changes_nothing", "0.2: a plunge at a loop's own entry changes nothing"),
        TestRow(
            "tests/benchmarks/test_survey_exact_sweeps.py::test_an_incident_loop_through_its_own_cut_annulus_engages_nothing_at_any_station",
            "0.2: incident loops read no engagement",
        ),
        TestRow("tests/benchmarks/test_gate_rules.py::test_a_link_at_depth_through_cleared_stock_is_a_transit_and_not_redundant", "0.3 rule (a)"),
        TestRow("tests/benchmarks/test_gate_rules.py::test_a_second_plunge_into_virgin_stock_is_counted_and_gated", "0.3 rule (b)"),
        TestRow("tests/benchmarks/test_gate_rules.py::test_a_straight_cut_under_the_cap_is_not_slotting", "0.3 rule (c)"),
        TestRow("tests/benchmarks/test_gate_rules.py::test_a_degenerate_loop_under_the_cap_is_reported_but_not_a_violation", "0.3 rule (d)"),
        TestRow("tests/benchmarks/test_gate_rules.py::test_a_plunge_that_lands_away_from_the_next_cut_is_misplaced", "0.3 rule (e)"),
        TestRow("tests/benchmarks/test_gate_rules.py::test_the_engagement_step_is_read_at_the_junction_and_not_from_peaks", "0.3 rule (f)"),
        TestRow("tests/tools/test_gate_ratchet.py::test_a_synthetic_regression_fails_the_ratchet", "0.3: the ratchet fails a synthetic regression"),
        TestRow("tests/tools/test_gate_table.py::test_the_committed_table_is_rendered_from_the_committed_vectors", "0.3: docs numbers are read from the vectors"),
        VectorRow(ALL_CASES, "misplaced_plunges", "==", "0", "0.3 rule (e): every frozen-generator plunge lands on its next motion's start"),
        VectorRow(ALL_CASES, "unsafe_rapids", "==", "0", "Stage 0 keeps every rapid off the cutting plane"),
        VectorRow(ALL_CASES, "zero_length_motions", "==", "0", "Stage 0 measures no zero-length motion on the gate cases"),
    ),
}
