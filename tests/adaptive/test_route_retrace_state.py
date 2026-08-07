"""Authoritative physical-chronology contracts for route retrace."""

import hashlib
import math
from dataclasses import replace

import pytest

from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.errors import CandidateStateMismatchError
from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.retrace_segment_trial import evaluate_retrace_segment_trial
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.task13f_fixture import task13f_retrace_decision
from tests.adaptive.task13f_fixture import task13f_route_one_terminal
from tests.adaptive.test_transaction import _state_fixture


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build one authenticated route-one retrace boundary."""
    return Task13FFixture.build()


def _evaluated_retrace(
    task13f: Task13FFixture,
) -> tuple[
    GenerationState,
    AdvanceSegmentOperation,
    RetraceSegmentOperation,
    Stock2Area,
    CoverageLedger,
]:
    """Evaluate the exact route-one reverse on independent physical owners."""
    parent, terminal, commits = task13f_route_one_terminal(task13f)
    source = parent.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    operation = RetraceSegmentOperation.build(
        source_operation=source,
        decision=task13f_retrace_decision(
            physical=parent,
            terminal=terminal,
            source_commit=commits[-1],
        ),
    )
    stock = parent.fork_stock()
    coverage = parent.clone_coverage()
    evaluate_retrace_segment_trial(
        containment_authority=GougeContainment.build(
            task13f.identity.reachable_domain,
        ),
        stock=stock,
        coverage=coverage,
        operation_index=len(parent.operations),
        operation=operation,
        tool_radius=task13f.identity.tool_radius,
        user_cap=task13f.identity.user_cap,
        effective_cap=task13f.identity.user_cap,
        depletion_policy=task13f.identity.depletion_policy,
    )
    return parent, source, operation, stock, coverage


@pytest.fixture(scope="module")
def evaluated_retrace(
    task13f: Task13FFixture,
) -> tuple[
    GenerationState,
    AdvanceSegmentOperation,
    RetraceSegmentOperation,
    Stock2Area,
    CoverageLedger,
]:
    """Share one expensive exact reverse while every test forks its owners."""
    return _evaluated_retrace(task13f)


def _unchecked_retrace(
    *,
    motion: ExactSegmentMotion,
    cut_z: CutZ,
    effective_cap_decision: FullCapDecision,
    decision: RouteRetraceDecision,
) -> RetraceSegmentOperation:
    """Forge one exact outer shell for state-boundary corruption tests."""
    operation = object.__new__(RetraceSegmentOperation)
    object.__setattr__(operation, "motion", motion)
    object.__setattr__(operation, "cut_z", cut_z)
    object.__setattr__(
        operation,
        "effective_cap_decision",
        effective_cap_decision,
    )
    object.__setattr__(operation, "decision", decision)
    return operation


def _build_child(
    *,
    parent: GenerationState,
    operation: RetraceSegmentOperation,
    stock: Stock2Area,
    coverage: CoverageLedger,
    phase_point: Point2[WorldXY] | None = None,
    operations: tuple[CanonicalOperation, ...] | None = None,
    passages: tuple[NeckPassage, ...] | None = None,
) -> GenerationState:
    """Submit one retrace child at the authoritative state boundary."""
    return GenerationState.build(
        stock=stock,
        coverage=coverage,
        tool_radius=parent.tool_radius,
        phase_point=operation.motion.end if phase_point is None else phase_point,
        traversal=parent.traversal,
        passages=parent.passages if passages is None else passages,
        operations=parent.operations + (operation,) if operations is None else operations,
    )


class _RetraceSubclass(RetraceSegmentOperation):
    """Foreign subtype used to prove exact evaluator variant closure."""


def test_retrace_state_changes_only_physical_chronology(
    task13f: Task13FFixture,
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Return through removed stock without manufacturing MAT progress.

    Args:
        task13f: Exact policy authority for candidate-parent validation.
        evaluated_retrace: Real route-one parent and mutated trial owners.
    """
    parent, source, operation, stock, coverage = evaluated_retrace
    child = _build_child(
        parent=parent,
        operation=operation,
        stock=stock.fork(),
        coverage=coverage.clone(),
    )

    assert child.phase_point == source.motion.start
    assert child.traversal.canonical_bytes == parent.traversal.canonical_bytes
    assert child.passages == parent.passages
    assert len(child.fork_stock().lineage) == len(parent.fork_stock().lineage) + 1
    assert len(child.clone_coverage().lineage) == len(parent.clone_coverage().lineage) + 1
    assert child.operations[-1] == operation
    task13f.evaluator._validate_physical_parent(child)


def test_generation_state_rejects_hollow_exact_retrace_with_named_error(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Close the state boundary before common lateral fields are read."""
    parent, _, operation, stock, coverage = evaluated_retrace
    hollow = object.__new__(RetraceSegmentOperation)

    with pytest.raises(InvalidGenerationStateError):
        GenerationState.build(
            stock=stock.fork(),
            coverage=coverage.clone(),
            tool_radius=parent.tool_radius,
            phase_point=operation.motion.end,
            traversal=parent.traversal,
            passages=parent.passages,
            operations=parent.operations + (hollow,),
        )


def test_candidate_evaluator_rejects_hollow_exact_retrace_with_named_error(
    task13f: Task13FFixture,
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Close evaluator cap validation before a hollow retrace is read."""
    parent, _, _, _, _ = evaluated_retrace
    hollow = object.__new__(RetraceSegmentOperation)
    state = object.__new__(GenerationState)
    object.__setattr__(
        state,
        "operations",
        parent.operations + (hollow,),
    )

    with pytest.raises(CandidateStateMismatchError):
        task13f.evaluator._validate_parent_operations(state)


@pytest.mark.parametrize(
    "foreign_operation",
    (
        pytest.param(object(), id="foreign"),
        pytest.param(object.__new__(_RetraceSubclass), id="subclass"),
    ),
)
def test_candidate_evaluator_rejects_nonexact_parent_operation_variant(
    task13f: Task13FFixture,
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
    foreign_operation: object,
) -> None:
    """Keep foreign and subtype operations outside the exact parent union."""
    parent, _, _, _, _ = evaluated_retrace
    state = object.__new__(GenerationState)
    object.__setattr__(
        state,
        "operations",
        parent.operations + (foreign_operation,),  # type: ignore[arg-type]
    )

    with pytest.raises(CandidateStateMismatchError):
        task13f.evaluator._validate_parent_operations(state)


def test_retrace_state_rejects_non_advance_predecessor(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a retrace naming an immediately preceding retrace as source."""
    parent, _, operation, stock, coverage = evaluated_retrace
    second = _unchecked_retrace(
        motion=operation.motion,
        cut_z=operation.cut_z,
        effective_cap_decision=operation.effective_cap_decision,
        decision=replace(
            operation.decision,
            source_operation_index=len(parent.operations),
            source_operation_digest=IdentityDigest(
                hashlib.sha256(operation.canonical_bytes).digest(),
            ),
        ),
    )
    second_stock = stock.fork()
    second_stock.deplete(second.motion, parent.tool_radius, parent.fork_stock().lineage[-1].policy)
    second_coverage = coverage.clone()
    second_coverage.add_sweep(second.motion, parent.tool_radius)

    with pytest.raises(InvalidGenerationStateError, match="source is not an advancing segment"):
        _build_child(
            parent=parent,
            operation=second,
            stock=second_stock,
            coverage=second_coverage,
            operations=parent.operations + (operation, second),
        )


def test_retrace_state_rejects_non_immediate_source_ordinal(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a source ordinal that does not name the immediate predecessor."""
    parent, _, operation, stock, coverage = evaluated_retrace
    forged = _unchecked_retrace(
        motion=operation.motion,
        cut_z=operation.cut_z,
        effective_cap_decision=operation.effective_cap_decision,
        decision=replace(
            operation.decision,
            source_operation_index=operation.decision.source_operation_index - 1,
        ),
    )

    with pytest.raises(InvalidGenerationStateError, match="immediately follow"):
        _build_child(
            parent=parent,
            operation=forged,
            stock=stock.fork(),
            coverage=coverage.clone(),
        )


def test_retrace_state_rejects_source_digest_corruption(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a retrace whose decision does not hash its stored source."""
    parent, _, operation, stock, coverage = evaluated_retrace
    forged = _unchecked_retrace(
        motion=operation.motion,
        cut_z=operation.cut_z,
        effective_cap_decision=operation.effective_cap_decision,
        decision=replace(
            operation.decision,
            source_operation_digest=IdentityDigest(
                hashlib.sha256(b"wrong-state-source").digest(),
            ),
        ),
    )

    with pytest.raises(InvalidGenerationStateError, match="contradicts its source or physical phase"):
        _build_child(
            parent=parent,
            operation=forged,
            stock=stock.fork(),
            coverage=coverage.clone(),
        )


def test_retrace_state_rejects_forged_non_reverse_motion(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a stored retrace shell that repeats the source direction."""
    parent, source, operation, _, _ = evaluated_retrace
    forged = _unchecked_retrace(
        motion=source.motion,
        cut_z=operation.cut_z,
        effective_cap_decision=operation.effective_cap_decision,
        decision=operation.decision,
    )
    stock = parent.fork_stock()
    stock.deplete(forged.motion, parent.tool_radius, parent.fork_stock().lineage[-1].policy)
    coverage = parent.clone_coverage()
    coverage.add_sweep(forged.motion, parent.tool_radius)

    with pytest.raises(InvalidGenerationStateError, match="contradicts its source or physical phase"):
        _build_child(
            parent=parent,
            operation=forged,
            stock=stock,
            coverage=coverage,
        )


def test_retrace_state_rejects_changed_cut_depth(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a retrace shell moved off its qualified cut plane."""
    parent, _, operation, stock, coverage = evaluated_retrace
    forged = _unchecked_retrace(
        motion=operation.motion,
        cut_z=CutZ.build(-3.0),
        effective_cap_decision=operation.effective_cap_decision,
        decision=operation.decision,
    )

    with pytest.raises(InvalidGenerationStateError, match="cut depth"):
        _build_child(
            parent=parent,
            operation=forged,
            stock=stock.fork(),
            coverage=coverage.clone(),
        )


def test_retrace_state_rejects_changed_full_cap(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a retrace shell that changes its source full-cap authority."""
    parent, _, operation, stock, coverage = evaluated_retrace
    foreign_cap = EngagementCap.build(math.pi / 2.0)
    forged = _unchecked_retrace(
        motion=operation.motion,
        cut_z=operation.cut_z,
        effective_cap_decision=FullCapDecision.build(
            user_cap=foreign_cap,
            effective_cap=foreign_cap,
        ),
        decision=operation.decision,
    )

    with pytest.raises(InvalidGenerationStateError, match="contradicts its source or physical phase"):
        _build_child(
            parent=parent,
            operation=forged,
            stock=stock.fork(),
            coverage=coverage.clone(),
        )


def test_retrace_state_rejects_changed_start_phase(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a reverse ending correctly but starting off the physical phase."""
    parent, source, operation, _, _ = evaluated_retrace
    changed_start = Point2[WorldXY].build(
        (source.motion.start.x + source.motion.end.x) / 2.0,
        source.motion.start.y,
    )
    forged = _unchecked_retrace(
        motion=ExactSegmentMotion.build(changed_start, source.motion.start),
        cut_z=operation.cut_z,
        effective_cap_decision=operation.effective_cap_decision,
        decision=operation.decision,
    )
    stock = parent.fork_stock()
    stock.deplete(forged.motion, parent.tool_radius, parent.fork_stock().lineage[-1].policy)
    coverage = parent.clone_coverage()
    coverage.add_sweep(forged.motion, parent.tool_radius)

    with pytest.raises(InvalidGenerationStateError, match="contradicts its source or physical phase"):
        _build_child(
            parent=parent,
            operation=forged,
            stock=stock,
            coverage=coverage,
        )


def test_retrace_state_rejects_pending_hold_link(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject retrace as a consumer of an unpaired hold link."""
    parent, _, operation, _, _ = evaluated_retrace
    link = LinkSegmentOperation.build(
        motion=operation.motion,
        cut_z=operation.cut_z,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=operation.effective_cap_decision,
        traversal_decision=HoldTraversalDecision.build(
            component_id=parent.traversal.component_id,
            edge_id=parent.traversal.edge_id,
            branch_id=parent.traversal.branch_id,
            cursor=parent.traversal.cursor,
        ),
    )
    stock = parent.fork_stock()
    coverage = parent.clone_coverage()
    for motion in (link.motion, operation.motion):
        stock.deplete(motion, parent.tool_radius, parent.fork_stock().lineage[-1].policy)
        coverage.add_sweep(motion, parent.tool_radius)

    with pytest.raises(InvalidGenerationStateError, match="pending hold link"):
        _build_child(
            parent=parent,
            operation=operation,
            stock=stock,
            coverage=coverage,
            operations=parent.operations + (link, operation),
        )


def test_retrace_state_rejects_changed_local_cursor(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Keep retrace from manufacturing a new local MAT cursor."""
    parent, _, operation, stock, coverage = evaluated_retrace
    changed_traversal = replace(
        parent.traversal,
        cursor=CursorIdentity(b"forged-retrace-cursor"),
    )

    with pytest.raises(InvalidGenerationStateError, match="traversal"):
        GenerationState.build(
            stock=stock.fork(),
            coverage=coverage.clone(),
            tool_radius=parent.tool_radius,
            phase_point=operation.motion.end,
            traversal=changed_traversal,
            passages=parent.passages,
            operations=parent.operations + (operation,),
        )


def test_retrace_state_rejects_changed_passage(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Keep retrace from manufacturing an unrelated neck transition."""
    parent, _, operation, stock, coverage = evaluated_retrace
    passage = _state_fixture().state.passages[0]
    changed_passage = replace(
        passage,
        state=PassageState.FIRST_PASS_COMPLETE,
    )

    with pytest.raises(InvalidGenerationStateError, match="passage state"):
        _build_child(
            parent=parent,
            operation=operation,
            stock=stock.fork(),
            coverage=coverage.clone(),
            passages=(changed_passage,),
        )


def test_retrace_state_rejects_missing_stock_lineage(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a retrace operation without its exact depletion event."""
    parent, _, operation, _, coverage = evaluated_retrace

    with pytest.raises(InvalidGenerationStateError, match="different lateral lengths"):
        _build_child(
            parent=parent,
            operation=operation,
            stock=parent.fork_stock(),
            coverage=coverage.clone(),
        )


def test_retrace_state_rejects_missing_coverage_lineage(
    evaluated_retrace: tuple[
        GenerationState,
        AdvanceSegmentOperation,
        RetraceSegmentOperation,
        Stock2Area,
        CoverageLedger,
    ],
) -> None:
    """Reject a retrace operation without its exact coverage event."""
    parent, _, operation, stock, _ = evaluated_retrace

    with pytest.raises(InvalidGenerationStateError, match="different lateral lengths"):
        _build_child(
            parent=parent,
            operation=operation,
            stock=stock.fork(),
            coverage=parent.clone_coverage(),
        )
