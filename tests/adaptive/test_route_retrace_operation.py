import hashlib
from dataclasses import replace
from inspect import signature
from typing import get_args

import pytest

from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import InvalidRetraceSegmentOperationError
from compas_cgal.adaptive.errors import InvalidRouteRetraceDecisionError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckOwnerId
from compas_cgal.adaptive.operation import NeckTraversalOrientation
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.operation import RouteNodeId
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.operation import WidthClassId
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.traversal import MatTraversalState
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.task13f_fixture import task13f_retrace_decision
from tests.adaptive.task13f_fixture import task13f_route_one_terminal


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build one immutable production fixture for all retrace attacks."""
    return Task13FFixture.build()


@pytest.fixture(scope="module")
def established_prefix(
    task13f: Task13FFixture,
) -> tuple[
    GenerationState,
    MatTraversalState,
    tuple[TraversalCommit, TraversalCommit],
]:
    """Reuse the accepted Task 13F prefix at the route-one terminal."""
    return task13f_route_one_terminal(task13f)


@pytest.fixture
def decision_fields() -> dict[str, object]:
    """Provide distinct preimages for one authenticated route boundary."""
    return {
        "completed_route_index": 1,
        "activated_route_index": 2,
        "completed_exit_node_id": RouteNodeId(b"completed-exit"),
        "activated_entry_node_id": RouteNodeId(b"activated-entry"),
        "terminal_traversal_digest": IdentityDigest(
            hashlib.sha256(b"terminal").digest(),
        ),
        "activated_traversal_digest": IdentityDigest(
            hashlib.sha256(b"activated").digest(),
        ),
        "source_commit_digest": IdentityDigest(
            hashlib.sha256(b"source-commit").digest(),
        ),
        "source_transaction_digest": IdentityDigest(
            hashlib.sha256(b"source-transaction").digest(),
        ),
        "source_operation_index": 5,
        "source_operation_digest": IdentityDigest(
            hashlib.sha256(b"source-operation").digest(),
        ),
    }


def test_route_retrace_decision_is_content_addressed(
    decision_fields: dict[str, object],
) -> None:
    """Bind a nonincident route switch and its accepted source by exact identity."""
    decision = RouteRetraceDecision.build(**decision_fields)
    rebuilt = RouteRetraceDecision.build(
        **dict(reversed(tuple(decision_fields.items()))),
    )

    assert decision.activated_route_index == decision.completed_route_index + 1
    assert decision.completed_exit_node_id != decision.activated_entry_node_id
    assert rebuilt.canonical_bytes == decision.canonical_bytes
    assert require_canonical_record(decision.canonical_bytes) == (decision.canonical_bytes)
    assert decision.digest == IdentityDigest(
        hashlib.sha256(decision.canonical_bytes).digest(),
    )


@pytest.mark.parametrize(
    ("field", "replacement"),
    (
        ("completed_route_index", -1),
        ("completed_route_index", True),
        ("activated_route_index", -1),
        ("activated_route_index", 4),
        ("activated_entry_node_id", RouteNodeId(b"")),
        ("completed_exit_node_id", RouteNodeId(b"same-node")),
        ("terminal_traversal_digest", IdentityDigest(b"short")),
        ("activated_traversal_digest", IdentityDigest(b"short")),
        ("source_commit_digest", IdentityDigest(b"short")),
        ("source_transaction_digest", IdentityDigest(b"short")),
        ("source_operation_digest", IdentityDigest(b"short")),
        ("source_operation_index", 1),
        ("source_operation_index", True),
    ),
)
def test_route_retrace_decision_rejects_malformed_identity(
    decision_fields: dict[str, object],
    field: str,
    replacement: object,
) -> None:
    """Reject a decision that cannot name one adjacent authenticated boundary."""
    changed = dict(decision_fields)
    changed[field] = replacement
    if field == "completed_exit_node_id":
        changed["activated_entry_node_id"] = replacement

    with pytest.raises(InvalidRouteRetraceDecisionError):
        RouteRetraceDecision.build(**changed)


def test_route_retrace_decision_rejects_foreign_strategy(
    decision_fields: dict[str, object],
) -> None:
    """Reject a decision whose source-reversal strategy changed after sealing."""
    decision = RouteRetraceDecision.build(**decision_fields)

    with pytest.raises(InvalidRouteRetraceDecisionError):
        replace(decision, strategy_identity=b"foreign-retrace-strategy-v1")


def _task13f_neck_cap(task13f: Task13FFixture) -> NeckCapDecision:
    """Build a valid neck-owned cap that retrace admission must reject."""
    return NeckCapDecision.build(
        neck_evidence_digest=IdentityDigest(
            hashlib.sha256(b"foreign-neck-evidence").digest(),
        ),
        width_class_id=WidthClassId.build(0),
        passage_before=PassageState.UNVISITED,
        passage_after=PassageState.FIRST_PASS_COMPLETE,
        user_cap=task13f.identity.user_cap,
        effective_cap=task13f.identity.user_cap,
    )


def _unchecked_non_full_cap_source(
    source: AdvanceSegmentOperation,
    cap: NeckCapDecision,
) -> AdvanceSegmentOperation:
    """Forge one internally cross-wired source to exercise fail-closed admission."""
    operation = object.__new__(AdvanceSegmentOperation)
    object.__setattr__(operation, "motion", source.motion)
    object.__setattr__(operation, "cut_z", source.cut_z)
    object.__setattr__(operation, "neck_scope", source.neck_scope)
    object.__setattr__(operation, "effective_cap_decision", cap)
    object.__setattr__(operation, "traversal_decision", source.traversal_decision)
    return operation


def test_retrace_operation_derives_only_the_exact_source_reverse(
    established_prefix: tuple[
        GenerationState,
        MatTraversalState,
        tuple[TraversalCommit, TraversalCommit],
    ],
) -> None:
    """Reverse the accepted horizontal leaf without accepting caller geometry."""
    physical, terminal, commits = established_prefix
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=commits[-1],
    )
    operation = RetraceSegmentOperation.build(
        source_operation=source,
        decision=decision,
    )
    rebuilt = RetraceSegmentOperation.build(
        source_operation=source,
        decision=decision,
    )

    assert operation.motion.start == source.motion.end
    assert operation.motion.end == source.motion.start
    assert operation.cut_z == source.cut_z
    assert operation.effective_cap_decision == source.effective_cap_decision
    assert not hasattr(operation, "traversal_decision")
    assert not hasattr(operation, "neck_scope")
    assert rebuilt.canonical_bytes == operation.canonical_bytes
    assert require_canonical_record(operation.canonical_bytes) == (operation.canonical_bytes)
    assert b"retrace-segment-v1" in operation.canonical_bytes
    assert RetraceSegmentOperation in get_args(CanonicalOperation)
    assert tuple(signature(RetraceSegmentOperation.build).parameters) == (
        "source_operation",
        "decision",
    )


@pytest.mark.parametrize(
    "source_variant",
    (
        "neck-owning",
        "non-full-cap",
        "nonterminal",
        "wrong-source-digest",
        "circle",
        "hold-link",
    ),
)
def test_retrace_operation_rejects_unadmitted_source(
    task13f: Task13FFixture,
    established_prefix: tuple[
        GenerationState,
        MatTraversalState,
        tuple[TraversalCommit, TraversalCommit],
    ],
    source_variant: str,
) -> None:
    """Reject every locally decidable corruption of the admitted source."""
    physical, terminal, commits = established_prefix
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=commits[-1],
    )
    replacement: object = source
    if source_variant == "neck-owning":
        replacement = AdvanceSegmentOperation.build(
            motion=source.motion,
            cut_z=source.cut_z,
            neck_scope=OrientedNeckScope.build(
                neck_owner_id=NeckOwnerId(b"foreign-neck"),
                orientation=NeckTraversalOrientation.FORWARD,
            ),
            effective_cap_decision=_task13f_neck_cap(task13f),
            traversal_decision=source.traversal_decision,
        )
    elif source_variant == "non-full-cap":
        replacement = _unchecked_non_full_cap_source(
            source,
            _task13f_neck_cap(task13f),
        )
    elif source_variant == "nonterminal":
        traversal = source.traversal_decision
        replacement = AdvanceSegmentOperation.build(
            motion=source.motion,
            cut_z=source.cut_z,
            neck_scope=source.neck_scope,
            effective_cap_decision=source.effective_cap_decision,
            traversal_decision=AdvanceTraversalDecision.build(
                component_id=traversal.component_id,
                edge_id=traversal.edge_id,
                branch_id=traversal.branch_id,
                cursor_before=traversal.cursor_before,
                cursor_after=traversal.cursor_after,
                makes_cursor_terminal=False,
            ),
        )
    elif source_variant == "wrong-source-digest":
        decision = replace(
            decision,
            source_operation_digest=IdentityDigest(
                hashlib.sha256(b"wrong-source-operation").digest(),
            ),
        )
    elif source_variant == "circle":
        replacement = next(operation for operation in reversed(physical.operations) if type(operation) is CutFullCircleOperation)
    elif source_variant == "hold-link":
        replacement = next(operation for operation in reversed(physical.operations) if type(operation) is LinkSegmentOperation)
    else:
        raise AssertionError(f"unhandled source variant: {source_variant}")

    with pytest.raises(InvalidRetraceSegmentOperationError):
        RetraceSegmentOperation.build(
            source_operation=replacement,  # type: ignore[arg-type]
            decision=decision,
        )


def test_retrace_operation_rejects_caller_supplied_motion(
    established_prefix: tuple[
        GenerationState,
        MatTraversalState,
        tuple[TraversalCommit, TraversalCommit],
    ],
) -> None:
    """Keep forward, shortened, offset, and degenerate geometry outside the API."""
    physical, terminal, commits = established_prefix
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=commits[-1],
    )

    with pytest.raises(TypeError, match="unexpected keyword argument 'motion'"):
        RetraceSegmentOperation.build(
            source_operation=source,
            decision=decision,
            motion=ExactSegmentMotion.build(  # type: ignore[call-arg]
                source.motion.start,
                source.motion.end,
            ),
        )
