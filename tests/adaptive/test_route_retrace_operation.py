import hashlib
import math
from dataclasses import replace
from inspect import signature
from typing import get_args

import pytest

from compas_cgal.adaptive.canonical import canonical_cut_z_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import InvalidRetraceSegmentOperationError
from compas_cgal.adaptive.errors import InvalidRouteRetraceDecisionError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.identity import OPERATION_SCHEMA_VERSION
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
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
from compas_cgal.adaptive.units import CutZ
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.task13f_fixture import task13f_retrace_decision
from tests.adaptive.task13f_fixture import task13f_route_one_terminal


class _ForeignTerminalTraversalOwner:
    makes_cursor_terminal = True

    def __init__(self) -> None:
        self.canonical_read_count = 0

    @property
    def canonical_bytes(self) -> bytes:
        self.canonical_read_count += 1
        return b"foreign-terminal-traversal-owner"


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
    expected = encode_tagged_union(
        b"route-retrace-decision-v1",
        encode_component_map(
            {
                b"activated-entry-node-id": encode_bytes(
                    bytes(decision.activated_entry_node_id),
                ),
                b"activated-route-index": encode_integer(
                    decision.activated_route_index,
                ),
                b"activated-traversal-digest": encode_bytes(
                    bytes(decision.activated_traversal_digest),
                ),
                b"completed-exit-node-id": encode_bytes(
                    bytes(decision.completed_exit_node_id),
                ),
                b"completed-route-index": encode_integer(
                    decision.completed_route_index,
                ),
                b"source-commit-digest": encode_bytes(
                    bytes(decision.source_commit_digest),
                ),
                b"source-operation-digest": encode_bytes(
                    bytes(decision.source_operation_digest),
                ),
                b"source-operation-index": encode_integer(
                    decision.source_operation_index,
                ),
                b"source-transaction-digest": encode_bytes(
                    bytes(decision.source_transaction_digest),
                ),
                b"strategy-identity": encode_bytes(decision.strategy_identity),
                b"terminal-traversal-digest": encode_bytes(
                    bytes(decision.terminal_traversal_digest),
                ),
            }
        ),
    )

    assert decision.activated_route_index == decision.completed_route_index + 1
    assert decision.completed_exit_node_id != decision.activated_entry_node_id
    assert rebuilt.canonical_bytes == decision.canonical_bytes == expected
    assert require_canonical_record(decision.canonical_bytes) == (decision.canonical_bytes)
    assert decision.digest == IdentityDigest(
        hashlib.sha256(expected).digest(),
    )


def test_route_retrace_decision_mutation_matrix_changes_bytes_and_digest(
    decision_fields: dict[str, object],
) -> None:
    """Make every independently variable causal preimage identity-bearing."""
    changes = (
        {"completed_route_index": 2, "activated_route_index": 3},
        {"completed_exit_node_id": RouteNodeId(b"completed-exit-other")},
        {"activated_entry_node_id": RouteNodeId(b"activated-entry-other")},
        {
            "terminal_traversal_digest": IdentityDigest(
                hashlib.sha256(b"terminal-other").digest(),
            )
        },
        {
            "activated_traversal_digest": IdentityDigest(
                hashlib.sha256(b"activated-other").digest(),
            )
        },
        {
            "source_commit_digest": IdentityDigest(
                hashlib.sha256(b"source-commit-other").digest(),
            )
        },
        {
            "source_transaction_digest": IdentityDigest(
                hashlib.sha256(b"source-transaction-other").digest(),
            )
        },
        {"source_operation_index": 6},
        {
            "source_operation_digest": IdentityDigest(
                hashlib.sha256(b"source-operation-other").digest(),
            )
        },
    )
    variants = [RouteRetraceDecision.build(**decision_fields)]
    for change in changes:
        variants.append(
            RouteRetraceDecision.build(**(decision_fields | change)),
        )

    assert len({decision.canonical_bytes for decision in variants}) == len(variants)
    assert len({decision.digest for decision in variants}) == len(variants)


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


def _unchecked_foreign_traversal_source(
    source: AdvanceSegmentOperation,
) -> tuple[AdvanceSegmentOperation, _ForeignTerminalTraversalOwner]:
    """Forge an exact outer shell with a foreign terminal authority."""
    operation = object.__new__(AdvanceSegmentOperation)
    traversal_owner = _ForeignTerminalTraversalOwner()
    object.__setattr__(operation, "motion", source.motion)
    object.__setattr__(operation, "cut_z", source.cut_z)
    object.__setattr__(operation, "neck_scope", source.neck_scope)
    object.__setattr__(
        operation,
        "effective_cap_decision",
        source.effective_cap_decision,
    )
    object.__setattr__(
        operation,
        "traversal_decision",
        traversal_owner,
    )
    return operation, traversal_owner


def _unchecked_hollow_full_cap_source(
    source: AdvanceSegmentOperation,
) -> AdvanceSegmentOperation:
    """Forge an exact full-cap shell with no required nested state."""
    operation = object.__new__(AdvanceSegmentOperation)
    object.__setattr__(operation, "motion", source.motion)
    object.__setattr__(operation, "cut_z", source.cut_z)
    object.__setattr__(operation, "neck_scope", source.neck_scope)
    object.__setattr__(
        operation,
        "effective_cap_decision",
        object.__new__(FullCapDecision),
    )
    object.__setattr__(operation, "traversal_decision", source.traversal_decision)
    return operation


def _unchecked_retrace_operation(
    *,
    motion: ExactSegmentMotion,
    cut_z: CutZ,
    effective_cap_decision: FullCapDecision,
    decision: RouteRetraceDecision,
) -> RetraceSegmentOperation:
    """Assemble one closed-type shell for canonical mutation testing only."""
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


@pytest.mark.parametrize(
    ("args", "kwargs"),
    (
        pytest.param((), {}, id="zero-argument"),
        pytest.param((object(),), {}, id="positional"),
        pytest.param((), {"motion": object()}, id="keyword"),
    ),
)
def test_retrace_operation_rejects_raw_construction(
    args: tuple[object, ...],
    kwargs: dict[str, object],
) -> None:
    """Keep every retrace instance behind source-derived admission."""
    with pytest.raises(
        InvalidRetraceSegmentOperationError,
        match="must be built from an admitted source",
    ):
        RetraceSegmentOperation(*args, **kwargs)


def test_retrace_operation_rejects_foreign_traversal_owner(
    established_prefix: tuple[
        GenerationState,
        MatTraversalState,
        tuple[TraversalCommit, TraversalCommit],
    ],
) -> None:
    """Reject a forged terminal authority before authenticating its bytes."""
    physical, terminal, commits = established_prefix
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    forged_source, traversal_owner = _unchecked_foreign_traversal_source(source)
    forged_source_bytes = forged_source.canonical_bytes
    assert traversal_owner.canonical_read_count == 1
    decision = replace(
        task13f_retrace_decision(
            physical=physical,
            terminal=terminal,
            source_commit=commits[-1],
        ),
        source_operation_digest=IdentityDigest(
            hashlib.sha256(forged_source_bytes).digest(),
        ),
    )

    with pytest.raises(InvalidRetraceSegmentOperationError):
        RetraceSegmentOperation.build(
            source_operation=forged_source,
            decision=decision,
        )
    assert traversal_owner.canonical_read_count == 1


def test_retrace_operation_rejects_hollow_full_cap_source(
    established_prefix: tuple[
        GenerationState,
        MatTraversalState,
        tuple[TraversalCommit, TraversalCommit],
    ],
) -> None:
    """Translate an exact hollow nested record through the retrace error model."""
    physical, terminal, commits = established_prefix
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=commits[-1],
    )

    with pytest.raises(
        InvalidRetraceSegmentOperationError,
        match="valid exact advancing segment",
    ):
        RetraceSegmentOperation.build(
            source_operation=_unchecked_hollow_full_cap_source(source),
            decision=decision,
        )


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
    expected = encode_tagged_union(
        OPERATION_SCHEMA_VERSION,
        encode_tagged_union(
            b"retrace-segment-v1",
            encode_component_map(
                {
                    b"cut-z": canonical_cut_z_bytes(operation.cut_z),
                    b"decision": operation.decision.canonical_bytes,
                    b"effective-cap-decision": operation.effective_cap_decision.canonical_bytes,
                    b"motion": canonical_task1_bytes(operation.motion),
                }
            ),
        ),
    )

    assert operation.motion.start == source.motion.end
    assert operation.motion.end == source.motion.start
    assert operation.cut_z == source.cut_z
    assert operation.effective_cap_decision == source.effective_cap_decision
    assert not hasattr(operation, "traversal_decision")
    assert not hasattr(operation, "neck_scope")
    assert rebuilt.canonical_bytes == operation.canonical_bytes == expected
    assert require_canonical_record(operation.canonical_bytes) == (operation.canonical_bytes)
    assert b"retrace-segment-v1" in operation.canonical_bytes
    assert RetraceSegmentOperation in get_args(CanonicalOperation)
    assert tuple(signature(RetraceSegmentOperation.build).parameters) == (
        "source_operation",
        "decision",
    )


def test_retrace_operation_mutation_matrix_binds_every_component(
    decision_fields: dict[str, object],
    established_prefix: tuple[
        GenerationState,
        MatTraversalState,
        tuple[TraversalCommit, TraversalCommit],
    ],
) -> None:
    """Bind motion, depth, cap, and causal decision independently."""
    physical, terminal, commits = established_prefix
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    operation = RetraceSegmentOperation.build(
        source_operation=source,
        decision=task13f_retrace_decision(
            physical=physical,
            terminal=terminal,
            source_commit=commits[-1],
        ),
    )
    changed_cap_value = EngagementCap.build(math.pi / 3.0)
    changed_cap = FullCapDecision.build(
        user_cap=changed_cap_value,
        effective_cap=changed_cap_value,
    )
    changed_decision = RouteRetraceDecision.build(**decision_fields)
    variants = (
        operation,
        _unchecked_retrace_operation(
            motion=source.motion,
            cut_z=operation.cut_z,
            effective_cap_decision=operation.effective_cap_decision,
            decision=operation.decision,
        ),
        _unchecked_retrace_operation(
            motion=operation.motion,
            cut_z=CutZ.build(operation.cut_z.value - 1.0),
            effective_cap_decision=operation.effective_cap_decision,
            decision=operation.decision,
        ),
        _unchecked_retrace_operation(
            motion=operation.motion,
            cut_z=operation.cut_z,
            effective_cap_decision=changed_cap,
            decision=operation.decision,
        ),
        _unchecked_retrace_operation(
            motion=operation.motion,
            cut_z=operation.cut_z,
            effective_cap_decision=operation.effective_cap_decision,
            decision=changed_decision,
        ),
    )

    assert len({variant.canonical_bytes for variant in variants}) == len(variants)


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
