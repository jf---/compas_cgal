"""Exact cross-axis route-retrace orchestration contracts."""

import hashlib
from dataclasses import dataclass
from dataclasses import fields
from typing import TypeVar
from typing import cast

import pytest

from compas_cgal.adaptive.errors import InvalidRouteRetraceDecisionError
from compas_cgal.adaptive.errors import InvalidRouteRetraceCommitError
from compas_cgal.adaptive.errors import UnsupportedRouteRetraceError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generator import _derive_route_retrace_decision
from compas_cgal.adaptive.generator import _route_retrace_required
from compas_cgal.adaptive.generator import GenerationContinuation
from compas_cgal.adaptive.generator import materialize_active_candidate_family
from compas_cgal.adaptive.generator import TraversalCommit
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_MOTION_STRATA
from compas_cgal.adaptive.motion_certificate import SweptPrefixMotionWitness
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.retrace_transaction import RouteRetraceEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.transaction import ZeroGuideLinkTransaction
from compas_cgal.adaptive.traversal import MatTraversalState
from tests.adaptive.task13f_fixture import TASK13F_ROUTE_ONE_COMMIT_DIGEST
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.task13f_fixture import task13f_retrace_continuation
from tests.adaptive.task13f_fixture import task13f_route_one_terminal

_T = TypeVar("_T")


def _identity(label: bytes) -> IdentityDigest:
    """Return one deterministic foreign identity for causal fault injection."""
    return IdentityDigest(hashlib.sha256(label).digest())


def _raw_copy(value: _T, **changes: object) -> _T:
    """Forge one dataclass shell without invoking the nested owner."""
    forged = object.__new__(type(value))
    for item in fields(value):
        object.__setattr__(
            forged,
            item.name,
            changes.get(item.name, getattr(value, item.name)),
        )
    return cast(_T, forged)


@dataclass(frozen=True)
class _RouteBoundary:
    """Established physical/global prefix at the first nonincident boundary."""

    physical: GenerationState
    terminal: MatTraversalState
    activated: MatTraversalState
    commits: tuple[TraversalCommit, TraversalCommit]


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build one authenticated Task 13F policy and physical authority."""
    return Task13FFixture.build()


@pytest.fixture(scope="module")
def boundary(task13f: Task13FFixture) -> _RouteBoundary:
    """Recover the unchanged accepted prefix before any retrace decision."""
    physical, terminal, commits = task13f_route_one_terminal(task13f)
    return _RouteBoundary(
        physical=physical,
        terminal=terminal,
        activated=terminal.activate_next(),
        commits=commits,
    )


def test_task13f_route_trigger_distinguishes_incident_and_nonincident_edges(
    boundary: _RouteBoundary,
) -> None:
    """Use stable MAT node identity, never coordinates, for route transport.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    route_zero_terminal = boundary.commits[0].traversal_after
    route_one_active = route_zero_terminal.activate_next()

    assert route_zero_terminal.active_cursor.route_step.exit_node_id == route_one_active.active_cursor.route_step.entry_node_id
    assert not _route_retrace_required(
        route_zero_terminal,
        route_one_active,
    )
    assert boundary.terminal.active_cursor.route_step.exit_node_id != boundary.activated.active_cursor.route_step.entry_node_id
    assert _route_retrace_required(
        boundary.terminal,
        boundary.activated,
    )


@pytest.mark.parametrize(
    ("terminal_variant", "activated_variant"),
    (
        (object(), None),
        (None, object()),
        ("nonterminal", None),
        (None, "same-state"),
        (None, "wrong-route-index"),
    ),
)
def test_route_trigger_rejects_foreign_or_noncausal_transitions(
    boundary: _RouteBoundary,
    terminal_variant: object,
    activated_variant: object,
) -> None:
    """Reject any trigger input not produced by one exact activation step.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
        terminal_variant: Fault injected into the completed-route state.
        activated_variant: Fault injected into the activated-route state.
    """
    terminal: object = boundary.terminal
    activated: object = boundary.activated
    if terminal_variant == "nonterminal":
        terminal = boundary.activated
    elif terminal_variant is not None:
        terminal = terminal_variant
    if activated_variant == "same-state":
        activated = boundary.terminal
    elif activated_variant == "wrong-route-index":
        activated = _raw_copy(
            boundary.activated,
            active_route_index=0,
        )
    elif activated_variant is not None:
        activated = activated_variant

    with pytest.raises(InvalidRouteRetraceDecisionError):
        _route_retrace_required(
            cast(MatTraversalState, terminal),
            cast(MatTraversalState, activated),
        )


def test_route_retrace_decision_binds_the_final_exact_source(
    boundary: _RouteBoundary,
) -> None:
    """Seal only the accepted final zero-guide advance as causal transport.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    source_commit = boundary.commits[-1]
    source_operation = boundary.physical.operations[-1]
    assert type(source_operation) is AdvanceSegmentOperation

    decision = _derive_route_retrace_decision(
        physical=boundary.physical,
        terminal=boundary.terminal,
        activated=boundary.activated,
        source_commit=source_commit,
    )

    assert type(decision) is RouteRetraceDecision
    assert decision.completed_route_index == 1
    assert decision.activated_route_index == 2
    assert decision.completed_exit_node_id == (boundary.terminal.active_cursor.route_step.exit_node_id)
    assert decision.activated_entry_node_id == (boundary.activated.active_cursor.route_step.entry_node_id)
    assert decision.terminal_traversal_digest == boundary.terminal.digest
    assert decision.activated_traversal_digest == boundary.activated.digest
    assert decision.source_commit_digest.hex() == (TASK13F_ROUTE_ONE_COMMIT_DIGEST)
    assert decision.source_transaction_digest == source_commit.transaction.digest
    assert decision.source_operation_index == len(boundary.physical.operations) - 1
    assert decision.source_operation_digest == IdentityDigest(
        hashlib.sha256(source_operation.canonical_bytes).digest(),
    )


def test_route_retrace_derivation_rejects_nonfinal_physical_or_global_source(
    boundary: _RouteBoundary,
) -> None:
    """Refuse a commit that is not the current physical and global child.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    source_commit = boundary.commits[-1]
    variants = (
        _raw_copy(
            source_commit,
            physical_child_digest=_identity(b"foreign-physical-child"),
        ),
        _raw_copy(
            source_commit,
            traversal_after=boundary.activated,
        ),
    )

    for variant in variants:
        with pytest.raises(InvalidRouteRetraceDecisionError):
            _derive_route_retrace_decision(
                physical=boundary.physical,
                terminal=boundary.terminal,
                activated=boundary.activated,
                source_commit=variant,
            )


def test_route_retrace_derivation_rejects_foreign_transaction_and_witness(
    boundary: _RouteBoundary,
) -> None:
    """Stop foreign source grammar at the route-decision owner boundary.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    source_commit = boundary.commits[-1]
    assert type(boundary.commits[0].transaction) is CandidateTransaction
    foreign_transaction = _raw_copy(
        source_commit,
        transaction=boundary.commits[0].transaction,
    )
    with pytest.raises(UnsupportedRouteRetraceError):
        _derive_route_retrace_decision(
            physical=boundary.physical,
            terminal=boundary.terminal,
            activated=boundary.activated,
            source_commit=foreign_transaction,
        )

    transaction = source_commit.transaction
    assert type(transaction) is ZeroGuideLinkTransaction
    malformed_transaction = _raw_copy(
        transaction,
        segment_witness=object(),
    )
    malformed_commit = _raw_copy(
        source_commit,
        transaction=malformed_transaction,
    )
    with pytest.raises(InvalidRouteRetraceDecisionError):
        _derive_route_retrace_decision(
            physical=boundary.physical,
            terminal=boundary.terminal,
            activated=boundary.activated,
            source_commit=malformed_commit,
        )


def test_route_retrace_derivation_rejects_cross_wired_source_operation(
    boundary: _RouteBoundary,
) -> None:
    """Bind witness ordinal, operation bytes, and physical endpoint together.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    source_commit = boundary.commits[-1]
    transaction = source_commit.transaction
    assert type(transaction) is ZeroGuideLinkTransaction
    witness = transaction.segment_witness
    wrong_ordinal = _raw_copy(
        witness,
        operation_index=witness.operation_index - 1,
    )
    wrong_operation = _raw_copy(
        witness,
        operation=boundary.physical.operations[-2],
    )
    for changed_witness in (wrong_ordinal, wrong_operation):
        changed_transaction = _raw_copy(
            transaction,
            segment_witness=changed_witness,
        )
        changed_commit = _raw_copy(
            source_commit,
            transaction=changed_transaction,
        )
        with pytest.raises(InvalidRouteRetraceDecisionError):
            _derive_route_retrace_decision(
                physical=boundary.physical,
                terminal=boundary.terminal,
                activated=boundary.activated,
                source_commit=changed_commit,
            )

    source = boundary.physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    wrong_phase = _raw_copy(
        boundary.physical,
        phase_point=source.motion.start,
    )
    with pytest.raises(InvalidRouteRetraceDecisionError):
        _derive_route_retrace_decision(
            physical=wrong_phase,
            terminal=boundary.terminal,
            activated=boundary.activated,
            source_commit=source_commit,
        )


@pytest.mark.parametrize(
    "source_field",
    ("type", "neck", "cap", "terminal"),
)
def test_route_retrace_derivation_rejects_unsupported_source_scope(
    boundary: _RouteBoundary,
    source_field: str,
) -> None:
    """Admit only one no-neck full-cap route-terminal advancing segment.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
        source_field: Independent source-grammar fault to inject.
    """
    source_commit = boundary.commits[-1]
    transaction = source_commit.transaction
    assert type(transaction) is ZeroGuideLinkTransaction
    source = boundary.physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    if source_field == "type":
        changed_source: object = boundary.physical.operations[-2]
    elif source_field == "neck":
        changed_source = _raw_copy(source, neck_scope=object())
    elif source_field == "cap":
        changed_source = _raw_copy(
            source,
            effective_cap_decision=object(),
        )
    else:
        changed_decision = _raw_copy(
            source.traversal_decision,
            makes_cursor_terminal=False,
        )
        changed_source = _raw_copy(
            source,
            traversal_decision=changed_decision,
        )
    changed_witness = _raw_copy(
        transaction.segment_witness,
        operation=changed_source,
    )
    changed_transaction = _raw_copy(
        transaction,
        segment_witness=changed_witness,
    )
    changed_commit = _raw_copy(
        source_commit,
        transaction=changed_transaction,
    )
    changed_physical = _raw_copy(
        boundary.physical,
        operations=(*boundary.physical.operations[:-1], changed_source),
    )

    with pytest.raises(UnsupportedRouteRetraceError):
        _derive_route_retrace_decision(
            physical=changed_physical,
            terminal=boundary.terminal,
            activated=boundary.activated,
            source_commit=changed_commit,
        )


def test_route_retrace_derivation_rejects_tampered_route_preimages(
    boundary: _RouteBoundary,
) -> None:
    """Reject altered activation, route ordinal, node, and terminal evidence.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    source_commit = boundary.commits[-1]
    nonterminal_cursor = _raw_copy(
        boundary.terminal.active_cursor,
        terminal=False,
    )
    terminal_cursors = list(boundary.terminal.cursors)
    terminal_cursors[boundary.terminal.active_route_index] = nonterminal_cursor
    nonterminal = _raw_copy(
        boundary.terminal,
        cursors=tuple(terminal_cursors),
    )
    wrong_index = _raw_copy(
        boundary.activated,
        active_route_index=0,
    )
    route_zero_terminal = boundary.commits[0].traversal_after
    incident_activation = route_zero_terminal.activate_next()
    variants = (
        (nonterminal, boundary.activated),
        (boundary.terminal, wrong_index),
        (boundary.terminal, incident_activation),
        (boundary.terminal, boundary.terminal),
    )

    for terminal, activated in variants:
        with pytest.raises(InvalidRouteRetraceDecisionError):
            _derive_route_retrace_decision(
                physical=boundary.physical,
                terminal=terminal,
                activated=activated,
                source_commit=source_commit,
            )


def test_route_retrace_decision_requires_no_neck_and_full_cap(
    boundary: _RouteBoundary,
) -> None:
    """Characterize the exact admitted Task 13F source policy.

    Args:
        boundary: Accepted route-one terminal physical/global prefix.
    """
    transaction = boundary.commits[-1].transaction
    assert type(transaction) is ZeroGuideLinkTransaction
    source = transaction.segment_witness.operation

    assert type(source) is AdvanceSegmentOperation
    assert type(source.neck_scope) is NoNeckScope
    assert type(source.effective_cap_decision) is FullCapDecision
    assert source.traversal_decision.makes_cursor_terminal is True
    assert transaction.passage_after is None


def test_continuation_rejects_missing_retrace_commit(
    task13f: Task13FFixture,
) -> None:
    """Require physical return adjacent to the activation it authorizes.

    Args:
        task13f: Authenticated Task 13F policy and physical authority.
    """
    valid = task13f_retrace_continuation(task13f)
    route_zero, route_one, _, route_two = valid.commits

    with pytest.raises(InvalidRouteRetraceCommitError):
        GenerationContinuation.build(
            launch_transaction=valid.launch_transaction,
            physical=valid.physical,
            traversal=valid.traversal,
            commits=(route_zero, route_one, route_two),
        )


def test_diagnostic_restored_route_two_rank_41_accepts_swept_prefix_theorem(
    task13f: Task13FFixture,
    boundary: _RouteBoundary,
) -> None:
    """Probe the dedicated two-stratum theorem for the first contained pair.

    Args:
        task13f: Authenticated Task 13F policy and physical authority.
        boundary: Accepted route-one terminal physical/global prefix.
    """
    decision = _derive_route_retrace_decision(
        physical=boundary.physical,
        terminal=boundary.terminal,
        activated=boundary.activated,
        source_commit=boundary.commits[-1],
    )
    retrace_evaluator = RouteRetraceEvaluator.build(
        evaluator=task13f.evaluator,
    )
    retrace_transaction = retrace_evaluator.evaluate(
        boundary.physical,
        decision,
    )
    restored = retrace_evaluator.commit(
        boundary.physical,
        retrace_transaction,
    )
    family = materialize_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=restored,
        traversal=boundary.activated,
    )
    candidate = family[41]
    cap_decision = candidate.effective_cap_decision
    assert type(candidate.neck_scope) is NoNeckScope
    assert type(cap_decision) is FullCapDecision
    assert cap_decision.user_cap_bytes == cap_decision.effective_cap_bytes == task13f.evaluator.user_cap.chord_ratio_bytes
    circle = candidate.motion
    link_motion = ExactSegmentMotion.build(
        restored.phase_point,
        type(restored.phase_point).build(
            circle.center.x + circle.phase_vector.x,
            circle.center.y + circle.phase_vector.y,
        ),
    )

    containment = task13f.evaluator._containment.certify_segment(
        link_motion,
        task13f.evaluator.tool_radius,
    )
    witness = MotionCertifier.build(
        stock=restored.fork_stock(),
        tool_radius=task13f.evaluator.tool_radius,
    ).certify_swept_prefix_segment(
        operation_index=len(restored.operations),
        motion=link_motion,
        user_cap=task13f.evaluator.user_cap,
        effective_cap=task13f.evaluator.user_cap,
    )

    assert containment.motion == link_motion
    assert type(witness) is SweptPrefixMotionWitness
    assert witness.event_cell_count == SWEPT_PREFIX_MOTION_STRATA == 2
    assert witness.unresolved_count == 0
