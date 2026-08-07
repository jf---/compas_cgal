"""Atomic evaluation and commit contracts for exact route retrace."""

import hashlib
import math
from dataclasses import FrozenInstanceError
from dataclasses import dataclass
from dataclasses import fields
from dataclasses import replace
from typing import NoReturn
from typing import TypeVar
from typing import cast

import pytest

import compas_cgal.adaptive.retrace_transaction as retrace_transaction_module
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.coverage import SweepWitness
from compas_cgal.adaptive.errors import CandidateStateMismatchError
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidCoverageSweepError
from compas_cgal.adaptive.errors import InvalidDepletionPolicyError
from compas_cgal.adaptive.errors import InvalidDepletionWitnessError
from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.errors import InvalidRouteRetraceTransactionError
from compas_cgal.adaptive.errors import InvalidStockAreaError
from compas_cgal.adaptive.errors import StaleRouteRetraceTransactionError
from compas_cgal.adaptive.errors import UnresolvedMotionEventError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_MOTION_STRATA
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_STRATEGY_VERSION
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_THEOREM_VERSION
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import SweptPrefixMotionWitness
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.retrace_transaction import RouteRetraceEvaluator
from compas_cgal.adaptive.retrace_transaction import RouteRetraceTransaction
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ToolRadius
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.task13f_fixture import task13f_retrace_decision
from tests.adaptive.task13f_fixture import task13f_route_one_terminal

_T = TypeVar("_T")


def _identity(label: bytes) -> IdentityDigest:
    """Return one deterministic foreign identity for a corruption test."""
    return IdentityDigest(hashlib.sha256(label).digest())


def _raw_copy(value: _T, **changes: object) -> _T:
    """Forge one dataclass shell without invoking its owner validation."""
    forged = object.__new__(type(value))
    for item in fields(value):
        object.__setattr__(
            forged,
            item.name,
            changes.get(item.name, getattr(value, item.name)),
        )
    return cast(_T, forged)


@dataclass(frozen=True)
class _RawSweptPrefixAudit:
    """Minimal corrupt audit surface read by replay witness validation."""

    strategy_version: bytes
    theorem_version: bytes
    motion_stratum_count: int


@dataclass(frozen=True)
class _RouteRetraceFixture:
    """Accepted route-one prefix and its exact retrace authority."""

    parent: GenerationState
    decision: RouteRetraceDecision
    evaluator: RouteRetraceEvaluator

    @classmethod
    def build(cls, task13f: Task13FFixture) -> "_RouteRetraceFixture":
        """Build the complete accepted prefix before any fault injection."""
        parent, terminal, commits = task13f_route_one_terminal(task13f)
        return cls(
            parent=parent,
            decision=task13f_retrace_decision(
                physical=parent,
                terminal=terminal,
                source_commit=commits[-1],
            ),
            evaluator=RouteRetraceEvaluator.build(
                evaluator=task13f.evaluator,
            ),
        )


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build one authenticated Task 13F policy authority."""
    return Task13FFixture.build()


@pytest.fixture(scope="module")
def retrace(task13f: Task13FFixture) -> _RouteRetraceFixture:
    """Build one accepted route-one retrace boundary."""
    return _RouteRetraceFixture.build(task13f)


@pytest.fixture(scope="module")
def accepted_retrace(
    retrace: _RouteRetraceFixture,
) -> tuple[RouteRetraceTransaction, GenerationState]:
    """Evaluate and independently commit one real retrace transaction."""
    transaction = retrace.evaluator.evaluate(
        retrace.parent,
        retrace.decision,
    )
    child = retrace.evaluator.commit(retrace.parent, transaction)
    return transaction, child


def _parent_identity(state: GenerationState) -> tuple[object, ...]:
    """Capture every authoritative parent owner and immutable chronology."""
    stock = state.fork_stock()
    coverage = state.clone_coverage()
    return (
        state.digest,
        state.stock_boundary_digest,
        state.stock_lineage_digest,
        tuple(witness.digest for witness in stock.lineage),
        coverage.certificate.digest,
        tuple(witness.digest for witness in coverage.lineage),
        state.phase_point,
        state.traversal.canonical_bytes,
        tuple(passage.canonical_bytes for passage in state.passages),
        tuple(operation.canonical_bytes for operation in state.operations),
    )


def _transaction_with_witness(
    transaction: RouteRetraceTransaction,
    witness: ReplayLateralWitness,
) -> RouteRetraceTransaction:
    """Submit one raw witness through the transaction owner boundary."""
    return RouteRetraceTransaction.build(
        parent_state_digest=transaction.parent_state_digest,
        decision=transaction.decision,
        segment_witness=witness,
        result_state_digest=transaction.result_state_digest,
    )


def test_route_retrace_transaction_binds_one_physical_suffix(
    retrace: _RouteRetraceFixture,
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Bind the exact reverse witness to one immutable parent and child.

    Args:
        retrace: Accepted physical parent and exact evaluator authority.
        accepted_retrace: Real transaction and independently replayed child.
    """
    transaction, child = accepted_retrace

    assert transaction.parent_state_digest == retrace.parent.digest
    assert transaction.decision == transaction.segment_witness.operation.decision
    assert transaction.result_state_digest == child.digest
    assert child.operations == (retrace.parent.operations + (transaction.segment_witness.operation,))
    assert child.phase_point == transaction.segment_witness.operation.motion.end
    assert child.traversal.canonical_bytes == retrace.parent.traversal.canonical_bytes
    assert child.passages == retrace.parent.passages


def test_route_retrace_transaction_has_one_frozen_four_field_identity(
    retrace: _RouteRetraceFixture,
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Encode and distinguish the parent, decision, proof, and child fields.

    Args:
        retrace: Accepted parent supplying the decision's owned source.
        accepted_retrace: Real accepted transaction and its reproduced child.
    """
    transaction, _ = accepted_retrace
    expected = encode_tagged_union(
        b"route-retrace-transaction-v1",
        encode_component_map(
            {
                b"decision": transaction.decision.canonical_bytes,
                b"parent-state": encode_bytes(
                    bytes(transaction.parent_state_digest),
                ),
                b"result-state": encode_bytes(
                    bytes(transaction.result_state_digest),
                ),
                b"segment-witness": transaction.segment_witness.canonical_bytes,
            }
        ),
    )

    assert transaction.canonical_bytes == expected
    assert transaction.digest == IdentityDigest(hashlib.sha256(expected).digest())

    changed_decision = replace(
        transaction.decision,
        source_commit_digest=_identity(b"changed-source-commit"),
    )
    source = retrace.parent.operations[transaction.decision.source_operation_index]
    assert type(source) is AdvanceSegmentOperation
    changed_operation = RetraceSegmentOperation.build(
        source_operation=source,
        decision=changed_decision,
    )
    changed_decision_witness = replace(
        transaction.segment_witness,
        operation=changed_operation,
    )
    sweep_parents = transaction.segment_witness.sweep_witness.parent_lineage
    assert sweep_parents
    changed_sweep = replace(
        transaction.segment_witness.sweep_witness,
        parent_lineage=(
            bytes(_identity(b"changed-opaque-coverage-parent")),
            *sweep_parents[1:],
        ),
    )
    changed_witness = replace(
        transaction.segment_witness,
        sweep_witness=changed_sweep,
    )
    variants = (
        replace(
            transaction,
            parent_state_digest=_identity(b"changed-parent-state"),
        ),
        RouteRetraceTransaction.build(
            parent_state_digest=transaction.parent_state_digest,
            decision=changed_decision,
            segment_witness=changed_decision_witness,
            result_state_digest=transaction.result_state_digest,
        ),
        replace(transaction, segment_witness=changed_witness),
        replace(
            transaction,
            result_state_digest=_identity(b"changed-result-state"),
        ),
    )
    records = (transaction.canonical_bytes,) + tuple(variant.canonical_bytes for variant in variants)
    digests = (transaction.digest,) + tuple(variant.digest for variant in variants)

    assert len(set(records)) == 5
    assert len(set(digests)) == 5
    for variant in variants:
        assert variant.digest == IdentityDigest(
            hashlib.sha256(variant.canonical_bytes).digest(),
        )
    with pytest.raises(FrozenInstanceError):
        setattr(transaction, "parent_state_digest", _identity(b"mutable-parent"))


def test_route_retrace_transaction_rejects_raw_semantic_corruptions(
    retrace: _RouteRetraceFixture,
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Reject every independently mutable operation and theorem identity.

    Args:
        retrace: Accepted parent used to obtain a foreign source operation.
        accepted_retrace: Real transaction whose nested records are corrupted.
    """
    transaction, _ = accepted_retrace
    witness = transaction.segment_witness
    operation = witness.operation
    motion_witness = witness.motion_witness
    assert type(operation) is RetraceSegmentOperation
    assert type(motion_witness) is SweptPrefixMotionWitness

    foreign_decision = replace(
        transaction.decision,
        source_commit_digest=_identity(b"foreign-source-commit"),
    )
    wrong_decision_operation = _raw_copy(
        operation,
        decision=foreign_decision,
    )
    foreign_operation = retrace.parent.operations[transaction.decision.source_operation_index]
    assert type(foreign_operation) is AdvanceSegmentOperation

    wrong_strategy = _raw_copy(
        motion_witness,
        native_audit=_RawSweptPrefixAudit(
            strategy_version=b"foreign-retrace-strategy-v1",
            theorem_version=SWEPT_PREFIX_THEOREM_VERSION,
            motion_stratum_count=SWEPT_PREFIX_MOTION_STRATA,
        ),
    )
    wrong_theorem = _raw_copy(
        motion_witness,
        native_audit=_RawSweptPrefixAudit(
            strategy_version=SWEPT_PREFIX_STRATEGY_VERSION,
            theorem_version=b"foreign-retrace-theorem-v1",
            motion_stratum_count=SWEPT_PREFIX_MOTION_STRATA,
        ),
    )
    wrong_strata = _raw_copy(
        motion_witness,
        native_audit=_RawSweptPrefixAudit(
            strategy_version=SWEPT_PREFIX_STRATEGY_VERSION,
            theorem_version=SWEPT_PREFIX_THEOREM_VERSION,
            motion_stratum_count=SWEPT_PREFIX_MOTION_STRATA + 1,
        ),
    )
    foreign_boundary = bytes(_identity(b"foreign-stock-boundary"))
    wrong_boundary_motion = _raw_copy(
        motion_witness,
        stock_boundary_digest=foreign_boundary,
    )
    wrong_boundary_witness = _raw_copy(
        witness,
        stock_boundary_digest=foreign_boundary,
        motion_witness=wrong_boundary_motion,
    )

    depletion_parents = witness.depletion_witness.parent_lineage
    assert depletion_parents
    wrong_depletion = _raw_copy(
        witness.depletion_witness,
        parent_lineage=(
            _identity(b"foreign-stock-parent"),
            *depletion_parents[1:],
        ),
    )
    wrong_coverage = _raw_copy(
        witness.sweep_witness,
        parent_lineage=(
            *witness.sweep_witness.parent_lineage,
            bytes(_identity(b"extra-coverage-parent")),
        ),
    )
    wrong_index_motion = _raw_copy(
        motion_witness,
        operation_index=motion_witness.operation_index + 1,
    )
    wrong_index_witness = _raw_copy(
        witness,
        operation_index=witness.operation_index + 1,
        motion_witness=wrong_index_motion,
    )

    corruptions = (
        _raw_copy(witness, operation=foreign_operation),
        _raw_copy(witness, operation=wrong_decision_operation),
        _raw_copy(witness, motion_witness=wrong_strategy),
        _raw_copy(witness, motion_witness=wrong_theorem),
        _raw_copy(witness, motion_witness=wrong_strata),
        wrong_boundary_witness,
        _raw_copy(witness, depletion_witness=wrong_depletion),
        _raw_copy(witness, sweep_witness=wrong_coverage),
        wrong_index_witness,
    )
    for corrupted in corruptions:
        with pytest.raises(
            InvalidRouteRetraceTransactionError,
            match="retrace|transaction|witness|chronology|parent",
        ):
            _transaction_with_witness(transaction, corrupted)


@pytest.mark.parametrize(
    ("field_name", "malformed"),
    (
        pytest.param("parent_state_digest", b"short", id="parent-short"),
        pytest.param("parent_state_digest", bytearray(32), id="parent-foreign-type"),
        pytest.param("result_state_digest", b"short", id="result-short"),
        pytest.param("result_state_digest", bytearray(32), id="result-foreign-type"),
    ),
)
def test_route_retrace_transaction_rejects_malformed_state_digest(
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
    field_name: str,
    malformed: object,
) -> None:
    """Reject malformed parent and result identities through the named owner.

    Args:
        accepted_retrace: Real transaction supplying all unmodified fields.
        field_name: Parent or result digest field under corruption.
        malformed: Wrong-size or wrong-runtime-type digest value.
    """
    transaction, _ = accepted_retrace

    with pytest.raises(InvalidRouteRetraceTransactionError, match="SHA-256"):
        replace(transaction, **{field_name: malformed})


def test_route_retrace_transaction_translates_hollow_nested_records(
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Translate every hollow exact nested shell to the transaction error.

    Args:
        accepted_retrace: Real transaction supplying complete sibling records.
    """
    transaction, _ = accepted_retrace
    witness = transaction.segment_witness
    hollows = (
        {
            "decision": object.__new__(RouteRetraceDecision),
            "segment_witness": witness,
        },
        {
            "decision": transaction.decision,
            "segment_witness": object.__new__(ReplayLateralWitness),
        },
        {
            "decision": transaction.decision,
            "segment_witness": _raw_copy(
                witness,
                operation=object.__new__(RetraceSegmentOperation),
            ),
        },
        {
            "decision": transaction.decision,
            "segment_witness": _raw_copy(
                witness,
                motion_witness=object.__new__(SweptPrefixMotionWitness),
            ),
        },
        {
            "decision": transaction.decision,
            "segment_witness": _raw_copy(
                witness,
                containment_certificate=object.__new__(
                    SegmentContainmentCertificate,
                ),
            ),
        },
        {
            "decision": transaction.decision,
            "segment_witness": _raw_copy(
                witness,
                depletion_witness=object.__new__(DepletionWitness),
            ),
        },
        {
            "decision": transaction.decision,
            "segment_witness": _raw_copy(
                witness,
                sweep_witness=object.__new__(SweepWitness),
            ),
        },
    )
    for nested in hollows:
        with pytest.raises(InvalidRouteRetraceTransactionError):
            RouteRetraceTransaction.build(
                parent_state_digest=transaction.parent_state_digest,
                decision=nested["decision"],  # type: ignore[arg-type]
                segment_witness=nested["segment_witness"],  # type: ignore[arg-type]
                result_state_digest=transaction.result_state_digest,
            )


def test_route_retrace_transaction_translates_malformed_depletion_policy(
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Translate a malformed exact nested policy through its witness owner.

    Args:
        accepted_retrace: Real transaction supplying all sibling evidence.
    """
    transaction, _ = accepted_retrace
    witness = transaction.segment_witness
    policy = witness.depletion_witness.policy
    assert type(policy) is DepletionPolicy
    malformed_policy = _raw_copy(policy, center_count_limit="4096")
    malformed_depletion = _raw_copy(
        witness.depletion_witness,
        policy=malformed_policy,
    )
    malformed_witness = _raw_copy(
        witness,
        depletion_witness=malformed_depletion,
    )

    with pytest.raises(
        InvalidRouteRetraceTransactionError,
        match="malformed nested exact state",
    ) as captured:
        _transaction_with_witness(transaction, malformed_witness)

    assert type(captured.value.__cause__) is InvalidDepletionWitnessError
    assert type(captured.value.__cause__.__cause__) is InvalidDepletionPolicyError


def test_route_retrace_evaluation_rejects_hollow_exact_state(
    retrace: _RouteRetraceFixture,
) -> None:
    """Translate an unsealed exact state before physical-policy delegation.

    Args:
        retrace: Accepted decision and exact evaluator authority.
    """
    hollow = object.__new__(GenerationState)

    with pytest.raises(
        InvalidRouteRetraceTransactionError,
        match="evaluation parent contains incomplete exact state",
    ):
        retrace.evaluator.evaluate(hollow, retrace.decision)


def test_route_retrace_commit_validates_transaction_before_hollow_state(
    retrace: _RouteRetraceFixture,
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Validate submitted evidence before reading the exact parent shell.

    Args:
        retrace: Exact evaluator authority for both commit attempts.
        accepted_retrace: Real accepted transaction for the state control.
    """
    transaction, _ = accepted_retrace
    hollow_state = object.__new__(GenerationState)
    hollow_transaction = object.__new__(RouteRetraceTransaction)

    with pytest.raises(
        InvalidRouteRetraceTransactionError,
        match="transaction contains incomplete nested exact state",
    ):
        retrace.evaluator.commit(hollow_state, hollow_transaction)
    with pytest.raises(
        InvalidRouteRetraceTransactionError,
        match="commit parent contains incomplete exact state",
    ):
        retrace.evaluator.commit(hollow_state, transaction)


def test_failed_retrace_containment_leaves_parent_byte_identical(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep every authoritative parent owner immutable after containment failure.

    Args:
        retrace: Accepted prefix constructed before fault injection.
        monkeypatch: Scoped containment failure injection.
    """
    before = _parent_identity(retrace.parent)

    def reject_containment(
        _owner: GougeContainment,
        _motion: ExactSegmentMotion,
        _tool_radius: ToolRadius,
    ) -> NoReturn:
        raise GougeContainmentError("injected retrace containment failure")

    monkeypatch.setattr(
        GougeContainment,
        "certify_segment",
        reject_containment,
    )
    with pytest.raises(
        GougeContainmentError,
        match="injected retrace containment failure",
    ):
        retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert _parent_identity(retrace.parent) == before


def test_failed_retrace_swept_prefix_leaves_parent_byte_identical(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep every authoritative parent owner immutable after theorem failure.

    Args:
        retrace: Accepted prefix constructed before fault injection.
        monkeypatch: Scoped swept-prefix theorem failure injection.
    """
    before = _parent_identity(retrace.parent)

    def reject_swept_prefix(
        _owner: MotionCertifier,
        **_kwargs: object,
    ) -> NoReturn:
        raise UnresolvedMotionEventError("injected retrace theorem failure")

    monkeypatch.setattr(
        MotionCertifier,
        "certify_swept_prefix_segment",
        reject_swept_prefix,
    )
    with pytest.raises(
        UnresolvedMotionEventError,
        match="injected retrace theorem failure",
    ):
        retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert _parent_identity(retrace.parent) == before


def test_failed_retrace_depletion_leaves_parent_byte_identical(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep every authoritative parent owner immutable after depletion failure.

    Args:
        retrace: Accepted prefix constructed before fault injection.
        monkeypatch: Scoped stock-depletion failure injection.
    """
    before = _parent_identity(retrace.parent)

    def reject_depletion(
        _owner: Stock2Area,
        _motion: ExactSegmentMotion,
        _tool_radius: ToolRadius,
        _policy: DepletionPolicy,
    ) -> NoReturn:
        raise InvalidStockAreaError("injected retrace depletion failure")

    monkeypatch.setattr(Stock2Area, "deplete", reject_depletion)
    with pytest.raises(
        InvalidStockAreaError,
        match="injected retrace depletion failure",
    ):
        retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert _parent_identity(retrace.parent) == before


def test_failed_retrace_coverage_leaves_parent_byte_identical(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep every authoritative parent owner immutable after coverage failure.

    Args:
        retrace: Accepted prefix constructed before fault injection.
        monkeypatch: Scoped coverage-ledger failure injection.
    """
    before = _parent_identity(retrace.parent)

    def reject_coverage(
        _owner: CoverageLedger,
        _motion: ExactSegmentMotion,
        _tool_radius: ToolRadius,
    ) -> NoReturn:
        raise InvalidCoverageSweepError("injected retrace coverage failure")

    monkeypatch.setattr(CoverageLedger, "add_sweep", reject_coverage)
    with pytest.raises(
        InvalidCoverageSweepError,
        match="injected retrace coverage failure",
    ):
        retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert _parent_identity(retrace.parent) == before


def test_failed_retrace_child_construction_leaves_parent_byte_identical(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep every authoritative parent owner immutable after child rejection.

    Args:
        retrace: Accepted prefix constructed before fault injection.
        monkeypatch: Scoped generation-state construction failure injection.
    """
    before = _parent_identity(retrace.parent)

    def reject_child(
        _cls: type[GenerationState],
        **_kwargs: object,
    ) -> NoReturn:
        raise InvalidGenerationStateError("injected retrace child failure")

    monkeypatch.setattr(
        GenerationState,
        "build",
        classmethod(reject_child),
    )
    with pytest.raises(
        InvalidGenerationStateError,
        match="injected retrace child failure",
    ):
        retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert _parent_identity(retrace.parent) == before


def test_route_retrace_commit_rejects_stale_parent(
    retrace: _RouteRetraceFixture,
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Reject a transaction after its authoritative physical parent advances.

    Args:
        retrace: Evaluator authority used for the stale commit attempt.
        accepted_retrace: Original transaction and its distinct child state.
    """
    transaction, child = accepted_retrace

    with pytest.raises(
        StaleRouteRetraceTransactionError,
        match="parent",
    ):
        retrace.evaluator.commit(child, transaction)


def test_route_retrace_evaluation_rejects_changed_depletion_policy(
    retrace: _RouteRetraceFixture,
) -> None:
    """Keep the candidate evaluator as the sole physical-policy authority.

    Args:
        retrace: Accepted parent and its original evaluator authority.
    """
    authority = retrace.evaluator.evaluator
    changed = CandidateEvaluator.build(
        reachable_domain=authority.reachable_domain,
        tool_radius=authority.tool_radius,
        user_cap=authority.user_cap,
        candidate_policy=authority.candidate_policy,
        neck_policy=authority.neck_policy,
        depletion_policy=DepletionPolicy.build(
            chord_bound=ChordBound.build(0.25),
            center_count_limit=4096,
        ),
        cut_direction_policy=authority.cut_direction_policy,
        cut_z=authority.cut_z,
        material_side=authority.material_side,
    )

    with pytest.raises(CandidateStateMismatchError, match="depletion policy"):
        RouteRetraceEvaluator.build(evaluator=changed).evaluate(
            retrace.parent,
            retrace.decision,
        )


def test_route_retrace_evaluation_rejects_wrong_full_cap_parent(
    retrace: _RouteRetraceFixture,
) -> None:
    """Reject a final source whose full cap differs from evaluator authority.

    Args:
        retrace: Accepted parent whose final source is raw-corrupted.
    """
    source = retrace.parent.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    foreign_cap = EngagementCap.build(math.pi / 2.0)
    wrong_source = replace(
        source,
        effective_cap_decision=FullCapDecision.build(
            user_cap=foreign_cap,
            effective_cap=foreign_cap,
        ),
    )
    wrong_parent = _raw_copy(
        retrace.parent,
        operations=retrace.parent.operations[:-1] + (wrong_source,),
    )

    with pytest.raises(CandidateStateMismatchError, match="cap decision"):
        retrace.evaluator.evaluate(wrong_parent, retrace.decision)


def test_route_retrace_evaluation_rejects_cross_wired_source_lineage(
    retrace: _RouteRetraceFixture,
) -> None:
    """Reject a valid decision that names a nonfinal physical source ordinal.

    Args:
        retrace: Accepted operation prefix supplying a foreign earlier source.
    """
    foreign_index = len(retrace.parent.operations) - 2
    foreign_source = retrace.parent.operations[foreign_index]
    cross_wired = replace(
        retrace.decision,
        source_operation_index=foreign_index,
        source_operation_digest=IdentityDigest(
            hashlib.sha256(foreign_source.canonical_bytes).digest(),
        ),
    )

    with pytest.raises(
        InvalidRouteRetraceTransactionError,
        match="final source",
    ):
        retrace.evaluator.evaluate(retrace.parent, cross_wired)


def test_route_retrace_evaluation_uses_one_physical_parent_authority(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Delegate physical-parent policy validation exactly once.

    Args:
        retrace: Accepted parent and short-lived retrace wrapper.
        monkeypatch: Scoped call-count instrumentation of the sole authority.
    """
    calls = 0
    real_validate = CandidateEvaluator._validate_physical_parent

    def track_validate(
        owner: CandidateEvaluator,
        state: GenerationState,
    ) -> None:
        nonlocal calls
        calls += 1
        real_validate(owner, state)

    monkeypatch.setattr(
        CandidateEvaluator,
        "_validate_physical_parent",
        track_validate,
    )
    retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert calls == 1


def test_route_retrace_commit_reproduces_transaction_and_child_bytes(
    retrace: _RouteRetraceFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Repeat the bounded physical proof before publishing its child.

    Args:
        retrace: Accepted parent and short-lived exact evaluator wrapper.
        monkeypatch: Scoped proof-path call-count instrumentation.
    """
    calls = 0
    real_trial = retrace_transaction_module.evaluate_retrace_segment_trial

    def track_trial(**kwargs: object) -> ReplayLateralWitness:
        nonlocal calls
        calls += 1
        return real_trial(**kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(
        retrace_transaction_module,
        "evaluate_retrace_segment_trial",
        track_trial,
    )
    transaction = retrace.evaluator.evaluate(retrace.parent, retrace.decision)
    child = retrace.evaluator.commit(retrace.parent, transaction)
    reproduced = retrace.evaluator.evaluate(retrace.parent, retrace.decision)

    assert calls == 3
    assert reproduced.canonical_bytes == transaction.canonical_bytes
    assert child.digest == transaction.result_state_digest


def test_route_retrace_commit_rejects_opaque_replay_cross_wires(
    retrace: _RouteRetraceFixture,
    accepted_retrace: tuple[RouteRetraceTransaction, GenerationState],
) -> None:
    """Reject structurally valid opaque evidence by independent byte replay.

    Args:
        retrace: Accepted parent and independent commit authority.
        accepted_retrace: Real transaction supplying structurally valid proof.
    """
    transaction, _ = accepted_retrace
    witness = transaction.segment_witness
    coverage_parents = witness.sweep_witness.parent_lineage
    assert coverage_parents
    cross_wired_sweep = _raw_copy(
        witness.sweep_witness,
        parent_lineage=(
            bytes(_identity(b"foreign-opaque-coverage-parent")),
            *coverage_parents[1:],
        ),
    )
    cross_wired_witness = _raw_copy(
        witness,
        sweep_witness=cross_wired_sweep,
    )
    cross_wired = _transaction_with_witness(
        transaction,
        cross_wired_witness,
    )
    foreign_child = replace(
        transaction,
        result_state_digest=_identity(b"foreign-retrace-child"),
    )

    for corruption in (cross_wired, foreign_child):
        with pytest.raises(
            InvalidRouteRetraceTransactionError,
            match="independent commit replay",
        ):
            retrace.evaluator.commit(retrace.parent, corruption)


def test_route_retrace_commit_rejects_hollow_exact_transaction(
    retrace: _RouteRetraceFixture,
) -> None:
    """Translate a hollow exact transaction before commit reads its fields.

    Args:
        retrace: Accepted parent and independent commit authority.
    """
    hollow = object.__new__(RouteRetraceTransaction)

    with pytest.raises(InvalidRouteRetraceTransactionError):
        retrace.evaluator.commit(retrace.parent, hollow)
