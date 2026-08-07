"""Atomic exact evaluation and independent commit for one route retrace."""

import hashlib
from dataclasses import dataclass
from typing import Self

from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import SweepWitness
from compas_cgal.adaptive.errors import InvalidContainmentCertificateError
from compas_cgal.adaptive.errors import InvalidCoverageSweepError
from compas_cgal.adaptive.errors import InvalidDepletionWitnessError
from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.errors import InvalidMotionCertificateError
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.errors import InvalidRetraceSegmentOperationError
from compas_cgal.adaptive.errors import InvalidRouteRetraceDecisionError
from compas_cgal.adaptive.errors import InvalidRouteRetraceTransactionError
from compas_cgal.adaptive.errors import StaleRouteRetraceTransactionError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_MOTION_STRATA
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_STRATEGY_VERSION
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_THEOREM_VERSION
from compas_cgal.adaptive.motion_certificate import SweptPrefixMotionWitness
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.retrace_segment_trial import evaluate_retrace_segment_trial
from compas_cgal.adaptive.transaction import CandidateEvaluator

_DIGEST_SIZE = hashlib.sha256().digest_size


def _require_digest(value: object, name: str) -> bytes:
    """Validate one opaque state identity at the transaction boundary."""
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise InvalidRouteRetraceTransactionError(
            f"{name} must be one exact SHA-256 digest.",
        )
    return value


def _stock_lineage_digest(
    lineage: tuple[IdentityDigest, ...],
) -> bytes:
    """Hash the ordered depletion prefix observed before certification."""
    if type(lineage) is not tuple or any(type(digest) is not bytes or len(digest) != _DIGEST_SIZE for digest in lineage):
        raise InvalidRouteRetraceTransactionError(
            "retrace stock parent lineage must contain exact SHA-256 digests.",
        )
    return hashlib.sha256(
        encode_tagged_union(
            b"exact-stock-lineage-v1",
            encode_sequence(tuple(bytes(digest) for digest in lineage)),
        )
    ).digest()


@dataclass(frozen=True)
class RouteRetraceTransaction:
    """Immutable physical proof binding one parent to one retrace child."""

    parent_state_digest: IdentityDigest
    decision: RouteRetraceDecision
    segment_witness: ReplayLateralWitness
    result_state_digest: IdentityDigest

    def __post_init__(self) -> None:
        RouteRetraceTransaction.validate(self)

    @classmethod
    def validate(
        cls,
        transaction: "RouteRetraceTransaction",
    ) -> None:
        """Validate every owned nested record and ordered proof invariant.

        Args:
            transaction: Exact transaction record to validate.

        Raises:
            InvalidRouteRetraceTransactionError: If a field or nested owner is
                malformed, cross-wired, or outside the fixed proof contract.
        """
        if type(transaction) is not cls:
            raise InvalidRouteRetraceTransactionError(
                "route retrace transaction must use the exact owned type.",
            )
        try:
            _require_digest(
                transaction.parent_state_digest,
                "retrace transaction parent state",
            )
            _require_digest(
                transaction.result_state_digest,
                "retrace transaction result state",
            )
            decision = transaction.decision
            witness = transaction.segment_witness
            if type(decision) is not RouteRetraceDecision:
                raise InvalidRouteRetraceTransactionError(
                    "retrace transaction requires one exact route decision.",
                )
            decision._validate()
            if type(witness) is not ReplayLateralWitness:
                raise InvalidRouteRetraceTransactionError(
                    "retrace transaction requires one exact segment witness.",
                )
            ReplayLateralWitness.validate(witness)
            operation = witness.operation
            if type(operation) is not RetraceSegmentOperation:
                raise InvalidRouteRetraceTransactionError(
                    "retrace transaction requires one exact retrace operation.",
                )
            operation._validate()
            if operation.decision != decision:
                raise InvalidRouteRetraceTransactionError(
                    "retrace operation decision contradicts its transaction.",
                )

            motion_witness = witness.motion_witness
            if (
                type(motion_witness) is not SweptPrefixMotionWitness
                or motion_witness.strategy_identity != SWEPT_PREFIX_STRATEGY_VERSION
                or motion_witness.theorem_identity != SWEPT_PREFIX_THEOREM_VERSION
                or motion_witness.event_cell_count != SWEPT_PREFIX_MOTION_STRATA
            ):
                raise InvalidRouteRetraceTransactionError(
                    "retrace transaction requires the fixed exact two-stratum theorem.",
                )
            containment = witness.containment_certificate
            if type(containment) is not SegmentContainmentCertificate:
                raise InvalidRouteRetraceTransactionError(
                    "retrace transaction requires one exact segment containment proof.",
                )
            containment.__post_init__()
            motion_witness.__post_init__()
            witness.depletion_witness.__post_init__()
            SweepWitness.validate(witness.sweep_witness)

            depletion_parents = witness.depletion_witness.parent_lineage
            if motion_witness.stock_lineage_digest != _stock_lineage_digest(
                depletion_parents,
            ):
                raise InvalidRouteRetraceTransactionError(
                    "retrace transaction violates certify-before-deplete chronology.",
                )
            coverage_parents = witness.sweep_witness.parent_lineage
            if len(depletion_parents) != len(coverage_parents) + 1 or witness.operation_index != len(depletion_parents) + 1:
                raise InvalidRouteRetraceTransactionError(
                    "retrace stock and coverage updates do not share one parent prefix.",
                )
        except AttributeError as error:
            raise InvalidRouteRetraceTransactionError(
                "retrace transaction contains incomplete nested exact state.",
            ) from error
        except (
            InvalidContainmentCertificateError,
            InvalidCoverageSweepError,
            InvalidDepletionWitnessError,
            InvalidMotionCertificateError,
            InvalidReplayTraceError,
            InvalidRetraceSegmentOperationError,
            InvalidRouteRetraceDecisionError,
        ) as error:
            raise InvalidRouteRetraceTransactionError(
                "retrace transaction contains malformed nested exact state.",
            ) from error

    @classmethod
    def build(
        cls,
        *,
        parent_state_digest: IdentityDigest,
        decision: RouteRetraceDecision,
        segment_witness: ReplayLateralWitness,
        result_state_digest: IdentityDigest,
    ) -> Self:
        """Build one validated, content-addressed retrace transaction.

        Args:
            parent_state_digest: Complete authoritative pre-state identity.
            decision: Authenticated route boundary and source lineage.
            segment_witness: Ordered exact proof for the source reversal.
            result_state_digest: Complete identity of the evaluated child.

        Returns:
            Frozen transaction binding parent, decision, proof, and child.

        Raises:
            InvalidRouteRetraceTransactionError: If any field or nested proof
                is malformed, foreign, cross-wired, or out of order.
        """
        return cls(
            parent_state_digest,
            decision,
            segment_witness,
            result_state_digest,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete four-field route-retrace acceptance record."""
        RouteRetraceTransaction.validate(self)
        return encode_tagged_union(
            b"route-retrace-transaction-v1",
            encode_component_map(
                {
                    b"decision": self.decision.canonical_bytes,
                    b"parent-state": encode_bytes(
                        bytes(self.parent_state_digest),
                    ),
                    b"result-state": encode_bytes(
                        bytes(self.result_state_digest),
                    ),
                    b"segment-witness": self.segment_witness.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class RouteRetraceEvaluator:
    """Short-lived wrapper around one authoritative physical-policy owner."""

    evaluator: CandidateEvaluator

    def __post_init__(self) -> None:
        if type(self) is not RouteRetraceEvaluator or type(self.evaluator) is not CandidateEvaluator:
            raise InvalidRouteRetraceTransactionError(
                "retrace evaluator requires one exact candidate authority.",
            )

    @classmethod
    def build(
        cls,
        *,
        evaluator: CandidateEvaluator,
    ) -> Self:
        """Bind retrace trials to the sole candidate physical authority.

        Args:
            evaluator: Exact physical-policy authority already used by the
                accepted operation prefix.

        Returns:
            Frozen short-lived retrace evaluator.

        Raises:
            InvalidRouteRetraceTransactionError: If the authority is foreign.
        """
        return cls(evaluator)

    def evaluate(
        self,
        state: GenerationState,
        decision: RouteRetraceDecision,
    ) -> RouteRetraceTransaction:
        """Evaluate one route retrace without mutating its physical parent.

        Args:
            state: Authoritative physical parent snapshot.
            decision: Authenticated final-source retrace decision.

        Returns:
            Immutable accepted proof transaction.
        """
        transaction, _ = self._evaluate_trial(state, decision)
        return transaction

    def commit(
        self,
        state: GenerationState,
        transaction: RouteRetraceTransaction,
    ) -> GenerationState:
        """Independently reproduce one retrace before returning its child.

        Args:
            state: Authoritative parent named by `transaction`.
            transaction: Previously accepted immutable retrace evidence.

        Returns:
            Independently reproduced physical child.

        Raises:
            InvalidRouteRetraceTransactionError: If types, nested evidence,
                or independent canonical bytes disagree.
            StaleRouteRetraceTransactionError: If the authoritative parent has
                changed since evaluation.
        """
        if type(state) is not GenerationState:
            raise InvalidRouteRetraceTransactionError(
                "retrace commit requires one exact generation state.",
            )
        RouteRetraceTransaction.validate(transaction)
        try:
            parent_digest = state.digest
        except (AttributeError, InvalidGenerationStateError) as error:
            raise InvalidRouteRetraceTransactionError(
                "retrace commit parent contains incomplete exact state.",
            ) from error
        if transaction.parent_state_digest != parent_digest:
            raise StaleRouteRetraceTransactionError(
                "retrace transaction parent is no longer authoritative.",
            )
        reproduced, child = self._evaluate_trial(
            state,
            transaction.decision,
        )
        if reproduced.canonical_bytes != transaction.canonical_bytes:
            raise InvalidRouteRetraceTransactionError(
                "retrace transaction differs from independent commit replay.",
            )
        return child

    def _evaluate_trial(
        self,
        state: GenerationState,
        decision: RouteRetraceDecision,
    ) -> tuple[RouteRetraceTransaction, GenerationState]:
        if type(state) is not GenerationState or type(decision) is not RouteRetraceDecision:
            raise InvalidRouteRetraceTransactionError(
                "retrace evaluation requires exact state and decision types.",
            )
        try:
            decision._validate()
        except AttributeError as error:
            raise InvalidRouteRetraceTransactionError(
                "retrace evaluation decision contains incomplete exact state.",
            ) from error
        except InvalidRouteRetraceDecisionError as error:
            raise InvalidRouteRetraceTransactionError(
                "retrace evaluation decision is malformed.",
            ) from error

        self.evaluator._validate_physical_parent(state)
        if decision.source_operation_index != len(state.operations) - 1:
            raise InvalidRouteRetraceTransactionError(
                "retrace decision must name the final source operation.",
            )
        source = state.operations[-1]
        if type(source) is not AdvanceSegmentOperation:
            raise InvalidRouteRetraceTransactionError(
                "retrace final source must be one exact advancing segment.",
            )
        try:
            operation = RetraceSegmentOperation.build(
                source_operation=source,
                decision=decision,
            )
        except InvalidRetraceSegmentOperationError as error:
            raise InvalidRouteRetraceTransactionError(
                "retrace decision contradicts its owned final source.",
            ) from error

        stock = state.fork_stock()
        coverage = state.clone_coverage()
        segment_witness = evaluate_retrace_segment_trial(
            containment_authority=self.evaluator._containment,
            stock=stock,
            coverage=coverage,
            operation_index=len(state.operations),
            operation=operation,
            tool_radius=self.evaluator.tool_radius,
            user_cap=self.evaluator.user_cap,
            effective_cap=self.evaluator.user_cap,
            depletion_policy=self.evaluator.depletion_policy,
        )
        child = GenerationState.build(
            stock=stock,
            coverage=coverage,
            tool_radius=self.evaluator.tool_radius,
            phase_point=operation.motion.end,
            traversal=state.traversal,
            passages=state.passages,
            operations=state.operations + (operation,),
        )
        transaction = RouteRetraceTransaction.build(
            parent_state_digest=state.digest,
            decision=decision,
            segment_witness=segment_witness,
            result_state_digest=child.digest,
        )
        return transaction, child
