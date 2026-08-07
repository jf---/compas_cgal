"""Atomic exact evaluation and commit for one finite-lattice candidate."""

import hashlib
from dataclasses import dataclass
from dataclasses import field
from typing import Self
from typing import TypeAlias
from typing import cast

from compas_cgal.adaptive.advancing_segment_trial import evaluate_advancing_segment_trial
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import CircleContainmentCertificate
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.errors import CandidateSelectionError
from compas_cgal.adaptive.errors import CandidateStateMismatchError
from compas_cgal.adaptive.errors import InvalidCandidateTransactionError
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.errors import InvalidZeroGuideTransactionError
from compas_cgal.adaptive.errors import StaleCandidateTransactionError
from compas_cgal.adaptive.errors import StaleZeroGuideTransactionError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_MOTION_STRATA
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_STRATEGY_VERSION
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_THEOREM_VERSION
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import SweptPrefixMotionWitness
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.toolpath import OperationType

_DIGEST_SIZE = hashlib.sha256().digest_size


def _digest(
    value: object,
    name: str,
    error_type: type[ValueError] = InvalidCandidateTransactionError,
) -> bytes:
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise error_type(f"{name} must be one exact SHA-256 digest.")
    return value


def _stock_lineage_digest(
    lineage: tuple[IdentityDigest, ...],
    error_type: type[ValueError] = InvalidCandidateTransactionError,
) -> bytes:
    if type(lineage) is not tuple or any(type(digest) is not bytes or len(digest) != _DIGEST_SIZE for digest in lineage):
        raise error_type("stock lineage must contain exact SHA-256 identities.")
    return hashlib.sha256(
        encode_tagged_union(
            b"exact-stock-lineage-v1",
            encode_sequence(tuple(bytes(digest) for digest in lineage)),
        )
    ).digest()


def _phase_point(candidate: MiddleCurveCandidate) -> Point2[WorldXY]:
    motion = candidate.motion
    return Point2[WorldXY].build(
        motion.center.x + motion.phase_vector.x,
        motion.center.y + motion.phase_vector.y,
    )


def _optional_passage_bytes(passage: NeckPassage | None) -> bytes:
    if passage is None:
        return encode_tagged_union(b"no-neck-passage-v1", b"")
    return encode_tagged_union(
        b"oriented-neck-passage-v1",
        passage.canonical_bytes,
    )


@dataclass(frozen=True)
class CandidateTransaction:
    """Immutable accepted evidence for one direct link and full circle."""

    parent_state_digest: IdentityDigest
    candidate: MiddleCurveCandidate
    link_witness: ReplayLateralWitness
    circle_witness: ReplayLateralWitness
    traversal_after: TraversalCursorState
    passage_after: NeckPassage | None
    result_state_digest: IdentityDigest

    def __post_init__(self) -> None:
        _digest(self.parent_state_digest, "transaction parent state")
        _digest(self.result_state_digest, "transaction result state")
        if type(self.candidate) is not MiddleCurveCandidate:
            raise InvalidCandidateTransactionError("transaction requires one exact finite-lattice candidate.")
        if type(self.link_witness) is not ReplayLateralWitness or type(self.circle_witness) is not ReplayLateralWitness:
            raise InvalidCandidateTransactionError("transaction requires one exact link and circle witness.")
        if type(self.traversal_after) is not TraversalCursorState:
            raise InvalidCandidateTransactionError("transaction requires one exact resulting traversal state.")
        if self.passage_after is not None and type(self.passage_after) is not NeckPassage:
            raise InvalidCandidateTransactionError("transaction passage result must use the closed exact type.")
        self._validate_pair()
        self._validate_chronology()

    @classmethod
    def build(
        cls,
        *,
        parent_state_digest: IdentityDigest,
        candidate: MiddleCurveCandidate,
        link_witness: ReplayLateralWitness,
        circle_witness: ReplayLateralWitness,
        traversal_after: TraversalCursorState,
        passage_after: NeckPassage | None,
        result_state_digest: IdentityDigest,
    ) -> Self:
        """Build and cross-validate one accepted candidate transaction.

        Args:
            parent_state_digest: Complete authoritative pre-state identity.
            candidate: Finite-lattice circle proposal and decisions.
            link_witness: Ordered proof bundle for the direct phase link.
            circle_witness: Ordered proof bundle for the candidate circle.
            traversal_after: Cursor state after the circle advance.
            passage_after: Oriented passage result, or `None` without a neck.
            result_state_digest: Complete identity of the evaluated child.

        Returns:
            Frozen content-addressed acceptance evidence.

        Raises:
            InvalidCandidateTransactionError: If any proof is cross-wired,
                reordered, or inconsistent with the candidate.
        """
        return cls(
            parent_state_digest,
            candidate,
            link_witness,
            circle_witness,
            traversal_after,
            passage_after,
            result_state_digest,
        )

    def _validate_pair(self) -> None:
        link = self.link_witness
        circle = self.circle_witness
        if type(link.operation) is not LinkSegmentOperation or type(circle.operation) is not CutFullCircleOperation:
            raise InvalidCandidateTransactionError("transaction requires canonical link then circle proof bundles.")
        if circle.operation_index != link.operation_index + 1:
            raise InvalidCandidateTransactionError("transaction link and circle operation indices are not consecutive.")
        candidate = self.candidate
        operation = circle.operation
        if (
            operation.motion != candidate.motion
            or operation.neck_scope != candidate.neck_scope
            or operation.effective_cap_decision != candidate.effective_cap_decision
            or operation.traversal_decision != candidate.traversal_decision
        ):
            raise InvalidCandidateTransactionError("transaction circle proof does not belong to its candidate.")
        traversal = candidate.traversal_decision
        expected_hold = HoldTraversalDecision.build(
            component_id=traversal.component_id,
            edge_id=traversal.edge_id,
            branch_id=traversal.branch_id,
            cursor=traversal.cursor_before,
        )
        if (
            link.operation.neck_scope != candidate.neck_scope
            or link.operation.effective_cap_decision != candidate.effective_cap_decision
            or link.operation.traversal_decision != expected_hold
            or link.operation.motion.end != _phase_point(candidate)
        ):
            raise InvalidCandidateTransactionError("transaction direct link does not belong to its candidate circle.")
        expected_traversal = TraversalCursorState.before(traversal).advance(traversal)
        if self.traversal_after != expected_traversal:
            raise InvalidCandidateTransactionError("transaction result traversal contradicts its candidate.")

        decision = candidate.effective_cap_decision
        if type(candidate.neck_scope) is NoNeckScope:
            if self.passage_after is not None or type(decision) is not FullCapDecision:
                raise InvalidCandidateTransactionError("no-neck transaction cannot advance an oriented passage.")
            return
        if type(candidate.neck_scope) is not OrientedNeckScope or type(decision) is not NeckCapDecision:
            raise InvalidCandidateTransactionError("transaction candidate uses a foreign neck decision.")
        passage = self.passage_after
        if (
            type(passage) is not NeckPassage
            or passage.scope != candidate.neck_scope
            or passage.neck.evidence_digest != decision.neck_evidence_digest
            or passage.neck.width_class_id != decision.width_class_id
            or passage.state is not decision.passage_after
        ):
            raise InvalidCandidateTransactionError("transaction oriented passage result contradicts its candidate.")

    def _validate_chronology(self) -> None:
        link = self.link_witness
        circle = self.circle_witness
        link_parents = link.depletion_witness.parent_lineage
        if link.motion_witness.stock_lineage_digest != _stock_lineage_digest(link_parents):
            raise InvalidCandidateTransactionError("transaction violates certify-before-deplete link chronology.")
        expected_circle_parents = link_parents + (link.depletion_witness.digest,)
        if circle.depletion_witness.parent_lineage != expected_circle_parents:
            raise InvalidCandidateTransactionError("transaction circle depletion does not follow its link.")
        if circle.motion_witness.stock_lineage_digest != _stock_lineage_digest(expected_circle_parents):
            raise InvalidCandidateTransactionError("transaction violates certify-before-deplete circle chronology.")
        expected_circle_sweeps = link.sweep_witness.parent_lineage + (link.sweep_witness.digest,)
        if circle.sweep_witness.parent_lineage != expected_circle_sweeps:
            raise InvalidCandidateTransactionError("transaction circle coverage does not follow its link.")

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete versioned acceptance record.

        Returns:
            Canonical CCAN bytes binding parent, proofs, and child identity.
        """
        return encode_tagged_union(
            b"candidate-transaction-v1",
            encode_component_map(
                {
                    b"candidate": self.candidate.canonical_bytes,
                    b"circle-witness": self.circle_witness.canonical_bytes,
                    b"link-witness": self.link_witness.canonical_bytes,
                    b"parent-state": encode_bytes(bytes(self.parent_state_digest)),
                    b"passage-after": _optional_passage_bytes(self.passage_after),
                    b"result-state": encode_bytes(bytes(self.result_state_digest)),
                    b"traversal-after": self.traversal_after.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`.

        Returns:
            Content identity of the accepted transaction.
        """
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class ZeroGuideLinkTransaction:
    """Immutable accepted evidence for one advancing zero-guide segment."""

    parent_state_digest: IdentityDigest
    candidate: ZeroGuideLinkCandidate
    segment_witness: ReplayLateralWitness
    traversal_after: TraversalCursorState
    passage_after: NeckPassage | None
    result_state_digest: IdentityDigest

    def __post_init__(self) -> None:
        _digest(
            self.parent_state_digest,
            "zero-guide transaction parent state",
            InvalidZeroGuideTransactionError,
        )
        _digest(
            self.result_state_digest,
            "zero-guide transaction result state",
            InvalidZeroGuideTransactionError,
        )
        if type(self.candidate) is not ZeroGuideLinkCandidate:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction requires one exact spatial-only candidate.",
            )
        if type(self.segment_witness) is not ReplayLateralWitness:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction requires one exact segment witness.",
            )
        if type(self.traversal_after) is not TraversalCursorState:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction requires one exact resulting traversal state.",
            )
        if self.passage_after is not None and type(self.passage_after) is not NeckPassage:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction passage result must use the closed exact type.",
            )
        self._validate_candidate_result()
        self._validate_chronology()

    @classmethod
    def build(
        cls,
        *,
        parent_state_digest: IdentityDigest,
        candidate: ZeroGuideLinkCandidate,
        segment_witness: ReplayLateralWitness,
        traversal_after: TraversalCursorState,
        passage_after: NeckPassage | None,
        result_state_digest: IdentityDigest,
    ) -> Self:
        """Build one content-addressed link-only acceptance record.

        Args:
            parent_state_digest: Complete authoritative pre-state identity.
            candidate: Spatial-only proposal retaining its native zero-guide
                proof.
            segment_witness: Ordered proof bundle for the advancing segment.
            traversal_after: Cursor state after the owned MAT advance.
            passage_after: Oriented passage result, or `None` without a neck.
            result_state_digest: Complete identity of the evaluated child.

        Returns:
            Frozen one-witness acceptance evidence.

        Raises:
            InvalidZeroGuideTransactionError: If any candidate field, proof,
                lineage, traversal result, or passage result is cross-wired.
        """
        return cls(
            parent_state_digest,
            candidate,
            segment_witness,
            traversal_after,
            passage_after,
            result_state_digest,
        )

    def _validate_candidate_result(self) -> None:
        candidate = self.candidate
        witness = self.segment_witness
        operation = witness.operation
        if type(operation) is not AdvanceSegmentOperation:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction requires one advancing segment operation.",
            )
        try:
            ReplayLateralWitness.validate(witness)
        except InvalidReplayTraceError as error:
            raise InvalidZeroGuideTransactionError(
                f"zero-guide swept-prefix replay witness is invalid ({error}).",
            ) from error
        motion_witness = witness.motion_witness
        if (
            type(motion_witness) is not SweptPrefixMotionWitness
            or motion_witness.strategy_identity != SWEPT_PREFIX_STRATEGY_VERSION
            or motion_witness.theorem_identity != SWEPT_PREFIX_THEOREM_VERSION
            or motion_witness.event_cell_count != SWEPT_PREFIX_MOTION_STRATA
        ):
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction requires the exact swept-prefix theorem witness.",
            )
        if (
            operation.motion.end != candidate.target
            or operation.neck_scope != candidate.neck_scope
            or operation.effective_cap_decision != candidate.effective_cap_decision
            or operation.traversal_decision != candidate.traversal_decision
        ):
            raise InvalidZeroGuideTransactionError(
                "zero-guide segment proof does not belong to its candidate.",
            )
        expected_traversal = TraversalCursorState.before(
            candidate.traversal_decision,
        ).advance(candidate.traversal_decision)
        if self.traversal_after != expected_traversal:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction result traversal contradicts its candidate.",
            )

        decision = candidate.effective_cap_decision
        if type(candidate.neck_scope) is NoNeckScope:
            if self.passage_after is not None or type(decision) is not FullCapDecision:
                raise InvalidZeroGuideTransactionError(
                    "no-neck zero-guide transaction cannot advance an oriented passage.",
                )
            return
        if type(candidate.neck_scope) is not OrientedNeckScope or type(decision) is not NeckCapDecision:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction candidate uses a foreign neck decision.",
            )
        passage = self.passage_after
        if (
            type(passage) is not NeckPassage
            or passage.scope != candidate.neck_scope
            or passage.neck.evidence_digest != decision.neck_evidence_digest
            or passage.neck.width_class_id != decision.width_class_id
            or passage.state is not decision.passage_after
        ):
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction passage result contradicts its candidate.",
            )

    def _validate_chronology(self) -> None:
        witness = self.segment_witness
        depletion_parents = witness.depletion_witness.parent_lineage
        sweep_parents = witness.sweep_witness.parent_lineage
        if witness.motion_witness.stock_lineage_digest != _stock_lineage_digest(
            depletion_parents,
            InvalidZeroGuideTransactionError,
        ):
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction violates certify-before-deplete chronology.",
            )
        if len(depletion_parents) != len(sweep_parents) + 1 or witness.operation_index != len(depletion_parents) + 1:
            raise InvalidZeroGuideTransactionError(
                "zero-guide stock and coverage updates do not share one parent prefix.",
            )

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete versioned one-witness acceptance record.

        Returns:
            Canonical CCAN bytes binding parent, candidate, proof, and child.
        """
        return encode_tagged_union(
            b"zero-guide-link-transaction-v1",
            encode_component_map(
                {
                    b"candidate": self.candidate.canonical_bytes,
                    b"parent-state": encode_bytes(bytes(self.parent_state_digest)),
                    b"passage-after": _optional_passage_bytes(self.passage_after),
                    b"result-state": encode_bytes(bytes(self.result_state_digest)),
                    b"segment-witness": self.segment_witness.canonical_bytes,
                    b"traversal-after": self.traversal_after.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`.

        Returns:
            Content identity of the accepted zero-guide transaction.
        """
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


AcceptedCandidateTransaction: TypeAlias = CandidateTransaction | ZeroGuideLinkTransaction


@dataclass(frozen=True)
class CandidateEvaluator:
    """Short-lived exact authority for isolated candidate trials."""

    reachable_domain: ReachableDomain
    tool_radius: ToolRadius
    user_cap: EngagementCap
    candidate_policy: CandidatePolicy
    neck_policy: NeckPolicy
    depletion_policy: DepletionPolicy
    cut_direction_policy: CutDirectionPolicy
    cut_z: CutZ
    material_side: MaterialSide
    _containment: GougeContainment = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if type(self.reachable_domain) is not ReachableDomain:
            raise CandidateStateMismatchError("candidate evaluator requires one exact reachable domain.")
        if type(self.tool_radius) is not ToolRadius:
            raise CandidateStateMismatchError("candidate evaluator requires one exact tool radius.")
        if self.reachable_domain.certificate.tool_radius != self.tool_radius:
            raise CandidateStateMismatchError("candidate evaluator tool radius contradicts its domain.")
        if type(self.user_cap) is not EngagementCap:
            raise CandidateStateMismatchError("candidate evaluator requires one exact user cap.")
        if type(self.candidate_policy) is not CandidatePolicy:
            raise CandidateStateMismatchError("candidate evaluator requires one exact candidate policy.")
        if type(self.neck_policy) is not NeckPolicy or self.neck_policy.user_cap != self.user_cap:
            raise CandidateStateMismatchError("candidate evaluator neck policy contradicts its user cap.")
        if type(self.depletion_policy) is not DepletionPolicy:
            raise CandidateStateMismatchError("candidate evaluator requires one exact depletion policy.")
        if type(self.cut_direction_policy) is not CutDirectionPolicy:
            raise CandidateStateMismatchError("candidate evaluator requires one exact cut-direction policy.")
        if type(self.cut_z) is not CutZ:
            raise CandidateStateMismatchError("candidate evaluator requires one exact cut depth.")
        if type(self.material_side) is not MaterialSide:
            raise CandidateStateMismatchError("candidate evaluator requires one exact material side.")
        object.__setattr__(
            self,
            "_containment",
            GougeContainment.build(self.reachable_domain),
        )

    @classmethod
    def build(
        cls,
        *,
        reachable_domain: ReachableDomain,
        tool_radius: ToolRadius,
        user_cap: EngagementCap,
        candidate_policy: CandidatePolicy,
        neck_policy: NeckPolicy,
        depletion_policy: DepletionPolicy,
        cut_direction_policy: CutDirectionPolicy,
        cut_z: CutZ,
        material_side: MaterialSide,
    ) -> Self:
        """Build an evaluator bound to all invariant proof authorities.

        Args:
            reachable_domain: Exact containment authority.
            tool_radius: Cutter radius bound to domain and state.
            user_cap: User engagement limit.
            candidate_policy: Authoritative finite-lattice policy.
            neck_policy: Exact state-dependent effective-cap mapping.
            depletion_policy: Exact stock-construction bounds.
            cut_direction_policy: Cut intent mapping material side to circle
                orientation.
            cut_z: Cut-plane depth for emitted operations.
            material_side: Radial side causing circle orientation.

        Returns:
            Validated short-lived evaluator.

        Raises:
            CandidateStateMismatchError: If invariant authorities disagree.
        """
        return cls(
            reachable_domain,
            tool_radius,
            user_cap,
            candidate_policy,
            neck_policy,
            depletion_policy,
            cut_direction_policy,
            cut_z,
            material_side,
        )

    def evaluate(
        self,
        state: GenerationState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        """Evaluate one candidate without changing authoritative state.

        Args:
            state: Immutable parent snapshot.
            candidate: Exact finite-lattice proposal.

        Returns:
            Content-addressed joint link-and-circle evidence.

        Raises:
            CandidateStateMismatchError: If state or candidate identity does
                not begin at the evaluator's exact authority.
            GougeContainmentError: If either cutter sweep leaves the pocket.
            EngagementCapExceededError: If either exact motion exceeds cap.
            UnresolvedMotionEventError: If either event proof is incomplete.
        """
        return self.evaluate_from_cursor(
            state,
            state.traversal,
            candidate,
        )

    def evaluate_from_cursor(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> CandidateTransaction:
        """Evaluate through one explicitly authenticated physical cursor.

        Args:
            state: Immutable physical parent snapshot.
            traversal_before: Active MAT cursor authorized by the caller's
                global traversal state.
            candidate: Exact finite-lattice proposal beginning at that cursor.

        Returns:
            Content-addressed joint link-and-circle evidence.

        Raises:
            CandidateStateMismatchError: If physical state, cursor, candidate,
                or evaluator authority disagree.
            GougeContainmentError: If either cutter sweep leaves the pocket.
            EngagementCapExceededError: If either exact motion exceeds cap.
            UnresolvedMotionEventError: If either event proof is incomplete.
        """
        transaction, _ = self._evaluate_trial(
            state,
            traversal_before,
            candidate,
        )
        return transaction

    def evaluate_zero_guide_from_cursor(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: ZeroGuideLinkCandidate,
    ) -> ZeroGuideLinkTransaction:
        """Evaluate one proved link-only candidate without mutating its parent.

        Args:
            state: Immutable physical parent snapshot.
            traversal_before: Active MAT cursor authenticated by global
                traversal.
            candidate: Spatial-only candidate beginning at that cursor.

        Returns:
            Content-addressed one-segment acceptance evidence.

        Raises:
            CandidateStateMismatchError: If physical state, cursor, candidate,
                or evaluator authority disagree.
            GougeContainmentError: If the segment sweep leaves the pocket.
            UnresolvedMotionEventError: If the clear-start exact-pi theorem
                cannot certify the advancing segment.
        """
        transaction, _ = self._evaluate_zero_guide_trial(
            state,
            traversal_before,
            candidate,
        )
        return transaction

    def commit(
        self,
        state: GenerationState,
        transaction: CandidateTransaction,
    ) -> GenerationState:
        """Independently replay and functionally commit one winner.

        Args:
            state: Authoritative parent named by `transaction`.
            transaction: Previously accepted immutable evidence.

        Returns:
            New authoritative child snapshot.

        Raises:
            StaleCandidateTransactionError: If the parent digest changed.
            InvalidCandidateTransactionError: If independent evaluation does
                not reproduce byte-identical evidence.
        """
        return self.commit_from_cursor(
            state,
            state.traversal,
            transaction,
        )

    def commit_from_cursor(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        transaction: CandidateTransaction,
    ) -> GenerationState:
        """Commit one winner from an explicitly authenticated MAT cursor.

        Args:
            state: Authoritative physical parent named by `transaction`.
            traversal_before: Same active cursor used during evaluation.
            transaction: Previously accepted immutable evidence.

        Returns:
            New authoritative physical child snapshot.

        Raises:
            StaleCandidateTransactionError: If the physical parent changed.
            InvalidCandidateTransactionError: If types are foreign or
                independent replay differs.
            CandidateStateMismatchError: If the supplied cursor is stale.
        """
        if type(state) is not GenerationState or type(traversal_before) is not TraversalCursorState or type(transaction) is not CandidateTransaction:
            raise InvalidCandidateTransactionError(
                "candidate commit requires exact state, cursor, and transaction types.",
            )
        if transaction.parent_state_digest != state.digest:
            raise StaleCandidateTransactionError("candidate transaction parent no longer names authoritative state.")
        reproduced, next_state = self._evaluate_trial(
            state,
            traversal_before,
            transaction.candidate,
        )
        if reproduced.canonical_bytes != transaction.canonical_bytes:
            raise InvalidCandidateTransactionError("candidate transaction differs from independent commit replay.")
        return next_state

    def commit_zero_guide_from_cursor(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        transaction: ZeroGuideLinkTransaction,
    ) -> GenerationState:
        """Independently replay and commit one zero-guide winner.

        Args:
            state: Authoritative physical parent named by `transaction`.
            traversal_before: Same authenticated cursor used for evaluation.
            transaction: Previously accepted one-segment evidence.

        Returns:
            New authoritative physical child snapshot.

        Raises:
            InvalidZeroGuideTransactionError: If types are foreign or replay
                does not reproduce byte-identical evidence.
            StaleZeroGuideTransactionError: If the physical parent or local
                traversal cursor changed after evaluation.
        """
        if type(state) is not GenerationState or type(traversal_before) is not TraversalCursorState or type(transaction) is not ZeroGuideLinkTransaction:
            raise InvalidZeroGuideTransactionError(
                "zero-guide commit requires exact state, cursor, and transaction types.",
            )
        if transaction.parent_state_digest != state.digest:
            raise StaleZeroGuideTransactionError(
                "zero-guide transaction parent no longer names authoritative state.",
            )
        expected_cursor = TraversalCursorState.before(
            transaction.candidate.traversal_decision,
        )
        if traversal_before != expected_cursor:
            raise StaleZeroGuideTransactionError(
                "zero-guide transaction cursor no longer names its evaluated route position.",
            )
        reproduced, next_state = self._evaluate_zero_guide_trial(
            state,
            traversal_before,
            transaction.candidate,
        )
        if reproduced.canonical_bytes != transaction.canonical_bytes:
            raise InvalidZeroGuideTransactionError(
                "zero-guide transaction differs from independent commit replay.",
            )
        return next_state

    def _evaluate_trial(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> tuple[CandidateTransaction, GenerationState]:
        self._validate_parent(
            state,
            traversal_before,
            candidate,
        )
        traversal_after = traversal_before.advance(
            candidate.traversal_decision,
        )
        effective_cap, passage_after = self._effective_cap(state, candidate)
        phase_point = _phase_point(candidate)
        link_motion = ExactSegmentMotion.build(
            state.phase_point,
            phase_point,
        )
        traversal = candidate.traversal_decision
        link_operation = LinkSegmentOperation.build(
            motion=link_motion,
            cut_z=self.cut_z,
            neck_scope=candidate.neck_scope,
            effective_cap_decision=candidate.effective_cap_decision,
            traversal_decision=HoldTraversalDecision.build(
                component_id=traversal.component_id,
                edge_id=traversal.edge_id,
                branch_id=traversal.branch_id,
                cursor=traversal.cursor_before,
            ),
        )
        circle_operation = CutFullCircleOperation.build(
            motion=candidate.motion,
            cut_z=self.cut_z,
            material_side=self.material_side,
            neck_scope=candidate.neck_scope,
            effective_cap_decision=candidate.effective_cap_decision,
            traversal_decision=traversal,
        )

        stock = state.fork_stock()
        coverage = state.clone_coverage()
        link_index = len(state.operations)
        link_witness = self._evaluate_link(
            stock=stock,
            coverage=coverage,
            operation_index=link_index,
            operation=link_operation,
            effective_cap=effective_cap,
        )
        circle_witness = self._evaluate_circle(
            stock=stock,
            coverage=coverage,
            operation_index=link_index + 1,
            operation=circle_operation,
            effective_cap=effective_cap,
        )
        passages = self._next_passages(
            state.passages,
            passage_after,
        )
        operations: tuple[CanonicalOperation, ...] = state.operations + (
            link_operation,
            circle_operation,
        )
        next_state = GenerationState.build(
            stock=stock,
            coverage=coverage,
            tool_radius=self.tool_radius,
            phase_point=phase_point,
            traversal=traversal_after,
            passages=passages,
            operations=operations,
        )
        transaction = CandidateTransaction.build(
            parent_state_digest=state.digest,
            candidate=candidate,
            link_witness=link_witness,
            circle_witness=circle_witness,
            traversal_after=traversal_after,
            passage_after=passage_after,
            result_state_digest=next_state.digest,
        )
        return transaction, next_state

    def _evaluate_zero_guide_trial(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: ZeroGuideLinkCandidate,
    ) -> tuple[ZeroGuideLinkTransaction, GenerationState]:
        self._validate_zero_guide_parent(
            state,
            traversal_before,
            candidate,
        )
        traversal_after = traversal_before.advance(
            candidate.traversal_decision,
        )
        effective_cap, passage_after = self._effective_cap(state, candidate)
        operation = AdvanceSegmentOperation.build(
            motion=ExactSegmentMotion.build(
                state.phase_point,
                candidate.target,
            ),
            cut_z=self.cut_z,
            neck_scope=candidate.neck_scope,
            effective_cap_decision=candidate.effective_cap_decision,
            traversal_decision=candidate.traversal_decision,
        )
        stock = state.fork_stock()
        coverage = state.clone_coverage()
        segment_witness = evaluate_advancing_segment_trial(
            containment_authority=self._containment,
            stock=stock,
            coverage=coverage,
            operation_index=len(state.operations),
            operation=operation,
            tool_radius=self.tool_radius,
            user_cap=self.user_cap,
            effective_cap=effective_cap,
            depletion_policy=self.depletion_policy,
        )
        passages = self._next_passages(
            state.passages,
            passage_after,
        )
        next_state = GenerationState.build(
            stock=stock,
            coverage=coverage,
            tool_radius=self.tool_radius,
            phase_point=candidate.target,
            traversal=traversal_after,
            passages=passages,
            operations=state.operations + (operation,),
        )
        transaction = ZeroGuideLinkTransaction.build(
            parent_state_digest=state.digest,
            candidate=candidate,
            segment_witness=segment_witness,
            traversal_after=traversal_after,
            passage_after=passage_after,
            result_state_digest=next_state.digest,
        )
        return transaction, next_state

    def _validate_parent(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: MiddleCurveCandidate,
    ) -> None:
        if type(state) is not GenerationState or type(traversal_before) is not TraversalCursorState or type(candidate) is not MiddleCurveCandidate:
            raise CandidateStateMismatchError(
                "candidate evaluation requires exact state, cursor, and candidate types.",
            )
        if candidate.policy != self.candidate_policy:
            raise CandidateStateMismatchError("candidate policy contradicts its evaluator.")
        expected_orientation = self.cut_direction_policy.circle_orientation(
            self.material_side,
        )
        if candidate.proposal.circle_orientation is not expected_orientation:
            raise CandidateStateMismatchError("candidate circle contradicts evaluator cut direction.")
        self._validate_physical_parent(state)
        traversal_before.advance(candidate.traversal_decision)

    def _validate_zero_guide_parent(
        self,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: ZeroGuideLinkCandidate,
    ) -> None:
        if type(state) is not GenerationState or type(traversal_before) is not TraversalCursorState or type(candidate) is not ZeroGuideLinkCandidate:
            raise CandidateStateMismatchError(
                "zero-guide evaluation requires exact state, cursor, and candidate types.",
            )
        if candidate.policy != self.candidate_policy:
            raise CandidateStateMismatchError(
                "zero-guide candidate policy contradicts its evaluator.",
            )
        if bytes(candidate.zero_guide_run.edge_id) != bytes(
            candidate.traversal_decision.edge_id,
        ):
            raise CandidateStateMismatchError(
                "zero-guide native proof contradicts its traversal edge.",
            )
        self._validate_physical_parent(state)
        traversal_before.advance(candidate.traversal_decision)

    def _validate_physical_parent(
        self,
        state: GenerationState,
    ) -> None:
        if state.tool_radius != self.tool_radius:
            raise CandidateStateMismatchError("candidate state tool radius contradicts its evaluator.")
        stock = state.fork_stock()
        if not stock.lineage or type(stock.lineage[0].motion) is not PreclearedEntry:
            raise CandidateStateMismatchError("candidate state has no exact precleared-entry root.")
        entry = stock.lineage[0].motion
        if entry.reachable_domain.certificate.digest != self.reachable_domain.certificate.digest:
            raise CandidateStateMismatchError("candidate state reachable domain contradicts its evaluator.")
        if entry.cut_plane.cut_z != self.cut_z:
            raise CandidateStateMismatchError("candidate evaluator cut depth contradicts its qualified entry.")
        if any(witness.policy != self.depletion_policy for witness in stock.lineage[1:]):
            raise CandidateStateMismatchError("candidate parent depletion policy contradicts its evaluator.")
        self._validate_parent_operations(state)

    def _validate_parent_operations(
        self,
        state: GenerationState,
    ) -> None:
        expected_full_cap = FullCapDecision.build(
            user_cap=self.user_cap,
            effective_cap=self.user_cap,
        )
        for raw_operation in state.operations[2:]:
            if type(raw_operation) is RetraceSegmentOperation:
                if raw_operation.effective_cap_decision != expected_full_cap:
                    raise CandidateStateMismatchError(
                        "candidate parent retrace contradicts evaluator full cap.",
                    )
                continue
            if type(raw_operation) not in (
                LinkSegmentOperation,
                CutFullCircleOperation,
                AdvanceSegmentOperation,
            ):
                raise CandidateStateMismatchError("candidate parent contains a foreign lateral operation.")
            operation = cast(
                LinkSegmentOperation | CutFullCircleOperation | AdvanceSegmentOperation,
                raw_operation,
            )
            decision = operation.effective_cap_decision
            if type(operation.neck_scope) is NoNeckScope:
                expected_decision: FullCapDecision | NeckCapDecision = expected_full_cap
            else:
                if type(operation.neck_scope) is not OrientedNeckScope or type(decision) is not NeckCapDecision:
                    raise CandidateStateMismatchError("candidate parent cap decision is outside the exact grammar.")
                passage = state.passage(operation.neck_scope)
                effective_cap = self.neck_policy.effective_cap(
                    passage.neck.width_class_id.value,
                    decision.passage_before,
                )
                expected_decision = NeckCapDecision.build(
                    neck_evidence_digest=passage.neck.evidence_digest,
                    width_class_id=passage.neck.width_class_id,
                    passage_before=decision.passage_before,
                    passage_after=decision.passage_after,
                    user_cap=self.user_cap,
                    effective_cap=effective_cap,
                )
            if decision != expected_decision:
                raise CandidateStateMismatchError("candidate parent cap decision contradicts evaluator policy.")
            if type(operation) is CutFullCircleOperation:
                expected_orientation = self.cut_direction_policy.circle_orientation(
                    operation.material_side,
                )
                if operation.motion.clockwise != expected_orientation.clockwise:
                    raise CandidateStateMismatchError("candidate parent circle contradicts evaluator cut direction.")

    def _effective_cap(
        self,
        state: GenerationState,
        candidate: MiddleCurveCandidate | ZeroGuideLinkCandidate,
    ) -> tuple[EngagementCap, NeckPassage | None]:
        decision = candidate.effective_cap_decision
        if type(candidate.neck_scope) is NoNeckScope:
            expected_full_cap = FullCapDecision.build(
                user_cap=self.user_cap,
                effective_cap=self.user_cap,
            )
            if type(decision) is not FullCapDecision or decision != expected_full_cap:
                raise CandidateStateMismatchError("candidate full-cap decision contradicts evaluator policy.")
            return self.user_cap, None
        if type(candidate.neck_scope) is not OrientedNeckScope:
            raise CandidateStateMismatchError("candidate neck scope is outside the closed exact grammar.")
        passage = state.passage(candidate.neck_scope)
        expected_neck_cap = passage.propose_cap_decision(self.neck_policy)
        if type(decision) is not NeckCapDecision or decision != expected_neck_cap:
            raise CandidateStateMismatchError("candidate neck-cap decision contradicts authoritative passage.")
        effective_cap = self.neck_policy.effective_cap(
            passage.neck.width_class_id.value,
            passage.state,
        )
        return effective_cap, passage.advance(expected_neck_cap)

    def _evaluate_link(
        self,
        *,
        stock: Stock2Area,
        coverage: CoverageLedger,
        operation_index: int,
        operation: LinkSegmentOperation,
        effective_cap: EngagementCap,
    ) -> ReplayLateralWitness:
        return self._evaluate_segment(
            stock=stock,
            coverage=coverage,
            operation_index=operation_index,
            operation=operation,
            effective_cap=effective_cap,
        )

    def _evaluate_segment(
        self,
        *,
        stock: Stock2Area,
        coverage: CoverageLedger,
        operation_index: int,
        operation: LinkSegmentOperation,
        effective_cap: EngagementCap,
    ) -> ReplayLateralWitness:
        containment: SegmentContainmentCertificate = self._containment.certify_segment(
            operation.motion,
            self.tool_radius,
        )
        certifier = MotionCertifier.build(
            stock=stock,
            tool_radius=self.tool_radius,
        )
        motion_witness = certifier.certify(
            operation_index=operation_index,
            operation_kind=OperationType.LINK,
            motion=operation.motion,
            user_cap=self.user_cap,
            effective_cap=effective_cap,
        )
        depletion = stock.deplete(
            operation.motion,
            self.tool_radius,
            self.depletion_policy,
        )
        sweep = coverage.add_sweep(
            operation.motion,
            self.tool_radius,
        )
        return ReplayLateralWitness(
            operation_index=operation_index,
            operation=operation,
            effective_cap_decision=operation.effective_cap_decision,
            stock_boundary_digest=certifier.canonical_boundary_digest,
            containment_certificate=containment,
            motion_witness=motion_witness,
            depletion_witness=depletion,
            sweep_witness=sweep,
        )

    def _evaluate_circle(
        self,
        *,
        stock: Stock2Area,
        coverage: CoverageLedger,
        operation_index: int,
        operation: CutFullCircleOperation,
        effective_cap: EngagementCap,
    ) -> ReplayLateralWitness:
        containment: CircleContainmentCertificate = self._containment.certify_full_circle(
            operation.motion,
            self.tool_radius,
        )
        certifier = MotionCertifier.build(
            stock=stock,
            tool_radius=self.tool_radius,
        )
        motion_witness = certifier.certify(
            operation_index=operation_index,
            operation_kind=OperationType.CUT,
            motion=operation.motion,
            user_cap=self.user_cap,
            effective_cap=effective_cap,
        )
        depletion = stock.deplete(
            operation.motion,
            self.tool_radius,
            self.depletion_policy,
        )
        sweep = coverage.add_sweep(
            operation.motion,
            self.tool_radius,
        )
        return ReplayLateralWitness(
            operation_index=operation_index,
            operation=operation,
            effective_cap_decision=operation.effective_cap_decision,
            stock_boundary_digest=certifier.canonical_boundary_digest,
            containment_certificate=containment,
            motion_witness=motion_witness,
            depletion_witness=depletion,
            sweep_witness=sweep,
        )

    @staticmethod
    def _next_passages(
        passages: tuple[NeckPassage, ...],
        passage_after: NeckPassage | None,
    ) -> tuple[NeckPassage, ...]:
        if passage_after is None:
            return passages
        found = False
        advanced: list[NeckPassage] = []
        for passage in passages:
            if passage.scope == passage_after.scope:
                advanced.append(passage_after)
                found = True
            else:
                advanced.append(passage)
        if not found:
            raise CandidateStateMismatchError("candidate passage result has no authoritative parent.")
        return tuple(advanced)


def select_candidate_transaction(
    transactions: tuple[CandidateTransaction, ...],
) -> CandidateTransaction:
    """Select one accepted transaction by invariant candidate policy.

    Args:
        transactions: Nonempty accepted set rooted at one parent state.

    Returns:
        Deterministic furthest/largest/canonical winner.

    Raises:
        CandidateSelectionError: If the set is empty, foreign, duplicated, or
            spans different parents or candidate policies.
    """
    if type(transactions) is not tuple or not transactions:
        raise CandidateSelectionError("candidate selection requires one nonempty transaction tuple.")
    if any(type(transaction) is not CandidateTransaction for transaction in transactions):
        raise CandidateSelectionError("candidate selection contains a foreign transaction type.")
    parent_digests = {transaction.parent_state_digest for transaction in transactions}
    if len(parent_digests) != 1:
        raise CandidateSelectionError("candidate selection transactions must share one parent state.")
    policies = {transaction.candidate.policy for transaction in transactions}
    if len(policies) != 1:
        raise CandidateSelectionError("candidate selection transactions must share one candidate policy.")
    digests = tuple(transaction.digest for transaction in transactions)
    if len(digests) != len(set(digests)):
        raise CandidateSelectionError("candidate selection transaction identities must be unique.")
    policy = next(iter(policies))
    return policy.order_candidates(
        transactions,
        key=lambda transaction: transaction.candidate.order_key,
    )[0]
