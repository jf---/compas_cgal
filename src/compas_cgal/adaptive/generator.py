"""Cross-bound exact continuation over physical and global MAT state."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Self
from typing import TypeAlias
from typing import cast

from compas_cgal.adaptive.bootstrap import InitialCandidateEvaluator
from compas_cgal.adaptive.bootstrap import InitialCandidateTransaction
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import TraversalCandidate
from compas_cgal.adaptive.candidates import ZeroGuideLinkCandidate
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
from compas_cgal.adaptive.candidates import enumerate_zero_guide_link_candidates
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import DegenerateSegmentMotionError
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import EngagementCapInfeasibleError
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidCandidateFamilyError
from compas_cgal.adaptive.errors import InvalidCandidatePolicyError
from compas_cgal.adaptive.errors import InvalidRouteRetraceCommitError
from compas_cgal.adaptive.errors import InvalidRouteRetraceDecisionError
from compas_cgal.adaptive.errors import InvalidTraversalCommitError
from compas_cgal.adaptive.errors import NeckTooTightError
from compas_cgal.adaptive.errors import NoFeasibleCandidateError
from compas_cgal.adaptive.errors import StaleTraversalCursorError
from compas_cgal.adaptive.errors import TerminalTraversalCursorError
from compas_cgal.adaptive.errors import UnsupportedRouteRetraceError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MatZeroGuideRun
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import EffectiveCapDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import RouteNodeId
from compas_cgal.adaptive.operation import RouteRetraceDecision
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.retrace_transaction import RouteRetraceEvaluator
from compas_cgal.adaptive.retrace_transaction import RouteRetraceTransaction
from compas_cgal.adaptive.transaction import AcceptedCandidateTransaction
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.transaction import ZeroGuideLinkTransaction
from compas_cgal.adaptive.traversal import MatTraversalState

_DIGEST_SIZE = hashlib.sha256().digest_size


def _digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise InvalidTraversalCommitError(
            f"{name} must be one exact SHA-256 digest.",
        )
    return value


def _route_retrace_digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise InvalidRouteRetraceCommitError(
            f"{name} must be one exact SHA-256 digest.",
        )
    return value


def _advance_global(
    traversal: MatTraversalState,
    candidate: TraversalCandidate,
) -> MatTraversalState:
    if type(traversal) is not MatTraversalState:
        raise InvalidTraversalCommitError(
            "global parent must be one exact MAT traversal state.",
        )
    if type(candidate) not in (MiddleCurveCandidate, ZeroGuideLinkCandidate):
        raise InvalidTraversalCommitError(
            "global continuation requires one exact candidate.",
        )
    try:
        return traversal.advance(candidate)
    except (
        StaleTraversalCursorError,
        TerminalTraversalCursorError,
    ) as error:
        raise InvalidTraversalCommitError(
            "candidate does not begin at the authoritative global parent.",
        ) from error


def _validate_evaluator_authority(
    evaluator: CandidateEvaluator,
    traversal: MatTraversalState,
) -> None:
    if type(evaluator) is not CandidateEvaluator:
        raise InvalidTraversalCommitError(
            "global continuation requires one exact candidate evaluator.",
        )
    if type(traversal) is not MatTraversalState:
        raise InvalidTraversalCommitError(
            "global continuation requires one exact MAT traversal state.",
        )
    certificate = evaluator.reachable_domain.certificate
    axis = traversal.authority.axis
    if (
        axis.design_boundary != certificate.design_boundary
        or axis.holes != certificate.holes
        or axis.tool_radius != evaluator.tool_radius
        or traversal.authority.inventory.policy != evaluator.neck_policy
    ):
        raise InvalidTraversalCommitError(
            "candidate evaluator contradicts global MAT authority.",
        )


def _active_zero_guide_run(
    traversal: MatTraversalState,
    error_type: type[ValueError],
) -> MatZeroGuideRun | None:
    """Return the active edge's inventory-authenticated zero-guide proof.

    Args:
        traversal: Authoritative nonterminal MAT traversal state.
        error_type: Named public-boundary error to raise on contradiction.

    Returns:
        Exact owned proof run, or `None` for the ordinary circle family.

    Raises:
        ValueError: Through `error_type` if the active edge or projected record
            contradicts the immutable MAT proof inventory.
    """
    if type(traversal) is not MatTraversalState or traversal.active_route_index is None:
        raise error_type(
            "candidate variant selection requires one active exact MAT route.",
        )
    axis = traversal.authority.axis
    edge_id = traversal.active_cursor.route_step.edge_id
    run = axis.zero_guide_run_by_edge_id.get(edge_id)
    if run is None:
        return None
    inventory_matches = tuple(owned for owned in axis.zero_guide_inventory.runs if bytes(owned.edge_id) == bytes(edge_id))
    if type(run) is not MatZeroGuideRun or len(inventory_matches) != 1 or inventory_matches[0] != run or inventory_matches[0].native_certificate != run.native_certificate:
        raise error_type(
            "active zero-guide record contradicts its MAT proof inventory.",
        )
    return run


def _validate_candidate_variant(
    candidate: TraversalCandidate,
    zero_guide_run: MatZeroGuideRun | None,
    error_type: type[ValueError],
) -> None:
    if zero_guide_run is None:
        if type(candidate) is MiddleCurveCandidate:
            return
        raise error_type(
            "candidate contradicts its active MAT proof variant.",
        )
    if type(candidate) is not ZeroGuideLinkCandidate:
        raise error_type(
            "candidate contradicts its active MAT proof variant.",
        )
    if candidate.zero_guide_run != zero_guide_run or candidate.zero_guide_run.native_certificate != zero_guide_run.native_certificate:
        raise error_type(
            "zero-guide candidate contradicts its owned native proof bytes.",
        )


def _validate_candidate_family(
    *,
    evaluator: CandidateEvaluator,
    traversal: MatTraversalState,
    candidates: tuple[TraversalCandidate, ...],
) -> None:
    _validate_evaluator_authority(
        evaluator,
        traversal,
    )
    if traversal.active_route_index is None:
        raise InvalidCandidateFamilyError(
            "finite search requires one active global MAT cursor.",
        )
    if type(candidates) is not tuple or any(type(candidate) not in (MiddleCurveCandidate, ZeroGuideLinkCandidate) for candidate in candidates):
        raise InvalidCandidateFamilyError(
            "finite search requires one immutable exact candidate tuple.",
        )
    zero_guide_run = _active_zero_guide_run(
        traversal,
        InvalidCandidateFamilyError,
    )
    for candidate in candidates:
        _validate_candidate_variant(
            candidate,
            zero_guide_run,
            InvalidCandidateFamilyError,
        )
    if any(candidate.policy != evaluator.candidate_policy for candidate in candidates):
        raise InvalidCandidateFamilyError(
            "finite candidate family contradicts evaluator policy.",
        )
    try:
        invariant_order = evaluator.candidate_policy.order_candidates(
            candidates,
            key=lambda candidate: candidate.order_key,
        )
    except InvalidCandidatePolicyError as error:
        raise InvalidCandidateFamilyError(
            "finite candidate family identities are not unique.",
        ) from error
    if candidates != invariant_order:
        raise InvalidCandidateFamilyError(
            "finite candidate family is not in invariant order.",
        )
    for candidate in candidates:
        try:
            _advance_global(
                traversal,
                candidate,
            )
        except InvalidTraversalCommitError as error:
            raise InvalidCandidateFamilyError(
                "finite candidate family contains a foreign global cursor.",
            ) from error


def _exhaustion_summary(
    traversal: MatTraversalState,
    *,
    attempts: int,
    cap: int,
    gouge: int,
    degenerate_link: int,
) -> str:
    cursor = bytes(traversal.active_cursor.cursor_identity).hex()
    return f"finite candidate family exhausted at cursor={cursor}; attempts={attempts}; cap={cap}; gouge={gouge}; degenerate-link={degenerate_link}."


def _route_retrace_required(
    terminal: MatTraversalState,
    activated: MatTraversalState,
) -> bool:
    """Decide route transport from stable MAT node identity.

    Args:
        terminal: Exact completed-route state before activation.
        activated: Exact state produced by one `activate_next()` call.

    Returns:
        Whether completed exit and activated entry nodes are nonincident.

    Raises:
        InvalidRouteRetraceDecisionError: If the pair is not one exact
            terminal-to-active transition.
    """
    if (
        type(terminal) is not MatTraversalState
        or type(activated) is not MatTraversalState
        or terminal.active_route_index is None
        or not terminal.active_cursor.terminal
        or activated != terminal.activate_next()
        or activated.active_route_index is None
    ):
        raise InvalidRouteRetraceDecisionError(
            "route trigger requires one exact terminal-to-active transition.",
        )
    completed = terminal.active_cursor.route_step
    next_step = activated.active_cursor.route_step
    return completed.exit_node_id != next_step.entry_node_id


def _activate_completed_incident_routes(
    traversal: MatTraversalState,
) -> MatTraversalState:
    current = traversal
    while current.active_route_index is not None and current.active_cursor.terminal:
        activated = current.activate_next()
        if activated.active_route_index is None:
            return activated
        if _route_retrace_required(current, activated):
            return current
        current = activated
    return current


@dataclass(frozen=True)
class TraversalCommit:
    """Immutable atomic binding for one physical and global continuation."""

    physical_parent_digest: IdentityDigest
    traversal_before: MatTraversalState
    transaction: AcceptedCandidateTransaction
    physical_child_digest: IdentityDigest
    traversal_after: MatTraversalState

    def __post_init__(self) -> None:
        _digest(
            self.physical_parent_digest,
            "traversal commit physical parent",
        )
        _digest(
            self.physical_child_digest,
            "traversal commit physical child",
        )
        if type(self.traversal_before) is not MatTraversalState or type(self.traversal_after) is not MatTraversalState:
            raise InvalidTraversalCommitError(
                "traversal commit requires exact global parent and child states.",
            )
        if type(self.transaction) not in (CandidateTransaction, ZeroGuideLinkTransaction):
            raise InvalidTraversalCommitError(
                "traversal commit requires one exact physical transaction.",
            )
        if self.transaction.parent_state_digest != self.physical_parent_digest:
            raise InvalidTraversalCommitError(
                "transaction contradicts traversal commit physical parent.",
            )
        if self.transaction.result_state_digest != self.physical_child_digest:
            raise InvalidTraversalCommitError(
                "transaction contradicts traversal commit physical child.",
            )
        _validate_candidate_variant(
            self.transaction.candidate,
            _active_zero_guide_run(
                self.traversal_before,
                InvalidTraversalCommitError,
            ),
            InvalidTraversalCommitError,
        )
        expected_after = _advance_global(
            self.traversal_before,
            self.transaction.candidate,
        )
        if expected_after != self.traversal_after:
            raise InvalidTraversalCommitError(
                "traversal commit global child contradicts its candidate.",
            )
        changed = tuple(
            index
            for index, (before, after) in enumerate(
                zip(
                    self.traversal_before.cursors,
                    self.traversal_after.cursors,
                    strict=True,
                )
            )
            if before.canonical_bytes != after.canonical_bytes
        )
        if changed != (self.traversal_before.active_route_index,):
            raise InvalidTraversalCommitError(
                "traversal commit must advance exactly one active global cursor.",
            )
        candidate = self.transaction.candidate
        if candidate.neck_scope != self.traversal_before.neck_scope:
            raise InvalidTraversalCommitError(
                "traversal commit candidate contradicts causal global scope.",
            )
        if self.traversal_before.pending_transit is None and type(candidate.neck_scope) is not NoNeckScope:
            raise InvalidTraversalCommitError(
                "traversal commit manufactures an unauthenticated neck transit.",
            )
        if self.traversal_before.pending_transit is not None and self.transaction.passage_after is None:
            raise InvalidTraversalCommitError(
                "traversal commit omits its causal passage result.",
            )

    @classmethod
    def build(
        cls,
        *,
        physical_before: GenerationState,
        traversal_before: MatTraversalState,
        transaction: AcceptedCandidateTransaction,
        physical_after: GenerationState,
        traversal_after: MatTraversalState,
    ) -> Self:
        """Cross-bind one independently reproduced continuation.

        Args:
            physical_before: Authoritative stock/coverage parent.
            traversal_before: Authoritative global graph parent.
            transaction: Accepted circle or advancing-segment evidence.
            physical_after: Independently reproduced physical child.
            traversal_after: Candidate-derived global child.

        Returns:
            Content-addressed commit binding both state axes.

        Raises:
            InvalidTraversalCommitError: If a type, digest, cursor, witness, or
                causal transition is stale or cross-wired.
        """
        if type(physical_before) is not GenerationState or type(physical_after) is not GenerationState:
            raise InvalidTraversalCommitError(
                "traversal commit requires exact physical parent and child states.",
            )
        if type(transaction) not in (CandidateTransaction, ZeroGuideLinkTransaction):
            raise InvalidTraversalCommitError(
                "traversal commit requires one exact physical transaction.",
            )
        if transaction.parent_state_digest != physical_before.digest:
            raise InvalidTraversalCommitError(
                "transaction contradicts authoritative physical parent.",
            )
        if transaction.result_state_digest != physical_after.digest or transaction.traversal_after != physical_after.traversal:
            raise InvalidTraversalCommitError(
                "transaction contradicts authoritative physical child.",
            )
        transaction_suffix: tuple[CanonicalOperation, ...]
        if type(transaction) is CandidateTransaction:
            transaction_suffix = (
                transaction.link_witness.operation,
                transaction.circle_witness.operation,
            )
        elif type(transaction) is ZeroGuideLinkTransaction:
            transaction_suffix = (transaction.segment_witness.operation,)
        else:
            raise InvalidTraversalCommitError(
                "traversal commit received a foreign physical transaction.",
            )
        if physical_after.operations != physical_before.operations + transaction_suffix:
            raise InvalidTraversalCommitError(
                "physical child must append the exact transaction suffix.",
            )
        return cls(
            physical_before.digest,
            traversal_before,
            transaction,
            physical_after.digest,
            traversal_after,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete cross-axis commit record.

        Returns:
            Canonical CCAN bytes binding both parent and child identities.
        """
        return encode_tagged_union(
            b"traversal-commit-v1",
            encode_component_map(
                {
                    b"physical-child": encode_bytes(
                        bytes(self.physical_child_digest),
                    ),
                    b"physical-parent": encode_bytes(
                        bytes(self.physical_parent_digest),
                    ),
                    b"transaction": self.transaction.canonical_bytes,
                    b"traversal-after": self.traversal_after.canonical_bytes,
                    b"traversal-before": self.traversal_before.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
        return IdentityDigest(
            hashlib.sha256(self.canonical_bytes).digest(),
        )


def _derive_route_retrace_decision(
    *,
    physical: GenerationState,
    terminal: MatTraversalState,
    activated: MatTraversalState,
    source_commit: TraversalCommit,
) -> RouteRetraceDecision:
    """Authenticate the sole admitted source for a nonincident route switch.

    Args:
        physical: Current physical child of `source_commit`.
        terminal: Current terminal global child of `source_commit`.
        activated: Exact next-route activation derived from `terminal`.
        source_commit: Immediately preceding accepted traversal commit.

    Returns:
        Content-addressed decision naming every causal preimage.

    Raises:
        InvalidRouteRetraceDecisionError: If physical, global, witness, or
            source identities are incomplete, stale, or cross-wired.
        UnsupportedRouteRetraceError: If the exact boundary has no admitted
            no-neck full-cap advancing-segment source.
    """
    if (
        type(physical) is not GenerationState
        or type(terminal) is not MatTraversalState
        or type(activated) is not MatTraversalState
        or type(source_commit) is not TraversalCommit
    ):
        raise InvalidRouteRetraceDecisionError(
            "route retrace derivation requires exact causal preimages.",
        )
    if not _route_retrace_required(terminal, activated):
        raise InvalidRouteRetraceDecisionError(
            "incident route activation requires no physical retrace.",
        )
    if (
        source_commit.physical_child_digest != physical.digest
        or source_commit.traversal_after != terminal
    ):
        raise InvalidRouteRetraceDecisionError(
            "route retrace source is not the final physical/global commit.",
        )

    source_transaction = source_commit.transaction
    if type(source_transaction) is not ZeroGuideLinkTransaction:
        raise UnsupportedRouteRetraceError(
            "nonincident route requires a final zero-guide source.",
        )
    if (
        source_commit.physical_parent_digest
        != source_transaction.parent_state_digest
        or source_transaction.result_state_digest != physical.digest
        or source_transaction.traversal_after != physical.traversal
    ):
        raise InvalidRouteRetraceDecisionError(
            "route retrace transaction is not the final physical source.",
        )
    source_witness = source_transaction.segment_witness
    if type(source_witness) is not ReplayLateralWitness:
        raise InvalidRouteRetraceDecisionError(
            "route retrace source requires one exact segment witness.",
        )
    source_index = len(physical.operations) - 1
    source_operation = physical.operations[source_index]
    if type(source_operation) is not AdvanceSegmentOperation:
        raise UnsupportedRouteRetraceError(
            "nonincident route requires a final advancing segment.",
        )
    if (
        type(source_witness.operation) is not AdvanceSegmentOperation
        or source_witness.operation != source_operation
        or source_witness.operation_index != source_index
        or source_operation.motion.end != physical.phase_point
        or source_operation.traversal_decision
        != source_transaction.candidate.traversal_decision
        or source_operation.motion.end != source_transaction.candidate.target
    ):
        raise InvalidRouteRetraceDecisionError(
            "route retrace witness, ordinal, or endpoint is cross-wired.",
        )
    if (
        type(source_operation.neck_scope) is not NoNeckScope
        or type(source_operation.effective_cap_decision) is not FullCapDecision
        or source_operation.traversal_decision.makes_cursor_terminal is not True
        or source_transaction.passage_after is not None
        or type(terminal.neck_scope) is not NoNeckScope
    ):
        raise UnsupportedRouteRetraceError(
            "route retrace source must be no-neck, full-cap, and route-terminal.",
        )

    return RouteRetraceDecision.build(
        completed_route_index=terminal.active_route_index,
        activated_route_index=activated.active_route_index,
        completed_exit_node_id=RouteNodeId(
            bytes(terminal.active_cursor.route_step.exit_node_id),
        ),
        activated_entry_node_id=RouteNodeId(
            bytes(activated.active_cursor.route_step.entry_node_id),
        ),
        terminal_traversal_digest=terminal.digest,
        activated_traversal_digest=activated.digest,
        source_commit_digest=source_commit.digest,
        source_transaction_digest=source_transaction.digest,
        source_operation_index=source_index,
        source_operation_digest=IdentityDigest(
            hashlib.sha256(source_operation.canonical_bytes).digest(),
        ),
    )


@dataclass(frozen=True)
class RouteRetraceCommit:
    """Atomic physical return and exact global route activation."""

    physical_parent_digest: IdentityDigest
    traversal_before: MatTraversalState
    transaction: RouteRetraceTransaction
    physical_child_digest: IdentityDigest
    traversal_after: MatTraversalState
    source_commit_digest: IdentityDigest

    def __post_init__(self) -> None:
        if type(self) is not RouteRetraceCommit:
            raise InvalidRouteRetraceCommitError(
                "route retrace commit must use the exact owned type.",
            )
        _route_retrace_digest(
            self.physical_parent_digest,
            "route retrace commit physical parent",
        )
        _route_retrace_digest(
            self.physical_child_digest,
            "route retrace commit physical child",
        )
        _route_retrace_digest(
            self.source_commit_digest,
            "route retrace commit source commit",
        )
        if (
            type(self.traversal_before) is not MatTraversalState
            or type(self.traversal_after) is not MatTraversalState
        ):
            raise InvalidRouteRetraceCommitError(
                "route retrace commit requires exact global parent and child states.",
            )
        if type(self.transaction) is not RouteRetraceTransaction:
            raise InvalidRouteRetraceCommitError(
                "route retrace commit requires one exact physical transaction.",
            )
        RouteRetraceTransaction.validate(self.transaction)
        if (
            self.transaction.parent_state_digest
            != self.physical_parent_digest
            or self.transaction.result_state_digest
            != self.physical_child_digest
        ):
            raise InvalidRouteRetraceCommitError(
                "route retrace transaction contradicts physical commit lineage.",
            )
        activated = self.traversal_before.activate_next()
        if (
            activated != self.traversal_after
            or self.traversal_after.active_route_index is None
            or not _route_retrace_required(
                self.traversal_before,
                self.traversal_after,
            )
        ):
            raise InvalidRouteRetraceCommitError(
                "route retrace commit requires one nonincident route activation.",
            )
        decision = self.transaction.decision
        if (
            decision.terminal_traversal_digest
            != self.traversal_before.digest
            or decision.activated_traversal_digest
            != self.traversal_after.digest
            or decision.completed_route_index
            != self.traversal_before.active_route_index
            or decision.activated_route_index
            != self.traversal_after.active_route_index
            or bytes(decision.completed_exit_node_id)
            != bytes(
                self.traversal_before.active_cursor.route_step.exit_node_id,
            )
            or bytes(decision.activated_entry_node_id)
            != bytes(
                self.traversal_after.active_cursor.route_step.entry_node_id,
            )
            or decision.source_commit_digest != self.source_commit_digest
        ):
            raise InvalidRouteRetraceCommitError(
                "route retrace decision contradicts its cross-axis commit.",
            )

    @classmethod
    def build(
        cls,
        *,
        physical_before: GenerationState,
        traversal_before: MatTraversalState,
        source_commit: TraversalCommit,
        transaction: RouteRetraceTransaction,
        physical_after: GenerationState,
        traversal_after: MatTraversalState,
    ) -> Self:
        """Cross-bind independently committed retrace and route activation.

        Args:
            physical_before: Authoritative physical source child.
            traversal_before: Terminal global source child.
            source_commit: Immediately preceding traversal commit.
            transaction: Independently reproducible physical retrace proof.
            physical_after: Independently committed physical retrace child.
            traversal_after: Exact activated next-route state.

        Returns:
            Content-addressed cross-axis retrace commit.

        Raises:
            InvalidRouteRetraceCommitError: If physical and global lineages,
                the source decision, or held state disagree.
        """
        if (
            type(physical_before) is not GenerationState
            or type(physical_after) is not GenerationState
            or type(traversal_before) is not MatTraversalState
            or type(traversal_after) is not MatTraversalState
            or type(source_commit) is not TraversalCommit
            or type(transaction) is not RouteRetraceTransaction
        ):
            raise InvalidRouteRetraceCommitError(
                "route retrace commit requires exact causal preimages.",
            )
        RouteRetraceTransaction.validate(transaction)
        decision = _derive_route_retrace_decision(
            physical=physical_before,
            terminal=traversal_before,
            activated=traversal_after,
            source_commit=source_commit,
        )
        if transaction.decision != decision:
            raise InvalidRouteRetraceCommitError(
                "retrace transaction contradicts route activation.",
            )
        if (
            transaction.parent_state_digest != physical_before.digest
            or transaction.result_state_digest != physical_after.digest
            or physical_after.operations
            != physical_before.operations
            + (transaction.segment_witness.operation,)
            or physical_after.phase_point
            != transaction.segment_witness.operation.motion.end
            or physical_after.traversal != physical_before.traversal
            or physical_after.passages != physical_before.passages
        ):
            raise InvalidRouteRetraceCommitError(
                "retrace commit breaks physical or held-state lineage.",
            )
        return cls(
            physical_before.digest,
            traversal_before,
            transaction,
            physical_after.digest,
            traversal_after,
            source_commit.digest,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete physical/global retrace commit record."""
        self.__post_init__()
        return encode_tagged_union(
            b"route-retrace-commit-v1",
            encode_component_map(
                {
                    b"physical-child": encode_bytes(
                        bytes(self.physical_child_digest),
                    ),
                    b"physical-parent": encode_bytes(
                        bytes(self.physical_parent_digest),
                    ),
                    b"source-commit": encode_bytes(
                        bytes(self.source_commit_digest),
                    ),
                    b"transaction": self.transaction.canonical_bytes,
                    b"traversal-after": self.traversal_after.canonical_bytes,
                    b"traversal-before": self.traversal_before.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
        return IdentityDigest(
            hashlib.sha256(self.canonical_bytes).digest(),
        )


ContinuationCommit: TypeAlias = TraversalCommit | RouteRetraceCommit


@dataclass(frozen=True)
class GenerationContinuation:
    """Content-addressed launch-rooted traversal prefix before coverage seal."""

    launch_transaction: InitialCandidateTransaction
    physical: GenerationState
    traversal: MatTraversalState
    commits: tuple[ContinuationCommit, ...]

    def __post_init__(self) -> None:
        if type(self.launch_transaction) is not InitialCandidateTransaction:
            raise InvalidTraversalCommitError(
                "generation continuation requires one exact launch transaction.",
            )
        if type(self.physical) is not GenerationState or type(self.traversal) is not MatTraversalState:
            raise InvalidTraversalCommitError(
                "generation continuation requires exact physical and global states.",
            )
        if type(self.commits) is not tuple or any(
            type(commit) not in (TraversalCommit, RouteRetraceCommit)
            for commit in self.commits
        ):
            raise InvalidTraversalCommitError(
                "generation continuation requires one immutable commit tuple.",
            )
        expected_physical_digest = self.launch_transaction.result_state_digest
        expected_traversal = _activate_completed_incident_routes(
            self.launch_transaction.traversal_after,
        )
        preceding_traversal_commit_digest: IdentityDigest | None = None
        for commit in self.commits:
            retrace_pending = (
                expected_traversal.active_route_index is not None
                and expected_traversal.active_cursor.terminal
            )
            if type(commit) is TraversalCommit:
                if retrace_pending:
                    raise InvalidRouteRetraceCommitError(
                        "nonincident route requires an adjacent retrace commit.",
                    )
                if commit.physical_parent_digest != expected_physical_digest:
                    raise InvalidTraversalCommitError(
                        "generation continuation breaks physical commit lineage.",
                    )
                if commit.traversal_before != expected_traversal:
                    raise InvalidTraversalCommitError(
                        "generation continuation breaks global commit lineage.",
                    )
                expected_physical_digest = commit.physical_child_digest
                expected_traversal = _activate_completed_incident_routes(
                    commit.traversal_after,
                )
                preceding_traversal_commit_digest = commit.digest
                continue
            if not retrace_pending:
                raise InvalidRouteRetraceCommitError(
                    "route retrace commit is unnecessary or out of order.",
                )
            if preceding_traversal_commit_digest is None:
                raise InvalidRouteRetraceCommitError(
                    "route retrace commit has no adjacent traversal source.",
                )
            if commit.physical_parent_digest != expected_physical_digest:
                raise InvalidRouteRetraceCommitError(
                    "generation continuation breaks retrace physical lineage.",
                )
            if commit.traversal_before != expected_traversal:
                raise InvalidRouteRetraceCommitError(
                    "generation continuation breaks retrace global lineage.",
                )
            if (
                commit.source_commit_digest
                != preceding_traversal_commit_digest
            ):
                raise InvalidRouteRetraceCommitError(
                    "route retrace does not name the immediately preceding commit.",
                )
            expected_physical_digest = commit.physical_child_digest
            expected_traversal = _activate_completed_incident_routes(
                commit.traversal_after,
            )
        if (
            expected_traversal.active_route_index is not None
            and expected_traversal.active_cursor.terminal
        ):
            raise InvalidRouteRetraceCommitError(
                "generation continuation ends before its required route retrace.",
            )
        if self.physical.digest != expected_physical_digest:
            raise InvalidTraversalCommitError(
                "generation continuation physical child contradicts commit lineage.",
            )
        if self.traversal != expected_traversal:
            raise InvalidTraversalCommitError(
                "generation continuation global child contradicts commit lineage.",
            )

    @classmethod
    def build(
        cls,
        *,
        launch_transaction: InitialCandidateTransaction,
        physical: GenerationState,
        traversal: MatTraversalState,
        commits: tuple[ContinuationCommit, ...],
    ) -> Self:
        """Build and validate one launch-rooted traversal prefix.

        Args:
            launch_transaction: Independently replayable entry-circle root.
            physical: Current stock/coverage state.
            traversal: Current normalized global MAT state.
            commits: Ordered post-launch cross-axis commits.

        Returns:
            Immutable content-addressed continuation.

        Raises:
            InvalidTraversalCommitError: If either lineage is incomplete,
                reordered, stale, or cross-wired.
        """
        return cls(
            launch_transaction,
            physical,
            traversal,
            commits,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return complete launch, commit, and current-state evidence."""
        return encode_tagged_union(
            b"generation-continuation-v1",
            encode_component_map(
                {
                    b"commits": encode_sequence(
                        tuple(commit.canonical_bytes for commit in self.commits),
                    ),
                    b"launch": self.launch_transaction.canonical_bytes,
                    b"physical": encode_bytes(bytes(self.physical.digest)),
                    b"traversal": self.traversal.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
        return IdentityDigest(
            hashlib.sha256(self.canonical_bytes).digest(),
        )


def evaluate_traversal_candidate(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
    candidate: TraversalCandidate,
) -> AcceptedCandidateTransaction:
    """Evaluate one globally authenticated candidate through Task 12.

    Args:
        evaluator: Exact invariant proof authority.
        physical: Current stock/coverage parent.
        traversal: Current global MAT parent.
        candidate: Candidate beginning at `traversal.active_cursor`.

    Returns:
        Accepted circle or advancing-segment transaction.

    Raises:
        InvalidTraversalCommitError: If global or evaluator authority differs.
        GougeContainmentError: If either cutter sweep leaves the pocket.
        EngagementCapExceededError: If either exact motion exceeds cap.
        UnresolvedMotionEventError: If exact certification is incomplete.
    """
    if type(physical) is not GenerationState:
        raise InvalidTraversalCommitError(
            "global continuation requires one exact physical parent.",
        )
    _validate_evaluator_authority(
        evaluator,
        traversal,
    )
    _validate_candidate_variant(
        candidate,
        _active_zero_guide_run(
            traversal,
            InvalidTraversalCommitError,
        ),
        InvalidTraversalCommitError,
    )
    _advance_global(
        traversal,
        candidate,
    )
    local_cursor = TraversalCursorState.before(
        candidate.traversal_decision,
    )
    if type(candidate) is MiddleCurveCandidate:
        return evaluator.evaluate_from_cursor(
            physical,
            local_cursor,
            candidate,
        )
    if type(candidate) is ZeroGuideLinkCandidate:
        return evaluator.evaluate_zero_guide_from_cursor(
            physical,
            local_cursor,
            candidate,
        )
    raise InvalidTraversalCommitError(
        "global continuation received a foreign candidate variant.",
    )


def evaluate_first_feasible_candidate(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
    candidates: tuple[TraversalCandidate, ...],
) -> AcceptedCandidateTransaction:
    """Evaluate one materialized family in invariant candidate order.

    Args:
        evaluator: Exact invariant proof authority.
        physical: Authoritative stock/coverage parent.
        traversal: Authoritative global MAT parent.
        candidates: Complete finite family for the active forward span.

    Returns:
        First accepted transaction in invariant policy order.

    Raises:
        InvalidCandidateFamilyError: If the family is foreign, duplicated, or
            not already in invariant order.
        EngagementCapInfeasibleError: If every non-neck trial exceeds cap.
        NeckTooTightError: If every causal neck trial exceeds cap.
        NoFeasibleCandidateError: If the exact rejection set is mixed or empty.
        UnresolvedMotionEventError: If any event proof is incomplete.
    """
    if type(physical) is not GenerationState:
        raise InvalidCandidateFamilyError(
            "finite search requires one exact physical parent.",
        )
    _validate_candidate_family(
        evaluator=evaluator,
        traversal=traversal,
        candidates=candidates,
    )
    cap_count = 0
    gouge_count = 0
    degenerate_link_count = 0
    for candidate in candidates:
        try:
            transaction = evaluate_traversal_candidate(
                evaluator=evaluator,
                physical=physical,
                traversal=traversal,
                candidate=candidate,
            )
        except EngagementCapExceededError:
            cap_count += 1
            continue
        except GougeContainmentError:
            gouge_count += 1
            continue
        except DegenerateSegmentMotionError:
            degenerate_link_count += 1
            continue
        if transaction.candidate != candidate:
            raise InvalidTraversalCommitError(
                "candidate proof engine returned a cross-wired transaction.",
            )
        return transaction

    attempts = len(candidates)
    summary = _exhaustion_summary(
        traversal,
        attempts=attempts,
        cap=cap_count,
        gouge=gouge_count,
        degenerate_link=degenerate_link_count,
    )
    if attempts > 0 and cap_count == attempts:
        if type(traversal.neck_scope) is OrientedNeckScope:
            raise NeckTooTightError(summary)
        raise EngagementCapInfeasibleError(summary)
    raise NoFeasibleCandidateError(summary)


def commit_traversal_candidate(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
    transaction: AcceptedCandidateTransaction,
) -> tuple[GenerationState, MatTraversalState, TraversalCommit]:
    """Independently commit one candidate on both state axes.

    Args:
        evaluator: Same exact proof authority used during evaluation.
        physical: Authoritative physical parent.
        traversal: Authoritative global parent.
        transaction: Previously accepted Task 12 evidence.

    Returns:
        Physical child, global child, and their atomic commit.

    Raises:
        InvalidTraversalCommitError: If either parent or authority is stale.
        InvalidCandidateTransactionError: If independent replay differs.
    """
    if type(transaction) not in (CandidateTransaction, ZeroGuideLinkTransaction):
        raise InvalidTraversalCommitError(
            "global continuation commit requires one exact transaction.",
        )
    _validate_evaluator_authority(
        evaluator,
        traversal,
    )
    _validate_candidate_variant(
        transaction.candidate,
        _active_zero_guide_run(
            traversal,
            InvalidTraversalCommitError,
        ),
        InvalidTraversalCommitError,
    )
    traversal_after = _advance_global(
        traversal,
        transaction.candidate,
    )
    local_cursor = TraversalCursorState.before(
        transaction.candidate.traversal_decision,
    )
    if type(transaction) is CandidateTransaction:
        physical_after = evaluator.commit_from_cursor(
            physical,
            local_cursor,
            transaction,
        )
    elif type(transaction) is ZeroGuideLinkTransaction:
        physical_after = evaluator.commit_zero_guide_from_cursor(
            physical,
            local_cursor,
            transaction,
        )
    else:
        raise InvalidTraversalCommitError(
            "global continuation commit received a foreign transaction variant.",
        )
    commit = TraversalCommit.build(
        physical_before=physical,
        traversal_before=traversal,
        transaction=transaction,
        physical_after=physical_after,
        traversal_after=traversal_after,
    )
    return physical_after, traversal_after, commit


def _active_effective_cap(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
) -> EffectiveCapDecision:
    transit = traversal.pending_transit
    if transit is None:
        return FullCapDecision.build(
            user_cap=evaluator.user_cap,
            effective_cap=evaluator.user_cap,
        )
    passage = transit.passage(physical.passages)
    return passage.propose_cap_decision(evaluator.neck_policy)


def _active_forward_limits(
    traversal: MatTraversalState,
) -> tuple[MatSample, ...]:
    if type(traversal) is not MatTraversalState:
        raise InvalidCandidateFamilyError(
            "forward window requires one exact MAT traversal state.",
        )
    if traversal.active_route_index is None:
        raise InvalidCandidateFamilyError(
            "terminal traversal has no active forward window.",
        )
    active = traversal.active_cursor
    if active.terminal:
        raise InvalidCandidateFamilyError(
            "terminal active cursor must be activated before materialization.",
        )
    sample_index = traversal.authority.sample_index
    samples = sample_index.samples_by_edge[active.route_step.edge_id]
    initial = sample_index.sample_by_cursor_id[active.route_step.initial_cursor_id]
    terminal = sample_index.sample_by_cursor_id[active.route_step.terminal_cursor_id]
    ordinal_step = 1 if terminal.ordinal_on_edge > initial.ordinal_on_edge else -1
    cursor = active.cursor
    if type(cursor) is MatSample:
        if ordinal_step == 1:
            eligible = tuple(sample for sample in samples if sample.ordinal_on_edge > cursor.ordinal_on_edge)
        else:
            eligible = tuple(sample for sample in reversed(samples) if sample.ordinal_on_edge < cursor.ordinal_on_edge)
    elif type(cursor) is DerivedCandidateCursor:
        if cursor.ordinal_step != ordinal_step:
            raise InvalidCandidateFamilyError(
                "derived cursor direction contradicts its global route.",
            )
        if ordinal_step == 1:
            eligible = tuple(sample for sample in samples if sample.ordinal_on_edge >= cursor.next_limit_ordinal)
        else:
            eligible = tuple(sample for sample in reversed(samples) if sample.ordinal_on_edge <= cursor.next_limit_ordinal)
    else:
        raise InvalidCandidateFamilyError(
            "active forward window requires native or derived cursor lineage.",
        )
    limits = eligible[: traversal.authority.policy.forward_window]
    if not limits:
        raise InvalidCandidateFamilyError(
            "nonterminal active cursor has no directed native limit.",
        )
    return limits


def materialize_active_candidate_family(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
) -> tuple[TraversalCandidate, ...]:
    """Materialize each active forward span exactly once.

    Args:
        evaluator: Exact continuation proof authority.
        physical: Current stock/coverage and passage state.
        traversal: Current normalized global MAT state.

    Returns:
        Complete family in invariant candidate-policy order.

    Raises:
        InvalidCandidateFamilyError: If authority, direction, or family
            identity is foreign or structurally incomplete.
        InvalidCausalNeckTransitError: If a pending transit has no unique
            physical passage state.
        TerminalNeckPassageError: If the oriented passage is already terminal.
    """
    if type(physical) is not GenerationState:
        raise InvalidCandidateFamilyError(
            "candidate materialization requires one exact physical parent.",
        )
    _validate_evaluator_authority(
        evaluator,
        traversal,
    )
    limits = _active_forward_limits(traversal)
    active = traversal.active_cursor
    effective_cap = _active_effective_cap(
        evaluator=evaluator,
        physical=physical,
        traversal=traversal,
    )
    zero_guide_run = _active_zero_guide_run(
        traversal,
        InvalidCandidateFamilyError,
    )
    candidates: list[TraversalCandidate] = []
    cursor_before = active.cursor
    if type(cursor_before) not in (MatSample, DerivedCandidateCursor):
        raise InvalidCandidateFamilyError(
            "active candidate family cannot begin at exhausted cursor lineage.",
        )
    for limit in limits:
        span = MiddleCurveSpan.build(
            axis=traversal.authority.axis,
            cursor_before=cast(
                MatSample | DerivedCandidateCursor,
                cursor_before,
            ),
            cursor_limit=limit,
        )
        if zero_guide_run is None:
            candidates.extend(
                enumerate_middle_curve_candidates(
                    span=span,
                    policy=evaluator.candidate_policy,
                    circle_orientation=(
                        evaluator.cut_direction_policy.circle_orientation(
                            evaluator.material_side,
                        )
                    ),
                    neck_scope=traversal.neck_scope,
                    effective_cap_decision=effective_cap,
                    makes_cursor_terminal_at_limit=(limit == active.terminal_cursor),
                )
            )
            continue
        if bytes(span.edge.identity) != bytes(zero_guide_run.edge_id):
            raise InvalidCandidateFamilyError(
                "active zero-guide run contradicts its directed forward span.",
            )
        candidates.extend(
            enumerate_zero_guide_link_candidates(
                span=span,
                policy=evaluator.candidate_policy,
                neck_scope=traversal.neck_scope,
                effective_cap_decision=effective_cap,
                makes_cursor_terminal_at_limit=(limit == active.terminal_cursor),
            )
        )
    try:
        family = evaluator.candidate_policy.order_candidates(
            tuple(candidates),
            key=lambda candidate: candidate.order_key,
        )
    except InvalidCandidatePolicyError as error:
        raise InvalidCandidateFamilyError(
            "materialized forward spans contain duplicate candidate identities.",
        ) from error
    return family


def advance_active_candidate_family(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
) -> tuple[GenerationState, MatTraversalState, TraversalCommit]:
    """Select and independently commit one active forward-window winner.

    Args:
        evaluator: Exact continuation proof authority.
        physical: Authoritative physical parent.
        traversal: Authoritative normalized global parent.

    Returns:
        Physical child, unactivated global child, and atomic cross-axis commit.
    """
    candidates = materialize_active_candidate_family(
        evaluator=evaluator,
        physical=physical,
        traversal=traversal,
    )
    transaction = evaluate_first_feasible_candidate(
        evaluator=evaluator,
        physical=physical,
        traversal=traversal,
        candidates=candidates,
    )
    return commit_traversal_candidate(
        evaluator=evaluator,
        physical=physical,
        traversal=traversal,
        transaction=transaction,
    )


def _validate_pipeline_authority(
    *,
    initial_evaluator: InitialCandidateEvaluator,
    evaluator: CandidateEvaluator,
) -> None:
    if type(initial_evaluator) is not InitialCandidateEvaluator:
        raise InvalidTraversalCommitError(
            "generation requires one exact initial candidate evaluator.",
        )
    if type(evaluator) is not CandidateEvaluator:
        raise InvalidTraversalCommitError(
            "generation requires one exact continuation evaluator.",
        )
    identity = initial_evaluator.input_identity
    if (
        evaluator.reachable_domain.certificate.digest != identity.reachable_domain_digest
        or evaluator.tool_radius != identity.tool_radius
        or evaluator.user_cap != identity.user_cap
        or evaluator.candidate_policy != identity.candidate_policy
        or evaluator.neck_policy != identity.neck_policy
        or evaluator.depletion_policy != identity.depletion_policy
        or evaluator.cut_direction_policy != identity.cut_direction_policy
        or evaluator.cut_z != identity.cut_plane.cut_z
        or evaluator.material_side is not initial_evaluator.material_side
    ):
        raise InvalidTraversalCommitError(
            "initial and continuation evaluators do not share one input authority.",
        )


def generate_exact_adaptive_continuation(
    *,
    initial_evaluator: InitialCandidateEvaluator,
    evaluator: CandidateEvaluator,
    seeded_traversal: MatTraversalState,
    launch_transaction: InitialCandidateTransaction,
) -> GenerationContinuation:
    """Generate one terminal MAT traversal prefix without sealing coverage.

    Args:
        initial_evaluator: Exact entry-circle proof authority.
        evaluator: Exact post-launch link-and-circle proof authority.
        seeded_traversal: Untouched global MAT root.
        launch_transaction: Previously accepted entry launch evidence.

    Returns:
        Terminal global traversal and its content-addressed physical lineage.

    Raises:
        InvalidTraversalCommitError: If launch and continuation authorities
            differ.
        NoFeasibleCandidateError: If a finite family has mixed exact rejection.
        EngagementCapInfeasibleError: If a non-neck family is cap-infeasible.
        NeckTooTightError: If a causal neck family is cap-infeasible.
        UnresolvedMotionEventError: If exact event certification is incomplete.

    Note:
        This artifact proves traversal exhaustion, not complete pocket
        coverage. `GenerationResult` adds exact residual emptiness and fresh
        replay in the terminal-seal stage.
    """
    _validate_pipeline_authority(
        initial_evaluator=initial_evaluator,
        evaluator=evaluator,
    )
    retrace_evaluator = RouteRetraceEvaluator.build(
        evaluator=evaluator,
    )
    physical, traversal = initial_evaluator.commit(
        seeded_traversal,
        launch_transaction,
    )
    traversal = _activate_completed_incident_routes(traversal)
    commits: list[ContinuationCommit] = []
    while traversal.active_route_index is not None:
        if traversal.active_cursor.terminal:
            raise UnsupportedRouteRetraceError(
                "nonincident launch boundary has no traversal source commit.",
            )
        physical, traversal, commit = advance_active_candidate_family(
            evaluator=evaluator,
            physical=physical,
            traversal=traversal,
        )
        commits.append(commit)
        traversal = _activate_completed_incident_routes(traversal)
        if (
            traversal.active_route_index is not None
            and traversal.active_cursor.terminal
        ):
            activated = traversal.activate_next()
            decision = _derive_route_retrace_decision(
                physical=physical,
                terminal=traversal,
                activated=activated,
                source_commit=commit,
            )
            retrace_transaction = retrace_evaluator.evaluate(
                physical,
                decision,
            )
            physical_after = retrace_evaluator.commit(
                physical,
                retrace_transaction,
            )
            retrace_commit = RouteRetraceCommit.build(
                physical_before=physical,
                traversal_before=traversal,
                source_commit=commit,
                transaction=retrace_transaction,
                physical_after=physical_after,
                traversal_after=activated,
            )
            physical = physical_after
            traversal = activated
            commits.append(retrace_commit)
    traversal.require_terminal()
    return GenerationContinuation.build(
        launch_transaction=launch_transaction,
        physical=physical,
        traversal=traversal,
        commits=tuple(commits),
    )
