"""Cross-bound exact continuation over physical and global MAT state."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Self
from typing import cast

from compas_cgal.adaptive.bootstrap import InitialCandidateEvaluator
from compas_cgal.adaptive.bootstrap import InitialCandidateTransaction
from compas_cgal.adaptive.candidates import DerivedCandidateCursor
from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.candidates import MiddleCurveSpan
from compas_cgal.adaptive.candidates import enumerate_middle_curve_candidates
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
from compas_cgal.adaptive.errors import InvalidTraversalCommitError
from compas_cgal.adaptive.errors import NeckTooTightError
from compas_cgal.adaptive.errors import NoFeasibleCandidateError
from compas_cgal.adaptive.errors import StaleTraversalCursorError
from compas_cgal.adaptive.errors import TerminalTraversalCursorError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.operation import EffectiveCapDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.transaction import CandidateEvaluator
from compas_cgal.adaptive.transaction import CandidateTransaction
from compas_cgal.adaptive.traversal import MatTraversalState

_DIGEST_SIZE = hashlib.sha256().digest_size


def _digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise InvalidTraversalCommitError(
            f"{name} must be one exact SHA-256 digest.",
        )
    return value


def _advance_global(
    traversal: MatTraversalState,
    candidate: MiddleCurveCandidate,
) -> MatTraversalState:
    if type(traversal) is not MatTraversalState:
        raise InvalidTraversalCommitError(
            "global parent must be one exact MAT traversal state.",
        )
    if type(candidate) is not MiddleCurveCandidate:
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


def _validate_candidate_family(
    *,
    evaluator: CandidateEvaluator,
    traversal: MatTraversalState,
    candidates: tuple[MiddleCurveCandidate, ...],
) -> None:
    _validate_evaluator_authority(
        evaluator,
        traversal,
    )
    if traversal.active_route_index is None:
        raise InvalidCandidateFamilyError(
            "finite search requires one active global MAT cursor.",
        )
    if type(candidates) is not tuple or any(type(candidate) is not MiddleCurveCandidate for candidate in candidates):
        raise InvalidCandidateFamilyError(
            "finite search requires one immutable exact candidate tuple.",
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


def _activate_completed_routes(
    traversal: MatTraversalState,
) -> MatTraversalState:
    current = traversal
    while current.active_route_index is not None and current.active_cursor.terminal:
        current = current.activate_next()
    return current


@dataclass(frozen=True)
class TraversalCommit:
    """Immutable atomic binding for one physical and global continuation."""

    physical_parent_digest: IdentityDigest
    traversal_before: MatTraversalState
    transaction: CandidateTransaction
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
        if type(self.transaction) is not CandidateTransaction:
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
        transaction: CandidateTransaction,
        physical_after: GenerationState,
        traversal_after: MatTraversalState,
    ) -> Self:
        """Cross-bind one independently reproduced continuation.

        Args:
            physical_before: Authoritative stock/coverage parent.
            traversal_before: Authoritative global graph parent.
            transaction: Accepted Task 12 link-and-circle evidence.
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
        if type(transaction) is not CandidateTransaction:
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
        if physical_after.operations[-2:] != (
            transaction.link_witness.operation,
            transaction.circle_witness.operation,
        ):
            raise InvalidTraversalCommitError(
                "physical child omits the transaction link-and-circle suffix.",
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


@dataclass(frozen=True)
class GenerationContinuation:
    """Content-addressed launch-rooted traversal prefix before coverage seal."""

    launch_transaction: InitialCandidateTransaction
    physical: GenerationState
    traversal: MatTraversalState
    commits: tuple[TraversalCommit, ...]

    def __post_init__(self) -> None:
        if type(self.launch_transaction) is not InitialCandidateTransaction:
            raise InvalidTraversalCommitError(
                "generation continuation requires one exact launch transaction.",
            )
        if type(self.physical) is not GenerationState or type(self.traversal) is not MatTraversalState:
            raise InvalidTraversalCommitError(
                "generation continuation requires exact physical and global states.",
            )
        if type(self.commits) is not tuple or any(type(commit) is not TraversalCommit for commit in self.commits):
            raise InvalidTraversalCommitError(
                "generation continuation requires one immutable commit tuple.",
            )
        expected_physical_digest = self.launch_transaction.result_state_digest
        expected_traversal = _activate_completed_routes(
            self.launch_transaction.traversal_after,
        )
        for commit in self.commits:
            if commit.physical_parent_digest != expected_physical_digest:
                raise InvalidTraversalCommitError(
                    "generation continuation breaks physical commit lineage.",
                )
            if commit.traversal_before != expected_traversal:
                raise InvalidTraversalCommitError(
                    "generation continuation breaks global commit lineage.",
                )
            expected_physical_digest = commit.physical_child_digest
            expected_traversal = _activate_completed_routes(
                commit.traversal_after,
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
        commits: tuple[TraversalCommit, ...],
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
    candidate: MiddleCurveCandidate,
) -> CandidateTransaction:
    """Evaluate one globally authenticated candidate through Task 12.

    Args:
        evaluator: Exact invariant proof authority.
        physical: Current stock/coverage parent.
        traversal: Current global MAT parent.
        candidate: Candidate beginning at `traversal.active_cursor`.

    Returns:
        Accepted link-and-circle transaction.

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
    _advance_global(
        traversal,
        candidate,
    )
    return evaluator.evaluate_from_cursor(
        physical,
        TraversalCursorState.before(candidate.traversal_decision),
        candidate,
    )


def evaluate_first_feasible_candidate(
    *,
    evaluator: CandidateEvaluator,
    physical: GenerationState,
    traversal: MatTraversalState,
    candidates: tuple[MiddleCurveCandidate, ...],
) -> CandidateTransaction:
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
    transaction: CandidateTransaction,
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
    if type(transaction) is not CandidateTransaction:
        raise InvalidTraversalCommitError(
            "global continuation commit requires one exact transaction.",
        )
    _validate_evaluator_authority(
        evaluator,
        traversal,
    )
    traversal_after = _advance_global(
        traversal,
        transaction.candidate,
    )
    local_cursor = TraversalCursorState.before(
        transaction.candidate.traversal_decision,
    )
    physical_after = evaluator.commit_from_cursor(
        physical,
        local_cursor,
        transaction,
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
) -> tuple[MiddleCurveCandidate, ...]:
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
    candidates: list[MiddleCurveCandidate] = []
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
    try:
        return evaluator.candidate_policy.order_candidates(
            tuple(candidates),
            key=lambda candidate: candidate.order_key,
        )
    except InvalidCandidatePolicyError as error:
        raise InvalidCandidateFamilyError(
            "materialized forward spans contain duplicate candidate identities.",
        ) from error


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
    physical, traversal = initial_evaluator.commit(
        seeded_traversal,
        launch_transaction,
    )
    traversal = _activate_completed_routes(traversal)
    commits: list[TraversalCommit] = []
    while traversal.active_route_index is not None:
        physical, traversal, commit = advance_active_candidate_family(
            evaluator=evaluator,
            physical=physical,
            traversal=traversal,
        )
        commits.append(commit)
        traversal = _activate_completed_routes(traversal)
    traversal.require_terminal()
    return GenerationContinuation.build(
        launch_transaction=launch_transaction,
        physical=physical,
        traversal=traversal,
        commits=tuple(commits),
    )
