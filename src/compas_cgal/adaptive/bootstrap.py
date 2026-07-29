"""Atomic exact evaluation and commit for the first entry-launched circle."""

import hashlib
from dataclasses import dataclass
from dataclasses import field
from typing import Self

from compas.geometry import Polygon  # type: ignore[import-untyped]

from compas_cgal.adaptive.candidates import MiddleCurveCandidate
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import CircleContainmentCertificate
from compas_cgal.adaptive.containment import CircleInEntryCertificate
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.coverage import CoverageCertificate
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.errors import InitialCandidatePhaseError
from compas_cgal.adaptive.errors import InitialCandidateStateMismatchError
from compas_cgal.adaptive.errors import InvalidCoverageCertificateError
from compas_cgal.adaptive.errors import InvalidInitialCandidateEvaluatorError
from compas_cgal.adaptive.errors import InvalidInitialCandidateTransactionError
from compas_cgal.adaptive.errors import StaleInitialCandidateTransactionError
from compas_cgal.adaptive.errors import StaleTraversalCursorError
from compas_cgal.adaptive.errors import TerminalTraversalCursorError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.replay_trace import ReplayLateralWitness
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.traversal import MatTraversalState
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType

_DIGEST_SIZE = hashlib.sha256().digest_size


def _digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise InvalidInitialCandidateTransactionError(
            f"{name} must be one exact SHA-256 digest.",
        )
    return value


def _stock_lineage_digest(
    witnesses: tuple[DepletionWitness, ...],
) -> bytes:
    return hashlib.sha256(
        encode_tagged_union(
            b"exact-stock-lineage-v1",
            encode_sequence(
                tuple(bytes(witness.digest) for witness in witnesses),
            ),
        )
    ).digest()


def _phase_point(
    candidate: MiddleCurveCandidate,
) -> Point2[WorldXY]:
    motion = candidate.motion
    return Point2[WorldXY].build(
        motion.center.x + motion.phase_vector.x,
        motion.center.y + motion.phase_vector.y,
    )


def _polygon(
    vertices: tuple[Point2[WorldXY], ...],
) -> Polygon:
    return Polygon(
        [[vertex.x, vertex.y, 0.0] for vertex in vertices],
    )


def _pristine_stock(identity: InputIdentity) -> Stock2Area:
    return Stock2Area.build(
        Stock(
            _polygon(identity.design_boundary.vertices),
            [_polygon(hole.vertices) for hole in identity.holes],
        )
    )


def _initial_passages(
    traversal: MatTraversalState,
) -> tuple[NeckPassage, ...]:
    return tuple(passage for neck in traversal.authority.inventory.necks for passage in (neck.forward, neck.reverse))


@dataclass(frozen=True)
class InitialCandidateTransaction:
    """Immutable proof bundle for one entry-launched full circle."""

    input_identity_digest: IdentityDigest
    traversal_before: MatTraversalState
    candidate: MiddleCurveCandidate
    pristine_stock_boundary_digest: bytes
    entry_depletion_witness: DepletionWitness
    initial_coverage_certificate: CoverageCertificate
    entry_circle_certificate: CircleInEntryCertificate
    circle_witness: ReplayLateralWitness
    traversal_after: MatTraversalState
    result_state_digest: IdentityDigest

    def __post_init__(self) -> None:
        if type(self) is not InitialCandidateTransaction:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction must use the exact owned type.",
            )
        _digest(
            self.input_identity_digest,
            "initial transaction input identity",
        )
        _digest(
            self.pristine_stock_boundary_digest,
            "initial transaction pristine stock boundary",
        )
        _digest(
            self.result_state_digest,
            "initial transaction result state",
        )
        if type(self.traversal_before) is not MatTraversalState or type(self.traversal_after) is not MatTraversalState:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires exact global traversal parents and children.",
            )
        if type(self.candidate) is not MiddleCurveCandidate:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires one exact finite-lattice candidate.",
            )
        if type(self.entry_depletion_witness) is not DepletionWitness:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires one exact entry depletion witness.",
            )
        if type(self.initial_coverage_certificate) is not CoverageCertificate:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires one exact initial coverage certificate.",
            )
        if type(self.entry_circle_certificate) is not CircleInEntryCertificate:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires one exact entry-circle certificate.",
            )
        if type(self.circle_witness) is not ReplayLateralWitness:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires one exact circle proof bundle.",
            )
        self._validate_traversal()
        self._validate_entry()
        self._validate_circle()

    @classmethod
    def build(
        cls,
        *,
        input_identity_digest: IdentityDigest,
        traversal_before: MatTraversalState,
        candidate: MiddleCurveCandidate,
        pristine_stock_boundary_digest: bytes,
        entry_depletion_witness: DepletionWitness,
        initial_coverage_certificate: CoverageCertificate,
        entry_circle_certificate: CircleInEntryCertificate,
        circle_witness: ReplayLateralWitness,
        traversal_after: MatTraversalState,
        result_state_digest: IdentityDigest,
    ) -> Self:
        """Build and cross-validate one accepted entry launch.

        Args:
            input_identity_digest: Complete planning-input root identity.
            traversal_before: Seeded global graph state.
            candidate: Exact entry-compatible finite-lattice circle.
            pristine_stock_boundary_digest: Exact stock boundary before entry.
            entry_depletion_witness: First-and-only qualified bore depletion.
            initial_coverage_certificate: Empty lateral coverage state.
            entry_circle_certificate: Complete sweep-in-entry proof.
            circle_witness: Ordered containment, motion, depletion, and coverage
                proof for the first full circle.
            traversal_after: Global state after exactly one cursor advance.
            result_state_digest: Content identity of the physical child.

        Returns:
            Frozen, content-addressed launch evidence.

        Raises:
            InvalidInitialCandidateTransactionError: If identities, proof
                chronology, operation ownership, or traversal disagree.
        """
        return cls(
            input_identity_digest,
            traversal_before,
            candidate,
            pristine_stock_boundary_digest,
            entry_depletion_witness,
            initial_coverage_certificate,
            entry_circle_certificate,
            circle_witness,
            traversal_after,
            result_state_digest,
        )

    def _validate_traversal(self) -> None:
        before = self.traversal_before
        if before.active_route_index != 0 or before.pending_transit is not None or any(cursor.accepted_candidate_count != 0 or cursor.terminal for cursor in before.cursors):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction parent must be one untouched seeded traversal.",
            )
        if type(self.candidate.neck_scope) is not NoNeckScope or type(self.candidate.effective_cap_decision) is not FullCapDecision:
            raise InvalidInitialCandidateTransactionError(
                "seeded entry launch cannot manufacture a causal neck passage.",
            )
        try:
            expected_after = before.advance(self.candidate)
        except (
            StaleTraversalCursorError,
            TerminalTraversalCursorError,
        ) as error:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction candidate does not advance its traversal parent.",
            ) from error
        if len(before.cursors) != len(self.traversal_after.cursors):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction must advance exactly one global cursor.",
            )
        changed = sum(
            parent.canonical_bytes != child.canonical_bytes
            for parent, child in zip(
                before.cursors,
                self.traversal_after.cursors,
                strict=True,
            )
        )
        if changed != 1:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction must advance exactly one global cursor.",
            )
        if self.traversal_after != expected_after:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction traversal child contradicts its candidate.",
            )

    def _validate_entry(self) -> None:
        witness = self.entry_depletion_witness
        entry = witness.motion
        if type(entry) is not PreclearedEntry or witness.parent_lineage or witness.tool_radius != entry.tool_radius:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction must begin with one unparented qualified entry depletion.",
            )
        coverage = self.initial_coverage_certificate
        try:
            CoverageCertificate.validate(coverage)
        except InvalidCoverageCertificateError as error:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction coverage certificate is invalid.",
            ) from error
        if (
            coverage.ordered_sweep_records
            or coverage.reachable_material_digest != entry.reachable_domain.certificate.reachable_material_digest
            or coverage.unreachable_residual_digest != entry.reachable_domain.certificate.unreachable_residual_digest
            or coverage.precleared_center != entry.center
            or coverage.precleared_radius != entry.radius
        ):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction coverage does not begin at its qualified entry.",
            )
        certificate = self.entry_circle_certificate
        if (
            certificate.entry_digest != entry.digest
            or certificate.entry_center != entry.center
            or certificate.entry_radius != entry.radius
            or certificate.tool_radius != entry.tool_radius
            or certificate.motion != self.candidate.motion
        ):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction entry-circle proof is cross-wired.",
            )
        if _phase_point(self.candidate) != entry.center:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction candidate phase differs from its entry.",
            )

    def _validate_circle(self) -> None:
        circle = self.circle_witness
        operation = circle.operation
        entry = self.entry_depletion_witness.motion
        if type(entry) is not PreclearedEntry:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction entry witness lost its exact owner.",
            )
        if circle.operation_index != 2 or type(operation) is not CutFullCircleOperation:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction requires one circle at canonical operation index two.",
            )
        if (
            operation.motion != self.candidate.motion
            or operation.neck_scope != self.candidate.neck_scope
            or operation.effective_cap_decision != self.candidate.effective_cap_decision
            or operation.traversal_decision != self.candidate.traversal_decision
        ):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction circle witness contradicts its candidate.",
            )
        if operation.cut_z != entry.cut_plane.cut_z:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction circle witness has a foreign cut depth.",
            )
        containment = circle.containment_certificate
        if (
            type(containment) is not CircleContainmentCertificate
            or containment.domain_digest != entry.reachable_domain.certificate.digest
            or containment.tool_radius != entry.tool_radius
        ):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction circle containment has a foreign domain or tool.",
            )
        expected_lineage = (self.entry_depletion_witness.digest,)
        if (
            circle.motion_witness.stock_lineage_digest
            != _stock_lineage_digest(
                (self.entry_depletion_witness,),
            )
            or circle.depletion_witness.parent_lineage != expected_lineage
        ):
            raise InvalidInitialCandidateTransactionError(
                "initial transaction violates certify-before-deplete circle chronology.",
            )
        if circle.sweep_witness.parent_lineage:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction circle must be the first lateral coverage sweep.",
            )

    @property
    def canonical_bytes(self) -> bytes:
        """Return complete versioned launch evidence.

        Returns:
            Canonical CCAN bytes binding input, proofs, and both state axes.
        """
        return encode_tagged_union(
            b"initial-candidate-transaction-v1",
            encode_component_map(
                {
                    b"candidate": self.candidate.canonical_bytes,
                    b"circle-witness": self.circle_witness.canonical_bytes,
                    b"entry-circle-certificate": self.entry_circle_certificate.canonical_bytes,
                    b"entry-depletion-witness": self.entry_depletion_witness.canonical_bytes,
                    b"initial-coverage-certificate": self.initial_coverage_certificate.canonical_bytes,
                    b"input-identity": encode_bytes(
                        bytes(self.input_identity_digest),
                    ),
                    b"pristine-stock-boundary": encode_bytes(
                        self.pristine_stock_boundary_digest,
                    ),
                    b"result-state": encode_bytes(
                        bytes(self.result_state_digest),
                    ),
                    b"traversal-after": self.traversal_after.canonical_bytes,
                    b"traversal-before": self.traversal_before.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`.

        Returns:
            Content identity of the accepted launch transaction.
        """
        return IdentityDigest(
            hashlib.sha256(self.canonical_bytes).digest(),
        )


@dataclass(frozen=True)
class InitialCandidateEvaluator:
    """Run isolated exact trials from one authenticated planning root."""

    input_identity: InputIdentity
    material_side: MaterialSide
    _pristine_stock: Stock2Area = field(init=False, repr=False)
    _containment: GougeContainment = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if type(self.input_identity) is not InputIdentity:
            raise InvalidInitialCandidateEvaluatorError(
                "initial evaluator requires one exact planning-input identity.",
            )
        if type(self.material_side) is not MaterialSide:
            raise InvalidInitialCandidateEvaluatorError(
                "initial evaluator requires one exact material side.",
            )
        object.__setattr__(
            self,
            "_pristine_stock",
            _pristine_stock(self.input_identity),
        )
        object.__setattr__(
            self,
            "_containment",
            GougeContainment.build(
                self.input_identity.reachable_domain,
            ),
        )

    @classmethod
    def build(
        cls,
        *,
        input_identity: InputIdentity,
        material_side: MaterialSide,
    ) -> Self:
        """Build one launch evaluator bound to complete root authority.

        Args:
            input_identity: Content root for geometry, entry, and all policies.
            material_side: Radial material side deciding cutter rotation.

        Returns:
            Validated short-lived entry evaluator.

        Raises:
            InvalidInitialCandidateEvaluatorError: If either authority has a
                foreign exact type.
        """
        return cls(
            input_identity,
            material_side,
        )

    def evaluate(
        self,
        traversal: MatTraversalState,
        candidate: MiddleCurveCandidate,
    ) -> InitialCandidateTransaction:
        """Evaluate one entry-compatible candidate without state mutation.

        Args:
            traversal: Untouched seeded global MAT state.
            candidate: Exact finite-lattice launch circle.

        Returns:
            Content-addressed launch evidence.

        Raises:
            InitialCandidateStateMismatchError: If graph, policies, cursor, or
                candidate authority differ from `input_identity`.
            InitialCandidatePhaseError: If the candidate phase is not the
                qualified entry center.
            GougeContainmentError: If the circle leaves the entry or design.
            EngagementCapExceededError: If the exact motion exceeds cap.
            UnresolvedMotionEventError: If the exact event proof is incomplete.
        """
        transaction, _, _ = self._evaluate_trial(
            traversal,
            candidate,
        )
        return transaction

    def commit(
        self,
        traversal: MatTraversalState,
        transaction: InitialCandidateTransaction,
    ) -> tuple[GenerationState, MatTraversalState]:
        """Independently reproduce and functionally commit one launch.

        Args:
            traversal: Authoritative global parent named by the transaction.
            transaction: Previously accepted immutable launch evidence.

        Returns:
            Physical generation child and global traversal child.

        Raises:
            StaleInitialCandidateTransactionError: If `traversal` changed.
            InvalidInitialCandidateTransactionError: If the transaction has a
                foreign root or differs from independent replay.
        """
        if type(traversal) is not MatTraversalState or type(transaction) is not InitialCandidateTransaction:
            raise InvalidInitialCandidateTransactionError(
                "initial commit requires exact traversal and transaction types.",
            )
        if transaction.input_identity_digest != self.input_identity.digest:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction belongs to another planning-input root.",
            )
        if transaction.traversal_before.digest != traversal.digest or transaction.traversal_before.canonical_bytes != traversal.canonical_bytes:
            raise StaleInitialCandidateTransactionError(
                "initial transaction parent no longer names authoritative traversal state.",
            )
        reproduced, state, traversal_after = self._evaluate_trial(
            traversal,
            transaction.candidate,
        )
        if reproduced.canonical_bytes != transaction.canonical_bytes:
            raise InvalidInitialCandidateTransactionError(
                "initial transaction differs from independent commit replay.",
            )
        return state, traversal_after

    def _evaluate_trial(
        self,
        traversal: MatTraversalState,
        candidate: MiddleCurveCandidate,
    ) -> tuple[
        InitialCandidateTransaction,
        GenerationState,
        MatTraversalState,
    ]:
        traversal_after = self._validate_parent(
            traversal,
            candidate,
        )
        identity = self.input_identity
        entry = identity.entry
        entry_circle_certificate = entry.certify_first_circle(
            candidate.motion,
        )
        containment: CircleContainmentCertificate = self._containment.certify_full_circle(
            candidate.motion,
            identity.tool_radius,
        )
        operation = CutFullCircleOperation.build(
            motion=candidate.motion,
            cut_z=identity.cut_plane.cut_z,
            material_side=self.material_side,
            neck_scope=candidate.neck_scope,
            effective_cap_decision=candidate.effective_cap_decision,
            traversal_decision=candidate.traversal_decision,
        )

        stock = self._pristine_stock.fork()
        pristine_stock_boundary_digest = MotionCertifier.build(
            stock=stock,
            tool_radius=identity.tool_radius,
        ).canonical_boundary_digest
        entry_depletion_witness = stock.deplete(entry)
        coverage = CoverageLedger.build(
            reachable_domain=identity.reachable_domain,
            precleared_center=entry.center,
            precleared_radius=entry.radius,
        )
        initial_coverage_certificate = coverage.certificate
        certifier = MotionCertifier.build(
            stock=stock,
            tool_radius=identity.tool_radius,
        )
        motion_witness = certifier.certify(
            operation_index=2,
            operation_kind=OperationType.CUT,
            motion=candidate.motion,
            user_cap=identity.user_cap,
            effective_cap=identity.user_cap,
        )
        depletion_witness = stock.deplete(
            candidate.motion,
            identity.tool_radius,
            identity.depletion_policy,
        )
        sweep_witness = coverage.add_sweep(
            candidate.motion,
            identity.tool_radius,
        )
        circle_witness = ReplayLateralWitness(
            operation_index=2,
            operation=operation,
            effective_cap_decision=candidate.effective_cap_decision,
            stock_boundary_digest=certifier.canonical_boundary_digest,
            containment_certificate=containment,
            motion_witness=motion_witness,
            depletion_witness=depletion_witness,
            sweep_witness=sweep_witness,
        )
        local_traversal = TraversalCursorState.before(
            candidate.traversal_decision,
        ).advance(
            candidate.traversal_decision,
        )
        state = GenerationState.build(
            stock=stock,
            coverage=coverage,
            tool_radius=identity.tool_radius,
            phase_point=entry.center,
            traversal=local_traversal,
            passages=_initial_passages(traversal),
            operations=(
                entry.approach,
                entry.plunge,
                operation,
            ),
        )
        transaction = InitialCandidateTransaction.build(
            input_identity_digest=identity.digest,
            traversal_before=traversal,
            candidate=candidate,
            pristine_stock_boundary_digest=pristine_stock_boundary_digest,
            entry_depletion_witness=entry_depletion_witness,
            initial_coverage_certificate=initial_coverage_certificate,
            entry_circle_certificate=entry_circle_certificate,
            circle_witness=circle_witness,
            traversal_after=traversal_after,
            result_state_digest=state.digest,
        )
        return transaction, state, traversal_after

    def _validate_parent(
        self,
        traversal: MatTraversalState,
        candidate: MiddleCurveCandidate,
    ) -> MatTraversalState:
        if type(traversal) is not MatTraversalState or type(candidate) is not MiddleCurveCandidate:
            raise InitialCandidateStateMismatchError(
                "initial evaluation requires exact traversal and candidate types.",
            )
        identity = self.input_identity
        authority = traversal.authority
        axis = authority.axis
        if axis.design_boundary != identity.design_boundary or axis.holes != identity.holes or axis.tool_radius != identity.tool_radius:
            raise InitialCandidateStateMismatchError(
                "initial traversal has foreign design, holes, or tool authority.",
            )
        sampling = identity.mat_sampling_policy
        if axis.station_spacing != sampling.station_spacing or axis.max_sagitta != sampling.max_sagitta or axis.max_refinement_depth != sampling.max_refinement_depth:
            raise InitialCandidateStateMismatchError(
                "initial traversal sampling contradicts its input identity.",
            )
        if authority.inventory.policy != identity.neck_policy or authority.policy != identity.traversal_policy:
            raise InitialCandidateStateMismatchError(
                "initial traversal neck or route policy contradicts its input identity.",
            )
        if (
            traversal.active_route_index != 0
            or traversal.pending_transit is not None
            or any(cursor.accepted_candidate_count != 0 or cursor.terminal for cursor in traversal.cursors)
        ):
            raise InitialCandidateStateMismatchError(
                "initial traversal must be one untouched seeded graph state.",
            )
        if candidate.policy != identity.candidate_policy:
            raise InitialCandidateStateMismatchError(
                "initial candidate policy contradicts its input identity.",
            )
        expected_orientation = identity.cut_direction_policy.circle_orientation(
            self.material_side,
        )
        if candidate.proposal.circle_orientation is not expected_orientation:
            raise InitialCandidateStateMismatchError(
                "initial candidate circle contradicts cut direction.",
            )
        expected_cap = FullCapDecision.build(
            user_cap=identity.user_cap,
            effective_cap=identity.user_cap,
        )
        if (
            type(candidate.neck_scope) is not NoNeckScope
            or type(candidate.effective_cap_decision) is not FullCapDecision
            or candidate.effective_cap_decision != expected_cap
            or candidate.neck_scope != traversal.neck_scope
        ):
            raise InitialCandidateStateMismatchError(
                "initial candidate scope or cap contradicts seeded causal state.",
            )
        try:
            traversal_after = traversal.advance(candidate)
        except (
            StaleTraversalCursorError,
            TerminalTraversalCursorError,
        ) as error:
            raise InitialCandidateStateMismatchError(
                "initial candidate does not begin at the active global cursor.",
            ) from error
        if _phase_point(candidate) != identity.entry.center:
            raise InitialCandidatePhaseError(
                "initial candidate phase must equal the qualified entry center.",
            )
        return traversal_after
