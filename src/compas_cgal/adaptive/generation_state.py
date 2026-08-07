"""Immutable authoritative state for atomic candidate evaluation."""

import hashlib
from dataclasses import dataclass
from dataclasses import field
from typing import Self
from typing import cast

from compas_cgal.adaptive.canonical import canonical_point2_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.errors import CandidateStateMismatchError
from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.errors import InvalidRetraceSegmentOperationError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.neck import NeckPassage
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import EdgeId
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

_DIGEST_SIZE = hashlib.sha256().digest_size


def _identity(value: object, name: str) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidGenerationStateError(f"{name} must be nonempty exact bytes.")
    return value


def _phase_point(motion: ExactCircleMotion) -> Point2[WorldXY]:
    return Point2[WorldXY].build(
        motion.center.x + motion.phase_vector.x,
        motion.center.y + motion.phase_vector.y,
    )


@dataclass(frozen=True)
class TraversalCursorState:
    """One active MAT cursor at an accepted candidate boundary."""

    component_id: ComponentId
    edge_id: EdgeId
    branch_id: BranchId
    cursor: CursorIdentity
    terminal: bool

    def __post_init__(self) -> None:
        _identity(self.component_id, "traversal component identity")
        _identity(self.edge_id, "traversal edge identity")
        _identity(self.branch_id, "traversal branch identity")
        _identity(self.cursor, "traversal cursor identity")
        if type(self.terminal) is not bool:
            raise InvalidGenerationStateError("traversal terminal state must be explicitly boolean.")

    @classmethod
    def before(cls, decision: AdvanceTraversalDecision) -> Self:
        """Build the cursor state immediately before one exact advance.

        Args:
            decision: Candidate-owned traversal advance.

        Returns:
            Nonterminal state at `decision.cursor_before`.

        Raises:
            InvalidGenerationStateError: If `decision` has a foreign type.
        """
        if type(decision) is not AdvanceTraversalDecision:
            raise InvalidGenerationStateError("traversal cursor factory requires one exact advance decision.")
        return cls(
            decision.component_id,
            decision.edge_id,
            decision.branch_id,
            decision.cursor_before,
            False,
        )

    def advance(self, decision: AdvanceTraversalDecision) -> Self:
        """Apply one candidate-owned cursor advance.

        Args:
            decision: Exact advance beginning at this state.

        Returns:
            New cursor state at `decision.cursor_after`.

        Raises:
            CandidateStateMismatchError: If this state is terminal or the
                decision begins at another component, edge, branch, or cursor.
        """
        if type(decision) is not AdvanceTraversalDecision:
            raise CandidateStateMismatchError("candidate traversal must use one exact advance decision.")
        if self.terminal:
            raise CandidateStateMismatchError("candidate cannot advance a terminal traversal cursor.")
        if decision.component_id != self.component_id or decision.edge_id != self.edge_id or decision.branch_id != self.branch_id or decision.cursor_before != self.cursor:
            raise CandidateStateMismatchError("candidate traversal does not begin at authoritative state.")
        return type(self)(
            self.component_id,
            self.edge_id,
            self.branch_id,
            decision.cursor_after,
            decision.makes_cursor_terminal,
        )

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete active-cursor identity.

        Returns:
            Canonical CCAN bytes for state hashing.
        """
        return encode_tagged_union(
            b"generation-traversal-cursor-v1",
            encode_component_map(
                {
                    b"branch-id": encode_bytes(bytes(self.branch_id)),
                    b"component-id": encode_bytes(bytes(self.component_id)),
                    b"cursor": encode_bytes(bytes(self.cursor)),
                    b"edge-id": encode_bytes(bytes(self.edge_id)),
                    b"terminal": encode_boolean(self.terminal),
                }
            ),
        )


@dataclass(frozen=True)
class GenerationState:
    """Alias-safe authoritative state at one accepted lateral boundary."""

    _stock: Stock2Area = field(repr=False)
    _coverage: CoverageLedger = field(repr=False)
    tool_radius: ToolRadius
    phase_point: Point2[WorldXY]
    traversal: TraversalCursorState
    passages: tuple[NeckPassage, ...]
    operations: tuple[CanonicalOperation, ...]
    _stock_boundary_digest: bytes = field(init=False, repr=False)
    _stock_lineage_digest: bytes = field(init=False, repr=False)
    _digest: IdentityDigest = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if type(self._stock) is not Stock2Area:
            raise InvalidGenerationStateError("generation state requires one exact stock owner.")
        if type(self._coverage) is not CoverageLedger:
            raise InvalidGenerationStateError("generation state requires one exact coverage owner.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidGenerationStateError("generation state requires one exact tool radius.")
        if type(self.phase_point) is not Point2:
            raise InvalidGenerationStateError("generation state phase must be one world-XY point.")
        if type(self.traversal) is not TraversalCursorState:
            raise InvalidGenerationStateError("generation state requires one exact traversal cursor.")
        if type(self.passages) is not tuple or any(type(passage) is not NeckPassage for passage in self.passages):
            raise InvalidGenerationStateError("generation state passages must use one exact immutable tuple.")
        if type(self.operations) is not tuple or any(
            type(operation)
            not in (
                ApproachOperation,
                PlungeOperation,
                LinkSegmentOperation,
                CutFullCircleOperation,
                AdvanceSegmentOperation,
                RetraceSegmentOperation,
            )
            for operation in self.operations
        ):
            raise InvalidGenerationStateError("generation state operations must use the closed exact grammar.")

        stock = self._stock.fork()
        coverage = self._coverage.clone()
        passages = tuple(sorted(self.passages, key=lambda passage: passage.scope.canonical_bytes))
        if len(passages) != len({passage.scope.canonical_bytes for passage in passages}):
            raise InvalidGenerationStateError("generation state passage scopes must be unique.")
        object.__setattr__(self, "_stock", stock)
        object.__setattr__(self, "_coverage", coverage)
        object.__setattr__(self, "passages", passages)

        entry = self._validate_chronology()
        self._validate_passages()
        certifier = self._certifier()
        object.__setattr__(
            self,
            "_stock_boundary_digest",
            certifier.canonical_boundary_digest,
        )
        object.__setattr__(
            self,
            "_stock_lineage_digest",
            certifier.stock_lineage_digest,
        )
        state_bytes = encode_tagged_union(
            b"generation-state-v1",
            encode_component_map(
                {
                    b"coverage": coverage.certificate.canonical_bytes,
                    b"entry": encode_bytes(entry.digest),
                    b"operations": encode_sequence(tuple(operation.canonical_bytes for operation in self.operations)),
                    b"passages": encode_sequence(tuple(passage.canonical_bytes for passage in self.passages)),
                    b"phase-point": canonical_point2_bytes(self.phase_point),
                    b"stock-boundary": encode_bytes(certifier.canonical_boundary_digest),
                    b"stock-lineage": encode_bytes(certifier.stock_lineage_digest),
                    b"tool-radius": encode_binary64(float(self.tool_radius.value)),
                    b"traversal": self.traversal.canonical_bytes,
                }
            ),
        )
        object.__setattr__(
            self,
            "_digest",
            IdentityDigest(hashlib.sha256(state_bytes).digest()),
        )

    @classmethod
    def build(
        cls,
        *,
        stock: Stock2Area,
        coverage: CoverageLedger,
        tool_radius: ToolRadius,
        phase_point: Point2[WorldXY],
        traversal: TraversalCursorState,
        passages: tuple[NeckPassage, ...],
        operations: tuple[CanonicalOperation, ...],
    ) -> Self:
        """Build and cross-validate one authoritative state snapshot.

        Args:
            stock: Exact stock owner after the operation prefix.
            coverage: Exact coverage owner after the same lateral prefix.
            tool_radius: Cutter radius bound by every witness.
            phase_point: Exact start point of the next lateral motion.
            traversal: Active MAT cursor after the final advancing operation.
            passages: Current oriented-neck passage states.
            operations: Complete accepted operation prefix.

        Returns:
            Alias-safe, content-addressed generation state.

        Raises:
            InvalidGenerationStateError: If owners, lineages, operations,
                phase, cursor, or passage states disagree.
        """
        return cls(
            stock,
            coverage,
            tool_radius,
            phase_point,
            traversal,
            passages,
            operations,
        )

    @classmethod
    def validate(
        cls,
        state: "GenerationState",
    ) -> None:
        """Revalidate one complete state without mutating its owners.

        Args:
            state: Exact immutable state whose owners and cached identities
                must reproduce through the canonical constructor.

        Raises:
            InvalidGenerationStateError: If the state is foreign, incomplete,
                internally invalid, or contradicts its cached identities.
        """
        if type(state) is not cls:
            raise InvalidGenerationStateError(
                "generation state validation requires the exact owned type.",
            )
        try:
            reproduced = cls(
                state._stock,
                state._coverage,
                state.tool_radius,
                state.phase_point,
                state.traversal,
                state.passages,
                state.operations,
            )
            if (
                reproduced.stock_boundary_digest != state.stock_boundary_digest
                or reproduced.stock_lineage_digest != state.stock_lineage_digest
                or reproduced.digest != state.digest
            ):
                raise InvalidGenerationStateError(
                    "generation state cached identities contradict its owners.",
                )
        except AttributeError as error:
            raise InvalidGenerationStateError(
                "generation state contains incomplete exact state.",
            ) from error

    def _certifier(self) -> MotionCertifier:
        return MotionCertifier.build(
            stock=self._stock,
            tool_radius=self.tool_radius,
        )

    def _validate_chronology(self) -> PreclearedEntry:
        stock_lineage = self._stock.lineage
        coverage_lineage = self._coverage.lineage
        if not stock_lineage or type(stock_lineage[0].motion) is not PreclearedEntry:
            raise InvalidGenerationStateError("generation stock chronology must begin with one precleared entry.")
        entry = stock_lineage[0].motion
        if any(witness.tool_radius != self.tool_radius for witness in stock_lineage) or any(witness.tool_radius != self.tool_radius for witness in coverage_lineage):
            raise InvalidGenerationStateError("generation tool radius contradicts stock or coverage lineage.")
        if len(stock_lineage) != len(coverage_lineage) + 1:
            raise InvalidGenerationStateError("stock and coverage chronology have different lateral lengths.")
        if tuple(witness.motion for witness in stock_lineage[1:]) != tuple(witness.motion for witness in coverage_lineage):
            raise InvalidGenerationStateError("stock and coverage chronology name different exact motions.")
        for index, witness in enumerate(stock_lineage):
            if witness.parent_lineage != tuple(parent.digest for parent in stock_lineage[:index]):
                raise InvalidGenerationStateError("stock depletion chronology has a broken parent lineage.")
        for index, coverage_witness in enumerate(coverage_lineage):
            if coverage_witness.parent_lineage != tuple(parent.digest for parent in coverage_lineage[:index]):
                raise InvalidGenerationStateError("coverage sweep chronology has a broken parent lineage.")

        if len(self.operations) < 3 or self.operations[0] != entry.approach or self.operations[1] != entry.plunge:
            raise InvalidGenerationStateError("generation operation chronology must begin with its exact entry.")
        raw_lateral = self.operations[2:]
        if any(
            type(operation)
            not in (
                LinkSegmentOperation,
                CutFullCircleOperation,
                AdvanceSegmentOperation,
                RetraceSegmentOperation,
            )
            for operation in raw_lateral
        ):
            raise InvalidGenerationStateError("generation operation chronology contains a foreign lateral type.")
        lateral = tuple(
            cast(
                LinkSegmentOperation | CutFullCircleOperation | AdvanceSegmentOperation | RetraceSegmentOperation,
                operation,
            )
            for operation in raw_lateral
        )
        for operation in lateral:
            if type(operation) is RetraceSegmentOperation:
                try:
                    operation._validate()
                except (InvalidRetraceSegmentOperationError, AttributeError) as error:
                    raise InvalidGenerationStateError(
                        "generation state contains an invalid retrace operation.",
                    ) from error
        if any(operation.cut_z != entry.cut_plane.cut_z for operation in lateral):
            raise InvalidGenerationStateError("generation lateral cut depth contradicts its qualified entry.")
        if tuple(operation.motion for operation in lateral) != tuple(witness.motion for witness in coverage_lineage):
            raise InvalidGenerationStateError("operation, stock, and coverage chronology name different motions.")

        if type(lateral[0]) is not CutFullCircleOperation:
            raise InvalidGenerationStateError("generation lateral chronology must begin with its entry circle.")

        current_phase = entry.center
        pending_link: LinkSegmentOperation | None = None
        last_advance: AdvanceTraversalDecision | None = None
        for operation_index, operation in enumerate(lateral, start=2):
            if type(operation) is LinkSegmentOperation:
                if pending_link is not None or operation.motion.start != current_phase:
                    raise InvalidGenerationStateError("link chronology does not begin at the accepted phase.")
                current_phase = operation.motion.end
                pending_link = operation
                continue

            if type(operation) is CutFullCircleOperation:
                if _phase_point(operation.motion) != current_phase:
                    raise InvalidGenerationStateError(
                        "circle phase contradicts the accepted operation chronology.",
                    )
                if operation_index > 2 and pending_link is None:
                    raise InvalidGenerationStateError(
                        "generation circle after entry requires one paired hold link.",
                    )
                if pending_link is not None:
                    hold = pending_link.traversal_decision
                    advance = operation.traversal_decision
                    if (
                        pending_link.neck_scope != operation.neck_scope
                        or pending_link.effective_cap_decision != operation.effective_cap_decision
                        or hold.component_id != advance.component_id
                        or hold.edge_id != advance.edge_id
                        or hold.branch_id != advance.branch_id
                        or hold.cursor_before != advance.cursor_before
                    ):
                        raise InvalidGenerationStateError(
                            "generation link and circle histories are cross-wired.",
                        )
                pending_link = None
                last_advance = operation.traversal_decision
                continue

            if type(operation) is RetraceSegmentOperation:
                if pending_link is not None:
                    raise InvalidGenerationStateError(
                        "route retrace cannot consume a pending hold link.",
                    )
                source_index = operation.decision.source_operation_index
                if source_index != operation_index - 1:
                    raise InvalidGenerationStateError(
                        "route retrace must immediately follow its named source.",
                    )
                source = self.operations[source_index]
                if type(source) is not AdvanceSegmentOperation:
                    raise InvalidGenerationStateError(
                        "route retrace source is not an advancing segment.",
                    )
                source_digest = IdentityDigest(
                    hashlib.sha256(source.canonical_bytes).digest(),
                )
                if (
                    source_digest != operation.decision.source_operation_digest
                    or operation.motion.start != current_phase
                    or operation.motion.start != source.motion.end
                    or operation.motion.end != source.motion.start
                    or operation.cut_z != source.cut_z
                    or operation.effective_cap_decision != source.effective_cap_decision
                ):
                    raise InvalidGenerationStateError(
                        "route retrace contradicts its source or physical phase.",
                    )
                current_phase = operation.motion.end
                continue

            if type(operation) is AdvanceSegmentOperation:
                if pending_link is not None:
                    raise InvalidGenerationStateError(
                        "hold link must be paired with a full circle, not an advancing segment.",
                    )
                if operation.motion.start != current_phase:
                    raise InvalidGenerationStateError(
                        "advancing segment does not begin at the accepted phase.",
                    )
                current_phase = operation.motion.end
                last_advance = operation.traversal_decision
                continue

            raise InvalidGenerationStateError(
                "operation chronology contains an unknown lateral type.",
            )

        if pending_link is not None:
            raise InvalidGenerationStateError(
                "generation state cannot end with a dangling hold link.",
            )
        if last_advance is None:
            raise InvalidGenerationStateError(
                "generation state requires one accepted traversal advance.",
            )
        if current_phase != self.phase_point:
            raise InvalidGenerationStateError(
                "generation phase contradicts its final advancing operation.",
            )
        decision = last_advance
        if (
            self.traversal.component_id != decision.component_id
            or self.traversal.edge_id != decision.edge_id
            or self.traversal.branch_id != decision.branch_id
            or self.traversal.cursor != decision.cursor_after
            or self.traversal.terminal != decision.makes_cursor_terminal
        ):
            raise InvalidGenerationStateError(
                "generation traversal contradicts its final advancing operation.",
            )

        certificate = self._coverage.certificate
        if (
            certificate.reachable_material_digest != entry.reachable_domain.certificate.reachable_material_digest
            or certificate.unreachable_residual_digest != entry.reachable_domain.certificate.unreachable_residual_digest
            or certificate.precleared_center != entry.center
            or certificate.precleared_radius != entry.radius
        ):
            raise InvalidGenerationStateError("coverage chronology contradicts the precleared entry.")
        return entry

    def _validate_passages(self) -> None:
        by_scope = {passage.scope: passage for passage in self.passages}
        expected_states = {scope: PassageState.UNVISITED for scope in by_scope}
        for operation in self.operations[2:]:
            if type(operation) is RetraceSegmentOperation:
                continue
            advancing_operation: CutFullCircleOperation | AdvanceSegmentOperation
            if type(operation) is CutFullCircleOperation:
                advancing_operation = operation
            elif type(operation) is AdvanceSegmentOperation:
                advancing_operation = operation
            else:
                continue
            if type(advancing_operation.neck_scope) is OrientedNeckScope:
                decision = advancing_operation.effective_cap_decision
                if type(decision) is not NeckCapDecision:
                    raise InvalidGenerationStateError("oriented operation omits its exact neck-cap transition.")
                scope = advancing_operation.neck_scope
                passage = by_scope.get(scope)
                expected_before = expected_states.get(scope)
                if (
                    passage is None
                    or expected_before is None
                    or passage.neck.evidence_digest != decision.neck_evidence_digest
                    or passage.neck.width_class_id != decision.width_class_id
                    or decision.passage_before is not expected_before
                ):
                    raise InvalidGenerationStateError("generation passage transition breaks its accepted chronology.")
                expected_states[scope] = decision.passage_after
        for scope, passage in by_scope.items():
            if passage.state is not expected_states[scope]:
                raise InvalidGenerationStateError("generation passage state contradicts operation chronology.")

    def fork_stock(self) -> Stock2Area:
        """Return an independent stock fork for one candidate trial.

        Returns:
            Alias-safe exact stock owner.
        """
        return self._stock.fork()

    def clone_coverage(self) -> CoverageLedger:
        """Return an independent coverage fork for one candidate trial.

        Returns:
            Alias-safe exact coverage owner.
        """
        return self._coverage.clone()

    def passage(self, scope: OrientedNeckScope) -> NeckPassage:
        """Resolve one exact oriented passage by owner and direction.

        Args:
            scope: Candidate-supplied oriented neck identity.

        Returns:
            Current immutable passage state.

        Raises:
            CandidateStateMismatchError: If no passage owns `scope`.
        """
        if type(scope) is not OrientedNeckScope:
            raise CandidateStateMismatchError("candidate passage lookup requires one oriented neck scope.")
        for passage in self.passages:
            if passage.scope == scope:
                return passage
        raise CandidateStateMismatchError("candidate oriented neck scope has no authoritative passage.")

    @property
    def stock_boundary_digest(self) -> bytes:
        """Return the exact native stock-boundary identity.

        Returns:
            SHA-256 boundary digest.
        """
        return self._stock_boundary_digest

    @property
    def stock_lineage_digest(self) -> bytes:
        """Return the accepted depletion-lineage identity.

        Returns:
            SHA-256 ordered-lineage digest.
        """
        return self._stock_lineage_digest

    @property
    def digest(self) -> IdentityDigest:
        """Return the complete authoritative state identity.

        Returns:
            SHA-256 digest binding every generation-state input.
        """
        if type(self._digest) is not bytes or len(self._digest) != _DIGEST_SIZE:
            raise InvalidGenerationStateError("generation state digest is not one exact SHA-256 identity.")
        return self._digest
