import hashlib
from dataclasses import dataclass

from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import CircleContainmentCertificate
from compas_cgal.adaptive.containment import CircleInEntryCertificate
from compas_cgal.adaptive.containment import SegmentContainmentCertificate
from compas_cgal.adaptive.coverage import CoverageCertificate
from compas_cgal.adaptive.coverage import SweepWitness
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.errors import InvalidReplayTraceError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion_certificate import MotionWitness
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import EffectiveCapDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.toolpath import OperationType


def _require_replay_digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != hashlib.sha256().digest_size:
        raise InvalidReplayTraceError(f"{name} must be one exact SHA-256 digest.")
    return value


def _stock_lineage_digest(
    witnesses: tuple[DepletionWitness, ...],
) -> bytes:
    return hashlib.sha256(
        encode_tagged_union(
            b"exact-stock-lineage-v1",
            encode_sequence(tuple(bytes(witness.digest) for witness in witnesses)),
        )
    ).digest()


@dataclass(frozen=True)
class ReplayLateralWitness:
    """One ordered, cross-validated proof bundle for a lateral operation."""

    operation_index: int
    operation: LinkSegmentOperation | CutFullCircleOperation
    effective_cap_decision: EffectiveCapDecision
    stock_boundary_digest: bytes
    containment_certificate: SegmentContainmentCertificate | CircleContainmentCertificate
    motion_witness: MotionWitness
    depletion_witness: DepletionWitness
    sweep_witness: SweepWitness

    def __post_init__(self) -> None:
        type(self).validate(self)

    @classmethod
    def validate(
        cls,
        witness: "ReplayLateralWitness",
    ) -> None:
        if type(witness) is not cls:
            raise InvalidReplayTraceError("lateral replay witness must use the exact owned type.")
        if type(witness.operation_index) is not int or witness.operation_index < 2:
            raise InvalidReplayTraceError("lateral replay witness requires one canonical operation index.")
        if type(witness.operation) not in (
            LinkSegmentOperation,
            CutFullCircleOperation,
        ):
            raise InvalidReplayTraceError("lateral replay witness requires one exact lateral operation.")
        if type(witness.effective_cap_decision) not in (
            FullCapDecision,
            NeckCapDecision,
        ):
            raise InvalidReplayTraceError("lateral replay witness requires one exact recomputed cap decision.")
        if witness.effective_cap_decision != witness.operation.effective_cap_decision:
            raise InvalidReplayTraceError("recomputed cap decision contradicts its recorded operation.")
        _require_replay_digest(
            witness.stock_boundary_digest,
            "lateral pre-motion stock boundary digest",
        )
        if type(witness.motion_witness) is not MotionWitness:
            raise InvalidReplayTraceError("lateral replay witness requires one exact motion witness.")
        if type(witness.depletion_witness) is not DepletionWitness:
            raise InvalidReplayTraceError("lateral replay witness requires one exact depletion witness.")
        if type(witness.sweep_witness) is not SweepWitness:
            raise InvalidReplayTraceError("lateral replay witness requires one exact coverage witness.")

        operation_kind = OperationType.LINK
        containment_type: type[SegmentContainmentCertificate] | type[CircleContainmentCertificate] = SegmentContainmentCertificate
        if type(witness.operation) is CutFullCircleOperation:
            operation_kind = OperationType.CUT
            containment_type = CircleContainmentCertificate
        if type(witness.containment_certificate) is not containment_type:
            raise InvalidReplayTraceError("lateral containment proof contradicts its operation kind.")
        if witness.motion_witness.operation_index != witness.operation_index or witness.motion_witness.operation_kind is not operation_kind:
            raise InvalidReplayTraceError("motion witness contradicts its canonical operation position.")
        motion = witness.operation.motion
        if not (witness.containment_certificate.motion == witness.motion_witness.motion == witness.depletion_witness.motion == witness.sweep_witness.motion == motion):
            raise InvalidReplayTraceError("lateral proof bundle contains cross-wired exact motions.")
        if not (witness.containment_certificate.tool_radius == witness.depletion_witness.tool_radius == witness.sweep_witness.tool_radius):
            raise InvalidReplayTraceError("lateral proof bundle contains cross-wired tool radii.")
        if (
            witness.motion_witness.user_cap_bytes != witness.effective_cap_decision.user_cap_bytes
            or witness.motion_witness.effective_cap_bytes != witness.effective_cap_decision.effective_cap_bytes
        ):
            raise InvalidReplayTraceError("motion witness contradicts the recomputed effective-cap decision.")

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete versioned per-operation proof record.

        Returns:
            Canonical CCAN bytes preserving all accepted fresh evidence.
        """
        type(self).validate(self)
        return encode_tagged_union(
            b"fresh-replay-lateral-witness-v1",
            encode_component_map(
                {
                    b"containment-certificate": self.containment_certificate.canonical_bytes,
                    b"depletion-witness": self.depletion_witness.canonical_bytes,
                    b"effective-cap-decision": self.effective_cap_decision.canonical_bytes,
                    b"motion-witness": self.motion_witness.canonical_bytes,
                    b"operation": self.operation.canonical_bytes,
                    b"operation-index": encode_integer(self.operation_index),
                    b"stock-boundary-digest": encode_bytes(self.stock_boundary_digest),
                    b"sweep-witness": self.sweep_witness.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`.

        Returns:
            Immutable per-operation replay identity.
        """
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class FreshReplayTrace:
    """Immutable evidence captured before terminal replay gates are applied."""

    pristine_stock_boundary_digest: bytes
    post_entry_stock_boundary_digest: bytes
    post_entry_stock_lineage_digest: bytes
    entry_depletion_witness: DepletionWitness
    first_circle_entry_certificate: CircleInEntryCertificate
    initial_coverage_certificate: CoverageCertificate
    lateral_witnesses: tuple[ReplayLateralWitness, ...]
    terminal_stock_boundary_digest: bytes
    terminal_stock_lineage_digest: bytes
    terminal_coverage_certificate: CoverageCertificate

    def __post_init__(self) -> None:
        type(self).validate(self)

    @classmethod
    def validate(
        cls,
        trace: "FreshReplayTrace",
    ) -> None:
        if type(trace) is not cls:
            raise InvalidReplayTraceError("fresh replay trace must use the exact owned type.")
        for name, digest in (
            (
                "pristine stock boundary digest",
                trace.pristine_stock_boundary_digest,
            ),
            (
                "post-entry stock boundary digest",
                trace.post_entry_stock_boundary_digest,
            ),
            (
                "post-entry stock lineage digest",
                trace.post_entry_stock_lineage_digest,
            ),
            (
                "terminal stock boundary digest",
                trace.terminal_stock_boundary_digest,
            ),
            (
                "terminal stock lineage digest",
                trace.terminal_stock_lineage_digest,
            ),
        ):
            _require_replay_digest(digest, name)
        if type(trace.entry_depletion_witness) is not DepletionWitness:
            raise InvalidReplayTraceError("fresh replay trace requires one exact entry depletion witness.")
        entry = trace.entry_depletion_witness.motion
        if type(entry) is not PreclearedEntry or trace.entry_depletion_witness.parent_lineage:
            raise InvalidReplayTraceError("fresh replay trace must begin with one unparented precleared entry.")
        if type(trace.first_circle_entry_certificate) is not CircleInEntryCertificate:
            raise InvalidReplayTraceError("fresh replay trace requires one exact first-circle entry proof.")
        CoverageCertificate.validate(trace.initial_coverage_certificate)
        CoverageCertificate.validate(trace.terminal_coverage_certificate)
        if trace.initial_coverage_certificate.ordered_sweep_records:
            raise InvalidReplayTraceError("initial replay coverage must precede every lateral sweep.")
        if type(trace.lateral_witnesses) is not tuple or not trace.lateral_witnesses:
            raise InvalidReplayTraceError("fresh replay trace requires one or more lateral proof bundles.")
        if any(type(witness) is not ReplayLateralWitness for witness in trace.lateral_witnesses):
            raise InvalidReplayTraceError("fresh replay trace contains a foreign lateral witness.")
        expected_indices = tuple(range(2, 2 + len(trace.lateral_witnesses)))
        if tuple(witness.operation_index for witness in trace.lateral_witnesses) != expected_indices:
            raise InvalidReplayTraceError("lateral replay witness operation chronology is not canonical.")

        first_circle = next(
            (witness for witness in trace.lateral_witnesses if type(witness.operation) is CutFullCircleOperation),
            None,
        )
        if first_circle is None:
            raise InvalidReplayTraceError("fresh replay trace has no accepted full circle.")
        entry_certificate = trace.first_circle_entry_certificate
        if (
            entry_certificate.entry_digest != entry.digest
            or entry_certificate.entry_center != entry.center
            or entry_certificate.entry_radius != entry.radius
            or entry_certificate.tool_radius != entry.tool_radius
            or entry_certificate.motion != first_circle.operation.motion
        ):
            raise InvalidReplayTraceError("first-circle entry proof contradicts the fresh entry or operation.")

        initial_coverage = trace.initial_coverage_certificate
        terminal_coverage = trace.terminal_coverage_certificate
        if (
            initial_coverage.reachable_material_digest != entry.reachable_domain.certificate.reachable_material_digest
            or initial_coverage.unreachable_residual_digest != entry.reachable_domain.certificate.unreachable_residual_digest
            or initial_coverage.precleared_center != entry.center
            or initial_coverage.precleared_radius != entry.radius
        ):
            raise InvalidReplayTraceError("initial coverage contradicts the authenticated precleared entry.")
        if (
            terminal_coverage.reachable_material_digest != initial_coverage.reachable_material_digest
            or terminal_coverage.unreachable_residual_digest != initial_coverage.unreachable_residual_digest
            or terminal_coverage.precleared_center != initial_coverage.precleared_center
            or terminal_coverage.precleared_radius != initial_coverage.precleared_radius
            or terminal_coverage.native_strategy_version != initial_coverage.native_strategy_version
        ):
            raise InvalidReplayTraceError("terminal coverage changed its fresh initial-state authority.")

        depletion_lineage: tuple[DepletionWitness, ...] = (trace.entry_depletion_witness,)
        expected_post_entry_lineage = _stock_lineage_digest(depletion_lineage)
        if trace.post_entry_stock_lineage_digest != expected_post_entry_lineage:
            raise InvalidReplayTraceError("post-entry stock lineage contradicts its depletion witness.")
        if trace.lateral_witnesses[0].stock_boundary_digest != trace.post_entry_stock_boundary_digest:
            raise InvalidReplayTraceError("first lateral proof does not observe the post-entry stock boundary.")
        sweep_lineage: tuple[SweepWitness, ...] = ()
        domain_digest = entry.reachable_domain.certificate.digest
        for witness in trace.lateral_witnesses:
            ReplayLateralWitness.validate(witness)
            if witness.containment_certificate.domain_digest != domain_digest:
                raise InvalidReplayTraceError("lateral containment proof has a foreign reachable-domain owner.")
            if witness.motion_witness.stock_lineage_digest != _stock_lineage_digest(depletion_lineage):
                raise InvalidReplayTraceError("motion witness observed the wrong pre-depletion stock lineage.")
            if witness.depletion_witness.parent_lineage != tuple(parent.digest for parent in depletion_lineage):
                raise InvalidReplayTraceError("depletion witness parent lineage breaks replay chronology.")
            if witness.sweep_witness.parent_lineage != tuple(parent.digest for parent in sweep_lineage):
                raise InvalidReplayTraceError("coverage witness parent lineage breaks replay chronology.")
            depletion_lineage += (witness.depletion_witness,)
            sweep_lineage += (witness.sweep_witness,)
        if trace.terminal_stock_lineage_digest != _stock_lineage_digest(depletion_lineage):
            raise InvalidReplayTraceError("terminal stock lineage contradicts the ordered depletion witnesses.")
        if terminal_coverage.ordered_sweep_records != tuple(witness.canonical_bytes for witness in sweep_lineage):
            raise InvalidReplayTraceError("terminal coverage contradicts the ordered sweep witnesses.")

    @property
    def canonical_bytes(self) -> bytes:
        """Return the versioned fresh replay evidence record.

        Returns:
            Canonical CCAN bytes for deterministic replay handoff.
        """
        type(self).validate(self)
        return encode_tagged_union(
            b"fresh-replay-trace-v1",
            encode_component_map(
                {
                    b"entry-depletion-witness": self.entry_depletion_witness.canonical_bytes,
                    b"first-circle-entry-certificate": self.first_circle_entry_certificate.canonical_bytes,
                    b"initial-coverage-certificate": self.initial_coverage_certificate.canonical_bytes,
                    b"lateral-witnesses": encode_sequence(tuple(witness.canonical_bytes for witness in self.lateral_witnesses)),
                    b"post-entry-stock-boundary-digest": encode_bytes(self.post_entry_stock_boundary_digest),
                    b"post-entry-stock-lineage-digest": encode_bytes(self.post_entry_stock_lineage_digest),
                    b"pristine-stock-boundary-digest": encode_bytes(self.pristine_stock_boundary_digest),
                    b"terminal-coverage-certificate": self.terminal_coverage_certificate.canonical_bytes,
                    b"terminal-stock-boundary-digest": encode_bytes(self.terminal_stock_boundary_digest),
                    b"terminal-stock-lineage-digest": encode_bytes(self.terminal_stock_lineage_digest),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`.

        Returns:
            Immutable fresh replay trace identity.
        """
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())
