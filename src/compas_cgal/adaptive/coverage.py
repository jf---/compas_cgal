import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Self

from compas_cgal import _coverage_2
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import IncompletePocketCoverageError
from compas_cgal.adaptive.errors import InvalidCoverageCertificateError
from compas_cgal.adaptive.errors import InvalidCoverageSweepError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _exact_bytes(value: object, name: str) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidCoverageCertificateError(f"{name} must be nonempty exact bytes.")
    return value


def _exact_bytes_tuple(values: object, name: str) -> tuple[bytes, ...]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise InvalidCoverageCertificateError(f"{name} must be a finite native byte sequence.") from None
    validated: list[bytes] = []
    for value in values:
        if type(value) is not bytes or not value:
            raise InvalidCoverageCertificateError(f"{name} must contain nonempty exact bytes.")
        validated.append(value)
    return tuple(validated)


@dataclass(frozen=True)
class SweepWitness:
    motion: ExactSegmentMotion | ExactCircleMotion
    tool_radius: ToolRadius
    native_strategy_version: bytes
    native_structural_record: bytes
    parent_lineage: tuple[bytes, ...]

    def __post_init__(self) -> None:
        type(self).validate(self)

    @classmethod
    def validate(cls, witness: "SweepWitness") -> None:
        if type(witness) is not cls:
            raise InvalidCoverageSweepError("coverage witness must use the exact owned type.")
        if type(witness.motion) not in (ExactSegmentMotion, ExactCircleMotion):
            raise InvalidCoverageSweepError("coverage witness requires one exact motion.")
        if type(witness.tool_radius) is not ToolRadius:
            raise InvalidCoverageSweepError("coverage witness requires one exact typed tool radius.")
        if type(witness.native_strategy_version) is not bytes or not witness.native_strategy_version:
            raise InvalidCoverageSweepError("coverage witness native strategy must be nonempty bytes.")
        if type(witness.native_structural_record) is not bytes or not witness.native_structural_record:
            raise InvalidCoverageSweepError("coverage witness native structural record must be nonempty bytes.")
        if type(witness.parent_lineage) is not tuple or any(type(digest) is not bytes or len(digest) != 32 for digest in witness.parent_lineage):
            raise InvalidCoverageSweepError("coverage witness parent lineage must contain exact SHA-256 digests.")

    @property
    def canonical_bytes(self) -> bytes:
        type(self).validate(self)
        return encode_tagged_union(
            b"coverage-sweep-witness-v1",
            encode_component_map(
                {
                    b"motion": canonical_task1_bytes(self.motion),
                    b"native-structural-record": encode_bytes(self.native_structural_record),
                    b"native-strategy-version": encode_bytes(self.native_strategy_version),
                    b"parent-lineage": encode_sequence(tuple(encode_bytes(digest) for digest in self.parent_lineage)),
                    b"tool-radius": encode_binary64(float(self.tool_radius.value)),
                }
            ),
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()


@dataclass(frozen=True)
class CoverageCertificate:
    reachable_material_digest: bytes
    unreachable_residual_digest: bytes
    precleared_center: Point2[WorldXY]
    precleared_radius: EntryRadius
    native_strategy_version: bytes
    ordered_sweep_records: tuple[bytes, ...]
    residual_component_records: tuple[bytes, ...]
    residual_component_count: int
    exact_residual_relation: bool
    exact_residual_empty: bool

    @classmethod
    def validate(cls, certificate: "CoverageCertificate") -> None:
        if type(certificate) is not cls:
            raise InvalidCoverageCertificateError("coverage certificate must use the exact owned type.")
        for name, digest in (
            ("reachable material digest", certificate.reachable_material_digest),
            ("unreachable residual digest", certificate.unreachable_residual_digest),
        ):
            if type(digest) is not bytes or len(digest) != 32:
                raise InvalidCoverageCertificateError(f"{name} must be one exact SHA-256 digest.")
        if type(certificate.precleared_center) is not Point2:
            raise InvalidCoverageCertificateError("coverage certificate preclear center must be an exact world-XY point.")
        if type(certificate.precleared_radius) is not EntryRadius:
            raise InvalidCoverageCertificateError("coverage certificate preclear radius must be one exact entry radius.")
        _exact_bytes(certificate.native_strategy_version, "coverage strategy")
        if type(certificate.ordered_sweep_records) is not tuple or any(type(record) is not bytes or not record for record in certificate.ordered_sweep_records):
            raise InvalidCoverageCertificateError("ordered sweep records must be an exact byte tuple.")
        if type(certificate.residual_component_records) is not tuple or any(type(record) is not bytes or not record for record in certificate.residual_component_records):
            raise InvalidCoverageCertificateError("residual component records must be an exact byte tuple.")
        if certificate.residual_component_records != tuple(sorted(certificate.residual_component_records)):
            raise InvalidCoverageCertificateError("residual component records must be canonically sorted.")
        if (
            type(certificate.residual_component_count) is not int
            or certificate.residual_component_count < 0
            or certificate.residual_component_count != len(certificate.residual_component_records)
        ):
            raise InvalidCoverageCertificateError("residual component count must match its exact records.")
        if type(certificate.exact_residual_relation) is not bool or not certificate.exact_residual_relation:
            raise InvalidCoverageCertificateError("coverage certificate residual relation is not exact.")
        if type(certificate.exact_residual_empty) is not bool:
            raise InvalidCoverageCertificateError("coverage certificate empty state must be exact bool.")
        if certificate.exact_residual_empty != (certificate.residual_component_count == 0):
            raise InvalidCoverageCertificateError("coverage certificate empty state contradicts residual components.")

    @property
    def canonical_bytes(self) -> bytes:
        type(self).validate(self)
        return encode_tagged_union(
            b"coverage-certificate-v2",
            encode_component_map(
                {
                    b"exact-residual-empty": encode_boolean(self.exact_residual_empty),
                    b"exact-residual-relation": encode_boolean(self.exact_residual_relation),
                    b"native-strategy-version": encode_bytes(self.native_strategy_version),
                    b"ordered-sweep-records": encode_sequence(self.ordered_sweep_records),
                    b"precleared-center": encode_sequence(
                        (
                            encode_binary64(float(self.precleared_center.x)),
                            encode_binary64(float(self.precleared_center.y)),
                        )
                    ),
                    b"precleared-radius": encode_binary64(float(self.precleared_radius.value)),
                    b"reachable-material-digest": encode_bytes(self.reachable_material_digest),
                    b"residual-component-count": encode_integer(self.residual_component_count),
                    b"residual-component-records": encode_sequence(self.residual_component_records),
                    b"unreachable-residual-digest": encode_bytes(self.unreachable_residual_digest),
                }
            ),
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()


class CoverageLedger:
    def __init__(
        self,
        *,
        native: _coverage_2.Coverage2,
        reachable_material_digest: bytes,
        unreachable_residual_digest: bytes,
        precleared_center: Point2[WorldXY],
        precleared_radius: EntryRadius,
        lineage: tuple[SweepWitness, ...],
    ) -> None:
        self._native = native
        self._reachable_material_digest = reachable_material_digest
        self._unreachable_residual_digest = unreachable_residual_digest
        self._precleared_center = precleared_center
        self._precleared_radius = precleared_radius
        self._lineage = lineage

    @classmethod
    def build(
        cls,
        *,
        reachable_domain: ReachableDomain,
        precleared_center: Point2[WorldXY],
        precleared_radius: EntryRadius,
    ) -> Self:
        if type(reachable_domain) is not ReachableDomain:
            raise InvalidCoverageCertificateError("coverage ledger requires one exact reachable-domain owner.")
        if type(precleared_center) is not Point2:
            raise InvalidCoverageCertificateError("coverage ledger preclear center must be an exact world-XY point.")
        if type(precleared_radius) is not EntryRadius:
            raise InvalidCoverageCertificateError("coverage ledger preclear radius must be one exact entry radius.")
        native = _coverage_2.Coverage2(
            reachable_domain.reachable_material,
            precleared_center.x,
            precleared_center.y,
            precleared_radius.value,
        )
        ledger = cls(
            native=native,
            reachable_material_digest=reachable_domain.certificate.reachable_material_digest,
            unreachable_residual_digest=reachable_domain.certificate.unreachable_residual_digest,
            precleared_center=precleared_center,
            precleared_radius=precleared_radius,
            lineage=(),
        )
        CoverageCertificate.validate(ledger.certificate)
        return ledger

    @property
    def lineage(self) -> tuple[SweepWitness, ...]:
        return self._lineage

    @property
    def residual(self) -> _coverage_2.ExactRegion2:
        return self._native.residual()

    @property
    def residual_component_records(self) -> tuple[bytes, ...]:
        return _exact_bytes_tuple(
            self._native.residual_component_records,
            "residual component records",
        )

    @property
    def certificate(self) -> CoverageCertificate:
        records = self.residual_component_records
        certificate = CoverageCertificate(
            reachable_material_digest=self._reachable_material_digest,
            unreachable_residual_digest=self._unreachable_residual_digest,
            precleared_center=self._precleared_center,
            precleared_radius=self._precleared_radius,
            native_strategy_version=_exact_bytes(
                self._native.strategy_version,
                "coverage strategy",
            ),
            ordered_sweep_records=tuple(witness.canonical_bytes for witness in self._lineage),
            residual_component_records=records,
            residual_component_count=self._native.residual_component_count(),
            exact_residual_relation=self._native.exact_residual_relation(),
            exact_residual_empty=self._native.residual_is_empty(),
        )
        CoverageCertificate.validate(certificate)
        return certificate

    def clone(self) -> Self:
        return type(self)(
            native=self._native.clone(),
            reachable_material_digest=self._reachable_material_digest,
            unreachable_residual_digest=self._unreachable_residual_digest,
            precleared_center=self._precleared_center,
            precleared_radius=self._precleared_radius,
            lineage=self._lineage,
        )

    def add_sweep(
        self,
        motion: ExactSegmentMotion | ExactCircleMotion,
        tool_radius: ToolRadius,
    ) -> SweepWitness:
        if type(motion) not in (ExactSegmentMotion, ExactCircleMotion):
            raise InvalidCoverageSweepError("coverage sweep requires one exact motion.")
        if type(tool_radius) is not ToolRadius:
            raise InvalidCoverageSweepError("coverage sweep requires one exact typed tool radius.")
        trial = self._native.clone()
        if isinstance(motion, ExactSegmentMotion):
            native_record = trial.add_segment_sweep(
                motion.start.x,
                motion.start.y,
                motion.end.x,
                motion.end.y,
                tool_radius.value,
            )
            matches = native_record.matches_exact_segment(
                motion.start.x,
                motion.start.y,
                motion.end.x,
                motion.end.y,
                tool_radius.value,
            )
        else:
            native_record = trial.add_full_circle_sweep(
                motion.center.x,
                motion.center.y,
                motion.phase_vector.x,
                motion.phase_vector.y,
                tool_radius.value,
            )
            matches = native_record.matches_exact_full_circle(
                motion.center.x,
                motion.center.y,
                motion.phase_vector.x,
                motion.phase_vector.y,
                tool_radius.value,
            )
        strategy = _exact_bytes(
            native_record.strategy_version,
            "sweep strategy version",
        )
        structural_record = _exact_bytes(
            native_record.structural_record,
            "sweep structural record",
        )
        if not matches or strategy != trial.strategy_version:
            raise InvalidCoverageSweepError("native coverage sweep does not bind the exact motion inputs.")
        expected_native_records = tuple(witness.native_structural_record for witness in self._lineage) + (structural_record,)
        if _exact_bytes_tuple(trial.sweep_records, "native sweep records") != expected_native_records:
            raise InvalidCoverageSweepError("native coverage sweep history diverged from owned lineage.")
        witness = SweepWitness(
            motion=motion,
            tool_radius=tool_radius,
            native_strategy_version=strategy,
            native_structural_record=structural_record,
            parent_lineage=tuple(parent.digest for parent in self._lineage),
        )
        next_lineage = self._lineage + (witness,)
        trial_ledger = type(self)(
            native=trial,
            reachable_material_digest=self._reachable_material_digest,
            unreachable_residual_digest=self._unreachable_residual_digest,
            precleared_center=self._precleared_center,
            precleared_radius=self._precleared_radius,
            lineage=next_lineage,
        )
        CoverageCertificate.validate(trial_ledger.certificate)
        self._native = trial
        self._lineage = next_lineage
        return witness

    def residual_is_empty(self) -> bool:
        return self._native.residual_is_empty()

    def residual_component_count(self) -> int:
        return self._native.residual_component_count()

    def require_complete(self) -> None:
        if not self.residual_is_empty():
            count = self.residual_component_count()
            raise IncompletePocketCoverageError(f"exact pocket coverage leaves {count} residual component(s).")
