import hashlib
from dataclasses import dataclass
from typing import Self

from compas_cgal import _containment_2
from compas_cgal.adaptive.canonical import canonical_point2_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidContainmentCertificateError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

_SHA256_BYTES = hashlib.sha256().digest_size


def _digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _SHA256_BYTES:
        raise InvalidContainmentCertificateError(f"{name} must be one SHA-256 digest.")
    return value


def _record(
    value: object,
) -> _containment_2.ContainmentRecord2:
    if type(value) is not _containment_2.ContainmentRecord2:
        raise InvalidContainmentCertificateError("containment certificate requires one exact native record.")
    if value.strategy_version != _containment_2.containment_strategy_version():
        raise InvalidContainmentCertificateError("containment certificate uses an unknown native strategy.")
    if type(value.structural_record) is not bytes or not value.structural_record:
        raise InvalidContainmentCertificateError("containment certificate has no native structural record.")
    if type(value.contained) is not bool:
        raise InvalidContainmentCertificateError("containment certificate verdict must be exact bool.")
    return value


def _common_certificate_bytes(
    *,
    tag: bytes,
    domain_or_entry_digest: bytes,
    motion_or_entry: bytes,
    tool_radius: ToolRadius,
    native_record: _containment_2.ContainmentRecord2,
) -> bytes:
    return encode_tagged_union(
        tag,
        encode_component_map(
            {
                b"authority-digest": encode_bytes(domain_or_entry_digest),
                b"geometry": motion_or_entry,
                b"native-structural-record": encode_bytes(native_record.structural_record),
                b"native-strategy-version": encode_bytes(native_record.strategy_version),
                b"tool-radius": encode_binary64(float(tool_radius.value)),
            }
        ),
    )


@dataclass(frozen=True)
class SegmentContainmentCertificate:
    """Replay-bound proof that one exact segment cutter sweep lies in `D`."""

    domain_digest: bytes
    motion: ExactSegmentMotion
    tool_radius: ToolRadius
    native_record: _containment_2.ContainmentRecord2

    def __post_init__(self) -> None:
        _digest(self.domain_digest, "containment domain digest")
        if type(self.motion) is not ExactSegmentMotion:
            raise InvalidContainmentCertificateError("segment certificate requires one exact segment motion.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidContainmentCertificateError("segment certificate requires one typed tool radius.")
        record = _record(self.native_record)
        if not record.contained or not record.guide_anchor_in_center_domain:
            raise InvalidContainmentCertificateError("segment certificate does not contain one accepted exact verdict.")
        if not record.matches_exact_segment(
            self.domain_digest,
            self.motion.start.x,
            self.motion.start.y,
            self.motion.end.x,
            self.motion.end.y,
            self.tool_radius.value,
        ):
            raise InvalidContainmentCertificateError("segment native record contradicts its exact motion.")

    @property
    def endpoints_in_center_domain(self) -> bool:
        return self.native_record.guide_anchor_in_center_domain

    @property
    def native_strategy_version(self) -> bytes:
        return self.native_record.strategy_version

    @property
    def native_structural_record(self) -> bytes:
        return self.native_record.structural_record

    @property
    def canonical_bytes(self) -> bytes:
        return _common_certificate_bytes(
            tag=b"exact-segment-containment-certificate-v1",
            domain_or_entry_digest=self.domain_digest,
            motion_or_entry=canonical_task1_bytes(self.motion),
            tool_radius=self.tool_radius,
            native_record=self.native_record,
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()


@dataclass(frozen=True)
class CircleContainmentCertificate:
    """Replay-bound proof that one full-circle cutter sweep lies in `D`."""

    domain_digest: bytes
    motion: ExactCircleMotion
    tool_radius: ToolRadius
    native_record: _containment_2.ContainmentRecord2

    def __post_init__(self) -> None:
        _digest(self.domain_digest, "containment domain digest")
        if type(self.motion) is not ExactCircleMotion:
            raise InvalidContainmentCertificateError("circle certificate requires one exact full-circle motion.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidContainmentCertificateError("circle certificate requires one typed tool radius.")
        record = _record(self.native_record)
        if not record.contained or not record.guide_anchor_in_center_domain:
            raise InvalidContainmentCertificateError("circle certificate does not contain one accepted exact verdict.")
        if not record.matches_exact_full_circle(
            self.domain_digest,
            self.motion.center.x,
            self.motion.center.y,
            self.motion.phase_vector.x,
            self.motion.phase_vector.y,
            self.tool_radius.value,
        ):
            raise InvalidContainmentCertificateError("circle native record contradicts its exact motion.")

    @property
    def guide_anchor_in_center_domain(self) -> bool:
        return self.native_record.guide_anchor_in_center_domain

    @property
    def outer_disk_contained(self) -> bool:
        return self.native_record.outer_disk_contained

    @property
    def disk_sweep(self) -> bool:
        return self.native_record.disk_sweep

    @property
    def native_strategy_version(self) -> bytes:
        return self.native_record.strategy_version

    @property
    def native_structural_record(self) -> bytes:
        return self.native_record.structural_record

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"exact-full-circle-containment-certificate-v1",
            encode_component_map(
                {
                    b"common": _common_certificate_bytes(
                        tag=b"exact-full-circle-containment-common-v1",
                        domain_or_entry_digest=self.domain_digest,
                        motion_or_entry=canonical_task1_bytes(self.motion),
                        tool_radius=self.tool_radius,
                        native_record=self.native_record,
                    ),
                    b"disk-sweep": encode_boolean(self.disk_sweep),
                    b"guide-anchor-in-center-domain": encode_boolean(self.guide_anchor_in_center_domain),
                    b"outer-disk-contained": encode_boolean(self.outer_disk_contained),
                }
            ),
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()


@dataclass(frozen=True)
class EntryContainmentCertificate:
    """Replay-bound proof that one complete entry disk lies in `D`."""

    domain_digest: bytes
    center: Point2[WorldXY]
    radius: EntryRadius
    tool_radius: ToolRadius
    native_record: _containment_2.ContainmentRecord2

    def __post_init__(self) -> None:
        _digest(self.domain_digest, "entry containment domain digest")
        if type(self.center) is not Point2:
            raise InvalidContainmentCertificateError("entry certificate requires one world-XY center.")
        if type(self.radius) is not EntryRadius or type(self.tool_radius) is not ToolRadius:
            raise InvalidContainmentCertificateError("entry certificate requires typed entry and tool radii.")
        record = _record(self.native_record)
        if not record.contained or not record.guide_anchor_in_center_domain:
            raise InvalidContainmentCertificateError("entry certificate does not contain one accepted exact verdict.")
        if not record.matches_exact_entry_disk(
            self.domain_digest,
            self.center.x,
            self.center.y,
            self.radius.value,
            self.tool_radius.value,
        ):
            raise InvalidContainmentCertificateError("entry native record contradicts its exact disk.")

    @property
    def center_in_center_domain(self) -> bool:
        return self.native_record.guide_anchor_in_center_domain

    @property
    def disk_contained(self) -> bool:
        return self.native_record.contained

    @property
    def native_strategy_version(self) -> bytes:
        return self.native_record.strategy_version

    @property
    def native_structural_record(self) -> bytes:
        return self.native_record.structural_record

    @property
    def canonical_bytes(self) -> bytes:
        return _common_certificate_bytes(
            tag=b"exact-entry-disk-containment-certificate-v1",
            domain_or_entry_digest=self.domain_digest,
            motion_or_entry=encode_tagged_union(
                b"entry-disk-v1",
                encode_component_map(
                    {
                        b"center": canonical_point2_bytes(self.center),
                        b"radius": encode_binary64(float(self.radius.value)),
                    }
                ),
            ),
            tool_radius=self.tool_radius,
            native_record=self.native_record,
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()


@dataclass(frozen=True)
class CircleInEntryCertificate:
    """Replay-bound proof that one circle sweep lies in a declared entry disk."""

    entry_digest: bytes
    entry_center: Point2[WorldXY]
    entry_radius: EntryRadius
    motion: ExactCircleMotion
    tool_radius: ToolRadius
    native_record: _containment_2.ContainmentRecord2

    def __post_init__(self) -> None:
        _digest(self.entry_digest, "entry identity digest")
        if type(self.entry_center) is not Point2 or type(self.entry_radius) is not EntryRadius:
            raise InvalidContainmentCertificateError("entry-circle certificate requires typed entry geometry.")
        if type(self.motion) is not ExactCircleMotion or type(self.tool_radius) is not ToolRadius:
            raise InvalidContainmentCertificateError("entry-circle certificate requires one exact circle and tool radius.")
        record = _record(self.native_record)
        if not record.contained:
            raise InvalidContainmentCertificateError("entry-circle certificate does not contain one accepted exact verdict.")
        if not record.matches_exact_full_circle_in_disk(
            self.entry_digest,
            self.entry_center.x,
            self.entry_center.y,
            self.entry_radius.value,
            self.motion.center.x,
            self.motion.center.y,
            self.motion.phase_vector.x,
            self.motion.phase_vector.y,
            self.tool_radius.value,
        ):
            raise InvalidContainmentCertificateError("entry-circle native record contradicts its exact geometry.")

    @property
    def native_strategy_version(self) -> bytes:
        return self.native_record.strategy_version

    @property
    def native_structural_record(self) -> bytes:
        return self.native_record.structural_record

    @property
    def canonical_bytes(self) -> bytes:
        return _common_certificate_bytes(
            tag=b"exact-full-circle-in-entry-certificate-v1",
            domain_or_entry_digest=self.entry_digest,
            motion_or_entry=encode_tagged_union(
                b"entry-and-circle-v1",
                encode_component_map(
                    {
                        b"entry-center": canonical_point2_bytes(self.entry_center),
                        b"entry-radius": encode_binary64(float(self.entry_radius.value)),
                        b"motion": canonical_task1_bytes(self.motion),
                    }
                ),
            ),
            tool_radius=self.tool_radius,
            native_record=self.native_record,
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()


class GougeContainment:
    """Single exact authority for segment, circle, and entry containment."""

    def __init__(self, reachable_domain: ReachableDomain) -> None:
        if type(reachable_domain) is not ReachableDomain:
            raise InvalidContainmentCertificateError("gouge containment requires one exact reachable-domain owner.")
        self._reachable_domain = reachable_domain

    @classmethod
    def build(
        cls,
        reachable_domain: ReachableDomain,
    ) -> Self:
        """Build a certifier bound to one exact reachable-domain owner.

        Args:
            reachable_domain: Exact owner of design and center-domain regions.

        Returns:
            Certifier that cannot accept another domain or tool radius.

        Raises:
            InvalidContainmentCertificateError: If the owner type is invalid.
        """
        return cls(reachable_domain)

    def _validate_tool_radius(self, tool_radius: ToolRadius) -> None:
        if type(tool_radius) is not ToolRadius:
            raise InvalidContainmentCertificateError("gouge containment requires one typed tool radius.")
        if tool_radius != self._reachable_domain.certificate.tool_radius:
            raise InvalidContainmentCertificateError("gouge containment tool radius contradicts its center domain.")

    def certify_segment(
        self,
        motion: ExactSegmentMotion,
        tool_radius: ToolRadius,
    ) -> SegmentContainmentCertificate:
        """Certify exact capsule containment for a segment center motion.

        Args:
            motion: Exact nondegenerate segment center motion.
            tool_radius: Radius owned by the reachable domain.

        Returns:
            Native-replay-bound segment containment certificate.

        Raises:
            InvalidContainmentCertificateError: If inputs are cross-wired.
            GougeContainmentError: If the capsule is not a subset of `D`.
        """
        if type(motion) is not ExactSegmentMotion:
            raise InvalidContainmentCertificateError("segment containment requires one exact segment motion.")
        self._validate_tool_radius(tool_radius)
        native_record = _containment_2.evaluate_exact_segment_containment(
            self._reachable_domain.design_region,
            self._reachable_domain.center_domain,
            self._reachable_domain.certificate.digest,
            motion.start.x,
            motion.start.y,
            motion.end.x,
            motion.end.y,
            tool_radius.value,
        )
        if not native_record.contained:
            raise GougeContainmentError("segment cutter sweep is not exactly contained in the design pocket.")
        return SegmentContainmentCertificate(
            self._reachable_domain.certificate.digest,
            motion,
            tool_radius,
            native_record,
        )

    def certify_full_circle(
        self,
        motion: ExactCircleMotion,
        tool_radius: ToolRadius,
    ) -> CircleContainmentCertificate:
        """Certify exact annular or disk containment for one full circle.

        Args:
            motion: Exact nondegenerate full-circle center motion.
            tool_radius: Radius owned by the reachable domain.

        Returns:
            Native-replay-bound full-circle containment certificate.

        Raises:
            InvalidContainmentCertificateError: If inputs are cross-wired.
            GougeContainmentError: If the sweep is not a subset of `D`.
        """
        if type(motion) is not ExactCircleMotion:
            raise InvalidContainmentCertificateError("circle containment requires one exact full-circle motion.")
        self._validate_tool_radius(tool_radius)
        native_record = _containment_2.evaluate_exact_full_circle_containment(
            self._reachable_domain.design_region,
            self._reachable_domain.center_domain,
            self._reachable_domain.certificate.digest,
            motion.center.x,
            motion.center.y,
            motion.phase_vector.x,
            motion.phase_vector.y,
            tool_radius.value,
        )
        if not native_record.contained:
            raise GougeContainmentError("full-circle cutter sweep is not exactly contained in the design pocket.")
        return CircleContainmentCertificate(
            self._reachable_domain.certificate.digest,
            motion,
            tool_radius,
            native_record,
        )

    def certify_entry_disk(
        self,
        *,
        center: Point2[WorldXY],
        radius: EntryRadius,
        tool_radius: ToolRadius,
    ) -> EntryContainmentCertificate:
        """Certify a complete precleared disk in the design pocket.

        Args:
            center: World-XY center of the declared disk.
            radius: Typed entry-disk radius.
            tool_radius: Radius owned by the reachable domain.

        Returns:
            Native-replay-bound entry containment certificate.

        Raises:
            InvalidContainmentCertificateError: If inputs are cross-wired.
            GougeContainmentError: If the disk is not a subset of `D`.
        """
        if type(center) is not Point2 or type(radius) is not EntryRadius:
            raise InvalidContainmentCertificateError("entry containment requires typed center and radius.")
        self._validate_tool_radius(tool_radius)
        native_record = _containment_2.evaluate_exact_entry_disk_containment(
            self._reachable_domain.design_region,
            self._reachable_domain.center_domain,
            self._reachable_domain.certificate.digest,
            center.x,
            center.y,
            radius.value,
            tool_radius.value,
        )
        if not native_record.contained:
            raise GougeContainmentError("entry disk is not exactly contained in the design pocket.")
        return EntryContainmentCertificate(
            self._reachable_domain.certificate.digest,
            center,
            radius,
            tool_radius,
            native_record,
        )

    def certify_full_circle_in_disk(
        self,
        *,
        entry_digest: bytes,
        entry_center: Point2[WorldXY],
        entry_radius: EntryRadius,
        motion: ExactCircleMotion,
        tool_radius: ToolRadius,
    ) -> CircleInEntryCertificate:
        """Certify a complete full-circle cutter sweep in an entry disk.

        Args:
            entry_digest: SHA-256 identity of the owning precleared entry.
            entry_center: World-XY center of the precleared disk.
            entry_radius: Radius of the precleared disk.
            motion: First exact full-circle center motion.
            tool_radius: Radius owned by the reachable domain.

        Returns:
            Certificate bound to the entry authority and circle geometry.

        Raises:
            InvalidContainmentCertificateError: If inputs are cross-wired.
            GougeContainmentError: If the sweep leaves the entry disk.
        """
        _digest(entry_digest, "entry identity digest")
        if type(entry_center) is not Point2 or type(entry_radius) is not EntryRadius:
            raise InvalidContainmentCertificateError("entry-circle containment requires typed entry geometry.")
        if type(motion) is not ExactCircleMotion:
            raise InvalidContainmentCertificateError("entry-circle containment requires one exact full-circle motion.")
        self._validate_tool_radius(tool_radius)
        native_record = _containment_2.evaluate_exact_full_circle_in_disk(
            entry_digest,
            entry_center.x,
            entry_center.y,
            entry_radius.value,
            motion.center.x,
            motion.center.y,
            motion.phase_vector.x,
            motion.phase_vector.y,
            tool_radius.value,
        )
        if not native_record.contained:
            raise GougeContainmentError("full-circle cutter sweep is not exactly contained in the entry disk.")
        return CircleInEntryCertificate(
            entry_digest,
            entry_center,
            entry_radius,
            motion,
            tool_radius,
            native_record,
        )
