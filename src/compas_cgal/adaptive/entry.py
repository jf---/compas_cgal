import hashlib
from dataclasses import dataclass
from fractions import Fraction
from typing import Final
from typing import NewType
from typing import Self

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import canonical_point2_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.containment import CircleInEntryCertificate
from compas_cgal.adaptive.containment import EntryContainmentCertificate
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.errors import InvalidPreclearedEntryError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

ENTRY_SCHEMA_VERSION: Final[bytes] = b"precleared-entry-schema-v1"
ENTRY_DEPLETION_STRATEGY_VERSION: Final[bytes] = _stock_2.exact_entry_depletion_strategy_version()

BoreProcessIdentity = NewType("BoreProcessIdentity", bytes)
BoreEvidenceDigest = NewType("BoreEvidenceDigest", bytes)

_SHA256_BYTES = hashlib.sha256().digest_size


def _nonempty_bytes(value: object, name: str) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidPreclearedEntryError(f"{name} must be nonempty exact bytes.")
    return value


@dataclass(frozen=True)
class QualifiedBore:
    """Process evidence for one bore over a complete vertical cut interval."""

    cut_plane: CutPlane
    process_identity: BoreProcessIdentity
    evidence_bytes: bytes
    evidence_digest: BoreEvidenceDigest

    def __post_init__(self) -> None:
        if type(self.cut_plane) is not CutPlane:
            raise InvalidPreclearedEntryError("qualified bore requires one exact cut-plane interval.")
        _nonempty_bytes(self.process_identity, "bore process identity")
        _nonempty_bytes(self.evidence_bytes, "bore qualification evidence")
        if type(self.evidence_digest) is not bytes or len(self.evidence_digest) != _SHA256_BYTES:
            raise InvalidPreclearedEntryError("bore evidence digest must be one SHA-256 digest.")
        if self.evidence_digest != hashlib.sha256(self.evidence_bytes).digest():
            raise InvalidPreclearedEntryError("bore evidence digest contradicts its qualification bytes.")

    @classmethod
    def build(
        cls,
        *,
        cut_plane: CutPlane,
        process_identity: BoreProcessIdentity,
        evidence_bytes: bytes,
    ) -> Self:
        """Build a bore qualification and derive its evidence digest.

        Args:
            cut_plane: Complete clearance-to-cut vertical interval.
            process_identity: Nonempty identity of the qualified bore process.
            evidence_bytes: Immutable qualification or metrology evidence.

        Returns:
            Validated bore qualification with SHA-256-addressed evidence.

        Raises:
            InvalidPreclearedEntryError: If any qualification input is invalid.
        """
        if type(cut_plane) is not CutPlane:
            raise InvalidPreclearedEntryError("bore qualification factory requires one exact cut plane.")
        process = BoreProcessIdentity(_nonempty_bytes(process_identity, "bore process identity"))
        evidence = _nonempty_bytes(evidence_bytes, "bore qualification evidence")
        return cls(
            cut_plane,
            process,
            evidence,
            BoreEvidenceDigest(hashlib.sha256(evidence).digest()),
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"qualified-bore-v1",
            encode_component_map(
                {
                    b"cut-plane": canonical_task1_bytes(self.cut_plane),
                    b"evidence": encode_bytes(self.evidence_bytes),
                    b"evidence-digest": encode_bytes(bytes(self.evidence_digest)),
                    b"process-identity": encode_bytes(bytes(self.process_identity)),
                }
            ),
        )


@dataclass(frozen=True)
class PreclearedEntry:
    """Exact precleared disk and physical bore qualification for path entry."""

    reachable_domain: ReachableDomain
    center: Point2[WorldXY]
    radius: EntryRadius
    tool_radius: ToolRadius
    cut_plane: CutPlane
    qualified_bore: QualifiedBore
    containment: EntryContainmentCertificate

    def __post_init__(self) -> None:
        if type(self.reachable_domain) is not ReachableDomain:
            raise InvalidPreclearedEntryError("precleared entry requires one exact reachable-domain owner.")
        if type(self.center) is not Point2:
            raise InvalidPreclearedEntryError("precleared entry requires one world-XY center.")
        if type(self.radius) is not EntryRadius or type(self.tool_radius) is not ToolRadius:
            raise InvalidPreclearedEntryError("precleared entry requires typed entry and tool radii.")
        if Fraction.from_float(self.radius.value) <= Fraction.from_float(self.tool_radius.value):
            raise InvalidPreclearedEntryError("precleared entry radius must exceed the tool radius to launch nondegenerate lateral motion.")
        if type(self.cut_plane) is not CutPlane or type(self.qualified_bore) is not QualifiedBore:
            raise InvalidPreclearedEntryError("precleared entry requires one cut plane and qualified bore.")
        if self.qualified_bore.cut_plane != self.cut_plane:
            raise InvalidPreclearedEntryError("qualified interval must equal the complete clearance-to-cut interval.")
        if self.reachable_domain.certificate.tool_radius != self.tool_radius:
            raise InvalidPreclearedEntryError("entry tool radius contradicts its reachable domain.")
        if type(self.containment) is not EntryContainmentCertificate:
            raise InvalidPreclearedEntryError("precleared entry requires one exact disk-containment certificate.")
        if (
            self.containment.domain_digest != self.reachable_domain.certificate.digest
            or self.containment.center != self.center
            or self.containment.radius != self.radius
            or self.containment.tool_radius != self.tool_radius
        ):
            raise InvalidPreclearedEntryError("entry containment certificate contradicts its owned input.")

    @classmethod
    def build(
        cls,
        *,
        reachable_domain: ReachableDomain,
        center: Point2[WorldXY],
        radius: EntryRadius,
        tool_radius: ToolRadius,
        cut_plane: CutPlane,
        qualified_bore: QualifiedBore,
    ) -> Self:
        """Build and exactly certify one precleared entry disk.

        Args:
            reachable_domain: Exact owner of `D` and admissible centers `C_r`.
            center: World-XY center of the qualified bore.
            radius: Radius of the complete precleared disk.
            tool_radius: Tool radius bound to the reachable domain.
            cut_plane: Complete clearance-to-cut vertical interval.
            qualified_bore: Process evidence for that same interval.

        Returns:
            Entry with an exact disk-containment certificate.

        Raises:
            InvalidPreclearedEntryError: If entry inputs are cross-wired or
                cannot launch nondegenerate lateral motion.
            GougeContainmentError: If the declared disk is not contained in
                the design pocket.
        """
        if type(reachable_domain) is not ReachableDomain:
            raise InvalidPreclearedEntryError("entry factory requires one exact reachable domain.")
        if type(center) is not Point2:
            raise InvalidPreclearedEntryError("entry factory requires one world-XY center.")
        if type(radius) is not EntryRadius or type(tool_radius) is not ToolRadius:
            raise InvalidPreclearedEntryError("entry factory requires typed entry and tool radii.")
        if Fraction.from_float(radius.value) <= Fraction.from_float(tool_radius.value):
            raise InvalidPreclearedEntryError("entry cannot launch from a one-tool-radius or smaller seed.")
        if type(cut_plane) is not CutPlane or type(qualified_bore) is not QualifiedBore:
            raise InvalidPreclearedEntryError("entry factory requires one cut plane and qualified bore.")
        if qualified_bore.cut_plane != cut_plane:
            raise InvalidPreclearedEntryError("qualified interval must equal the complete clearance-to-cut interval.")
        containment = GougeContainment.build(reachable_domain).certify_entry_disk(
            center=center,
            radius=radius,
            tool_radius=tool_radius,
        )
        return cls(
            reachable_domain,
            center,
            radius,
            tool_radius,
            cut_plane,
            qualified_bore,
            containment,
        )

    @property
    def approach(self) -> ApproachOperation:
        return ApproachOperation.build(
            position=self.center,
            clearance_z=self.cut_plane.clearance_z,
        )

    @property
    def plunge(self) -> PlungeOperation:
        return PlungeOperation.build(
            position=self.center,
            clearance_z=self.cut_plane.clearance_z,
            cut_z=self.cut_plane.cut_z,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            ENTRY_SCHEMA_VERSION,
            encode_component_map(
                {
                    b"center": canonical_point2_bytes(self.center),
                    b"containment-certificate": self.containment.canonical_bytes,
                    b"entry-radius": encode_binary64(float(self.radius.value)),
                    b"qualified-bore": self.qualified_bore.canonical_bytes,
                    b"reachable-domain-digest": encode_bytes(self.reachable_domain.certificate.digest),
                    b"tool-radius": encode_binary64(float(self.tool_radius.value)),
                }
            ),
        )

    @property
    def digest(self) -> bytes:
        return hashlib.sha256(self.canonical_bytes).digest()

    def certify_first_circle(
        self,
        motion: ExactCircleMotion,
    ) -> CircleInEntryCertificate:
        """Certify the first circle's complete cutter sweep in the entry void.

        Args:
            motion: Exact full-circle center motion to launch.

        Returns:
            Certificate bound to this entry digest and circle geometry.

        Raises:
            InvalidPreclearedEntryError: If `motion` is not an exact circle.
            GougeContainmentError: If the sweep leaves the precleared disk.
        """
        if type(motion) is not ExactCircleMotion:
            raise InvalidPreclearedEntryError("entry launch requires one exact full-circle motion.")
        return GougeContainment.build(self.reachable_domain).certify_full_circle_in_disk(
            entry_digest=self.digest,
            entry_center=self.center,
            entry_radius=self.radius,
            motion=motion,
            tool_radius=self.tool_radius,
        )
