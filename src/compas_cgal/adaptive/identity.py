from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from math import gcd
from typing import TYPE_CHECKING
from typing import Final
from typing import NewType
from typing import Self
from typing import TypeAlias

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import CANONICAL_ENCODING_VERSION
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import canonical_polygon_bytes
from compas_cgal.adaptive.canonical import canonical_record_kind
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import CanonicalEncodingError
from compas_cgal.adaptive.errors import InvalidBoundaryVertexIdentityError
from compas_cgal.adaptive.errors import InvalidComponentIdentityError
from compas_cgal.adaptive.errors import InvalidInputIdentityError

if TYPE_CHECKING:
    from compas_cgal.adaptive.entry import PreclearedEntry
    from compas_cgal.adaptive.motion import EngagementCap
    from compas_cgal.adaptive.policy import CandidatePolicy
    from compas_cgal.adaptive.policy import CutDirectionPolicy
    from compas_cgal.adaptive.policy import DepletionPolicy
    from compas_cgal.adaptive.policy import NeckPolicy
    from compas_cgal.adaptive.policy import TraversalPolicy
    from compas_cgal.adaptive.reachable_domain import ReachableDomain
    from compas_cgal.adaptive.units import CutPlane
    from compas_cgal.adaptive.units import ToolRadius

INPUT_SCHEMA_VERSION: Final[bytes] = b"adaptive-input-schema-v1"
OPERATION_SCHEMA_VERSION: Final[bytes] = b"adaptive-operation-schema-v1"
COMPONENT_IDENTITY_VERSION: Final[bytes] = b"component-identity-v1"
BOUNDARY_VERTEX_ID_VERSION: Final[bytes] = b"boundary-vertex-id-v1"
WORLD_XY_MILLIMETRE_FRAME_VERSION: Final[bytes] = b"world-xy-millimetre-frame-v1"

IdentityDigest = NewType("IdentityDigest", bytes)
FrameIdentity = NewType("FrameIdentity", bytes)
ComponentDomainTag = NewType("ComponentDomainTag", bytes)
StrategyVersion = NewType("StrategyVersion", bytes)
SourceRevision = NewType("SourceRevision", bytes)
NativeSourceTreeDigest = NewType("NativeSourceTreeDigest", bytes)


def _nonempty_bytes(value: object, name: str) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidComponentIdentityError(f"{name} must be nonempty bytes.")
    return value


@dataclass(frozen=True)
class ComponentIdentity:
    component_domain: ComponentDomainTag
    strategy_version: StrategyVersion
    source_revision: SourceRevision
    native_source_tree_digest: NativeSourceTreeDigest
    canonical_parameter_bytes: bytes

    def __post_init__(self) -> None:
        _nonempty_bytes(self.component_domain, "component domain")
        _nonempty_bytes(self.strategy_version, "strategy version")
        _nonempty_bytes(self.source_revision, "source revision")
        if type(self.native_source_tree_digest) is not bytes or len(self.native_source_tree_digest) != hashlib.sha256().digest_size:
            raise InvalidComponentIdentityError("native source-tree digest must be exactly one 32-byte SHA-256 digest.")
        try:
            require_canonical_record(self.canonical_parameter_bytes)
        except CanonicalEncodingError:
            raise InvalidComponentIdentityError("canonical parameter bytes must contain one complete versioned CCAN record.") from None
        if canonical_record_kind(self.canonical_parameter_bytes) != b"T":
            raise InvalidComponentIdentityError("canonical parameter bytes must be a tagged domain record.")

    @classmethod
    def build(
        cls,
        *,
        component_domain: ComponentDomainTag,
        strategy_version: StrategyVersion,
        source_revision: SourceRevision,
        native_source_tree_digest: NativeSourceTreeDigest,
        canonical_parameter_bytes: bytes,
    ) -> Self:
        return cls(
            component_domain,
            strategy_version,
            source_revision,
            native_source_tree_digest,
            canonical_parameter_bytes,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not ComponentIdentity:
            raise InvalidComponentIdentityError("component identity must be exact ComponentIdentity, not a subclass.")
        return encode_tagged_union(
            COMPONENT_IDENTITY_VERSION,
            encode_component_map(
                {
                    b"canonical-parameters": self.canonical_parameter_bytes,
                    b"component-domain": bytes(self.component_domain),
                    b"native-source-tree-digest": bytes(self.native_source_tree_digest),
                    b"source-revision": bytes(self.source_revision),
                    b"strategy-version": bytes(self.strategy_version),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class ComponentVersionBinding:
    """Canonical name/version pair for one active planning component."""

    component: bytes
    version: bytes

    def __post_init__(self) -> None:
        if type(self.component) is not bytes or not self.component:
            raise InvalidInputIdentityError("component-version binding requires one nonempty component name.")
        if type(self.version) is not bytes or not self.version:
            raise InvalidInputIdentityError("component-version binding requires one nonempty version.")

    @classmethod
    def build(
        cls,
        *,
        component: bytes,
        version: bytes,
    ) -> Self:
        """Build one nonempty component-version binding.

        Args:
            component: Stable component domain name.
            version: Active schema or strategy version.

        Returns:
            Validated immutable binding.

        Raises:
            InvalidInputIdentityError: If either byte string is empty.
        """
        return cls(component, version)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not ComponentVersionBinding:
            raise InvalidInputIdentityError("component version must be exact ComponentVersionBinding.")
        return encode_tagged_union(
            b"component-version-binding-v1",
            encode_component_map(
                {
                    b"component": encode_bytes(self.component),
                    b"version": encode_bytes(self.version),
                }
            ),
        )


def _current_component_versions(
    reachable_domain: ReachableDomain,
    entry: PreclearedEntry,
) -> tuple[ComponentVersionBinding, ...]:
    from compas_cgal.adaptive.entry import ENTRY_DEPLETION_STRATEGY_VERSION
    from compas_cgal.adaptive.entry import ENTRY_SCHEMA_VERSION
    from compas_cgal.adaptive.motion_certificate import MOTION_CERTIFICATE_SCHEMA_VERSION

    bindings = (
        ComponentVersionBinding.build(
            component=b"canonical-encoding",
            version=CANONICAL_ENCODING_VERSION,
        ),
        ComponentVersionBinding.build(
            component=b"input-schema",
            version=INPUT_SCHEMA_VERSION,
        ),
        ComponentVersionBinding.build(
            component=b"operation-schema",
            version=OPERATION_SCHEMA_VERSION,
        ),
        ComponentVersionBinding.build(
            component=b"reachable-domain",
            version=reachable_domain.certificate.strategy_version,
        ),
        ComponentVersionBinding.build(
            component=b"gouge-containment",
            version=entry.containment.native_strategy_version,
        ),
        ComponentVersionBinding.build(
            component=b"precleared-entry-schema",
            version=ENTRY_SCHEMA_VERSION,
        ),
        ComponentVersionBinding.build(
            component=b"precleared-entry-depletion",
            version=ENTRY_DEPLETION_STRATEGY_VERSION,
        ),
        ComponentVersionBinding.build(
            component=b"motion-depletion",
            version=_stock_2.exact_depletion_strategy_version(),
        ),
        ComponentVersionBinding.build(
            component=b"motion-certificate-schema",
            version=MOTION_CERTIFICATE_SCHEMA_VERSION,
        ),
        ComponentVersionBinding.build(
            component=b"event-exact-motion-oracle",
            version=_continuous_tea_2.event_oracle_component_version(),
        ),
    )
    return tuple(sorted(bindings, key=lambda binding: binding.component))


@dataclass(frozen=True)
class InputIdentity:
    """Content-addressed root of all validated Phase-1 planning inputs."""

    design_boundary: CanonicalRingV1
    holes: tuple[CanonicalRingV1, ...]
    frame_identity: FrameIdentity
    cut_plane: CutPlane
    tool_radius: ToolRadius
    reachable_domain: ReachableDomain
    reachable_domain_digest: bytes
    entry: PreclearedEntry
    user_cap: EngagementCap
    candidate_policy: CandidatePolicy
    neck_policy: NeckPolicy
    depletion_policy: DepletionPolicy
    traversal_policy: TraversalPolicy
    cut_direction_policy: CutDirectionPolicy
    component_versions: tuple[ComponentVersionBinding, ...]

    def __post_init__(self) -> None:
        from compas_cgal.adaptive.entry import PreclearedEntry
        from compas_cgal.adaptive.motion import EngagementCap
        from compas_cgal.adaptive.policy import CandidatePolicy
        from compas_cgal.adaptive.policy import CutDirectionPolicy
        from compas_cgal.adaptive.policy import DepletionPolicy
        from compas_cgal.adaptive.policy import NeckPolicy
        from compas_cgal.adaptive.policy import TraversalPolicy
        from compas_cgal.adaptive.reachable_domain import ReachableDomain
        from compas_cgal.adaptive.units import CutPlane
        from compas_cgal.adaptive.units import ToolRadius

        if type(self.design_boundary) is not CanonicalRingV1 or not self.design_boundary.is_outer:
            raise InvalidInputIdentityError("input identity requires one canonical outer design boundary.")
        if type(self.holes) is not tuple or any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in self.holes):
            raise InvalidInputIdentityError("input identity holes must be one canonical hole tuple.")
        hole_bytes = tuple(hole.canonical_bytes for hole in self.holes)
        if hole_bytes != tuple(sorted(hole_bytes)) or len(hole_bytes) != len(set(hole_bytes)):
            raise InvalidInputIdentityError("input identity holes must be sorted and unique.")
        if type(self.frame_identity) is not bytes or self.frame_identity != WORLD_XY_MILLIMETRE_FRAME_VERSION:
            raise InvalidInputIdentityError("input identity requires the world-XY millimetre frame.")
        if type(self.cut_plane) is not CutPlane:
            raise InvalidInputIdentityError("input identity requires one typed cut plane.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidInputIdentityError("input identity requires one typed tool radius.")
        if type(self.reachable_domain) is not ReachableDomain:
            raise InvalidInputIdentityError("input identity requires one exact reachable domain.")
        certificate = self.reachable_domain.certificate
        if certificate.design_boundary != self.design_boundary or certificate.holes != self.holes or certificate.tool_radius != self.tool_radius:
            raise InvalidInputIdentityError("reachable domain contradicts the owned design geometry or tool.")
        if (
            type(self.reachable_domain_digest) is not bytes
            or len(self.reachable_domain_digest) != hashlib.sha256().digest_size
            or self.reachable_domain_digest != certificate.digest
        ):
            raise InvalidInputIdentityError("reachable-domain digest contradicts its exact certificate.")
        if type(self.entry) is not PreclearedEntry:
            raise InvalidInputIdentityError("input identity entry must be exact PreclearedEntry.")
        if self.entry.reachable_domain is not self.reachable_domain or self.entry.cut_plane != self.cut_plane or self.entry.tool_radius != self.tool_radius:
            raise InvalidInputIdentityError("PreclearedEntry contradicts its input owner, cut plane, or tool.")
        if type(self.user_cap) is not EngagementCap:
            raise InvalidInputIdentityError("input identity requires one exact engagement cap.")
        for value, expected, name in (
            (self.candidate_policy, CandidatePolicy, "candidate"),
            (self.neck_policy, NeckPolicy, "neck"),
            (self.depletion_policy, DepletionPolicy, "depletion"),
            (self.traversal_policy, TraversalPolicy, "traversal"),
            (
                self.cut_direction_policy,
                CutDirectionPolicy,
                "cut-direction",
            ),
        ):
            if type(value) is not expected:
                raise InvalidInputIdentityError(f"input identity {name} policy has the wrong exact type.")
        if self.neck_policy.user_cap != self.user_cap:
            raise InvalidInputIdentityError("neck policy contradicts the owned user engagement cap.")
        if type(self.component_versions) is not tuple or any(type(binding) is not ComponentVersionBinding for binding in self.component_versions):
            raise InvalidInputIdentityError("input identity requires exact component-version bindings.")
        expected_versions = _current_component_versions(
            self.reachable_domain,
            self.entry,
        )
        if self.component_versions != expected_versions:
            raise InvalidInputIdentityError("component versions do not match the active exact pipeline.")
        try:
            self.canonical_bytes
        except CanonicalEncodingError as error:
            raise InvalidInputIdentityError(str(error)) from None

    @classmethod
    def build(
        cls,
        *,
        design_boundary: CanonicalRingV1,
        holes: Sequence[CanonicalRingV1],
        cut_plane: CutPlane,
        tool_radius: ToolRadius,
        reachable_domain: ReachableDomain,
        entry: PreclearedEntry,
        user_cap: EngagementCap,
        candidate_policy: CandidatePolicy,
        neck_policy: NeckPolicy,
        depletion_policy: DepletionPolicy,
        traversal_policy: TraversalPolicy,
        cut_direction_policy: CutDirectionPolicy,
    ) -> Self:
        """Build the complete canonical identity root for a planning run.

        Args:
            design_boundary: Canonical outer design ring.
            holes: Canonical island rings in any input order.
            cut_plane: Typed clearance-to-cut interval.
            tool_radius: Tool radius owned by the reachable domain.
            reachable_domain: Exact design/center-domain proof owner.
            entry: Validated precleared entry owned by that domain.
            user_cap: Maximum lateral engagement cap.
            candidate_policy: Complete finite proposal-lattice policy.
            neck_policy: Exact neck classification and effective-cap policy.
            depletion_policy: Exact motion-depletion construction policy.
            traversal_policy: Finite traversal lookahead policy.
            cut_direction_policy: Explicit climb/conventional intent.

        Returns:
            Identity with canonical holes, current component versions, and a
            derived reachable-domain digest.

        Raises:
            InvalidInputIdentityError: If any owner or policy is malformed or
                cross-wired.
        """
        from compas_cgal.adaptive.entry import PreclearedEntry
        from compas_cgal.adaptive.reachable_domain import ReachableDomain

        try:
            canonical_holes = tuple(holes)
        except TypeError:
            raise InvalidInputIdentityError("input identity holes must be a finite canonical sequence.") from None
        if any(type(hole) is not CanonicalRingV1 or hole.is_outer for hole in canonical_holes):
            raise InvalidInputIdentityError("input identity holes must be canonical hole rings.")
        canonical_holes = tuple(sorted(canonical_holes, key=lambda hole: hole.canonical_bytes))
        if len(canonical_holes) != len({hole.canonical_bytes for hole in canonical_holes}):
            raise InvalidInputIdentityError("input identity holes must be canonically unique.")
        if type(reachable_domain) is not ReachableDomain or type(entry) is not PreclearedEntry:
            raise InvalidInputIdentityError("input identity factory requires exact reachable and entry owners.")
        versions = _current_component_versions(
            reachable_domain,
            entry,
        )
        return cls(
            design_boundary,
            canonical_holes,
            FrameIdentity(WORLD_XY_MILLIMETRE_FRAME_VERSION),
            cut_plane,
            tool_radius,
            reachable_domain,
            reachable_domain.certificate.digest,
            entry,
            user_cap,
            candidate_policy,
            neck_policy,
            depletion_policy,
            traversal_policy,
            cut_direction_policy,
            versions,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not InputIdentity:
            raise InvalidInputIdentityError("input identity must be exact InputIdentity, not a subclass.")
        return encode_tagged_union(
            INPUT_SCHEMA_VERSION,
            encode_component_map(
                {
                    b"candidate-policy": canonical_task1_bytes(self.candidate_policy),
                    b"component-versions": encode_sequence(tuple(binding.canonical_bytes for binding in self.component_versions)),
                    b"cut-direction-policy": canonical_task1_bytes(self.cut_direction_policy),
                    b"cut-plane": canonical_task1_bytes(self.cut_plane),
                    b"depletion-policy": canonical_task1_bytes(self.depletion_policy),
                    b"design": canonical_polygon_bytes(
                        self.design_boundary.vertices,
                        tuple(hole.vertices for hole in self.holes),
                    ),
                    b"entry": self.entry.canonical_bytes,
                    b"frame": encode_bytes(bytes(self.frame_identity)),
                    b"neck-policy": canonical_task1_bytes(self.neck_policy),
                    b"operation-schema-version": encode_bytes(OPERATION_SCHEMA_VERSION),
                    b"reachable-domain-digest": encode_bytes(self.reachable_domain_digest),
                    b"tool-radius": encode_tagged_union(
                        b"tool-radius-mm-v1",
                        encode_binary64(float(self.tool_radius.value)),
                    ),
                    b"traversal-policy": canonical_task1_bytes(self.traversal_policy),
                    b"user-cap": canonical_task1_bytes(self.user_cap),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


class SupportKind(Enum):
    LINE = b"line-support-v1"
    CIRCLE = b"circle-support-v1"


class TrimIncidenceOrientation(Enum):
    ENTERING = b"trim-entering-v1"
    LEAVING = b"trim-leaving-v1"


def _least_common_multiple(left: int, right: int) -> int:
    return abs(left * right) // gcd(left, right)


def _primitive_coefficients(coefficients: tuple[ExactRationalV1, ...]) -> tuple[int, ...]:
    denominator_lcm = 1
    for exact_coefficient in coefficients:
        denominator_lcm = _least_common_multiple(denominator_lcm, exact_coefficient.denominator)
    integers = tuple(coefficient.numerator * (denominator_lcm // coefficient.denominator) for coefficient in coefficients)
    common_factor = 0
    for integer_coefficient in integers:
        common_factor = gcd(common_factor, abs(integer_coefficient))
    if common_factor == 0:
        raise InvalidBoundaryVertexIdentityError("incident support coefficients must include a nonzero exact value.")
    primitive = tuple(coefficient // common_factor for coefficient in integers)
    first_nonzero = next(coefficient for coefficient in primitive if coefficient != 0)
    if first_nonzero < 0:
        primitive = tuple(-coefficient for coefficient in primitive)
    return primitive


@dataclass(frozen=True)
class IncidentSupportIdV1:
    support_kind: SupportKind
    primitive_coefficients: tuple[int, ...]

    def __post_init__(self) -> None:
        if type(self.support_kind) is not SupportKind:
            raise InvalidBoundaryVertexIdentityError("incident support kind must be exact SupportKind.")
        if type(self.primitive_coefficients) is not tuple or any(type(value) is not int for value in self.primitive_coefficients):
            raise InvalidBoundaryVertexIdentityError("incident support coefficients must be exact primitive integers.")
        expected_count = 3 if self.support_kind is SupportKind.LINE else 4
        if len(self.primitive_coefficients) != expected_count:
            expected_name = "three" if expected_count == 3 else "four"
            raise InvalidBoundaryVertexIdentityError(f"{self.support_kind.name.lower()} support requires exactly {expected_name} coefficients.")
        normalized = _primitive_coefficients(tuple(ExactRationalV1(value, 1) for value in self.primitive_coefficients))
        if normalized != self.primitive_coefficients:
            raise InvalidBoundaryVertexIdentityError("incident support coefficients must already be primitive with positive leading sign.")
        if self.support_kind is SupportKind.LINE and self.primitive_coefficients[:2] == (0, 0):
            raise InvalidBoundaryVertexIdentityError("line support requires a nonzero normal.")
        if self.support_kind is SupportKind.CIRCLE and self.primitive_coefficients[0] == 0:
            raise InvalidBoundaryVertexIdentityError("circle support requires a nonzero quadratic coefficient.")

    @classmethod
    def build(
        cls,
        *,
        support_kind: SupportKind,
        normalized_coefficients: Sequence[ExactRationalV1],
    ) -> Self:
        try:
            coefficients = tuple(normalized_coefficients)
        except TypeError:
            raise InvalidBoundaryVertexIdentityError("incident support coefficients must be finite.") from None
        if type(support_kind) is not SupportKind:
            raise InvalidBoundaryVertexIdentityError("incident support kind must be exact SupportKind.")
        if any(type(coefficient) is not ExactRationalV1 for coefficient in coefficients):
            raise InvalidBoundaryVertexIdentityError("incident support coefficients must be exact ExactRationalV1 values.")
        expected_count = 3 if support_kind is SupportKind.LINE else 4 if support_kind is SupportKind.CIRCLE else 0
        if len(coefficients) != expected_count:
            expected_name = "three" if expected_count == 3 else "four"
            raise InvalidBoundaryVertexIdentityError(f"incident support requires exactly {expected_name} exact coefficients.")
        return cls(support_kind, _primitive_coefficients(coefficients))

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not IncidentSupportIdV1:
            raise InvalidBoundaryVertexIdentityError("incident support identity must be exact IncidentSupportIdV1, not a subclass.")
        return encode_tagged_union(
            b"incident-support-id-v1",
            encode_component_map(
                {
                    b"coefficients": encode_sequence(tuple(encode_integer(coefficient) for coefficient in self.primitive_coefficients)),
                    b"kind": encode_tagged_union(bytes(self.support_kind.value), b""),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class IncidentSupport:
    support_id: IncidentSupportIdV1
    trim_incidence_orientation: TrimIncidenceOrientation

    def __post_init__(self) -> None:
        if type(self.support_id) is not IncidentSupportIdV1:
            raise InvalidBoundaryVertexIdentityError("incident support requires exact IncidentSupportIdV1.")
        if type(self.trim_incidence_orientation) is not TrimIncidenceOrientation:
            raise InvalidBoundaryVertexIdentityError("incident support requires exact TrimIncidenceOrientation.")

    @classmethod
    def build(
        cls,
        *,
        support_id: IncidentSupportIdV1,
        trim_incidence_orientation: TrimIncidenceOrientation,
    ) -> Self:
        return cls(support_id, trim_incidence_orientation)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not IncidentSupport:
            raise InvalidBoundaryVertexIdentityError("incident support must be exact IncidentSupport, not a subclass.")
        return encode_tagged_union(
            b"incident-support-v1",
            encode_component_map(
                {
                    b"orientation": encode_tagged_union(
                        bytes(self.trim_incidence_orientation.value),
                        b"",
                    ),
                    b"support-id": self.support_id.canonical_bytes,
                }
            ),
        )


@dataclass(frozen=True)
class InputRingVertexIdV1:
    canonical_ring: CanonicalRingV1
    vertex_ordinal: int

    def __post_init__(self) -> None:
        if type(self.canonical_ring) is not CanonicalRingV1:
            raise InvalidBoundaryVertexIdentityError("input-ring vertex requires exact CanonicalRingV1.")
        if type(self.vertex_ordinal) is not int or self.vertex_ordinal < 0:
            raise InvalidBoundaryVertexIdentityError("input-ring vertex ordinal must be an exact nonnegative integer.")
        if self.vertex_ordinal >= self.canonical_ring.vertex_count:
            raise InvalidBoundaryVertexIdentityError("input-ring vertex ordinal must be below the ring vertex count.")

    @classmethod
    def build(cls, *, canonical_ring: CanonicalRingV1, vertex_ordinal: int) -> Self:
        return cls(canonical_ring, vertex_ordinal)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not InputRingVertexIdV1:
            raise InvalidBoundaryVertexIdentityError("input-ring vertex identity must be exact InputRingVertexIdV1, not a subclass.")
        return encode_tagged_union(
            BOUNDARY_VERTEX_ID_VERSION,
            encode_tagged_union(
                b"input-ring-vertex-v1",
                encode_component_map(
                    {
                        b"ring": self.canonical_ring.canonical_bytes,
                        b"vertex-ordinal": encode_integer(self.vertex_ordinal),
                    }
                ),
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class IntersectionBoundaryVertexIdV1:
    incident_supports: tuple[IncidentSupport, ...]
    intersection_ordinal: int

    def __post_init__(self) -> None:
        if type(self.incident_supports) is not tuple or len(self.incident_supports) < 2:
            raise InvalidBoundaryVertexIdentityError("intersection boundary vertex requires at least two incident supports.")
        if any(type(incident) is not IncidentSupport for incident in self.incident_supports):
            raise InvalidBoundaryVertexIdentityError("intersection incident supports must be exact IncidentSupport.")
        canonical = tuple(sorted(self.incident_supports, key=lambda incident: incident.support_id.canonical_bytes))
        if canonical != self.incident_supports:
            raise InvalidBoundaryVertexIdentityError("intersection incident supports must use canonical pair order.")
        support_bytes = tuple(incident.support_id.canonical_bytes for incident in self.incident_supports)
        if len(support_bytes) != len(set(support_bytes)):
            raise InvalidBoundaryVertexIdentityError("intersection incident support identities must be distinct.")
        if type(self.intersection_ordinal) is not int or self.intersection_ordinal < 0:
            raise InvalidBoundaryVertexIdentityError("intersection ordinal must be an exact nonnegative integer.")

    @classmethod
    def build(
        cls,
        *,
        incident_supports: Sequence[IncidentSupport],
        intersection_ordinal: int,
    ) -> Self:
        try:
            supports = tuple(incident_supports)
        except TypeError:
            raise InvalidBoundaryVertexIdentityError("intersection incident supports must be finite.") from None
        if any(type(incident) is not IncidentSupport for incident in supports):
            raise InvalidBoundaryVertexIdentityError("intersection incident supports must be exact IncidentSupport.")
        supports = tuple(sorted(supports, key=lambda incident: incident.support_id.canonical_bytes))
        return cls(supports, intersection_ordinal)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not IntersectionBoundaryVertexIdV1:
            raise InvalidBoundaryVertexIdentityError("intersection identity must be exact IntersectionBoundaryVertexIdV1, not a subclass.")
        return encode_tagged_union(
            BOUNDARY_VERTEX_ID_VERSION,
            encode_tagged_union(
                b"support-intersection-v1",
                encode_component_map(
                    {
                        b"incident-supports": encode_sequence(tuple(incident.canonical_bytes for incident in self.incident_supports)),
                        b"intersection-ordinal": encode_integer(self.intersection_ordinal),
                    }
                ),
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class MultiIncidenceIntersectionBoundaryVertexIdV1:
    incident_supports: tuple[IncidentSupport, ...]
    intersection_ordinal: int

    def __post_init__(self) -> None:
        if type(self.incident_supports) is not tuple or len(self.incident_supports) < 3:
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection requires at least three incidents.")
        if any(type(incident) is not IncidentSupport for incident in self.incident_supports):
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection supports must be exact IncidentSupport.")
        canonical = tuple(
            sorted(
                self.incident_supports,
                key=lambda incident: (
                    incident.support_id.canonical_bytes,
                    incident.canonical_bytes,
                ),
            )
        )
        if canonical != self.incident_supports:
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection supports must use canonical order.")
        support_bytes = tuple(incident.support_id.canonical_bytes for incident in self.incident_supports)
        distinct_supports = set(support_bytes)
        if len(distinct_supports) < 2:
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection requires at least two distinct supports.")
        if len(distinct_supports) == len(support_bytes):
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection requires a repeated support.")
        if type(self.intersection_ordinal) is not int or self.intersection_ordinal < 0:
            raise InvalidBoundaryVertexIdentityError("intersection ordinal must be an exact nonnegative integer.")

    @classmethod
    def build(
        cls,
        *,
        incident_supports: Sequence[IncidentSupport],
        intersection_ordinal: int,
    ) -> Self:
        try:
            supports = tuple(incident_supports)
        except TypeError:
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection supports must be finite.") from None
        if any(type(incident) is not IncidentSupport for incident in supports):
            raise InvalidBoundaryVertexIdentityError("multi-incidence intersection supports must be exact IncidentSupport.")
        canonical = tuple(
            sorted(
                supports,
                key=lambda incident: (
                    incident.support_id.canonical_bytes,
                    incident.canonical_bytes,
                ),
            )
        )
        return cls(canonical, intersection_ordinal)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not MultiIncidenceIntersectionBoundaryVertexIdV1:
            raise InvalidBoundaryVertexIdentityError("multi-incidence identity must be exact MultiIncidenceIntersectionBoundaryVertexIdV1, not a subclass.")
        return encode_tagged_union(
            BOUNDARY_VERTEX_ID_VERSION,
            encode_tagged_union(
                b"multi-incidence-intersection-v1",
                encode_component_map(
                    {
                        b"incident-supports": encode_sequence(tuple(incident.canonical_bytes for incident in self.incident_supports)),
                        b"intersection-ordinal": encode_integer(self.intersection_ordinal),
                    }
                ),
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


@dataclass(frozen=True)
class ParameterSeamBoundaryVertexIdV1:
    incident_supports: tuple[IncidentSupport, IncidentSupport]
    seam_ordinal: int

    def __post_init__(self) -> None:
        if type(self.incident_supports) is not tuple or len(self.incident_supports) != 2:
            raise InvalidBoundaryVertexIdentityError("parameter seam requires exactly two incident supports.")
        if any(type(incident) is not IncidentSupport for incident in self.incident_supports):
            raise InvalidBoundaryVertexIdentityError("parameter seam incident supports must be exact IncidentSupport.")
        canonical = tuple(sorted(self.incident_supports, key=lambda incident: incident.canonical_bytes))
        if canonical != self.incident_supports:
            raise InvalidBoundaryVertexIdentityError("parameter seam incident supports must use canonical order.")
        support_bytes = {incident.support_id.canonical_bytes for incident in self.incident_supports}
        if len(support_bytes) != 1:
            raise InvalidBoundaryVertexIdentityError("parameter seam incidences must use the same support.")
        orientations = {incident.trim_incidence_orientation for incident in self.incident_supports}
        if orientations != {TrimIncidenceOrientation.ENTERING, TrimIncidenceOrientation.LEAVING}:
            raise InvalidBoundaryVertexIdentityError("parameter seam requires one entering and leaving incidence.")
        if type(self.seam_ordinal) is not int or self.seam_ordinal < 0:
            raise InvalidBoundaryVertexIdentityError("parameter seam ordinal must be an exact nonnegative integer.")

    @classmethod
    def build(
        cls,
        *,
        incident_supports: Sequence[IncidentSupport],
        seam_ordinal: int,
    ) -> Self:
        try:
            supports = tuple(incident_supports)
        except TypeError:
            raise InvalidBoundaryVertexIdentityError("parameter seam incident supports must be finite.") from None
        if len(supports) != 2 or any(type(incident) is not IncidentSupport for incident in supports):
            raise InvalidBoundaryVertexIdentityError("parameter seam requires exactly two exact IncidentSupport values.")
        canonical = tuple(sorted(supports, key=lambda incident: incident.canonical_bytes))
        return cls((canonical[0], canonical[1]), seam_ordinal)

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not ParameterSeamBoundaryVertexIdV1:
            raise InvalidBoundaryVertexIdentityError("parameter seam identity must be exact ParameterSeamBoundaryVertexIdV1, not a subclass.")
        return encode_tagged_union(
            BOUNDARY_VERTEX_ID_VERSION,
            encode_tagged_union(
                b"parameter-seam-v1",
                encode_component_map(
                    {
                        b"incident-supports": encode_sequence(tuple(incident.canonical_bytes for incident in self.incident_supports)),
                        b"seam-ordinal": encode_integer(self.seam_ordinal),
                    }
                ),
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


BoundaryVertexIdV1: TypeAlias = InputRingVertexIdV1 | IntersectionBoundaryVertexIdV1 | MultiIncidenceIntersectionBoundaryVertexIdV1 | ParameterSeamBoundaryVertexIdV1
