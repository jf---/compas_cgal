import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from math import gcd
from typing import Final
from typing import NewType
from typing import Self
from typing import TypeAlias

from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import canonical_record_kind
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import ArtifactIdentityError
from compas_cgal.adaptive.errors import CanonicalEncodingError

INPUT_SCHEMA_VERSION: Final[bytes] = b"adaptive-input-schema-v1"
OPERATION_SCHEMA_VERSION: Final[bytes] = b"adaptive-operation-schema-v1"
COMPONENT_IDENTITY_VERSION: Final[bytes] = b"component-identity-v1"
BOUNDARY_VERTEX_ID_VERSION: Final[bytes] = b"boundary-vertex-id-v1"

IdentityDigest = NewType("IdentityDigest", bytes)
ComponentDomainTag = NewType("ComponentDomainTag", bytes)
StrategyVersion = NewType("StrategyVersion", bytes)
SourceRevision = NewType("SourceRevision", bytes)
NativeSourceTreeDigest = NewType("NativeSourceTreeDigest", bytes)


def _nonempty_bytes(value: object, name: str) -> bytes:
    if type(value) is not bytes or not value:
        raise ArtifactIdentityError(f"{name} must be nonempty bytes.")
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
            raise ArtifactIdentityError("native source-tree digest must be exactly one 32-byte SHA-256 digest.")
        try:
            require_canonical_record(self.canonical_parameter_bytes)
        except CanonicalEncodingError:
            raise ArtifactIdentityError("canonical parameter bytes must contain one complete versioned CCAN record.") from None
        if canonical_record_kind(self.canonical_parameter_bytes) != b"T":
            raise ArtifactIdentityError("canonical parameter bytes must be a tagged domain record.")

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
        raise ArtifactIdentityError("incident support coefficients must include a nonzero exact value.")
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
        if not isinstance(self.support_kind, SupportKind):
            raise ArtifactIdentityError("incident support kind must be line or circle.")
        if not isinstance(self.primitive_coefficients, tuple) or any(type(value) is not int for value in self.primitive_coefficients):
            raise ArtifactIdentityError("incident support coefficients must be exact primitive integers.")
        expected_count = 3 if self.support_kind is SupportKind.LINE else 4
        if len(self.primitive_coefficients) != expected_count:
            expected_name = "three" if expected_count == 3 else "four"
            raise ArtifactIdentityError(f"{self.support_kind.name.lower()} support requires exactly {expected_name} coefficients.")
        normalized = _primitive_coefficients(tuple(ExactRationalV1(value, 1) for value in self.primitive_coefficients))
        if normalized != self.primitive_coefficients:
            raise ArtifactIdentityError("incident support coefficients must already be primitive with positive leading sign.")
        if self.support_kind is SupportKind.LINE and self.primitive_coefficients[:2] == (0, 0):
            raise ArtifactIdentityError("line support requires a nonzero normal.")
        if self.support_kind is SupportKind.CIRCLE and self.primitive_coefficients[0] == 0:
            raise ArtifactIdentityError("circle support requires a nonzero quadratic coefficient.")

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
            raise ArtifactIdentityError("incident support coefficients must be finite.") from None
        if any(not isinstance(coefficient, ExactRationalV1) for coefficient in coefficients):
            raise ArtifactIdentityError("incident support coefficients must be exact rationals.")
        expected_count = 3 if support_kind is SupportKind.LINE else 4 if support_kind is SupportKind.CIRCLE else 0
        if len(coefficients) != expected_count:
            expected_name = "three" if expected_count == 3 else "four"
            raise ArtifactIdentityError(f"incident support requires exactly {expected_name} exact coefficients.")
        return cls(support_kind, _primitive_coefficients(coefficients))

    @property
    def canonical_bytes(self) -> bytes:
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
        if not isinstance(self.support_id, IncidentSupportIdV1):
            raise ArtifactIdentityError("incident support requires a normalized support identity.")
        if not isinstance(self.trim_incidence_orientation, TrimIncidenceOrientation):
            raise ArtifactIdentityError("incident support requires an explicit trim-incidence orientation.")

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
    canonical_ring_bytes: bytes
    vertex_ordinal: int

    def __post_init__(self) -> None:
        try:
            require_canonical_record(self.canonical_ring_bytes)
        except CanonicalEncodingError:
            raise ArtifactIdentityError("input-ring vertex requires one complete versioned canonical ring record.") from None
        if canonical_record_kind(self.canonical_ring_bytes) != b"T":
            raise ArtifactIdentityError("input-ring vertex requires a tagged canonical ring record.")
        if type(self.vertex_ordinal) is not int or self.vertex_ordinal < 0:
            raise ArtifactIdentityError("input-ring vertex ordinal must be an exact nonnegative integer.")

    @classmethod
    def build(cls, *, canonical_ring_bytes: bytes, vertex_ordinal: int) -> Self:
        return cls(canonical_ring_bytes, vertex_ordinal)

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            BOUNDARY_VERTEX_ID_VERSION,
            encode_tagged_union(
                b"input-ring-vertex-v1",
                encode_component_map(
                    {
                        b"ring": self.canonical_ring_bytes,
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
        if not isinstance(self.incident_supports, tuple) or len(self.incident_supports) < 2:
            raise ArtifactIdentityError("intersection boundary vertex requires at least two incident supports.")
        if any(not isinstance(incident, IncidentSupport) for incident in self.incident_supports):
            raise ArtifactIdentityError("intersection incident supports must be typed.")
        canonical = tuple(sorted(self.incident_supports, key=lambda incident: incident.canonical_bytes))
        if canonical != self.incident_supports:
            raise ArtifactIdentityError("intersection incident supports must use canonical pair order.")
        support_bytes = tuple(incident.support_id.canonical_bytes for incident in self.incident_supports)
        if len(support_bytes) != len(set(support_bytes)):
            raise ArtifactIdentityError("intersection incident support identities must be distinct.")
        if type(self.intersection_ordinal) is not int or self.intersection_ordinal < 0:
            raise ArtifactIdentityError("intersection ordinal must be an exact nonnegative integer.")

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
            raise ArtifactIdentityError("intersection incident supports must be finite.") from None
        if any(not isinstance(incident, IncidentSupport) for incident in supports):
            raise ArtifactIdentityError("intersection incident supports must be typed.")
        supports = tuple(sorted(supports, key=lambda incident: incident.canonical_bytes))
        return cls(supports, intersection_ordinal)

    @property
    def canonical_bytes(self) -> bytes:
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


BoundaryVertexIdV1: TypeAlias = InputRingVertexIdV1 | IntersectionBoundaryVertexIdV1
