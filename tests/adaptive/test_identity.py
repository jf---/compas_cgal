import hashlib
from dataclasses import fields
from fractions import Fraction

import pytest

from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import ArtifactIdentityError
from compas_cgal.adaptive.identity import BOUNDARY_VERTEX_ID_VERSION
from compas_cgal.adaptive.identity import COMPONENT_IDENTITY_VERSION
from compas_cgal.adaptive.identity import INPUT_SCHEMA_VERSION
from compas_cgal.adaptive.identity import OPERATION_SCHEMA_VERSION
from compas_cgal.adaptive.identity import BoundaryVertexIdV1
from compas_cgal.adaptive.identity import ComponentIdentity
from compas_cgal.adaptive.identity import ComponentDomainTag
from compas_cgal.adaptive.identity import IncidentSupport
from compas_cgal.adaptive.identity import IncidentSupportIdV1
from compas_cgal.adaptive.identity import InputRingVertexIdV1
from compas_cgal.adaptive.identity import IntersectionBoundaryVertexIdV1
from compas_cgal.adaptive.identity import NativeSourceTreeDigest
from compas_cgal.adaptive.identity import SourceRevision
from compas_cgal.adaptive.identity import StrategyVersion
from compas_cgal.adaptive.identity import SupportKind
from compas_cgal.adaptive.identity import TrimIncidenceOrientation


def _component_identity(
    *,
    component_domain: bytes = b"motion-certifier",
    strategy_version: bytes = b"event-exact-v1",
    source_revision: bytes = b"0123456789abcdef",
    native_digest: bytes = b"\x10" * 32,
    parameters: bytes = encode_tagged_union(b"motion-parameters-v1", b"parameters"),
) -> ComponentIdentity:
    return ComponentIdentity.build(
        component_domain=ComponentDomainTag(component_domain),
        strategy_version=StrategyVersion(strategy_version),
        source_revision=SourceRevision(source_revision),
        native_source_tree_digest=NativeSourceTreeDigest(native_digest),
        canonical_parameter_bytes=parameters,
    )


def _line_support(a: int, b: int, c: int) -> IncidentSupportIdV1:
    return IncidentSupportIdV1.build(
        support_kind=SupportKind.LINE,
        normalized_coefficients=(
            ExactRationalV1.build(a, 1),
            ExactRationalV1.build(b, 1),
            ExactRationalV1.build(c, 1),
        ),
    )


def test_schema_and_component_versions_are_frozen_nonempty_bytes() -> None:
    versions = (
        INPUT_SCHEMA_VERSION,
        OPERATION_SCHEMA_VERSION,
        COMPONENT_IDENTITY_VERSION,
        BOUNDARY_VERTEX_ID_VERSION,
    )

    assert all(type(version) is bytes and version for version in versions)
    assert len(set(versions)) == len(versions)


def test_component_identity_binds_every_declared_input() -> None:
    baseline = _component_identity()

    variants = (
        _component_identity(component_domain=b"depletion"),
        _component_identity(strategy_version=b"event-exact-v0"),
        _component_identity(source_revision=b"0123456789abcdee"),
        _component_identity(native_digest=b"\x10" * 31 + b"\x11"),
        _component_identity(parameters=encode_tagged_union(b"motion-parameters-v1", b"changed")),
    )

    assert all(variant.digest != baseline.digest for variant in variants)


def test_component_identity_hashes_only_its_canonical_inputs() -> None:
    identity = _component_identity()

    assert "digest" not in {field.name for field in fields(identity)}
    assert identity.digest == hashlib.sha256(identity.canonical_bytes).digest()
    assert identity.digest not in identity.canonical_bytes


def test_component_identity_rejects_untyped_or_malformed_inputs() -> None:
    with pytest.raises(ArtifactIdentityError, match="strategy version"):
        ComponentIdentity.build(
            strategy_version=StrategyVersion(b""),
            component_domain=ComponentDomainTag(b"motion-certifier"),
            source_revision=SourceRevision(b"revision"),
            native_source_tree_digest=NativeSourceTreeDigest(b"\x00" * 32),
            canonical_parameter_bytes=b"parameters",
        )
    with pytest.raises(ArtifactIdentityError, match="32-byte"):
        _component_identity(native_digest=b"short")
    with pytest.raises(ArtifactIdentityError, match="versioned"):
        _component_identity(parameters=b"parameters")
    with pytest.raises(ArtifactIdentityError, match="tagged"):
        _component_identity(parameters=encode_bytes(b"parameters"))
    with pytest.raises(ArtifactIdentityError, match="one complete"):
        _component_identity(parameters=encode_tagged_union(b"motion-parameters-v1", b"parameters") + b"trailing")


def test_incident_support_identity_normalizes_rational_scale_and_sign() -> None:
    primitive = _line_support(2, 4, 6)
    scaled = IncidentSupportIdV1.build(
        support_kind=SupportKind.LINE,
        normalized_coefficients=(
            ExactRationalV1.build(-1, 2),
            ExactRationalV1.build(-1, 1),
            ExactRationalV1.build(-3, 2),
        ),
    )

    assert primitive.canonical_bytes == scaled.canonical_bytes
    assert primitive.primitive_coefficients == (1, 2, 3)


def test_input_ring_boundary_vertex_binds_ring_and_exact_ordinal() -> None:
    ring_bytes = encode_tagged_union(b"input-ring-v1", b"outer-ring")
    vertex = InputRingVertexIdV1.build(canonical_ring_bytes=ring_bytes, vertex_ordinal=3)
    changed = InputRingVertexIdV1.build(canonical_ring_bytes=ring_bytes, vertex_ordinal=4)

    boundary_vertex: BoundaryVertexIdV1 = vertex

    assert boundary_vertex.canonical_bytes != changed.canonical_bytes


def test_input_ring_boundary_vertex_rejects_incomplete_canonical_record() -> None:
    incomplete = encode_tagged_union(b"input-ring-v1", b"outer-ring")[:-1]

    with pytest.raises(ArtifactIdentityError, match="one complete"):
        InputRingVertexIdV1.build(canonical_ring_bytes=incomplete, vertex_ordinal=0)
    with pytest.raises(ArtifactIdentityError, match="tagged"):
        InputRingVertexIdV1.build(
            canonical_ring_bytes=encode_bytes(b"outer-ring"),
            vertex_ordinal=0,
        )


def test_intersection_boundary_vertex_sorts_incident_support_ids() -> None:
    first = _line_support(1, 0, -1)
    second = _line_support(0, 1, -2)
    first_incidence = IncidentSupport.build(
        support_id=first,
        trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
    )
    second_incidence = IncidentSupport.build(
        support_id=second,
        trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
    )
    forward = IntersectionBoundaryVertexIdV1.build(
        incident_supports=(first_incidence, second_incidence),
        intersection_ordinal=0,
    )
    reversed_supports = IntersectionBoundaryVertexIdV1.build(
        incident_supports=(second_incidence, first_incidence),
        intersection_ordinal=0,
    )

    assert forward == reversed_supports
    assert forward.canonical_bytes == reversed_supports.canonical_bytes


def test_intersection_boundary_vertex_binds_ordinal_and_trim_orientation() -> None:
    first = _line_support(1, 0, -1)
    second = _line_support(0, 1, -2)
    supports = (
        IncidentSupport.build(
            support_id=first,
            trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
        ),
        IncidentSupport.build(
            support_id=second,
            trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
        ),
    )
    baseline = IntersectionBoundaryVertexIdV1.build(
        incident_supports=supports,
        intersection_ordinal=0,
    )
    changed_ordinal = IntersectionBoundaryVertexIdV1.build(
        incident_supports=supports,
        intersection_ordinal=1,
    )
    changed_support_orientation = IncidentSupport.build(
        support_id=first,
        trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
    )
    changed_orientation = IntersectionBoundaryVertexIdV1.build(
        incident_supports=(changed_support_orientation, supports[1]),
        intersection_ordinal=0,
    )

    assert len({baseline.digest, changed_ordinal.digest, changed_orientation.digest}) == 3
    assert "root" not in {field.name for field in fields(baseline)}
    assert "coordinate" not in {field.name for field in fields(baseline)}


def test_intersection_boundary_vertex_rejects_duplicate_supports() -> None:
    support = _line_support(1, 0, -1)
    incidence = IncidentSupport.build(
        support_id=support,
        trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
    )

    with pytest.raises(ArtifactIdentityError, match="distinct"):
        IntersectionBoundaryVertexIdV1.build(
            incident_supports=(incidence, incidence),
            intersection_ordinal=0,
        )


def test_intersection_boundary_vertex_rejects_malformed_support_before_sorting() -> None:
    with pytest.raises(ArtifactIdentityError, match="typed"):
        IntersectionBoundaryVertexIdV1.build(
            incident_supports=(object(), object()),  # type: ignore[arg-type]
            intersection_ordinal=0,
        )


def test_incident_support_rejects_degenerate_exact_coefficients() -> None:
    zero = ExactRationalV1.build(0, 1)

    with pytest.raises(ArtifactIdentityError, match="nonzero"):
        IncidentSupportIdV1.build(
            support_kind=SupportKind.LINE,
            normalized_coefficients=(zero, zero, zero),
        )
    with pytest.raises(ArtifactIdentityError, match="four"):
        IncidentSupportIdV1.build(
            support_kind=SupportKind.CIRCLE,
            normalized_coefficients=(ExactRationalV1.build(Fraction(1, 2).numerator, 2), zero),
        )
