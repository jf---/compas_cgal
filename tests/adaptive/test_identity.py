import hashlib
import math
from dataclasses import dataclass
from dataclasses import fields
from dataclasses import replace
from fractions import Fraction

import pytest

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import InvalidBoundaryVertexIdentityError
from compas_cgal.adaptive.errors import InvalidComponentIdentityError
from compas_cgal.adaptive.errors import InvalidInputIdentityError
from compas_cgal.adaptive.identity import BOUNDARY_VERTEX_ID_VERSION
from compas_cgal.adaptive.identity import COMPONENT_IDENTITY_VERSION
from compas_cgal.adaptive.identity import INPUT_SCHEMA_VERSION
from compas_cgal.adaptive.identity import OPERATION_SCHEMA_VERSION
from compas_cgal.adaptive.identity import BoundaryVertexIdV1
from compas_cgal.adaptive.identity import ComponentIdentity
from compas_cgal.adaptive.identity import ComponentDomainTag
from compas_cgal.adaptive.identity import FrameIdentity
from compas_cgal.adaptive.identity import IncidentSupport
from compas_cgal.adaptive.identity import IncidentSupportIdV1
from compas_cgal.adaptive.identity import InputIdentity
from compas_cgal.adaptive.identity import InputRingVertexIdV1
from compas_cgal.adaptive.identity import IntersectionBoundaryVertexIdV1
from compas_cgal.adaptive.identity import MultiIncidenceIntersectionBoundaryVertexIdV1
from compas_cgal.adaptive.identity import NativeSourceTreeDigest
from compas_cgal.adaptive.identity import ParameterSeamBoundaryVertexIdV1
from compas_cgal.adaptive.identity import SourceRevision
from compas_cgal.adaptive.identity import StrategyVersion
from compas_cgal.adaptive.identity import SupportKind
from compas_cgal.adaptive.identity import TrimIncidenceOrientation
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


@dataclass(frozen=True)
class SemanticComponentIdentity(ComponentIdentity):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticIncidentSupport(IncidentSupport):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticInputRingVertex(InputRingVertexIdV1):
    provenance_bit: bool


def _outer_ring() -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(
        (
            Point2[WorldXY].build(0.0, 0.0),
            Point2[WorldXY].build(4.0, 0.0),
            Point2[WorldXY].build(4.0, 4.0),
            Point2[WorldXY].build(0.0, 4.0),
        )
    )


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
    with pytest.raises(InvalidComponentIdentityError, match="strategy version"):
        ComponentIdentity.build(
            strategy_version=StrategyVersion(b""),
            component_domain=ComponentDomainTag(b"motion-certifier"),
            source_revision=SourceRevision(b"revision"),
            native_source_tree_digest=NativeSourceTreeDigest(b"\x00" * 32),
            canonical_parameter_bytes=b"parameters",
        )
    with pytest.raises(InvalidComponentIdentityError, match="32-byte"):
        _component_identity(native_digest=b"short")
    with pytest.raises(InvalidComponentIdentityError, match="versioned"):
        _component_identity(parameters=b"parameters")
    with pytest.raises(InvalidComponentIdentityError, match="tagged"):
        _component_identity(parameters=encode_bytes(b"parameters"))
    with pytest.raises(InvalidComponentIdentityError, match="one complete"):
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
    ring = _outer_ring()
    vertex = InputRingVertexIdV1.build(canonical_ring=ring, vertex_ordinal=2)
    changed = InputRingVertexIdV1.build(canonical_ring=ring, vertex_ordinal=3)

    boundary_vertex: BoundaryVertexIdV1 = vertex

    assert boundary_vertex.canonical_bytes != changed.canonical_bytes


def test_input_ring_boundary_vertex_requires_typed_ring_and_bounded_ordinal() -> None:
    ring = _outer_ring()

    with pytest.raises(InvalidBoundaryVertexIdentityError, match="CanonicalRingV1"):
        InputRingVertexIdV1.build(
            canonical_ring=encode_tagged_union(b"not-a-ring-v1", b"payload"),  # type: ignore[arg-type]
            vertex_ordinal=0,
        )
    with pytest.raises(InvalidBoundaryVertexIdentityError, match="vertex count"):
        InputRingVertexIdV1.build(canonical_ring=ring, vertex_ordinal=ring.vertex_count)


def _oppositely_ordered_incident_supports() -> tuple[IncidentSupport, IncidentSupport]:
    support_low, support_high = sorted(
        (_line_support(1, 0, -1), _line_support(0, 1, -2)),
        key=lambda support: support.canonical_bytes,
    )
    incidence_low = IncidentSupport.build(
        support_id=support_low,
        trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
    )
    incidence_high = IncidentSupport.build(
        support_id=support_high,
        trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
    )

    assert incidence_high.canonical_bytes < incidence_low.canonical_bytes
    return incidence_low, incidence_high


def test_intersection_boundary_vertex_build_sorts_only_by_support_id() -> None:
    incidence_low, incidence_high = _oppositely_ordered_incident_supports()
    canonical_order = (incidence_low, incidence_high)
    forward = IntersectionBoundaryVertexIdV1.build(
        incident_supports=canonical_order,
        intersection_ordinal=0,
    )
    reversed_supports = IntersectionBoundaryVertexIdV1.build(
        incident_supports=tuple(reversed(canonical_order)),
        intersection_ordinal=0,
    )

    assert forward.incident_supports == canonical_order
    assert forward == reversed_supports
    assert forward.canonical_bytes == reversed_supports.canonical_bytes


def test_intersection_boundary_vertex_raw_order_is_by_support_id() -> None:
    incidence_low, incidence_high = _oppositely_ordered_incident_supports()
    canonical_order = (incidence_low, incidence_high)

    raw = IntersectionBoundaryVertexIdV1(canonical_order, 0)

    assert raw.incident_supports == canonical_order
    with pytest.raises(InvalidBoundaryVertexIdentityError, match="canonical pair order"):
        IntersectionBoundaryVertexIdV1(tuple(reversed(canonical_order)), 0)


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

    with pytest.raises(InvalidBoundaryVertexIdentityError, match="distinct"):
        IntersectionBoundaryVertexIdV1.build(
            incident_supports=(incidence, incidence),
            intersection_ordinal=0,
        )


def test_parameter_seam_vertex_preserves_both_oriented_incidences() -> None:
    support = _line_support(0, 1, 0)
    entering = IncidentSupport.build(
        support_id=support,
        trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
    )
    leaving = IncidentSupport.build(
        support_id=support,
        trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
    )

    forward = ParameterSeamBoundaryVertexIdV1.build(
        incident_supports=(entering, leaving),
        seam_ordinal=0,
    )
    reverse = ParameterSeamBoundaryVertexIdV1.build(
        incident_supports=(leaving, entering),
        seam_ordinal=0,
    )
    changed_ordinal = ParameterSeamBoundaryVertexIdV1.build(
        incident_supports=(entering, leaving),
        seam_ordinal=1,
    )
    boundary_vertex: BoundaryVertexIdV1 = forward

    assert forward == reverse
    assert boundary_vertex.canonical_bytes != changed_ordinal.canonical_bytes


def test_parameter_seam_vertex_rejects_non_seam_incidences() -> None:
    first = IncidentSupport.build(
        support_id=_line_support(0, 1, 0),
        trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
    )
    second_support = IncidentSupport.build(
        support_id=_line_support(1, 0, 0),
        trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
    )
    same_orientation = IncidentSupport.build(
        support_id=first.support_id,
        trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
    )

    with pytest.raises(InvalidBoundaryVertexIdentityError, match="same support"):
        ParameterSeamBoundaryVertexIdV1.build(
            incident_supports=(first, second_support),
            seam_ordinal=0,
        )
    with pytest.raises(InvalidBoundaryVertexIdentityError, match="entering and leaving"):
        ParameterSeamBoundaryVertexIdV1.build(
            incident_supports=(first, same_orientation),
            seam_ordinal=0,
        )


def test_multi_incidence_intersection_preserves_every_orientation() -> None:
    first_support = _line_support(0, 1, 0)
    second_support = _line_support(1, 0, 0)
    incidents = tuple(
        IncidentSupport.build(
            support_id=support,
            trim_incidence_orientation=orientation,
        )
        for support in (first_support, second_support)
        for orientation in (
            TrimIncidenceOrientation.ENTERING,
            TrimIncidenceOrientation.LEAVING,
        )
    )

    forward = MultiIncidenceIntersectionBoundaryVertexIdV1.build(
        incident_supports=incidents,
        intersection_ordinal=0,
    )
    reverse = MultiIncidenceIntersectionBoundaryVertexIdV1.build(
        incident_supports=reversed(incidents),
        intersection_ordinal=0,
    )
    boundary_vertex: BoundaryVertexIdV1 = forward

    assert forward == reverse
    assert len(boundary_vertex.incident_supports) == 4


def test_multi_incidence_intersection_requires_repeated_support() -> None:
    incidents = (
        IncidentSupport.build(
            support_id=_line_support(0, 1, 0),
            trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
        ),
        IncidentSupport.build(
            support_id=_line_support(1, 0, 0),
            trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
        ),
        IncidentSupport.build(
            support_id=_line_support(1, 1, 0),
            trim_incidence_orientation=TrimIncidenceOrientation.ENTERING,
        ),
    )

    with pytest.raises(InvalidBoundaryVertexIdentityError, match="repeated support"):
        MultiIncidenceIntersectionBoundaryVertexIdV1.build(
            incident_supports=incidents,
            intersection_ordinal=0,
        )


def test_intersection_boundary_vertex_rejects_malformed_support_before_sorting() -> None:
    with pytest.raises(InvalidBoundaryVertexIdentityError, match="exact"):
        IntersectionBoundaryVertexIdV1.build(
            incident_supports=(object(), object()),  # type: ignore[arg-type]
            intersection_ordinal=0,
        )


def test_incident_support_rejects_degenerate_exact_coefficients() -> None:
    zero = ExactRationalV1.build(0, 1)

    with pytest.raises(InvalidBoundaryVertexIdentityError, match="nonzero"):
        IncidentSupportIdV1.build(
            support_kind=SupportKind.LINE,
            normalized_coefficients=(zero, zero, zero),
        )
    with pytest.raises(InvalidBoundaryVertexIdentityError, match="four"):
        IncidentSupportIdV1.build(
            support_kind=SupportKind.CIRCLE,
            normalized_coefficients=(ExactRationalV1.build(Fraction(1, 2).numerator, 2), zero),
        )


def test_component_identity_subclass_cannot_drop_added_semantics() -> None:
    baseline = _component_identity()
    variants = (
        SemanticComponentIdentity(
            baseline.component_domain,
            baseline.strategy_version,
            baseline.source_revision,
            baseline.native_source_tree_digest,
            baseline.canonical_parameter_bytes,
            False,
        ),
        SemanticComponentIdentity(
            baseline.component_domain,
            baseline.strategy_version,
            baseline.source_revision,
            baseline.native_source_tree_digest,
            baseline.canonical_parameter_bytes,
            True,
        ),
    )

    assert variants[0] != variants[1]
    for variant in variants:
        with pytest.raises(InvalidComponentIdentityError, match="exact ComponentIdentity"):
            _ = variant.canonical_bytes


def test_boundary_identity_subclasses_fail_before_canonical_bytes() -> None:
    support = _line_support(1, 0, -1)
    semantic_incident = SemanticIncidentSupport(
        support,
        TrimIncidenceOrientation.ENTERING,
        True,
    )
    ring = _outer_ring()
    semantic_vertex = SemanticInputRingVertex(ring, 0, True)

    with pytest.raises(InvalidBoundaryVertexIdentityError, match="exact IncidentSupport"):
        IntersectionBoundaryVertexIdV1.build(
            incident_supports=(
                semantic_incident,
                IncidentSupport.build(
                    support_id=_line_support(0, 1, -1),
                    trim_incidence_orientation=TrimIncidenceOrientation.LEAVING,
                ),
            ),
            intersection_ordinal=0,
        )
    with pytest.raises(InvalidBoundaryVertexIdentityError, match="exact InputRingVertexIdV1"):
        _ = semantic_vertex.canonical_bytes


def _input_identity(
    *,
    bore_evidence: bytes = b"qualified-bore-metrology-v1",
    spatial_resolution: float = 0.5,
    cut_intent: CutIntent = CutIntent.CLIMB,
) -> InputIdentity:
    design_boundary = _outer_ring()
    tool_radius = ToolRadius.build(0.5)
    reachable_domain = ReachableDomain.build(
        design_boundary=design_boundary,
        holes=(),
        tool_radius=tool_radius,
    )
    cut_plane = CutPlane.build(
        CutZ.build(-2.0),
        ClearanceZ.build(5.0),
    )
    entry = PreclearedEntry.build(
        reachable_domain=reachable_domain,
        center=Point2[WorldXY].build(2.0, 2.0),
        radius=EntryRadius.build(1.0),
        tool_radius=tool_radius,
        cut_plane=cut_plane,
        qualified_bore=QualifiedBore.build(
            cut_plane=cut_plane,
            process_identity=BoreProcessIdentity(b"predrill-cycle-revision-7"),
            evidence_bytes=bore_evidence,
        ),
    )
    user_cap = EngagementCap.build(math.radians(120.0))
    neck_policy = NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
        effective_caps={
            (neck_class, passage_state): EngagementCap.build(
                math.radians(90.0 - 10.0 * passage_state.rank),
            )
            for neck_class in range(2)
            for passage_state in ACTIVE_PASSAGE_STATES
        },
    )
    return InputIdentity.build(
        design_boundary=design_boundary,
        holes=(),
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        reachable_domain=reachable_domain,
        entry=entry,
        user_cap=user_cap,
        candidate_policy=CandidatePolicy.build(
            spatial_resolution=Spacing.build(spatial_resolution),
            spatial_refinement_levels=2,
            radius_resolution=Spacing.build(0.125),
            radius_refinement_levels=2,
            phase_count=4,
            minimum_guide_radius=GuideRadius.build(0.125),
            minimum_progress=Spacing.build(0.25),
        ),
        neck_policy=neck_policy,
        depletion_policy=DepletionPolicy.build(
            chord_bound=ChordBound.build(0.125),
            center_count_limit=4096,
        ),
        traversal_policy=TraversalPolicy.build(forward_window=4),
        cut_direction_policy=CutDirectionPolicy.build(cut_intent),
    )


def test_input_identity_binds_entry_domain_policies_schemas_and_components() -> None:
    identity = _input_identity()
    canonical = identity.canonical_bytes

    assert identity.reachable_domain_digest == identity.reachable_domain.certificate.digest
    assert identity.entry.reachable_domain is identity.reachable_domain
    assert INPUT_SCHEMA_VERSION in canonical
    assert OPERATION_SCHEMA_VERSION in canonical
    assert identity.entry.canonical_bytes in canonical
    assert identity.reachable_domain_digest in canonical
    assert identity.component_versions
    assert all(binding.version in canonical for binding in identity.component_versions)
    assert "digest" not in {field.name for field in fields(identity)}
    assert identity.digest == hashlib.sha256(canonical).digest()
    assert identity.digest not in canonical

    variants = (
        _input_identity(bore_evidence=b"qualified-bore-metrology-v2"),
        _input_identity(spatial_resolution=0.25),
        _input_identity(cut_intent=CutIntent.CONVENTIONAL),
    )
    assert all(variant.digest != identity.digest for variant in variants)


def test_input_identity_accepts_only_its_validated_precleared_entry() -> None:
    identity = _input_identity()

    with pytest.raises(InvalidInputIdentityError, match="PreclearedEntry"):
        replace(identity, entry=object())

    equivalent_but_distinct = _input_identity()
    with pytest.raises(InvalidInputIdentityError, match="PreclearedEntry"):
        replace(
            identity,
            reachable_domain=equivalent_but_distinct.reachable_domain,
            reachable_domain_digest=(equivalent_but_distinct.reachable_domain_digest),
        )


def test_input_identity_rejects_frame_and_component_version_relabelling() -> None:
    identity = _input_identity()
    changed_binding = replace(
        identity.component_versions[0],
        version=b"unowned-component-version-v0",
    )

    with pytest.raises(InvalidInputIdentityError, match="world-XY"):
        replace(
            identity,
            frame_identity=FrameIdentity(b"machine-xy-millimetre-frame-v1"),
        )
    with pytest.raises(InvalidInputIdentityError, match="component versions"):
        replace(
            identity,
            component_versions=(
                changed_binding,
                *identity.component_versions[1:],
            ),
        )
