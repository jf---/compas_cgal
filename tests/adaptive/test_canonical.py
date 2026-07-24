import math
import struct
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations

import pytest

from compas_cgal.adaptive.canonical import CANONICAL_ENCODING_VERSION
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import canonical_cut_z_bytes
from compas_cgal.adaptive.canonical import canonical_point2_bytes
from compas_cgal.adaptive.canonical import canonical_polygon_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import canonical_vector2_bytes
from compas_cgal.adaptive.canonical import encode_boolean
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import CanonicalEncodingError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


@dataclass(frozen=True)
class SemanticCutZ(CutZ):
    rapid: bool


@dataclass(frozen=True)
class SemanticPoint(Point2[WorldXY]):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticVector(Vector2[WorldXY]):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticSegment(ExactSegmentMotion):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticCircle(ExactCircleMotion):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticCandidatePolicy(CandidatePolicy):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticSpacing(Spacing):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticTraversalPolicy(TraversalPolicy):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticRational(ExactRationalV1):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticRing(CanonicalRingV1):
    provenance_bit: bool


def _candidate_policy() -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.5),
        spatial_refinement_levels=3,
        radius_resolution=Spacing.build(0.25),
        radius_refinement_levels=2,
        phase_count=8,
        minimum_guide_radius=GuideRadius.build(1.0),
        minimum_progress=Spacing.build(0.125),
    )


def _neck_policy() -> NeckPolicy:
    user_cap = EngagementCap.build(math.pi / 2.0)
    effective_cap = EngagementCap.build(math.pi / 3.0)
    return NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(Fraction(25, 1),),
        effective_caps={(neck_class, passage_state): effective_cap for neck_class in range(2) for passage_state in ACTIVE_PASSAGE_STATES},
    )


def test_primitives_are_versioned_big_endian_and_domain_separated() -> None:
    encoded_integer = encode_integer(0x0102)

    assert encoded_integer.startswith(CANONICAL_ENCODING_VERSION)
    assert encoded_integer.endswith(b"\x01\x02")
    assert encoded_integer != encode_binary64(258.0)
    assert encode_integer(-258) != encoded_integer
    with pytest.raises(CanonicalEncodingError, match="not bool"):
        encode_integer(True)


@pytest.mark.parametrize(
    ("encoded", "expected_hex"),
    (
        (encode_integer(0), "4343414e000149000000000000000100"),
        (encode_integer(-258), "4343414e0001490000000000000003010102"),
        (encode_binary64(1.5), "4343414e00014400000000000000083ff8000000000000"),
        (encode_boolean(True), "4343414e00013f000000000000000101"),
        (
            ExactRationalV1.build(1, 2).canonical_bytes,
            "4343414e00015200000000000000224343414e000149000000000000000200014343414e00014900000000000000020002",
        ),
        (
            encode_sequence((b"a", b"bc")),
            "4343414e000153000000000000002900000000000000024343414e0001420000000000000001614343414e00014200000000000000026263",
        ),
        (
            encode_tagged_union(b"kind-v1", b"payload"),
            "4343414e000154000000000000002c4343414e00014200000000000000076b696e642d76314343414e00014200000000000000077061796c6f6164",
        ),
        (
            encode_component_map({b"b": b"2", b"a": b"1"}),
            "4343414e00014d0000000000000048"
            "0000000000000002"
            "4343414e000142000000000000000161"
            "4343414e000142000000000000000131"
            "4343414e000142000000000000000162"
            "4343414e000142000000000000000132",
        ),
    ),
)
def test_frozen_grammar_golden_vectors(encoded: bytes, expected_hex: str) -> None:
    assert encoded == bytes.fromhex(expected_hex)


def test_binary64_normalizes_signed_zero_and_rejects_nonfinite_values() -> None:
    assert encode_binary64(-0.0) == encode_binary64(0.0)
    assert encode_binary64(1.0).endswith(struct.pack(">d", 1.0))

    for nonfinite in (math.nan, math.inf, -math.inf):
        with pytest.raises(CanonicalEncodingError, match="finite"):
            encode_binary64(nonfinite)


def test_exact_rational_normalizes_sign_and_common_factors() -> None:
    reduced = ExactRationalV1.build(1, 2)

    assert ExactRationalV1.build(2, 4) == reduced
    assert ExactRationalV1.build(-2, -4) == reduced
    assert reduced.numerator == 1
    assert reduced.denominator == 2
    assert reduced.canonical_bytes.startswith(CANONICAL_ENCODING_VERSION)


def test_sequences_and_tagged_unions_are_collision_unambiguous() -> None:
    assert encode_sequence((b"a", b"bc")) != encode_sequence((b"ab", b"c"))
    assert encode_tagged_union(b"a", b"bc") != encode_tagged_union(b"ab", b"c")
    assert encode_sequence(()) != encode_sequence((b"",))


def test_canonical_record_validator_rejects_unknown_or_malformed_nodes() -> None:
    unknown_kind = CANONICAL_ENCODING_VERSION + b"X" + struct.pack(">Q", 0)
    invalid_integer_sign = CANONICAL_ENCODING_VERSION + b"I" + struct.pack(">Q", 1) + b"\x02"
    truncated_binary64 = CANONICAL_ENCODING_VERSION + b"D" + struct.pack(">Q", 7) + b"\x00" * 7

    for malformed in (unknown_kind, invalid_integer_sign, truncated_binary64):
        with pytest.raises(CanonicalEncodingError):
            require_canonical_record(malformed)


def test_canonical_record_validator_rejects_noncanonical_structured_payloads() -> None:
    def node(kind: bytes, payload: bytes) -> bytes:
        return CANONICAL_ENCODING_VERSION + kind + struct.pack(">Q", len(payload)) + payload

    nonminimal_integer = node(b"I", b"\x00\x00")
    negative_zero = node(b"I", b"\x01")
    negative_zero_binary64 = node(b"D", struct.pack(">d", -0.0))
    nonfinite_binary64 = node(b"D", struct.pack(">d", math.inf))
    invalid_boolean = node(b"?", b"\x02")
    nonreduced_rational = node(b"R", encode_integer(2) + encode_integer(4))
    negative_denominator = node(b"R", encode_integer(1) + encode_integer(-2))
    wrong_sequence_count = node(b"S", struct.pack(">Q", 2) + encode_bytes(b"only-one"))
    wrong_union_count = node(b"T", encode_bytes(b"only-tag"))
    empty_union_tag = node(b"T", encode_bytes(b"") + encode_bytes(b"payload"))
    incomplete_map = node(b"M", struct.pack(">Q", 1) + encode_bytes(b"key"))
    empty_map_key = node(b"M", struct.pack(">Q", 1) + encode_bytes(b"") + encode_bytes(b"value"))
    unordered_map = node(
        b"M",
        struct.pack(">Q", 2) + encode_bytes(b"z") + encode_bytes(b"1") + encode_bytes(b"a") + encode_bytes(b"2"),
    )
    duplicate_map = node(
        b"M",
        struct.pack(">Q", 2) + encode_bytes(b"a") + encode_bytes(b"1") + encode_bytes(b"a") + encode_bytes(b"2"),
    )

    for malformed in (
        nonminimal_integer,
        negative_zero,
        negative_zero_binary64,
        nonfinite_binary64,
        invalid_boolean,
        nonreduced_rational,
        negative_denominator,
        wrong_sequence_count,
        wrong_union_count,
        empty_union_tag,
        incomplete_map,
        empty_map_key,
        unordered_map,
        duplicate_map,
        encode_integer(1) + b"trailing",
    ):
        with pytest.raises(CanonicalEncodingError):
            require_canonical_record(malformed)


def test_component_maps_ignore_insertion_order_but_bind_names() -> None:
    left = encode_component_map({b"plane": b"a", b"policy": b"b"})
    right = encode_component_map({b"policy": b"b", b"plane": b"a"})

    assert left == right
    assert left != encode_component_map({b"plane": b"b", b"policy": b"a"})


def test_polygon_encoding_normalizes_outer_start_winding_and_hole_order() -> None:
    outer = (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(8.0, 0.0),
        Point2[WorldXY].build(8.0, 8.0),
        Point2[WorldXY].build(0.0, 8.0),
    )
    hole_a = (
        Point2[WorldXY].build(1.0, 1.0),
        Point2[WorldXY].build(2.0, 1.0),
        Point2[WorldXY].build(2.0, 2.0),
        Point2[WorldXY].build(1.0, 2.0),
    )
    hole_b = (
        Point2[WorldXY].build(5.0, 5.0),
        Point2[WorldXY].build(7.0, 5.0),
        Point2[WorldXY].build(7.0, 7.0),
        Point2[WorldXY].build(5.0, 7.0),
    )

    canonical = canonical_polygon_bytes(outer, (hole_a, hole_b))
    reordered = canonical_polygon_bytes(
        tuple(reversed(outer[1:] + outer[:1])),
        (
            tuple(reversed(hole_b[2:] + hole_b[:2])),
            hole_a[3:] + hole_a[:3],
        ),
    )

    assert reordered == canonical


def test_canonical_ring_exhausts_rotation_reversal_closure_and_signed_zero() -> None:
    outer = (
        Point2[WorldXY].build(-0.0, 0.0),
        Point2[WorldXY].build(8.0, -0.0),
        Point2[WorldXY].build(8.0, 8.0),
        Point2[WorldXY].build(0.0, 8.0),
    )
    canonical = CanonicalRingV1.build_outer(outer)
    positive_zero = (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(8.0, 0.0),
        Point2[WorldXY].build(8.0, 8.0),
        Point2[WorldXY].build(0.0, 8.0),
    )

    assert CanonicalRingV1.build_outer(positive_zero) == canonical
    assert canonical.vertex_count == 4
    for source in (outer, tuple(reversed(outer))):
        for offset in range(len(source)):
            rotated = source[offset:] + source[:offset]
            assert CanonicalRingV1.build_outer(rotated) == canonical
            assert CanonicalRingV1.build_outer(rotated + (rotated[0],)) == canonical

    hole = (
        Point2[WorldXY].build(-0.0, 0.0),
        Point2[WorldXY].build(2.0, 0.0),
        Point2[WorldXY].build(2.0, 2.0),
        Point2[WorldXY].build(0.0, 2.0),
    )
    canonical_hole = CanonicalRingV1.build_hole(hole)
    assert canonical_hole.canonical_bytes != CanonicalRingV1.build_outer(hole).canonical_bytes
    for source in (hole, tuple(reversed(hole))):
        for offset in range(len(source)):
            rotated = source[offset:] + source[:offset]
            assert CanonicalRingV1.build_hole(rotated) == canonical_hole
            assert CanonicalRingV1.build_hole(rotated + (rotated[0],)) == canonical_hole


def test_canonical_polygon_exhausts_hole_permutations() -> None:
    outer = (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(10.0, 0.0),
        Point2[WorldXY].build(10.0, 10.0),
        Point2[WorldXY].build(0.0, 10.0),
    )
    holes = (
        (
            Point2[WorldXY].build(1.0, 1.0),
            Point2[WorldXY].build(2.0, 1.0),
            Point2[WorldXY].build(2.0, 2.0),
            Point2[WorldXY].build(1.0, 2.0),
        ),
        (
            Point2[WorldXY].build(4.0, 1.0),
            Point2[WorldXY].build(5.0, 1.0),
            Point2[WorldXY].build(5.0, 2.0),
            Point2[WorldXY].build(4.0, 2.0),
        ),
        (
            Point2[WorldXY].build(7.0, 1.0),
            Point2[WorldXY].build(8.0, 1.0),
            Point2[WorldXY].build(8.0, 2.0),
            Point2[WorldXY].build(7.0, 2.0),
        ),
    )
    canonical = canonical_polygon_bytes(outer, holes)

    assert all(canonical_polygon_bytes(outer, permutation) == canonical for permutation in permutations(holes))


def test_canonical_ring_raw_constructor_accepts_only_canonical_state() -> None:
    canonical = CanonicalRingV1.build_outer(
        (
            Point2[WorldXY].build(0.0, 0.0),
            Point2[WorldXY].build(2.0, 0.0),
            Point2[WorldXY].build(2.0, 2.0),
            Point2[WorldXY].build(0.0, 2.0),
        )
    )

    assert CanonicalRingV1(canonical.vertices, True) == canonical
    with pytest.raises(CanonicalEncodingError, match="canonical rotation"):
        CanonicalRingV1(canonical.vertices[1:] + canonical.vertices[:1], True)
    with pytest.raises(CanonicalEncodingError, match="exact bool"):
        CanonicalRingV1(canonical.vertices, 1)  # type: ignore[arg-type]
    with pytest.raises(CanonicalEncodingError, match="exact CanonicalRingV1"):
        SemanticRing(canonical.vertices, True, False)


def test_polygon_encoding_rejects_degenerate_rings() -> None:
    collinear = (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
        Point2[WorldXY].build(2.0, 0.0),
    )

    with pytest.raises(CanonicalEncodingError, match="nonzero exact area"):
        canonical_polygon_bytes(collinear, ())


def test_polygon_encoding_rejects_repeated_vertices_and_duplicate_holes() -> None:
    repeated = (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(4.0, 0.0),
        Point2[WorldXY].build(4.0, 4.0),
        Point2[WorldXY].build(0.0, 4.0),
        Point2[WorldXY].build(4.0, 0.0),
    )
    outer = (
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(8.0, 0.0),
        Point2[WorldXY].build(8.0, 8.0),
        Point2[WorldXY].build(0.0, 8.0),
    )
    hole = (
        Point2[WorldXY].build(1.0, 1.0),
        Point2[WorldXY].build(2.0, 1.0),
        Point2[WorldXY].build(2.0, 2.0),
        Point2[WorldXY].build(1.0, 2.0),
    )

    with pytest.raises(CanonicalEncodingError, match="repeated"):
        canonical_polygon_bytes(repeated, ())
    with pytest.raises(CanonicalEncodingError, match="duplicates"):
        canonical_polygon_bytes(outer, (hole, tuple(reversed(hole))))


def test_every_task1_motion_plane_cap_and_policy_has_canonical_bytes() -> None:
    start = Point2[WorldXY].build(1.0, 2.0)
    end = Point2[WorldXY].build(3.0, 4.0)
    cap = EngagementCap.build(math.pi / 2.0)
    plane = CutPlane.build(CutZ.build(-2.0), ClearanceZ.build(5.0))
    values = (
        plane,
        ExactSegmentMotion.build(start, end),
        ExactCircleMotion.build(start, Vector2[WorldXY].build(2.0, 0.0), False),
        cap,
        _candidate_policy(),
        _neck_policy(),
        DepletionPolicy.build(chord_bound=ChordBound.build(0.01), center_count_limit=256),
        TraversalPolicy.build(forward_window=6),
        CutDirectionPolicy.build(CutIntent.CLIMB),
    )

    encodings = tuple(canonical_task1_bytes(value) for value in values)

    assert all(encoded.startswith(CANONICAL_ENCODING_VERSION) for encoded in encodings)
    assert len(set(encodings)) == len(values)
    assert b"world-xy-v1" in canonical_task1_bytes(ExactSegmentMotion.build(start, end))


def test_task1_encoding_binds_each_policy_decision() -> None:
    policy = _candidate_policy()
    changed = CandidatePolicy.build(
        spatial_resolution=policy.spatial_resolution,
        spatial_refinement_levels=policy.spatial_refinement_levels,
        radius_resolution=policy.radius_resolution,
        radius_refinement_levels=policy.radius_refinement_levels,
        phase_count=policy.phase_count + 1,
        minimum_guide_radius=policy.minimum_guide_radius,
        minimum_progress=policy.minimum_progress,
    )

    assert canonical_task1_bytes(policy) != canonical_task1_bytes(changed)
    assert canonical_task1_bytes(TraversalPolicy.build(forward_window=1)) != canonical_task1_bytes(TraversalPolicy.build(forward_window=2))
    assert canonical_task1_bytes(CutDirectionPolicy.build(CutIntent.CLIMB)) != canonical_task1_bytes(CutDirectionPolicy.build(CutIntent.CONVENTIONAL))


def test_canonical_encoder_rejects_unsupported_values() -> None:
    with pytest.raises(CanonicalEncodingError, match="Task 1"):
        canonical_task1_bytes("climb")


def test_closed_canonical_variants_reject_semantic_subclasses() -> None:
    start = Point2[WorldXY].build(0.0, 0.0)
    end = Point2[WorldXY].build(1.0, 0.0)
    semantic_point = SemanticPoint(Millimetre(0.0), Millimetre(0.0), False)
    semantic_vector = SemanticVector(Millimetre(1.0), Millimetre(0.0), False)
    semantic_segment = SemanticSegment(start, end, False)
    semantic_circle = SemanticCircle(start, Vector2[WorldXY].build(1.0, 0.0), False, False)
    base_policy = _candidate_policy()
    semantic_policy = SemanticCandidatePolicy(
        base_policy.spatial_resolution,
        base_policy.spatial_refinement_levels,
        base_policy.radius_resolution,
        base_policy.radius_refinement_levels,
        base_policy.phase_count,
        base_policy.minimum_guide_radius,
        base_policy.minimum_progress,
        False,
    )

    with pytest.raises(CanonicalEncodingError, match="exact Point2"):
        canonical_point2_bytes(semantic_point)
    with pytest.raises(CanonicalEncodingError, match="exact Vector2"):
        canonical_vector2_bytes(semantic_vector)
    with pytest.raises(CanonicalEncodingError, match="exact CutZ"):
        canonical_cut_z_bytes(SemanticCutZ(Millimetre(-2.0), False))
    with pytest.raises(CanonicalEncodingError, match="exact ExactRationalV1"):
        _ = SemanticRational(1, 2, False).canonical_bytes
    with pytest.raises(CanonicalEncodingError, match="Task 1"):
        canonical_task1_bytes(semantic_segment)
    with pytest.raises(CanonicalEncodingError, match="Task 1"):
        canonical_task1_bytes(semantic_circle)
    with pytest.raises(CanonicalEncodingError, match="Task 1"):
        canonical_task1_bytes(semantic_policy)
    with pytest.raises(CanonicalEncodingError, match="Task 1"):
        canonical_task1_bytes(SemanticTraversalPolicy(4, False))


def test_nested_task1_fields_reject_semantic_subclasses() -> None:
    semantic_cut_z = SemanticCutZ(Millimetre(-2.0), False)
    cut_plane = CutPlane.build(semantic_cut_z, ClearanceZ.build(5.0))
    semantic_point = SemanticPoint(Millimetre(0.0), Millimetre(0.0), False)
    semantic_vector = SemanticVector(Millimetre(1.0), Millimetre(0.0), False)
    segment = ExactSegmentMotion.build(semantic_point, Point2[WorldXY].build(1.0, 0.0))
    circle = ExactCircleMotion.build(Point2[WorldXY].build(0.0, 0.0), semantic_vector, False)
    base_policy = _candidate_policy()
    candidate_policy = CandidatePolicy.build(
        spatial_resolution=SemanticSpacing(Millimetre(0.5), False),
        spatial_refinement_levels=base_policy.spatial_refinement_levels,
        radius_resolution=base_policy.radius_resolution,
        radius_refinement_levels=base_policy.radius_refinement_levels,
        phase_count=base_policy.phase_count,
        minimum_guide_radius=base_policy.minimum_guide_radius,
        minimum_progress=base_policy.minimum_progress,
    )

    for value in (cut_plane, segment, circle, candidate_policy):
        with pytest.raises(CanonicalEncodingError, match="exact"):
            canonical_task1_bytes(value)
