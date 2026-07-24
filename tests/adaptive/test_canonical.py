import math
import struct
from fractions import Fraction

import pytest

from compas_cgal.adaptive.canonical import CANONICAL_ENCODING_VERSION
from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import canonical_polygon_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
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
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


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
