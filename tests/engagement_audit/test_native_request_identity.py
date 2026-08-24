from __future__ import annotations

import math
import hashlib

import numpy as np
import pytest

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


def _boundary(*, extent: float = 10.0) -> np.ndarray:
    return np.array(
        ((0.0, 0.0), (extent, 0.0), (extent, extent), (0.0, extent)),
        dtype=np.float64,
    )


def _hole(*, offset: float = 2.0) -> np.ndarray:
    return np.array(
        (
            (offset, offset),
            (offset, offset + 1.0),
            (offset + 1.0, offset + 1.0),
            (offset + 1.0, offset),
        ),
        dtype=np.float64,
    )


def _policy(
    *,
    cap: float = math.pi / 2.0,
    chord_bound: float = 0.0625,
    center_limit: int = 4096,
) -> _stock_2.AuditPolicy2:
    return _stock_2.build_audit_policy(
        2.0,
        cap,
        _stock_2.cap_chord_ratio(cap),
        chord_bound,
        center_limit,
    )


def _segment(*, end_x: float = 4.0) -> _stock_2.AuditSegmentMotion2:
    motion = _stock_2.classify_audit_line(
        (1.0, 1.0, 0.0),
        (end_x, 1.0, 0.0),
        0.0,
        5.0,
        "cut",
    )
    assert isinstance(motion, _stock_2.AuditSegmentMotion2)
    return motion


def test_policy_seals_the_authored_cap_observation() -> None:
    cap = 0.7
    policy = _stock_2.build_audit_policy(
        2.0,
        cap,
        _stock_2.cap_chord_ratio(cap),
        0.0625,
        4096,
    )

    assert len(policy.digest) == 32
    with pytest.raises(TypeError):
        _stock_2.AuditPolicy2()
    with pytest.raises(_stock_2.AuditPolicyCapSurrogateMismatchError):
        _stock_2.build_audit_policy(2.0, cap, math.nextafter(_stock_2.cap_chord_ratio(cap), math.inf), 0.0625, 4096)


@pytest.mark.parametrize("cap", [math.nextafter(math.pi, 0.0), math.pi])
def test_policy_accepts_upper_cap_boundary_and_adjacent_interior(cap: float) -> None:
    assert (
        len(
            _stock_2.build_audit_policy(
                2.0,
                cap,
                _stock_2.cap_chord_ratio(cap),
                0.0625,
                4096,
            ).digest
        )
        == 32
    )


@pytest.mark.parametrize("cap", [0.0, math.nextafter(math.pi, math.inf)])
def test_policy_rejects_cap_outside_closed_authored_domain(cap: float) -> None:
    with pytest.raises(_stock_2.AuditPolicyEngagementCapRangeError):
        _stock_2.build_audit_policy(2.0, cap, 1.0, 0.0625, 4096)


def test_policy_names_unrepresentable_positive_cap_surrogate() -> None:
    with pytest.raises(_stock_2.AuditPolicyCapSurrogateMismatchError):
        _stock_2.build_audit_policy(
            2.0,
            math.nextafter(0.0, math.inf),
            0.0,
            0.0625,
            4096,
        )


@pytest.mark.parametrize(
    ("tool_radius", "cap", "ratio", "chord_bound", "center_limit", "error"),
    (
        (math.nan, 0.7, 1.0, 0.0625, 4096, "AuditPolicyNonFiniteInputError"),
        (2.0, math.nan, math.nan, 0.0625, 4096, "AuditPolicyNonFiniteInputError"),
        (2.0, 0.7, math.inf, 0.0625, 4096, "AuditPolicyNonFiniteInputError"),
        (2.0, 0.7, 1.0, math.inf, 4096, "AuditPolicyNonFiniteInputError"),
        (0.0, 0.7, 1.0, 0.0625, 4096, "AuditPolicyToolRadiusError"),
        (-1.0, 0.7, 1.0, 0.0625, 4096, "AuditPolicyToolRadiusError"),
        (2.0, 0.7, 1.0, 0.0, 4096, "AuditPolicyDepletionChordBoundError"),
        (2.0, 0.7, 1.0, -0.0625, 4096, "AuditPolicyDepletionChordBoundError"),
        (2.0, 0.7, 1.0, 2.0, 4096, "AuditPolicyDepletionChordBoundError"),
        (2.0, 0.7, 1.0, 0.0625, 0, "AuditPolicyCenterCountLimitError"),
        (2.0, 0.7, 1.0, 0.0625, -1, "AuditPolicyCenterCountLimitError"),
        (2.0, 0.7, 1.0, 0.0625, True, "AuditPolicyCenterCountLimitError"),
    ),
)
def test_policy_rejects_each_invalid_input_with_named_error(
    tool_radius: float,
    cap: float,
    ratio: float,
    chord_bound: float,
    center_limit: int | bool,
    error: str,
) -> None:
    supplied_ratio = _stock_2.cap_chord_ratio(cap) if ratio == 1.0 and math.isfinite(cap) else ratio
    with pytest.raises(getattr(_stock_2, error)):
        _stock_2.build_audit_policy(
            tool_radius,
            cap,
            supplied_ratio,
            chord_bound,
            center_limit,
        )


def test_native_request_recomputes_stock_policy_motion_and_order_identity() -> None:
    first = _segment(end_x=4.0)
    second = _segment(end_x=5.0)
    baseline = _stock_2.build_audit_native_request_identity(
        _boundary(),
        [_hole()],
        _policy(),
        (first, second),
    )
    mutations = (
        _stock_2.build_audit_native_request_identity(_boundary(extent=11.0), [_hole()], _policy(), (first, second)),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole(offset=3.0)], _policy(), (first, second)),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(chord_bound=0.03125), (first, second)),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(cap=0.7), (first, second)),
        _stock_2.build_audit_native_request_identity(
            _boundary(),
            [_hole()],
            _stock_2.build_audit_policy(2.5, math.pi / 2.0, _stock_2.cap_chord_ratio(math.pi / 2.0), 0.0625, 4096),
            (first, second),
        ),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(center_limit=2048), (first, second)),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(), (_segment(end_x=6.0), second)),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(), (first, _segment(end_x=6.0))),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(), (first,)),
        _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], _policy(), (second, first)),
    )

    assert len(baseline.digest) == 32
    assert len({baseline.digest, *(mutation.digest for mutation in mutations)}) == 1 + len(mutations)
    with pytest.raises(TypeError):
        _stock_2.AuditNativeRequestIdentity2()


def test_native_request_rejects_empty_list_and_foreign_motion_sequence() -> None:
    with pytest.raises(_stock_2.AuditNativeRequestMotionError):
        _stock_2.build_audit_native_request_identity(_boundary(), [], _policy(), ())
    with pytest.raises(_stock_2.AuditNativeRequestMotionError):
        _stock_2.build_audit_native_request_identity(_boundary(), [], _policy(), (_segment(), object()))
    with pytest.raises(TypeError):
        _stock_2.build_audit_native_request_identity(_boundary(), [], _policy(), [_segment()])


def test_native_request_ccan_binds_stock_policy_and_ordered_motion_digests() -> None:
    first = _segment(end_x=4.0)
    second = _segment(end_x=5.0)
    stock = _stock_2.build_audit_native_stock_identity(_boundary(), [_hole()])
    policy = _policy()
    request = _stock_2.build_audit_native_request_identity(_boundary(), [_hole()], policy, (first, second))
    expected = encode_tagged_union(
        b"audit-native-request-v1",
        encode_component_map(
            {
                b"native-motion-digests": encode_sequence((first.digest, second.digest)),
                b"policy-digest": policy.digest,
                b"stock-digest": stock.digest,
            }
        ),
    )

    assert request.canonical_bytes == expected
    assert request.digest == hashlib.sha256(expected).digest()


def test_native_stock_identity_normalizes_ring_and_hole_order() -> None:
    boundary = _boundary()
    first_hole = _hole(offset=2.0)
    second_hole = _hole(offset=5.0)
    motion = _segment()

    baseline = _stock_2.build_audit_native_request_identity(
        boundary,
        [first_hole, second_hole],
        _policy(),
        (motion,),
    )
    equivalent = _stock_2.build_audit_native_request_identity(
        np.roll(boundary[::-1], 1, axis=0),
        [np.roll(second_hole[::-1], 2, axis=0), np.roll(first_hole[::-1], 1, axis=0)],
        _policy(),
        (motion,),
    )

    assert baseline.digest == equivalent.digest


def _canonical_ring(points: np.ndarray, *, outer: bool) -> CanonicalRingV1:
    typed = tuple(Point2[WorldXY].build(float(x), float(y)) for x, y in points)
    return CanonicalRingV1.build_outer(typed) if outer else CanonicalRingV1.build_hole(typed)


def test_native_stock_identity_matches_cross_language_ccan_golden() -> None:
    boundary = _boundary()
    hole = _hole()
    identity = _stock_2.build_audit_native_stock_identity(boundary, [hole])
    outer = _canonical_ring(boundary, outer=True)
    inner = _canonical_ring(hole, outer=False)
    expected = encode_tagged_union(
        b"audit-native-stock-v1",
        encode_component_map(
            {
                b"boundary": outer.canonical_bytes,
                b"holes": encode_sequence((inner.canonical_bytes,)),
            }
        ),
    )

    assert identity.canonical_bytes == expected
    assert identity.digest == hashlib.sha256(expected).digest()
    with pytest.raises(TypeError):
        _stock_2.AuditNativeStockIdentity2()


def test_native_stock_identity_normalizes_signed_zero() -> None:
    positive = _stock_2.build_audit_native_stock_identity(_boundary(), [])
    negative_zero = _boundary()
    negative_zero[0, 0] = -0.0

    assert _stock_2.build_audit_native_stock_identity(negative_zero, []).digest == positive.digest


def test_native_stock_identity_normalizes_one_repeated_closing_vertex() -> None:
    open_ring = _boundary()
    closed_ring = np.vstack((open_ring, open_ring[0]))

    assert _stock_2.build_audit_native_stock_identity(open_ring, []).digest == _stock_2.build_audit_native_stock_identity(closed_ring, []).digest


@pytest.mark.parametrize(
    "invalid",
    (
        np.array(((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)), dtype=np.float64),
        np.array(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (1.0, 0.0)), dtype=np.float64),
        np.array(((0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)), dtype=np.float64),
    ),
)
def test_native_stock_identity_rejects_degenerate_repeated_or_self_intersecting_ring(
    invalid: np.ndarray,
) -> None:
    with pytest.raises(_stock_2.AuditNativeStockRingError):
        _stock_2.build_audit_native_stock_identity(invalid, [])


def test_native_stock_identity_rejects_nonfinite_shape_and_duplicate_holes() -> None:
    nonfinite = _boundary()
    nonfinite[1, 0] = math.inf
    wrong_columns = np.zeros((4, 3), dtype=np.float64)
    hole = _hole()

    with pytest.raises(_stock_2.AuditNativeStockNonFiniteInputError):
        _stock_2.build_audit_native_stock_identity(nonfinite, [])
    with pytest.raises(_stock_2.AuditNativeStockShapeError):
        _stock_2.build_audit_native_stock_identity(wrong_columns, [])
    with pytest.raises(_stock_2.AuditNativeStockDuplicateHoleError):
        _stock_2.build_audit_native_stock_identity(_boundary(), [hole, hole.copy()])
