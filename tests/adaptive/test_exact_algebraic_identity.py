from __future__ import annotations

from compas_cgal import _continuous_tea_2
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union

PRIMITIVE_FACTOR = ("1", "0", "-5", "0", "6")
EQUIVALENT_POLYNOMIALS = (
    PRIMITIVE_FACTOR,
    ("-1", "0", "5", "0", "-6"),
    ("7", "0", "-35", "0", "42"),
)


def _root_id(coefficients: tuple[int, ...], ordinal: int) -> bytes:
    return encode_tagged_union(
        b"algebraic-root-id-v1",
        encode_component_map(
            {
                b"coefficients": encode_sequence(tuple(encode_integer(value) for value in coefficients)),
                b"root-ordinal": encode_integer(ordinal),
            }
        ),
    )


def _event(suffix: bytes) -> _continuous_tea_2.PartitionEvent2:
    return _continuous_tea_2.PartitionEvent2(
        "tangency",
        b"feature-" + suffix,
        b"support-" + suffix,
        b"trim-" + suffix,
        b"vertex-" + suffix,
        b"branch-" + suffix,
        "accepted",
    )


def _projection(
    identifier: str,
    coefficients: tuple[str, ...],
    event: _continuous_tea_2.PartitionEvent2,
) -> _continuous_tea_2.ProjectionInput2:
    return _continuous_tea_2.ProjectionInput2(
        identifier,
        coefficients,
        (event,),
    )


def _root_signature(
    certificate: _continuous_tea_2.EventPartitionCertificate2,
) -> tuple[tuple[bytes, tuple[str, ...], int], ...]:
    return tuple(
        (
            root.root_id,
            root.factor_coefficients,
            root.root_ordinal,
        )
        for root in certificate.roots
    )


def _cell_signature(
    certificate: _continuous_tea_2.EventPartitionCertificate2,
) -> tuple[tuple[bytes, bytes, str, str], ...]:
    return tuple(
        (
            cell.lower_root_id,
            cell.upper_root_id,
            cell.witness_numerator,
            cell.witness_denominator,
        )
        for cell in certificate.cells
    )


def _fibre_signature(
    certificate: _continuous_tea_2.EventPartitionCertificate2,
) -> tuple[tuple[bytes, tuple[bytes, ...]], ...]:
    return tuple(
        (
            fibre.root_id,
            tuple(event.feature_id for event in fibre.events),
        )
        for fibre in certificate.fibres
    )


def test_scaled_and_sign_flipped_polynomials_share_exact_identity() -> None:
    certificates = tuple(_continuous_tea_2.partition_projections((_projection("equivalent", coefficients, _event(b"same")),)) for coefficients in EQUIVALENT_POLYNOMIALS)
    expected_roots = (
        (_root_id((1, 0, -5, 0, 6), 2), PRIMITIVE_FACTOR, 2),
        (_root_id((1, 0, -5, 0, 6), 3), PRIMITIVE_FACTOR, 3),
    )
    expected_cells = (
        (b"", expected_roots[0][0], "1", "2"),
        (expected_roots[0][0], expected_roots[1][0], "5", "8"),
        (expected_roots[1][0], b"", "3", "4"),
    )

    assert all(certificate.projections[0].factor_coefficients == (PRIMITIVE_FACTOR,) for certificate in certificates)
    assert all(_root_signature(certificate) == expected_roots for certificate in certificates)
    assert all(_cell_signature(certificate) == expected_cells for certificate in certificates)


def test_equivalent_projection_merge_is_input_order_deterministic() -> None:
    projections = (
        _projection("positive", PRIMITIVE_FACTOR, _event(b"positive")),
        _projection(
            "negative",
            EQUIVALENT_POLYNOMIALS[1],
            _event(b"negative"),
        ),
        _projection(
            "scaled",
            EQUIVALENT_POLYNOMIALS[2],
            _event(b"scaled"),
        ),
    )
    forward = _continuous_tea_2.partition_projections(projections)
    reverse = _continuous_tea_2.partition_projections(tuple(reversed(projections)))

    assert _root_signature(forward) == _root_signature(reverse)
    assert _cell_signature(forward) == _cell_signature(reverse)
    assert _fibre_signature(forward) == _fibre_signature(reverse)
    assert all(
        feature_ids
        == (
            b"feature-negative",
            b"feature-positive",
            b"feature-scaled",
        )
        for _, feature_ids in _fibre_signature(forward)
    )
