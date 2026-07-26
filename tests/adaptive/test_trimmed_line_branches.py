from __future__ import annotations

from fractions import Fraction

import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union


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


def test_native_trim_filter_retains_only_closed_branch_domains() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        ("3/5", "0"),
        ("3/5", "1"),
        ("0", "0", "0", "2"),
        "1",
        "rim-half-0-v1",
    )

    assert [
        (
            branch.rim_parameter,
            branch.motion_domain_low,
            branch.motion_domain_high,
            branch.trim_disposition,
        )
        for branch in branches
    ] == [
        ("1/2", "0", "1/10", "accepted"),
        ("-1/2", "2/5", "9/10", "accepted"),
    ]
    assert all(branch.rejected_outside_closed_domain for branch in branches)
    assert all(branch.feature_id and branch.trim_id and branch.branch_id for branch in branches)
    assert [
        (
            branch.rim_root.root_id,
            branch.rim_root.factor_coefficients,
            branch.rim_root.root_ordinal,
            branch.rim_root.multiplicity,
        )
        for branch in branches
    ] == [
        (_root_id((-1, 0, 4), 1), ("-1", "0", "4"), 1, 1),
        (_root_id((-1, 0, 4), 0), ("-1", "0", "4"), 0, 1),
    ]
    assert all(branch.lower_trim_predicate_rows == (("0", "2", "0"), ("2", "0", "2")) for branch in branches)
    assert all(branch.upper_trim_predicate_rows == (("1", "-2", "1"), ("-2", "0", "-2")) for branch in branches)
    assert all(branch.rational_convenience_available for branch in branches)
    assert all(branch.domain_nonempty_rechecked for branch in branches)
    assert all(branch.complementarity_rechecked for branch in branches)


def test_native_trim_filter_mirrors_exact_domains_on_second_rim_chart() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "3"),
        ("-3/5", "0"),
        ("-3/5", "1"),
        ("0", "0", "0", "2"),
        "1",
        "rim-half-1-v1",
    )

    assert [
        (
            branch.rim_parameter,
            branch.motion_domain_low,
            branch.motion_domain_high,
        )
        for branch in branches
    ] == [
        ("1/2", "2/5", "9/10"),
        ("-1/2", "0", "1/10"),
    ]
    assert [branch.rim_root.root_ordinal for branch in branches] == [1, 0]
    assert all(branch.domain_nonempty_rechecked for branch in branches)
    assert all(branch.complementarity_rechecked for branch in branches)


def test_rational_trim_convenience_matches_legacy_fields() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        ("3/5", "0"),
        ("3/5", "1"),
        ("0", "0", "0", "2"),
        "1",
        "rim-half-0-v1",
    )
    expected = (
        ("1/2", "0", "1/10"),
        ("-1/2", "2/5", "9/10"),
    )

    for branch, values in zip(branches, expected, strict=True):
        convenience = branch.rational_convenience
        assert convenience is not None
        assert (
            convenience.rim_parameter,
            convenience.motion_domain_low,
            convenience.motion_domain_high,
        ) == values
        assert (
            branch.rim_parameter,
            branch.motion_domain_low,
            branch.motion_domain_high,
        ) == values
        assert branch.rational_convenience_available


def test_native_trim_filter_retains_nonsquare_root_and_signed_predicates() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("1", "1", "0"),
        ("0", "0"),
        ("1", "-1"),
        ("0", "0", "1", "-1"),
        "1",
        "rim-half-0-v1",
    )

    assert len(branches) == 1
    branch = branches[0]
    assert branch.rim_parameter == ""
    assert branch.motion_domain_low == ""
    assert branch.motion_domain_high == ""
    assert branch.rim_root.root_id == _root_id((-1, -2, 1), 0)
    assert branch.rim_root.factor_coefficients == ("-1", "-2", "1")
    assert branch.rim_root.root_ordinal == 0
    assert branch.rim_root.multiplicity == 1
    interval_low = Fraction(branch.rim_root.interval_low)
    interval_high = Fraction(branch.rim_root.interval_high)
    assert interval_low < interval_high < 1
    factor_at_low = interval_low**2 - 2 * interval_low - 1
    factor_at_high = interval_high**2 - 2 * interval_high - 1
    assert factor_at_low * factor_at_high < 0
    assert branch.lower_trim_predicate_rows == (
        ("1", "0", "-1"),
        ("1", "0", "1"),
    )
    assert branch.upper_trim_predicate_rows == (
        ("0", "0", "2"),
        ("-1", "0", "-1"),
    )
    assert not branch.rational_convenience_available
    assert branch.domain_nonempty_rechecked
    assert branch.complementarity_rechecked


def test_irrational_trim_convenience_is_absent_with_empty_legacy_fields() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("1", "1", "0"),
        ("0", "0"),
        ("1", "-1"),
        ("0", "0", "1", "-1"),
        "1",
        "rim-half-0-v1",
    )

    assert len(branches) == 1
    branch = branches[0]
    assert branch.rational_convenience is None
    assert (
        branch.rim_parameter,
        branch.motion_domain_low,
        branch.motion_domain_high,
    ) == ("", "", "")
    assert not branch.rational_convenience_available


def test_native_trim_filter_rechecks_oblique_contact_quadratic_coefficient() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("1", "1", "-7/5"),
        ("0", "7/5"),
        ("1", "2/5"),
        ("0", "0", "1", "-1"),
        "1",
        "rim-half-0-v1",
    )

    branch = next(branch for branch in branches if branch.rim_parameter == "1/2")
    assert branch.trim_disposition == "accepted"
    assert branch.domain_nonempty_rechecked
    assert branch.complementarity_rechecked


def test_native_trim_filter_returns_empty_when_support_equation_has_no_real_root() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("1", "0", "-2"),
        ("2", "0"),
        ("2", "1"),
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )

    assert not branches


def test_native_trim_filter_rejects_zero_radius_before_overlap_classification() -> None:
    with pytest.raises(
        _continuous_tea_2.NonPositiveToolRadiusError,
        match="positive",
    ):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("0", "0"),
            ("0", "1"),
            ("0", "0", "0", "1"),
            "0",
            "rim-half-0-v1",
        )


def test_native_trim_filter_includes_chart_boundaries_and_excludes_outside_roots() -> None:
    boundary_branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("1", "0", "0"),
        ("0", "-2"),
        ("0", "3"),
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )

    assert [branch.rim_parameter for branch in boundary_branches] == ["1", "-1"]
    assert all(branch.rational_convenience_available for branch in boundary_branches)
    negative, positive = sorted(
        boundary_branches,
        key=lambda branch: branch.rim_root.root_ordinal,
    )
    assert Fraction(negative.rim_root.interval_high) < 1
    assert Fraction(positive.rim_root.interval_low) > -1

    outside = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "3"),
        ("-3/5", "-2"),
        ("-3/5", "3"),
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )
    assert not outside


def test_native_trim_identity_is_scale_invariant_and_collision_complete() -> None:
    trim_start = ("3/5", "-3")
    trim_end = ("3/5", "3")
    canonical = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        trim_start,
        trim_end,
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )
    scaled = _continuous_tea_2.solve_trimmed_line_branches(
        ("-10", "0", "6"),
        trim_start,
        trim_end,
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )

    assert [branch.feature_id for branch in canonical] == [branch.feature_id for branch in scaled]
    assert all(branch.feature_id == branch.local_support_id for branch in (*canonical, *scaled))
    assert [branch.local_support_id for branch in canonical] == [branch.local_support_id for branch in scaled]
    assert [branch.local_trimmed_feature_id for branch in canonical] == [branch.local_trimmed_feature_id for branch in scaled]
    assert [branch.branch_id for branch in canonical] == [branch.branch_id for branch in scaled]

    different_trim = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        ("3/5", "-2"),
        trim_end,
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )
    reversed_trim = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        trim_end,
        trim_start,
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )
    shifted_motion = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        trim_start,
        trim_end,
        ("0", "1", "0", "2"),
        "1",
        "rim-half-0-v1",
    )
    larger_radius = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        trim_start,
        trim_end,
        ("-3/5", "0", "-3/5", "1"),
        "2",
        "rim-half-0-v1",
    )

    canonical_ids = {branch.rim_root.root_id: branch.branch_id for branch in canonical}
    shifted_ids = {branch.rim_root.root_id: branch.branch_id for branch in shifted_motion}
    larger_radius_ids = {branch.rim_root.root_id: branch.branch_id for branch in larger_radius}
    different_trim_ids = {branch.rim_root.root_id: branch.branch_id for branch in different_trim}
    reversed_trim_ids = {branch.rim_root.root_id: branch.branch_id for branch in reversed_trim}
    canonical_local_ids = {branch.rim_root.root_id: branch.local_trimmed_feature_id for branch in canonical}
    different_trim_local_ids = {branch.rim_root.root_id: branch.local_trimmed_feature_id for branch in different_trim}
    reversed_trim_local_ids = {branch.rim_root.root_id: branch.local_trimmed_feature_id for branch in reversed_trim}
    assert canonical_ids.keys() == shifted_ids.keys() == larger_radius_ids.keys()
    assert canonical_ids.keys() == different_trim_ids.keys() == reversed_trim_ids.keys()
    assert {branch.local_support_id for branch in canonical} == {branch.local_support_id for branch in different_trim} == {branch.local_support_id for branch in reversed_trim}
    assert all(canonical_ids[root_id] != shifted_ids[root_id] for root_id in canonical_ids)
    assert all(canonical_ids[root_id] != larger_radius_ids[root_id] for root_id in canonical_ids)
    assert all(canonical_local_ids[root_id] != different_trim_local_ids[root_id] and canonical_ids[root_id] != different_trim_ids[root_id] for root_id in canonical_ids)
    assert all(canonical_local_ids[root_id] != reversed_trim_local_ids[root_id] and canonical_ids[root_id] != reversed_trim_ids[root_id] for root_id in canonical_ids)


@pytest.mark.parametrize("radius", ("0", "-1"))
def test_native_trim_filter_rejects_nonpositive_radius(radius: str) -> None:
    with pytest.raises(_continuous_tea_2.NonPositiveToolRadiusError):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("0", "0"),
            ("0", "1"),
            ("0", "0", "0", "1"),
            radius,
            "rim-half-0-v1",
        )


def test_native_trim_filter_rejects_zero_line_normal() -> None:
    with pytest.raises(_continuous_tea_2.DegenerateBoundarySupportError):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("0", "0", "1"),
            ("0", "0"),
            ("0", "1"),
            ("0", "0", "0", "1"),
            "1",
            "rim-half-0-v1",
        )


def test_native_trim_filter_rejects_unknown_rim_chart() -> None:
    with pytest.raises(_continuous_tea_2.ChartCoverageError, match="chart"):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("0", "0"),
            ("0", "1"),
            ("1", "0", "1", "1"),
            "1",
            "unknown-chart",
        )


@pytest.mark.parametrize(
    (
        "line_support",
        "trim_start",
        "trim_end",
        "segment_motion",
        "cutter_radius",
    ),
    (
        (("1", "0"), ("0", "0"), ("0", "1"), ("1", "0", "1", "1"), "1"),
        (("1", "0", "0"), ("0",), ("0", "1"), ("1", "0", "1", "1"), "1"),
        (("1", "0", "0"), ("0", "0"), ("0", "1"), ("1", "0", "1"), "1"),
        (("1", "0", "0"), ("0", "0"), ("0", "1"), ("1", "0", "1", "1"), "1/0"),
        (("1", "0", "0"), ("0", "0"), ("0", "1//2"), ("1", "0", "1", "1"), "1"),
    ),
)
def test_native_trim_filter_rejects_malformed_exact_inputs(
    line_support: tuple[str, ...],
    trim_start: tuple[str, ...],
    trim_end: tuple[str, ...],
    segment_motion: tuple[str, ...],
    cutter_radius: str,
) -> None:
    with pytest.raises(_continuous_tea_2.TrimFilterError):
        _continuous_tea_2.solve_trimmed_line_branches(
            line_support,
            trim_start,
            trim_end,
            segment_motion,
            cutter_radius,
            "rim-half-0-v1",
        )


def test_native_trim_filter_rejects_endpoint_off_exact_support() -> None:
    assert issubclass(
        _continuous_tea_2.TrimEndpointOffSupportError,
        _continuous_tea_2.TrimFilterError,
    )
    with pytest.raises(
        _continuous_tea_2.TrimEndpointOffSupportError,
        match="endpoint",
    ):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("2", "0"),
            ("2", "1"),
            ("5", "0", "5", "1"),
            "1",
            "rim-half-0-v1",
        )


def test_native_trim_identity_normalizes_rational_aliases() -> None:
    canonical = _continuous_tea_2.solve_trimmed_line_branches(
        ("5", "0", "-3"),
        ("3/5", "-3"),
        ("3/5", "3"),
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )
    aliases = _continuous_tea_2.solve_trimmed_line_branches(
        ("10/2", "0/1", "-6/2"),
        ("6/10", "-6/2"),
        ("6/10", "6/2"),
        ("0/2", "0/7", "0/9", "2/2"),
        "2/2",
        "rim-half-0-v1",
    )

    def identities(
        branches: list[_continuous_tea_2.TrimmedLineBranch2],
    ) -> dict[bytes, tuple[bytes, bytes, bytes]]:
        return {
            branch.rim_root.root_id: (
                branch.local_support_id,
                branch.local_trimmed_feature_id,
                branch.branch_id,
            )
            for branch in branches
        }

    assert identities(canonical) == identities(aliases)


def test_native_trim_filter_rejects_zero_length_motion() -> None:
    with pytest.raises(
        _continuous_tea_2.ZeroLengthSegmentMotionError,
        match="distinct",
    ):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("0", "0"),
            ("0", "1"),
            ("0", "0", "0", "0"),
            "1",
            "rim-half-0-v1",
        )


def test_native_trim_filter_rejects_motion_not_parallel_to_support() -> None:
    assert issubclass(
        _continuous_tea_2.UnsupportedLineMotionError,
        _continuous_tea_2.TrimFilterError,
    )
    with pytest.raises(
        _continuous_tea_2.UnsupportedLineMotionError,
        match="parallel to support",
    ):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("0", "0"),
            ("0", "1"),
            ("0", "0", "1", "0"),
            "1",
            "rim-half-0-v1",
        )


def test_native_trim_filter_rejects_coincident_trim_before_no_root_return() -> None:
    assert issubclass(
        _continuous_tea_2.DegenerateTrimError,
        _continuous_tea_2.TrimFilterError,
    )
    with pytest.raises(_continuous_tea_2.DegenerateTrimError, match="distinct"):
        _continuous_tea_2.solve_trimmed_line_branches(
            ("1", "0", "0"),
            ("0", "0"),
            ("0", "0"),
            ("1", "0", "1", "1"),
            "1",
            "rim-half-0-v1",
        )


def test_native_trim_filter_preserves_repeated_tangency_multiplicity() -> None:
    branches = _continuous_tea_2.solve_trimmed_line_branches(
        ("1", "0", "-1"),
        ("1", "-2"),
        ("1", "3"),
        ("0", "0", "0", "1"),
        "1",
        "rim-half-0-v1",
    )

    assert len(branches) == 1
    branch = branches[0]
    assert branch.rim_root.factor_coefficients == ("0", "1")
    assert branch.rim_root.root_ordinal == 0
    assert branch.rim_root.multiplicity == 2
    assert branch.motion_domain_low == "0"
    assert branch.motion_domain_high == "1"
