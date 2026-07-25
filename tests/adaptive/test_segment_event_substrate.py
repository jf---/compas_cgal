from __future__ import annotations

import hashlib
import math

import pytest

from compas_cgal import _continuous_tea_2


def test_segment_source_exact_lifts_each_binary64_once() -> None:
    source = _continuous_tea_2.SegmentEventSource2.from_binary64(
        0.1,
        -0.5,
        1.25,
        2.0,
        0.25,
        4.0,
    )

    assert (source.x0.numerator, source.x0.denominator) == (
        "3602879701896397",
        "36028797018963968",
    )
    assert (source.y0.numerator, source.y0.denominator) == ("-1", "2")
    assert (source.x1.numerator, source.x1.denominator) == ("5", "4")
    assert (source.y1.numerator, source.y1.denominator) == ("2", "1")
    assert (source.tool_radius.numerator, source.tool_radius.denominator) == ("1", "4")
    assert (source.cap_chord_ratio.numerator, source.cap_chord_ratio.denominator) == ("4", "1")
    assert source.motion_data == (
        "3602879701896397/36028797018963968",
        "-1/2",
        "5/4",
        "2",
    )
    assert source.canonical_digest == hashlib.sha256(source.canonical_bytes).digest()


@pytest.mark.parametrize(
    ("values", "error"),
    (
        ((math.nan, 0.0, 1.0, 0.0, 0.25, 1.0), "NonFiniteSegmentInputError"),
        ((0.0, 0.0, math.inf, 0.0, 0.25, 1.0), "NonFiniteSegmentInputError"),
        ((0.0, 0.0, 0.0, 0.0, 0.25, 1.0), "ZeroLengthSegmentMotionError"),
        ((0.0, 0.0, 1.0, 0.0, 0.0, 1.0), "NonPositiveToolRadiusError"),
        ((0.0, 0.0, 1.0, 0.0, -0.25, 1.0), "NonPositiveToolRadiusError"),
        ((0.0, 0.0, 1.0, 0.0, 0.25, 0.0), "InvalidCapChordRatioError"),
        ((0.0, 0.0, 1.0, 0.0, 0.25, 4.000000000000001), "InvalidCapChordRatioError"),
    ),
)
def test_segment_source_rejects_each_invalid_boundary_with_named_error(
    values: tuple[float, float, float, float, float, float],
    error: str,
) -> None:
    with pytest.raises(getattr(_continuous_tea_2, error)):
        _continuous_tea_2.SegmentEventSource2.from_binary64(*values)
