"""Pin the binary64 -> exact-rational contract the attestation view depends on.

`fractions.Fraction(d)` is an exact, independent oracle: a binary64 IS a dyadic
rational, so its exact value is representable with no approximation.
"""

from fractions import Fraction

from hypothesis import given, settings
from hypothesis import strategies as st

from compas_cgal._continuous_tea_2 import SegmentEventSource2

# Values a uniform float strategy will essentially never generate, and where a
# bit-decomposition bug would actually live.
EDGE_DOUBLES = [
    5e-324,  # smallest positive subnormal
    2.2250738585072014e-308,  # smallest positive normal
    1.0,
    0.5,
    2.0**52,
    2.0**53,
    2.0**53 + 2.0,  # first representable integer above 2**53
    0.0,  # +0.0; -0.0 is covered by test_signed_zero_lifts_to_positive_zero
    0.1,  # exactly representable in decimal, not in binary64
    1e308,
]


def _source(value: float) -> SegmentEventSource2:
    """Build a source whose x0 carries `value` and whose other fields are valid."""
    return SegmentEventSource2.from_binary64(value, 0.0, value + 1.0, 1.0, 1.0, 1.0)


def _assert_matches_oracle(value: float) -> None:
    exact = Fraction(value)
    rational = _source(value).x0
    assert int(rational.numerator) == exact.numerator
    assert int(rational.denominator) > 0
    assert int(rational.denominator) == exact.denominator


@given(
    st.floats(
        allow_nan=False,
        allow_infinity=False,
        min_value=-1e6,
        max_value=1e6,
    )
)
@settings(max_examples=400)
def test_binary64_lift_is_exact(value: float) -> None:
    _assert_matches_oracle(value)


def test_binary64_lift_is_exact_at_edges() -> None:
    for value in EDGE_DOUBLES:
        _assert_matches_oracle(value)
        _assert_matches_oracle(-value)


def test_signed_zero_lifts_to_positive_zero() -> None:
    """-0.0 and +0.0 denote the same rational; the lift must not carry the sign bit."""
    for value in (0.0, -0.0):
        rational = _source(value).x0
        assert int(rational.numerator) == 0
        assert int(rational.denominator) == 1
        assert not rational.numerator.startswith("-")


def test_binary64_lift_denominator_is_a_power_of_two() -> None:
    """Every binary64 is dyadic, so the reduced denominator is a power of two."""
    for value in EDGE_DOUBLES:
        denominator = int(_source(value).x0.denominator)
        assert denominator > 0
        assert denominator & (denominator - 1) == 0
