"""Fail-loud contracts of the coverage certificate's exact search loops (review C1, I1, q17).

- C1: the circle-sweep start-witness search is bounded by its theorem bound
  1 + 2 (2 |counted| + |target circles| + |lines|); an uncut coincident support
  (every witness on a support) raises ``StartWitnessSearchError`` instead of
  hanging.
- I1a: every exact sign whose certified enclosure excludes 0 must agree with
  that enclosure; a wrong sign rule raises ``SignRuleContractError`` at once.
- I1b: the P3 chart bisection is bounded by a root-separation bound derived
  from the operands' bit sizes, not a fixed 2^16; a non-separable pair (a
  contract breach: lo == hi) raises ``SeparationBoundExceededError``.
- q17: a vertex disk never carries an extension over a zero radicand or a zero
  coefficient.
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _coverage_certificate_2 as native

Text = tuple[str, str]


def _text(value: Fraction) -> Text:
    return str(value.numerator), str(value.denominator)


def _sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def _sign_plus_root(a: Fraction, b: Fraction, d: Fraction) -> int:
    """Exact sign of a + b sqrt(d), d >= 0, by squaring (independent oracle)."""
    sa, sb = _sign(a), _sign(b) if d > 0 else 0
    if sb == 0:
        return sa
    if sa == 0 or sa == sb:
        return sb
    return sa if a * a > b * b * d else (sb if a * a < b * b * d else 0)


# ---------------------------------------------------------------------------
# C1: start-witness search
# ---------------------------------------------------------------------------


def test_native_error_hierarchy_names_every_contract_breach() -> None:
    for name in ("StartWitnessSearchError", "SignRuleContractError", "SeparationBoundExceededError"):
        assert issubclass(getattr(native, name), native.CoverageCertificateNativeError)


def test_start_witness_search_with_an_uncut_coincident_support_raises_not_hangs() -> None:
    with pytest.raises(native.StartWitnessSearchError, match="bound 3"):
        native.contract_probe_start_witness(2.0, -1.0, 1.5, True)


def test_start_witness_search_returns_within_the_bound_on_a_crossing_support() -> None:
    attempts, bound = native.contract_probe_start_witness(2.0, -1.0, 1.5, False)
    assert bound == 3
    assert 1 <= attempts <= bound


# ---------------------------------------------------------------------------
# I1a: sign rules cross-checked against certified enclosures
# ---------------------------------------------------------------------------


def test_a_sign_contradicting_a_zero_free_enclosure_raises() -> None:
    with pytest.raises(native.SignRuleContractError):
        native.contract_probe_sign_check(-1, 1.0, 2.0)
    with pytest.raises(native.SignRuleContractError):
        native.contract_probe_sign_check(0, -2.0, -1.0)


def test_a_sign_inside_a_zero_straddling_enclosure_is_not_judged() -> None:
    for claimed in (-1, 0, 1):
        native.contract_probe_sign_check(claimed, -1.0, 1.0)
    native.contract_probe_sign_check(1, 1.0, 2.0)


_SMALL = st.fractions(min_value=-50, max_value=50, max_denominator=64)
_RADICAND = st.fractions(min_value=0, max_value=50, max_denominator=64)


@settings(max_examples=300, deadline=None)
@given(_SMALL, _SMALL, _RADICAND, _SMALL, _RADICAND)
def test_two_radical_sign_rule_matches_an_independent_oracle(u: Fraction, v: Fraction, p: Fraction, w: Fraction, q: Fraction) -> None:
    # sign(u + v sqrt(p) + w sqrt(q)): oracle by the same mixed-radical rule in Fractions.
    first = _sign_plus_root(u, v, p)
    sw = _sign(w) if q > 0 else 0
    if sw == 0:
        expected = first
    elif first == 0 or first == sw:
        expected = first or sw
    else:
        square = _sign_plus_root(u * u + v * v * p - w * w * q, 2 * u * v, p)
        expected = first if square > 0 else (sw if square < 0 else 0)
    assert native.contract_probe_sign_rational_two(*map(_text, (u, v, p, w, q))) == expected
    approx = float(u) + float(v) * math.sqrt(float(p)) + float(w) * math.sqrt(float(q))
    if abs(approx) > 1e-9:
        assert expected == (1 if approx > 0 else -1)


# ---------------------------------------------------------------------------
# I1b: P3 bisection bounded by a root-separation bound
# ---------------------------------------------------------------------------


def test_chart_bisection_separates_values_closer_than_two_to_the_minus_300() -> None:
    gap = Fraction(1, 2**300)
    q = Fraction(*map(int, native.contract_probe_rational_between(_text(Fraction(0)), _text(Fraction(1)), _text(Fraction(2)), _text(gap), _text(Fraction(1)), _text(Fraction(2)))))
    # lo = sqrt(2) < q < sqrt(2) + gap = hi.
    assert q * q > 2
    assert _sign_plus_root(gap - q, Fraction(1), Fraction(2)) > 0


def test_chart_bisection_of_a_non_separable_pair_raises_its_named_bound() -> None:
    lo = (_text(Fraction(0)), _text(Fraction(1)), _text(Fraction(2)))
    with pytest.raises(native.SeparationBoundExceededError, match="root-separation bound"):
        native.contract_probe_rational_between(*lo, *lo)


# ---------------------------------------------------------------------------
# q17: vertex disks fold zero radicands and zero coefficients
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("parts", "extended"),
    [
        ((Fraction(1), Fraction(3), Fraction(2), Fraction(5), Fraction(0)), (False, False)),
        ((Fraction(1), Fraction(0), Fraction(2), Fraction(5), Fraction(7)), (False, True)),
        ((Fraction(1), Fraction(4), Fraction(2), Fraction(0), Fraction(7)), (True, False)),
        ((Fraction(1), Fraction(1), Fraction(2), Fraction(1), Fraction(4)), (False, False)),
    ],
)
def test_vertex_disk_never_extends_over_a_zero_radicand_or_coefficient(parts: tuple[Fraction, ...], extended: tuple[bool, bool]) -> None:
    assert native.contract_probe_vertex_disk(*map(_text, parts)) == extended
