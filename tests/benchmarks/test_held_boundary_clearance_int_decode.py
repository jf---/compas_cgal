"""The integer decode of the native clearance records equals the Fraction decode it replaces.

The oracle is the pre-change ``_site`` (HEAD 2445a48a), kept verbatim: it
materialises every value as a reduced ``Fraction``. The integer decode keeps
the native integers and derives each value on demand. Required, record by
record: the same segment and parameter, the same exact footpoint and squared
distance, doubles equal to ``float`` of the oracle's Fractions (both are the
correctly rounded quotient of the same rational), and the exact
``squared_distance > bound`` verdict equal to the Fraction comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.held_boundary_clearance import BoundaryClearanceSite
from compas_cgal.held_boundary_clearance import HeldBoundaryClearanceIndex
from compas_cgal.held_boundary_clearance import _site
from compas_cgal.units import Millimetre
from compas_cgal.units import SquaredMillimetre

Record = tuple[int, str, str, str, str, str, int]


@dataclass(frozen=True)
class _OracleSite:
    """The pre-change ``BoundaryClearanceSite`` fields."""

    segment: int
    parameter: Fraction
    footpoint: tuple[Fraction, Fraction]
    squared_distance: Fraction


def _oracle_scaled(numerator: int, denominator: int, exponent: int) -> Fraction:
    """Verbatim pre-change ``_scaled``."""
    if exponent >= 0:
        return Fraction(numerator << exponent, denominator)
    return Fraction(numerator, denominator << -exponent)


def _oracle_site(record: Record) -> _OracleSite:
    """Verbatim pre-change ``_site``."""
    segment, along, length_squared, foot_x, foot_y, distance, exponent = record
    denominator = int(length_squared, 16)
    return _OracleSite(
        segment,
        Fraction(int(along, 16), denominator),
        (
            _oracle_scaled(int(foot_x, 16), denominator, exponent),
            _oracle_scaled(int(foot_y, 16), denominator, exponent),
        ),
        _oracle_scaled(int(distance, 16), denominator, 2 * exponent),
    )


def _float_outcome(value: Fraction) -> float | type[Exception]:
    try:
        return float(value)
    except OverflowError:
        return OverflowError


def _double_outcome(site: BoundaryClearanceSite, field: str) -> object:
    try:
        return getattr(site, field)
    except OverflowError:
        return OverflowError


def _assert_matches_oracle(record: Record, bounds: tuple[Fraction, ...]) -> None:
    site = _site(record)
    oracle = _oracle_site(record)

    assert site.segment == oracle.segment
    assert site.parameter == oracle.parameter
    assert type(site.parameter) is Fraction
    assert site.footpoint == oracle.footpoint
    assert site.squared_distance == oracle.squared_distance
    footpoint_outcomes = tuple(_float_outcome(value) for value in oracle.footpoint)
    expected_footpoint = OverflowError if OverflowError in footpoint_outcomes else footpoint_outcomes
    assert _double_outcome(site, "footpoint_doubles") == expected_footpoint, record
    assert _double_outcome(site, "squared_distance_double") == _float_outcome(oracle.squared_distance), record
    for bound in bounds:
        assert site.squared_distance_exceeds(SquaredMillimetre(bound)) == (oracle.squared_distance > bound), (record, bound)


def _boundary_bounds(value: Fraction, extra: tuple[Fraction, ...] = ()) -> tuple[Fraction, ...]:
    """The value itself and its tightest neighbours: the verdict flips exactly there."""
    nudge = Fraction(1, value.denominator * 2**64)
    return (value, value - nudge, value + nudge, *extra)


def test_int_decode_matches_fraction_decode_on_every_figure5_record() -> None:
    case = load_held_reference_case("figure5")
    guide = build_held_reference_raw_guide(case)
    index = HeldBoundaryClearanceIndex.build(case.projection.points)
    stations = tuple(station.middle_point for run in guide.runs for station in run.stations)
    gap = Millimetre(float(guide.site_budget.admissible_distance_gap))
    tool_radius = Fraction.from_float(float(case.tool_radius.value))
    rows = index._native.distance_admissible_sites([(float(point.x), float(point.y)) for point in stations], float(gap))

    checked = 0
    for row in rows:
        for record in row:
            _assert_matches_oracle(record, (tool_radius * tool_radius,))
            checked += 1

    # LB's decode probe counted 31,742 Figure 5 records; at least one per station.
    assert checked > guide.station_count


def test_wrapper_batch_returns_the_int_decoded_sites() -> None:
    case = load_held_reference_case("figure5")
    index = HeldBoundaryClearanceIndex.build(case.projection.points)
    station = case.projection.points[0]

    (row,) = index.distance_admissible_sites((station,), Millimetre(0.0))

    assert row
    assert all(type(site) is BoundaryClearanceSite for site in row)
    assert all(site.squared_distance == 0 for site in row)


_HEX_MAGNITUDE = st.one_of(st.integers(min_value=0, max_value=2**16), st.integers(min_value=0, max_value=2**320))
_SIGNED = st.one_of(_HEX_MAGNITUDE, _HEX_MAGNITUDE.map(lambda value: -value))
# Native exponents put dyadic inputs on a common integer grid; the extremes
# overflow or underflow binary64, which both decodes must report identically.
_EXPONENT = st.one_of(st.integers(min_value=-80, max_value=80), st.integers(min_value=-1200, max_value=1200))


def _hex(value: int) -> str:
    return f"-{-value:x}" if value < 0 else f"{value:x}"


@st.composite
def _records(draw: st.DrawFn) -> Record:
    length_squared = draw(st.integers(min_value=1, max_value=2**320))
    along = draw(st.integers(min_value=0, max_value=length_squared - 1))
    return (
        draw(st.integers(min_value=0, max_value=10_000)),
        _hex(along),
        _hex(length_squared),
        _hex(draw(_SIGNED)),
        _hex(draw(_SIGNED)),
        _hex(draw(_HEX_MAGNITUDE)),
        draw(_EXPONENT),
    )


@settings(max_examples=1500, deadline=None)
@given(record=_records(), bound=st.fractions(min_value=0))
@example(record=(0, "0", "1", "0", "0", "0", 0), bound=Fraction(0))
@example(record=(3, "1", "3", "-5", "7", "9", -3), bound=Fraction(9, 3 * 2**6))
@example(record=(1, "0", "1", "1", "1", "1", 1100), bound=Fraction(1))
@example(record=(1, "0", "3", "1", "-1", "1", -1100), bound=Fraction(0))
def test_int_decode_matches_fraction_decode_on_random_records(record: Record, bound: Fraction) -> None:
    oracle = _oracle_site(record).squared_distance
    _assert_matches_oracle(record, _boundary_bounds(oracle, (bound,)))


@settings(max_examples=300, deadline=None)
@given(record=_records(), factor=st.integers(min_value=1, max_value=2**64))
def test_equality_is_exact_value_equality(record: Record, factor: int) -> None:
    # The same rational on a scaled integer grid is the same site.
    segment, along, length_squared, foot_x, foot_y, distance, exponent = record
    scaled = (
        segment,
        _hex(int(along, 16) * factor),
        _hex(int(length_squared, 16) * factor),
        _hex(int(foot_x, 16) * factor),
        _hex(int(foot_y, 16) * factor),
        _hex(int(distance, 16) * factor),
        exponent,
    )

    assert _site(scaled) == _site(record)
    assert hash(_site(scaled)) == hash(_site(record))
    moved = (segment, along, length_squared, foot_x, foot_y, _hex(int(distance, 16) + 1), exponent)
    assert _site(moved) != _site(record)
