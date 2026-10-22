"""Exact ownership on one common scale, and the Ω* re-phase that reuses the first exact sort.

Both are speed changes that must decide exactly what the code before them
decided: ownership against the pre-change dyadic implementation (verbatim
below), the restricted order against ``sorted`` with the same comparator, and
the subset re-phase against :func:`benchmarks.held_omega_path.phase_from_start`
on real Figure 5 Ω* circles.
"""

from __future__ import annotations

import functools
import struct
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import benchmarks.held_omega_path as omega_path
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_figure5_path import _common_scale_squared_distances
from benchmarks.held_figure5_path import _cyclic_boundary_neighborhoods
from benchmarks.held_figure5_path import _dyadic_compare
from benchmarks.held_figure5_path import _dyadic_minimum
from benchmarks.held_figure5_path import _dyadic_squared_distance
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_omega_path import AmbiguousOmegaStartError
from benchmarks.held_omega_path import _junction_bridged_indices
from benchmarks.held_omega_path import _owner_indices
from benchmarks.held_omega_path import _phase_subset_from_start
from benchmarks.held_omega_path import _restricted_order
from benchmarks.held_omega_path import phase_from_start
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

_PIECES = 6
_LATTICE = st.integers(min_value=-8, max_value=8).map(lambda value: value / 4)
_SPECIAL = st.sampled_from((0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1e-300, -1e300, 1.7976931348623157e308))
_SUBNORMAL = st.integers(min_value=1, max_value=(1 << 52) - 1).map(lambda bits: struct.unpack("<d", struct.pack("<Q", bits))[0])
_COORDINATE = st.one_of(_LATTICE, _SPECIAL, _SUBNORMAL, st.floats(allow_nan=False, allow_infinity=False))
_POINT = st.builds(lambda x, y: Point2[WorldXY].build(x, y), _COORDINATE, _COORDINATE)
_KINDS = st.lists(st.sampled_from(tuple(PieceKind)), min_size=_PIECES, max_size=_PIECES).filter(lambda kinds: any(kind is not PieceKind.JUNCTION_ARC for kind in kinds))


def _exact(first: Point2[WorldXY], second: Point2[WorldXY]) -> Fraction:
    delta_x = Fraction(float(second.x)) - Fraction(float(first.x))
    delta_y = Fraction(float(second.y)) - Fraction(float(first.y))
    return delta_x * delta_x + delta_y * delta_y


def _old_owner_indices(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    pieces: tuple[int, ...],
    piece_kinds: tuple[PieceKind, ...],
    published_start_center: Point2[WorldXY],
) -> tuple[int, ...]:
    """The pre-change ``held_omega_path._owner_indices``, verbatim."""
    bridged, piece_count = _junction_bridged_indices(piece_kinds)
    pieces = tuple(bridged[piece] for piece in pieces)
    by_station: dict[tuple[GuideRunId, GuideRunStationOrdinal], list[int]] = {}
    for index, hypothesis in enumerate(hypotheses):
        by_station.setdefault((hypothesis.run_id, hypothesis.station_ordinal), []).append(index)
    owners: list[int] = []
    for station_indices in by_station.values():
        for neighborhood in _cyclic_boundary_neighborhoods({pieces[index] for index in station_indices}, piece_count):
            candidates = tuple(index for index in station_indices if pieces[index] in neighborhood)
            distances = tuple(_dyadic_squared_distance(hypotheses[index].middle_point, hypotheses[index].boundary_footpoint) for index in candidates)
            minimum = _dyadic_minimum(distances)
            owners.extend(index for index, distance in zip(candidates, distances, strict=True) if _dyadic_compare(distance, minimum) == 0)
    start_distances = tuple(_dyadic_squared_distance(published_start_center, hypothesis.contact_point) for hypothesis in hypotheses)
    start_minimum = _dyadic_minimum(start_distances)
    for index, distance in enumerate(start_distances):
        if _dyadic_compare(distance, start_minimum) == 0 and hypotheses[index] not in (hypotheses[owner] for owner in owners):
            owners.append(index)
    return tuple(owners)


@settings(max_examples=1500, deadline=None)
@given(pairs=st.lists(st.tuples(_POINT, _POINT), max_size=6))
def test_common_scale_numerators_are_the_exact_squared_distances(pairs: list[tuple[Point2[WorldXY], Point2[WorldXY]]]) -> None:
    numerators, exponent = _common_scale_squared_distances(pairs)
    assert len(numerators) == len(pairs)
    for numerator, (first, second) in zip(numerators, pairs, strict=True):
        assert Fraction(numerator, 1 << (2 * exponent)) == _exact(first, second)


@st.composite
def _hypothesis(draw: st.DrawFn) -> ProjectionAdmissibleBoundaryHypothesis:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(1, 2), _PIECES)
    contact = draw(_POINT)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(draw(st.integers(min_value=0, max_value=1))),
        GuideRunStationOrdinal(draw(st.integers(min_value=0, max_value=2))),
        site,
        (site,),
        Millimetre(0.0),
        draw(_POINT),
        draw(_POINT),
        contact,
        contact,
        GuideRadius.build(1.0),
    )


@settings(max_examples=1500, deadline=None)
@given(
    pool=st.lists(st.tuples(_hypothesis(), st.integers(min_value=0, max_value=_PIECES - 1)), min_size=1, max_size=8),
    picks=st.lists(st.integers(min_value=0, max_value=7), min_size=1, max_size=24),
    kinds=_KINDS,
    start=_POINT,
)
def test_omega_owners_are_the_dyadic_owners(
    pool: list[tuple[ProjectionAdmissibleBoundaryHypothesis, int]],
    picks: list[int],
    kinds: list[PieceKind],
    start: Point2[WorldXY],
) -> None:
    # Re-picking from a small pool makes equal duplicates and exact distance ties.
    chosen = [pool[pick % len(pool)] for pick in picks]
    hypotheses = tuple(hypothesis for hypothesis, _ in chosen)
    pieces = tuple(piece for _, piece in chosen)
    assert _owner_indices(hypotheses, pieces, tuple(kinds), start) == _old_owner_indices(hypotheses, pieces, tuple(kinds), start)


def _coarse(first: tuple[int, int], second: tuple[int, int]) -> int:
    """A strict weak order with many ties: only the first field decides."""
    return (first[0] > second[0]) - (first[0] < second[0])


@settings(max_examples=2000, deadline=None)
@given(data=st.data(), full=st.lists(st.tuples(st.integers(min_value=0, max_value=4), st.integers()), max_size=30))
def test_restricted_order_is_sorted_with_the_same_comparator(data: st.DataObject, full: list[tuple[int, int]]) -> None:
    items = [list(item) for item in full]  # distinct objects even for equal values
    key = functools.cmp_to_key(lambda first, second: _coarse(tuple(first), tuple(second)))
    ordered = sorted(items, key=key)
    position = {id(item): index for index, item in enumerate(ordered)}
    subset = data.draw(st.permutations(items)) if items else []
    subset = subset[: data.draw(st.integers(min_value=0, max_value=len(subset)))]
    expected = sorted(subset, key=key)
    restricted = _restricted_order(tuple(subset), position, lambda first, second: _coarse(tuple(first), tuple(second)))
    assert len(restricted) == len(expected)
    assert all(first is second for first, second in zip(restricted, expected, strict=True))


@pytest.fixture(scope="module")
def figure5_omega_circles() -> tuple[tuple[object, ...], Point2[WorldXY], Millimetre]:
    """The placed Ω* circles of Figure 5 that the first exact sort orders, with the start evidence."""
    captured: list[tuple[object, ...]] = []
    original = omega_path._start_order

    def record(circles, start):  # type: ignore[no-untyped-def]
        captured.append((circles, start))
        return original(circles, start)

    omega_path._start_order = record
    try:
        from benchmarks.held_omega_toolpath import generate_omega

        generate_omega("figure5_refined")
    finally:
        omega_path._start_order = original
    circles, start = captured[0]
    return circles, start, Millimetre(1.0e6)


@settings(max_examples=40, deadline=None)
@given(data=st.data())
def test_subset_phase_is_phase_from_start(data: st.DataObject, figure5_omega_circles: tuple[tuple[object, ...], Point2[WorldXY], Millimetre]) -> None:
    circles, start, bound = figure5_omega_circles
    order = omega_path._start_order(circles, start)
    indices = data.draw(st.lists(st.integers(min_value=0, max_value=len(circles) - 1), min_size=1, max_size=200, unique=True))
    subset = tuple(circles[index] for index in indices)
    expected = phase_from_start(subset, start, bound)
    reused = _phase_subset_from_start(order, subset, bound)
    assert len(reused) == len(expected)
    assert all(first is second for first, second in zip(reused, expected, strict=True))


def test_subset_phase_refuses_as_phase_from_start(figure5_omega_circles: tuple[tuple[object, ...], Point2[WorldXY], Millimetre]) -> None:
    circles, start, _ = figure5_omega_circles
    order = omega_path._start_order(circles, start)
    far = Millimetre(0.0)
    subset = circles[len(circles) // 2 :]
    with pytest.raises(AmbiguousOmegaStartError):
        phase_from_start(subset, start, far)
    with pytest.raises(AmbiguousOmegaStartError):
        _phase_subset_from_start(order, subset, far)
