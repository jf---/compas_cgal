"""Held path bookkeeping decides exactly what the Fraction implementation decided.

The oracles below are the pre-change ``benchmarks/held_figure5_path.py``
functions, verbatim. Every property compares outputs with ``==`` and object
identity, so they pin bit-identical behaviour, not approximate agreement.
"""

from __future__ import annotations

import math
import struct
from collections.abc import Callable
from fractions import Fraction

from hypothesis import event
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import tests.benchmarks.held_figure5_path_pre_lever_oracle as oracle
from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_path import _CanonicalHypothesis
from benchmarks.held_figure5_path import _boundary_neighborhood_owners
from benchmarks.held_figure5_path import _canonicalize
from benchmarks.held_figure5_path import _circle_geometry_key
from benchmarks.held_figure5_path import _cyclic_boundary_neighborhoods
from benchmarks.held_figure5_path import _dyadic_squared_distance
from benchmarks.held_figure5_path import build_hypothesis_figure5_path
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal.motion import EngagementCap
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

_SIDES = 6

# Lattice values force exact ties; signed zeros, subnormals and the extremes of
# the finite range exercise the dyadic exponent bookkeeping.
_LATTICE = st.integers(min_value=-8, max_value=8).map(lambda value: value / 4)
_SPECIAL = st.sampled_from((0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1e-300, -1e300, 1.7976931348623157e308))
_SUBNORMAL = st.integers(min_value=1, max_value=(1 << 52) - 1).map(lambda bits: struct.unpack("<d", struct.pack("<Q", bits))[0])
_GENERIC = st.floats(allow_nan=False, allow_infinity=False)
_COORDINATE = st.one_of(_LATTICE, _SPECIAL, _SUBNORMAL, _GENERIC)


def _point(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


# --- oracles: the pre-change implementations, verbatim ----------------------


def _old_exact_squared_distance(
    first: Point2[WorldXY],
    second: Point2[WorldXY],
) -> Fraction:
    delta_x = Fraction.from_float(float(second.x)) - Fraction.from_float(float(first.x))
    delta_y = Fraction.from_float(float(second.y)) - Fraction.from_float(float(first.y))
    return delta_x * delta_x + delta_y * delta_y


def _old_boundary_neighborhood_owners(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    segment_count: int,
    published_start_center: Point2[WorldXY],
) -> tuple[ProjectionAdmissibleBoundaryHypothesis, ...]:
    by_station: dict[
        tuple[GuideRunId, GuideRunStationOrdinal],
        list[ProjectionAdmissibleBoundaryHypothesis],
    ] = {}
    for hypothesis in hypotheses:
        by_station.setdefault(
            (hypothesis.run_id, hypothesis.station_ordinal),
            [],
        ).append(hypothesis)

    owners: list[ProjectionAdmissibleBoundaryHypothesis] = []
    for station_hypotheses in by_station.values():
        neighborhoods = _cyclic_boundary_neighborhoods(
            {int(hypothesis.projected_boundary_site.segment_id) for hypothesis in station_hypotheses},
            segment_count,
        )
        for neighborhood in neighborhoods:
            candidates = tuple(hypothesis for hypothesis in station_hypotheses if int(hypothesis.projected_boundary_site.segment_id) in neighborhood)
            # Parallel tuples, not hypothesis-keyed dicts: equal hypotheses carry
            # equal distances, so the dict only cost a dataclass hash per entry.
            distances = tuple(
                _old_exact_squared_distance(
                    hypothesis.middle_point,
                    hypothesis.boundary_footpoint,
                )
                for hypothesis in candidates
            )
            minimum = min(distances)
            owners.extend(hypothesis for hypothesis, distance in zip(candidates, distances, strict=True) if distance == minimum)

    start_distances = tuple(
        _old_exact_squared_distance(
            published_start_center,
            hypothesis.contact_point,
        )
        for hypothesis in hypotheses
    )
    start_minimum = min(start_distances)
    # Only exact start-distance ties reach the membership scan (one or two per
    # workload), so a list scan beats hashing every owner into a set.
    for hypothesis, distance in zip(hypotheses, start_distances, strict=True):
        if distance == start_minimum and hypothesis not in owners:
            owners.append(hypothesis)
    return tuple(owners)


def _old_circle_geometry_key(
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
    return (
        Fraction.from_float(float(hypothesis.center.x)),
        Fraction.from_float(float(hypothesis.center.y)),
        Fraction.from_float(float(hypothesis.guide_radius.value)),
        Fraction.from_float(float(hypothesis.contact_point.x)),
        Fraction.from_float(float(hypothesis.contact_point.y)),
    )


def _old_hypothesis_order_key(
    hypothesis: ProjectionAdmissibleBoundaryHypothesis,
) -> tuple[tuple[Fraction, Fraction, Fraction, Fraction, Fraction], int, Fraction, int, int]:
    return (
        _old_circle_geometry_key(hypothesis),
        int(hypothesis.projected_boundary_site.segment_id),
        hypothesis.projected_boundary_site.parameter,
        int(hypothesis.run_id),
        int(hypothesis.station_ordinal),
    )


def _old_canonicalize(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
) -> tuple[_CanonicalHypothesis, ...]:
    grouped: dict[
        tuple[Fraction, Fraction, Fraction, Fraction, Fraction],
        list[ProjectionAdmissibleBoundaryHypothesis],
    ] = {}
    for hypothesis in hypotheses:
        grouped.setdefault(_old_circle_geometry_key(hypothesis), []).append(hypothesis)
    return tuple(
        _CanonicalHypothesis(
            min(group, key=_old_hypothesis_order_key),
            tuple(sorted({candidate.run_id for candidate in group}, key=int)),
        )
        for _, group in sorted(grouped.items())
    )


# --- P1: ownership by exact dyadic integers ----------------------------------


@settings(max_examples=2000, deadline=None)
@given(first=st.tuples(_COORDINATE, _COORDINATE), second=st.tuples(_COORDINATE, _COORDINATE))
def test_dyadic_squared_distance_is_the_exact_fraction(first: tuple[float, float], second: tuple[float, float]) -> None:
    numerator, exponent = _dyadic_squared_distance(_point(*first), _point(*second))
    assert Fraction(numerator, 1 << (2 * exponent)) == _old_exact_squared_distance(_point(*first), _point(*second))


@st.composite
def _ownership_hypothesis(draw: st.DrawFn) -> ProjectionAdmissibleBoundaryHypothesis:
    side = draw(st.integers(min_value=0, max_value=_SIDES - 1))
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), draw(st.sampled_from((Fraction(1, 4), Fraction(1, 2)))), _SIDES)
    contact = _point(draw(_COORDINATE), draw(_COORDINATE))
    middle = _point(draw(_COORDINATE), draw(_COORDINATE))
    footpoint = _point(draw(_COORDINATE), draw(_COORDINATE))
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(draw(st.integers(min_value=0, max_value=2))),
        GuideRunStationOrdinal(draw(st.integers(min_value=0, max_value=2))),
        site,
        (site,),
        Millimetre(0.0),
        middle,
        footpoint,
        contact,
        contact,
        GuideRadius.build(1.0),
    )


@settings(max_examples=1500, deadline=None)
@given(
    pool=st.lists(_ownership_hypothesis(), min_size=1, max_size=8),
    picks=st.lists(st.integers(min_value=0, max_value=7), min_size=1, max_size=24),
    start=st.tuples(_COORDINATE, _COORDINATE),
)
def test_dyadic_owners_are_the_fraction_owners_object_for_object(
    pool: list[ProjectionAdmissibleBoundaryHypothesis],
    picks: list[int],
    start: tuple[float, float],
) -> None:
    # Re-picking from a small pool makes equal and identical duplicates.
    hypotheses = tuple(pool[pick % len(pool)] for pick in picks)
    center = _point(*start)
    new = _boundary_neighborhood_owners(hypotheses, _SIDES, center)
    old = _old_boundary_neighborhood_owners(hypotheses, _SIDES, center)
    assert len(new) == len(old)
    assert all(first is second for first, second in zip(new, old, strict=True))


# --- P2: canonical keys on the raw doubles -----------------------------------

_POSITIVE = st.one_of(st.sampled_from((0.25, 0.5, 1.0, 5e-324, 2.2250738585072014e-308)), _SUBNORMAL, st.floats(min_value=5e-324, allow_infinity=False))
# A small value pool makes equal geometry (grouping) and ±0.0 aliasing common.
_KEY_COORDINATE = st.one_of(st.sampled_from((0.0, -0.0, 0.5, -0.5, 5e-324, -5e-324)), _COORDINATE)


@st.composite
def _canonical_hypothesis(draw: st.DrawFn) -> ProjectionAdmissibleBoundaryHypothesis:
    side = draw(st.integers(min_value=0, max_value=_SIDES - 1))
    site = ProjectionBoundarySite.build(
        ProjectionBoundarySegmentId(side),
        draw(st.sampled_from((Fraction(0), Fraction(1, 4), Fraction(1, 3), Fraction(1, 2)))),
        _SIDES,
    )
    contact = _point(draw(_KEY_COORDINATE), draw(_KEY_COORDINATE))
    center = _point(draw(_KEY_COORDINATE), draw(_KEY_COORDINATE))
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(draw(st.integers(min_value=0, max_value=3))),
        GuideRunStationOrdinal(draw(st.integers(min_value=0, max_value=3))),
        site,
        (site,),
        Millimetre(0.0),
        center,
        center,
        contact,
        center,
        GuideRadius.build(draw(_POSITIVE)),
    )


@settings(max_examples=1500, deadline=None)
@given(first=_canonical_hypothesis(), second=_canonical_hypothesis())
def test_double_geometry_key_orders_equates_and_hashes_like_the_fraction_key(
    first: ProjectionAdmissibleBoundaryHypothesis,
    second: ProjectionAdmissibleBoundaryHypothesis,
) -> None:
    new_first, new_second = _circle_geometry_key(first), _circle_geometry_key(second)
    old_first, old_second = _old_circle_geometry_key(first), _old_circle_geometry_key(second)
    assert all(type(value) is float for value in new_first)
    assert (new_first == new_second) == (old_first == old_second)
    assert (new_first < new_second) == (old_first < old_second)
    if new_first == new_second:
        assert hash(new_first) == hash(new_second)


@settings(max_examples=1500, deadline=None)
@given(
    pool=st.lists(_canonical_hypothesis(), min_size=1, max_size=8),
    picks=st.lists(st.integers(min_value=0, max_value=7), min_size=1, max_size=24),
)
def test_double_keyed_canonicalization_is_the_fraction_one_object_for_object(
    pool: list[ProjectionAdmissibleBoundaryHypothesis],
    picks: list[int],
) -> None:
    hypotheses = tuple(pool[pick % len(pool)] for pick in picks)
    new = _canonicalize(hypotheses)
    old = _old_canonicalize(hypotheses)
    assert new == old
    assert all(first.hypothesis is second.hypothesis for first, second in zip(new, old, strict=True))


# --- P6: offset sites as a tuple parallel to the hypotheses ------------------

_SQUARE = (_point(0.0, 0.0), _point(4.0, 0.0), _point(4.0, 4.0), _point(0.0, 4.0))
_INWARD = ((0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0))
_TOOL_RADIUS = ToolRadius.build(1.0)
_CAP = EngagementCap.build(math.radians(80.0))


@st.composite
def _square_hypotheses(draw: st.DrawFn) -> tuple[tuple[ProjectionAdmissibleBoundaryHypothesis, ...], tuple[ProjectionBoundarySite, ...]]:
    """Runs walking the square's sides, each with an offset site near its own contact."""
    hypotheses = []
    sites = []
    for run in range(draw(st.integers(min_value=1, max_value=3))):
        side = draw(st.integers(min_value=0, max_value=3))
        ticks = sorted(set(draw(st.lists(st.integers(min_value=1, max_value=63), min_size=1, max_size=10))))
        radius = draw(st.sampled_from((1.0, 1.0, 0.75, 1.25)))
        for ordinal, tick in enumerate(ticks):
            parameter = Fraction(tick, 64)
            start, end = _SQUARE[side], _SQUARE[(side + 1) % 4]
            contact = _point(
                float(start.x) + float(parameter) * (float(end.x) - float(start.x)),
                float(start.y) + float(parameter) * (float(end.y) - float(start.y)),
            )
            center = _point(float(contact.x) + radius * _INWARD[side][0], float(contact.y) + radius * _INWARD[side][1])
            site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), parameter, 4)
            hypotheses.append(
                ProjectionAdmissibleBoundaryHypothesis(
                    GuideRunId(run), GuideRunStationOrdinal(ordinal), site, (site,), Millimetre(0.0), center, contact, contact, center, GuideRadius.build(radius)
                )
            )
            offset_tick = min(64, max(0, tick + draw(st.integers(min_value=-1, max_value=1))))
            sites.append(ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), Fraction(offset_tick, 64), 4))
    order = draw(st.permutations(range(len(hypotheses))))
    return tuple(hypotheses[index] for index in order), tuple(sites[index] for index in order)


def _outcome(call: Callable[[], object]) -> str:
    # Typed path failures are outcomes to compare; any other exception
    # (KeyError, TypeError) is a defect and propagates.
    try:
        result = repr(call())
    except (BenchmarkError, ValueError, RuntimeError) as error:
        event(type(error).__name__)
        return f"{type(error).__name__}: {error}"
    event("path")
    return result


@settings(max_examples=400, deadline=None)
@given(case=_square_hypotheses(), preserve=st.booleans())
def test_parallel_offset_sites_build_the_dict_keyed_path_byte_for_byte(
    case: tuple[tuple[ProjectionAdmissibleBoundaryHypothesis, ...], tuple[ProjectionBoundarySite, ...]],
    preserve: bool,
) -> None:
    hypotheses, sites = case
    common = {
        "inward_components": (_SQUARE,),
        "hypotheses": hypotheses,
        "reached_run_ids": frozenset(hypothesis.run_id for hypothesis in hypotheses),
        "published_start_center": hypotheses[0].contact_point,
        "published_start_radius": Millimetre(1.0),
        "boundary_evidence_bound": Millimetre(0.25),
        "start_evidence_bound": Millimetre(0.25),
        "tool_radius": _TOOL_RADIUS,
        "cap": _CAP,
        "preserve_source_runs": preserve,
    }
    new = _outcome(lambda: build_hypothesis_figure5_path(offset_sites=sites, **common))
    old = _outcome(lambda: oracle.build_hypothesis_figure5_path(offset_sites=dict(zip(hypotheses, sites, strict=True)), **common))
    assert new == old
