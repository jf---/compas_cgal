"""Ownership pins at one unit of the common scale, and both ownership oracles rerun on a lattice.

The common-scale ownership decides every squared distance as an integer
numerator over one shared ``4**exponent``. A one-unit slack (``<= minimum + 1``
instead of ``== minimum``) on either minimum is invisible to the oracle tests
that draw arbitrary doubles: one non-lattice coordinate pushes the shared
exponent up to 1074, so two distances one unit apart are practically never
drawn. The pins below place the competing distances exactly one unit apart
(coordinates on the quarter lattice, so the unit is 1/16 mm²), and the slack
mutants are run against them to show they are killed. The lattice reruns draw
every coordinate, the start included, from the quarter lattice, where near
ties are as common as ties.
"""

from __future__ import annotations

from collections.abc import Callable
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import benchmarks.held_figure5_path as figure5_path
import benchmarks.held_omega_path as omega_path
import tests.benchmarks.test_held_order_reuse_and_common_scale as omega_oracle
import tests.benchmarks.test_held_path_exact_bookkeeping as polygon_oracle
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_figure5_path import _boundary_neighborhood_owner_indices
from benchmarks.held_figure5_path import _boundary_neighborhood_owners
from benchmarks.held_figure5_path import _common_scale_squared_distances
from benchmarks.held_figure5_path import _cyclic_boundary_neighborhoods
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_omega_path import _owner_indices
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Ring / piece count of the fixtures (any count >= 1 works; six matches the oracle files).
_SEGMENTS = 6
_PIECE_KINDS = (PieceKind.OFFSET_LINE,) * _SEGMENTS
# One unit of slack on the integer numerators: the smallest step the common scale can express.
_ONE_UNIT = 1
_START = (0.0, 0.0)
_EXPECTED_OWNERS = (0,)
# Quarter lattice in [-1, 1] mm: squared distances are integers over 16, so near ties are frequent.
_LATTICE = st.integers(min_value=-4, max_value=4).map(lambda value: value / 4)
_LATTICE_POINT = st.builds(lambda x, y: Point2[WorldXY].build(x, y), _LATTICE, _LATTICE)

OwnerIndices = Callable[[tuple[ProjectionAdmissibleBoundaryHypothesis, ...], Point2[WorldXY]], tuple[int, ...]]


def _point(xy: tuple[float, float]) -> Point2[WorldXY]:
    return Point2[WorldXY].build(*xy)


def _hypothesis(foot: tuple[float, float], contact: tuple[float, float]) -> ProjectionAdmissibleBoundaryHypothesis:
    """One hypothesis of station (0, 0) on segment 0 with middle point (0, 0)."""
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(1, 2), _SEGMENTS)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(0),
        GuideRunStationOrdinal(0),
        site,
        (site,),
        Millimetre(0.0),
        _point((0.0, 0.0)),
        _point(foot),
        _point(contact),
        _point(contact),
        GuideRadius.build(1.0),
    )


# A owns the neighbourhood (foot distance 1 vs 4); B's contact is one unit
# (1/16 mm²) farther from the start than A's: start numerators 1 and 2 at exponent 2.
_START_NEAR_TIE = (_hypothesis((0.0, 1.0), (0.25, 0.0)), _hypothesis((0.0, 2.0), (0.25, 0.25)))
# B's foot is one unit (1/16 mm²) farther from the middle point than A's: numerators 16 and 17.
_NEIGHBORHOOD_NEAR_TIE = (_hypothesis((0.0, 1.0), (0.25, 0.0)), _hypothesis((0.25, 1.0), (0.25, 0.25)))


def _polygon_owners(hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...], start: Point2[WorldXY]) -> tuple[int, ...]:
    return _boundary_neighborhood_owner_indices(hypotheses, _SEGMENTS, start)


def _omega_owners(hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...], start: Point2[WorldXY]) -> tuple[int, ...]:
    pieces = tuple(int(hypothesis.projected_boundary_site.segment_id) for hypothesis in hypotheses)
    return _owner_indices(hypotheses, pieces, _PIECE_KINDS, start)


_ENTRY_POINTS: dict[str, OwnerIndices] = {"polygon": _polygon_owners, "omega": _omega_owners}
_FIXTURES = {"start_near_tie": _START_NEAR_TIE, "neighborhood_near_tie": _NEIGHBORHOOD_NEAR_TIE}


def _common_scale_numerators(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    first: Callable[[ProjectionAdmissibleBoundaryHypothesis], Point2[WorldXY]],
    second: Callable[[ProjectionAdmissibleBoundaryHypothesis], Point2[WorldXY]],
) -> tuple[int, ...]:
    return _common_scale_squared_distances((first(hypothesis), second(hypothesis)) for hypothesis in hypotheses).numerators


def _slack_start_owners(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    owners: list[int],
    published_start_center: Point2[WorldXY],
) -> tuple[int, ...]:
    """``_with_start_owners`` with one unit of slack on the start minimum (the mutant)."""
    start_distances = _common_scale_numerators(hypotheses, lambda _: published_start_center, lambda hypothesis: hypothesis.contact_point)
    start_minimum = min(start_distances)
    for index, distance in enumerate(start_distances):
        if distance <= start_minimum + _ONE_UNIT and hypotheses[index] not in (hypotheses[owner] for owner in owners):
            owners.append(index)
    return tuple(owners)


def _slack_neighborhood_owners(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    by_station: dict[tuple[GuideRunId, GuideRunStationOrdinal], list[int]],
    segments: tuple[int, ...],
    segment_count: int,
) -> list[int]:
    """``_station_neighborhood_owners`` with one unit of slack on each neighbourhood minimum (the mutant)."""
    distances = _common_scale_numerators(hypotheses, lambda hypothesis: hypothesis.middle_point, lambda hypothesis: hypothesis.boundary_footpoint)
    owners: list[int] = []
    for station_indices in by_station.values():
        for neighborhood in _cyclic_boundary_neighborhoods({segments[index] for index in station_indices}, segment_count):
            candidates = [index for index in station_indices if segments[index] in neighborhood]
            minimum = min(distances[index] for index in candidates)
            owners.extend(index for index in candidates if distances[index] <= minimum + _ONE_UNIT)
    return owners


def test_fixtures_sit_one_unit_apart() -> None:
    start = _point(_START)
    assert _common_scale_numerators(_START_NEAR_TIE, lambda _: start, lambda hypothesis: hypothesis.contact_point) == (1, 2)
    assert _common_scale_numerators(_NEIGHBORHOOD_NEAR_TIE, lambda hypothesis: hypothesis.middle_point, lambda hypothesis: hypothesis.boundary_footpoint) == (16, 17)


@pytest.mark.parametrize("fixture", sorted(_FIXTURES))
@pytest.mark.parametrize("entry", sorted(_ENTRY_POINTS))
def test_one_unit_farther_hypothesis_does_not_own(entry: str, fixture: str) -> None:
    assert _ENTRY_POINTS[entry](_FIXTURES[fixture], _point(_START)) == _EXPECTED_OWNERS


@pytest.mark.parametrize(
    ("fixture", "target", "mutant"),
    [
        ("start_near_tie", "_with_start_owners", _slack_start_owners),
        ("neighborhood_near_tie", "_station_neighborhood_owners", _slack_neighborhood_owners),
    ],
)
@pytest.mark.parametrize("entry", sorted(_ENTRY_POINTS))
def test_one_unit_slack_mutant_is_killed(monkeypatch: pytest.MonkeyPatch, entry: str, fixture: str, target: str, mutant: Callable[..., object]) -> None:
    # Both entry points resolve the helpers through their own module globals.
    monkeypatch.setattr(figure5_path, target, mutant)
    monkeypatch.setattr(omega_path, target, mutant)
    assert _ENTRY_POINTS[entry](_FIXTURES[fixture], _point(_START)) != _EXPECTED_OWNERS


@st.composite
def _lattice_omega_hypothesis(draw: st.DrawFn) -> ProjectionAdmissibleBoundaryHypothesis:
    """``omega_oracle._hypothesis`` with every point on the quarter lattice."""
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(1, 2), omega_oracle._PIECES)
    contact = draw(_LATTICE_POINT)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(draw(st.integers(min_value=0, max_value=1))),
        GuideRunStationOrdinal(draw(st.integers(min_value=0, max_value=2))),
        site,
        (site,),
        Millimetre(0.0),
        draw(_LATTICE_POINT),
        draw(_LATTICE_POINT),
        contact,
        contact,
        GuideRadius.build(1.0),
    )


@settings(max_examples=1500, deadline=None)
@given(
    pool=st.lists(st.tuples(_lattice_omega_hypothesis(), st.integers(min_value=0, max_value=omega_oracle._PIECES - 1)), min_size=1, max_size=8),
    picks=st.lists(st.integers(min_value=0, max_value=7), min_size=1, max_size=24),
    kinds=omega_oracle._KINDS,
    start=_LATTICE_POINT,
)
def test_omega_owners_are_the_dyadic_owners_on_the_lattice(
    pool: list[tuple[ProjectionAdmissibleBoundaryHypothesis, int]],
    picks: list[int],
    kinds: list[PieceKind],
    start: Point2[WorldXY],
) -> None:
    chosen = [pool[pick % len(pool)] for pick in picks]
    hypotheses = tuple(hypothesis for hypothesis, _ in chosen)
    pieces = tuple(piece for _, piece in chosen)
    assert _owner_indices(hypotheses, pieces, tuple(kinds), start) == omega_oracle._old_owner_indices(hypotheses, pieces, tuple(kinds), start)


@st.composite
def _lattice_polygon_hypothesis(draw: st.DrawFn) -> ProjectionAdmissibleBoundaryHypothesis:
    """``polygon_oracle._ownership_hypothesis`` with every point on the quarter lattice."""
    side = draw(st.integers(min_value=0, max_value=polygon_oracle._SIDES - 1))
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), draw(st.sampled_from((Fraction(1, 4), Fraction(1, 2)))), polygon_oracle._SIDES)
    contact = draw(_LATTICE_POINT)
    middle = draw(_LATTICE_POINT)
    footpoint = draw(_LATTICE_POINT)
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
    pool=st.lists(_lattice_polygon_hypothesis(), min_size=1, max_size=8),
    picks=st.lists(st.integers(min_value=0, max_value=7), min_size=1, max_size=24),
    start=_LATTICE_POINT,
)
def test_polygon_owners_are_the_fraction_owners_on_the_lattice(
    pool: list[ProjectionAdmissibleBoundaryHypothesis],
    picks: list[int],
    start: Point2[WorldXY],
) -> None:
    hypotheses = tuple(pool[pick % len(pool)] for pick in picks)
    new = _boundary_neighborhood_owners(hypotheses, polygon_oracle._SIDES, start)
    old = polygon_oracle._old_boundary_neighborhood_owners(hypotheses, polygon_oracle._SIDES, start)
    assert len(new) == len(old)
    assert all(first is second for first, second in zip(new, old, strict=True))
