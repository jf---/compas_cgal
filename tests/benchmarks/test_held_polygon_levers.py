"""Polygon-generation levers P1, P2, P5, P6: identity of the fast forms with the forms they replace."""

from __future__ import annotations

from collections.abc import Callable
from fractions import Fraction

import pytest
from compas.geometry import Point
from compas.geometry import Polygon
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.errors import BenchmarkError
from benchmarks.errors import DegeneratePocketError
from benchmarks.errors import InvalidCapError
from benchmarks.errors import NonPositiveToolError
from benchmarks.errors import PocketNotSimpleError
from benchmarks.held_common_scale import common_scale_squared_distances
from benchmarks.held_figure5_boundary_path import MisalignedFigure5OffsetSitesError
from benchmarks.held_figure5_boundary_path import _squared_distance
from benchmarks.held_figure5_boundary_path import build_figure5_boundary_path
from benchmarks.held_figure5_boundary_path import build_indexed_figure5_boundary_path
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import _native_guide_rows
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.spec import PocketSpec
from compas_cgal.toolpath import _polygon_to_ccw_vertices
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

SQUARE = (Point2.build(0.0, 0.0), Point2.build(4.0, 0.0), Point2.build(4.0, 4.0), Point2.build(0.0, 4.0))
TOOL_RADIUS = ToolRadius.build(1.0)
# Lattice of 1/64 mm along each side, so exact ties of the start distance occur.
TICKS = 64
# Self-intersecting (edge 2 crosses edge 0 at (4.8, 0)) with area 8, so only simplicity rejects it.
BOWTIE = Polygon([[0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [8.0, 8.0, 0.0], [4.0, -2.0, 0.0]])


def _candidate(run: int, side: int, tick: int, radius: float) -> ProjectionAdmissibleBoundaryHypothesis:
    start, end = SQUARE[side], SQUARE[(side + 1) % 4]
    fraction = Fraction(tick, TICKS)
    contact: Point2[WorldXY] = Point2.build(
        float(start.x) + float(fraction) * (float(end.x) - float(start.x)),
        float(start.y) + float(fraction) * (float(end.y) - float(start.y)),
    )
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), fraction, 4)
    centre: Point2[WorldXY] = Point2.build(2.0, 2.0)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(run), GuideRunStationOrdinal(0), site, (site,), Millimetre(0.0), centre, contact, contact, centre, GuideRadius.build(radius)
    )


@st.composite
def _candidates(draw: st.DrawFn) -> tuple[tuple[ProjectionAdmissibleBoundaryHypothesis, ...], tuple[ProjectionBoundarySite, ...], Point2[WorldXY]]:
    count = draw(st.integers(min_value=1, max_value=8))
    candidates = []
    sites = []
    for run in range(count):
        side = draw(st.integers(min_value=0, max_value=3))
        tick = draw(st.integers(min_value=1, max_value=TICKS - 1))
        candidates.append(_candidate(run, side, tick, draw(st.sampled_from((0.5, 0.75, 1.0)))))
        offset_tick = draw(st.integers(min_value=1, max_value=TICKS - 1))
        sites.append(ProjectionBoundarySite.build(ProjectionBoundarySegmentId(side), Fraction(offset_tick, TICKS), 4))
    start: Point2[WorldXY] = Point2.build(draw(st.sampled_from((0.0, 2.0, 4.0))), draw(st.sampled_from((0.0, 2.0, 4.0))))
    return tuple(candidates), tuple(sites), start


def _common(candidates: tuple[ProjectionAdmissibleBoundaryHypothesis, ...], start: Point2[WorldXY], bound: float) -> dict[str, object]:
    return {
        "inward_components": (SQUARE,),
        "side_candidates": candidates,
        "reached_run_ids": frozenset(candidate.run_id for candidate in candidates),
        "published_start_center": start,
        "published_start_radius": Millimetre(1.0),
        "tool_radius": TOOL_RADIUS,
        "boundary_evidence_bound": Millimetre(0.25),
        "start_evidence_bound": Millimetre(bound),
    }


def _outcome(call: Callable[[], object]) -> str:
    # Typed path failures are outcomes to compare; any other exception is a defect and propagates.
    try:
        return repr(call())
    except (BenchmarkError, ValueError) as error:
        return f"{type(error).__name__}: {error}"


@settings(max_examples=300, deadline=None)
@given(case=_candidates(), bound=st.sampled_from((0.25, 2.0, 2.0 * 2.0**0.5, 8.0)))
def test_indexed_parallel_path_is_the_dict_keyed_path_and_names_its_sources(case, bound) -> None:  # type: ignore[no-untyped-def]
    """P1 + P5: parallel sites give the dict-keyed path; each circle names the candidate it came from."""
    candidates, sites, start = case
    common = _common(candidates, start, bound)
    keyed = _outcome(lambda: build_figure5_boundary_path(offset_sites=dict(zip(candidates, sites)), **common))
    indexed = _outcome(lambda: build_indexed_figure5_boundary_path(offset_sites=sites, **common).path)
    assert indexed == keyed
    if keyed.startswith("Figure5BoundaryPath("):
        result = build_indexed_figure5_boundary_path(offset_sites=sites, **common)
        assert sorted(result.source_indices) == list(range(len(candidates)))
        for circle, index in zip(result.path.circles, result.source_indices):
            assert circle.run_id == candidates[index].run_id
            assert circle.boundary_site == sites[index]


def test_indexed_path_refuses_misaligned_offset_sites() -> None:
    candidates = (_candidate(0, 0, 16, 1.0), _candidate(1, 1, 16, 1.0))
    with pytest.raises(MisalignedFigure5OffsetSitesError):
        build_indexed_figure5_boundary_path(offset_sites=(candidates[0].projected_boundary_site,), **_common(candidates, SQUARE[0], 8.0))


finite = st.floats(allow_nan=False, allow_infinity=False, min_value=-1e6, max_value=1e6)
lattice = st.integers(min_value=-256, max_value=256).map(lambda n: n / 64.0)


@settings(max_examples=1500, deadline=None)
@given(
    pairs=st.lists(st.tuples(st.one_of(finite, lattice), st.one_of(finite, lattice), st.one_of(finite, lattice), st.one_of(finite, lattice)), min_size=1, max_size=6),
    bound=st.one_of(st.floats(min_value=2.0**-30, max_value=1e6, allow_nan=False), st.integers(min_value=1, max_value=512).map(lambda n: n / 64.0)),
)
def test_common_scale_start_seed_decides_as_the_fraction_seed(pairs, bound) -> None:  # type: ignore[no-untyped-def]
    """P2: the integer argmin set and the bound verdict equal the ``Fraction`` ones, exact ties included."""
    points = [(Point2.build(ax, ay), Point2.build(bx, by)) for ax, ay, bx, by in pairs]
    fractions = [_squared_distance(first, second) for first, second in points]
    scaled = common_scale_squared_distances(points)
    minimum = min(scaled.numerators)
    assert [i for i, n in enumerate(scaled.numerators) if n == minimum] == [i for i, f in enumerate(fractions) if f == min(fractions)]
    exact_bound = Fraction(bound)
    assert scaled.exceeds_squared(minimum, bound) == (min(fractions) > exact_bound * exact_bound)


def test_exceeds_squared_holds_exactly_at_the_bound() -> None:
    first: Point2[WorldXY] = Point2.build(0.0, 0.0)
    second: Point2[WorldXY] = Point2.build(0.1, 0.0)
    scaled = common_scale_squared_distances([(first, second)])
    assert scaled.exceeds_squared(scaled.numerators[0], 0.1) is False
    assert scaled.exceeds_squared(scaled.numerators[0], 0.09999999999999999) is True


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_guide_pocket_spec_is_the_validated_spec(name: str) -> None:
    """P6: the guide's spec carries the same fields as the fully validated one."""
    case = load_held_reference_case(name)
    assert repr(case.guide_pocket_spec()) == repr(case.pocket_spec())


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_native_guide_sees_exactly_the_spec_ring(name: str) -> None:
    """P6: the native exact simplicity check is decided on the spec's own doubles (no ring point dropped or moved)."""
    spec = load_held_reference_case(name).guide_pocket_spec()
    assert spec.holes == ()
    vertices = _polygon_to_ccw_vertices(spec.polygon)
    assert [(float(x), float(y)) for x, y, _ in vertices.tolist()] == [(float(p[0]), float(p[1])) for p in spec.polygon.points]


def test_native_guide_refuses_a_non_simple_ring() -> None:
    """P6: the exact native check is live on the guide request (CGAL ``Polygon_2::is_simple``)."""
    spec = PocketSpec.build_natively_simple("bowtie", "test", BOWTIE, 1.0, 80.0)
    with pytest.raises(ValueError, match="simple"):
        _native_guide_rows(spec)


def test_natively_simple_factory_defers_only_simplicity() -> None:
    with pytest.raises(PocketNotSimpleError):
        PocketSpec.build("bowtie", "test", BOWTIE, 1.0, 80.0)
    assert PocketSpec.build_natively_simple("bowtie", "test", BOWTIE, 1.0, 80.0).polygon is BOWTIE
    square = Polygon([Point(0.0, 0.0, 0.0), Point(4.0, 0.0, 0.0), Point(4.0, 4.0, 0.0), Point(0.0, 4.0, 0.0)])
    with pytest.raises(NonPositiveToolError):
        PocketSpec.build_natively_simple("s", "test", square, 0.0, 80.0)
    with pytest.raises(InvalidCapError):
        PocketSpec.build_natively_simple("s", "test", square, 1.0, 181.0)
    with pytest.raises(DegeneratePocketError):
        PocketSpec.build_natively_simple("s", "test", square, 5.0, 80.0)
