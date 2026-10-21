"""Contract tests for the declared (Pythagorean) Held boundary."""

from __future__ import annotations

import functools
import math
from dataclasses import replace
from fractions import Fraction
from functools import cache

import pytest

from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.held_declared_boundary import T_GRID_DENOMINATOR
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import DegenerateArcChordError
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import InvalidDeclaredProvenanceError
from benchmarks.held_declared_boundary import InvalidPythagoreanDirectionError
from benchmarks.held_declared_boundary import OpenDeclaredCycleError
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_declared_boundary_certificate import DeclaredBoundaryCertificate
from benchmarks.held_declared_boundary_certificate import DeclaredProvenanceMismatchError
from benchmarks.held_declared_boundary_certificate import OpenLemma
from benchmarks.held_declared_boundary_certificate import arc_move_bound
from benchmarks.held_declared_boundary_certificate import measure_moves
from benchmarks.held_declared_boundary_certificate import sqrt_lower
from benchmarks.held_declared_boundary_certificate import sqrt_upper
from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_declared_boundary_emission import reemit_primitives
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferencePrimitive
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

# The design prototype's measured maximum displacement (arc centre move,
# crossed skis): docs/held_exact_centre_domain.md, "Boundary contract".
MEASURED_DISPLACEMENT_BOUND = Fraction(23, 10**12)
# Exact junction census (exact G1, left, right) of the declared cycle. The
# exact-G1 counts equal the design prototype's
# (build/held-exact-omega/pythagorean_reemission_tgrid56.json). Its left/right
# split differs on upper (37/33), crossed skis (25/33) and Monstera (140/171):
# the prototype reduced the quarter index mod 4 before subtracting it from the
# angle, so negative-quadrant residuals were off by 2π and rounded to other
# grid points (24 of 70 upper directions). The ~1e-16 rad kinks that split
# then flip sign; G1 is not imposed, so either split is a valid cycle.
DECLARED_JUNCTION_CENSUS = {
    "figure5": (1, 15, 15),
    "figure8_upper": (0, 34, 36),
    "figure8_crossed_skis": (0, 24, 34),
    "figure8_monstera": (6, 136, 175),
}


@cache
def _case(name: str) -> HeldReferenceCase:
    return load_held_reference_case(name)


@cache
def _declared(name: str) -> DeclaredBoundary:
    return declared_boundary_from_reference(name, _case(name).boundary)


def _anchor(primitives: tuple[ReferencePrimitive, ...]) -> int:
    return next(index for index in range(len(primitives)) if isinstance(primitives[index - 1], ReferenceArc))


def _nudge_arc_end(primitives: tuple[ReferencePrimitive, ...], index: int) -> tuple[ReferencePrimitive, ...]:
    arc = primitives[index]
    assert isinstance(arc, ReferenceArc)
    nudged = Point2[WorldXY].build(math.nextafter(float(arc.end.x), math.inf), float(arc.end.y))
    return (*primitives[:index], ReferenceArc.build(arc.start, nudged, arc.centre, arc.sweep), *primitives[index + 1 :])


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_declared_boundary_closes_exactly(name: str) -> None:
    declared = _declared(name)
    primitives = declared.primitives
    assert len(primitives) == len(_case(name).boundary.primitives)
    assert [primitive.source_index for primitive in primitives] == list(range(len(primitives)))
    for current, following in zip(primitives, (*primitives[1:], primitives[0])):
        assert current.end == following.start


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_every_direction_centre_and_radius_is_exact_and_rational(name: str) -> None:
    for primitive in _declared(name).primitives:
        start_tangent = primitive.start_tangent
        end_tangent = primitive.end_tangent
        assert start_tangent[0] ** 2 + start_tangent[1] ** 2 == 1
        assert end_tangent[0] ** 2 + end_tangent[1] ** 2 == 1
        if isinstance(primitive, DeclaredLine):
            assert T_GRID_DENOMINATOR % primitive.direction.t.denominator == 0
            assert Fraction(float(primitive.length)) == primitive.length
            continue
        assert T_GRID_DENOMINATOR % primitive.start_tangent_direction.t.denominator == 0
        centre = primitive.centre
        radius = primitive.radius
        assert isinstance(radius, Fraction) and radius > 0
        for point in (primitive.start, primitive.end):
            radial = point.minus(centre)
            assert radial[0] ** 2 + radial[1] ** 2 == primitive.radius_squared
        # (ρ ± r)² is rational because ρ is: the offset circle is a rational circle.
        tool_radius = Fraction(float(_case(name).tool_radius.value))
        assert isinstance((radius + tool_radius) ** 2, Fraction)


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_junction_census_is_exact_and_imposes_no_g1(name: str) -> None:
    census = _declared(name).junction_census()
    assert (census.exact_g1, census.left_turns, census.right_turns) == DECLARED_JUNCTION_CENSUS[name]
    assert census.reversals == 0


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_move_from_stored_is_within_the_measured_displacement(name: str) -> None:
    moves = measure_moves(_case(name).boundary.primitives, _declared(name))
    maximum = max(move.bound for move in moves)
    print(f"{name}: certified max move {float(maximum):.3e} mm = {maximum.numerator.bit_length()}/{maximum.denominator.bit_length()}-bit rational")
    assert 0 < maximum <= MEASURED_DISPLACEMENT_BOUND


@pytest.mark.parametrize("name", CANONICAL_CASE_NAMES)
def test_declared_boundary_is_recertified_against_the_publisher_evidence(name: str) -> None:
    case = _case(name)
    certificate = DeclaredBoundaryCertificate.build(case, _declared(name))
    # Review I2.2: B enters as the next double above the recorded one (API
    # semantic), so a record rounded to nearest below the true bound is covered.
    recorded = float(case.reconstruction.deviation_upper_bound)
    stored_bound = Fraction(math.nextafter(recorded, math.inf))
    assert certificate.stored_bound == stored_bound > Fraction(recorded)
    assert certificate.declared_bound == stored_bound + certificate.maximum_move
    assert certificate.declared_bound <= stored_bound + MEASURED_DISPLACEMENT_BOUND
    assert certificate.declared_bound <= certificate.reconstruction_limit
    assert certificate.boundary == _declared(name).identity


def test_composed_bound_names_the_open_biarc_lemma_in_its_type_and_identity() -> None:
    # Review I2.1: B + δ is a composed bound conditional on an unproved lemma
    # (B is a two-sided point-set bound for biarc spans); the type carries it.
    certificate = DeclaredBoundaryCertificate.build(_case("figure5"), _declared("figure5"))
    assert certificate.conditional_on == (OpenLemma.BIARC_POINT_SET_BOUND,)
    assert replace(certificate, conditional_on=()).identity_sha256 != certificate.identity_sha256


def test_derived_values_never_write_the_frozen_declared_types() -> None:
    # Review I1.1: no cached_property on frozen types; memos are per-instance
    # objects created in __init__, so reads never touch the frozen __dict__.
    declared = _declared("figure5")
    for kind in (PythagoreanDirection, DeclaredLine, DeclaredArc, DeclaredBoundary):
        assert not any(isinstance(value, functools.cached_property) for value in vars(kind).values()), kind
    snapshot = {id(item): dict(vars(item)) for item in (declared, *declared.primitives)}
    identity = declared.identity
    for primitive in declared.primitives:
        if isinstance(primitive, DeclaredArc):
            assert primitive.centre is primitive.centre and primitive.radius is primitive.radius
        else:
            assert primitive.end is primitive.end and primitive.direction.vector is primitive.direction.vector
    assert {id(item): dict(vars(item)) for item in (declared, *declared.primitives)} == snapshot
    assert declared.identity is identity
    twin = replace(declared)
    assert twin == declared and twin._identity is not declared._identity and twin._identity.value is None
    assert twin.identity == identity


def test_identity_is_content_addressed() -> None:
    stored = _case("figure5").boundary.primitives
    first = reemit_primitives("figure5", stored)
    again = reemit_primitives("figure5", stored)
    assert first.identity == again.identity
    corrupted = reemit_primitives("figure5", _nudge_arc_end(stored, next(i for i, p in enumerate(stored) if isinstance(p, ReferenceArc) and i != _anchor(stored) - 1)))
    assert corrupted.identity.declared_sha256 != first.identity.declared_sha256
    assert corrupted.identity.source_sha256 != first.identity.source_sha256


def test_certificate_refuses_a_boundary_from_other_primitives() -> None:
    stored = _case("figure5").boundary.primitives
    index = next(i for i, p in enumerate(stored) if isinstance(p, ReferenceArc) and i != _anchor(stored) - 1)
    with pytest.raises(DeclaredProvenanceMismatchError):
        DeclaredBoundaryCertificate.build(_case("figure5"), reemit_primitives("figure5", _nudge_arc_end(stored, index)))


def test_one_ulp_corrupted_arc_end_is_repaired_and_the_stored_closure_check_raises() -> None:
    case = _case("figure5")
    stored = case.boundary.primitives
    index = next(i for i, p in enumerate(stored) if isinstance(p, ReferenceArc) and i != _anchor(stored) - 1)
    corrupted = _nudge_arc_end(stored, index)
    with pytest.raises(DisconnectedPublishedBoundaryError):
        ReferenceBoundary.build(corrupted, case.boundary.tool_radius, case.boundary.boundary_stroke_width)
    repaired = reemit_primitives("figure5", corrupted)
    following = repaired.primitives[(index + 1) % len(stored)]
    assert repaired.primitives[index].end == following.start
    moves = measure_moves(corrupted, repaired)
    assert moves[(index + 1) % len(stored)].start_gap > 0
    assert max(move.bound for move in moves) <= MEASURED_DISPLACEMENT_BOUND


def test_corrupting_the_closing_arc_end_leaves_the_cycle_open() -> None:
    stored = _case("figure5").boundary.primitives
    with pytest.raises(OpenDeclaredCycleError):
        reemit_primitives("figure5", _nudge_arc_end(stored, _anchor(stored) - 1))


def _p(x: int, y: int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(Fraction(x), Fraction(y))


def test_pythagorean_three_four_five_arc_has_exact_centre_and_radius() -> None:
    # Leaves (0, 0) along +x, counter-clockwise, through (8, 4): s = 80 / 8 = 10.
    arc = DeclaredArc.build(_p(0, 0), _p(8, 4), PythagoreanDirection.build(Fraction(0), 0), 1, 0)
    assert arc.radius == 10
    assert arc.centre == _p(0, 10)
    assert arc.end_tangent == (Fraction(3, 5), Fraction(4, 5))


def test_direction_contract_is_enforced() -> None:
    with pytest.raises(InvalidPythagoreanDirectionError):
        PythagoreanDirection.build(Fraction(1, 3), 0)
    with pytest.raises(InvalidPythagoreanDirectionError):
        PythagoreanDirection.build(Fraction(3, 4), 0)
    with pytest.raises(InvalidPythagoreanDirectionError):
        PythagoreanDirection.build(Fraction(0), 4)
    assert PythagoreanDirection.build(Fraction(1, 2), 1).vector == (Fraction(-4, 5), Fraction(3, 5))


def test_chord_away_from_the_centre_is_degenerate() -> None:
    with pytest.raises(DegenerateArcChordError):
        DeclaredArc.build(_p(0, 0), _p(8, -4), PythagoreanDirection.build(Fraction(0), 0), 1, 0)


def test_open_or_misordered_cycles_are_refused() -> None:
    east = PythagoreanDirection.build(Fraction(0), 0)
    north = PythagoreanDirection.build(Fraction(0), 1)
    west = PythagoreanDirection.build(Fraction(0), 2)
    south = PythagoreanDirection.build(Fraction(0), 3)
    square = tuple(
        DeclaredLine.build(start, direction, Fraction(1), index)
        for index, (start, direction) in enumerate(((_p(0, 0), east), (_p(1, 0), north), (_p(1, 1), west), (_p(0, 1), south)))
    )
    assert DeclaredBoundary.build("square", square, "").junction_census().left_turns == 4
    with pytest.raises(OpenDeclaredCycleError):
        DeclaredBoundary.build("square", (*square[:3], replace(square[3], length=Fraction(2))), "")
    with pytest.raises(InvalidDeclaredProvenanceError):
        DeclaredBoundary.build("square", (square[1], square[2], square[3], square[0]), "")


def test_square_root_enclosures_bracket_the_root() -> None:
    for value in (Fraction(2), Fraction(1, 3), Fraction(10**-30), Fraction(9, 4)):
        assert sqrt_lower(value) ** 2 <= value <= sqrt_upper(value) ** 2
    assert sqrt_lower(Fraction(9, 4)) == sqrt_upper(Fraction(9, 4)) == Fraction(3, 2)


def _stored_arc(start: tuple[float, float], end: tuple[float, float], centre: tuple[float, float], sweep: float) -> ReferenceArc:
    return ReferenceArc.build(Point2[WorldXY].build(start), Point2[WorldXY].build(end), Point2[WorldXY].build(centre), Radian(sweep))


def _distance_to_arc(point: tuple[float, float], centre: tuple[float, float], radius: float, start_angle: float, sweep: float) -> float:
    angle = math.atan2(point[1] - centre[1], point[0] - centre[0])
    offset = (angle - start_angle) % math.tau if sweep > 0 else (start_angle - angle) % math.tau
    if offset <= abs(sweep):
        return abs(math.dist(point, centre) - radius)
    ends = [(centre[0] + radius * math.cos(start_angle + a), centre[1] + radius * math.sin(start_angle + a)) for a in (0.0, sweep)]
    return min(math.dist(point, end) for end in ends)


def test_identical_arc_moves_by_exactly_zero() -> None:
    stored = _stored_arc((5.0, 0.0), (3.0, 4.0), (0.0, 0.0), math.atan2(4.0, 3.0))
    declared = DeclaredArc.build(_p(5, 0), _p(3, 4), PythagoreanDirection.build(Fraction(0), 1), 1, 0)
    assert declared.centre == _p(0, 0)
    assert arc_move_bound(stored, declared).bound == 0


@pytest.mark.parametrize("start_offset", [0.0, 2.0**-30])
def test_arc_move_bound_dominates_the_sampled_hausdorff_distance(start_offset: float) -> None:
    # Declared tangent rotated by 2^-19 rad: the declared circle differs from the stored one
    # by ~1e-6 mm over the arc, so float sampling resolves the true distance to ~1e-15.
    stored = _stored_arc((5.0, 0.0), (3.0, 4.0), (0.0, 0.0), math.atan2(4.0, 3.0))
    start = ExactPoint2[WorldXY].from_binary64(5.0, start_offset)
    declared = DeclaredArc.build(start, _p(3, 4), PythagoreanDirection.build(Fraction(1, 2**20), 1), 1, 0)
    bound = float(arc_move_bound(stored, declared).bound)
    centre = declared.centre.as_floats()
    radius = float(declared.radius)
    start_angle = math.atan2(start_offset - centre[1], 5.0 - centre[0])
    end_angle = math.atan2(4.0 - centre[1], 3.0 - centre[0])
    sweep = end_angle - start_angle
    samples = [index / 2000 for index in range(2001)]
    declared_points = [(centre[0] + radius * math.cos(start_angle + s * sweep), centre[1] + radius * math.sin(start_angle + s * sweep)) for s in samples]
    stored_sweep = math.atan2(4.0, 3.0)
    stored_points = [(5.0 * math.cos(s * stored_sweep), 5.0 * math.sin(s * stored_sweep)) for s in samples]
    sampled = max(
        max(_distance_to_arc(point, (0.0, 0.0), 5.0, 0.0, stored_sweep) for point in declared_points),
        max(_distance_to_arc(point, centre, radius, start_angle, sweep) for point in stored_points),
    )
    print(f"start offset {start_offset}: sampled {sampled:.6e} mm, bound {bound:.6e} mm")
    assert sampled > 0
    assert sampled <= bound
    assert bound <= 4 * sampled + 4 * start_offset
