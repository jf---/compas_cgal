"""Contract tests of the native clearance wire: exact text round trips, owner extraction, binary64 lanes refuse non-doubles."""

from __future__ import annotations

from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import strategies as st

from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_declared_clearance_wire import ARC_KIND
from benchmarks.held_declared_clearance_wire import LINE_KIND
from benchmarks.held_declared_clearance_wire import NonBinary64ClearanceInputError
from benchmarks.held_declared_clearance_wire import binary64_circle
from benchmarks.held_declared_clearance_wire import element_wire
from benchmarks.held_declared_clearance_wire import offset_owners
from benchmarks.held_declared_clearance_wire import primitive_wires
from benchmarks.held_declared_clearance_wire import rational_wire
from benchmarks.held_declared_clearance_wire import root_wire
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket

F = Fraction


def _rational(wire: tuple[str, str]) -> Fraction:
    return Fraction(int(wire[0]), int(wire[1]))


@given(st.fractions())
def test_rational_wire_round_trips_exactly(value: Fraction) -> None:
    numerator, denominator = rational_wire(value)
    assert int(denominator) > 0
    assert _rational((numerator, denominator)) == value


@given(st.fractions(), st.fractions(), st.fractions(min_value=0))
def test_root_wire_carries_the_folded_representation(a: Fraction, b: Fraction, c: Fraction) -> None:
    value = OneRootNumber.build(a, b, c)
    wire = root_wire(value)
    assert (_rational(wire[0]), _rational(wire[1]), _rational(wire[2])) == (value.rational, value.coefficient, value.radicand)


def test_primitive_wires_encode_lines_and_arcs_with_their_exact_data() -> None:
    boundary = concave_arc_pocket()
    wires = primitive_wires(boundary)
    assert len(wires) == len(boundary.primitives)
    for wire, primitive in zip(wires, boundary.primitives, strict=True):
        kind, sense, values = wire
        assert len(values) == 7
        assert (_rational(values[0]), _rational(values[1])) == (primitive.start.x, primitive.start.y)
        assert (_rational(values[2]), _rational(values[3])) == (primitive.end.x, primitive.end.y)
        if isinstance(primitive, DeclaredArc):
            assert kind == ARC_KIND
            assert sense == primitive.sense
            assert (_rational(values[4]), _rational(values[5]), _rational(values[6])) == (primitive.centre.x, primitive.centre.y, primitive.radius)
        else:
            assert kind == LINE_KIND
            assert (_rational(values[4]), _rational(values[5])) == primitive.direction.vector
            assert _rational(values[6]) == primitive.length


def test_element_wires_carry_owners_and_support_data() -> None:
    start = OneRootPoint2.build(OneRootNumber.build(F(1), F(1), F(2)), OneRootNumber.of(F(3)))
    end = OneRootPoint2.rational(F(5), F(3))
    segment = SegmentElement(start, end, (F(1), F(0)), (F(0), F(1)), F(3))
    kind, sense, owners, start_wire, _, values = element_wire(segment, (4, 7))
    assert (kind, sense, owners) == (LINE_KIND, 0, [4, 7])
    assert start_wire[0] == root_wire(start.x)
    assert [_rational(value) for value in values] == [F(1), F(0), F(0), F(1), F(3)]
    arc = ArcElement(ExactPoint2[WorldXY].build(F(0), F(0)), F(5), -1, OneRootPoint2.rational(F(3), F(4)), OneRootPoint2.rational(F(5), F(0)))
    kind, sense, owners, _, _, values = element_wire(arc, ())
    assert (kind, sense, owners) == (ARC_KIND, -1, [])
    assert [_rational(value) for value in values] == [F(0), F(0), F(5)]


def test_offset_owners_keep_line_and_arc_features_only() -> None:
    owners = frozenset({BoundaryFeature(BoundaryFeatureKind.ARC, 5), BoundaryFeature(BoundaryFeatureKind.JUNCTION, 2), BoundaryFeature(BoundaryFeatureKind.LINE, 1)})
    assert offset_owners(owners) == (1, 5)


def test_binary64_lane_refuses_a_value_that_is_not_a_double() -> None:
    centre = ExactPoint2[WorldXY].build(F(1, 3), F(0))
    with pytest.raises(NonBinary64ClearanceInputError):
        binary64_circle(CircleElement(centre, F(1)))
    assert binary64_circle(CircleElement(ExactPoint2[WorldXY].from_binary64(0.1, 2.5), F(0.75))) == (0.1, 2.5, 0.75)
