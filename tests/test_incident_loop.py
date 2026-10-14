"""Loop entries authored EXACTLY on their guide circle, with an exact certificate."""

import math
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import strategies as st

from compas_cgal.errors import InvalidUnitValueError
from compas_cgal.incident_loop import INCIDENCE_CENTRE_BOUND_MM
from compas_cgal.incident_loop import INCIDENCE_CHART_DENOMINATOR
from compas_cgal.incident_loop import INCIDENCE_COORDINATE_BOUND_MM
from compas_cgal.incident_loop import INCIDENCE_PHASE_BOUND_RAD
from compas_cgal.incident_loop import INCIDENCE_RADIUS_BOUND_MM
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.incident_loop import IncidentLoopRadiusBelowQuantumError
from compas_cgal.incident_loop import IncidentLoopRangeError
from compas_cgal.incident_loop import NonFiniteIncidentLoopPhaseError
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

coordinates = st.floats(-INCIDENCE_COORDINATE_BOUND_MM, INCIDENCE_COORDINATE_BOUND_MM, allow_nan=False)
radii = st.floats(1e-3, INCIDENCE_COORDINATE_BOUND_MM, allow_nan=False)
phases = st.floats(-4 * math.pi, 4 * math.pi, allow_nan=False)


def _angle_gap(a: float, b: float) -> float:
    return abs(math.remainder(a - b, math.tau))


def _build(cx: float, cy: float, radius: float, phase: float) -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(cx, cy), GuideRadius.build(radius), Radian(phase))


def test_the_quantisation_constants_meet_the_ruling() -> None:
    assert INCIDENCE_CHART_DENOMINATOR >= 2**10
    assert INCIDENCE_PHASE_BOUND_RAD <= 1.0 / 2**10
    assert INCIDENCE_RADIUS_BOUND_MM <= 1e-5
    assert INCIDENCE_COORDINATE_BOUND_MM >= 2.0**10


@given(cx=coordinates, cy=coordinates, radius=radii, phase=phases)
def test_an_authored_entry_lies_exactly_on_its_guide(cx: float, cy: float, radius: float, phase: float) -> None:
    loop = _build(cx, cy, radius, phase)
    dx = Fraction(loop.entry.x) - Fraction(loop.centre.x)
    dy = Fraction(loop.entry.y) - Fraction(loop.centre.y)
    assert dx * dx + dy * dy == Fraction(loop.radius.value) ** 2
    assert loop.certifies_incidence()


@given(cx=coordinates, cy=coordinates, radius=radii, phase=phases)
def test_authoring_moves_the_loop_by_at_most_the_declared_quanta(cx: float, cy: float, radius: float, phase: float) -> None:
    loop = _build(cx, cy, radius, phase)
    assert abs(loop.radius.value - radius) <= INCIDENCE_RADIUS_BOUND_MM
    assert _angle_gap(loop.phase, phase) <= INCIDENCE_PHASE_BOUND_RAD
    assert abs(loop.centre.x - cx) <= INCIDENCE_CENTRE_BOUND_MM
    assert abs(loop.centre.y - cy) <= INCIDENCE_CENTRE_BOUND_MM
    # The realised phase is the direction of the exact entry, not a label.
    assert _angle_gap(math.atan2(loop.entry.y - loop.centre.y, loop.entry.x - loop.centre.x), loop.phase) <= 1e-9


def test_axis_phases_author_the_axis_entry() -> None:
    loop = _build(10.0, 6.0, 1.7, 0.0)
    assert loop.entry.y == loop.centre.y
    assert loop.entry.x == loop.centre.x + loop.radius.value


def test_the_circle_carries_the_authored_radius_and_centre() -> None:
    loop = _build(3.3, 4.1, 1.7, 1.0)
    circle = loop.circle()
    assert float(circle.radius) == loop.radius.value
    assert (float(circle.frame.point[0]), float(circle.frame.point[1])) == (loop.centre.x, loop.centre.y)


@pytest.mark.parametrize(
    ("centre", "radius", "phase", "error"),
    [
        ((0.0, 0.0), 1.0, float("nan"), NonFiniteIncidentLoopPhaseError),
        ((0.0, 0.0), 1.0, float("inf"), NonFiniteIncidentLoopPhaseError),
        ((2.0**11, 0.0), 1.0, 0.0, IncidentLoopRangeError),
        ((0.0, 0.0), 2.0**11, 0.0, IncidentLoopRangeError),
        ((0.0, 0.0), 1e-15, 0.5, IncidentLoopRadiusBelowQuantumError),
    ],
)
def test_unsupported_input_is_refused_by_name(centre: tuple[float, float], radius: float, phase: float, error: type[Exception]) -> None:
    with pytest.raises(error):
        _build(centre[0], centre[1], radius, phase)


def test_the_unit_types_refuse_non_finite_and_non_positive_input() -> None:
    with pytest.raises(InvalidUnitValueError):
        _build(float("nan"), 0.0, 1.0, 0.0)
    with pytest.raises(InvalidUnitValueError):
        _build(0.0, 0.0, 0.0, 0.0)
    with pytest.raises(InvalidUnitValueError):
        _build(0.0, 0.0, -1.0, 0.0)
