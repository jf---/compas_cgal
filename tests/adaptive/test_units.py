import math
from collections.abc import Callable
from fractions import Fraction
from typing import assert_type

import pytest

import compas_cgal.adaptive as adaptive
from compas_cgal.adaptive.errors import InvalidUnitValueError
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import Clearance
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


def test_point_factory_accepts_scalar_and_sequence_forms() -> None:
    scalar = Point2[WorldXY].build(1.25, -2.5)
    sequence = Point2[WorldXY].build([1.25, -2.5])

    assert scalar == sequence
    assert_type(scalar, Point2[WorldXY])


def test_vector_factory_accepts_scalar_and_sequence_forms() -> None:
    scalar = Vector2[WorldXY].build(3.0, 4.0)
    sequence = Vector2[WorldXY].build((3.0, 4.0))

    assert scalar == sequence
    assert_type(scalar, Vector2[WorldXY])


@pytest.mark.parametrize(
    ("factory", "value"),
    [
        (ToolRadius.build, 0.0),
        (ToolRadius.build, -1.0),
        (EntryRadius.build, 0.0),
        (Spacing.build, 0.0),
        (ChordBound.build, -1.0),
        (GuideRadius.build, 0.0),
        (Clearance.build, -1.0),
        (CutZ.build, math.inf),
        (ClearanceZ.build, math.nan),
        (ToolRadius.build, True),
    ],
)
def test_unit_factories_reject_invalid_values(factory: Callable[[float], object], value: float) -> None:
    with pytest.raises(InvalidUnitValueError):
        factory(value)


def test_cut_plane_requires_clearance_above_cut() -> None:
    cut_z = CutZ.build(-3.0)

    with pytest.raises(InvalidUnitValueError):
        CutPlane.build(cut_z, ClearanceZ.build(-3.0))
    with pytest.raises(InvalidUnitValueError):
        CutPlane.build(cut_z, ClearanceZ.build(-4.0))

    plane = CutPlane.build(cut_z, ClearanceZ.build(5.0))
    assert plane.cut_z == cut_z
    assert plane.clearance_z == ClearanceZ.build(5.0)


def test_coordinates_reject_non_finite_components() -> None:
    with pytest.raises(InvalidUnitValueError):
        Point2[WorldXY].build(math.nan, 0.0)
    with pytest.raises(InvalidUnitValueError):
        Vector2[WorldXY].build([0.0, math.inf])
    with pytest.raises(InvalidUnitValueError):
        Point2[WorldXY].build(True, 0.0)


@pytest.mark.parametrize("components", [[], [1.0], [1.0, 2.0, 3.0]])
def test_coordinate_sequences_have_exactly_two_components(components: list[float]) -> None:
    with pytest.raises(InvalidUnitValueError):
        Point2[WorldXY].build(components)
    with pytest.raises(InvalidUnitValueError):
        Vector2[WorldXY].build(components)


def test_raw_constructors_cannot_bypass_unit_invariants() -> None:
    with pytest.raises(InvalidUnitValueError):
        ToolRadius(Millimetre(0.0))
    with pytest.raises(InvalidUnitValueError):
        Point2[WorldXY](Millimetre(math.nan), Millimetre(0.0))

    vector = Vector2[WorldXY](Millimetre(3), Millimetre(4))  # type: ignore[arg-type]
    assert type(vector.x) is float
    assert type(vector.y) is float


def test_squared_millimetre_retains_exact_rational_value() -> None:
    exact = SquaredMillimetre(Fraction(9, 4))
    assert exact == Fraction(9, 4)


def test_adaptive_package_does_not_publish_an_all_manifest() -> None:
    assert not hasattr(adaptive, "__all__")
