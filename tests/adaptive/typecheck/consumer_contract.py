from typing import assert_type

from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


class MachineXY:
    """Distinct machine-coordinate frame used only by this compile contract."""


world_point_scalar = Point2[WorldXY].build(1.0, 2.0)
world_point_sequence = Point2[WorldXY].build([1.0, 2.0])
world_vector_scalar = Vector2[WorldXY].build(1.0, 2.0)
world_vector_sequence = Vector2[WorldXY].build((1.0, 2.0))

assert_type(world_point_scalar, Point2[WorldXY])
assert_type(world_point_sequence, Point2[WorldXY])
assert_type(world_vector_scalar, Vector2[WorldXY])
assert_type(world_vector_sequence, Vector2[WorldXY])


def accepts_world_point(point: Point2[WorldXY]) -> None:
    pass


machine_point = Point2[MachineXY].build(1.0, 2.0)
accepts_world_point(machine_point)  # type: ignore[arg-type]
accepts_world_point(world_vector_scalar)  # type: ignore[arg-type]
