from typing import assert_never
from typing import assert_type

from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.operation import TraversalDecision
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutZ
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

approach_operation = ApproachOperation.build(
    position=world_point_scalar,
    clearance_z=ClearanceZ.build(5.0),
)
plunge_operation = PlungeOperation.build(
    position=world_point_scalar,
    clearance_z=ClearanceZ.build(5.0),
    cut_z=CutZ.build(-2.0),
)
assert_type(approach_operation, ApproachOperation)
assert_type(plunge_operation, PlungeOperation)


def accepts_world_point(point: Point2[WorldXY]) -> None:
    pass


machine_point = Point2[MachineXY].build(1.0, 2.0)
accepts_world_point(machine_point)  # type: ignore[arg-type]
accepts_world_point(world_vector_scalar)  # type: ignore[arg-type]


def operation_kind(operation: CanonicalOperation) -> str:
    if isinstance(operation, ApproachOperation):
        return "approach"
    if isinstance(operation, PlungeOperation):
        return "plunge"
    if isinstance(operation, LinkSegmentOperation):
        return "link"
    if isinstance(operation, CutFullCircleOperation):
        return "circle"
    assert_never(operation)


def cap_kind(decision: FullCapDecision | NeckCapDecision) -> str:
    if isinstance(decision, FullCapDecision):
        return "full"
    if isinstance(decision, NeckCapDecision):
        return "neck"
    assert_never(decision)


def traversal_kind(decision: TraversalDecision) -> str:
    if isinstance(decision, HoldTraversalDecision):
        return "hold"
    if isinstance(decision, AdvanceTraversalDecision):
        return "advance"
    assert_never(decision)
