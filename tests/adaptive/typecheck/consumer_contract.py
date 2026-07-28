from typing import assert_never
from typing import assert_type

import numpy as np

from compas_cgal import _medial_axis_2
from compas_cgal._medial_axis_2 import MedialAxisResult
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
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock


class MachineXY:
    """Distinct machine-coordinate frame used only by this compile contract."""


mat_result = _medial_axis_2.segment_site_medial_axis(
    np.asarray(
        (
            (0.0, 0.0, 0.0),
            (6.0, 0.0, 0.0),
            (6.0, 2.0, 0.0),
            (2.0, 2.0, 0.0),
            (2.0, 6.0, 0.0),
            (0.0, 6.0, 0.0),
        ),
        dtype=np.float64,
    ),
    (),
    0.5,
    0.75,
    0.02,
    32,
)
assert_type(mat_result, MedialAxisResult)

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

stock_area = Stock2Area.build(Stock.__new__(Stock))
depletion_policy = DepletionPolicy.build(
    chord_bound=ChordBound.build(0.25),
    center_count_limit=256,
)
segment_motion = ExactSegmentMotion.build(
    Point2[WorldXY].build(0.0, 0.0),
    Point2[WorldXY].build(1.0, 0.0),
)
circle_motion = ExactCircleMotion.build(
    Point2[WorldXY].build(0.0, 0.0),
    Vector2[WorldXY].build(1.0, 0.0),
    False,
)
segment_witness = stock_area.deplete(segment_motion, ToolRadius.build(0.5), depletion_policy)
circle_witness = stock_area.deplete(circle_motion, ToolRadius.build(0.5), depletion_policy)
assert_type(segment_witness, DepletionWitness)
assert_type(circle_witness, DepletionWitness)


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
