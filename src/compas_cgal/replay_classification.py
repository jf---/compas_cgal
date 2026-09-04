"""Shared cut-plane replay classification for production and evidence consumers."""

import math
from typing import Iterable
from typing import Literal
from typing import TypeAlias

from compas.geometry import Line
from compas.tolerance import TOL

from compas_cgal.adaptive.units import Millimetre
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

ReplayCategory: TypeAlias = Literal["motion", "rapid", "plunge", "retract"]

AUDIT_ENGAGED = frozenset(
    {
        OperationType.CUT,
        OperationType.LEAD_IN,
        OperationType.LEAD_OUT,
        OperationType.LINK,
    }
)


class CutPlaneRampError(ValueError):
    """A differing-Z line also travels in XY outside the cut-plane model."""


class OffPlaneReplayCurveError(ValueError):
    """An engaged arc or circle lies outside the inferred cutting plane."""


def minimum_cut_height(heights: Iterable[Millimetre]) -> Millimetre:
    """Return the lowest supplied motion height, or zero for no motions."""
    return min(heights, default=Millimetre(0.0))


def classify_line_replay(
    operation: OperationType,
    *,
    start_z: Millimetre,
    end_z: Millimetre,
    xy_travel: Millimetre,
    cut_z: Millimetre,
) -> ReplayCategory:
    """Classify a line under the single-plane replay contract."""
    if operation is OperationType.RETRACT:
        return "retract"
    if abs(start_z - end_z) > TOL.absolute:
        if xy_travel > TOL.absolute:
            raise CutPlaneRampError(f"differing-z line has nonzero XY travel (len={xy_travel:.3e})")
        return "plunge" if end_z < start_z else "rapid"
    if start_z > cut_z + TOL.absolute:
        return "rapid"
    if operation not in AUDIT_ENGAGED:
        return "rapid"
    return "motion"


def classify_planar_replay(operation: OperationType, *, motion_z: Millimetre, cut_z: Millimetre) -> ReplayCategory:
    """Classify a planar arc or circle under the replay contract."""
    if operation is OperationType.RETRACT:
        return "retract"
    if operation not in AUDIT_ENGAGED:
        return "rapid"
    if operation is OperationType.LINK and motion_z > cut_z + TOL.absolute:
        return "rapid"
    if not TOL.is_between(motion_z, cut_z, cut_z, atol=TOL.absolute):
        raise OffPlaneReplayCurveError("engaged curve lies outside the inferred cut plane")
    return "motion"


def classify_operation_replay(operation: ToolpathOperation, cut_z: Millimetre) -> ReplayCategory:
    """Classify one typed toolpath operation without stock side effects."""
    geometry = operation.geometry
    if isinstance(geometry, Line):
        return classify_line_replay(
            operation.operation,
            start_z=Millimetre(float(geometry.start[2])),
            end_z=Millimetre(float(geometry.end[2])),
            xy_travel=Millimetre(
                math.hypot(
                    float(geometry.end[0]) - float(geometry.start[0]),
                    float(geometry.end[1]) - float(geometry.start[1]),
                )
            ),
            cut_z=cut_z,
        )
    return classify_planar_replay(
        operation.operation,
        motion_z=Millimetre(float(geometry.frame.point[2])),
        cut_z=cut_z,
    )
