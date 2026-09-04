"""Behavioral contract for the shared cut-plane replay decision."""

import pytest
from compas.tolerance import TOL

from compas_cgal.adaptive.units import Millimetre
from compas_cgal.replay_classification import CutPlaneRampError
from compas_cgal.replay_classification import OffPlaneReplayCurveError
from compas_cgal.replay_classification import classify_line_replay
from compas_cgal.replay_classification import classify_planar_replay
from compas_cgal.toolpath import OperationType

CUT_HEIGHT = Millimetre(0.0)
AT_TOLERANCE = Millimetre(TOL.absolute)
ABOVE_TOLERANCE = Millimetre(2.0 * TOL.absolute)


@pytest.mark.parametrize(
    ("operation", "start_z", "end_z", "xy_travel", "expected"),
    [
        (OperationType.CUT, CUT_HEIGHT, CUT_HEIGHT, Millimetre(1.0), "motion"),
        (OperationType.CUT, AT_TOLERANCE, AT_TOLERANCE, Millimetre(1.0), "motion"),
        (OperationType.CUT, ABOVE_TOLERANCE, ABOVE_TOLERANCE, Millimetre(1.0), "rapid"),
        (OperationType.CUT, ABOVE_TOLERANCE, CUT_HEIGHT, AT_TOLERANCE, "plunge"),
        (OperationType.CUT, CUT_HEIGHT, ABOVE_TOLERANCE, AT_TOLERANCE, "rapid"),
        (OperationType.PLUNGE, CUT_HEIGHT, CUT_HEIGHT, Millimetre(1.0), "rapid"),
        (OperationType.RETRACT, ABOVE_TOLERANCE, CUT_HEIGHT, ABOVE_TOLERANCE, "retract"),
    ],
)
def test_line_replay_boundaries(
    operation: OperationType,
    start_z: Millimetre,
    end_z: Millimetre,
    xy_travel: Millimetre,
    expected: str,
) -> None:
    assert (
        classify_line_replay(
            operation,
            start_z=start_z,
            end_z=end_z,
            xy_travel=xy_travel,
            cut_z=CUT_HEIGHT,
        )
        == expected
    )


def test_line_replay_refuses_xy_travel_just_above_tolerance() -> None:
    with pytest.raises(CutPlaneRampError):
        classify_line_replay(
            OperationType.CUT,
            start_z=ABOVE_TOLERANCE,
            end_z=CUT_HEIGHT,
            xy_travel=ABOVE_TOLERANCE,
            cut_z=CUT_HEIGHT,
        )


@pytest.mark.parametrize(
    ("operation", "motion_z", "expected"),
    [
        (OperationType.CUT, CUT_HEIGHT, "motion"),
        (OperationType.CUT, AT_TOLERANCE, "motion"),
        (OperationType.LINK, ABOVE_TOLERANCE, "rapid"),
        (OperationType.PLUNGE, ABOVE_TOLERANCE, "rapid"),
        (OperationType.RETRACT, ABOVE_TOLERANCE, "retract"),
    ],
)
def test_planar_replay_boundaries(operation: OperationType, motion_z: Millimetre, expected: str) -> None:
    assert classify_planar_replay(operation, motion_z=motion_z, cut_z=CUT_HEIGHT) == expected


def test_planar_replay_refuses_engaged_curve_off_plane() -> None:
    with pytest.raises(OffPlaneReplayCurveError):
        classify_planar_replay(OperationType.CUT, motion_z=ABOVE_TOLERANCE, cut_z=CUT_HEIGHT)
