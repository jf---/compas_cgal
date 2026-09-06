"""Exact residual checks at the emitted-operation consumer boundary."""

import numpy as np
import pytest
from compas.geometry import Arc, Circle, Frame, Line, Polygon

from benchmarks.held_exact_motion_coverage import IncompleteMotionCoverageError
from benchmarks.spec import PocketSpec
from benchmarks.toolpath_coverage import UnsupportedCoverageMotionError, require_toolpath_coverage, replay_toolpath_coverage
from compas_cgal.toolpath import OperationType, ToolpathOperation, ToolpathResult


def _spec() -> PocketSpec:
    return PocketSpec.build("coverage", "test", Polygon([[-2, -2, 0], [2, -2, 0], [2, 2, 0], [-2, 2, 0]]), 2.0, 80.0)


def _path(*geometries: Line | Circle | Arc) -> ToolpathResult:
    return ToolpathResult([ToolpathOperation(g, OperationType.CUT, 0) for g in geometries], np.zeros((0, 3)))


def test_full_circle_preserves_uncut_centre_without_actual_plunge() -> None:
    result = _path(Circle(1.5))
    residual = replay_toolpath_coverage(_spec(), result)
    assert residual.contains(0.0, 0.0)
    result.operations.insert(0, ToolpathOperation(Line([0, 0, 5], [0, 0, 0]), OperationType.PLUNGE, 0))
    assert not replay_toolpath_coverage(_spec(), result).contains(0.0, 0.0)


def test_actual_cutting_segment_removes_capsule_but_rapid_does_not() -> None:
    result = _path(Line([-2, 0, 0], [2, 0, 0]))
    residual = replay_toolpath_coverage(_spec(), result)
    assert not residual.contains(0.0, 0.5)
    assert residual.contains(0.0, 1.5)
    result.operations[0] = ToolpathOperation(Line([-2, 0, 5], [2, 0, 5]), OperationType.LINK, 0)
    assert replay_toolpath_coverage(_spec(), result).contains(0.0, 0.5)


def test_incomplete_path_cannot_be_accepted_as_length_evidence() -> None:
    with pytest.raises(IncompleteMotionCoverageError, match="residual"):
        require_toolpath_coverage(_spec(), _path(Circle(1.5)))
    require_toolpath_coverage(_spec(), _path(Circle(0.5), Line([-2, -2, 0], [2, -2, 0]), Line([-2, 0, 0], [2, 0, 0]), Line([-2, 2, 0], [2, 2, 0])))


def test_partial_arc_and_tilted_circle_fail_loudly() -> None:
    with pytest.raises(UnsupportedCoverageMotionError, match="partial arc"):
        replay_toolpath_coverage(_spec(), _path(Arc(1.0, 0.0, 1.0)))
    with pytest.raises(UnsupportedCoverageMotionError, match="world XY"):
        replay_toolpath_coverage(_spec(), _path(Circle(1.0, frame=Frame.worldYZ())))


@pytest.mark.parametrize("consumer", ["repository", "controlled", "constant_spacing"])
def test_figure6_consumers_reject_incomplete_cutting_paths(monkeypatch: pytest.MonkeyPatch, consumer: str) -> None:
    from benchmarks import figure6, held_figure6_comparison, mathsm

    incomplete = _path(Circle(1.5))
    monkeypatch.setattr(figure6, "controlled_path", lambda *args: incomplete)
    monkeypatch.setattr(held_figure6_comparison, "controlled_path", lambda *args: incomplete)
    monkeypatch.setattr(held_figure6_comparison, "reference_pocket", _spec)
    monkeypatch.setattr(mathsm, "constant_spacing_path", lambda *args: incomplete)
    with pytest.raises(IncompleteMotionCoverageError, match="residual"):
        if consumer == "repository":
            held_figure6_comparison.measure_repository_figure6()
        elif consumer == "controlled":
            figure6.figure6_points(_spec(), [80.0], [])
        else:
            mathsm.measure_constant_spacing(_spec(), 0.1)
