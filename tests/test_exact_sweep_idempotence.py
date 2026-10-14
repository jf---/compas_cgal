"""Replay idempotence: a swept region re-swept, or entered at its own entry, changes nothing.

The replay (`compas_cgal.engagement._subtract_operation`) removes a loop as the
exact annulus of its AUTHORED radius and a straight cut as one exact-or-certified
line sweep, so re-cutting what is already cut is an exact no-op -- not a sliver.
"""

import math

import numpy as np
import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

from compas_cgal.engagement import _subtract_operation
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

TOOL_RADIUS = 1.0
GUIDE_RADIUS = 2.5
CENTRE = (20.0, 20.0)
TRIALS = 12
# The arrangement after four 16 mm wall runs: the disk chain built 8,014 vertices.
FOUR_RUN_VERTEX_CEILING = 200


def _square() -> Stock:
    return Stock(Polygon([[0, 0, 0], [40, 0, 0], [40, 40, 0], [0, 40, 0]]))


def _rect(width: float, height: float) -> Stock:
    return Stock(Polygon([[0, 0, 0], [width, 0, 0], [width, height, 0], [0, height, 0]]))


def _loop(centre: tuple[float, float], radius: float, phase: float) -> ToolpathOperation:
    ux, uy = math.cos(phase), math.sin(phase)
    frame = Frame([centre[0], centre[1], 0.0], [ux, uy, 0.0], [-uy, ux, 0.0])
    return ToolpathOperation(Circle(radius, frame=frame), OperationType.CUT, 0, True, np.zeros(3), np.zeros(3))


def _cut(x0: float, y0: float, x1: float, y1: float) -> ToolpathOperation:
    return ToolpathOperation(Line([x0, y0, 0.0], [x1, y1, 0.0]), OperationType.CUT, 0, False, np.zeros(3), np.zeros(3))


@pytest.mark.parametrize("k", range(TRIALS))
def test_a_loop_re_swept_from_any_start_changes_nothing(k: int) -> None:
    stock = _square()
    _subtract_operation(stock, _loop(CENTRE, GUIDE_RADIUS, 2 * math.pi * k / TRIALS), TOOL_RADIUS)
    again = stock.clone()
    _subtract_operation(again, _loop(CENTRE, GUIDE_RADIUS, 2 * math.pi * k / TRIALS + 1.0), TOOL_RADIUS)
    assert again.exactly_equals(stock)


@pytest.mark.parametrize("k", range(TRIALS))
def test_a_plunge_at_a_loops_own_entry_changes_nothing(k: int) -> None:
    loop = IncidentLoop.build(Point2[WorldXY].build(*CENTRE), GuideRadius.build(GUIDE_RADIUS), Radian(2 * math.pi * k / TRIALS))
    stock = _square()
    _subtract_operation(stock, _loop((loop.centre.x, loop.centre.y), loop.radius.value, loop.phase), TOOL_RADIUS)
    plunged = stock.clone()
    plunged.subtract_disk(loop.entry.x, loop.entry.y, TOOL_RADIUS)
    assert plunged.exactly_equals(stock)


def test_four_wall_runs_leave_a_small_arrangement() -> None:
    stock = _rect(20.0, 12.0)
    for y in (1.5, 3.5, 8.5, 10.5):
        _subtract_operation(stock, _cut(2.0, y, 18.0, y), TOOL_RADIUS)
    vertices, _halfedges, _faces = stock.arrangement_stats()
    assert vertices < FOUR_RUN_VERTEX_CEILING
    assert stock.representation_is_valid()


@pytest.mark.parametrize(
    ("outer", "inner"),
    [
        pytest.param((2.0, 4.0, 8.0, 4.0), (3.0, 4.0, 7.0, 4.0), id="axis-exact-capsule"),
        pytest.param((2.0, 2.0, 8.0, 10.0), (3.5, 4.0, 5.0, 6.0), id="pythagorean-3-4-5"),
        pytest.param((2.0, 2.0, 8.0, 10.0), (5.0, 6.0, 3.5, 4.0), id="pythagorean-reversed"),
    ],
)
def test_an_exact_cut_inside_an_earlier_collinear_cut_changes_nothing(outer: tuple[float, float, float, float], inner: tuple[float, float, float, float]) -> None:
    stock = _rect(12.0, 12.0)
    _subtract_operation(stock, _cut(*outer), TOOL_RADIUS)
    again = stock.clone()
    _subtract_operation(again, _cut(*inner), TOOL_RADIUS)
    assert again.exactly_equals(stock)


@pytest.mark.parametrize(
    ("outer", "inner"),
    [
        pytest.param((2.0, 3.0, 8.0, 5.0), (3.5, 3.5, 6.5, 4.5), id="oblique"),
        pytest.param((2.0, 3.0, 8.0, 5.0), (6.5, 4.5, 3.5, 3.5), id="oblique-reversed"),
        pytest.param((2.0, 3.0, 8.0, 5.0), (2.0, 3.0, 8.0, 5.0), id="oblique-identical"),
    ],
)
def test_an_oblique_collinear_recut_removes_only_inside_the_outer_over_region(outer: tuple[float, float, float, float], inner: tuple[float, float, float, float]) -> None:
    """The exact property of the end-disk quad.

    A collinear re-cut's full-radius end disks may poke into the earlier cut's
    side band, never past its certified OVER region. Decided with exact set
    operations: removed = before - after lies inside OVER iff
    (before - OVER) == (after - OVER).
    """
    before = _rect(12.0, 12.0)
    _subtract_operation(before, _cut(*outer), TOOL_RADIUS)
    after = before.clone()
    _subtract_operation(after, _cut(*inner), TOOL_RADIUS)
    assert after.is_subset_of(before)
    before_outside = before.clone()
    before_outside.raw.subtract_line_sweep_over(*outer, TOOL_RADIUS)
    after_outside = after.clone()
    after_outside.raw.subtract_line_sweep_over(*outer, TOOL_RADIUS)
    assert after_outside.exactly_equals(before_outside)


@pytest.mark.parametrize(
    "bridge",
    [
        pytest.param((2.0, 3.0, 8.0, 5.0), id="oblique"),
        pytest.param((1.3, 8.7, 8.9, 1.1), id="oblique-steep"),
        pytest.param((4.9728933988224115, 2.8516680046863603, 4.962616443776326, 2.7917695211359552), id="gate-bridge"),
    ],
)
def test_a_full_radius_disk_at_a_non_pythagorean_bridge_end_changes_nothing(bridge: tuple[float, float, float, float]) -> None:
    """The junction case: a plunge, or the next motion's start, at either end.

    The ends are read back from the compas `Line`, which recomputes its end from
    start + direction * length (so (8.9, 1.1) comes back as 1.0999999999999996):
    the junction is with the end the motion actually sweeps.
    """
    stock = _rect(12.0, 12.0)
    operation = _cut(*bridge)
    _subtract_operation(stock, operation, TOOL_RADIUS)
    line = operation.geometry
    assert isinstance(line, Line)
    for x, y in ((float(line.start[0]), float(line.start[1])), (float(line.end[0]), float(line.end[1]))):
        plunged = stock.clone()
        plunged.subtract_disk(x, y, TOOL_RADIUS)
        assert plunged.exactly_equals(stock)


def test_a_pythagorean_line_sweep_is_the_exact_capsule() -> None:
    # Exact: a point at r(1 - 1e-6) from the axis is removed interior. The quad
    # (half-width (1 - 5e-5) r) and the chain (slack 1e-4 r) would both keep it.
    stock = _rect(12.0, 12.0)
    stock.subtract_line_sweep(2.0, 6.0, 10.0, 6.0, 1.0)
    assert not stock.contains(6.0, 6.0 + 1.0 - 1e-6)
    assert not stock.contains(6.0, 6.0 - 1.0 + 1e-6)
    assert stock.contains(6.0, 6.0 + 1.0 + 1e-6)
