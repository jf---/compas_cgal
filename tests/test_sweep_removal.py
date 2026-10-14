"""The replay's redundancy verdict: one exact interior test, three honest answers.

`compas_cgal.engagement._operation_removal` brackets what a motion removes
between the region the replay subtracts (UNDER) and a certified region that
contains the true swept area (OVER):

* REMOVES      -- the replay's own region meets the stock's interior;
* NOTHING      -- even the over-region misses it, so the tool cuts no material;
* WITHIN_BAND  -- neither: material survives only in the certified band.

For EXACT sweeps (circles, Pythagorean lines) under and over coincide, so the
verdict is two-valued and must agree with clone + subtract + compare.
"""

import math

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.engagement import UnexpectedToolpathGeometryError
from compas_cgal.engagement import _operation_removal
from compas_cgal.engagement import _subtract_operation
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation

TOOL_RADIUS = 1.0
LATTICE = st.integers(0, 40).map(lambda q: q / 4.0)
PYTHAGOREAN_STEPS = ((3, 4), (4, -3), (1, 0), (0, 1), (5, 12))


def _rect() -> Stock:
    return Stock(Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]))


def _line(x0: float, y0: float, x1: float, y1: float) -> ToolpathOperation:
    return ToolpathOperation(Line([x0, y0, 0.0], [x1, y1, 0.0]), OperationType.CUT, 0, False, np.zeros(3), np.zeros(3))


def _circle(cx: float, cy: float, radius: float, phase: float = 0.0) -> ToolpathOperation:
    frame = Frame([cx, cy, 0.0], [math.cos(phase), math.sin(phase), 0.0], [-math.sin(phase), math.cos(phase), 0.0])
    return ToolpathOperation(Circle(radius, frame=frame), OperationType.CUT, 0, True, np.zeros(3), np.zeros(3))


@st.composite
def _exact_operations(draw: st.DrawFn) -> ToolpathOperation:
    x0, y0 = draw(LATTICE), draw(LATTICE)
    if draw(st.booleans()):
        return _circle(x0, y0, draw(st.sampled_from((0.5, 1.0, 1.5, 2.5))), draw(st.floats(0.0, 6.0)))
    dx, dy = draw(st.sampled_from(PYTHAGOREAN_STEPS))
    scale = draw(st.sampled_from((0.25, 0.5, 1.0)))
    return _line(x0, y0, x0 + scale * dx, y0 + scale * dy)


@settings(max_examples=60, deadline=None)
@given(operations=st.lists(_exact_operations(), min_size=1, max_size=8))
def test_an_exact_verdict_agrees_with_clone_subtract_compare(operations: list[ToolpathOperation]) -> None:
    stock = _rect()
    for operation in operations + operations[:1]:
        probe = stock.clone()
        _subtract_operation(probe, operation, TOOL_RADIUS)
        changes = not probe.exactly_equals(stock)
        verdict = _operation_removal(stock, operation, TOOL_RADIUS)
        assert verdict is not SweepRemoval.WITHIN_BAND
        assert (verdict is SweepRemoval.REMOVES) is changes
        stock = probe


def test_engagement_cannot_decide_redundancy_in_either_direction() -> None:
    """The two motions of the gate's reference test, through the verdict itself."""
    stock = Stock(Polygon([[-10, -10, 0], [10, -10, 0], [10, 10, 0], [-10, 10, 0]]))
    stock.subtract_disk(3.0, 0.0, TOOL_RADIUS)
    _subtract_operation(stock, _circle(0.0, 0.0, 1.5), TOOL_RADIUS)
    clears_the_core = _line(-0.1, 0.0, 0.1, 0.0)
    assert _operation_removal(stock, clears_the_core, TOOL_RADIUS) is SweepRemoval.REMOVES
    _subtract_operation(stock, clears_the_core, TOOL_RADIUS)
    assert _operation_removal(stock, _line(-0.1, 0.0, 0.1, 0.0), TOOL_RADIUS) is SweepRemoval.NOTHING


def test_a_non_pythagorean_recut_is_within_band() -> None:
    stock = _rect()
    oblique = _line(2.0, 2.0, 8.0, 8.0)
    _subtract_operation(stock, oblique, TOOL_RADIUS)
    assert _operation_removal(stock, oblique, TOOL_RADIUS) is SweepRemoval.WITHIN_BAND


def test_an_arc_verdict_reads_its_chain_and_its_authored_annulus() -> None:
    stock = _rect()
    arc = ToolpathOperation(Arc(2.0, 0.0, math.pi / 2, frame=Frame([5.0, 5.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])), OperationType.CUT, 0, False, np.zeros(3), np.zeros(3))
    assert _operation_removal(stock, arc, TOOL_RADIUS) is SweepRemoval.REMOVES
    far = ToolpathOperation(Arc(2.0, 0.0, math.pi / 2, frame=Frame([30.0, 30.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])), OperationType.CUT, 0, False, np.zeros(3), np.zeros(3))
    assert _operation_removal(stock, far, TOOL_RADIUS) is SweepRemoval.NOTHING
    _subtract_operation(stock, arc, TOOL_RADIUS)
    assert _operation_removal(stock, arc, TOOL_RADIUS) is SweepRemoval.WITHIN_BAND


def test_an_unsupported_geometry_is_refused_by_name() -> None:
    bad = ToolpathOperation(Polygon([[0, 0, 0], [1, 0, 0], [1, 1, 0]]), OperationType.CUT, 0, False, np.zeros(3), np.zeros(3))
    with pytest.raises(UnexpectedToolpathGeometryError):
        _operation_removal(_rect(), bad, TOOL_RADIUS)
