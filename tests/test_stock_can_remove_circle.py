"""Contract of `Stock.can_remove_circle`: a removal it authorises never loses coverage.

When the test says True, replaying the neighbours and the preserved connectors
WITHOUT the current loop must leave exactly the stock that replaying them WITH
it leaves. The connectors are swept with the replay's own line sweep, which the
test itself uses, so "empty lost region => removable" keeps its direction.
"""

import numpy as np
from compas.geometry import Polygon
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.stock import Stock

TOOL_RADIUS = 1.0
LATTICE = st.integers(0, 40).map(lambda q: q / 4.0)
GUIDES = st.sampled_from((0.0, 0.5, 1.0, 1.5, 2.0, 3.0))


def _rect() -> Stock:
    return Stock(Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]))


def _replay(stock: Stock, circles: list[tuple[float, float, float]], connector: np.ndarray) -> Stock:
    result = stock.clone()
    for cx, cy, rho in circles:
        result.subtract_circle_sweep(cx, cy, rho, TOOL_RADIUS)
    for (x0, y0), (x1, y1) in zip(connector[:-1], connector[1:], strict=True):
        result.subtract_line_sweep(float(x0), float(y0), float(x1), float(y1), TOOL_RADIUS)
    return result


@st.composite
def _cases(draw: st.DrawFn) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float], np.ndarray]:
    def circle() -> tuple[float, float, float]:
        return (draw(LATTICE), draw(LATTICE), draw(GUIDES))

    previous, current, following = circle(), circle(), circle()
    if draw(st.booleans()):
        current = (previous[0], previous[1], previous[2])
    points = draw(st.lists(st.tuples(LATTICE, LATTICE), min_size=2, max_size=4))
    return previous, current, following, np.array(points, dtype=np.float64)


@settings(max_examples=60, deadline=None)
@given(case=_cases())
def test_an_authorised_removal_loses_no_coverage(case: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float], np.ndarray]) -> None:
    previous, current, following, connector = case
    stock = _rect()
    reference = stock.clone()
    if stock.can_remove_circle(previous, current, following, connector, TOOL_RADIUS):
        without = _replay(stock, [previous, following], connector)
        with_it = _replay(stock, [previous, current, following], connector)
        assert without.exactly_equals(with_it)
    assert stock.exactly_equals(reference)


def test_an_identical_neighbour_authorises_removal_through_the_wrapper() -> None:
    stock = _rect()
    assert stock.can_remove_circle((5.0, 5.0, 2.0), (5.0, 5.0, 2.0), (1.0, 1.0, 0.5), np.zeros((2, 2)), TOOL_RADIUS)


def test_a_unique_loop_is_retained_through_the_wrapper() -> None:
    stock = _rect()
    assert not stock.can_remove_circle((2.0, 5.0, 1.0), (5.0, 5.0, 2.0), (8.0, 5.0, 1.0), np.array([[2.0, 5.0], [8.0, 5.0]]), TOOL_RADIUS)
