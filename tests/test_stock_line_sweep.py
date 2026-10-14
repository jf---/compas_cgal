"""Exact contract of the replay's line sweep (`Stock.subtract_line_sweep`).

One model, picked by an EXACT test on the injected endpoints:

* ``dx^2 + dy^2`` a square in Q (axis-aligned, 3-4-5, ...): the swept capsule
  itself, equality, because its side lines are then rational.
* otherwise: the canonical-direction QUAD, the exact full-radius end disks
  joined with a rectangle of half-width ``h`` whose band ``(1 - f) r <= h <= r``
  is checked exactly. The direction class fixes ``h``, so collinear cuts share
  one rectangle width; only the side band is inexact.

SAFETY and COVERAGE are decided per sample over `fractions.Fraction`, exactly as
for the quad (tests/test_stock_capsule_quad.py), whose oracle and sampler are
reused verbatim.
"""

from fractions import Fraction

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from tests.test_stock_capsule_quad import MIN_DISCRIMINATING_SAMPLES
from tests.test_stock_capsule_quad import SEGMENTS
from tests.test_stock_capsule_quad import SQUARE
from tests.test_stock_capsule_quad import exact_squared_distance_to_segment
from tests.test_stock_capsule_quad import sample_points
from tests.test_stock_capsule_quad import slack_half_width

PYTHAGOREAN_SEGMENTS = [
    (2.0, 5.0, 8.0, 5.0, 0.5),
    (5.0, 8.0, 5.0, 2.0, 0.75),
    (2.0, 2.0, 8.0, 10.0, 1.0),  # 6-8-10
    (1.0, 1.0, 4.0, 5.0, 0.6),  # 3-4-5
    (9.0, 1.0, 4.0, 13.0, 0.6),  # 5-12-13, leaves the square
    (5.0, 5.0, 5.75, 6.0, 0.4),  # 3/4-1-5/4
]


def _stock() -> Stock:
    return Stock(Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]))


def _model(segment: tuple[float, float, float, float, float]) -> "_stock_2.LineSweepModel":
    return _stock_2.line_sweep_model(*segment)


def test_the_model_is_decided_by_the_exact_pythagorean_test() -> None:
    assert _model((2.0, 5.0, 8.0, 5.0, 0.5)) == _stock_2.LineSweepModel.EXACT_CAPSULE
    assert _model((1.0, 1.0, 4.0, 5.0, 0.6)) == _stock_2.LineSweepModel.EXACT_CAPSULE
    assert _model((5.0, 5.0, 5.75, 6.0, 0.4)) == _stock_2.LineSweepModel.EXACT_CAPSULE
    assert _model((2.0, 2.0, 8.0, 8.0, 1.0)) == _stock_2.LineSweepModel.CANONICAL_QUAD
    assert _model((1.0, 4.0, 9.0, 4.6, 0.6)) == _stock_2.LineSweepModel.CANONICAL_QUAD
    assert _model((3.0, 3.0, 3.0, 3.0, 0.6)) == _stock_2.LineSweepModel.DISK


@pytest.mark.parametrize("segment", SEGMENTS + PYTHAGOREAN_SEGMENTS)
def test_line_sweep_removes_nothing_outside_the_true_capsule(segment: tuple[float, float, float, float, float]) -> None:
    x0, y0, x1, y1, radius = segment
    before = _stock_2.Stock2(SQUARE, [])
    after = _stock_2.Stock2(SQUARE, [])
    after.subtract_line_sweep(x0, y0, x1, y1, radius)
    radius_squared = Fraction(radius) ** 2
    removed = 0
    violations = []
    for px, py in sample_points(*segment):
        if after.contains(px, py) or not before.contains(px, py):
            continue
        removed += 1
        if exact_squared_distance_to_segment(px, py, x0, y0, x1, y1) > radius_squared:
            violations.append((px, py))
    assert violations == []
    assert removed >= MIN_DISCRIMINATING_SAMPLES


@pytest.mark.parametrize("segment", SEGMENTS + PYTHAGOREAN_SEGMENTS)
def test_line_sweep_covers_the_documented_slack_band(segment: tuple[float, float, float, float, float]) -> None:
    x0, y0, x1, y1, radius = segment
    before = _stock_2.Stock2(SQUARE, [])
    after = _stock_2.Stock2(SQUARE, [])
    after.subtract_line_sweep(x0, y0, x1, y1, radius)
    band_squared = slack_half_width(radius) ** 2
    in_band = 0
    survivors = []
    for px, py in sample_points(*segment):
        if not before.contains(px, py):
            continue
        if exact_squared_distance_to_segment(px, py, x0, y0, x1, y1) > band_squared:
            continue
        in_band += 1
        if after.contains(px, py):
            survivors.append((px, py))
    assert survivors == []
    assert in_band >= MIN_DISCRIMINATING_SAMPLES


@pytest.mark.parametrize("segment", PYTHAGOREAN_SEGMENTS)
def test_a_pythagorean_sweep_removes_everything_strictly_inside_the_capsule(segment: tuple[float, float, float, float, float]) -> None:
    """EQUALITY, not a band: every sample strictly within r of the segment is gone."""
    x0, y0, x1, y1, radius = segment
    before = _stock_2.Stock2(SQUARE, [])
    after = _stock_2.Stock2(SQUARE, [])
    after.subtract_line_sweep(x0, y0, x1, y1, radius)
    radius_squared = Fraction(radius) ** 2
    inside = [(px, py) for px, py in sample_points(*segment) if before.contains(px, py) and exact_squared_distance_to_segment(px, py, x0, y0, x1, y1) < radius_squared]
    assert [point for point in inside if after.contains(*point)] == []
    assert len(inside) >= MIN_DISCRIMINATING_SAMPLES


@pytest.mark.parametrize("segment", SEGMENTS + PYTHAGOREAN_SEGMENTS)
def test_line_sweep_keeps_the_boolean_set_canonical(segment: tuple[float, float, float, float, float]) -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_line_sweep(*segment)
    assert stock.representation_is_valid()


@pytest.mark.parametrize("segment", SEGMENTS + PYTHAGOREAN_SEGMENTS)
def test_the_removal_verdict_on_virgin_stock_is_removes(segment: tuple[float, float, float, float, float]) -> None:
    stock = _stock()
    assert stock.line_sweep_removal(*segment) is SweepRemoval.REMOVES


@pytest.mark.parametrize("segment", PYTHAGOREAN_SEGMENTS)
def test_an_exact_sweep_already_cut_reads_nothing(segment: tuple[float, float, float, float, float]) -> None:
    stock = _stock()
    stock.subtract_line_sweep(*segment)
    assert stock.line_sweep_removal(*segment) is SweepRemoval.NOTHING


@pytest.mark.parametrize("segment", [s for s in SEGMENTS if _stock_2.line_sweep_model(*s) == _stock_2.LineSweepModel.CANONICAL_QUAD])
def test_a_quad_already_cut_reads_within_band_never_nothing(segment: tuple[float, float, float, float, float]) -> None:
    """The model left the side band ``h < d < r`` standing; the tool did not.

    So the verdict cannot be NOTHING (the true capsule still meets material),
    and it is not REMOVES (the model would change nothing): WITHIN_BAND.
    """
    stock = _stock()
    stock.subtract_line_sweep(*segment)
    assert stock.line_sweep_removal(*segment) is SweepRemoval.WITHIN_BAND


# Side offsets of a thin strip of stock parallel to the oblique segment
# (2, 2) -> (8, 8), r = 1, as multiples of r. The OVER rectangle is certified in
# [r, (1 + f) r] and aimed at (1 + f/2) r, f = 1e-4; UNDER in [(1 - f) r, r].
BAND_STRIPS = [
    pytest.param((1.0 - 3e-4, 1.0 - 2e-4), SweepRemoval.REMOVES, id="inside-under"),
    pytest.param((1.0 - 1e-6, 1.0 + 1e-6), SweepRemoval.WITHIN_BAND, id="straddles-r"),
    pytest.param((1.0 + 1e-6, 1.0 + 2e-6), SweepRemoval.WITHIN_BAND, id="inside-over-beyond-r"),
    pytest.param((1.0 + 3e-4, 1.0 + 4e-4), SweepRemoval.NOTHING, id="beyond-certified-over"),
]


@pytest.mark.parametrize(("offsets", "expected"), BAND_STRIPS)
def test_the_over_band_edges_are_pinned(offsets: tuple[float, float], expected: SweepRemoval) -> None:
    """Material only at a chosen distance from the axis of a canonical quad."""
    root_half = 0.5**0.5
    ux, uy, nx, ny = root_half, root_half, -root_half, root_half
    mx, my = 5.0, 5.0
    near, far = offsets
    corners = [
        [mx - 0.5 * ux + near * nx, my - 0.5 * uy + near * ny, 0.0],
        [mx + 0.5 * ux + near * nx, my + 0.5 * uy + near * ny, 0.0],
        [mx + 0.5 * ux + far * nx, my + 0.5 * uy + far * ny, 0.0],
        [mx - 0.5 * ux + far * nx, my - 0.5 * uy + far * ny, 0.0],
    ]
    # The native constructor, not compas `Polygon`: the strip is 1e-6 wide, below
    # compas's point-merging tolerance.
    strip = Stock._from_raw(_stock_2.Stock2(np.array(corners, dtype=np.float64), []))
    assert _stock_2.line_sweep_model(2.0, 2.0, 8.0, 8.0, 1.0) == _stock_2.LineSweepModel.CANONICAL_QUAD
    assert strip.line_sweep_removal(2.0, 2.0, 8.0, 8.0, 1.0) is expected


def test_a_sweep_far_from_material_reads_nothing() -> None:
    stock = _stock()
    assert stock.line_sweep_removal(20.0, 20.0, 27.0, 23.0, 1.0) is SweepRemoval.NOTHING
    assert stock.line_sweep_removal(20.0, 20.0, 26.0, 28.0, 1.0) is SweepRemoval.NOTHING


def test_the_verdict_never_mutates_the_stock() -> None:
    stock = _stock()
    reference = stock.clone()
    stock.line_sweep_removal(2.0, 2.0, 8.0, 8.0, 1.0)
    stock.circle_sweep_removal(5.0, 5.0, 2.0, 1.0)
    assert stock.exactly_equals(reference)


def test_a_circle_sweep_verdict_is_exact_both_ways() -> None:
    stock = _stock()
    assert stock.circle_sweep_removal(5.0, 5.0, 2.0, 1.0) is SweepRemoval.REMOVES
    stock.subtract_circle_sweep(5.0, 5.0, 2.0, 1.0)
    assert stock.circle_sweep_removal(5.0, 5.0, 2.0, 1.0) is SweepRemoval.NOTHING
    assert stock.circle_sweep_removal(5.0, 5.0, 2.5, 0.5) is SweepRemoval.NOTHING
    assert stock.circle_sweep_removal(5.0, 5.0, 0.5, 1.0) is SweepRemoval.REMOVES


@pytest.mark.parametrize(
    ("segment", "error"),
    [
        ((0.0, 0.0, float("nan"), 1.0, 1.0), _stock_2.NonFiniteCapsuleInputError),
        ((0.0, 0.0, 1.0, 1.0, float("inf")), _stock_2.NonFiniteCapsuleInputError),
        ((0.0, 0.0, 1.0, 1.0, 0.0), _stock_2.InvalidLineSweepRadiusError),
        ((0.0, 0.0, 1.0, 1.0, -1.0), _stock_2.InvalidLineSweepRadiusError),
    ],
)
def test_malformed_line_sweeps_are_refused_by_name(segment: tuple[float, float, float, float, float], error: type[Exception]) -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    for method in (stock.subtract_line_sweep, stock.subtract_line_sweep_local, stock.line_sweep_removal):
        with pytest.raises(error):
            method(*segment)
    assert stock.exactly_equals(_stock_2.Stock2(SQUARE, []))


def test_the_line_sweep_errors_are_value_errors() -> None:
    assert issubclass(_stock_2.InvalidLineSweepRadiusError, ValueError)
    assert issubclass(_stock_2.NonFiniteCapsuleInputError, ValueError)
