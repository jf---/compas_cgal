from __future__ import annotations

from benchmarks.families.analytic import rectangle
from benchmarks.instrument import probe_digits
from benchmarks.instrument import probe_size
from compas_cgal.stock import Stock

TOOL_RADIUS = 0.5

# Six disks along the rectangle's spine, spaced two radii apart so each removal
# meets the previous one and the boolean chain actually compounds.
SPINE_DISKS = 6
SPINE_START_X = -6.0
SPINE_STEP_X = 2.0

# A staircase of overlapping capsules walking across the pocket: each one starts
# where the previous ended and drifts in y, so every removal cuts into coordinates
# the previous removal constructed. That chaining is what grows the rationals.
CHAIN_LINKS = 10
CHAIN_START_X = -5.0
CHAIN_STEP_X = 1.0
CHAIN_LENGTH_X = 1.0
CHAIN_DRIFT_Y = 0.1
CHAIN_RISE_Y = 0.3


def _rect_stock() -> Stock:
    """A virgin rectangular stock, large enough to absorb the removals below."""
    spec = rectangle(width=20.0, height=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    return Stock(spec.polygon, list(spec.holes))


def test_arrangement_size_grows_with_depletion() -> None:
    stock = _rect_stock()
    before = probe_size(stock)
    assert before.vertices > 0
    assert before.faces > 0
    for i in range(SPINE_DISKS):
        stock.subtract_disk(SPINE_START_X + SPINE_STEP_X * i, 0.0, TOOL_RADIUS)
    after = probe_size(stock)
    assert after.vertices > before.vertices
    assert after.halfedges > before.halfedges


def test_coordinate_digits_are_reported_and_nonzero() -> None:
    stock = _rect_stock()
    stock.subtract_disk(0.0, 0.0, TOOL_RADIUS)
    digits = probe_digits(stock)
    assert digits.sampled > 0
    assert digits.max_digits >= digits.mean_digits > 0.0


def test_coordinate_digits_grow_under_chained_subtraction() -> None:
    stock = _rect_stock()
    stock.subtract_disk(0.0, 0.0, TOOL_RADIUS)
    first = probe_digits(stock).max_digits
    for i in range(CHAIN_LINKS):
        x0 = CHAIN_START_X + i * CHAIN_STEP_X
        y0 = CHAIN_DRIFT_Y * i
        stock.subtract_capsule(x0, y0, x0 + CHAIN_LENGTH_X, y0 + CHAIN_RISE_Y, TOOL_RADIUS)
    later = probe_digits(stock).max_digits
    assert later >= first
