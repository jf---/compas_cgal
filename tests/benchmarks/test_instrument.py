from __future__ import annotations

from benchmarks.families.analytic import rectangle
from benchmarks.instrument import probe_digits, probe_size
from compas_cgal.stock import Stock


def _rect_stock() -> Stock:
    spec = rectangle(width=20.0, height=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    return Stock(spec.polygon, [])


def test_arrangement_size_grows_with_depletion() -> None:
    stock = _rect_stock()
    before = probe_size(stock)
    assert before.vertices > 0
    assert before.faces > 0
    for i in range(6):
        stock.subtract_disk(-6.0 + 2.0 * i, 0.0, 0.5)
    after = probe_size(stock)
    assert after.vertices > before.vertices
    assert after.halfedges > before.halfedges


def test_coordinate_digits_are_reported_and_nonzero() -> None:
    stock = _rect_stock()
    stock.subtract_disk(0.0, 0.0, 0.5)
    digits = probe_digits(stock)
    assert digits.sampled > 0
    assert digits.max_digits >= digits.mean_digits > 0.0


def test_coordinate_digits_grow_under_chained_subtraction() -> None:
    stock = _rect_stock()
    stock.subtract_disk(0.0, 0.0, 0.5)
    first = probe_digits(stock).max_digits
    for i in range(10):
        stock.subtract_capsule(-5.0 + i * 1.0, 0.1 * i, -4.0 + i * 1.0, 0.1 * i + 0.3, 0.5)
    later = probe_digits(stock).max_digits
    assert later >= first
