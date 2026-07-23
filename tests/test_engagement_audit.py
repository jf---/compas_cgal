"""Stock wrapper and toolpath engagement-audit tests (Task 6).

Task 6a covers the typed `Stock` wrapper roundtrip; the engagement-audit tests
are appended by Task 6b once `engagement.py` lands.
"""

from compas.geometry import Polygon

from compas_cgal.stock import Stock

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


def test_stock_wrapper_roundtrip():
    stock = Stock(SQUARE)
    assert stock.contains(5.0, 5.0)
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not stock.contains(5.0, 5.0)
