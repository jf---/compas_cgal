"""A Stock2 clone must never reach geometry through its parent's freed traits.

CGAL's `Gps_on_surface_base_2` copy constructor allocates a fresh `Traits_2` for
the copy but builds the copy's arrangement with `Aos_2(*(ps.m_arr))`, and
`Arrangement_on_surface_2::assign` propagates a BORROWED traits pointer verbatim
(`m_geom_traits = arr.m_own_traits ? new Traits_adaptor_2 : arr.m_geom_traits`).
Every `Gps` arrangement borrows, so a clone's arrangement points at the ROOT
`Gps`'s traits object -- which the root's destructor deletes.

`Stock2::contains` is the one method that reads it: it builds an
`Arr_trapezoid_ric_point_location`, whose constructor copy-constructs a
`Td_traits` from `arr.geometry_traits()` and copies
`Arr_circle_segment_traits_2::inter_map` (a `std::map`) out of that memory. Every
other method reaches geometry through the `Gps`'s own traits and is unaffected.

The failure is a use-after-free, so it is NOT reliably a crash: when the freed
32 bytes still parse as a plausible `std::map`, `contains` returns a containment
answer computed against a corrupted intersection cache with nothing reported.
These tests therefore assert the ANSWERS as well as clean termination, and run
in a child process so a hard fault is one failed test rather than a dead
pytest worker.
"""

import subprocess
import sys

# 10 x 10 square with a 4..6 square island, so point location has to separate a
# contained face from an enclosed void rather than answer from a single face.
_STOCK_SETUP = """
import numpy as np
from compas_cgal import _stock_2

SQUARE = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=np.float64)
ISLAND = np.array([[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]], dtype=np.float64)


def probe(stock):
    return (
        stock.contains(2.0, 2.0),
        stock.contains(5.0, 5.0),
        stock.contains(-1.0, 5.0),
        stock.contains(11.0, 5.0),
    )


EXPECTED = (True, False, False, False)
"""

CLONE_FROM_TEMPORARY = (
    _STOCK_SETUP
    + """
# The parent is a temporary: it is destroyed on this line, before any query.
survivor = _stock_2.Stock2(SQUARE, [ISLAND]).clone()
assert probe(survivor) == EXPECTED, probe(survivor)
"""
)

CLONE_FAMILY_TWO_DEEP = (
    _STOCK_SETUP
    + """
# The shipped path clones twice (Stock2Area.__init__ and Stock2Area.raw), so the
# retained snapshot outlives BOTH of its ancestors.
root = _stock_2.Stock2(SQUARE, [ISLAND])
child = root.clone()
grandchild = child.clone()
del child
del root
assert probe(grandchild) == EXPECTED, probe(grandchild)
"""
)

CLONE_UNDER_ALLOCATION_CHURN = (
    _STOCK_SETUP
    + """
# Freed traits are only detected when something reuses the block, so churn the
# allocator between the parent's death and the first query, and repeat: a clone
# family that reads freed memory fails on some iteration even if the first one
# happens to find the bytes intact.
for _ in range(16):
    survivor = _stock_2.Stock2(SQUARE, [ISLAND]).clone()
    ballast = [bytearray(32) for _ in range(4096)]
    assert probe(survivor) == EXPECTED, probe(survivor)
    del ballast, survivor
"""
)


def _run_isolated(code: str) -> None:
    """Run `code` in a fresh interpreter; a fault becomes a named assertion failure."""
    finished = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    if finished.returncode < 0:
        raise AssertionError(f"child process died on signal {-finished.returncode} (use-after-free on the borrowed Gps traits)\n{finished.stderr}")
    assert finished.returncode == 0, finished.stderr


def test_clone_of_a_destroyed_parent_locates_points() -> None:
    _run_isolated(CLONE_FROM_TEMPORARY)


def test_clone_family_two_deep_locates_points() -> None:
    _run_isolated(CLONE_FAMILY_TWO_DEEP)


def test_clone_locates_points_under_allocation_churn() -> None:
    _run_isolated(CLONE_UNDER_ALLOCATION_CHURN)
