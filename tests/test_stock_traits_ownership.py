"""Every `Stock2` must keep alive the object owning the traits its arrangement reads.

CGAL's `Gps_on_surface_base_2` copy constructor allocates a fresh `Traits_2` for
the copy but builds the copy's arrangement as `Aos_2(*(ps.m_arr))`, and
`Arrangement_on_surface_2::assign` propagates a BORROWED traits pointer verbatim
(`m_geom_traits = arr.m_own_traits ? new Traits_adaptor_2 : arr.m_geom_traits`).
Every `Gps` arrangement borrows, so a copied set reads the traits of the set it
was copied from and its own freshly allocated traits stays unused -- until the
first two-operand boolean operation, which rebuilds the arrangement on the copy's
OWN traits (`_difference(const Aos_2&)` does `new Aos_2(m_traits)`).

That second half is why a `shared_ptr<const GpsTraits>` family member cannot fix
this: CGAL reassigns the pointer the arrangement reads on the first operation, so
one traits object cannot stand for a whole clone family. The owner has to be
resolved from the object graph wherever the arrangement changes.

These tests assert the OWNERSHIP INVARIANT, not a symptom. The four shapes below
are undefined behaviour whether or not the allocator happens to be kind:

* an operated clone cloned again (`Stock2::clone`), which SIGSEGVs in
  `Stock2::contains` -- `Arr_trapezoid_ric_point_location`'s constructor
  copy-constructs a `Td_traits` out of the arrangement's traits;
* a clone taken after `replace_set` (the `subtract_exact_*` depletion path),
  where even the root has left the family object;
* `replace_set` onto an ALREADY EMPTY stock, where `_difference` early-returns
  and the replacement is left reading traits `replace_set` then frees -- no
  clone involved, and it faults about one run in twelve;
* the sweep oracle's two sets, latent today because its consumers reach geometry
  only through the overlay result and `is_empty()`.

Ancestors are deliberately kept alive here: freeing them is what turns the first
shape into a hard fault, and a dead pytest worker is worse evidence than a failed
assertion. `tests/test_stock_clone_lifetime.py` runs the fault arm in child
processes.
"""

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

# 10 x 10 square with a 4..6 square island, so point location has to separate a
# contained face from an enclosed void rather than answer from a single face.
SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)
ISLAND = np.array(
    [[4.0, 4.0, 0.0], [6.0, 4.0, 0.0], [6.0, 6.0, 0.0], [4.0, 6.0, 0.0]],
    dtype=np.float64,
)
# Inside the ring, inside the island void, outside on two sides, and on a corner.
# The corner answers False: `Stock2::contains` locates a FACE and asks whether it
# is contained, so a point ON the boundary is located as a vertex and is not
# material -- unlike `ExactRegion2::contains`, which uses `oriented_side`.
PROBES = ((2.0, 2.0), (5.0, 5.0), (-1.0, 5.0), (11.0, 5.0), (0.0, 0.0))
UNCUT = (True, False, False, False, False)

# A small square a single exact segment depletion can empty, so the
# replace-set-onto-empty shape is reachable in two calls.
SMALL_SQUARE = np.array(
    [[-2.0, -2.0, 0.0], [2.0, -2.0, 0.0], [2.0, 2.0, 0.0], [-2.0, 2.0, 0.0]],
    dtype=np.float64,
)
FLOODING_SEGMENT = (-1.0, 0.0, 1.0, 0.0)
FLOODING_TOOL_RADIUS = 5.0
COARSE_CHORD = 10.0
CENTER_COUNT_LIMIT = 4096


def owns_traits(stock: _stock_2.Stock2) -> bool:
    return stock.arrangement_traits_are_owned_for_audit()


def probe(stock: _stock_2.Stock2) -> tuple:
    return tuple(stock.contains(x, y) for x, y in PROBES)


def raw_stock() -> _stock_2.Stock2:
    return _stock_2.Stock2(SQUARE, [ISLAND])


def flooded_stock() -> _stock_2.Stock2:
    """A stock emptied through `replace_set`, so its set owns its own traits."""
    stock = _stock_2.Stock2(SMALL_SQUARE, [])
    stock.subtract_exact_segment(
        *FLOODING_SEGMENT,
        FLOODING_TOOL_RADIUS,
        COARSE_CHORD,
        CENTER_COUNT_LIMIT,
    )
    return stock


def test_root_owns_its_traits() -> None:
    root = raw_stock()
    assert owns_traits(root)
    assert probe(root) == UNCUT


def test_unoperated_clone_owns_its_traits() -> None:
    root = raw_stock()
    clone = root.clone()
    assert owns_traits(clone)
    assert probe(clone) == UNCUT


def test_operated_stock_owns_its_traits() -> None:
    """The operated clone itself is fine -- the operation moved it onto its own."""
    root = raw_stock()
    child = root.clone()
    child.subtract_disk(3.0, 3.0, 1.0)
    assert owns_traits(child)


def test_clone_of_an_operated_clone_owns_its_traits() -> None:
    root = raw_stock()
    child = root.clone()
    child.subtract_disk(3.0, 3.0, 1.0)
    grandchild = child.clone()
    assert owns_traits(grandchild)
    assert grandchild.contains(8.0, 8.0)
    assert not grandchild.contains(3.0, 3.0)


def test_clone_after_replace_set_owns_its_traits() -> None:
    """`subtract_exact_*` installs a trial that owns its traits; the root leaves the family."""
    stock = _stock_2.Stock2(SQUARE, [ISLAND])
    stock.subtract_exact_segment(1.0, 1.0, 3.0, 1.0, 0.5, 0.5, CENTER_COUNT_LIMIT)
    assert owns_traits(stock)
    clone = stock.clone()
    assert owns_traits(clone)
    assert clone.contains(8.0, 8.0)


def test_replace_set_onto_an_empty_stock_owns_its_traits() -> None:
    """No clone involved: `_difference` early-returns on an empty set, so the
    replacement never rebuilds and is left reading traits `replace_set` frees."""
    stock = flooded_stock()
    assert stock.is_empty()
    assert owns_traits(stock)
    stock.subtract_exact_segment(
        *FLOODING_SEGMENT,
        FLOODING_TOOL_RADIUS,
        COARSE_CHORD,
        CENTER_COUNT_LIMIT,
    )
    assert stock.is_empty()
    assert owns_traits(stock)


def test_segment_sweep_oracle_sets_own_their_traits() -> None:
    assert _stock_2.exact_segment_sweep_oracle_traits_are_owned_for_audit(
        -3.0,
        -1.0,
        3.0,
        7.0,
        10.0,
        0.75,
        0.5,
        CENTER_COUNT_LIMIT,
    )


def test_full_circle_sweep_oracle_sets_own_their_traits() -> None:
    assert _stock_2.exact_full_circle_sweep_oracle_traits_are_owned_for_audit(
        0.0,
        0.0,
        5.0,
        0.0,
        5.0,
        0.75,
        0.5,
        CENTER_COUNT_LIMIT,
    )


@pytest.fixture
def public_stock() -> Stock:
    return Stock(
        Polygon([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]),
        [Polygon([[4.0, 4.0, 0.0], [6.0, 4.0, 0.0], [6.0, 6.0, 0.0], [4.0, 6.0, 0.0]])],
    )


def test_public_stock_root_and_clone_own_their_traits(public_stock: Stock) -> None:
    assert owns_traits(public_stock.raw)
    clone = public_stock.clone()
    assert owns_traits(clone.raw)


def test_public_stock_clone_of_an_operated_clone_owns_its_traits(public_stock: Stock) -> None:
    """`Stock.clone()` is the shipped surface: two clones with any depletion between them."""
    child = public_stock.clone()
    child.subtract_disk(3.0, 3.0, 1.0)
    grandchild = child.clone()
    assert owns_traits(grandchild.raw)
    assert grandchild.contains(8.0, 8.0)
    assert not grandchild.contains(3.0, 3.0)
