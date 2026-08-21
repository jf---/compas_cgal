"""Equivalence gate for the LOCAL depletion path.

A depleted stock is the input to every later engagement decision, so a local
update whose face containment is subtly wrong produces a stock that looks fine
and silently corrupts every downstream verdict. The gate is therefore exact
point-set EQUALITY against the global `General_polygon_set_2::difference` path --
never containment, never a subset relation -- plus the boolean engine's own
representation invariant, which point-set equality alone cannot see.
"""

import math

import numpy as np
import pytest

from compas_cgal import _stock_2

SQUARE = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=np.float64)
ISLAND = np.array([[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]], dtype=np.float64)
L_SHAPE = np.array(
    [[0, 0, 0], [12, 0, 0], [12, 5, 0], [5, 5, 0], [5, 10, 0], [0, 10, 0]],
    dtype=np.float64,
)

# ("disk", cx, cy, radius) or ("annulus", cx, cy, inner, outer).
#
# 26 depletions covering, in order: disjoint interior disks; overlapping disks;
# a disk repeated exactly (idempotent removal); disks straddling the outer
# boundary on three sides and one corner; a disk landing entirely OUTSIDE the
# stock; annuli fully interior, overlapping each other, degenerate (inner == 0,
# a plain disk), straddling the boundary, concentric-nested, one whose hole
# encloses previously removed material, and one whose ring exactly re-traces an
# earlier annulus radius (a coincident boundary circle, the case where the
# "edge lies on the region boundary" test has to be exact).
MIXED_DEPLETIONS = [
    ("disk", 2.0, 2.0, 0.75),
    ("disk", 8.0, 2.0, 0.75),
    ("disk", 2.5, 2.5, 0.75),
    ("disk", 2.5, 2.5, 0.75),
    ("disk", 3.0, 2.0, 1.25),
    ("disk", 0.0, 5.0, 1.5),
    ("disk", 10.0, 7.0, 1.5),
    ("disk", 5.0, 0.0, 1.25),
    ("disk", 0.0, 0.0, 2.0),
    ("disk", 15.0, 15.0, 1.0),
    ("disk", 9.5, 9.5, 0.4),
    ("annulus", 5.0, 5.0, 1.0, 2.0),
    ("annulus", 5.0, 5.0, 2.5, 3.0),
    ("annulus", 5.0, 5.0, 0.0, 0.5),
    ("annulus", 6.0, 5.5, 1.0, 2.0),
    ("annulus", 3.0, 7.0, 0.5, 1.5),
    ("annulus", 3.6, 7.4, 0.5, 1.5),
    ("annulus", 5.0, 5.0, 1.0, 2.0),
    ("annulus", 0.0, 10.0, 1.0, 2.5),
    ("annulus", 10.0, 0.0, 0.75, 2.25),
    ("annulus", 7.5, 3.0, 0.25, 0.75),
    ("annulus", 7.5, 3.0, 1.25, 1.75),
    ("annulus", 2.0, 2.0, 2.0, 3.0),
    ("annulus", 14.0, 14.0, 0.5, 1.0),
    ("disk", 5.0, 9.0, 1.1),
    ("annulus", 5.0, 9.0, 1.1, 2.2),
]

assert len(MIXED_DEPLETIONS) >= 20


def apply_global(stock, op):
    """Remove one operation through the global `Gps::difference` path."""
    if op[0] == "disk":
        stock.subtract_disk(op[1], op[2], op[3])
    else:
        stock.subtract_annulus(op[1], op[2], op[3], op[4])


def apply_local(stock, op):
    """Remove one operation through the local arrangement-surgery path."""
    if op[0] == "disk":
        stock.subtract_disk_local(op[1], op[2], op[3])
    else:
        stock.subtract_annulus_local(op[1], op[2], op[3], op[4])


@pytest.mark.parametrize(
    "boundary,holes",
    [(SQUARE, []), (SQUARE, [ISLAND]), (L_SHAPE, [])],
    ids=["square", "square_with_island", "L_shape"],
)
def test_local_depletion_matches_global_step_by_step(boundary, holes):
    """After EVERY one of the 26 mixed depletions the two stocks are equal.

    Checked after each step rather than once at the end, so a divergence names
    the operation that introduced it instead of surfacing many removals later.
    """
    reference = _stock_2.Stock2(boundary, holes)
    local = _stock_2.Stock2(boundary, holes)
    for index, op in enumerate(MIXED_DEPLETIONS):
        apply_global(reference, op)
        apply_local(local, op)
        assert local.exactly_equals(reference), f"step {index} diverged: {op}"
        assert reference.exactly_equals(local), f"step {index} diverged: {op}"
        assert local.representation_is_valid(), f"step {index} broke the Gps invariant: {op}"


def test_local_depletion_matches_global_in_reverse_order():
    """The same 26 removals applied back to front still agree exactly.

    Different arrival order builds a different intermediate arrangement, so this
    exercises the same removals against face layouts the forward run never
    produces.
    """
    reference = _stock_2.Stock2(SQUARE, [ISLAND])
    local = _stock_2.Stock2(SQUARE, [ISLAND])
    for index, op in enumerate(reversed(MIXED_DEPLETIONS)):
        apply_global(reference, op)
        apply_local(local, op)
        assert local.exactly_equals(reference), f"reverse step {index} diverged: {op}"
        assert local.representation_is_valid(), f"reverse step {index} broke the Gps invariant: {op}"


def test_local_depletion_matches_global_on_a_random_sequence():
    """40 pseudo-random disks and annuli, fixed seed, exact equality throughout.

    Random placement reaches tangency, near-tangency and containment
    configurations a hand-written list does not enumerate; the seed is fixed so a
    failure is reproducible.
    """
    rng = np.random.default_rng(20260820)
    reference = _stock_2.Stock2(SQUARE, [ISLAND])
    local = _stock_2.Stock2(SQUARE, [ISLAND])
    for index in range(40):
        cx = float(rng.uniform(-1.0, 11.0))
        cy = float(rng.uniform(-1.0, 11.0))
        outer = float(rng.uniform(0.3, 2.5))
        if index % 2:
            op = ("disk", cx, cy, outer)
        else:
            op = ("annulus", cx, cy, float(rng.uniform(0.0, 0.9)) * outer, outer)
        apply_global(reference, op)
        apply_local(local, op)
        assert local.exactly_equals(reference), f"random step {index} diverged: {op}"
        assert local.representation_is_valid(), f"random step {index} broke the Gps invariant: {op}"


def test_local_arc_sweep_matches_global_full_turn():
    """A full-turn arc sweep removes the identical annulus on both paths.

    This is the call the generator actually makes, so the equality has to hold at
    the arc-sweep entry point and not only at the annulus core it dispatches to.
    """
    reference = _stock_2.Stock2(SQUARE, [])
    local = _stock_2.Stock2(SQUARE, [])
    for index, (cx, cy, guide, angle) in enumerate(
        [
            (5.0, 5.0, 2.0, 0.0),
            (3.0, 7.0, 1.0, 0.75),
            (8.0, 3.0, 0.4, 2.1),  # guide narrower than the tool: a filled disk
            (0.5, 0.5, 1.5, 4.0),  # straddles two boundary edges
            (5.0, 5.0, 2.0, 1.3),  # exactly retraces the first annulus
        ]
    ):
        ex = cx + guide * math.cos(angle)
        ey = cy + guide * math.sin(angle)
        reference.subtract_arc_sweep(cx, cy, ex, ey, ex, ey, True, 0.5)
        local.subtract_arc_sweep_local(cx, cy, ex, ey, ex, ey, True, 0.5)
        assert local.exactly_equals(reference), f"full turn {index} diverged"
        assert local.representation_is_valid(), f"full turn {index} broke the Gps invariant"


def test_local_arc_sweep_defers_partial_arcs_to_the_global_chain():
    """A partial arc is out of the local path's scope and must still be exact."""
    reference = _stock_2.Stock2(SQUARE, [])
    local = _stock_2.Stock2(SQUARE, [])
    reference.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 5.0, 7.0, False, 0.5)
    local.subtract_arc_sweep_local(5.0, 5.0, 7.0, 5.0, 5.0, 7.0, False, 0.5)
    assert local.exactly_equals(reference)


def test_local_depletion_on_an_already_empty_stock():
    """Removing from nothing leaves nothing, and leaves a valid representation."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk_local(5.0, 5.0, 20.0)
    assert stock.is_empty()
    stock.subtract_annulus_local(5.0, 5.0, 1.0, 2.0)
    assert stock.is_empty()
    assert stock.representation_is_valid()


def test_local_depletion_rejects_the_same_arguments_the_global_path_rejects():
    """The two paths share one argument contract, exception types included."""
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(ValueError):
        stock.subtract_disk_local(5.0, 5.0, 0.0)
    with pytest.raises(_stock_2.InvalidAnnulusRadiiError):
        stock.subtract_annulus_local(5.0, 5.0, 2.0, 1.0)
    with pytest.raises(_stock_2.InvalidAnnulusRadiiError):
        stock.subtract_annulus_local(5.0, 5.0, -1.0, 1.0)
    with pytest.raises(_stock_2.NonFiniteAnnulusInputError):
        stock.subtract_annulus_local(float("nan"), 5.0, 1.0, 2.0)


def test_local_depletion_keeps_the_arrangement_the_same_size_as_the_global_path():
    """Equal point sets are not enough: the arrangements must not diverge in size.

    A local update that leaves extra vertices or edges behind would still pass
    `exactly_equals` while making every later zone query and removal slower, one
    depletion at a time.
    """
    reference = _stock_2.Stock2(SQUARE, [ISLAND])
    local = _stock_2.Stock2(SQUARE, [ISLAND])
    for op in MIXED_DEPLETIONS:
        apply_global(reference, op)
        apply_local(local, op)
    assert local.arrangement_stats() == reference.arrangement_stats()
