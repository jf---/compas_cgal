"""Witness attribution of the Omega* coverage tool: exact nearest declared piece, deterministic ties.

Attribution is a reporting label, but the choice of piece is an exact decision:
every distance is a closed form ``|√a + k|`` over the declared rationals and the
minimum is taken with exact comparisons, ties going to the lowest primitive
index. No proximity window around junctions may hide a piece's own distance.

The tool's independent sweep oracle registers each sweep in every integer cell
its rational bounding box meets (``math.floor`` on Fractions), so a point's cell
lists every sweep that can contain it; membership is Fraction arithmetic. No
binary64 distance and no slack decides which sweeps are tested.
"""

from __future__ import annotations

import hashlib
import random
from fractions import Fraction

from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_one_root import OneRootNumber
from compas_cgal.units import WorldXY
from tools.coverage_certificate_omega_star import SweptOracle
from tools.coverage_certificate_omega_star import WitnessAttributor
from tools.coverage_certificate_omega_star import attribute_witness

F = Fraction
# A binary64 hairline offset from a wall, mm: 2^-46 ~ 1.4e-14 (the CC defect's depth scale).
HAIRLINE_OFFSET = F(1, 2**46)


def _point(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(F(x), F(y))


class _Pocket:
    """Sequential builder of a declared CCW cycle."""

    def __init__(self, start: ExactPoint2[WorldXY]) -> None:
        self.cursor = start
        self.primitives: list[DeclaredLine | DeclaredArc] = []

    def line(self, quarter: int, length: Fraction | int) -> _Pocket:
        line = DeclaredLine.build(self.cursor, PythagoreanDirection.build(F(0), quarter), F(length), len(self.primitives))
        self.primitives.append(line)
        self.cursor = line.end
        return self

    def arc(self, end: ExactPoint2[WorldXY], t: Fraction, quarter: int, sense: int) -> _Pocket:
        arc = DeclaredArc.build(self.cursor, end, PythagoreanDirection.build(t, quarter), 1 if sense == 1 else -1, len(self.primitives))
        self.primitives.append(arc)
        self.cursor = arc.end
        return self

    def close(self, name: str) -> DeclaredBoundary:
        return DeclaredBoundary.build(name, tuple(self.primitives), hashlib.sha256(name.encode()).hexdigest())


def _concave_arc() -> DeclaredBoundary:
    """20 x 20 square whose bottom edge is a CW arc about (10, -7.5), rho = 12.5; line 1 is x = 20 from (20, 0)."""
    return _Pocket(_point(0, 0)).arc(_point(20, 0), F(1, 2), 0, -1).line(1, 20).line(2, 20).line(3, 20).close("concave_arc")


def _l_shape() -> DeclaredBoundary:
    """L-shape with a reflex corner at (10, 10), the junction of lines 2 and 3."""
    return _Pocket(_point(0, 0)).line(0, 20).line(1, 10).line(2, 10).line(1, 10).line(2, 10).line(3, 20).close("l_shape")


def test_a_witness_near_a_line_start_is_attributed_to_its_own_line() -> None:
    """The CC defect: 0.01 mm from line 1's start, 2^-46 mm from line 1, ~6e-3 mm from the concave arc."""
    boundary = _concave_arc()
    witness = (F(20) - HAIRLINE_OFFSET, F(1, 100))
    attribution = attribute_witness(witness, boundary.primitives)
    assert attribution.primitive_index == 1
    assert attribution.feature == "line"
    assert attribution.squared_depth.compare(HAIRLINE_OFFSET * HAIRLINE_OFFSET) == 0


def test_a_witness_on_the_concave_arc_interior_is_attributed_to_the_arc() -> None:
    """Above the arc apex (10, 5) by 2^-46: the arc wins exactly with depth 2^-46."""
    boundary = _concave_arc()
    attribution = attribute_witness((F(10), F(5) + HAIRLINE_OFFSET), boundary.primitives)
    assert attribution.primitive_index == 0
    assert attribution.feature == "concave arc"
    assert attribution.squared_depth.compare(HAIRLINE_OFFSET * HAIRLINE_OFFSET) == 0


def test_a_junction_nearest_witness_ties_to_the_lower_index_and_is_labelled_junction() -> None:
    """(9, 9) is sqrt(2) from the reflex corner (10, 10), the end of line 2 and the start of line 3."""
    boundary = _l_shape()
    attribution = attribute_witness((F(9), F(9)), boundary.primitives)
    assert attribution.primitive_index == 2
    assert attribution.feature == "junction"
    assert attribution.squared_depth.compare(OneRootNumber.of(2)) == 0


def test_an_interior_tie_between_two_walls_goes_to_the_lower_index() -> None:
    """(1, 1) in the L-shape is exactly 1 from line 0 (y = 0) and line 5 (x = 0)."""
    boundary = _l_shape()
    attribution = attribute_witness((F(1), F(1)), boundary.primitives)
    assert attribution.primitive_index == 0
    assert attribution.feature == "line"
    assert attribution.squared_depth.compare(OneRootNumber.of(1)) == 0


def test_a_sub_ulp_difference_between_walls_is_decided_exactly() -> None:
    """(1 + 2^-80, 1): line 0 (y = 0) is nearer than line 5 (x = 0) by 2^-80, below any binary64 resolution at 1 mm."""
    boundary = _l_shape()
    attribution = attribute_witness((F(1) + F(1, 2**80), F(1)), boundary.primitives)
    assert attribution.primitive_index == 0
    witness = (F(1), F(1) + F(1, 2**80))
    attribution = attribute_witness(witness, boundary.primitives)
    assert attribution.primitive_index == 5


def test_box_pruning_keeps_every_nearest_piece() -> None:
    """Attribution within the box candidates equals attribution over all pieces (seeded points, both pockets)."""
    rng = random.Random(20261022)
    for boundary in (_concave_arc(), _l_shape()):
        attributor = WitnessAttributor.build(boundary.primitives)
        for _ in range(40):
            x, y = F(rng.randint(0, 2**20), 2**16), F(rng.randint(0, 2**20), 2**16)
            points = [(x + F(rng.randint(0, 2**10), 2**12), y + F(rng.randint(0, 2**10), 2**12)) for _ in range(4)]
            low = (min(p[0] for p in points), min(p[1] for p in points))
            high = (max(p[0] for p in points), max(p[1] for p in points))
            within = attributor.candidates(low, high)
            for point in points:
                pruned, full = attributor.attribute(point, within), attributor.attribute(point)
                assert (pruned.primitive_index, pruned.feature) == (full.primitive_index, full.feature)
                assert pruned.squared_depth.compare(full.squared_depth) == 0


def test_a_point_on_the_circle_beyond_the_arc_span_measures_to_the_arc_end() -> None:
    """(22.5, -7.5) is on the arc's circle but outside its span: nearest is the junction (20, 0), squared distance 62.5."""
    boundary = _concave_arc()
    attribution = attribute_witness((F(45, 2), F(-15, 2)), boundary.primitives)
    assert attribution.primitive_index == 0
    assert attribution.feature == "junction"
    assert attribution.squared_depth.compare(OneRootNumber.of(F(125, 2))) == 0


# ---------------------------------------------------------------------------
# Sweep oracle
# ---------------------------------------------------------------------------

TOOL = F(1)
# Sub-ulp offset at the 10 mm scale: 2^-60 mm.
TINY = F(1, 2**60)


def _oracle() -> SweptOracle:
    # Annulus about (10, 10), path radius 3: closed band 2 <= |p - c| <= 4.
    # Capsule from (20, 0) to (25, 0), tool 1: closed stadium.
    return SweptOracle.build(TOOL, [(F(10), F(10), F(3))], [((F(20), F(0)), (F(25), F(0)))])


def test_annulus_boundaries_are_closed_and_decided_exactly() -> None:
    oracle = _oracle()
    assert oracle.covers((F(14), F(10)))  # outer circle, exactly R + r
    assert oracle.covers((F(12), F(10)))  # inner circle, exactly R - r
    assert not oracle.covers((F(14) + TINY, F(10)))
    assert not oracle.covers((F(12) - TINY, F(10)))
    assert not oracle.covers((F(10), F(10)))  # the hole


def test_capsule_boundaries_are_closed_and_decided_exactly() -> None:
    oracle = _oracle()
    assert oracle.covers((F(22), F(1)))  # flank, exactly r
    assert oracle.covers((F(26), F(0)))  # cap apex, exactly r
    assert not oracle.covers((F(22), F(1) + TINY))
    assert not oracle.covers((F(26) + TINY, F(0)))


def test_points_on_cell_boundaries_find_their_sweeps() -> None:
    """Box corners and integer grid lines: (6, 10) is the annulus box's left edge and an integer cell edge."""
    oracle = _oracle()
    assert oracle.covers((F(6), F(10)))
    assert oracle.covers((F(19), F(0)))  # capsule box left edge, cap apex


def test_the_index_agrees_with_a_brute_force_scan() -> None:
    """Seeded random annuli, capsules and points: covers == any exact membership over all sweeps."""
    rng = random.Random(20261022)
    annuli = [(F(rng.randint(0, 2**12), 2**8), F(rng.randint(0, 2**12), 2**8), F(rng.randint(2**6, 2**10), 2**8)) for _ in range(30)]
    capsules = []
    for _ in range(30):
        a = (F(rng.randint(0, 2**12), 2**8), F(rng.randint(0, 2**12), 2**8))
        capsules.append((a, (a[0] + F(rng.randint(-(2**9), 2**9), 2**8), a[1] + F(rng.randint(-(2**9), 2**9), 2**8))))
    oracle = SweptOracle.build(TOOL, annuli, capsules)
    for _ in range(400):
        point = (F(rng.randint(-(2**8), 2**12 + 2**8), 2**8), F(rng.randint(-(2**8), 2**12 + 2**8), 2**8))
        assert oracle.covers(point) == oracle.member_of_any(point)
