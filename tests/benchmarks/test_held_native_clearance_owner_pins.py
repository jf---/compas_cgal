"""Pins on the native oracle's line-owner lemma: an owner claim is a hint, never a verdict.

``line_owner_touches`` may turn a (segment, line) pair into a TOUCH only when
both segment ends lie exactly on the owner's offset ``n·(p - S) = r`` inside
the owner's span. Every pair whose hypotheses fail is decided in full, so the
verdict with any owner list must equal the verdict with none, and the Python
oracle's verdict. The fixtures put segments next to the L-shape's bottom line
(primitive 0, y = 0, x in [0, 20]) at heights in (0, 2] with x in [2, 8]: no
other primitive is within reach, so the line-0 decision is the verdict, and a
lemma that accepted an owner without its hypotheses would report TOUCH for
every segment dipping below y = r.
"""

from __future__ import annotations

from fractions import Fraction

from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_declared_clearance_native import NativeDeclaredClearance
from benchmarks.held_declared_clearance_native import last_counters
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import WorldXY
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket

F = Fraction
R = F(1)
BOTTOM_LINE = 0
_ORACLES: list[tuple[DeclaredBoundary, DeclaredClearance, NativeDeclaredClearance]] = []


def _oracles() -> tuple[DeclaredBoundary, DeclaredClearance, NativeDeclaredClearance]:
    if not _ORACLES:
        boundary = l_shape_pocket()
        _ORACLES.append((boundary, DeclaredClearance.build(boundary, R), NativeDeclaredClearance.build(boundary, R)))
    return _ORACLES[0]


def _decision(verdict: ClearanceVerdict) -> tuple[Clearance, int | None]:
    return verdict.clearance, verdict.primitive


def _horizontal(x0: Fraction, x1: Fraction, y: Fraction) -> SegmentElement:
    """The segment y = level from x0 to x1 (an offset of the bottom line when y == r)."""
    return SegmentElement(OneRootPoint2.rational(x0, y), OneRootPoint2.rational(x1, y), (F(1), F(0)), (F(0), F(1)), y)


def test_a_parallel_segment_below_the_offset_gouges_despite_its_owner_claim() -> None:
    """y = 1/2 over x in [5, 8], owner line 0: n·(p - S) = 1/2 != r, so the lemma refuses and line 0 gouges."""
    _, python, native = _oracles()
    element = _horizontal(F(5), F(8), F(1, 2))
    (claimed,) = native.elements([element], [(BOTTOM_LINE,)])
    assert last_counters().owner_touches == 0
    assert _decision(claimed) == (Clearance.GOUGE, BOTTOM_LINE)
    assert _decision(claimed) == _decision(python.element(element))


def test_a_segment_on_the_offset_touches_through_its_owner() -> None:
    """Positive control: y = r over x in [5, 8] satisfies both hypotheses, one lemma touch."""
    _, python, native = _oracles()
    element = _horizontal(F(5), F(8), R)
    (claimed,) = native.elements([element], [(BOTTOM_LINE,)])
    assert last_counters().owner_touches == 1
    assert _decision(claimed) == (Clearance.TOUCH, BOTTOM_LINE)
    assert _decision(claimed) == _decision(python.element(element))


def test_a_slanted_segment_from_the_offset_into_the_wall_is_decided_in_full() -> None:
    """One end on y = r, the other at y = 1/4: only one hypothesis holds, the pair gouges."""
    _, python, native = _oracles()
    start, end = ExactPoint2[WorldXY].build(F(3), R), ExactPoint2[WorldXY].build(F(7), F(1, 4))
    element = SegmentElement.between(start, end)
    (claimed,) = native.elements([element], [(BOTTOM_LINE,)])
    assert last_counters().owner_touches == 0
    assert _decision(claimed) == (Clearance.GOUGE, BOTTOM_LINE)
    assert _decision(claimed) == _decision(python.element(element))


_X = st.fractions(min_value=2, max_value=8, max_denominator=64)
_HEIGHT = st.fractions(min_value=F(1, 64), max_value=2, max_denominator=64)


@settings(max_examples=300, deadline=None)
@given(_X, _HEIGHT, _X, _HEIGHT, st.booleans())
def test_owner_claims_on_every_line_never_change_a_verdict(x0: Fraction, y0: Fraction, x1: Fraction, y1: Fraction, level: bool) -> None:
    """Every primitive claimed as owner of a segment near line 0: the verdict equals the unclaimed and the Python verdicts.

    ``level`` forces a horizontal segment at height ``y0`` (the lemma's own
    shape, offset by ``y0 - r``); otherwise the ends are independent.
    """
    if level:
        y1 = y0
    if (x0, y0) == (x1, y1):
        return
    boundary, python, native = _oracles()
    element = SegmentElement.between(ExactPoint2[WorldXY].build(x0, y0), ExactPoint2[WorldXY].build(x1, y1))
    everyone = tuple(range(len(boundary.primitives)))
    claimed, unclaimed = native.elements([element, element], [everyone, ()])
    expected = _decision(python.element(element))
    assert _decision(unclaimed) == expected
    assert _decision(claimed) == expected
    if min(y0, y1) < R:
        assert expected == (Clearance.GOUGE, BOTTOM_LINE)
