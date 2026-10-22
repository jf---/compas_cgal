"""Exact values of immutable sweep data are computed once per object and decide exactly as before.

``ArcSweep.full`` is decided once at construction (filtered, then exact);
``PointBox.of`` reuses the enclosure of a recently boxed point object; the
connector continuity and junction checks read ``points_apart``. None of these
changes a decision: each test pins the exact answer on a near tie the interval
filter cannot decide, and pins that the cached value is not recomputed.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from benchmarks.held_connector_engagement_arcs import ArcSweep
from benchmarks.held_connector_engagement_exact import ExactConnector
from benchmarks.held_connector_engagement_exact import ExactConnectorPath
from benchmarks.held_connector_engagement_lines import LineSweep
from benchmarks.held_connector_engagement_range import DiscontinuousConnectorError
from benchmarks.held_connector_engagement_range import PointBox
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import points_apart
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_dyadic_interval import FRACTION_BITS
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import WorldXY

ORIGIN = ExactPoint2[WorldXY].build(Fraction(0), Fraction(0))
# Half-angle tangent step far below the interval filter's resolution
# (2^-FRACTION_BITS mm), so the two points' enclosures overlap and only the
# exact comparison can tell them apart.
SUB_FILTER_STEP = Fraction(1, 1 << (FRACTION_BITS + 24))


def _on_unit_circle(t: Fraction) -> OneRootPoint2:
    """The rational point of the unit circle with half-angle tangent ``t``."""
    return OneRootPoint2.rational((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))


def _one_root_on_unit_circle() -> OneRootPoint2:
    """``(1/2, √3/2)`` on the unit circle: a point in ``Q(√3)``, built twice as distinct equal objects."""
    return OneRootPoint2.build(OneRootNumber.of(Fraction(1, 2)), OneRootNumber.build(0, Fraction(1, 2), 3))


class _ExactComparisonCalled(AssertionError):
    """The exact point comparison ran when the cached value should have answered."""


def _forbid_exact_point_comparison(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(self: OneRootPoint2, other: OneRootPoint2) -> int:
        raise _ExactComparisonCalled(f"{self!r} vs {other!r}")

    monkeypatch.setattr(OneRootPoint2, "compare_xy", refuse)


def test_full_is_exact_on_a_near_tie_the_filter_cannot_decide() -> None:
    start, end = _on_unit_circle(Fraction(1, 3)), _on_unit_circle(Fraction(1, 3) + SUB_FILTER_STEP)
    box_start, box_end = PointBox.of(start), PointBox.of(end)
    assert (box_start.x - box_end.x).sign() is None and (box_start.y - box_end.y).sign() is None
    assert not ArcSweep.build(ORIGIN, Fraction(1), 1, start, end).full
    assert ArcSweep.build(ORIGIN, Fraction(1), 1, start, _on_unit_circle(Fraction(1, 3))).full


def test_full_turn_between_distinct_equal_one_root_objects() -> None:
    first, second = _one_root_on_unit_circle(), _one_root_on_unit_circle()
    assert first is not second
    assert ArcSweep.build(ORIGIN, Fraction(1), -1, first, second).full


def test_full_is_decided_once_at_construction(monkeypatch: pytest.MonkeyPatch) -> None:
    point = _one_root_on_unit_circle()
    sweep = ArcSweep.build(ORIGIN, Fraction(1), 1, point, _one_root_on_unit_circle())
    near_tie = ArcSweep.build(ORIGIN, Fraction(1), 1, _on_unit_circle(Fraction(1, 3)), _on_unit_circle(Fraction(1, 3) + SUB_FILTER_STEP))
    _forbid_exact_point_comparison(monkeypatch)
    circle = PriorCircle(ExactPoint2[WorldXY].build(Fraction(3), Fraction(1)), Fraction(2))
    for _ in range(3):
        assert sweep.full
        assert not near_tie.full
        sweep.squared_distances(circle)
        near_tie.squared_distances(circle)


def test_point_box_is_reused_for_the_same_object_and_fresh_for_an_equal_one() -> None:
    first, second = _one_root_on_unit_circle(), _one_root_on_unit_circle()
    box = PointBox.of(first)
    assert PointBox.of(first) is box
    other = PointBox.of(second)
    assert other is not box
    assert other == box


def test_points_apart_is_exact_on_the_near_tie_and_reflexive_without_comparing(monkeypatch: pytest.MonkeyPatch) -> None:
    start, end = _on_unit_circle(Fraction(1, 3)), _on_unit_circle(Fraction(1, 3) + SUB_FILTER_STEP)
    assert points_apart(start, end)
    assert not points_apart(start, _on_unit_circle(Fraction(1, 3)))
    _forbid_exact_point_comparison(monkeypatch)
    assert not points_apart(start, start)


def test_continuity_of_shared_endpoint_objects_needs_no_exact_comparison(monkeypatch: pytest.MonkeyPatch) -> None:
    a, b, c = OneRootPoint2.rational(0, 0), _one_root_on_unit_circle(), OneRootPoint2.rational(3, 1)
    first, second = LineSweep.build(a, b), LineSweep.build(b, c)
    _forbid_exact_point_comparison(monkeypatch)
    connector = ExactConnector.build((first, second))
    circles = (PriorCircle(ORIGIN, Fraction(1)), PriorCircle(ExactPoint2[WorldXY].build(Fraction(3), Fraction(0)), Fraction(1)))
    assert ExactConnectorPath.build(circles, (connector,), (a, c)).connectors == (connector,)


def test_continuity_still_decides_exactly_across_distinct_objects() -> None:
    a, c = OneRootPoint2.rational(0, 0), OneRootPoint2.rational(3, 1)
    assert ExactConnector.build((LineSweep.build(a, _one_root_on_unit_circle()), LineSweep.build(_one_root_on_unit_circle(), c))).sweeps
    near = _on_unit_circle(Fraction(1, 3))
    with pytest.raises(DiscontinuousConnectorError):
        ExactConnector.build((LineSweep.build(a, near), LineSweep.build(_on_unit_circle(Fraction(1, 3) + SUB_FILTER_STEP), c)))
