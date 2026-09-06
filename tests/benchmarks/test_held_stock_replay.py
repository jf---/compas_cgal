"""Accumulated circle sweeps preserve history and uncut centre islands."""

from fractions import Fraction

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_stock_replay import replay_circle_stock
from compas_cgal import _stock_2
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((-5, -5), (8, -5), (8, 5), (-5, 5)))
TOOL = ToolRadius.build(0.5)


def _circle(x: float) -> Figure5CounterclockwiseCircle:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(0), 4)
    return Figure5CounterclockwiseCircle(GuideRunId(0), GuideRunStationOrdinal(0), site, site, Point2[WorldXY].build(x, 0), Point2[WorldXY].build(x + 1, 0), GuideRadius.build(1))


def test_replay_retains_cuts_older_than_the_immediate_predecessor() -> None:
    earlier, latest = _circle(0), _circle(4)
    accumulated = replay_circle_stock(BOUNDARY, (earlier, latest), TOOL)
    predecessor_only = replay_circle_stock(BOUNDARY, (latest,), TOOL)
    assert predecessor_only.contains(1, 0)
    assert not accumulated.contains(1, 0)
    assert accumulated.is_subset_of(predecessor_only)
    # At this earlier-cut location the whole smaller cutter rim is clear.
    # The exact native decision therefore differs from predecessor-only stock.
    cap_squared_chord_ratio = 2.0  # Exactly the 90-degree cap surrogate.
    assert not _stock_2.engagement_at(accumulated, 1, 0, 0.25, cap_squared_chord_ratio)[2]
    assert _stock_2.engagement_at(predecessor_only, 1, 0, 0.25, cap_squared_chord_ratio)[2]


def test_replay_does_not_invent_center_clearing_or_connectors() -> None:
    stock = replay_circle_stock(BOUNDARY, (_circle(0), _circle(4)), TOOL)
    assert stock.contains(0, 0)
    assert stock.contains(2, 0)
    assert not stock.contains(1, 0)


def test_empty_prefix_leaves_the_pocket_uncut() -> None:
    stock = replay_circle_stock(BOUNDARY, (), TOOL)
    assert stock.contains(0, 0)
    assert not stock.contains(20, 0)
