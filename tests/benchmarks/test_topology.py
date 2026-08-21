from __future__ import annotations

import pytest

from benchmarks.errors import CrowdedIslandGridError, InvalidIslandCountError
from benchmarks.families.topology import ISLAND_SIDE, island_grid, island_sweep


def test_island_grid_produces_the_requested_hole_count() -> None:
    spec = island_grid(rows=2, cols=3, tool_diameter=1.0, tea_cap_deg=120.0)
    assert len(spec.holes) == 6
    assert spec.params["islands"] == 6.0


def test_island_sweep_is_ordered_by_count() -> None:
    specs = island_sweep(counts=((1, 1), (2, 2), (3, 3)), tool_diameter=1.0, tea_cap_deg=120.0)
    assert [len(s.holes) for s in specs] == [1, 4, 9]
    assert all(s.family == "topology" for s in specs)


def test_islands_lie_strictly_inside_the_outer_boundary() -> None:
    spec = island_grid(rows=2, cols=2, tool_diameter=1.0, tea_cap_deg=120.0)
    xs = [p[0] for p in spec.polygon.points]
    ys = [p[1] for p in spec.polygon.points]
    for hole in spec.holes:
        for p in hole.points:
            assert min(xs) < p[0] < max(xs)
            assert min(ys) < p[1] < max(ys)


def test_islands_are_pairwise_disjoint() -> None:
    # The exact stock model rejects touching holes outright, so a grid dense
    # enough to make two islands meet is a broken instance, not a hard one.
    spec = island_grid(rows=4, cols=5, tool_diameter=1.0, tea_cap_deg=120.0)
    boxes = [(min(p[0] for p in h.points), min(p[1] for p in h.points), max(p[0] for p in h.points), max(p[1] for p in h.points)) for h in spec.holes]
    for i, (ax0, ay0, ax1, ay1) in enumerate(boxes):
        for bx0, by0, bx1, by1 in boxes[i + 1 :]:
            assert ax1 < bx0 or bx1 < ax0 or ay1 < by0 or by1 < ay0


def test_every_island_is_a_square_of_the_declared_side() -> None:
    spec = island_grid(rows=2, cols=2, tool_diameter=1.0, tea_cap_deg=120.0)
    for hole in spec.holes:
        width = max(p[0] for p in hole.points) - min(p[0] for p in hole.points)
        height = max(p[1] for p in hole.points) - min(p[1] for p in hole.points)
        assert width == pytest.approx(ISLAND_SIDE)
        assert height == pytest.approx(ISLAND_SIDE)


def test_grid_too_dense_for_the_tool_is_refused() -> None:
    # At this density the channel between neighbouring islands is narrower than
    # the tool needs, so the pocket between them cannot be machined at all.
    with pytest.raises(CrowdedIslandGridError):
        island_grid(rows=9, cols=9, tool_diameter=1.0, tea_cap_deg=120.0)


def test_a_tool_too_wide_for_the_channels_is_refused_at_a_workable_count() -> None:
    island_grid(rows=3, cols=3, tool_diameter=1.0, tea_cap_deg=120.0)
    with pytest.raises(CrowdedIslandGridError):
        island_grid(rows=3, cols=3, tool_diameter=4.0, tea_cap_deg=120.0)


def test_grid_rejects_a_non_positive_count() -> None:
    with pytest.raises(InvalidIslandCountError):
        island_grid(rows=0, cols=3, tool_diameter=1.0, tea_cap_deg=120.0)
    with pytest.raises(InvalidIslandCountError):
        island_grid(rows=3, cols=-1, tool_diameter=1.0, tea_cap_deg=120.0)
