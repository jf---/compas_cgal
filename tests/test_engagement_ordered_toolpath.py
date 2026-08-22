"""Tests for chain ordering: the link predicate, the orientation choice, the entry count."""

from __future__ import annotations

from typing import List

import pytest
from compas.geometry import Circle
from compas.geometry import Polygon

from compas_cgal import _coverage_2
from compas_cgal import engagement_ordered_toolpath as ordered_module
from compas_cgal.engagement_ordered_toolpath import LINK_PROBE_SPACING_TOOL_RADII
from compas_cgal.engagement_ordered_toolpath import MIN_LINK_PROBES
from compas_cgal.engagement_ordered_toolpath import STRAIGHT_MOVE_CAP_FRACTION
from compas_cgal.engagement_ordered_toolpath import _largest_admissible_advance
from compas_cgal.engagement_ordered_toolpath import _orientations
from compas_cgal.engagement_ordered_toolpath import _straight_move_is_admissible
from compas_cgal.engagement_ordered_toolpath import _straight_move_ratio
from compas_cgal.engagement_ordered_toolpath import _straight_probe_positions
from compas_cgal.engagement_ordered_toolpath import chain_ordered_toolpath
from compas_cgal.engagement_rho_toolpath import _rho_stations
from compas_cgal.engagement_toolpath import UnavoidableEngagementWarning
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.engagement_toolpath import _Regulation
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import _polygon_to_ccw_vertices

TOOL_DIAMETER = 2.0
TOOL_RADIUS = 1.0
CAP_DEG = 120.0
SYNTHETIC_SPACING = 0.05

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])


def _regulation() -> _Regulation:
    return _Regulation.build(
        tool_diameter=TOOL_DIAMETER,
        tea_cap_deg=CAP_DEG,
        guide_step_tool_diameters=0.025,
        max_advance_tool_diameters=1.0,
        radial_clearance=None,
        cut_z=0.0,
        clearance_z=None,
    )


def _centre_domain(polygon: Polygon = RECT_12X8):
    return _coverage_2.ReachableDomain2(_polygon_to_ccw_vertices(polygon), [], TOOL_RADIUS).center_domain()


def _straight_chain(count: int, slope: float, first_radius: float) -> List[_GuideStation]:
    return [
        _GuideStation(
            cx=index * SYNTHETIC_SPACING,
            cy=0.0,
            radius=first_radius + slope * index * SYNTHETIC_SPACING,
            clockwise=True,
            tx=1.0,
            ty=0.0,
        )
        for index in range(count)
    ]


def _result():
    with pytest.warns(UnavoidableEngagementWarning):
        return chain_ordered_toolpath(RECT_12X8, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)


# ---------------------------------------------------------------------------
# The link predicate. Two conditions, and containment is the one an
# engagement-only test misses -- it let a gouging link through on L_shape.
# ---------------------------------------------------------------------------


def test_a_link_that_leaves_the_centre_domain_is_refused_however_little_it_engages() -> None:
    domain = _centre_domain()
    stock = Stock(RECT_12X8, holes=None)
    regulation = _regulation()
    # The pocket spans x in [-6, 6]; a tool of radius 1 may not be centred past
    # x = 5. Target x = 5.9 is inside the pocket and outside the centre domain.
    assert domain.contains(0.0, 0.0)
    assert not domain.contains(5.9, 0.0)
    assert not _straight_move_is_admissible(stock, (0.0, 0.0), (5.9, 0.0), regulation, domain)


def test_a_straight_move_is_held_to_a_stricter_threshold_than_a_machining_circle() -> None:
    regulation = _regulation()
    assert STRAIGHT_MOVE_CAP_FRACTION < 1.0
    assert _straight_move_ratio(regulation) < regulation.cap_ratio, "a straight move must be judged more strictly than a circle, never less"


def test_a_link_through_solid_stock_is_refused() -> None:
    domain = _centre_domain()
    stock = Stock(RECT_12X8, holes=None)
    assert not _straight_move_is_admissible(stock, (-4.0, 0.0), (4.0, 0.0), _regulation(), domain), "a full-length cut through virgin stock is a slot, not a link"


def test_probe_spacing_puts_a_tool_diameter_of_overlap_between_evaluated_positions() -> None:
    positions = _straight_probe_positions((0.0, 0.0), (10.0, 0.0), TOOL_RADIUS)
    assert len(positions) >= MIN_LINK_PROBES
    assert positions[0] == (0.0, 0.0)
    assert positions[-1] == pytest.approx((10.0, 0.0))
    gaps = [abs(later[0] - earlier[0]) for earlier, later in zip(positions, positions[1:])]
    assert max(gaps) <= LINK_PROBE_SPACING_TOOL_RADII * TOOL_RADIUS + 1e-9


def test_a_move_shorter_than_one_tool_radius_still_evaluates_both_ends() -> None:
    positions = _straight_probe_positions((0.0, 0.0), (0.1, 0.0), TOOL_RADIUS)
    assert len(positions) == MIN_LINK_PROBES


# ---------------------------------------------------------------------------
# Orientation.
# ---------------------------------------------------------------------------


def test_a_chain_is_offered_in_both_walk_directions() -> None:
    chain = _straight_chain(8, slope=0.0, first_radius=3.0)
    orientations = _orientations(chain, origin=2, tool_radius=TOOL_RADIUS)
    assert len(orientations) == 2
    assert {orientation.reversed_walk for orientation in orientations} == {False, True}
    assert all(orientation.origin == 2 for orientation in orientations)
    forward, backward = orientations
    assert forward.entry_point != backward.entry_point, "the two directions must begin at opposite ends"


def test_a_chain_with_no_trochoidal_station_is_offered_in_no_direction() -> None:
    chain = _straight_chain(6, slope=0.0, first_radius=0.5 * TOOL_RADIUS)
    assert _orientations(chain, origin=0, tool_radius=TOOL_RADIUS) == []


def test_reversing_a_chain_reverses_its_tangents_and_keeps_its_radii() -> None:
    chain = _straight_chain(5, slope=0.0, first_radius=3.0)
    reversed_chain = ordered_module._reversed_guide_chain(chain)
    assert [station.radius for station in reversed_chain] == [station.radius for station in reversed(chain)]
    assert [(station.tx, station.ty) for station in reversed_chain] == [(-1.0, 0.0)] * len(chain)
    assert [station.clockwise for station in reversed_chain] == [station.clockwise for station in chain]


# ---------------------------------------------------------------------------
# The advance scan keeps its direction, now with the bridge condition.
# ---------------------------------------------------------------------------


def test_the_advance_scan_still_returns_the_furthest_admissible_station(monkeypatch) -> None:
    stations = _rho_stations(_straight_chain(8, slope=0.0, first_radius=2.0))
    admissible = {2, 5}
    seen: List[int] = []

    def fake_station(_stock, guide_station, _advance, _tool_radius, _cap_ratio) -> bool:
        index = round(guide_station.cx / SYNTHETIC_SPACING)
        seen.append(index)
        return index in admissible

    monkeypatch.setattr(ordered_module, "_station_is_admissible", fake_station)
    monkeypatch.setattr(ordered_module, "_straight_move_is_admissible", lambda *_a, **_k: True)
    index, forced = _largest_admissible_advance(None, stations, 0, 6, _regulation(), None)
    assert index == 5, f"expected the furthest admissible station 5, got {index}"
    assert forced is False
    assert seen[0] == 6, "the scan must begin at the furthest candidate"


def test_an_inadmissible_bridge_refuses_a_candidate_whose_circle_complies(monkeypatch) -> None:
    """The bridge condition is load-bearing: without it a tip-ward walk slots."""
    stations = _rho_stations(_straight_chain(8, slope=0.0, first_radius=2.0))
    monkeypatch.setattr(ordered_module, "_station_is_admissible", lambda *_a, **_k: True)
    monkeypatch.setattr(ordered_module, "_straight_move_is_admissible", lambda *_a, **_k: False)
    _index, forced = _largest_admissible_advance(None, stations, 0, 6, _regulation(), None)
    assert forced is True, "every candidate's circle complies; only the bridge refuses, and that must still refuse"


# ---------------------------------------------------------------------------
# End to end: the point of the module.
# ---------------------------------------------------------------------------


def test_the_rectangle_is_entered_once_instead_of_once_per_chain() -> None:
    result = _result()
    plunges = [op for op in result.operations if op.operation is OperationType.PLUNGE]
    assert len(plunges) == 1, f"the four corner chains are reachable through the spine's cut stock; got {len(plunges)} plunges"


def test_the_one_unavoidable_entry_is_named_and_explained() -> None:
    result = _result()
    assert len(result.forced_entries) == 1
    isolated = result.forced_entries[0]
    assert isolated.path_index == 0
    assert "first cut" in isolated.reason
    assert isolated.best_link_length == 0.0, "the first chain has nothing to link from"


def test_declined_regions_survive_the_reordering() -> None:
    result = _result()
    assert result.declined_regions, "the corners are still below the degeneracy floor whatever the order"


def test_no_emitted_machining_circle_is_a_bore() -> None:
    result = _result()
    radii = [float(op.geometry.radius) for op in result.operations if isinstance(op.geometry, Circle)]
    assert radii and min(radii) > TOOL_RADIUS


# ---------------------------------------------------------------------------
# Chain numbering. A chain with no trochoidal station anywhere never enters the
# machining sequence, so it has no machining index to borrow -- and borrowing one
# makes `path_index` name two different chains at once.
# ---------------------------------------------------------------------------

# Body twelve by eight with a four-wide arm on top. The arm is narrower than a
# trochoid needs at this tool -- its clearance is 2.0 against a floor of
# 1.0 + r + radial clearance -- so the guide's two arm-spine chains carry no
# machining circle at all, while the eight body chains do. That is the shape the
# numbering bug needs and the corpus does not have. AXIS-ALIGNED deliberately:
# `ReachableDomain2.center_domain()` is milliseconds on axis-parallel input and
# seconds to minutes on oblique input (`docs/oblique_edge_cost.md`), and a test
# is not the place to pay that.
T_ARM4 = Polygon(
    [
        [0.0, 0.0, 0.0],
        [12.0, 0.0, 0.0],
        [12.0, 8.0, 0.0],
        [8.0, 8.0, 0.0],
        [8.0, 12.0, 0.0],
        [4.0, 12.0, 0.0],
        [4.0, 8.0, 0.0],
        [0.0, 8.0, 0.0],
    ]
)

# The narrow arm starts here. Every wholly declined run on `T_ARM4` lies above it
# and every machined loop below it, so "which chain is this region on" needs no
# tolerance to answer.
T_ARM4_ARM_FLOOR_Y = 9.0

# Chains of `T_ARM4` whose every station is below the degeneracy floor.
T_ARM4_FULLY_DECLINED_CHAINS = 2


def _arm_result():
    with pytest.warns(UnavoidableEngagementWarning):
        return chain_ordered_toolpath(T_ARM4, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)


def test_a_chain_declined_in_full_is_numbered_after_every_machined_chain() -> None:
    """`path_index` must name ONE chain.

    A chain with no trochoidal station emits no operations, so it cannot be
    stamped with a machining index without that index also naming a chain that
    does. Numbering the fully declined chains after the machined ones keeps the
    field meaning one thing.
    """
    result = _arm_result()
    machined = {op.path_index for op in result.operations}
    declined = {region.path_index for region in result.declined_regions}
    declined_only = declined - machined
    assert len(declined_only) == T_ARM4_FULLY_DECLINED_CHAINS, (
        f"expected {T_ARM4_FULLY_DECLINED_CHAINS} chains numbered outside the machining sequence, got {sorted(declined_only)}"
    )
    assert min(declined_only) > max(machined), "a chain that emits no operations must not borrow a machined chain's index"


def test_a_declined_region_lies_on_the_chain_its_index_names() -> None:
    """The consequence a consumer would actually hit: joining regions to operations by index."""
    result = _arm_result()
    loops_by_index: dict = {}
    for op in result.operations:
        if isinstance(op.geometry, Circle):
            loops_by_index.setdefault(op.path_index, []).append(float(op.geometry.frame.point[1]))
    for region in result.declined_regions:
        loop_ys = loops_by_index.get(region.path_index)
        if loop_ys is None:
            continue
        in_arm = region.first_center[1] > T_ARM4_ARM_FLOOR_Y and region.last_center[1] > T_ARM4_ARM_FLOOR_Y
        chain_in_arm = min(loop_ys) > T_ARM4_ARM_FLOOR_Y
        assert in_arm == chain_in_arm, (
            f"region at y={region.first_center[1]:.2f} is stamped path_index {region.path_index}, whose loops sit at y={min(loop_ys):.2f}..{max(loop_ys):.2f}"
        )


def test_the_fully_declined_arm_is_still_reported_not_dropped() -> None:
    """Renumbering must not lose the regions -- the arm material is still declared."""
    result = _arm_result()
    arm = [region for region in result.declined_regions if region.first_center[1] > T_ARM4_ARM_FLOOR_Y]
    assert len(arm) == T_ARM4_FULLY_DECLINED_CHAINS
    assert all(region.largest_gouge_free_radius <= TOOL_RADIUS for region in arm)


def test_a_pocket_with_no_fully_declined_chain_numbers_every_region_onto_a_machined_chain() -> None:
    """The corpus case, pinned so the renumbering cannot start inventing indices."""
    result = _result()
    machined = {op.path_index for op in result.operations}
    declined = {region.path_index for region in result.declined_regions}
    assert declined <= machined, "with every chain machined, every declined run belongs to one of them"
