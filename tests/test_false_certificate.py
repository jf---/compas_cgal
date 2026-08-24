"""Independent P0 falsifiers retained as decision-side negative controls.

These probes never certify. They use the exact station predicate only to prove
that a proposed witness is live or dead before the native three-way adapter is
asked to interpret proof state.
"""

from __future__ import annotations

import math

from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

TOOL_RADIUS = 0.5
CAP_RADIANS = math.pi / 2.0
CAP_CHORD_RATIO = _stock_2.cap_chord_ratio(CAP_RADIANS)
RIB_THICKNESS = 0.004
RIB_FACETS = 256
SHORT_MOTION = 0.025
DEPLETION_CHORD_BOUND = 0.02
DEPLETION_CENTER_LIMIT = 4096
SPIRAL_GROWTH = 0.005
SPIRAL_HALF_TURN = 2.4
SPIRAL_STEPS = round(RIB_FACETS * 2.0 * SPIRAL_HALF_TURN / (2.0 * math.pi))
SPIRAL_PROBE_CENTRE = (0.0, 0.0034)
SPIRAL_DIRECTIONS = (
    (1.0, 0.0),
    (12.0 / 13.0, 5.0 / 13.0),
    (4.0 / 5.0, 3.0 / 5.0),
    (3.0 / 5.0, 4.0 / 5.0),
    (0.0, 1.0),
    (-3.0 / 5.0, 4.0 / 5.0),
)
SECTOR_EXTENT = 0.75 * math.pi
SECTOR_FACETS = round(RIB_FACETS * SECTOR_EXTENT / (2.0 * math.pi))


def _ngon(radius: float) -> Polygon:
    return Polygon(
        [
            (
                radius * math.cos(2.0 * math.pi * index / RIB_FACETS),
                radius * math.sin(2.0 * math.pi * index / RIB_FACETS),
                0.0,
            )
            for index in range(RIB_FACETS)
        ]
    )


def _rib_stock() -> Stock:
    return Stock(
        _ngon(TOOL_RADIUS + 0.5 * RIB_THICKNESS),
        holes=[_ngon(TOOL_RADIUS - 0.5 * RIB_THICKNESS)],
    )


def _machined_rib_stock() -> Stock:
    """Leave the same rib by a bore and Task 3 exact full-arc depletion."""
    stock = Stock(Polygon([(-6, -6, 0), (6, -6, 0), (6, 6, 0), (-6, 6, 0)]))
    stock.subtract_disk(0.0, 0.0, TOOL_RADIUS - 0.5 * RIB_THICKNESS)
    guide_radius = 2.0 * TOOL_RADIUS + 0.5 * RIB_THICKNESS
    motion = _stock_2.classify_audit_arc(
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        guide_radius,
        0.0,
        2.0 * math.pi,
        False,
        0.0,
        5.0,
        "cut",
    )
    assert isinstance(motion, _stock_2.AuditArcMotion2)
    stock.raw.subtract_exact_arc(
        motion,
        TOOL_RADIUS,
        DEPLETION_CHORD_BOUND,
        DEPLETION_CENTER_LIMIT,
    )
    return stock


def _sector_rib_stock() -> Stock:
    outer = TOOL_RADIUS + 0.5 * RIB_THICKNESS
    inner = TOOL_RADIUS - 0.5 * RIB_THICKNESS
    points = [
        (
            outer * math.cos(-0.5 * SECTOR_EXTENT + SECTOR_EXTENT * index / SECTOR_FACETS),
            outer * math.sin(-0.5 * SECTOR_EXTENT + SECTOR_EXTENT * index / SECTOR_FACETS),
            0.0,
        )
        for index in range(SECTOR_FACETS + 1)
    ]
    points.extend(
        (
            inner * math.cos(-0.5 * SECTOR_EXTENT + SECTOR_EXTENT * index / SECTOR_FACETS),
            inner * math.sin(-0.5 * SECTOR_EXTENT + SECTOR_EXTENT * index / SECTOR_FACETS),
            0.0,
        )
        for index in reversed(range(SECTOR_FACETS + 1))
    )
    return Stock(Polygon(points))


def _spiral_rib_stock() -> Stock:
    angles = [-SPIRAL_HALF_TURN + 2.0 * SPIRAL_HALF_TURN * index / SPIRAL_STEPS for index in range(SPIRAL_STEPS + 1)]

    def wall(theta: float, offset: float) -> tuple[float, float, float]:
        radius = TOOL_RADIUS + SPIRAL_GROWTH * theta + offset
        return radius * math.cos(theta), radius * math.sin(theta), 0.0

    points = [wall(angle, 0.5 * RIB_THICKNESS) for angle in angles]
    points.extend(wall(angle, -0.5 * RIB_THICKNESS) for angle in reversed(angles))
    return Stock(Polygon(points))


def _station_exceeds(stock: Stock, point: tuple[float, float]) -> bool:
    return bool(
        _stock_2.engagement_at(
            stock.raw,
            point[0],
            point[1],
            TOOL_RADIUS,
            CAP_CHORD_RATIO,
            0.0,
        )[2]
    )


def _arc_point(
    centre: tuple[float, float],
    guide_radius: float,
    start_angle: float,
    sweep: float,
    fraction: float,
) -> tuple[float, float]:
    angle = start_angle + fraction * sweep
    return (
        centre[0] + guide_radius * math.cos(angle),
        centre[1] + guide_radius * math.sin(angle),
    )


# Production mutation caught: treating two station-green endpoints as complete
# coverage certifies a segment whose exact interior station violates the cap.
def test_annular_rib_witness_is_live_between_safe_stations() -> None:
    stock = _rib_stock()
    half = 0.5 * SHORT_MOTION

    assert not _station_exceeds(stock, (-half, 0.0))
    assert _station_exceeds(stock, (0.0, 0.0))
    assert not _station_exceeds(stock, (half, 0.0))


# Production mutation caught: allowing the P0 falsifier to depend on the
# retired trigonometric arc sweep hides divergence from the Task 3 surrogate.
def test_machined_rib_from_exact_arc_has_the_same_live_interior_witness() -> None:
    stock = _machined_rib_stock()
    half = 0.5 * SHORT_MOTION

    assert not _station_exceeds(stock, (-half, 0.0))
    assert _station_exceeds(stock, (0.0, 0.0))
    assert not _station_exceeds(stock, (half, 0.0))


# Production mutation caught: deriving liveness from a certifier verdict lets a
# vanished or misplaced rib turn a negative control green without a witness.
def test_live_and_dead_witnesses_are_distinguished_without_the_certifier() -> None:
    stock = _rib_stock()

    assert _station_exceeds(stock, (0.0, 0.0))
    assert not any(_station_exceeds(stock, (x, 0.0)) for x in (0.20, 0.225, 0.25))


# Production mutation caught: specializing the repair to one constant-radius
# annulus misses a positive-area spiral-rib family and oblique motion directions.
def test_spiral_rib_witness_is_live_in_every_control_direction() -> None:
    stock = _spiral_rib_stock()
    px, py = SPIRAL_PROBE_CENTRE
    half = 0.5 * SHORT_MOTION

    assert _station_exceeds(stock, SPIRAL_PROBE_CENTRE)
    for cosine, sine in SPIRAL_DIRECTIONS:
        dx = half * cosine
        dy = half * sine
        assert not _station_exceeds(stock, (px - dx, py - dy))
        assert not _station_exceeds(stock, (px + dx, py + dy))


# Production mutation caught: interpreting zero contact at both guarded stations
# as evidence of clearance misses an open sector engaged above the authored cap.
def test_sector_rib_has_zero_contact_stations_and_live_interior() -> None:
    stock = _sector_rib_stock()
    half = 0.5 * SHORT_MOTION

    for x in (-half, half):
        total, run, exceeded = _stock_2.engagement_at(
            stock.raw,
            x,
            0.0,
            TOOL_RADIUS,
            CAP_CHORD_RATIO,
            0.0,
        )
        assert total == run == 0.0
        assert not exceeded
    assert _station_exceeds(stock, (0.0, 0.0))


# Predecessor-evidence mutation caught: an endpoint-only trigonometric research
# adapter misses the known rib witness. Native Task 4 tests own all decisions.
def test_predecessor_trig_arc_has_safe_endpoints_and_live_midpoint() -> None:
    stock = _rib_stock()
    for guide_radius in (0.1, 0.5, 1.0):
        sweep = SHORT_MOTION / guide_radius
        centre = (0.0, -guide_radius)
        start = 0.5 * math.pi - 0.5 * sweep
        assert not _station_exceeds(stock, _arc_point(centre, guide_radius, start, sweep, 0.0))
        assert _station_exceeds(stock, _arc_point(centre, guide_radius, start, sweep, 0.5))
        assert not _station_exceeds(stock, _arc_point(centre, guide_radius, start, sweep, 1.0))
