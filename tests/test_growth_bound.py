"""P0 theorem falsifiers retained without importing predecessor authority."""

from __future__ import annotations

import math

import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

TOOL_RADIUS = 0.5
CAP_90 = math.pi / 2.0
RIB_THICKNESS = 0.004
RIB_FACETS = 256


def _cap_ratio(radians: float) -> float:
    return _stock_2.cap_chord_ratio(radians)


def _retired_growth_bound(travel: float, radius: float) -> float:
    """The falsified predecessor lemma, retained only as test evidence."""
    return 4.0 * math.asin(travel / (2.0 * radius)) + 2.0 * math.acos(1.0 - travel / radius)


def _block() -> Polygon:
    return Polygon([(-4, -4, 0), (4, -4, 0), (4, 4, 0), (-4, 4, 0)])


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


def _integer_similarity_polygon(
    points: tuple[tuple[float, float], ...],
    a: int,
    b: int,
) -> Polygon:
    return Polygon(
        [
            (
                a * x - b * y,
                b * x + a * y,
                0.0,
            )
            for x, y in points
        ]
    )


def _dyadic_square_diamond_rib(a: int, b: int) -> Stock:
    outer_half_width = 65.0 / 128.0
    inner_axis_radius = 63.0 / 128.0
    outer = (
        (-outer_half_width, -outer_half_width),
        (outer_half_width, -outer_half_width),
        (outer_half_width, outer_half_width),
        (-outer_half_width, outer_half_width),
    )
    inner = (
        (inner_axis_radius, 0.0),
        (0.0, inner_axis_radius),
        (-inner_axis_radius, 0.0),
        (0.0, -inner_axis_radius),
    )
    return Stock(
        _integer_similarity_polygon(outer, a, b),
        holes=[_integer_similarity_polygon(inner, a, b)],
    )


def _integer_similarity_point(x: float, y: float, a: int, b: int) -> tuple[float, float]:
    return a * x - b * y, b * x + a * y


# Production mutation caught: reintroducing the factor-one growth lemma as a
# proof bound certifies a cap that the exact station predicate already exceeds.
@pytest.mark.parametrize("travel", [1e-2, 1e-3, 1e-4])
def test_retired_growth_lemma_is_falsified_by_a_tool_sized_void(travel: float) -> None:
    stock = Stock(_block())
    stock.subtract_disk(0.0, 0.0, TOOL_RADIUS)
    claimed_bound = _retired_growth_bound(travel, TOOL_RADIUS)

    base = _stock_2.engagement_at(
        stock.raw,
        0.0,
        0.0,
        TOOL_RADIUS,
        _cap_ratio(claimed_bound),
        0.0,
    )
    displaced = _stock_2.engagement_at(
        stock.raw,
        travel,
        0.0,
        TOOL_RADIUS,
        _cap_ratio(claimed_bound),
        0.0,
    )

    assert base[0] == base[1] == 0.0
    assert displaced[2]


# Production mutation caught: station-only closure misses stock in the swept
# annulus even when both endpoints pass the exact unguarded cap predicate.
def test_swept_annulus_control_contains_material_between_safe_stations() -> None:
    stock = _rib_stock()
    half_spacing = 0.0125
    ratio = _cap_ratio(CAP_90)

    left = _stock_2.engagement_at(stock.raw, -half_spacing, 0.0, TOOL_RADIUS, ratio, 0.0)
    middle = _stock_2.engagement_at(stock.raw, 0.0, 0.0, TOOL_RADIUS, ratio, 0.0)
    right = _stock_2.engagement_at(stock.raw, half_spacing, 0.0, TOOL_RADIUS, ratio, 0.0)

    assert not left[2]
    assert middle[2]
    assert not right[2]


# Production mutation caught: an axis-aligned between-station theorem misses
# the live midpoint of an exact dyadic rib after integer similarity.
@pytest.mark.parametrize(
    ("a", "b", "scale"),
    [(1, 0, 1), (3, 4, 5), (12, 5, 13), (0, 1, 1)],
)
def test_dyadic_square_diamond_rib_is_integer_similarity_invariant(
    a: int,
    b: int,
    scale: int,
) -> None:
    assert a * a + b * b == scale * scale
    stock = _dyadic_square_diamond_rib(a, b)
    tool_radius = scale * TOOL_RADIUS
    start = _integer_similarity_point(-3.0 / 32.0, 0.0, a, b)
    end = _integer_similarity_point(3.0 / 32.0, 0.0, a, b)
    midpoint = (
        0.5 * (start[0] + end[0]),
        0.5 * (start[1] + end[1]),
    )

    def exceeds(point: tuple[float, float]) -> bool:
        return bool(
            _stock_2.engagement_at(
                stock.raw,
                *point,
                tool_radius,
                _cap_ratio(CAP_90),
                0.0,
            )[2]
        )

    assert not exceeds(start)
    assert exceeds(midpoint)
    assert not exceeds(end)
