from __future__ import annotations

import pytest

from benchmarks.errors import DegeneratePocketError, InvalidArcRatioError, InvalidSideCountError
from benchmarks.families.complexity import arc_fraction, arc_fraction_sweep, ngon_sweep, regular_ngon


def test_ngon_area_is_held_constant_across_k() -> None:
    for k in (3, 8, 64, 256):
        spec = regular_ngon(k=k, area=100.0, tool_diameter=1.0, tea_cap_deg=120.0)
        assert abs(spec.polygon.area) == pytest.approx(100.0, rel=1e-9)


def test_ngon_vertex_count_is_exactly_the_requested_k() -> None:
    for k in (3, 16, 128):
        spec = regular_ngon(k=k, area=100.0, tool_diameter=1.0, tea_cap_deg=120.0)
        assert len(spec.polygon.points) == k
        assert int(spec.params["vertices"]) == k


def test_ngon_sweep_covers_the_requested_ks_in_order() -> None:
    specs = ngon_sweep(ks=(3, 16, 128), area=100.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [int(s.params["k"]) for s in specs] == [3, 16, 128]
    assert all(s.family == "complexity" for s in specs)


def test_ngon_rejects_fewer_than_three_sides() -> None:
    with pytest.raises(InvalidSideCountError):
        regular_ngon(k=2, area=100.0, tool_diameter=1.0, tea_cap_deg=120.0)


def test_ngon_rejects_a_non_positive_area() -> None:
    with pytest.raises(DegeneratePocketError):
        regular_ngon(k=6, area=0.0, tool_diameter=1.0, tea_cap_deg=120.0)


def test_arc_fraction_zero_is_all_straight_and_one_is_all_arc() -> None:
    straight = arc_fraction(n=16, arc_ratio=0.0, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    curved = arc_fraction(n=16, arc_ratio=1.0, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert straight.params["arc_ratio"] == 0.0
    assert curved.params["arc_ratio"] == 1.0
    # A boundary with bulged arc sides encloses more area than its chord polygon.
    assert abs(curved.polygon.area) > abs(straight.polygon.area)


def test_arc_fraction_records_the_vertex_count_that_co_varies_with_the_ratio() -> None:
    # A side cannot become curved without gaining tessellation vertices, so the
    # ratio and the element count move together. The covariate is recorded on
    # every instance rather than left for a reader of the report to infer.
    low = arc_fraction(n=16, arc_ratio=0.0, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    high = arc_fraction(n=16, arc_ratio=1.0, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert int(low.params["vertices"]) == len(low.polygon.points)
    assert int(high.params["vertices"]) == len(high.polygon.points)
    assert int(high.params["vertices"]) > int(low.params["vertices"])


def test_arc_fraction_sweep_is_monotone_in_ratio() -> None:
    specs = arc_fraction_sweep(n=16, ratios=(0.0, 0.25, 0.5, 1.0), radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    ratios = [s.params["arc_ratio"] for s in specs]
    assert ratios == sorted(ratios)


def test_arc_fraction_rejects_a_ratio_outside_the_unit_interval() -> None:
    with pytest.raises(InvalidArcRatioError):
        arc_fraction(n=16, arc_ratio=1.5, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    with pytest.raises(InvalidArcRatioError):
        arc_fraction(n=16, arc_ratio=-0.1, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)


def test_arc_fraction_rejects_fewer_than_three_sides() -> None:
    with pytest.raises(InvalidSideCountError):
        arc_fraction(n=2, arc_ratio=0.5, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
