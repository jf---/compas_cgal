from __future__ import annotations

import pytest

from benchmarks.families.analytic import arc_channel, axis_clearance, disk, rectangle, stadium


def test_disk_area_matches_closed_form() -> None:
    spec = disk(radius=5.0, tool_diameter=1.0, tea_cap_deg=120.0)
    # Inscribed regular polygon under-approximates the disk; 2% is ample slack
    # for the default tessellation and still catches a wrong radius.
    assert abs(spec.polygon.area) == pytest.approx(3.14159265358979 * 25.0, rel=0.02)


def test_stadium_axis_clearance_is_constant() -> None:
    spec = stadium(straight_length=20.0, half_width=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    samples = [axis_clearance(spec, x, 0.0) for x in (-8.0, -4.0, 0.0, 4.0, 8.0)]
    for value in samples:
        assert value == pytest.approx(3.0, abs=1e-9)


def test_rectangle_axis_clearance_is_half_height_on_the_spine() -> None:
    spec = rectangle(width=20.0, height=6.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert axis_clearance(spec, 0.0, 0.0) == pytest.approx(3.0, abs=1e-9)


def test_arc_channel_axis_clearance_is_constant() -> None:
    spec = arc_channel(guide_radius=10.0, half_width=2.0, sweep_deg=90.0, tool_diameter=1.0, tea_cap_deg=120.0)
    for angle_deg in (10.0, 45.0, 80.0):
        import math

        a = math.radians(angle_deg)
        x, y = 10.0 * math.cos(a), 10.0 * math.sin(a)
        assert axis_clearance(spec, x, y) == pytest.approx(2.0, abs=1e-9)


def test_families_reject_tool_wider_than_channel() -> None:
    from benchmarks.errors import DegeneratePocketError

    with pytest.raises(DegeneratePocketError):
        stadium(straight_length=20.0, half_width=0.4, tool_diameter=1.0, tea_cap_deg=120.0)
    with pytest.raises(DegeneratePocketError):
        arc_channel(guide_radius=10.0, half_width=0.4, sweep_deg=90.0, tool_diameter=1.0, tea_cap_deg=120.0)
