from __future__ import annotations

import pytest
from compas.geometry import Polygon

from benchmarks.errors import InvalidCapError, NonPositiveToolError, PocketNotSimpleError
from benchmarks.spec import PocketSpec

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


def test_build_accepts_valid_pocket() -> None:
    spec = PocketSpec.build(name="square", family="analytic", polygon=SQUARE, tool_diameter=1.0, tea_cap_deg=120.0)
    assert spec.tool_radius == pytest.approx(0.5)
    assert spec.holes == ()


def test_build_rejects_non_positive_tool() -> None:
    with pytest.raises(NonPositiveToolError):
        PocketSpec.build(name="x", family="analytic", polygon=SQUARE, tool_diameter=0.0, tea_cap_deg=120.0)


def test_build_rejects_cap_outside_range() -> None:
    with pytest.raises(InvalidCapError):
        PocketSpec.build(name="x", family="analytic", polygon=SQUARE, tool_diameter=1.0, tea_cap_deg=181.0)


def test_build_rejects_self_intersecting_boundary() -> None:
    bowtie = Polygon([[0, 0, 0], [10, 10, 0], [10, 0, 0], [0, 10, 0]])
    with pytest.raises(PocketNotSimpleError):
        PocketSpec.build(name="x", family="analytic", polygon=bowtie, tool_diameter=1.0, tea_cap_deg=120.0)


def test_spec_is_frozen() -> None:
    spec = PocketSpec.build(name="square", family="analytic", polygon=SQUARE, tool_diameter=1.0, tea_cap_deg=120.0)
    with pytest.raises(AttributeError):
        spec.tool_diameter = 2.0  # type: ignore[misc]
