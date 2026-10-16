"""The pipeline's guide: the C++ straight-skeleton stations, read once, typed, reversible."""

from __future__ import annotations

import math

import pytest
from compas.geometry import Polygon

from benchmarks.gate import gate_pocket
from compas_cgal.pipeline.errors import EmptyPipelineGuideError
from compas_cgal.pipeline.errors import GuideTruncatedError
from compas_cgal.pipeline.guide import Guide
from compas_cgal.pipeline.guide import pipeline_spacing
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.units import Clearance
from compas_cgal.units import ClearanceZ
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius

CAP_120 = Radian(math.radians(120.0))


def _guide(name: str, max_passes: int = 1000) -> Guide:
    spec = gate_pocket(name)
    return Guide.build(
        spec.polygon,
        list(spec.holes),
        ToolRadius.build(spec.tool_radius),
        pipeline_spacing(spec.tool_diameter, CAP_120),
        Clearance.build(RADIAL_CLEARANCE_FRACTION * spec.tool_diameter),
        climb=True,
        max_passes=max_passes,
        clearance_z=ClearanceZ.build(2.0 * spec.tool_diameter),
    )


def test_the_rectangle_guide_is_a_spine_and_four_corner_chains() -> None:
    chains = _guide("rect_12x8").chains
    assert len(chains) == 5
    assert [chain.origin for chain in chains] == list(range(5))
    assert not any(chain.reversed_walk for chain in chains)


def test_reversing_a_chain_reverses_tangents_and_keeps_radii_and_turn() -> None:
    chain = _guide("rect_12x8").chains[1]
    back = chain.reversed()
    assert back.origin == chain.origin and back.reversed_walk and not chain.reversed_walk
    assert [s.radius for s in back.stations] == [s.radius for s in reversed(chain.stations)]
    assert [s.centre for s in back.stations] == [s.centre for s in reversed(chain.stations)]
    assert all((b.tangent.x, b.tangent.y) == (-f.tangent.x, -f.tangent.y) for b, f in zip(back.stations, reversed(chain.stations)))
    assert all(b.clockwise == f.clockwise for b, f in zip(back.stations, reversed(chain.stations)))
    assert back.reversed() == chain


def test_a_truncated_guide_raises_instead_of_warning() -> None:
    # Confirmed path: max_passes=1 on rect_20x12 leaves 4 skeleton chains unmachined
    # (the native `skipped` count), which the compas wrapper only warns about.
    with pytest.raises(GuideTruncatedError):
        _guide("rect_20x12", max_passes=1)


def test_a_pocket_too_small_for_the_tool_has_no_guide() -> None:
    square = Polygon([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    with pytest.raises(EmptyPipelineGuideError):
        Guide.build(square, [], ToolRadius.build(1.0), pipeline_spacing(2.0, CAP_120), Clearance.build(0.002), climb=True, max_passes=1000, clearance_z=ClearanceZ.build(4.0))
