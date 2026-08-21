from __future__ import annotations

import pytest
from compas.geometry import Polygon

from benchmarks.degeneracy import DEGENERACY_NAMES, degeneracy_corpus

TOOL = 1.0
CAP = 120.0


def _by_name(tool_diameter: float = TOOL) -> dict:
    return {spec.name: spec for spec in degeneracy_corpus(tool_diameter=tool_diameter, tea_cap_deg=CAP)}


def test_corpus_names_are_unique_and_stable() -> None:
    names = [s.name for s in degeneracy_corpus(tool_diameter=TOOL, tea_cap_deg=CAP)]
    assert len(names) == len(set(names))
    assert set(names) == set(DEGENERACY_NAMES)
    assert {"tangent_island", "collinear_run", "half_turn_arm", "pinch_exactly_tool", "cocircular_square"} <= set(names)


def test_every_degenerate_instance_is_still_a_valid_spec() -> None:
    for spec in degeneracy_corpus(tool_diameter=TOOL, tea_cap_deg=CAP):
        assert spec.family == "degeneracy"
        assert len(spec.polygon.points) >= 3
        assert spec.tool_radius > 0.0


def test_tangent_island_leaves_exactly_one_tool_diameter_to_the_wall() -> None:
    # The name's whole claim: a cutter of radius r centred at x = r touches the
    # wall at x = 0 and the island at x = 2r at the same instant.
    for tool in (1.0, 2.5):
        spec = _by_name(tool)["tangent_island"]
        wall_x = min(p[0] for p in spec.polygon.points)
        island_x = min(p[0] for hole in spec.holes for p in hole.points)
        assert island_x - wall_x == pytest.approx(tool, abs=1e-12)
        assert spec.params["wall_gap"] == pytest.approx(tool)


def test_collinear_run_is_exactly_collinear() -> None:
    spec = _by_name()["collinear_run"]
    a, b, c = (spec.polygon.points[i] for i in (0, 1, 2))
    cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    # Exact zero, not near-zero: these three vertices are constructed on one line
    # so CGAL::orientation returns COLLINEAR rather than a sign.
    assert cross == 0.0


def test_half_turn_arm_is_exactly_one_tool_diameter_wide() -> None:
    for tool in (1.0, 2.5):
        spec = _by_name(tool)["half_turn_arm"]
        far_x = max(p[0] for p in spec.polygon.points)
        ys = sorted({p[1] for p in spec.polygon.points if p[0] == far_x})
        assert ys[-1] - ys[0] == pytest.approx(tool, abs=1e-12)
        assert spec.params["arm_width"] == pytest.approx(tool)


def test_pinch_is_exactly_the_tool_diameter() -> None:
    for tool in (1.0, 2.5):
        spec = _by_name(tool)["pinch_exactly_tool"]
        assert spec.params["neck_width"] == pytest.approx(tool)


def test_cocircular_square_has_four_equidistant_walls() -> None:
    spec = _by_name()["cocircular_square"]
    xs = [p[0] for p in spec.polygon.points]
    ys = [p[1] for p in spec.polygon.points]
    side_x, side_y = max(xs) - min(xs), max(ys) - min(ys)
    assert side_x == pytest.approx(side_y, abs=1e-12)
    assert spec.params["inradius"] == pytest.approx(0.5 * side_x)


def test_a_repeated_vertex_cannot_be_expressed_and_so_is_not_in_the_corpus() -> None:
    # `compas.geometry.Polygon` COLLAPSES a repeated vertex on construction, so a
    # "duplicate vertex" instance would be a silent copy of the plain box under a
    # misleading name. This pins the reason the corpus omits one; if compas ever
    # stops collapsing, this fails and the instance becomes worth authoring.
    repeated = [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 12.0, 0.0], [0.0, 12.0, 0.0]]
    assert len(Polygon(repeated).points) == len(repeated) - 1
    assert "duplicate_vertex" not in DEGENERACY_NAMES
