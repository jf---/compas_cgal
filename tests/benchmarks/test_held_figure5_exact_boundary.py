from __future__ import annotations

import numpy as np

from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import _coverage_2


def test_figure5_exact_center_boundary_preserves_order_kind_and_lineage() -> None:
    case = load_held_reference_case("figure5")
    spec = case.pocket_spec()
    boundary = np.asarray(
        [[point.x, point.y, point.z] for point in spec.polygon.points],
        dtype=np.float64,
    )

    cycle = _coverage_2.build_center_boundary_cycle(
        boundary,
        [],
        float(spec.tool_radius),
    )

    assert cycle.counterclockwise
    assert len(cycle.primitives) == 82
    assert {primitive.kind for primitive in cycle.primitives} == {
        "arc",
        "line",
    }
    assert all(primitive.source_piece_records for primitive in cycle.primitives)
    assert all(
        first.end_mm == second.start_mm
        for first, second in zip(
            cycle.primitives,
            (*cycle.primitives[1:], cycle.primitives[0]),
        )
    )
    assert all(primitive.arc_radius_mm > 0.0 for primitive in cycle.primitives if primitive.kind == "arc")
