"""Prepared figure boundaries retain their authored native endpoint chain."""

import pytest

from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from compas_cgal import _coverage_2 as native


@pytest.mark.parametrize("case_name", CANONICAL_CASE_NAMES)
def test_prepared_curve_import_preserves_exact_endpoint_chain(case_name: str) -> None:
    from benchmarks.held_native_curve_import import import_held_boundary

    case = load_held_reference_case(case_name)
    boundary = import_held_boundary(case)
    curves = tuple(boundary.curves)
    assert len(curves) == len(case.boundary.primitives)
    assert all(first.end == second.start for first, second in zip(curves, curves[1:] + curves[:1]))
    for source, curve in zip(case.boundary.primitives, curves):
        endpoint_witness = native.NativeBoundaryCurve2.line(
            (float(source.start.x), float(source.start.y)),
            (float(source.end.x), float(source.end.y)),
        )
        assert curve.start == endpoint_witness.start
        assert curve.end == endpoint_witness.end
    assert sum(piece.kind == "arc" for piece in boundary.cycle.primitives) >= sum(isinstance(source, ReferenceArc) for source in case.boundary.primitives)


def test_record_stops_keeps_a_geometry_stop_as_a_recorded_outcome(tmp_path) -> None:
    # Upper stops at a clearance-below-tool event; in record_stops mode the
    # driver renders, records the stop in the report, and returns instead of raising.
    from benchmarks.held_native_curved_circle_plot import render_native_curved_circles
    import json

    output = tmp_path / "upper.png"
    render_native_curved_circles("figure8_upper", output=output, record_stops=True)
    report = json.loads(output.with_suffix(".json").read_text())
    assert report["failure"]["type"] == "NoPositiveNativeBoundaryCircleError"
    assert report["completed_queries"] == 76
