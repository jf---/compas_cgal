"""Exact C0 line/arc construction; source fitting and G1 are separate contracts."""

import math

import pytest

from compas_cgal import _coverage_2 as native


def test_semicircle_chain_preserves_exact_endpoints_and_sampling() -> None:
    upper = native.NativeBoundaryCurve2.arc((1.0, 0.0), (-1.0, 0.0), (0.0, 0.0), True)
    bottom = native.NativeBoundaryCurve2.line((-1.0, 0.0), (1.0, 0.0))
    boundary = native.NativeBoundary2([upper, bottom])
    assert upper.end == bottom.start
    assert bottom.end == upper.start
    assert boundary.cycle.counterclockwise
    arc = boundary.cycle.primitives[0]
    assert arc.kind == "arc"
    assert arc.sample(0.0) == upper.start
    assert arc.sample(1.0) == upper.end
    assert arc.sample(0.5).reporting_xy_mm == pytest.approx((0.0, 1.0))
    assert arc.sample(0.25).reporting_xy_mm[0] > 0.0
    assert arc.sample(0.75).reporting_xy_mm[0] < 0.0
    transition = boundary.cycle.ccw_transition(arc.sample(0.25), arc.sample(0.75))
    assert transition[0].start == arc.sample(0.25)
    assert transition[-1].end == arc.sample(0.75)
    assert boundary.design_region().contains(0.0, 0.5)
    assert not boundary.design_region().contains(0.0, -0.5)


def test_center_projection_retains_authored_join_despite_radius_mismatch() -> None:
    # Use a normal sized one-ulp displacement for the native exact arithmetic.
    center_x = math.nextafter(1.0, 2.0) - 1.0
    arc = native.NativeBoundaryCurve2.arc((1.0, 0.0), (-1.0, 0.0), (center_x, 0.0), True)
    assert arc.center_mm == pytest.approx((0.0, 0.0))
    assert arc.center_adjustment_mm == pytest.approx(center_x, rel=0.0, abs=math.ulp(center_x))
    closing = native.NativeBoundaryCurve2.line((-1.0, 0.0), (1.0, 0.0))
    chain = native.NativeBoundary2([arc, closing])
    assert chain.cycle.primitives[0].end == closing.start


def test_clockwise_arc_orientation_and_major_arc_splitting() -> None:
    arc = native.NativeBoundaryCurve2.arc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0), False)
    closing = native.NativeBoundaryCurve2.line((0.0, 1.0), (1.0, 0.0))
    chain = native.NativeBoundary2([arc, closing])
    assert not chain.cycle.counterclockwise
    arcs = [piece for piece in chain.cycle.primitives if piece.kind == "arc"]
    assert len(arcs) == 2
    assert all(not piece.arc_counterclockwise for piece in arcs)
    assert arcs[0].start == arc.start
    assert arcs[-1].end == arc.end
    assert arcs[0].end == arcs[1].start
    assert arcs[0].sample(0.5).reporting_xy_mm == pytest.approx((0.0, -1.0))


@pytest.mark.parametrize("bad", [(math.nan, 0.0), (0.0, math.inf)])
def test_nonfinite_input_rejected(bad: tuple[float, float]) -> None:
    with pytest.raises(native.InvalidNativeBoundaryCurveError):
        native.NativeBoundaryCurve2.line(bad, (1.0, 0.0))
    with pytest.raises(native.InvalidNativeBoundaryCurveError):
        native.NativeBoundaryCurve2.arc((1.0, 0.0), (-1.0, 0.0), bad, True)


def test_degenerate_curve_and_disconnected_chain_rejected() -> None:
    with pytest.raises(native.InvalidNativeBoundaryCurveError):
        native.NativeBoundaryCurve2.line((0.0, 0.0), (0.0, 0.0))
    with pytest.raises(native.InvalidNativeBoundaryCurveError):
        native.NativeBoundaryCurve2.arc((0.0, 0.0), (0.0, 0.0), (1.0, 1.0), True)
    with pytest.raises(native.InvalidNativeBoundaryChainError):
        native.NativeBoundary2([])
    with pytest.raises(native.InvalidNativeBoundaryChainError):
        native.NativeBoundary2([native.NativeBoundaryCurve2.line((0.0, 0.0), (1.0, 0.0))])


def test_self_intersecting_and_retraced_chains_rejected() -> None:
    points = [(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]
    curves = [native.NativeBoundaryCurve2.line(a, b) for a, b in zip(points, points[1:] + points[:1])]
    with pytest.raises(native.InvalidNativeBoundaryChainError):
        native.NativeBoundary2(curves)
    curves = [native.NativeBoundaryCurve2.line((0.0, 0.0), (1.0, 0.0)), native.NativeBoundaryCurve2.line((1.0, 0.0), (0.0, 0.0))]
    with pytest.raises(native.InvalidNativeBoundaryChainError):
        native.NativeBoundary2(curves)
