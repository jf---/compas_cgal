"""Analytic curved-boundary medial contacts, including finite-arc ownership."""

import subprocess
import sys

import pytest

from compas_cgal import _coverage_2 as native


def _disk() -> native.NativeBoundary2:
    return native.NativeBoundary2(
        [
            native.NativeBoundaryCurve2.arc((5.0, 0.0), (-5.0, 0.0), (0.0, 0.0), True),
            native.NativeBoundaryCurve2.arc((-5.0, 0.0), (5.0, 0.0), (0.0, 0.0), True),
        ]
    )


def test_disk_source_focal_event_constructs_coupled_circle() -> None:
    circle = _disk().circle_on_piece(0, 0.5, 1.0)
    assert circle.p_mm == pytest.approx((0.0, 5.0))
    assert circle.m_mm == pytest.approx((0.0, 0.0))
    assert circle.q_mm == pytest.approx((0.0, 4.0))
    assert circle.center_mm == pytest.approx((0.0, 2.0))
    assert circle.clearance_mm == pytest.approx(5.0)
    assert circle.guide_radius_mm == pytest.approx(2.0)
    assert circle.competing_arc_indices
    assert not circle.is_stationary


def test_capsule_line_and_arc_have_same_medial_clearance() -> None:
    owner = native.NativeBoundary2(
        [
            native.NativeBoundaryCurve2.line((-1.0, -2.0), (1.0, -2.0)),
            native.NativeBoundaryCurve2.arc((1.0, -2.0), (1.0, 2.0), (1.0, 0.0), True),
            native.NativeBoundaryCurve2.line((1.0, 2.0), (-1.0, 2.0)),
            native.NativeBoundaryCurve2.arc((-1.0, 2.0), (-1.0, -2.0), (-1.0, 0.0), True),
        ]
    )
    # Inner halves of cap support circles would tie before clearance 2;
    # their radial feet lie outside the authored semicircles and must not win.
    line = owner.circle_on_piece(0, 0.5, 1.0)
    assert line.m_mm == pytest.approx((0.0, 0.0))
    assert line.clearance_mm == pytest.approx(2.0)
    assert line.guide_radius_mm == pytest.approx(0.5)
    assert line.competing_segment_indices == [2]
    arc = owner.circle_on_piece(1, 1.0, 1.0)
    assert arc.p_mm == pytest.approx((3.0, 0.0))
    assert arc.m_mm == pytest.approx((1.0, 0.0))
    assert arc.guide_radius_mm == pytest.approx(0.5)


def test_crescent_competing_concave_arc_beats_source_focal_event() -> None:
    owner = native.NativeBoundary2(
        [
            native.NativeBoundaryCurve2.arc((-4.0, -3.0), (-4.0, 3.0), (0.0, 0.0), True),
            native.NativeBoundaryCurve2.arc((-4.0, 3.0), (-4.0, -3.0), (-4.0, 0.0), False),
        ]
    )
    circle = owner.circle_on_piece(0, 1.0, 1.0)
    assert circle.p_mm == pytest.approx((5.0, 0.0))
    assert circle.m_mm == pytest.approx((2.0, 0.0))
    assert circle.q_mm == pytest.approx((4.0, 0.0))
    assert circle.center_mm == pytest.approx((3.0, 0.0))
    assert circle.clearance_mm == pytest.approx(3.0)
    assert circle.competing_arc_indices == [1]
    # Concave source points away from its support center; no positive focal shortcut.
    circle = owner.circle_on_piece(2, 1.0, 1.0)
    assert circle.p_mm == pytest.approx((-1.0, 0.0))
    assert circle.m_mm == pytest.approx((2.0, 0.0))
    assert circle.competing_arc_indices == [0]


def test_equal_tool_radius_is_stationary_and_larger_tool_rejected() -> None:
    owner = _disk()
    assert owner.circle_on_piece(0, 0.5, 5.0).is_stationary
    with pytest.raises(native.NoPositiveNativeBoundaryCircleError):
        owner.circle_on_piece(0, 0.5, 6.0)


@pytest.mark.parametrize(
    "index,parameter,radius", [(-1, 0.5, 1.0), (3, 0.5, 1.0), (0, -0.1, 1.0), (0, 1.1, 1.0), (0, float("nan"), 1.0), (0, 0.5, 0.0), (0, 0.5, -1.0), (0, 0.5, float("inf"))]
)
def test_invalid_queries_fail_named(index: int, parameter: float, radius: float) -> None:
    with pytest.raises(native.InvalidNativeBoundaryMedialInputError):
        _disk().circle_on_piece(index, parameter, radius)


def test_disk_one_sided_piece_endpoint_has_native_focal_contact() -> None:
    proposal = _disk().circle_on_piece(0, 0.0, 1.0)
    assert proposal.p_mm == pytest.approx((5.0, 0.0))
    assert proposal.m_mm == pytest.approx((0.0, 0.0))
    assert proposal.q_mm == pytest.approx((4.0, 0.0))


def test_clockwise_disk_uses_inward_one_sided_normal() -> None:
    owner = native.NativeBoundary2(
        [
            native.NativeBoundaryCurve2.arc((5.0, 0.0), (-5.0, 0.0), (0.0, 0.0), False),
            native.NativeBoundaryCurve2.arc((-5.0, 0.0), (5.0, 0.0), (0.0, 0.0), False),
        ]
    )
    proposal = owner.circle_on_piece(0, 0.5, 1.0)
    assert proposal.p_mm == pytest.approx((0.0, -5.0))
    assert proposal.m_mm == pytest.approx((0.0, 0.0))
    assert proposal.q_mm == pytest.approx((0.0, -4.0))
    assert proposal.guide_radius_mm == pytest.approx(2.0)


def test_coverage_only_fresh_process_registers_proposal_return_type() -> None:
    # A combined collection imports the circle module elsewhere and can mask
    # a missing cross-module return-type registration in coverage consumers.
    code = """from compas_cgal import _coverage_2 as native
owner = native.NativeBoundary2([
    native.NativeBoundaryCurve2.arc((5.,0.),(-5.,0.),(0.,0.),True),
    native.NativeBoundaryCurve2.arc((-5.,0.),(5.,0.),(0.,0.),True),
])
assert not owner.circle_on_piece(0,0.5,1.).is_stationary
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True, timeout=30)
