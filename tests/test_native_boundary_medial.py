"""Analytic curved-boundary medial contacts, including finite-arc ownership."""

import math
import subprocess
import sys
import time

import pytest
from hypothesis import Phase
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import _circle_geometry_2
from compas_cgal import _coverage_2 as native

# Wall-clock budget for ONE query on generic double coordinates. Measured on
# 2026-09-07: the identity-deciding implementation took 6.5 s (rounded
# rectangle, line source) to 19.6 s (two-arc disk) per query because CORE
# refined identically-zero decisions to their root-separation bound; the
# scratch benchmark of the same decision without identities ran in 23 us to
# 131 us. Two seconds sits three orders of magnitude from both regimes, so a
# loaded CI box passes and the failure mode cannot return unnoticed.
GENERIC_DOUBLE_QUERY_BUDGET_S = 2.0


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


def _timed_query(
    owner: native.NativeBoundary2, piece: int, parameter: float, tool: float
) -> tuple[_circle_geometry_2.BoundaryNormalCircleProposal2, float]:
    started = time.perf_counter()
    proposal = owner.circle_on_piece(piece, parameter, tool)
    return proposal, time.perf_counter() - started


def _assert_coupled(proposal: _circle_geometry_2.BoundaryNormalCircleProposal2, tool: float) -> None:
    """The q-m diameter invariants, checked on reporting doubles, never re-decided in CORE."""
    p, m, q, c = proposal.p_mm, proposal.m_mm, proposal.q_mm, proposal.center_mm
    assert math.dist(m, p) == pytest.approx(proposal.clearance_mm)
    assert math.dist(q, p) == pytest.approx(tool)
    assert math.dist(c, q) == pytest.approx(proposal.guide_radius_mm)
    assert math.dist(c, m) == pytest.approx(proposal.guide_radius_mm)


def _generic_disk(cx: float, cy: float, r: float) -> native.NativeBoundary2:
    arc = native.NativeBoundaryCurve2.arc
    return native.NativeBoundary2([arc((cx + r, cy), (cx, cy + r), (cx, cy), True), arc((cx, cy + r), (cx + r, cy), (cx, cy), True)])


def _rounded_rectangle(
    x0: float, y0: float, width: float, height: float, rx: float, ry: float | None = None
) -> native.NativeBoundary2:
    """CCW rounded rectangle; unequal rx/ry gives corner arcs that are not tangent to the edges."""
    ry = rx if ry is None else ry
    arc, line = native.NativeBoundaryCurve2.arc, native.NativeBoundaryCurve2.line
    # Every shared coordinate is computed once so chain endpoints are bit-identical.
    xr, xw, xwr, yr, yh, yhr = x0 + rx, x0 + width, x0 + width - rx, y0 + ry, y0 + height, y0 + height - ry
    return native.NativeBoundary2(
        [
            line((xr, y0), (xwr, y0)),
            arc((xwr, y0), (xw, yr), (xwr, yr), True),
            line((xw, yr), (xw, yhr)),
            arc((xw, yhr), (xwr, yh), (xwr, yhr), True),
            line((xwr, yh), (xr, yh)),
            arc((xr, yh), (x0, yhr), (xr, yhr), True),
            line((x0, yhr), (x0, yr)),
            arc((x0, yr), (xr, y0), (xr, yr), True),
        ]
    )


def test_generic_double_disk_focal_query_is_fast_and_coupled() -> None:
    cx, cy, r = 0.37, -1.21, 4.123  # irrational radius: sqrt(r*r) has no exact double
    proposal, seconds = _timed_query(_generic_disk(cx, cy, r), 0, 0.5, 1.0)
    assert seconds < GENERIC_DOUBLE_QUERY_BUDGET_S
    assert proposal.m_mm == pytest.approx((cx, cy))
    assert proposal.clearance_mm == pytest.approx(r)
    assert proposal.guide_radius_mm == pytest.approx((r - 1.0) / 2)
    assert proposal.competing_arc_indices == [0, 1]
    _assert_coupled(proposal, 1.0)


def test_integer_fillet_arc_midpoint_hits_focal_event_with_join_vertices() -> None:
    # Exact fillets: both adjacent edges are tangent to the arc's circle, so the
    # focal disk touches them at the join vertices. Those ties are decided on
    # rational quantities; only the arc's own support and its two join
    # vertices may be reported, never the edges as interior contacts.
    owner = _rounded_rectangle(0.0, 0.0, 12.0, 8.0, 2.0)
    proposal, seconds = _timed_query(owner, 1, 0.5, 1.0)
    assert seconds < GENERIC_DOUBLE_QUERY_BUDGET_S
    assert proposal.m_mm == pytest.approx((10.0, 2.0))
    assert proposal.clearance_mm == pytest.approx(2.0)
    assert proposal.competing_arc_indices == [1]
    assert proposal.competing_vertex_indices == [1, 2]
    assert not proposal.competing_segment_indices
    _assert_coupled(proposal, 1.0)


def test_generic_double_fillet_arc_midpoint_is_fast() -> None:
    # Fitted centres make these fillets only nearly tangent; the adjacent join
    # vertices still lie exactly on the arc's circle.
    x0, y0, width, height, r = 0.13, -0.29, 12.37, 7.91, 1.37
    owner = _rounded_rectangle(x0, y0, width, height, r)
    arc_piece = [primitive.kind for primitive in owner.cycle.primitives].index("arc")
    proposal, seconds = _timed_query(owner, arc_piece, 0.5, 0.5)
    assert seconds < GENERIC_DOUBLE_QUERY_BUDGET_S
    assert proposal.clearance_mm == pytest.approx(r, rel=1e-6)
    _assert_coupled(proposal, 0.5)


def test_generic_double_rounded_rectangle_line_query_is_fast_and_analytic() -> None:
    x0, y0, width, height, r = 0.13, -0.29, 12.37, 7.91, 1.37
    proposal, seconds = _timed_query(_rounded_rectangle(x0, y0, width, height, r), 0, 0.5, 1.0)
    assert seconds < GENERIC_DOUBLE_QUERY_BUDGET_S
    # Bottom-edge midpoint: the opposite (parallel) edge is the first contact at half the height.
    assert proposal.m_mm == pytest.approx((x0 + width / 2, y0 + height / 2))
    assert proposal.clearance_mm == pytest.approx(height / 2)
    assert proposal.competing_segment_indices == [4]
    assert not proposal.competing_arc_indices
    assert not proposal.competing_vertex_indices
    _assert_coupled(proposal, 1.0)


# Generic-position pockets: coordinates of physical magnitude (no 1e-42 offsets
# that create near-tangencies far below double resolution) and corner arcs with
# unequal x/y extents, so no arc is tangent to an edge. Exact and near ties are
# a separate, documented cost and are witnessed by the deterministic tests.
_physical = st.floats(-50.0, 50.0).filter(lambda v: abs(v) >= 1e-3)


@settings(max_examples=25, deadline=None, phases=(Phase.explicit, Phase.reuse, Phase.generate))
@given(
    x0=_physical,
    y0=_physical,
    height=st.floats(3.0, 30.0),
    aspect=st.floats(1.1, 3.0),
    fillet=st.floats(0.05, 0.3),
    squash=st.floats(1.1, 1.5),
    tool_fraction=st.floats(0.05, 0.9),
)
def test_generic_double_rounded_rectangles_stay_fast_and_coupled(
    x0: float, y0: float, height: float, aspect: float, fillet: float, squash: float, tool_fraction: float
) -> None:
    width, rx = height * aspect, height * fillet
    ry = rx * squash
    tool = rx * tool_fraction * 0.25  # below every clearance queried here
    owner = _rounded_rectangle(x0, y0, width, height, rx, ry)
    kinds = [primitive.kind for primitive in owner.cycle.primitives]
    line_piece, arc_piece = kinds.index("line"), kinds.index("arc")
    proposal, seconds = _timed_query(owner, line_piece, 0.5, tool)
    assert seconds < GENERIC_DOUBLE_QUERY_BUDGET_S
    assert proposal.clearance_mm == pytest.approx(height / 2)
    _assert_coupled(proposal, tool)
    proposal, seconds = _timed_query(owner, arc_piece, 0.5, tool)
    assert seconds < GENERIC_DOUBLE_QUERY_BUDGET_S
    assert proposal.clearance_mm >= tool
    assert proposal.competing_vertex_indices or proposal.competing_segment_indices or proposal.competing_arc_indices
    _assert_coupled(proposal, tool)


@pytest.mark.parametrize("parameter", [0.1, 0.3, 0.5, 0.7, 0.9])
def test_arc_samples_stay_on_their_supporting_circle(parameter: float) -> None:
    # The sampler constructs its point on the circle and no longer re-decides
    # that incidence exactly; witness it here on reporting values instead.
    for cx, cy, r in ((0.0, 0.0, 5.0), (0.37, -1.21, 4.123)):
        owner = _generic_disk(cx, cy, r)
        for primitive in owner.cycle.primitives:
            assert math.dist(primitive.sample(parameter).reporting_xy_mm, (cx, cy)) == pytest.approx(r)
