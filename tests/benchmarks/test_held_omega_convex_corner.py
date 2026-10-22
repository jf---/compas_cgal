"""The Ω* convex-corner refusal claims a true violation only when its hypothesis and precondition hold, exactly.

Each hypothesis check has a fixture where it alone fails (so a refusal there
would be a false "no certificate exists"), next to a control where the same
corner IS refused. The arc-adjacent scope fails loud by name through the
depth budget, and a connector that is not the first arrival at the vertex is
a named contract breach.
"""

from __future__ import annotations

import math
from fractions import Fraction
from functools import cache

import pytest

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_connector_repair import OmegaConnectorRepairDepthError
from benchmarks.held_omega_convex_corner import OmegaCornerFirstArrivalError
from benchmarks.held_omega_convex_corner import _Wedge
from benchmarks.held_omega_convex_corner import sharp_vertex
from benchmarks.held_omega_sites import OmegaSite
from tests.benchmarks.test_held_exact_centre_domain import _point
from tests.benchmarks.test_held_exact_centre_domain import _Pocket
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_omega_connector_repair import CAP
from tests.benchmarks.test_held_omega_connector_repair import GAP
from tests.benchmarks.test_held_omega_connector_repair import _Fixture
from tests.benchmarks.test_held_omega_connector_repair import _provenance

F = Fraction
TOOL_RADIUS = F(1)
COS_SQ_HALF_CAP = RangeEngagementLemma.build(TOOL_RADIUS, CAP.chord_ratio).cos_sq_half_cap
# The L-shape's Ω* corner V = (19, 1) joins the bottom line y = 1 to the right line x = 19 (piece 6 starts at V).
L_CORNER_PIECE = 6
# Concave-arc pocket: the Ω* offset arc about (10, -7.5) has ρ' = 13.5; it meets the line x = 19 at a convex vertex.
OFFSET_ARC_CENTRE = (10.0, -7.5)
OFFSET_ARC_RADIUS = 13.5
# Depth for the arc-adjacent repair: enough to show the midpoints never certify, small enough to stay fast.
ARC_REPAIR_DEPTH = 4


def step_pocket() -> DeclaredBoundary:
    """The L-shape's lower-right corner with a lower lobe left of x = 18: its reflex corner (18, 0) puts a ∂Ω* arc inside B(V, 2r) outside T."""
    return _Pocket(_point(0, -5)).line(0, 18).line(1, 5).line(0, 2).line(1, 10).line(2, 20).line(3, 15).close("step")


@cache
def l_corner_fixture() -> _Fixture:
    """L-shape: sites (18.5, 1) and (19, 3) around V = (19, 1); the control the step pocket is compared with."""
    return _Fixture(l_shape_pocket(), [(18.5, 1.0), (19.0, 3.0)], (0.5, 1.0))


@cache
def step_corner_fixture() -> _Fixture:
    """The same two sites and radii, on the step pocket."""
    return _Fixture(step_pocket(), [(18.5, 1.0), (19.0, 3.0)], (0.5, 1.0))


def _sharp(fixture: _Fixture, emitted: tuple[OmegaCircle, ...], target: OmegaCircle) -> int | None:
    return sharp_vertex(fixture.loop, emitted, target, COS_SQ_HALF_CAP, TOOL_RADIUS)


def test_the_control_corner_is_refused_with_its_outgoing_piece() -> None:
    fixture = l_corner_fixture()
    first, second = fixture.circles
    assert _sharp(fixture, (first,), second) == second.site.index == L_CORNER_PIECE


def test_an_earlier_circle_whose_disk_crosses_a_wedge_line_near_the_vertex_voids_the_refusal() -> None:
    fixture = l_corner_fixture()
    first, second = fixture.circles
    # Before ``first`` on the lap: site (16.5, 1), radius 3, centre (16.5, 4): 2.5 from x = 19 (< ρ), so its
    # disk leaves T, and |c - V|² = 15.25 < (2r + ρ)² = 25, so it is not far. Earlier cuts may then have
    # cleared stock at V, and the lemma does not apply.
    site = OmegaSite.decided(fixture.loop, fixture.domain.decide_contacts([(16.5, 1.0)], GAP).contacts[0])
    crossing = place_local(site, 3.0, _provenance(9))
    centre_x, radius = F(float(crossing.center.x)), F(float(crossing.radius.value))
    assert 19 - centre_x < radius
    assert _sharp(fixture, (crossing, first), second) is None


def test_a_domain_piece_inside_the_ball_but_outside_the_wedge_voids_the_refusal() -> None:
    step, control = step_corner_fixture(), l_corner_fixture()
    arcs = [index for index, piece in enumerate(step.loop.pieces) if piece.kind is PieceKind.JUNCTION_ARC]
    assert len(arcs) == 1  # the arc about the reflex corner (18, 0), from (17, 0) to (18, 1): y < 1 near V
    first, second = step.circles
    assert [(float(c.center.x), float(c.center.y), float(c.radius.value)) for c in step.circles] == [
        (float(c.center.x), float(c.center.y), float(c.radius.value)) for c in control.circles
    ]
    assert _sharp(step, (first,), second) is None


def test_an_arc_adjacent_convex_vertex_is_outside_the_lemma_and_fails_loud_by_depth() -> None:
    angle_x = 18.0
    on_arc = (angle_x, OFFSET_ARC_CENTRE[1] + math.sqrt(OFFSET_ARC_RADIUS**2 - (angle_x - OFFSET_ARC_CENTRE[0]) ** 2))
    fixture = _Fixture(concave_arc_pocket(), [on_arc, (19.0, 5.0)], (1.0, 1.0))
    first, second = fixture.circles
    assert fixture.loop.pieces[first.site.index].kind is PieceKind.OFFSET_ARC
    assert isinstance(fixture.loop.pieces[second.site.index].support, OffsetLineSupport)
    assert second.site.index == (first.site.index + 1) % len(fixture.loop.pieces)
    assert 0 in fixture.assessment(fixture.circles).refused
    assert _sharp(fixture, (first,), second) is None
    with pytest.raises(OmegaConnectorRepairDepthError) as raised:
        fixture.repair(max_depth=ARC_REPAIR_DEPTH)
    assert (raised.value.connector, raised.value.seam) == (0, False)


def test_a_connector_after_an_earlier_one_has_passed_the_vertex_is_a_named_contract_breach() -> None:
    # Lap from (18, 1): (19, 5) is past V, so connector 0 crosses V; (17, 1) wraps; then (17, 1) -> (19, 3) crosses V again.
    fixture = _Fixture(l_shape_pocket(), [(18.0, 1.0), (19.0, 5.0), (17.0, 1.0), (19.0, 3.0)], (1.0, 1.0, 1.0, 1.0))
    with pytest.raises(OmegaCornerFirstArrivalError, match="already reaches") as raised:
        _sharp(fixture, fixture.circles[:3], fixture.circles[3])
    assert raised.value.piece == L_CORNER_PIECE


def test_emitted_sites_out_of_lap_order_are_a_named_contract_breach() -> None:
    # Lap from (16, 1): (19, 5), then (17, 1) steps back on the lap.
    fixture = _Fixture(l_shape_pocket(), [(16.0, 1.0), (19.0, 5.0), (17.0, 1.0), (19.0, 3.0)], (1.0, 1.0, 1.0, 1.0))
    with pytest.raises(OmegaCornerFirstArrivalError, match="lies before") as raised:
        _sharp(fixture, fixture.circles[:3], fixture.circles[3])
    assert raised.value.piece == L_CORNER_PIECE


def test_in_order_emitted_sites_before_the_vertex_satisfy_the_precondition() -> None:
    fixture = _Fixture(l_shape_pocket(), [(16.0, 1.0), (17.0, 1.0), (19.0, 3.0)], (1.0, 1.0, 1.0))
    assert _sharp(fixture, fixture.circles[:2], fixture.circles[2]) == L_CORNER_PIECE


def test_the_hull_far_test_is_exact_at_the_boundary() -> None:
    fixture = l_corner_fixture()
    lines = tuple(fixture.loop.pieces[index].support for index in (L_CORNER_PIECE - 1, L_CORNER_PIECE))
    wedge = _Wedge(lines[0], lines[1], (F(0), F(0)), TOOL_RADIUS)
    # Hull radius max(1, √(9/4)) = 3/2: far iff |c|² > (2 + 3/2)² = 49/4, decided without a square root.
    squared_reach = F(9, 4)
    boundary = F(7, 2)
    tiny = F(1, 2**80)
    assert not wedge.hull_far((boundary, F(0)), F(1), squared_reach)
    assert wedge.hull_far((boundary + tiny, F(0)), F(1), squared_reach)
    # Radius branch: max(2, 3/2) = 2, far iff |c| > 4.
    assert not wedge.hull_far((F(4), F(0)), F(2), squared_reach)
    assert wedge.hull_far((F(4) + tiny, F(0)), F(2), squared_reach)
