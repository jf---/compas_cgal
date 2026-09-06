"""Exact proposal contacts enter the existing boundary transition consumer."""

import numpy as np
import pytest

from compas_cgal import _circle_geometry_2 as circles
from compas_cgal import _coverage_2 as coverage

RECTANGLE = [(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)]
L_SHAPE = [(0.0, 0.0), (6.0, 0.0), (6.0, 2.0), (2.0, 2.0), (2.0, 6.0), (0.0, 6.0)]


def _cycle(boundary: list[tuple[float, float]], radius: float) -> coverage.ReachableBoundaryCycle2:
    return coverage.build_center_boundary_cycle(np.array([(x, y, 0.0) for x, y in boundary]), [], radius)


def test_segment_proposals_keep_exact_contacts_in_line_transition() -> None:
    owner = circles.BoundaryNormalCircle2(RECTANGLE)
    cycle = _cycle(RECTANGLE, 1.0)
    start = coverage.boundary_circle_contact(owner.query(0, 0.3, 1.0))
    end = coverage.boundary_circle_contact(owner.query(0, 0.6, 1.0))
    transition = cycle.ccw_transition(start, end)
    assert len(transition) == 1
    assert transition[0].kind == "line"
    assert transition[0].start == start
    assert transition[0].end == end


def test_reflex_contacts_keep_exact_irrational_arc_and_sector_joins() -> None:
    owner = circles.BoundaryNormalCircle2(L_SHAPE)
    cycle = _cycle(L_SHAPE, 0.25)
    points = [
        coverage.boundary_circle_contact(owner.query(2, 0.875, 0.25)),
        coverage.boundary_circle_contact(owner.query_vertex(3, (0.0, -1.0), 0.25)),
        coverage.boundary_circle_contact(owner.query_vertex(3, (-1.0, -1.0), 0.25)),
        coverage.boundary_circle_contact(owner.query_vertex(3, (-1.0, 0.0), 0.25)),
        coverage.boundary_circle_contact(owner.query(3, 0.125, 0.25)),
    ]
    for start, end, kind in zip(points, points[1:], ["line", "arc", "arc", "line"]):
        transition = cycle.ccw_transition(start, end)
        assert transition[0].start == start
        assert transition[-1].end == end
        assert all(piece.kind == kind for piece in transition)
        assert all(a.end == b.start for a, b in zip(transition, transition[1:]))
        if kind == "arc":
            assert all(not piece.arc_counterclockwise for piece in transition)
    arc = next(piece for piece in cycle.primitives if piece.kind == "arc")
    assert points[1] == arc.start
    assert points[3] == arc.end


def test_foreign_proposal_contact_is_rejected_by_existing_geometry_check() -> None:
    cycle = _cycle(RECTANGLE, 1.0)
    owner = circles.BoundaryNormalCircle2([(x + 100.0, y) for x, y in RECTANGLE])
    foreign = coverage.boundary_circle_contact(owner.query(0, 0.5, 1.0))
    with pytest.raises(coverage.ReachableArrangementTopologyError, match="does not lie"):
        cycle.ccw_transition(cycle.primitives[0].start, foreign)


def test_bridge_does_not_accept_reporting_coordinates() -> None:
    proposal = circles.BoundaryNormalCircle2(RECTANGLE).query(0, 0.5, 1.0)
    with pytest.raises(TypeError):
        coverage.boundary_circle_contact(proposal.q_mm)  # type: ignore[arg-type]
