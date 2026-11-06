"""The part polygon of a centre-domain fixture: one tool radius outside it, every side offset and every vertex mitred."""

from __future__ import annotations

import math

from benchmarks.held_material_domain import PartPolygon
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# A part far from every fixture: no rim reaches it, so the stock decision accepts. It isolates
# a test of the pair model or the connector gate from the third conjunct of the refinement.
NO_MATERIAL = PartPolygon.build(tuple(Point2[WorldXY].build(x, y) for x, y in ((1e6, 1e6), (1e6 + 1.0, 1e6), (1e6, 1e6 + 1.0))))


def part_around(boundary: tuple[Point2[WorldXY], ...], tool_radius: ToolRadius) -> PartPolygon:
    """The material domain whose one-radius inward offset is the CCW loop ``boundary``."""
    count = len(boundary)
    radius = float(tool_radius.value)
    corners = []
    for k in range(count):
        ax, ay = float(boundary[k - 1].x), float(boundary[k - 1].y)
        bx, by = float(boundary[k].x), float(boundary[k].y)
        cx, cy = float(boundary[(k + 1) % count].x), float(boundary[(k + 1) % count].y)
        first, second = math.hypot(bx - ax, by - ay), math.hypot(cx - bx, cy - by)
        # Unit outward normals of the two sides: the CCW edge direction turned clockwise.
        n0x, n0y = (by - ay) / first, (ax - bx) / first
        n1x, n1y = (cy - by) / second, (bx - cx) / second
        scale = radius / (1.0 + n0x * n1x + n0y * n1y)
        corners.append(Point2[WorldXY].build(bx + scale * (n0x + n1x), by + scale * (n0y + n1y)))
    return PartPolygon.build(tuple(corners))
