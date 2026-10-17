"""P-D: a pencil is authored only when the root's swept disk lies in the design pocket (Task 1.2, review I3).

Every member is internally tangent to the root at ``P`` with ``rho_k <= rho_root``,
so its swept disk ``disk(c_k, rho_k + r)`` lies in ``disk(c_root, rho_root + r)``.
One exact containment certificate of that disk therefore covers the whole
pencil; an island inside it is refused by name, never silently cut.
"""

from __future__ import annotations

import pytest
from compas.geometry import Polygon

from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.errors import HoledDomainPencilError
from compas_cgal.pipeline.pencil import TangentPencil
from compas_cgal.pipeline.pencil import entry_containment
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)
# Root at the origin, radius 2.5: its swept disk has radius 3.5 and sits inside the 12 x 8 rectangle.
ROOT = IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(2.5), Radian(0.0))
# An island at the root's centre: inside the root disk, outside the root loop's own annulus.
ISLAND_IN_PENCIL = Polygon([[-0.25, -0.25, 0.0], [-0.25, 0.25, 0.0], [0.25, 0.25, 0.0], [0.25, -0.25, 0.0]])
# An island in a corner, far from the root disk.
ISLAND_ELSEWHERE = Polygon([[4.5, 2.5, 0.0], [4.5, 3.0, 0.0], [5.0, 3.0, 0.0], [5.0, 2.5, 0.0]])


def test_a_pencil_in_a_simply_connected_pocket_is_authored() -> None:
    pencil = TangentPencil.build(ROOT, TOOL, entry_containment(RECT_12X8, [], TOOL))
    assert pencil.member(pencil.root_steps) == ROOT


def test_a_pencil_in_a_holed_domain_is_refused_by_name() -> None:
    with pytest.raises(HoledDomainPencilError):
        TangentPencil.build(ROOT, TOOL, entry_containment(RECT_12X8, [ISLAND_IN_PENCIL], TOOL))


def test_an_island_outside_the_root_disk_does_not_refuse_the_pencil() -> None:
    pencil = TangentPencil.build(ROOT, TOOL, entry_containment(RECT_12X8, [ISLAND_ELSEWHERE], TOOL))
    assert pencil.root_steps > 0
