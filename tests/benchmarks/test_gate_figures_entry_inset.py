"""The entry inset of the gate figures (Task 1.2, brief §6).

A panel whose path enters by a HELIX gets a zoomed entry inset beside the row:
a plan view of the helix (its footprint), the bore disk it clears, the pencil
members and the root loop, with the entry P marked, and a side elevation of Z
against plan arc length showing the ramp and the flat closing turn. A row with
no helix gets no inset. The paths here are built by the real emitters
(`emit.helix_operation`, `emit.loop_operation`).
"""

from __future__ import annotations


import matplotlib
import pytest

matplotlib.use("Agg")

from benchmarks.gate_figures import ENTRY_ELEVATION_LABEL  # noqa: E402
from benchmarks.gate_figures import ENTRY_PLAN_LABEL  # noqa: E402
from benchmarks.gate_figures import HelixWithoutRootLoopError  # noqa: E402
from benchmarks.gate_figures import draw_row  # noqa: E402
from benchmarks.gate_figures import entry_group  # noqa: E402
from benchmarks.gate_figures import entry_inset_bounds  # noqa: E402
from benchmarks.pathgeometry import drawable_motions  # noqa: E402
from compas_cgal.incident_loop import IncidentLoop  # noqa: E402
from compas_cgal.pipeline.emit import loop_operation  # noqa: E402
from compas_cgal.toolpath import ToolpathOperation  # noqa: E402
from compas_cgal.toolpath import ToolpathResult  # noqa: E402
from compas_cgal.units import CutZ  # noqa: E402
from compas_cgal.units import GuideRadius  # noqa: E402
from compas_cgal.units import Point2  # noqa: E402
from compas_cgal.units import Radian  # noqa: E402
from compas_cgal.units import WorldXY  # noqa: E402
from tests.benchmarks.helix_fixtures import CUT_Z  # noqa: E402
from tests.benchmarks.helix_fixtures import POCKET  # noqa: E402
from tests.benchmarks.helix_fixtures import TOOL_RADIUS  # noqa: E402
from tests.benchmarks.helix_fixtures import emitted  # noqa: E402
from tests.benchmarks.helix_fixtures import first_member  # noqa: E402
from tests.benchmarks.helix_fixtures import make_helix  # noqa: E402
from tests.benchmarks.helix_fixtures import result  # noqa: E402

# The root loop: tangent to the footprint and the member at P = (1, 0), centre on the ray from P through c_h.
ROOT_CENTRE = (-1.5, 0.0)
ROOT_RADIUS = 2.5
# The walk's next loop after the root: same centre, larger, its entry off P.
NEXT_RADIUS = 3.0


def _loop(centre: tuple[float, float], radius: float, clockwise: bool) -> ToolpathOperation:
    loop = IncidentLoop.build(Point2[WorldXY].build(*centre), GuideRadius.build(radius), Radian(0.0))
    return loop_operation(loop, CutZ.build(CUT_Z), 0, clockwise=clockwise)


def _helix_entry_path(clockwise: bool = True) -> ToolpathResult:
    helix = make_helix(clockwise)
    return result([emitted(helix), first_member(helix), _loop(ROOT_CENTRE, ROOT_RADIUS, clockwise), _loop(ROOT_CENTRE, NEXT_RADIUS, clockwise)])


def _plain_path() -> ToolpathResult:
    return result([_loop(ROOT_CENTRE, ROOT_RADIUS, True), _loop(ROOT_CENTRE, NEXT_RADIUS, True)])


def _row(panels: dict[str, ToolpathResult]) -> matplotlib.figure.Figure:
    return draw_row(
        panels,
        {name: name for name in panels},
        boundary=POCKET.polygon,
        holes=[],
        tool_diameter=2.0 * TOOL_RADIUS,
        title="entry inset",
        provenance="test",
    )


def _labels(figure: matplotlib.figure.Figure) -> list[str]:
    return [str(axis.get_label()) for axis in figure.axes]


def test_a_helix_panel_gets_a_plan_and_an_elevation_inset() -> None:
    figure = _row({"frozen": _plain_path(), "pipeline": _helix_entry_path()})
    labels = _labels(figure)
    assert len(figure.axes) == 2 + 2
    assert labels.count(ENTRY_PLAN_LABEL.format(name="pipeline")) == 1
    assert labels.count(ENTRY_ELEVATION_LABEL.format(name="pipeline")) == 1
    assert ENTRY_PLAN_LABEL.format(name="frozen") not in labels


def test_a_row_without_a_helix_has_no_inset() -> None:
    figure = _row({"frozen": _plain_path(), "other": _plain_path()})
    assert len(figure.axes) == 2
    assert not any(label.startswith(("entry-plan", "entry-elevation")) for label in _labels(figure))


def test_the_entry_group_is_the_helix_its_members_and_the_root_loop() -> None:
    path = _helix_entry_path()
    group = entry_group(path)
    assert group is not None
    assert (group.helix_index, group.member_indices, group.root_index) == (0, (1,), 2)
    assert entry_group(_plain_path()) is None


@pytest.mark.parametrize("clockwise", [True, False])
def test_the_inset_view_contains_the_bore_and_the_root_loop(clockwise: bool) -> None:
    path = _helix_entry_path(clockwise)
    group = entry_group(path)
    assert group is not None
    helix = group.helix
    bore = float(helix.footprint.radius.value) + TOOL_RADIUS
    cx, cy = float(helix.footprint.centre.x), float(helix.footprint.centre.y)
    root = drawable_motions(path)[group.root_index]
    xmin, xmax, ymin, ymax = entry_inset_bounds(group, drawable_motions(path), TOOL_RADIUS)
    assert xmin < cx - bore and xmax > cx + bore and ymin < cy - bore and ymax > cy + bore
    assert all(xmin < x < xmax and ymin < y < ymax for x, y in root.points)

    figure = _row({"pipeline": path})
    (plan,) = [axis for axis in figure.axes if axis.get_label() == ENTRY_PLAN_LABEL.format(name="pipeline")]
    (lo_x, hi_x), (lo_y, hi_y) = plan.get_xlim(), plan.get_ylim()
    assert lo_x <= xmin and hi_x >= xmax and lo_y <= ymin and hi_y >= ymax


def test_a_helix_followed_by_no_loop_at_its_entry_is_refused_by_name() -> None:
    helix = make_helix(True)
    with pytest.raises(HelixWithoutRootLoopError):
        entry_group(result([emitted(helix), _loop(ROOT_CENTRE, NEXT_RADIUS, True)]))
