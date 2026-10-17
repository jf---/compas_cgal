"""The entry inset of the gate figures (Task 1.2, brief §6).

A panel whose path enters by a HELIX gets a zoomed entry inset beside the row:
a plan view of the helix (its footprint), the bore disk it clears, the pencil
members and the root loop, with the entry P marked, and a side elevation of Z
against plan arc length showing the ramp and the flat closing turn. A row with
no helix gets no inset. The paths here are built by the real emitters
(`emit.helix_operation`, `emit.loop_operation`), and the entry's structure is
read from the generator's own `RampOutRecord` (brief N2), never re-derived
from coincident entries.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from benchmarks.gate_figures import ENTRY_ELEVATION_LABEL  # noqa: E402
from benchmarks.gate_figures import ENTRY_PLAN_LABEL  # noqa: E402
from benchmarks.gate_figures import REPO_ROOT  # noqa: E402
from benchmarks.gate_figures import UnrecordedHelixEntryError  # noqa: E402
from benchmarks.gate_figures import draw_row  # noqa: E402
from benchmarks.gate_figures import entry_group  # noqa: E402
from benchmarks.gate_figures import entry_inset_bounds  # noqa: E402
from benchmarks.pathgeometry import drawable_motions  # noqa: E402
from benchmarks.plotting import save_figure  # noqa: E402
from compas_cgal.incident_loop import IncidentLoop  # noqa: E402
from compas_cgal.pipeline.emit import loop_operation  # noqa: E402
from compas_cgal.pipeline.result import EntryKind  # noqa: E402
from compas_cgal.pipeline.result import ForcedDecisionLedger  # noqa: E402
from compas_cgal.pipeline.result import ForcedEntry  # noqa: E402
from compas_cgal.pipeline.result import ForcedEntryReason  # noqa: E402
from compas_cgal.pipeline.result import PipelineResult  # noqa: E402
from compas_cgal.pipeline.result import RampOutRecord  # noqa: E402
from compas_cgal.helix import Helix  # noqa: E402
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
# A walk loop after the root that is ALSO entered at P: tangent there, larger
# than the root. Coincidence would make it a member; the record says it is not.
TANGENT_CENTRE = (-2.0, 0.0)
TANGENT_RADIUS = 3.0

# The commit whose `benchmarks/gate_figures.py` predates the entry inset: a row
# with no helix must draw the bytes that module draws.
PRE_INSET_REV = "7b001d47"

INSET_LEGEND = ("bore (footprint + r)", "pencil member", "root loop", "P (helix entry)")


def _loop(centre: tuple[float, float], radius: float, clockwise: bool) -> ToolpathOperation:
    loop = IncidentLoop.build(Point2[WorldXY].build(*centre), GuideRadius.build(radius), Radian(0.0))
    return loop_operation(loop, CutZ.build(CUT_Z), 0, clockwise=clockwise)


def _pipeline(operations: list[ToolpathOperation], root: int) -> PipelineResult:
    """A pipeline result whose record puts the helix at 0, the members between it and *root*."""
    helix_op = operations[0]
    assert isinstance(helix_op.geometry, Helix)
    first_cut = ForcedEntry(path_index=0, entry=helix_op.geometry.entry_xy, reason=ForcedEntryReason.FIRST_CUT, into_material=True, entry_kind=EntryKind.HELIX)
    return PipelineResult(
        operations=operations,
        polyline=np.zeros((0, 3), dtype=float),
        ledger=ForcedDecisionLedger.build([first_cut]),
        omitted=(),
        links=(),
        ramp_out=RampOutRecord.build(operations, 0, 0, range(1, root), root),
    )


def _helix_entry_path(clockwise: bool = True, *, with_member: bool = True) -> PipelineResult:
    helix = make_helix(clockwise)
    members = [first_member(helix)] if with_member else []
    operations = [emitted(helix), *members, _loop(ROOT_CENTRE, ROOT_RADIUS, clockwise), _loop(ROOT_CENTRE, NEXT_RADIUS, clockwise)]
    return _pipeline(operations, root=1 + len(members))


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


def _legend(figure: matplotlib.figure.Figure) -> list[str]:
    (legend,) = figure.legends
    return [text.get_text() for text in legend.get_texts()]


def _pre_inset_module(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """`benchmarks/gate_figures.py` as of `PRE_INSET_REV`, loaded under a private name."""
    shown = subprocess.run(["git", "show", f"{PRE_INSET_REV}:benchmarks/gate_figures.py"], cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    name = "pre_inset_gate_figures"
    spec = importlib.util.spec_from_loader(name, loader=None)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    # It locates the repository from its own path, which is the same file's.
    module.__file__ = str(REPO_ROOT / "benchmarks" / "gate_figures.py")
    # dataclasses resolve their module through sys.modules while the source executes.
    monkeypatch.setitem(sys.modules, name, module)
    exec(compile(shown.stdout, f"{name}.py", "exec"), module.__dict__)
    return module


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


def test_the_entry_group_is_the_generators_ramp_out_record() -> None:
    path = _helix_entry_path()
    group = entry_group(path)
    assert group is not None
    record = path.ramp_out
    assert (group.helix_index, group.member_indices, group.root_index) == (record.helix_operation, tuple(record.members), record.root_operation)
    assert entry_group(_plain_path()) is None


def test_a_walk_loop_entered_at_p_after_the_root_is_not_a_member() -> None:
    helix = make_helix(True)
    operations = [emitted(helix), first_member(helix), _loop(ROOT_CENTRE, ROOT_RADIUS, True), _loop(TANGENT_CENTRE, TANGENT_RADIUS, True)]
    assert operations[3].entry_xy == operations[0].entry_xy
    group = entry_group(_pipeline(operations, root=2))
    assert group is not None
    assert (group.member_indices, group.root_index) == ((1,), 2)


def test_an_entry_with_no_members_draws_the_root_and_no_member_legend() -> None:
    path = _helix_entry_path(with_member=False)
    group = entry_group(path)
    assert group is not None
    assert (group.member_indices, group.root_index) == ((), 1)
    legend = _legend(_row({"pipeline": path}))
    assert "pencil member" not in legend
    assert [label for label in INSET_LEGEND if label != "pencil member"] == [label for label in legend if label in INSET_LEGEND]


def test_the_row_legend_carries_the_inset_entries_once() -> None:
    legend = _legend(_row({"frozen": _plain_path(), "pipeline": _helix_entry_path()}))
    assert [label for label in legend if label in INSET_LEGEND] == list(INSET_LEGEND)
    assert len(legend) == len(set(legend))
    assert not set(INSET_LEGEND) & set(_legend(_row({"frozen": _plain_path()})))


def test_a_row_without_a_helix_draws_the_same_bytes_as_before_the_inset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    before = _pre_inset_module(monkeypatch)
    panels = {"frozen": _plain_path(), "other": _plain_path()}
    arguments = {"boundary": POCKET.polygon, "holes": [], "tool_diameter": 2.0 * TOOL_RADIUS, "title": "entry inset", "provenance": "test"}
    old = save_figure(before.draw_row(panels, {name: name for name in panels}, **arguments), tmp_path / "before.png")
    new = save_figure(_row(panels), tmp_path / "after.png")
    assert new.read_bytes() == old.read_bytes()


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


def test_a_helix_without_a_ramp_out_record_is_refused_by_name() -> None:
    helix = make_helix(True)
    with pytest.raises(UnrecordedHelixEntryError):
        entry_group(result([emitted(helix), first_member(helix), _loop(ROOT_CENTRE, ROOT_RADIUS, True)]))
