from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import List

import matplotlib

from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line

# The recipe module draws, so it must draw without a display here too.
matplotlib.use("Agg")

from benchmarks.families.analytic import rectangle  # noqa: E402
from benchmarks.figures import CONTROLLED_CAP_DEG  # noqa: E402
from benchmarks.figures import TOOLPATHS_FIGURE_NAME  # noqa: E402
from benchmarks.figures import UNREGULATED_SPACING_TOOL_DIAMETERS  # noqa: E402
from benchmarks.figures import CONTROLLED_LABEL  # noqa: E402
from benchmarks.figures import UNREGULATED_LABEL  # noqa: E402
from benchmarks.figures import FigurePanels  # noqa: E402
from benchmarks.figures import draw_toolpath_figure  # noqa: E402
from benchmarks.figures import figure_subtitle  # noqa: E402
from benchmarks.figures import panel_subtitle  # noqa: E402
from benchmarks.figures import DEFAULT_FIGURE_FORMATS  # noqa: E402
from benchmarks.figures import _PUBLISHED_FIGURES  # noqa: E402
from benchmarks.figures import regenerate_all  # noqa: E402
from benchmarks.figures import save_figure  # noqa: E402
from benchmarks.figures import themed_name  # noqa: E402
from benchmarks.palette import Theme  # noqa: E402

# The published pocket, restated here so a change to the figure's own constants
# shows up as a failing assertion about the caption rather than silently.
POCKET_WIDTH = 20.0
POCKET_HEIGHT = 12.0
TOOL_DIAMETER = 2.0
CAP_DEG = 120.0

# A machining circle whose length is exactly 2*pi: one radius unit, so the
# expected subtitle can be written down rather than recomputed by the test.
UNIT_RADIUS = 1.0

# The end-to-end regeneration runs the real generators, so it runs them on the
# smoke pocket rather than the published one: what is under test is that the
# registry is walked and files land, not the size of the pocket.
SMOKE_WIDTH = 8.0
SMOKE_HEIGHT = 6.0


@dataclass(frozen=True)
class FakeOperation:
    """A toolpath operation with only the fields a drawing and a length read."""

    geometry: Any
    operation: str
    path_index: int


@dataclass(frozen=True)
class FakeResult:
    """A toolpath result carrying only its operation stream."""

    operations: List[FakeOperation]


def _spec() -> Any:
    """The pocket the published figure is drawn on."""
    return rectangle(width=POCKET_WIDTH, height=POCKET_HEIGHT, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)


def _result(path_index: int = 0) -> FakeResult:
    """A three-operation path: plunge, one unit circle, retract."""
    centre = [0.0, 0.0, 0.0]
    return FakeResult(
        operations=[
            FakeOperation(Line([0.0, 0.0, 4.0], [0.0, 0.0, 0.0]), "plunge", path_index),
            FakeOperation(Circle(radius=UNIT_RADIUS, frame=Frame(centre, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])), "cut", path_index),
            FakeOperation(Line([UNIT_RADIUS, 0.0, 0.0], [UNIT_RADIUS, 0.0, 4.0]), "retract", path_index),
        ]
    )


def _panels() -> FigurePanels:
    """Both panels of the published comparison, without generating anything."""
    return FigurePanels(spec=_spec(), results={UNREGULATED_LABEL: _result(0), CONTROLLED_LABEL: _result(1)})


def test_a_panel_subtitle_is_measured_from_the_path_it_sits_under() -> None:
    """The caption is formatted from the operations, so it cannot go stale."""
    assert panel_subtitle(_result()) == "length 14  ·  3 ops"


def test_the_figure_subtitle_names_the_pocket_from_the_spec_it_was_drawn_on() -> None:
    assert figure_subtitle(_spec()) == "rect 20×12, tool ⌀2 — tool-centre path, every cutting and linking move"


def test_the_baseline_panel_says_what_makes_it_unregulated() -> None:
    """A reader must be able to tell what the comparison's other side actually was."""
    assert f"{UNREGULATED_SPACING_TOOL_DIAMETERS:g}" in UNREGULATED_LABEL
    assert f"{CONTROLLED_CAP_DEG:.0f}" in CONTROLLED_LABEL


def test_the_published_figure_is_two_labelled_panels_over_one_pocket() -> None:
    drawing = draw_toolpath_figure(_panels())
    assert len(drawing.axes) == 2
    assert drawing.axes[0].get_xlim() == drawing.axes[1].get_xlim()
    labels = {text.get_text() for axis in drawing.axes for text in axis.texts}
    assert UNREGULATED_LABEL in labels
    assert panel_subtitle(_result()) in labels
    assert drawing.legend_labels == ("cut", "plunge", "retract")


def test_saving_writes_one_file_per_format_under_the_published_name(tmp_path: Path) -> None:
    written = save_figure(draw_toolpath_figure(_panels()), tmp_path, TOOLPATHS_FIGURE_NAME, ("svg", "pdf"))
    assert [path.name for path in written] == [f"{TOOLPATHS_FIGURE_NAME}.svg", f"{TOOLPATHS_FIGURE_NAME}.pdf"]
    assert written[0].read_bytes().startswith(b"<?xml")
    assert written[1].read_bytes().startswith(b"%PDF")


def test_saving_creates_a_directory_that_is_not_there_yet(tmp_path: Path) -> None:
    target = tmp_path / "images"
    written = save_figure(draw_toolpath_figure(_panels()), target, TOOLPATHS_FIGURE_NAME)
    assert written[0].parent == target


def test_the_dark_variant_gets_its_own_file_rather_than_overwriting_the_published_one() -> None:
    """The docs link the light figure by path; a dark redraw must not land on it."""
    assert themed_name(TOOLPATHS_FIGURE_NAME, Theme.LIGHT) == TOOLPATHS_FIGURE_NAME
    assert themed_name(TOOLPATHS_FIGURE_NAME, Theme.DARK) == f"{TOOLPATHS_FIGURE_NAME}_dark"


def test_regenerating_everything_walks_the_registry_and_writes_real_files(tmp_path: Path) -> None:
    """One entry point: a figure absent from the registry cannot be regenerated at all."""
    smoke = rectangle(width=SMOKE_WIDTH, height=SMOKE_HEIGHT, tool_diameter=TOOL_DIAMETER, tea_cap_deg=CAP_DEG)
    written = regenerate_all(tmp_path, spec=smoke)
    assert len(written) == len(_PUBLISHED_FIGURES) * len(DEFAULT_FIGURE_FORMATS)
    assert {path.name for path in written} == {f"{TOOLPATHS_FIGURE_NAME}.{suffix}" for suffix in DEFAULT_FIGURE_FORMATS}
    for path in written:
        assert path.read_bytes().startswith(b"<?xml")


def test_the_registry_is_not_empty() -> None:
    """An empty registry would make `regenerate_all` a silent no-op."""
    assert _PUBLISHED_FIGURES


def test_redrawing_the_same_paths_rewrites_the_same_bytes(tmp_path: Path) -> None:
    """A published figure is tracked; if its bytes churn, nobody reads its diff."""
    first = save_figure(draw_toolpath_figure(_panels()), tmp_path / "a", TOOLPATHS_FIGURE_NAME, ("svg", "pdf"))
    second = save_figure(draw_toolpath_figure(_panels()), tmp_path / "b", TOOLPATHS_FIGURE_NAME, ("svg", "pdf"))
    for one, two in zip(first, second):
        assert one.read_bytes() == two.read_bytes(), one.suffix
