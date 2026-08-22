"""The tool-path figures the docs publish, and the one place their recipe lives.

A figure is a claim, so the numbers under it are measured at the moment it is
drawn, never typed into a caption: the panel subtitles below are formatted from
`benchmarks.pathmetrics` and from the operation stream itself, so a figure whose
generator changed cannot keep an old caption. Which pocket, which generators and
which settings are named constants here for the same reason -- a figure produced
by a throwaway script cannot be re-drawn after the script is gone, and cannot be
compared against the one that came before it.

`benchmarks.plotting` draws; this module decides WHAT is drawn.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from benchmarks.figure6 import controlled_path
from benchmarks.figure6 import reference_pocket
from benchmarks.mathsm import constant_spacing_path
from benchmarks.palette import Theme
from benchmarks.pathmetrics import path_length
from benchmarks.plotting import ColourBy
from benchmarks.plotting import ToolpathDrawing
from benchmarks.plotting import draw_comparison
from benchmarks.spec import PocketSpec
from compas_cgal.toolpath import ToolpathResult

# Where a published figure is written, and under what name. The docs link these
# by path, so the name is part of the interface.
DEFAULT_FIGURES_OUT = Path("docs/assets/images")
TOOLPATHS_FIGURE_NAME = "fig6_toolpaths"

# Vector, always: these are line drawings of tens of thousands of segments, and a
# raster of one is both larger and worse.
DEFAULT_FIGURE_FORMATS: Tuple[str, ...] = ("svg",)

# The baseline panel: constant machining-circle spacing with NO engagement input
# of any kind, which is what "unregulated" means here. Four tenths of a tool
# diameter is what the generator picks for itself on this pocket when `stepover`
# is left unset -- both produce 187 operations of length 1,718 -- so the panel
# shows what a user gets by not asking for engagement control. It is pinned to a
# number rather than left to the generator's default because a figure whose
# baseline moves with a default cannot be compared against the one before it.
UNREGULATED_SPACING_TOOL_DIAMETERS = 0.4

# The controlled panel's cap. The same 120 degrees the corpus defaults to, so the
# figure and the reports describe one operating point.
CONTROLLED_CAP_DEG = 120.0

UNREGULATED_LABEL = f"unregulated, constant {UNREGULATED_SPACING_TOOL_DIAMETERS:g}⌀ spacing"
CONTROLLED_LABEL = f"engagement-controlled, {CONTROLLED_CAP_DEG:.0f}° cap"

TOOLPATHS_TITLE = "Trochoidal tool paths on the same pocket"


@dataclass(frozen=True)
class FigurePanels:
    """The generated paths one figure compares.

    Attributes:
        spec: The pocket and tool every panel was generated on.
        results: Panel label to generated path, in the order the panels stack.
    """

    spec: PocketSpec
    results: Dict[str, ToolpathResult]


def toolpath_panels(spec: PocketSpec) -> FigurePanels:
    """Generate both panels of the tool-path comparison.

    Args:
        spec: The pocket and tool to generate on.

    Returns:
        The panels, baseline first.

    Warns:
        UnavoidableEngagementWarning: Raised through from the controlled
            generator where the cap could not be honoured. Never suppressed: it
            is the generator's own report about the path in the figure.
    """
    return FigurePanels(
        spec=spec,
        results={
            UNREGULATED_LABEL: constant_spacing_path(spec, UNREGULATED_SPACING_TOOL_DIAMETERS),
            CONTROLLED_LABEL: controlled_path(spec, CONTROLLED_CAP_DEG),
        },
    )


def panel_subtitle(result: ToolpathResult) -> str:
    """What one panel costs, measured from the path rather than remembered.

    Args:
        result: The generated path.

    Returns:
        A line naming its analytic length and its operation count.
    """
    return f"length {path_length(result):,.0f}  ·  {len(result.operations)} ops"


def figure_subtitle(spec: PocketSpec) -> str:
    """The line under the title that says what the reader is looking at.

    The pocket is named from the spec's own sweep coordinates rather than from
    this module's constants, so the caption follows the geometry when a caller
    draws the same figure on another pocket.

    Args:
        spec: The pocket and tool.

    Returns:
        The subtitle.
    """
    width, height = spec.params.get("width"), spec.params.get("height")
    shape = f"rect {width:g}×{height:g}" if width is not None and height is not None else spec.name.replace("_", " ")
    return f"{shape}, tool ⌀{spec.tool_diameter:g} — tool-centre path, every cutting and linking move"


def draw_toolpath_figure(panels: FigurePanels, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """Draw the published tool-path comparison.

    Coloured by operation, because the question this figure answers is where the
    cutting is and where the machine is only moving; the traversal and engagement
    views of the same paths are one keyword away for a reader who wants them.

    Args:
        panels: The generated paths.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.
    """
    return draw_comparison(
        panels.results,
        boundary=panels.spec.polygon,
        holes=list(panels.spec.holes),
        colour_by=ColourBy.OPERATION,
        tool_diameter=panels.spec.tool_diameter,
        panel_subtitles={label: panel_subtitle(result) for label, result in panels.results.items()},
        title=TOOLPATHS_TITLE,
        subtitle=figure_subtitle(panels.spec),
        theme=theme,
    )


def themed_name(name: str, theme: Theme) -> str:
    """The file stem a figure gets under one theme.

    The light variant keeps the bare name because the docs link it by path, and
    every other theme takes a suffix -- otherwise regenerating the dark variant
    would silently overwrite the published figure with one that reads as a hole
    on a white page.

    Args:
        name: The figure's base name.
        theme: Which surface it was drawn for.

    Returns:
        The stem, without a suffix.
    """
    return name if theme is Theme.LIGHT else f"{name}_{theme.value}"


def save_figure(drawing: ToolpathDrawing, out_dir: Path, name: str, formats: Sequence[str] = DEFAULT_FIGURE_FORMATS) -> Tuple[Path, ...]:
    """Write one drawing to disk, once per format.

    Args:
        drawing: The drawing to write.
        out_dir: Directory to write into; created if absent.
        name: File stem, without a suffix.
        formats: File extensions to write, one file each.

    Returns:
        The paths written, in the order the formats were given.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for suffix in formats:
        path = out_dir / f"{name}.{suffix}"
        drawing.figure.savefig(path, format=suffix, facecolor=drawing.figure.get_facecolor())
        written.append(path)
    return tuple(written)


def write_toolpath_figure(
    out_dir: Path = DEFAULT_FIGURES_OUT,
    *,
    spec: Optional[PocketSpec] = None,
    formats: Sequence[str] = DEFAULT_FIGURE_FORMATS,
    theme: Theme = Theme.LIGHT,
) -> Tuple[Path, ...]:
    """Generate, draw and write the published tool-path comparison.

    Args:
        out_dir: Directory the figure is written to; created if absent.
        spec: The pocket and tool; defaults to the Figure 6 reference pocket.
        formats: File extensions to write, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written, in the order the formats were given.

    Warns:
        UnavoidableEngagementWarning: Raised through from the controlled
            generator where the cap could not be honoured.
    """
    panels = toolpath_panels(reference_pocket() if spec is None else spec)
    return save_figure(draw_toolpath_figure(panels, theme=theme), out_dir, themed_name(TOOLPATHS_FIGURE_NAME, theme), formats)
