"""The nine machining-quality figures, and the one place their recipe lives.

Every figure here is a CLAIM about a generated tool path, so every number in
every caption is measured at the moment the figure is drawn -- from
`benchmarks.survey`, `benchmarks.quality` and `benchmarks.coverage` -- and never
typed in. A caption that cannot go stale is the whole reason these live in a
committed module instead of a throwaway script.

`benchmarks.plotting` draws tool paths; `benchmarks.palette` owns every colour;
this module decides WHAT is drawn and what the reader should conclude from it.
Chart panels that are not tool paths (a curve, a histogram, a scatter) are built
against matplotlib directly, but they take every colour from the same validated
palette and paint the same surface, so the page reads as one system.

COLOUR, BY JOB. Magnitude -- engagement, curvature -- takes the single-hue
ordinal ramp. Identity -- which path, which panel -- takes the fixed categorical
slots in `Palette.cut`, `.link`, `.lead`, in that order and never cycled. Text is
always an ink token and never a series colour, so a label is never mistaken for a
datum.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Callable
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

import numpy as np
from matplotlib import rc_context
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.patheffects import withStroke

from benchmarks.coverage import COVERAGE_GRID_SAMPLES
from benchmarks.coverage import _grid_axes
from benchmarks.coverage import measure_coverage
from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import gate_pocket
from benchmarks.models import MachineModel
from benchmarks.models import MaterialModel
from benchmarks.palette import Palette
from benchmarks.palette import Theme
from benchmarks.palette import palette_for
from benchmarks.pathmetrics import entry_cut_indices
from benchmarks.plotting import LABEL_HALO_ALPHA
from benchmarks.plotting import LABEL_HALO_PAD_PT
from benchmarks.plotting import Z_LABEL
from benchmarks.plotting import ColourBy
from benchmarks.plotting import ToolpathDrawing
from benchmarks.plotting import draw_toolpath
from benchmarks.plotting import save_drawing
from benchmarks.quality import FULL_TURN_DEG
from benchmarks.quality import chip_thickness_ratio_from_rim
from benchmarks.quality import measure_quality
from benchmarks.spec import PocketSpec
from benchmarks.survey import MotionQuality
from benchmarks.survey import PathSurvey
from benchmarks.survey import survey_path
from compas_cgal import _coverage_2
from compas_cgal.stock import Stock
from compas_cgal.stock import _polygon_to_ccw_vertices
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

# The pocket every measured figure is drawn on, and the one the corner defect was
# found on by eye. Named here rather than passed in: these figures annotate
# specific motions, so a figure drawn on another pocket would carry captions
# about geometry it is not showing.
FIGURE_POCKET = "rect_20x12"

# The second pocket, used only where a figure compares two instances.
FIGURE_SMALL_POCKET = "rect_12x8"

# The generator the measured figures show. The radius-regulated generator emits
# an identical motion stream at this cap, so drawing both would be two copies of
# one picture; the docs say so in prose instead.
FIGURE_GENERATOR = "engagement_controlled"

# Stepover, in tool diameters, for the deliberately under-covered path in the
# residual figure. Nine tenths of a diameter leaves the loops barely touching, so
# the residue is a visible lattice rather than a hairline -- the point of that
# panel is to show what the metric detects, not to show a good path.
UNDERCOVERED_STEPOVER_TOOL_DIAMETERS = 0.9

# Margin around the annotated motions in the corner zoom, in tool diameters. Wide
# enough to carry the pocket wall and the neighbouring proper trochoid, tight
# enough that the corner still reads as a corner.
CORNER_ZOOM_MARGIN_TOOL_DIAMETERS = 0.75

# The corner figure's own geometry. The panel takes the left half and the labels
# a right-hand column, because every label on the first draft of this figure
# landed on top of the path it was pointing at.
CORNER_FIGURE_HEIGHT_IN = 3.9
CORNER_PANEL_BOX = (0.03, 0.12, 0.50, 0.70)
CORNER_LABEL_X_FRACTION = 1.08
CORNER_LABEL_FRACTIONS: Tuple[float, ...] = (0.90, 0.52, 0.14)

# The inset, in axes fractions of the main panel, and how many loop radii it
# spans. It sits in the strip BEYOND the pocket wall, which is empty, and low
# enough that it does not cover the corner it exists to explain. Three radii puts the 0.02-radius loop at a third of the inset's width,
# which is the point at which a reader can see it is a CLOSED LOOP and not a
# dot; at sixteen it was still a speck.
CORNER_INSET_BOX = (0.64, 0.05, 0.34, 0.34)
CORNER_INSET_LOOP_RADII = 3.0

# Samples across the chip-thinning curve. One per half degree over 0-180 renders
# the plateau corner at 90 degrees without a visible facet.
CHIP_CURVE_SAMPLES = 361

# Where the chip reaches full thickness, in the kernel's ENGAGED RIM ARC: the
# textbook plateau is at a quarter turn of the entry-to-exit angle, and the rim
# arc is twice that (`benchmarks.quality.textbook_engagement_deg`).
CHIP_PLATEAU_RIM_DEG = 180.0

# The two formulas the figures carry on their own faces, so a reader never has to
# take an axis label on trust.
CHIP_FORMULA = "$h_{ex} = f_z\\,\\sin(\\min(θ_{rim}/2,\\,90°))$"
FEED_FORMULA = "$v \\leq \\sqrt{a_{max}/κ}$"

# The engagement range this project cannot yet deliver, in degrees. Held's
# comparison sweeps caps from 20 upward, and no generator here honours a cap
# under about 60 on a pocket entered from solid stock.
UNREACHABLE_CAP_LOW_DEG = 20.0
UNREACHABLE_CAP_HIGH_DEG = 60.0

# Figure geometry, in inches. The width matches `benchmarks.plotting`'s so a
# chart and a tool-path panel sit at the same measure on the page.
CHART_WIDTH_IN = 6.5
CHART_HEIGHT_IN = 3.6
STACKED_CHART_HEIGHT_IN = 5.0
FIGURE_DPI = 200.0

# Type sizes, in points, matching `benchmarks.plotting`.
TITLE_PT = 10.5
SUBTITLE_PT = 8.5
LABEL_PT = 8.5
TICK_PT = 8.0
ANNOTATION_PT = 7.5
LEGEND_PT = 8.0

# Stroke weights, in points.
CURVE_WIDTH_PT = 1.8
REFERENCE_WIDTH_PT = 1.1
HAIRLINE_PT = 0.8

# Alpha for a shaded band. Light enough that a curve crossing it stays readable.
BAND_ALPHA = 0.16
FLOOR_ALPHA = 0.28

# Halo behind an annotation, in points, so a label over a dense path stays legible.
LABEL_HALO_PT = 2.0

# Marker size for a point event, in points.
EVENT_MARKER_PT = 5.5

# Residual shading marker size, in points squared, for the coverage scatter.
RESIDUAL_MARKER_AREA_PT2 = 1.4

# The residual figure's height. Two 20x12 panels stacked need less than the
# generic stacked height, which left a band of empty surface under the second.
RESIDUAL_FIGURE_HEIGHT_IN = 4.4

# Every figure name this module publishes. The docs link these by path, so a name
# is part of the interface.
CORNER_DEFECT_NAME = "quality_corner_defect"
COVERAGE_RESIDUAL_NAME = "quality_coverage_residual"
ENGAGEMENT_ALONG_PATH_NAME = "quality_engagement_along_path"
CHIP_THINNING_NAME = "quality_chip_thinning"
ENGAGEMENT_MAP_NAME = "quality_engagement_map"
CURVATURE_FEED_NAME = "quality_curvature_feed"
LENGTH_VS_TIME_NAME = "quality_length_vs_time"
ENGAGEMENT_HISTOGRAM_NAME = "quality_engagement_histogram"
MATERIAL_ENTRIES_NAME = "quality_material_entries"


@dataclass(frozen=True)
class MeasuredPath:
    """One generated path and everything the figures read off it.

    Generating and surveying a path costs seconds, and six figures want the same
    one, so it is done once and passed around.

    Attributes:
        spec: The pocket, tool and cap.
        result: The generated tool path.
        survey: The per-motion measurement.
    """

    spec: PocketSpec
    result: ToolpathResult
    survey: PathSurvey


def measured_path(pocket: str = FIGURE_POCKET, generator: str = FIGURE_GENERATOR) -> MeasuredPath:
    """Generate and survey one gate path.

    Args:
        pocket: A `benchmarks.gate` pocket name.
        generator: A `benchmarks.gate` generator name.

    Returns:
        The path and its survey.

    Warns:
        UnavoidableEngagementWarning: Raised through from the generator where the
            cap could not be honoured. Never suppressed: it is the generator's
            own report about the path in the figure.
    """
    spec = gate_pocket(pocket)
    result = GATE_GENERATORS[generator](spec)
    return MeasuredPath(spec=spec, result=result, survey=survey_path(spec, result))


# ---------------------------------------------------------------------------
# Chart scaffolding shared by the panels that are not tool paths.
# ---------------------------------------------------------------------------


def _chart_figure(palette: Palette, *, panels: int = 1, height_in: float = CHART_HEIGHT_IN, share_x: bool = False) -> Tuple[Figure, Tuple[Any, ...]]:
    """A figure painted with the theme's surface, with an Agg canvas attached.

    Args:
        palette: The theme's tokens.
        panels: How many panels to stack.
        height_in: Total figure height.
        share_x: Whether stacked panels share one x axis.

    Returns:
        The figure and its panels, top to bottom.
    """
    figure = Figure(figsize=(CHART_WIDTH_IN, height_in), dpi=FIGURE_DPI, facecolor=palette.surface)
    FigureCanvasAgg(figure)
    axes = figure.subplots(panels, 1, sharex=share_x, squeeze=False)[:, 0]
    for axis in axes:
        _style_axes(axis, palette)
    return figure, tuple(axes)


def _style_axes(axis: Any, palette: Palette) -> None:
    """Recessive grid, recessive spines, ink-coloured ticks.

    Args:
        axis: The axes to style.
        palette: The theme's tokens.
    """
    axis.set_facecolor(palette.surface)
    axis.grid(True, color=palette.grid, linewidth=HAIRLINE_PT, zorder=0)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color(palette.grid)
        axis.spines[side].set_linewidth(HAIRLINE_PT)
    axis.tick_params(colors=palette.secondary, labelsize=TICK_PT, length=0.0)


def _title(figure: Figure, palette: Palette, title: str, subtitle: str) -> None:
    """Put a title and one subtitle line at the top of a chart figure.

    Args:
        figure: The figure.
        palette: The theme's tokens.
        title: Primary-ink heading.
        subtitle: Secondary-ink line beneath it.
    """
    figure.text(0.012, 0.975, title, ha="left", va="top", fontsize=TITLE_PT, color=palette.ink)
    figure.text(0.012, 0.918, subtitle, ha="left", va="top", fontsize=SUBTITLE_PT, color=palette.secondary)


def _legend(axis: Any, palette: Palette, handles: Sequence[Any], labels: Sequence[str], location: str = "best") -> None:
    """A legend in ink, on the surface, with no frame competing with the data.

    Args:
        axis: The axes to attach it to.
        palette: The theme's tokens.
        handles: Proxy artists, in legend order.
        labels: Their labels, in the same order.
        location: A matplotlib legend location; ``best`` unless a panel has a
            known empty corner, in which case naming it beats letting the
            solver drop the box on a bar.
    """
    legend = axis.legend(handles, labels, loc=location, fontsize=LEGEND_PT, frameon=True, facecolor=palette.surface, edgecolor=palette.grid, framealpha=0.9)
    for text in legend.get_texts():
        text.set_color(palette.secondary)


def _annotate(axis: Any, palette: Palette, text: str, xy: Tuple[float, float], xytext: Tuple[float, float], colour: Optional[str] = None) -> None:
    """A direct label with a leader line and a surface halo.

    Args:
        axis: The axes.
        palette: The theme's tokens.
        text: The label.
        xy: The point being labelled, in data coordinates.
        xytext: Where the text sits, in data coordinates.
        colour: Leader-line colour; defaults to the muted ink.
    """
    axis.annotate(
        text,
        xy=xy,
        xytext=xytext,
        fontsize=ANNOTATION_PT,
        color=palette.ink,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "-", "color": colour or palette.muted, "linewidth": HAIRLINE_PT, "shrinkA": 0.0, "shrinkB": 2.0},
        path_effects=[withStroke(linewidth=LABEL_HALO_PT, foreground=palette.surface)],
        zorder=10,
    )


def _drawing(figure: Figure, axes: Sequence[Any], labels: Sequence[str] = ()) -> ToolpathDrawing:
    """Wrap a hand-built chart as a `ToolpathDrawing` so `save_drawing` accepts it.

    Args:
        figure: The figure.
        axes: Its panels.
        labels: Legend labels actually drawn.

    Returns:
        The drawing.
    """
    return ToolpathDrawing(figure=figure, axes=tuple(axes), legend_labels=tuple(labels))


def _write(drawing: ToolpathDrawing, out_dir: Path, name: str, formats: Sequence[str], theme: Theme) -> Tuple[Path, ...]:
    """Write one drawing, with LIVE TEXT in the SVG.

    matplotlib's default `svg.fonttype` is ``path``, which converts every glyph
    to outlines: the text is then unselectable, unsearchable, and immune to a
    reader's font settings. ``none`` keeps real `<text>` elements. It is set here
    in a context rather than as a global rcParam so this module cannot change how
    anything else in the process draws.

    Args:
        drawing: The drawing to write.
        out_dir: Directory to write into; created if absent.
        name: Base file stem, before the theme suffix.
        formats: File extensions, one file each.
        theme: The theme it was drawn for.

    Returns:
        The paths written.
    """
    from benchmarks.figures import themed_name

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = themed_name(name, theme)
    with rc_context({"svg.fonttype": "none"}):
        return tuple(save_drawing(drawing, out_dir / f"{stem}.{suffix}") for suffix in formats)


def _cumulative_samples(survey: PathSurvey) -> Tuple[List[float], List[float]]:
    """Engagement against cumulative arc length over every cut motion.

    Args:
        survey: The measured path.

    Returns:
        ``(distance, engagement_deg)``, in travel order.
    """
    distances: List[float] = []
    engagements: List[float] = []
    travelled = 0.0
    for motion in survey.motions:
        for sample in motion.samples:
            distances.append(travelled + sample.distance)
            engagements.append(sample.engagement_deg)
        travelled += motion.length
    return distances, engagements


def _motion_starts(survey: PathSurvey) -> List[float]:
    """Cumulative arc length at which each cut motion begins.

    Args:
        survey: The measured path.

    Returns:
        One distance per motion, in travel order.
    """
    starts: List[float] = []
    travelled = 0.0
    for motion in survey.motions:
        starts.append(travelled)
        travelled += motion.length
    return starts


# ---------------------------------------------------------------------------
# 1. The corner defect.
# ---------------------------------------------------------------------------


def corner_motions(survey: PathSurvey, spec: PocketSpec) -> Tuple[MotionQuality, MotionQuality, MotionQuality]:
    """The degenerate loop, the slotting link, and the proper trochoid beside them.

    Found by SEARCH rather than by index, so the figure annotates whatever the
    current generator emits at that corner instead of three operation numbers
    that a regeneration can invalidate.

    Args:
        survey: The measured path.
        spec: The pocket, whose top-right corner is the subject.

    Returns:
        ``(degenerate loop, following link, following loop)``.

    Raises:
        NoCornerDefectError: The path emits no degenerate loop near that corner,
            which would mean the defect this figure exists to show is gone.
    """
    from benchmarks.errors import NoCornerDefectError

    corner_x = max(float(point[0]) for point in spec.polygon.points)
    corner_y = max(float(point[1]) for point in spec.polygon.points)
    degenerate = [
        motion
        for motion in survey.motions
        if motion.loop_radius is not None
        and motion.loop_radius <= spec.tool_radius
        and math.hypot(motion.start[0] - corner_x, motion.start[1] - corner_y) < spec.tool_diameter * CORNER_ZOOM_MARGIN_TOOL_DIAMETERS
    ]
    if not degenerate:
        raise NoCornerDefectError(
            f"{spec.name}: no degenerate loop within {CORNER_ZOOM_MARGIN_TOOL_DIAMETERS:g} tool diameters of the top-right corner; this figure has nothing to annotate."
        )
    loop = degenerate[0]
    following = [motion for motion in survey.motions if motion.index > loop.index]
    link = next(motion for motion in following if motion.loop_radius is None)
    trochoid = next(motion for motion in following if motion.loop_radius is not None and motion.loop_radius > spec.tool_radius)
    return loop, link, trochoid


def draw_corner_defect(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """Zoom on the corner where the generator plunges and then slots out.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.

    Raises:
        NoCornerDefectError: The corner carries no degenerate loop to annotate.
    """
    palette = palette_for(theme)
    spec = measured.spec
    loop, link, trochoid = corner_motions(measured.survey, spec)

    # A square window, so an equal aspect can be satisfied by the BOX and the
    # limits are honoured. Asking for an equal aspect on a non-square window with
    # `adjustable="datalim"` silently widens the view instead, which is how the
    # first draft of this figure quietly showed a whole quadrant.
    radius = trochoid.loop_radius or 0.0
    interest_x = (loop.start[0], link.start[0], link.end[0], trochoid.start[0] - radius, trochoid.start[0] + radius)
    interest_y = (loop.start[1], link.start[1], link.end[1], trochoid.start[1] - radius, trochoid.start[1] + radius)
    centre_x = 0.5 * (min(interest_x) + max(interest_x))
    centre_y = 0.5 * (min(interest_y) + max(interest_y))
    half = 0.5 * max(max(interest_x) - min(interest_x), max(interest_y) - min(interest_y)) + CORNER_ZOOM_MARGIN_TOOL_DIAMETERS * spec.tool_diameter

    figure = Figure(figsize=(CHART_WIDTH_IN, CORNER_FIGURE_HEIGHT_IN), dpi=FIGURE_DPI, facecolor=palette.surface)
    FigureCanvasAgg(figure)
    axis = figure.add_axes(CORNER_PANEL_BOX)
    draw_toolpath(
        measured.result,
        boundary=spec.polygon,
        holes=list(spec.holes),
        colour_by=ColourBy.OPERATION,
        tool_diameter=spec.tool_diameter,
        ax=axis,
        theme=theme,
    )
    axis.set_xlim(centre_x - half, centre_x + half)
    axis.set_ylim(centre_y - half, centre_y + half)
    axis.set_aspect("equal", adjustable="box")
    axis.set_axis_off()

    # The degenerate loop is 2% of the tool radius across: at the scale that
    # shows its neighbours it is one pixel. The inset is the only way to show
    # that it IS a loop and not a point, and it sits outside the pocket wall
    # where the main panel has nothing to occlude.
    inset = axis.inset_axes(CORNER_INSET_BOX)
    draw_toolpath(
        measured.result,
        boundary=spec.polygon,
        holes=list(spec.holes),
        colour_by=ColourBy.OPERATION,
        tool_diameter=spec.tool_diameter,
        ax=inset,
        theme=theme,
    )
    inset_half = CORNER_INSET_LOOP_RADII * (loop.loop_radius or spec.tool_radius)
    inset.set_xlim(loop.start[0] - inset_half, loop.start[0] + inset_half)
    inset.set_ylim(loop.start[1] - inset_half, loop.start[1] + inset_half)
    inset.set_aspect("equal", adjustable="box")
    # `draw_toolpath` leaves every panel with `set_axis_off()`, which clears the
    # axes' `axison` flag and suppresses the patch and the spines WHATEVER their
    # own visibility says -- so the inset drew as a transparent hole with the
    # main panel showing through. Turn the axis back on, then hide only the
    # ticks.
    inset.set_axis_on()
    inset.set_xticks([])
    inset.set_yticks([])
    inset.patch.set_visible(True)
    inset.set_facecolor(palette.surface)
    inset.set_zorder(9)
    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_color(palette.grid)
        spine.set_linewidth(HAIRLINE_PT)
    # A halo in the surface colour: the pocket wall passes behind the inset, and
    # a magnification label sitting on a boundary line reads as neither.
    inset.text(
        0.04,
        0.96,
        f"×{half / inset_half:.0f}",
        transform=inset.transAxes,
        fontsize=ANNOTATION_PT,
        color=palette.secondary,
        ha="left",
        va="top",
        zorder=Z_LABEL,
        bbox={"facecolor": palette.surface, "edgecolor": "none", "pad": LABEL_HALO_PAD_PT, "alpha": LABEL_HALO_ALPHA},
    )

    # `draw_toolpath` gives every axes it is handed its own legend, and this
    # figure hands it two. Two legends stacked on one panel is worse than none,
    # so both are lifted into a single one at the foot of the figure.
    handles, labels = _strip_axes_legends((axis, inset))
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(labels),
        fontsize=LEGEND_PT,
        frameon=False,
        bbox_to_anchor=(0.5, 0.0),
        labelcolor=palette.secondary,
    )

    loop_radius = loop.loop_radius or 0.0
    trochoid_radius = trochoid.loop_radius or 0.0
    callouts: Tuple[Tuple[Tuple[float, float], str, str], ...] = (
        (
            loop.start,
            palette.secondary,
            f"degenerate loop\nρ={loop_radius:.4f}, {loop_radius / spec.tool_radius:.1%} of tool radius\nlength {loop.length:.3f}, {loop.peak_engagement_deg:.0f}° — a plunge",
        ),
        (
            (0.5 * (link.start[0] + link.end[0]), 0.5 * (link.start[1] + link.end[1])),
            palette.secondary,
            f"straight {link.operation.value} between loops\nlength {link.length:.3f}, {link.peak_engagement_deg:.0f}°\ncutting into the corner — a slot",
        ),
        (
            trochoid.start,
            palette.secondary,
            f"proper trochoid\nρ={trochoid_radius:.4f}, {trochoid_radius / spec.tool_radius:.2f}× tool radius\nlength {trochoid.length:.3f}, {trochoid.peak_engagement_deg:.0f}°",
        ),
    )
    for (target, colour, text), fraction in zip(callouts, CORNER_LABEL_FRACTIONS):
        axis.annotate(
            text,
            xy=target,
            xycoords="data",
            xytext=(CORNER_LABEL_X_FRACTION, fraction),
            textcoords="axes fraction",
            fontsize=ANNOTATION_PT,
            color=palette.ink,
            ha="left",
            va="center",
            arrowprops={"arrowstyle": "-", "color": colour, "linewidth": HAIRLINE_PT, "shrinkA": 0.0, "shrinkB": 1.0},
            annotation_clip=False,
            zorder=10,
        )

    _title(
        figure,
        palette,
        "A plunge and a slot, in the corner of a 20×12 pocket",
        f"{spec.name}, tool ⌀{spec.tool_diameter:g}, {spec.tea_cap_deg:.0f}° cap — tool-centre path, top-right corner",
    )
    return _drawing(figure, (axis,))


def write_corner_defect(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the corner-defect figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_corner_defect(measured_path(), theme=theme), out_dir, CORNER_DEFECT_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 2. Coverage and residue.
# ---------------------------------------------------------------------------


def residual_points(spec: PocketSpec, stock: Stock, *, grid: int = COVERAGE_GRID_SAMPLES) -> Tuple[List[float], List[float]]:
    """Grid samples inside the reachable material that the path did not remove.

    Uses `benchmarks.coverage`'s own grid, so the picture shades exactly the
    samples the metric counts. A figure drawn on a different grid from the number
    beside it is two measurements pretending to be one.

    Args:
        spec: The pocket and tool.
        stock: The depleted stock the path left.
        grid: Samples along the pocket's longer bounding-box side.

    Returns:
        ``(xs, ys)`` of the uncut reachable samples.
    """
    xs_axis, ys_axis, _cell = _grid_axes(spec, grid)
    region = _coverage_2.ReachableDomain2(
        _polygon_to_ccw_vertices(spec.polygon),
        [_polygon_to_ccw_vertices(hole) for hole in spec.holes],
        spec.tool_radius,
    ).reachable_material()
    xs: List[float] = []
    ys: List[float] = []
    for x in xs_axis:
        for y in ys_axis:
            if region.contains(x, y) and stock.contains(x, y):
                xs.append(x)
                ys.append(y)
    return xs, ys


def draw_coverage_residual(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """The residue the metric detects, on a good path and a deliberately bad one.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.
    """
    palette = palette_for(theme)
    spec = measured.spec
    coarse = trochoidal_mat_toolpath_circular(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        stepover=UNDERCOVERED_STEPOVER_TOOL_DIAMETERS * spec.tool_diameter,
        holes=list(spec.holes),
        clearance_z=2.0 * spec.tool_diameter,
    )
    coarse_survey = survey_path(spec, coarse)
    panels = (
        (f"engagement-controlled, {spec.tea_cap_deg:.0f}° cap", measured.survey),
        (f"deliberately under-covered, constant {UNDERCOVERED_STEPOVER_TOOL_DIAMETERS:g}⌀ stepover", coarse_survey),
    )

    figure = Figure(figsize=(CHART_WIDTH_IN, RESIDUAL_FIGURE_HEIGHT_IN), dpi=FIGURE_DPI, facecolor=palette.surface)
    FigureCanvasAgg(figure)
    axes = figure.subplots(2, 1, squeeze=False)[:, 0]
    ring = [(float(p[0]), float(p[1])) for p in spec.polygon.points]
    ring.append(ring[0])
    for axis, (label, survey) in zip(axes, panels):
        estimate = measure_coverage(spec, survey.final_stock)
        xs, ys = residual_points(spec, survey.final_stock)
        axis.set_facecolor(palette.surface)
        axis.plot([p[0] for p in ring], [p[1] for p in ring], color=palette.ink, linewidth=REFERENCE_WIDTH_PT, zorder=2)
        axis.scatter(xs, ys, s=RESIDUAL_MARKER_AREA_PT2, c=palette.link, marker="s", linewidths=0.0, zorder=3)
        axis.set_aspect("equal", adjustable="datalim")
        axis.set_axis_off()
        axis.set_title(f"{label}   ·   {estimate.uncut_fraction:.3%} of the reachable material left", fontsize=SUBTITLE_PT, color=palette.secondary, loc="left", pad=4.0)
    _title(
        figure,
        palette,
        "What the coverage metric detects",
        f"{spec.name}, tool ⌀{spec.tool_diameter:g} — every square is a grid sample inside the reachable region that still holds material",
    )
    figure.subplots_adjust(top=0.82, bottom=0.02, left=0.02, right=0.98, hspace=0.22)
    return _drawing(figure, tuple(axes), ("uncut reachable sample",))


def write_coverage_residual(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the coverage-residual figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_coverage_residual(measured_path(), theme=theme), out_dir, COVERAGE_RESIDUAL_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 3. Engagement along the path.
# ---------------------------------------------------------------------------


def draw_engagement_along_path(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """Engagement against arc length, with the cap and every exceedance marked.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.
    """
    palette = palette_for(theme)
    cap = measured.spec.tea_cap_deg
    distances, engagements = _cumulative_samples(measured.survey)
    starts = _motion_starts(measured.survey)
    over = [(start + 0.5 * motion.length, motion.peak_engagement_deg) for start, motion in zip(starts, measured.survey.motions) if motion.cap_exceeded]

    figure, axes = _chart_figure(palette)
    axis = axes[0]
    axis.plot(distances, engagements, color=palette.ramp[2], linewidth=HAIRLINE_PT, zorder=3)
    axis.axhline(cap, color=palette.ink, linewidth=REFERENCE_WIDTH_PT, linestyle=(0, (5, 3)), zorder=4)
    axis.scatter([d for d, _ in over], [e for _, e in over], s=18.0, facecolors="none", edgecolors=palette.link, linewidths=1.1, zorder=5)
    axis.set_xlabel("cumulative cut length", fontsize=LABEL_PT, color=palette.secondary)
    axis.set_ylabel("engagement (°)", fontsize=LABEL_PT, color=palette.secondary)
    axis.set_xlim(0.0, max(distances) if distances else 1.0)
    axis.set_ylim(0.0, max(max(engagements, default=cap), cap) * 1.06)
    handles = [
        Line2D([], [], color=palette.ramp[2], linewidth=CURVE_WIDTH_PT),
        Line2D([], [], color=palette.ink, linewidth=REFERENCE_WIDTH_PT, linestyle=(0, (5, 3))),
        Line2D([], [], color=palette.link, marker="o", markerfacecolor="none", linestyle="none", markersize=5.0),
    ]
    _legend(axis, palette, handles, ["engagement at each probe", f"{cap:.0f}° cap", "motion over the cap"])
    _title(
        figure,
        palette,
        "Engagement along the path, and where it breaks the cap",
        f"{measured.spec.name}, tool ⌀{measured.spec.tool_diameter:g} — {len(over)} of {len(measured.survey.motions)} cut motions exceed the cap",
    )
    figure.subplots_adjust(top=0.80, bottom=0.14, left=0.10, right=0.98)
    return _drawing(figure, axes, ("engagement at each probe", f"{cap:.0f}° cap", "motion over the cap"))


def write_engagement_along_path(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the engagement-trace figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_engagement_along_path(measured_path(), theme=theme), out_dir, ENGAGEMENT_ALONG_PATH_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 4. Chip thinning.
# ---------------------------------------------------------------------------


def operating_band(survey: PathSurvey) -> Tuple[float, float]:
    """Where the generator's ENGAGED cutting actually sits, in degrees.

    Measured as the interquartile range of the engaged probes weighted by the cut
    length they stand for, so the band is where the tool spends its time rather
    than the extremes it touches.

    Args:
        survey: The measured path.

    Returns:
        ``(low_deg, high_deg)``; ``(0, 0)`` when nothing is engaged.
    """
    weighted = sorted(
        (sample.engagement_deg, motion.length / len(motion.samples)) for motion in survey.motions if motion.samples for sample in motion.samples if sample.engagement_deg > 0.0
    )
    total = sum(weight for _, weight in weighted)
    if total <= 0.0:
        return (0.0, 0.0)
    bounds: List[float] = []
    for target in (0.25 * total, 0.75 * total):
        cumulative = 0.0
        for value, weight in weighted:
            cumulative += weight
            if cumulative >= target:
                bounds.append(value)
                break
    return (bounds[0], bounds[-1]) if len(bounds) == 2 else (0.0, 0.0)


def draw_chip_thinning(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT, material: Optional[MaterialModel] = None) -> ToolpathDrawing:
    """Why engagement angle is a proxy and chip thickness is the load.

    Args:
        measured: The generated and surveyed path, for the operating band.
        theme: Which surface the figure is drawn for.
        material: The coefficients the rubbing floor comes from; the documented
            default when omitted, which the figure says on its face.

    Returns:
        The drawing.
    """
    palette = palette_for(theme)
    model = MaterialModel.build() if material is None else material
    floor = model.min_chip_thickness_mm / model.feed_per_tooth_mm
    low, high = operating_band(measured.survey)

    angles = [FULL_TURN_DEG * index / (CHIP_CURVE_SAMPLES - 1) for index in range(CHIP_CURVE_SAMPLES)]
    ratios = [chip_thickness_ratio_from_rim(angle) for angle in angles]

    figure, axes = _chart_figure(palette, height_in=CHART_HEIGHT_IN + 0.25)
    axis = axes[0]
    axis.axhspan(0.0, floor, color=palette.link, alpha=FLOOR_ALPHA, linewidth=0.0, zorder=1)
    axis.axvspan(UNREACHABLE_CAP_LOW_DEG, UNREACHABLE_CAP_HIGH_DEG, color=palette.muted, alpha=BAND_ALPHA, linewidth=0.0, zorder=1)
    if high > low:
        axis.axvspan(low, high, color=palette.cut, alpha=BAND_ALPHA, linewidth=0.0, zorder=1)
    axis.plot(angles, ratios, color=palette.ramp[3], linewidth=CURVE_WIDTH_PT, zorder=4)
    axis.axvline(CHIP_PLATEAU_RIM_DEG, color=palette.muted, linewidth=HAIRLINE_PT, linestyle=(0, (2, 3)), zorder=2)

    axis.set_xlim(0.0, FULL_TURN_DEG)
    axis.set_ylim(0.0, 1.08)
    axis.set_xticks([0, 60, 120, 180, 240, 300, 360])
    axis.set_xlabel("engaged rim arc θ (°) — the angle this kernel reports; a full turn is a surrounded cutter", fontsize=LABEL_PT, color=palette.secondary)
    axis.set_ylabel("$h_{ex}/f_z$  (chip thickness ÷ feed per tooth)", fontsize=LABEL_PT, color=palette.secondary)

    # A reparameterisation of the SAME x axis -- radial immersion a_e/D --
    # so it is a second scale on one quantity and not a second quantity.
    top = axis.secondary_xaxis("top", functions=(_immersion_from_angle, _angle_from_immersion))
    top.set_xlabel("radial immersion $a_e/D$", fontsize=LABEL_PT, color=palette.secondary)
    top.tick_params(colors=palette.secondary, labelsize=TICK_PT, length=0.0)
    for side in top.spines:
        top.spines[side].set_visible(False)

    _annotate(axis, palette, "plateau: past 180° of rim the chip\nis already at full thickness", xy=(280.0, 1.0), xytext=(196.0, 0.66))
    _annotate(axis, palette, f"rubbing floor $h_{{min}}/f_z$ = {floor:.3g}\nbelow this the edge ploughs", xy=(330.0, floor), xytext=(210.0, 0.18))
    axis.text(0.5 * (UNREACHABLE_CAP_LOW_DEG + UNREACHABLE_CAP_HIGH_DEG), 1.03, "caps we\ncannot deliver", fontsize=ANNOTATION_PT, color=palette.secondary, ha="center", va="top")
    if high > low:
        axis.text(
            0.5 * (low + high), 0.06, f"where this generator\nactually cuts\n{low:.0f}–{high:.0f}°", fontsize=ANNOTATION_PT, color=palette.secondary, ha="center", va="bottom"
        )

    _title(
        figure,
        palette,
        "Chip thinning: why the engagement cap is a proxy for the load",
        f"{CHIP_FORMULA} — $f_z$={model.feed_per_tooth_mm:g} mm, $h_{{min}}$={model.min_chip_thickness_mm:g} mm: documented DEFAULTS",
    )
    figure.subplots_adjust(top=0.74, bottom=0.13, left=0.11, right=0.97)
    return _drawing(figure, axes)


def _immersion_from_angle(degrees: Any) -> Any:
    """Radial immersion ``a_e/D = (1 - cos(θ_rim/2))/2`` from an engaged rim arc.

    Vectorised, because matplotlib hands a secondary axis' transform whole
    arrays of tick positions and `math.cos` accepts only scalars.

    Args:
        degrees: Engagement angle, scalar or array, in degrees.

    Returns:
        The immersion ratio, same shape as the input.
    """
    return 0.5 * (1.0 - np.cos(np.radians(np.clip(np.asarray(degrees, dtype=float), 0.0, FULL_TURN_DEG) / 2.0)))


def _angle_from_immersion(ratio: Any) -> Any:
    """The inverse of `_immersion_from_angle`, clamped to the domain of ``arccos``.

    Args:
        ratio: Radial immersion, scalar or array.

    Returns:
        The engagement angle in degrees, same shape as the input.
    """
    return 2.0 * np.degrees(np.arccos(np.clip(1.0 - 2.0 * np.asarray(ratio, dtype=float), -1.0, 1.0)))


def write_chip_thinning(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the chip-thinning figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_chip_thinning(measured_path(FIGURE_SMALL_POCKET), theme=theme), out_dir, CHIP_THINNING_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 5. The path coloured by engagement.
# ---------------------------------------------------------------------------


def per_operation_engagement(measured: MeasuredPath) -> List[Optional[float]]:
    """One peak engagement per operation, ``None`` where nothing was measured.

    Args:
        measured: The generated and surveyed path.

    Returns:
        A list as long as the operation stream.
    """
    peaks: List[Optional[float]] = [None] * len(measured.result.operations)
    for motion in measured.survey.motions:
        peaks[motion.index] = motion.peak_engagement_deg
    return peaks


def draw_engagement_map(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """Where on the pocket the load actually peaks.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.
    """
    peaks = [value for value in per_operation_engagement(measured) if value is not None]
    return draw_toolpath(
        measured.result,
        boundary=measured.spec.polygon,
        holes=list(measured.spec.holes),
        colour_by=ColourBy.ENGAGEMENT,
        engagement_deg=per_operation_engagement(measured),
        tool_diameter=measured.spec.tool_diameter,
        title="Where the load is, not just how big it gets",
        subtitle=f"{measured.spec.name}, tool ⌀{measured.spec.tool_diameter:g} — engagement {min(peaks):.0f}° to {max(peaks):.0f}° over {len(peaks)} cut motions",
        theme=theme,
    )


def write_engagement_map(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the engagement-map figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_engagement_map(measured_path(), theme=theme), out_dir, ENGAGEMENT_MAP_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 6. Curvature and the feed ceiling.
# ---------------------------------------------------------------------------


def draw_curvature_feed(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT, machine: Optional[MachineModel] = None) -> ToolpathDrawing:
    """Curvature along the path, and the feed the machine can hold over it.

    Two stacked panels sharing one x axis, never a dual y axis: curvature and
    feed are different quantities in different units, and overlaying them on one
    frame would invite a reader to compare two arbitrary scalings.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.
        machine: The limits the ceiling comes from; the documented default when
            omitted.

    Returns:
        The drawing.
    """
    palette = palette_for(theme)
    model = MachineModel.build() if machine is None else machine
    starts = _motion_starts(measured.survey)

    steps_x: List[float] = []
    curvature: List[float] = []
    ceiling: List[float] = []
    for start, motion in zip(starts, measured.survey.motions):
        for position in (start, start + motion.length):
            steps_x.append(position)
            curvature.append(motion.curvature)
            ceiling.append(model.curvature_limited_feed_mm_per_s(motion.curvature))

    figure, axes = _chart_figure(palette, panels=2, height_in=STACKED_CHART_HEIGHT_IN, share_x=True)
    top, bottom = axes
    top.plot(steps_x, curvature, color=palette.ramp[3], linewidth=HAIRLINE_PT, zorder=3)
    top.set_ylabel("curvature κ (1/mm)", fontsize=LABEL_PT, color=palette.secondary)
    top.set_yscale("symlog", linthresh=0.1)

    bottom.plot(steps_x, ceiling, color=palette.cut, linewidth=HAIRLINE_PT, zorder=3)
    bottom.axhline(model.feed_mm_per_s, color=palette.ink, linewidth=REFERENCE_WIDTH_PT, linestyle=(0, (5, 3)), zorder=4)
    bottom.set_ylabel("feed ceiling (mm/s)", fontsize=LABEL_PT, color=palette.secondary)
    bottom.set_xlabel("cumulative cut length", fontsize=LABEL_PT, color=palette.secondary)
    bottom.set_xlim(0.0, max(steps_x) if steps_x else 1.0)
    _legend(
        bottom,
        palette,
        [Line2D([], [], color=palette.cut, linewidth=CURVE_WIDTH_PT), Line2D([], [], color=palette.ink, linewidth=REFERENCE_WIDTH_PT, linestyle=(0, (5, 3)))],
        ["curvature-limited ceiling", f"programmed feed {model.feed_mm_per_s:.0f} mm/s"],
    )

    worst = min(ceiling) if ceiling else model.feed_mm_per_s
    _title(
        figure,
        palette,
        "Curvature sets a feed ceiling the machine cannot exceed",
        f"{FEED_FORMULA}, $a_{{max}}$={model.max_acceleration_mm_per_s2:g} mm/s² — tightest corner {worst:.1f} mm/s, {model.feed_mm_per_s / worst:.0f}× under feed",
    )
    figure.subplots_adjust(top=0.84, bottom=0.10, left=0.12, right=0.97, hspace=0.12)
    return _drawing(figure, axes, ("feed ceiling", "programmed feed"))


def write_curvature_feed(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the curvature/feed figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_curvature_feed(measured_path(), theme=theme), out_dir, CURVATURE_FEED_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 7. Length against feed-limited time.
# ---------------------------------------------------------------------------


def draw_length_vs_time(paths: Sequence[MeasuredPath], *, theme: Theme = Theme.LIGHT, machine: Optional[MachineModel] = None) -> ToolpathDrawing:
    """Path length beside the time a machine would actually take.

    Args:
        paths: The measured paths to compare, in plotting order.
        theme: Which surface the figure is drawn for.
        machine: The limits the time comes from; the documented default when
            omitted.

    Returns:
        The drawing.
    """
    from benchmarks.quality import machine_outcome

    palette = palette_for(theme)
    model = MachineModel.build() if machine is None else machine
    labels = [path.spec.name for path in paths]
    lengths = [path.survey.cut_length for path in paths]
    outcomes = [machine_outcome(path.survey, model) for path in paths]
    seconds = [outcome.cutting_seconds for outcome in outcomes]
    ideal = [length / model.feed_mm_per_s for length in lengths]

    figure, axes = _chart_figure(palette, panels=2, height_in=STACKED_CHART_HEIGHT_IN - 0.6)
    top, bottom = axes
    positions = list(range(len(labels)))
    top.barh(positions, lengths, color=palette.cut, height=0.55, zorder=3)
    top.set_yticks(positions)
    top.set_yticklabels(labels, fontsize=TICK_PT, color=palette.secondary)
    top.set_xlabel("cut length", fontsize=LABEL_PT, color=palette.secondary)
    top.invert_yaxis()
    for position, value in zip(positions, lengths):
        top.text(value, position, f" {value:,.0f}", fontsize=ANNOTATION_PT, color=palette.secondary, va="center", ha="left")

    width = 0.34
    # The ideal time is a REFERENCE, not a second series, so it is an outlined
    # bar rather than a second hue. `palette.muted` would have been the obvious
    # choice and is wrong: it is an ink token, it reads grey, and the palette
    # validator fails it on the chroma floor for exactly that reason.
    bottom.barh([p - width / 2 for p in positions], ideal, facecolor="none", edgecolor=palette.secondary, linewidth=HAIRLINE_PT, height=width, zorder=3)
    bottom.barh([p + width / 2 for p in positions], seconds, color=palette.cut, height=width, zorder=3)
    bottom.set_yticks(positions)
    bottom.set_yticklabels(labels, fontsize=TICK_PT, color=palette.secondary)
    bottom.set_xlabel("cutting time (s)", fontsize=LABEL_PT, color=palette.secondary)
    bottom.invert_yaxis()
    for position, value in zip(positions, seconds):
        bottom.text(value, position + width / 2, f" {value:.1f} s", fontsize=ANNOTATION_PT, color=palette.secondary, va="center", ha="left")
    handles = [Patch(facecolor="none", edgecolor=palette.secondary, linewidth=HAIRLINE_PT), Patch(facecolor=palette.cut)]
    _legend(bottom, palette, handles, ["length ÷ programmed feed", "feed-limited, curvature bound applied"], location="upper right")

    utilisation = [outcome.feed_utilisation for outcome in outcomes]
    band = f"{min(utilisation):.1%}–{max(utilisation):.1%}"
    _title(
        figure,
        palette,
        "Path length is a proxy; time under the machine's limits is the cost",
        f"$a_{{max}}$={model.max_acceleration_mm_per_s2:g} mm/s², feed {model.feed_mm_per_s:.0f} mm/s — utilisation {band}: the ceiling barely binds",
    )
    figure.subplots_adjust(top=0.80, bottom=0.12, left=0.17, right=0.94, hspace=0.55)
    return _drawing(figure, axes, ("length ÷ programmed feed", "feed-limited, curvature bound applied"))


def write_length_vs_time(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the length-against-time figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    paths = [measured_path(FIGURE_SMALL_POCKET), measured_path(FIGURE_POCKET)]
    return _write(draw_length_vs_time(paths, theme=theme), out_dir, LENGTH_VS_TIME_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 8. Time at engagement.
# ---------------------------------------------------------------------------


def draw_engagement_histogram(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """How much cutting happens at each load, not just the worst moment.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.
    """
    palette = palette_for(theme)
    quality = measure_quality(measured.spec, measured.result)
    bands = quality.longevity.engagement_length_histogram
    cap = measured.spec.tea_cap_deg

    figure, axes = _chart_figure(palette)
    axis = axes[0]
    positions = list(range(len(bands)))
    # The ramp encodes the band's magnitude; bands past the ramp's validated
    # length reuse its top step rather than inventing a sixth colour.
    colours = [palette.ramp[min(index, len(palette.ramp) - 1)] for index in positions]
    lengths = [length for _low, _high, length in bands]
    axis.bar(positions, lengths, color=colours, width=0.72, zorder=3)
    axis.set_xticks(positions)
    axis.set_xticklabels([f"{low:.0f}–{high:.0f}" for low, high, _length in bands], fontsize=TICK_PT, color=palette.secondary)
    axis.set_xlabel("engagement band (°)", fontsize=LABEL_PT, color=palette.secondary)
    axis.set_ylabel("cut length in band", fontsize=LABEL_PT, color=palette.secondary)
    for position, length in zip(positions, lengths):
        if length > 0.0:
            # A band holding a tenth of a unit is not "0": rounding it to the
            # nearest whole number would print a zero over a bar that exists,
            # and those are exactly the over-cap bands a reader is looking for.
            printed = f"{length:,.0f}" if length >= 1.0 else f"{length:.2f}"
            axis.text(position, length, printed, fontsize=ANNOTATION_PT, color=palette.secondary, ha="center", va="bottom")

    over_cap = sum(length for low, _high, length in bands if low >= cap)
    total = sum(lengths)
    _title(
        figure,
        palette,
        "Cumulative damage: how much cutting happens at each load",
        f"{measured.spec.name}, tool ⌀{measured.spec.tool_diameter:g} — {over_cap:,.1f} of {total:,.0f} units ({over_cap / total:.2%}) at or above the {cap:.0f}° cap",
    )
    figure.subplots_adjust(top=0.80, bottom=0.15, left=0.11, right=0.98)
    return _drawing(figure, axes)


def write_engagement_histogram(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the time-at-engagement figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_engagement_histogram(measured_path(), theme=theme), out_dir, ENGAGEMENT_HISTOGRAM_NAME, formats, theme)


# ---------------------------------------------------------------------------
# 9. Entries into material.
# ---------------------------------------------------------------------------


def draw_material_entries(measured: MeasuredPath, *, theme: Theme = Theme.LIGHT) -> ToolpathDrawing:
    """Every impact the tool takes, marked where it happens.

    Args:
        measured: The generated and surveyed path.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.
    """
    palette = palette_for(theme)
    entries = entry_cut_indices(measured.result)
    marked = [motion for motion in measured.survey.motions if motion.index in entries]
    worst_entry = max((motion.peak_engagement_deg for motion in marked), default=0.0)
    drawing = draw_toolpath(
        measured.result,
        boundary=measured.spec.polygon,
        holes=list(measured.spec.holes),
        colour_by=ColourBy.OPERATION,
        tool_diameter=measured.spec.tool_diameter,
        title="Every entry into material is an impact",
        subtitle=f"{measured.spec.name}, tool ⌀{measured.spec.tool_diameter:g} — {len(marked)} entries, worst {worst_entry:.0f}° engagement",
        theme=theme,
    )
    axis = drawing.axis
    axis.scatter(
        [motion.start[0] for motion in marked],
        [motion.start[1] for motion in marked],
        s=EVENT_MARKER_PT**2,
        facecolors="none",
        edgecolors=palette.link,
        linewidths=1.3,
        zorder=8,
    )
    # The ring is this figure's whole subject, so it belongs in the legend beside
    # the motion classes rather than only in the caption.
    handles, labels = _strip_figure_legends(drawing.figure)
    handles.append(Line2D([], [], color=palette.link, marker="o", markerfacecolor="none", linestyle="none", markersize=EVENT_MARKER_PT))
    labels.append("entry into material")
    drawing.figure.legend(handles, labels, loc="lower center", ncol=len(labels), fontsize=LEGEND_PT, frameon=False, bbox_to_anchor=(0.5, 0.0), labelcolor=palette.secondary)
    return drawing


def write_material_entries(out_dir: Path, *, formats: Sequence[str] = ("svg",), theme: Theme = Theme.LIGHT) -> Tuple[Path, ...]:
    """Draw and write the material-entries figure.

    Args:
        out_dir: Directory to write into.
        formats: File extensions, one file each.
        theme: Which surface the figure is drawn for.

    Returns:
        The paths written.
    """
    return _write(draw_material_entries(measured_path(), theme=theme), out_dir, MATERIAL_ENTRIES_NAME, formats, theme)


# Every quality figure, in publication order. A figure is published by adding its
# writer here and nowhere else, so one cannot reach the docs without being
# regenerable.
QualityWriter = Callable[..., Tuple[Path, ...]]

QUALITY_FIGURES: Tuple[Tuple[str, QualityWriter], ...] = (
    (CORNER_DEFECT_NAME, write_corner_defect),
    (COVERAGE_RESIDUAL_NAME, write_coverage_residual),
    (ENGAGEMENT_ALONG_PATH_NAME, write_engagement_along_path),
    (CHIP_THINNING_NAME, write_chip_thinning),
    (ENGAGEMENT_MAP_NAME, write_engagement_map),
    (CURVATURE_FEED_NAME, write_curvature_feed),
    (LENGTH_VS_TIME_NAME, write_length_vs_time),
    (ENGAGEMENT_HISTOGRAM_NAME, write_engagement_histogram),
    (MATERIAL_ENTRIES_NAME, write_material_entries),
)


def write_all_quality_figures(out_dir: Path, *, formats: Sequence[str] = ("svg",), themes: Sequence[Theme] = (Theme.LIGHT, Theme.DARK)) -> Tuple[Path, ...]:
    """Redraw every machining-quality figure, in both themes.

    Args:
        out_dir: Directory the figures are written to; created if absent.
        formats: File extensions, one file each per figure.
        themes: Which surfaces to draw for.

    Returns:
        Every path written, figure by figure in publication order.

    Warns:
        UnavoidableEngagementWarning: Raised through from the generators where
            the cap could not be honoured.
    """
    written: List[Path] = []
    for _name, write in QUALITY_FIGURES:
        for theme in themes:
            written.extend(write(out_dir, formats=formats, theme=theme))
    return tuple(written)


def _strip_figure_legends(figure: Figure) -> Tuple[List[Any], List[str]]:
    """Remove every FIGURE-level legend and return the entries of the first.

    `draw_toolpath` attaches its legend to the FIGURE when it builds its own,
    and to the AXES when it is handed one. A caller that wants to extend the
    legend has to find whichever it got, or it ends up drawing a second legend
    on top of the first.

    Args:
        figure: The figure to strip.

    Returns:
        ``(handles, labels)`` from the first legend removed; empty when none.
    """
    handles: List[Any] = []
    labels: List[str] = []
    for legend in list(figure.legends):
        if not handles:
            handles = list(legend.legend_handles)
            labels = [text.get_text() for text in legend.get_texts()]
        legend.remove()
    return handles, labels


def _strip_axes_legends(axes: Sequence[Any]) -> Tuple[List[Any], List[str]]:
    """Remove every per-axes legend and return the entries of the first one found.

    Args:
        axes: The panels to strip.

    Returns:
        ``(handles, labels)`` from the first legend removed; empty when none.
    """
    handles: List[Any] = []
    labels: List[str] = []
    for axis in axes:
        legend = axis.get_legend()
        if legend is None:
            continue
        if not handles:
            handles = list(legend.legend_handles)
            labels = [text.get_text() for text in legend.get_texts()]
        legend.remove()
    return handles, labels
