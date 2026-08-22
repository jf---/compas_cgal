"""Draw a generated tool path as a figure, on this corpus's validated palette.

This is the public drawing API: `draw_toolpath` for one path over its pocket,
`draw_comparison` for several over the same pocket as stacked equal-aspect
panels. Both return a `ToolpathDrawing` rather than a bare axes, so a caller gets
the legend it drew along with the figure it drew it on.

Three modules stand behind it, and the split is the point. `benchmarks.palette`
holds the design tokens and the evidence they passed; `benchmarks.pathgeometry`
turns a result into polylines, sampling every arc from its own parametrisation
rather than from the generator's tessellation; `benchmarks.marks` decides what
each motion looks like under the active colour mode. Only this module imports
matplotlib, and it builds its own Agg canvas, so a figure can be drawn and saved
with no display and without touching global backend state.

`ColourBy.ENGAGEMENT` draws measurements the caller supplies. Nothing in the
drawing path imports the kernel, so a figure cannot disagree with the audit that
produced its numbers, and a path that was replayed or deserialised draws exactly
like one just generated.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any
from typing import Dict
from typing import List
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Tuple

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from benchmarks.errors import AmbiguousPanelError
from benchmarks.errors import EmptyComparisonError
from benchmarks.errors import MissingEngagementDataError
from benchmarks.errors import MissingToolDiameterError
from benchmarks.marks import POINT_MARKER_PT

# Re-exported under its public spelling: the colour mode is chosen by callers of
# this module, and `benchmarks.marks` is where it is defined rather than where it
# is used. The `as` form is what marks a re-export as deliberate.
from benchmarks.marks import ColourBy as ColourBy
from benchmarks.marks import Mark
from benchmarks.marks import encode
from benchmarks.marks import legend_entries
from benchmarks.palette import Palette
from benchmarks.palette import Theme
from benchmarks.palette import palette_for
from benchmarks.pathgeometry import Motion
from benchmarks.pathgeometry import closed_rings
from benchmarks.pathgeometry import drawable_motions
from benchmarks.pathgeometry import view_bounds

# Stroke weight of the pocket boundary, in points.
BOUNDARY_WIDTH_PT = 1.1

# Layout, in inches. The figure is sized for a one-column page; a panel's height
# follows the pocket's own aspect ratio, so an equal-aspect drawing does not sit
# in a band of empty surface, and is clamped so that a very long or very square
# pocket still produces a usable page.
FIGURE_WIDTH_IN = 6.5
SIDE_MARGIN_IN = 0.12
TOP_MARGIN_IN = 0.10
BOTTOM_MARGIN_IN = 0.10
PANEL_HEADER_IN = 0.42
PANEL_GAP_IN = 0.18
LEGEND_ROW_IN = 0.22
LEGEND_PAD_IN = 0.12
SUPTITLE_HEIGHT_IN = 0.52
MIN_PANEL_HEIGHT_IN = 1.0
MAX_PANEL_HEIGHT_IN = 5.5

# Type sizes in points, and the resolution a raster export defaults to.
TITLE_PT = 10.5
SUBTITLE_PT = 8.5
PANEL_LABEL_PT = 9.5
PANEL_SUBTITLE_PT = 8.0
LEGEND_PT = 8.0
TRAVERSAL_LABEL_PT = 7.5
FIGURE_DPI = 200.0

# Line spacing as a multiple of type size, for stacking the title and subtitle.
TEXT_LEADING = 1.55

# Legend metrics. The column count is estimated from the widest label, in
# characters: `LEGEND_GLYPH_EM` is the mean advance of a sans glyph as a fraction
# of the type size, and `LEGEND_HANDLE_CHARS` charges the sample line and the
# gutter to the same character budget.
LEGEND_GLYPH_EM = 0.62
LEGEND_HANDLE_CHARS = 8
LEGEND_HANDLE_LEN_EM = 2.4
LEGEND_COLUMN_SPACING_EM = 1.4

# A traversal label sits on top of the strokes it names, so it carries a halo of
# the surface colour, kept slightly translucent so it dims rather than erases.
LABEL_HALO_PAD_PT = 1.0
LABEL_HALO_ALPHA = 0.85

# Points per inch, for converting a width in data units into a stroke weight.
POINTS_PER_INCH = 72.0

# Draw order. The swept envelope is a footprint under everything; the pocket
# boundary sits above it; travel below cutting; point events and labels on top.
Z_ENVELOPE = 1
Z_BOUNDARY = 2
Z_TRAVEL = 3
Z_CUT = 4
Z_POINT = 5
Z_LABEL = 6


@dataclass(frozen=True)
class ToolpathDrawing:
    """A drawn figure, and what a reader can be told about it.

    Attributes:
        figure: The matplotlib figure. It owns an Agg canvas unless the caller
            supplied their own axes, so `figure.savefig` works with no display
            and without touching global backend state.
        axes: One axes per geometry panel, in draw order. Always a tuple, even
            for a single panel, so a caller never has to test which shape it got.
        legend_labels: Every legend entry actually drawn, in legend order.
    """

    figure: Any
    axes: Tuple[Any, ...]
    legend_labels: Tuple[str, ...]

    @property
    def axis(self) -> Any:
        """The one geometry panel, for the single-panel case.

        Returns:
            The single axes.

        Raises:
            AmbiguousPanelError: The drawing has any number of panels but one.
        """
        if len(self.axes) != 1:
            raise AmbiguousPanelError(f"This drawing has {len(self.axes)} panels; read `.axes` and say which one.")
        return self.axes[0]


def draw_toolpath(
    result: Any,
    *,
    boundary: Any,
    holes: Sequence[Any] = (),
    colour_by: ColourBy = ColourBy.OPERATION,
    engagement_deg: Optional[Sequence[Optional[float]]] = None,
    annotate_traversals: bool = False,
    tool_diameter: Optional[float] = None,
    show_tool_envelope: bool = False,
    ax: Optional[Any] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    theme: Theme = Theme.LIGHT,
) -> ToolpathDrawing:
    """Draw one tool path over its pocket.

    Args:
        result: Anything carrying `.operations`; each operation is read for
            `.geometry`, `.operation` and `.path_index`.
        boundary: The pocket's outer boundary, as a compas polygon or any
            sequence of points.
        holes: Island boundaries, in the same form.
        colour_by: What the stroke colour encodes.
        engagement_deg: One engagement angle per operation, or None for an
            operation carrying no measurement. Required by, and used only by,
            `ColourBy.ENGAGEMENT`.
        annotate_traversals: Label each chain at its first cutting move with its
            path index.
        tool_diameter: Cutter diameter, in the geometry's own units. Required by
            `show_tool_envelope`, and otherwise used only to pad the view.
        show_tool_envelope: Draw the swept width of the tool under the cutting
            moves, rather than the centre path alone.
        ax: Draw into these axes instead of building a figure.
        title: Figure title, in primary ink.
        subtitle: One line under the title, in secondary ink.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing.

    Raises:
        EmptyToolpathError: The result carries no operations.
        MissingEngagementDataError: `ColourBy.ENGAGEMENT` without measurements.
        EngagementLengthMismatchError: Not exactly one measurement per operation.
        MissingToolDiameterError: An envelope was asked for without a diameter.
        UnknownColourModeError: The colour mode is not one this module draws.
        UnknownOperationClassError: An operation names an unknown class.
        UnplottableGeometryError: An operation carries a primitive with no path.
        UnplottableBoundaryError: A boundary has fewer than three points.
    """
    palette = palette_for(theme)
    motions = drawable_motions(result)
    marks = encode(motions, colour_by, palette, engagement_deg)
    entries = legend_entries(marks)
    rings = closed_rings(boundary, holes)
    bounds = view_bounds(rings, motions, tool_diameter if show_tool_envelope else None)

    if ax is None:
        figure, axes = _new_figure(
            panels=1,
            bounds=bounds,
            palette=palette,
            has_panel_header=False,
            has_suptitle=title is not None or subtitle is not None,
            legend_rows=_legend_layout(tuple(entries))[1],
        )
        panel = axes[0]
    else:
        figure, panel = ax.figure, ax
        _ensure_canvas(figure)

    _draw_panel(
        panel,
        motions=motions,
        marks=marks,
        rings=rings,
        bounds=bounds,
        palette=palette,
        annotate_traversals=annotate_traversals,
        tool_diameter=tool_diameter,
        show_tool_envelope=show_tool_envelope,
    )
    if ax is None:
        _draw_headings(figure, title, subtitle, palette)
        labels = _draw_figure_legend(figure, entries, palette)
    else:
        _draw_panel_heading(panel, title, subtitle, palette)
        labels = _draw_axes_legend(panel, entries, palette)
    return ToolpathDrawing(figure=figure, axes=(panel,), legend_labels=labels)


def draw_comparison(
    panels: Mapping[str, Any],
    *,
    boundary: Any,
    holes: Sequence[Any] = (),
    colour_by: ColourBy = ColourBy.OPERATION,
    engagement_deg: Optional[Mapping[str, Sequence[Optional[float]]]] = None,
    annotate_traversals: bool = False,
    tool_diameter: Optional[float] = None,
    show_tool_envelope: bool = False,
    panel_subtitles: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    theme: Theme = Theme.LIGHT,
) -> ToolpathDrawing:
    """Draw several tool paths over the same pocket, as stacked equal-aspect panels.

    Every panel shares one view and one legend, which is the whole point: two
    paths drawn at two scales are not comparable however they are labelled. The
    panels are stacked rather than placed side by side so that the pocket is as
    wide as the page in each of them.

    Args:
        panels: Label to result, in the order the panels are stacked.
        boundary: The pocket's outer boundary, shared by every panel.
        holes: Island boundaries.
        colour_by: What the stroke colour encodes, in every panel.
        engagement_deg: Panel label to that panel's per-operation measurements.
            Required by, and used only by, `ColourBy.ENGAGEMENT`.
        annotate_traversals: Label each chain at its first cutting move.
        tool_diameter: Cutter diameter, in the geometry's own units.
        show_tool_envelope: Draw the swept width of the tool.
        panel_subtitles: Label to one line of secondary text beside that panel's
            label; a panel with no entry gets none.
        title: Figure title, in primary ink.
        subtitle: One line under the title, in secondary ink.
        theme: Which surface the figure is drawn for.

    Returns:
        The drawing, whose `axes` holds one panel per entry in `panels`.

    Raises:
        EmptyComparisonError: No panels were given.
        EmptyToolpathError: A result carries no operations.
        MissingEngagementDataError: `ColourBy.ENGAGEMENT` without a panel's
            measurements.
        EngagementLengthMismatchError: Not exactly one measurement per operation.
        MissingToolDiameterError: An envelope was asked for without a diameter.
        UnknownColourModeError: The colour mode is not one this module draws.
        UnknownOperationClassError: An operation names an unknown class.
        UnplottableGeometryError: An operation carries a primitive with no path.
        UnplottableBoundaryError: A boundary has fewer than three points.
    """
    if not panels:
        raise EmptyComparisonError("A comparison needs at least one panel.")

    palette = palette_for(theme)
    rings = closed_rings(boundary, holes)
    per_panel = {label: drawable_motions(result) for label, result in panels.items()}
    per_panel_marks = {label: encode(motions, colour_by, palette, _panel_engagement(colour_by, label, engagement_deg)) for label, motions in per_panel.items()}
    entries: Dict[str, Mark] = {}
    for marks in per_panel_marks.values():
        entries.update({label: mark for label, mark in legend_entries(marks).items() if label not in entries})

    every_motion = [motion for motions in per_panel.values() for motion in motions]
    bounds = view_bounds(rings, every_motion, tool_diameter if show_tool_envelope else None)
    figure, axes = _new_figure(
        panels=len(panels),
        bounds=bounds,
        palette=palette,
        has_panel_header=True,
        has_suptitle=title is not None or subtitle is not None,
        legend_rows=_legend_layout(tuple(entries))[1],
    )

    for axis, (label, motions) in zip(axes, per_panel.items()):
        _draw_panel(
            axis,
            motions=motions,
            marks=per_panel_marks[label],
            rings=rings,
            bounds=bounds,
            palette=palette,
            annotate_traversals=annotate_traversals,
            tool_diameter=tool_diameter,
            show_tool_envelope=show_tool_envelope,
        )
        _draw_panel_heading(axis, label, None if panel_subtitles is None else panel_subtitles.get(label), palette)

    _draw_headings(figure, title, subtitle, palette)
    labels = _draw_figure_legend(figure, entries, palette)
    return ToolpathDrawing(figure=figure, axes=tuple(axes), legend_labels=labels)


def _panel_engagement(
    colour_by: ColourBy,
    label: str,
    engagement_deg: Optional[Mapping[str, Sequence[Optional[float]]]],
) -> Optional[Sequence[Optional[float]]]:
    """One panel's measurements, checked for presence before the panel is drawn.

    Args:
        colour_by: The active colour mode.
        label: The panel's label.
        engagement_deg: The mapping the caller supplied, if any.

    Returns:
        The panel's measurements, or None when the mode does not need them.

    Raises:
        MissingEngagementDataError: The mode needs them and this panel has none.
    """
    if colour_by is not ColourBy.ENGAGEMENT:
        return None
    if engagement_deg is None or label not in engagement_deg:
        raise MissingEngagementDataError(f"ColourBy.ENGAGEMENT needs measurements for every panel; panel {label!r} has none.")
    return engagement_deg[label]


def _draw_panel(
    axis: Any,
    *,
    motions: Sequence[Motion],
    marks: Mapping[int, Mark],
    rings: Sequence[Sequence[Tuple[float, float]]],
    bounds: Tuple[float, float, float, float],
    palette: Palette,
    annotate_traversals: bool,
    tool_diameter: Optional[float],
    show_tool_envelope: bool,
) -> None:
    """Draw one panel.

    Args:
        axis: The axes to draw into.
        motions: The motions to draw.
        marks: Their marks.
        rings: The pocket rings.
        bounds: The shared view.
        palette: The theme's tokens.
        annotate_traversals: Label each chain at its first cutting move.
        tool_diameter: Cutter diameter, in the geometry's own units.
        show_tool_envelope: Draw the swept width of the tool.

    Raises:
        MissingToolDiameterError: An envelope was asked for without a diameter.
    """
    _prepare_axes(axis, bounds, palette)
    if show_tool_envelope:
        _draw_envelope(axis, motions, tool_diameter, palette)
    for ring in rings:
        axis.plot(
            [x for x, _ in ring],
            [y for _, y in ring],
            color=palette.ink,
            linewidth=BOUNDARY_WIDTH_PT,
            solid_joinstyle="round",
            zorder=Z_BOUNDARY,
        )

    for motion in motions:
        mark = marks[motion.index]
        axis.plot(
            [x for x, _ in motion.points],
            [y for _, y in motion.points],
            color=mark.colour,
            linewidth=mark.width_pt,
            linestyle=mark.linestyle,
            marker=mark.marker or "none",
            markersize=POINT_MARKER_PT,
            markeredgewidth=0.0,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=_motion_zorder(motion),
        )
    if annotate_traversals:
        _draw_traversal_labels(axis, motions, palette)


def _motion_zorder(motion: Motion) -> int:
    """Where a motion sits in the draw order.

    Args:
        motion: The motion.

    Returns:
        Its z order: point events above cutting, cutting above travel.
    """
    if motion.is_point:
        return Z_POINT
    return Z_CUT if motion.is_cutting else Z_TRAVEL


def _prepare_axes(axis: Any, bounds: Tuple[float, float, float, float], palette: Palette) -> None:
    """Give a panel an equal aspect, the shared view, and no axis furniture.

    The aspect is applied immediately rather than left to the first draw. An
    equal aspect on a fixed box works by widening the data limits, so until it
    has been applied the data transform is wrong -- and the swept envelope, whose
    width is a length in data units, would be sized from that wrong transform and
    only corrected on a second render.

    Args:
        axis: The axes.
        bounds: The shared view.
        palette: The theme's tokens.
    """
    xmin, xmax, ymin, ymax = bounds
    axis.set_facecolor(palette.surface)
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.set_aspect("equal", adjustable="datalim")
    axis.apply_aspect()
    axis.set_axis_off()


def _draw_envelope(axis: Any, motions: Sequence[Motion], tool_diameter: Optional[float], palette: Palette) -> None:
    """Draw the swept width of the tool under the moves that cut at depth.

    The width is a length in the geometry's units, not a stroke weight, so it is
    recomputed from the axes transform on every draw and stays true through a
    resize or an export at another size.

    Args:
        axis: The axes.
        motions: The motions. Only cutting moves sweep material at depth, and a
            plunge bores the tool disk at its foot.
        tool_diameter: Cutter diameter, in the geometry's own units.
        palette: The theme's tokens.

    Raises:
        MissingToolDiameterError: No diameter was supplied.
    """
    if tool_diameter is None:
        raise MissingToolDiameterError("show_tool_envelope draws the swept width of the tool; pass tool_diameter=... so that there is a width to draw.")

    strokes: List[Any] = []
    disks: List[Any] = []
    for motion in motions:
        bores = motion.name == "plunge"
        if not (motion.is_cutting or bores) or (motion.is_point and not bores):
            continue
        artists = axis.plot(
            [x for x, _ in motion.points],
            [y for _, y in motion.points],
            color=palette.grid,
            linestyle="none" if motion.is_point else "-",
            marker="o" if motion.is_point else "none",
            markeredgewidth=0.0,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=Z_ENVELOPE,
        )
        (disks if motion.is_point else strokes).extend(artists)
    _bind_data_width(axis, strokes, disks, tool_diameter)


def _bind_data_width(axis: Any, strokes: Sequence[Any], disks: Sequence[Any], width: float) -> None:
    """Keep artists as wide as *width* in DATA units, across every redraw.

    Args:
        axis: The axes the artists live in.
        strokes: Artists whose linewidth carries the width.
        disks: Artists whose marker size carries the width.
        width: The width in data units.
    """

    def _resize(_event: Any = None) -> None:
        points = _data_units_in_points(axis, width)
        for artist in strokes:
            artist.set_linewidth(points)
        for artist in disks:
            artist.set_markersize(points)

    figure = axis.figure
    _ensure_canvas(figure)
    figure.canvas.mpl_connect("draw_event", _resize)
    _resize()


def _data_units_in_points(axis: Any, width: float) -> float:
    """Convert a length in data units into a stroke weight in points.

    Args:
        axis: The axes whose transform defines the scale.
        width: The length in data units.

    Returns:
        The equivalent number of points at the figure's current size.
    """
    origin = axis.transData.transform((0.0, 0.0))
    offset = axis.transData.transform((width, 0.0))
    pixels = math.hypot(float(offset[0]) - float(origin[0]), float(offset[1]) - float(origin[1]))
    return pixels * POINTS_PER_INCH / float(axis.figure.dpi)


def _draw_traversal_labels(axis: Any, motions: Sequence[Motion], palette: Palette) -> None:
    """Label every chain at its first cutting move.

    Args:
        axis: The axes.
        motions: The motions, in stream order.
        palette: The theme's tokens.
    """
    labelled: List[int] = []
    for motion in motions:
        if not motion.is_cutting or motion.is_point or motion.path_index in labelled:
            continue
        labelled.append(motion.path_index)
        axis.text(
            motion.points[0][0],
            motion.points[0][1],
            str(motion.path_index),
            color=palette.ink,
            fontsize=TRAVERSAL_LABEL_PT,
            ha="center",
            va="center",
            zorder=Z_LABEL,
            bbox={"facecolor": palette.surface, "edgecolor": "none", "pad": LABEL_HALO_PAD_PT, "alpha": LABEL_HALO_ALPHA},
        )


# --- figure furniture ---------------------------------------------------------


def _new_figure(
    *,
    panels: int,
    bounds: Tuple[float, float, float, float],
    palette: Palette,
    has_panel_header: bool,
    has_suptitle: bool,
    legend_rows: int,
) -> Tuple[Any, Tuple[Any, ...]]:
    """Build a figure with one axes per panel, sized from the pocket's aspect.

    Panel boxes are placed in explicit figure fractions rather than by a layout
    engine, so that every panel of a comparison gets exactly the same box and the
    drawings stacked above one another are at the same scale.

    Args:
        panels: How many panels to stack.
        bounds: The shared view, whose aspect sets the panel height.
        palette: The theme's tokens.
        has_panel_header: Reserve a header strip above every panel.
        has_suptitle: Reserve a strip at the top of the figure.
        legend_rows: How many rows the legend will wrap to, so the strip reserved
            for it is the strip it needs.

    Returns:
        The figure and its panels, top to bottom.
    """
    xmin, xmax, ymin, ymax = bounds
    usable_width = FIGURE_WIDTH_IN - 2.0 * SIDE_MARGIN_IN
    aspect = (ymax - ymin) / (xmax - xmin) if xmax > xmin else 1.0
    panel_height = min(MAX_PANEL_HEIGHT_IN, max(MIN_PANEL_HEIGHT_IN, usable_width * aspect))
    header = PANEL_HEADER_IN if has_panel_header else 0.0
    suptitle = SUPTITLE_HEIGHT_IN if has_suptitle else 0.0
    legend = legend_rows * LEGEND_ROW_IN + LEGEND_PAD_IN

    figure_height = TOP_MARGIN_IN + suptitle + panels * (panel_height + header) + (panels - 1) * PANEL_GAP_IN + legend + BOTTOM_MARGIN_IN
    figure = Figure(figsize=(FIGURE_WIDTH_IN, figure_height), dpi=FIGURE_DPI, facecolor=palette.surface)
    FigureCanvasAgg(figure)

    axes: List[Any] = []
    cursor = figure_height - TOP_MARGIN_IN - suptitle
    for _ in range(panels):
        cursor -= header + panel_height
        axes.append(
            figure.add_axes(
                (
                    SIDE_MARGIN_IN / FIGURE_WIDTH_IN,
                    cursor / figure_height,
                    usable_width / FIGURE_WIDTH_IN,
                    panel_height / figure_height,
                )
            )
        )
        cursor -= PANEL_GAP_IN
    return figure, tuple(axes)


def _ensure_canvas(figure: Any) -> None:
    """Give a figure an Agg canvas if it has none, so that it draws with no display.

    Args:
        figure: The figure.
    """
    if getattr(figure, "canvas", None) is None:
        FigureCanvasAgg(figure)


def _draw_panel_heading(axis: Any, label: Optional[str], subtitle: Optional[str], palette: Palette) -> None:
    """Put a panel's label, and any subtitle, above its box.

    Args:
        axis: The panel.
        label: The panel's name, in primary ink, or None.
        subtitle: One line of secondary text, or None.
        palette: The theme's tokens.
    """
    if label is not None:
        axis.text(0.0, 1.0, label, transform=axis.transAxes, ha="left", va="bottom", color=palette.ink, fontsize=PANEL_LABEL_PT)
    if subtitle is not None:
        axis.text(1.0, 1.0, subtitle, transform=axis.transAxes, ha="right", va="bottom", color=palette.secondary, fontsize=PANEL_SUBTITLE_PT)


def _draw_headings(figure: Any, title: Optional[str], subtitle: Optional[str], palette: Palette) -> None:
    """Put the figure's title and subtitle at the top left.

    Args:
        figure: The figure.
        title: The title, in primary ink.
        subtitle: One line under it, in secondary ink.
        palette: The theme's tokens.
    """
    height = figure.get_figheight()
    left = SIDE_MARGIN_IN / FIGURE_WIDTH_IN
    if title is not None:
        figure.text(left, 1.0 - TOP_MARGIN_IN / height, title, ha="left", va="top", color=palette.ink, fontsize=TITLE_PT)
    if subtitle is not None:
        below = TOP_MARGIN_IN + (TITLE_PT * TEXT_LEADING) / POINTS_PER_INCH
        figure.text(left, 1.0 - below / height, subtitle, ha="left", va="top", color=palette.secondary, fontsize=SUBTITLE_PT)


def _legend_handles(entries: Mapping[str, Mark]) -> Tuple[List[Any], Tuple[str, ...]]:
    """Legend handles and labels, in legend order.

    Args:
        entries: Label to mark, as drawn.

    Returns:
        The handles and their labels.
    """
    ordered = sorted(entries.values(), key=lambda mark: (mark.order, mark.label))
    handles = [
        Line2D(
            [],
            [],
            color=mark.colour,
            linewidth=mark.width_pt,
            linestyle=mark.linestyle,
            marker=mark.marker or "none",
            markersize=POINT_MARKER_PT,
            markeredgewidth=0.0,
        )
        for mark in ordered
    ]
    return handles, tuple(mark.label for mark in ordered)


def _style_legend(legend: Any, palette: Palette) -> None:
    """Put legend text in ink, never in a series colour.

    Args:
        legend: The legend.
        palette: The theme's tokens.
    """
    for text in legend.get_texts():
        text.set_color(palette.ink)


def _draw_figure_legend(figure: Any, entries: Mapping[str, Mark], palette: Palette) -> Tuple[str, ...]:
    """Draw one legend for the whole figure, in the strip reserved for it.

    Args:
        figure: The figure.
        entries: Label to mark, as drawn.
        palette: The theme's tokens.

    Returns:
        The labels, in legend order.
    """
    handles, labels = _legend_handles(entries)
    legend = figure.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=_legend_layout(labels)[0],
        frameon=False,
        fontsize=LEGEND_PT,
        handlelength=LEGEND_HANDLE_LEN_EM,
        columnspacing=LEGEND_COLUMN_SPACING_EM,
    )
    _style_legend(legend, palette)
    return labels


def _draw_axes_legend(axis: Any, entries: Mapping[str, Mark], palette: Palette) -> Tuple[str, ...]:
    """Draw the legend inside a caller-supplied panel.

    Args:
        axis: The panel.
        entries: Label to mark, as drawn.
        palette: The theme's tokens.

    Returns:
        The labels, in legend order.
    """
    handles, labels = _legend_handles(entries)
    legend = axis.legend(handles, labels, loc="best", frameon=False, fontsize=LEGEND_PT, handlelength=LEGEND_HANDLE_LEN_EM)
    _style_legend(legend, palette)
    return labels


def _legend_layout(labels: Sequence[str]) -> Tuple[int, int]:
    """How the legend wraps: columns across the page, and the rows that implies.

    Estimated from the widest label rather than measured, because the figure must
    be sized before there is a canvas to measure text on. The estimate is used
    for both the reserved strip and the drawn legend, so the two always agree.

    Args:
        labels: The legend labels.

    Returns:
        ``(columns, rows)``, each at least one.
    """
    if not labels:
        return (1, 1)
    widest = max(len(label) for label in labels)
    budget = FIGURE_WIDTH_IN * POINTS_PER_INCH / (LEGEND_PT * LEGEND_GLYPH_EM * (widest + LEGEND_HANDLE_CHARS))
    columns = max(1, min(len(labels), int(budget)))
    return (columns, int(math.ceil(len(labels) / columns)))
