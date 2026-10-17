r"""Plot the gate's toolpaths, one figure per pocket, generators side by side in one view.

Every Stage 1 commit is judged by these figures: the ``pipeline`` path beside
the frozen generators on the three gate pockets, each panel captioned with that
case's gate cells. The path drawn is exactly the gate's -- the generator is
looked up in `benchmarks.gate.GATE_GENERATORS` and fed `benchmarks.gate.gate_pocket`
at the gate's own cap -- and the caption is read from the pinned
`benchmarks/results/gate_vectors.json`, never typed. ``--check-fresh`` re-measures
each plotted case with `tools.gate_measure.measure_case` and refuses to draw a
caption the current code no longer reproduces.

One row per pocket, one panel per generator, side by side: never overlaid,
never stacked. Each panel is drawn by `benchmarks.plotting.draw_toolpath` into
its own axes (operation-class encoding: cut / lead / link strokes, plunge and
retract as point markers, arcs resampled from their own parametrisation); the
row then gets one equal-aspect view shared by every panel, so they compare
directly, and one figure legend. A panel whose path enters by a HELIX also
gets an entry inset in its own column right of the row (`draw_entry_inset`):
a zoomed plan of the helix footprint, its bore disk, the pencil members, the
root loop and the entry P, over a side elevation of Z against plan arc length
showing the ramp turns and the flat closing turn. ``--overview`` also stacks each cap's three
pocket rows into one sheet (not published).

Run (rebuilds the editable extension and re-measures every plotted case
first)::

    pixi run gate-figures

Ad hoc: ``python -m benchmarks.gate_figures [generators ...] [--caps 120 40]
[--out DIR] [--check-fresh] [--overview]``. Generators default to
`benchmarks.gate.GATE_GENERATOR_NAMES` (gate order), caps to every
`benchmarks.gate.GATE_CAP_CASES` value.

Before/after page: add ``--before <ref>`` where ``<ref>`` is a gate-vectors JSON
file, a directory holding ``gate_vectors.json``, or a git ref whose
``benchmarks/results/gate_vectors.json`` is read with ``git show``.

Writes ``<out>/<cap>/<pocket>.png`` (default out:
``docs/assets/images/pipeline``) through `benchmarks.plotting.save_figure`, so
redrawing the same paths rewrites the same bytes.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
from typing import Sequence

from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle as DiskPatch
from PIL import Image

from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_GENERATOR_NAMES
from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket
from benchmarks.marks import CUT_WIDTH_PT
from benchmarks.marks import LINK_WIDTH_PT
from benchmarks.marks import POINT_MARKER_PT
from benchmarks.marks import Mark
from benchmarks.marks import encode
from benchmarks.marks import legend_entries
from benchmarks.palette import Palette
from benchmarks.palette import Theme
from benchmarks.palette import palette_for
from benchmarks.pathgeometry import MARGIN_FRACTION
from benchmarks.pathgeometry import Motion
from benchmarks.pathgeometry import closed_rings
from benchmarks.pathgeometry import drawable_motions
from benchmarks.pathgeometry import helix_elevation
from benchmarks.pathgeometry import view_bounds
from benchmarks.plotting import FIGURE_DPI
from benchmarks.plotting import ColourBy
from benchmarks.plotting import draw_toolpath
from benchmarks.plotting import save_figure
from compas_cgal.helix import HELIX_CLOSING_TURNS
from compas_cgal.helix import Helix
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.toolpath import ToolpathResult
from tools.gate_measure import case_key
from tools.gate_measure import measure_case
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import GateVectors
from tools.gate_vectors import load_gate_vectors
from tools.gate_vectors import parse_gate_vectors

# Repository root: this file is benchmarks/gate_figures.py.
REPO_ROOT = Path(__file__).resolve().parents[1]
# Published location, relative to the working directory (the repository root
# under `pixi run gate-figures`), as `benchmarks.cli figures --out` does.
DEFAULT_OUT = Path("docs/assets/images/pipeline")
HEAD_VECTORS = REPO_ROOT / GATE_VECTORS_PATH

# The file a directory given to --before must hold.
BEFORE_DIR_FILE = "gate_vectors.json"

OVERVIEW_NAME = "overview.png"

# The gate cells every panel is captioned with, as (criterion name, caption
# label, unit suffix), in caption order. Names are checked against the loaded
# document's criteria, so a renamed criterion fails here rather than captioning
# a blank.
CAPTION_CELLS: tuple[tuple[str, str, str], ...] = (
    ("uncut_fraction", "uncut", ""),
    ("material_plunges", "plunges", ""),
    ("redundant_operations", "redundant", ""),
    ("tangent_breaks", "tangent breaks", ""),
    ("max_engagement_step_deg", "max engagement step", "°"),
)

# Caption cells per text line. A panel is `PANEL_WIDTH_IN` wide; two
# before->after cells run to ~2.3 in at `CAPTION_PT`, three to ~3.5 in, so two
# per line is the widest wrap that stays inside its own panel.
CELLS_PER_LINE = 2

# Row geometry, in inches. Panel width fixes the scale; the panel height follows
# from the shared view's aspect.
PANEL_WIDTH_IN = 3.6
PANEL_GAP_IN = 0.25
SIDE_MARGIN_IN = 0.15
# Header above the panels: figure title + provenance line.
FIGURE_HEADER_IN = 0.55
# Per-panel caption strip: generator name + up to four caption lines
# (new-case note + five cells at two per line) at `CAPTION_PT` with 1.3 leading.
PANEL_CAPTION_IN = 0.75
# Strip under the panels for the single figure legend (one row).
LEGEND_STRIP_IN = 0.35

TITLE_PT = 10.5
PROVENANCE_PT = 8.0
CAPTION_PT = 7.5
LEGEND_PT = 8.0

# Entry inset (brief §6): one column right of the row per panel whose path
# enters by a HELIX. The plan zoom sits over a short side elevation; both are
# labelled (`Axes.get_label`) so a test can find them. Inches.
ENTRY_PLAN_LABEL = "entry-plan:{name}"
ENTRY_ELEVATION_LABEL = "entry-elevation:{name}"
INSET_WIDTH_IN = 2.6
# Elevation block: the axes plus room under it for tick labels and the axis label.
ELEVATION_AXES_IN = 0.55
ELEVATION_LABELS_IN = 0.4
ELEVATION_GAP_IN = 0.15
# The plan zoom never shrinks below this, so a flat pocket's row grows instead.
INSET_MIN_PLAN_IN = 1.4
# Left room inside the inset column for the elevation's Z tick labels.
ELEVATION_LEFT_IN = 0.4
TICK_PT = 6.5
# The bore's legend swatch is a thick stroke in the bore's fill colour.
BORE_SWATCH_PT = 6.0
# Legend order of the inset entries: after every motion mark of `benchmarks.marks`.
_ORDER_BORE = 400
_ORDER_MEMBER = 401
_ORDER_ROOT = 402
_ORDER_ENTRY_POINT = 403

CELL_SEPARATOR = "  ·  "
CHANGE_ARROW = " → "


class UnregisteredGeneratorError(Exception):
    """A generator named on the command line is not in `GATE_GENERATORS`."""


class UnknownCapError(Exception):
    """The cap asked for is not one of `GATE_CAP_CASES`."""


class UnknownCaptionCriterionError(Exception):
    """A caption cell names a criterion the gate-vectors document does not carry."""


class MissingGateCaseError(Exception):
    """A plotted case has no cells in the gate-vectors document."""


class ExcludedGateCaseError(Exception):
    """A plotted case is excluded by an admission: its generator raises, so there is no path."""


class BeforeSourceError(Exception):
    """The --before argument is neither a gate-vectors file, a directory holding one, nor a git ref carrying one."""


class StaleGateCellsError(Exception):
    """A fresh measurement of a plotted case disagrees with the pinned cells."""


@dataclass(frozen=True)
class VectorsSource:
    """A gate-vectors document and where it came from.

    Attributes:
        vectors: The validated document.
        label: Provenance for the figure subtitle (path or short commit).
    """

    vectors: GateVectors
    label: str


class UnrecordedHelixEntryError(Exception):
    """A path carries a HELIX op that no generator `RampOutRecord` accounts for.

    The entry inset reads the helix, the pencil members and the root loop from
    the generator's record (brief N2); it never re-derives them from the
    stream, so a helix without a record has no drawable entry.
    """


@dataclass(frozen=True)
class EntryGroup:
    """The helical entry of one path: the HELIX op, the pencil members and the root loop.

    Attributes:
        helix: The helix.
        helix_index: Its position in the operation stream.
        member_indices: The pencil members, in stream order; empty when the
            first scan admitted the root.
        root_index: The root loop.
    """

    helix: Helix
    helix_index: int
    member_indices: tuple[int, ...]
    root_index: int


class MissingPanelLegendError(RuntimeError):
    """A panel came back without the axes legend `draw_toolpath` draws when given an axes."""


def resolve_generators(names: Sequence[str]) -> tuple[str, ...]:
    """Check every name against the gate's registry, keeping the given order.

    Args:
        names: Generator names from the command line.

    Returns:
        The names, unchanged.

    Raises:
        UnregisteredGeneratorError: A name is not in `GATE_GENERATORS`.
    """
    unknown = [name for name in names if name not in GATE_GENERATORS]
    if unknown:
        raise UnregisteredGeneratorError(f"{unknown} are not registered in benchmarks.gate.GATE_GENERATORS; registered: {sorted(GATE_GENERATORS)}.")
    return tuple(names)


def resolve_cap(label: str) -> tuple[str, GateCapDegrees]:
    """The gate cap case whose degrees print as *label*.

    Args:
        label: Cap in degrees as `f"{deg:g}"`, e.g. ``"120"``.

    Returns:
        The cap id and degrees.

    Raises:
        UnknownCapError: No gate cap case matches.
    """
    for cap_id, cap in GATE_CAP_CASES:
        if f"{cap:g}" == label:
            return cap_id, cap
    raise UnknownCapError(f"cap {label!r} is not a gate cap; gate caps: {[f'{cap:g}' for _, cap in GATE_CAP_CASES]}.")


def load_before(argument: str) -> VectorsSource:
    """Load the document a before/after page compares against.

    Args:
        argument: A gate-vectors file, a directory holding `BEFORE_DIR_FILE`,
            or a git ref of this repository.

    Returns:
        The document and its provenance.

    Raises:
        BeforeSourceError: The argument resolves to none of the three.
    """
    path = Path(argument)
    if path.is_file():
        return VectorsSource(load_gate_vectors(path), str(path))
    if path.is_dir():
        inner = path / BEFORE_DIR_FILE
        if not inner.is_file():
            raise BeforeSourceError(f"--before directory {path} holds no {BEFORE_DIR_FILE}.")
        return VectorsSource(load_gate_vectors(inner), str(inner))
    commit = _git(("rev-parse", "--short", "--verify", f"{argument}^{{commit}}"), argument).strip()
    text = _git(("show", f"{commit}:{GATE_VECTORS_PATH.as_posix()}"), argument)
    return VectorsSource(parse_gate_vectors(text), f"{argument} ({commit})")


def _git(args: Sequence[str], argument: str) -> str:
    """Run a read-only git query in the worktree.

    Args:
        args: The git arguments.
        argument: The --before value, for the error message.

    Returns:
        Standard output.

    Raises:
        BeforeSourceError: git exited non-zero.
    """
    done = subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise BeforeSourceError(f"--before {argument!r} is not a file, a directory, or a git ref carrying {GATE_VECTORS_PATH}: git {' '.join(args)}: {done.stderr.strip()}")
    return done.stdout


def check_caption_criteria(vectors: GateVectors, label: str) -> None:
    """Require every caption criterion in *vectors*.

    Args:
        vectors: The document.
        label: Its provenance, for the error message.

    Raises:
        UnknownCaptionCriterionError: A caption criterion is absent.
    """
    known = {criterion.name for criterion in vectors.criteria}
    missing = [name for name, _, _ in CAPTION_CELLS if name not in known]
    if missing:
        raise UnknownCaptionCriterionError(f"{label} carries no criteria {missing}.")


def head_cells(vectors: GateVectors, key: str) -> Mapping[str, str]:
    """The pinned cells of one plotted case.

    Args:
        vectors: The head document.
        key: The case key.

    Returns:
        Criterion name to stored decimal string.

    Raises:
        ExcludedGateCaseError: An admission excludes the case.
        MissingGateCaseError: The document has no such case.
    """
    excluded = vectors.excluded
    if key in excluded:
        raise ExcludedGateCaseError(f"{key} is excluded (its generator raises {excluded[key]}); there is no path to plot. Drop the generator or admit the case.")
    if key not in vectors.cases:
        raise MissingGateCaseError(f"{key} has no cells in {HEAD_VECTORS}; re-measure the gate (pixi run gate-vectors) after registering the generator.")
    return vectors.cases[key]


def check_fresh(case: tuple[str, float, str, str], pinned: Mapping[str, str], vectors: GateVectors) -> None:
    """Re-measure one case and require the pinned freshness-checked cells.

    Args:
        case: ``(cap id, cap degrees, generator, pocket)``.
        pinned: The pinned cells.
        vectors: The head document, for its freshness-checked criteria.

    Raises:
        StaleGateCellsError: A freshness-checked cell differs.
    """
    key, fresh = measure_case(case)
    stale = {spec.name: (pinned[spec.name], fresh[spec.name]) for spec in vectors.freshness_checked if pinned[spec.name] != fresh[spec.name]}
    if stale:
        raise StaleGateCellsError(f"{key}: pinned cells are stale against the current code (pinned, fresh): {stale}.")


def caption(new: Mapping[str, str], old: Mapping[str, str] | None, before: VectorsSource | None) -> str:
    """The panel subtitle: the caption cells, before -> after when a before document is given.

    Args:
        new: The case's head cells.
        old: The case's cells in the before document, or None when it has none.
        before: The before document, or None when no comparison was asked for.

    Returns:
        The caption, wrapped `CELLS_PER_LINE` cells per line.
    """
    parts: list[str] = [] if before is None or old is not None else [f"new case (absent at {before.label})"]
    for name, label, unit in CAPTION_CELLS:
        if old is None:
            parts.append(f"{label} {new[name]}{unit}")
        elif old[name] == new[name]:
            parts.append(f"{label} {new[name]}{unit} (=)")
        else:
            parts.append(f"{label} {old[name]}{CHANGE_ARROW}{new[name]}{unit}")
    lines = [CELL_SEPARATOR.join(parts[i : i + CELLS_PER_LINE]) for i in range(0, len(parts), CELLS_PER_LINE)]
    return "\n".join(lines)


def compose_pocket(
    pocket: str,
    cap_id: str,
    cap: GateCapDegrees,
    generators: Sequence[str],
    head: VectorsSource,
    before: VectorsSource | None,
    fresh: bool,
    out_dir: Path,
) -> Path:
    """Generate, caption and draw one pocket's figure.

    Args:
        pocket: One of `GATE_POCKET_NAMES`.
        cap_id: The gate cap id.
        cap: The cap in degrees.
        generators: Registered generator names, in panel order.
        head: The pinned document.
        before: The comparison document, or None.
        fresh: Re-measure every case against its pinned cells first.
        out_dir: Where the PNG goes.

    Returns:
        The PNG written.
    """
    spec = gate_pocket(pocket, tea_cap_deg=cap)
    panels: dict[str, ToolpathResult] = {}
    subtitles: dict[str, str] = {}
    for generator in generators:
        key = case_key(cap_id, generator, pocket)
        cells = head_cells(head.vectors, key)
        if fresh:
            check_fresh((cap_id, float(cap), generator, pocket), cells, head.vectors)
        started = time.perf_counter()
        panels[generator] = GATE_GENERATORS[generator](spec)
        print(f"{key}: generated in {time.perf_counter() - started:.1f} s", flush=True)
        old = None if before is None else before.vectors.cases.get(key)
        subtitles[generator] = caption(cells, old, before)
    provenance = f"gate cells: {head.label} (rules {head.vectors.rules_version})"
    if before is not None:
        provenance += f"; before: {before.label} (rules {before.vectors.rules_version})"
    figure = draw_row(
        panels,
        subtitles,
        boundary=spec.polygon,
        holes=list(spec.holes),
        tool_diameter=spec.tool_diameter,
        title=f"{pocket}  ·  cap {cap:g}°  ·  tool ⌀{spec.tool_diameter:g}",
        provenance=provenance,
    )
    return save_figure(figure, out_dir / f"{pocket}.png")


def draw_row(
    panels: Mapping[str, ToolpathResult],
    captions: Mapping[str, str],
    *,
    boundary: object,
    holes: Sequence[object],
    tool_diameter: float,
    title: str,
    provenance: str,
) -> Figure:
    """Draw one pocket's generators side by side, one panel each, in one shared view.

    Every panel is drawn by `draw_toolpath` into its own axes; its per-axes
    legend is removed and replaced by one figure legend built from the same
    marks. The view is the union over every panel, so equal aspect plus equal
    limits puts every panel at one scale.

    Args:
        panels: Generator name to its path, in panel order.
        captions: Generator name to its caption.
        boundary: The pocket boundary.
        holes: The pocket islands.
        tool_diameter: Cutter diameter, passed through for view padding.
        title: Figure title.
        provenance: Where the captions' cells came from.

    Returns:
        The figure, on its own Agg canvas.
    """
    palette = palette_for(Theme.LIGHT)
    rings = closed_rings(boundary, holes)
    motions = {name: drawable_motions(result) for name, result in panels.items()}
    entries_of: dict[str, EntryGroup] = {}
    for name, result in panels.items():
        group = entry_group(result)
        if group is not None:
            entries_of[name] = group
    xmin, xmax, ymin, ymax = view_bounds(rings, [motion for row in motions.values() for motion in row], None)
    panel_height = PANEL_WIDTH_IN * (ymax - ymin) / (xmax - xmin)
    elevation_block = ELEVATION_GAP_IN + ELEVATION_AXES_IN + ELEVATION_LABELS_IN
    body = max(panel_height, INSET_MIN_PLAN_IN + elevation_block) if entries_of else panel_height
    count = len(panels)
    width = 2.0 * SIDE_MARGIN_IN + count * PANEL_WIDTH_IN + (count - 1) * PANEL_GAP_IN + len(entries_of) * (PANEL_GAP_IN + INSET_WIDTH_IN)
    height = FIGURE_HEADER_IN + PANEL_CAPTION_IN + body + LEGEND_STRIP_IN
    figure = Figure(figsize=(width, height), dpi=FIGURE_DPI, facecolor=palette.surface)
    FigureCanvasAgg(figure)

    entries: dict[str, Mark] = {}
    bottom = LEGEND_STRIP_IN / height
    for slot, (name, result) in enumerate(panels.items()):
        left = (SIDE_MARGIN_IN + slot * (PANEL_WIDTH_IN + PANEL_GAP_IN)) / width
        axis = figure.add_axes((left, bottom, PANEL_WIDTH_IN / width, panel_height / height))
        draw_toolpath(result, boundary=boundary, holes=holes, colour_by=ColourBy.OPERATION, tool_diameter=tool_diameter, ax=axis, theme=Theme.LIGHT)
        panel_legend = axis.get_legend()
        if panel_legend is None:
            raise MissingPanelLegendError(f"draw_toolpath drew no axes legend for {name}; the row legend is built from the panels' entries.")
        panel_legend.remove()
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)
        axis.set_aspect("equal", adjustable="datalim")
        axis.apply_aspect()
        axis.set_title(f"{name}\n{captions[name]}", loc="left", fontsize=CAPTION_PT, color=palette.ink, linespacing=1.3)
        for label, mark in legend_entries(encode(motions[name], ColourBy.OPERATION, palette, None)).items():
            entries.setdefault(label, mark)

    top = (LEGEND_STRIP_IN + body) / height
    for slot, (name, group) in enumerate(entries_of.items()):
        left_in = SIDE_MARGIN_IN + count * (PANEL_WIDTH_IN + PANEL_GAP_IN) + slot * (INSET_WIDTH_IN + PANEL_GAP_IN)
        plan_height = body - elevation_block
        plan = figure.add_axes((left_in / width, top - plan_height / height, INSET_WIDTH_IN / width, plan_height / height))
        elevation = figure.add_axes(
            (
                (left_in + ELEVATION_LEFT_IN) / width,
                bottom + ELEVATION_LABELS_IN / height,
                (INSET_WIDTH_IN - ELEVATION_LEFT_IN) / width,
                ELEVATION_AXES_IN / height,
            )
        )
        draw_entry_inset(plan, elevation, name=name, group=group, motions=motions[name], tool_radius=0.5 * tool_diameter, palette=palette)
        for label, mark in entry_inset_marks(palette, group).items():
            entries.setdefault(label, mark)

    figure.text(SIDE_MARGIN_IN / width, 1.0 - 0.08 / height, title, ha="left", va="top", color=palette.ink, fontsize=TITLE_PT)
    figure.text(SIDE_MARGIN_IN / width, 1.0 - 0.32 / height, provenance, ha="left", va="top", color=palette.secondary, fontsize=PROVENANCE_PT)
    ordered = sorted(entries.values(), key=lambda mark: (mark.order, mark.label))
    handles = [
        Line2D([], [], color=mark.colour, linewidth=mark.width_pt, linestyle=mark.linestyle, marker=mark.marker or "none", markersize=POINT_MARKER_PT, markeredgewidth=0.0)
        for mark in ordered
    ]
    legend = figure.legend(handles, [mark.label for mark in ordered], loc="lower center", ncol=len(ordered), frameon=False, fontsize=LEGEND_PT)
    for text in legend.get_texts():
        text.set_color(palette.ink)
    return figure


def entry_group(result: ToolpathResult) -> EntryGroup | None:
    """The helical entry of *result*, as its generator recorded it.

    The helix, the pencil members and the root loop are read from the
    pipeline's `RampOutRecord` (brief N2), which the generator validated
    against the stream when it emitted the entry. Nothing is inferred from
    coincident entries: a walk loop entered at P after the root, or a root the
    walk skipped, would make any such inference misclassify silently.

    Args:
        result: The path; a `PipelineResult` carries the record.

    Returns:
        The entry, or None when the path carries no HELIX op.

    Raises:
        UnrecordedHelixEntryError: A HELIX op is not the one a record names
            (the path carries no record, or a second helix).
    """
    operations = list(result.operations)
    helices = [index for index, operation in enumerate(operations) if isinstance(operation.geometry, Helix)]
    record = result.ramp_out if isinstance(result, PipelineResult) else None
    if record is None:
        if helices:
            raise UnrecordedHelixEntryError(f"HELIX ops at {helices} on a {type(result).__name__}, which carries no RampOutRecord to read the entry from.")
        return None
    if helices != [record.helix_operation]:
        raise UnrecordedHelixEntryError(f"HELIX ops at {helices}; the ramp-out record names only {record.helix_operation}.")
    helix = operations[record.helix_operation].geometry
    assert isinstance(helix, Helix)
    return EntryGroup(helix=helix, helix_index=record.helix_operation, member_indices=tuple(record.members), root_index=record.root_operation)


def entry_inset_bounds(group: EntryGroup, motions: Sequence[Motion], tool_radius: float) -> tuple[float, float, float, float]:
    """The entry inset's view: the bore disk and every entry loop, padded.

    The padding is `MARGIN_FRACTION` of the longer side plus one tool radius,
    so the root loop's swept width would fit inside the view too.

    Args:
        group: The entry.
        motions: The path's motions, indexed as its operations.
        tool_radius: Cutter radius r; the bore is ``disk(c_h, rho_h + r)``.

    Returns:
        ``(xmin, xmax, ymin, ymax)``.
    """
    cx, cy = float(group.helix.footprint.centre.x), float(group.helix.footprint.centre.y)
    bore = float(group.helix.footprint.radius.value) + tool_radius
    points = [point for index in (*group.member_indices, group.root_index) for point in motions[index].points]
    xs = [cx - bore, cx + bore, *(x for x, _ in points)]
    ys = [cy - bore, cy + bore, *(y for _, y in points)]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    margin = MARGIN_FRACTION * max(xmax - xmin, ymax - ymin) + tool_radius
    return (xmin - margin, xmax + margin, ymin - margin, ymax + margin)


def entry_inset_marks(palette: Palette, group: EntryGroup) -> dict[str, Mark]:
    """The marks the entry inset draws, by label: exactly those it draws, so the legend names nothing absent.

    Args:
        palette: The theme's tokens.
        group: The entry; the pencil-member mark is present only when it has members.

    Returns:
        Bore, pencil member (when any), root loop and entry point marks.
    """
    marks = [Mark(label="bore (footprint + r)", order=_ORDER_BORE, colour=palette.grid, width_pt=BORE_SWATCH_PT, linestyle="-", marker="")]
    if group.member_indices:
        marks.append(Mark(label="pencil member", order=_ORDER_MEMBER, colour=palette.cut, width_pt=LINK_WIDTH_PT, linestyle="-", marker=""))
    marks.append(Mark(label="root loop", order=_ORDER_ROOT, colour=palette.cut, width_pt=CUT_WIDTH_PT, linestyle="-", marker=""))
    marks.append(Mark(label="P (helix entry)", order=_ORDER_ENTRY_POINT, colour=palette.ink, width_pt=0.0, linestyle="none", marker="o"))
    return {mark.label: mark for mark in marks}


def draw_entry_inset(
    plan: Axes,
    elevation: Axes,
    *,
    name: str,
    group: EntryGroup,
    motions: Sequence[Motion],
    tool_radius: float,
    palette: Palette,
) -> None:
    """Draw one path's entry: a plan zoom and a Z-against-arc-length elevation.

    The plan draws the bore disk, the helix as its plan projection (its own
    mark, `benchmarks.marks`), the pencil members, the root loop and P, inside
    a full zoom-window frame; the elevation draws the helix op's ramp and its flat
    closing turn (`benchmarks.pathgeometry.helix_elevation`). REPORTING only.

    Args:
        plan: The plan axes.
        elevation: The elevation axes.
        name: The panel's generator, for titles and labels.
        group: The entry.
        motions: The path's motions, indexed as its operations.
        tool_radius: Cutter radius r.
        palette: The theme's tokens.
    """
    marks = entry_inset_marks(palette, group)
    helix = group.helix
    xmin, xmax, ymin, ymax = entry_inset_bounds(group, motions, tool_radius)
    plan.set_label(ENTRY_PLAN_LABEL.format(name=name))
    plan.set_facecolor(palette.surface)
    plan.set_xlim(xmin, xmax)
    plan.set_ylim(ymin, ymax)
    # "box", not the row's "datalim": the zoom's limits are the contract (they
    # hold the bore and the root loop), so the axes box shrinks to fit instead.
    plan.set_aspect("equal", adjustable="box", anchor="N")
    plan.apply_aspect()
    # A zoom window, framed on all four sides: the pocket rings are not drawn,
    # because the zoom clips them to partial walls that read as a cut frame.
    plan.set_xticks([])
    plan.set_yticks([])
    for spine in plan.spines.values():
        spine.set_color(palette.secondary)
        spine.set_linewidth(LINK_WIDTH_PT)
    centre = (float(helix.footprint.centre.x), float(helix.footprint.centre.y))
    plan.add_patch(DiskPatch(centre, float(helix.footprint.radius.value) + tool_radius, facecolor=palette.grid, edgecolor="none", zorder=1))
    helix_mark = encode([motions[group.helix_index]], ColourBy.OPERATION, palette, None)[group.helix_index]
    strokes = [(group.helix_index, helix_mark)]
    strokes += [(index, marks["pencil member"]) for index in group.member_indices]
    strokes.append((group.root_index, marks["root loop"]))
    for index, mark in strokes:
        points = motions[index].points
        plan.plot([x for x, _ in points], [y for _, y in points], color=mark.colour, linewidth=mark.width_pt, linestyle=mark.linestyle, solid_joinstyle="round", zorder=3)
    point = marks["P (helix entry)"]
    plan.plot([helix.entry_xy.x], [helix.entry_xy.y], color=point.colour, linestyle="none", marker=point.marker, markersize=POINT_MARKER_PT, markeredgewidth=0.0, zorder=4)
    plan.set_title(
        f"{name}  ·  entry (zoom)\nhelix {helix.turns} + {HELIX_CLOSING_TURNS} flat turns  ·  {len(group.member_indices)} members",
        loc="left",
        fontsize=CAPTION_PT,
        color=palette.ink,
        linespacing=1.3,
    )

    profile = helix_elevation(helix)
    elevation.set_label(ENTRY_ELEVATION_LABEL.format(name=name))
    elevation.set_facecolor(palette.surface)
    ramp, flat = profile[:2], profile[1:]
    # The ramp in the helix mark; the flat closing turn solid in the same ink, so
    # the one op's two parts read apart even when the flat turn is a short tail.
    elevation.plot([s for s, _ in ramp], [z for _, z in ramp], color=helix_mark.colour, linewidth=helix_mark.width_pt, linestyle=helix_mark.linestyle, zorder=3)
    elevation.plot([s for s, _ in flat], [z for _, z in flat], color=helix_mark.colour, linewidth=helix_mark.width_pt, linestyle="-", solid_capstyle="butt", zorder=3)
    ramp_end, closed = profile[1][0], profile[-1][0]
    top, bottom = profile[0][1], profile[-1][1]
    elevation.text(0.5 * ramp_end, 0.5 * (top + bottom), f"  ramp ×{helix.turns}", ha="left", va="bottom", fontsize=TICK_PT, color=palette.secondary)
    elevation.text(closed, top, f"flat closing turn ×{HELIX_CLOSING_TURNS} at bottom Z", ha="right", va="top", fontsize=TICK_PT, color=palette.secondary)
    elevation.set_xlim(0.0, closed)
    elevation.set_ylim(bottom - 0.1 * (top - bottom), top + 0.1 * (top - bottom))
    elevation.set_xlabel("plan arc length (mm)", fontsize=TICK_PT, color=palette.secondary, labelpad=1.0)
    elevation.set_ylabel("Z (mm)", fontsize=TICK_PT, color=palette.secondary, labelpad=1.0)
    elevation.tick_params(labelsize=TICK_PT, colors=palette.secondary, length=2.0, pad=1.0)
    for side in ("top", "right"):
        elevation.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        elevation.spines[side].set_color(palette.secondary)


def overview(paths: Sequence[Path], out_dir: Path) -> Path:
    """Stack the pocket figures into one sheet.

    Args:
        paths: The pocket PNGs, top to bottom.
        out_dir: Where the sheet goes.

    Returns:
        The sheet written.
    """
    images = [Image.open(path).convert("RGB") for path in paths]
    sheet_width = images[0].width
    thumbs = [image.resize((sheet_width, round(image.height * sheet_width / image.width))) for image in images]
    sheet = Image.new("RGB", (sheet_width, sum(thumb.height for thumb in thumbs)), "white")
    y = 0
    for thumb in thumbs:
        sheet.paste(thumb, (0, y))
        y += thumb.height
    out = out_dir / OVERVIEW_NAME
    sheet.save(out)
    return out


def main(argv: Sequence[str] | None = None) -> int:
    """Draw every gate pocket for the named generators at every named cap.

    Args:
        argv: Command-line arguments; None reads `sys.argv`.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(prog="python -m benchmarks.gate_figures", description="Plot gate toolpaths per pocket, captioned with their gate cells.")
    parser.add_argument(
        "generators",
        nargs="*",
        default=list(GATE_GENERATOR_NAMES),
        help="Generator names registered in benchmarks.gate.GATE_GENERATORS, in panel order (default: GATE_GENERATOR_NAMES).",
    )
    parser.add_argument("--caps", nargs="+", default=[f"{cap:g}" for _, cap in GATE_CAP_CASES], help="Gate caps in degrees (GATE_CAP_CASES values; default: all).")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help=f"Output root; figures go to <out>/<cap>/<pocket>.png (default: {DEFAULT_OUT}).")
    parser.add_argument("--before", default=None, help="Gate-vectors file, directory holding gate_vectors.json, or git ref, for before -> after captions.")
    parser.add_argument("--check-fresh", action="store_true", help="Re-measure each case and refuse captions the current code no longer reproduces.")
    parser.add_argument("--overview", action="store_true", help=f"Also stack each cap's pocket figures into <out>/<cap>/{OVERVIEW_NAME}.")
    args = parser.parse_args(argv)

    generators = resolve_generators(args.generators)
    caps = [resolve_cap(label) for label in args.caps]
    head = VectorsSource(load_gate_vectors(HEAD_VECTORS), GATE_VECTORS_PATH.as_posix())
    check_caption_criteria(head.vectors, head.label)
    before = None if args.before is None else load_before(args.before)
    if before is not None:
        check_caption_criteria(before.vectors, before.label)

    for cap_id, cap in caps:
        out_dir = args.out / f"{cap:g}"
        out_dir.mkdir(parents=True, exist_ok=True)
        written = [compose_pocket(pocket, cap_id, cap, generators, head, before, args.check_fresh, out_dir) for pocket in GATE_POCKET_NAMES]
        if args.overview:
            written.append(overview(written, out_dir))
        for path in written:
            print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
