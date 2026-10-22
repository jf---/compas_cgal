"""Draw a Held toolpath as the material each motion sweeps, painted in emission order.

The tool (radius r) following a centre circle of radius R sweeps the annulus
``max(R - r, 0) <= |p - c| <= R + r``; following a connector polyline it sweeps
the capsule ``polyline (+) D_r``. Each swept region is painted opaque on a
sequential ramp, in emission order (circle i, then the connector i -> i+1),
later over earlier, so the crescent left of a region is the material only that
motion reached. Those are exactly the regions the coverage certificate
(``compas_cgal.coverage_certificate``) decides on; its verdict and exact
witnesses are drawn over the picture. This is a display: every double here is
a reported value, no geometric decision is made.

Display invariant: band paint lies inside ``path (+) D_r`` and fills it, so
page colour inside the pocket means uncut, up to the raster (one pixel of
antialiasing) and the turn polygon's chord sagitta ``R (1 - cos(pi / N))``,
N = ``TURN_DISPLAY_SAMPLES``. The largest emitted R on the published pockets is
7.08 mm (repaired Monstera), sagitta 3.8e-3 mm against a pixel of about
0.06 mm at the figure's scale. To keep it:

- each band's separating rim lies inside the band: the rim stroke is exactly
  2r wide and the colour stroke 2r minus both rims on top of it (the rim is
  dropped where it would leave the colour narrower than one rim), and the rim
  is a darker shade of the band's own colour, never near the page;
- the pocket outline is clipped to the outside of the pocket, so it covers no
  interior pixel.

Markers drawn on top (start, seam centre line, over-cap crosses, witness
rings) are annotations outside this invariant.

The figure is drawn from a :class:`SequencePlotData` of plain doubles; the
generator runs and the cached JSON round trip live in ``held_sequence_cli``,
the certificate run in ``held_sequence_coverage``.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Colormap
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.image import imread
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_motion_coverage import connector_pairs
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import maximum_predecessor_engagement

DATA_VERSION = 3
# Display resolution of one full turn; no geometric decisions. Chord sagitta
# R (1 - cos(pi / 96)) = 5.4e-4 R: 3.8e-3 mm at R = 7.08 mm, 4.9e-3 mm at the
# Figure 5 maximum R = 9.02 mm, both under 0.13 px at the page scale.
TURN_DISPLAY_SAMPLES = 96
# Ramp window of matplotlib's perceptually uniform 'magma': the lowest 4 % is
# near-black (loses the stroke against the grid) and the top 12 % is pale cream
# that vanishes on the light page; [0.04, 0.88] keeps lightness monotone and
# every stroke visible.
RAMP_LOW = 0.04
RAMP_HIGH = 0.88
RAMP_STEPS = 256
PAGE_COLOR = "#faf9f5"
BOUNDARY_COLOR = "#2f8f3a"
SEAM_COLOR = "#0b6fa4"
CAP_COLOR = "#c34c32"
COVERED_COLOR = "#2f8f3a"
WITNESS_COLOR = "#1f6fd1"
# Each swept region wears a thin rim, inside its own band, so the crescent it
# leaves behind is separated from the next region's crescent. The rim is the
# band's colour scaled toward black by RIM_SHADE (a darker edge of the same hue:
# a crescent thinner than the rim still reads as paint, never as page). Width in
# points per side, small next to a band's 2r.
DISK_EDGE_PT = 0.3
RIM_SHADE = 0.72
# The rim is dropped for a member whose colour stroke (2r minus both rims)
# would be narrower than this many rim widths: the band never vanishes.
MIN_COLOUR_RIMS = 1.0
# The pocket outline is drawn centred on the wall and clipped to the outside of
# the pocket: half of this width (points) shows, outside, and nothing inside.
BOUNDARY_LINE_PT = 2.4
FIGURE_WIDTH_IN = 13.0
# Page layout in inches. The pocket panel's height follows the pocket's aspect
# (an elongated pocket does not leave a band of empty page), within bounds.
PATH_PANEL_WIDTH_IN = 10.5
PATH_PANEL_MIN_IN = 2.5
PATH_PANEL_MAX_IN = 8.5
ENGAGEMENT_PANEL_IN = 2.3
PANEL_GAP_IN = 1.0
HEADER_IN = 1.05
FOOTER_IN = 0.85
NOTE_LINE_IN = 0.2  # One 9 pt footer line plus leading.
DPI = 180
# Agg's clip of unfilled paths: the canvas grown by this many pixels on each side.
CANVAS_CLIP_PAD_PX = 1.0
WITNESS_LIST_LIMIT = 4  # Uncut sites spelled out in the caption; the rest are counted.
VERDICT_LINE_IN = 0.3  # Header line for the certificate verdict.


class SequencePlotError(BenchmarkError):
    """Base error of the sequence plot."""


class SequenceDataVersionError(SequencePlotError):
    """A cached sequence file does not carry the version this module reads."""


class EmptySequenceError(SequencePlotError):
    """A sequence plot needs at least one circle."""


class InvalidCoverageVerdictError(SequencePlotError):
    """A coverage verdict's witnesses contradict its outcome."""


class StrokeClippedAtCanvasError(SequencePlotError):
    """Agg would cut a band's centre path at the canvas edge close enough to the axes to leave page colour inside a band."""


@dataclass(frozen=True)
class SequenceCircle:
    """One emitted circle as reported doubles (millimetres, phase in radians) and its guide run."""

    x: float
    y: float
    radius: float
    phase: float
    run: int


@dataclass(frozen=True)
class CoverageVerdict:
    """The coverage certificate's outcome on the drawn path, as reported doubles.

    ``witnesses`` are the certificate's exact uncovered points rounded to
    binary64 (millimetres), with the kind of the exposed piece beside each.
    """

    covered: bool
    witnesses: tuple[tuple[float, float], ...]
    piece_kinds: tuple[str, ...]

    @classmethod
    def build(cls, *, covered: bool, witnesses: tuple[tuple[float, float], ...], piece_kinds: tuple[str, ...]) -> CoverageVerdict:
        """Validate that a covered verdict has no witness and an uncovered one has at least one.

        Raises:
            InvalidCoverageVerdictError: the witnesses contradict ``covered``,
                or there is not one piece kind per witness.
        """
        if covered == bool(witnesses):
            raise InvalidCoverageVerdictError(f"covered={covered} with {len(witnesses)} witnesses.")
        if len(piece_kinds) != len(witnesses):
            raise InvalidCoverageVerdictError(f"{len(piece_kinds)} piece kinds for {len(witnesses)} witnesses.")
        return cls(covered, witnesses, piece_kinds)

    def sites(self, join_mm: float) -> tuple[tuple[float, float, int], ...]:
        """Witnesses grouped for display: single linkage at ``join_mm``, each site (mean x, mean y, count), ordered by (y, x).

        One uncut island is bounded by several exposed edges, each with its own
        witness; a site is one island as the eye sees it. Reporting only.
        """
        points = np.array(self.witnesses, dtype=np.float64).reshape(-1, 2)
        group = list(range(len(points)))

        def root(index: int) -> int:
            while group[index] != index:
                index = group[index]
            return index

        for first in range(len(points)):
            for second in range(first + 1, len(points)):
                if float(np.hypot(*(points[first] - points[second]))) <= join_mm:
                    group[root(second)] = root(first)
        members: dict[int, list[int]] = {}
        for index in range(len(points)):
            members.setdefault(root(index), []).append(index)
        sites = [(float(points[m, 0].mean()), float(points[m, 1].mean()), len(m)) for m in members.values()]
        return tuple(sorted(sites, key=lambda site: (site[1], site[0])))

    def caption(self, join_mm: float) -> str:
        """One line: the verdict and the uncut sites, labelled A, B, ... as on the plot."""
        if self.covered:
            return "Coverage certificate: CERTIFIED COVERED (exact decision on these swept regions, closing seam included)."
        sites = self.sites(join_mm)
        listed = ", ".join(f"{witness_label(index)} ({x:.2f}, {y:.2f}) \u00d7{count}" for index, (x, y, count) in enumerate(sites[:WITNESS_LIST_LIMIT]))
        more = f" and {len(sites) - WITNESS_LIST_LIMIT} more" if len(sites) > WITNESS_LIST_LIMIT else ""
        return f"Coverage certificate: NOT COVERED \u00b7 {len(self.witnesses)} exact witnesses of uncut material at {len(sites)} site(s): {listed}{more} (mm)."


@dataclass(frozen=True)
class SequencePlotData:
    """Everything the figure draws, in emission order, as plain doubles.

    ``engagement_deg[i]`` is the engagement value plotted for ``circles[i + 1]``
    after ``circles[i]`` (NaN: unresolved model). ``connectors[i]`` is the
    centre polyline from ``circles[i]`` to ``circles[i + 1]``;
    ``seam_connector`` is the closing motion last -> first, if one is emitted.
    ``tool_radius`` is r: every centre path is drawn with the band it sweeps.
    """

    workload: str
    title: str
    boundary: tuple[tuple[float, float], ...]
    circles: tuple[SequenceCircle, ...]
    connectors: tuple[tuple[tuple[float, float], ...], ...]
    seam_connector: tuple[tuple[float, float], ...] | None
    engagement_deg: tuple[float, ...]
    cap_deg: float
    tool_radius: float
    coverage: CoverageVerdict | None

    @classmethod
    def build(
        cls,
        *,
        workload: str,
        title: str,
        boundary: tuple[tuple[float, float], ...],
        circles: tuple[SequenceCircle, ...],
        connectors: tuple[tuple[tuple[float, float], ...], ...],
        seam_connector: tuple[tuple[float, float], ...] | None,
        engagement_deg: tuple[float, ...],
        cap_deg: float,
        tool_radius: float,
        coverage: CoverageVerdict | None,
    ) -> SequencePlotData:
        """Validate the order invariants.

        Raises:
            EmptySequenceError: no circles.
            SequencePlotError: engagement count is not one per successor pair,
                there is more than one connector per successor pair, or the
                tool radius is not positive.
        """
        if not circles:
            raise EmptySequenceError("A sequence plot needs at least one circle.")
        if len(engagement_deg) != len(circles) - 1:
            raise SequencePlotError(f"{len(engagement_deg)} engagements for {len(circles)} circles; expected one per successor pair.")
        if len(connectors) > len(circles) - 1:
            raise SequencePlotError(f"{len(connectors)} connectors for {len(circles)} circles; at most one per successor pair.")
        if not tool_radius > 0.0:
            raise SequencePlotError(f"tool radius {tool_radius!r} mm is not positive.")
        return cls(workload, title, boundary, circles, connectors, seam_connector, engagement_deg, cap_deg, tool_radius, coverage)

    def to_json(self) -> str:
        coverage = self.coverage
        return json.dumps(
            {
                "version": DATA_VERSION,
                "workload": self.workload,
                "title": self.title,
                "boundary": self.boundary,
                "circles": [(c.x, c.y, c.radius, c.phase, c.run) for c in self.circles],
                "connectors": self.connectors,
                "seam_connector": self.seam_connector,
                # JSON has no NaN: an unresolved engagement is stored as null.
                "engagement_deg": [None if math.isnan(angle) else angle for angle in self.engagement_deg],
                "cap_deg": self.cap_deg,
                "tool_radius": self.tool_radius,
                "coverage": None if coverage is None else {"covered": coverage.covered, "witnesses": coverage.witnesses, "piece_kinds": coverage.piece_kinds},
            }
        )

    @classmethod
    def from_json(cls, text: str) -> SequencePlotData:
        """Read a cached sequence.

        Raises:
            SequenceDataVersionError: the file's version is not ``DATA_VERSION``.
        """
        raw = json.loads(text)
        if raw.get("version") != DATA_VERSION:
            raise SequenceDataVersionError(f"version {raw.get('version')!r} is not {DATA_VERSION}.")
        seam = raw["seam_connector"]
        coverage = raw["coverage"]
        return cls.build(
            workload=raw["workload"],
            title=raw["title"],
            boundary=tuple((float(x), float(y)) for x, y in raw["boundary"]),
            circles=tuple(SequenceCircle(float(x), float(y), float(r), float(p), int(run)) for x, y, r, p, run in raw["circles"]),
            connectors=tuple(tuple((float(x), float(y)) for x, y in line) for line in raw["connectors"]),
            seam_connector=None if seam is None else tuple((float(x), float(y)) for x, y in seam),
            engagement_deg=tuple(math.nan if angle is None else float(angle) for angle in raw["engagement_deg"]),
            cap_deg=float(raw["cap_deg"]),
            tool_radius=float(raw["tool_radius"]),
            coverage=None
            if coverage is None
            else CoverageVerdict.build(
                covered=bool(coverage["covered"]),
                witnesses=tuple((float(x), float(y)) for x, y in coverage["witnesses"]),
                piece_kinds=tuple(str(kind) for kind in coverage["piece_kinds"]),
            ),
        )


def emitted_engagement(
    case: HeldReferenceCase,
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
) -> tuple[tuple[float, ...], float | None]:
    """Standard predecessor engagement in emitted order, degrees: (successor pairs, seam pair or None).

    NaN marks an unresolved model (missing evidence is shown as a gap, never a
    passing angle). The seam pair exists for a closed path of two or more
    circles.
    """
    pairs = connector_pairs(circles, transitions)
    # A lone circle's closing lap returns to itself: no seam pair to report.
    seam_pair = len(circles) > 1 and len(pairs) == len(circles)
    measured: list[float] = []
    for first, second, _ in pairs[: len(circles) - 1 + seam_pair]:
        try:
            angle = math.degrees(float(maximum_predecessor_engagement(_paper_candidate(first), _paper_candidate(second), case.tool_radius)))
        except StandardPlacementFragmentationError:
            angle = math.nan
        measured.append(angle)
    return tuple(measured[: len(circles) - 1]), (measured[-1] if seam_pair else None)


def plot_data(
    *,
    workload: str,
    title: str,
    case: HeldReferenceCase,
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    transitions: tuple[Figure5BoundaryTransition, ...],
    engagement_deg: tuple[float, ...],
    coverage: CoverageVerdict | None,
) -> SequencePlotData:
    """Report an emitted path as plot doubles; a closed path's last transition is the seam connector."""
    closed = len(circles) > 1 and len(transitions) == len(circles)
    lines = tuple(tuple((float(p.x), float(p.y)) for p in transition.samples) for transition in transitions)
    return SequencePlotData.build(
        workload=workload,
        title=title,
        boundary=tuple((float(p.x), float(p.y)) for p in case.projection.points),
        circles=tuple(
            SequenceCircle(
                float(c.center.x),
                float(c.center.y),
                float(c.radius.value),
                math.atan2(float(c.contact_point.y) - float(c.center.y), float(c.contact_point.x) - float(c.center.x)),
                int(c.run_id),
            )
            for c in circles
        ),
        connectors=lines[:-1] if closed else lines,
        seam_connector=lines[-1] if closed else None,
        engagement_deg=engagement_deg,
        cap_deg=float(case.tea_cap),
        tool_radius=float(case.tool_radius.value),
        coverage=coverage,
    )


def order_colormap() -> Colormap:
    """The sequential order ramp: a lightness-monotone window of 'magma'."""
    return ListedColormap(colormaps["magma"](np.linspace(RAMP_LOW, RAMP_HIGH, RAMP_STEPS)), name="emission_order")


class ColourMode(Enum):
    """What the order ramp spans.

    ORDER: the whole path, circle 0 darkest, the last circle lightest.
    RUN: each guide run on its own, its first circle darkest and its last
    lightest, so every fan reads as a dark-to-light sweep (Held–Pfeiffer
    Figure 5(b) look). Runs are maximal blocks of consecutive circles that
    share a ``run`` id.
    """

    ORDER = "order"
    RUN = "run"


COLOURBAR_LABELS = {
    ColourMode.ORDER: "emission order (circle index)",
    ColourMode.RUN: "position within guide run (first → last circle)",
}


def order_positions(circles: tuple[SequenceCircle, ...], mode: ColourMode) -> npt.NDArray[np.float64]:
    """Ramp position in ``[0, 1]`` per circle, non-decreasing in emission index (within a run for RUN)."""
    count = len(circles)
    if mode is ColourMode.ORDER:
        return np.arange(count, dtype=np.float64) / max(count - 1, 1)
    positions = np.zeros(count, dtype=np.float64)
    start = 0
    for index in range(1, count + 1):
        if index == count or circles[index].run != circles[start].run:
            length = index - start
            positions[start:index] = np.arange(length, dtype=np.float64) / max(length - 1, 1)
            start = index
    return positions


def order_colors(circles: tuple[SequenceCircle, ...], mode: ColourMode) -> npt.NDArray[np.float64]:
    """RGBA per circle on the order ramp, dark (first) to light (last)."""
    return np.asarray(order_colormap()(order_positions(circles, mode)), dtype=np.float64)


class SweptStrokes(LineCollection):
    """Centre paths stroked to a millimetre width: a round-capped, round-joined stroke of width 2r is path (+) D_r.

    ``half_widths_mm`` is r per path; ``inset_pt`` narrows a member by a
    display width in points (both rims of a colour stroke), so the member
    stays inside its band. A member whose inset would leave it narrower than
    ``MIN_COLOUR_RIMS`` rims keeps the full band width instead (its rim
    stroke is then hidden under it). Widths are resolved at draw time from the
    axes' data scale, so the band is r in millimetres whatever the layout.
    Members are drawn in sequence, later over earlier; no member is wider
    than 2r.
    """

    def __init__(self, paths: Sequence[npt.NDArray[np.float64]], half_widths_mm: npt.NDArray[np.float64], inset_pt: npt.NDArray[np.float64], colors: Any) -> None:
        super().__init__(list(paths), colors=colors, capstyle="round", joinstyle="round")
        self._widths_mm = 2.0 * half_widths_mm
        self._inset_pt = inset_pt

    def draw(self, renderer: RendererBase) -> None:
        axes = self.axes
        if axes is None:
            raise SequencePlotError("SweptStrokes draws only inside an Axes.")
        origin, unit = axes.transData.transform([(0.0, 0.0), (1.0, 0.0)])
        pixels_per_mm = abs(float(unit[0] - origin[0]))
        points_per_pixel = 1.0 / float(np.asarray(renderer.points_to_pixels(1.0), dtype=np.float64))
        band_pt = self._widths_mm * pixels_per_mm * points_per_pixel
        rim_pt = self._inset_pt / 2.0
        narrowed = band_pt - self._inset_pt
        self._require_unclipped_bands(renderer, 0.5 * float(self._widths_mm.max()) * pixels_per_mm)
        self.set_linewidth(list(np.where(narrowed >= MIN_COLOUR_RIMS * rim_pt, narrowed, band_pt)))
        super().draw(renderer)

    def _require_unclipped_bands(self, renderer: RendererBase, half_band_px: float) -> None:
        """Fail if Agg's canvas clip could remove band paint inside the axes.

        Agg cuts an unfilled path at the canvas edge (padded by 1 px) before
        stroking it, so a centre path leaving the canvas is stroked only up to
        the cut; the missing paint lies within half a band of that edge. It is
        invisible when every path stays on the canvas or the axes keep at least
        half a band from every canvas edge; otherwise a zoom would show page
        colour inside a band. Display geometry only; no geometric decision.

        Raises:
            StrokeClippedAtCanvasError: a path leaves the canvas and the axes
                lie within half a band of a canvas edge.
        """
        axes = self.axes
        if axes is None:
            raise SequencePlotError("SweptStrokes draws only inside an Axes.")
        width, height = renderer.get_canvas_width_height()
        vertices = axes.transData.transform(np.concatenate([path.vertices for path in self.get_paths()]))
        leaves = bool(np.any(vertices < -CANVAS_CLIP_PAD_PX) or np.any(vertices[:, 0] > width + CANVAS_CLIP_PAD_PX) or np.any(vertices[:, 1] > height + CANVAS_CLIP_PAD_PX))
        box = axes.bbox
        gap = min(box.x0, box.y0, width - box.x1, height - box.y1)
        if leaves and gap < half_band_px:
            raise StrokeClippedAtCanvasError(
                f"band half-width {half_band_px:.0f} px exceeds the {gap:.0f} px between the axes and the canvas edge while centre paths leave the canvas; "
                "render the full view and crop instead of zooming the axes."
            )


def swept_strokes(data: SequencePlotData, colors: npt.NDArray[np.float64]) -> SweptStrokes:
    """Every motion in emission order (circle i, connector i -> i+1, ..., seam), each a rim stroke then a colour stroke.

    The rim stroke is exactly the band (width 2r) in a darker shade of the
    band's colour; the colour stroke on top is 2r minus both rims, so the rim
    shows only as the band's own edge and nothing is painted outside it.
    """
    paths: list[npt.NDArray[np.float64]] = []
    rgba: list[npt.NDArray[np.float64]] = []

    def emit(path: npt.NDArray[np.float64], color: npt.NDArray[np.float64]) -> None:
        rim = np.array((*(RIM_SHADE * color[:3]), color[3]))
        paths.extend((path, path))
        rgba.extend((rim, color))

    for index, turn in enumerate(_turns(data.circles)):
        emit(turn, colors[index])
        if index < len(data.connectors):
            emit(np.asarray(data.connectors[index], dtype=np.float64), colors[index])
    if data.seam_connector is not None:
        emit(np.asarray(data.seam_connector, dtype=np.float64), colors[-1])
    count = len(paths)
    inset = np.tile(np.array((0.0, 2.0 * DISK_EDGE_PT)), count // 2)
    return SweptStrokes(paths, np.full(count, data.tool_radius), inset, np.array(rgba))


def render_sequence(
    data: SequencePlotData,
    output: Path,
    mode: ColourMode = ColourMode.ORDER,
    *,
    subtitle: str | None = None,
    notes: tuple[str, ...] = (),
) -> None:
    """Draw the order-coloured swept regions over the order-coloured engagement trace.

    Args:
        subtitle: The metrics line under the title; default: counts from ``data``.
        notes: Caption lines above the certificate verdict and colour legend.
    """
    count = len(data.circles)
    colors = order_colors(data.circles, mode)
    extent = np.ptp(np.array(data.boundary), axis=0)
    path_height = float(np.clip(PATH_PANEL_WIDTH_IN * extent[1] / extent[0], PATH_PANEL_MIN_IN, PATH_PANEL_MAX_IN))
    lines = caption_lines(data, mode, notes)
    footer = FOOTER_IN + NOTE_LINE_IN * (len(lines) - 1)
    header = HEADER_IN + (VERDICT_LINE_IN if data.coverage is not None else 0.0)
    height = header + path_height + PANEL_GAP_IN + ENGAGEMENT_PANEL_IN + footer
    figure = Figure(figsize=(FIGURE_WIDTH_IN, height), facecolor=PAGE_COLOR, dpi=DPI)
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 2, height_ratios=(path_height, ENGAGEMENT_PANEL_IN), width_ratios=(1, 0.018), wspace=0.03)
    path_axis = figure.add_subplot(grid[0, 0])
    engagement_axis = figure.add_subplot(grid[1, 0])
    colorbar_axis = figure.add_subplot(grid[:, 1])
    draw_sweep(path_axis, data, colors)
    draw_engagement(engagement_axis, data, colors)
    add_order_colorbar(figure, colorbar_axis, mode, count)
    finite = [angle for angle in data.engagement_deg if math.isfinite(angle)]
    over = sum(1 for angle in finite if angle > data.cap_deg)
    unresolved = len(data.engagement_deg) - len(finite)
    seam = "seam connector emitted" if data.seam_connector is not None else "no seam connector"
    figure.suptitle(data.title, x=0.06, y=1 - 0.12 / height, ha="left", va="top", fontsize=19, fontweight="bold")
    if subtitle is None:
        subtitle = (
            f"{count:,} circles in emission order · {_run_count(data.circles):,} guide runs · {over:,} over cap · {unresolved} unresolved"
            f" · max known {max(finite, default=0.0):.2f}° · {seam}"
        )
    figure.text(0.06, 1 - 0.55 / height, subtitle)
    draw_verdict(figure, data.coverage, x=0.06, y=1 - 0.82 / height, join_mm=data.tool_radius)
    for row, line in enumerate(reversed(lines)):
        figure.text(0.06, (0.15 + NOTE_LINE_IN * row) / height, line, fontsize=9)
    figure.subplots_adjust(left=0.06, right=0.93, top=1 - header / height, bottom=footer / height, hspace=PANEL_GAP_IN / ((path_height + ENGAGEMENT_PANEL_IN) / 2))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=figure.dpi, facecolor=figure.get_facecolor())


def caption_lines(data: SequencePlotData, mode: ColourMode, notes: tuple[str, ...]) -> tuple[str, ...]:
    """Footer lines top to bottom: notes, the certificate verdict (if any), what colour means."""
    return (*notes, *(() if data.coverage is None else (data.coverage.caption(data.tool_radius),)), colour_legend_text(mode))


def draw_verdict(figure: Figure, coverage: CoverageVerdict | None, *, x: float, y: float, join_mm: float) -> None:
    """The certificate verdict in the header, left-aligned at (x, y) in figure coordinates."""
    if coverage is None:
        return
    figure.text(
        x,
        y,
        "CERTIFIED COVERED" if coverage.covered else f"NOT COVERED · {len(coverage.sites(join_mm))} uncut site(s), {len(coverage.witnesses)} exact witnesses",
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=COVERED_COLOR if coverage.covered else CAP_COLOR,
    )


def colour_legend_text(mode: ColourMode) -> str:
    """The caption line that says what colour and paint order mean."""
    meaning = "emission order over the whole path" if mode is ColourMode.ORDER else "position within each guide run"
    return f"Colour = {meaning}. Swept bands (circle: annulus R ± r; connector: capsule) painted in emission order, later on top; white inside the pocket = uncut."


def add_order_colorbar(figure: Figure, axis: Axes, mode: ColourMode, count: int) -> None:
    """The order ramp's colourbar on ``axis``."""
    upper = max(count - 1, 1) if mode is ColourMode.ORDER else 1
    colorbar = figure.colorbar(ScalarMappable(norm=Normalize(vmin=0, vmax=upper), cmap=order_colormap()), cax=axis)
    colorbar.set_label(COLOURBAR_LABELS[mode])
    if mode is ColourMode.RUN:
        colorbar.set_ticks([0, 1], labels=["first", "last"])
    colorbar.outline.set_visible(False)


def witness_label(index: int) -> str:
    """Plot label of the index-th witness: A, B, ..., then W27, W28, ..."""
    return chr(ord("A") + index) if index < 26 else f"W{index + 1}"


def _run_count(circles: tuple[SequenceCircle, ...]) -> int:
    return 1 + sum(1 for first, second in zip(circles, circles[1:]) if first.run != second.run)


def _turns(circles: tuple[SequenceCircle, ...]) -> list[npt.NDArray[np.float64]]:
    # Each turn starts at the circle's contact (its phase), as emitted.
    parameters = np.linspace(0.0, math.tau, TURN_DISPLAY_SAMPLES + 1)
    return [np.array((circle.x, circle.y)) + circle.radius * np.column_stack((np.cos(parameters + circle.phase), np.sin(parameters + circle.phase))) for circle in circles]


def _outside_of(ring: npt.NDArray[np.float64]) -> MplPath:
    """A compound path covering a large frame around ``ring`` minus the ring's interior (opposite orientation: a hole)."""
    low = ring.min(axis=0)
    high = ring.max(axis=0)
    pad = float(np.ptp(ring, axis=0).max())
    frame = np.array(((low[0] - pad, low[1] - pad), (high[0] + pad, low[1] - pad), (high[0] + pad, high[1] + pad), (low[0] - pad, high[1] + pad)))
    hole = ring[:-1]
    signed_area = 0.5 * float(np.sum(hole[:, 0] * np.roll(hole[:, 1], -1) - np.roll(hole[:, 0], -1) * hole[:, 1]))
    if signed_area > 0.0:
        hole = hole[::-1]
    vertices = np.concatenate((frame, frame[:1], hole, hole[:1]))
    codes = [MplPath.MOVETO, *([MplPath.LINETO] * 3), MplPath.CLOSEPOLY, MplPath.MOVETO, *([MplPath.LINETO] * (len(hole) - 1)), MplPath.CLOSEPOLY]
    return MplPath(vertices, codes)


def draw_sweep(axis: Axes, data: SequencePlotData, colors: npt.NDArray[np.float64], *, legend: bool = True) -> None:
    """The swept regions in emission order, the pocket outline, start, seam and certificate witnesses."""
    boundary = np.array(data.boundary + data.boundary[:1])
    axis.set_facecolor(PAGE_COLOR)
    axis.set(aspect="equal", xlabel="World X / mm", ylabel="World Y / mm")
    axis.grid(color="#9aa09c", linestyle=":", linewidth=0.6, alpha=0.6)
    axis.set_axisbelow(True)
    margin = 0.03 * float(np.ptp(boundary, axis=0).max())
    axis.set_xlim(boundary[:, 0].min() - margin, boundary[:, 0].max() + margin)
    axis.set_ylim(boundary[:, 1].min() - margin, boundary[:, 1].max() + margin)
    # A collection draws its members in sequence: each later band covers the
    # earlier ones, so what is left of a band is the crescent only it reached.
    sweeps = swept_strokes(data, colors)
    sweeps.set_zorder(2)
    axis.add_collection(sweeps)
    (outline,) = axis.plot(boundary[:, 0], boundary[:, 1], color=BOUNDARY_COLOR, linewidth=BOUNDARY_LINE_PT, solid_joinstyle="round", zorder=3, label="Pocket boundary")
    outline.set_clip_path(PathPatch(_outside_of(boundary), transform=axis.transData))
    if data.seam_connector is not None:
        # The seam is often a short arc; its end marker keeps it findable.
        seam = np.array(data.seam_connector)
        axis.plot(seam[:, 0], seam[:, 1], color=SEAM_COLOR, linewidth=1.4, zorder=4, label="Seam connector centre (last to first)")
        axis.scatter(seam[0, 0], seam[0, 1], marker="D", s=40, color=SEAM_COLOR, edgecolor="white", linewidth=0.8, zorder=6, label="Seam start (last circle)")
    angles = np.array(data.engagement_deg)
    over = np.flatnonzero(np.isfinite(angles) & (angles > data.cap_deg)) + 1
    if len(over):
        # The successor circle of each over-cap pair, marked at its centre.
        centres = np.array([(data.circles[index].x, data.circles[index].y) for index in over])
        axis.scatter(centres[:, 0], centres[:, 1], marker="x", s=22, color=CAP_COLOR, linewidth=1.1, zorder=6, label="Over-cap successor (centre)")
    first = data.circles[0]
    start = (first.x + first.radius * math.cos(first.phase), first.y + first.radius * math.sin(first.phase))
    axis.scatter(*start, s=80, color="white", edgecolor="#1c1c1c", linewidth=1.2, zorder=5, label="Start")
    if data.coverage is not None and data.coverage.witnesses:
        witnesses = np.array(data.coverage.witnesses)
        axis.scatter(witnesses[:, 0], witnesses[:, 1], s=6, color=WITNESS_COLOR, zorder=7)
        sites = np.array([site[:2] for site in data.coverage.sites(data.tool_radius)])
        axis.scatter(sites[:, 0], sites[:, 1], s=150, facecolor="none", edgecolor=WITNESS_COLOR, linewidth=1.6, zorder=7, label="Uncut site (exact witnesses)")
        for index, (x, y) in enumerate(sites):
            axis.annotate(witness_label(index), xy=(x, y), xytext=(8, 6), textcoords="offset points", color=WITNESS_COLOR, fontsize=11, fontweight="bold", zorder=7)
    if legend:
        axis.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=5, frameon=False, fontsize=9, borderaxespad=0.2)


def draw_engagement(
    axis: Axes,
    data: SequencePlotData,
    colors: npt.NDArray[np.float64],
    *,
    ylabel: str = "Engagement / degrees",
    ylim: tuple[float, float] = (0.0, 190.0),
) -> None:
    """The engagement trace on the emission-order axis, each segment coloured like the circle it ends at."""
    successors = np.arange(1, len(data.circles))
    angles = np.array(data.engagement_deg)
    if len(angles):
        segments = np.stack((np.column_stack((successors[:-1], angles[:-1])), np.column_stack((successors[1:], angles[1:]))), axis=1)
        # Segment i ends at successor i + 1: coloured like that circle.
        axis.add_collection(LineCollection(list(segments), colors=colors[2:], linewidths=0.9))
        unresolved = successors[np.isnan(angles)]
        if len(unresolved):
            axis.scatter(unresolved, np.full(len(unresolved), ylim[1] - 5.0), marker="x", color="#6e39a6", s=24, label="Unresolved model")
        over = np.isfinite(angles) & (angles > data.cap_deg)
        if over.any():
            axis.scatter(successors[over], angles[over], color=CAP_COLOR, s=9, zorder=3, label="Over cap")
    axis.axhline(data.cap_deg, color=CAP_COLOR, linestyle="--", linewidth=1, label=f"{data.cap_deg:g}° cap")
    axis.set(xlabel="Emission order (circle index)", ylabel=ylabel, ylim=ylim, xlim=(0, max(len(data.circles) - 1, 1)))
    axis.set_facecolor(PAGE_COLOR)
    axis.grid(color="#9aa09c", linestyle=":", linewidth=0.6, alpha=0.6)
    # Above the panel, so no marker of the trace is hidden.
    axis.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, frameon=False, fontsize=9, borderaxespad=0.2)


def render_side_by_side(old: Path, new: Path, output: Path) -> None:
    """Place two rendered PNGs next to each other for review."""
    figure = Figure(figsize=(26, 11), facecolor=PAGE_COLOR)
    FigureCanvasAgg(figure)
    for index, (path, label) in enumerate(((old, f"old · {old.name}"), (new, f"new · {new.name}"))):
        axis = figure.add_subplot(1, 2, index + 1)
        axis.imshow(imread(str(path)))
        axis.set_title(label, fontsize=12)
        axis.axis("off")
    figure.subplots_adjust(left=0.005, right=0.995, top=0.95, bottom=0.005, wspace=0.01)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=110, facecolor=figure.get_facecolor())
