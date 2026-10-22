"""Qualification overlay of a published Held toolpath: what the tool removes, what it leaves, where it gouges.

The published Held figures draw tool-centre circles only, so they can neither
show nor hide uncut stock or gouges. This renderer draws the material view:

- the analytic pocket boundary (lines and arcs of the reference case);
- the swept area, light fill: the annulus of every full circle and the capsule
  of every connector piece. RENDERED GEOMETRY, a raster for reporting only;
- uncut material in red, split into material inside the machinable target
  T = Omega (+) D_r (a placement defect) and material outside it (pocket tips
  and chord slivers that no disk centred in the declared Omega reaches);
- exact uncovered witnesses (red crosses) from the coverage certificate;
- exact gouges from :mod:`benchmarks.held_gouge_audit` in the second accent,
  with reported depths on the worst few;
- zoom insets on the largest witness clusters and the deepest gouge;
- a caption that states what is EXACT (decided) and what is DIAGNOSTIC.

Decisions vs reporting. Two things here are exact decisions, and both are made
by their owning modules, never here: coverage, by
:meth:`CoverageCertificate.build` on the whole closed path (every full circle
and every connector, the seam connector included), which either certifies
T = Omega (+) D_r covered or raises with one exact witness per exposed edge;
and the gouge verdict of every element, by the gouge audit. Everything else,
the swept fill, the uncut regions, the inset framing, the gouge locations and
every printed number, is diagnostic and never decides.

Run as ``pixi run held-qualification-figures`` (generates both paths live).
"""

from __future__ import annotations

import argparse
import math
import textwrap
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
from matplotlib.path import Path as MplPath

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_gouge_audit import REPORTING_SPLIT_MM
from benchmarks.held_gouge_audit import AuditCircle
from benchmarks.held_gouge_audit import AuditConnector
from benchmarks.held_gouge_audit import ElementClearance
from benchmarks.held_gouge_audit import HeldGougeAudit
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_motion_coverage import connector_pairs
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_toolpath_benchmark import GeneratedToolpath
from benchmarks.held_toolpath_benchmark import generate
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import UncoveredWitness
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tools.performance_style import FIGURE_DPI
from tools.performance_style import FOOTER_FONT_SIZE
from tools.performance_style import INK
from tools.performance_style import MUTED
from tools.performance_style import SURFACE
from tools.performance_style import TITLE_FONT_SIZE

ROOT = Path(__file__).resolve().parents[2]
IMAGES = ROOT / "docs" / "assets" / "images"

# Accents on the shared surface. Swept fill: step 100 of the sequential blue
# ramp (tools.performance_style.RAMP). Uncut: red. Gouge: violet, the second
# accent, distinct from red and from the blue fill under colour-vision deficits.
SWEPT = "#cde2fb"
UNCUT = "#c8312a"
UNREACHABLE = "#e9b8b2"
GOUGE = "#7a3db8"
OMEGA = "#9a9a94"

# Diagnostic raster pitch of the main view, mm: 1/20 of the 1 mm Held tool
# radius, so the 0.03 mm^2 Monstera slivers (about 0.3 mm across) span ~6 px.
MAIN_PIXEL_MM = 0.05
# Inset rasters are this many pixels across, whatever their window.
INSET_PIXELS = 480
# Smallest inset half-width, mm: frames a 0.3 mm sliver or a 0.05 mm-deep gouge
# lens with context of about one tool radius on each side of the wall.
INSET_MIN_HALF_MM = 0.6
# Gouge inset half-width, mm: the Held gouges are 0.02-0.05 mm deep, so a
# 0.1 mm half-window shows the gouged lens at 10-25 % of the inset height.
GOUGE_INSET_HALF_MM = 0.1
# Inset half-width = this factor times the cluster's half-extent (padding).
INSET_PAD_FACTOR = 2.5
UNCUT_INSETS = 2
# Single-linkage distance grouping exact witnesses into clusters for inset
# framing, mm (reporting only): one tool radius; the Monstera slivers are
# ~0.3 mm across and > 10 mm apart, so the grouping is unambiguous there.
WITNESS_CLUSTER_LINK_MM = 1.0
# Gouge depth labels on the main view: the worst few, at least this far apart
# (mm) so neighbouring elements of one gouge cluster do not stack labels.
LABELLED_GOUGES = 3
GOUGE_LABEL_SEPARATION_MM = 8.0
# Angular samples along a gouging circle when locating its deepest point
# (reporting only): 2 pi R / 4096 <= 0.011 mm along the largest Held circle.
GOUGE_LOCATION_SAMPLES = 4096
GOUGE_SEGMENT_SAMPLES = 512
# Pocket-boundary rendering: an arc is drawn as a polyline of this many
# segments per radian (chord sagitta <= rho / 8 * (1/64)^2 ~ 0.002 mm at rho = 7).
ARC_SEGMENTS_PER_RAD = 64
# Figure layout, inches (matplotlib's unit): main-view height, page margin,
# gap before the inset column, inset side, note block under an inset, title
# band, legend band, footer line pitch at FOOTER_FONT_SIZE, gap above the footer.
MAIN_HEIGHT_IN = 7.0
MARGIN_IN = 0.25
COLUMN_GAP_IN = 0.35
INSET_SIDE_IN = 2.2
INSET_NOTE_IN = 0.55
INSET_NOTE_PAD_IN = 0.06
TITLE_IN = 0.6
LEGEND_IN = 0.45
FOOTER_LINE_IN = 0.15
FOOTER_GAP_IN = 0.05
# Characters of FOOTER_FONT_SIZE (7.5 pt) DejaVu Sans per inch of line, about
# 72 / (0.55 * 7.5); wrapping is a layout choice, measured on the rendered PNGs.
FOOTER_CHARS_PER_INCH = 16.5
# Inset scale bars: the largest of these lengths (mm) that fits in this
# fraction of the inset window.
SCALE_BAR_STEPS_MM = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
SCALE_BAR_FRACTION = 0.3
# Gouge depth labels sit this fraction of the plot span off their wall point.
LABEL_OFFSET_X = 0.07
LABEL_OFFSET_Y = 0.05
# Bytes per RGBA pixel channel range.
CHANNEL_MAX = 255


class QualificationOverlayError(BenchmarkError):
    """Base error of the qualification overlay."""


class UnknownOverlayWorkloadError(QualificationOverlayError):
    """The requested workload has no overlay definition."""


class InconsistentConnectorFanError(QualificationOverlayError):
    """A connector's within-fan flag disagrees with the exact sites of the circles it joins."""


class CoverageStatus(Enum):
    """The coverage certificate's exact verdict on one closed path."""

    CERTIFIED_COVERED = "certified_covered"
    NOT_COVERED = "not_covered"


@dataclass(frozen=True)
class OverlayCircle:
    """A full-circle cut: the tool centre runs once round (centre, radius)."""

    centre: Point2[WorldXY]
    radius: GuideRadius


@dataclass(frozen=True)
class OverlayInput:
    """Everything the overlay draws or decides on, independent of the generator.

    Attributes:
        name: Short workload name, used in the output file name.
        title: Figure title.
        pocket: Analytic pocket boundary primitives, one CCW cycle.
        domain: Declared centre domain Omega (the connectors follow its outer loop).
        tool_radius: Tool radius r.
        circles: Full-circle cuts in emission order.
        connectors: Emitted connector polylines (tool-centre samples) in
            emission order, the seam connector last when the path closes with
            one. They are the motion the gouge audit decides on and the raster
            draws.
        connector_sweeps: The same connectors as exact CCW site intervals of
            Omega, the regions the coverage certificate decides on; a
            zero-progress fan connector sweeps nothing new and has none.
    """

    name: str
    title: str
    pocket: tuple[ReferenceLine | ReferenceArc, ...]
    domain: CentreDomain
    tool_radius: ToolRadius
    circles: tuple[OverlayCircle, ...]
    connectors: tuple[tuple[Point2[WorldXY], ...], ...]
    connector_sweeps: tuple[BoundaryConnectorSweep, ...]

    @property
    def omega(self) -> tuple[Point2[WorldXY], ...]:
        """The outer loop of Omega."""
        return self.domain.outer

    @classmethod
    def from_generated(cls, path: GeneratedToolpath, name: str, title: str) -> OverlayInput:
        """Adapt a live generator output: Omega is the loop its connectors follow, the seam connector included.

        Raises:
            InvalidCoverageMotionError: the connector count is neither open nor
                closed, or a connector does not join its circles' contacts.
            InconsistentConnectorFanError: a connector is flagged within-fan
                but joins different exact sites, or the converse.
        """
        domain = CentreDomain.build(list(path.boundary))
        sweeps: list[BoundaryConnectorSweep] = []
        for first, second, transition in connector_pairs(path.circles, path.transitions):
            same_site = first.boundary_site == second.boundary_site and first is not second
            if transition.within_fan != same_site:
                raise InconsistentConnectorFanError(f"connector within_fan={transition.within_fan} joins sites {first.boundary_site} and {second.boundary_site}")
            if same_site:
                continue
            sweeps.append(BoundaryConnectorSweep.build(domain, _site(first.boundary_site), _site(second.boundary_site)))
        return cls(
            name=name,
            title=title,
            pocket=tuple(path.case.boundary.primitives),
            domain=domain,
            tool_radius=path.case.tool_radius,
            circles=tuple(OverlayCircle(c.center, c.radius) for c in path.circles),
            connectors=tuple(t.samples for t in path.transitions),
            connector_sweeps=tuple(sweeps),
        )


def _site(site: ProjectionBoundarySite) -> CentreBoundarySite:
    """The generator's exact contact site as a certificate site on loop 0 (the outer loop)."""
    return CentreBoundarySite.build(loop=0, segment=int(site.segment_id), parameter=site.parameter)


@dataclass(frozen=True)
class CoverageQualification:
    """The certificate's exact verdict on the whole closed path.

    Attributes:
        status: CERTIFIED_COVERED when ``CoverageCertificate.build`` returns,
            NOT_COVERED when it raises ``UncoveredPocketError``.
        witnesses: One exact uncovered point per exposed edge (empty when covered).
        undecided_edge_count: Exposed edges the certificate could not decide
            (zero on a polygonal Omega without vertex disks).
        regions: Swept regions the certificate decided on (circles + connector sweeps).
    """

    status: CoverageStatus
    witnesses: tuple[UncoveredWitness, ...]
    undecided_edge_count: int
    regions: int

    @property
    def witnesses_mm(self) -> tuple[tuple[float, float], ...]:
        """Witness points rounded to binary64, for drawing only."""
        return tuple((float(x), float(y)) for x, y in (w.uncovered_point_mm for w in self.witnesses))


@dataclass(frozen=True)
class GougeMark:
    """One gouging element with its reported depth and the wall point it cuts (reporting only)."""

    verdict: ElementClearance
    location_mm: tuple[float, float]


@dataclass(frozen=True)
class Inset:
    """One zoom window: a label, its centre and half-width in mm, and why it was chosen."""

    label: str
    centre: tuple[float, float]
    half_width: float
    note: str


FloatArray = npt.NDArray[np.float64]
BoolArray = npt.NDArray[np.bool_]


# ---------------------------------------------------------------------------
# Exact decisions
# ---------------------------------------------------------------------------


def _segments(connectors: Sequence[Sequence[Point2[WorldXY]]]) -> FloatArray:
    """Connector pieces as rows (ax, ay, bx, by); a one-sample connector is a zero-length piece."""
    rows: list[tuple[float, float, float, float]] = []
    for samples in connectors:
        if len(samples) == 1:
            rows.append((samples[0].x, samples[0].y, samples[0].x, samples[0].y))
        rows.extend((a.x, a.y, b.x, b.y) for a, b in zip(samples, samples[1:]))
    return np.array(rows, dtype=np.float64).reshape(-1, 4)


def qualify_coverage(overlay: OverlayInput) -> CoverageQualification:
    """The coverage certificate's exact verdict on every circle and every connector sweep.

    No decision is made here: ``CoverageCertificate.build`` either certifies
    T = Omega (+) D_r covered or raises ``UncoveredPocketError`` with its exact
    witnesses, and this only records which.
    """
    cuts = [FullCircleSweep.build(c.centre, c.radius) for c in overlay.circles]
    regions = len(cuts) + len(overlay.connector_sweeps)
    try:
        CoverageCertificate.build(overlay.domain, overlay.tool_radius, cuts, overlay.connector_sweeps)
    except UncoveredPocketError as uncovered:
        return CoverageQualification(CoverageStatus.NOT_COVERED, uncovered.witnesses, uncovered.undecided_edge_count, regions)
    return CoverageQualification(CoverageStatus.CERTIFIED_COVERED, (), 0, regions)


def gouge_audit(overlay: OverlayInput) -> HeldGougeAuditReport:
    """Exact gouge verdict of every element against the analytic pocket."""
    audit = HeldGougeAudit.build(overlay.pocket, overlay.tool_radius.value)
    circles = [AuditCircle(c.centre, c.radius.value) for c in overlay.circles]
    return audit.audit(circles, [AuditConnector(samples) for samples in overlay.connectors])


# ---------------------------------------------------------------------------
# Reporting geometry (floats, never decides)
# ---------------------------------------------------------------------------


def _arc_frame(arc: ReferenceArc) -> tuple[float, float, float, float, float]:
    """(ox, oy, rho, start angle, sweep) of an arc, rho = |start - centre|."""
    ox, oy = arc.centre.x, arc.centre.y
    return ox, oy, math.hypot(arc.start.x - ox, arc.start.y - oy), math.atan2(arc.start.y - oy, arc.start.x - ox), float(arc.sweep)


def pocket_polyline(pocket: Sequence[ReferenceLine | ReferenceArc]) -> FloatArray:
    """The analytic boundary sampled for drawing and for the diagnostic pocket mask."""
    points: list[tuple[float, float]] = []
    for primitive in pocket:
        if isinstance(primitive, ReferenceLine):
            points.append((primitive.start.x, primitive.start.y))
            continue
        ox, oy, rho, a0, sweep = _arc_frame(primitive)
        count = max(2, math.ceil(abs(sweep) * ARC_SEGMENTS_PER_RAD))
        angles = a0 + sweep * np.arange(count) / count
        points.extend(zip((ox + rho * np.cos(angles)).tolist(), (oy + rho * np.sin(angles)).tolist()))
    points.append(points[0])
    return np.array(points, dtype=np.float64)


def _nearest_on_primitive(primitive: ReferenceLine | ReferenceArc, p: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Distance from each row of ``p`` to the primitive and the nearest primitive point."""
    if isinstance(primitive, ReferenceLine):
        a = np.array([primitive.start.x, primitive.start.y])
        d = np.array([primitive.end.x, primitive.end.y]) - a
        t = np.clip(((p - a) @ d) / float(d @ d), 0.0, 1.0)
        q = a + t[:, None] * d
        return np.hypot(*(p - q).T), q
    ox, oy, rho, a0, sweep = _arc_frame(primitive)
    o = np.array([ox, oy])
    v = p - o
    phi = np.arctan2(v[:, 1], v[:, 0])
    rel = np.mod(phi - a0, math.tau) if sweep > 0 else np.mod(a0 - phi, math.tau)
    inside = rel <= abs(sweep)
    norm = np.hypot(v[:, 0], v[:, 1])
    radial = o + rho * v / norm[:, None]
    ends = np.array([[ox + rho * math.cos(a0), oy + rho * math.sin(a0)], [ox + rho * math.cos(a0 + sweep), oy + rho * math.sin(a0 + sweep)]])
    to_ends = np.hypot(p[:, None, 0] - ends[None, :, 0], p[:, None, 1] - ends[None, :, 1])
    nearest_end = ends[np.argmin(to_ends, axis=1)]
    q = np.where(inside[:, None], radial, nearest_end)
    return np.hypot(*(p - q).T), q


def gouge_marks(overlay: OverlayInput, report: HeldGougeAuditReport) -> tuple[GougeMark, ...]:
    """Locate each exact gouge at the wall point nearest its element (reporting only), deepest first."""
    marks: list[GougeMark] = []
    for item in report.gouges:
        primitive = overlay.pocket[item.primitive]
        if item.element == "circle":
            circle = overlay.circles[item.index]
            angles = np.linspace(0.0, math.tau, GOUGE_LOCATION_SAMPLES, endpoint=False)
            path = np.column_stack((circle.centre.x + circle.radius.value * np.cos(angles), circle.centre.y + circle.radius.value * np.sin(angles)))
        else:
            samples = overlay.connectors[item.index]
            piece = item.segment if item.segment is not None else 0
            a, b = samples[piece], samples[min(piece + 1, len(samples) - 1)]
            t = np.linspace(0.0, 1.0, GOUGE_SEGMENT_SAMPLES)
            path = np.column_stack((a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)))
        distance, nearest = _nearest_on_primitive(primitive, path)
        q = nearest[int(np.argmin(distance))]
        marks.append(GougeMark(item, (float(q[0]), float(q[1]))))
    return tuple(sorted(marks, key=lambda mark: -mark.verdict.depth_mm))


@dataclass(frozen=True)
class MaterialRaster:
    """Diagnostic raster of one window: pocket, machinable target, swept area (reporting only)."""

    extent: tuple[float, float, float, float]
    pixel_mm: float
    pocket: BoolArray
    target: BoolArray
    swept: BoolArray

    @property
    def uncut_in_target(self) -> BoolArray:
        return self.pocket & self.target & ~self.swept

    @property
    def uncut_unreachable(self) -> BoolArray:
        return self.pocket & ~self.target & ~self.swept

    @property
    def gouged(self) -> BoolArray:
        return self.swept & ~self.pocket

    def rgba(self) -> npt.NDArray[np.uint8]:
        """Surface, swept fill, uncut in T, uncut outside T, gouged; row 0 at the bottom."""
        image = np.empty((*self.pocket.shape, 4), dtype=np.uint8)
        image[...] = _rgba(SURFACE)
        image[self.pocket & self.swept] = _rgba(SWEPT)
        image[self.uncut_unreachable] = _rgba(UNREACHABLE)
        image[self.uncut_in_target] = _rgba(UNCUT)
        image[self.gouged] = _rgba(GOUGE)
        return image


def _rgba(colour: str) -> tuple[int, int, int, int]:
    rgb = matplotlib.colors.to_rgb(colour)
    return (round(rgb[0] * CHANNEL_MAX), round(rgb[1] * CHANNEL_MAX), round(rgb[2] * CHANNEL_MAX), CHANNEL_MAX)


def _paint_capsules(mask: BoolArray, xs: FloatArray, ys: FloatArray, segments: FloatArray, radius: float) -> None:
    pixel = xs[1] - xs[0] if xs.size > 1 else 1.0
    for ax, ay, bx, by in segments:
        i0, i1 = np.searchsorted(xs, [min(ax, bx) - radius - pixel, max(ax, bx) + radius + pixel])
        j0, j1 = np.searchsorted(ys, [min(ay, by) - radius - pixel, max(ay, by) + radius + pixel])
        if i0 >= i1 or j0 >= j1:
            continue
        gx, gy = np.meshgrid(xs[i0:i1] - ax, ys[j0:j1] - ay)
        dx, dy = bx - ax, by - ay
        length2 = dx * dx + dy * dy
        t = np.zeros_like(gx) if length2 == 0.0 else np.clip((gx * dx + gy * dy) / length2, 0.0, 1.0)
        mask[j0:j1, i0:i1] |= (gx - t * dx) ** 2 + (gy - t * dy) ** 2 <= radius * radius


def material_raster(overlay: OverlayInput, extent: tuple[float, float, float, float], pixel_mm: float) -> MaterialRaster:
    """Rasterise one window at pixel centres (diagnostic, reporting only)."""
    x0, x1, y0, y1 = extent
    xs = (x0 + pixel_mm * (np.arange(math.ceil((x1 - x0) / pixel_mm)) + 0.5)).astype(np.float64)
    ys = (y0 + pixel_mm * (np.arange(math.ceil((y1 - y0) / pixel_mm)) + 0.5)).astype(np.float64)
    grid = np.column_stack([g.ravel() for g in np.meshgrid(xs, ys)])
    shape = (ys.size, xs.size)
    radius = float(overlay.tool_radius.value)
    pocket = MplPath(pocket_polyline(overlay.pocket)).contains_points(grid).reshape(shape)
    omega = np.array([(p.x, p.y) for p in overlay.omega] + [(overlay.omega[0].x, overlay.omega[0].y)])
    target = MplPath(omega).contains_points(grid).reshape(shape)
    _paint_capsules(target, xs, ys, np.column_stack((omega[:-1], omega[1:])), radius)
    swept = np.zeros(shape, dtype=np.bool_)
    for circle in overlay.circles:
        cx, cy, big_r = circle.centre.x, circle.centre.y, circle.radius.value
        outer, inner = big_r + radius, max(big_r - radius, 0.0)
        i0, i1 = np.searchsorted(xs, [cx - outer - pixel_mm, cx + outer + pixel_mm])
        j0, j1 = np.searchsorted(ys, [cy - outer - pixel_mm, cy + outer + pixel_mm])
        if i0 >= i1 or j0 >= j1:
            continue
        gx, gy = np.meshgrid(xs[i0:i1] - cx, ys[j0:j1] - cy)
        d2 = gx * gx + gy * gy
        swept[j0:j1, i0:i1] |= (d2 <= outer * outer) & (d2 >= inner * inner)
    _paint_capsules(swept, xs, ys, _segments(overlay.connectors), radius)
    return MaterialRaster(extent, pixel_mm, pocket, target, swept)


def _bounds(overlay: OverlayInput) -> tuple[float, float, float, float]:
    poly = pocket_polyline(overlay.pocket)
    margin = 2.0 * float(overlay.tool_radius.value)
    return float(poly[:, 0].min() - margin), float(poly[:, 0].max() + margin), float(poly[:, 1].min() - margin), float(poly[:, 1].max() + margin)


def witness_clusters(points: Sequence[tuple[float, float]]) -> list[list[tuple[float, float]]]:
    """Single-linkage groups of witness points at WITNESS_CLUSTER_LINK_MM, largest first (reporting only)."""
    parent = list(range(len(points)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, p in enumerate(points):
        for j in range(i):
            if math.dist(p, points[j]) <= WITNESS_CLUSTER_LINK_MM:
                parent[root(i)] = root(j)
    groups: dict[int, list[tuple[float, float]]] = {}
    for i, p in enumerate(points):
        groups.setdefault(root(i), []).append(p)
    return sorted(groups.values(), key=lambda group: (-len(group), group[0]))


def choose_insets(overlay: OverlayInput, coverage: CoverageQualification, marks: Sequence[GougeMark]) -> tuple[Inset, ...]:
    """The exact-witness clusters with the most witnesses (ties: most uncut-in-T area at inset resolution), then the deepest gouge.

    A rendering choice, never a decision: every uncut inset frames exact
    witnesses, the ranking only orders them.
    """
    framed: list[tuple[int, int, tuple[float, float], float]] = []
    for group in witness_clusters(coverage.witnesses_mm):
        xs, ys = [p[0] for p in group], [p[1] for p in group]
        cx, cy = 0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys))
        half = max(INSET_MIN_HALF_MM, INSET_PAD_FACTOR * 0.5 * max(max(xs) - min(xs), max(ys) - min(ys)))
        local = material_raster(overlay, (cx - half, cx + half, cy - half, cy + half), 2 * half / INSET_PIXELS)
        framed.append((len(group), int(np.count_nonzero(local.uncut_in_target)), (cx, cy), half))
    framed.sort(key=lambda item: (-item[0], -item[1], item[2]))
    insets = [Inset(chr(ord("A") + rank), centre, half, "uncut") for rank, (_, _, centre, half) in enumerate(framed[:UNCUT_INSETS])]
    if marks:
        insets.append(Inset("G", marks[0].location_mm, GOUGE_INSET_HALF_MM, "gouge"))
    return tuple(insets)


# ---------------------------------------------------------------------------
# Caption and drawing
# ---------------------------------------------------------------------------


def caption(coverage: CoverageQualification, report: HeldGougeAuditReport) -> tuple[str, str]:
    """(EXACT line, DIAGNOSTIC line) of the figure footer."""
    geometric = sum(1 for item in report.gouges if item.depth_mm > REPORTING_SPLIT_MM)
    gouge = (
        f"gouge audit: {len(report.gouging_circles)} circles + {len(report.gouging_connectors)} connectors cut outside the analytic pocket "
        f"(dist < r decided exactly; {geometric} deeper than {REPORTING_SPLIT_MM:g} mm, max {max(report.max_depth_mm('circle'), report.max_depth_mm('connector')):.3f} mm)"
        if report.gouges
        else "gouge audit: no element cuts outside the analytic pocket (dist >= r decided exactly)"
    )
    regions = f"{coverage.regions} swept regions (every circle and connector, the seam included)"
    if coverage.status is CoverageStatus.CERTIFIED_COVERED:
        cover = f"coverage certificate: CERTIFIED COVERED, T = Ω ⊕ D_r lies inside the union of {regions}"
    else:
        undecided = f", {coverage.undecided_edge_count} exposed edge(s) undecided" if coverage.undecided_edge_count else ""
        noun = "witness" if len(coverage.witnesses) == 1 else "witnesses"
        cover = f"coverage certificate: NOT COVERED, {len(coverage.witnesses)} exact {noun} (x) in T outside all {regions}{undecided}"
    certified = f"EXACT - {cover}; {gouge}."
    diagnostic = (
        f"DIAGNOSTIC (raster {MAIN_PIXEL_MM:g} mm/px, insets finer; never a decision) - swept fill, red uncut regions, "
        "pale red = uncut outside T = Ω ⊕ D_r, gouge locations, inset framing, every printed area and depth."
    )
    return certified, diagnostic


def _draw_layers(ax: Axes, overlay: OverlayInput, raster: MaterialRaster) -> None:
    x0, x1, y0, y1 = raster.extent
    ax.imshow(raster.rgba(), origin="lower", extent=(x0, x1, y0, y1), interpolation="nearest", zorder=1)
    omega = np.array([(p.x, p.y) for p in overlay.omega] + [(overlay.omega[0].x, overlay.omega[0].y)])
    ax.plot(omega[:, 0], omega[:, 1], color=OMEGA, lw=0.5, ls=(0, (3, 2)), zorder=2)
    pocket = pocket_polyline(overlay.pocket)
    ax.plot(pocket[:, 0], pocket[:, 1], color=INK, lw=0.9, zorder=3)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_facecolor(SURFACE)


def _style_inset_axes(ax: Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color(MUTED)
        spine.set_linewidth(0.6)


def _scale_bar(ax: Axes, inset: Inset) -> None:
    length = max(step for step in SCALE_BAR_STEPS_MM if step <= SCALE_BAR_FRACTION * 2 * inset.half_width)
    cx, cy = inset.centre
    x = cx - 0.9 * inset.half_width
    y = cy - 0.9 * inset.half_width
    ax.plot([x, x + length], [y, y], color=INK, lw=1.4, solid_capstyle="butt", zorder=6)
    ax.text(
        x + 0.5 * length,
        y + 0.04 * inset.half_width,
        f"{length:g} mm",
        color=INK,
        fontsize=6.5,
        ha="center",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.1", "fc": SURFACE, "ec": "none"},
        zorder=6,
    )


def _in_window(point: tuple[float, float], inset: Inset) -> bool:
    return abs(point[0] - inset.centre[0]) < inset.half_width and abs(point[1] - inset.centre[1]) < inset.half_width


def _inset_note(inset: Inset, raster: MaterialRaster, coverage: CoverageQualification, marks: Sequence[GougeMark]) -> tuple[str, str]:
    """(heading, detail) under one inset; areas and depths are reported, never decided."""
    if inset.note == "uncut":
        count = sum(1 for point in coverage.witnesses_mm if _in_window(point, inset))
        area = float(raster.uncut_in_target.sum()) * raster.pixel_mm**2
        noun = "witness" if count == 1 else "witnesses"
        return f"{inset.label}  {count} exact {noun}", f"({inset.centre[0]:.2f}, {inset.centre[1]:.2f}), raster {area:.3f} mm² uncut"
    deepest = marks[0].verdict
    return f"{inset.label}  deepest gouge, {deepest.depth_mm:.3f} mm", f"{deepest.element} {deepest.index} against {deepest.primitive_kind} {deepest.primitive}"


def _labelled_gouges(marks: Sequence[GougeMark]) -> list[GougeMark]:
    chosen: list[GougeMark] = []
    for mark in marks:
        if len(chosen) == LABELLED_GOUGES or mark.verdict.depth_mm <= REPORTING_SPLIT_MM:
            break
        if all(math.dist(mark.location_mm, other.location_mm) >= GOUGE_LABEL_SEPARATION_MM for other in chosen):
            chosen.append(mark)
    return chosen


def _label_offset(location: tuple[float, float], extent: tuple[float, float, float, float]) -> tuple[float, float]:
    """Push a label towards the nearer outside edge of the plot, off the pocket interior."""
    x0, x1, y0, y1 = extent
    span = max(x1 - x0, y1 - y0)
    dx = -LABEL_OFFSET_X * span if location[0] - x0 < x1 - location[0] else LABEL_OFFSET_X * span
    dy = -LABEL_OFFSET_Y * span if location[1] - y0 < y1 - location[1] else LABEL_OFFSET_Y * span
    return dx, dy


def _footer_lines(text: str, width_in: float) -> list[str]:
    return textwrap.wrap(text, width=int(width_in * FOOTER_CHARS_PER_INCH))


def _legend_handles() -> list[Artist]:
    return [
        Patch(fc=SWEPT, ec="none", label="swept by the tool"),
        Patch(fc=UNCUT, ec="none", label="uncut, inside T"),
        Patch(fc=UNREACHABLE, ec="none", label="uncut, outside T"),
        Patch(fc=GOUGE, ec="none", label="gouged (outside pocket)"),
        Line2D([], [], color=INK, lw=0.9, label="analytic pocket"),
        Line2D([], [], color=OMEGA, lw=0.6, ls=(0, (3, 2)), label="centre domain Ω"),
        Line2D([], [], color=UNCUT, marker="x", lw=0, ms=6, mew=1.4, label="proven uncovered (exact)"),
        Line2D([], [], color=GOUGE, marker="o", lw=0, ms=4, label=f"gouge > {REPORTING_SPLIT_MM:g} mm (exact)"),
        Line2D([], [], color=GOUGE, marker="o", mfc="none", lw=0, ms=3, mew=0.6, label=f"gouge ≤ {REPORTING_SPLIT_MM:g} mm (exact)"),
    ]


def render(
    overlay: OverlayInput,
    coverage: CoverageQualification,
    report: HeldGougeAuditReport,
    output: Path,
    *,
    main_pixel_mm: float = MAIN_PIXEL_MM,
) -> tuple[Inset, ...]:
    """Draw the qualification overlay to ``output``; returns the insets it drew."""
    marks = gouge_marks(overlay, report)
    extent = _bounds(overlay)
    main = material_raster(overlay, extent, main_pixel_mm)
    insets = choose_insets(overlay, coverage, marks)

    # Layout in inches: [margin | main view | gap | inset column | margin], footer below.
    main_w = MAIN_HEIGHT_IN * (extent[1] - extent[0]) / (extent[3] - extent[2])
    fig_w = MARGIN_IN + main_w + COLUMN_GAP_IN + INSET_SIDE_IN + MARGIN_IN
    certified, diagnostic = caption(coverage, report)
    footer = _footer_lines(certified, fig_w - 2 * MARGIN_IN), _footer_lines(f"{diagnostic} Regenerate: pixi run held-qualification-figures", fig_w - 2 * MARGIN_IN)
    footer_h = FOOTER_LINE_IN * (len(footer[0]) + len(footer[1])) + FOOTER_GAP_IN
    bottom = MARGIN_IN + footer_h + LEGEND_IN
    fig_h = bottom + MAIN_HEIGHT_IN + TITLE_IN
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=SURFACE)

    ax = fig.add_axes((MARGIN_IN / fig_w, bottom / fig_h, main_w / fig_w, MAIN_HEIGHT_IN / fig_h))
    _draw_layers(ax, overlay, main)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    for mark in marks:
        big = mark.verdict.depth_mm > REPORTING_SPLIT_MM
        ax.plot(*mark.location_mm, marker="o", ms=3.0 if big else 1.6, mfc=GOUGE if big else "none", mec=GOUGE, mew=0.5, zorder=4)
    for mark in _labelled_gouges(marks):
        dx, dy = _label_offset(mark.location_mm, extent)
        ax.annotate(
            f"{mark.verdict.depth_mm:.3f} mm",
            xy=mark.location_mm,
            xytext=(mark.location_mm[0] + dx, mark.location_mm[1] + dy),
            color=GOUGE,
            fontsize=7.5,
            fontweight="bold",
            ha="center",
            va="center",
            arrowprops={"arrowstyle": "-", "color": GOUGE, "lw": 0.6, "shrinkB": 2},
            bbox={"boxstyle": "round,pad=0.2", "fc": SURFACE, "ec": GOUGE, "lw": 0.5},
            zorder=7,
        )
    for x, y in coverage.witnesses_mm:
        ax.plot(x, y, marker="x", ms=5.5, mew=1.4, color=UNCUT, zorder=5)

    column_x = (MARGIN_IN + main_w + COLUMN_GAP_IN) / fig_w
    # The inset column never runs below the main view (the footer is there).
    side_in = min(INSET_SIDE_IN, MAIN_HEIGHT_IN / max(len(insets), 1) - INSET_NOTE_IN)
    for k, inset in enumerate(insets):
        cx, cy = inset.centre
        half = inset.half_width
        ax.add_patch(Rectangle((cx - half, cy - half), 2 * half, 2 * half, fill=False, ec=INK, lw=0.8, zorder=6))
        ax.annotate(
            inset.label, xy=(cx - half, cy + half), xytext=(-3, 3), textcoords="offset points", color=INK, fontsize=8, fontweight="bold", ha="right", va="bottom", zorder=7
        )
        top_in = bottom + MAIN_HEIGHT_IN - k * (side_in + INSET_NOTE_IN)
        sub = fig.add_axes((column_x, (top_in - side_in) / fig_h, side_in / fig_w, side_in / fig_h))
        local = material_raster(overlay, (cx - half, cx + half, cy - half, cy + half), 2 * half / INSET_PIXELS)
        _draw_layers(sub, overlay, local)
        for point in coverage.witnesses_mm:
            if _in_window(point, inset):
                sub.plot(*point, marker="x", ms=7, mew=1.6, color=INK, zorder=5)
        for mark in marks:
            if _in_window(mark.location_mm, inset):
                sub.plot(*mark.location_mm, marker="o", ms=4, mfc="none", mec=GOUGE, mew=0.9, zorder=5)
        _scale_bar(sub, inset)
        _style_inset_axes(sub)
        heading, detail = _inset_note(inset, local, coverage, marks)
        note_y = (top_in - side_in - INSET_NOTE_PAD_IN) / fig_h
        fig.text(column_x, note_y, heading, fontsize=8, fontweight="bold", color=INK, ha="left", va="top")
        fig.text(column_x, note_y - FOOTER_LINE_IN * 1.15 / fig_h, detail, fontsize=7.5, color=MUTED, ha="left", va="top")

    fig.legend(
        handles=_legend_handles(),
        loc="lower left",
        bbox_to_anchor=(MARGIN_IN / fig_w, (MARGIN_IN + footer_h) / fig_h),
        ncol=5,
        frameon=False,
        fontsize=7.5,
        labelcolor=INK,
        handlelength=1.4,
        columnspacing=1.2,
        borderaxespad=0.0,
    )
    fig.text(
        MARGIN_IN / fig_w,
        1 - MARGIN_IN / fig_h,
        f"{overlay.title} — qualification overlay (swept material, not tool centres)",
        fontsize=TITLE_FONT_SIZE,
        color=INK,
        ha="left",
        va="top",
    )
    y = MARGIN_IN + footer_h - FOOTER_GAP_IN
    for colour, lines in ((INK, footer[0]), (MUTED, footer[1])):
        for line in lines:
            fig.text(MARGIN_IN / fig_w, y / fig_h, line, fontsize=FOOTER_FONT_SIZE, color=colour, ha="left", va="top")
            y -= FOOTER_LINE_IN
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=FIGURE_DPI, facecolor=SURFACE)
    plt.close(fig)
    return insets


@dataclass(frozen=True)
class OverlayWorkload:
    """A published workload and how its overlay is named."""

    workload: str
    name: str
    title: str


WORKLOADS: dict[str, OverlayWorkload] = {
    "figure5": OverlayWorkload("figure5_refined", "figure5", "Held Figure 5"),
    "monstera": OverlayWorkload("figure8_monstera_boundary_qualified_corners", "monstera", "Held Figure 8 Monstera"),
}


def output_path(name: str) -> Path:
    """Where the overlay of one workload is published."""
    return IMAGES / f"held_{name}_qualification.png"


def qualify_and_render(overlay: OverlayInput, output: Path, *, main_pixel_mm: float = MAIN_PIXEL_MM) -> tuple[CoverageQualification, HeldGougeAuditReport, tuple[Inset, ...]]:
    """Decide (certificate witnesses, gouge audit), then draw."""
    coverage = qualify_coverage(overlay)
    report = gouge_audit(overlay)
    insets = render(overlay, coverage, report, output, main_pixel_mm=main_pixel_mm)
    return coverage, report, insets


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("names", nargs="*", default=list(WORKLOADS), help=f"overlays to render (default: {' '.join(WORKLOADS)})")
    args = parser.parse_args(argv)
    for name in args.names:
        if name not in WORKLOADS:
            raise UnknownOverlayWorkloadError(f"unknown overlay {name!r}; known: {list(WORKLOADS)}")
        spec = WORKLOADS[name]
        overlay = OverlayInput.from_generated(generate(spec.workload), spec.name, spec.title)
        coverage, report, insets = qualify_and_render(overlay, output_path(spec.name))
        print(f"{spec.name}: {coverage.status.value}, {len(coverage.witnesses)} exact witness(es) over {coverage.regions} regions; {len(report.gouges)} gouging elements")
        for witness in coverage.witnesses:
            print(f"  witness {witness.piece_kind.name} region {witness.region_index} at {witness.uncovered_point_mm}")
        for inset in insets:
            print(f"  inset {inset.label} ({inset.note}) at ({inset.centre[0]:.2f}, {inset.centre[1]:.2f}), half-width {inset.half_width:.2f} mm")
        print(f"  wrote {output_path(spec.name).relative_to(ROOT)}")


if __name__ == "__main__":
    main()
