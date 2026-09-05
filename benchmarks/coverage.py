"""What the path left behind, measured on a grid of exact membership queries.

Coverage is the one family of quality questions the exact kernel cannot answer
outright: `ExactRegion2` decides membership exactly but exposes no area, so a
FRACTION has to be counted rather than integrated. Every individual query here is
exact -- an exact point location in an exact reachable region and in the exact
depleted stock -- and only the aggregation is sampled. A residue thinner than one
cell is invisible, so every number this module returns is a LOWER BOUND on
residue and never a certificate of coverage.

THE DENOMINATOR IS THE REACHABLE MATERIAL, never the pocket. A sharp corner a
round tool can never enter is not residue; counting it would put a floor under
every path and bury the residue that is real.

`CoarseCoverageGridError` is the discipline that keeps the bound meaningful: a
grid too coarse to resolve residue would silently report LESS residue, so a
coarse grid is refused rather than answered.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List
from typing import Sequence
from typing import Tuple

from benchmarks.errors import CoarseCoverageGridError
from benchmarks.errors import EmptyReachableRegionError
from benchmarks.errors import InvalidGridResolutionError
from benchmarks.spec import PocketSpec
from compas_cgal import _coverage_2
from compas_cgal.stock import Stock
from compas_cgal.stock import _polygon_to_ccw_vertices

# Grid samples along the LONGER side of the pocket's bounding box; the shorter
# side is scaled to keep cells square, so every sample carries the same area and
# the count is an unweighted area estimate. A fixed count rather than a fixed
# cell size bounds the cost on a large pocket, and `MAX_CELL_TOOL_RADIUS_FRACTION`
# is what stops that bound from silently coarsening the measurement.
COVERAGE_GRID_SAMPLES = 200

# The coarsest cell the coverage measurement accepts, as a fraction of the tool
# radius. At a tenth of the tool radius at least two samples land across any
# residue a fifth of a tool radius wide, which is the smallest leftover a
# finishing pass at this tool size would still have to remove. ENGINEERING
# JUDGEMENT on the feature size; the CONSEQUENCE is not a judgement -- a coarser
# grid under-reports, so it raises instead of answering.
MAX_CELL_TOOL_RADIUS_FRACTION = 0.1

# How close to the wall a residue sample must be, in tool diameters, to count
# toward the wall scallop. One tool diameter is the widest band a single pass
# against the wall can leave, so residue further in than that belongs to the
# interior of the pocket and is reported by `uncut_fraction` instead.
WALL_BAND_TOOL_DIAMETERS = 1.0


def minimum_coverage_grid(spec: PocketSpec) -> int:
    """Return the smallest longer-axis count accepted for ``spec``."""
    maximum_cell = MAX_CELL_TOOL_RADIUS_FRACTION * spec.tool_radius
    width, height = _bounding_box_size(spec)
    grid = max(1, math.ceil(max(width, height) / maximum_cell))
    while _grid_shape(width, height, grid)[2] > maximum_cell:
        grid += 1
    return grid


@dataclass(frozen=True)
class CoverageEstimate:
    """What one grid pass found.

    Attributes:
        nx: Samples along x.
        ny: Samples along y.
        cell_area: Area each sample stands for.
        reachable_samples: Samples inside the tool-reachable material.
        uncut_reachable_samples: Reachable samples the path did not remove.
        remaining_samples: Samples still inside the stock anywhere in the
            pocket, reachable or not.
        wall_scallop_height: Deepest uncut reachable sample within
            `WALL_BAND_TOOL_DIAMETERS` of the boundary, measured perpendicular
            from the boundary. A residue DEPTH near the wall, not the classical
            two-pass cusp formula.
    """

    nx: int
    ny: int
    cell_area: float
    reachable_samples: int
    uncut_reachable_samples: int
    remaining_samples: int
    wall_scallop_height: float

    @property
    def uncut_fraction(self) -> float:
        """Uncut reachable samples over reachable samples, in ``[0, 1]``."""
        return self.uncut_reachable_samples / self.reachable_samples

    @property
    def remaining_area(self) -> float:
        """Estimated area of material still standing anywhere in the pocket."""
        return self.remaining_samples * self.cell_area


def measure_coverage(spec: PocketSpec, stock: Stock, *, grid: int = COVERAGE_GRID_SAMPLES) -> CoverageEstimate:
    """Count what *stock* still holds, against what the tool could have reached.

    Args:
        spec: The pocket and tool the reachable region is derived from.
        stock: The depleted stock the path left.
        grid: Samples along the pocket's longer bounding-box side.

    Returns:
        The counts and the wall residue depth.

    Raises:
        InvalidGridResolutionError: *grid* is below one.
        CoarseCoverageGridError: The resulting cell is coarser than
            `MAX_CELL_TOOL_RADIUS_FRACTION` of the tool radius.
        EmptyReachableRegionError: No sample landed in the reachable region.
        ReachableDomainConstructionError: The kernel could not build the
            reachable domain for this pocket and tool.
    """
    xs, ys, cell = _grid_axes(spec, grid)
    region = _coverage_2.ReachableMaterialPredicate2.build(
        _polygon_to_ccw_vertices(spec.polygon),
        [_polygon_to_ccw_vertices(hole) for hole in spec.holes],
        spec.tool_radius,
    )
    edges = _boundary_edges(spec)
    wall_band = WALL_BAND_TOOL_DIAMETERS * spec.tool_diameter

    reachable = 0
    uncut = 0
    remaining = 0
    scallop = 0.0
    for x in xs:
        for y in ys:
            # The depleted stock is the more expensive arrangement, but
            # `remaining_samples` needs it everywhere in the pocket, so it is
            # queried first and the reachable test narrows the numerator.
            in_stock = stock.contains(x, y)
            if in_stock:
                remaining += 1
            if not region.contains(x, y):
                continue
            reachable += 1
            if not in_stock:
                continue
            uncut += 1
            depth = _distance_to_boundary(x, y, edges)
            if depth <= wall_band and depth > scallop:
                scallop = depth
    if reachable == 0:
        raise EmptyReachableRegionError(f"{spec.name}: no sample of a {len(xs)}x{len(ys)} grid landed in the tool-reachable region, so no coverage fraction is defined.")
    return CoverageEstimate(
        nx=len(xs),
        ny=len(ys),
        cell_area=cell * cell,
        reachable_samples=reachable,
        uncut_reachable_samples=uncut,
        remaining_samples=remaining,
        wall_scallop_height=scallop,
    )


def _grid_axes(spec: PocketSpec, grid: int) -> Tuple[List[float], List[float], float]:
    """Cell-centre coordinates of the coverage grid over *spec*'s bounding box.

    Cell CENTRES, never the box edges: a sample sitting exactly on an
    axis-aligned wall is an on-boundary query the exact predicate answers
    deterministically and identically for a whole row of samples, which would
    bias the count by one arbitrary convention. Half a cell in from the edge, no
    sample of an axis-aligned pocket lands on its boundary at all.

    Args:
        spec: The pocket whose outer boundary sets the box and whose tool radius
            sets the coarsest acceptable cell.
        grid: Samples along the longer side.

    Returns:
        ``(xs, ys, cell)`` -- the sample coordinates and the cell size.

    Raises:
        InvalidGridResolutionError: *grid* is below one.
        CoarseCoverageGridError: The cell is coarser than
            `MAX_CELL_TOOL_RADIUS_FRACTION` of the tool radius.
    """
    if grid < 1:
        raise InvalidGridResolutionError(f"grid must be at least 1 sample per axis, got {grid!r}.")
    x_min, y_min, width, height = _bounding_box(spec)
    nx, ny, cell = _grid_shape(width, height, grid)
    coarsest = MAX_CELL_TOOL_RADIUS_FRACTION * spec.tool_radius
    if cell > coarsest:
        raise CoarseCoverageGridError(
            f"{spec.name}: a {nx}x{ny} grid over {width:g}x{height:g} gives a {cell:.4g} cell, coarser than the {coarsest:.4g} "
            f"({MAX_CELL_TOOL_RADIUS_FRACTION:g} of the {spec.tool_radius:g} tool radius) at which residue is still resolved; raise `grid`."
        )
    return (
        [x_min + (index + 0.5) * width / nx for index in range(nx)],
        [y_min + (index + 0.5) * height / ny for index in range(ny)],
        cell,
    )


def _bounding_box(spec: PocketSpec) -> tuple[float, float, float, float]:
    xs = [float(point[0]) for point in spec.polygon.points]
    ys = [float(point[1]) for point in spec.polygon.points]
    x_min, y_min = min(xs), min(ys)
    return x_min, y_min, max(xs) - x_min, max(ys) - y_min


def _bounding_box_size(spec: PocketSpec) -> tuple[float, float]:
    _, _, width, height = _bounding_box(spec)
    return width, height


def _grid_shape(width: float, height: float, grid: int) -> tuple[int, int, float]:
    if width >= height:
        nx, ny = grid, max(1, round(grid * height / width))
    else:
        nx, ny = max(1, round(grid * width / height)), grid
    return nx, ny, max(width / nx, height / ny)


def _boundary_edges(spec: PocketSpec) -> Tuple[Tuple[float, float, float, float], ...]:
    """Every wall segment of the pocket, outer boundary and islands alike.

    Args:
        spec: The instance.

    Returns:
        One ``(x0, y0, x1, y1)`` per edge.
    """
    edges: List[Tuple[float, float, float, float]] = []
    for ring in (spec.polygon,) + tuple(spec.holes):
        points = [(float(p[0]), float(p[1])) for p in ring.points]
        for index, (x0, y0) in enumerate(points):
            x1, y1 = points[(index + 1) % len(points)]
            edges.append((x0, y0, x1, y1))
    return tuple(edges)


def _distance_to_boundary(x: float, y: float, edges: Sequence[Tuple[float, float, float, float]]) -> float:
    """Perpendicular distance from ``(x, y)`` to the nearest wall.

    Args:
        x: Sample x.
        y: Sample y.
        edges: Wall segments from `_boundary_edges`.

    Returns:
        The distance; ``inf`` when the pocket has no edges, which
        `PocketSpec.build` already rules out.
    """
    best = math.inf
    for x0, y0, x1, y1 in edges:
        dx, dy = x1 - x0, y1 - y0
        span = dx * dx + dy * dy
        if span == 0.0:
            best = min(best, math.hypot(x - x0, y - y0))
            continue
        t = max(0.0, min(1.0, ((x - x0) * dx + (y - y0) * dy) / span))
        best = min(best, math.hypot(x - (x0 + t * dx), y - (y0 + t * dy)))
    return best
