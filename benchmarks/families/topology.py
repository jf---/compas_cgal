"""Island-count sweep: boolean topology at fixed outer boundary.

The outer pocket and the island size are both fixed, so the only thing the sweep
varies is how many disconnected boundary components the exact stock arrangement
carries -- and, with them, how many separate channels the toolpath must thread.

The grid is refused rather than crowded. `Stock` rejects touching holes outright,
and a channel narrower than the tool is not a hard instance but an unmachinable
one, so `island_grid` fails loudly at the density where either would happen
instead of quietly producing an instance that is no longer the same experiment.

MEASURED (2026-08-21, tool 4.0, cap 120 deg, via `benchmarks.runner.run_spec`):

    islands                1       4       8
    generate (s)      0.0057  0.0134  0.0258
    certify (s)        12.21   24.15   33.49
    operations           316     703    1253
    uncertified          146     368     538
    arrangement verts   3521    3768    4372

Eight times the islands costs 4.5x the generation and 2.7x the certification, so
the axis moves what it claims to move. Note that certification grows SUBLINEARLY
in operations here (4.0x the operations, 2.7x the time) -- islands add channels
faster than they add per-operation difficulty.
"""

from __future__ import annotations

from compas.geometry import Polygon

from benchmarks.errors import CrowdedIslandGridError
from benchmarks.errors import InvalidIslandCountError
from benchmarks.spec import PocketSpec

# Outer pocket size. Fixed so the sweep varies hole count only.
OUTER_WIDTH = 40.0
OUTER_HEIGHT = 30.0

# Each island is a square of this side length.
ISLAND_SIDE = 3.0

# A channel between two islands must exceed the tool diameter by this factor
# before the instance counts as machinable rather than merely traversable, on the
# same reasoning as `analytic.MIN_CHANNEL_CLEARANCE_FACTOR`: at exactly the tool
# diameter the generator has zero working room and emits no operations there, so
# the island would silently stop contributing to the measurement.
MIN_ISLAND_GAP_FACTOR = 1.2


def island_grid(rows: int, cols: int, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A rectangular pocket containing a ``rows x cols`` grid of square islands.

    Args:
        rows: Island rows, at least one.
        cols: Island columns, at least one.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        InvalidIslandCountError: Either count is below one.
        CrowdedIslandGridError: The grid leaves a channel the tool cannot work in.
    """
    if rows < 1 or cols < 1:
        raise InvalidIslandCountError(f"island_grid: rows and cols must both be at least 1, got rows={rows}, cols={cols}.")
    pitch_x = OUTER_WIDTH / (cols + 1)
    pitch_y = OUTER_HEIGHT / (rows + 1)
    needed = MIN_ISLAND_GAP_FACTOR * tool_diameter
    gap = min(pitch_x, pitch_y) - ISLAND_SIDE
    if gap < needed:
        raise CrowdedIslandGridError(
            f"island_grid: a {rows}x{cols} grid of {ISLAND_SIDE:g}-sided islands leaves a {gap:g} channel, below the {needed:g} a tool of diameter {tool_diameter:g} can work in."
        )
    outer = Polygon([[0.0, 0.0, 0.0], [OUTER_WIDTH, 0.0, 0.0], [OUTER_WIDTH, OUTER_HEIGHT, 0.0], [0.0, OUTER_HEIGHT, 0.0]])
    holes: list[Polygon] = []
    half = 0.5 * ISLAND_SIDE
    for r in range(rows):
        for c in range(cols):
            cx = pitch_x * (c + 1)
            cy = pitch_y * (r + 1)
            # Islands are CW so they read as holes to the exact stock model.
            holes.append(Polygon([[cx - half, cy - half, 0.0], [cx - half, cy + half, 0.0], [cx + half, cy + half, 0.0], [cx + half, cy - half, 0.0]]))
    return PocketSpec.build(
        name=f"islands_{rows}x{cols}",
        family="topology",
        polygon=outer,
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        holes=tuple(holes),
        params={"rows": float(rows), "cols": float(cols), "islands": float(rows * cols)},
    )


def island_sweep(counts: tuple[tuple[int, int], ...], tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep the island count.

    Args:
        counts: ``(rows, cols)`` pairs to generate.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per pair, in input order.

    Raises:
        InvalidIslandCountError: A requested count is below one.
        CrowdedIslandGridError: A requested grid leaves an unmachinable channel.
    """
    return [island_grid(rows=r, cols=c, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for r, c in counts]
