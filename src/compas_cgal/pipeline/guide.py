"""The guide: ordered straight-skeleton stations, read once from the C++ generator and typed.

`compas_cgal.toolpath.trochoidal_mat_toolpath_circular` already extracts the
interior straight skeleton, splits it into chains, orders them, and derives each
machining-circle radius from an exact clearance query at its own centre. Its
circle rows, grouped by runs of ``path_index``, are exactly the ordered station
sequence the pipeline walks. This module reads them through the native call
directly, because the wrapper only WARNS when ``max_passes`` truncates the guide,
and a truncated guide is a silent under-cut: here it raises `GuideTruncatedError`.

The read reproduces `compas_cgal.engagement_toolpath._guide_chains` bit for bit
(same native arguments, same tangent estimate): held against rho's public loops
by `tests/pipeline/test_against_rho.py`, and bit for bit against the private
read by a logged scratch probe, not in the repository (Ruling 1.1-B-oracle). This module imports no frozen generator.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Self

import numpy as np
from compas.geometry import Polygon

from compas_cgal import _toolpath  # type: ignore[attr-defined]
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import EmptyPipelineGuideError
from compas_cgal.pipeline.errors import GuideStationOrderError
from compas_cgal.pipeline.errors import GuideTruncatedError
from compas_cgal.toolpath import _holes_to_vertex_arrays
from compas_cgal.toolpath import _polygon_to_ccw_vertices
from compas_cgal.units import Clearance
from compas_cgal.units import ClearanceZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import Spacing
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# Guide station spacing in tool diameters: one station is D/40 = r/20 of radial
# depth. Same derivation and value as the advance-regulated generators' guide
# (engagement_toolpath.GUIDE_STEP_TOOL_DIAMETERS, held equal by test_against_rho):
# at half-radius immersion dTEA/dae ~ 2.31/r, so one station is worth at most
# ~0.12 rad (~7 deg) of engagement, the cap headroom the quantisation leaves.
GUIDE_STEP_TOOL_DIAMETERS = 0.025

# The guide's own stepover is this multiple of its pitch, so the pitch governs the
# station spacing and the grid is uniform. Two suffices: the clearance along a
# chain is 1-Lipschitz, so advance + radius growth per station <= 2 * pitch.
GUIDE_STEPOVER_PITCH_MULTIPLE = 2.0

# Native guide arguments that switch off everything but the stations: no lower or
# upper radius clamp, unscaled clearance, no leads, no links, no final retract.
GUIDE_MIN_TROCHOID_RADIUS_MM = 0.0
GUIDE_MAX_TROCHOID_RADIUS_MM = math.inf
GUIDE_MAT_SCALE = 1.0
GUIDE_LEAD_MM = 0.0
GUIDE_CUT_Z_MM = 0.0

# Polyline densities of the discarded native tessellation: the wrapper's default
# samples per arc primitive, and the coarsest samples-per-radian the API takes.
GUIDE_SAMPLES_PER_CYCLE = 20
GUIDE_SAMPLES_PER_RADIAN = 1.0

# Native primitive-type code of a full machining circle in the meta matrix
# (0 line, 1 arc, 2 circle; `compas_cgal.toolpath._matrices_to_operations`).
NATIVE_CIRCLE_PRIMITIVE = 2

# A meta row's turn flag is stored as 0.0 / 1.0; any value above one half reads
# clockwise, the wrapper's own decoding.
NATIVE_FLAG_THRESHOLD = 0.5

# Tangent of a single-station chain, where the central difference has no
# direction: +X, deterministic. The tangent only orients the entry around a
# circle that is traversed in full, so it moves where the loop starts, never
# what it removes.
SINGLE_STATION_TANGENT = (1.0, 0.0)


def pipeline_spacing(diameter: float, cap: Radian) -> Spacing:
    """Guide station spacing for a tool of *diameter* at engagement cap *cap*.

    In 1.1 the spacing is `GUIDE_STEP_TOOL_DIAMETERS` tool diameters at every
    cap; the cap is in the signature so that a regulated spacing (Task 1.5) has
    one named place to land.

    Args:
        diameter: Tool diameter, mm.
        cap: Engagement cap (unused by the 1.1 rule).

    Returns:
        The station spacing.
    """
    del cap
    return Spacing.build(GUIDE_STEP_TOOL_DIAMETERS * diameter)


@dataclass(frozen=True)
class GuideStation:
    """One ordered station on a skeleton chain.

    Attributes:
        centre: Machining-circle centre, on the guide.
        radius: Gouge-free machining-circle radius the C++ guide derived from
            the exact clearance at *centre*, less tool radius and clearance.
        tangent: Unit guide tangent in the walk direction (a central-difference
            estimate; it orients the entry, it decides nothing).
        clockwise: Turn direction of the machining circle.
    """

    centre: Point2[WorldXY]
    radius: GuideRadius
    tangent: Direction2[WorldXY]
    clockwise: bool


@dataclass(frozen=True)
class GuideChain:
    """One skeleton chain in one walk direction.

    Attributes:
        stations: Stations in walk order.
        origin: Index of the chain in the guide's emission order, shared by both
            walk directions.
        reversed_walk: Whether the stations run against the guide's order.
    """

    stations: tuple[GuideStation, ...]
    origin: int
    reversed_walk: bool = False

    def reversed(self) -> GuideChain:
        """The same chain walked the other way.

        The tangent is a direction and reverses with the walk; radius, centre
        and turn direction belong to the station and do not.

        Returns:
            The reversed chain; reversing twice returns an equal chain.
        """
        stations = tuple(
            GuideStation(
                centre=station.centre,
                radius=station.radius,
                tangent=station.tangent.reversed(),
                clockwise=station.clockwise,
            )
            for station in reversed(self.stations)
        )
        return GuideChain(stations=stations, origin=self.origin, reversed_walk=not self.reversed_walk)


@dataclass(frozen=True)
class Guide:
    """The ordered skeleton chains of one pocket, in the C++ guide's emission order.

    Attributes:
        chains: One chain per skeleton chain, forward walk direction.
    """

    chains: tuple[GuideChain, ...]

    @classmethod
    def build(
        cls,
        polygon: Polygon,
        holes: list[Polygon],
        tool: ToolRadius,
        step: Spacing,
        radial_clearance: Clearance,
        climb: bool,
        max_passes: int,
        clearance_z: ClearanceZ,
    ) -> Self:
        """Read the C++ guide's stations for one pocket.

        Args:
            polygon: Outer pocket boundary in world XY.
            holes: Island boundaries strictly inside *polygon*.
            tool: Tool radius.
            step: Station spacing along the guide.
            radial_clearance: Safety clearance subtracted from each radius.
            climb: True for climb milling (clockwise loops).
            max_passes: Most skeleton chains the guide may emit.
            clearance_z: Rapid plane handed to the native guide. It emits no
                links, so the plane moves no station.

        Returns:
            The guide.

        Raises:
            GuideTruncatedError: *max_passes* left skeleton chains unread.
            EmptyPipelineGuideError: The guide emitted no machining station.
            GuideStationOrderError: A chain's station run was interrupted.
        """
        diameter = 2.0 * tool.value
        meta, _starts, _ends, centres, radii, _polyline, _start_tangents, _end_tangents, skipped = _toolpath.trochoidal_mat_toolpath_circular(
            _polygon_to_ccw_vertices(polygon),
            _holes_to_vertex_arrays(holes),
            float(diameter),
            float(GUIDE_STEPOVER_PITCH_MULTIPLE * step.value),
            float(step.value),
            GUIDE_MIN_TROCHOID_RADIUS_MM,
            GUIDE_MAX_TROCHOID_RADIUS_MM,
            GUIDE_MAT_SCALE,
            float(radial_clearance.value),
            GUIDE_SAMPLES_PER_CYCLE,
            int(max_passes),
            bool(climb),
            GUIDE_LEAD_MM,
            GUIDE_LEAD_MM,
            False,
            True,
            GUIDE_CUT_Z_MM,
            float(clearance_z.value),
            True,
            False,
            GUIDE_SAMPLES_PER_RADIAN,
        )
        if int(skipped) > 0:
            raise GuideTruncatedError(f"max_passes={max_passes} left {int(skipped)} skeleton chain(s) unread; the guide does not cover the pocket.")
        rows = _circle_rows(np.asarray(meta, dtype=np.float64), np.asarray(centres, dtype=np.float64), np.asarray(radii, dtype=np.float64).reshape(-1))
        if not rows:
            raise EmptyPipelineGuideError(f"The straight-skeleton guide produced no machining station for tool radius {tool.value!r}; the pocket admits no gouge-free loop.")
        return cls(chains=tuple(_chain(run, origin) for origin, run in enumerate(rows)))


_Row = tuple[float, float, float, bool]


def _circle_rows(meta: np.ndarray, centres: np.ndarray, radii: np.ndarray) -> list[list[_Row]]:
    """Circle rows grouped by runs of equal ``path_index``, in emission order.

    Args:
        meta: Native meta matrix (path index, primitive type, turn flag, op code).
        centres: Native centre rows.
        radii: Native radii.

    Returns:
        One run of ``(cx, cy, radius, clockwise)`` per chain.

    Raises:
        GuideStationOrderError: A path index reappears after its run ended.
    """
    runs: list[list[_Row]] = []
    seen: set[int] = set()
    current: int | None = None
    for row in range(meta.shape[0]):
        if int(meta[row, 1]) != NATIVE_CIRCLE_PRIMITIVE:
            continue
        path_index = int(meta[row, 0])
        if path_index != current:
            if path_index in seen:
                raise GuideStationOrderError(f"guide chain {path_index} resumes after another chain's stations.")
            seen.add(path_index)
            runs.append([])
            current = path_index
        runs[-1].append((float(centres[row, 0]), float(centres[row, 1]), float(radii[row]), bool(meta[row, 2] > NATIVE_FLAG_THRESHOLD)))
    return runs


def _unit(ax: float, ay: float, bx: float, by: float) -> tuple[float, float]:
    """Unit vector from ``a`` to ``b``; the zero vector when they coincide (construction guard)."""
    dx, dy = bx - ax, by - ay
    length = math.hypot(dx, dy)
    if length == 0.0:
        return 0.0, 0.0
    return dx / length, dy / length


def _chain(run: list[_Row], origin: int) -> GuideChain:
    """Type one run of circle rows as a forward chain.

    The tangent is the central difference, one-sided at the chain ends: a local
    estimate that only orients the entry.

    Args:
        run: The chain's ``(cx, cy, radius, clockwise)`` rows in walk order.
        origin: The chain's emission index.

    Returns:
        The chain.
    """
    count = len(run)
    stations: list[GuideStation] = []
    for index, (cx, cy, radius, clockwise) in enumerate(run):
        before = run[max(0, index - 1)]
        after = run[min(count - 1, index + 1)]
        tx, ty = _unit(before[0], before[1], after[0], after[1])
        if (tx, ty) == (0.0, 0.0):
            tx, ty = SINGLE_STATION_TANGENT
        stations.append(
            GuideStation(
                centre=Point2[WorldXY](Millimetre(cx), Millimetre(cy)),
                radius=GuideRadius(Millimetre(radius)),
                tangent=Direction2[WorldXY].build(tx, ty),
                clockwise=clockwise,
            )
        )
    return GuideChain(stations=tuple(stations), origin=origin)
