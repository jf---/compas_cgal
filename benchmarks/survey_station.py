"""Station-level primitives shared by the motion survey and the helix observation.

One probed cutter position (`EngagementSample`), the full-turn seam query that
reads a closed turn's hand-over station, the closed-form swept area of a
primitive, and the unit travel direction on a circle. `benchmarks.survey` and
`benchmarks.helix_observation` both import these; neither imports the other's
internals through here, so there is no cycle.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from benchmarks.errors import UnmeasurableOperationLengthError
from compas_cgal import _coverage_2
from compas_cgal import _stock_2
from compas_cgal.stock import Stock
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# The engagement query is asked WITHOUT gap-closure pessimism, matching
# `benchmarks.exceedance.NO_GAP_CLOSURE`: this is a question about the material
# actually engaged at a position, not about what a certifier must assume could
# happen before the next station.
NO_GAP_CLOSURE = 0.0


# The prefix parameter of a closed loop's own end: the full turn.
FULL_TURN_PREFIX = 1.0


@dataclass(frozen=True)
class EngagementSample:
    """One probed cutter position on a motion.

    Attributes:
        distance: Arc length from the motion's start to this position.
        position: Cutter-centre position in the world XY frame.
        engagement_deg: The largest engaged run of the cutter rim, in degrees.
            REPORTING ONLY: this value never decides cap exceedance.
        cap_exceeded: The exact cap predicate's verdict at this position.
        inside_centre_domain: Whether the cutter centre lies where a tool of this
            radius may legally be. False is a demonstrated gouge.

    Each retained ``cap_exceeded`` verdict is exact. Absence of exceedance across
    this finite station set remains sampled-negative evidence, not a certificate
    for the continuous motion.
    """

    distance: float
    position: Point2[WorldXY]
    engagement_deg: float
    cap_exceeded: bool
    inside_centre_domain: bool


def seam_station(
    stock: Stock,
    circle: Circle,
    seam: Tuple[float, float],
    clockwise: bool,
    distance: float,
    tool_radius: float,
    cap_ratio: float,
    centre_domain: "_coverage_2.CutterCentreDomain2",
) -> EngagementSample:
    """The seam of a closed turn, read after the whole turn has been swept.

    Args:
        stock: The stock the turn started from; read, never mutated.
        circle: The guide circle (centre and authored radius).
        seam: Where the turn starts and ends.
        clockwise: Turn sense in world XY.
        distance: Path length travelled to the seam, REPORTING.
        tool_radius: Tool radius.
        cap_ratio: Exact squared-chord surrogate of the cap.
        centre_domain: Where a cutter of this radius may legally be centred.

    Returns:
        The seam station with the full-turn prefix.
    """
    cx, cy = float(circle.frame.point[0]), float(circle.frame.point[1])
    sx, sy = seam
    ((_total, max_run, exceeded),) = _stock_2.engagement_circle_prefixes(
        stock.raw, cx, cy, float(circle.radius), sx, sy, clockwise, [FULL_TURN_PREFIX], tool_radius, cap_ratio, NO_GAP_CLOSURE
    )
    return EngagementSample(
        distance=distance,
        position=Point2[WorldXY].build(sx, sy),
        engagement_deg=math.degrees(max_run),
        cap_exceeded=exceeded,
        inside_centre_domain=centre_domain.contains(sx, sy),
    )


def swept_area(geometry: object, tool_radius: float) -> float:
    """Closed-form area swept by a tool of *tool_radius* along *geometry*.

    A segment of length L sweeps a capsule, ``2 r L + pi r^2``. A guide circle of
    radius rho sweeps an annulus of area ``4 pi rho r`` when ``rho > r``, and a
    filled disk of area ``pi (rho + r)^2`` when it does not -- which is the same
    boundary `benchmarks.quality.degenerate_loops` tests, arrived at from the
    area rather than from the hole. An arc of swept angle alpha sweeps the
    corresponding annular sector plus one full end disk.

    Args:
        geometry: The motion's primitive.
        tool_radius: Tool radius.

    Returns:
        The swept area.

    Raises:
        UnmeasurableOperationLengthError: The primitive is not a line, arc, or
            circle.
    """
    if isinstance(geometry, Circle):
        rho = float(geometry.radius)
        if rho > tool_radius:
            return 4.0 * math.pi * rho * tool_radius
        return math.pi * (rho + tool_radius) ** 2
    if isinstance(geometry, Arc):
        rho = float(geometry.radius)
        alpha = float(geometry.angle)
        inner = max(0.0, rho - tool_radius)
        return 0.5 * alpha * ((rho + tool_radius) ** 2 - inner**2) + math.pi * tool_radius**2
    if isinstance(geometry, Line):
        return 2.0 * tool_radius * float(geometry.length) + math.pi * tool_radius**2
    raise UnmeasurableOperationLengthError(f"Geometry {type(geometry).__name__!r} has no closed-form swept area.")


def perpendicular(dx: float, dy: float, turn: float) -> Tuple[float, float]:
    """The unit travel direction at a point on a circle, given its turn sense.

    Args:
        dx: X offset from the centre to the point.
        dy: Y offset from the centre to the point.
        turn: ``+1`` counterclockwise, ``-1`` clockwise.

    Returns:
        The unit tangent, or ``(0, 0)`` at the centre.
    """
    return unit(-turn * dy, turn * dx)


def unit(dx: float, dy: float) -> Tuple[float, float]:
    """Normalise a vector, returning ``(0, 0)`` when it has no length.

    A zero return is a fact about the primitive -- a zero-length move has no
    direction -- and `benchmarks.quality` counts those separately rather than
    treating them as a direction of (1, 0).

    Args:
        dx: X component.
        dy: Y component.

    Returns:
        The unit vector, or ``(0.0, 0.0)``.
    """
    norm = math.hypot(dx, dy)
    if norm == 0.0:
        return (0.0, 0.0)
    return (dx / norm, dy / norm)
