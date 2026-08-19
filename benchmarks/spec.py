from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field
from typing import Mapping

from compas.geometry import Polygon

from benchmarks.errors import DegeneratePocketError
from benchmarks.errors import InvalidCapError
from benchmarks.errors import NonPositiveToolError
from benchmarks.errors import PocketNotSimpleError

# The exact kernel's documented cap contract is theta in (0, pi]: one engaged run
# subtends at most a half turn before ">pi" is an exact orientation verdict, so a
# cap above 180 degrees is meaningless (src/engagement_2.cpp::certify_segment_tea).
MAX_CAP_DEG = 180.0

# A pocket must admit at least one full tool circle; below this multiple of the
# tool radius the pocket is a slot the generator cannot enter, which is a corpus
# authoring error rather than an interesting instance.
MIN_AREA_TOOL_RADII_SQ = 4.0


@dataclass(frozen=True)
class PocketSpec:
    """One reference machining problem: a pocket, a tool, and an engagement cap.

    Attributes:
        name: Unique instance name; becomes the row key in every report.
        family: Sweep family that produced this instance (e.g. ``"necks"``).
        polygon: Outer pocket boundary in the world XY plane, CCW.
        holes: Island boundaries, possibly empty.
        tool_diameter: Cutter diameter in the same units as the polygon.
        tea_cap_deg: Maximum tool-engagement angle in degrees.
        params: The sweep coordinates that produced this instance, for plotting.
    """

    name: str
    family: str
    polygon: Polygon
    holes: tuple[Polygon, ...] = ()
    tool_diameter: float = 1.0
    tea_cap_deg: float = 120.0
    params: Mapping[str, float] = field(default_factory=dict)

    @property
    def tool_radius(self) -> float:
        """Half the tool diameter."""
        return 0.5 * self.tool_diameter

    @property
    def tea_cap_rad(self) -> float:
        """The engagement cap in radians, as the kernel consumes it."""
        return math.radians(self.tea_cap_deg)

    @classmethod
    def build(
        cls,
        name: str,
        family: str,
        polygon: Polygon,
        tool_diameter: float,
        tea_cap_deg: float,
        holes: tuple[Polygon, ...] = (),
        params: Mapping[str, float] | None = None,
    ) -> "PocketSpec":
        """Validate every invariant and return a frozen spec.

        Args:
            name: Unique instance name.
            family: Sweep family name.
            polygon: Outer boundary.
            tool_diameter: Cutter diameter; must be finite and positive.
            tea_cap_deg: Engagement cap; must lie in (0, 180].
            holes: Island boundaries.
            params: Sweep coordinates recorded for plotting.

        Returns:
            The validated spec.

        Raises:
            NonPositiveToolError: The diameter is NaN, zero, or negative.
            InvalidCapError: The cap is NaN or outside (0, 180].
            PocketNotSimpleError: The boundary or a hole self-intersects.
            DegeneratePocketError: The pocket has fewer than three vertices or
                too little area to admit the tool.
        """
        if not math.isfinite(tool_diameter) or tool_diameter <= 0.0:
            raise NonPositiveToolError(f"tool_diameter must be finite and positive, got {tool_diameter!r}.")
        if not math.isfinite(tea_cap_deg) or not (0.0 < tea_cap_deg <= MAX_CAP_DEG):
            raise InvalidCapError(f"tea_cap_deg must lie in (0, {MAX_CAP_DEG}], got {tea_cap_deg!r}.")
        for ring in (polygon,) + tuple(holes):
            if len(ring.points) < 3:
                raise DegeneratePocketError(f"{name}: a ring has {len(ring.points)} vertices; at least three are required.")
            if not ring.is_convex and _self_intersects(ring):
                raise PocketNotSimpleError(f"{name}: a ring self-intersects and is not a valid general polygon.")
        radius = 0.5 * tool_diameter
        if abs(polygon.area) < MIN_AREA_TOOL_RADII_SQ * radius * radius:
            raise DegeneratePocketError(f"{name}: area {polygon.area:.6g} cannot admit a tool of radius {radius:.6g}.")
        return cls(
            name=name,
            family=family,
            polygon=polygon,
            holes=tuple(holes),
            tool_diameter=tool_diameter,
            tea_cap_deg=tea_cap_deg,
            params=dict(params or {}),
        )


def _self_intersects(ring: Polygon) -> bool:
    """Report whether any two non-adjacent edges of *ring* cross.

    Args:
        ring: The closed boundary to test.

    Returns:
        True when a crossing exists.
    """
    from compas.geometry import intersection_segment_segment_xy

    points = [(p[0], p[1], 0.0) for p in ring.points]
    n = len(points)
    edges = [(points[i], points[(i + 1) % n]) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if j == i or (j + 1) % n == i or (i + 1) % n == j:
                continue
            if intersection_segment_segment_xy(edges[i], edges[j]) is not None:
                return True
    return False
