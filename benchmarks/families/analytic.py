"""Pocket families whose medial axis and clearance field are known in closed form.

These instances are the corpus's correctness net: because the exact clearance at
any axis point is derivable, a certifier disagreeing with it has a bug, and the
disagreement is attributable without a second implementation.
"""

from __future__ import annotations

import math

from compas.geometry import Polygon

from benchmarks.errors import DegeneratePocketError
from benchmarks.spec import PocketSpec

# Segments used to tessellate each semicircular cap or full circle. 64 keeps the
# inscribed-polygon sagitta below 0.13% of the radius, well under the 2% slack the
# area assertions allow, while keeping arrangement sizes small enough to run fast.
ANALYTIC_SEGMENTS_PER_CAP = 64

# A channel must be wider than the tool by at least this multiple of the tool
# radius, or the generator has no room to place a machining circle at all.
MIN_CHANNEL_CLEARANCE_FACTOR = 1.2


def disk(radius: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A circular pocket: the medial axis is a single point.

    Args:
        radius: Pocket radius.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    n = 2 * ANALYTIC_SEGMENTS_PER_CAP
    points = [[radius * math.cos(2.0 * math.pi * i / n), radius * math.sin(2.0 * math.pi * i / n), 0.0] for i in range(n)]
    return PocketSpec.build(
        name=f"disk_r{radius:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"radius": radius},
    )


def rectangle(width: float, height: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """An axis-aligned rectangular pocket centred on the origin.

    Args:
        width: Extent along x.
        height: Extent along y.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The rectangle is narrower than the tool needs.
    """
    _require_channel(0.5 * height, tool_diameter, "rectangle")
    hw, hh = 0.5 * width, 0.5 * height
    points = [[-hw, -hh, 0.0], [hw, -hh, 0.0], [hw, hh, 0.0], [-hw, hh, 0.0]]
    return PocketSpec.build(
        name=f"rect_{width:g}x{height:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"width": width, "height": height},
    )


def stadium(straight_length: float, half_width: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A constant-width channel: rectangle of *straight_length* with semicircular caps.

    Its medial axis is the segment ``y = 0, |x| <= straight_length / 2`` and the
    clearance is exactly *half_width* everywhere on it — the conservation law the
    congruence tests in Task 3 exploit.

    Args:
        straight_length: Length of the straight section.
        half_width: Channel half-width, equal to the cap radius.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The channel is narrower than the tool needs.
    """
    _require_channel(half_width, tool_diameter, "stadium")
    hl = 0.5 * straight_length
    points: list[list[float]] = []
    for i in range(ANALYTIC_SEGMENTS_PER_CAP + 1):  # right cap, -90 deg -> +90 deg
        a = -0.5 * math.pi + math.pi * i / ANALYTIC_SEGMENTS_PER_CAP
        points.append([hl + half_width * math.cos(a), half_width * math.sin(a), 0.0])
    for i in range(ANALYTIC_SEGMENTS_PER_CAP + 1):  # left cap, +90 deg -> +270 deg
        a = 0.5 * math.pi + math.pi * i / ANALYTIC_SEGMENTS_PER_CAP
        points.append([-hl + half_width * math.cos(a), half_width * math.sin(a), 0.0])
    return PocketSpec.build(
        name=f"stadium_l{straight_length:g}_w{half_width:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"straight_length": straight_length, "half_width": half_width},
    )


def arc_channel(guide_radius: float, half_width: float, sweep_deg: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A constant-width channel whose spine is a circular arc about the origin.

    Args:
        guide_radius: Radius of the spine arc.
        half_width: Channel half-width.
        sweep_deg: Arc sweep in degrees, starting at the +x axis.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The channel is narrower than the tool needs.
    """
    _require_channel(half_width, tool_diameter, "arc_channel")
    sweep = math.radians(sweep_deg)
    steps = ANALYTIC_SEGMENTS_PER_CAP
    outer = [[(guide_radius + half_width) * math.cos(sweep * i / steps), (guide_radius + half_width) * math.sin(sweep * i / steps), 0.0] for i in range(steps + 1)]
    inner = [[(guide_radius - half_width) * math.cos(sweep * i / steps), (guide_radius - half_width) * math.sin(sweep * i / steps), 0.0] for i in range(steps, -1, -1)]
    end_cap = [
        [
            guide_radius * math.cos(sweep) + half_width * math.cos(sweep + math.pi * j / steps),
            guide_radius * math.sin(sweep) + half_width * math.sin(sweep + math.pi * j / steps),
            0.0,
        ]
        for j in range(1, steps)
    ]
    start_cap = [[guide_radius + half_width * math.cos(math.pi + math.pi * j / steps), half_width * math.sin(math.pi + math.pi * j / steps), 0.0] for j in range(1, steps)]
    points = outer + end_cap + inner + start_cap
    return PocketSpec.build(
        name=f"arcchan_r{guide_radius:g}_w{half_width:g}_s{sweep_deg:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"guide_radius": guide_radius, "half_width": half_width, "sweep_deg": sweep_deg},
    )


def axis_clearance(spec: PocketSpec, x: float, y: float) -> float:
    """Closed-form clearance at a medial-axis point of an analytic pocket.

    Args:
        spec: An instance produced by this module.
        x: Axis-point x coordinate.
        y: Axis-point y coordinate.

    Returns:
        The exact distance from the point to the pocket boundary.

    Raises:
        DegeneratePocketError: The spec did not come from this module.
    """
    params = spec.params
    if "half_width" in params:
        return float(params["half_width"])
    if "radius" in params:
        return float(params["radius"]) - math.hypot(x, y)
    if "height" in params:
        return 0.5 * float(params["height"]) - abs(y)
    raise DegeneratePocketError(f"{spec.name}: not an analytic family instance.")


def _require_channel(half_width: float, tool_diameter: float, family: str) -> None:
    """Raise when a channel cannot admit the tool with working room.

    Args:
        half_width: Channel half-width.
        tool_diameter: Cutter diameter.
        family: Family name, for the error message.

    Raises:
        DegeneratePocketError: The channel is too narrow.
    """
    needed = MIN_CHANNEL_CLEARANCE_FACTOR * 0.5 * tool_diameter
    if half_width < needed:
        raise DegeneratePocketError(f"{family}: half_width {half_width:g} is below the {needed:g} required for a tool of diameter {tool_diameter:g}.")
