"""Boundary-complexity sweeps: vertex count, and straight/curved composition.

Two axes, deliberately separated:

* `regular_ngon` / `ngon_sweep` vary the vertex count at FIXED area, so the only
  thing changing is how many boundary elements the arrangement and the medial
  axis must carry -- not how much material there is to clear.
* `arc_fraction` / `arc_fraction_sweep` vary how much of the boundary is curved.

The second axis cannot be isolated the way the first is: a side only becomes
curved by acquiring tessellation vertices, so the arc ratio and the element count
move together. Rather than hide that, every arc-fraction instance records its
vertex count in `params["vertices"]`, and the k-gon sweep is the vertex-count-only
control a reader subtracts to attribute the difference. Holding the count fixed
instead -- by giving straight sides the same interior points, lying on the chord
-- was rejected: it would seed every instance with collinear vertices, replacing
the covariate with a degeneracy that `benchmarks.degeneracy` exists to isolate.

MEASURED (2026-08-21, tool 4.0, cap 120 deg, via `benchmarks.runner.run_spec`)
-- both axes move the numbers they exist to move, so neither sweep is decorative:

    k-gon, area held at 100:   k 3 -> 64   generate 0.0013 -> 0.0975 s (75x),
                               certify 0.35 -> 4.66 s (13x), operations 47 -> 767
    arc fraction, n=8, r=10:   ratio 0 -> 1  generate 0.0072 -> 0.159 s (22x),
                               certify 0.75 -> 10.29 s (14x), operations 191 -> 1191

Note what the k-gon sweep does NOT hold fixed: equal area does not mean equal
work. A 64-gon's medial axis carries an order of magnitude more branches than a
triangle's, so the toolpath grows even though the material removed does not. Area
is held constant to remove volume as an explanation, not to predict a flat cost.
"""

from __future__ import annotations

import math

from compas.geometry import Polygon

from benchmarks.errors import DegeneratePocketError
from benchmarks.errors import InvalidArcRatioError
from benchmarks.errors import InvalidSideCountError
from benchmarks.spec import PocketSpec

# A closed ring needs three sides before it encloses anything.
MIN_POLYGON_SIDES = 3

# Chord segments used to tessellate each bulged (arc) side. At ARC_BULGE_RATIO the
# bulge subtends about 106 degrees, and 12 chords keep the tessellation within
# 0.6% of that arc, while leaving the vertex count dominated by the requested n
# rather than by the tessellation.
ARC_SIDE_SEGMENTS = 12

# Bulge height of an "arc" side as a fraction of that side's chord length. A
# circular arc of sweep theta has sagitta/chord = tan(theta/4)/2, so 0.25 is the
# 106.3-degree arc: visibly curved, and its outward bulge meets its neighbours at
# a reflex corner -- exactly the medial-axis feature a straight side does not
# produce -- without folding the ring back on itself at small n.
ARC_BULGE_RATIO = 0.25


def regular_ngon(k: int, area: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A regular k-gon of the requested area.

    Holding area fixed isolates boundary complexity: only the element count
    changes, not how much material there is to clear.

    Args:
        k: Number of sides, at least three.
        area: Target enclosed area; must be finite and positive.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        InvalidSideCountError: Fewer than three sides were requested.
        DegeneratePocketError: The target area is NaN, zero, or negative.
    """
    _require_sides(k, "regular_ngon")
    if not math.isfinite(area) or area <= 0.0:
        raise DegeneratePocketError(f"regular_ngon: area must be finite and positive, got {area!r}.")
    # area = 0.5 * k * R^2 * sin(2*pi/k)  =>  R = sqrt(2*area / (k*sin(2*pi/k)))
    radius = math.sqrt(2.0 * area / (k * math.sin(2.0 * math.pi / k)))
    points = [[radius * math.cos(2.0 * math.pi * i / k), radius * math.sin(2.0 * math.pi * i / k), 0.0] for i in range(k)]
    return PocketSpec.build(
        name=f"ngon_k{k}",
        family="complexity",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"k": float(k), "area": area, "vertices": float(k)},
    )


def ngon_sweep(ks: tuple[int, ...], area: float, tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep the vertex count at fixed area.

    Args:
        ks: Vertex counts to generate.
        area: Target enclosed area.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per requested k, in input order.

    Raises:
        InvalidSideCountError: A requested count is below three.
        DegeneratePocketError: The target area is NaN, zero, or negative.
    """
    return [regular_ngon(k=k, area=area, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for k in ks]


def arc_fraction(n: int, arc_ratio: float, radius: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """An n-sided pocket where a fraction of the sides bulge outward as arcs.

    The first ``round(arc_ratio * n)`` sides are replaced by an outward bulge; the
    rest stay straight chords. The resulting vertex count is recorded in
    ``params["vertices"]`` because it co-varies with the ratio by construction.

    Args:
        n: Number of sides, at least three.
        arc_ratio: Fraction of sides rendered as outward arcs, in [0, 1].
        radius: Circumradius of the underlying n-gon; must be finite and positive.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        InvalidSideCountError: Fewer than three sides were requested.
        InvalidArcRatioError: The ratio is NaN or outside [0, 1].
        DegeneratePocketError: The circumradius is NaN, zero, or negative.
    """
    _require_sides(n, "arc_fraction")
    if not math.isfinite(arc_ratio) or not (0.0 <= arc_ratio <= 1.0):
        raise InvalidArcRatioError(f"arc_fraction: arc_ratio must lie in [0, 1], got {arc_ratio!r}.")
    if not math.isfinite(radius) or radius <= 0.0:
        raise DegeneratePocketError(f"arc_fraction: radius must be finite and positive, got {radius!r}.")
    corners = [(radius * math.cos(2.0 * math.pi * i / n), radius * math.sin(2.0 * math.pi * i / n)) for i in range(n)]
    arc_sides = int(round(arc_ratio * n))
    points: list[list[float]] = []
    for i in range(n):
        ax, ay = corners[i]
        bx, by = corners[(i + 1) % n]
        points.append([ax, ay, 0.0])
        if i < arc_sides:
            points.extend(_bulged_side(ax, ay, bx, by))
    return PocketSpec.build(
        name=f"arcfrac_n{n}_r{arc_ratio:g}",
        family="complexity",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"n": float(n), "arc_ratio": arc_ratio, "radius": radius, "vertices": float(len(points))},
    )


def arc_fraction_sweep(n: int, ratios: tuple[float, ...], radius: float, tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep the straight/curved composition of the boundary.

    Args:
        n: Number of sides.
        ratios: Arc fractions to generate, each in [0, 1].
        radius: Circumradius of the underlying n-gon.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per ratio, in input order.

    Raises:
        InvalidSideCountError: Fewer than three sides were requested.
        InvalidArcRatioError: A requested ratio is NaN or outside [0, 1].
        DegeneratePocketError: The circumradius is NaN, zero, or negative.
    """
    return [arc_fraction(n=n, arc_ratio=r, radius=radius, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for r in ratios]


def _bulged_side(ax: float, ay: float, bx: float, by: float) -> list[list[float]]:
    """Interior points of an outward bulge from ``(ax, ay)`` to ``(bx, by)``.

    The chord is strictly positive by construction -- every caller has already
    established ``n >= 3`` and ``radius > 0``, so no two corners coincide -- which
    is why the normalisation below needs no floor.

    Args:
        ax: Start x.
        ay: Start y.
        bx: End x.
        by: End y.

    Returns:
        The interior tessellation points, excluding both endpoints.
    """
    dx, dy = bx - ax, by - ay
    chord = math.hypot(dx, dy)
    # Outward normal for a CCW ring is the clockwise perpendicular of the edge.
    nx, ny = dy / chord, -dx / chord
    sagitta = ARC_BULGE_RATIO * chord
    out: list[list[float]] = []
    for j in range(1, ARC_SIDE_SEGMENTS):
        t = j / ARC_SIDE_SEGMENTS
        bulge = sagitta * math.sin(math.pi * t)  # zero at both ends, peak at mid-side
        out.append([ax + dx * t + nx * bulge, ay + dy * t + ny * bulge, 0.0])
    return out


def _require_sides(sides: int, family: str) -> None:
    """Raise when a family is asked for a ring that cannot close.

    Args:
        sides: Requested side count.
        family: Family name, for the error message.

    Raises:
        InvalidSideCountError: The count is below three.
    """
    if sides < MIN_POLYGON_SIDES:
        raise InvalidSideCountError(f"{family}: at least {MIN_POLYGON_SIDES} sides are required, got {sides}.")
