"""Coordinate-precision sweeps: magnitude and significant decimals.

In an exact-constructions kernel, wall time depends on how often the lazy
interval filter fails and how many bits the fallback rationals carry. Both grow
with input precision and with construction depth. This family varies precision
alone, holding shape fixed, so growth attributable to bit length can be separated
from growth attributable to feature count.

`scale_sweep` rounds AFTER scaling, which is a deliberate choice and the only one
faithful to CAD: a part machined at 1000 mm carries the same four decimal places
as one at 1 mm, so growing the part grows the total significant digits. Scaling a
once-rounded unit shape instead -- the alternative that would make the instances
exactly similar across scales -- multiplies every exact rational by a constant and
so cannot move bit length by more than that constant's own width, i.e. it would be
a sweep guaranteed in advance not to measure anything. The price is that instances
at different scales are similar only to within the rounding, which is why the
tests derive their tolerance from `decimals` rather than asserting exact scaling.

MEASURED (2026-08-21, k=12, seed 1, cap 120 deg) -- READ THIS BEFORE READING A
REPORT ROW FROM THIS FAMILY. Neither axis moves cost on this branch, and the
digit axis barely moves even the quantity it was built to instrument.

Decimal-digit width of the INJECTED rationals, on a virgin stock at radius 10
(`benchmarks.instrument.probe_digits`):

    decimals        0      1      2      4      6      9     12
    max digits      3     34     34     34     34     34     34
    mean digits  1.58  19.75  26.08  27.17  26.96  26.83  27.04

The width is FLAT from one decimal upward. Coordinates enter the exact kernel as
IEEE-754 doubles, and every double is a dyadic rational p/2^n; a decimal that is
not itself dyadic -- 9.4 as much as 9.428571429 -- lands on a full 53-bit
mantissa and so carries full width no matter how few digits it was written with.
The one real step is 0 -> 1, where exactly-representable integers give an 11x
shorter rational. `decimals` is therefore a DYADIC-OR-NOT switch, not a
bit-length dial, which is why `DIGIT_SWEEP_DEFAULT` starts at 0: drop that entry
and the sweep is a flat line.

End to end, radius 10, tool 4.0:

    decimals        0      1      4      9
    certify (s)  1.313  1.214  1.191  1.186
    operations     317    303    301    301
    uncertified    137    131    130    130
    max digits     135    141    139    139

Even the d=0 column differs only because rounding to integers at radius 10 MOVES
THE VERTICES by up to 5% of the circumradius; that is a shape difference, not a
precision effect. After depletion all four sit at 135-141 digits, so whatever
head start exact integers give is gone by the first boolean construction.

Scale, decimals 4, tool = 0.1 x scale:

    scale          0.1     10    1000
    certify (s)  12.12  12.22   12.14
    operations    1145   1145    1145
    uncertified    552    552     552
    truly_exceed    47     47      47
    max digits     149    139     134

Operation counts and BOTH cap columns are bit-identical across four orders of
magnitude: the pipeline is scale-invariant, and certify time varies by 0.8%,
which is noise. Injected width does move monotonically -- 37 digits at scale 0.01
down to 31 at 10000 -- and in the direction opposite to intuition, because a
small value rounded to four decimals needs more FRACTIONAL bits (a wider 2^n
denominator) while a large one is dominated by its integer part.

The conclusion this family exists to support, then, is a negative one, and it is
worth more than a curve: exact-kernel cost here is driven by construction depth,
not by input precision or magnitude. The plan's premise -- that tidy benchmark
coordinates understate the cost of nine-significant-digit CAD input -- does not
survive measurement past the integer boundary.
"""

from __future__ import annotations

import math
import random

from compas.geometry import Polygon

from benchmarks.errors import DegeneratePocketError
from benchmarks.errors import InvalidDecimalsError
from benchmarks.errors import InvalidSideCountError
from benchmarks.spec import PocketSpec

# A closed ring needs three sides before it encloses anything.
MIN_POLYGON_SIDES = 3

# Radial jitter as a fraction of the circumradius. Large enough to force distinct
# non-round coordinates on every vertex, small enough to keep the ring simple.
PRECISION_JITTER_RATIO = 0.05

# Tool diameter as a fraction of the circumradius, so the tool scales with the
# geometry and the scale sweep varies magnitude ONLY, never relative feature size.
TOOL_TO_RADIUS_RATIO = 0.1

# Decimal counts spanning the only step this axis can resolve. 0 is load-bearing:
# it is the exactly-representable-integer case, and it is the ONLY entry that
# differs from the others in injected rational width (see the module docstring).
# The counts above it are kept so a future kernel that does become
# precision-sensitive shows up here rather than needing a new sweep.
DIGIT_SWEEP_DEFAULT: tuple[int, ...] = (0, 1, 2, 4, 6, 9, 12)

# Circumradii spanning six orders of magnitude, from a watchmaking part to a
# gantry-scale one, with the tool scaled to match so the machining problem is the
# same at every entry and only the magnitude of the numbers changes.
SCALE_SWEEP_DEFAULT: tuple[float, ...] = (0.01, 1.0, 100.0, 10000.0)


def perturbed_ngon(k: int, radius: float, decimals: int, seed: int, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A k-gon whose vertices are jittered and rounded to *decimals* places.

    Args:
        k: Number of vertices, at least three.
        radius: Circumradius before jitter; must be finite and positive.
        decimals: Decimal places retained on every coordinate; must not be
            negative.
        seed: Seed for the deterministic jitter.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        InvalidSideCountError: Fewer than three vertices were requested.
        InvalidDecimalsError: The decimal count is negative.
        DegeneratePocketError: The circumradius is NaN, zero, or negative, or
            rounding collapsed the ring below the tool's working area.
    """
    if k < MIN_POLYGON_SIDES:
        raise InvalidSideCountError(f"perturbed_ngon: at least {MIN_POLYGON_SIDES} vertices are required, got {k}.")
    if decimals < 0:
        raise InvalidDecimalsError(f"perturbed_ngon: decimals must not be negative, got {decimals}.")
    if not math.isfinite(radius) or radius <= 0.0:
        raise DegeneratePocketError(f"perturbed_ngon: radius must be finite and positive, got {radius!r}.")
    rng = random.Random(seed)
    points: list[list[float]] = []
    for i in range(k):
        angle = 2.0 * math.pi * i / k
        r = radius * (1.0 + PRECISION_JITTER_RATIO * (2.0 * rng.random() - 1.0))
        points.append([round(r * math.cos(angle), decimals), round(r * math.sin(angle), decimals), 0.0])
    return PocketSpec.build(
        name=f"prec_k{k}_r{radius:g}_d{decimals}",
        family="precision",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"k": float(k), "radius": radius, "decimals": float(decimals), "seed": float(seed)},
    )


def scale_sweep(k: int, scales: tuple[float, ...], decimals: int, seed: int, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep coordinate magnitude at fixed shape and fixed relative tool size.

    The tool scales with the geometry, so every instance is the same machining
    problem expressed in different units -- only the magnitude of the numbers the
    exact kernel injects changes.

    Args:
        k: Number of vertices.
        scales: Circumradii to generate.
        decimals: Decimal places retained on every coordinate, applied AFTER
            scaling; see the module docstring for why.
        seed: Seed for the deterministic jitter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per scale, in input order.

    Raises:
        InvalidSideCountError: Fewer than three vertices were requested.
        InvalidDecimalsError: The decimal count is negative.
        DegeneratePocketError: A scale is not positive, or rounding collapsed the
            ring at that scale.
    """
    return [perturbed_ngon(k=k, radius=s, decimals=decimals, seed=seed, tool_diameter=TOOL_TO_RADIUS_RATIO * s, tea_cap_deg=tea_cap_deg) for s in scales]


def digit_sweep(k: int, radius: float, decimal_counts: tuple[int, ...], seed: int, tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep significant decimals at fixed shape and magnitude.

    Args:
        k: Number of vertices.
        radius: Circumradius.
        decimal_counts: Decimal-place counts to generate.
        seed: Seed for the deterministic jitter.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per decimal count, in input order.

    Raises:
        InvalidSideCountError: Fewer than three vertices were requested.
        InvalidDecimalsError: A requested decimal count is negative.
        DegeneratePocketError: The circumradius is NaN, zero, or negative, or
            rounding collapsed the ring.
    """
    return [perturbed_ngon(k=k, radius=radius, decimals=d, seed=seed, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for d in decimal_counts]
