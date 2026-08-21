"""The MATHSM baseline: constant machining-circle spacing, cap met by brute-force search.

Held & Pfeiffer compare their engagement-controlled paths against MATHSM, whose
machining circles are spaced at a CONSTANT step along the medial axis. The
relationship between that step and the maximum engagement it produces has no
closed form, so the protocol is empirical: vary the spacing, measure the maximum
engagement each spacing gives, and for a required cap take the shortest path among
those that met it.

This module reproduces the PROTOCOL, not Held's implementation or his pocket. The
constant-spacing generator is this project's own
`trochoidal_mat_toolpath_circular` driven by its `stepover` knob, and the
engagement measurement is this project's exact audit rather than the paper's
admitted discretisation. Nothing here reproduces a published number; what it
reproduces is a like-for-like curve a reader can regenerate.

THE RELATION IS NOT MONOTONE, which is why the search is brute force and the
selection is a minimum over all compliant trials rather than a bisection. Measured
on the 20x12 pocket at tool 2.0, the after-entry maximum engagement runs 131.14,
110.11, 120.99, 98.73, 119.42, 131.49, 117.79, 128.84, 138.38, 166.54, 189.56,
205.73 degrees across the twelve default spacings -- it falls and rises three
times before settling into a monotone climb, and its minimum sits in the middle of
the sweep rather than at its fine end. Bisecting on spacing would have found a
compliant trial and missed shorter ones.

Spacings are expressed in TOOL DIAMETERS so a sweep means the same thing at every
scale, and engagement is measured through `benchmarks.pathmetrics` so the entry
cut -- a full slot for any generator, at any spacing -- cannot pin the axis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from benchmarks.pathmetrics import PathMetrics
from benchmarks.pathmetrics import measure_path
from benchmarks.spec import PocketSpec
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

# Trial spacings in tool diameters. The lower end is past the turn where finer
# spacing stops helping -- on the reference pocket 0.025 measures 131.14 degrees
# against 98.73 at 0.1 -- and is kept so the sweep shows that turn rather than
# implying engagement falls without limit. The upper end is where the path is a
# rough slot at any cap worth setting. Twelve trials sweep the reference pocket in
# about three minutes, dominated by the audit rather than by generation.
SPACING_SWEEP_TOOL_DIAMETERS: Tuple[float, ...] = (0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6)


@dataclass(frozen=True)
class MathsmPoint:
    """One constant-spacing trial.

    Attributes:
        spacing_tool_diameters: The stepover used, as a multiple of the tool
            diameter.
        metrics: What that spacing cost and how hard it cut.
    """

    spacing_tool_diameters: float
    metrics: PathMetrics


def constant_spacing_path(spec: PocketSpec, spacing_tool_diameters: float) -> ToolpathResult:
    """Generate one constant-spacing path.

    The single place the baseline fixes its generation parameters, so a trial, a
    test, and a figure all describe the same path.

    Args:
        spec: The pocket and tool. Its engagement cap is NOT consulted: this
            generator has no cap input, which is exactly what makes it the
            baseline.
        spacing_tool_diameters: Stepover as a multiple of the tool diameter.

    Returns:
        The generated toolpath.
    """
    return trochoidal_mat_toolpath_circular(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        stepover=spacing_tool_diameters * spec.tool_diameter,
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )


def measure_constant_spacing(spec: PocketSpec, spacing_tool_diameters: float) -> MathsmPoint:
    """Generate and measure one trial.

    Args:
        spec: The pocket and tool.
        spacing_tool_diameters: Stepover as a multiple of the tool diameter.

    Returns:
        The trial.
    """
    return MathsmPoint(spacing_tool_diameters=spacing_tool_diameters, metrics=measure_path(spec, constant_spacing_path(spec, spacing_tool_diameters)))


def sweep_spacing(spec: PocketSpec, spacings_tool_diameters: Sequence[float]) -> List[MathsmPoint]:
    """Measure every requested spacing.

    Args:
        spec: The pocket and tool.
        spacings_tool_diameters: Stepovers to trial, as multiples of the tool
            diameter, in input order.

    Returns:
        One trial per spacing, in input order.
    """
    return [measure_constant_spacing(spec, spacing) for spacing in spacings_tool_diameters]


def shortest_within_cap(points: Sequence[MathsmPoint], cap_deg: float) -> Optional[MathsmPoint]:
    """The shortest trial whose measured engagement respects *cap_deg*.

    Held records the shorter path whenever several spacings give roughly the same
    maximum engagement; minimising length over every compliant trial is that rule
    stated without the "roughly".

    Compliance is read off `max_tea_after_entry_deg`, never the raw maximum. The
    raw maximum is the full slot the entry cut takes, so filtering on it would
    reject every trial at every cap and the baseline curve would silently not
    exist. This is a MEASURED comparison reproducing the paper's protocol, not a
    certificate: it says no audited station exceeded the cap away from the entry,
    which is a stronger statement than the paper's discretisation makes and a
    weaker one than a continuous proof.

    Args:
        points: Trials from `sweep_spacing`.
        cap_deg: The engagement cap in degrees, met when the measurement is at or
            below it.

    Returns:
        The shortest compliant trial, or None when none comply.
    """
    compliant = [p for p in points if p.metrics.max_tea_after_entry_deg <= cap_deg]
    return min(compliant, key=lambda p: p.metrics.length) if compliant else None
