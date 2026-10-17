"""The pockets and generators the machining-quality gate is asserted on.

Fixed here rather than in the test so that the gate, a report table, and every
figure describe the SAME six paths. Three pockets across two generators is the
smallest product in which a defect can be attributed: a finding on every pocket
is the generator's, a finding on every generator is the pocket's, and a finding
on one cell is neither.

* `rect_12x8` -- the convex reference. Its medial axis is two degree-3 nodes and
  a spine, so nothing about it is hard except the corners.
* `rect_20x12` -- the SAME shape at 2.5 times the area, and the instance on
  which the corner defect was first seen by eye: a loop of radius 0.0215 against
  a tool radius of 1.0, followed by a 1.956-long straight link cutting into the
  corner at 106.5 degrees. A criterion that never fires on the corpus is worth
  less than one that does, and this is the pocket that fires it. Whether the
  12x8 rectangle -- which has corners too -- fires it as well is itself a
  measurement worth reporting.
* `L_shape` -- adds the one feature neither rectangle has: a reflex corner,
  where the guide must turn and the skeleton carries more chains.

`rect_12x8` and `L_shape` have the same area, 96 square units, so their path
lengths and run times compare without normalising anything.
"""

from __future__ import annotations

import math
from typing import Callable
from typing import Dict
from typing import List
from typing import NewType
from typing import Tuple

from compas.geometry import Polygon

from benchmarks.families.analytic import rectangle
from benchmarks.models import GATE_MACHINE
from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from benchmarks.spec import PocketSpec
from compas_cgal.engagement_radial_toolpath import radius_regulated_toolpath
from compas_cgal.engagement_rho_toolpath import rho_regulated_toolpath
from compas_cgal.engagement_toolpath import engagement_controlled_toolpath
from compas_cgal.helix_units import EntryRampLength
from compas_cgal.pipeline.generate import generate_narrow
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import ClearanceZ
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius

# Cutter diameter for every gate instance. Matches `benchmarks.figure6`, where it
# was chosen so the skeleton carries several chains -- and therefore bridge cuts
# between them, which is where the worst engagement lives -- while a full run
# stays in seconds. It is also the diameter the corner defect was found at.
GATE_TOOL_DIAMETER = 2.0

# Gate-only engagement caps, in degrees. The 120-degree cap remains the corpus
# default and reporting snapshot. Forty degrees is the existing tight-cap case
# and the sole value in the approved 40--100 sweep where every pocket's two
# generators diverged; at 50--100 degrees L_shape stayed saturated (453/453 at
# the original 60-degree stop).
GateCapDegrees = NewType("GateCapDegrees", float)

GATE_CAP_DEG = GateCapDegrees(120.0)
GATE_ATTRIBUTION_CAP_DEG = GateCapDegrees(40.0)
GATE_DEFAULT_CAP_ID = f"cap-{GATE_CAP_DEG:g}-default"
GATE_ATTRIBUTION_CAP_ID = f"cap-{GATE_ATTRIBUTION_CAP_DEG:g}-attribution"
GATE_CAP_CASES: Tuple[Tuple[str, GateCapDegrees], ...] = (
    (GATE_DEFAULT_CAP_ID, GATE_CAP_DEG),
    (GATE_ATTRIBUTION_CAP_ID, GATE_ATTRIBUTION_CAP_DEG),
)

# Six by four tool diameters at `GATE_TOOL_DIAMETER`. Wide enough that the
# largest spine loop comes out at three times the tool radius -- well clear of
# the degeneracy boundary the L's arm width is derived against below -- and small
# enough that generation and measurement together stay in single-digit seconds.
GATE_RECT_WIDTH = 12.0
GATE_RECT_HEIGHT = 8.0

# The pocket the corner defect was first observed on. Ten by six tool diameters.
GATE_LARGE_RECT_WIDTH = 20.0
GATE_LARGE_RECT_HEIGHT = 12.0

# The L's arm width, in tool diameters, and the one number in this module that
# had to be DERIVED rather than picked. Along a channel of width W the clearance
# is W/2, so a loop centred on the spine has radius at most W/2 - r; a trochoid
# needs radius > r (see `benchmarks.quality.DEGENERATE_LOOP_RATIO`), hence
# W > 4r, hence an arm STRICTLY WIDER THAN TWO TOOL DIAMETERS. Exactly two tool
# diameters is the open degeneracy boundary. The prior empirical anecdote cannot
# be reproduced without its generator, circle-selection, entry-treatment, and
# operation-enumeration inputs, so it is deleted. Three tool diameters permits a
# spine-centred radius up to 2r, which leaves room for a trochoid.
L_ARM_TOOL_DIAMETERS = 3.0

# The L's two arm lengths, in tool diameters. Chosen so the instance's area comes
# out at 96 square units, equal to `rect_12x8`'s.
L_LONG_ARM_TOOL_DIAMETERS = 6.0
L_SHORT_ARM_TOOL_DIAMETERS = 5.0

# Instance names in emission order, so a report or a test can assert the gate's
# membership without rebuilding its geometry.
GATE_POCKET_NAMES: Tuple[str, ...] = ("rect_12x8", "rect_20x12", "L_shape")

# Generator names in the order the gate reports them: the two frozen columns,
# the rho-regulated reference column, and the product pipeline (Task 1.1).
GATE_GENERATOR_NAMES: Tuple[str, ...] = ("engagement_controlled", "radius_regulated", "rho_regulated", "pipeline")

# The two frozen generators the attribution witnesses compare (Ruling
# 1.1-C3-nonref): fixed by name, so registering further columns never changes
# which pair a two-generator attribution reads.
FROZEN_PAIR: Tuple[str, str] = ("engagement_controlled", "radius_regulated")

Generator = Callable[[PocketSpec], ToolpathResult]


def gate_pockets(
    tool_diameter: float = GATE_TOOL_DIAMETER,
    tea_cap_deg: GateCapDegrees = GATE_CAP_DEG,
) -> List[PocketSpec]:
    """The pockets the quality gate is asserted on, in `GATE_POCKET_NAMES` order.

    Args:
        tool_diameter: Cutter diameter. The L-shape's arm width is expressed in
            it, so the corpus keeps its meaning at any tool size.
        tea_cap_deg: Engagement cap in degrees. Defaults to the 120-degree
            reporting cap; the quality product also supplies the 40-degree
            attribution cap explicitly.

    Returns:
        The three instances.

    Raises:
        DegeneratePocketError: A tool too large for the fixed geometry.
        NonPositiveToolError: The diameter is NaN, zero, or negative.
        InvalidCapError: The cap is NaN or outside (0, 180].
    """
    arm = L_ARM_TOOL_DIAMETERS * tool_diameter
    long_arm = L_LONG_ARM_TOOL_DIAMETERS * tool_diameter
    short_arm = L_SHORT_ARM_TOOL_DIAMETERS * tool_diameter
    l_shape = Polygon([[0.0, 0.0, 0.0], [long_arm, 0.0, 0.0], [long_arm, arm, 0.0], [arm, arm, 0.0], [arm, short_arm, 0.0], [0.0, short_arm, 0.0]])
    return [
        rectangle(width=GATE_RECT_WIDTH, height=GATE_RECT_HEIGHT, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
        rectangle(width=GATE_LARGE_RECT_WIDTH, height=GATE_LARGE_RECT_HEIGHT, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
        PocketSpec.build(
            name="L_shape",
            family="analytic",
            polygon=l_shape,
            tool_diameter=tool_diameter,
            tea_cap_deg=tea_cap_deg,
            params={"arm_width": arm},
        ),
    ]


def gate_pocket(
    name: str,
    tool_diameter: float = GATE_TOOL_DIAMETER,
    tea_cap_deg: GateCapDegrees = GATE_CAP_DEG,
) -> PocketSpec:
    """One gate pocket by name.

    Args:
        name: One of `GATE_POCKET_NAMES`.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees. Defaults to the 120-degree
            reporting cap; callers may supply the 40-degree attribution cap.

    Returns:
        The instance.

    Raises:
        UnknownCorpusError: No gate pocket answers to that name.
    """
    from benchmarks.errors import UnknownCorpusError

    for spec in gate_pockets(tool_diameter, tea_cap_deg):
        if spec.name == name:
            return spec
    raise UnknownCorpusError(f"Unknown gate pocket {name!r}; expected one of {GATE_POCKET_NAMES}.")


def engagement_controlled(spec: PocketSpec) -> ToolpathResult:
    """Generate *spec*'s path with the advance-regulated generator.

    Args:
        spec: The instance to machine.

    Returns:
        The generated toolpath.
    """
    return engagement_controlled_toolpath(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=spec.tea_cap_deg,
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )


def radius_regulated(spec: PocketSpec) -> ToolpathResult:
    """Generate *spec*'s path with the advance- and radius-regulated generator.

    Args:
        spec: The instance to machine.

    Returns:
        The generated toolpath.
    """
    return radius_regulated_toolpath(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=spec.tea_cap_deg,
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )


def rho_regulated(spec: PocketSpec) -> ToolpathResult:
    """Generate *spec*'s path with the frozen rho-regulated generator, the pipeline's reference column.

    The arguments have the same shape as the frozen columns' adapters,
    `engagement_controlled` and `radius_regulated`: the pocket polygon and
    holes, the tool diameter, the cap and the clearance height.

    Args:
        spec: The instance to machine.

    Returns:
        The generated toolpath.
    """
    return rho_regulated_toolpath(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=spec.tea_cap_deg,
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )


# The gate's declared entry ramp for the pipeline's helical entry (USER 2026-10-16
# Q4): one tool diameter, the placeholder declared intent of
# `EntryRampLength.default`, stated here so the adapter passes it explicitly.
GATE_ENTRY_RAMP = EntryRampLength.default()


def pipeline(spec: PocketSpec) -> ToolpathResult:
    """Generate *spec*'s path with the product pipeline's narrow regime (`compas_cgal.pipeline`).

    The adapter from a benchmark `PocketSpec` to the pipeline's typed arguments
    (Ruling 1.1-gen-args): the pipeline imports nothing from `benchmarks`.
    The helix pitch is `GATE_MACHINE`'s own surrogate, the one the survey
    judges the ramp by, and the entry ramp is `GATE_ENTRY_RAMP`: both are
    passed explicitly, so generator and gate cannot disagree.

    Args:
        spec: The instance to machine.

    Returns:
        The generated toolpath, a `compas_cgal.pipeline.result.PipelineResult`.
    """
    return generate_narrow(
        spec.polygon,
        list(spec.holes),
        ToolRadius.build(spec.tool_radius),
        Radian(math.radians(spec.tea_cap_deg)),
        ClearanceZ.build(CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter),
        helix_pitch_per_radius=GATE_MACHINE.helix_pitch_per_radius,
        entry_ramp=GATE_ENTRY_RAMP,
    )


GATE_GENERATORS: Dict[str, Generator] = {
    "engagement_controlled": engagement_controlled,
    "radius_regulated": radius_regulated,
    "rho_regulated": rho_regulated,
    "pipeline": pipeline,
}
