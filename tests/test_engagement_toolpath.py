"""Engagement-controlled trochoidal generator tests.

The point of `compas_cgal.engagement_toolpath` is regulation: it must beat the
unregulated `trochoidal_mat_toolpath_circular` on the engagement the audit
measures, on the same pocket at the same cap. These tests are sized so that a
full generate-plus-audit round trip stays tractable -- the audit replays every
motion against a depleting exact stock, which is the dominant cost.
"""

import math

import pytest
from compas.geometry import Circle, Line, Polygon

from compas_cgal.engagement import _subtract_operation, audit_toolpath_engagement
from compas_cgal.engagement_toolpath import (
    InvalidEngagementCapDegreesError,
    InvalidGuideResolutionError,
    NonPositiveToolDiameterError,
    engagement_controlled_toolpath,
)
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType, ToolpathResult, trochoidal_mat_toolpath_circular

# Small pocket, deliberately: generation and audit both replay against a depleting
# exact arrangement, so cost grows with the op count. 6x4 with a 2 mm tool keeps
# the structural, monotonicity, and determinism checks cheap while still
# exercising the full topology (a central skeleton chain plus four corner spokes).
POCKET = Polygon([[0, 0, 0], [6, 0, 0], [6, 4, 0], [0, 4, 0]])
TOOL_DIAMETER = 2.0
CAP_DEG = 120.0

# The engagement comparison needs a pocket in the genuine trochoidal regime. On
# 6x4 the clearance-derived loop radius (0.998) is below the tool radius, so every
# machining circle sweeps a full disk rather than an annulus and four of the five
# skeleton chains run through air the central chain already cleared -- a regime
# where neither generator's advance is the binding constraint. 10x6 gives a loop
# radius of 1.998 and a central chain that dominates the path, which is where a
# stepover proxy and a measured cap actually differ.
COMPARISON_POCKET = Polygon([[0, 0, 0], [10, 0, 0], [10, 6, 0], [0, 6, 0]])

# The unregulated generator's stepover chosen to match the regulated cap in the
# textbook radial-immersion relation TEA = 2*acos(1 - ae/r): ae = r*(1 - cos(cap/2))
# = 1.0*(1 - cos 60 deg) = 0.5. That is the stepover a CAM programmer would dial in
# aiming for a 120 deg engagement, which is exactly the comparison of interest --
# a geometric proxy for the cap versus a measured cap.
COMPARABLE_STEPOVER = 0.5 * TOOL_DIAMETER * (1.0 - math.cos(0.5 * math.radians(CAP_DEG)))


def _cut_length(result: ToolpathResult) -> float:
    """Total XY length of the cutting motions (circles plus bridges)."""
    total = 0.0
    for op in result.operations:
        if op.operation is not OperationType.CUT:
            continue
        geometry = op.geometry
        if isinstance(geometry, Circle):
            total += 2.0 * math.pi * float(geometry.radius)
        elif isinstance(geometry, Line):
            total += math.hypot(
                float(geometry.end[0]) - float(geometry.start[0]),
                float(geometry.end[1]) - float(geometry.start[1]),
            )
    return total


def _op_signature(result: ToolpathResult):
    """Hashable, exact-float signature of the emitted operation sequence."""
    signature = []
    for op in result.operations:
        geometry = op.geometry
        if isinstance(geometry, Circle):
            center = geometry.frame.point
            start = geometry.point_at(0.0)
            key = ("circle", float(geometry.radius), tuple(float(v) for v in center), tuple(float(v) for v in start))
        else:
            key = ("line", tuple(float(v) for v in geometry.start), tuple(float(v) for v in geometry.end))
        signature.append((op.operation.value, op.path_index, op.clockwise, key))
    return tuple(signature)


def test_rectangular_pocket_yields_cut_operations():
    result = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, CAP_DEG)

    assert isinstance(result, ToolpathResult)
    assert result.operations, "generator emitted no operations at all"

    cuts = [op for op in result.operations if op.operation is OperationType.CUT]
    loops = [op for op in cuts if isinstance(op.geometry, Circle)]
    bridges = [op for op in cuts if isinstance(op.geometry, Line)]
    assert loops, "no machining circles emitted"
    assert bridges, "no bridge cuts emitted between machining circles"

    # Every chain is entered by a plunge and left by a retract, so the vertical
    # moves come in pairs and the chain count is consistent.
    plunges = [op for op in result.operations if op.operation is OperationType.PLUNGE]
    retracts = [op for op in result.operations if op.operation is OperationType.RETRACT]
    assert len(plunges) == len(retracts) == len({op.path_index for op in cuts})

    assert result.polyline.ndim == 2 and result.polyline.shape[1] == 3
    assert result.polyline.shape[0] > 0


def _measured_exceedances(report, cap_rad: float) -> int:
    """Operations whose audited peak engagement is above the cap.

    Distinct from ``report.cap_violations``, which counts operations the audit
    could not CERTIFY. At a 120 deg cap the audit's fixed circular-station density
    (`AUDIT_ARC_STEP_FRACTION`) leaves no positive guarded cap, so it declines to
    certify every circular motion regardless of what it measured -- that count
    therefore tracks the number of circles, not their engagement. The peak TEA it
    reports is a measurement, and that is what regulation has to move.
    """
    return sum(1 for e in report.operations if e.max_tea > cap_rad)


def test_regulated_path_beats_unregulated_on_measured_engagement():
    cap_rad = math.radians(CAP_DEG)
    pocket = COMPARISON_POCKET

    regulated = engagement_controlled_toolpath(pocket, TOOL_DIAMETER, CAP_DEG)
    unregulated = trochoidal_mat_toolpath_circular(pocket, tool_diameter=TOOL_DIAMETER, stepover=COMPARABLE_STEPOVER)

    regulated_report = audit_toolpath_engagement(pocket, regulated, TOOL_DIAMETER, cap_rad)
    unregulated_report = audit_toolpath_engagement(pocket, unregulated, TOOL_DIAMETER, cap_rad)

    assert regulated_report.engaged_ops > 0, "regulated path never engages material"
    assert unregulated_report.cap_violations > 0, "baseline is unexpectedly already under the cap"

    # The whole point of the module: regulating the advance on the exact
    # per-position predicate must cut how many motions the audit measures above
    # the cap, and by a wide margin rather than by a hair.
    regulated_over = _measured_exceedances(regulated_report, cap_rad)
    unregulated_over = _measured_exceedances(unregulated_report, cap_rad)
    assert unregulated_over > 0, "baseline is unexpectedly already under the cap"
    assert regulated_over * 2 <= unregulated_over, f"regulated {regulated_over} vs unregulated {unregulated_over} measured exceedances"

    # And the audit's own certifiability count must not get worse either.
    assert regulated_report.cap_violations < unregulated_report.cap_violations


def test_only_chain_entry_loops_are_measured_above_the_cap():
    """The 4-probe accept rule holds up against the audit's 20-station sweep.

    The generator evaluates four positions per machining circle; the audit
    re-measures twenty. This pins the gap between them: on this pocket the only
    motions the audit finds above the cap are the chain-entry loops, which the
    generator already refuses and reports through `UnavoidableEngagementWarning`
    (the tool meeting virgin stock cuts a full slot, whatever the advance). It is a
    measurement on one pocket, not a proof that four probes suffice in general.
    """
    cap_rad = math.radians(CAP_DEG)
    result = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, CAP_DEG)
    report = audit_toolpath_engagement(POCKET, result, TOOL_DIAMETER, cap_rad)

    entry_indices = set()
    for index, op in enumerate(result.operations):
        if op.operation is OperationType.PLUNGE:
            entry_indices.add(index + 1)  # the machining circle the plunge lands in

    above_cap = {e.op_index: math.degrees(e.max_tea) for e in report.operations if e.max_tea > cap_rad}
    assert set(above_cap) == entry_indices, f"above cap {above_cap}, chain-entry loops at {sorted(entry_indices)}"


def test_regulated_path_clears_everything_the_tool_can_reach():
    """The advance search must not buy low engagement by skipping material.

    Nothing in the accept/reject rule mentions coverage: a candidate whose probes
    all sit in void is accepted regardless of what lies between. What keeps the
    walk honest is `MAX_ADVANCE_TOOL_DIAMETERS`, which bounds the advance at the
    annulus width so consecutive machining circles always overlap. This checks the
    consequence directly -- replay the path and confirm nothing survives that a
    round tool of this radius could have reached. Residue nearer the wall than the
    tool radius is the round tool's inherent corner residue, not a generator defect.
    """
    tool_radius = 0.5 * TOOL_DIAMETER
    result = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, CAP_DEG)

    stock = Stock(POCKET)
    for op in result.operations:
        geometry = op.geometry
        if isinstance(geometry, Line):
            z_start, z_end = float(geometry.start[2]), float(geometry.end[2])
            if z_start != z_end:
                if z_end < z_start:  # plunge bores a full disk at the cutting plane
                    stock.subtract_disk(float(geometry.end[0]), float(geometry.end[1]), tool_radius)
                continue
            if z_start != 0.0:  # clearance-height rapid, cuts nothing
                continue
        _subtract_operation(stock, op, tool_radius)

    width, height = 6.0, 4.0
    columns, rows = 60, 40
    survivors = []
    for row in range(rows):
        y = (row + 0.5) * height / rows
        for column in range(columns):
            x = (column + 0.5) * width / columns
            wall_distance = min(x, width - x, y, height - y)
            if wall_distance > tool_radius and stock.contains(x, y):
                survivors.append((round(x, 3), round(y, 3), round(wall_distance, 3)))
    assert not survivors, f"material left further than the tool radius from the wall: {survivors[:10]}"


def test_tighter_cap_is_not_shorter():
    tight = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, 90.0)
    loose = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, 170.0)
    # A tighter cap forces shorter advances, hence more machining circles over the
    # same guide, hence at least as much cutting travel.
    assert _cut_length(tight) >= _cut_length(loose)


def test_generation_is_deterministic():
    first = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, CAP_DEG)
    second = engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, CAP_DEG)
    assert _op_signature(first) == _op_signature(second)


@pytest.mark.parametrize("bad_cap_deg", [0.0, -1.0, 180.001, 360.0, float("nan")])
def test_rejects_out_of_range_cap(bad_cap_deg):
    with pytest.raises(InvalidEngagementCapDegreesError):
        engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, bad_cap_deg)


@pytest.mark.parametrize("bad_diameter", [0.0, -2.0])
def test_rejects_non_positive_tool_diameter(bad_diameter):
    with pytest.raises(NonPositiveToolDiameterError):
        engagement_controlled_toolpath(POCKET, bad_diameter, CAP_DEG)


def test_rejects_guide_resolution_coarser_than_the_advance_bound():
    # A guide step at or above the maximum admissible advance leaves the search no
    # bracket to bisect, so it is rejected at the boundary instead of silently
    # degenerating into a fixed-stepover walk.
    with pytest.raises(InvalidGuideResolutionError):
        engagement_controlled_toolpath(POCKET, TOOL_DIAMETER, CAP_DEG, guide_step_tool_diameters=1.0)
