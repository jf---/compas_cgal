"""The helical entry on the six gate cases (Task 1.2 generator core).

Through the gate's own adapter (`GATE_GENERATORS["pipeline"]`) at both caps on
every gate pocket: exactly one HELIX, at the FIRST_CUT root's exact entry,
followed by the pencil members and the root loop; every pencil member is
re-decided here exactly on the replayed stock (the junction at the declared
rise limit, then the full cap), the survey samples no cap exceedance, and no
PLUNGE meets material.
"""

from __future__ import annotations

import functools
import math
from fractions import Fraction

import pytest
from compas.geometry import Circle

from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import gate_pocket
from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from benchmarks.survey import PathSurvey
from benchmarks.survey import survey_path
from compas_cgal.engagement import _subtract_operation
from compas_cgal.helix import Helix
from compas_cgal.helix_units import EntryRampLength
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import entry_rise_ratio
from compas_cgal.pipeline.generate import DEFAULT_MAX_PASSES
from compas_cgal.pipeline.guide import Guide
from compas_cgal.pipeline.guide import pipeline_spacing
from compas_cgal.pipeline.narrow import orientations
from compas_cgal.pipeline.result import ForcedEntryReason
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.stock import Stock
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import Clearance
from compas_cgal.units import ClearanceZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

CASES = [
    pytest.param(cap, pocket, id=f"{cap_id}-{pocket}", marks=pytest.mark.xdist_group(f"helix-entry-{cap_id}-{pocket}"))
    for cap_id, cap in GATE_CAP_CASES
    for pocket in GATE_POCKET_NAMES
]


@functools.cache
def _generated(cap: float, pocket: str) -> tuple[PipelineResult, PathSurvey]:
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
    result = GATE_GENERATORS["pipeline"](spec)
    assert isinstance(result, PipelineResult)
    return result, survey_path(spec, result)


def _first_cut_root(cap: float, pocket: str) -> IncidentLoop:
    """The 1.1 FIRST_CUT root, derived from the guide as `generate_narrow` derives it: the first pending chain's first loop."""
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
    tool = ToolRadius.build(spec.tool_radius)
    guide = Guide.build(
        spec.polygon,
        list(spec.holes),
        tool,
        pipeline_spacing(spec.tool_diameter, Radian(math.radians(cap))),
        Clearance.build(RADIAL_CLEARANCE_FRACTION * spec.tool_diameter),
        climb=True,
        max_passes=DEFAULT_MAX_PASSES,
        clearance_z=ClearanceZ.build(CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter),
    )
    first = next(offered for offered in (orientations(chain, tool) for chain in guide.chains) if offered)[0]
    root = first.loops[first.first]
    assert root is not None
    return root


def _is_loop(op: ToolpathOperation, loop: IncidentLoop) -> bool:
    return isinstance(op.geometry, Circle) and (float(op.geometry.frame.point[0]), float(op.geometry.frame.point[1]), float(op.geometry.radius)) == (
        loop.centre.x,
        loop.centre.y,
        loop.radius.value,
    )


def _pencil(cap: float, pocket: str) -> tuple[ToolpathOperation, list[ToolpathOperation]]:
    """The helix (operation 0) and every operation from it up to and including the root loop's one emission."""
    result, _survey = _generated(cap, pocket)
    helix_op = result.operations[0]
    assert helix_op.operation is OperationType.HELIX and isinstance(helix_op.geometry, Helix)
    root = _first_cut_root(cap, pocket)
    emissions = [index for index, op in enumerate(result.operations) if _is_loop(op, root)]
    assert len(emissions) == 1, f"the root loop is emitted {len(emissions)} times"
    return helix_op, list(result.operations[1 : emissions[0] + 1])


def _loop(op: ToolpathOperation) -> IncidentLoop:
    """The exactly incident loop an emitted pencil circle carries."""
    centre = op.geometry.frame.point
    assert op.entry_xy is not None
    loop = IncidentLoop.build(
        Point2[WorldXY](Millimetre(float(centre[0])), Millimetre(float(centre[1]))),
        GuideRadius(Millimetre(float(op.geometry.radius))),
        Radian(math.atan2(op.entry_xy[1] - float(centre[1]), op.entry_xy[0] - float(centre[0]))),
    )
    assert (loop.entry.x, loop.entry.y) == op.entry_xy
    return loop


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_exactly_one_helix_and_it_opens_the_path(cap: float, pocket: str) -> None:
    result, _survey = _generated(cap, pocket)
    assert [index for index, op in enumerate(result.operations) if op.operation is OperationType.HELIX] == [0]


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_the_ledger_has_one_first_cut_at_the_helix_entry(cap: float, pocket: str) -> None:
    result, _survey = _generated(cap, pocket)
    helix_op, members = _pencil(cap, pocket)
    first_cuts = [entry for entry in result.ledger.entries if entry.reason is ForcedEntryReason.FIRST_CUT]
    assert len(first_cuts) == 1
    (first,) = first_cuts
    assert first.into_material and (first.entry.x, first.entry.y) == helix_op.entry_xy
    assert members, "the root loop follows the helix"


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_the_pencil_is_cut_loops_entered_at_the_root_entry_rising_to_the_root(cap: float, pocket: str) -> None:
    """Between the helix and the root loop's one emission: only CUT circles, each entered exactly at P, radii strictly rising; the root is last."""
    helix_op, members = _pencil(cap, pocket)
    root = _first_cut_root(cap, pocket)
    assert helix_op.entry_xy == (root.entry.x, root.entry.y)
    assert all(op.operation is OperationType.CUT and isinstance(op.geometry, Circle) and op.entry_xy == helix_op.entry_xy for op in members)
    radii = [Fraction(float(op.geometry.radius)) for op in members]
    assert all(a < b for a, b in zip(radii, radii[1:]))
    assert _is_loop(members[-1], root)


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_every_pencil_member_passes_the_junction_and_the_full_cap_exactly(cap: float, pocket: str) -> None:
    """Re-decided on the replayed stock: the helix bore, then each member after the one before; the root loop is the last member."""
    result, _survey = _generated(cap, pocket)
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
    tool = ToolRadius.build(spec.tool_radius)
    cap_rad = Radian(math.radians(cap))
    admissibility = Admissibility.build(
        tool, cap_ratio(cap_rad), centre_domain(spec.polygon, list(spec.holes), tool), entry_rise_ratio=entry_rise_ratio(cap_rad, EntryRampLength.default())
    )
    helix_op, members = _pencil(cap, pocket)
    stock = Stock(spec.polygon, list(spec.holes))
    _subtract_operation(stock, helix_op, spec.tool_radius)
    previous = helix_op.geometry.footprint
    for op in members:
        loop = _loop(op)
        assert Fraction(loop.radius.value) - Fraction(previous.radius.value) <= Fraction(spec.tool_radius)
        assert admissibility.loop_entry_within_rise(stock, loop, op.clockwise)
        assert admissibility.loop_within_cap(stock, loop, op.clockwise)
        _subtract_operation(stock, op, spec.tool_radius)
        previous = loop


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_the_survey_samples_no_cap_exceedance(cap: float, pocket: str) -> None:
    _result, survey = _generated(cap, pocket)
    assert [int(motion.index) for motion in survey.motions if motion.cap_exceeded] == []


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_no_plunge_meets_material(cap: float, pocket: str) -> None:
    result, survey = _generated(cap, pocket)
    plunges = {index for index, op in enumerate(result.operations) if op.operation is OperationType.PLUNGE}
    assert not plunges & {int(index) for index in survey.material_plunge_indices}
    assert all(not entry.into_material for entry in result.ledger.entries if entry.reason is not ForcedEntryReason.FIRST_CUT)
