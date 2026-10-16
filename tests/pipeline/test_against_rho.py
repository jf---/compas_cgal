"""Differential oracle against the frozen rho generator, on PUBLIC outputs only (Rulings 1.1-C7, 1.1-B-oracle).

No frozen module private is imported here. The port's bit-parity with rho's
privates is a logged scratch probe, not in the repository;
this file holds the port to rho's public constants, its returned toolpath and
declined regions, and -- the reference rule -- the measured gate cells.
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path

import pytest
from compas.geometry import Circle

from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import gate_pocket
from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from compas_cgal import engagement_rho_toolpath as rho
from compas_cgal import engagement_toolpath as ec
from compas_cgal.pipeline import admissibility
from compas_cgal.pipeline import guide as pipeline_guide
from compas_cgal.pipeline import narrow
from compas_cgal.pipeline.guide import Guide
from compas_cgal.pipeline.guide import pipeline_spacing
from compas_cgal.pipeline.result import DeclinedRun
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.units import Clearance
from compas_cgal.units import ClearanceZ
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import load_gate_vectors

ROOT = Path(__file__).resolve().parents[2]
CAP_120 = Radian(math.radians(120.0))


def test_the_ported_constants_equal_rho_and_the_guide() -> None:
    assert pipeline_guide.GUIDE_STEP_TOOL_DIAMETERS == ec.GUIDE_STEP_TOOL_DIAMETERS
    assert pipeline_guide.GUIDE_STEPOVER_PITCH_MULTIPLE == ec.GUIDE_STEPOVER_ADVANCE_MULTIPLE
    assert pipeline_guide.GUIDE_SAMPLES_PER_RADIAN == ec.GUIDE_SAMPLES_PER_RADIAN
    assert pipeline_guide.SINGLE_STATION_TANGENT == ec.SINGLE_STATION_TANGENT
    assert narrow.MAX_ADVANCE_TOOL_DIAMETERS == ec.MAX_ADVANCE_TOOL_DIAMETERS
    assert narrow.DEGENERATE_LOOP_TOOL_RADII == rho.DEGENERATE_LOOP_TOOL_RADII
    assert narrow.MAX_LOOP_RADIUS_STEP_TOOL_RADII == rho.MAX_LOOP_RADIUS_STEP_TOOL_RADII
    assert narrow.MAX_CLEARANCE_SLOPE == rho.MAX_CLEARANCE_SLOPE
    assert admissibility.LOOP_PROBE_COUNT_AT_REST == ec.LOOP_PROBE_COUNT
    ring = tuple(admissibility.FULL_TURN_DEG * index / admissibility.LOOP_PROBE_COUNT_AT_REST for index in range(admissibility.LOOP_PROBE_COUNT_AT_REST))
    assert ring == ec.LOOP_PROBE_ANGLES_DEG


def _guide(name: str) -> Guide:
    spec = gate_pocket(name)
    return Guide.build(
        spec.polygon,
        list(spec.holes),
        ToolRadius.build(spec.tool_radius),
        pipeline_spacing(spec.tool_diameter, CAP_120),
        Clearance.build(RADIAL_CLEARANCE_FRACTION * spec.tool_diameter),
        climb=True,
        max_passes=1000,
        clearance_z=ClearanceZ.build(CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter),
    )


def _rho_result(name: str) -> rho.RhoToolpathResult:
    spec = gate_pocket(name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", rho.UnavoidableEngagementWarning)
        return rho.rho_regulated_toolpath(
            spec.polygon,
            tool_diameter=spec.tool_diameter,
            tea_cap_deg=spec.tea_cap_deg,
            holes=list(spec.holes),
            clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
        )


@pytest.mark.parametrize("name", GATE_POCKET_NAMES)
def test_every_rho_loop_sits_on_a_guide_station_with_its_radius(name: str) -> None:
    """rho's loops are guide stations above the floor: the port's guide carries each, centre and radius to the bit."""
    stations = {(s.centre.x, s.centre.y): s.radius.value for chain in _guide(name).chains for s in chain.stations}
    tool = ToolRadius.build(gate_pocket(name).tool_radius)
    loops = [op.geometry for op in _rho_result(name).operations if isinstance(op.geometry, Circle)]
    assert loops
    for loop in loops:
        centre = (float(loop.frame.point[0]), float(loop.frame.point[1]))
        assert stations[centre] == float(loop.radius)
        assert narrow.clears_floor(float(loop.radius), tool)


@pytest.mark.parametrize("name", GATE_POCKET_NAMES)
def test_the_declined_runs_equal_rhos_declined_regions(name: str) -> None:
    """The floor declines the station runs rho's public `declined_regions` reports (centres and counts)."""
    tool = ToolRadius.build(gate_pocket(name).tool_radius)
    ours: list[tuple[tuple[float, float], tuple[float, float], int]] = []
    for chain in _guide(name).chains:
        oriented = narrow.OrientedChain.build(chain, tool)
        runs = (DeclinedRun(chain=chain.origin, first_station=0, last_station=len(chain.stations) - 1),) if oriented is None else oriented.declined
        for run in runs:
            first, last = chain.stations[run.first_station].centre, chain.stations[run.last_station].centre
            ours.append(((first.x, first.y), (last.x, last.y), run.last_station - run.first_station + 1))
    theirs = [(region.first_center, region.last_center, region.station_count) for region in _rho_result(name).declined_regions]
    assert sorted(ours) == sorted(theirs)


@pytest.mark.parametrize(("cap_id", "pocket"), [(cap_id, pocket) for cap_id, _cap in GATE_CAP_CASES for pocket in GATE_POCKET_NAMES])
def test_no_gated_cell_is_worse_than_rho_regulated(cap_id: str, pocket: str) -> None:
    """The reference rule (USER 2026-10-15), read from the committed vectors: pipeline <= rho_regulated on every gated cell."""
    vectors = load_gate_vectors(ROOT / GATE_VECTORS_PATH)
    ours, reference = f"{cap_id}/pipeline/{pocket}", f"{cap_id}/rho_regulated/{pocket}"
    assert ours in vectors.cases, f"{ours} is not measured"
    assert reference in vectors.cases, f"{reference} is not measured"
    worse = {
        criterion.name: (vectors.value(ours, criterion.name), vectors.value(reference, criterion.name))
        for criterion in vectors.gated
        if vectors.value(ours, criterion.name) > vectors.value(reference, criterion.name)
    }
    assert worse == {}
