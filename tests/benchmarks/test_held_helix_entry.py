"""The Held helical entry: a HELIX on circle 0's footprint (or its pencil's), ramp-out loops, hand-off into the published circle 0.

Rulings 2026-11-04: declared cut depth (1 mm = r), the D1 ramp default (2 deg),
an exact injection of circle 0 (published circles and connectors byte-identical),
a centre-cutting cutter (fail closed), Q4's theta_e at every junction of the
entry, ending with the hand-off into circle 0.
"""

from __future__ import annotations

import dataclasses
import json
import math
from fractions import Fraction

import pytest
from compas.geometry import Polygon

from benchmarks.held_helix_entry import HELD_STOCK_TOP
from benchmarks.held_helix_entry import HeldCircleFootprint
from benchmarks.held_helix_entry import HeldHelixEntry
from benchmarks.held_helix_entry import entry_admissibility
from benchmarks.held_helix_entry import entry_digest
from benchmarks.held_helix_entry import part_polygon
from benchmarks.held_reference_cases import DEFAULT_HELD_CUT_DEPTH
from benchmarks.held_toolpath_benchmark import GENERATORS
from benchmarks.held_toolpath_benchmark import GeneratedToolpath
from benchmarks.held_toolpath_benchmark import generate
from compas_cgal.engagement import _subtract_operation
from compas_cgal.helix import Helix
from compas_cgal.helix import ramp_turns
from compas_cgal.helix_units import DEFAULT_HELIX_PITCH_PER_RADIUS
from compas_cgal.helix_units import DEFAULT_MAX_RAMP_ANGLE_DEG
from compas_cgal.helix_units import EntryRampLength
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.admissibility import entry_rise_ratio
from compas_cgal.post.errors import NonCentreCuttingHelixError
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import Radian
from tools.held_path_digest import DIGEST_PATH
from tools.held_path_digest import DIGEST_WORKLOADS
from tools.held_path_digest import HeldPathDigests
from tools.held_path_digest import measure_workload

WORKLOADS = tuple(GENERATORS)
_PATHS: dict[str, GeneratedToolpath] = {}


def _path(workload: str) -> GeneratedToolpath:
    if workload not in _PATHS:
        _PATHS[workload] = generate(workload)
    return _PATHS[workload]


def _rho0_within_tool(path: GeneratedToolpath) -> bool:
    return Fraction(path.circles[0].radius.value) <= Fraction(path.case.tool_radius.value)


def _stock(path: GeneratedToolpath) -> Stock:
    polygon: Polygon = part_polygon(path.case)
    return Stock(polygon)


@pytest.mark.parametrize("workload", WORKLOADS)
def test_the_first_operation_is_a_helix_within_the_rulings(workload: str) -> None:
    path = _path(workload)
    entry = path.entry
    operations = entry.operations()
    assert operations[0].operation is OperationType.HELIX
    helix = operations[0].geometry
    assert isinstance(helix, Helix)
    assert helix is entry.helix
    # Declared cut depth, 1 mm = r, from the stock top.
    assert path.case.cut_depth == DEFAULT_HELD_CUT_DEPTH
    assert float(DEFAULT_HELD_CUT_DEPTH.value) == float(path.case.tool_radius.value) == 1.0
    assert helix.top_z.value == HELD_STOCK_TOP.value
    assert helix.descent == Fraction(path.case.cut_depth.value)
    assert helix.clockwise is path.circles[0].clockwise
    # The footprint is no wider than the tool: no core pin in the bore.
    assert Fraction(helix.footprint.radius.value) <= Fraction(path.case.tool_radius.value)
    assert helix.footprint.certifies_incidence()
    # Fewest turns the D1 default ramp admits.
    assert DEFAULT_MAX_RAMP_ANGLE_DEG == 2.0
    assert helix.turns == ramp_turns(helix.descent, helix.footprint.radius, DEFAULT_HELIX_PITCH_PER_RADIUS)
    # Every later operation is at the cut plane.
    assert all(op.operation is not OperationType.HELIX for op in operations[1:])


@pytest.mark.parametrize("workload", WORKLOADS)
def test_the_pitch_is_within_the_k_bound_and_minimal(workload: str) -> None:
    helix = _path(workload).entry.helix
    k = DEFAULT_HELIX_PITCH_PER_RADIUS
    assert k.value == 2.0 * math.pi * math.tan(math.radians(2.0))
    assert round(k.value, 4) == 0.2194
    bound = k.exact * Fraction(helix.footprint.radius.value)
    assert helix.pitch <= bound
    assert helix.ramp_within(k)
    if helix.turns > 1:
        assert ExactMillimetre(helix.descent / (helix.turns - 1)) > bound


@pytest.mark.parametrize("workload", [w for w in WORKLOADS if w != "figure8_crossed_skis_boundary" and w != "figure8_monstera_boundary_qualified_corners"])
def test_the_bore_equals_circle_zero_swept_disk(workload: str) -> None:
    path = _path(workload)
    assert _rho0_within_tool(path)
    circle0 = path.circles[0]
    footprint = path.entry.helix.footprint
    assert isinstance(footprint, HeldCircleFootprint)
    # The exact injection: circle 0's own doubles, not re-authored.
    assert (footprint.centre, footprint.radius, footprint.entry) == (circle0.center, circle0.radius, circle0.contact_point)
    assert path.entry.members == ()
    assert path.entry.hand_off is None
    r = float(path.case.tool_radius.value)
    bored = _stock(path)
    _subtract_operation(bored, path.entry.operations()[0], r)
    swept = _stock(path)
    swept.subtract_circle_sweep(circle0.center.x, circle0.center.y, circle0.radius.value, r)
    assert bored.exactly_equals(swept)
    assert bored.circle_sweep_removal(circle0.center.x, circle0.center.y, circle0.radius.value, r) is SweepRemoval.NOTHING


@pytest.mark.parametrize("workload", ["figure8_crossed_skis_boundary", "figure8_monstera_boundary_qualified_corners"])
def test_a_circle_zero_wider_than_the_tool_enters_on_its_pencil(workload: str) -> None:
    path = _path(workload)
    assert not _rho0_within_tool(path)
    entry = path.entry
    assert isinstance(entry.helix.footprint, IncidentLoop)
    assert entry.root is not None
    c0 = path.circles[0]
    # The pencil root is inscribed in circle 0's swept disk, exactly.
    gap = Fraction(c0.radius.value) - Fraction(entry.root.radius.value)
    dx = Fraction(entry.root.centre.x) - Fraction(c0.center.x)
    dy = Fraction(entry.root.centre.y) - Fraction(c0.center.y)
    assert gap >= 0 and gap * gap >= dx * dx + dy * dy
    # Every ramp-out member is on the root's pencil, strictly below the root.
    for member in entry.members:
        assert member.entry == entry.root.entry
        assert Fraction(member.radius.value) < Fraction(entry.root.radius.value)


@pytest.mark.parametrize("workload", WORKLOADS)
def test_every_junction_of_the_entry_is_within_theta_e(workload: str) -> None:
    """Independent replay: stock after each entry operation, Q4's theta_e at station 1 of every next loop, circle 0 last."""
    path = _path(workload)
    case, entry = path.case, path.entry
    cap = math.radians(float(case.tea_cap))
    theta_e = math.degrees(2.0 * math.asin(math.sqrt(float(entry_rise_ratio(Radian(cap), EntryRampLength.default()))) / 2.0))
    assert theta_e == pytest.approx(float(case.tea_cap) * math.pi / 45.0, rel=1e-12)
    admissibility = entry_admissibility(case)
    stock = _stock(path)
    r = float(case.tool_radius.value)
    operations = entry.operations()
    links = 0 if entry.hand_off is None else 1
    assert len(operations) == 1 + len(entry.members) + links
    _subtract_operation(stock, operations[0], r)
    for member, op in zip(entry.members, operations[1 : 1 + len(entry.members)], strict=True):
        assert op.entry_xy == (member.entry.x, member.entry.y)
        assert admissibility.loop_entry_within_rise(stock, member, False), f"{workload}: member {member.radius.value} over theta_e"
        assert admissibility.loop_within_cap(stock, member, False)
        _subtract_operation(stock, op, r)
    if entry.hand_off is not None:
        start, end = entry.hand_off
        assert end == path.circles[0].contact_point
        assert admissibility.straight_within_cap(stock, start, end)
        _subtract_operation(stock, operations[-1], r)
    circle0 = HeldCircleFootprint.build(path.circles[0])
    assert admissibility.loop_entry_within_rise(stock, circle0, False), f"{workload}: hand-off into circle 0 over theta_e"
    assert admissibility.loop_within_cap(stock, circle0, False)


@pytest.mark.parametrize("workload", ["figure5_refined"])
def test_a_cutter_that_is_not_centre_cutting_is_refused(workload: str) -> None:
    path = _path(workload)
    case = dataclasses.replace(path.case, centre_cutting=False)
    with pytest.raises(NonCentreCuttingHelixError):
        HeldHelixEntry.build(case, path.circles)


@pytest.mark.parametrize("workload", ["figure5_refined", "figure8_upper_boundary"])
def test_the_entry_digest_is_deterministic(workload: str) -> None:
    path = _path(workload)
    rebuilt = HeldHelixEntry.build(path.case, path.circles)
    assert rebuilt == path.entry
    assert entry_digest(rebuilt) == entry_digest(path.entry)


@pytest.mark.parametrize("workload", ["figure5_refined"])
def test_the_published_digests_are_unchanged_and_the_entry_digest_is_recorded(workload: str) -> None:
    """In-process regeneration of the fastest digest workload; `tools.held_path_digest --check` proves all four."""
    recorded = HeldPathDigests.from_json(DIGEST_PATH.read_text(encoding="utf-8")).workloads[workload]
    measured = measure_workload(workload).digest
    for field in ("initial", "circles", "connectors", "circle_count", "connector_count", "seam"):
        assert getattr(measured, field) == getattr(recorded, field), field
    assert recorded.entry is not None
    assert measured.entry == recorded.entry


def test_the_digest_record_carries_an_entry_per_workload() -> None:
    raw = json.loads(DIGEST_PATH.read_text(encoding="utf-8"))
    assert all("entry" in raw["workloads"][name] for name in DIGEST_WORKLOADS)
