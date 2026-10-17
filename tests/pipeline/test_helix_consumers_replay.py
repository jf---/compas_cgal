"""Replay consumers of a HELIX operation (Task 1.2 consumer set A: replay and removal).

A helix replays as its exact bore: the full-turn sweep of its footprint, which
is the disk of radius ``rho_h + r`` because the footprint satisfies
``rho_h <= r``. Both replays (the engagement audit and the depletion replay)
remove exactly that disk; neither skips it. Like a plunge, the helix is not
cap-sampled (E8): its record is the plunge's unmeasured-by-design record, and
the motion after it is measured against the bored stock. Every category and
kind dispatch is exhaustive: an unknown category raises a named error.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

import benchmarks.depletion as depletion
import compas_cgal.engagement as engagement
from benchmarks.depletion import CutMotion
from benchmarks.depletion import HelixEvent
from benchmarks.depletion import ReplayKind
from benchmarks.depletion import replay_cuts
from benchmarks.depletion import replay_depletion
from benchmarks.depletion import replay_events
from benchmarks.errors import UnreplayableOperationError
from benchmarks.families.analytic import rectangle
from benchmarks.spec import PocketSpec
from compas_cgal.engagement import UnhandledReplayCategoryError
from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

EMPTY_POLYLINE = np.zeros((0, 3), dtype=np.float64)
TOOL_DIAMETER = 2.0
CUT_Z = 0.0
CLEARANCE_Z = 4.0
# Footprint radius below the tool radius (1 mm), so the full-turn sweep is a disk.
FOOTPRINT_RADIUS = 0.5
# The boundary footprint rho_h == r (tool radius 1 mm): the inner sweep radius is 0.
BOUNDARY_FOOTPRINT_RADIUS = 1.0
# A loop whose swept disk (0.25 + 1 = 1.25 mm) lies strictly inside the bore (1.5 mm).
INNER_LOOP_RADIUS = 0.25
RAMP_DEG = 2.0
TEA_CAP_RAD = 2.0
BOGUS_CATEGORY = "bogus"


def _spec() -> PocketSpec:
    return rectangle(width=8.0, height=6.0, tool_diameter=TOOL_DIAMETER, tea_cap_deg=120.0)


def _helix_op(footprint_radius: float = FOOTPRINT_RADIUS) -> ToolpathOperation:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(0.0), Millimetre(0.0)), GuideRadius(Millimetre(footprint_radius)), Radian(0.0))
    plane = CutPlane.build(CutZ.build(CUT_Z), ClearanceZ.build(CLEARANCE_Z))
    helix = Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(RAMP_DEG)), clockwise=False)
    entry = helix.entry_xy
    return ToolpathOperation(geometry=helix, operation=OperationType.HELIX, path_index=0, clockwise=False, entry_xy=(entry.x, entry.y))


def _inner_loop() -> ToolpathOperation:
    circle = Circle(INNER_LOOP_RADIUS, frame=Frame([0.0, 0.0, CUT_Z], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]))
    return ToolpathOperation(geometry=circle, operation=OperationType.CUT, path_index=0)


def _plunge() -> ToolpathOperation:
    return ToolpathOperation(geometry=Line([0.0, 0.0, CLEARANCE_Z], [0.0, 0.0, CUT_Z]), operation=OperationType.PLUNGE, path_index=0)


def _bored(spec: PocketSpec, helix: Helix) -> Stock:
    """The virgin stock minus the bore disk of radius ``rho_h + r``, built independently."""
    bore_radius = Fraction(helix.footprint.radius.value) + Fraction(spec.tool_radius)
    assert Fraction(float(bore_radius)) == bore_radius, "bore radius is not an exact double"
    stock = Stock(spec.polygon, list(spec.holes))
    stock.subtract_disk(helix.footprint.centre.x, helix.footprint.centre.y, float(bore_radius))
    return stock


def _result(*operations: ToolpathOperation) -> ToolpathResult:
    return ToolpathResult(operations=list(operations), polyline=EMPTY_POLYLINE)


def test_engagement_replay_removes_exactly_the_bore_disk() -> None:
    spec = _spec()
    op = _helix_op()
    assert isinstance(op.geometry, Helix)
    stock = Stock(spec.polygon, list(spec.holes))
    engagement._replay_operation(stock, 0, op, spec.tool_radius, TEA_CAP_RAD, CUT_Z)
    assert stock.exactly_equals(_bored(spec, op.geometry))


def test_audit_toolpath_engagement_removes_the_helix_bore() -> None:
    # A loop inside the bore meets no material only if the bore was removed.
    spec = _spec()
    report = audit_toolpath_engagement(spec.polygon, _result(_helix_op(), _inner_loop()), TOOL_DIAMETER, TEA_CAP_RAD)
    loop = report.operations[1]
    assert loop.stations > 0
    assert loop.max_tea == 0.0


def test_the_helix_is_not_cap_sampled_and_records_as_a_plunge_does() -> None:
    spec = _spec()
    helix_record = audit_toolpath_engagement(spec.polygon, _result(_helix_op()), TOOL_DIAMETER, TEA_CAP_RAD).operations[0]
    plunge_record = audit_toolpath_engagement(spec.polygon, _result(_plunge()), TOOL_DIAMETER, TEA_CAP_RAD).operations[0]
    assert helix_record.operation is OperationType.HELIX
    assert (helix_record.max_tea, helix_record.cap_certified, helix_record.stations) == (plunge_record.max_tea, plunge_record.cap_certified, plunge_record.stations)
    assert helix_record.stations == 0


def test_engagement_replay_refuses_an_unknown_category(monkeypatch: pytest.MonkeyPatch) -> None:
    spec = _spec()
    monkeypatch.setattr(engagement, "classify_operation_replay", lambda op, cut_z: BOGUS_CATEGORY)
    stock = Stock(spec.polygon, list(spec.holes))
    with pytest.raises(UnhandledReplayCategoryError):
        engagement._replay_operation(stock, 0, _inner_loop(), spec.tool_radius, TEA_CAP_RAD, CUT_Z)


def test_depletion_bores_exactly_the_helix_disk() -> None:
    spec = _spec()
    op = _helix_op()
    assert isinstance(op.geometry, Helix)
    assert depletion._replay_kind(0, op, CUT_Z) is ReplayKind.HELIX
    assert replay_depletion(spec, _result(op)).exactly_equals(_bored(spec, op.geometry))


def test_depletion_reports_the_helix_bore_before_it_is_removed() -> None:
    spec = _spec()
    op = _helix_op()
    stock = Stock(spec.polygon, list(spec.holes))
    events = []
    for event in replay_events(spec, _result(op, _inner_loop()), stock):
        if isinstance(event, HelixEvent):
            # The stock is the pre-bore state while the event is held.
            assert event.stock_before.exactly_equals(Stock(spec.polygon, list(spec.holes)))
            assert (event.centre.x, event.centre.y, event.guide_radius.value) == (0.0, 0.0, FOOTPRINT_RADIUS)
        events.append(event)
    assert [type(event) for event in events] == [HelixEvent, CutMotion]
    helix_event = events[0]
    assert isinstance(helix_event, HelixEvent)
    assert helix_event.index == 0
    assert helix_event.into_material


def test_the_helix_is_not_a_cut_motion_of_the_depletion_replay() -> None:
    spec = _spec()
    stock = Stock(spec.polygon, list(spec.holes))
    assert [motion.index for motion in replay_cuts(spec, _result(_helix_op(), _inner_loop()), stock)] == [1]


def test_depletion_refuses_an_unknown_category(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(depletion, "classify_planar_replay", lambda operation, motion_z, cut_z: BOGUS_CATEGORY)
    with pytest.raises(UnhandledReplayCategoryError):
        depletion._replay_kind(0, _inner_loop(), CUT_Z)


def test_depletion_refuses_an_unknown_kind(monkeypatch: pytest.MonkeyPatch) -> None:
    spec = _spec()
    monkeypatch.setattr(depletion, "_replay_kind", lambda index, op, cut_z: BOGUS_CATEGORY)
    with pytest.raises(depletion.UnhandledReplayKindError):
        list(replay_events(spec, _result(_inner_loop()), Stock(spec.polygon, list(spec.holes))))


def test_a_helix_into_an_already_bored_stock_meets_no_material() -> None:
    spec = _spec()
    stock = Stock(spec.polygon, list(spec.holes))
    verdicts = [event.into_material for event in replay_events(spec, _result(_helix_op(), _helix_op()), stock) if isinstance(event, HelixEvent)]
    assert verdicts == [True, False]


def test_the_boundary_footprint_rho_equal_to_r_bores_the_disk_in_both_replays() -> None:
    spec = _spec()
    op = _helix_op(BOUNDARY_FOOTPRINT_RADIUS)
    assert isinstance(op.geometry, Helix)
    assert op.geometry.footprint.radius.value == spec.tool_radius
    bored = _bored(spec, op.geometry)
    stock = Stock(spec.polygon, list(spec.holes))
    engagement._replay_operation(stock, 0, op, spec.tool_radius, TEA_CAP_RAD, CUT_Z)
    assert stock.exactly_equals(bored)
    assert replay_depletion(spec, _result(op)).exactly_equals(bored)


def test_depletion_refuses_an_unsupported_geometry_by_name() -> None:
    bad = ToolpathOperation(Polygon([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]), OperationType.CUT, 0)
    with pytest.raises(UnreplayableOperationError):
        depletion._replay_kind(0, bad, CUT_Z)
