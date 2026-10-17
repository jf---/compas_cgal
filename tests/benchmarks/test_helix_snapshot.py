"""HELIX at the snapshot, evidence and identity boundaries (Task 1.2, consumer set C).

A helix is recorded faithfully (footprint centre/radius/entry, top/bottom Z,
turns, exact pitch), round-trips to the core `Helix`, and is refused loudly
when malformed. Figure 5 evidence names it explicitly instead of falling
through. Its operation identity is content-derived and distinct from a plunge
at the same point.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

import numpy as np
import pytest
from compas.geometry import Line

from benchmarks import held_path_evidence
from benchmarks.errors import InvalidHeldOperationSnapshotError
from benchmarks.errors import InvalidHeldPathEvidenceError
from benchmarks.errors import MutatedHeldToolpathError
from benchmarks.held_path_snapshot import HeldHelixSnapshot
from benchmarks.held_path_snapshot import InvalidHeldHelixSnapshotError
from benchmarks.held_path_snapshot import assert_toolpath_matches_snapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.quality_observations import HelixHasNoMotionKindError
from benchmarks.quality_observations import _snapshot_cut_height
from benchmarks.quality_observations import snapshot_motion_kind
from benchmarks.quality_observations import _snapshot_replay_category
from benchmarks.survey import MotionKind
from benchmarks.units import OperationIndex
from compas_cgal.engagement_audit.errors import InvalidAuditOperationError
from compas_cgal.engagement_audit.operation_identity import AuditHelixGeometryError
from compas_cgal.engagement_audit.operation_identity import HelixOperationSnapshot
from compas_cgal.engagement_audit.operation_identity import UnsupportedAuditHelixError
from compas_cgal.engagement_audit.operation_identity import canonical_toolpath_operation_bytes
from compas_cgal.engagement_audit.operation_identity import operation_digest
from compas_cgal.engagement_audit.operation_identity import snapshot_identity_operation
from compas_cgal.engagement_audit.operation_identity import snapshot_toolpath_operation
from compas_cgal.helix import Helix
from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.helix_units import RampAngleDegrees
from compas_cgal.incident_loop import INCIDENCE_CHART_DENOMINATOR
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import helix_operation
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import Direction3
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY
from compas_cgal.units import WorldXYZ

BOTTOM_Z = 0.0
TOP_Z = 4.0
RAMP_DEG = 2.0
CENTRE_XY = (1.0, 2.0)
FOOTPRINT_RADIUS = 1.0


def _helix(*, radius: float = FOOTPRINT_RADIUS, top_z: float = TOP_Z, ramp_deg: float = RAMP_DEG, clockwise: bool = False, bottom_z: float = BOTTOM_Z) -> Helix:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(CENTRE_XY[0]), Millimetre(CENTRE_XY[1])), GuideRadius(Millimetre(radius)), Radian(0.0))
    plane = CutPlane.build(CutZ.build(bottom_z), ClearanceZ.build(top_z))
    return Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(ramp_deg)), clockwise=clockwise)


def _helix_operation(helix: Helix, *, operation: OperationType = OperationType.HELIX, clockwise: bool | None = None) -> ToolpathOperation:
    entry = helix.entry_xy
    return ToolpathOperation(
        geometry=helix,
        operation=operation,
        path_index=0,
        clockwise=helix.clockwise if clockwise is None else clockwise,
        entry_xy=(entry.x, entry.y),
    )


def _result(*operations: ToolpathOperation) -> ToolpathResult:
    return ToolpathResult(operations=list(operations), polyline=np.zeros((0, 3)))


def _plunge_at_entry(helix: Helix) -> ToolpathOperation:
    entry = helix.entry_xy
    line = Line([entry.x, entry.y, TOP_Z], [entry.x, entry.y, BOTTOM_Z])
    return ToolpathOperation(geometry=line, operation=OperationType.PLUNGE, path_index=0)


def _snapshot_fields(helix: Helix) -> dict[str, Any]:
    footprint = helix.footprint
    return {
        "ordinal": OperationIndex(0),
        "operation": OperationType.HELIX,
        "path_index": 0,
        "clockwise": helix.clockwise,
        "centre": footprint.centre,
        "entry": footprint.entry,
        "radius": footprint.radius.value,
        "phase": footprint.phase,
        "chart": footprint.chart,
        "chart_numerator": footprint.chart_numerator,
        "top_z": helix.top_z.value,
        "bottom_z": helix.bottom_z.value,
        "turns": helix.turns,
        "pitch": helix.pitch,
        "start_tangent": None,
        "end_tangent": None,
    }


# --- 1. snapshot -----------------------------------------------------------------


def test_a_helix_operation_snapshots_faithfully() -> None:
    helix = _helix()
    (record,) = snapshot_toolpath(_result(_helix_operation(helix)))
    assert type(record) is HeldHelixSnapshot
    assert record.operation is OperationType.HELIX
    assert record.centre == helix.footprint.centre
    assert record.entry == helix.entry_xy
    assert record.radius == helix.footprint.radius.value
    assert record.top_z == TOP_Z
    assert record.bottom_z == BOTTOM_Z
    assert record.turns == helix.turns == 19
    assert type(record.pitch) is Fraction
    assert record.pitch == Fraction(4) / 19 == helix.pitch
    assert record.clockwise is False


@pytest.mark.parametrize("clockwise", [False, True])
def test_a_helix_snapshot_round_trips_to_the_core_helix(clockwise: bool) -> None:
    helix = _helix(clockwise=clockwise)
    result = _result(_helix_operation(helix))
    snapshot = snapshot_toolpath(result)
    (record,) = snapshot
    assert isinstance(record, HeldHelixSnapshot)
    assert record.helix() == helix
    assert HeldHelixSnapshot.build(**_snapshot_fields(helix)) == record
    assert_toolpath_matches_snapshot(result, snapshot)


def test_a_changed_helix_is_a_mutated_toolpath() -> None:
    result = _result(_helix_operation(_helix()))
    snapshot = snapshot_toolpath(result)
    steeper = _helix(ramp_deg=3.0)
    result.operations[0] = _helix_operation(steeper)
    with pytest.raises(MutatedHeldToolpathError):
        assert_toolpath_matches_snapshot(result, snapshot)


def _malformed(field: str, value: object) -> dict[str, Any]:
    fields = _snapshot_fields(_helix())
    fields[field] = value
    return fields


MALFORMED = {
    "pitch-off-by-a-turn": ("pitch", ExactMillimetre(Fraction(4, 18))),
    "pitch-as-float": ("pitch", 4.0 / 19.0),
    "zero-turns": ("turns", 0),
    "bool-turns": ("turns", True),
    "flat-descent": ("top_z", Millimetre(BOTTOM_Z)),
    "entry-off-circle": ("entry", Point2[WorldXY](Millimetre(2.0), Millimetre(2.5))),
    "nonpositive-radius": ("radius", Millimetre(0.0)),
    "chart-out-of-range": ("chart", 4),
    "role-is-plunge": ("operation", OperationType.PLUNGE),
    "untyped-centre": ("centre", (1.0, 2.0)),
    "numerator-above-denominator": ("chart_numerator", INCIDENCE_CHART_DENOMINATOR + 1),
    "chart-contradicts-entry": ("chart", 1),
    "numerator-contradicts-entry": ("chart_numerator", 1),
    "phase-contradicts-entry": ("phase", Radian(0.5)),
}


@pytest.mark.parametrize("case", sorted(MALFORMED))
def test_a_malformed_helix_snapshot_is_refused_by_name(case: str) -> None:
    field, value = MALFORMED[case]
    with pytest.raises(InvalidHeldHelixSnapshotError):
        HeldHelixSnapshot.build(**_malformed(field, value))


def test_the_named_helix_error_is_a_snapshot_error() -> None:
    assert issubclass(InvalidHeldHelixSnapshotError, InvalidHeldOperationSnapshotError)


@pytest.mark.parametrize("clockwise", [False, True])
def test_a_helix_tangent_is_the_unit_xy_travel_direction_with_z_zero(clockwise: bool) -> None:
    """The emitter convention: the slope is exact in ``pitch``, never in the tangent."""
    helix = _helix(clockwise=clockwise)
    fields = _snapshot_fields(helix)
    record = HeldHelixSnapshot.build(**fields)
    tangent = record.analytic_tangent()
    assert tangent.z == 0.0
    assert tangent == Direction3[WorldXYZ].build(0.0, -1.0 if clockwise else 1.0, 0.0)  # entry at phase 0
    fields["start_tangent"] = tangent
    fields["end_tangent"] = tangent
    HeldHelixSnapshot.build(**fields)
    for wrong in (
        Direction3[WorldXYZ].build(-tangent.x, -tangent.y, 0.0),
        Direction3[WorldXYZ].build(1.0, 0.0, 0.0),
    ):
        fields["start_tangent"] = wrong
        with pytest.raises(InvalidHeldHelixSnapshotError):
            HeldHelixSnapshot.build(**fields)


@pytest.mark.parametrize("clockwise", [False, True])
@pytest.mark.parametrize("phase", [0.0, 0.7, 2.5, -1.9])
def test_an_emitted_helix_operation_snapshots_and_round_trips(clockwise: bool, phase: float) -> None:
    footprint = IncidentLoop.build(Point2[WorldXY](Millimetre(CENTRE_XY[0]), Millimetre(CENTRE_XY[1])), GuideRadius(Millimetre(FOOTPRINT_RADIUS)), Radian(phase))
    plane = CutPlane.build(CutZ.build(BOTTOM_Z), ClearanceZ.build(TOP_Z))
    helix = Helix.build(footprint, plane, HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees.build(RAMP_DEG)), clockwise=clockwise)
    result = _result(helix_operation(helix, 0))
    snapshot = snapshot_toolpath(result)
    (record,) = snapshot
    assert isinstance(record, HeldHelixSnapshot)
    assert record.start_tangent is not None and record.end_tangent is not None
    assert record.helix() == helix
    assert_toolpath_matches_snapshot(result, snapshot)


def test_helix_geometry_under_another_role_is_refused() -> None:
    with pytest.raises(InvalidHeldHelixSnapshotError):
        snapshot_toolpath(_result(_helix_operation(_helix(), operation=OperationType.CUT)))


def test_a_helix_role_with_line_geometry_is_refused() -> None:
    plunge = _plunge_at_entry(_helix())
    plunge.operation = OperationType.HELIX
    with pytest.raises(InvalidHeldHelixSnapshotError):
        snapshot_toolpath(_result(plunge))


def test_a_helix_operation_must_carry_its_exact_entry_and_turn_sense() -> None:
    helix = _helix()
    missing_entry = _helix_operation(helix)
    missing_entry.entry_xy = None
    shifted_entry = _helix_operation(helix)
    shifted_entry.entry_xy = (helix.entry_xy.x, float(np.nextafter(helix.entry_xy.y, 3.0)))
    flipped = _helix_operation(helix, clockwise=True)
    for operation in (missing_entry, shifted_entry, flipped):
        with pytest.raises(InvalidHeldHelixSnapshotError):
            snapshot_toolpath(_result(operation))


# --- 2. evidence -----------------------------------------------------------------


def test_held_path_evidence_kinds_a_helix_explicitly() -> None:
    """Evidence and quality observations share one kind dispatch and one named helix error."""
    assert vars(held_path_evidence)["snapshot_motion_kind"] is snapshot_motion_kind
    (record,) = snapshot_toolpath(_result(_helix_operation(_helix())))
    with pytest.raises(HelixHasNoMotionKindError):
        held_path_evidence._validate_snapshot((record,))
    (plunge,) = snapshot_toolpath(_result(_plunge_at_entry(_helix())))
    assert held_path_evidence._validate_snapshot((plunge,)) == 1


# --- 3. identity -----------------------------------------------------------------


def test_helix_identity_is_content_derived_and_stable() -> None:
    first = canonical_toolpath_operation_bytes(_helix_operation(_helix()))
    again = canonical_toolpath_operation_bytes(_helix_operation(_helix()))
    assert first == again
    assert operation_digest(first) == operation_digest(again)
    assert type(snapshot_identity_operation(_helix_operation(_helix()))) is HelixOperationSnapshot


@pytest.mark.parametrize(
    "variant",
    [
        {"ramp_deg": 3.0},
        {"top_z": 5.0},
        {"radius": 0.5},
        {"clockwise": True},
    ],
)
def test_helix_identity_moves_with_every_defining_value(variant: dict[str, Any]) -> None:
    base = canonical_toolpath_operation_bytes(_helix_operation(_helix()))
    assert canonical_toolpath_operation_bytes(_helix_operation(_helix(**variant))) != base


def test_helix_identity_is_distinct_from_a_plunge_at_the_same_point() -> None:
    helix = _helix()
    helix_bytes = canonical_toolpath_operation_bytes(_helix_operation(helix))
    plunge_bytes = canonical_toolpath_operation_bytes(_plunge_at_entry(helix))
    assert helix_bytes != plunge_bytes
    assert operation_digest(helix_bytes) != operation_digest(plunge_bytes)


def test_the_native_audit_ingress_refuses_a_helix_by_name() -> None:
    with pytest.raises(UnsupportedAuditHelixError):
        snapshot_toolpath_operation(_helix_operation(_helix()))
    assert issubclass(UnsupportedAuditHelixError, InvalidAuditOperationError)


def test_identity_refuses_a_helix_role_with_line_geometry() -> None:
    plunge = _plunge_at_entry(_helix())
    plunge.operation = OperationType.HELIX
    with pytest.raises(AuditHelixGeometryError):
        snapshot_identity_operation(plunge)


def test_identity_refuses_helix_geometry_under_a_cut_role() -> None:
    with pytest.raises(AuditHelixGeometryError):
        snapshot_identity_operation(_helix_operation(_helix(), operation=OperationType.CUT))


# --- 4. quality observations (Q-1..Q-3, round 1) ------------------------------------

SUNK_BOTTOM_Z = -1.0  # A cut plane below z = 0, so a skipped anchor (default 0.0) is visible.


def _sunk_helix_record() -> HeldHelixSnapshot:
    (record,) = snapshot_toolpath(_result(_helix_operation(_helix(bottom_z=SUNK_BOTTOM_Z))))
    assert isinstance(record, HeldHelixSnapshot)
    return record


def test_snapshot_cut_height_anchors_on_the_helix_bottom() -> None:
    assert _snapshot_cut_height((_sunk_helix_record(),)) == SUNK_BOTTOM_Z


def test_snapshot_replay_category_of_a_helix_is_helix() -> None:
    assert _snapshot_replay_category(_sunk_helix_record(), SUNK_BOTTOM_Z) == "helix"


def test_an_off_plane_helix_is_unreplayable_evidence() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        _snapshot_replay_category(_sunk_helix_record(), BOTTOM_Z)


def test_quality_observations_kinds_a_helix_explicitly() -> None:
    with pytest.raises(HelixHasNoMotionKindError):
        snapshot_motion_kind(_sunk_helix_record())
    assert issubclass(HelixHasNoMotionKindError, InvalidHeldPathEvidenceError)
    (plunge,) = snapshot_toolpath(_result(_plunge_at_entry(_helix())))
    assert snapshot_motion_kind(plunge) is MotionKind.LINE
