import math
import struct
from dataclasses import fields

import pytest

from compas_cgal.adaptive.errors import ArtifactIdentityError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import EdgeId
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckOwnerId
from compas_cgal.adaptive.operation import NeckTraversalOrientation
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import TraversalDecision
from compas_cgal.adaptive.operation import WidthClassId
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


def _plane() -> CutPlane:
    return CutPlane.build(CutZ.build(-2.0), ClearanceZ.build(5.0))


def _advance(*, terminal: bool = False) -> AdvanceTraversalDecision:
    return AdvanceTraversalDecision.build(
        component_id=ComponentId(b"component-a"),
        edge_id=EdgeId(b"edge-a"),
        branch_id=BranchId(b"branch-a"),
        cursor_before=CursorIdentity(b"cursor-0"),
        cursor_after=CursorIdentity(b"cursor-1"),
        makes_cursor_terminal=terminal,
    )


def _hold() -> HoldTraversalDecision:
    return HoldTraversalDecision.build(
        component_id=ComponentId(b"component-a"),
        edge_id=EdgeId(b"edge-a"),
        branch_id=BranchId(b"branch-a"),
        cursor=CursorIdentity(b"cursor-0"),
    )


def _full_cap() -> FullCapDecision:
    cap = EngagementCap.build(math.pi / 2.0)
    return FullCapDecision.build(user_cap=cap, effective_cap=cap)


def _neck_cap() -> NeckCapDecision:
    return NeckCapDecision.build(
        neck_evidence_digest=IdentityDigest(b"\x22" * 32),
        width_class_id=WidthClassId.build(0),
        passage_before=PassageState.UNVISITED,
        passage_after=PassageState.FIRST_PASS_COMPLETE,
        user_cap=EngagementCap.build(math.pi / 2.0),
        effective_cap=EngagementCap.build(math.pi / 3.0),
    )


def test_full_cap_requires_equal_user_and_effective_surrogate_bytes() -> None:
    user_cap = EngagementCap.build(math.pi / 2.0)

    with pytest.raises(ArtifactIdentityError, match="equal"):
        FullCapDecision.build(
            user_cap=user_cap,
            effective_cap=EngagementCap.build(math.pi / 3.0),
        )


def test_neck_cap_binds_evidence_class_passage_and_both_caps() -> None:
    baseline = _neck_cap()
    changed = NeckCapDecision.build(
        neck_evidence_digest=IdentityDigest(b"\x23" * 32),
        width_class_id=baseline.width_class_id,
        passage_before=baseline.passage_before,
        passage_after=baseline.passage_after,
        user_cap=EngagementCap.build(math.pi / 2.0),
        effective_cap=EngagementCap.build(math.pi / 3.0),
    )

    assert baseline.user_cap_bytes != baseline.effective_cap_bytes
    assert len(baseline.user_cap_bytes) == struct.calcsize(">d")
    assert len(baseline.effective_cap_bytes) == struct.calcsize(">d")
    assert baseline.canonical_bytes != changed.canonical_bytes


def test_neck_cap_rejects_nonadvancing_passage_or_excessive_effective_cap() -> None:
    with pytest.raises(ArtifactIdentityError, match="advance"):
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.UNVISITED,
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi / 3.0),
        )
    with pytest.raises(ArtifactIdentityError, match="exceeds"):
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.FIRST_PASS_COMPLETE,
            user_cap=EngagementCap.build(math.pi / 3.0),
            effective_cap=EngagementCap.build(math.pi / 2.0),
        )


def test_neck_cap_raw_constructor_rejects_untyped_passage_state() -> None:
    cap = EngagementCap.build(math.pi / 2.0).chord_ratio_bytes

    with pytest.raises(ArtifactIdentityError, match="typed passage"):
        NeckCapDecision(
            IdentityDigest(b"\x22" * 32),
            WidthClassId.build(0),
            "unvisited",  # type: ignore[arg-type]
            PassageState.FIRST_PASS_COMPLETE,
            cap,
            cap,
        )


def test_advance_traversal_binds_all_ids_cursor_transition_and_terminal_bit() -> None:
    baseline = _advance()
    terminal = _advance(terminal=True)
    changed_edge = AdvanceTraversalDecision.build(
        component_id=baseline.component_id,
        edge_id=EdgeId(b"edge-b"),
        branch_id=baseline.branch_id,
        cursor_before=baseline.cursor_before,
        cursor_after=baseline.cursor_after,
        makes_cursor_terminal=False,
    )

    assert len({baseline.canonical_bytes, terminal.canonical_bytes, changed_edge.canonical_bytes}) == 3


def test_advance_traversal_requires_real_progress_and_typed_terminal_state() -> None:
    with pytest.raises(ArtifactIdentityError, match="advance"):
        AdvanceTraversalDecision.build(
            component_id=ComponentId(b"component-a"),
            edge_id=EdgeId(b"edge-a"),
            branch_id=BranchId(b"branch-a"),
            cursor_before=CursorIdentity(b"cursor"),
            cursor_after=CursorIdentity(b"cursor"),
            makes_cursor_terminal=False,
        )
    with pytest.raises(ArtifactIdentityError, match="terminal"):
        AdvanceTraversalDecision.build(
            component_id=ComponentId(b"component-a"),
            edge_id=EdgeId(b"edge-a"),
            branch_id=BranchId(b"branch-a"),
            cursor_before=CursorIdentity(b"cursor-0"),
            cursor_after=CursorIdentity(b"cursor-1"),
            makes_cursor_terminal=1,  # type: ignore[arg-type]
        )


def test_hold_traversal_binds_equal_cursor_identities_and_cannot_be_terminal() -> None:
    hold = _hold()
    traversal: TraversalDecision = hold

    assert traversal.cursor_before == traversal.cursor_after
    assert traversal.makes_cursor_terminal is False
    assert hold.canonical_bytes != _advance().canonical_bytes


def test_approach_and_plunge_bind_the_same_xy_and_complete_cut_plane() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)
    approach = ApproachOperation.build(position=point, cut_plane=_plane())
    plunge = PlungeOperation.build(position=point, cut_plane=_plane())
    changed_plane = ApproachOperation.build(
        position=point,
        cut_plane=CutPlane.build(CutZ.build(-3.0), ClearanceZ.build(5.0)),
    )

    operations: tuple[CanonicalOperation, ...] = (approach, plunge)

    assert approach.canonical_bytes != plunge.canonical_bytes
    assert approach.canonical_bytes != changed_plane.canonical_bytes
    assert len(operations) == 2


def test_link_segment_binds_motion_neck_scope_cap_and_hold() -> None:
    motion = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )
    baseline = LinkSegmentOperation.build(
        motion=motion,
        cut_z=CutZ.build(-2.0),
        neck_scope=OrientedNeckScope.build(
            neck_owner_id=NeckOwnerId(b"neck-a"),
            orientation=NeckTraversalOrientation.FORWARD,
        ),
        effective_cap_decision=_neck_cap(),
        traversal_decision=_hold(),
    )
    changed_scope = LinkSegmentOperation.build(
        motion=motion,
        cut_z=baseline.cut_z,
        neck_scope=OrientedNeckScope.build(
            neck_owner_id=NeckOwnerId(b"neck-a"),
            orientation=NeckTraversalOrientation.REVERSE,
        ),
        effective_cap_decision=baseline.effective_cap_decision,
        traversal_decision=baseline.traversal_decision,
    )

    assert baseline.canonical_bytes != changed_scope.canonical_bytes


def test_cut_full_circle_binds_phase_orientation_scope_cap_and_advance() -> None:
    center = Point2[WorldXY].build(2.0, 3.0)
    baseline = CutFullCircleOperation.build(
        motion=ExactCircleMotion.build(center, Vector2[WorldXY].build(1.0, 0.0), False),
        cut_z=CutZ.build(-2.0),
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        traversal_decision=_advance(),
    )
    reversed_orientation = CutFullCircleOperation.build(
        motion=ExactCircleMotion.build(center, Vector2[WorldXY].build(1.0, 0.0), True),
        cut_z=baseline.cut_z,
        neck_scope=baseline.neck_scope,
        effective_cap_decision=baseline.effective_cap_decision,
        traversal_decision=baseline.traversal_decision,
    )

    assert baseline.canonical_bytes != reversed_orientation.canonical_bytes


def test_lateral_operation_scope_and_cap_decision_must_agree() -> None:
    motion = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )

    with pytest.raises(ArtifactIdentityError, match="neck scope"):
        LinkSegmentOperation.build(
            motion=motion,
            cut_z=CutZ.build(-2.0),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_neck_cap(),
            traversal_decision=_hold(),
        )


def test_link_and_circle_cannot_duplicate_or_exchange_the_cursor_advance() -> None:
    segment = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )
    circle = ExactCircleMotion.build(
        Point2[WorldXY].build(1.0, 0.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )

    with pytest.raises(ArtifactIdentityError, match="hold traversal"):
        LinkSegmentOperation(  # type: ignore[arg-type]
            segment,
            CutZ.build(-2.0),
            NoNeckScope.build(),
            _full_cap(),
            _advance(),
        )
    with pytest.raises(ArtifactIdentityError, match="own the one accepted"):
        CutFullCircleOperation(  # type: ignore[arg-type]
            circle,
            CutZ.build(-2.0),
            NoNeckScope.build(),
            _full_cap(),
            _hold(),
        )


def test_canonical_operations_do_not_embed_or_trust_witnesses() -> None:
    operation_types = (
        ApproachOperation,
        PlungeOperation,
        LinkSegmentOperation,
        CutFullCircleOperation,
    )

    for operation_type in operation_types:
        field_names = {field.name for field in fields(operation_type)}
        assert not any("witness" in field_name for field_name in field_names)
        assert "digest" not in field_names
