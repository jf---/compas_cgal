import math
import struct
from dataclasses import dataclass
from dataclasses import fields

import pytest

from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import InvalidAdvanceSegmentOperationError
from compas_cgal.adaptive.errors import InvalidEffectiveCapDecisionError
from compas_cgal.adaptive.errors import InvalidOperationIdentityError
from compas_cgal.adaptive.errors import InvalidTraversalDecisionError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
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
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY


@dataclass(frozen=True)
class SemanticApproach(ApproachOperation):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticPlunge(PlungeOperation):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticLink(LinkSegmentOperation):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticCircle(CutFullCircleOperation):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticPoint(Point2[WorldXY]):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticVector(Vector2[WorldXY]):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticFullCap(FullCapDecision):
    provenance_bit: bool


@dataclass(frozen=True)
class SemanticHold(HoldTraversalDecision):
    provenance_bit: bool


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

    with pytest.raises(InvalidEffectiveCapDecisionError, match="equal"):
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
    with pytest.raises(InvalidEffectiveCapDecisionError, match="advance"):
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.UNVISITED,
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi / 3.0),
        )
    with pytest.raises(InvalidEffectiveCapDecisionError, match="exceeds"):
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

    with pytest.raises(InvalidEffectiveCapDecisionError, match="typed passage"):
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
    with pytest.raises(InvalidTraversalDecisionError, match="advance"):
        AdvanceTraversalDecision.build(
            component_id=ComponentId(b"component-a"),
            edge_id=EdgeId(b"edge-a"),
            branch_id=BranchId(b"branch-a"),
            cursor_before=CursorIdentity(b"cursor"),
            cursor_after=CursorIdentity(b"cursor"),
            makes_cursor_terminal=False,
        )
    with pytest.raises(InvalidTraversalDecisionError, match="terminal"):
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


def test_approach_binds_only_target_xy_and_clearance_z() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)
    approach = ApproachOperation.build(position=point, clearance_z=ClearanceZ.build(5.0))
    changed_clearance = ApproachOperation.build(
        position=point,
        clearance_z=ClearanceZ.build(6.0),
    )
    unrelated_cut_z_a = CutZ.build(-2.0)
    unrelated_cut_z_b = CutZ.build(-3.0)

    assert unrelated_cut_z_a != unrelated_cut_z_b
    assert {field.name for field in fields(ApproachOperation)} == {"position", "clearance_z"}
    assert approach.canonical_bytes != changed_clearance.canonical_bytes


def test_plunge_binds_same_xy_clearance_z_and_cut_z() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)
    approach = ApproachOperation.build(position=point, clearance_z=ClearanceZ.build(5.0))
    plunge = PlungeOperation.build(
        position=point,
        clearance_z=ClearanceZ.build(5.0),
        cut_z=CutZ.build(-2.0),
    )
    changed_cut = PlungeOperation.build(
        position=point,
        clearance_z=plunge.clearance_z,
        cut_z=CutZ.build(-3.0),
    )
    operations: tuple[CanonicalOperation, ...] = (approach, plunge)

    assert approach.position == plunge.position
    assert approach.clearance_z == plunge.clearance_z
    assert plunge.canonical_bytes != changed_cut.canonical_bytes
    assert approach.canonical_bytes != plunge.canonical_bytes
    assert len(operations) == 2


@pytest.mark.parametrize(
    ("clearance_z", "cut_z"),
    (
        (ClearanceZ.build(-2.0), CutZ.build(-2.0)),
        (ClearanceZ.build(-3.0), CutZ.build(-2.0)),
    ),
)
def test_plunge_requires_clearance_strictly_above_cut(
    clearance_z: ClearanceZ,
    cut_z: CutZ,
) -> None:
    with pytest.raises(InvalidOperationIdentityError, match="strictly above"):
        PlungeOperation.build(
            position=Point2[WorldXY].build(1.0, 2.0),
            clearance_z=clearance_z,
            cut_z=cut_z,
        )


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


def test_advance_segment_is_a_distinct_canonical_traversal_operation() -> None:
    """Separate link-only progress from a hold link awaiting one circle.

    A mode bit on `LinkSegmentOperation` would make consumers interpret one
    record through two chronologies. The advancing variant instead owns its
    cursor transition and changes identity when its physical endpoint changes.
    """
    segment = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )
    scope = OrientedNeckScope.build(
        neck_owner_id=NeckOwnerId(b"neck-a"),
        orientation=NeckTraversalOrientation.FORWARD,
    )
    operation = AdvanceSegmentOperation.build(
        motion=segment,
        cut_z=CutZ.build(-2.0),
        neck_scope=scope,
        effective_cap_decision=_neck_cap(),
        traversal_decision=_advance(),
    )
    hold_link = LinkSegmentOperation.build(
        motion=segment,
        cut_z=operation.cut_z,
        neck_scope=operation.neck_scope,
        effective_cap_decision=operation.effective_cap_decision,
        traversal_decision=_hold(),
    )
    variants = (
        operation,
        AdvanceSegmentOperation.build(
            motion=ExactSegmentMotion.build(
                segment.start,
                Point2[WorldXY].build(2.0, 0.0),
            ),
            cut_z=operation.cut_z,
            neck_scope=operation.neck_scope,
            effective_cap_decision=operation.effective_cap_decision,
            traversal_decision=operation.traversal_decision,
        ),
        AdvanceSegmentOperation.build(
            motion=operation.motion,
            cut_z=CutZ.build(-3.0),
            neck_scope=operation.neck_scope,
            effective_cap_decision=operation.effective_cap_decision,
            traversal_decision=operation.traversal_decision,
        ),
        AdvanceSegmentOperation.build(
            motion=operation.motion,
            cut_z=operation.cut_z,
            neck_scope=OrientedNeckScope.build(
                neck_owner_id=NeckOwnerId(b"neck-a"),
                orientation=NeckTraversalOrientation.REVERSE,
            ),
            effective_cap_decision=operation.effective_cap_decision,
            traversal_decision=operation.traversal_decision,
        ),
        AdvanceSegmentOperation.build(
            motion=operation.motion,
            cut_z=operation.cut_z,
            neck_scope=operation.neck_scope,
            effective_cap_decision=NeckCapDecision.build(
                neck_evidence_digest=IdentityDigest(b"\x23" * 32),
                width_class_id=WidthClassId.build(0),
                passage_before=PassageState.UNVISITED,
                passage_after=PassageState.FIRST_PASS_COMPLETE,
                user_cap=EngagementCap.build(math.pi / 2.0),
                effective_cap=EngagementCap.build(math.pi / 3.0),
            ),
            traversal_decision=operation.traversal_decision,
        ),
        AdvanceSegmentOperation.build(
            motion=operation.motion,
            cut_z=operation.cut_z,
            neck_scope=operation.neck_scope,
            effective_cap_decision=operation.effective_cap_decision,
            traversal_decision=_advance(terminal=True),
        ),
    )

    assert require_canonical_record(operation.canonical_bytes) == (operation.canonical_bytes)
    assert b"advance-segment-v1" in operation.canonical_bytes
    assert operation != hold_link
    assert len({variant.canonical_bytes for variant in variants}) == len(variants)


def test_advance_segment_fails_through_its_named_error_boundary() -> None:
    """Keep every malformed advancing record distinguishable from old links.

    Evaluation and replay need to report corruption of the new one-motion
    chronology directly. Letting nested generic operation errors escape would
    erase which closed-union variant failed authentication.
    """
    segment = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )
    cap = _full_cap()
    advance = _advance()
    with pytest.raises(InvalidAdvanceSegmentOperationError, match="advance traversal"):
        AdvanceSegmentOperation(  # type: ignore[arg-type]
            segment,
            CutZ.build(-2.0),
            NoNeckScope.build(),
            cap,
            _hold(),
        )
    with pytest.raises(InvalidAdvanceSegmentOperationError, match="neck scope"):
        AdvanceSegmentOperation.build(
            motion=segment,
            cut_z=CutZ.build(-2.0),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_neck_cap(),
            traversal_decision=advance,
        )
    with pytest.raises(InvalidAdvanceSegmentOperationError, match="ExactSegmentMotion"):
        AdvanceSegmentOperation(  # type: ignore[arg-type]
            object(),
            CutZ.build(-2.0),
            NoNeckScope.build(),
            cap,
            advance,
        )
    with pytest.raises(InvalidAdvanceSegmentOperationError, match="exact Point2"):
        AdvanceSegmentOperation.build(
            motion=ExactSegmentMotion.build(
                SemanticPoint(0.0, 0.0, True),
                segment.end,
            ),
            cut_z=CutZ.build(-2.0),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=cap,
            traversal_decision=advance,
        )

    @dataclass(frozen=True)
    class SemanticAdvance(AdvanceSegmentOperation):
        provenance_bit: bool

    semantic = SemanticAdvance(
        segment,
        CutZ.build(-2.0),
        NoNeckScope.build(),
        cap,
        advance,
        True,
    )
    with pytest.raises(
        InvalidAdvanceSegmentOperationError,
        match="exact AdvanceSegmentOperation",
    ):
        _ = semantic.canonical_bytes


def test_cut_full_circle_binds_phase_orientation_scope_cap_and_advance() -> None:
    center = Point2[WorldXY].build(2.0, 3.0)
    baseline = CutFullCircleOperation.build(
        motion=ExactCircleMotion.build(center, Vector2[WorldXY].build(1.0, 0.0), False),
        cut_z=CutZ.build(-2.0),
        material_side=MaterialSide.OUTSIDE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        traversal_decision=_advance(),
    )
    reversed_orientation = CutFullCircleOperation.build(
        motion=ExactCircleMotion.build(center, Vector2[WorldXY].build(1.0, 0.0), True),
        cut_z=baseline.cut_z,
        material_side=baseline.material_side,
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

    with pytest.raises(InvalidOperationIdentityError, match="neck scope"):
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

    with pytest.raises(InvalidOperationIdentityError, match="hold traversal"):
        LinkSegmentOperation(  # type: ignore[arg-type]
            segment,
            CutZ.build(-2.0),
            NoNeckScope.build(),
            _full_cap(),
            _advance(),
        )
    with pytest.raises(InvalidOperationIdentityError, match="own the one accepted"):
        CutFullCircleOperation(  # type: ignore[arg-type]
            circle,
            CutZ.build(-2.0),
            MaterialSide.OUTSIDE,
            NoNeckScope.build(),
            _full_cap(),
            _hold(),
        )


def test_operation_and_nested_decision_subclasses_cannot_drop_semantics() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)
    approaches = (
        SemanticApproach(point, ClearanceZ.build(5.0), False),
        SemanticApproach(point, ClearanceZ.build(5.0), True),
    )
    plunge = SemanticPlunge(point, ClearanceZ.build(5.0), CutZ.build(-2.0), True)
    cap = _full_cap()
    semantic_cap = SemanticFullCap(cap.user_cap_bytes, cap.effective_cap_bytes, True)
    hold = _hold()
    semantic_hold = SemanticHold(
        hold.component_id,
        hold.edge_id,
        hold.branch_id,
        hold.cursor_before,
        hold.cursor_after,
        hold.makes_cursor_terminal,
        True,
    )
    link_motion = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )
    link = SemanticLink(
        link_motion,
        CutZ.build(-2.0),
        NoNeckScope.build(),
        cap,
        hold,
        True,
    )
    circle = SemanticCircle(
        ExactCircleMotion.build(
            Point2[WorldXY].build(1.0, 0.0),
            Vector2[WorldXY].build(1.0, 0.0),
            False,
        ),
        CutZ.build(-2.0),
        MaterialSide.OUTSIDE,
        NoNeckScope.build(),
        cap,
        _advance(),
        True,
    )

    assert approaches[0] != approaches[1]
    for approach in approaches:
        with pytest.raises(InvalidOperationIdentityError, match="exact ApproachOperation"):
            _ = approach.canonical_bytes
    with pytest.raises(InvalidOperationIdentityError, match="exact PlungeOperation"):
        _ = plunge.canonical_bytes
    with pytest.raises(InvalidOperationIdentityError, match="exact LinkSegmentOperation"):
        _ = link.canonical_bytes
    with pytest.raises(InvalidOperationIdentityError, match="exact CutFullCircleOperation"):
        _ = circle.canonical_bytes
    with pytest.raises(InvalidEffectiveCapDecisionError, match="exact FullCapDecision"):
        _ = semantic_cap.canonical_bytes
    with pytest.raises(InvalidTraversalDecisionError, match="exact HoldTraversalDecision"):
        _ = semantic_hold.canonical_bytes


def test_lateral_operation_rejects_subclassed_nested_decisions() -> None:
    cap = _full_cap()
    semantic_cap = SemanticFullCap(cap.user_cap_bytes, cap.effective_cap_bytes, True)
    hold = _hold()
    semantic_hold = SemanticHold(
        hold.component_id,
        hold.edge_id,
        hold.branch_id,
        hold.cursor_before,
        hold.cursor_after,
        hold.makes_cursor_terminal,
        True,
    )
    motion = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )

    with pytest.raises(InvalidOperationIdentityError, match="exact effective-cap"):
        LinkSegmentOperation.build(
            motion=motion,
            cut_z=CutZ.build(-2.0),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=semantic_cap,
            traversal_decision=hold,
        )
    with pytest.raises(InvalidOperationIdentityError, match="exact HoldTraversalDecision"):
        LinkSegmentOperation.build(
            motion=motion,
            cut_z=CutZ.build(-2.0),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=cap,
            traversal_decision=semantic_hold,
        )


def test_lateral_raw_constructor_rejects_motion_with_subclassed_geometry() -> None:
    semantic_point = SemanticPoint(0.0, 0.0, True)
    segment = ExactSegmentMotion.build(
        semantic_point,
        Point2[WorldXY].build(1.0, 0.0),
    )
    semantic_vector = SemanticVector(1.0, 0.0, True)
    circle = ExactCircleMotion.build(
        Point2[WorldXY].build(1.0, 0.0),
        semantic_vector,
        False,
    )

    with pytest.raises(InvalidOperationIdentityError, match="exact Point2"):
        LinkSegmentOperation.build(
            motion=segment,
            cut_z=CutZ.build(-2.0),
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            traversal_decision=_hold(),
        )
    with pytest.raises(InvalidOperationIdentityError, match="exact Vector2"):
        CutFullCircleOperation.build(
            motion=circle,
            cut_z=CutZ.build(-2.0),
            material_side=MaterialSide.OUTSIDE,
            neck_scope=NoNeckScope.build(),
            effective_cap_decision=_full_cap(),
            traversal_decision=_advance(),
        )


def test_cap_decision_mutation_matrix_binds_every_field() -> None:
    full_cap = _full_cap()
    changed_full_cap_value = EngagementCap.build(math.pi / 3.0)
    full_variants = (
        full_cap,
        FullCapDecision.build(
            user_cap=changed_full_cap_value,
            effective_cap=changed_full_cap_value,
        ),
    )
    neck_variants = (
        _neck_cap(),
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x23" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.FIRST_PASS_COMPLETE,
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi / 3.0),
        ),
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(1),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.FIRST_PASS_COMPLETE,
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi / 3.0),
        ),
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.FIRST_PASS_COMPLETE,
            passage_after=PassageState.SECOND_PASS_COMPLETE,
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi / 3.0),
        ),
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.FIRST_PASS_COMPLETE,
            user_cap=EngagementCap.build(3.0 * math.pi / 4.0),
            effective_cap=EngagementCap.build(math.pi / 3.0),
        ),
        NeckCapDecision.build(
            neck_evidence_digest=IdentityDigest(b"\x22" * 32),
            width_class_id=WidthClassId.build(0),
            passage_before=PassageState.UNVISITED,
            passage_after=PassageState.FIRST_PASS_COMPLETE,
            user_cap=EngagementCap.build(math.pi / 2.0),
            effective_cap=EngagementCap.build(math.pi / 4.0),
        ),
    )

    assert len({decision.canonical_bytes for decision in full_variants}) == len(full_variants)
    assert len({decision.canonical_bytes for decision in neck_variants}) == len(neck_variants)


def test_full_circle_binds_the_material_side_that_caused_its_orientation() -> None:
    motion = ExactCircleMotion.build(
        Point2[WorldXY].build(1.0, 0.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )
    outside = CutFullCircleOperation.build(
        motion=motion,
        cut_z=CutZ.build(-2.0),
        material_side=MaterialSide.OUTSIDE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        traversal_decision=_advance(),
    )
    inside = CutFullCircleOperation.build(
        motion=motion,
        cut_z=CutZ.build(-2.0),
        material_side=MaterialSide.INSIDE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        traversal_decision=_advance(),
    )

    assert outside.material_side is MaterialSide.OUTSIDE
    assert outside.canonical_bytes != inside.canonical_bytes


def test_traversal_decision_mutation_matrix_binds_every_field() -> None:
    hold = _hold()
    hold_variants = (
        hold,
        HoldTraversalDecision.build(
            component_id=ComponentId(b"component-b"),
            edge_id=hold.edge_id,
            branch_id=hold.branch_id,
            cursor=hold.cursor_before,
        ),
        HoldTraversalDecision.build(
            component_id=hold.component_id,
            edge_id=EdgeId(b"edge-b"),
            branch_id=hold.branch_id,
            cursor=hold.cursor_before,
        ),
        HoldTraversalDecision.build(
            component_id=hold.component_id,
            edge_id=hold.edge_id,
            branch_id=BranchId(b"branch-b"),
            cursor=hold.cursor_before,
        ),
        HoldTraversalDecision.build(
            component_id=hold.component_id,
            edge_id=hold.edge_id,
            branch_id=hold.branch_id,
            cursor=CursorIdentity(b"cursor-other"),
        ),
    )
    advance = _advance()
    advance_variants = (
        advance,
        AdvanceTraversalDecision.build(
            component_id=ComponentId(b"component-b"),
            edge_id=advance.edge_id,
            branch_id=advance.branch_id,
            cursor_before=advance.cursor_before,
            cursor_after=advance.cursor_after,
            makes_cursor_terminal=advance.makes_cursor_terminal,
        ),
        AdvanceTraversalDecision.build(
            component_id=advance.component_id,
            edge_id=EdgeId(b"edge-b"),
            branch_id=advance.branch_id,
            cursor_before=advance.cursor_before,
            cursor_after=advance.cursor_after,
            makes_cursor_terminal=advance.makes_cursor_terminal,
        ),
        AdvanceTraversalDecision.build(
            component_id=advance.component_id,
            edge_id=advance.edge_id,
            branch_id=BranchId(b"branch-b"),
            cursor_before=advance.cursor_before,
            cursor_after=advance.cursor_after,
            makes_cursor_terminal=advance.makes_cursor_terminal,
        ),
        AdvanceTraversalDecision.build(
            component_id=advance.component_id,
            edge_id=advance.edge_id,
            branch_id=advance.branch_id,
            cursor_before=CursorIdentity(b"cursor-before-other"),
            cursor_after=advance.cursor_after,
            makes_cursor_terminal=advance.makes_cursor_terminal,
        ),
        AdvanceTraversalDecision.build(
            component_id=advance.component_id,
            edge_id=advance.edge_id,
            branch_id=advance.branch_id,
            cursor_before=advance.cursor_before,
            cursor_after=CursorIdentity(b"cursor-after-other"),
            makes_cursor_terminal=advance.makes_cursor_terminal,
        ),
        _advance(terminal=True),
    )

    assert len({decision.canonical_bytes for decision in hold_variants}) == len(hold_variants)
    assert len({decision.canonical_bytes for decision in advance_variants}) == len(advance_variants)


def test_operation_mutation_matrix_binds_every_field() -> None:
    point = Point2[WorldXY].build(1.0, 2.0)
    approach = ApproachOperation.build(position=point, clearance_z=ClearanceZ.build(5.0))
    approach_variants = (
        approach,
        ApproachOperation.build(position=Point2[WorldXY].build(2.0, 2.0), clearance_z=approach.clearance_z),
        ApproachOperation.build(position=point, clearance_z=ClearanceZ.build(6.0)),
    )
    plunge = PlungeOperation.build(
        position=point,
        clearance_z=ClearanceZ.build(5.0),
        cut_z=CutZ.build(-2.0),
    )
    plunge_variants = (
        plunge,
        PlungeOperation.build(
            position=Point2[WorldXY].build(2.0, 2.0),
            clearance_z=plunge.clearance_z,
            cut_z=plunge.cut_z,
        ),
        PlungeOperation.build(position=point, clearance_z=ClearanceZ.build(6.0), cut_z=plunge.cut_z),
        PlungeOperation.build(position=point, clearance_z=plunge.clearance_z, cut_z=CutZ.build(-3.0)),
    )
    link_motion = ExactSegmentMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Point2[WorldXY].build(1.0, 0.0),
    )
    link = LinkSegmentOperation.build(
        motion=link_motion,
        cut_z=CutZ.build(-2.0),
        neck_scope=OrientedNeckScope.build(
            neck_owner_id=NeckOwnerId(b"neck-a"),
            orientation=NeckTraversalOrientation.FORWARD,
        ),
        effective_cap_decision=_neck_cap(),
        traversal_decision=_hold(),
    )
    link_variants = (
        link,
        LinkSegmentOperation.build(
            motion=ExactSegmentMotion.build(
                Point2[WorldXY].build(0.0, 0.0),
                Point2[WorldXY].build(2.0, 0.0),
            ),
            cut_z=link.cut_z,
            neck_scope=link.neck_scope,
            effective_cap_decision=link.effective_cap_decision,
            traversal_decision=link.traversal_decision,
        ),
        LinkSegmentOperation.build(
            motion=link.motion,
            cut_z=CutZ.build(-3.0),
            neck_scope=link.neck_scope,
            effective_cap_decision=link.effective_cap_decision,
            traversal_decision=link.traversal_decision,
        ),
        LinkSegmentOperation.build(
            motion=link.motion,
            cut_z=link.cut_z,
            neck_scope=OrientedNeckScope.build(
                neck_owner_id=NeckOwnerId(b"neck-b"),
                orientation=NeckTraversalOrientation.FORWARD,
            ),
            effective_cap_decision=link.effective_cap_decision,
            traversal_decision=link.traversal_decision,
        ),
        LinkSegmentOperation.build(
            motion=link.motion,
            cut_z=link.cut_z,
            neck_scope=link.neck_scope,
            effective_cap_decision=NeckCapDecision.build(
                neck_evidence_digest=IdentityDigest(b"\x23" * 32),
                width_class_id=WidthClassId.build(0),
                passage_before=PassageState.UNVISITED,
                passage_after=PassageState.FIRST_PASS_COMPLETE,
                user_cap=EngagementCap.build(math.pi / 2.0),
                effective_cap=EngagementCap.build(math.pi / 3.0),
            ),
            traversal_decision=link.traversal_decision,
        ),
        LinkSegmentOperation.build(
            motion=link.motion,
            cut_z=link.cut_z,
            neck_scope=link.neck_scope,
            effective_cap_decision=link.effective_cap_decision,
            traversal_decision=HoldTraversalDecision.build(
                component_id=link.traversal_decision.component_id,
                edge_id=link.traversal_decision.edge_id,
                branch_id=link.traversal_decision.branch_id,
                cursor=CursorIdentity(b"cursor-other"),
            ),
        ),
    )
    circle_motion = ExactCircleMotion.build(
        Point2[WorldXY].build(1.0, 0.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )
    circle = CutFullCircleOperation.build(
        motion=circle_motion,
        cut_z=CutZ.build(-2.0),
        material_side=MaterialSide.OUTSIDE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=_full_cap(),
        traversal_decision=_advance(),
    )
    changed_cap_value = EngagementCap.build(math.pi / 3.0)
    circle_variants = (
        circle,
        CutFullCircleOperation.build(
            motion=ExactCircleMotion.build(circle.motion.center, circle.motion.phase_vector, True),
            cut_z=circle.cut_z,
            material_side=circle.material_side,
            neck_scope=circle.neck_scope,
            effective_cap_decision=circle.effective_cap_decision,
            traversal_decision=circle.traversal_decision,
        ),
        CutFullCircleOperation.build(
            motion=circle.motion,
            cut_z=CutZ.build(-3.0),
            material_side=circle.material_side,
            neck_scope=circle.neck_scope,
            effective_cap_decision=circle.effective_cap_decision,
            traversal_decision=circle.traversal_decision,
        ),
        CutFullCircleOperation.build(
            motion=circle.motion,
            cut_z=circle.cut_z,
            material_side=circle.material_side,
            neck_scope=OrientedNeckScope.build(
                neck_owner_id=NeckOwnerId(b"neck-a"),
                orientation=NeckTraversalOrientation.FORWARD,
            ),
            effective_cap_decision=_neck_cap(),
            traversal_decision=circle.traversal_decision,
        ),
        CutFullCircleOperation.build(
            motion=circle.motion,
            cut_z=circle.cut_z,
            material_side=circle.material_side,
            neck_scope=circle.neck_scope,
            effective_cap_decision=FullCapDecision.build(
                user_cap=changed_cap_value,
                effective_cap=changed_cap_value,
            ),
            traversal_decision=circle.traversal_decision,
        ),
        CutFullCircleOperation.build(
            motion=circle.motion,
            cut_z=circle.cut_z,
            material_side=circle.material_side,
            neck_scope=circle.neck_scope,
            effective_cap_decision=circle.effective_cap_decision,
            traversal_decision=_advance(terminal=True),
        ),
    )

    for variants in (approach_variants, plunge_variants, link_variants, circle_variants):
        assert len({operation.canonical_bytes for operation in variants}) == len(variants)


def test_canonical_operations_do_not_embed_or_trust_witnesses() -> None:
    operation_types = (
        ApproachOperation,
        PlungeOperation,
        LinkSegmentOperation,
        CutFullCircleOperation,
        AdvanceSegmentOperation,
    )

    for operation_type in operation_types:
        field_names = {field.name for field in fields(operation_type)}
        assert not any("witness" in field_name for field_name in field_names)
        assert "digest" not in field_names
