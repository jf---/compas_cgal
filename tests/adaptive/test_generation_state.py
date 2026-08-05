"""Closed operation-chronology contracts for authoritative generation state."""

from collections.abc import Sequence

import pytest

from compas_cgal.adaptive.errors import InvalidGenerationStateError
from compas_cgal.adaptive.generation_state import GenerationState
from compas_cgal.adaptive.generation_state import TraversalCursorState
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CanonicalOperation
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import CutFullCircleOperation
from compas_cgal.adaptive.operation import HoldTraversalDecision
from compas_cgal.adaptive.operation import LinkSegmentOperation
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.test_transaction import _StateFixture
from tests.adaptive.test_transaction import _state_fixture


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build one real exact entry-circle state shared by chronology cases."""
    return Task13FFixture.build()


@pytest.fixture(scope="module")
def oriented_state() -> _StateFixture:
    """Build one real state whose first neck passage is already accepted."""
    return _state_fixture(oriented=True)


def _point(x: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, 1.0625)


def _advance_decision(
    fixture: Task13FFixture,
    *,
    cursor_before: CursorIdentity,
    cursor_after: CursorIdentity,
    terminal: bool = False,
) -> AdvanceTraversalDecision:
    traversal = fixture.physical.traversal
    return AdvanceTraversalDecision.build(
        component_id=traversal.component_id,
        edge_id=traversal.edge_id,
        branch_id=traversal.branch_id,
        cursor_before=cursor_before,
        cursor_after=cursor_after,
        makes_cursor_terminal=terminal,
    )


def _hold_decision(
    fixture: Task13FFixture,
    cursor: CursorIdentity,
) -> HoldTraversalDecision:
    traversal = fixture.physical.traversal
    return HoldTraversalDecision.build(
        component_id=traversal.component_id,
        edge_id=traversal.edge_id,
        branch_id=traversal.branch_id,
        cursor=cursor,
    )


def _advance_segment(
    fixture: Task13FFixture,
    *,
    start: Point2[WorldXY],
    end: Point2[WorldXY],
    cursor_before: CursorIdentity,
    cursor_after: CursorIdentity,
) -> AdvanceSegmentOperation:
    final = fixture.physical.operations[-1]
    assert type(final) is CutFullCircleOperation
    return AdvanceSegmentOperation.build(
        motion=ExactSegmentMotion.build(start, end),
        cut_z=fixture.identity.cut_plane.cut_z,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=final.effective_cap_decision,
        traversal_decision=_advance_decision(
            fixture,
            cursor_before=cursor_before,
            cursor_after=cursor_after,
        ),
    )


def _hold_link(
    fixture: Task13FFixture,
    *,
    start: Point2[WorldXY],
    end: Point2[WorldXY],
    cursor: CursorIdentity,
) -> LinkSegmentOperation:
    final = fixture.physical.operations[-1]
    assert type(final) is CutFullCircleOperation
    return LinkSegmentOperation.build(
        motion=ExactSegmentMotion.build(start, end),
        cut_z=fixture.identity.cut_plane.cut_z,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=final.effective_cap_decision,
        traversal_decision=_hold_decision(fixture, cursor),
    )


def _circle(
    fixture: Task13FFixture,
    *,
    phase: Point2[WorldXY],
    cursor_before: CursorIdentity,
    cursor_after: CursorIdentity,
) -> CutFullCircleOperation:
    final = fixture.physical.operations[-1]
    assert type(final) is CutFullCircleOperation
    phase_vector = Vector2[WorldXY].build(0.0625, 0.0)
    return CutFullCircleOperation.build(
        motion=ExactCircleMotion.build(
            Point2[WorldXY].build(
                phase.x - phase_vector.x,
                phase.y,
            ),
            phase_vector,
            False,
        ),
        cut_z=fixture.identity.cut_plane.cut_z,
        material_side=MaterialSide.OUTSIDE,
        neck_scope=NoNeckScope.build(),
        effective_cap_decision=final.effective_cap_decision,
        traversal_decision=_advance_decision(
            fixture,
            cursor_before=cursor_before,
            cursor_after=cursor_after,
        ),
    )


def _build_state(
    fixture: Task13FFixture,
    *,
    lateral: Sequence[CanonicalOperation],
    phase_point: Point2[WorldXY],
    traversal: TraversalCursorState,
) -> GenerationState:
    stock = fixture.physical.fork_stock()
    coverage = fixture.physical.clone_coverage()
    for operation in lateral:
        assert type(operation) in (
            LinkSegmentOperation,
            CutFullCircleOperation,
            AdvanceSegmentOperation,
        )
        stock.deplete(
            operation.motion,
            fixture.identity.tool_radius,
            fixture.identity.depletion_policy,
        )
        coverage.add_sweep(
            operation.motion,
            fixture.identity.tool_radius,
        )
    return GenerationState.build(
        stock=stock,
        coverage=coverage,
        tool_radius=fixture.identity.tool_radius,
        phase_point=phase_point,
        traversal=traversal,
        passages=fixture.physical.passages,
        operations=fixture.physical.operations + tuple(lateral),
    )


def test_generation_state_accepts_entry_circle_then_advancing_segment(
    task13f: Task13FFixture,
) -> None:
    """Allow a zero-guide motion to own progress without a successor circle."""
    cursor_before = task13f.physical.traversal.cursor
    decision_after = CursorIdentity(b"zero-guide-cursor-1")
    segment = _advance_segment(
        task13f,
        start=task13f.physical.phase_point,
        end=_point(1.625),
        cursor_before=cursor_before,
        cursor_after=decision_after,
    )
    assert type(segment) is AdvanceSegmentOperation

    state = _build_state(
        task13f,
        lateral=(segment,),
        phase_point=segment.motion.end,
        traversal=task13f.physical.traversal.advance(
            segment.traversal_decision,
        ),
    )

    assert state.phase_point == segment.motion.end
    assert state.traversal.cursor == decision_after
    assert state.operations[-1] == segment


def test_generation_state_advances_neck_passage_on_advancing_segment(
    oriented_state: _StateFixture,
) -> None:
    """Replay neck-cap causality for either operation that owns an advance.

    Ignoring an advancing segment here would leave the immutable passage at
    its old state even though the operation canonically owns the next cap
    decision. That would split traversal and neck authority during replay.
    """
    scope = oriented_state.second.neck_scope
    decision = oriented_state.second.effective_cap_decision
    assert type(scope) is OrientedNeckScope
    assert type(decision) is NeckCapDecision
    passage_before = oriented_state.state.passage(scope)
    passage_after = passage_before.advance(decision)
    segment = AdvanceSegmentOperation.build(
        motion=ExactSegmentMotion.build(
            oriented_state.state.phase_point,
            oriented_state.second.middle_point,
        ),
        cut_z=oriented_state.identity.cut_plane.cut_z,
        neck_scope=scope,
        effective_cap_decision=decision,
        traversal_decision=oriented_state.second.traversal_decision,
    )
    stock = oriented_state.state.fork_stock()
    stock.deplete(
        segment.motion,
        oriented_state.identity.tool_radius,
        oriented_state.identity.depletion_policy,
    )
    coverage = oriented_state.state.clone_coverage()
    coverage.add_sweep(
        segment.motion,
        oriented_state.identity.tool_radius,
    )

    state = GenerationState.build(
        stock=stock,
        coverage=coverage,
        tool_radius=oriented_state.identity.tool_radius,
        phase_point=segment.motion.end,
        traversal=oriented_state.state.traversal.advance(
            segment.traversal_decision,
        ),
        passages=tuple(passage_after if passage.scope == scope else passage for passage in oriented_state.state.passages),
        operations=oriented_state.state.operations + (segment,),
    )

    assert state.passage(scope) == passage_after


def test_generation_state_accepts_circle_pair_then_advancing_segment(
    task13f: Task13FFixture,
) -> None:
    """Compose the old link-circle pair with the new one-motion chronology."""
    cursor_1 = task13f.physical.traversal.cursor
    cursor_2 = CursorIdentity(b"circle-cursor-2")
    cursor_3 = CursorIdentity(b"zero-guide-cursor-3")
    phase_1 = _point(1.625)
    phase_2 = _point(1.75)
    link = _hold_link(
        task13f,
        start=task13f.physical.phase_point,
        end=phase_1,
        cursor=cursor_1,
    )
    circle = _circle(
        task13f,
        phase=phase_1,
        cursor_before=cursor_1,
        cursor_after=cursor_2,
    )
    segment = _advance_segment(
        task13f,
        start=phase_1,
        end=phase_2,
        cursor_before=cursor_2,
        cursor_after=cursor_3,
    )
    traversal = task13f.physical.traversal.advance(
        circle.traversal_decision,
    ).advance(segment.traversal_decision)

    state = _build_state(
        task13f,
        lateral=(link, circle, segment),
        phase_point=phase_2,
        traversal=traversal,
    )

    assert state.operations[-3:] == (link, circle, segment)
    assert state.phase_point == phase_2
    assert state.traversal == traversal


def test_generation_state_accepts_advancing_segment_then_circle_pair(
    task13f: Task13FFixture,
) -> None:
    """Resume ordinary trochoids after one exact zero-guide MAT advance."""
    cursor_1 = task13f.physical.traversal.cursor
    cursor_2 = CursorIdentity(b"zero-guide-cursor-2")
    cursor_3 = CursorIdentity(b"circle-cursor-3")
    phase_1 = _point(1.625)
    phase_2 = _point(1.75)
    segment = _advance_segment(
        task13f,
        start=task13f.physical.phase_point,
        end=phase_1,
        cursor_before=cursor_1,
        cursor_after=cursor_2,
    )
    link = _hold_link(
        task13f,
        start=phase_1,
        end=phase_2,
        cursor=cursor_2,
    )
    circle = _circle(
        task13f,
        phase=phase_2,
        cursor_before=cursor_2,
        cursor_after=cursor_3,
    )
    traversal = task13f.physical.traversal.advance(
        segment.traversal_decision,
    ).advance(circle.traversal_decision)

    state = _build_state(
        task13f,
        lateral=(segment, link, circle),
        phase_point=phase_2,
        traversal=traversal,
    )

    assert state.operations[-3:] == (segment, link, circle)
    assert state.phase_point == phase_2
    assert state.traversal == traversal


def test_generation_state_rejects_dangling_hold_link(
    task13f: Task13FFixture,
) -> None:
    """Prevent a hold operation from ending state without its owned circle."""
    link = _hold_link(
        task13f,
        start=task13f.physical.phase_point,
        end=_point(1.625),
        cursor=task13f.physical.traversal.cursor,
    )

    with pytest.raises(InvalidGenerationStateError, match="hold link"):
        _build_state(
            task13f,
            lateral=(link,),
            phase_point=link.motion.end,
            traversal=task13f.physical.traversal,
        )


def test_generation_state_rejects_hold_link_followed_by_advancing_segment(
    task13f: Task13FFixture,
) -> None:
    """Forbid one hold from being consumed by the wrong advancing variant."""
    cursor_1 = task13f.physical.traversal.cursor
    cursor_2 = CursorIdentity(b"zero-guide-cursor-2")
    phase_1 = _point(1.625)
    link = _hold_link(
        task13f,
        start=task13f.physical.phase_point,
        end=phase_1,
        cursor=cursor_1,
    )
    segment = _advance_segment(
        task13f,
        start=phase_1,
        end=_point(1.75),
        cursor_before=cursor_1,
        cursor_after=cursor_2,
    )

    with pytest.raises(InvalidGenerationStateError, match="hold link"):
        _build_state(
            task13f,
            lateral=(link, segment),
            phase_point=segment.motion.end,
            traversal=task13f.physical.traversal.advance(
                segment.traversal_decision,
            ),
        )


def test_generation_state_rejects_circle_directly_after_advancing_segment(
    task13f: Task13FFixture,
) -> None:
    """Stop a circle from claiming an advance already owned by a segment."""
    cursor_1 = task13f.physical.traversal.cursor
    cursor_2 = CursorIdentity(b"shared-cursor-2")
    phase_1 = _point(1.625)
    segment = _advance_segment(
        task13f,
        start=task13f.physical.phase_point,
        end=phase_1,
        cursor_before=cursor_1,
        cursor_after=cursor_2,
    )
    circle = _circle(
        task13f,
        phase=phase_1,
        cursor_before=cursor_1,
        cursor_after=cursor_2,
    )

    with pytest.raises(InvalidGenerationStateError, match="circle.*link"):
        _build_state(
            task13f,
            lateral=(segment, circle),
            phase_point=phase_1,
            traversal=task13f.physical.traversal.advance(
                segment.traversal_decision,
            ),
        )


@pytest.mark.parametrize("mismatch", ("phase", "traversal"))
def test_generation_state_rejects_final_advancing_segment_state_mismatch(
    task13f: Task13FFixture,
    mismatch: str,
) -> None:
    """Bind both physical phase and graph cursor to the final segment."""
    cursor_before = task13f.physical.traversal.cursor
    segment = _advance_segment(
        task13f,
        start=task13f.physical.phase_point,
        end=_point(1.625),
        cursor_before=cursor_before,
        cursor_after=CursorIdentity(b"zero-guide-cursor-2"),
    )
    expected_traversal = task13f.physical.traversal.advance(
        segment.traversal_decision,
    )

    with pytest.raises(
        InvalidGenerationStateError,
        match=mismatch,
    ):
        _build_state(
            task13f,
            lateral=(segment,),
            phase_point=(task13f.physical.phase_point if mismatch == "phase" else segment.motion.end),
            traversal=(task13f.physical.traversal if mismatch == "traversal" else expected_traversal),
        )
