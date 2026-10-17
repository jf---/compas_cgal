"""The helical entry's structured record on the result (Task 1.2, brief N2; read by the entry figure).

`enter_first_cut` builds a `RampOutRecord` from what it actually emitted: the
helix's operation index, the pencil members' operation range (possibly
empty) and the root loop's operation index. `RampOutRecord.build` checks the
record against the stream and raises by name on any contradiction. The FIRST_CUT
ledger entry is recorded as a helix (`EntryKind.HELIX`); re-entries stay plunges.
"""

from __future__ import annotations

import math

import pytest
from compas.geometry import Circle
from compas.geometry import Polygon

import compas_cgal.pipeline.generate as generate
from compas_cgal.helix_units import DEFAULT_HELIX_PITCH_PER_RADIUS
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import EntryRiseChordRatio
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import RampOutRecordMismatchError
from compas_cgal.pipeline.errors import RootLoopNotEmittedError
from compas_cgal.pipeline.generate import HelicalEntry
from compas_cgal.pipeline.generate import PathState
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.pencil import entry_containment
from compas_cgal.pipeline.result import EntryKind
from compas_cgal.pipeline.result import ForcedEntryReason
from compas_cgal.pipeline.result import RampOutRecord
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tests.pipeline.test_gate_helix_entry import CASES
from tests.pipeline.test_gate_helix_entry import _first_cut_root
from tests.pipeline.test_gate_helix_entry import _generated
from tests.pipeline.test_gate_helix_entry import _is_loop

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)
# A root within one tool radius of the footprint (rho_root - rho_h <= r): the first scan's furthest rung is the root.
ONE_SCAN_ROOT_RADIUS = 1.5


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_the_record_matches_the_emitted_operations(cap: float, pocket: str) -> None:
    result, _survey = _generated(cap, pocket)
    record = result.ramp_out
    operations = result.operations
    assert record.helix_operation == 0 and operations[0].operation is OperationType.HELIX
    assert record.members.start == record.helix_operation + 1 and record.members.stop == record.root_operation
    assert all(operations[index].operation is OperationType.CUT and isinstance(operations[index].geometry, Circle) for index in record.members)
    assert _is_loop(operations[record.root_operation], _first_cut_root(cap, pocket))
    assert record == RampOutRecord.build(operations, record.path_index, record.helix_operation, record.members, record.root_operation)


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_the_first_cut_is_ledgered_as_a_helix_and_re_entries_as_plunges(cap: float, pocket: str) -> None:
    result, _survey = _generated(cap, pocket)
    kinds = [(entry.reason is ForcedEntryReason.FIRST_CUT, entry.entry_kind) for entry in result.ledger.entries]
    assert all(kind is (EntryKind.HELIX if first else EntryKind.PLUNGE) for first, kind in kinds)


def _one_scan_path() -> tuple[OrientedChain, PathState]:
    stations = tuple(
        GuideStation(
            centre=Point2[WorldXY].build(-1.0 + 1.2 * index, 0.0),
            radius=GuideRadius.build(ONE_SCAN_ROOT_RADIUS),
            tangent=Direction2[WorldXY].build(1.0, 0.0),
            clockwise=True,
        )
        for index in range(2)
    )
    chain = OrientedChain.build(GuideChain(stations=stations, origin=0), TOOL)
    assert chain is not None
    # Half-turn cap and rise limit: every rung passes, so the furthest one (the root) is admitted at once.
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(math.pi)), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=EntryRiseChordRatio(cap_ratio(Radian(math.pi))))
    helical_entry = HelicalEntry(pitch_per_radius=DEFAULT_HELIX_PITCH_PER_RADIUS, pocket=entry_containment(RECT_12X8, [], TOOL))
    return chain, PathState(stock=Stock(RECT_12X8), admissibility=admissibility, cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0), helical_entry=helical_entry)


def test_a_root_within_one_tool_radius_of_the_footprint_has_no_members() -> None:
    chain, path = _one_scan_path()
    record = generate.enter_first_cut(path, chain, chain.first)
    assert record.members == range(1, 1) and record.root_operation == 1
    assert [op.operation for op in path.operations] == [OperationType.HELIX, OperationType.CUT]


def test_a_root_loop_skipped_as_removing_nothing_raises_by_name(monkeypatch: pytest.MonkeyPatch) -> None:
    """Unreachable by geometry (the root's swept disk minus the last member's is a lune of positive area in the pocket); fail loud if it ever fires."""
    chain, path = _one_scan_path()
    monkeypatch.setattr(generate, "skip_loop_removing_nothing", lambda *_args: True)
    with pytest.raises(RootLoopNotEmittedError):
        generate.enter_first_cut(path, chain, chain.first)


def test_a_record_that_contradicts_the_stream_raises_by_name() -> None:
    chain, path = _one_scan_path()
    record = generate.enter_first_cut(path, chain, chain.first)
    with pytest.raises(RampOutRecordMismatchError):
        RampOutRecord.build(path.operations, record.path_index, record.helix_operation, range(1, 2), 1)
    with pytest.raises(RampOutRecordMismatchError):
        RampOutRecord.build(path.operations, record.path_index, 1, record.members, record.root_operation)
