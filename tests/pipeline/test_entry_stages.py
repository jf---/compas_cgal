"""The entry stages are split by reason (Task 1.2 amendment M2): the first cut is its own stage.

The FIRST_CUT stage is the helix, its pencil ramp-out and the root loop; the
two re-entry reasons (NO_ADMISSIBLE_LINK, RADIUS_STEP_SPLIT) keep the
retract, rapid and plunge.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pytest
from compas.geometry import Polygon

import compas_cgal.pipeline.generate as generate
from compas_cgal.helix_units import DEFAULT_HELIX_PITCH_PER_RADIUS
from compas_cgal.helix_units import EntryRampLength
from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import gate_pocket
from benchmarks.gate import pipeline
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import entry_rise_ratio
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import FirstCutAfterEmissionError
from compas_cgal.pipeline.errors import FirstCutReEntryError
from compas_cgal.pipeline.errors import MissingHelicalEntryError
from compas_cgal.pipeline.generate import HelicalEntry
from compas_cgal.pipeline.generate import PathState
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.pencil import entry_containment
from compas_cgal.pipeline.result import ForcedEntry
from compas_cgal.pipeline.result import ForcedEntryReason
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)
RE_ENTRY_REASONS = (ForcedEntryReason.NO_ADMISSIBLE_LINK, ForcedEntryReason.RADIUS_STEP_SPLIT)
# A permissive cap (radians, below pi) so only the stage split is under test.
PERMISSIVE_CAP = Radian(3.0)
# First station x, mm: the root's swept disk (guide radius 2 + tool radius 1) stays
# inside the 12 x 8 rectangle, so the pencil is certified in the pocket (P-D).
STATION_X0 = -2.4


def _chain() -> OrientedChain:
    stations = tuple(
        GuideStation(centre=Point2[WorldXY].build(STATION_X0 + 1.2 * index, 0.0), radius=GuideRadius.build(2.0), tangent=Direction2[WorldXY].build(1.0, 0.0), clockwise=True)
        for index in range(3)
    )
    chain = OrientedChain.build(GuideChain(stations=stations, origin=0), TOOL)
    assert chain is not None
    return chain


def _path() -> PathState:
    admissibility = Admissibility.build(
        TOOL, cap_ratio(PERMISSIVE_CAP), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=entry_rise_ratio(PERMISSIVE_CAP, EntryRampLength.default())
    )
    helical_entry = HelicalEntry(pitch_per_radius=DEFAULT_HELIX_PITCH_PER_RADIUS, pocket=entry_containment(RECT_12X8, [], TOOL))
    return PathState(stock=Stock(RECT_12X8), admissibility=admissibility, cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0), helical_entry=helical_entry)


def test_the_first_cut_without_a_helical_entry_declaration_raises_by_name() -> None:
    """No silent placeholder: a path with no declared helix pitch and pocket cannot enter, and nothing is emitted."""
    chain, path = _chain(), _path()
    bare = PathState(stock=path.stock, admissibility=path.admissibility, cut=path.cut, clearance=path.clearance)
    with pytest.raises(MissingHelicalEntryError):
        generate.enter_first_cut(bare, chain, chain.first)
    assert bare.operations == [] and bare.records == []


def test_the_first_cut_is_its_own_stage() -> None:
    """`enter_first_cut` emits the helix, then CUT loops entered at the root's entry, the root last; one FIRST_CUT record."""
    chain, path = _chain(), _path()
    generate.enter_first_cut(path, chain, chain.first)
    root = chain.loops[chain.first]
    assert root is not None
    entry = (root.entry.x, root.entry.y)
    operations = [op.operation for op in path.operations]
    assert operations[0] is OperationType.HELIX and set(operations[1:]) == {OperationType.CUT}
    assert all(op.entry_xy == entry for op in path.operations)
    assert path.operations[-1].geometry.radius == root.radius.value
    assert path.records == [ForcedEntry(path_index=0, entry=root.entry, reason=ForcedEntryReason.FIRST_CUT, into_material=True)]


def test_the_re_entry_stage_refuses_a_first_cut_by_name() -> None:
    """A re-entry is never the first cut: `re_enter` raises before touching the state."""
    chain, path = _chain(), _path()
    root = chain.loops[chain.first]
    assert root is not None
    with pytest.raises(FirstCutReEntryError):
        generate.re_enter(path, chain, chain.first, ForcedEntryReason.FIRST_CUT, root.entry)  # type: ignore[arg-type]
    assert path.operations == [] and path.records == []


def test_the_re_entry_stage_retracts_rapids_and_plunges() -> None:
    """The re-entry keeps today's sequence: retract, rapid at clearance, plunge, root loop."""
    chain, path = _chain(), _path()
    root = chain.loops[chain.first]
    assert root is not None
    origin = Point2[WorldXY].build(-5.0, 0.0)
    generate.re_enter(path, chain, chain.first, ForcedEntryReason.NO_ADMISSIBLE_LINK, origin)
    assert [op.operation for op in path.operations] == [OperationType.RETRACT, OperationType.LINK, OperationType.PLUNGE, OperationType.CUT]
    assert path.records == [ForcedEntry(path_index=0, entry=root.entry, reason=ForcedEntryReason.NO_ADMISSIBLE_LINK, into_material=True)]


def test_the_first_cut_refuses_a_path_already_emitted() -> None:
    """The first cut starts at the clearance plane with nothing emitted; a second first cut raises before touching the state."""
    chain, path = _chain(), _path()
    generate.enter_first_cut(path, chain, chain.first)
    operations, records = list(path.operations), list(path.records)
    with pytest.raises(FirstCutAfterEmissionError):
        generate.enter_first_cut(path, chain, chain.first + 1)
    assert path.operations == operations and path.records == records


def _spy(monkeypatch: pytest.MonkeyPatch, name: str, calls: list[tuple[str, ForcedEntryReason]]) -> None:
    original: Callable[..., None] = getattr(generate, name)

    def spy(path: PathState, *args: Any) -> None:
        before = len(path.records)
        original(path, *args)
        for record in path.records[before:]:
            if isinstance(record, ForcedEntry):
                calls.append((name, record.reason))

    monkeypatch.setattr(generate, name, spy)


@pytest.mark.parametrize("cap", [cap for _cap_id, cap in GATE_CAP_CASES])
def test_only_the_first_cut_stage_records_first_cut_and_re_entries_still_plunge(monkeypatch: pytest.MonkeyPatch, cap: float) -> None:
    """On the L at both caps: one FIRST_CUT, a helix from `enter_first_cut`; every other entry is a re-entry PLUNGE from `re_enter`."""
    calls: list[tuple[str, ForcedEntryReason]] = []
    _spy(monkeypatch, "enter_first_cut", calls)
    _spy(monkeypatch, "re_enter", calls)
    result = pipeline(gate_pocket("L_shape", tea_cap_deg=cap))  # type: ignore[arg-type]
    entries = [record for record in result.ledger.records if isinstance(record, ForcedEntry)]
    assert calls == [("enter_first_cut", ForcedEntryReason.FIRST_CUT)] + [("re_enter", entry.reason) for entry in entries[1:]]
    assert entries[0].reason is ForcedEntryReason.FIRST_CUT
    assert len(entries) > 1 and all(entry.reason in RE_ENTRY_REASONS for entry in entries[1:])
    plunges = [op for op in result.operations if op.operation is OperationType.PLUNGE]
    assert len(plunges) == len(entries) - 1
    assert [op.operation for op in result.operations].count(OperationType.HELIX) == 1
