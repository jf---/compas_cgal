"""The narrow pipeline on the gate pockets, judged by the gate's own survey and assessment (Task 1.1 acceptance).

Every case is generated through `benchmarks.gate.pipeline`, the gate adapter,
so a number asserted here is the number the gate vectors pin. One plunge into
material is asserted through the gate's ``material_plunges`` cell, never by
counting PLUNGE operations (Ruling 1.1-G13).
"""

from __future__ import annotations

import functools
import math
from collections import Counter
from collections.abc import Callable
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import pytest
from compas.geometry import Circle
from compas.geometry import Polygon

from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_POCKET_NAMES
from benchmarks.gate import gate_pocket
from benchmarks.gate import pipeline
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.quality_observations import QualityEvidence
from benchmarks.quality_observations import reduce_quality_evidence
from benchmarks.survey import MotionKind
from benchmarks.survey import PathSurvey
from benchmarks.survey import survey_path
from compas_cgal import _stock_2
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import loop_parameters
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.errors import CoincidentStationsError
from compas_cgal.pipeline.errors import LedgerInvariantError
from compas_cgal.pipeline.generate import LoopConditions
from compas_cgal.pipeline.generate import PathState
from compas_cgal.pipeline.generate import enter_first_cut
from compas_cgal.pipeline.generate import walk_chain
from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from compas_cgal.pipeline.generate import DEFAULT_MAX_PASSES
from compas_cgal.pipeline.guide import Guide
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import pipeline_spacing
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.units import Clearance
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.narrow import Advance
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.narrow import advance_ladder
from compas_cgal.pipeline.result import DeclineCause
from compas_cgal.pipeline.result import DeclinedRun
from compas_cgal.pipeline.result import ForcedDecisionLedger
from compas_cgal.pipeline.result import ForcedEntry
from compas_cgal.pipeline.result import ForcedEntryReason
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import load_gate_vectors

# One xdist group per case, so a case's tests share one generation on one worker.
CASES = [
    pytest.param(cap, pocket, id=f"{cap_id}-{pocket}", marks=pytest.mark.xdist_group(f"pipeline-{cap_id}-{pocket}"))
    for cap_id, cap in GATE_CAP_CASES
    for pocket in GATE_POCKET_NAMES
]
CASES_AT_120 = [case for case in CASES if case.values[0] == GATE_CAP_CASES[0][1]]

# rho_regulated's measured 120-degree tangent_breaks cells (P1 at a88cdab0): the
# acceptance ceiling is rho's own column, never a tuned number.
EXPECTED_TANGENT_CEILING = {"rect_12x8": 0, "rect_20x12": 0, "L_shape": 2}


@functools.cache
def _generated(cap: float, pocket: str) -> tuple[PipelineResult, PathSurvey, QualityEvidence]:
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
    result = pipeline(spec)
    assert isinstance(result, PipelineResult)
    survey = survey_path(spec, result)
    return result, survey, reduce_quality_evidence(spec, snapshot_toolpath(result), survey)


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_one_plunge_into_material_and_one_first_cut(cap: float, pocket: str) -> None:
    result, survey, evidence = _generated(cap, pocket)
    assert int(evidence.assessment.material_plunges.measured) == 1
    entries = result.ledger.entries
    assert [entry.reason for entry in entries].count(ForcedEntryReason.FIRST_CUT) == 1
    plunge_feet = {(op.geometry.end[0], op.geometry.end[1]) for op in (result.operations[int(i)] for i in survey.material_plunge_indices)}
    assert {(entry.entry.x, entry.entry.y) for entry in entries if entry.into_material} == plunge_feet


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_stage_0_rows_hold_on_the_pipeline(cap: float, pocket: str) -> None:
    """The stage-0 all-case rows bind the pipeline at both caps, with no exemption."""
    _result, _survey, evidence = _generated(cap, pocket)
    assessment = evidence.assessment
    assert (int(assessment.misplaced_plunges.measured), int(assessment.unsafe_rapids.measured), int(assessment.zero_length_motions.measured)) == (0, 0, 0)


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_every_link_at_depth_removes_nothing(cap: float, pocket: str) -> None:
    """C11: the gate does not flag a LINK that removes, so the pipeline certifies every cut-plane LINK is NOTHING."""
    _result, survey, _evidence = _generated(cap, pocket)
    links = [motion for motion in survey.motions if motion.operation is OperationType.LINK]
    assert links, "the pipeline links at depth on every gate pocket"
    assert all(motion.removal is SweepRemoval.NOTHING for motion in links)


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_no_within_band_motion_is_labelled_link(cap: float, pocket: str) -> None:
    """C17: a motion whose bracket is WITHIN_BAND (or REMOVES) is a CUT, never a LINK."""
    _result, survey, _evidence = _generated(cap, pocket)
    assert all(motion.operation is OperationType.CUT for motion in survey.motions if motion.removal is not SweepRemoval.NOTHING)


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_every_emitted_loop_is_authored_exactly_incident(cap: float, pocket: str) -> None:
    result, _survey, _evidence = _generated(cap, pocket)
    loops = [op for op in result.operations if isinstance(op.geometry, Circle)]
    assert loops and all(op.entry_xy is not None for op in loops)
    for op in loops:
        assert op.entry_xy is not None
        centre = op.geometry.frame.point
        dx, dy = Fraction(op.entry_xy[0]) - Fraction(float(centre[0])), Fraction(op.entry_xy[1]) - Fraction(float(centre[1]))
        assert dx * dx + dy * dy == Fraction(float(op.geometry.radius)) ** 2


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_no_loop_that_removes_nothing_is_emitted(cap: float, pocket: str) -> None:
    """The NOTHING-loop skip: every emitted loop removes material; the skipped ones are reported, not hidden."""
    result, survey, _evidence = _generated(cap, pocket)
    assert all(motion.removes_material for motion in survey.motions if motion.kind is MotionKind.LOOP)
    assert result.omitted and all(loop.verdict is SweepRemoval.NOTHING for loop in result.omitted)


@pytest.mark.parametrize(("cap", "pocket"), CASES_AT_120)
def test_the_narrow_acceptance_at_120(cap: float, pocket: str) -> None:
    """Task 1.1 at 120 deg: redundant 0, degenerate 0, tangent breaks within rho's own cells."""
    _result, _survey, evidence = _generated(cap, pocket)
    assert int(evidence.assessment.redundant_operations.measured) == 0
    assert evidence.path_quality.elementary.degenerate_loops == 0
    assert int(evidence.assessment.tangent_breaks.measured) <= EXPECTED_TANGENT_CEILING[pocket]


def test_ledger_has_exactly_one_first_cut() -> None:
    entry = Point2[WorldXY].build(0.0, 0.0)
    first = ForcedEntry(path_index=0, entry=entry, reason=ForcedEntryReason.FIRST_CUT, into_material=True)
    retract = ForcedEntry(path_index=1, entry=entry, reason=ForcedEntryReason.NO_ADMISSIBLE_LINK, into_material=False)
    assert ForcedDecisionLedger.build([first, retract]).entries == (first, retract)
    with pytest.raises(LedgerInvariantError):
        ForcedDecisionLedger.build([retract])
    with pytest.raises(LedgerInvariantError):
        ForcedDecisionLedger.build([first, first])


RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)


def test_a_chain_split_at_an_unbridgeable_step_is_re_entered_as_a_forced_entry() -> None:
    """B-unbridgeable: the walk records the gap and re-enters past it; it never raises, never bridges the bound."""
    radii = [2.0, 2.0, 0.5, 0.5, 0.5, 0.5, 2.0, 2.0]
    stations = tuple(
        GuideStation(centre=Point2[WorldXY].build(-4.2 + 1.2 * index, 0.0), radius=GuideRadius.build(value), tangent=Direction2[WorldXY].build(1.0, 0.0), clockwise=True)
        for index, value in enumerate(radii)
    )
    chain = OrientedChain.build(GuideChain(stations=stations, origin=0), TOOL)
    assert chain is not None
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(3.0)), centre_domain(RECT_12X8, [], TOOL))
    path = PathState(stock=Stock(RECT_12X8), admissibility=admissibility, cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0))
    enter_first_cut(path, chain, chain.first)
    end = walk_chain(path, chain, chain.first, 3, LoopConditions(admissibility), emit_start_loop=False)
    assert DeclinedRun(chain=0, first_station=1, last_station=6, cause=DeclineCause.RADIUS_STEP) in path.records
    split = [record for record in path.records if isinstance(record, ForcedEntry) and record.reason is ForcedEntryReason.RADIUS_STEP_SPLIT]
    resumed = chain.loops[6]
    assert isinstance(resumed, IncidentLoop)
    assert [(record.entry, record.into_material) for record in split] == [(resumed.entry, True)]
    assert end.loop == chain.loops[chain.last]


def _synthetic_chain(radii: list[float]) -> OrientedChain:
    stations = tuple(
        GuideStation(centre=Point2[WorldXY].build(-4.2 + 1.2 * index, 0.0), radius=GuideRadius.build(value), tangent=Direction2[WorldXY].build(1.0, 0.0), clockwise=True)
        for index, value in enumerate(radii)
    )
    chain = OrientedChain.build(GuideChain(stations=stations, origin=0), TOOL)
    assert chain is not None
    return chain


def _synthetic_path() -> tuple[PathState, Admissibility]:
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(3.0)), centre_domain(RECT_12X8, [], TOOL))
    return PathState(stock=Stock(RECT_12X8), admissibility=admissibility, cut=CutZ.build(0.0), clearance=ClearanceZ.build(4.0)), admissibility


def test_the_entry_stage_is_the_plunge_and_the_root_loop() -> None:
    """The entry stage is one stage: `enter` plunges and cuts the root loop; the walk after it does not cut it again."""
    chain = _synthetic_chain([2.0, 2.0, 2.0])
    path, admissibility = _synthetic_path()
    enter_first_cut(path, chain, chain.first)
    root = chain.loops[chain.first]
    assert root is not None
    assert [op.operation for op in path.operations] == [OperationType.PLUNGE, OperationType.CUT]
    assert path.operations[-1].entry_xy == (root.entry.x, root.entry.y)
    walk_chain(path, chain, chain.first, 3, LoopConditions(admissibility), emit_start_loop=False)
    assert [op.entry_xy for op in path.operations].count((root.entry.x, root.entry.y)) == 1


def test_coincident_consecutive_stations_raise_by_name() -> None:
    """R1-loud: no advance direction exists between coincident centres; refuse loudly, never divide by zero."""
    chain = _synthetic_chain([2.0, 2.0])
    _path, admissibility = _synthetic_path()
    candidate = chain.candidate(0)
    with pytest.raises(CoincidentStationsError):
        LoopConditions(admissibility).loop_admissible_at_rest(Stock(RECT_12X8), candidate, candidate)


def test_the_at_rest_ring_short_circuits_before_the_in_process_bridge(monkeypatch: pytest.MonkeyPatch) -> None:
    """Evaluation order: a candidate the cheap at-rest ring refuses never pays for the 46-station in-process bridge.

    The ladder takes the conjunction, so the order decides nothing; it decides
    cost only (Task 1.1 perf, docs/pipeline.md).
    """
    chain = _synthetic_chain([2.0, 2.0, 2.0, 2.0, 2.0])
    path, admissibility = _synthetic_path()
    calls: Counter[str] = Counter()

    def refuse(_self: Admissibility, *_args: object) -> bool:
        calls["loop_within_cap_at_rest"] += 1
        return False

    def counted(name: str) -> Callable[..., bool]:
        original = getattr(Admissibility, name)

        def wrapper(self: Admissibility, *args: object) -> bool:
            calls[name] += 1
            return bool(original(self, *args))

        return wrapper

    monkeypatch.setattr(Admissibility, "loop_within_cap_at_rest", refuse)
    for name in ("straight_in_domain", "straight_within_cap", "arc_within_cap", "loop_within_cap"):
        monkeypatch.setattr(Admissibility, name, counted(name))
    advance = advance_ladder(path.stock, chain, chain.first, chain.last, LoopConditions(admissibility).ordered)
    assert isinstance(advance, Advance) and advance.forced is not None
    assert calls == Counter({"loop_within_cap_at_rest": chain.last - chain.first})


def test_the_at_rest_ring_refuses_a_candidate_the_in_process_checks_admit() -> None:
    """I3 (User ruling (2026-10-16), decision 1 = A): the at-rest ring decides on its own; dropping it changes the advance.

    After the root loop, the next station's bridge and full loop are under the
    cap in process, but rho's at-rest probe ring is not: the ladder must refuse
    the candidate and record a forced advance.
    """
    chain = _synthetic_chain([2.0, 2.0, 2.0])
    path, admissibility = _synthetic_path()
    enter_first_cut(path, chain, chain.first)
    conditions = LoopConditions(admissibility)
    origin, candidate = chain.candidate(chain.first), chain.candidate(chain.first + 1)
    assert conditions.bridge_admissible(path.stock, origin, candidate)
    assert conditions.loop_admissible(path.stock, origin, candidate)
    assert not conditions.loop_admissible_at_rest(path.stock, origin, candidate)
    advance = advance_ladder(path.stock, chain, chain.first, chain.first + 1, conditions.ordered)
    assert isinstance(advance, Advance) and advance.forced is not None


def test_the_ladders_loop_condition_leaves_station_zero_to_the_ring(monkeypatch: pytest.MonkeyPatch) -> None:
    """Perf round 2: the ring's entry probe IS the loop's station 0 (same point, same stock), so the conjunction asks it once."""
    chain = _synthetic_chain([2.0, 2.0, 2.0])
    path, admissibility = _synthetic_path()
    enter_first_cut(path, chain, chain.first)
    asked: list[list[float]] = []
    native = _stock_2.engagement_circle_exceeded

    def recorded(*args: object) -> bool:
        asked.append(list(args[7]))  # type: ignore[call-overload]
        return bool(native(*args))  # type: ignore[arg-type]

    monkeypatch.setattr(_stock_2, "engagement_circle_exceeded", recorded)
    # Let the ring admit, so the scan reaches the loop condition (the ring's verdict is not what is pinned here).
    monkeypatch.setattr(Admissibility, "loop_within_cap_at_rest", lambda *_args: True)
    conditions = LoopConditions(admissibility)
    assert conditions.loop_admissible_at_rest in conditions.ordered
    advance_ladder(path.stock, chain, chain.first, chain.last, conditions.ordered)
    assert asked and all(0.0 not in parameters for parameters in asked)
    assert all(parameters == loop_parameters()[1:] for parameters in asked)


def test_every_ladder_condition_reads_one_unchanged_stock(monkeypatch: pytest.MonkeyPatch) -> None:
    """The station-0 reuse never crosses a stock mutation: one conjunction, one stock object, never mutated during the scan."""
    chain = _synthetic_chain([2.0, 2.0, 2.0, 2.0])
    path, admissibility = _synthetic_path()
    enter_first_cut(path, chain, chain.first)
    before = path.stock.clone()
    seen: list[tuple[int, bool]] = []

    def watched(condition: Callable[..., bool]) -> Callable[..., bool]:
        def wrapper(stock: Stock, origin: object, candidate: object) -> bool:
            seen.append((id(stock), stock.exactly_equals(before)))
            return condition(stock, origin, candidate)

        return wrapper

    ordered = tuple(watched(condition) for condition in LoopConditions(admissibility).ordered)
    advance_ladder(path.stock, chain, chain.first, chain.last, ordered)
    assert seen and all(stock_id == id(path.stock) and unchanged for stock_id, unchanged in seen)
    assert path.stock.exactly_equals(before)


ROOT = Path(__file__).resolve().parents[2]
CAP_IDS = {cap: cap_id for cap_id, cap in GATE_CAP_CASES}
ENGAGEMENT_STEP = "max_engagement_step_deg"


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_the_engagement_step_is_within_rhos_committed_cell(cap: float, pocket: str) -> None:
    """I3: freshly generated, the pipeline's engagement step stays within rho_regulated's committed cell.

    This is the cell the at-rest ring closed (docs/pipeline.md); it is held on
    a fresh generation, not only on the committed vectors.
    """
    _result, _survey, evidence = _generated(cap, pocket)
    vectors = load_gate_vectors(ROOT / GATE_VECTORS_PATH)
    criterion = next(spec for spec in vectors.criteria if spec.name == ENGAGEMENT_STEP)
    ours = Decimal(criterion.format(float(evidence.assessment.max_engagement_step.measured)))
    assert ours <= vectors.value(f"{CAP_IDS[cap]}/rho_regulated/{pocket}", ENGAGEMENT_STEP)


def _forced_advance_centres(cap: float, pocket: str, result: PipelineResult) -> set[tuple[float, float]]:
    """Centres of the stations forced advances landed on (the ledger names them in guide order)."""
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
    guide = Guide.build(
        spec.polygon,
        list(spec.holes),
        ToolRadius.build(spec.tool_radius),
        pipeline_spacing(spec.tool_diameter, Radian(math.radians(cap))),
        Clearance.build(RADIAL_CLEARANCE_FRACTION * spec.tool_diameter),
        climb=True,
        max_passes=DEFAULT_MAX_PASSES,
        clearance_z=ClearanceZ.build(CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter),
    )
    centres = set()
    for advance in result.ledger.advances:
        centre = guide.chains[advance.chain].stations[advance.station].centre
        centres.add((centre.x, centre.y))
    return centres


@pytest.mark.parametrize(("cap", "pocket"), CASES)
def test_no_decided_motion_exceeds_the_cap_at_a_gate_sample(cap: float, pocket: str) -> None:
    """I4: only forced motions may exceed the cap at the gate's own samples.

    Forced are the loop entered from a plunge, and a forced advance's bridge
    and loop (the ledger's `ForcedAdvance`); every other motion was admitted
    by an in-process check at the gate's own sample stations, on the gate's
    stock model. A link's pieces and the root loop it arrives at are judged on
    the pre-link stock, which holds at least the material the gate replays
    them against, so that check refuses more, never admits more.
    """
    result, survey, _evidence = _generated(cap, pocket)
    operations = result.operations
    forced = {index + 1 for index, op in enumerate(operations) if op.operation is OperationType.PLUNGE}
    landing = _forced_advance_centres(cap, pocket, result)
    for index, op in enumerate(operations):
        if isinstance(op.geometry, Circle) and (float(op.geometry.frame.point[0]), float(op.geometry.frame.point[1])) in landing:
            forced |= {index, index - 1}
    exceeded = {int(motion.index) for motion in survey.motions if motion.cap_exceeded}
    assert exceeded <= forced
