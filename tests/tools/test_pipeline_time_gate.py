"""The Task 1.2 generation-time gate (D6): generated from three records, a verdict only from quiet-machine evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from tests.tools.test_run_provenance import provenance

from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import LoadAverage
from tools.held_generation_clock import MachineNotQuietError
from tools.pipeline_pencil_cost_results import ENTER_FIRST_CUT
from tools.pipeline_pencil_cost_results import RAMP_OUT_BUILD
from tools.pipeline_pencil_cost_results import RAMP_OUT_NEXT_MEMBER
from tools.pipeline_pencil_cost_results import TANGENT_PENCIL_BUILD
from tools.pipeline_pencil_cost_results import CasePencilCost
from tools.pipeline_pencil_cost_results import PencilRun
from tools.pipeline_pencil_cost_results import PIPELINE_PENCIL_COST_PATH
from tools.pipeline_pencil_cost_results import PipelinePencilCosts
from tools.pipeline_pencil_cost_results import load_pipeline_pencil_costs
from tools.pipeline_pencil_cost_results import pencil_case_keys
from tools.pipeline_time_gate import COMMITTED_1_1_REV
from tools.pipeline_time_gate import PIPELINE_TIME_GATE_TABLE_PATH
from tools.pipeline_time_gate import BACK_TO_BACK_WINDOW_S
from tools.pipeline_time_gate import BaseRevisionMismatchError
from tools.pipeline_time_gate import CommittedBaseError
from tools.pipeline_time_gate import BaseToolMismatchError
from tools.pipeline_time_gate import DirtyBaseError
from tools.pipeline_time_gate import EnvironmentMismatchError
from tools.pipeline_time_gate import ForeignBaseLibraryError
from tools.pipeline_time_gate import HostMismatchError
from tools.pipeline_time_gate import LockMismatchError
from tools.pipeline_time_gate import NativeLibraryMismatchError
from tools.pipeline_time_gate import NotBackToBackError
from tools.pipeline_time_gate import PencilTreeMismatchError
from tools.pipeline_time_gate import SessionBaseOutsideRepositoryError
from tools.pipeline_time_gate import IncompleteTimeGateInputError
from tools.pipeline_time_gate import PipelineTimeGate
from tools.pipeline_time_gate import StalePencilCostError
from tools.pipeline_time_gate import TimeGateBase
from tools.pipeline_time_gate import TimeGateMode
from tools.pipeline_time_gate import UnrecordedLoadError
from tools.pipeline_time_gate import committed_time_gate_base
from tools.pipeline_time_gate import main
from tools.pipeline_time_gate import render_time_gate_table
from tools.pipeline_time_gate import session_time_gate_base
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import CaseTiming
from tools.pipeline_timing_results import PipelineTimings
from tools.run_provenance import RunProvenance

ROOT = Path(__file__).resolve().parents[2]
KEYS = pencil_case_keys()
# Quiet at 1 minute, loaded over the 15 minutes before the end of the run.
LOADED_MID_RUN: LoadAverage = (1.0, 1.5, QUIET_LOAD_AVERAGE_MAX + 0.5)
SMOKE_KEY = "cap-120-default/pipeline/rect_12x8"
DIGEST = "0123456789abcdef"
OTHER_DIGEST = "fedcba9876543210"
TIMING_SHA256 = hashlib.sha256(b"the 1.2 timing document").hexdigest()
BASE_SHA256 = "1" * 64
OPERATIONS = 7
MEMBERS = 4
CALLS = {TANGENT_PENCIL_BUILD: 1, RAMP_OUT_BUILD: 1, RAMP_OUT_NEXT_MEMBER: MEMBERS + 1, ENTER_FIRST_CUT: 1}
# Every value below is exact in binary, so the equality edge is an exact equality.
UNSORTED_DECIDE = (0.375, 0.125, 0.25)
C_PENCIL = 0.25
BASE_RUNS = (1.0, 1.5, 1.25)
BASE_SLOWEST = 1.5
BOUND = BASE_SLOWEST + C_PENCIL
EDGE_MEDIAN = 1.75
HOLDING_MEDIAN = 1.625
FAILING_MEDIAN = 1.875
FIRST_CUT_S = 0.5
GENERATION_S = 1.0
QUIET: LoadAverage = (0.5, 0.75, 1.25)
LOADED: LoadAverage = (QUIET_LOAD_AVERAGE_MAX + 1.0, 3.5, 4.5)
COMMITTED_1_1_SLOWEST_120_RECT_12X8 = 0.6331
BASE_REV = "7b001d473d6e94b79fc8ff3cb852ecde1ed67733"
HEAD_REV = "b" * 40
BASE_TOOL_SHA = "9" * 64
# The 1.1 build's extension digest; the 1.2 fixtures carry provenance()'s default "a" * 64.
BASE_STOCK_2 = "f" * 64
PATCH = "c" * 64
# The 1.1 tree's native sources digest; the 1.2 fixtures carry provenance()'s default "a" * 64.
BASE_NATIVE_SOURCES = "d" * 64
# Back to back: base 18:00-18:10, 1.2 timing 18:12-18:30, pencil 18:31-18:40 (UTC).
BASE_INTERVAL = ("2026-10-22T18:00:00+00:00", "2026-10-22T18:10:00+00:00")
TIMING_INTERVAL = ("2026-10-22T18:12:00+00:00", "2026-10-22T18:30:00+00:00")
PENCIL_INTERVAL = ("2026-10-22T18:31:00+00:00", "2026-10-22T18:40:00+00:00")


def _base_provenance(**overrides: object) -> RunProvenance:
    fields: dict[str, object] = {
        "rev": BASE_REV,
        "started_utc": BASE_INTERVAL[0],
        "ended_utc": BASE_INTERVAL[1],
        "timing_tool_sha256": BASE_TOOL_SHA,
        "stock_2_sha256": BASE_STOCK_2,
        "native_sha256": BASE_STOCK_2,
        "native_sources_sha256": BASE_NATIVE_SOURCES,
    }
    fields.update(overrides)
    return provenance(**fields)


def _head_provenance(interval: tuple[str, str], **overrides: object) -> RunProvenance:
    fields: dict[str, object] = {"rev": HEAD_REV, "dirty": True, "patch_sha256": PATCH, "started_utc": interval[0], "ended_utc": interval[1]}
    fields.update(overrides)
    return provenance(**fields)


def _costs(*, load: LoadAverage = QUIET, indicative: bool = False, origin: RunProvenance | None = None) -> PipelinePencilCosts:
    runs = tuple(PencilRun(generation_s=GENERATION_S, decide_s=decide, first_cut_s=FIRST_CUT_S) for decide in UNSORTED_DECIDE)
    cases = {
        key: CasePencilCost(runs=runs, calls=dict(CALLS), ramp_out_members=MEMBERS, operations=OPERATIONS, digest=DIGEST, load_before=load, load_after=load) for key in KEYS[:-1]
    }
    return PipelinePencilCosts(timing_sha256=TIMING_SHA256, indicative=indicative, cases=cases, excluded=KEYS[-1:], provenance=origin or _head_provenance(PENCIL_INTERVAL))


def _timings(median_s: float, *, load: LoadAverage | None = QUIET, digest: str = DIGEST, origin: RunProvenance | None = None) -> PipelineTimings:
    cases = {key: CaseTiming(runs_s=(median_s,) * RUNS_PER_CASE, operations=OPERATIONS, digest=digest) for key in KEYS}
    recorded = origin or _head_provenance(TIMING_INTERVAL)
    return PipelineTimings(cases=cases, excluded=(), load_start=load, load_end=load, provenance=None if load is None else recorded)


def _base(*, same_session: bool = True, load: LoadAverage | None = QUIET, origin: RunProvenance | None = None) -> TimeGateBase:
    cases = {key: CaseTiming(runs_s=BASE_RUNS, operations=OPERATIONS, digest=OTHER_DIGEST) for key in KEYS}
    recorded = origin or _base_provenance()
    timings = PipelineTimings(cases=cases, excluded=(), load_start=load, load_end=load, provenance=None if load is None else recorded)
    source = "same-session measurement `synthetic`" if same_session else "committed record `synthetic`"
    return TimeGateBase(source=source, sha256=BASE_SHA256, timings=timings, same_session=same_session)


def _gate(median_s: float) -> PipelineTimeGate:
    return PipelineTimeGate.build(_base(), _timings(median_s), TIMING_SHA256, _costs(), expected_base_rev=BASE_REV, expected_base_tool_sha256=BASE_TOOL_SHA, indicative=False)


def test_the_fixture_is_exact_at_the_edge() -> None:
    assert _costs().cases[KEYS[0]].median_decide_s == C_PENCIL
    assert BOUND == EDGE_MEDIAN


def test_the_gate_holds_where_the_median_is_within_the_slowest_1_1_run_plus_c_pencil() -> None:
    gate = _gate(HOLDING_MEDIAN)
    assert gate.mode is TimeGateMode.VERDICT
    assert [case.key for case in gate.cases] == list(KEYS[:-1])
    for case in gate.cases:
        assert case.base_slowest_s == BASE_SLOWEST
        assert case.c_pencil_s == C_PENCIL
        assert case.holds is True


def test_a_median_exactly_at_the_bound_holds() -> None:
    gate = _gate(EDGE_MEDIAN)
    assert all(case.bound_s == case.median_s and case.holds is True for case in gate.cases)


def test_the_gate_fails_where_the_median_exceeds_the_bound() -> None:
    gate = _gate(FAILING_MEDIAN)
    assert all(case.holds is False for case in gate.cases)
    assert "fails" in render_time_gate_table(gate)


def test_a_loaded_pencil_record_yields_no_verdict() -> None:
    with pytest.raises(MachineNotQuietError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(load=LOADED),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_an_indicative_pencil_record_yields_no_verdict() -> None:
    with pytest.raises(MachineNotQuietError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(indicative=True),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_loaded_timing_record_yields_no_verdict() -> None:
    with pytest.raises(MachineNotQuietError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN, load=LOADED),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_timing_record_without_load_yields_no_verdict() -> None:
    with pytest.raises(UnrecordedLoadError) as refused:
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN, load=None),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )
    assert "python -m tools.pipeline_time_gate --indicative" in str(refused.value)


def test_an_indicative_gate_is_rendered_but_labelled_not_a_verdict_and_says_why() -> None:
    gate = PipelineTimeGate.build(
        _base(),
        _timings(HOLDING_MEDIAN, load=None),
        TIMING_SHA256,
        _costs(load=LOADED),
        expected_base_rev=BASE_REV,
        expected_base_tool_sha256=BASE_TOOL_SHA,
        indicative=True,
    )
    assert gate.mode is TimeGateMode.INDICATIVE
    text = render_time_gate_table(gate)
    assert "INDICATIVE, not a verdict" in text
    assert "the 1.2 timing record holds no load or provenance" in text
    assert "the pencil record's 1-minute load reached 3.00" in text
    assert "holds (indicative)" in text


def test_a_pencil_record_of_another_timing_document_is_refused() -> None:
    with pytest.raises(StalePencilCostError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            BASE_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_pencil_record_whose_digest_is_not_the_timing_digest_is_refused() -> None:
    with pytest.raises(StalePencilCostError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN, digest=OTHER_DIGEST),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_base_without_a_case_is_refused() -> None:
    base = _base()
    partial = TimeGateBase(
        source=base.source,
        sha256=base.sha256,
        timings=PipelineTimings(cases={KEYS[0]: base.timings.cases[KEYS[0]]}, excluded=(), load_start=QUIET, load_end=QUIET, provenance=_base_provenance()),
        same_session=True,
    )
    with pytest.raises(IncompleteTimeGateInputError):
        PipelineTimeGate.build(
            partial,
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_without_a_pencil_record_the_gate_is_pending_and_says_why() -> None:
    gate = PipelineTimeGate.build(
        _base(),
        _timings(HOLDING_MEDIAN, load=None),
        TIMING_SHA256,
        None,
        expected_base_rev=BASE_REV,
        expected_base_tool_sha256=BASE_TOOL_SHA,
        indicative=False,
    )
    assert gate.mode is TimeGateMode.PENDING
    assert all(case.holds is None and case.c_pencil_s is None for case in gate.cases)
    text = render_time_gate_table(gate)
    assert "Verdict pending" in text
    assert "1.2 timing load: unrecorded" in text
    assert "holds" not in text and "fails" not in text


def test_the_gate_rendering_states_its_inputs_both_loads_and_what_c_pencil_excuses() -> None:
    text = render_time_gate_table(_gate(HOLDING_MEDIAN))
    assert text == render_time_gate_table(_gate(HOLDING_MEDIAN))
    assert text.startswith("<!-- Generated by tools.pipeline_time_gate")
    assert "python -m tools.pipeline_time_gate" in text
    for phrase in ("same-session measurement `synthetic`", BASE_SHA256, TIMING_SHA256, "slowest", "holds", "0.50 / 0.75 / 1.25 at the start", "credited"):
        assert phrase in text
    assert "nothing else running" not in text


def test_the_committed_1_1_base_is_the_default_and_read_from_git() -> None:
    base = committed_time_gate_base()
    assert COMMITTED_1_1_REV == "7b001d47"
    assert base.source == f"committed record `{COMMITTED_1_1_REV}:{PIPELINE_TIMING_PATH.as_posix()}`"
    assert max(base.timings.cases[SMOKE_KEY].runs_s) == COMMITTED_1_1_SLOWEST_120_RECT_12X8
    assert base.timings.load_start is None


def test_the_cli_uses_the_committed_1_1_base_unless_another_is_named(tmp_path: Path) -> None:
    out = tmp_path / "gate.md"
    assert main(["--pending", "--write", str(out)]) == 0
    assert f"committed record `{COMMITTED_1_1_REV}:" in out.read_text()


def test_a_session_base_names_its_file_and_content_hash(tmp_path: Path) -> None:
    path = tmp_path / "timing-1.1.json"
    committed = committed_time_gate_base()
    path.write_text(committed.timings.dumps())
    base = session_time_gate_base(path, root=tmp_path)
    assert base.timings == committed.timings
    assert base.source == "same-session measurement `timing-1.1.json`"
    assert base.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()


def test_the_committed_1_1_base_yields_no_verdict() -> None:
    with pytest.raises(CommittedBaseError):
        PipelineTimeGate.build(
            _base(same_session=False, load=None),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_with_the_committed_base_the_gate_renders_indicative_and_says_its_load_is_unrecorded() -> None:
    gate = PipelineTimeGate.build(
        _base(same_session=False, load=None),
        _timings(HOLDING_MEDIAN),
        TIMING_SHA256,
        _costs(),
        expected_base_rev=BASE_REV,
        expected_base_tool_sha256=BASE_TOOL_SHA,
        indicative=True,
    )
    text = render_time_gate_table(gate)
    assert gate.mode is TimeGateMode.INDICATIVE
    assert "not a same-session re-timing" in text
    assert "1.1 base load: unrecorded" in text


def test_a_same_session_base_without_load_yields_no_verdict() -> None:
    with pytest.raises(UnrecordedLoadError):
        PipelineTimeGate.build(
            _base(load=None),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_loaded_same_session_base_yields_no_verdict() -> None:
    with pytest.raises(MachineNotQuietError):
        PipelineTimeGate.build(
            _base(load=LOADED),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_timing_record_loaded_mid_run_yields_no_verdict() -> None:
    with pytest.raises(MachineNotQuietError) as refused:
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN, load=LOADED_MID_RUN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )
    assert "15-minute load at the end reached 2.50" in str(refused.value)


def test_the_committed_gate_page_is_the_render_of_the_committed_records() -> None:
    """The page as `python -m tools.pipeline_time_gate` writes it from the default committed base: pending without a pencil record, indicative with one."""
    pencil = ROOT / PIPELINE_PENCIL_COST_PATH
    costs = load_pipeline_pencil_costs(pencil) if pencil.exists() else None
    timing_bytes = (ROOT / PIPELINE_TIMING_PATH).read_bytes()
    timings = PipelineTimings.build(json.loads(timing_bytes))
    gate = PipelineTimeGate.build(
        committed_time_gate_base(),
        timings,
        hashlib.sha256(timing_bytes).hexdigest(),
        costs,
        expected_base_rev=BASE_REV,
        expected_base_tool_sha256=BASE_TOOL_SHA,
        indicative=costs is not None,
    )
    assert (ROOT / PIPELINE_TIME_GATE_TABLE_PATH).read_text() == render_time_gate_table(gate)


def test_the_back_to_back_window_is_fifteen_minutes() -> None:
    assert BACK_TO_BACK_WINDOW_S == 15 * 60


def test_a_base_of_another_commit_yields_no_verdict() -> None:
    with pytest.raises(BaseRevisionMismatchError):
        PipelineTimeGate.build(
            _base(origin=_base_provenance(rev=HEAD_REV)),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_modified_base_yields_no_verdict() -> None:
    with pytest.raises(DirtyBaseError):
        PipelineTimeGate.build(
            _base(origin=_base_provenance(dirty=True, patch_sha256=PATCH)),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


@pytest.mark.parametrize("field", ["host", "cpu_model"])
def test_records_from_another_machine_yield_no_verdict(field: str) -> None:
    other = _head_provenance(PENCIL_INTERVAL, **{field: "elsewhere"})
    with pytest.raises(HostMismatchError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(origin=other),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_records_from_another_lock_yield_no_verdict() -> None:
    with pytest.raises(LockMismatchError):
        PipelineTimeGate.build(
            _base(origin=_base_provenance(pixi_lock_sha256="d" * 64)),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_records_pair_only_on_an_equal_commit_dirty_flag_and_digest() -> None:
    """The pairing rule: (rev, dirty, digest) all equal. Two export records of one content whose
    containing repository committed between them (same digest, another parent rev) are refused, the
    safe direction: re-measure back to back. A worktree and an export record cannot pair, since their
    digests are domain-separated (``tools.run_provenance.EXPORT_DIGEST_DOMAIN``)."""
    moved = _head_provenance(PENCIL_INTERVAL, rev="d" * 40)
    with pytest.raises(PencilTreeMismatchError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(origin=moved),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_a_pencil_record_of_another_tree_than_the_1_2_timing_yields_no_verdict() -> None:
    other = _head_provenance(PENCIL_INTERVAL, patch_sha256="e" * 64)
    with pytest.raises(PencilTreeMismatchError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(origin=other),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


@pytest.mark.parametrize(
    ("timing_interval", "pencil_interval"),
    [
        (("2026-10-22T18:05:00+00:00", "2026-10-22T18:30:00+00:00"), PENCIL_INTERVAL),
        (("2026-10-22T18:25:01+00:00", "2026-10-22T18:30:00+00:00"), PENCIL_INTERVAL),
        (TIMING_INTERVAL, ("2026-10-22T18:29:00+00:00", "2026-10-22T18:40:00+00:00")),
        (TIMING_INTERVAL, ("2026-10-22T18:45:01+00:00", "2026-10-22T18:50:00+00:00")),
    ],
)
def test_records_not_back_to_back_yield_no_verdict(timing_interval: tuple[str, str], pencil_interval: tuple[str, str]) -> None:
    timings = _timings(HOLDING_MEDIAN, origin=_head_provenance(timing_interval))
    with pytest.raises(NotBackToBackError):
        PipelineTimeGate.build(
            _base(),
            timings,
            TIMING_SHA256,
            _costs(origin=_head_provenance(pencil_interval)),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_gaps_of_exactly_the_window_are_back_to_back() -> None:
    timings = _timings(HOLDING_MEDIAN, origin=_head_provenance(("2026-10-22T18:25:00+00:00", "2026-10-22T18:30:00+00:00")))
    costs = _costs(origin=_head_provenance(("2026-10-22T18:45:00+00:00", "2026-10-22T18:50:00+00:00")))
    assert (
        PipelineTimeGate.build(_base(), timings, TIMING_SHA256, costs, expected_base_rev=BASE_REV, expected_base_tool_sha256=BASE_TOOL_SHA, indicative=False).mode
        is TimeGateMode.VERDICT
    )


def test_a_pencil_record_loaded_over_15_minutes_yields_no_verdict() -> None:
    with pytest.raises(MachineNotQuietError):
        PipelineTimeGate.build(
            _base(),
            _timings(HOLDING_MEDIAN),
            TIMING_SHA256,
            _costs(load=LOADED_MID_RUN),
            expected_base_rev=BASE_REV,
            expected_base_tool_sha256=BASE_TOOL_SHA,
            indicative=False,
        )


def test_an_indicative_render_names_the_identity_breaches() -> None:
    other = _head_provenance(PENCIL_INTERVAL, host="elsewhere")
    gate = PipelineTimeGate.build(
        _base(),
        _timings(HOLDING_MEDIAN),
        TIMING_SHA256,
        _costs(origin=other),
        expected_base_rev=BASE_REV,
        expected_base_tool_sha256=BASE_TOOL_SHA,
        indicative=True,
    )
    text = render_time_gate_table(gate)
    assert "INDICATIVE, not a verdict" in text and "elsewhere" in text


def test_the_verdict_page_states_each_records_provenance() -> None:
    text = render_time_gate_table(_gate(HOLDING_MEDIAN))
    for phrase in (BASE_REV, HEAD_REV, BASE_INTERVAL[0], TIMING_INTERVAL[0], PENCIL_INTERVAL[0], "host.local"):
        assert phrase in text


def test_a_session_base_outside_the_repository_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "timing-1.1.json"
    path.write_text(committed_time_gate_base().timings.dumps())
    with pytest.raises(SessionBaseOutsideRepositoryError):
        session_time_gate_base(path, root=tmp_path / "repository")


def _build(base: TimeGateBase, costs: PipelinePencilCosts, *, tool: str = BASE_TOOL_SHA) -> PipelineTimeGate:
    return PipelineTimeGate.build(base, _timings(HOLDING_MEDIAN), TIMING_SHA256, costs, expected_base_rev=BASE_REV, expected_base_tool_sha256=tool, indicative=False)


def test_a_base_built_from_the_1_2_library_yields_no_verdict() -> None:
    """The measured native sources differ (records' ``native_sources_sha256``), yet the base ran the 1.2 extensions."""
    with pytest.raises(ForeignBaseLibraryError):
        _build(_base(origin=_base_provenance(native_sha256="a" * 64)), _costs())


def test_an_equal_extension_digest_is_admitted_when_the_measured_native_sources_are_equal() -> None:
    same = _base_provenance(native_sha256="a" * 64, native_sources_sha256="a" * 64)
    assert _build(_base(origin=same), _costs()).mode is TimeGateMode.VERDICT


def test_a_pencil_record_of_another_native_build_than_the_1_2_timing_yields_no_verdict() -> None:
    """The 1.2 timing and pencil records must have run one native build; the base runs its own and is not compared with them."""
    with pytest.raises(NativeLibraryMismatchError):
        _build(_base(), _costs(origin=_head_provenance(PENCIL_INTERVAL, native_sha256="e" * 64)))


@pytest.mark.parametrize("field", ["environment", "python_version"])
def test_records_from_another_environment_yield_no_verdict(field: str) -> None:
    with pytest.raises(EnvironmentMismatchError):
        _build(_base(), _costs(origin=_head_provenance(PENCIL_INTERVAL, **{field: "other"})))


def test_a_base_from_another_harness_yields_no_verdict() -> None:
    with pytest.raises(BaseToolMismatchError):
        _build(_base(), _costs(), tool="8" * 64)
