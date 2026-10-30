"""The pencil decision-cost tool.

One fresh process per run, runs agree, the digest equals the uninstrumented
one, the wrappers decide nothing, every wrapped callable is counted on its own,
and the load is recorded; a loaded machine measures nothing unless indicative.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Iterator

import pytest

from tests.tools.test_run_provenance import provenance

from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket
from compas_cgal.pipeline import generate
from compas_cgal.pipeline.pencil import RampOut
from compas_cgal.pipeline.pencil import TangentPencil
from compas_cgal.pipeline.result import PipelineResult
from tools.held_native_identity import NativeBuildChangedError
from tools.run_native_identity import NativeFinisher
from tools.run_native_identity import NativeIdentity
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import LoadAverage
from tools.held_generation_clock import MachineNotQuietError
from tools.pipeline_pencil_cost import RESULT_LINE_PREFIX
from tools.pipeline_pencil_cost import InstrumentationNotReachedError
from tools.pipeline_pencil_cost import MissingTimingReferenceError
from tools.pipeline_pencil_cost import PencilCaseRun
from tools.pipeline_pencil_cost import PencilDigestMismatchError
from tools.pipeline_pencil_cost import UnknownPencilCaseError
from tools.pipeline_pencil_cost import launch_fresh_process
from tools.pipeline_pencil_cost import measure_pipeline_pencil_costs
from tools.pipeline_pencil_cost import pencil_instrumentation
from tools.pipeline_pencil_cost import result_payload
from tools.pipeline_pencil_cost import rooted
from tools.pipeline_pencil_cost import run_case_here
from tools.pipeline_pencil_cost_results import DECIDE_CALLABLES
from tools.pipeline_pencil_cost_results import ENTER_FIRST_CUT
from tools.pipeline_pencil_cost_results import RAMP_OUT_BUILD
from tools.pipeline_pencil_cost_results import RAMP_OUT_NEXT_MEMBER
from tools.pipeline_pencil_cost_results import RENDERED_SECONDS_DECIMALS
from tools.pipeline_pencil_cost_results import TANGENT_PENCIL_BUILD
from tools.pipeline_pencil_cost_results import WRAPPED_CALLABLES
from tools.pipeline_pencil_cost_results import CallCompositionError
from tools.pipeline_pencil_cost_results import CasePencilCost
from tools.pipeline_pencil_cost_results import MalformedPencilCostError
from tools.pipeline_pencil_cost_results import PencilRun
from tools.pipeline_pencil_cost_results import PipelinePencilCosts
from tools.pipeline_pencil_cost_results import load_pipeline_pencil_costs
from tools.pipeline_pencil_cost_results import pencil_case_keys
from tools.pipeline_pencil_cost_results import render_pipeline_pencil_cost_table
from tools.pipeline_timing import CaseProcessFailedError
from tools.pipeline_timing import MissingCaseResultError
from tools.pipeline_timing import NondeterministicCaseError
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import load_pipeline_timings
from tools.run_provenance import RunProvenance
from tools.run_provenance import TreeChangedDuringRunError
from tools.run_provenance import TreeIdentity

ROOT = Path(__file__).resolve().parents[2]
KEYS = pencil_case_keys()
SMOKE_KEY = "cap-120-default/pipeline/rect_12x8"
SMOKE_CAP = GateCapDegrees(120.0)
SMOKE_POCKET = "rect_12x8"
NON_PIPELINE_KEY = "cap-120-default/rho_regulated/rect_12x8"
DIGEST = "0123456789abcdef"
OTHER_DIGEST = "fedcba9876543210"
TIMING_SHA256 = "0" * 64
OPERATIONS = 7
MEMBERS = 4
CALLS = {TANGENT_PENCIL_BUILD: 1, RAMP_OUT_BUILD: 1, RAMP_OUT_NEXT_MEMBER: MEMBERS + 1, ENTER_FIRST_CUT: 1}
UNSORTED_DECIDE = (0.3, 0.1, 0.2)
MEDIAN_OF_UNSORTED_DECIDE = 0.2
FIRST_CUT_S = 0.5
GENERATION_S = 1.0
QUIET: LoadAverage = (0.5, 0.75, 1.25)
LOADED: LoadAverage = (QUIET_LOAD_AVERAGE_MAX + 1.0, 3.5, 4.5)


def _run(decide_s: float = 0.1, *, digest: str = DIGEST, operations: int = OPERATIONS, calls: dict[str, int] | None = None) -> PencilCaseRun:
    return PencilCaseRun(
        run=PencilRun(generation_s=GENERATION_S, decide_s=decide_s, first_cut_s=FIRST_CUT_S),
        calls=dict(CALLS if calls is None else calls),
        ramp_out_members=MEMBERS,
        operations=operations,
        digest=digest,
    )


def _reference(keys: tuple[str, ...] = KEYS) -> dict[str, str]:
    return {key: DIGEST for key in keys}


TREE = ("b" * 40, True, "c" * 64)


def _tree() -> TreeIdentity:
    return TREE


NATIVE = NativeIdentity("1" * 64, "2" * 64)


def _track_native() -> NativeFinisher:
    return lambda: NATIVE


def _native_build_changed() -> NativeFinisher:
    def finish() -> NativeIdentity:
        raise NativeBuildChangedError("a native extension changed on disk during the run")

    return finish


def _describe(tree: TreeIdentity, native: NativeIdentity, started_utc: str, ended_utc: str) -> RunProvenance:
    return provenance(
        rev=tree[0],
        dirty=tree[1],
        patch_sha256=tree[2],
        native_sha256=native.native_sha256,
        native_sources_sha256=native.native_sources_sha256,
        started_utc=started_utc,
        ended_utc=ended_utc,
    )


def _quiet() -> LoadAverage:
    return QUIET


def _loads(*loads: LoadAverage) -> Iterator[LoadAverage]:
    return iter(loads)


def _synthetic(*, load: LoadAverage = QUIET, indicative: bool = False) -> PipelinePencilCosts:
    runs = tuple(PencilRun(generation_s=GENERATION_S, decide_s=decide, first_cut_s=FIRST_CUT_S) for decide in UNSORTED_DECIDE)
    cases = {
        key: CasePencilCost(runs=runs, calls=dict(CALLS), ramp_out_members=MEMBERS, operations=OPERATIONS, digest=DIGEST, load_before=load, load_after=load) for key in KEYS[:-1]
    }
    return PipelinePencilCosts(timing_sha256=TIMING_SHA256, indicative=indicative, cases=cases, excluded=KEYS[-1:], provenance=provenance())


# --- measurement: serial, fresh process, agreement, digest -------------------------------------------


def test_the_pencil_cases_are_the_pipeline_generator_cases_in_gate_order() -> None:
    assert SMOKE_KEY in KEYS
    assert all(key.split("/")[1] == "pipeline" for key in KEYS)
    assert NON_PIPELINE_KEY not in KEYS


def test_every_case_is_launched_once_per_run_serially_case_by_case() -> None:
    launched: list[str] = []

    def launch(key: str) -> PencilCaseRun:
        launched.append(key)
        return _run()

    costs = measure_pipeline_pencil_costs(
        KEYS[:-1],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=KEYS[-1:],
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=launch,
        load_reader=_quiet,
    )
    assert launched == [key for key in KEYS[:-1] for _run_index in range(RUNS_PER_CASE)]
    assert tuple(costs.cases) == KEYS[:-1]
    assert costs.excluded == KEYS[-1:]
    assert PipelinePencilCosts.build(costs.to_json()) == costs


def test_the_low_median_of_three_decision_costs_is_reported() -> None:
    runs = iter(UNSORTED_DECIDE)
    costs = measure_pipeline_pencil_costs(
        KEYS[:1],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=(),
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=lambda _key: _run(next(runs)),
        load_reader=_quiet,
    )
    (case,) = costs.cases.values()
    assert tuple(run.decide_s for run in case.runs) == UNSORTED_DECIDE
    assert case.median_decide_s == MEDIAN_OF_UNSORTED_DECIDE


def test_runs_that_disagree_on_the_digest_are_refused() -> None:
    digests = iter((DIGEST, DIGEST, OTHER_DIGEST))
    with pytest.raises(NondeterministicCaseError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(digest=next(digests)),
            load_reader=_quiet,
        )


def test_runs_that_disagree_on_one_callables_count_are_refused() -> None:
    other = {**CALLS, RAMP_OUT_NEXT_MEMBER: CALLS[RAMP_OUT_NEXT_MEMBER] + 1}
    calls = iter((CALLS, other, CALLS))
    with pytest.raises(NondeterministicCaseError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(calls=next(calls)),
            load_reader=_quiet,
        )


def test_agreeing_runs_off_the_uninstrumented_digest_are_refused() -> None:
    with pytest.raises(PencilDigestMismatchError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(digest=OTHER_DIGEST),
            load_reader=_quiet,
        )


def test_a_case_without_a_reference_digest_is_refused_before_any_launch() -> None:
    def launch(_key: str) -> PencilCaseRun:
        raise AssertionError("launched without a reference digest")

    with pytest.raises(MissingTimingReferenceError):
        measure_pipeline_pencil_costs(
            KEYS[:1], describe=_describe, track_native=_track_native, read_tree=_tree, excluded=(), reference={}, timing_sha256=TIMING_SHA256, launch=launch, load_reader=_quiet
        )


def test_the_default_launcher_starts_a_fresh_interpreter_per_run(monkeypatch: pytest.MonkeyPatch) -> None:
    commands: list[list[str]] = []
    payload = json.dumps(result_payload(_run()))

    def run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout=f"noise\n{RESULT_LINE_PREFIX}{payload}\n", stderr="")

    monkeypatch.setattr(subprocess, "run", run)
    costs = measure_pipeline_pencil_costs(
        KEYS[:2],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=(),
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=launch_fresh_process,
        load_reader=_quiet,
    )
    assert commands == [[sys.executable, "-m", "tools.pipeline_pencil_cost", "--case", key] for key in KEYS[:2] for _run_index in range(RUNS_PER_CASE)]
    assert all(case.digest == DIGEST and dict(case.calls) == CALLS for case in costs.cases.values())


def test_a_child_without_a_result_line_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subprocess, "run", lambda command, **_kwargs: subprocess.CompletedProcess(command, 0, stdout="noise\n", stderr=""))
    with pytest.raises(MissingCaseResultError):
        launch_fresh_process(KEYS[0])


def test_a_failing_child_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subprocess, "run", lambda command, **_kwargs: subprocess.CompletedProcess(command, 1, stdout="", stderr="boom"))
    with pytest.raises(CaseProcessFailedError):
        launch_fresh_process(KEYS[0])


def test_a_non_pipeline_case_is_refused() -> None:
    with pytest.raises(UnknownPencilCaseError):
        run_case_here(NON_PIPELINE_KEY)


def test_relative_output_paths_resolve_against_the_repository_root() -> None:
    relative = Path("benchmarks") / "results" / "x.json"
    assert rooted(relative) == ROOT / relative
    assert rooted(ROOT / relative) == ROOT / relative


# --- load: recorded, and a loaded machine measures nothing unless indicative -------------------------


def test_a_loaded_machine_is_refused_before_any_launch() -> None:
    def launch(_key: str) -> PencilCaseRun:
        raise AssertionError("launched on a loaded machine")

    with pytest.raises(MachineNotQuietError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=launch,
            load_reader=lambda: LOADED,
        )


def test_the_loaded_machine_refusal_names_this_tools_indicative_mode() -> None:
    with pytest.raises(MachineNotQuietError) as refused:
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(),
            load_reader=lambda: LOADED,
        )
    assert "python -m tools.pipeline_pencil_cost --indicative" in str(refused.value)
    assert "held_generation_clock" not in str(refused.value)


def test_a_case_that_ends_on_a_loaded_machine_is_refused() -> None:
    loads = _loads(QUIET, LOADED)
    with pytest.raises(MachineNotQuietError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(),
            load_reader=lambda: next(loads),
        )


def test_the_load_before_and_after_each_case_is_recorded() -> None:
    after = (1.0, 1.5, 2.0)
    loads = _loads(QUIET, after)
    costs = measure_pipeline_pencil_costs(
        KEYS[:1],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=(),
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=lambda _key: _run(),
        load_reader=lambda: next(loads),
    )
    (case,) = costs.cases.values()
    assert (case.load_before, case.load_after) == (QUIET, after)
    assert costs.indicative is False


def test_an_indicative_run_measures_under_load_and_says_so() -> None:
    costs = measure_pipeline_pencil_costs(
        KEYS[:1],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=(),
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=lambda _key: _run(),
        load_reader=lambda: LOADED,
        indicative=True,
    )
    (case,) = costs.cases.values()
    assert case.load_before == LOADED
    assert costs.indicative is True


# --- the record --------------------------------------------------------------------------------------


def test_the_document_round_trips_through_json_loads_and_the_loader(tmp_path: Path) -> None:
    costs = _synthetic()
    path = tmp_path / "pencil_cost.json"
    path.write_text(costs.dumps())
    assert json.loads(path.read_text()) == costs.to_json()
    assert load_pipeline_pencil_costs(path) == costs


def test_a_document_with_a_missing_case_is_refused() -> None:
    raw = _synthetic().to_json()
    assert isinstance(raw["cases"], dict)
    del raw["cases"][KEYS[0]]
    with pytest.raises(MalformedPencilCostError):
        PipelinePencilCosts.build(raw)


def test_a_document_whose_decision_cost_exceeds_its_first_cut_is_refused() -> None:
    raw = _synthetic().to_json()
    assert isinstance(raw["cases"], dict)
    raw["cases"][KEYS[0]]["runs"][0]["decide_s"] = FIRST_CUT_S + 1
    with pytest.raises(MalformedPencilCostError):
        PipelinePencilCosts.build(raw)


def test_a_measured_record_carries_its_provenance() -> None:
    costs = measure_pipeline_pencil_costs(
        KEYS[:1],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=(),
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=lambda _key: _run(),
        load_reader=_quiet,
    )
    assert costs.provenance.started <= costs.provenance.ended
    assert RunProvenance.build(costs.provenance.to_json()) == costs.provenance


def test_a_document_without_provenance_is_refused() -> None:
    raw = _synthetic().to_json()
    del raw["provenance"]
    with pytest.raises(MalformedPencilCostError):
        PipelinePencilCosts.build(raw)


def test_the_highest_recorded_15_minute_load_is_reported() -> None:
    assert _synthetic(load=(0.5, 0.75, 3.0)).max_fifteen_minute_load == 3.0


def test_a_document_without_a_load_is_refused() -> None:
    raw = _synthetic().to_json()
    assert isinstance(raw["cases"], dict)
    del raw["cases"][KEYS[0]]["load_before"]
    with pytest.raises(MalformedPencilCostError):
        PipelinePencilCosts.build(raw)


@pytest.mark.parametrize("callable_name", WRAPPED_CALLABLES)
def test_a_document_missing_one_callables_count_is_refused(callable_name: str) -> None:
    raw = _synthetic().to_json()
    assert isinstance(raw["cases"], dict)
    del raw["cases"][KEYS[0]]["calls"][callable_name]
    with pytest.raises(MalformedPencilCostError):
        PipelinePencilCosts.build(raw)


def test_a_document_whose_next_member_calls_are_not_one_per_member_plus_the_root_is_refused() -> None:
    raw = _synthetic().to_json()
    assert isinstance(raw["cases"], dict)
    raw["cases"][KEYS[0]]["calls"][RAMP_OUT_NEXT_MEMBER] = MEMBERS
    with pytest.raises(CallCompositionError):
        PipelinePencilCosts.build(raw)


def test_a_document_with_two_pencil_builds_per_first_cut_is_refused() -> None:
    raw = _synthetic().to_json()
    assert isinstance(raw["cases"], dict)
    raw["cases"][KEYS[0]]["calls"][TANGENT_PENCIL_BUILD] = 2
    with pytest.raises(CallCompositionError):
        PipelinePencilCosts.build(raw)


def test_rendering_is_deterministic_and_states_the_method_and_the_measured_load() -> None:
    text = render_pipeline_pencil_cost_table(_synthetic())
    assert text == render_pipeline_pencil_cost_table(_synthetic())
    assert text.startswith("<!-- Generated by tools.pipeline_pencil_cost")
    assert "pixi run pipeline-pencil-cost" in text
    for phrase in ("C_pencil", "fresh process", f"{RUNS_PER_CASE} generations", "uninstrumented digest", TIMING_SHA256, "same session", "instrumented"):
        assert phrase in text
    for name in WRAPPED_CALLABLES:
        assert name in text
    assert "nothing else running" not in text
    assert "0.50 / 0.75 / 1.25" in text
    assert f"{MEDIAN_OF_UNSORTED_DECIDE:.{RENDERED_SECONDS_DECIMALS}f}" in text
    assert f"`{KEYS[-1]}`" in text
    assert "INDICATIVE" not in text


def test_an_indicative_record_is_labelled_indicative() -> None:
    text = render_pipeline_pencil_cost_table(_synthetic(load=LOADED, indicative=True))
    assert "INDICATIVE" in text
    assert "3.00 / 3.50 / 4.50" in text


# --- instrumentation ---------------------------------------------------------------------------------


def test_the_instrumentation_restores_every_wrapped_attribute_even_on_error() -> None:
    before = (TangentPencil.__dict__["build"], RampOut.__dict__["build"], RampOut.__dict__["next_member"], generate.__dict__["enter_first_cut"])
    with pytest.raises(InstrumentationNotReachedError):
        with pencil_instrumentation():
            during = (TangentPencil.__dict__["build"], RampOut.__dict__["build"], RampOut.__dict__["next_member"], generate.__dict__["enter_first_cut"])
            assert all(now is not then for now, then in zip(during, before))
            raise InstrumentationNotReachedError("probe")
    assert (TangentPencil.__dict__["build"], RampOut.__dict__["build"], RampOut.__dict__["next_member"], generate.__dict__["enter_first_cut"]) == before


def test_smoke_the_instrumented_case_reproduces_the_uninstrumented_digest_and_counts_each_callable() -> None:
    """Real generation of the fast cap-120 rect_12x8 case: each callable counted on its own, nested clocks, the committed digest unchanged."""
    reference = load_pipeline_timings(ROOT / PIPELINE_TIMING_PATH).cases[SMOKE_KEY]
    uninstrumented = GATE_GENERATORS["pipeline"](gate_pocket(SMOKE_POCKET, tea_cap_deg=SMOKE_CAP))
    assert isinstance(uninstrumented, PipelineResult)
    members = len(uninstrumented.ramp_out.members)
    run = run_case_here(SMOKE_KEY)
    assert run.digest == reference.digest
    assert run.operations == reference.operations
    assert run.ramp_out_members == members
    assert set(run.calls) == set(WRAPPED_CALLABLES)
    assert run.calls[TANGENT_PENCIL_BUILD] == 1
    assert run.calls[RAMP_OUT_BUILD] == 1
    assert run.calls[ENTER_FIRST_CUT] == 1
    assert run.calls[RAMP_OUT_NEXT_MEMBER] == members + 1
    assert set(DECIDE_CALLABLES) == {TANGENT_PENCIL_BUILD, RAMP_OUT_BUILD, RAMP_OUT_NEXT_MEMBER}
    assert 0 < run.run.decide_s <= run.run.first_cut_s <= run.run.generation_s


def test_smoke_a_real_fresh_process_round_trip() -> None:
    """The child mode end to end: a fresh interpreter generates the smoke case and its result line parses to the committed digest."""
    reference = load_pipeline_timings(ROOT / PIPELINE_TIMING_PATH).cases[SMOKE_KEY]
    run = launch_fresh_process(SMOKE_KEY)
    assert run.digest == reference.digest
    assert run.calls[RAMP_OUT_NEXT_MEMBER] == run.ramp_out_members + 1


def test_a_tree_that_changes_during_the_run_is_refused() -> None:
    trees = iter((TREE, ("b" * 40, True, "d" * 64)))
    with pytest.raises(TreeChangedDuringRunError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_track_native,
            read_tree=lambda: next(trees),
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(),
            load_reader=_quiet,
        )


def test_a_native_build_that_changes_during_the_pencil_run_is_refused() -> None:
    with pytest.raises(NativeBuildChangedError):
        measure_pipeline_pencil_costs(
            KEYS[:1],
            describe=_describe,
            track_native=_native_build_changed,
            read_tree=_tree,
            excluded=(),
            reference=_reference(),
            timing_sha256=TIMING_SHA256,
            launch=lambda _key: _run(),
            load_reader=_quiet,
        )


def test_a_measured_pencil_record_carries_the_native_identity_of_the_run() -> None:
    costs = measure_pipeline_pencil_costs(
        KEYS[:1],
        describe=_describe,
        track_native=_track_native,
        read_tree=_tree,
        excluded=(),
        reference=_reference(),
        timing_sha256=TIMING_SHA256,
        launch=lambda _key: _run(),
        load_reader=_quiet,
    )
    assert (costs.provenance.native_sha256, costs.provenance.native_sources_sha256) == NATIVE
