"""The CPU/wall clock of the Held generation budget: what it counts and what it leaves out."""

from __future__ import annotations

import dataclasses
import inspect
import time
from dataclasses import dataclass

import pytest

from benchmarks import held_omega_toolpath as omega_toolpath
from benchmarks import held_toolpath_benchmark as toolpath_benchmark
from benchmarks.held_omega_toolpath import PathDomain
from tools.held_generation_clock import OMEGA_CONNECTOR_GATE_STAGE
from tools.held_generation_clock import OMEGA_STAGE_CLOCK_ATTRIBUTE
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.held_generation_clock import ClockedGeneration
from tools.held_generation_clock import CpuSummary
from tools.held_generation_clock import GenerationClock
from tools.held_generation_clock import MachineNotQuietError
from tools.held_generation_clock import NoClockedRunsError
from tools.held_generation_clock import StageClockCaptureError
from tools.held_generation_clock import clock_generation
from tools.held_generation_clock import clock_omega_generation
from tools.held_generation_clock import describe_indicative
from tools.held_generation_clock import measure_gated_runs
from tools.held_generation_clock import measure_indicative_runs
from tools.held_generation_clock import require_quiet_machine
from tools.held_generation_clock import worker_command
from tools.held_path_digest import WorkloadDigest

# CPU seconds each fake stage burns; far above process_time's resolution (about 1 µs on macOS).
BURN_S = 0.02
# CPU seconds the fake case load burns before the first stage. Four stage burns, so
# a leak of the case load into the counted total overshoots the bound below by
# CASE_LOAD_BURN_S / 2 = 2 stage burns, and a correct clock has the same headroom
# under it for the clock's own overhead.
CASE_LOAD_BURN_S = 4 * BURN_S


def _counted_total_bound(counted_burns: int) -> float:
    """Upper bound on a correct clock's CPU over ``counted_burns`` counted burns: half the case load above them."""
    return counted_burns * BURN_S + CASE_LOAD_BURN_S / 2


# A well-formed digest record for fake runs: the content is irrelevant to the clock.
_FAKE_DIGEST = {"initial": "0" * 64, "circles": "1" * 64, "connectors": "2" * 64, "circle_count": 3, "connector_count": 2, "seam": "closed"}


def _burn(seconds: float) -> None:
    started = time.process_time()
    while time.process_time() - started < seconds:
        pass


def _fake_generator() -> str:
    _burn(CASE_LOAD_BURN_S)  # case load: before the first stage, not counted
    stages = toolpath_benchmark._Stages()
    stages.start()
    _burn(BURN_S)
    stages.stop("guide")
    stages.start()
    _burn(BURN_S)
    stages.stop("initial")
    _burn(BURN_S)  # the exit collect of paused_cyclic_gc: after the last stage, counted
    return "path"


def test_cpu_counts_every_stage_and_the_exit_tail_but_not_the_case_load() -> None:
    result, clock = clock_generation(_fake_generator)
    assert result == "path"
    assert list(clock.stage_cpu_seconds) == ["guide", "initial"]
    assert all(seconds >= BURN_S for seconds in clock.stage_cpu_seconds.values())
    assert clock.exit_cpu_seconds >= BURN_S
    assert clock.cpu_seconds == pytest.approx(sum(clock.stage_cpu_seconds.values()) + clock.exit_cpu_seconds)
    # Two stages and the exit tail are counted; a leak of the case load would add CASE_LOAD_BURN_S.
    assert clock.cpu_seconds < _counted_total_bound(3)


def test_the_stage_clock_is_restored_after_the_run() -> None:
    original = toolpath_benchmark._Stages
    clock_generation(_fake_generator)
    assert toolpath_benchmark._Stages is original


def test_a_run_without_exactly_one_stage_clock_is_refused() -> None:
    original = toolpath_benchmark._Stages
    with pytest.raises(StageClockCaptureError):
        clock_generation(lambda: None)
    assert toolpath_benchmark._Stages is original


def test_the_clock_round_trips_through_json() -> None:
    _, clock = clock_generation(_fake_generator)
    assert GenerationClock.from_json(clock.to_json()) == clock


def test_the_gate_reads_the_minimum_and_reports_median_and_spread() -> None:
    summary = CpuSummary.build([2.5, 2.1, 3.0, 2.2, 2.4])
    assert summary.gated == 2.1
    assert summary.minimum == 2.1
    assert summary.median == 2.4
    assert summary.maximum == 3.0
    assert summary.spread == pytest.approx(0.9)


def test_an_empty_run_set_is_refused() -> None:
    with pytest.raises(NoClockedRunsError):
        CpuSummary.build([])


# Fake 1-minute load averages either side of the quiet-machine limit (user ruling D2, 2026-10-22: > 2 refuses).
LOADED = (2.5, 1.0, 1.0)
QUIET = (2.0, 9.0, 9.0)


def _never_measure(workload: str) -> ClockedGeneration:
    raise AssertionError(f"{workload}: measured on a machine that is not quiet")


def test_a_loaded_machine_is_refused_by_name_before_measuring() -> None:
    with pytest.raises(MachineNotQuietError, match="2.5"):
        measure_gated_runs("figure5_refined", 5, load_reader=lambda: LOADED, measure=_never_measure)


def test_the_quiet_limit_is_inclusive_and_reads_the_one_minute_load() -> None:
    assert QUIET_LOAD_AVERAGE_MAX == 2.0
    require_quiet_machine(lambda: QUIET)
    with pytest.raises(MachineNotQuietError):
        require_quiet_machine(lambda: LOADED)


def test_a_load_rising_mid_sequence_refuses_the_remaining_runs() -> None:
    loads = iter([QUIET, LOADED])
    measured: list[str] = []

    def measure(workload: str) -> ClockedGeneration:
        measured.append(workload)
        return _fake_run(workload)

    with pytest.raises(MachineNotQuietError):
        measure_gated_runs("figure5_refined", 2, load_reader=lambda: next(loads), measure=measure)
    assert measured == ["figure5_refined"]


def test_gated_runs_record_the_load_before_each_run() -> None:
    runs = measure_gated_runs("figure5_refined", 2, load_reader=lambda: QUIET, measure=_fake_run)
    assert [run.load_before for run in runs] == [QUIET, QUIET]
    assert all(not run.indicative for run in runs)


def test_indicative_runs_measure_under_load_and_are_marked_indicative() -> None:
    runs = measure_indicative_runs("figure5_refined", 2, load_reader=lambda: LOADED, measure=_fake_run)
    assert [run.load_before for run in runs] == [LOADED, LOADED]
    assert all(run.indicative for run in runs)


def test_the_indicative_report_never_states_a_budget_verdict() -> None:
    runs = measure_indicative_runs("figure5_refined", 2, load_reader=lambda: LOADED, measure=_fake_run)
    text = describe_indicative(runs)
    assert "INDICATIVE" in text
    assert "2.5" in text
    for verdict in ("pass", "fail", "budget"):
        assert verdict not in text.lower()


def _fake_run(workload: str) -> ClockedGeneration:
    _, clock = clock_generation(_fake_generator)
    digest = WorkloadDigest.from_raw(workload, _FAKE_DIGEST)
    return ClockedGeneration(workload=workload, digest=digest, clock=clock, worker_cpu_seconds=clock.cpu_seconds)


# --- the Ω* generator (benchmarks.held_omega_toolpath.generate_omega) -------------------------
# It times its stages by calling its module's ``time.perf_counter`` in start/stop pairs, so the
# fake below does the same through that module attribute, which the clock replaces for the run.


@dataclass(frozen=True)
class _FakeOmegaPath:
    stage_seconds: dict[str, float]


def _omega_read() -> float:
    """One stage-clock read through the Ω* module's clock attribute, as the generator makes it."""
    return float(vars(omega_toolpath)[OMEGA_STAGE_CLOCK_ATTRIBUTE].perf_counter())


def _omega_module_clock() -> object:
    return vars(omega_toolpath)[OMEGA_STAGE_CLOCK_ATTRIBUTE]


def _fake_omega_generator() -> _FakeOmegaPath:
    _burn(CASE_LOAD_BURN_S)  # case load: before the first stage, not counted
    seconds: dict[str, float] = {}
    started = _omega_read()
    _burn(BURN_S)
    seconds["guide"] = _omega_read() - started
    started = _omega_read()
    _burn(BURN_S)
    seconds["initial"] = _omega_read() - started
    started = _omega_read()
    _burn(BURN_S)
    seconds[OMEGA_CONNECTOR_GATE_STAGE] = _omega_read() - started
    _burn(BURN_S)  # the exit collect of paused_cyclic_gc: after the last stage, counted
    return _FakeOmegaPath(seconds)


def test_omega_cpu_counts_every_stage_and_the_exit_tail_but_not_the_case_load() -> None:
    path, clock = clock_omega_generation(_fake_omega_generator)
    assert list(clock.stage_cpu_seconds) == ["guide", "initial", OMEGA_CONNECTOR_GATE_STAGE]
    assert all(seconds >= BURN_S for seconds in clock.stage_cpu_seconds.values())
    assert clock.exit_cpu_seconds >= BURN_S
    assert clock.cpu_seconds == pytest.approx(sum(clock.stage_cpu_seconds.values()) + clock.exit_cpu_seconds)
    # Three stages and the exit tail are counted; a leak of the case load would add CASE_LOAD_BURN_S.
    assert clock.cpu_seconds < _counted_total_bound(4)
    assert clock.stage_wall_seconds == path.stage_seconds


def test_the_omega_connector_gate_stage_is_the_generators_connector_repair_stage() -> None:
    """Ω*'s connector gate and repair runs inside the generator; the clock counts it under the generator's own stage name."""
    assert OMEGA_CONNECTOR_GATE_STAGE == "connector_repair"
    assert f'seconds["{OMEGA_CONNECTOR_GATE_STAGE}"] = time.perf_counter() - started' in inspect.getsource(omega_toolpath._generate_omega)


def test_an_omega_run_without_its_connector_gate_stage_is_refused() -> None:
    def ungated() -> _FakeOmegaPath:
        _burn(CASE_LOAD_BURN_S)
        started = _omega_read()
        _burn(BURN_S)
        return _FakeOmegaPath({"guide": _omega_read() - started})

    with pytest.raises(StageClockCaptureError, match=OMEGA_CONNECTOR_GATE_STAGE):
        clock_omega_generation(ungated)


def test_the_omega_module_clock_is_restored_after_the_run_and_after_an_error() -> None:
    original = _omega_module_clock()
    clock_omega_generation(_fake_omega_generator)
    assert _omega_module_clock() is original

    def failing() -> _FakeOmegaPath:
        _omega_read()
        raise RuntimeError("generator failed")

    with pytest.raises(RuntimeError, match="generator failed"):
        clock_omega_generation(failing)
    assert _omega_module_clock() is original


def test_an_omega_run_whose_clock_reads_do_not_pair_with_its_stages_is_refused() -> None:
    def extra_read() -> _FakeOmegaPath:
        path = _fake_omega_generator()
        _omega_read()
        return path

    with pytest.raises(StageClockCaptureError, match="pair"):
        clock_omega_generation(extra_read)


def test_an_omega_run_whose_stage_seconds_are_not_its_clock_reads_is_refused() -> None:
    def forged() -> _FakeOmegaPath:
        path = _fake_omega_generator()
        return _FakeOmegaPath({name: seconds + 1.0 for name, seconds in path.stage_seconds.items()})

    with pytest.raises(StageClockCaptureError, match="guide"):
        clock_omega_generation(forged)


def test_an_omega_run_without_stages_is_refused() -> None:
    with pytest.raises(StageClockCaptureError):
        clock_omega_generation(lambda: _FakeOmegaPath({}))


def test_the_worker_command_names_the_domain() -> None:
    assert worker_command("figure5_refined", PathDomain.POLYGON)[-3:] == ["--domain", "polygon", "figure5_refined"]
    assert worker_command("figure5_refined", PathDomain.OMEGA)[-3:] == ["--domain", "omega", "figure5_refined"]


def test_an_omega_indicative_report_names_its_domain() -> None:
    def measure(workload: str) -> ClockedGeneration:
        return dataclasses.replace(_fake_run(workload), domain=PathDomain.OMEGA)

    runs = measure_indicative_runs("figure5_refined", 2, load_reader=lambda: LOADED, measure=measure)
    assert "omega" in describe_indicative(runs).splitlines()[0]
