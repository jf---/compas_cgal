"""Process CPU and wall clocks of one Held generation, in a fresh interpreter.

The Held generation budget gates process CPU time; wall time is reported beside
it and not gated (user ruling 2026-10-22). Wall time grows with the load other
sessions put on the machine; CPU time grows much less, but still some (shared
caches, core assignment). Interference only ever adds time, so the gate reads
the minimum CPU over N fresh-process runs, with median and spread reported
beside it (second user ruling 2026-10-22): ``CpuSummary``.

Budget verdicts count only from quiet-machine runs (user ruling D2,
2026-10-22): ``measure_gated_runs`` reads the 1-minute load average before
every run and raises ``MachineNotQuietError`` instead of measuring when it is
above ``QUIET_LOAD_AVERAGE_MAX``. Under load, ``measure_indicative_runs`` (CLI:
``python -m tools.held_generation_clock --indicative <workload>...``) measures
the same way and prints the numbers with their load, marked indicative; it
never compares them with a budget.

What is counted: the intervals of ``GeneratedToolpath.stage_seconds`` (guide,
inward offset, initial placement, refinement, NOTHING-circle removal,
connectors, connector post-gate), clocked with ``time.process_time`` over the
same start/stop calls, plus the tail from the end of the last stage to the
generator's return. That tail is the one ``gc.collect()`` that
``paused_cyclic_gc`` runs on exit. It is counted because the budgets were set
(lever plan, 2026-10-21) before the collector was paused, when the collector's
passes ran inside the stages: the pause moves that work out of the stage
clocks, it does not remove it. Case load, before the first stage, is excluded,
as in ``generation_seconds``. Work a generator does between one stage's stop
and the next stage's start is not counted either; today that is a few cheap
checks (the Ω* start-marker check and ``EngagementCap.build``), but work moved
into such a gap leaves the gated quantity silently.

The stage clock is replaced for the run by a subclass that reads both clocks;
the generator code is not touched, and the path it produces is digested by
``tools.held_path_digest.measure_workload`` unchanged.

The Ω* generator (``benchmarks.held_omega_toolpath.generate_omega``; CLI
``--domain omega``) is clocked the same way, over the same counted quantity:
its stages (guide, Ω* domain, initial placement, refinement, NOTHING-circle
removal, connector repair) plus the exit tail, case load excluded. The stage
sets differ at the connectors on purpose (orchestrator ruling on RV19 A-I2,
2026-10-22). The polygon path certifies its connectors in a post-gate stage
(``connector_certificate``). On Ω* that post-gate,
``held_toolpath_benchmark.require_omega_connector_certificate``, refuses every
moving connector by name: it is not Ω*'s gate and is not counted. Ω*'s
connector gate and repair is its in-generator stage ``connector_repair``
(``OMEGA_CONNECTOR_GATE_STAGE``), which is counted: a run whose stages lack it
is refused. It times its
stages with its module's ``time.perf_counter`` in start/stop pairs, not with
the polygon path's stage clock, so ``clock_omega_generation`` replaces that
module attribute for the run by a clock that returns the same wall reading and
records process CPU beside it. The reads must pair with ``stage_seconds``
exactly (two per stage, and each pair's wall difference equal to the stage's
recorded seconds), or the run is refused: a change to how the generator times
its stages fails loudly here instead of mis-attributing CPU. Again the
generator code is not touched.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import resource
import statistics
import subprocess
import sys
import time
from collections.abc import Callable
from collections.abc import Iterator
from collections.abc import Mapping
from collections.abc import Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Protocol
from typing import TypeVar

from benchmarks import held_omega_toolpath as omega_toolpath
from benchmarks import held_toolpath_benchmark as toolpath_benchmark
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import PathDomain
from benchmarks.held_omega_toolpath import generate_omega
from benchmarks.performance_record import REPOSITORY
from tools.held_path_digest import MissingWorkerResultError
from tools.held_path_digest import WorkerProcessError
from tools.held_path_digest import WorkloadDigest
from tools.held_path_digest import canonical_digest
from tools.held_path_digest import measure_workload

# Prefix of the one stdout line a worker process reports its clock on.
WORKER_MARKER = "HELD_GENERATION_CLOCK_RESULT "
# Highest 1-minute load average (runnable processes) at which a run counts toward
# a budget verdict; inclusive. User ruling D2, 2026-10-22. The gate's own serial
# worker plus its parent contribute about 1.
QUIET_LOAD_AVERAGE_MAX = 2.0
# The Ω* generator module's attribute whose ``perf_counter`` times its stages.
OMEGA_STAGE_CLOCK_ATTRIBUTE = "time"
# The Ω* generator's in-generator connector gate and repair stage; the Ω* clock
# refuses a run without it (ruling on RV19 A-I2, 2026-10-22).
OMEGA_CONNECTOR_GATE_STAGE = "connector_repair"
# Fresh-process runs per workload of the indicative CLI; the gate's count.
INDICATIVE_RUNS_PER_WORKLOAD = 5

LoadAverage = tuple[float, float, float]

R = TypeVar("R")


class _StageTimed(Protocol):
    @property
    def stage_seconds(self) -> Mapping[str, float]: ...


S = TypeVar("S", bound=_StageTimed)


class StageClockCaptureError(Exception):
    """The clocked run's stage reads cannot be attributed to its stages.

    Polygon: the run did not create exactly one generator stage clock. Ω*: the
    run recorded no stage or not its connector gate stage, its clock reads do
    not pair two per stage, or a pair's wall difference is not the stage's
    recorded seconds.
    """


class NoClockedRunsError(Exception):
    """A CPU summary was asked for over no runs."""


class MachineNotQuietError(Exception):
    """The 1-minute load average is above ``QUIET_LOAD_AVERAGE_MAX``: a budget verdict would not count."""


@dataclass(frozen=True)
class CpuSummary:
    """CPU seconds over N runs of one workload: the gated minimum, and the median and spread beside it.

    Attributes:
        minimum: The fastest run; the gated value, since interference only adds time.
        median: Reported, not gated.
        maximum: Reported, not gated.
    """

    minimum: float
    median: float
    maximum: float

    @classmethod
    def build(cls, cpu_seconds: Sequence[float]) -> CpuSummary:
        """Raises:
        NoClockedRunsError: ``cpu_seconds`` is empty.
        """
        if not cpu_seconds:
            raise NoClockedRunsError("a CPU summary needs at least one run.")
        return cls(min(cpu_seconds), statistics.median(cpu_seconds), max(cpu_seconds))

    @property
    def gated(self) -> float:
        """The value the budget is compared with."""
        return self.minimum

    @property
    def spread(self) -> float:
        """Maximum minus minimum."""
        return self.maximum - self.minimum


@dataclass(frozen=True)
class GenerationClock:
    """CPU and wall seconds of one generation, stage by stage, plus the exit tail.

    Attributes:
        stage_cpu_seconds: Process CPU seconds per generator stage.
        stage_wall_seconds: Wall seconds per generator stage (``GeneratedToolpath.stage_seconds``).
        exit_cpu_seconds: Process CPU seconds from the last stage's end to the generator's return.
        exit_wall_seconds: Wall seconds of the same tail.
    """

    stage_cpu_seconds: Mapping[str, float]
    stage_wall_seconds: Mapping[str, float]
    exit_cpu_seconds: float
    exit_wall_seconds: float

    @property
    def cpu_seconds(self) -> float:
        """The gated quantity: stage CPU plus the exit tail."""
        return sum(self.stage_cpu_seconds.values()) + self.exit_cpu_seconds

    @property
    def wall_seconds(self) -> float:
        """Reported beside the gate: stage wall plus the exit tail."""
        return sum(self.stage_wall_seconds.values()) + self.exit_wall_seconds

    def to_json(self) -> str:
        return json.dumps(dataclasses.asdict(self))

    @classmethod
    def from_json(cls, text: str) -> GenerationClock:
        raw = json.loads(text)
        return cls(
            stage_cpu_seconds={str(k): float(v) for k, v in raw["stage_cpu_seconds"].items()},
            stage_wall_seconds={str(k): float(v) for k, v in raw["stage_wall_seconds"].items()},
            exit_cpu_seconds=float(raw["exit_cpu_seconds"]),
            exit_wall_seconds=float(raw["exit_wall_seconds"]),
        )


class _ClockedStages(toolpath_benchmark._Stages):
    """The generator's stage clock, also reading process CPU over the same intervals."""

    def __init__(self) -> None:
        super().__init__()
        self.cpu_seconds: dict[str, float] = {}
        self._cpu_started = 0.0
        self.last_stop_cpu = 0.0
        self.last_stop_wall = 0.0

    def start(self) -> None:
        super().start()
        self._cpu_started = time.process_time()

    def stop(self, stage: str) -> None:
        self.last_stop_cpu = time.process_time()
        self.cpu_seconds[stage] = self.last_stop_cpu - self._cpu_started
        super().stop(stage)
        self.last_stop_wall = time.perf_counter()


@contextmanager
def _clocked_stages() -> Iterator[list[_ClockedStages]]:
    created: list[_ClockedStages] = []
    original = toolpath_benchmark._Stages

    class _Recorded(_ClockedStages):
        def __init__(self) -> None:
            super().__init__()
            created.append(self)

    toolpath_benchmark._Stages = _Recorded  # type: ignore[misc]
    try:
        yield created
    finally:
        toolpath_benchmark._Stages = original  # type: ignore[misc]


def clock_generation(run: Callable[[], R]) -> tuple[R, GenerationClock]:
    """Run ``run`` (one generator call) with the stage clock reading CPU too.

    Raises:
        StageClockCaptureError: ``run`` created no stage clock, or more than one.
    """
    with _clocked_stages() as created:
        result = run()
        end_cpu = time.process_time()
        end_wall = time.perf_counter()
    if len(created) != 1:
        raise StageClockCaptureError(f"expected one generator stage clock, the run created {len(created)}.")
    stages = created[0]
    clock = GenerationClock(
        stage_cpu_seconds=dict(stages.cpu_seconds),
        stage_wall_seconds=dict(stages.seconds),
        exit_cpu_seconds=end_cpu - stages.last_stop_cpu,
        exit_wall_seconds=end_wall - stages.last_stop_wall,
    )
    return result, clock


class _PairedReadClock:
    """Stands in for a generator module's ``time``: ``perf_counter`` returns the wall reading and records CPU beside it."""

    def __init__(self) -> None:
        self.reads: list[tuple[float, float]] = []

    def perf_counter(self) -> float:
        cpu = time.process_time()
        wall = time.perf_counter()
        self.reads.append((wall, cpu))
        return wall


@contextmanager
def _clocked_omega_module() -> Iterator[_PairedReadClock]:
    clock = _PairedReadClock()
    namespace = vars(omega_toolpath)
    original = namespace[OMEGA_STAGE_CLOCK_ATTRIBUTE]
    namespace[OMEGA_STAGE_CLOCK_ATTRIBUTE] = clock
    try:
        yield clock
    finally:
        namespace[OMEGA_STAGE_CLOCK_ATTRIBUTE] = original


def clock_omega_generation(run: Callable[[], S]) -> tuple[S, GenerationClock]:
    """Run ``run`` (one Ω* generator call) with its module's stage reads also reading CPU.

    Raises:
        StageClockCaptureError: the run recorded no stage or no
            ``OMEGA_CONNECTOR_GATE_STAGE``, its clock reads do not pair two per
            recorded stage, or a pair's wall difference is not the stage's
            recorded seconds.
    """
    with _clocked_omega_module() as clock:
        result = run()
        end_cpu = time.process_time()
        end_wall = time.perf_counter()
    stages = dict(result.stage_seconds)
    reads = clock.reads
    if not stages:
        raise StageClockCaptureError("the Ω* run recorded no stage.")
    if OMEGA_CONNECTOR_GATE_STAGE not in stages:
        raise StageClockCaptureError(f"the Ω* run recorded no connector gate stage {OMEGA_CONNECTOR_GATE_STAGE!r}; its stages: {list(stages)}.")
    if len(reads) != 2 * len(stages):
        raise StageClockCaptureError(f"the Ω* run made {len(reads)} clock reads; they must pair two per stage over its {len(stages)} stages {list(stages)}.")
    stage_cpu: dict[str, float] = {}
    for index, (stage, seconds) in enumerate(stages.items()):
        (start_wall, start_cpu), (stop_wall, stop_cpu) = reads[2 * index], reads[2 * index + 1]
        if stop_wall - start_wall != seconds:
            raise StageClockCaptureError(f"Ω* stage {stage!r}: recorded {seconds!r} s, but its clock read pair spans {stop_wall - start_wall!r} s.")
        stage_cpu[stage] = stop_cpu - start_cpu
    last_wall, last_cpu = reads[-1]
    generation = GenerationClock(
        stage_cpu_seconds=stage_cpu,
        stage_wall_seconds=stages,
        exit_cpu_seconds=end_cpu - last_cpu,
        exit_wall_seconds=end_wall - last_wall,
    )
    return result, generation


def omega_digest(path: OmegaGeneratedToolpath) -> WorkloadDigest:
    """Digests of one Ω* path: its initial placement, its emitted circles and connectors.

    Raises:
        NonCanonicalReprError: a generated value has no canonical ``repr``.
    """
    connectors = path.connectors
    return WorkloadDigest(
        initial=canonical_digest(path.initial),
        circles=canonical_digest(path.circles),
        connectors=canonical_digest(connectors.connectors),
        circle_count=len(path.circles),
        connector_count=len(connectors.connectors),
        seam=connectors.seam.value,
    )


def measure_clocked_omega(workload: str) -> tuple[WorkloadDigest, GenerationClock]:
    """Generate and digest ``workload`` on Ω* in this process with both clocks.

    Raises:
        UnknownOmegaWorkloadError: ``workload`` has no Ω* generator.
    """
    path, clock = clock_omega_generation(lambda: generate_omega(workload))
    return omega_digest(path), clock


@dataclass(frozen=True)
class ClockedGeneration:
    """One fresh-process generation: its path digest, its clocks, and the whole worker's CPU.

    Attributes:
        workload: The ``held_toolpath_benchmark.GENERATORS`` key.
        digest: The digests of the generated path.
        clock: CPU and wall seconds of the generation.
        worker_cpu_seconds: CPU of the whole worker process (imports, case load,
            digesting included), from the parent's ``RUSAGE_CHILDREN``; context only.
        domain: The centre domain the path was generated on.
    """

    workload: str
    digest: WorkloadDigest
    clock: GenerationClock
    worker_cpu_seconds: float
    domain: PathDomain = PathDomain.POLYGON


def measure_clocked(workload: str) -> tuple[WorkloadDigest, GenerationClock]:
    """Generate and digest ``workload`` in this process with both clocks."""
    clocks: list[GenerationClock] = []
    original = toolpath_benchmark.generate

    def clocked(workload: str) -> toolpath_benchmark.GeneratedToolpath:
        path, clock = clock_generation(lambda: original(workload))
        clocks.append(clock)
        return path

    toolpath_benchmark.generate = clocked
    try:
        measurement = measure_workload(workload)
    finally:
        toolpath_benchmark.generate = original
    if len(clocks) != 1:
        raise StageClockCaptureError(f"{workload}: expected one clocked generation, got {len(clocks)}.")
    return measurement.digest, clocks[0]


def _children_cpu_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def worker_command(workload: str, domain: PathDomain) -> list[str]:
    """The fresh-interpreter command that generates ``workload`` on ``domain`` and reports its clock."""
    return [sys.executable, "-m", "tools.held_generation_clock", "--domain", domain.value, workload]


def measure_clocked_in_fresh_process(workload: str, domain: PathDomain = PathDomain.POLYGON) -> ClockedGeneration:
    """Generate one workload on ``domain`` in a fresh interpreter (same executable and environment).

    Raises:
        WorkerProcessError: the worker exited non-zero.
        MissingWorkerResultError: the worker reported no clock.
    """
    before = _children_cpu_seconds()
    completed = subprocess.run(
        worker_command(workload, domain),
        cwd=REPOSITORY,
        capture_output=True,
        text=True,
        check=False,
    )
    worker_cpu = _children_cpu_seconds() - before
    if completed.returncode != 0:
        raise WorkerProcessError(f"{workload}: worker exited {completed.returncode}; stderr tail:\n{completed.stderr[-4000:]}")
    lines = [line for line in completed.stdout.splitlines() if line.startswith(WORKER_MARKER)]
    if len(lines) != 1:
        raise MissingWorkerResultError(f"expected one {WORKER_MARKER.strip()} line, found {len(lines)}; stdout tail:\n{completed.stdout[-2000:]}")
    raw = json.loads(lines[0][len(WORKER_MARKER) :])
    return ClockedGeneration(
        workload=workload,
        digest=WorkloadDigest.from_raw(workload, raw["digest"]),
        clock=GenerationClock.from_json(json.dumps(raw["clock"])),
        worker_cpu_seconds=worker_cpu,
        domain=domain,
    )


@dataclass(frozen=True)
class MeasuredRun:
    """One fresh-process generation with the load average read just before it.

    Attributes:
        generation: The clocked generation.
        load_before: ``os.getloadavg()`` (1, 5, 15 minutes) before the run.
        indicative: True when the run was measured without the quiet-machine
            check; such a run never counts toward a budget verdict.
    """

    generation: ClockedGeneration
    load_before: LoadAverage
    indicative: bool


def require_quiet_machine(load_reader: Callable[[], LoadAverage] = os.getloadavg) -> LoadAverage:
    """Return the load average, refusing a machine whose 1-minute load is above the quiet limit.

    Raises:
        MachineNotQuietError: the 1-minute load average is above ``QUIET_LOAD_AVERAGE_MAX``.
    """
    load = load_reader()
    if load[0] > QUIET_LOAD_AVERAGE_MAX:
        raise MachineNotQuietError(
            f"1-minute load average {load[0]:.2f} > {QUIET_LOAD_AVERAGE_MAX} (load {load[0]:.2f} / {load[1]:.2f} / {load[2]:.2f}): "
            "budget verdicts count only from quiet-machine runs (user ruling D2, 2026-10-22). "
            "Re-run on a quiet machine, or measure indicatively: python -m tools.held_generation_clock --indicative <workload>."
        )
    return load


def measure_gated_runs(
    workload: str,
    runs: int,
    load_reader: Callable[[], LoadAverage] = os.getloadavg,
    measure: Callable[[str], ClockedGeneration] = measure_clocked_in_fresh_process,
) -> list[MeasuredRun]:
    """Measure ``runs`` serial fresh-process generations, each only on a quiet machine.

    Raises:
        MachineNotQuietError: the load was above the quiet limit before a run; no further run is measured.
    """
    measured = []
    for _ in range(runs):
        load = require_quiet_machine(load_reader)
        measured.append(MeasuredRun(measure(workload), load, indicative=False))
    return measured


def measure_indicative_runs(
    workload: str,
    runs: int,
    load_reader: Callable[[], LoadAverage] = os.getloadavg,
    measure: Callable[[str], ClockedGeneration] = measure_clocked_in_fresh_process,
) -> list[MeasuredRun]:
    """Measure ``runs`` serial fresh-process generations at any load; every run is marked indicative."""
    measured = []
    for _ in range(runs):
        load = load_reader()
        measured.append(MeasuredRun(measure(workload), load, indicative=True))
    return measured


def describe_indicative(runs: Sequence[MeasuredRun]) -> str:
    """Report indicative runs of one workload: CPU min / median / max, wall, stages, load. No verdict."""
    if not runs:
        raise NoClockedRunsError("an indicative report needs at least one run.")
    workload = runs[0].generation.workload
    domain = runs[0].generation.domain
    cpu = CpuSummary.build([run.generation.clock.cpu_seconds for run in runs])
    walls = [run.generation.clock.wall_seconds for run in runs]
    one_minute = [run.load_before[0] for run in runs]
    digests = {run.generation.digest for run in runs}
    lines = [
        f"{workload} ({domain.value}): INDICATIVE (loaded machine, not a verdict): CPU min {cpu.minimum:.2f} / median {cpu.median:.2f} / max {cpu.maximum:.2f} s, "
        f"spread {cpu.spread:.2f} s; wall median {statistics.median(walls):.2f} s; "
        f"1-min load before runs {min(one_minute):.2f}-{max(one_minute):.2f}; "
        f"{len(digests)} distinct digest(s), circles {runs[0].generation.digest.circle_count} {runs[0].generation.digest.circles[:12]}"
    ]
    for index, run in enumerate(runs):
        clock = run.generation.clock
        stages = " ".join(f"{stage} {seconds:.2f}" for stage, seconds in clock.stage_cpu_seconds.items())
        load = run.load_before
        lines.append(
            f"  run {index + 1}: CPU {clock.cpu_seconds:.2f} s [{stages} exit-collect {clock.exit_cpu_seconds:.2f}]"
            f"  wall {clock.wall_seconds:.2f} s  load {load[0]:.2f} / {load[1]:.2f} / {load[2]:.2f}"
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workload", nargs="+", choices=list(toolpath_benchmark.GENERATORS))
    parser.add_argument(
        "--indicative",
        action="store_true",
        help="measure each workload in fresh serial processes at any load and print the numbers; never a budget verdict",
    )
    parser.add_argument("--runs", type=int, default=INDICATIVE_RUNS_PER_WORKLOAD, help="runs per workload with --indicative")
    parser.add_argument("--domain", choices=[domain.value for domain in PathDomain], default=PathDomain.POLYGON.value, help="centre domain the paths are generated on")
    args = parser.parse_args()
    domain = PathDomain(args.domain)
    if args.indicative:
        for workload in args.workload:
            runs = measure_indicative_runs(workload, args.runs, measure=lambda name: measure_clocked_in_fresh_process(name, domain))
            print(describe_indicative(runs), flush=True)
        return
    if len(args.workload) != 1:
        parser.error("the worker form takes exactly one workload; pass --indicative to measure several.")
    digest, clock = measure_clocked_omega(args.workload[0]) if domain is PathDomain.OMEGA else measure_clocked(args.workload[0])
    print(WORKER_MARKER + json.dumps({"digest": dataclasses.asdict(digest), "clock": json.loads(clock.to_json())}), flush=True)


if __name__ == "__main__":
    main()
