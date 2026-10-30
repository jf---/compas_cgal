"""Measure every gate generator's generation time serially, one fresh process per run.

The gate measurement (`tools.gate_measure`) runs its cases in a parallel pool,
so its ``generation_s`` is wall clock under contention. This tool is the
generation-time authority: it times every gate case (`tools.gate_measure.gate_cases`,
the same cap x generator x pocket construction) one case at a time, each of
`RUNS_PER_CASE` runs in a fresh interpreter (``python -m tools.pipeline_timing
--case KEY``), and records the median, the operation count and the ops digest.
Runs of one case that disagree on the digest or the count fail loud: the
generators are deterministic, so disagreement is a defect, never noise.

It records the 1/5/15-minute load average before the first case and after the
last, and the run's `tools.run_provenance.RunProvenance` (schema
``pipeline-timing/2``). It does not refuse a loaded machine; the
consumers that need quiet seconds (the generation-time gate,
`tools.pipeline_time_gate`) read the recorded load and refuse a verdict.

Run as ``pixi run pipeline-timing``; it writes
`benchmarks/results/pipeline_timing.json` and its rendering
`benchmarks/results/pipeline_timing.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Sequence

from benchmarks.gate import GATE_GENERATORS
from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket
from compas_cgal.helix import Helix
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from tools.gate_measure import case_key
from tools.gate_measure import gate_cases
from tools.gate_measure import partition_cases
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import load_gate_vectors
from tools.gate_vectors import replay_admissions
from tools.held_generation_clock import LoadAverage
from tools.pipeline_timing_results import DIGEST_HEX_CHARS
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import PIPELINE_TIMING_TABLE_PATH
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import STORED_SECONDS_DECIMALS
from tools.pipeline_timing_results import CaseTiming
from tools.pipeline_timing_results import PipelineTimings
from tools.pipeline_timing_results import stored_load
from tools.run_native_identity import NativeIdentity
from tools.run_native_identity import NativeTracker
from tools.run_native_identity import track_native_here
from tools.run_provenance import RunProvenance
from tools.run_provenance import TreeIdentity
from tools.run_provenance import TreeReader
from tools.run_provenance import file_set_sha256
from tools.run_provenance import interpreter_identity_here
from tools.run_provenance import machine_identity
from tools.run_provenance import require_unchanged_tree
from tools.run_provenance import sha256_file
from tools.run_provenance import tree_identity
from tools.run_provenance import utc_now

ROOT = Path(__file__).resolve().parents[1]
# The child prints its result on one stdout line after this prefix, so library
# output on stdout cannot be mistaken for it.
RESULT_LINE_PREFIX = "PIPELINE_TIMING_RESULT "


class UnknownTimingCaseError(Exception):
    """A case key that is not a gate case."""


class CaseProcessFailedError(Exception):
    """A run's child process exited non-zero."""


class MissingCaseResultError(Exception):
    """A run's child process printed no result line."""


class NondeterministicCaseError(Exception):
    """The runs of one case produced different operations."""


@dataclass(frozen=True)
class CaseRun:
    """One generation of one case, as its child process reports it.

    Attributes:
        generation_s: Wall-clock seconds of the generator call.
        operations: The number of operations generated.
        digest: The ops digest (`ops_digest`).
    """

    generation_s: float
    operations: int
    digest: str


Launcher = Callable[[str], CaseRun]
LoadReader = Callable[[], LoadAverage]
# Builds a run's provenance from its UTC start and end.
Describer = Callable[[TreeIdentity, NativeIdentity, str, str], RunProvenance]
# The sources of this timing tool, digested into the provenance.
TIMING_TOOL_SOURCES = (
    Path("tools") / "pipeline_timing.py",
    Path("tools") / "pipeline_timing_results.py",
    Path("tools") / "run_provenance.py",
    Path("tools") / "run_native_identity.py",
)
PIXI_LOCK = Path("pixi.lock")


def _helix_data(helix: Helix) -> dict[str, Any]:
    """Every defining field of *helix*, as plain values in a fixed key order.

    A `Helix` is no compas geometry and has no ``__data__``; its digest content
    is its exact definition: the footprint loop (centre, radius, exact entry,
    realised phase and chart parameter), both heights, the turn count and the
    turn sense. Every float is a double whose ``repr`` round-trips exactly.
    """
    footprint = helix.footprint
    return {
        "footprint": {
            "centre": (float(footprint.centre.x), float(footprint.centre.y)),
            "radius": float(footprint.radius.value),
            "entry": (float(footprint.entry.x), float(footprint.entry.y)),
            "phase": float(footprint.phase),
            "chart": footprint.chart,
            "chart_numerator": footprint.chart_numerator,
        },
        "top_z": float(helix.top_z.value),
        "bottom_z": float(helix.bottom_z.value),
        "turns": helix.turns,
        "clockwise": helix.clockwise,
    }


def _geometry_data(operation: ToolpathOperation) -> Any:
    """The digest content of *operation*'s geometry: compas ``__data__``, or `_helix_data` for a helix."""
    geometry = operation.geometry
    if isinstance(geometry, Helix):
        return _helix_data(geometry)
    return geometry.__data__


def ops_digest(result: ToolpathResult) -> str:
    """The truncated SHA-256 over every operation's kind, geometry type and geometry data, in order.

    Line, arc and circle data is compas ``__data__``; a helix's is its
    defining fields (`_helix_data`). Output without a helix digests exactly as
    before the helix existed.

    Args:
        result: A generated path.

    Returns:
        `DIGEST_HEX_CHARS` lowercase hex characters.
    """
    digest = hashlib.sha256()
    for operation in result.operations:
        digest.update(repr((operation.operation, type(operation.geometry).__name__, _geometry_data(operation))).encode())
    return digest.hexdigest()[:DIGEST_HEX_CHARS]


def run_case_here(key: str) -> CaseRun:
    """Generate *key*'s case once in this process and time the generator call alone.

    Args:
        key: A gate case key.

    Returns:
        The run.

    Raises:
        UnknownTimingCaseError: *key* is not a gate case.
    """
    cases = {case_key(cap_id, generator, pocket): (cap, generator, pocket) for cap_id, cap, generator, pocket in gate_cases()}
    if key not in cases:
        raise UnknownTimingCaseError(f"{key!r} is not a gate case.")
    cap, generator, pocket = cases[key]
    spec = gate_pocket(pocket, tea_cap_deg=GateCapDegrees(cap))
    started = time.perf_counter()
    result = GATE_GENERATORS[generator](spec)
    elapsed = time.perf_counter() - started
    return CaseRun(generation_s=elapsed, operations=len(result.operations), digest=ops_digest(result))


def launch_fresh_process(key: str) -> CaseRun:
    """Run *key*'s case once in a fresh interpreter and read its result line.

    The child inherits this process's environment (the editable build the pixi
    task selects) and runs from the repository root.

    Args:
        key: A gate case key.

    Returns:
        The child's run.

    Raises:
        CaseProcessFailedError: The child exited non-zero; its stderr is quoted.
        MissingCaseResultError: The child printed no result line.
    """
    command = [sys.executable, "-m", "tools.pipeline_timing", "--case", key]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise CaseProcessFailedError(f"{key}: {' '.join(command)} exited {completed.returncode}:\n{completed.stderr}")
    lines = [line for line in completed.stdout.splitlines() if line.startswith(RESULT_LINE_PREFIX)]
    if len(lines) != 1:
        raise MissingCaseResultError(f"{key}: the child printed {len(lines)} result lines, not one.")
    payload = json.loads(lines[0][len(RESULT_LINE_PREFIX) :])
    return CaseRun(generation_s=float(payload["generation_s"]), operations=int(payload["operations"]), digest=str(payload["digest"]))


def time_case(key: str, launch: Launcher) -> CaseTiming:
    """Run *key* `RUNS_PER_CASE` times through *launch*, consecutively, and require agreement.

    Args:
        key: A gate case key.
        launch: Runs the case once, in a fresh process for a real measurement.

    Returns:
        The case's timing.

    Raises:
        NondeterministicCaseError: Two runs disagree on the digest or the
            operation count.
    """
    runs = [launch(key) for _run in range(RUNS_PER_CASE)]
    first = runs[0]
    if any((run.digest, run.operations) != (first.digest, first.operations) for run in runs):
        produced = [(run.digest, run.operations) for run in runs]
        raise NondeterministicCaseError(f"{key}: the runs produced different operations (digest, count): {produced}.")
    return CaseTiming(runs_s=tuple(round(run.generation_s, STORED_SECONDS_DECIMALS) for run in runs), operations=first.operations, digest=first.digest)


def read_this_tree() -> TreeIdentity:
    """This tree's identity (`ROOT`)."""
    return tree_identity(ROOT)


def describe_this_run(tree: TreeIdentity, native: NativeIdentity, started_utc: str, ended_utc: str) -> RunProvenance:
    """The provenance of a timing run of this tree (`ROOT`) in this interpreter, between *started_utc* and *ended_utc*."""
    rev, dirty, patch = tree
    environment, python_version, stock_2 = interpreter_identity_here()
    host, cpu_model = machine_identity()
    return RunProvenance(
        rev=rev,
        dirty=dirty,
        patch_sha256=patch,
        timing_tool_sha256=file_set_sha256(ROOT, TIMING_TOOL_SOURCES),
        pixi_lock_sha256=sha256_file(ROOT / PIXI_LOCK),
        environment=environment,
        python_version=python_version,
        stock_2_sha256=stock_2,
        native_sha256=native.native_sha256,
        native_sources_sha256=native.native_sources_sha256,
        host=host,
        cpu_model=cpu_model,
        started_utc=started_utc,
        ended_utc=ended_utc,
    )


def measure_pipeline_timings(
    keys: Sequence[str],
    *,
    excluded: Sequence[str],
    launch: Launcher = launch_fresh_process,
    load_reader: LoadReader = os.getloadavg,
    describe: Describer = describe_this_run,
    read_tree: TreeReader = read_this_tree,
    track_native: NativeTracker = track_native_here(ROOT),
) -> PipelineTimings:
    """Time *keys* serially, one case at a time, in order, reading the load before the first case and after the last.

    Args:
        keys: The case keys to time.
        excluded: The case keys the admissions exclude, recorded untimed.
        launch: Runs one case once.
        load_reader: Reads the 1/5/15-minute load average.
        describe: Builds the run's provenance from its tree identity and UTC start and end.
        read_tree: Reads the measured tree's identity; read at the start and the end, which must agree.
        track_native: Starts the native bracket after the start read; its end, after the end read, yields the run's native identity.

    Returns:
        The timings (unvalidated; `PipelineTimings.build` checks coverage).

    Raises:
        TreeChangedDuringRunError: The tree's identity changed during the run.
        NativeIdentityError: The native build changed during the run, is stale, or was not compiled from this tree
            (``tools.held_native_identity``; ``NativeSourcesMismatchError`` from ``tools.run_native_identity``).
    """
    tree_at_start, started_utc, load_start = read_tree(), utc_now(), load_reader()
    finish_native = track_native()
    cases = {key: time_case(key, launch) for key in keys}
    load_end, ended_utc, tree_at_end = load_reader(), utc_now(), read_tree()
    require_unchanged_tree(tree_at_start, tree_at_end, "the 1.2 timing tree")
    native = finish_native()
    return PipelineTimings(
        cases=cases, excluded=tuple(excluded), load_start=stored_load(load_start), load_end=stored_load(load_end), provenance=describe(tree_at_end, native, started_utc, ended_utc)
    )


def gate_timing_keys() -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Every gate case key to time and the keys the committed admissions exclude, in gate-case order."""
    admissions = load_gate_vectors(ROOT / GATE_VECTORS_PATH).admissions
    to_time, excluded = partition_cases(gate_cases(), replay_admissions(admissions)[1])
    return (
        tuple(case_key(cap_id, generator, pocket) for cap_id, _cap, generator, pocket in to_time),
        tuple(case_key(cap_id, generator, pocket) for (cap_id, _cap, generator, pocket), _error in excluded),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry: time every case and write the document and its table, or run one case as a child.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    from tools.pipeline_timing_table import render_pipeline_timing_table

    parser = argparse.ArgumentParser(prog="python -m tools.pipeline_timing", description="Time every gate case serially, one fresh process per run.")
    parser.add_argument("--case", default=None, metavar="KEY", help="child mode: run one case once in this process and print its result line")
    parser.add_argument("--write", type=Path, default=PIPELINE_TIMING_PATH)
    parser.add_argument("--table", type=Path, default=PIPELINE_TIMING_TABLE_PATH)
    args = parser.parse_args(argv)
    if args.case is not None:
        run = run_case_here(args.case)
        print(RESULT_LINE_PREFIX + json.dumps({"generation_s": run.generation_s, "operations": run.operations, "digest": run.digest}))
        return 0
    keys, excluded = gate_timing_keys()
    timings = PipelineTimings.build(measure_pipeline_timings(keys, excluded=excluded).to_json())
    args.write.write_text(timings.dumps())
    args.table.write_text(render_pipeline_timing_table(timings))
    print(f"wrote {args.write} and {args.table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
