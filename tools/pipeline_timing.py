"""Measure every gate generator's generation time serially, one fresh process per run.

The gate measurement (`tools.gate_measure`) runs its cases in a parallel pool,
so its ``generation_s`` is wall clock under contention. This tool is the
generation-time authority: it times every gate case (`tools.gate_measure.gate_cases`,
the same cap x generator x pocket construction) one case at a time, each of
`RUNS_PER_CASE` runs in a fresh interpreter (``python -m tools.pipeline_timing
--case KEY``), and records the median, the operation count and the ops digest.
Runs of one case that disagree on the digest or the count fail loud: the
generators are deterministic, so disagreement is a defect, never noise.

Run as ``pixi run pipeline-timing``; it writes
`benchmarks/results/pipeline_timing.json` and its rendering
`benchmarks/results/pipeline_timing.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Sequence

from benchmarks.gate import GATE_GENERATORS
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
from tools.pipeline_timing_results import DIGEST_HEX_CHARS
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import PIPELINE_TIMING_TABLE_PATH
from tools.pipeline_timing_results import RUNS_PER_CASE
from tools.pipeline_timing_results import STORED_SECONDS_DECIMALS
from tools.pipeline_timing_results import CaseTiming
from tools.pipeline_timing_results import PipelineTimings

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
    spec = gate_pocket(pocket, tea_cap_deg=cap)  # type: ignore[arg-type]
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


def measure_pipeline_timings(keys: Sequence[str], *, excluded: Sequence[str], launch: Launcher = launch_fresh_process) -> PipelineTimings:
    """Time *keys* serially, one case at a time, in order.

    Args:
        keys: The case keys to time.
        excluded: The case keys the admissions exclude, recorded untimed.
        launch: Runs one case once.

    Returns:
        The timings (unvalidated; `PipelineTimings.build` checks coverage).
    """
    return PipelineTimings(cases={key: time_case(key, launch) for key in keys}, excluded=tuple(excluded))


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
