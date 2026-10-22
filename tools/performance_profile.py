"""Profile the performance workloads with cProfile, one fresh process per workload.

Each workload (`tools.performance_workloads.PERFORMANCE_WORKLOADS`) runs once
under cProfile in its own interpreter (``python -m tools.performance_profile
--workload KEY``), so no workload inherits another's caches or allocator
state. Only the workload call is profiled: imports and interpreter start-up
are outside it. The child writes ``build/performance/<key>.prof`` and a sidecar
``<key>.json`` with the commit, the date and the profiled wall time.

Run as ``pixi run performance-profile``.
"""

from __future__ import annotations

import argparse
import cProfile
import json
import subprocess
import sys
import time
from datetime import date
from typing import Callable
from typing import Sequence

from tools.performance_workloads import PERFORMANCE_WORKLOADS
from tools.performance_workloads import PROFILE_DIR
from tools.performance_workloads import ROOT
from tools.performance_workloads import WORKLOADS_BY_KEY
from tools.performance_workloads import Lane
from tools.performance_workloads import PerformanceWorkload


class UnknownPerformanceWorkloadError(Exception):
    """A workload key that is not in the registry."""


class ProfileProcessFailedError(Exception):
    """A workload's child process exited non-zero."""


class GitStateUnavailableError(Exception):
    """``git`` could not report the commit the profile was taken at."""


def git_state() -> tuple[str, bool]:
    """The current commit and whether the working tree has uncommitted changes.

    Returns:
        ``(commit, dirty)``.

    Raises:
        GitStateUnavailableError: ``git`` failed.
    """
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False)
    status = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, check=False)
    if head.returncode != 0 or status.returncode != 0:
        raise GitStateUnavailableError(f"git failed: {head.stderr}{status.stderr}")
    return head.stdout.strip(), bool(status.stdout.strip())


def _workload_call(workload: PerformanceWorkload) -> Callable[[], object]:
    if workload.lane is Lane.HELD:
        from benchmarks.held_toolpath_benchmark import WORKLOADS

        held = WORKLOADS[workload.target]
        return lambda: held()
    from tools.pipeline_timing import run_case_here

    return lambda: run_case_here(workload.target)


def profile_here(workload: PerformanceWorkload) -> None:
    """Run *workload* once under cProfile in this process and write its profile and sidecar.

    Args:
        workload: The workload.
    """
    call = _workload_call(workload)
    commit, dirty = git_state()
    profiler = cProfile.Profile()
    started = time.perf_counter()
    profiler.enable()
    call()
    profiler.disable()
    elapsed = time.perf_counter() - started
    PROFILE_DIR.mkdir(parents=True, exist_ok=True)
    profiler.dump_stats(str(workload.profile_path))
    meta = {"workload": workload.key, "target": workload.target, "commit": commit, "dirty": dirty, "date": date.today().isoformat(), "profiled_wall_s": elapsed}
    workload.meta_path.write_text(json.dumps(meta, indent=2) + "\n")


def launch_fresh_process(workload: PerformanceWorkload) -> None:
    """Profile *workload* in a fresh interpreter that inherits this environment.

    Args:
        workload: The workload.

    Raises:
        ProfileProcessFailedError: The child exited non-zero; its stderr is quoted.
    """
    command = [sys.executable, "-m", "tools.performance_profile", "--workload", workload.key]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise ProfileProcessFailedError(f"{workload.key}: {' '.join(command)} exited {completed.returncode}:\n{completed.stderr}")


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry: profile every workload serially, or one workload as a child.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.

    Raises:
        UnknownPerformanceWorkloadError: ``--workload`` names no registered workload.
    """
    parser = argparse.ArgumentParser(prog="python -m tools.performance_profile", description="Profile every performance workload, one fresh process each.")
    parser.add_argument("--workload", default=None, metavar="KEY", help="child mode: profile one workload in this process")
    args = parser.parse_args(argv)
    if args.workload is not None:
        if args.workload not in WORKLOADS_BY_KEY:
            raise UnknownPerformanceWorkloadError(f"{args.workload!r}; known: {sorted(WORKLOADS_BY_KEY)}")
        profile_here(WORKLOADS_BY_KEY[args.workload])
        return 0
    for workload in PERFORMANCE_WORKLOADS:
        started = time.perf_counter()
        launch_fresh_process(workload)
        print(f"{workload.key}: wrote {workload.profile_path.relative_to(ROOT)} in {time.perf_counter() - started:.1f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
