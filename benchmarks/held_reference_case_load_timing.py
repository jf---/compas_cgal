"""Cold and warm Held case-load timing against the reconstruction cache.

Every measurement runs in a fresh interpreter, so process-level memoization
(``functools.cache`` in the certification kernel, the code digest) never
leaks from one sample into the next. Cold: an empty cache directory, so the
load reconstructs, validates and stores. Warm: the directory the cold sample
populated. Serial, median of ``--repeats``.

Run: ``python -m benchmarks.held_reference_case_load_timing`` inside the pixi env.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_reconstruction_cache import ReconstructionCache
from benchmarks.held_reference_reconstruction_identity import reconstruction_code_digest

OUTPUT_DIRECTORY = Path(__file__).resolve().parents[1] / "build" / "held-case-cache-timing"
DEFAULT_REPEATS = 3


class CaseLoadTimingError(RuntimeError):
    """A timing sample subprocess failed or printed no measurement."""


def _single(case: str, cache_directory: Path) -> None:
    cache = ReconstructionCache.build(cache_directory)
    start = time.perf_counter()
    load_held_reference_case(case, reconstruction_cache=cache)
    print(json.dumps({"seconds": time.perf_counter() - start}))


def _sample(case: str, cache_directory: Path) -> float:
    completed = subprocess.run(
        [sys.executable, "-m", "benchmarks.held_reference_case_load_timing", "--single", case, "--cache-dir", str(cache_directory)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise CaseLoadTimingError(f"{case}: sample failed:\n{completed.stderr}")
    seconds = json.loads(completed.stdout.strip().splitlines()[-1])["seconds"]
    if not isinstance(seconds, float):
        raise CaseLoadTimingError(f"{case}: sample printed {completed.stdout!r}.")
    return seconds


def measure(repeats: int) -> dict[str, object]:
    """Time cold and warm loads of every canonical case, serially."""
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    cases: dict[str, dict[str, object]] = {}
    for case in CANONICAL_CASE_NAMES:
        cold: list[float] = []
        warm: list[float] = []
        for _ in range(repeats):
            with tempfile.TemporaryDirectory(dir=OUTPUT_DIRECTORY, prefix=f"{case}-") as directory:
                cold.append(_sample(case, Path(directory)))
                warm.append(_sample(case, Path(directory)))
        cases[case] = {
            "cold_seconds": cold,
            "warm_seconds": warm,
            "cold_median_seconds": statistics.median(cold),
            "warm_median_seconds": statistics.median(warm),
        }
        print(f"{case}: cold {statistics.median(cold):.3f} s, warm {statistics.median(warm):.3f} s (median of {repeats})", flush=True)
    return {
        "repeats": repeats,
        "platform": platform.platform(),
        "python": sys.version,
        "code_sha256": reconstruction_code_digest().hex(),
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--single", choices=CANONICAL_CASE_NAMES)
    parser.add_argument("--cache-dir", type=Path)
    arguments = parser.parse_args()
    if arguments.single is not None:
        if arguments.cache_dir is None:
            parser.error("--single requires --cache-dir")
        _single(arguments.single, arguments.cache_dir)
        return
    result = measure(arguments.repeats)
    path = OUTPUT_DIRECTORY / "timing.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
