"""Time the Fraction and native station-boundary analyses over every raw-guide station.

Serial, ``time.perf_counter``, median of three runs per path. Both paths analyse
the same prepared guide; the guide build and case load are outside the clock.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from collections.abc import Callable
from pathlib import Path

from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import _analyze_station_boundary
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_figure5_raw_guide_native import native_station_boundary_analyses
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case

RUNS_PER_PATH = 3
DEFAULT_CASES = ("figure5", "figure8_crossed_skis")


def _fraction_path(guide: Figure5RawGuide) -> int:
    return sum(len(_analyze_station_boundary(guide, station).distance_admissible_sites) for run in guide.runs for station in run.stations)


def _native_path(guide: Figure5RawGuide) -> int:
    return sum(len(analysis.distance_admissible_sites) for run in native_station_boundary_analyses(guide) for analysis in run)


def _median_seconds(path: Callable[[Figure5RawGuide], int], guide: Figure5RawGuide) -> tuple[float, list[float], int]:
    samples: list[float] = []
    sites = -1
    for _ in range(RUNS_PER_PATH):
        started = time.perf_counter()
        sites = path(guide)
        samples.append(time.perf_counter() - started)
    return statistics.median(samples), samples, sites


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", choices=CANONICAL_CASE_NAMES)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    results = {}
    for name in arguments.case or DEFAULT_CASES:
        guide = build_held_reference_raw_guide(load_held_reference_case(name))
        native_median, native_samples, native_sites = _median_seconds(_native_path, guide)
        fraction_median, fraction_samples, fraction_sites = _median_seconds(_fraction_path, guide)
        if native_sites != fraction_sites:
            raise RuntimeError(f"{name}: native admitted {native_sites} sites, Fraction path {fraction_sites}.")
        results[name] = {
            "stations": guide.station_count,
            "boundary_segments": guide.boundary_segment_count,
            "admissible_sites": fraction_sites,
            "fraction_median_s": fraction_median,
            "fraction_samples_s": fraction_samples,
            "native_median_s": native_median,
            "native_samples_s": native_samples,
            "speedup": fraction_median / native_median,
        }
        print(name, json.dumps(results[name]), flush=True)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
