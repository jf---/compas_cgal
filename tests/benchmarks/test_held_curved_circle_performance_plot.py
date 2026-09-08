"""The performance figure renders from a benchmark record, with a baseline overlay and stop markers."""

from __future__ import annotations

import pathlib

from benchmarks import held_curved_circle_benchmark as bench
from benchmarks import held_curved_circle_performance_plot as plot


def _record(seconds: dict[str, list[float]], cpu: str, stops: dict[str, tuple[int, str]] | None = None) -> bench.CurvedCircleBenchmark:
    stops = stops or {}
    cases = {}
    for name, values in seconds.items():
        stop = stops.get(name)
        pieces = tuple(bench.PieceTiming(i, "arc", s, (s,)) for i, s in enumerate(values))
        cases[name] = bench.CaseBenchmark(
            name, len(values) + (5 if stop else 0), len(values), stop[0] if stop else None, stop[1] if stop else None, 0.01, pieces
        )
    identity = bench.MachineIdentity("arm64", "Darwin", cpu, "3.12")
    build = bench.BuildIdentity("abc", False, "0.9.1", "2026-09-08T00:00:00+00:00")
    return bench.CurvedCircleBenchmark(bench.SCHEMA_VERSION, identity, build, 3, 0.5, cases)


def test_figure_renders_every_case_with_baseline_and_stop(tmp_path: pathlib.Path) -> None:
    current = _record({"figure5": [0.001, 0.002], "figure8_upper": [0.003, 0.004]}, "Apple M1 Max", {"figure8_upper": (2, "NoPositiveNativeBoundaryCircleError")})
    baseline = _record({"figure5": [0.001, 0.001], "figure8_upper": [0.002, 0.002]}, "Apple M1 Max", {"figure8_upper": (2, "NoPositiveNativeBoundaryCircleError")})
    output = tmp_path / "perf.png"
    summary = plot.render_performance_figure(current, baseline, output)
    assert output.is_file() and output.stat().st_size > 0
    assert summary.cases == ("figure5", "figure8_upper")
    assert summary.stops == {"figure8_upper": 2}
