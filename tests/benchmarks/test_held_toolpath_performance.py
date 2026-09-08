"""Regression gate: the Held toolpath generators against the recorded baseline, geometry and time."""

from __future__ import annotations

from benchmarks import held_toolpath_benchmark as tb


def test_baseline_is_recorded_for_every_workload() -> None:
    assert tb.BASELINE_PATH.is_file(), f"record it with `pixi run held-toolpath-baseline`: {tb.BASELINE_PATH}"
    baseline = tb.read_json(tb.BASELINE_PATH)
    assert set(baseline.runs) == set(tb.WORKLOADS)
    assert baseline.build.commit and baseline.identity.cpu_brand


def test_fast_generators_reproduce_the_baseline_geometry_no_slower() -> None:
    # Figure 5 refined and crossed skis, about a minute together; the slow
    # workloads gate in benchmarks/acceptance via `pixi run held-toolpath-gate`.
    assert tb.BASELINE_PATH.is_file(), f"record it with `pixi run held-toolpath-baseline`: {tb.BASELINE_PATH}"
    baseline = tb.read_json(tb.BASELINE_PATH)
    current = tb.measure(tb.FAST_WORKLOADS)
    findings = tb.compare(current, tb.ToolpathBenchmark(baseline.schema_version, baseline.identity, baseline.build, {name: baseline.runs[name] for name in tb.FAST_WORKLOADS}))
    assert findings == [], "\n".join(tb.summary_lines(current) + [finding.message for finding in findings])
