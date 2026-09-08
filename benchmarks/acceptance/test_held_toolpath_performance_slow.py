"""Slow-tier regression gate: the upper and Monstera generators against the recorded baseline.

Not collected by pytest; run with `pixi run held-toolpath-gate`. Together these
take about ten minutes on the recording machine.
"""

from __future__ import annotations

from benchmarks import held_toolpath_benchmark as tb


def test_slow_generators_reproduce_the_baseline_geometry_no_slower() -> None:
    assert tb.BASELINE_PATH.is_file(), f"record it with `pixi run held-toolpath-baseline`: {tb.BASELINE_PATH}"
    baseline = tb.read_json(tb.BASELINE_PATH)
    current = tb.measure(tb.SLOW_WORKLOADS)
    subset = tb.ToolpathBenchmark(baseline.schema_version, baseline.identity, baseline.build, {name: baseline.runs[name] for name in tb.SLOW_WORKLOADS})
    findings = tb.compare(current, subset)
    assert findings == [], "\n".join(tb.summary_lines(current) + [finding.message for finding in findings])
