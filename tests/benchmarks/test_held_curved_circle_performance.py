"""Regression gate: the curved-circle query on the prepared Held pockets against the recorded baseline."""

from __future__ import annotations

from benchmarks import held_curved_circle_benchmark as bench


def test_baseline_is_recorded_with_identity() -> None:
    assert bench.BASELINE_PATH.is_file(), f"record it with `pixi run held-curved-circle-baseline`: {bench.BASELINE_PATH}"
    baseline = bench.read_json(bench.BASELINE_PATH)
    assert baseline.schema_version == bench.SCHEMA_VERSION
    assert baseline.build.commit and baseline.identity.cpu_brand
    assert set(baseline.cases) == set(bench.CANONICAL_CASE_NAMES)


def test_current_build_is_not_slower_than_the_baseline_and_stops_where_it_stopped() -> None:
    assert bench.BASELINE_PATH.is_file(), f"record it with `pixi run held-curved-circle-baseline`: {bench.BASELINE_PATH}"
    baseline = bench.read_json(bench.BASELINE_PATH)
    current = bench.measure_all(repeats=bench.DEFAULT_REPEATS)
    findings = bench.compare(current, baseline)
    report = "\n".join(bench.summary_lines(current) + [finding.message for finding in findings])
    assert findings == [], report
