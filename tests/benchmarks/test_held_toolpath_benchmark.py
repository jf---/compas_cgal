"""Record and regression policy of the Held toolpath generator benchmark."""

from __future__ import annotations

import pytest

from benchmarks import held_toolpath_benchmark as tb
from benchmarks import performance_record as record


def _identity(cpu: str = "Apple M1 Max") -> record.MachineIdentity:
    return record.MachineIdentity(machine="arm64", system="Darwin", cpu_brand=cpu, python="3.12.13")


def _run(name: str, seconds: float, circles: int = 100, length: float = 1234.5, engagement: float = 79.99, rejections: int = 0) -> tb.WorkloadRun:
    # Held paths close: one connector per circle, the last runs back to the first.
    return tb.WorkloadRun(
        workload=name,
        metrics=tb.ToolpathMetrics(circles=circles, connectors=circles, path_length_mm=length, max_engagement_deg=engagement, containment_rejections=rejections),
        stage_seconds={"guide": seconds * 0.1, "offset": seconds * 0.1, "initial": seconds * 0.6, "refine": seconds * 0.1, "connectors": seconds * 0.1},
        generation_seconds=seconds,
    )


def _benchmark(runs: dict[str, tb.WorkloadRun], cpu: str = "Apple M1 Max") -> tb.ToolpathBenchmark:
    build = record.BuildIdentity(commit="abc123", dirty=False, package_version="0.9.1", recorded_at="2026-09-08T00:00:00+00:00")
    return tb.ToolpathBenchmark(schema_version=tb.SCHEMA_VERSION, identity=_identity(cpu), build=build, runs=runs)


CEILINGS = {"figure5_refined": 60.0}


def test_json_round_trip_preserves_every_field() -> None:
    original = _benchmark({"figure5_refined": _run("figure5_refined", 20.0)})
    assert tb.ToolpathBenchmark.from_json(original.to_json()) == original


def test_identical_geometry_and_faster_time_has_no_findings() -> None:
    baseline = _benchmark({"figure5_refined": _run("figure5_refined", 20.0)})
    current = _benchmark({"figure5_refined": _run("figure5_refined", 15.0)})
    assert tb.compare(current, baseline, ceilings=CEILINGS) == []


def test_any_change_in_circles_connectors_length_engagement_or_rejections_is_a_finding() -> None:
    baseline = _benchmark({"figure5_refined": _run("figure5_refined", 20.0)})
    changed = _benchmark({"figure5_refined": _run("figure5_refined", 20.0, circles=101, length=1234.6, engagement=80.1, rejections=1)})
    kinds = [f.kind for f in tb.compare(changed, baseline, ceilings=CEILINGS)]
    assert kinds == ["circles", "connectors", "path_length", "max_engagement", "containment_rejections"]


def test_same_machine_time_factor_and_absolute_ceiling() -> None:
    baseline = _benchmark({"figure5_refined": _run("figure5_refined", 20.0)})
    slower = _benchmark({"figure5_refined": _run("figure5_refined", 20.0 * record.SAME_MACHINE_FACTOR + 1.0)})
    kinds = [f.kind for f in tb.compare(slower, baseline, ceilings=CEILINGS)]
    assert kinds == ["generation_time", "absolute_ceiling"]


def test_cross_machine_uses_the_looser_factor() -> None:
    baseline = _benchmark({"figure5_refined": _run("figure5_refined", 2.0)}, cpu="Apple M1 Max")
    eight = _benchmark({"figure5_refined": _run("figure5_refined", 16.0)}, cpu="AMD EPYC 7763")
    twelve = _benchmark({"figure5_refined": _run("figure5_refined", 24.0)}, cpu="AMD EPYC 7763")
    assert tb.compare(eight, baseline, ceilings=CEILINGS) == []
    assert [f.kind for f in tb.compare(twelve, baseline, ceilings=CEILINGS)] == ["generation_time"]


def test_missing_workload_is_a_finding() -> None:
    baseline = _benchmark({"figure5_refined": _run("figure5_refined", 20.0)})
    assert [f.kind for f in tb.compare(_benchmark({}), baseline, ceilings=CEILINGS)] == ["missing_workload"]


def test_every_gated_workload_has_an_absolute_ceiling() -> None:
    assert set(tb.WORKLOADS) == set(tb.ABSOLUTE_CEILINGS_S)


def test_path_length_ignores_reporting_round_off_only() -> None:
    baseline = _benchmark({"figure5_refined": _run("figure5_refined", 20.0, length=1234.5)})
    same = _benchmark({"figure5_refined": _run("figure5_refined", 20.0, length=1234.5 + tb.REPORTED_LENGTH_TOLERANCE_MM / 2)})
    assert tb.compare(same, baseline, ceilings=CEILINGS) == []
    with pytest.raises(AssertionError):
        assert tb.compare(_benchmark({"figure5_refined": _run("figure5_refined", 20.0, length=1235.5)}), baseline, ceilings=CEILINGS) == []


def test_seam_engagement_round_trips_and_is_reported_not_gated() -> None:
    run = _run("figure5_refined", 20.0)
    seamed = tb.WorkloadRun(run.workload, tb.ToolpathMetrics(**{**run.metrics.__dict__, "seam_engagement_deg": 123.0}), run.stage_seconds, run.generation_seconds)
    benchmark = _benchmark({"figure5_refined": seamed})
    assert tb.ToolpathBenchmark.from_json(benchmark.to_json()) == benchmark
    assert tb.compare(benchmark, _benchmark({"figure5_refined": run}), ceilings=CEILINGS) == []
    assert "seam 123.0000 deg" in tb.summary_lines(benchmark)[1]
