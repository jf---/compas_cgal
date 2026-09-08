"""Pure contracts of the curved-circle performance benchmark: aggregation, JSON round trip, regression policy."""

from __future__ import annotations

import pytest

from benchmarks import held_curved_circle_benchmark as bench


def _identity(cpu: str = "Apple M1 Max") -> bench.MachineIdentity:
    return bench.MachineIdentity(machine="arm64", system="Darwin", cpu_brand=cpu, python="3.12.13")


def _case(name: str, seconds: list[float], stop: int | None = None, stop_type: str | None = None) -> bench.CaseBenchmark:
    pieces = tuple(bench.PieceTiming(piece_index=i, kind="arc", seconds=s, samples=(s, s * 1.1, s * 0.9)) for i, s in enumerate(seconds))
    return bench.CaseBenchmark(
        case=name,
        native_pieces=len(seconds) if stop is None else stop + 1 + 5,
        completed_queries=len(seconds),
        stop_piece_index=stop,
        stop_type=stop_type,
        import_seconds=0.05,
        pieces=pieces,
    )


def _benchmark(cases: dict[str, bench.CaseBenchmark], cpu: str = "Apple M1 Max") -> bench.CurvedCircleBenchmark:
    return bench.CurvedCircleBenchmark(
        schema_version=bench.SCHEMA_VERSION,
        identity=_identity(cpu),
        build=bench.BuildIdentity(commit="abc123", dirty=False, package_version="0.9.1", recorded_at="2026-09-08T00:00:00+00:00"),
        repeats=3,
        parameter=0.5,
        cases=cases,
    )


def test_case_aggregates_median_p95_and_max_over_piece_medians() -> None:
    case = _case("figure5", [0.001, 0.002, 0.003, 0.004, 0.010])
    assert case.median_seconds == pytest.approx(0.003)
    assert case.max_seconds == pytest.approx(0.010)
    assert case.p95_seconds == pytest.approx(0.010)


def test_json_round_trip_preserves_every_field() -> None:
    original = _benchmark({"figure5": _case("figure5", [0.001, 0.002]), "figure8_upper": _case("figure8_upper", [0.003], stop=1, stop_type="NoPositiveNativeBoundaryCircleError")})
    assert bench.CurvedCircleBenchmark.from_json(original.to_json()) == original


def test_unknown_schema_version_is_rejected() -> None:
    text = _benchmark({"figure5": _case("figure5", [0.001])}).to_json().replace(f'"schema_version": {bench.SCHEMA_VERSION}', '"schema_version": 999')
    with pytest.raises(bench.BaselineSchemaError):
        bench.CurvedCircleBenchmark.from_json(text)


def test_same_machine_within_factors_has_no_findings() -> None:
    baseline = _benchmark({"figure5": _case("figure5", [0.001, 0.002, 0.003])})
    current = _benchmark({"figure5": _case("figure5", [0.002, 0.004, 0.006])})  # 2x median, 2x max
    assert bench.compare(current, baseline) == []


def test_same_machine_median_beyond_factor_is_a_finding() -> None:
    baseline = _benchmark({"figure5": _case("figure5", [0.001, 0.001, 0.001])})
    current = _benchmark({"figure5": _case("figure5", [0.0035, 0.0035, 0.0035])})
    findings = bench.compare(current, baseline)
    assert [f.kind for f in findings] == ["median"]
    assert findings[0].bound == pytest.approx(0.001 * bench.SAME_MACHINE_MEDIAN_FACTOR)


def test_cross_machine_uses_the_looser_factor() -> None:
    baseline = _benchmark({"figure5": _case("figure5", [0.001, 0.001, 0.001])}, cpu="Apple M1 Max")
    eight_times = _benchmark({"figure5": _case("figure5", [0.008, 0.008, 0.008])}, cpu="AMD EPYC 7763")
    twelve_times = _benchmark({"figure5": _case("figure5", [0.012, 0.012, 0.012])}, cpu="AMD EPYC 7763")
    assert bench.compare(eight_times, baseline) == []
    assert [f.kind for f in bench.compare(twelve_times, baseline)] == ["median", "max"]


def test_absolute_ceilings_apply_on_every_machine() -> None:
    baseline = _benchmark({"figure5": _case("figure5", [0.1, 0.1, 0.1])})
    current = _benchmark({"figure5": _case("figure5", [0.1, 0.1, 0.3])}, cpu="AMD EPYC 7763")
    kinds = [f.kind for f in bench.compare(current, baseline)]
    assert "absolute_piece_ceiling" in kinds
    assert "absolute_median_ceiling" in kinds


def test_changed_stop_or_fewer_completed_pieces_is_a_finding() -> None:
    baseline = _benchmark({"figure8_upper": _case("figure8_upper", [0.001, 0.001], stop=2, stop_type="NoPositiveNativeBoundaryCircleError")})
    earlier_stop = _benchmark({"figure8_upper": _case("figure8_upper", [0.001], stop=1, stop_type="NoPositiveNativeBoundaryCircleError")})
    no_stop = _benchmark({"figure8_upper": _case("figure8_upper", [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])})
    assert {f.kind for f in bench.compare(earlier_stop, baseline)} == {"completed", "stop"}
    assert {f.kind for f in bench.compare(no_stop, baseline)} == {"completed", "stop"}


def test_missing_case_is_a_finding() -> None:
    baseline = _benchmark({"figure5": _case("figure5", [0.001])})
    assert [f.kind for f in bench.compare(_benchmark({}), baseline)] == ["missing_case"]
