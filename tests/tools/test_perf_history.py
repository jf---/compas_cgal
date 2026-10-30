"""The performance history: adapters per record family, the append-only store, the quiet baseline and the regression check."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from pathlib import PurePosixPath
from typing import Optional

import pytest
import tools.perf_history_adapters as perf_history_adapters
import tools.perf_history_ingest as perf_history_ingest
from tests.tools.test_coverage_certificate_omega_star_budget import LOADED as BUDGET_LOAD_LOADED
from tests.tools.test_coverage_certificate_omega_star_budget import _quiet_budget as quiet_omega_star_budget
from tests.tools.test_coverage_certificate_omega_star_budget import _record as omega_star_record
from tests.tools.test_coverage_certificate_omega_star_budget import _run as omega_star_budget_run
from tools.coverage_certificate_budget_gate import BUDGET_RUNS
from tools.coverage_certificate_budget_gate import PRODUCTION_WORKLOAD
from tools.coverage_certificate_budget_gate import ProductionBudget
from tools.coverage_certificate_omega_star import load_record as load_omega_star_record
from tools.held_generation_clock import QUIET_LOAD_AVERAGE_MAX
from tools.perf_history import REGRESSION_EXIT_STATUS
from tools.perf_history import main
from tools.perf_history_adapters import ADAPTER_VERSIONS
from tools.perf_history_adapters import FAMILIES
from tools.perf_history_adapters import NON_TIMING_RECORDS
from tools.perf_history_adapters import PENDING_FAMILIES
from tools.perf_history_adapters import NotATimingRecordError
from tools.perf_history_adapters import PendingRecordFamilyError
from tools.perf_history_adapters import UnknownRecordFamilyError
from tools.perf_history_adapters import UnknownRecordSchemaError
from tools.perf_history_adapters import adapt
from tools.perf_history_baseline import REGRESSION_RATIO_MAX
from tools.perf_history_baseline import ClockStatisticMismatchError
from tools.perf_history_baseline import GatedMetricMismatchError
from tools.perf_history_baseline import PerformanceRegressionError
from tools.perf_history_baseline import baselines
from tools.perf_history_baseline import check
from tools.perf_history_baseline import current_entries
from tools.perf_history_entry import ClockSummary
from tools.perf_history_entry import EntryDigestMismatchError
from tools.perf_history_entry import LoadAverage1m
from tools.perf_history_entry import Measurement
from tools.perf_history_entry import NaiveTimestampError
from tools.perf_history_entry import NoClockError
from tools.perf_history_entry import PerfEntry
from tools.perf_history_entry import RecordProvenance
from tools.perf_history_entry import RecordSource
from tools.perf_history_ingest import backfill_entries
from tools.perf_history_ingest import file_entries
from tools.perf_history_ingest import record_entries
from tools.perf_history_ingest import record_versions
from tools.perf_history_json import MalformedTimingRecordError
from tools.perf_history_rebaseline import Rebaseline
from tools.perf_history_rebaseline import RebaselineKeyMismatchError
from tools.perf_history_rebaseline import RebaselineReasonError
from tools.perf_history_rebaseline import RebaselineTargetError
from tools.perf_history_render import GENERATED_BY
from tools.perf_history_render import render_markdown
from tools.perf_history_store import AppendOnlyViolationError
from tools.perf_history_store import History
from tools.perf_history_store import append_entries
from tools.perf_history_store import append_rebaseline
from tools.perf_history_store import load_history

ROOT = Path(__file__).resolve().parents[2]
RESULTS = "benchmarks/results"
# Every record present under benchmarks/results (the committed tree, or an export of it).
COMMITTED_RECORDS = tuple(sorted(path.relative_to(ROOT).as_posix() for path in (ROOT / RESULTS).rglob("*.json")))
COMMITTED_TIMING_RECORDS = tuple(path for path in COMMITTED_RECORDS if PurePosixPath(path).name in FAMILIES)
MACHINE = "Apple M1 Max"
QUIET = LoadAverage1m(QUIET_LOAD_AVERAGE_MAX)
LOADED = LoadAverage1m(QUIET_LOAD_AVERAGE_MAX + 0.5)
BASE_SECONDS = 10.0
# A commit whose pipeline_timing.json is schema /1 (the last one before the PENCIL unit wrote /2).
PIPELINE_TIMING_SCHEMA_1_REVISION = "1fcda9c4"


def _no_sibling(name: str) -> bytes:
    raise AssertionError(f"no sibling expected, asked for {name}")


def _entry(
    seconds: float = BASE_SECONDS,
    *,
    load: Optional[LoadAverage1m] = QUIET,
    recorded_at: str = "2026-10-23T10:00:00+00:00",
    record_sha256: str = "a" * 64,
    cpu: bool = True,
    workload: str = "figure8_monstera_boundary_qualified_corners",
) -> PerfEntry:
    clock = ClockSummary.over_runs([seconds, seconds * 1.01, seconds * 1.02])
    measurement = Measurement(workload, "held_generation_clock/generation", clock if cpu else None, None if cpu else clock, load, load)
    source = RecordSource(family="test", record_path="benchmarks/results/test.json", record_sha256=record_sha256)
    return PerfEntry.build(source, RecordProvenance(commit="c" * 40, dirty=False, machine=MACHINE, recorded_at=recorded_at), measurement)


# --- classification ----------------------------------------------------------------------------------


def test_every_committed_record_is_a_timing_family_or_a_listed_non_timing_record() -> None:
    unclassified = [path for path in COMMITTED_RECORDS if PurePosixPath(path).name not in FAMILIES and PurePosixPath(path).name not in NON_TIMING_RECORDS]
    assert not unclassified, f"records with no adapter and no non-timing listing: {unclassified}"


def test_every_family_has_a_committed_record() -> None:
    committed = {PurePosixPath(path).name for path in COMMITTED_TIMING_RECORDS}
    assert set(FAMILIES) <= committed


def test_the_pipeline_families_are_registered_and_a_pending_family_raises_its_own_error(monkeypatch: pytest.MonkeyPatch) -> None:
    assert {"pipeline_pencil_cost.json", "pipeline_timing.json"} <= set(FAMILIES) and not PENDING_FAMILIES
    monkeypatch.setattr(perf_history_adapters, "PENDING_FAMILIES", {"future_timing.json": "its adapter is not registered yet"})
    with pytest.raises(PendingRecordFamilyError, match="not registered yet"):
        adapt(f"{RESULTS}/future_timing.json", b"{}", _no_sibling)


def test_pipeline_timing_reads_schema_1_and_2_and_refuses_another_schema() -> None:
    schema_1 = subprocess.run(["git", "show", f"{PIPELINE_TIMING_SCHEMA_1_REVISION}:{RESULTS}/pipeline_timing.json"], cwd=ROOT, capture_output=True, check=True).stdout
    current = (ROOT / RESULTS / "pipeline_timing.json").read_bytes()
    assert [json.loads(schema_1)["schema"], json.loads(current)["schema"]] == ["pipeline-timing/1", "pipeline-timing/2"]
    assert record_entries(f"{RESULTS}/pipeline_timing.json", schema_1, _no_sibling)
    assert record_entries(f"{RESULTS}/pipeline_timing.json", current, _no_sibling)
    with pytest.raises(MalformedTimingRecordError, match="pipeline_timing"):
        adapt(f"{RESULTS}/pipeline_timing.json", json.dumps({"schema": "pipeline-timing/3"}).encode(), _no_sibling)


def test_an_unknown_record_family_raises() -> None:
    with pytest.raises(UnknownRecordFamilyError):
        adapt("benchmarks/results/brand_new_timing.json", b"{}", _no_sibling)


def test_a_non_timing_record_raises_rather_than_ingesting_nothing() -> None:
    with pytest.raises(NotATimingRecordError):
        adapt("benchmarks/results/held_path_digests.json", (ROOT / RESULTS / "held_path_digests.json").read_bytes(), _no_sibling)


def test_an_unknown_schema_of_a_known_family_raises() -> None:
    record = json.loads((ROOT / RESULTS / "coverage_certificate_omega_star.json").read_text())
    record["schema"] = "coverage-certificate-omega-star/99"
    with pytest.raises(UnknownRecordSchemaError):
        adapt("benchmarks/results/coverage_certificate_omega_star.json", json.dumps(record).encode(), _no_sibling)


OMEGA_STAR = f"{RESULTS}/coverage_certificate_omega_star.json"


def _omega_star_3(tmp_path: Path, budget: ProductionBudget) -> bytes:
    """A /3 record as the tool writes it, read back through the tool's own loader (``load_record``)."""
    path = tmp_path / "coverage_certificate_omega_star.json"
    path.write_text(omega_star_record(budget).to_json(), encoding="utf-8")
    return load_omega_star_record(path).to_json().encode()


def test_omega_star_3_yields_the_production_budget_beside_the_diagnostic(tmp_path: Path) -> None:
    budget = quiet_omega_star_budget()
    entries = {(entry.measurement.workload, entry.measurement.stage): entry for entry in record_entries(OMEGA_STAR, _omega_star_3(tmp_path, budget), _no_sibling)}
    family = "coverage_certificate_omega_star"
    assert set(entries) == {
        ("figure8_monstera", f"{family}/certificate"),
        ("figure8_monstera", f"{family}/target_build"),
        ("figure8_monstera", f"{family}/attribution"),
        (PRODUCTION_WORKLOAD, f"{family}/budget"),
    }
    gated = entries[(PRODUCTION_WORKLOAD, f"{family}/budget")]
    cpu = [run.cpu_s for run in budget.runs]
    assert gated.measurement.cpu == ClockSummary.over_runs(cpu) and gated.measurement.cpu.minimum == min(cpu) == budget.cpu.gated
    assert gated.measurement.wall == ClockSummary.over_runs([run.wall_s for run in budget.runs])
    assert gated.quiet and gated.measurement.load_before == LoadAverage1m(budget.runs[0].load_1min)


def test_an_indicative_omega_star_3_budget_yields_indicative_entries(tmp_path: Path) -> None:
    budget = ProductionBudget.build(tuple(omega_star_budget_run(4.0, load=BUDGET_LOAD_LOADED) for _ in range(BUDGET_RUNS)), indicative=True)
    gated = next(entry for entry in record_entries(OMEGA_STAR, _omega_star_3(tmp_path, budget), _no_sibling) if entry.measurement.stage.endswith("/budget"))
    assert not gated.quiet and gated.measurement.record_indicative is not None


def test_an_omega_star_3_budget_the_gate_refuses_is_refused(tmp_path: Path) -> None:
    raw = json.loads(_omega_star_3(tmp_path, quiet_omega_star_budget()))
    raw["budget"]["runs"] = raw["budget"]["runs"][:-1]
    with pytest.raises(MalformedTimingRecordError, match="budget"):
        adapt(OMEGA_STAR, json.dumps(raw).encode(), _no_sibling)


# --- adapters round-trip every committed record ------------------------------------------------------


@pytest.mark.parametrize("record", COMMITTED_TIMING_RECORDS, ids=[PurePosixPath(r).name for r in COMMITTED_TIMING_RECORDS])
def test_a_committed_record_round_trips_through_the_history_line(record: str) -> None:
    entries = file_entries(ROOT / record)
    assert entries, f"{record} yielded no entry"
    for entry in entries:
        assert PerfEntry.from_json_line(entry.to_json_line()) == entry
        assert entry.source.record_path == record
        assert entry.gated_clock.minimum >= 0


def test_a_record_without_load_is_indicative_and_says_so() -> None:
    entries = file_entries(ROOT / RESULTS / "held_toolpath_baseline.json")
    assert all(not entry.quiet and entry.indicative_reason == "no load recorded" for entry in entries)


def test_the_450_second_monstera_generation_is_in_the_baseline_record() -> None:
    entries = file_entries(ROOT / RESULTS / "held_toolpath_baseline.json")
    generation = [e for e in entries if e.key[:2] == ("figure8_monstera_boundary_qualified_corners", "held_toolpath_baseline/generation")]
    assert len(generation) == 1 and 449 < generation[0].gated_clock.minimum < 450


def test_backfill_reads_every_committed_version_of_every_record() -> None:
    versions = record_versions()
    result = backfill_entries(versions)
    assert len(result.timing_versions) + len(result.non_timing_versions) == len(versions)
    assert {PurePosixPath(v.path).name for v in result.timing_versions} == {PurePosixPath(p).name for p in COMMITTED_TIMING_RECORDS}
    assert all(entry.source.committed_in is not None for entry in result.entries)


# --- quiet rule ---------------------------------------------------------------------------------------


def test_quiet_is_inclusive_at_the_threshold_and_needs_both_loads() -> None:
    assert _entry(load=QUIET).quiet
    assert not _entry(load=LOADED).quiet
    assert _entry(load=None).indicative_reason == "no load recorded"


def test_a_measurement_without_any_clock_is_refused() -> None:
    source = RecordSource(family="test", record_path="x.json", record_sha256="a" * 64)
    with pytest.raises(NoClockError):
        PerfEntry.build(source, RecordProvenance(), Measurement("w", "s", None, None, QUIET, QUIET))


# --- append-only store --------------------------------------------------------------------------------


def test_appending_the_same_entries_twice_dedupes(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    entries = (_entry(), _entry(11.0, record_sha256="b" * 64))
    first = append_entries(entries, history)
    before = history.read_bytes()
    second = append_entries(entries, history)
    assert len(first.added) == 2 and not first.duplicates
    assert not second.added and len(second.duplicates) == 2
    assert history.read_bytes() == before


def test_existing_lines_never_change_on_append(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    append_entries([_entry()], history)
    before = history.read_bytes()
    append_entries([_entry(12.0, record_sha256="b" * 64)], history)
    after = history.read_bytes()
    assert after.startswith(before) and len(after) > len(before)
    assert load_history(history).entries == (_entry(), _entry(12.0, record_sha256="b" * 64))


def test_an_edited_line_is_refused(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    append_entries([_entry()], history)
    history.write_text(history.read_text().replace('"minimum":10.0', '"minimum":9.0'))
    with pytest.raises(EntryDigestMismatchError):
        load_history(history)


def test_an_append_that_rewrites_the_existing_bytes_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    history = tmp_path / "history.jsonl"
    append_entries([_entry()], history)
    original_open = Path.open

    def truncating_open(self: Path, mode: str = "r", *args: object, **kwargs: object) -> object:
        return original_open(self, "w" if mode == "a" else mode, *args, **kwargs)  # type: ignore[call-overload]

    monkeypatch.setattr(Path, "open", truncating_open)
    with pytest.raises(AppendOnlyViolationError):
        append_entries([_entry(12.0, record_sha256="b" * 64)], history)


# --- baseline and check -------------------------------------------------------------------------------


def _history(*entries: PerfEntry, rebaselines: tuple[Rebaseline, ...] = ()) -> History:
    return History(entries=entries, rebaselines=rebaselines)


def _sha(index: int) -> str:
    return f"{index:064x}"


def test_indicative_entries_never_set_a_baseline() -> None:
    quiet = _entry(recorded_at="2026-10-01T00:00:00+00:00")
    later_loaded = _entry(5.0, load=LOADED, recorded_at="2026-10-20T00:00:00+00:00", record_sha256="b" * 64)
    assert baselines(_history(quiet, later_loaded))[quiet.key] == quiet
    assert baselines(_history(later_loaded)) == {}


def test_the_baseline_is_the_best_quiet_entry_not_the_latest() -> None:
    best = _entry(8.0, recorded_at="2026-10-01T00:00:00+00:00")
    later_slower = _entry(9.0, recorded_at="2026-10-20T00:00:00+00:00", record_sha256="b" * 64)
    assert baselines(_history(later_slower, best))[best.key] == best


def test_a_regression_already_ingested_still_fails() -> None:
    baseline = _entry()
    regressed = _entry(2 * BASE_SECONDS, recorded_at="2026-10-24T00:00:00+00:00", record_sha256="b" * 64)
    with pytest.raises(PerformanceRegressionError):
        check(_history(baseline, regressed), [regressed])


def test_the_entry_under_test_is_never_its_own_baseline() -> None:
    only = _entry()
    (comparison,) = check(_history(only), [only])
    assert comparison.baseline is None


def test_an_identical_rerun_passes() -> None:
    baseline = _entry()
    rerun = _entry(recorded_at="2026-10-24T00:00:00+00:00", record_sha256="b" * 64)
    (comparison,) = check(_history(baseline, rerun), [rerun])
    assert comparison.baseline == baseline and not comparison.over_threshold and comparison.ratio == 1.0


def test_creeping_nine_percent_steps_fail_by_the_second_step() -> None:
    history = _history(_entry(recorded_at="2026-10-01T00:00:00+00:00"))
    failed_at = None
    for step in range(1, 6):
        run = _entry(BASE_SECONDS * 1.09**step, recorded_at=f"2026-10-{1 + step:02d}T00:00:00+00:00", record_sha256=_sha(step))
        history = _history(*history.entries, run)  # ingested first: the order must not matter
        try:
            check(history, [run])
        except PerformanceRegressionError:
            failed_at = step
            break
    assert failed_at is not None and failed_at <= 2


def test_a_rebaseline_makes_the_slower_run_the_baseline() -> None:
    best = _entry(recorded_at="2026-10-01T00:00:00+00:00")
    slower = _entry(BASE_SECONDS * 1.3, recorded_at="2026-10-05T00:00:00+00:00", record_sha256=_sha(1))
    rerun = _entry(BASE_SECONDS * 1.31, recorded_at="2026-10-06T00:00:00+00:00", record_sha256=_sha(2))
    with pytest.raises(PerformanceRegressionError):
        check(_history(best, slower), [rerun])
    accepted = Rebaseline.build(_history(best, slower), slower.entry_sha256, "accepted: exact predicate replaces float filter", "2026-10-05T12:00:00+00:00")
    history = _history(best, slower, rebaselines=(accepted,))
    assert baselines(history)[best.key] == slower
    (comparison,) = check(history, [rerun])
    assert comparison.baseline == slower and not comparison.over_threshold


def test_a_rebaseline_must_name_a_quiet_entry_of_the_history() -> None:
    loaded = _entry(load=LOADED)
    with pytest.raises(RebaselineTargetError):
        Rebaseline.build(_history(loaded), loaded.entry_sha256, "reason", "2026-10-05T12:00:00+00:00")
    with pytest.raises(RebaselineTargetError):
        Rebaseline.build(_history(_entry()), "f" * 64, "reason", "2026-10-05T12:00:00+00:00")
    with pytest.raises(RebaselineReasonError):
        Rebaseline.build(_history(_entry()), _entry().entry_sha256, " ", "2026-10-05T12:00:00+00:00")


def test_a_rebaseline_round_trips_through_the_history_file(tmp_path: Path) -> None:
    path = tmp_path / "history.jsonl"
    slower = _entry(BASE_SECONDS * 1.3)
    append_entries([slower], path)
    accepted = Rebaseline.build(load_history(path), slower.entry_sha256, "accepted", "2026-10-24T00:00:00+00:00")
    before = path.read_bytes()
    append_rebaseline(accepted, path)
    assert path.read_bytes().startswith(before)
    assert load_history(path).rebaselines == (accepted,)


def test_times_order_as_instants_across_offsets() -> None:
    # 23:00 at -05:00 is 04:00Z on the 21st: later than 01:00Z on the 21st, though it sorts first as text.
    later = _entry(recorded_at="2026-10-20T23:00:00-05:00")
    earlier = _entry(recorded_at="2026-10-21T01:00:00+00:00", record_sha256="b" * 64)
    assert later.when is not None and earlier.when is not None and later.when > earlier.when
    rebaseline = Rebaseline.build(_history(earlier, later), earlier.entry_sha256, "accepted", "2026-10-21T02:00:00+02:00")
    # the window starts at the named entry: the later (by instant) entry is a candidate, and it is faster here
    faster_later = _entry(9.0, recorded_at="2026-10-20T23:00:00-05:00", record_sha256="c" * 64)
    assert baselines(_history(earlier, faster_later, rebaselines=(rebaseline,)))[earlier.key] == faster_later


def test_a_naive_recorded_at_is_refused_when_the_entry_is_built() -> None:
    measurement = _entry().measurement
    source = RecordSource(family="test", record_path="benchmarks/results/test.json", record_sha256="a" * 64)
    with pytest.raises(NaiveTimestampError):
        PerfEntry.build(source, RecordProvenance(machine=MACHINE, recorded_at="2026-10-21T01:00:00"), measurement)
    naive_commit = RecordSource(family="test", record_path="benchmarks/results/test.json", record_sha256="a" * 64, committed_at="2026-10-21T01:00:00")
    with pytest.raises(NaiveTimestampError):
        PerfEntry.build(naive_commit, RecordProvenance(machine=MACHINE), measurement)


def test_a_record_with_a_naive_recorded_at_appends_nothing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    record = json.loads((ROOT / RESULTS / "held_splice_ab.json").read_text())
    record["build"]["recorded_at"] = "2026-10-23T08:23:12"
    path = tmp_path / RESULTS / "held_splice_ab.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(record))
    history = tmp_path / "history.jsonl"
    monkeypatch.setattr(perf_history_ingest, "REPOSITORY", tmp_path)
    with pytest.raises(NaiveTimestampError):
        main(["--history", str(history), "ingest", str(path)])
    assert not history.exists()


def test_a_naive_rebaseline_time_is_refused() -> None:
    entry = _entry()
    with pytest.raises(NaiveTimestampError):
        Rebaseline.build(_history(entry), entry.entry_sha256, "accepted", "2026-10-24T00:00:00")


def test_a_rebaseline_whose_key_differs_from_its_entry_is_refused_on_load(tmp_path: Path) -> None:
    path = tmp_path / "history.jsonl"
    entry = _entry()
    append_entries([entry], path)
    forged = Rebaseline.from_json(json.loads(Rebaseline.build(_history(entry), entry.entry_sha256, "accepted", "2026-10-24T00:00:00+00:00").to_json_line()))
    other_key = json.loads(forged.to_json_line())
    other_key["workload"] = "figure5_refined"
    del other_key["rebaseline_sha256"]
    resealed = Rebaseline._sealed((other_key["workload"], other_key["stage"], other_key["machine"]), forged.entry_sha256, forged.reason, forged.recorded_at)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(resealed.to_json_line())
    with pytest.raises(RebaselineKeyMismatchError):
        load_history(path)


def test_check_passes_just_under_the_threshold() -> None:
    (comparison,) = check(_history(_entry()), [_entry(BASE_SECONDS * (REGRESSION_RATIO_MAX - 0.01), record_sha256="b" * 64)])
    assert not comparison.over_threshold


def test_check_fails_on_an_injected_quiet_regression() -> None:
    regressed = _entry(BASE_SECONDS * (REGRESSION_RATIO_MAX + 0.05), record_sha256="b" * 64)
    with pytest.raises(PerformanceRegressionError) as raised:
        check(_history(_entry()), [regressed])
    assert raised.value.regressions[0].entry == regressed


def test_an_indicative_entry_over_the_threshold_is_reported_but_never_fails() -> None:
    (comparison,) = check(_history(_entry()), [_entry(BASE_SECONDS * 3, load=LOADED, record_sha256="b" * 64)])
    assert comparison.over_threshold and not comparison.regressed


def test_check_refuses_comparing_cpu_with_wall() -> None:
    with pytest.raises(GatedMetricMismatchError):
        check(_history(_entry()), [_entry(cpu=False, record_sha256="b" * 64)])


def test_check_refuses_comparing_a_minimum_with_a_median() -> None:
    baseline = _entry()
    median = Measurement(baseline.measurement.workload, baseline.measurement.stage, ClockSummary.reported_median(BASE_SECONDS, 3), None, QUIET, QUIET)
    entry = PerfEntry.build(RecordSource(family="test", record_path="benchmarks/results/test.json", record_sha256="b" * 64), baseline.provenance, median)
    with pytest.raises(ClockStatisticMismatchError):
        check(_history(baseline), [entry])


# --- adapter versions: a corrected adapter supersedes, never rewrites -----------------------------------


def test_actual_stock_profile_generation_is_process_cpu() -> None:
    entries = file_entries(ROOT / RESULTS / "actual_stock_profile.json")
    assert entries and all(e.measurement.cpu is not None and e.measurement.wall is None for e in entries)
    assert all(e.source.adapter_version == ADAPTER_VERSIONS["actual_stock_profile"] >= 2 for e in entries)


def _wall_v1(entry: PerfEntry) -> PerfEntry:
    m = entry.measurement
    source = RecordSource(entry.source.family, entry.source.record_path, entry.source.record_sha256, adapter_version=1)
    return PerfEntry.build(source, entry.provenance, Measurement(m.workload, m.stage, None, m.cpu, m.load_before, m.load_after))


def test_a_corrected_adapter_version_supersedes_the_old_entry_without_rewriting_it(tmp_path: Path) -> None:
    path = tmp_path / "history.jsonl"
    corrected = _entry()
    corrected = PerfEntry.build(RecordSource("test", corrected.source.record_path, corrected.source.record_sha256, adapter_version=2), corrected.provenance, corrected.measurement)
    wrong = _wall_v1(corrected)
    append_entries([wrong], path)
    before = path.read_bytes()
    result = append_entries([corrected], path)
    assert result.added == (corrected,) and path.read_bytes().startswith(before)
    history = load_history(path)
    assert history.entries == (wrong, corrected)
    assert current_entries(history.entries) == (corrected,)
    assert not append_entries([wrong, corrected], path).added


def test_a_legacy_schema_1_line_still_loads_as_adapter_version_1() -> None:
    first = (ROOT / RESULTS / "perf_history.jsonl").read_text().splitlines()[0]
    assert json.loads(first)["schema"] == "perf-history-entry/1"
    assert PerfEntry.from_json_line(first).source.adapter_version == 1


def test_rebackfilling_unchanged_families_appends_nothing_for_them(tmp_path: Path) -> None:
    path = tmp_path / "history.jsonl"
    path.write_bytes((ROOT / RESULTS / "perf_history.jsonl").read_bytes())
    result = append_entries(backfill_entries(record_versions()).entries, path)
    assert all(entry.source.adapter_version > 1 for entry in result.added)


def test_the_check_command_exits_non_zero_on_a_regression(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    record = json.loads((ROOT / RESULTS / "held_splice_ab.json").read_text())
    record["load_1min_before_after"] = [1.0, 1.5]  # quiet
    path = tmp_path / RESULTS / "held_splice_ab.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(record))
    history = tmp_path / "history.jsonl"
    monkeypatch.setattr(perf_history_ingest, "REPOSITORY", tmp_path)
    assert main(["--history", str(history), "ingest", str(path)]) == 0
    assert main(["--history", str(history), "check", str(path)]) == 0
    record["full_refinement_seconds"] *= REGRESSION_RATIO_MAX + 0.5
    record["build"]["recorded_at"] = "2026-10-24T00:00:00+00:00"
    path.write_text(json.dumps(record))
    assert main(["--history", str(history), "check", str(path)]) == REGRESSION_EXIT_STATUS


# --- rendering ----------------------------------------------------------------------------------------


def test_the_rendering_names_its_generator_and_shows_the_baseline() -> None:
    text = render_markdown(_history(_entry(), _entry(12.0, load=LOADED, recorded_at="2026-10-24T00:00:00+00:00", record_sha256="b" * 64)))
    assert text.startswith(GENERATED_BY)
    assert "| figure8_monstera_boundary_qualified_corners | held_generation_clock/generation | Apple M1 Max |" in text
    assert "10.000 s" in text and "12.000 s" in text and "| 1.20 |" in text and "| 1 / 2 |" in text
