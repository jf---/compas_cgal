"""The ratio history: ingest of ``perf-ratio/1`` sessions into ``perf-ratio-entry/1`` lines, and its render."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
from typing import Callable

import pytest
from tools import perf_ratio_history
from tools import perf_ratio_record
from tools.perf_history_adapters import FAMILIES
from tools.perf_history_adapters import NON_TIMING_RECORDS
from tools.perf_history_store import AppendOnlyViolationError
from tools.perf_history_store import DuplicateHistoryEntryError
from tools.perf_history_store import UnterminatedHistoryError
from tools.perf_ratio_history import GENERATED_BY
from tools.perf_ratio_history import SEMANTICS
from tools.perf_ratio_history import InjectedRecordError
from tools.perf_ratio_history import MalformedRatioEntryError
from tools.perf_ratio_history import RatioEntry
from tools.perf_ratio_history import RatioEntryDigestMismatchError
from tools.perf_ratio_history import ingest
from tools.perf_ratio_history import load_record
from tools.perf_ratio_history import main
from tools.perf_ratio_history import render_markdown
from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import LoadOneMinute
from tools.perf_ratio_record import PairRun
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import Ratio
from tools.perf_ratio_record import WorkloadRatio

REF_COMMIT = "15fbed70ebb016bdd8ef49453e04606d03e0b1fe"
HEAD_COMMIT = "0123456789abcdef0123456789abcdef01234567"
HOST = "bench-host"
WORKLOADS = ("gt-cap120-case-1", "polygon-square")
# (ref, head) CPU seconds per pair.
PAIR_CPU = ((2.0, 3.0), (1.0, 4.0), (4.0, 2.0))
OLDER = ("2026-10-31T09:00:00+00:00", "2026-10-31T09:30:00+00:00")
NEWER = ("2026-11-01T10:00:00+00:00", "2026-11-01T10:30:00+00:00")


def _arm(rev: str, name: str, dirty: bool = False) -> ArmIdentity:
    return ArmIdentity(
        rev=rev,
        dirty=dirty,
        patch_sha256=f"patch-{name}" if dirty else None,
        native_sha256=f"native-{name}",
        native_sources_sha256=f"sources-{name}",
        environment="default",
        python_version="3.12.7",
        stock_2_sha256=f"stock-{name}",
        pixi_lock_sha256=f"lock-{name}",
        compas_cgal_file=f"/trees/{name}/src/compas_cgal/__init__.py",
        stock_2_origin=f"/trees/{name}/build/_stock_2.so",
        editable_build_dir=f"/trees/{name}/build",
        interpreter_prefix=f"/trees/{name}/.pixi/envs/default",
    )


def _pairs(pair_cpu: tuple[tuple[float, float], ...] = PAIR_CPU) -> tuple[PairRun, ...]:
    orders = ((ArmName.REF, ArmName.HEAD), (ArmName.HEAD, ArmName.REF))
    return tuple(
        PairRun.build(
            order=orders[index % 2],
            ref_cpu=CpuSeconds(ref),
            head_cpu=CpuSeconds(head),
            load_before=LoadOneMinute(100.0 + index),
            load_after=LoadOneMinute(110.0 + index),
        )
        for index, (ref, head) in enumerate(pair_cpu)
    )


def _session(
    *,
    workloads: tuple[str, ...] = WORKLOADS,
    when: tuple[str, str] = NEWER,
    head_dirty: bool = False,
    head_digest: str = "digest-a",
    injected: Any = None,
    pair_cpu: tuple[tuple[float, float], ...] = PAIR_CPU,
) -> PerfRatioSession:
    return PerfRatioSession(
        ref_commit=REF_COMMIT,
        arms={ArmName.REF: _arm(REF_COMMIT, "ref"), ArmName.HEAD: _arm(HEAD_COMMIT, "head", head_dirty)},
        workloads=tuple(WorkloadRatio(workload=name, family="gt", pairs=_pairs(pair_cpu), ref_digest="digest-a", head_digest=head_digest) for name in workloads),
        host=HOST,
        cpu_model="Apple M1 Max",
        started_utc=when[0],
        ended_utc=when[1],
        harness_sha256="harness-digest",
        injected=injected,
    )


def _write(path: Path, session: PerfRatioSession) -> Path:
    path.write_text(session.dumps(), encoding="utf-8")
    return path


def _restamp(document: dict[str, Any]) -> str:
    """A line whose ``entry_sha256`` is re-derived from the (edited) content, so only the shape check can refuse it."""
    content = {name: value for name, value in document.items() if name != "entry_sha256"}
    document = dict(content, entry_sha256=hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest())
    return json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"


def _one_line(tmp_path: Path) -> dict[str, Any]:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "session.json", _session(workloads=WORKLOADS[:1])), history)
    document: dict[str, Any] = json.loads(history.read_text(encoding="utf-8"))
    return document


# --- ingest ------------------------------------------------------------------------------------------


def test_ingest_appends_one_entry_per_workload_keyed_by_the_reference(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    added = ingest(_write(tmp_path / "session.json", _session()), history)
    assert len(history.read_text(encoding="utf-8").splitlines()) == 2
    assert [entry.stage for entry in added] == ["ratio-vs-15fbed70"] * 2
    assert [entry.key for entry in added] == [(name, "ratio-vs-15fbed70", HOST) for name in WORKLOADS]
    assert load_record(history) == added


def test_an_entry_names_both_arms_their_runs_the_harness_and_the_cpu(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    session = _session(head_dirty=True)
    added = ingest(_write(tmp_path / "session.json", session), history)
    loaded = load_record(history)
    assert loaded == added
    for entry, workload in zip(loaded, session.workloads):
        assert entry.arms == session.arms
        assert entry.run == workload
        assert entry.pairs == workload.pairs
        assert (entry.harness_sha256, entry.cpu_model) == (session.harness_sha256, session.cpu_model)
        assert (entry.ref_commit, entry.head_rev, entry.head_dirty) == (REF_COMMIT, HEAD_COMMIT, True)
        assert (entry.ratio, entry.pair_ratio_median, entry.same_output) == (workload.ratio, workload.pair_ratio_median, workload.same_output)
        assert (entry.load_min, entry.load_max) == workload.load_range
    document = json.loads(history.read_text(encoding="utf-8").splitlines()[0])
    for arm in ArmName:
        stored = document["arms"][arm.value]
        identity = session.arms[arm]
        assert (stored["rev"], stored["native_sha256"], stored["native_sources_sha256"]) == (identity.rev, identity.native_sha256, identity.native_sources_sha256)
        assert (stored["environment"], stored["python_version"], stored["stock_2_sha256"]) == (identity.environment, identity.python_version, identity.stock_2_sha256)
    assert [(pair["ref_cpu"], pair["head_cpu"]) for pair in document["run"]["pairs"]] == list(PAIR_CPU)
    assert [(pair["load_before"], pair["load_after"]) for pair in document["run"]["pairs"]] == [(100.0 + index, 110.0 + index) for index in range(len(PAIR_CPU))]
    assert (document["harness_sha256"], document["cpu_model"]) == (session.harness_sha256, session.cpu_model)


def test_the_stage_is_derived_by_the_record_definition(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # One source for the stage name: a longer reference prefix in the record must still load from the history.
    monkeypatch.setattr(perf_ratio_record, "REF_SHORT_CHARS", 12)
    history = tmp_path / "history.jsonl"
    added = ingest(_write(tmp_path / "session.json", _session()), history)
    assert [entry.stage for entry in load_record(history)] == [f"ratio-vs-{REF_COMMIT[:12]}"] * len(added)


def test_reingesting_a_record_appends_nothing(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    record = _write(tmp_path / "session.json", _session())
    ingest(record, history)
    before = history.read_bytes()
    assert ingest(record, history) == ()
    assert history.read_bytes() == before


def test_an_injected_record_is_refused_at_ingest(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    record = _write(tmp_path / "session.json", _session(injected=Injection(ArmName.HEAD, CpuSeconds(1.0))))
    with pytest.raises(InjectedRecordError):
        ingest(record, history)
    assert not history.exists()


def test_the_cli_refuses_an_injected_record_before_appending_any(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    history = tmp_path / "history.jsonl"
    monkeypatch.setattr(perf_ratio_history, "JSON_PATH", history)
    honest = _write(tmp_path / "honest.json", _session())
    injected = _write(tmp_path / "injected.json", _session(when=OLDER, injected=Injection(ArmName.REF, CpuSeconds(1.0))))
    with pytest.raises(InjectedRecordError):
        main(["ingest", str(honest), str(injected)])
    assert not history.exists()
    assert main(["ingest", str(honest)]) == 0
    assert main(["ingest", str(honest)]) == 0
    assert capsys.readouterr().out.splitlines() == ["appended 2 entries", "appended 0 entries"]


def test_the_same_record_named_twice_in_one_ingest_appends_its_entries_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    history = tmp_path / "history.jsonl"
    monkeypatch.setattr(perf_ratio_history, "JSON_PATH", history)
    record = _write(tmp_path / "session.json", _session())
    assert main(["ingest", str(record), str(record)]) == 0
    assert capsys.readouterr().out.splitlines() == [f"appended {len(WORKLOADS)} entries"]
    assert len(load_record(history)) == len(WORKLOADS)


def test_an_append_that_rewrites_the_existing_bytes_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "older.json", _session(when=OLDER)), history)
    original_open = Path.open

    def truncating_open(self: Path, mode: str = "r", *args: object, **kwargs: object) -> object:
        return original_open(self, "w" if mode == "a" else mode, *args, **kwargs)  # type: ignore[call-overload]

    monkeypatch.setattr(Path, "open", truncating_open)
    with pytest.raises(AppendOnlyViolationError):
        ingest(_write(tmp_path / "newer.json", _session()), history)


# --- load --------------------------------------------------------------------------------------------


def _edit_head_cpu(document: dict[str, Any]) -> None:
    document["run"]["pairs"][0]["head_cpu"] *= 0.5


def _edit_load(document: dict[str, Any]) -> None:
    document["run"]["pairs"][1]["load_after"] += 1.0


def _edit_native(document: dict[str, Any]) -> None:
    document["arms"]["head"]["native_sha256"] = "native-other"


def _edit_interpreter(document: dict[str, Any]) -> None:
    document["arms"]["ref"]["python_version"] = "3.9.0"


def _edit_harness(document: dict[str, Any]) -> None:
    document["harness_sha256"] = "harness-other"


def _edit_cpu_model(document: dict[str, Any]) -> None:
    document["cpu_model"] = "another cpu"


# Each edit changes a field the entry carries without re-stamping its digest.
EDITS: dict[str, Callable[[dict[str, Any]], None]] = {
    "a pair's head CPU": _edit_head_cpu,
    "a pair's load": _edit_load,
    "an arm's native digest": _edit_native,
    "an arm's interpreter": _edit_interpreter,
    "the harness digest": _edit_harness,
    "the cpu model": _edit_cpu_model,
}


@pytest.mark.parametrize("edit", sorted(EDITS))
def test_an_edited_history_line_is_refused(tmp_path: Path, edit: str) -> None:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "session.json", _session()), history)
    lines = history.read_text(encoding="utf-8").splitlines(keepends=True)
    document = json.loads(lines[0])
    EDITS[edit](document)
    history.write_text(json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n" + lines[1], encoding="utf-8")
    with pytest.raises(RatioEntryDigestMismatchError):
        load_record(history)


def test_an_unterminated_history_is_refused(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "session.json", _session()), history)
    history.write_text(history.read_text(encoding="utf-8").rstrip("\n"), encoding="utf-8")
    with pytest.raises(UnterminatedHistoryError):
        load_record(history)


def test_a_repeated_line_is_refused(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "session.json", _session()), history)
    first = history.read_text(encoding="utf-8").splitlines(keepends=True)[0]
    with history.open("a", encoding="utf-8") as handle:
        handle.write(first)
    with pytest.raises(DuplicateHistoryEntryError):
        load_record(history)


def test_a_missing_history_loads_as_empty(tmp_path: Path) -> None:
    assert load_record(tmp_path / "absent.jsonl") == ()


# Each edit is re-stamped with a matching digest, so only the named shape check stands between it and the history.
def _with_arm(document: dict[str, Any], arm: str, **fields: object) -> dict[str, Any]:
    return dict(document, arms=dict(document["arms"], **{arm: dict(document["arms"][arm], **fields)}))


def _without_arm_field(document: dict[str, Any], arm: str, name: str) -> dict[str, Any]:
    stripped = {key: value for key, value in document["arms"][arm].items() if key != name}
    return dict(document, arms=dict(document["arms"], **{arm: stripped}))


def _with_run(document: dict[str, Any], **fields: object) -> dict[str, Any]:
    return dict(document, run=dict(document["run"], **fields))


def _with_pair(document: dict[str, Any], **fields: object) -> dict[str, Any]:
    pairs = [dict(document["run"]["pairs"][0], **fields), *document["run"]["pairs"][1:]]
    return _with_run(document, pairs=pairs)


MALFORMED_LINES: dict[str, Callable[[dict[str, Any]], str]] = {
    "not json": lambda document: "{not json\n",
    "not an object": lambda document: "7\n",
    "another schema": lambda document: _restamp(dict(document, schema="perf-ratio-entry/0")),
    "a missing key": lambda document: _restamp({name: value for name, value in document.items() if name != "cpu_model"}),
    "an extra key": lambda document: _restamp(dict(document, note="hand-added")),
    "a text field holding a number": lambda document: _restamp(dict(document, machine=7)),
    "a time without a UTC offset": lambda document: _restamp(dict(document, recorded_at="2026-11-01T10:30:00")),
    "a time that is not ISO 8601": lambda document: _restamp(dict(document, recorded_at="yesterday")),
    "a stage that is not the reference's": lambda document: _restamp(dict(document, stage="ratio-vs-deadbeef")),
    "arms that are not an object": lambda document: _restamp(dict(document, arms=[])),
    "a missing arm": lambda document: _restamp(dict(document, arms={"ref": document["arms"]["ref"]})),
    "an arm without its native digest": lambda document: _restamp(_without_arm_field(document, "head", "native_sha256")),
    "an arm whose interpreter is a number": lambda document: _restamp(_with_arm(document, "ref", python_version=3.12)),
    "a dirty arm without a patch digest": lambda document: _restamp(_with_arm(document, "head", dirty=True)),
    "a run that is not an object": lambda document: _restamp(dict(document, run="gt")),
    "a run without pairs": lambda document: _restamp(_with_run(document, pairs=[])),
    "a stored ratio its pairs do not give": lambda document: _restamp(_with_run(document, ratio=document["run"]["ratio"] * 2.0)),
    "a stored same_output its digests do not give": lambda document: _restamp(_with_run(document, same_output=not document["run"]["same_output"])),
    "a pair CPU that is not positive": lambda document: _restamp(_with_pair(document, ref_cpu=0.0)),
    "a pair load holding text": lambda document: _restamp(_with_pair(document, load_before="100")),
    "a CPU written as an integer": lambda document: _restamp(_with_pair(document, ref_cpu=int(document["run"]["pairs"][0]["ref_cpu"]))),
}


@pytest.mark.parametrize("edit", sorted(MALFORMED_LINES))
def test_a_malformed_history_line_is_refused_by_name(tmp_path: Path, edit: str) -> None:
    line = MALFORMED_LINES[edit](_one_line(tmp_path))
    with pytest.raises(MalformedRatioEntryError):
        RatioEntry.from_json_line(line)


def test_a_restamped_line_round_trips(tmp_path: Path) -> None:
    document = _one_line(tmp_path)
    assert RatioEntry.from_json_line(_restamp(document)).to_json_line() == _restamp(document)


# --- render ------------------------------------------------------------------------------------------


def test_render_states_the_semantics_and_shows_only_the_latest_entry_per_key(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "older.json", _session(when=OLDER, pair_cpu=((1.0, 1.9871),))), history)
    ingest(_write(tmp_path / "newer.json", _session(workloads=WORKLOADS[:1], when=NEWER, pair_cpu=((1.0, 1.2341),))), history)
    entries = load_record(history)
    assert len(entries) == 3
    assert entries[-1].ratio == Ratio(1.2341)
    markdown = render_markdown(entries)
    assert markdown.startswith(GENERATED_BY + "\n\n" + SEMANTICS)
    rows = [line for line in markdown.splitlines() if line.startswith("| ") and not line.startswith(("| Workload ", "| --- "))]
    # One row per key, keys sorted: the newer entry replaces the older one of its key; the other key keeps its only entry.
    assert [row.split("|")[1].strip() for row in rows] == sorted(WORKLOADS)
    first = next(row for row in rows if row.startswith(f"| {WORKLOADS[0]} "))
    assert "2026-11-01" in first and "1.234×" in first
    assert "2026-10-31" not in first and "1.987×" not in first
    second = next(row for row in rows if row.startswith(f"| {WORKLOADS[1]} "))
    assert "2026-10-31" in second and "1.987×" in second
    assert render_markdown(tuple(reversed(entries))) == markdown


def test_render_marks_a_dirty_head_and_different_output(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    ingest(_write(tmp_path / "session.json", _session(workloads=WORKLOADS[:1], head_dirty=True, head_digest="digest-b")), history)
    row = next(line for line in render_markdown(load_record(history)).splitlines() if line.startswith(f"| {WORKLOADS[0]} "))
    cells = [cell.strip() for cell in row.strip("|").split("|")]
    header = [cell.strip() for cell in render_markdown(()).splitlines()[4].strip("|").split("|")]
    assert header == ["Workload", "Reference", "Machine", "Date", "HEAD", "Ratio", "Pair-ratio median", "Same output", "1-min load"]
    assert cells[header.index("HEAD")] == f"`{HEAD_COMMIT[:8]}`+dirty"
    assert cells[header.index("Same output")] == "no"


def test_write_markdown_renders_the_history(tmp_path: Path) -> None:
    history = tmp_path / "history.jsonl"
    markdown = tmp_path / "history.md"
    ingest(_write(tmp_path / "session.json", _session()), history)
    text = perf_ratio_history.write_markdown(history, markdown)
    assert markdown.read_text(encoding="utf-8") == text == render_markdown(load_record(history))


# --- ownership ---------------------------------------------------------------------------------------


def test_the_module_follows_the_rendered_record_protocol() -> None:
    assert all(hasattr(perf_ratio_history, name) for name in ("SCHEMA", "JSON_PATH", "MD_PATH", "load_record", "render_markdown"))


def test_the_session_record_is_owned_outside_the_absolute_history() -> None:
    assert "perf_ratio_session.json" in NON_TIMING_RECORDS
    assert "perf_ratio_session.json" not in FAMILIES
