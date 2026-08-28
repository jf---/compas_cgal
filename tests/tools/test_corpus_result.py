from __future__ import annotations

import datetime
import hashlib
import importlib
import importlib.util
import json
import pathlib
import subprocess
from typing import Any
from typing import Optional

import pytest

from benchmarks.cli import AGGREGATE_CORPUS_NAMES
from benchmarks.cli import DEFAULT_CAP_DEG
from benchmarks.cli import DEFAULT_TOOL_DIAMETER
from benchmarks.measurement import MeasurementRecord
from benchmarks.report import render_markdown


UTC = datetime.timezone.utc


def _artifact() -> Any:
    return importlib.import_module("tools.measurement_artifact")


def _corpus() -> Any:
    return importlib.import_module("tools.corpus_result")


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    return subprocess.run(["git", "-C", str(repository), *arguments], check=True, capture_output=True).stdout


def _repository(tmp_path: pathlib.Path) -> pathlib.Path:
    repository = tmp_path / "repository"
    repository.mkdir()
    _git(repository, "init", "-q")
    (repository / "pixi.lock").write_bytes(b"corpus lock\n")
    _git(repository, "add", "pixi.lock")
    _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "lock")
    return repository.resolve()


def _record() -> MeasurementRecord:
    return MeasurementRecord(
        name="rectangle",
        family="analytic",
        params={"width": 20.0},
        tool_diameter=DEFAULT_TOOL_DIAMETER,
        tea_cap_deg=DEFAULT_CAP_DEG,
        generate_seconds=1.0,
        certify_seconds=2.0,
        operations=3,
        cut_operations=2,
        stations=4,
        max_tea_deg=90.0,
        uncertified=0,
        truly_exceeding=0,
        unresolved=0,
        arrangement_vertices_final=5,
        max_coordinate_digits=6,
    )


def _input_payload() -> dict[str, object]:
    return {
        "name": "all",
        "tool_diameter": DEFAULT_TOOL_DIAMETER,
        "cap_deg": DEFAULT_CAP_DEG,
        "collect_digits": True,
        "aggregate_corpora": list(AGGREGATE_CORPUS_NAMES),
    }


def _source_for_commit(module: Any, repository: pathlib.Path, commit: str) -> Any:
    lock_bytes = _git(repository, "show", f"{commit}:pixi.lock")
    return module.SourceSnapshot.build(repository=repository, commit=commit, pixi_lock_sha256=hashlib.sha256(lock_bytes).hexdigest())


def _write_result(
    repository: pathlib.Path,
    *,
    commit: Optional[str] = None,
    started: Optional[datetime.datetime] = None,
    finished: Optional[datetime.datetime] = None,
    markdown: Optional[bytes] = None,
    json_payload: Optional[bytes] = None,
    input_payload: Optional[dict[str, object]] = None,
    argv: Optional[tuple[str, ...]] = None,
    results_root: Optional[pathlib.Path] = None,
) -> pathlib.Path:
    artifact = _artifact()
    source_commit = commit or _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    source = _source_for_commit(artifact, repository, source_commit)
    actual_started = started or datetime.datetime(2026, 8, 28, 12, 0, tzinfo=UTC)
    actual_finished = finished or actual_started + datetime.timedelta(seconds=1)
    logical_name = f"{actual_started.date().isoformat()}-{source_commit[:12]}"
    result = (results_root or repository / "benchmarks" / "results") / logical_name
    result.mkdir(parents=True)
    records = [_record()]
    payloads = {
        "benchmark_report.md": markdown if markdown is not None else render_markdown(records).encode("utf-8"),
        "benchmark_report.json": json_payload if json_payload is not None else (json.dumps([record.to_dict() for record in records], indent=2) + "\n").encode("utf-8"),
    }
    for name, payload in payloads.items():
        (result / name).write_bytes(payload)
    hidden_stage = f"benchmarks/results/.{logical_name}.stage-fixture"
    envelope = artifact.build_envelope(
        artifact_kind=artifact.ArtifactKind("benchmark-corpus-result/v1"),
        source=source,
        started=actual_started,
        finished=actual_finished,
        argv=argv or ("/python", "-m", "benchmarks.cli", "corpus", "--name", "all", "--out", hidden_stage),
        input_version=artifact.IdentityVersion("benchmark-corpus-input/v1"),
        input_payload=input_payload or _input_payload(),
        result_version=artifact.IdentityVersion("benchmark-corpus-result/v1"),
        payloads=payloads,
    )
    artifact.write_envelope(result, envelope)
    return result


def test_corpus_result_module_exists() -> None:
    assert importlib.util.find_spec("tools.corpus_result") is not None


def test_validate_result_accepts_authenticated_exact_column_reports(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    result = _write_result(repository)
    assert _corpus().validate_result(result) == datetime.datetime(2026, 8, 28, 12, 0, 1, tzinfo=UTC)


@pytest.mark.parametrize(
    ("damage", "markdown", "json_payload"),
    [
        ("markdown", b"not a corpus report\n", None),
        ("empty-markdown", b"", None),
        ("invalid-utf8", b"\xff", None),
        ("json-container", None, b"{}\n"),
        ("json-columns", None, b'[{"unknown": 1}]\n'),
        ("empty-json", None, b"[]\n"),
    ],
)
def test_validate_result_rejects_corpus_contract_damage(
    tmp_path: pathlib.Path,
    damage: str,
    markdown: Optional[bytes],
    json_payload: Optional[bytes],
) -> None:
    del damage
    repository = _repository(tmp_path)
    result = _write_result(repository, markdown=markdown, json_payload=json_payload)
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().validate_result(result)


def test_validate_result_rejects_authenticated_wrong_directory_name(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    result = _write_result(repository)
    renamed = result.with_name("wrong-name")
    result.rename(renamed)
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().validate_result(renamed)


def test_validate_result_accepts_authenticated_failed_record(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    failed = MeasurementRecord.failed("bad", "analytic", {}, DEFAULT_TOOL_DIAMETER, DEFAULT_CAP_DEG, "measured failure")
    records = [failed]
    result = _write_result(
        repository,
        markdown=render_markdown(records).encode("utf-8"),
        json_payload=(json.dumps([failed.to_dict()]) + "\n").encode("utf-8"),
    )
    assert _corpus().validate_result(result).tzinfo == UTC


def test_validate_result_rejects_symlink_payload(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    outside = repository / "outside-report.md"
    outside.write_text(render_markdown([_record()]), encoding="utf-8")
    result = _write_result(repository, markdown=outside.read_bytes())
    (result / "benchmark_report.md").unlink()
    (result / "benchmark_report.md").symlink_to(outside)
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().validate_result(result)


def test_validate_result_rejects_result_directory_symlink_laundering(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    result = _write_result(repository)
    real = result.with_name(".real-result")
    result.rename(real)
    result.symlink_to(real, target_is_directory=True)
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().validate_result(result)


@pytest.mark.parametrize(
    "input_payload",
    [
        {**_input_payload(), "name": "smoke"},
        {**_input_payload(), "tool_diameter": 2.0},
        {**_input_payload(), "cap_deg": 60.0},
        {**_input_payload(), "collect_digits": False},
        {**_input_payload(), "aggregate_corpora": list(reversed(AGGREGATE_CORPUS_NAMES))},
    ],
)
def test_validate_result_rejects_authenticated_wrong_effective_input(tmp_path: pathlib.Path, input_payload: dict[str, object]) -> None:
    repository = _repository(tmp_path)
    result = _write_result(repository, input_payload=input_payload)
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().validate_result(result)


def test_validate_result_accepts_input_object_in_any_key_order(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    expected = _input_payload()
    reordered = {key: expected[key] for key in reversed(tuple(expected))}
    result = _write_result(repository, input_payload=reordered)
    assert _corpus().validate_result(result).tzinfo == UTC


@pytest.mark.parametrize(
    "argv_tail",
    [
        ("-m", "benchmarks.cli", "corpus", "--name", "all"),
        ("-m", "benchmarks.cli", "corpus", "--name", "all", "--out", "one", "--out", "two"),
        ("-m", "benchmarks.cli", "corpus", "--name", "all", "--out", "/absolute/stage"),
        ("-m", "benchmarks.cli", "corpus", "--name", "all", "--out", "build/.stage"),
        ("-m", "benchmarks.cli", "corpus", "--name", "all", "--out", "benchmarks/results/stage"),
        ("-m", "benchmarks.cli", "corpus", "--name", "all", "--out", "benchmarks/results/../.stage"),
    ],
)
def test_validate_result_rejects_authenticated_wrong_argv(tmp_path: pathlib.Path, argv_tail: tuple[str, ...]) -> None:
    repository = _repository(tmp_path)
    result = _write_result(repository, argv=("/python", *argv_tail))
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().validate_result(result)


@pytest.mark.parametrize("damage", ["visible-file", "visible-symlink"])
def test_latest_result_fails_loud_on_visible_root_damage(tmp_path: pathlib.Path, damage: str) -> None:
    repository = _repository(tmp_path)
    results = repository / "benchmarks" / "results"
    results.mkdir(parents=True)
    if damage == "visible-file":
        (results / "damage").write_text("bad", encoding="utf-8")
    else:
        (results / "damage").symlink_to(repository / "pixi.lock")
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().latest_result(results)


def test_latest_result_selects_greatest_finished_time_not_name(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    old_commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    (repository / "next.txt").write_text("next", encoding="utf-8")
    _git(repository, "add", "next.txt")
    _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "next")
    new_commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    lexically_later = _write_result(
        repository,
        commit=old_commit,
        started=datetime.datetime(2026, 8, 29, tzinfo=UTC),
        finished=datetime.datetime(2026, 8, 29, 1, tzinfo=UTC),
    )
    finished_later = _write_result(
        repository,
        commit=new_commit,
        started=datetime.datetime(2026, 8, 28, tzinfo=UTC),
        finished=datetime.datetime(2026, 8, 30, 1, tzinfo=UTC),
    )
    assert lexically_later.name > finished_later.name
    assert _corpus().latest_result(repository / "benchmarks" / "results") == finished_later


def test_latest_result_ignores_hidden_stage_and_rejects_empty_root(tmp_path: pathlib.Path) -> None:
    results = tmp_path / "results"
    results.mkdir()
    (results / ".interrupted.stage-1").mkdir()
    with pytest.raises(_corpus().NoMeasuredResultError):
        _corpus().latest_result(results)


def test_latest_result_rejects_missing_root_and_malformed_visible_directory(tmp_path: pathlib.Path) -> None:
    missing = tmp_path / "missing"
    with pytest.raises(_corpus().NoMeasuredResultError):
        _corpus().latest_result(missing)
    malformed = tmp_path / "results"
    (malformed / "visible").mkdir(parents=True)
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().latest_result(malformed)


def test_latest_result_validates_each_visible_child_once_and_breaks_finish_tie_by_name(
    tmp_path: pathlib.Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repository = _repository(tmp_path)
    first_commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    (repository / "next.txt").write_text("next", encoding="utf-8")
    _git(repository, "add", "next.txt")
    _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "next")
    second_commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    finish = datetime.datetime(2026, 8, 30, tzinfo=UTC)
    first = _write_result(repository, commit=first_commit, started=datetime.datetime(2026, 8, 28, tzinfo=UTC), finished=finish)
    second = _write_result(repository, commit=second_commit, started=datetime.datetime(2026, 8, 28, tzinfo=UTC), finished=finish)
    calls: list[pathlib.Path] = []
    original = _corpus()._validate_result

    def counting_validate(path: pathlib.Path, *, repository: pathlib.Path, logical_name: str) -> datetime.datetime:
        calls.append(path)
        return original(path, repository=repository, logical_name=logical_name)

    monkeypatch.setattr(_corpus(), "_validate_result", counting_validate)
    expected = max((first, second), key=lambda path: path.name)
    assert _corpus().latest_result(repository / "benchmarks" / "results") == expected
    assert sorted(calls) == sorted((first, second))


def test_latest_result_binds_all_children_to_results_root_repository(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    repository = _repository(tmp_path)
    result = _write_result(repository)
    calls: list[tuple[pathlib.Path, pathlib.Path, str]] = []
    original = _corpus()._validate_result

    def bound_validate(path: pathlib.Path, *, repository: pathlib.Path, logical_name: str) -> datetime.datetime:
        calls.append((path, repository, logical_name))
        return original(path, repository=repository, logical_name=logical_name)

    def forbidden_public(path: pathlib.Path) -> datetime.datetime:
        raise AssertionError(f"latest_result rediscovered repository from child {path}")

    monkeypatch.setattr(_corpus(), "_validate_result", bound_validate)
    monkeypatch.setattr(_corpus(), "validate_result", forbidden_public)
    assert _corpus().latest_result(repository / "benchmarks" / "results") == result
    assert calls == [(result, repository, result.name)]


def test_latest_result_rejects_visible_nested_git_repository(tmp_path: pathlib.Path) -> None:
    repository = _repository(tmp_path)
    candidate = repository / "benchmarks" / "results" / "2026-08-28-aaaaaaaaaaaa"
    candidate.mkdir(parents=True)
    _git(candidate, "init", "-q")
    (candidate / "pixi.lock").write_text("nested\n", encoding="utf-8")
    _git(candidate, "add", "pixi.lock")
    _git(candidate, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "nested")
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().latest_result(repository / "benchmarks" / "results")


def test_latest_result_cannot_be_rebound_by_git_repository_at_results_root(tmp_path: pathlib.Path) -> None:
    outer = _repository(tmp_path)
    results = outer / "benchmarks" / "results"
    results.mkdir(parents=True)
    _git(results, "init", "-q")
    (results / "pixi.lock").write_text("nested results lock\n", encoding="utf-8")
    _git(results, "add", "pixi.lock")
    _git(results, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "nested results")
    _write_result(results, results_root=results)
    (results / "pixi.lock").unlink()
    with pytest.raises(_corpus().InvalidMeasuredResultError):
        _corpus().latest_result(results)
