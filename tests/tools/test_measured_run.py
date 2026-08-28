from __future__ import annotations

import datetime
import ast
import importlib
import importlib.util
import json
import pathlib
import subprocess
import sys
from typing import Any

import pytest

from benchmarks.cli import AGGREGATE_CORPUS_NAMES
from benchmarks.cli import DEFAULT_CAP_DEG
from benchmarks.cli import DEFAULT_TOOL_DIAMETER
from benchmarks.measurement import MeasurementRecord
from benchmarks.report import write_report


UTC = datetime.timezone.utc


def _module() -> Any:
    return importlib.import_module("tools.measured_run")


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    return subprocess.run(["git", "-C", str(repository), *arguments], check=True, capture_output=True).stdout


def _repository(tmp_path: pathlib.Path) -> pathlib.Path:
    repository = tmp_path / "repository"
    repository.mkdir()
    _git(repository, "init", "-q")
    (repository / "pixi.lock").write_bytes(b"runner lock\n")
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


def test_measured_run_module_exists() -> None:
    assert importlib.util.find_spec("tools.measured_run") is not None


def test_plan_public_imports_are_direct_forwarding_objects() -> None:
    measured_run = _module()
    corpus_result = importlib.import_module("tools.corpus_result")
    artifact = importlib.import_module("tools.measurement_artifact")
    assert measured_run.latest_result is corpus_result.latest_result
    assert measured_run.validate_result is corpus_result.validate_result
    assert measured_run.NoMeasuredResultError is corpus_result.NoMeasuredResultError
    assert measured_run.InvalidMeasuredResultError is corpus_result.InvalidMeasuredResultError
    assert measured_run.STAMP_KEYS is artifact.ENVELOPE_KEYS
    assert measured_run.MeasuredRunError is artifact.MeasurementArtifactError
    assert measured_run.DirtyMeasuredRunError is artifact.DirtyMeasurementTreeError
    assert measured_run.MeasuredRunInputChangedError is artifact.MeasurementInputChangedError
    assert measured_run.MeasuredResultCollisionError is artifact.MeasurementArtifactCollisionError


def test_main_runs_fixed_all_corpus_from_repository_and_publishes(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    measured_run = _module()
    repository = _repository(tmp_path)
    calls: list[tuple[tuple[str, ...], pathlib.Path]] = []
    times = iter(
        (
            datetime.datetime(2026, 8, 28, 12, 0, tzinfo=UTC),
            datetime.datetime(2026, 8, 28, 12, 1, tzinfo=UTC),
        )
    )

    def run_child(command: tuple[str, ...], cwd: pathlib.Path) -> int:
        calls.append((command, cwd))
        out = repository / pathlib.Path(command[-1])
        write_report([_record()], out)
        return 0

    monkeypatch.setattr(measured_run, "_repository_root", lambda: repository)
    monkeypatch.setattr(measured_run, "_utc_now", lambda: next(times))
    monkeypatch.setattr(measured_run, "_run_child", run_child)
    assert measured_run.main([]) == 0
    commit = _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    final = repository / "benchmarks" / "results" / f"2026-08-28-{commit[:12]}"
    assert capsys.readouterr().out == f"benchmarks/results/{final.name}\n"
    assert len(calls) == 1
    command, cwd = calls[0]
    assert cwd == repository
    assert command[:-1] == (sys.executable, "-m", "benchmarks.cli", "corpus", "--name", "all", "--out")
    stage = pathlib.PurePosixPath(command[-1])
    assert stage.parent == pathlib.PurePosixPath("benchmarks/results")
    assert stage.name.startswith(f".{final.name}.stage-")
    stamp = json.loads((final / "stamp.json").read_text(encoding="utf-8"))
    assert stamp["input_identity"]["payload"] == {
        "name": "all",
        "tool_diameter": DEFAULT_TOOL_DIAMETER,
        "cap_deg": DEFAULT_CAP_DEG,
        "collect_digits": True,
        "aggregate_corpora": list(AGGREGATE_CORPUS_NAMES),
    }
    assert measured_run.validate_result(final) == datetime.datetime(2026, 8, 28, 12, 1, tzinfo=UTC)


def test_main_preserves_positive_child_exit_without_partial_output(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    measured_run = _module()
    repository = _repository(tmp_path)
    monkeypatch.setattr(measured_run, "_repository_root", lambda: repository)
    monkeypatch.setattr(measured_run, "_run_child", lambda command, cwd: 37)
    assert measured_run.main([]) == 37
    results = repository / "benchmarks" / "results"
    assert not results.exists() or list(results.iterdir()) == []


@pytest.mark.parametrize("damage", ["missing-report", "dirty-tree", "changed-head"])
def test_main_zero_exit_damage_fails_without_visible_result(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, damage: str) -> None:
    measured_run = _module()
    repository = _repository(tmp_path)

    def damaged_child(command: tuple[str, ...], cwd: pathlib.Path) -> int:
        stage = repository / pathlib.Path(command[-1])
        if damage != "missing-report":
            write_report([_record()], stage)
        else:
            stage.mkdir(parents=True, exist_ok=True)
            (stage / "benchmark_report.md").write_text("# Benchmark corpus result\n", encoding="utf-8")
        if damage == "dirty-tree":
            (repository / "unexpected.txt").write_text("dirty", encoding="utf-8")
        elif damage == "changed-head":
            (repository / "tracked.txt").write_text("changed", encoding="utf-8")
            _git(repository, "add", "tracked.txt")
            _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "changed")
        return 0

    monkeypatch.setattr(measured_run, "_repository_root", lambda: repository)
    monkeypatch.setattr(measured_run, "_run_child", damaged_child)
    with pytest.raises(measured_run.MeasuredRunError):
        measured_run.main([])
    results = repository / "benchmarks" / "results"
    assert not results.exists() or all(path.name.startswith(".") for path in results.iterdir())
    assert not results.exists() or not any(path.name[0].isdigit() for path in results.iterdir())


@pytest.mark.parametrize("outcome", ["signal", "spawn"])
def test_child_signal_and_spawn_failure_raise_named_error(monkeypatch: pytest.MonkeyPatch, outcome: str) -> None:
    measured_run = _module()

    def failed_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        if outcome == "spawn":
            raise OSError("cannot spawn")
        return subprocess.CompletedProcess(["child"], -9)

    monkeypatch.setattr(measured_run.subprocess, "run", failed_run)
    with pytest.raises(measured_run.MeasuredRunChildError):
        measured_run._run_child(("child",), pathlib.Path.cwd())


def test_docs_and_cli_explain_committed_measurement_contract() -> None:
    repository = pathlib.Path(__file__).parents[2]
    cli = (repository / "benchmarks" / "cli.py").read_text(encoding="utf-8")
    docs = (repository / "docs" / "benchmarks.md").read_text(encoding="utf-8")
    normalized_docs = " ".join(docs.split())
    assert "pixi run measured-run" in cli
    assert "benchmarks/results/" in cli
    assert "benchmarks/measurement_claim_results/" in normalized_docs
    assert "reporting evidence" in normalized_docs
    assert "single writer" in normalized_docs
    assert "corpus --name all" in normalized_docs
    assert "clean committed worktree" in normalized_docs
    assert "rechecks full HEAD and worktree cleanliness after the child, immediately before publication" in normalized_docs
    assert "build identity" in normalized_docs
    assert "input identity" in normalized_docs
    assert "result identity" in normalized_docs
    assert "not geometric or continuous-engagement certificates" in normalized_docs
    assert "empty destination in the final race window" in normalized_docs


def test_report_writer_declares_utf8_at_both_payload_boundaries() -> None:
    repository = pathlib.Path(__file__).parents[2]
    report = (repository / "benchmarks" / "report.py").read_text(encoding="utf-8")
    assert report.count("write_text(") == 2
    assert report.count('encoding="utf-8"') == 2


def test_measured_run_pixi_task_uses_editable_build_and_module_entry() -> None:
    repository = pathlib.Path(__file__).parents[2]
    manifest = (repository / "pyproject.toml").read_text(encoding="utf-8")
    assert "measured-run =" in manifest
    assert "python -m tools.measured_run" in manifest
    assert 'depends-on = ["_editable-rebuild"]' in manifest


def test_new_modules_preserve_python39_syntax_and_api_floor() -> None:
    repository = pathlib.Path(__file__).parents[2]
    paths = (
        repository / "tools" / "measurement_artifact.py",
        repository / "tools" / "corpus_result.py",
        repository / "tools" / "measured_run.py",
        repository / "tests" / "tools" / "test_measurement_artifact.py",
        repository / "tests" / "tools" / "test_corpus_result.py",
        repository / "tests" / "tools" / "test_measured_run.py",
    )
    forbidden_apis = tuple(
        left + right
        for left, right in (
            ("datetime", ".UTC"),
            ("Path", ".walk"),
            ("hashlib", ".file_digest"),
            ("tom", "llib"),
            ("typing", ".Self"),
            ("slots", "=True"),
            ("kw_only", "=True"),
        )
    )
    for path in paths:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path), feature_version=(3, 9))
        assert not any(isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr) for node in ast.walk(tree)), path
        for forbidden in forbidden_apis:
            assert forbidden not in source, f"{path}: Python 3.9 forbids {forbidden}"
