from __future__ import annotations

import contextlib
import datetime
import importlib
import importlib.util
import pathlib
import signal
from types import SimpleNamespace
from typing import Any

import pytest

from tools.measurement_artifact import DirtyMeasurementTreeError
from tools.measurement_artifact import MeasurementInputChangedError
from tools.measurement_artifact import SourceSnapshot
from tools.measurement_claim_errors import MeasurementClaimChildError


EXPECTED_CASE_ORDER = (
    "radial-station",
    "radial-subdivisions",
    "radial-floor",
    "radial-margin",
    "advance-placement",
    "advance-probe-count",
)


def _probes() -> Any:
    return importlib.import_module("tools.measurement_claim_probes")


def test_measurement_claim_probes_module_exists() -> None:
    assert importlib.util.find_spec("tools.measurement_claim_probes") is not None


def test_generator_case_order_is_fixed() -> None:
    assert _probes().GENERATOR_CASE_ORDER == EXPECTED_CASE_ORDER


def test_run_cases_dispatches_serially_in_requested_order(monkeypatch: pytest.MonkeyPatch) -> None:
    probes = _probes()
    calls: list[str] = []

    def radial(case: str) -> dict[str, str]:
        calls.append(case)
        return {"case": case}

    def advance(case: str) -> dict[str, str]:
        calls.append(case)
        return {"case": case}

    monkeypatch.setattr(probes, "run_radial_case", radial)
    monkeypatch.setattr(probes, "run_advance_case", advance)

    results = probes._run_cases(EXPECTED_CASE_ORDER)

    assert calls == list(EXPECTED_CASE_ORDER)
    assert [result["case"] for result in results] == list(EXPECTED_CASE_ORDER)


@pytest.mark.parametrize(
    "selected",
    (
        (),
        ("radial-station", "radial-station"),
        ("radial-station", "unknown-case"),
        ("advance-placement", "radial-station"),
    ),
)
def test_invalid_selection_fails_before_runner(
    monkeypatch: pytest.MonkeyPatch,
    selected: tuple[str, ...],
) -> None:
    probes = _probes()
    calls: list[str] = []
    monkeypatch.setattr(probes, "run_radial_case", lambda case: calls.append(case))
    monkeypatch.setattr(probes, "run_advance_case", lambda case: calls.append(case))

    with pytest.raises(probes.InvalidMeasurementClaimConfigError):
        probes._run_cases(selected)

    assert calls == []


def test_list_prints_only_fixed_case_order(capsys: pytest.CaptureFixture[str]) -> None:
    assert _probes().main(["list"]) == 0
    assert capsys.readouterr().out.splitlines() == list(EXPECTED_CASE_ORDER)


def test_benchmark_cli_dispatches_the_named_case(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    probes = _probes()
    monkeypatch.chdir(tmp_path)
    calls: list[tuple[str, pathlib.Path, pathlib.Path]] = []
    result = tmp_path / "benchmarks" / "measurement_claim_results" / "result"
    monkeypatch.setattr(
        probes,
        "_run_benchmark_case",
        lambda case, repository, results_root: calls.append((case, repository, results_root)) or result,
    )
    assert probes.main(["run-benchmark", "--case", "figure6-spacing"]) == 0
    assert calls == [
        (
            "figure6-spacing",
            tmp_path.resolve(),
            tmp_path.resolve() / "benchmarks" / "measurement_claim_results",
        )
    ]
    assert capsys.readouterr().out == "benchmarks/measurement_claim_results/result\n"


def test_generator_preflight_failure_stops_before_any_runner(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> None:
    probes = _probes()
    calls: list[str] = []

    def fail_preflight(repository: pathlib.Path) -> None:
        del repository
        raise DirtyMeasurementTreeError("dirty")

    monkeypatch.setattr(probes, "capture_clean_source", fail_preflight)
    monkeypatch.setattr(probes, "_run_cases", lambda cases: calls.extend(cases))

    with pytest.raises(DirtyMeasurementTreeError, match="dirty"):
        probes._produce_generator(tmp_path, tmp_path / "benchmarks" / "measurement_claim_results")

    assert calls == []


def test_generator_publication_validates_stage_then_final_and_renders_only_final_tuple(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    probes = _probes()
    repository = tmp_path
    results = repository / "benchmarks" / "measurement_claim_results"
    started = datetime.datetime(2026, 8, 29, 12, 30, 1, 123456, tzinfo=datetime.timezone.utc)
    finished = started + datetime.timedelta(seconds=2)
    source = SourceSnapshot.build(repository=repository, commit="a" * 40, pixi_lock_sha256="b" * 64)
    payload = {"source_commit": "a" * 40, "cases": []}
    final_tuple = ("final-payload", "final-envelope", "final-start", "final-path")
    stage_tuple = ("stage-payload", "stage-envelope", "stage-start", "final-path")
    events: list[object] = []
    stage_path: pathlib.Path | None = None

    monkeypatch.setattr(probes, "capture_clean_source", lambda path: source)
    monkeypatch.setattr(probes, "_now_utc", lambda: started if not events else finished)
    monkeypatch.setattr(probes, "_run_cases", lambda cases: events.append(("run", tuple(cases))) or [])
    monkeypatch.setattr(probes, "compose_generator_payload", lambda commit, cases: payload)
    monkeypatch.setattr(probes, "generator_semantic_input", lambda value: {"semantic": value["source_commit"]})
    monkeypatch.setattr(
        probes,
        "build_envelope",
        lambda **kwargs: {
            "input_identity": {"sha256": "c" * 64},
            "built": events.append(("envelope", kwargs["argv"], kwargs["input_payload"])),
        },
    )
    monkeypatch.setattr(probes, "write_envelope", lambda stage, envelope: events.append(("write-envelope", stage, envelope)))

    @contextlib.contextmanager
    def stage(root: pathlib.Path, logical_name: str) -> Any:
        nonlocal stage_path
        root.mkdir(parents=True)
        stage_path = root / f".{logical_name}.stage-owned"
        stage_path.mkdir()
        events.append(("stage", logical_name))
        yield stage_path

    monkeypatch.setattr(probes, "publication_stage", stage)

    def validate(path: pathlib.Path) -> tuple[str, str, str, str]:
        events.append(("validate", path))
        return stage_tuple if path == stage_path else final_tuple

    monkeypatch.setattr(probes, "validate_claim_artifact", validate)

    def publish(*, source: object, stage: pathlib.Path, final: pathlib.Path) -> pathlib.Path:
        del source
        events.append(("publish", stage, final))
        stage.rename(final)
        return final

    monkeypatch.setattr(probes, "publish_stage", publish)
    rendered: list[tuple[object, ...]] = []
    monkeypatch.setattr(
        probes,
        "render_ledger_evidence",
        lambda value, envelope, *, started, artifact_directory: rendered.append((value, envelope, started, artifact_directory)),
    )

    final = probes._produce_generator(repository, results)

    expected = results / "2026-08-29-aaaaaaaaaaaa-generator-cccccccccccc"
    assert final == expected
    assert [event for event in events if isinstance(event, tuple) and event[0] == "validate"] == [
        ("validate", stage_path),
        ("validate", expected),
    ]
    assert rendered == [final_tuple]
    assert events.index(("validate", stage_path)) < events.index(("publish", stage_path, expected)) < events.index(("validate", expected))


def test_benchmark_publication_reuses_precomputed_semantics_digest_and_executed_argv(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    probes = _probes()
    repository = tmp_path
    results = repository / "benchmarks" / "measurement_claim_results"
    started = datetime.datetime(2026, 8, 29, 12, 30, 1, 123456, tzinfo=datetime.timezone.utc)
    source = SourceSnapshot.build(repository=repository, commit="a" * 40, pixi_lock_sha256="b" * 64)
    semantic = {"source_commit": "a" * 40, "command": ["figure6"]}
    payload = {"source_commit": "a" * 40}
    executed = ("python", "-m", "benchmarks.cli", "figure6", "--out", "owned-stage")
    events: list[object] = []

    monkeypatch.setattr(probes, "capture_clean_source", lambda path: source)
    monkeypatch.setattr(probes, "_now_utc", lambda: started)
    monkeypatch.setattr(probes, "benchmark_semantic_input", lambda commit: semantic)
    monkeypatch.setattr(probes, "benchmark_payload_semantic_input", lambda candidate: semantic)
    monkeypatch.setattr(probes, "input_identity_sha256", lambda **kwargs: "c" * 64)

    @contextlib.contextmanager
    def stage(root: pathlib.Path, logical_name: str) -> Any:
        path = root / f".{logical_name}.stage-owned"
        path.mkdir(parents=True)
        events.append(("stage", logical_name))
        yield path

    monkeypatch.setattr(probes, "publication_stage", stage)

    def run_case(path: pathlib.Path, *, source_commit: str) -> tuple[dict[str, str], tuple[str, ...]]:
        assert source_commit == "a" * 40
        (path / "figure6.md").write_bytes(b"markdown")
        (path / "figure6.json").write_bytes(b"{}\n")
        events.append(("case", path))
        return payload, executed

    monkeypatch.setattr(probes, "run_figure6_case", run_case)

    def envelope(**kwargs: object) -> dict[str, object]:
        events.append(("envelope", kwargs["argv"], kwargs["input_payload"]))
        return {"input_identity": {"sha256": "c" * 64}}

    monkeypatch.setattr(probes, "build_envelope", envelope)
    monkeypatch.setattr(probes, "write_envelope", lambda path, value: events.append(("write", path, value)))
    validated = ({}, object(), object(), "canonical")
    monkeypatch.setattr(probes, "validate_benchmark_claim_artifact", lambda path: events.append(("validate", path)) or validated)

    def publish(*, source: object, stage: pathlib.Path, final: pathlib.Path) -> pathlib.Path:
        del source
        events.append(("publish", stage, final))
        stage.rename(final)
        return final

    monkeypatch.setattr(probes, "publish_stage", publish)
    final = probes._produce_benchmark(repository, results)
    expected = results / "2026-08-29-aaaaaaaaaaaa-benchmark-cccccccccccc"
    assert final == expected
    assert ("envelope", executed, semantic) in events
    assert [event[0] for event in events if isinstance(event, tuple)] == ["stage", "case", "envelope", "write", "validate", "publish", "validate"]


def test_benchmark_preflight_failure_stops_before_stage_or_child(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> None:
    probes = _probes()
    events: list[str] = []

    def fail_preflight(repository: pathlib.Path) -> None:
        del repository
        raise DirtyMeasurementTreeError("dirty")

    monkeypatch.setattr(probes, "capture_clean_source", fail_preflight)
    monkeypatch.setattr(probes, "publication_stage", lambda *args: pytest.fail("stage created after dirty preflight"))
    monkeypatch.setattr(probes, "run_figure6_case", lambda *args, **kwargs: events.append("child"))

    with pytest.raises(DirtyMeasurementTreeError, match="dirty"):
        probes._produce_benchmark(tmp_path, tmp_path / "benchmarks" / "measurement_claim_results")

    assert events == []


@pytest.mark.parametrize(
    "failure",
    [
        OSError("spawn denied"),
        -signal.SIGTERM,
        7,
    ],
)
def test_benchmark_child_failure_stops_before_serialization_and_publication(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
    failure: object,
) -> None:
    probes = _probes()
    figure6 = importlib.import_module("tools.measurement_claim_figure6")
    results = tmp_path / "benchmarks" / "measurement_claim_results"
    source = SourceSnapshot.build(repository=tmp_path, commit="a" * 40, pixi_lock_sha256="b" * 64)
    monkeypatch.setattr(probes, "capture_clean_source", lambda path: source)
    monkeypatch.setattr(probes, "_now_utc", lambda: datetime.datetime(2026, 8, 29, tzinfo=datetime.timezone.utc))
    monkeypatch.setattr(probes, "benchmark_semantic_input", lambda commit: {"source_commit": commit})
    monkeypatch.setattr(probes, "input_identity_sha256", lambda **kwargs: "c" * 64)

    @contextlib.contextmanager
    def stage(root: pathlib.Path, logical_name: str) -> Any:
        path = root / f".{logical_name}.stage-owned"
        path.mkdir(parents=True)
        try:
            yield path
        finally:
            path.rmdir()

    monkeypatch.setattr(probes, "publication_stage", stage)

    def run(argv: tuple[str, ...], *, check: bool) -> SimpleNamespace:
        del argv, check
        if isinstance(failure, OSError):
            raise failure
        return SimpleNamespace(returncode=failure)

    monkeypatch.setattr(figure6.subprocess, "run", run)
    monkeypatch.setattr(figure6, "_read_raw_outputs", lambda path: pytest.fail(f"parsed failed child output: {path}"))
    for name in (
        "benchmark_payload_semantic_input",
        "_payload_bytes",
        "build_envelope",
        "write_envelope",
        "validate_benchmark_claim_artifact",
        "publish_stage",
    ):
        monkeypatch.setattr(probes, name, lambda *args, _name=name, **kwargs: pytest.fail(f"{_name} called after child failure"))

    with pytest.raises(MeasurementClaimChildError) as caught:
        probes._produce_benchmark(tmp_path, results)

    message = str(caught.value)
    assert ("spawn denied" in message) if isinstance(failure, OSError) else (str(failure) in message)
    if failure == -signal.SIGTERM:
        assert str(signal.SIGTERM) in message
    assert not tuple(results.iterdir())


def test_benchmark_semantic_mutation_stops_before_write_or_stamp(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> None:
    probes = _probes()
    results = tmp_path / "benchmarks" / "measurement_claim_results"
    source = SourceSnapshot.build(repository=tmp_path, commit="a" * 40, pixi_lock_sha256="b" * 64)
    semantic = {
        "extraction_commit": "e" * 40,
        "source_commit": "a" * 40,
        "semantic_command": ["canonical"],
        "config": {"width": 20.0},
        "claim_sources": [],
    }
    mutated = dict(semantic)
    mutated["semantic_command"] = ["mutated"]
    mutated["claims"] = []
    monkeypatch.setattr(probes, "capture_clean_source", lambda path: source)
    monkeypatch.setattr(probes, "_now_utc", lambda: datetime.datetime(2026, 8, 29, tzinfo=datetime.timezone.utc))
    monkeypatch.setattr(probes, "benchmark_semantic_input", lambda commit: semantic)
    monkeypatch.setattr(probes, "input_identity_sha256", lambda **kwargs: "c" * 64)

    @contextlib.contextmanager
    def stage(root: pathlib.Path, logical_name: str) -> Any:
        path = root / f".{logical_name}.stage-owned"
        path.mkdir(parents=True)
        try:
            yield path
        finally:
            assert not tuple(path.iterdir())
            path.rmdir()

    monkeypatch.setattr(probes, "publication_stage", stage)
    monkeypatch.setattr(probes, "run_figure6_case", lambda path, *, source_commit: (mutated, ("python", "figure6")))
    for name in ("_payload_bytes", "build_envelope", "write_envelope", "validate_benchmark_claim_artifact", "publish_stage"):
        monkeypatch.setattr(probes, name, lambda *args, _name=name, **kwargs: pytest.fail(f"{_name} called after semantic mutation"))

    with pytest.raises(probes.InvalidMeasurementClaimPayloadError, match="semantic input differs"):
        probes._produce_benchmark(tmp_path, results)

    assert not tuple(results.iterdir())


@pytest.mark.parametrize(
    "failure",
    [DirtyMeasurementTreeError("changed outside stage"), MeasurementInputChangedError("HEAD changed")],
)
def test_benchmark_finish_change_leaves_no_publication(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
    failure: Exception,
) -> None:
    probes = _probes()
    results = tmp_path / "benchmarks" / "measurement_claim_results"
    source = SourceSnapshot.build(repository=tmp_path, commit="a" * 40, pixi_lock_sha256="b" * 64)
    semantic = {"source_commit": "a" * 40}
    payload = {"source_commit": "a" * 40}
    stage_path: pathlib.Path | None = None
    monkeypatch.setattr(probes, "capture_clean_source", lambda path: source)
    monkeypatch.setattr(probes, "_now_utc", lambda: datetime.datetime(2026, 8, 29, tzinfo=datetime.timezone.utc))
    monkeypatch.setattr(probes, "benchmark_semantic_input", lambda commit: semantic)
    monkeypatch.setattr(probes, "benchmark_payload_semantic_input", lambda candidate: semantic)
    monkeypatch.setattr(probes, "input_identity_sha256", lambda **kwargs: "c" * 64)

    @contextlib.contextmanager
    def stage(root: pathlib.Path, logical_name: str) -> Any:
        nonlocal stage_path
        stage_path = root / f".{logical_name}.stage-owned"
        stage_path.mkdir(parents=True)
        try:
            yield stage_path
        finally:
            for child in stage_path.iterdir():
                child.unlink()
            stage_path.rmdir()

    monkeypatch.setattr(probes, "publication_stage", stage)

    def run_case(path: pathlib.Path, *, source_commit: str) -> tuple[dict[str, str], tuple[str, ...]]:
        del source_commit
        (path / "figure6.md").write_bytes(b"markdown")
        (path / "figure6.json").write_bytes(b"{}\n")
        return payload, ("python", "figure6")

    monkeypatch.setattr(probes, "run_figure6_case", run_case)
    monkeypatch.setattr(probes, "build_envelope", lambda **kwargs: {"input_identity": {"sha256": "c" * 64}})
    monkeypatch.setattr(probes, "write_envelope", lambda path, envelope: None)
    monkeypatch.setattr(
        probes,
        "validate_benchmark_claim_artifact",
        lambda path: ({}, object(), object(), "stage") if path == stage_path else pytest.fail("final artifact validated after failed publish"),
    )

    def fail_publish(*, source: object, stage: pathlib.Path, final: pathlib.Path) -> None:
        del source, stage, final
        raise failure

    monkeypatch.setattr(probes, "publish_stage", fail_publish)

    with pytest.raises(type(failure), match=str(failure)):
        probes._produce_benchmark(tmp_path, results)

    assert not tuple(results.iterdir())


def test_validate_cli_prints_authenticated_canonical_path(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    probes = _probes()
    artifact = pathlib.Path("benchmarks/measurement_claim_results/example")
    canonical = pathlib.PurePosixPath("benchmarks/measurement_claim_results/canonical")
    monkeypatch.setattr(probes, "validate_claim_artifact", lambda path: ({}, object(), object(), canonical))

    assert probes.main(["validate", str(artifact)]) == 0
    assert capsys.readouterr().out == f"{canonical}\n"


def test_validate_ledger_cli_delegates_the_complete_consumer_boundary(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    probes = _probes()
    ledger = pathlib.Path("docs/measurement_claims.md")
    generator = pathlib.Path("benchmarks/measurement_claim_results/generator")
    benchmark = pathlib.Path("benchmarks/measurement_claim_results/benchmark")
    calls: list[tuple[pathlib.Path, tuple[pathlib.Path, pathlib.Path]]] = []
    monkeypatch.setattr(
        probes,
        "validate_ledger_evidence",
        lambda candidate, artifacts: calls.append((candidate, artifacts)),
    )
    monkeypatch.setattr(
        probes,
        "validate_claim_artifact",
        lambda path: pytest.fail(f"CLI bypassed ledger consumer boundary: {path}"),
    )

    assert probes.main(["validate-ledger", "--ledger", str(ledger), str(generator), str(benchmark)]) == 0

    assert calls == [(ledger, (generator, benchmark))]
    assert capsys.readouterr().out == f"{ledger}\n"


def test_validate_ledger_cli_rejects_non_pair() -> None:
    probes = _probes()
    with pytest.raises(SystemExit):
        probes.main(["validate-ledger", "artifact-only"])


def test_run_generator_does_not_advertise_an_unstampable_results_override(capsys: pytest.CaptureFixture[str]) -> None:
    probes = _probes()

    with pytest.raises(SystemExit):
        probes.main(["run-generator", "--all", "--results", "elsewhere"])

    assert "unrecognized arguments: --results elsewhere" in capsys.readouterr().err
