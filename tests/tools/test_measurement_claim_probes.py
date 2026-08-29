from __future__ import annotations

import contextlib
import datetime
import importlib
import importlib.util
import pathlib
from typing import Any

import pytest

from tools.measurement_artifact import DirtyMeasurementTreeError
from tools.measurement_artifact import SourceSnapshot


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


def test_validate_cli_prints_authenticated_canonical_path(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    probes = _probes()
    artifact = pathlib.Path("benchmarks/measurement_claim_results/example")
    canonical = pathlib.PurePosixPath("benchmarks/measurement_claim_results/canonical")
    monkeypatch.setattr(probes, "validate_claim_artifact", lambda path: ({}, object(), object(), canonical))

    assert probes.main(["validate", str(artifact)]) == 0
    assert capsys.readouterr().out == f"{canonical}\n"
