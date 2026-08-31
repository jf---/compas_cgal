from __future__ import annotations

import importlib
import pathlib
import subprocess
from types import SimpleNamespace
from typing import Any

import pytest


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
PINNED_BENCHMARK_LEDGER_COMMIT = "9fd38c674ad321617c6cf982c0fdd233959adfed"


def _module() -> Any:
    return importlib.import_module("tools.measurement_claim_ledger")


def _rows() -> tuple[dict[str, str], ...]:
    dispositions = ["corrected"] * 10 + ["not-a-claim", "deleted", "corrected", "not-a-claim"]
    return tuple(
        {
            "ordinal": f"{index:03d}",
            "claim_id": f"MC-{index:03d}",
            "extracted_location": "source.py:1",
            "stable_anchor": "anchor",
            "anchor_match": "1/1",
            "disposition": dispositions[index - 1],
            "evidence": f"evidence-{index}",
        }
        for index in range(1, 15)
    )


def _generator_payload() -> dict[str, object]:
    return {
        "source_commit": "a" * 40,
        "source_correction_commit": "b" * 40,
        "claims": [{"claim_id": f"MC-{index:03d}", "disposition": "corrected"} for index in range(1, 11)],
    }


def _benchmark_payload() -> dict[str, object]:
    dispositions = ("not-a-claim", "deleted", "corrected", "not-a-claim")
    return {
        "source_commit": "a" * 40,
        "claims": [{"claim_id": f"MC-{index:03d}", "disposition": dispositions[index - 11]} for index in range(11, 15)],
    }


def _artifact_paths(tmp_path: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    root = tmp_path / "repository" / "benchmarks" / "measurement_claim_results"
    return root / "generator", root / "benchmark"


@pytest.fixture(autouse=True)
def _trust_benchmark_history_for_joint_unit_tests(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(_module(), "validate_benchmark_history", lambda repository, payload: None)


def _install_joint_authentication(
    monkeypatch: pytest.MonkeyPatch,
    module: Any,
    rows: tuple[dict[str, str], ...],
    rendered: list[str],
) -> None:
    monkeypatch.setattr(module, "validate_ledger_structure", lambda ledger: rows)
    monkeypatch.setattr(
        module,
        "validate_claim_artifact",
        lambda path: (_generator_payload(), SimpleNamespace(commit="a" * 40), "generator-started", "generator-artifact"),
    )
    monkeypatch.setattr(module, "validate_task6_source_correction", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_task6_source_lineage", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        module,
        "validate_benchmark_claim_artifact",
        lambda path: (_benchmark_payload(), SimpleNamespace(commit="a" * 40), "benchmark-started", "benchmark-artifact"),
    )
    monkeypatch.setattr(
        module,
        "render_generator_ledger_evidence",
        lambda *args, **kwargs: rendered.append("generator") or {f"MC-{index:03d}": f"evidence-{index}" for index in range(1, 11)},
    )
    monkeypatch.setattr(
        module,
        "render_benchmark_ledger_evidence",
        lambda *args, **kwargs: rendered.append("benchmark") or {f"MC-{index:03d}": f"evidence-{index}" for index in range(11, 15)},
    )


@pytest.mark.parametrize(
    "artifacts",
    [[], (), (pathlib.Path("generator"),), (pathlib.Path("generator"), pathlib.Path("benchmark"), pathlib.Path("extra"))],
)
def test_joint_ledger_rejects_non_exact_tuple_before_rows_or_authentication(monkeypatch: pytest.MonkeyPatch, artifacts: object) -> None:
    module = _module()
    monkeypatch.setattr(module, "validate_ledger_structure", lambda ledger: pytest.fail(f"compared rows: {ledger}"))
    monkeypatch.setattr(module, "validate_claim_artifact", lambda path: pytest.fail(f"authenticated generator: {path}"))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="two-member tuple"):
        module.validate_ledger_evidence(pathlib.Path("ledger.md"), artifacts)


def test_generator_failure_prevents_benchmark_authentication(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    monkeypatch.setattr(module, "validate_ledger_structure", lambda ledger: _rows())
    monkeypatch.setattr(module, "validate_claim_artifact", lambda path: (_ for _ in ()).throw(module.InvalidMeasurementClaimLedgerError("generator failed")))
    monkeypatch.setattr(module, "validate_benchmark_claim_artifact", lambda path: pytest.fail(f"benchmark authenticated: {path}"))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="generator failed"):
        module.validate_ledger_evidence(pathlib.Path("ledger.md"), (pathlib.Path("generator"), pathlib.Path("benchmark")))


def test_benchmark_failure_prevents_row_comparison(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> None:
    module = _module()
    generator = tmp_path / "benchmarks" / "measurement_claim_results" / "generator"
    benchmark = tmp_path / "benchmarks" / "measurement_claim_results" / "benchmark"
    monkeypatch.setattr(module, "validate_ledger_structure", lambda ledger: _rows())
    monkeypatch.setattr(module, "validate_claim_artifact", lambda path: (_generator_payload(), SimpleNamespace(commit="a" * 40), object(), object()))
    monkeypatch.setattr(module, "validate_task6_source_correction", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_task6_source_lineage", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_benchmark_claim_artifact", lambda path: (_ for _ in ()).throw(module.InvalidMeasurementClaimLedgerError("benchmark failed")))
    monkeypatch.setattr(module, "_compare_rows", lambda *args: pytest.fail("rows compared"))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="benchmark failed"):
        module.validate_ledger_evidence(pathlib.Path("ledger.md"), (generator, benchmark))


def test_joint_ledger_requires_exact_disjoint_fourteen_claim_union(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> None:
    module = _module()
    generator = tmp_path / "benchmarks" / "measurement_claim_results" / "generator"
    benchmark = tmp_path / "benchmarks" / "measurement_claim_results" / "benchmark"
    payload = _benchmark_payload()
    payload["claims"][0]["claim_id"] = "MC-010"
    monkeypatch.setattr(module, "validate_ledger_structure", lambda ledger: _rows())
    monkeypatch.setattr(module, "validate_claim_artifact", lambda path: (_generator_payload(), SimpleNamespace(commit="a" * 40), object(), object()))
    monkeypatch.setattr(module, "validate_task6_source_correction", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_task6_source_lineage", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_benchmark_claim_artifact", lambda path: (payload, SimpleNamespace(commit="a" * 40), object(), object()))
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="MC-001.*MC-014|union"):
        module.validate_ledger_evidence(pathlib.Path("ledger.md"), (generator, benchmark))


@pytest.mark.parametrize(
    ("first", "second"),
    [("benchmark", "generator"), ("generator", "generator"), ("benchmark", "benchmark")],
)
def test_joint_ledger_rejects_reversed_or_wrong_family_pair_before_comparison(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
    first: str,
    second: str,
) -> None:
    module = _module()
    generator, benchmark = _artifact_paths(tmp_path)
    paths = {"generator": generator, "benchmark": benchmark}
    monkeypatch.setattr(module, "validate_ledger_structure", lambda ledger: _rows())

    def validate_generator(path: pathlib.Path) -> tuple[object, ...]:
        if path.name != "generator":
            raise module.InvalidMeasurementClaimLedgerError("wrong generator family")
        return _generator_payload(), SimpleNamespace(commit="a" * 40), object(), object()

    def validate_benchmark(path: pathlib.Path) -> tuple[object, ...]:
        if path.name != "benchmark":
            raise module.InvalidMeasurementClaimLedgerError("wrong benchmark family")
        return _benchmark_payload(), SimpleNamespace(commit="a" * 40), object(), object()

    monkeypatch.setattr(module, "validate_claim_artifact", validate_generator)
    monkeypatch.setattr(module, "validate_task6_source_correction", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_task6_source_lineage", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_benchmark_claim_artifact", validate_benchmark)
    monkeypatch.setattr(module, "_compare_rows", lambda *args: pytest.fail("rows compared for a wrong-family pair"))

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="wrong .* family"):
        module.validate_ledger_evidence(pathlib.Path("ledger.md"), (paths[first], paths[second]))


def test_joint_ledger_renders_both_families_and_accepts_fourteen_byte_equal_rows(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    rows = _rows()
    rendered: list[str] = []
    _install_joint_authentication(monkeypatch, module, rows, rendered)
    ledger = tmp_path / "ledger.md"
    ledger.write_text(f"{module.COMPLETE_STATUS}\n", encoding="utf-8")

    assert module.validate_ledger_evidence(ledger, _artifact_paths(tmp_path)) is None
    assert rendered == ["generator", "benchmark"]


def test_joint_ledger_rejects_truncated_rows_before_zip_comparison(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    rows = _rows()[:-1]
    _install_joint_authentication(monkeypatch, module, rows, [])
    ledger = tmp_path / "ledger.md"
    ledger.write_text(f"{module.COMPLETE_STATUS}\n", encoding="utf-8")
    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="exactly 14 rows"):
        module.validate_ledger_evidence(ledger, _artifact_paths(tmp_path))


def test_joint_ledger_invokes_benchmark_history_before_rendering(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
) -> None:
    module = _module()
    events: list[str] = []
    _install_joint_authentication(monkeypatch, module, _rows(), events)
    monkeypatch.setattr(module, "validate_benchmark_history", lambda repository, payload: events.append("history"))
    ledger = tmp_path / "ledger.md"
    ledger.write_text(f"{module.COMPLETE_STATUS}\n", encoding="utf-8")
    assert module.validate_ledger_evidence(ledger, _artifact_paths(tmp_path)) is None
    assert events == ["history", "generator", "benchmark"]


@pytest.mark.parametrize("damage", ["disposition", "evidence", "pending"])
def test_joint_ledger_rejects_row_damage(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: pathlib.Path,
    damage: str,
) -> None:
    module = _module()
    rows = tuple(dict(row) for row in _rows())
    if damage == "disposition":
        rows[10]["disposition"] = "corrected"
    elif damage == "evidence":
        rows[13]["evidence"] = "manually transcribed"
    elif damage == "pending":
        rows[11]["disposition"] = "pending"
        rows[11]["evidence"] = "—"
    _install_joint_authentication(monkeypatch, module, rows, [])
    ledger = tmp_path / "ledger.md"
    ledger.write_text(f"{module.COMPLETE_STATUS}\n", encoding="utf-8")

    with pytest.raises(module.InvalidMeasurementClaimLedgerError, match="disposition|evidence|pending"):
        module.validate_ledger_evidence(ledger, _artifact_paths(tmp_path))


def test_repository_ledger_accepts_exact_authenticated_generator_and_benchmark_artifacts() -> None:
    module = _module()
    repository = pathlib.Path(__file__).resolve().parents[2]
    artifact_root = repository / "benchmarks" / "measurement_claim_results"
    generators = tuple(sorted(artifact_root.glob("*-generator-*")))
    benchmarks = tuple(sorted(artifact_root.glob("*-benchmark-*")))
    assert len(generators) == 1
    assert len(benchmarks) == 1

    assert (
        module.validate_ledger_evidence(
            repository / "docs" / "measurement_claims.md",
            (generators[0], benchmarks[0]),
        )
        is None
    )


def test_pinned_task7_ledger_blob_accepts_the_committed_benchmark_artifact(tmp_path: pathlib.Path) -> None:
    module = _module()
    raw = subprocess.run(
        [
            "git",
            "--no-replace-objects",
            "show",
            f"{PINNED_BENCHMARK_LEDGER_COMMIT}:docs/measurement_claims.md",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    ledger = tmp_path / "measurement_claims.md"
    ledger.write_bytes(raw)
    artifact_root = PROJECT_ROOT / "benchmarks" / "measurement_claim_results"
    generators = tuple(sorted(path for path in artifact_root.glob("*-generator-*") if path.is_dir()))
    benchmarks = tuple(sorted(path for path in artifact_root.glob("*-benchmark-*") if path.is_dir()))
    assert len(generators) == 1
    assert len(benchmarks) == 1
    assert module.validate_ledger_evidence(ledger, (generators[0], benchmarks[0])) is None
