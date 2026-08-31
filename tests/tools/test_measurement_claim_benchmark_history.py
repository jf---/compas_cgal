from __future__ import annotations

import copy
import importlib
import pathlib

import pytest

from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS_ONLY_COMMIT = "9fd38c674ad321617c6cf982c0fdd233959adfed"


def _benchmark_payload() -> BenchmarkClaimPayload:
    result = importlib.import_module("tools.measurement_claim_result")
    artifact_root = PROJECT_ROOT / "benchmarks" / "measurement_claim_results"
    artifacts = tuple(path for path in artifact_root.glob("*-benchmark-*") if path.is_dir())
    assert len(artifacts) == 1
    payload, _, _, _ = result.validate_benchmark_claim_artifact(artifacts[0])
    return payload


def test_benchmark_history_commits_exist_and_touch_each_claim_source() -> None:
    history = importlib.import_module("tools.measurement_claim_benchmark_history")
    assert history.validate_benchmark_history(PROJECT_ROOT, _benchmark_payload()) is None


@pytest.mark.parametrize("damage", ["missing", "wrong-path"])
def test_benchmark_history_rejects_unverifiable_commit_identity(damage: str) -> None:
    history = importlib.import_module("tools.measurement_claim_benchmark_history")
    payload = copy.deepcopy(_benchmark_payload())
    payload["claims"][0]["history_commit"] = "a" * 40 if damage == "missing" else DOCS_ONLY_COMMIT
    with pytest.raises(history.InvalidMeasurementClaimPayloadError, match="history|commit|source"):
        history.validate_benchmark_history(PROJECT_ROOT, payload)
