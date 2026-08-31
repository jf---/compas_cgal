"""Authenticated ledger rendering for the benchmark claim family."""

from __future__ import annotations

import datetime
import pathlib
from typing import Dict
from typing import Mapping
from typing import cast

from tools.measurement_artifact import ValidatedEnvelope
from tools.measurement_claim_benchmark_identity import PAYLOAD_NAME
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_json import canonical_text
from tools.measurement_claim_markdown import require_markdown_safe
from tools.measurement_claim_schema import ValidatedArtifactDirectory
from tools.measurement_claim_schema import ValidatedArtifactStartedUtc


def _artifact_text(
    payload: BenchmarkClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> str:
    if type(started) is not datetime.datetime or started.tzinfo is None or started.utcoffset() != datetime.timedelta(0):
        raise InvalidMeasurementClaimPayloadError("benchmark artifact started must be one timezone-aware UTC datetime")
    if started > envelope.finished:
        raise InvalidMeasurementClaimPayloadError("benchmark artifact started must not follow the authenticated finish")
    if payload["source_commit"] != envelope.commit:
        raise InvalidMeasurementClaimPayloadError("benchmark payload source commit must equal the authenticated envelope commit")
    expected = pathlib.PurePosixPath(
        "benchmarks",
        "measurement_claim_results",
        f"{started.date().isoformat()}-{str(envelope.commit)[:12]}-benchmark-{str(envelope.input_sha256)[:12]}",
    )
    if type(artifact_directory) is not pathlib.PurePosixPath or artifact_directory != expected:
        raise InvalidMeasurementClaimPayloadError(f"benchmark artifact directory must equal the authenticated repository-relative path: {expected}")
    relative = require_markdown_safe(artifact_directory.as_posix(), "benchmark artifact directory", InvalidMeasurementClaimPayloadError)
    return f"{relative}/{PAYLOAD_NAME}@sha256:{envelope.payload_sha256[PAYLOAD_NAME]}"


def render_benchmark_ledger_evidence(
    payload: BenchmarkClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> Mapping[str, str]:
    """Render four ledger cells from values rebound to one envelope.

    Args:
        payload: Validated benchmark claim payload.
        envelope: Authenticated artifact envelope.
        started: Validated artifact start timestamp.
        artifact_directory: Validated repository-relative artifact path.

    Returns:
        Claim IDs mapped to Markdown-safe evidence cells.

    Raises:
        InvalidMeasurementClaimPayloadError: A caller value, per-claim source,
            JSON field, or rendered cell is not bound to the envelope.
    """
    artifact = _artifact_text(
        payload,
        envelope,
        started=started,
        artifact_directory=artifact_directory,
    )
    semantic_command = canonical_text(payload["semantic_command"], "benchmark semantic command", InvalidMeasurementClaimPayloadError)
    config = canonical_text(payload["config"], "benchmark config", InvalidMeasurementClaimPayloadError)
    values: Dict[str, str] = {}
    for claim in payload["claims"]:
        record = cast(Dict[str, object], claim)
        claim_id = cast(str, record["claim_id"])
        if record["source_commit"] != envelope.commit:
            raise InvalidMeasurementClaimPayloadError(f"{claim_id}: claim source commit must equal the authenticated envelope commit")
        cell = (
            f"artifact={artifact}; started={started.isoformat(timespec='microseconds')}; claim={claim_id}; "
            f"input={envelope.input_sha256}; result={envelope.result_sha256}; source={record['source_commit']}; "
            f"history={record['history_commit']}; disposition={record['disposition']}; reason={record['reason']}; "
            f"semantic_command={semantic_command}; config={config}; "
            f"missing_inputs={canonical_text(record['missing_inputs'], f'{claim_id} missing_inputs', InvalidMeasurementClaimPayloadError)}"
        )
        selected = record.get("selected_values")
        if selected is not None:
            cell += f"; selected_values={canonical_text(selected, f'{claim_id} selected_values', InvalidMeasurementClaimPayloadError)}"
        values[claim_id] = require_markdown_safe(cell, claim_id, InvalidMeasurementClaimPayloadError)
    return values
