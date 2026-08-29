"""Run one fixed Figure-6 child and interpret its raw outputs."""

from __future__ import annotations

import pathlib
import subprocess
import sys
from typing import Tuple

from tools import measurement_claim_json
from tools.measurement_artifact import GitObjectId
from tools.measurement_claim_benchmark_identity import MARKDOWN_NAME
from tools.measurement_claim_benchmark_identity import RAW_JSON_NAME
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SEMANTIC_COMMAND
from tools.measurement_claim_benchmark_validation import compose_benchmark_payload
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_errors import MeasurementClaimChildError


def _read_regular(path: pathlib.Path, field: str) -> bytes:
    if not path.is_file() or path.is_symlink():
        raise InvalidMeasurementClaimPayloadError(f"{field}: must be a regular non-symlink file")
    return path.read_bytes()


def _read_raw_outputs(stage: pathlib.Path) -> Tuple[object, bytes]:
    raw_json = _read_regular(stage / RAW_JSON_NAME, RAW_JSON_NAME)
    markdown = _read_regular(stage / MARKDOWN_NAME, MARKDOWN_NAME)
    decoded = measurement_claim_json.decode_strict(raw_json, RAW_JSON_NAME, InvalidMeasurementClaimPayloadError)
    return decoded, markdown


def _executed_argv(stage: pathlib.Path) -> Tuple[str, ...]:
    semantic = FIGURE6_SEMANTIC_COMMAND
    return (
        sys.executable,
        semantic[0],
        semantic[1],
        semantic[2],
        "--out",
        str(stage),
        *semantic[3:],
    )


def run_figure6_case(stage: pathlib.Path, *, source_commit: GitObjectId) -> Tuple[BenchmarkClaimPayload, Tuple[str, ...]]:
    """Run the fixed benchmark child and derive four claim records."""
    executed_argv = _executed_argv(stage)
    try:
        completed = subprocess.run(executed_argv, check=False)
    except OSError as exc:
        raise MeasurementClaimChildError(f"Figure-6 child could not start: argv={executed_argv!r}; os_error={exc}") from exc
    if completed.returncode < 0:
        raise MeasurementClaimChildError(f"Figure-6 child terminated by signal: argv={executed_argv!r}; returncode={completed.returncode}; signal={-completed.returncode}")
    if completed.returncode > 0:
        raise MeasurementClaimChildError(f"Figure-6 child failed: argv={executed_argv!r}; returncode={completed.returncode}")
    raw, markdown = _read_raw_outputs(stage)
    payload = compose_benchmark_payload(source_commit=source_commit, figure6_payload=raw, figure6_markdown=markdown)
    return payload, executed_argv
