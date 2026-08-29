"""Filesystem and envelope authentication for generator-v2 and benchmark-v1 claim artifacts."""

from __future__ import annotations

import datetime
import pathlib
import re
from typing import Tuple
from typing import cast

from tools import measurement_artifact
from tools import measurement_claim_json
from tools.measurement_artifact import ValidatedEnvelope
from tools.measurement_claim_benchmark_identity import ARTIFACT_KIND as BENCHMARK_ARTIFACT_KIND
from tools.measurement_claim_benchmark_identity import INPUT_VERSION as BENCHMARK_INPUT_VERSION
from tools.measurement_claim_benchmark_identity import MARKDOWN_NAME
from tools.measurement_claim_benchmark_identity import PAYLOAD_NAME as BENCHMARK_PAYLOAD_NAME
from tools.measurement_claim_benchmark_identity import RAW_JSON_NAME
from tools.measurement_claim_benchmark_identity import RESULT_VERSION as BENCHMARK_RESULT_VERSION
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_semantic_input import FIGURE6_SEMANTIC_COMMAND
from tools.measurement_claim_benchmark_semantic_input import benchmark_payload_semantic_input
from tools.measurement_claim_benchmark_validation import validate_benchmark_payload
from tools.measurement_claim_case_validation import fail_payload
from tools.measurement_claim_case_validation import validate_array
from tools.measurement_claim_case_validation import validate_literal
from tools.measurement_claim_case_validation import validate_object
from tools.measurement_claim_case_validation import validate_same
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_identity import ARTIFACT_KIND
from tools.measurement_claim_identity import INPUT_VERSION
from tools.measurement_claim_identity import PAYLOAD_NAME
from tools.measurement_claim_identity import RESULT_VERSION
from tools.measurement_claim_payload import validate_generator_payload
from tools.measurement_claim_schema import GeneratorClaimPayload
from tools.measurement_claim_schema import ValidatedArtifactDirectory
from tools.measurement_claim_schema import ValidatedArtifactStartedUtc
from tools.measurement_claim_semantic_input import generator_semantic_input

_FINAL_NAME = re.compile(r"\d{4}-\d{2}-\d{2}-[0-9a-f]{12}-generator-[0-9a-f]{12}\Z")
_BENCHMARK_FINAL_NAME = re.compile(r"\d{4}-\d{2}-\d{2}-[0-9a-f]{12}-benchmark-[0-9a-f]{12}\Z")


def _decode(data: bytes, field: str) -> object:
    return measurement_claim_json.decode_strict(data, field, InvalidMeasurementClaimPayloadError)


def _read(path: pathlib.Path, field: str) -> bytes:
    if not path.is_file() or path.is_symlink():
        raise InvalidMeasurementClaimPayloadError(f"{field}: must be a regular non-symlink file")
    return path.read_bytes()


def _logical_name(actual_name: str) -> str:
    if actual_name.startswith("."):
        marker = ".stage-"
        if marker not in actual_name[1:]:
            raise InvalidMeasurementClaimPayloadError("artifact name is not an owned hidden stage")
        logical, suffix = actual_name[1:].split(marker, 1)
        if not suffix or marker in suffix or _FINAL_NAME.fullmatch(logical) is None:
            raise InvalidMeasurementClaimPayloadError("artifact stage wrapper is malformed")
        return logical
    if _FINAL_NAME.fullmatch(actual_name) is None:
        raise InvalidMeasurementClaimPayloadError("artifact final name is malformed")
    return actual_name


def _benchmark_logical_name(actual_name: str) -> str:
    if actual_name.startswith("."):
        marker = ".stage-"
        if marker not in actual_name[1:]:
            raise InvalidMeasurementClaimPayloadError("benchmark artifact name is not an owned hidden stage")
        logical, suffix = actual_name[1:].split(marker, 1)
        if not suffix or marker in suffix or _BENCHMARK_FINAL_NAME.fullmatch(logical) is None:
            raise InvalidMeasurementClaimPayloadError("benchmark artifact stage wrapper is malformed")
        return logical
    if _BENCHMARK_FINAL_NAME.fullmatch(actual_name) is None:
        raise InvalidMeasurementClaimPayloadError("benchmark artifact final name is malformed")
    return actual_name


def validate_claim_artifact(
    result: pathlib.Path,
) -> Tuple[GeneratorClaimPayload, ValidatedEnvelope, ValidatedArtifactStartedUtc, ValidatedArtifactDirectory]:
    """Authenticate one canonical final or exact Task-4-owned hidden stage."""
    if result.is_symlink():
        raise InvalidMeasurementClaimPayloadError("artifact directory must not be a symlink")
    lexical = result.absolute()
    if lexical.parent.name != "measurement_claim_results" or lexical.parent.parent.name != "benchmarks":
        raise InvalidMeasurementClaimPayloadError("artifact must be directly below benchmarks/measurement_claim_results")
    logical_name = _logical_name(lexical.name)
    repository = lexical.parent.parent.parent
    stamp_path, payload_path = lexical / measurement_artifact.STAMP_NAME, lexical / PAYLOAD_NAME
    stamp_before = _read(stamp_path, "stamp")
    payload_before = _read(payload_path, PAYLOAD_NAME)
    envelope = measurement_artifact.validate_envelope(
        lexical,
        logical_name=logical_name,
        artifact_kind=measurement_artifact.ArtifactKind(ARTIFACT_KIND),
        repository=repository,
    )
    stamp_after = _read(stamp_path, "stamp")
    payload_after = _read(payload_path, PAYLOAD_NAME)
    if stamp_before != stamp_after or payload_before != payload_after:
        raise InvalidMeasurementClaimPayloadError("artifact bytes changed during validation")
    stamp = validate_object(_decode(stamp_before, "stamp"), measurement_artifact.ENVELOPE_KEYS, "stamp")
    payload = validate_generator_payload(_decode(payload_before, PAYLOAD_NAME))
    validate_literal(stamp["artifact_kind"], ARTIFACT_KIND, "stamp.artifact_kind")
    argv = validate_array(stamp["argv"], "stamp.argv")
    if len(argv) != 5 or type(argv[0]) is not str or not argv[0]:
        fail_payload("stamp.argv", "must contain one executable and the canonical command")
    validate_same(argv[1:], ["-m", "tools.measurement_claim_probes", "run-generator", "--all"], "stamp.argv[1:]")
    input_identity = validate_object(stamp["input_identity"], ("version", "payload", "sha256"), "stamp.input_identity")
    result_identity = validate_object(stamp["result_identity"], ("version", "payloads", "sha256"), "stamp.result_identity")
    validate_literal(input_identity["version"], INPUT_VERSION, "stamp.input_identity.version")
    validate_literal(result_identity["version"], RESULT_VERSION, "stamp.result_identity.version")
    payloads = validate_object(result_identity["payloads"], (PAYLOAD_NAME,), "stamp.result_identity.payloads")
    validate_object(payloads[PAYLOAD_NAME], ("sha256",), f"stamp.result_identity.payloads.{PAYLOAD_NAME}")
    validate_same(input_identity["payload"], generator_semantic_input(payload), "stamp.input_identity.payload")
    validate_literal(payload["source_commit"], str(envelope.commit), "source_commit")
    started_value = stamp["started"]
    if type(started_value) is not str:
        fail_payload("stamp.started", "must be a canonical UTC timestamp")
    try:
        started = datetime.datetime.fromisoformat(cast(str, started_value))
    except ValueError as exc:
        raise InvalidMeasurementClaimPayloadError("stamp.started: not ISO-8601") from exc
    if started.utcoffset() != datetime.timedelta(0) or started.isoformat(timespec="microseconds") != started_value:
        fail_payload("stamp.started", "must be canonical timezone-aware UTC with microseconds")
    expected_name = f"{started.date().isoformat()}-{str(envelope.commit)[:12]}-generator-{str(envelope.input_sha256)[:12]}"
    validate_literal(logical_name, expected_name, "artifact logical name")
    canonical = pathlib.PurePosixPath("benchmarks", "measurement_claim_results", expected_name)
    return payload, envelope, ValidatedArtifactStartedUtc(started), ValidatedArtifactDirectory(canonical)


def validate_benchmark_claim_artifact(
    result: pathlib.Path,
) -> Tuple[BenchmarkClaimPayload, ValidatedEnvelope, ValidatedArtifactStartedUtc, ValidatedArtifactDirectory]:
    """Authenticate one explicit benchmark-family artifact."""
    if result.is_symlink():
        raise InvalidMeasurementClaimPayloadError("benchmark artifact directory must not be a symlink")
    lexical = result.absolute()
    if lexical.parent.name != "measurement_claim_results" or lexical.parent.parent.name != "benchmarks":
        raise InvalidMeasurementClaimPayloadError("benchmark artifact must be directly below benchmarks/measurement_claim_results")
    logical_name = _benchmark_logical_name(lexical.name)
    repository = lexical.parent.parent.parent
    names = (measurement_artifact.STAMP_NAME, MARKDOWN_NAME, RAW_JSON_NAME, BENCHMARK_PAYLOAD_NAME)
    before = {name: _read(lexical / name, name) for name in names}
    envelope = measurement_artifact.validate_envelope(
        lexical,
        logical_name=logical_name,
        artifact_kind=measurement_artifact.ArtifactKind(BENCHMARK_ARTIFACT_KIND),
        repository=repository,
    )
    after = {name: _read(lexical / name, name) for name in names}
    if before != after:
        raise InvalidMeasurementClaimPayloadError("benchmark artifact bytes changed during validation")
    stamp = validate_object(_decode(before[measurement_artifact.STAMP_NAME], "stamp"), measurement_artifact.ENVELOPE_KEYS, "stamp")
    raw = _decode(before[RAW_JSON_NAME], RAW_JSON_NAME)
    payload = validate_benchmark_payload(
        _decode(before[BENCHMARK_PAYLOAD_NAME], BENCHMARK_PAYLOAD_NAME),
        figure6_payload=raw,
        figure6_markdown=before[MARKDOWN_NAME],
    )
    argv = validate_array(stamp["argv"], "stamp.argv")
    if len(argv) != len(FIGURE6_SEMANTIC_COMMAND) + 3 or type(argv[0]) is not str or not argv[0]:
        fail_payload("stamp.argv", "must contain one executable, one --out pair, and the canonical command")
    validate_same(argv[1:4], list(FIGURE6_SEMANTIC_COMMAND[:3]), "stamp.argv[1:4]")
    validate_literal(argv[4], "--out", "stamp.argv[4]")
    if type(argv[5]) is not str or not argv[5]:
        fail_payload("stamp.argv[5]", "must be the executed stage path")
    validate_same(argv[6:], list(FIGURE6_SEMANTIC_COMMAND[3:]), "stamp.argv[6:]")
    input_identity = validate_object(stamp["input_identity"], ("version", "payload", "sha256"), "stamp.input_identity")
    result_identity = validate_object(stamp["result_identity"], ("version", "payloads", "sha256"), "stamp.result_identity")
    validate_literal(input_identity["version"], BENCHMARK_INPUT_VERSION, "stamp.input_identity.version")
    validate_literal(result_identity["version"], BENCHMARK_RESULT_VERSION, "stamp.result_identity.version")
    payloads = validate_object(
        result_identity["payloads"],
        (MARKDOWN_NAME, RAW_JSON_NAME, BENCHMARK_PAYLOAD_NAME),
        "stamp.result_identity.payloads",
    )
    for name in (MARKDOWN_NAME, RAW_JSON_NAME, BENCHMARK_PAYLOAD_NAME):
        validate_object(payloads[name], ("sha256",), f"stamp.result_identity.payloads.{name}")
    validate_same(input_identity["payload"], benchmark_payload_semantic_input(payload), "stamp.input_identity.payload")
    validate_literal(payload["source_commit"], str(envelope.commit), "benchmark-claims.json.source_commit")
    started_value = stamp["started"]
    if type(started_value) is not str:
        fail_payload("stamp.started", "must be a canonical UTC timestamp")
    try:
        started = datetime.datetime.fromisoformat(cast(str, started_value))
    except ValueError as exc:
        raise InvalidMeasurementClaimPayloadError("stamp.started: not ISO-8601") from exc
    if started.utcoffset() != datetime.timedelta(0) or started.isoformat(timespec="microseconds") != started_value:
        fail_payload("stamp.started", "must be canonical timezone-aware UTC with microseconds")
    expected_name = f"{started.date().isoformat()}-{str(envelope.commit)[:12]}-benchmark-{str(envelope.input_sha256)[:12]}"
    validate_literal(logical_name, expected_name, "benchmark artifact logical name")
    canonical = pathlib.PurePosixPath("benchmarks", "measurement_claim_results", expected_name)
    return payload, envelope, ValidatedArtifactStartedUtc(started), ValidatedArtifactDirectory(canonical)
