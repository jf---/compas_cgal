"""Filesystem and envelope authentication for v2 claim artifacts."""

from __future__ import annotations

import datetime
import pathlib
import re
from typing import Tuple
from typing import cast

from tools import measurement_artifact
from tools import measurement_claim_json
from tools.measurement_artifact import ValidatedEnvelope
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
